import os, sys
import gzip
from cStringIO import StringIO
from contextlib import contextmanager
from xml.etree.cElementTree import fromstring, parse, ParseError, tostring

import requests

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import natural_keys, get_pdb_residues, izip_missing


class InvalidSIFTS(RuntimeError):
    pass

ns = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"
def process_residue(residue):
    try:
        pdb = int(residue.findall('.//{}crossRefDb[@dbSource="PDB"]'.format(ns))[0].attrib["dbResNum"])
    except (IndexError, ValueError) as e:
        pdb = None

    try:
        ncbi = int(residue.findall('.//{}crossRefDb[@dbSource="NCBI"]'.format(ns))[0].attrib["dbResNum"])
    except (IndexError, ValueError) as e:
        ncbi = None

    return pdb, ncbi

def mmdb_to_pdb_resi(pdb_name, chain, resi, replace_nulls=False, job=None):
    if len(resi) == 1 and isinstance(resi[0], str) and "," in resi[0]:
        resi = resi[0].split(",")

    complete_resi = ",".join(map(str,resi))

    resi = [natural_keys(r, use_int=True) for r in sorted(resi, key=natural_keys)]

    residue_mapping = {}
    #sys.stdout = None
    with get_sifts(pdb_name, job=job) as sifts_file:
        xml = gzip.GzipFile(None, "rb", 9, sifts_file).read()
        try:
            root = fromstring(xml)
        except ParseError:
            raise InvalidSIFTS("Error parsing: "+pdb_name)

        for entityNum, entity in enumerate(root.findall(".//{}entity".format(ns))):
            try:
                dbChainId = entity.find(".//{ns}segment/{ns}listResidue/{ns}residue/{ns}crossRefDb[@dbSource='PDB']".format(ns=ns)).attrib["dbChainId"]
            except (AttributeError, KeyError):
                raise InvalidSIFTS(pdb_name)
            if  dbChainId == chain or (chain == " " and entityNum == 0):
                residues = entity.findall(".//{ns}segment/{ns}listResidue/{ns}residue".format(ns=ns))
                for residue in residues:
                    pdb, ncbi = process_residue(residue)
                    residue_mapping[tuple(natural_keys(residue.attrib["dbResNum"], use_int=True))] = pdb
                break
    del xml
    del root

    sorted_map = sorted(residue_mapping.iteritems(), key=lambda x:natural_keys(x[0]))

    if len(sorted_map) == 0:
        raise InvalidSIFTS("No mapping for {}.{}".format(pdb_name, chain))

    for i, r in enumerate(resi):
        try:
            pdb_resnum = residue_mapping[tuple(r)]
        except KeyError:
            try:
                sorted_map[-1]
            except IndexError:
                print pdb_name, chain, r, sorted_map
            if r[1] > sorted_map[-1][0][1]:
                last = sorted_map[-1][1]
                i=1
                while last is None:
                    last = sorted_map[-i][1]
                    i+=1
                yield last
            elif r[1] < sorted_map[0][0][1]:
                first = sorted_map[0][1]
                i=1
                while first is None:
                    first = sorted_map[i][1]
                    i+=1
                yield first
            else:
                raise InvalidSIFTS("?? resnum {} does not exit {}".format(r, sorted_map))
            raise StopIteration

        if pdb_resnum is not None:
            yield pdb_resnum
        else:
            lookahead = False
            for sorted_resi, sorted_pdb in sorted_map:
                if sorted_resi == tuple(r):
                    lookahead = True
                    continue
                if lookahead and sorted_pdb is not None:
                    yield sorted_pdb
                    break
            else:
                lookbehind = False
                for sorted_resi, sorted_pdb in reversed(sorted_map):
                    if sorted_resi == tuple(r):
                        lookbehind = True
                        continue
                    if lookbehind and sorted_pdb is not None:
                        yield sorted_pdb
                        break
                else:
                    raise InvalidSIFTS("{} {} resnum {} does not exit {}".format(pdb_name, chain, r, sorted_map))
    del sorted_map
    del residue_mapping

@contextmanager
def get_sifts(pdb, job=None):
    if job is not None:
        work_dir = job.fileStore.getLocalTempDir()
        prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
        in_store = IOStore.get("{}:molmimic-sifts".format(prefix))
        sifts_prefix = "{}/{}.xml.gz".format(pdb[1:3].lower(), pdb.lower())
	sifts_path = os.path.join(work_dir, os.path.basename(sifts_prefix))        
        job.log("Saving {}:molmimic-sifts :: {} to {}".format(prefix, sifts_prefix, sifts_path))
        
        try:
            in_store.read_input_file(sifts_prefix, sifts_path)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise InvalidSIFTS("Cannot open {}".format(pdb))

        with open(sifts_path) as f:
	    yield f

        os.remove(sifts_path)
    else: 
        path = os.path.join(os.environ.get("PDB_SNAPSHOT", os.path.join(data_path_prefix, "pdb")), "sifts", pdb[1:3].lower(), "{}.xml.gz".format(pdb.lower()))
        try:
            with open(path) as f:
                yield file
        except IOError as e:
            url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}.xml.gz".format(pdb.lower())
            try:
                sifts = requests.get(url)
                file = StringIO(sifts.content)
                yield file
            except requests.exceptions.RequestException:
                file.close()
                raise InvalidSIFTS("Not found: "+url+" orgig error:"+str(e))
            finally:
                file.close()
        os.remove(sifts)

def comare_to_pdb(pdb_file, resi):
    if len(resi) == 1 and isinstance(resi[0], str) and "," in resi[0]:
        resi = resi[0].split(",")

    resi = (natural_keys(r, use_int=True) for r in sorted(resi, key=natural_keys))

    for test_resi, known_resi in izip_missing(iter(resi), get_pdb_residues(pdb_file)):
        yield test_resi
