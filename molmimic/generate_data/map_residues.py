import os, sys
import gzip
from io import StringIO
from contextlib import contextmanager
from xml.etree.cElementTree import fromstring, parse, ParseError, tostring
import numpy as np
import requests
import binascii
import zlib
import json
from pyasn1.codec.ber import decoder
from itertools import groupby

import numpy as np

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import natural_keys, get_pdb_residues, izip_missing

from toil.realtimeLogger import RealtimeLogger

from joblib import Memory
memory = Memory("/tmp", verbose=0)

class InvalidSIFTS(RuntimeError):
    pass

class ChangePDB(Exception):
    def __init__(self, gi, old_pdb, old_chain, new_pdb, new_chain, new_res, num_residues):
        self.gi = gi
        self.old_pdb  = old_pdb
        self.old_chain = old_chain
        self.new_pdb = new_pdb
        self.new_chain = new_chain
        self.new_res = new_res
        self.num_residues = num_residues

ns = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"
def process_residue(residue):
    try:
        pdb = natural_keys(residue.findall('.//{}crossRefDb[@dbSource="PDB"]'.format(ns))[0].attrib["dbResNum"], use_int=True)
    except (IndexError, ValueError, AssertionError) as e:
        pdb = None

    try:
        ncbi = natural_keys(residue.findall('.//{}crossRefDb[@dbSource="NCBI"]'.format(ns))[0].attrib["dbResNum"], use_int=True)
    except (IndexError, ValueError, AssertionError) as e:
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
                    residue_mapping[tuple(natural_keys(int(residue.attrib["dbResNum"]), use_int=True))] = pdb
                break
    del xml
    del root

    sorted_map = list(sorted(residue_mapping.items(), key=lambda x: x[0]))

    if len(sorted_map) == 0:
        raise InvalidSIFTS("No mapping for {}.{}".format(pdb_name, chain))

    for i, r in enumerate(resi):
        try:
            pdb_resnum = residue_mapping[tuple(r)]
        except KeyError:
            last = sorted_map[-1][1]
            i=1
            while last is None:
                last = sorted_map[-i][1]
                i+=1

            if r[1] > last:
                yield last
            else:
                first = sorted_map[0][1]
                i=1
                while first is None:
                    first = sorted_map[i][1]
                    i+=1

                if r[1] < first:
                    yield first
                else:
                    raise InvalidSIFTS("?? resnum {} does not exit {}".format(r, sorted_map))
            #raise StopIteration
            continue

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

def get_entrez(gi, work_dir):
    store = IOStore.get("aws:us-east-1:molmimic-mmdb-mappings")
    key = "{}.json".format(gi)

    if store.exists(key):
        entrez_file = os.path.join(work_dir, key)
        try:
            store.read_input_file(key, entrez_file)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise InvalidSIFTS("Unable to donwload Entrex MMDB mapping for {}".format(gi))

        with open(entrez_file) as f:
            mapping = json.load(f)

        try:
            os.remove(entrez_file)
        except OSError:
            pass

        return mapping["mmdb_indices"], mapping["seq"], mapping["num_residues"], \
            mapping["pdb"], mapping["chain"]


    from Bio import Entrez
    Entrez.email = "edraizen"

    mmdb_indices = []
    seq = ""
    num_residues = None

    last_pdb = None
    last_chain = None

    parse_gi = False
    parse_indices = False
    parse_seq = False

    try:
        handle = Entrez.efetch(db="protein", id=str(gi), retmode="text")
    except AttributeError:
        RealtimeLogger.info("Cannot download GI: {}".format(gi))
        raise InvalidSIFTS("Cannot download GI: {}".format(gi))

    for i, line in enumerate(handle):
        if "mol" in line and not parse_gi:
            RealtimeLogger.info("MOL LINE is '{}'".format(line.strip()))
            last_pdb = line.strip().split()[-1][1:-2]
            RealtimeLogger.info("MOL LINE is '{}'".format(last_pdb))


        if "chain" in line and not parse_gi:
            try:
                last_chain = chr(int(line.strip().split()[-1][:-1]))
            except ValueError:
                last_chain = ""

        if not parse_gi and "gi" in line and line.strip().split()[-1] == str(gi):
            RealtimeLogger.info("{} {}".format(i, line))
            parse_gi = True
            continue

        if not parse_gi:
            continue

        if not parse_indices and "num enum {" in line:
            parse_indices = True
            continue

        if parse_indices and "num" in line:
            num_residues = int(line.strip().split()[-1][:-1])
            continue
        elif parse_indices and "names" in line:
            continue
        elif parse_indices and "}" in line:
            parse_indices = False
            continue
        elif parse_indices:
            if "," in line:
                if line.strip()[1:-2] == "":
                    mmdb_indices.append(None)
                else:
                    mmdb_indices.append(natural_keys(line.strip()[1:-2], use_int=True))

            else:
                if line.strip()[1:-1] == "":
                    mmdb_indices.append(None)
                else:
                    mmdb_indices.append(natural_keys(line.strip()[1:-1], use_int=True))

                parse_indices = False
            # else:
            #     mmdb_indices.append(None)
            continue

        if not parse_seq and "seq-data ncbieaa" in line:
            seq += line.rstrip().split()[-1][1:]
            parse_seq = True
        elif parse_seq and not '"' in line:
            seq += line.strip()
        elif parse_seq and '"' in line:
            seq += line.strip()[:-1]
            parse_seq = False
            break

    if '"' in seq or "}" in seq:
        seq = seq.split('"')[0]

    mappings = {"mmdb_indices":mmdb_indices, "seq":seq, "num_residues":num_residues,
        "pdb":last_pdb, "chain": last_chain}

    entrez_file = os.path.join(work_dir, key)
    with open(entrez_file, "w") as f:
        json.dump(mappings, f)

    store.write_output_file(entrez_file, key)

    return mmdb_indices, seq, num_residues, last_pdb, last_chain

def mmdb_to_pdb_resi_obsolete(gi, pdb_name, chain, resi, shrink_or_expand=False, replace_nulls=False, return_gi_len=False, work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    from Bio import Entrez
    from Bio import SeqIO
    from Bio.PDB import PDBParser, MMCIFParser, Polypeptide


    from molmimic.generate_data.util import s3_download_pdb

    if len(resi) == 1 and isinstance(resi[0], str) and "," in resi[0]:
        resi = resi[0].split(",")

    resi = map(int, resi)

    # try:
    #     pdb_file, file_format = s3_download_pdb(pdb_name, work_dir=work_dir)
    # except (SystemExit, KeyboardInterrupt):
    #     raise
    # except IOError:
    #     RealtimeLogger.info("{} {} {} {}".format(gi, pdb_name, chain, resi))
    #     raise
    # except:
    #     raise
    #     raise InvalidSIFTS
    #
    # if pdb_file.endswith(".gz"):
    #     _pdb_file = pdb_file[:-3]
    #     with gzip.open(pdb_file, 'rt') as zipf, open(_pdb_file, "w") as pdbf:
    #         try:
    #             pdbf.write(zipf.read())
    #         except IOError:
    #             RealtimeLogger.info("PDB IS", open(pdb_file).read())
    #             raise
    #
    #     try:
    #         os.remove(pdb_file)
    #     except OSError:
    #         pass
    #
    #     pdb_file = _pdb_file
    #
    # parser = PDBParser() if file_format == "pdb" else MMCIFParser()
    # structure = parser.get_structure('dummy', pdb_file)

    # if chain in ["", " "]:
    #     chain = next(structure.get_chains()).get_id()
    #     RealtimeLogger.info("NEW CHAIN IS '{}'".format(chain))

    mmdb_indices, mmdb_seq, num_residues, mmdb_pdb, mmdb_chain = get_entrez(gi, work_dir)

    if chain in ["", " "]:
        chain = mmdb_chain
        #RealtimeLogger.info("NEW CHAIN IS '{}'".format(chain))

    if mmdb_pdb.lower() != pdb_name.lower() or  mmdb_chain != chain:
        #RealtimeLogger.info("CHAINGING PDBs {}.{}=>{}.{}".format(pdb_name, chain, mmdb_pdb, mmdb_chain))
        new_res = list(mmdb_to_pdb_resi_obsolete(gi, mmdb_pdb, mmdb_chain, resi, \
            shrink_or_expand=shrink_or_expand, return_gi_len=return_gi_len,
            replace_nulls=False, work_dir=work_dir, job=job))

        if return_gi_len:
            num_residues = new_res[-1]
            new_res = new_res[:-1]
        else:
            num_residues = None

        if len(new_res) > 0:
            raise ChangePDB(gi, pdb_name, chain, mmdb_pdb, mmdb_chain, new_res, num_residues)
        else:
            raise InvalidSIFTS("Cannot map old codes: {} {} {} {} {} {} {}".format(pdb_name, chain, mmdb_pdb, mmdb_chain, mmdb_pdb_map, len(mmdb_pdb_map)))

    # try:
    #     pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() \
    #         if Polypeptide.is_aa(r) and r.get_id()[0] == ' ']
    # except KeyError:
    #     RealtimeLogger.info("FAILED key error {} {} {}".format(pdb_name, chain, resi))
    #     RealtimeLogger.info("Models={}; CHAINS={}; INPUT-CHAIN={}.{} ({}.{})".format(
    #         [m.get_id() for m in structure.get_models()],
    #         [c.get_id() for c in structure.get_chains()],
    #         pdb_name, chain, mmdb_pdb, mmdb_chain))
    #     if replace_nulls:
    #         for _ in range(len(resi)):
    #             yield None
    #     raise StopIteration

    #RealtimeLogger.info("MMDB PDB is {} {}".format(mmdb_pdb, mmdb_chain))

    if len(mmdb_indices) != len(mmdb_seq) and len(mmdb_seq) != num_residues:
        #RealtimeLogger.info("{} {} {}".format(mmdb_indices, mmdb_seq, num_residues))
        raise InvalidSIFTS("Cannot map old codes")

    if len(mmdb_indices) == 0:
        if replace_nulls:
            RealtimeLogger.info("mmdb_indices DNE {}".format(mmdb_indices))
            yield None
            raise StopIteration
        else:
            raise InvalidSIFTS("mmdb_indices DNE {}".format(mmdb_indices))


    last = None
    last_index = len(mmdb_indices)-1
    while last is None and last_index >= 0:
        try:
            last = mmdb_indices[last_index]
        except IndexError:
            pass
        last_index -= 1

    first = None
    first_index = 0
    while first is None and first_index < len(mmdb_indices):
        try:
            first = mmdb_indices[first_index]
        except:
            pass
        first_index += 1

    if None is (first, last):
        if replace_nulls:
            RealtimeLogger.info("mmdb_indices DNE {}".format(mmdb_indices))
            yield None
            raise StopIteration
        else:
            raise InvalidSIFTS("mmdb_indices DNE {}".format(mmdb_indices))

    #mmdb_pdb_map = [i+1 for i, ind in enumerate(mmdb_indices) if ind is not None]

    #MMDB uses 0-based indexing
    mmdb_pdb_map = {i:ind for i, ind in enumerate(mmdb_indices) if ind if ind is not None}#[ind for i, ind in enumerate(mmdb_indices) if ind is not None]
    if False and len(mmdb_pdb_map) != len(pdb_residues):
        if mmdb_pdb.lower() != pdb_name.lower():
            #RealtimeLogger.info("CHAINGING PDBs {}.{}=>{}.{}".format(pdb_name, chain, mmdb_pdb, mmdb_chain))
            new_res = list(mmdb_to_pdb_resi_obsolete(gi, mmdb_pdb, mmdb_chain, resi, \
                replace_nulls=replace_nulls, work_dir=work_dir, job=job))
            if len(new_res) > 0:
                raise ChangePDB(gi, pdb_name, chain, mmdb_pdb, mmdb_chain, new_res, num_residues)
            else:
                raise InvalidSIFTS("Cannot map old codes: {} {} {} {} {} {} {}".format(pdb_name, mmdb_pdb, chain, mmdb_pdb_map, pdb_residues, len(mmdb_pdb_map), len(pdb_residues)))

        else:
            failed = False
            try:
                pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() \
                    if Polypeptide.is_aa(r) and r.get_id()[0] != "W"]
            except KeyError:
                failed = True

            if not failed and len(mmdb_pdb_map) != len(pdb_residues):
                try:
                    pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() if r.get_id()[0] != "W"]
                except KeyError:
                    failed = True

            if not failed and len(mmdb_pdb_map) != len(pdb_residues):
                try:
                    pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() if r.get_id()[0] == " "]
                except KeyError:
                    failed = True

            if not failed and len(mmdb_pdb_map) != len(pdb_residues):
                try:
                    pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() if Polypeptide.is_aa(r)]
                except KeyError:
                    failed = True

            if not failed and len(mmdb_pdb_map) != len(pdb_residues):
                try:
                    pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() \
                        if Polypeptide.is_aa(r) and not (not Polypeptide.is_aa(r, standard=True) and r.get_id()[0].startswith("H"))]
                except KeyError:
                    failed = True

            if not failed and len(mmdb_pdb_map) != len(pdb_residues):
                try:
                    pdb_residues = [r.get_id() for r in structure[0][chain].get_residues() \
                        if Polypeptide.is_aa(r) and not (not Polypeptide.is_aa(r, standard=True) and r.get_id()[0].startswith("H"))]
                    pdb_residues = [k for k, g in groupby(pdb_residues)]
                except KeyError:
                    failed = True

            if not failed and len(mmdb_pdb_map) != len(pdb_residues):
                failed = True

            if failed:
                #RealtimeLogger.info("FAILED matching lengths {} {} {}".format(pdb_name, chain, resi))
                if replace_nulls:
                    for _ in range(len(resi)):
                        yield None
                        raise StopIteration
                else:
                    #RealtimeLogger.info("RES are {}".format([r.get_id() for r in structure[0][chain].get_residues()]))
                    raise RuntimeError("Cannot map old codes: {} {} {} {} {} {} {}".format(pdb_name, mmdb_pdb, chain, mmdb_pdb_map, pdb_residues, len(mmdb_pdb_map), len(pdb_residues)))

    #RealtimeLogger.info("{}".format(mmdb_pdb_map))

    for i, r in enumerate(resi):
        try:
            yield mmdb_pdb_map[r] #.index(r)
            #RealtimeLogger.info("pdb_index {}".format(pdb_index))
        except KeyError:
            if shrink_or_expand:
                if r <= first_index:
                    yield first
                elif r >= last_index:
                    yield last
                elif first_index <= r <= last_index:
                    if len(resi) == 2:
                        if r>resi[1-i]:
                            #Expand to closest residue forward
                            resi_range = range(r+1, last_index+1)
                        else:
                            #Shrink to closest residue backward
                            resi_range = range(r-1, first_index+1, -1)
                        for next_r in resi_range:
                            try:
                                yield mmdb_pdb_map[next_r]
                                break
                            except (KeyError, ValueError):
                                continue
                        else:
                            if replace_nulls:
                                #RealtimeLogger.info("not found pdb_index for {} {} {} ".format(r, first, last))
                                yield None
                                continue
                            else:
                                raise InvalidSIFTS("{} {} resnum {} does not exit {}".format(pdb_name, chain, r, mmdb_pdb_map))
                    else:
                        choose_resi = []
                        for range_ in (range(r+1, last_index+1), range(r-1, first_index+1, -1)):
                            for j, next_r in enumerate(range_):
                                try:
                                    choose_resi.append((j, mmdb_pdb_map[next_r]))
                                    break
                                except (KeyError, ValueError):
                                    continue

                        if len(choose_resi)> 0:
                            yield min(choose_resi, key=lambda x: x[0])[1]
                        else:
                            if replace_nulls:
                                yield None
                                continue
                            else:
                                raise InvalidSIFTS("{} {} resnum {} does not exit {}".format(pdb_name, chain, r, mmdb_pdb_map))

                elif replace_nulls:
                    yield None
                    continue
                else:
                    raise InvalidSIFTS("{} {} resnum {} does not exit {}".format(pdb_name, chain, r, mmdb_pdb_map))

            elif replace_nulls:
                yield None
                continue
            else:
                raise InvalidSIFTS("{} {} resnum {} does not exit {}".format(pdb_name, chain, r, mmdb_pdb_map))

    if return_gi_len:
        yield num_residues

        # try:
        #     yield pdb_residues[pdb_index]
        # except IndexError:
        #     if replace_nulls:
        #         yield None
        #     else:
        #         raise InvalidSIFTS("{} {} resnum {} does not exit {}".format(pdb_name, chain, r, sorted_map))

    #
    #     """
    #     1->0->7
    #     2->0->7
    #     """
    # else:
    #     last = mmdb_indices[-1]
    #     i=1
    #     while last is None:
    #         last = mmdb_indices[-i][1]
    #         i+=1
    #
    #     if r > last:
    #         yield last
    #     else:
    #         first = mmdb_indices[0]
    #         i=1
    #         while first is None:
    #             first = mmdb_indices[i]
    #             i+=1
    #
    #         if r < first:
    #             yield first
    #         else:
    #             raise InvalidSIFTS("?? resnum {} does not exit {}".format(r, sorted_map))
    #     #raise StopIteration
    #     continue
    #
    #     RealtimeLogger.info("ALIGNMENT NEEDED {} {} {}".format(mmdb_indices, mmdb_seq, num_residues))
    #     from Bio.SubsMat.MatrixInfo import blosum62
    #     from Bio.Align import PairwiseAligner
    #
    #     ppb = PPBuilder()
    #     pdb_seq = "".join([Polypeptide.three_to_one(r.get_resname()) for r in \
    #         structure[0][chain].get_residues() if Polypeptide.is_aa(r)])
    #
    #     aligner = PairwiseAligner()
    #     aligner.mode = 'global'
    #     aligner.open_gap_score = -0.5
    #     aligner.substitution_matrix = blosum62
    #     alignment = max(aligner.align(mmdb_seq, pdb_seq), key=lambda a: a.score)
    #     RealtimeLogger.info("ALIGNMENT: {}".format(alignment))
    #     query, _, target = str(alignment).splitlines()
    #
    #     query_map = [i for i, a in enumerate(query) if a != "-"]
    #     target_map = [i for i, a in enumerate(target) if a != "-"]
    #
    #     for r in resi:
    #         try:
    #             yield pdb_residues[target_map.index(query_map[mmdb_indices.index(r)])]
    #         except IndexError:
    #             if not replace_nulls:
    #                 yield None

    # complete_resi = ",".join(map(str,resi))
    #
    # resi = [natural_keys(r, use_int=True) for r in sorted(resi, key=natural_keys)]
    # for r in resi:
    #     try:
    #         print("CONV RESI", r[1], pdb_residues[r[1]])
    #         yield pdb_residues[r[1]]
    #     except IndexError:
    #         continue
    #
    # try:
    #     os.remove(pdb_file)
    # except OSError:
    #     pass



@contextmanager
def get_sifts(pdb, job=None):
    if job is not None:
        work_dir = job.fileStore.getLocalTempDir()
        in_store = IOStore.get("aws:us-east-1:molmimic-sifts")
        sifts_prefix = "{}/{}.xml.gz".format(pdb[1:3].lower(), pdb.lower())
        sifts_path = os.path.join(work_dir, os.path.basename(sifts_prefix))

        try:
            in_store.read_input_file(sifts_prefix, sifts_path)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise InvalidSIFTS("Cannot open {}".format(pdb))

        with open(sifts_path) as f:
            yield f

        os.remove(sifts_path)
    # else:
    #     path = os.path.join(os.environ.get("PDB_SNAPSHOT", os.path.join(data_path_prefix, "pdb")), "sifts", pdb[1:3].lower(), "{}.xml.gz".format(pdb.lower()))
    #     try:
    #         with open(path) as f:
    #             yield file
    #     except IOError as e:
    #         url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}.xml.gz".format(pdb.lower())
    #         try:
    #             sifts = requests.get(url)
    #             file = StringIO(sifts.content)
    #             yield file
    #         except requests.exceptions.RequestException:
    #             file.close()
    #             raise InvalidSIFTS("Not found: "+url+" orgig error:"+str(e))
    #         finally:
    #             file.close()
    #     os.remove(sifts)

def comare_to_pdb(pdb_file, resi):
    if len(resi) == 1 and isinstance(resi[0], str) and "," in resi[0]:
        resi = resi[0].split(",")

    resi = (natural_keys(r, use_int=True) for r in sorted(resi, key=natural_keys))

    for test_resi, known_resi in izip_missing(iter(resi), get_pdb_residues(pdb_file)):
        yield test_resi

def decode_residues(job, gi, pdb, chain, res, work_dir=None, row=None):
    if not pdb or pdb == np.NaN or not isinstance(pdb, str):
        raise InvalidSIFTS

    if work_dir is None:
        work_dir = os.getcwd()

    residues = []

    if res.startswith("0x"):
        res = res[2:]
    try:
        res = binascii.unhexlify(res)
    except:
        pass

    try:
        code, rest = decoder.decode(zlib.decompress(res, 16 + zlib.MAX_WBITS))
    except Exception as e:
        if type(res, str) and "," in res:
            return res
        else:
            return ""

    for i in range(len(code)):
        c = code[i]
        range_from, range_to, gi = tuple([c[j] for j in range(len(c))])
        for x in range(range_from, range_to + 1):
            residues.append(x)

    # try:
    #     return ",".join(map(str, mmdb_to_pdb_resi(pdb, chain, residues, job=job)))
    # except Exception as error:
    #     RealtimeLogger.info("TRY OBSOLETE")
    try:
        residues = ["".join(map(str, r)).strip() for r in mmdb_to_pdb_resi_obsolete(
            gi, pdb, chain, residues, replace_nulls=True, work_dir=work_dir,
            job=job) if r is not None]
        return ",".join(map(str, residues))
    except Exception as error:
        raise
        RealtimeLogger.info("Error mapping mmdb for", pdb, chain, error, row)
        raise InvalidSIFTS
        residues.insert(0, "mmdb")
        return ",".join(map(str,residues))
