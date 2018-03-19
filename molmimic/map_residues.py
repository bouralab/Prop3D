import os
import gzip
import urllib
from xml.etree.cElementTree import parse, ParseError

class InvalidSIFTS(RuntimeError):
    pass

def map_residues(pdb_name, chain, resi):
    if len(resi) == 1 and isinstance(resi[0], str) and "," in resi[0]:
        resi = map(int, resi[0].split(","))
    comlete_resi = ",".join(map(str,resi))

    resi = iter(sorted(resi))
    parsing = False
    current_resi = resi.next()

    sifts_path = get_sifts(pdb_name)
    if sifts_path is None:
        raise InvalidSIFTS(pdb_name)

    ns = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"

    with gzip.open(sifts_path, 'rb') as xml:
        try:
            tree = parse(xml)
        except ParseError:
            raise InvalidSIFTS(pdb_name)

        elem = tree.getroot()
        for entityNum, entity in enumerate(elem.findall(".//{}entity".format(ns))):
            try:
                dbChainId = entity.find(".//{ns}segment/{ns}listResidue/{ns}residue/{ns}crossRefDb[@dbSource='PDB']".format(ns=ns)).attrib["dbChainId"]
            except (AttributeError, KeyError):
                raise InvalidSIFTS(pdb_name)
            if  dbChainId == chain:
                for residue in entity.findall(".//{ns}segment/{ns}listResidue/{ns}residue".format(ns=ns)):
                    if int(residue.attrib["dbResNum"]) == current_resi+1:
                        try:
                            pdb = int(residue.findall('.//{}crossRefDb[@dbSource="PDB"]'.format(ns))[0].attrib["dbResNum"])
                        except (IndexError, ValueError) as e:
                            pdb = None

                        try:
                            resn = residue.findall('.//{}crossRefDb[@dbSource="NCBI"]'.format(ns))[0].attrib["dbResName"]
                        except (IndexError, ValueError) as e:
                            resn = None

                        try:
                            ncbi = int(residue.findall('.//{}crossRefDb[@dbSource="NCBI"]'.format(ns))[0].attrib["dbResNum"])
                        except (IndexError, ValueError) as e:
                            ncbi = None

                        yield current_resi, resn, pdb, ncbi
                        current_resi = resi.next()
                break

def get_sifts(pdb):
    path = os.path.join(os.environ.get("PDB_SNAPSHOT", "/pdb"), "sifts", pdb[1:3].lower(), "{}.xml.gz".format(pdb.lower()))
    if os.path.isfile(path):
        return path
    else:
        import urllib2
        path = "{}.xml.gz".format(pdb.lower())
        sifts = urllib2.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}".format(path))
        with open(path, "w") as f:
            f.write(sifts.read())
        return path
