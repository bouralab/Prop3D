import os
import time
from tempfile import gettempdir

from Bio import PDB
from Bio.PDB import MMCIFParser

tempdir = gettempdir()

class SelectChain(PDB.Select):
    """ Only accept the specified chain and remove hydruogens and hetatms when saving. """
    def __init__(self, chain, model=0):
        self.model = model
        self.chain = chain

    def accept_model(self, model):
        return model.get_id() == self.model

    def accept_chain(self, chain):
        return chain.get_id() == self.chain

def mmcif2pdb(pdb, chain, model=0):
    parser = MMCIFParser()
    writer = PDB.PDBIO()

    remove = False
    if not os.path.isfile(pdb) and len(pdb) == 4:
        pdb_code = pdb
        pdb = os.path.join(tempdir, "{}.cif".format(pdb.lower()))
        attempts = 0
        while not os.path.isfile(pdb) and attempts < 5:
            try:
                pdb = PDB.PDBList().retrieve_pdb_file(pdb_code, pdir=tempdir, file_format="mmCif")
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                pass
            attempts += 1
            time.sleep(1)

        if not os.path.isfile(pdb):
            raise IOError("Cannot download file")
        remove = True

    name = os.path.splitext(pdb)[0]

    try:
        structure = parser.get_structure(name, pdb)
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        raise IOError("Unable to open pdb {}".format(pdb))

    try:
        writer.set_structure(structure)
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        raise IOError("Unable to save pdb {}".format(pdb))

    if not isinstance(chain, (list, tuple)):
        chain = [chain]

    for c in chain:
        if c is None:
            yield None
            continue
        
        #Make chain one character to fit in PDB format
        new_pdb = "{}_{}.pdb".format(name, c)
        new_chain = c[0]
        try:
            writer.save(new_pdb, select=SelectChain(new_chain, model))
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            raise IOError("Unable to save pdb {}".format(new_pdb))

        if not os.path.isfile(new_pdb):
            raise IOError("Cannot extract chain")

        yield new_pdb

    if remove:
        try:
            os.remove(pdb)
        except OSError:
            pass
