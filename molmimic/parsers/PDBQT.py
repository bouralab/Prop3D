from tempfile import mkdtemp

from joblib import Memory
import pybel

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

@memory.cache
def run_pdbqt(struct, modified=False):
    if modified:
        mol = pybel.readstring("pdb", struct.save_pdb(file_like=True).read())
    else:
        mol = next(pybel.readfile("pdb", struct.path))

    mol.addh()

    pdbqt = mol.write("pdbqt")
    autodock_features = {}
    for atom_index in range(mol.OBMol.NumAtoms()):
        a = mol.OBMol.GetAtom(atom_index + 1)

        if a.IsCarbon() and a.IsAromatic():
            element = 'A'
        elif a.IsOxygen():
            element = 'OA'
        elif a.IsNitrogen() and a.IsHbondAcceptor():
            element = 'NA'
        elif a.IsSulfur() and a.IsHbondAcceptor():
            element ='SA'
        else:
            element = "".join([c for c in a.GetType() if c.isalnum()])

        autodock_features[a.GetIdx()] = (element, a.IsHbondDonor())
    return autodock_features

def get_autodock_features(struct, atom):
    autodock_features = run_pdbqt(struct, struct.modified_pdb_file)
    try:
        return autodock_features[atom.serial_number]
    except KeyError:
        return "  ", False
