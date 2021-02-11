from molmimic.parsers.container import Container

atom_type_h_bond_donor = {
    "H":  False, # Non H-bonding Hydrogen
    "HD": False, # Donor 1 H-bond Hydrogen
    "HS": False, # Donor S Spherical Hydrogen
    "C":  False, # Non H-bonding Aliphatic Carbon
    "A":  False, # Non H-bonding Aromatic Carbon
    "N":  False, # Non H-bonding Nitrogen
    "NA": True,  # Acceptor 1 H-bond Nitrogen
    "NS": True,  # Acceptor S Spherical Nitrogen
    "OA": True,  # Acceptor 2 H-bonds Oxygen
    "OS": True,  # Acceptor S Spherical Oxygen
    "F":  False, # Non H-bonding Fluorine
    "Mg": False, # Non H-bonding Magnesium
    "MG": False, # Non H-bonding Magnesium
    "P":  False, # Non H-bonding Phosphorus
    "SA": True,  # Acceptor 2 H-bonds Sulphur
    "S":  False, # Non H-bonding Sulphur
    "Cl": False, # Non H-bonding Chlorine
    "CL": False, # Non H-bonding Chlorine
    "Ca": False, # Non H-bonding Calcium
    "CA": False, # Non H-bonding Calcium
    "Mn": False, # Non H-bonding Manganese
    "MN": False, # Non H-bonding Manganese
    "Fe": False, # Non H-bonding Iron
    "FE": False, # Non H-bonding Iron
    "Zn": False, # Non H-bonding Zinc
    "ZN": False, # Non H-bonding Zinc
    "Br": False, # Non H-bonding Bromine
    "BR": False, # Non H-bonding Bromine
    "I":  False # Non H-bonding Iodine
}

class OpenBabel(Container):
    IMAGE = 'docker://edraizen/openbabel:latest'
    LOCAL = None
    PARAMETERS = [
        ("in_format", "str", "i"),
        ("in_file", "path:in", ""),
        ("out_format", "str", "o"),
        ("out_file", "path:out", "O")]
    RETURN_FILES = True
    ARG_START="-"

    def get_autodock_features(self, pdb_path):
        out_file = pdb_path+".pdbqt"
        pdbqt_file = self(in_format="pdb", in_file=pdb_path, out_format="pdbqt",
            out_file=out_file)
        autodock_features = self._get_autodock_features_from_pdbqt(pdbqt_file)
        self.files_to_remove += [out_file, pdbqt_file]
        self.clean()
        return autodock_features

    def _get_autodock_features_from_pdbqt(self, pdbqt_file):
        autodock_features = {}
        with open(pdbqt_file) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_serial = int(line[6:11])
                    autodock_type = line[77:79].strip().upper()
                    h_bond_donor = atom_type_h_bond_donor.get(autodock_type, False)
                    autodock_features[atom_serial] = (autodock_type, h_bond_donor)
        return autodock_features

def run_pdbqt_python(pdb_path):
    import pybel
    mol = next(pybel.readfile("pdb", pdb_path))

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

# def run_open_babel(in_format, in_file, out_format, out_file, work_dir=None, docker=True, job=None):
#     """Run APBS. Calculates correct size using Psize and defualt from Chimera
#     """
#     if work_dir is None:
#         work_dir = os.getcwd()
#
#     if docker and apiDockerCall is not None and job is not None:
#         parameters = ["-i{}".format(in_format), os.path.basename(in_file),
#             "-o{}".format(out_format), "-O", os.path.basename(out_file)]
#
#         if not os.path.abspath(os.path.dirname(in_file)) == os.path.abspath(work_dir):
#             shutil.copy(in_file, os.path.join(work_dir, in_file))
#
#         try:
#             apiDockerCall(job,
#                           image='edraizen/openbabel:latest',
#                           working_dir="/data",
#                           volumes={work_dir:{"bind":"/data", "mode":"rw"}},
#                           parameters=parameters)
#         except (SystemExit, KeyboardInterrupt):
#             raise
#         except:
#             raise
#
#     else:
#         raise RuntimeError("Openbabel needs to be run in docker")
#
#     out_file = os.path.join(work_dir, os.path.basename(out_file))
#     assert os.path.isfile(out_file), "Outfile not found: {}".format(os.listdir(work_dir))
#     return out_file
