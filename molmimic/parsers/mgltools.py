from __future__ import print_function
import os
import tempfile

from molmimic.parsers.container import Container
from molmimic.util import safe_remove

from toil.realtimeLogger import RealtimeLogger

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

class MGLTools(Container):
    IMAGE = "docker://edraizen/mgltools"

class PrepareReceptor(MGLTools):
    PARAMETERS = [
        "/opt/mgltools/utilities/prepare_receptor4.py",
        ("receptor", "path:in", ["-r", "{}"]),
        ("out_file", "path:out", ["-o", "{}"]),
        (":v", "str"),
        (":A", "str"),
        (":C", "str", "c"),
        (":p", "str"),
        (":U", "str"),
        (":e", "str"),
        (":M", "str"),
        (":d", "str"),
        (":w", "str")
    ]
    RETURN_FILES = True
    ARG_START = "-"

    def convert_to_pdbqt(self, receptor, out_file=None, **kwds):
        if out_file is None:
            out_file = os.path.join(self.work_dir, os.path.basename(receptor)+".pdbqt")

        if "A" not in kwds:
            #Do not repair
            kwds["A"] = "None"
        if "U" not in kwds:
            #Do not cleanup -- we want to keep extra hydrogens
            kwds["U"] = "None"

        return self(receptor=receptor, out_file=out_file, **kwds)

    def get_autodock_atom_types(self, receptor, verify=False, **kwds):
        pdbqt_file = self.convert_to_pdbqt(receptor, **kwds)

        if verify:
            pdbqt_atoms = {}
            pdbqt_residues = {}
            for serial, atom_name, resid, charge, atom_type in self.parse_pdbqt(pdbqt_file):
                if atom_name[0].isdigit() and "H" in atom_name:
                    atom_name = f"{atom_name[1:]}{atom_name[0]}"
                if resid in pdbqt_residues:
                    pdbqt_residues[resid].add(atom_name)
                else:
                    pdbqt_residues[resid] = set([atom_name])
                pdbqt_atoms[serial] = (atom_name, resid)

            pdb_atoms = {}
            pdb_residues = {}
            for serial, atom_name, resid in self.parse_pdbqt(receptor, pdb=True):
                if resid in pdb_residues:
                    pdb_residues[resid].add(atom_name)
                else:
                    pdb_residues[resid] = set([atom_name])
                pdb_atoms[serial] = (atom_name, resid)

            missing = {resid: pdb_atom_names-pdbqt_residues[resid] for resid, pdb_atom_names in pdb_residues.items()}
            extra = {resid: pdbqt_residues[resid]-pdb_atom_names for resid, pdb_atom_names in pdb_residues.items()}

            assert len(set(pdb_residues.keys()).intersection(set(pdbqt_residues.keys()))) == len(pdb_residues.keys())
            assert sum([len(x) for x in missing.values()])==0, "Error residues in receptor do not match residues in the generated pdbqt file. (Missing: {}; Extra: {})".format(missing, extra)
            assert pdbqt_atoms==pdb_atoms

        autodock_atom_types = {atom[0]: (atom[-1], atom_type_h_bond_donor.get(atom[-1], False)) \
            for atom in self.parse_pdbqt(pdbqt_file)}
        return autodock_atom_types


    @staticmethod
    def parse_pdbqt(pdbqt_file, pdb=False):
        with open(pdbqt_file) as f:
            for i, line in enumerate(f):
                if line.startswith("ATOM"):
                    #if i<14: print(line)
                    serial = line[6:11]
                    atom_name = line[12:16]
                    resid = (line[21], line[17:20], line[22:27].strip()) #chain, resn, resi,
                    if not pdb:
                        charge = line[70:76]
                        atom_type = line[77:79]
                        yield serial, atom_name, resid, charge, atom_type
                    else:
                        yield serial, atom_name, resid
