import os
import gzip
import itertools as it
import numpy as np
from sklearn.decomposition import PCA
from Bio import PDB
from Bio import SeqIO
#from pdb2pqr import mainCommand
#:import freesasa

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

hydrophobicity_scales = { 
    "kd": {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
           'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
           'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
           'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 },
    "biological": {
           "A": -0.11,"C":  0.13,"D": -3.49,"E": -2.68,"F": 0.32,
           "G": -0.74,"H": -2.06,"I":  0.60,"K": -2.71,"L": 0.55,
           "M":  0.10,"N": -2.05,"P": -2.23,"Q": -2.36,"R": -2.58,
           "S": -0.84,"T": -0.52,"V":  0.31,"W": -0.30,"Y": -0.68},
    "octanal":{
           "A": -0.50, "C":  0.02, "D": -3.64, "E": -3.63, 
           "F":  1.71, "G": -1.15, "H": -0.11, "I":  1.12, 
           "K": -2.80, "L":  1.25, "M":  0.67, "N": -0.85, 
           "P": -0.14, "Q": -0.77, "R": -1.81, "S": -0.46, 
           "T": -0.25, "V":  0.46, "W":  2.09, "Y":  0.71,}
    }

vdw = {'Ru': 1.2, 'Re': 4.3, 'Ra': 2.57, 'Rb': 2.65, 'Rn': 2.3, 'Rh': 1.22, 'Be': 0.63, 'Ba': 2.41, 'Bi': 1.73, 'Bk': 1.64, 'Br': 1.85, 'D': 1.2, 'H': 1.2, 'P': 1.9, 'Os': 1.58, 'Es': 1.62, 'Hg': 1.55, 'Ge': 1.48, 'Gd': 1.69, 'Ga': 1.87, 'Pr': 1.62, 'Pt': 1.72, 'Pu': 1.67, 'C': 1.775, 'Pb': 2.02, 'Pa': 1.6, 'Pd': 1.63, 'Cd': 1.58, 'Po': 1.21, 'Pm': 1.76, 'Ho': 1.61, 'Hf': 1.4, 'K': 2.75, 'He': 1.4, 'Md': 1.6, 'Mg': 1.73, 'Mo': 1.75, 'Mn': 1.19, 'O': 1.45, 'S': 1.8, 'W': 1.26, 'Zn': 1.39, 'Eu': 1.96, 'Zr': 1.42, 'Er': 1.59, 'Ni': 1.63, 'No': 1.59, 'Na': 2.27, 'Nb': 1.33, 'Nd': 1.79, 'Ne': 1.54, 'Np': 1.71, 'Fr': 3.24, 'Fe': 1.26, 'Fm': 1.61, 'B': 1.75, 'F': 1.47, 'Sr': 2.02, 'N': 1.5, 'Kr': 2.02, 'Si': 2.1, 'Sn': 2.17, 'Sm': 1.74, 'V': 1.06, 'Sc': 1.32, 'Sb': 1.12, 'Se': 1.9, 'Co': 1.13, 'Cm': 1.65, 'Cl': 1.75, 'Ca': 1.95, 'Cf': 1.63, 'Ce': 1.86, 'Xe': 2.16, 'Lu': 1.53, 'Cs': 3.01, 'Cr': 1.13, 'Cu': 1.4, 'La': 1.83, 'Li': 1.82, 'Tl': 1.96, 'Tm': 1.57, 'Lr': 1.58, 'Th': 1.84, 'Ti': 1.95, 'Te': 1.26, 'Tb': 1.66, 'Tc': 2.0, 'Ta': 1.22, 'Yb': 1.54, 'Dy': 1.63, 'I': 1.98, 'U': 1.75, 'Y': 1.61, 'Ac': 2.12, 'Ag': 1.72, 'Ir': 1.22, 'Am': 1.66, 'Al': 1.5, 'As': 0.83, 'Ar': 1.88, 'Au': 1.66, 'At': 1.12, 'In': 1.93}

maxASA = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0, "K": 201.0, "L": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}

parser = PDB.PDBParser()
writer = PDB.PDBIO()

class SelectChain(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        return chain.get_id() == self.chain

class Structure(object):
    def __init__(self, pdb, chain, path=None, snapshot=True):
        if path is None and not os.path.isfile(pdb) and len(pdb) == 4:
            if snapshot:
                path = "{}/pdb/{}/pdb{}.ent.gz".format(os.environ.get("PDB_SNAPSHOT", "/pdb"), pdb[1:3].lower(), pdb.lower())
            else:
                path = download_pdb(pdb)
            pdb = pdb.lower()
        elif path is None:
            path = pdb

        print path

        if os.path.isfile(path):
            if path.endswith(".gz"):
                with gzip.open(path, 'rb') as f:
                    self.structure = parser.get_structure(pdb, f)
                path = path[:-2]
            else:
                self.structure = parser.get_structure(pdb, path)
        else:
            raise RuntimeError("Invalid PDB id or file")

        self.pdb = pdb
        self.path = path
        self.dssp = None
        self.pqr = None
        self.sasa = None
        self.sasa_classes = None
        self.mean_coord = np.zeros(3)
        self.mean_coord_updated = False
        self.starting_residue = None
        self.starting_index_seq = None
        self.starting_index_struc = None
        self.chain = chain

        try:
            aa = SeqIO.index("{}/sequences/pdb_seqres.txt".format(os.environ.get("PDB_SNAPSHOT", "/pdb")), "fasta")
            self.aa = str(aa["{}_{}".format(pdb, chain.upper())].seq)
        except KeyError:
            self.aa = None
            print "Canot read seq file"

    @staticmethod
    def features_from_string(pdb, chain, resi, input_shape=(96,96,96), rotations=200):
        """Get features

        Paramters
        ---------
        pdb : str
            PDB ID
        chain : str
            Chain to extract
        resi : str
            Residue numbers separated by whitespace corresponding to positionts in the full protein sequence
        input_shape : 3-tuple
        rotations : int
        """
        s = Structure(pdb, chain).extract_chain(chain)
        s.align_to_pai()
        binding_site = [s.align_seq_to_struc(int(r), return_residue=True) for r in resi.split(",")]
        for rotation_num in self.rotate(rotations):
            data_grid, truth_grid = self.get_features(input_shape=input_shape, residue_list=binding_site)
            yield data_grid, truth_grid

    def get_atoms(self, include_hetatms=False):
        for a in self.structure.get_atoms():
            hetflag, resseq, icode = a.get_parent().get_id()
            if not include_hetatms and hetflag is not ' ':
                continue
            yield a

    def _get_dssp(self):
        if self.dssp is None:
            self.dssp = PDB.DSSP(self.structure[0], self.path, dssp='mkdssp')
        return self.dssp

    def _get_sasa(self):
        if self.sasa is None:
            self.sasa, self.sasa_classes = freesasa.calcBioPDB(self.structure)
        return self.sasa

    def _get_pqr(self):
        if self.pqr is None:
            pqr_out = "{}.pqr".format(os.path.splitext(self.path)[0])
            mainCommand(["pdb2pqr.py", "--ff=amber", "--whitespace", self.path, pqr_out])

            self.pqr = {}
            with open(pqr_out) as pqr:
                for line in pqr:
                    if line.startswith("REMARK") or line.startswith("HETATM"): continue
                    
                    fields = line.rstrip().split()
                    if len(fields) == 11:
                        recordName, serial, atomName, residueName, chainID, residueNumber, X, Y, Z, charge, radius = fields
                    elif len(fields) == 10:
                        recordName, serial, atomName, residueName, residueNumber, X, Y, Z, charge, radius = fields
                    else:
                        print len(fields)
                        raise RuntimeError("Invalid PQR File")
                    key = ((' ', int(residueNumber), ' '), (atomName.strip(), ' '))
                    self.pqr[key] = float(charge)
        return self.pqr

    def _mean_coord(self):
        if not self.mean_coord_updated:
            self.mean_coord = np.mean(self.get_coords(), axis=0)
            self.mean_coord_updated = True
        return self.mean_coord

    def extract_chain(self, chain):
        chain = chain.upper()
        out_path = "{}_{}.pdb".format(os.path.splitext(self.path)[0], chain)
        if not os.path.isfile(out_path):
            writer.set_structure(self.structure)
            writer.save(out_path, select=SelectChain(chain))
        return Structure(self.pdb, chain, path=out_path)

    def get_coords(self):
        return np.array([a.get_coord() for a in self.structure.get_atoms()])

    def orient_to_pai(self, flip_axis=(0.2, 0.2, 0.2)):
        coords = PCA(n_components = 3).fit_transform(self.get_coords())
        coords = flip_around_axis(coords, axis=flip_axis)
        self.update_coords(coords)

    def rotate(self, num=1):
        """Rotate structure in randomly in place"""
        for r in xrange(num):
            M = rotation_matrix(random=True)
            coords = np.dot(self.get_coords(), M)
            self.update_coords(coords)
            yield r

    def update_coords(self, coords):
        for atom, coord in it.izip(self.structure.get_atoms(), coords):
            atom.set_coord(coord)

    def create_full_volume(self, input_shape=(96, 96, 96)):
        truth_grid = np.zeros(list(input_shape)+[1])
        for atom in self.get_atoms():
            grid = self.get_grid_coord(atom, vsize=input_shape[0])
            truth_grid[grid[0], grid[1], grid[2], 0] = 1
        return truth_grid

    def align_seq_to_struc(self, *seq_num, **kwds):
        return_residue=kwds.get("return_residue", False)
        if self.starting_index_struc is None and self.aa is not None:
            first_struct_aa = self.structure.get_residues().next()
            self.starting_index_seq = 0
            self.starting_index_struc = first_struct_aa.get_id()[1]
            srarting_resname_struc = PDB.Polypeptide.three_to_index(first_struct_aa.get_resname())

            for i, seq_aa in enumerate(self.aa):
                if seq_aa == srarting_resname_struc:
                    starting_index_seq = i
                    break


        if len(seq_num) == 1 and isinstance(seq_num[0], str) and "," in seq_num[0]:
            seq_num = map(int, seq_num[0].split(","))

        mapped_resdues = []
        for num in seq_num:
            if num < self.starting_index_struc:
                index = int(num)-self.starting_index_seq+self.starting_index_struc
            else:
                index = num
            print num, "->", index
            if return_residue:
                mapped_resdues.append(self.structure[0][self.structure.get_chains().next().get_id()][index])
            else:
                mapped_resdues(index)
        print
        return mapped_resdues

    def get_features(self, input_shape=(96, 96, 96), residue_list=None, return_data=True, return_truth=True):
        if residue_list is not None:
            atoms = [a for r in residue_list for a in r]
        else:
            atoms = self.get_atoms()

        data_grid = np.zeros(list(input_shape)+[52])
        truth_grid = np.zeros(list(input_shape)+[1])
        for atom in atoms:
            grid = self.get_grid_coord(atom, vsize=input_shape[0])
            if return_data:
                data_grid[grid[0], grid[1], grid[2], :] = self.get_features_for_atom(atom)
            if return_truth:
                truth_grid[grid[0], grid[1], grid[2], 0] = 1

        if return_data and return_truth:
            return data_grid, truth_grid
        elif return_data:
            return data_grid
        else:
            return truth_grid

    def get_features_for_atom(self, atom):
        """Calculate FEATUREs"""
        atom_type             = self.get_atom_type(atom)
        #partial_charge        = None
        element_type          = self.get_element_type(atom)
        #hydroxyl              = None
        #amide                 = None
        #amine                 = None
        #carbonyl              = None
        #ring_system           = None
        #peptide               = None
        vdw_volume            = self.get_vdw(atom)
        charge                = self.get_charge(atom)
        neg_charge            = int(charge < 0)
        pos_charge            = int(charge > 0)
        neutral_charge        = int(charge == 0)
        #charge_with_his       = None
        hydrophobicity        = self.get_hydrophobicity(atom)
        #mobility              = None
        solvent_accessibility = self.get_accessible_surface_area(atom)

        residue               = self.get_residue(atom, one_hot=True)
        # residue_class1        = self.get_residue_class2(atom, one_hot=True)
        # residue_class2        = self.get_residue_class2(atom, one_hot=True)
        ss_class             = self.get_ss(atom)

        return atom_type+element_type+[vdw_volume, charge, neg_charge, pos_charge, neutral_charge, hydrophobicity]+solvent_accessibility+residue+ss_class

    def get_atom_type(self, atom):
        """
        ATOM_TYPE_IS_C
        ATOM_TYPE_IS_CT
        ATOM_TYPE_IS_CA
        ATOM_TYPE_IS_N
        ATOM_TYPE_IS_N2
        ATOM_TYPE_IS_N3
        ATOM_TYPE_IS_NA
        ATOM_TYPE_IS_O
        ATOM_TYPE_IS_O2
        ATOM_TYPE_IS_OH
        ATOM_TYPE_IS_S
        ATOM_TYPE_IS_SH
        ATOM_TYPE_IS_OTHER"""
        atom_types = ["C", "CT", "CA", "N", "N2", "N3", "NA", "O", "O2", "OH", "S", "SH"]
        return [int(atom.get_name().strip() == a) for a in atom_types]+[int(atom.get_name().strip() not in atom_types)]

    def get_element_type(self, atom):
        """ELEMENT_IS_ANY
        ELEMENT_IS_C
        ELEMENT_IS_N
        ELEMENT_IS_O
        ELEMENT_IS_S
        ELEMENT_IS_OTHER"""
        elems = "CNOS"
        return [int(atom.element == e) for e in elems]+[int(atom.element not in elems)]

    def get_vdw(self, atom):
        return vdw.get(atom.get_name().strip()[0], 0)

    def get_charge(self, atom):
        pqr = self._get_pqr()
        return pqr[atom.get_full_id()[3:]]

    def get_hydrophobicity(self, atom, scale="kd"):
        assert scale in hydrophobicity_scales.keys()
        try:
            return hydrophobicity_scales[scale][PDB.Polypeptide.three_to_one(atom.get_parent().get_resname())]
        except KeyError:
            return 0

    def get_accessible_surface_area(self, atom):
        sasa = self._get_sasa()
        dssp = self._get_dssp()
        atom_area = sasa.atomArea(atom.serial_number-1)
        residue_area = dssp[atom.get_full_id()[2:4]][3]
        return [atom_area, residue_area, int(residue_area<0.2), int(residue_area>=0.2)]

    def get_residue(self, atom, one_hot=False):
        """ADD:
        RESIDUE_NAME_IS_HOH
        """
        if one_hot:
            residue = [0]*(len(PDB.Polypeptide.aa1)+1) #one hot residue
            try:
                residue[PDB.Polypeptide.three_to_index(atom.get_parent().get_resname())] = 1
            except ValueError:
                residue[-1] = 1
            return residue
        else:
            return atom.get_parent().get_resname()

    def get_ss(self, atom):
        dssp = self._get_dssp()
        atom_ss = dssp[atom.get_full_id()[2:4]][2]
        ss = [int(atom_ss in "GHIT"), int(atom_ss in "BE"), int(atom_ss in "GHITBE")]
        return ss

    def get_grid_coord(self, atom, vsize=96, max_radius=40):
        center = np.array(atom.coord) - self._mean_coord()
        adjusted = center*(vsize/2.-1)/float(max_radius)
        translated = adjusted + (vsize-1)/2. # Translate center
        rounded = translated.astype(int) # Round components
        return rounded

    def get_atoms_from_grids(self, grids, vsize=96, max_radius=40):
        for atom in self.get_atoms():
            if self.get_grid_coord(atom, vsize, max_radius) in grids:
                yield atom

def download_pdb(id):
    import urllib
    if not os.path.isfile("{}.pdb".format(id)):
        urllib.urlretrieve("http://files.rcsb.org/download/{}.pdb".format(id), "{}.pdb".format(id))
    return "{}.pdb".format(id)

def flip_around_axis(coords, axis = (0.2, 0.2, 0.2)):
    'Flips coordinates randomly w.r.t. each axis with its associated probability'
    for col in xrange(3):
        if np.random.binomial(1, axis[col]):
            coords[:,col] = np.negative(coords[:,col])
    return coords

def rotation_matrix(random = False, theta = 0, phi = 0, z = 0):
    'Creates a rotation matrix'
    # Adapted from: http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
    # Initialization
    if random == True:
        randnums = np.random.uniform(size=(3,))
        theta, phi, z = randnums
    theta = theta * 2.0*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0  # For magnitude of pole deflection.
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    st = np.sin(theta)
    ct = np.cos(theta)
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
 
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    M = (np.outer(V, V) - np.eye(3)).dot(R)
 
    return M

if __name__ == "__main__":
  import sys
  assert len(sys.argv) > 1
  structure = Structure(sys.argv[1], panfs=False).extract_chain("A")
  num_atoms = sum(1 for a in structure.structure.get_atoms() if a.get_parent().get_id()[0] == " ")
  for atom in structure.get_atoms():
    features = structure.get_features_for_atom(atom)
    print atom.get_full_id(), atom.serial_number, "of", num_atoms, features, len(features)
    
