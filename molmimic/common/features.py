import os

import numpy as np
from sklearn.gaussian_process.kernels import RBF
from Bio import PDB
import freesasa

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from toil.realtimeLogger import RealtimeLogger

try:
    import pybel
except ImportError:
    pybel = None

from molmimic.util import InvalidPDB, natural_keys, silence_stdout, silence_stderr
from molmimic.parsers.FreeSASA import run_freesasa_biopython
from molmimic.parsers.Electrostatics import run_pdb2pqr_APBS
from molmimic.parsers.CX import run_cx
from molmimic.parsers.Consurf import run_consurf
from molmimic.parsers.DSSP import run_dssp
from molmimic.parsers.PDBQT import get_autodock_features

from molmimic.common.Structure import Structure, angle_between, get_dihedral
from molmimic.common.ProteinTables import vdw_radii, vdw_aa_radii, surface_areas, \
    hydrophobicity_scales

def number_of_features(only_aa=False, only_atom=False, non_geom_features=False,
  use_deepsite_features=False, course_grained=False):
    if course_grained:
        if non_geom_features:
            return 35
        elif only_aa:
            return 21
        else:
            return 45
    else:
        if non_geom_features:
            return 13
        elif only_atom:
            return 5
        elif only_aa:
            return 21
        elif use_deepsite_features:
            return 9
        else:
            return 73

def get_feature_names(only_aa=False, only_atom=False, use_deepsite_features=False, course_grained=False):
    if only_atom:
        return ["C", "N", "O", "S", "Unk_element"]
    if only_aa:
        return PDB.Polypeptide.aa3
    feature_names = [
        "C", "CT", "CA", "N", "N2", "N3", "NA", "O", "O2", "OH", "S", "SH", "Unk_atom",
        "C", "N", "O", "S", "Unk_element",
        "vdw_volume", "charge", "neg_charge", "pos_charge", "neutral_charge",
        "electrostatic_potential", "is_electropositive", "is_electronegative"
        "cx", "is_concave", "is_convex", "is_concave_and_convex",
        "hydrophobicity", "is_hydrophbic", "is_hydrophilic", "hydrophobicity_is_0"
        "atom_asa", "atom_is_buried", "atom_exposed",
        "residue_asa", "residue_buried", "residue_exposed"]
    feature_names += PDB.Polypeptide.aa3
    feature_names += ["Unk_residue", "is_helix", "is_sheet", "Unk_SS"]
    feature_names += [ #DeepSite features
        "hydrophobic_atom",
        "aromatic_atom",
        "hbond_acceptor",
        "hbond_donor",
        "metal",
        "is_hydrogen",
    ]
    fearues_names += [
        "conservation_normalized",
        "conservation_scale",
        "is_conserved"
    ]
    return feature_names

class ProteinFeaturizer(Structure):
    def __init__(self, path, pdb, chain, sdi, domNo, job, work_dir,
      input_format="pdb", force_feature_calculation=False):
        super(self, ProteinVoxelizer).__init__(path, pdb, chain, sdi, domain,
            input_format=input_format, feature_mode="w+" if force_feature_calculation else "r")

    def calculate_flat_features(self, course_grained=False, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False):
        if course_grained:
            features = [self.calculate_features_for_residue(r, only_aa=only_aa,
                non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for r in self.structure.get_residues()]
            return features, self.residue_features_file
        else:
            features = [self.calculate_features_for_atom(atom, only_aa=only_aa,
                only_atom=only_atom, non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for atom in self.structure.get_atoms()]
            return features, self.atom_features_file

    def get_features_per_atom(residue_list):
        """Get features for eah atom, but not organized in grid"""
        features = [self.get_features_for_atom(a) for r in residue_list for a in r]
        return features

    def get_features_per_residue(self, residue_list):
        features = [self.get_features_for_residue(a) for r in residue_list]
        return features

    def calculate_features_for_atom(self, atom, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, warn_if_buried=False):
        if use_deepsite_features:
            features = self.get_deepsite_features(atom)
            if warn_if_buried:
                is_burried = self.get_accessible_surface_area(atom)[-2]
        elif only_atom:
            features = self.get_element_type(atom)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(atom)[-2]
        elif only_aa:
            features = self.get_residue(atom)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(atom)[-2]
        elif non_geom_features:
            features = np.zeros(13)
            features[0:5] = self.get_element_type(atom)
            features[5:9] = self.get_charge_and_electrostatics(atom)
            features[9:13] = self.get_hydrophobicity(atom)
            is_buried = self.get_accessible_surface_area(atom)[-2]
        else:
            features = np.empty(self.n_atom_features)

            features[0:13]  = self.get_atom_type(atom)
            features[13:18] = self.get_element_type(atom)
            features[18:19] = self.get_vdw(atom)
            features[19:26] = self.get_charge_and_electrostatics(atom)
            features[26:30] = self.get_concavity(atom)
            features[30:34] = self.get_hydrophobicity(atom)
            features[34:40] = self.get_accessible_surface_area(atom)
            features[40:61] = self.get_residue(atom)
            features[61:64] = self.get_ss(atom)
            features[64:70] = self.get_deepsite_features(atom, calc_charge=False,
                calc_conservation=False)
            features[70:73] = self.get_evolutionary_conservation_score(atom)

            is_buried = bool(features[35])

        RealtimeLogger.info("Finished atom {}".format(atom))

        self.atom_features[atom.serial_number-1] = features

        if warn_if_buried:
            return features, is_buried
        else:
            return features

    def calculate_features_for_residue(self, residue, only_aa=False, non_geom_features=False,
      use_deepsite_features=False):
        """Calculate FEATUREs"""
        if non_geom_features:
            features = np.concatenate((
                self.get_residue(residue),
                self.get_charge_and_electrostatics(residue),
                self.get_hydrophobicity(residue),
                self.get_evolutionary_conservation_score(residue)), axis=None)
        elif only_aa:
            features = self.get_residue(residue)
        else:
            features = np.concatenate((
                self.get_charge_and_electrostatics(residue),
                self.get_concavity(residue),
                self.get_hydrophobicity(residue),
                self.get_accessible_surface_area(residue),
                self.get_residue(residue),
                self.get_ss(residue),
                self.get_evolutionary_conservation_score(residue)), axis=None)

        self.residue_features[residue.get_id()[1]-1] = features

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
        atom_type = np.zeros(13)
        try:
            index = atom_types.index(atom.get_name().strip())
            atom_type[index] = 1.
        except ValueError:
            atom_type[12] = 1.
        return atom_type

    def get_element_type(self, atom):
        """ELEMENT_IS_ANY
        ELEMENT_IS_C
        ELEMENT_IS_N
        ELEMENT_IS_O
        ELEMENT_IS_S
        ELEMENT_IS_OTHER"""
        elems = "CNOS"
        elem_type = np.zeros(5)
        try:
            index = elems.index(atom.element)
            elem_type[index] = 1.
        except ValueError:
            elem_type[4] = 1.
        return elem_type

    def get_vdw(self, atom):
        return np.array([vdw_radii.get(atom.element.title(), 2.0)])

    def get_charge_and_electrostatics(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            atoms = [self.get_charge_and_electrostatics(a) for a in residue]
            charge_value = np.sum([a[0] for a in atoms])
            electrostatic_pot_value = np.sum([a[4] for a in atoms])
            charge = np.zeros(7)
            charge[0] = charge_value
            charge[1] = float(charge_value < 0)
            charge[2] = float(charge_value > 0)
            charge[3] = float(charge_value == 0)
            charge[4] = electrostatic_pot_value
            charge[5] = float(electrostatic_pot_value < 0)
            charge[6] = float(electrostatic_pot_value > 0)
            return charge
        else:
            raise RuntimeError("Input must be Atom or Residue: {}".format(type(atom_or_residue)))

        if not hasattr(self, "_pqr"):
            self._pqr = run_pdb2pqr_APBS(self.pdb_path, pdb2pqr_whitespace=True,
                work_dir=self.work_dir, job=self.job)
        atom_id = atom.get_full_id()[3:5]

        if atom_id[1][1] != " ":
            #pdb2pqr removes alternate conformations and only uses the first
            atom_id = (atom_id[0], (atom_id[1][0], " "))

        try:
            charge_value, electrostatic_pot_value = self._pqr[atom_id]
        except KeyError:
            charge_value, electrostatic_pot_value = np.NaN, np.NaN

        charge = np.zeros(7)
        charge[0] = charge_value
        #print "charge", charge_value, float(charge_value < 0)
        charge[1] = float(charge_value < 0)
        charge[2] = float(charge_value > 0)
        charge[3] = float(charge_value == 0)
        charge[4] = electrostatic_pot_value
        charge[5] = float(electrostatic_pot_value < 0)
        charge[6] = float(electrostatic_pot_value > 0)
        return charge

    def get_concavity(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            concavity_value = np.mean([self.get_concavity(a)[0] for a in residue])
            #concave, convex, or both
            concavity = np.zeros(4)
            concavity[0] = concavity_value
            concavity[1] = float(concavity_value <= 2)
            concavity[2] = float(concavity_value > 5)
            concavity[3] = float(2 < concavity_value <= 5)
            return concavity
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if not hasattr(self, "_cx"):
            self._cx = run_cx(self.pdb_path, work_dir=self.work_dir, job=self.job)

        concavity_value = self._cx.get(atom.serial_number, np.NaN)
        #concave, convex, or both
        concavity = np.zeros(4)
        concavity[0] = concavity_value
        concavity[1] = float(concavity_value <= 2)
        concavity[2] = float(concavity_value > 5)
        concavity[3] = float(2 < concavity_value <= 5)
        return concavity

    def get_hydrophobicity(self, atom_or_residue, scale="kd"):
        assert scale in hydrophobicity_scales.keys()
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        hydrophobicity = np.zeros(4)
        try:
            resname = PDB.Polypeptide.three_to_one(residue.get_resname())

            hydrophobicity_value = hydrophobicity_scales[scale][resname]
            hydrophobicity[0] = hydrophobicity_value
            hydrophobicity[1] = float(hydrophobicity_value < 0)
            hydrophobicity[2] = float(hydrophobicity_value > 0)
            hydrophobicity[3] = float(hydrophobicity_value == 0)
        except KeyError:
            hydrophobicity[0] = np.NaN
            hydrophobicity[1] = np.NaN
            hydrophobicity[2] = np.NaN
            hydrophobicity[3] = np.NaN

        return hydrophobicity

    def get_accessible_surface_area(self, atom_or_residue):
        """Returns the ASA value from freesasa (if inout is Atom) and the DSSP
        value (if input is Atom or Residue)

        Returns
        -------
        If input is residue a 3-vector is returned, otherwise a 4-vector is returned
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if not hasattr(self, "_sasa"):
            self._sasa = run_freesasa_biopython(self.pdb_path)

        sasa, sasa_struct = self._sasa

        if is_atom:
            total_area = surface_areas.get(atom.element.title(), 1.0)
            try:
                selection = "sele, chain {} and resi {} and name {}".format(self.chain, atom.get_parent().get_id()[1], atom.get_id()[0])
                with silence_stdout(), silence_stderr():
                    selections = freesasa.selectArea([selection], sasa_struct, sasa)
                    atom_area = selections["sele"]
                fraction = atom_area/total_area
            except (KeyError, AssertionError, AttributeError, TypeError):
                raise
                atom_area = np.NaN
                fraction = np.NaN

        if not hasattr(self, "_dssp"):
            self._dssp = run_dssp(self.structure, self.pdb_path,
                work_dir=self.work_dir, job=self.job)

        try:
            residue_area = self._dssp[residue.get_full_id()[2:]][3]
        except (KeyError, AssertionError, AttributeError, TypeError):
            try:
                #Remove HETATMs
                residue_area = self._dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][3]
                if residue_area == "NA":
                    residue_area = np.NaN
            except (KeyError, AssertionError, AttributeError, TypeError):
                residue_area = np.NaN

        if is_atom:
            asa = np.zeros(6)
            asa[0] = atom_area
            asa[1] = float(fraction <= 0.2) #buried
            asa[2] = float(fraction > 0.2) #exposed
            asa[3] = residue_area
            asa[4] = float(residue_area < 0.2)
            asa[5] = float(residue_area >= 0.2)
        else:
            asa = np.zeros(3)
            asa[0] = residue_area
            asa[1] = float(residue_area < 0.2)
            asa[2] = float(residue_area >= 0.2)

        return asa

    def get_residue(self, atom_or_residue):
        """
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        residues = [0.]*(len(PDB.Polypeptide.aa1)+1) #one hot residue

        try:
            residues[PDB.Polypeptide.three_to_index(residue.get_resname())] = 1.
        except (ValueError, KeyError) as e:
            residues[-1] = 1.
        return np.array(residues)

    def get_ss(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if not hasattr(self, "_dssp"):
            self._dssp = run_dssp(self.structure, self.pdb_path,
                work_dir=self.work_dir, job=self.job)

        try:
            atom_ss = self._dssp[residue.get_full_id()[2:]][2]
        except (KeyError, AssertionError, AttributeError, TypeError):
            try:
                #Remove HETATMs
                atom_ss = self._dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][2]
            except (KeyError, AssertionError, AttributeError, TypeError):
                atom_ss = "X"

        ss = np.zeros(3)
        ss[0] = float(atom_ss in "GHIT")
        ss[1] = float(atom_ss in "BE")
        ss[2] = float(atom_ss not in "GHITBE")
        return ss

    def get_deepsite_features(self, atom, calc_charge=True, calc_conservation=True):
        """Use DeepSites rules for autodock atom types
        """
        if not hasattr(self, "_autodock"):
            self._autodock = get_autodock_features(self.pdb_path, work_dir=self.work_dir, job=self.job)

        try:
            element = self._autodock[atom.serial_number]
        except KeyError:
            element = "  "

        nFeatures = 6
        if calc_charge:
            nFeatures += 2
        if calc_conservation:
            nFeatures += 1
        features = np.zeros(nFeatures, dtype=bool)

        #hydrophobic
        features[0] = (element == 'C') | (element == 'A')

        #aromatic
        features[1] = element == 'A'

        #hbond_acceptor
        features[2] = (element == 'NA') | (element == 'NS') | (element == 'OA') | \
                      (element == 'OS') | (element == 'SA')

        #hbond_donor
        features[3] = (element != 'HS') | (element != 'HD')

        #metal
        features[4] = (element == 'MG') | (element == 'ZN') | (element == 'MN') | \
                      (element == 'CA') | (element == 'FE')

        #occupancies / excluded volume
        features[5] = (element != 'H') & (element != 'HS') & (element != 'HD')

        if calc_charge:
            #positive_ionizable
            charge = self.get_charge(atom)[0]
            features[6] = charge > 0

            #negative_ionizable
            features[7] = charge < 0

        if calc_conservation:
            features[nFeatures-1] = self.get_evolutionary_conservation_score(atom)[-1]

        return features.astype(float)

    def get_evolutionary_conservation_score(self, atom_or_residue):
        if not hasattr(self, "_consurf"):
            self._consurf = run_consurf(self, self.pdb, self.chain)

        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        try:
            normalized_score, conservation_score = self._consurf[residue.get_id()]
        except KeyError:
            normalized_score, conservation_score = np.NaN, np.NaN

        features = np.zeros(3)
        features[0] = normalized_score
        features[1] = conservation_score
        features[2] = float(conservation_score>=6)

        return features

    def calculate_graph(self, d_cutoff=100.):
        import networkx as nx
        structure_graph = nx.Graph()
        for r1, r2 in calculate_neighbors(self.structure, d_cutoff=d_cutoff):
            r1_id, r2_id = r1.get_id()[1], r2.get_id()[1]
            # if r1_id not in structure_graph:
            #     structure_graph.add_node(r1_id)
            # if r2_id not in structure_graph:
            #     structure_graph.add_node(r2_id)
            structure_graph.add_edge(r1_id, r2_id, attr_dict=self.get_edge_features(r1, r2))

        edge_file = os.path.join(self.work_dir, "{}.edges.gz".format(self.id))
        nx.write_edgelist(structure_graph, edge_file)
        return edge_file

    def get_edge_features(self, r1, r2):
        r1_pos = np.array([a.get_coord() for a in r1]).mean(axis=0)
        r2_pos = np.array([a.get_coord() for a in r2]).mean(axis=0)

        distance = np.linalg.norm(r1_pos-r2_pos)
        angle = angle_between(r1_pos, r2_pos)

        #Add in features from gregarious paper
        #Left handed or right handed: direction of cross product

        a = np.array(r1["O"].get_coord())-np.array(r1["N"].get_coord())
        b = np.array(r2["O"].get_coord())-np.array(r2["N"].get_coord())

        omega = angle_between(a, b)

        c_a, c_b = np.array(r1["O"].get_coord()), np.array(r2["O"].get_coord())
        mid_a = np.array([r1["O"].get_coord(), r1["N"].get_coord()]).mean(axis=0)
        mid_b = np.array([r2["O"].get_coord(), r2["N"].get_coord()]).mean(axis=0)

        theta = get_dihedral(c_a, mid_a, mid_b, c_b)

        chirality = int(theta > 0)

        return {
            "distance":distance,
            "angle":angle,
            "omega":omega,
            "theta":theta,
            "chirality":chirality
        }

    ### Preload features will remove

    def get_features_for_atom(self, atom, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, warn_if_buried=False, preload=True):
        """Calculate FEATUREs"""
        if isinstance(atom, PDB.Atom.DisorderedAtom):
            #All altlocs have been removed so onlt one remains
            atom = atom.disordered_get_list()[0]

        if preload and not self.force_feature_calculation and self.precalc_features is not None:
            try:
                features = self.precalc_features[atom.serial_number-1]
                is_buried = bool(features[35]) #Residue asa #[self.precalc_features[a.serial_number-1][31] for a in atom.get_parent()]
                # if asa > 0.0:
                #     asa /= surface_areas.get(atom.element.title(), 1.0)
                #     is_buried = asa <= 0.2
                # else:
                #     is_buried = False

                if use_deepsite_features:
                    feats = np.concatenate((
                        features[64:70],
                        features[20:22],
                        features[72:]))
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                if only_atom:
                    feats = features[13:18]
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                elif only_aa:
                    feats = features[40:61]
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                elif non_geom_features:
                    feats = np.concatenate((
                        features[13:18],
                        features[19:23],
                        features[30:33], #1]))
                        np.array([float(is_buried)])))
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                else:
                    if warn_if_buried:
                        return features, is_buried
                    else:
                        return features
            except ValueError as e:
                # print e
                # pass
                raise

    def get_features_for_residue(self, residue, only_aa=False, non_geom_features=False, use_deepsite_features=False, preload=True):
        """Calculate FEATUREs"""
        print("Using precalc features")
        try:
            features = self.precalc_features[residue.get_id()[1]-1]
            if non_geom_features:
                return np.concatenate((
                    features[15:36],
                    features[0:4],
                    features[8:12],
                    ))
            elif only_aa:
                return features[15:36]
            else:
                return features[:self.nFeatures]
        except ValueError:
            pass
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
