import os, sys
import networkx as nx
import numpy as np
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_index
from scipy.interpolate import Rbf
from sklearn.gaussian_process.kernels import RBF


try:
    STRUCTURE_PATH = os.environ["IBIS_STRUCTURES"]
except KeyError:
    raise RuntimeError("Must set IBIS_STRUCTURES env variable")

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def calculate_neighbors(struct, d_cutoff=100.0, level="R"):
    """
    Calculates intermolecular contacts in a parsed struct object.
    Modified from haddocking/prodigy

    Parameters
    ----------
    struct : Bio.PDB.Structure
        The structure object
    d_cuttoff: float
        Distance to find neighbors

    Returns
    -------
    A list of lists of nearby elements at the specified level: [(a1,b2),]
    """
    atom_list = list(struct.get_atoms())
    ns = NeighborSearch(atom_list)
    all_list = ns.search_all(radius=d_cutoff, level=level)

    if not all_list:
        raise ValueError('No contacts found for selection')

    return all_list

def get_node_features(r):
    res = r.get_resname()
    if res == "MSE":
        res = "MET"
    try:
        index = three_to_index(res)
        mask = [int(i==index) for i in range(20)]
        return mask
    except (KeyError, IndexError):
        return [0]*20

def get_edge_features(r1, r2):
    r1_pos = np.array([a.get_coords() for a in r1]).mean()
    r2_pos = np.array([a.get_coords() for a in r2]).mean()

    angle = angle_between(r1_pos, r2_pos)
    distance = RBF(1./18)(np.linalg.norm(r1_pos-r2_pos))

    #Add in features from gregarious paper
    #Left handed or right handed: direction of cross product
    cross = np.cross(r1_pos, r2_pos)
    cross_angle = angle_between(cross, r1_pos)
    chirality = int(np.cross(r1_pos, r2_pos) > 0)

    return {
        "angle":angle,
        "distance":distance,
        "cross_angle":cross_angle,
        "chirality":chirality
    }

def get_protein_graph(sfam_id, pdb, chain, sdi, domNo, d_cuttoff=100.0):
    """Get a protein graph using IBIS information. All domains have been spit out
    from their full chain structures and can be accessed using sfam_id, pdb,
    chain, sdi, and domNo.

    Parameters
    ----------
    sfam_id : int
        Superfamily of protein (usually from mol_superfam_id column)
    pdb : str
        4 letter pdb code  (usually from mol_pdb column)
    chain : str
        chain ID  (usually from mol_chain column)
    sdi : int
        Structural Domain ID  (usually from mol_sdi_id column)
    domNo : int
        Domain number  (usually from mol_domNo column)
    """
    parser = PDBParser()
    structure_path = os.path.join(STRUCTURE_PATH, pdb[1:3].lower(),
        "{}_{}_sdi{}_d{}".format(pdb.upper(), chain, sdi, domNo))
    structure = parser.get_structure('dummy', structure_path)

    structure_graph = nx.Graph()
    for r1, r2 in calculate_neighbors(structure, d_cutoff=d_cutoff):
        if r1 not in structure_graph:
            structure_graph.add_node(r1, attr_dict=get_node_features(r1))
        if r2 not in structure_graph:
            structure_graph.add_node(r2, attr_dict=get_node_features(r2))
        structure_graph.add_edge(r1, r2, attr_dict=get_edge_features(r1, r2))

    return structure_graph

if __name__ == "__main__":
    get_protein_graph(*sys.argv[1:])
