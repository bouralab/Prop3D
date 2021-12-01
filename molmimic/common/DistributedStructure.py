import numpy as np
import numpy.lib.recfunctions
from sklearn import preprocessing
from sklearn.decomposition import PCA

from molmimic.common.features import default_atom_feature_np, default_residue_feature_np, \
    atom_features, residue_features

residue_columns = ["residue_id", "chain", "bfactor", "occupancy", "X", "Y", "Z"]
atom_columns = ["serial_number", "atom_name", "residue_id", "chain", "bfactor",
"X", "Y", "Z"]

class DistributedStructure(object):
    def __init__(self, path, key, cath_domain_dataset, coarse_grained=False, file_mode="r"):
        self.path = path
        self.key = key

        if cath_domain_dataset is None:
            if file_mode == "r":
                raise RuntimeError("A h5 dataset must be passed in with param cath_domain_dataset")
            f = h5pyd.File(path, cache=False)
            try:
                group = f[key]
            except KeyError:
                raise RuntimeError("Must create datasets first")

            if coarse_grained:
                rec_array = default_residue_feature_np() #Need regular Structure object...
            else:
                rec_aray = default_atom_feature_np()

            cath_domain_dataset = f.create_dataset(key, rec_arr.shape, data=rec_arr)

        self.cath_domain_dataset = cath_domain_dataset

        self.cath_domain = key
        self.pdb = key[:4]
        self.chain = key[4]
        self.domNo = key[5:]
        self.file_mode = file_mode
        self.coarse_grained = coarse_grained

        if coarse_grained:
            self.data = self.cath_domain_dataset["residue"][:]
            self.pdb_info = self.data[residue_columns]
            self.feature_names = [name for name in self.data.dtype.names if name not in residue_columns]
            self.features = self.data[self.feature_names]
        else:
            self.data = self.cath_domain_dataset["atom"][:]
            self.pdb_info = self.data[atom_columns]
            self.feature_names = [name for name in self.data.dtype.names if name not in atom_columns]
            self.features = self.data[self.feature_names]

        self.n = len(self.data)
        self.coords = None
        self.get_coords()

    def get_atoms(self, include_hetatms=False, exclude_atoms=None):
        for a in self.data:
            yield a

    def get_surface(self):
        surface = self.features[self.features["residue_buried"]==0]
        return surface

    def get_mean_coord(self):
        if not self.mean_coord_updated:
            self.mean_coord = np.around(np.nanmean(self.coords, axis=0), decimals=4)
            self.mean_coord_updated = True
        return self.mean_coord

    def get_max_coord(self):
        return np.nanmax(self.coords, axis=0)

    def get_min_coord(self):
        return np.nanmin(self.coords, axis=0)

    def get_max_length(self, buffer=0, pct_buffer=0):
        length = int(np.ceil(np.linalg.norm(self.get_max_coord()-self.get_min_coord())))
        if pct_buffer!=0:
            length += int(np.ceil(length*pct_buffer))
        else:
            length += buffer
        if length%2 == 1:
            length += 1
        return length

    def shift_coords(self, new_center=None, from_origin=True):
        """if new_center is None, it will shift to the origin"""
        coords = self.coords
        if from_origin or new_center is None:
            mean_coord = self.get_mean_coord()
            coords -= mean_coord
        if new_center is not None:
            coords += np.around(new_center, decimals=4)

        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.nanmean(coords, axis=0)

    def shift_coords_to_origin(self):
        return self.shift_coords()

    def get_coords(self, include_hetatms=False, exclude_atoms=None):
        if self.coords is None:
            self.coords = numpy.lib.recfunctions.structured_to_unstructured(
                self.pdb_info[["X", "Y", "Z"]]).round(decimals=4)
        return self.coords

    def orient_to_pai(self, random_flip=False, flip_axis=(0.2, 0.2, 0.2)):
        self.shift_coords_to_origin()

        coords = PCA(n_components = 3).fit_transform(self.get_coords())
        if random_flip:
            coords = flip_around_axis(coords, axis=flip_axis)

        self.update_coords(coords)

    def rotate(self, rvs=None, num=1, return_to=None):
        """Rotate structure in randomly in place"""
        for r in range(num):
            if rvs is None:
                M, theta, phi, z = rotation_matrix(random=True)
                #print(M, theta, phi, z)
            else:
                M=rvs
            self.shift_coords_to_origin()
            old_coords = self.get_coords()
            coords = np.dot(self.coords, M).round(decimals=4)
            self.update_coords(coords)
            # if rvs is None or rvs!=np.eye(3):
            #     assert not np.array_equal(old_coords, self.get_coords()), M
            #self.shift_coords_to_volume_center()

            self.shift_coords(self.get_mean_coord() if return_to is None else return_to)
            # if rvs is None or rvs!=np.eye(3):
            #     assert not np.array_equal(coords, self.get_coords()), M

            yield r, M

    def update_coords(self, coords):
        self.coords = coords
        self.mean_coord = None
        self.mean_coord_updated = False

    def update_bfactors(self, b_factors):
        self.pdb_info["bfactor"] = b_factors

    def calculate_neighbors(self, d_cutoff=100.0):
        """
        Calculates intermolecular contacts in a parsed struct object.

        Parameters
        ----------
        d_cuttoff: float
            Distance to find neighbors

        Returns
        -------
        A list of lists of nearby elements at the specified level: [(a1,b2),]
        """
        ns = spatial.cKDTree(self.coords)
        return ns.query_pairs(d_cutoff)

    def get_vdw(self, atom_or_residue):
        return atom_or_residue["vdw"]


def flip_around_axis(coords, axis = (0.2, 0.2, 0.2)):
    'Flips coordinates randomly w.r.t. each axis with its associated probability'
    for col in range(3):
        if np.random.binomial(1, axis[col]):
            coords[:,col] = np.negative(coords[:,col])
    return coords

def rotation_matrix(random = False, theta = 0, phi = 0, z = 0, uniform=True):
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

    return M, theta, phi, z

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

def get_dihedral(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

    return M, theta, phi, z
