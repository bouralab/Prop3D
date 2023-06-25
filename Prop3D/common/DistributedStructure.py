import os
import copy

import numpy as np
import numpy.lib.recfunctions
from sklearn import preprocessing
from sklearn.decomposition import PCA
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBIO

import h5pyd
import pandas as pd

from Prop3D.common.AbstractStructure import AbstractStructure
from Prop3D.common.features import default_atom_feature_np, default_residue_feature_np, \
    atom_features, residue_features

residue_columns = ["residue_id", "chain", "bfactor", "occupancy", "X", "Y", "Z"]
atom_columns = ["serial_number", "atom_name", "residue_id", "chain", "bfactor",
"X", "Y", "Z"]
entity_levels = ["A", "R", "C", "M", "S"]

class DistributedStructure(AbstractStructure):
    """A structure class to deal with structures originated from a distributed
     HSDS instance.

    Parameters:
    -----------
    path : str
        HSDS endpoint to access structures
    key : str
        Key to access speficic protein inside the HDF file
    cath_domain_dataset : str
        The CATH superfamily if endpoint is setup to use CATH (use '/' instead of '.')
    coarse_grained: boolean
        Use a residue only model instead of an all atom model. Defualt False. Warning, not fully implemented.
    """
    def __init__(self, path, key, cath_domain_dataset=None, coarse_grained=False):
        self.path = path
        self.key = key
        self.f = None
        self.cath_domain_dataset_key = cath_domain_dataset

        if cath_domain_dataset is None:
            #Full key given
            self.f = h5pyd.File(path, use_cache=False)
            try:
                self.full_key = key
                cath_domain_dataset = self.f[key]
            except KeyError:
                raise RuntimeError(f"Structure with key {key} does not exist in {path}")
        elif isinstance(cath_domain_dataset, str):
            #Name of domain
            self.f = h5pyd.File(path, use_cache=False)
            try:
                self.full_key = f"{key}/domains/{cath_domain_dataset}"
                cath_domain_dataset = self.f[f"{key}/domains/{cath_domain_dataset}"]
            except KeyError:
                raise RuntimeError(f"Structure with key {key}/domains/{cath_domain_dataset} does not exist in {path}")
        elif not isinstance(cath_domain_dataset, h5pyd.Group):
            raise RuntimeError("cath_domain_dataset must be None (key suppllied w/ previous argument), a domain name within the key, or a h5pyd.Group")
        else:
            #Is is already a group
            self.full_key = cath_domain_dataset.name

        self.cath_domain_dataset = cath_domain_dataset

        self.cath_domain = self.full_key[1:].rsplit('/', 1)[1]
        self.pdb = key[:4]
        self.chain = key[4]
        self.domNo = key[5:]
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

        super().__init__(f"{key}-{self.cath_domain}", coarse_grained=coarse_grained)

        if self.f is not None:
            self.f.close()

    def deep_copy_feature(self, feature_name, memo):
        """Deep copy a  specific feature

        Parameters:
        -----------
        feature_name: str
            Feature name to copy
        memo :
            objects to pass to deepcopy

        Raises
        ------
        NotImeplementedError if no method to handle feature
        """
        if feature_name == "data":
            return copy.deepcopy(self.data, memo)
        if feature_name == "features":
            print("copying features")
            return copy.deepcopy(self.features, memo)
        elif feature_name == "f":
            return None
        elif feature_name == "cath_domain_dataset":
            return None
        else:
            raise NotImplementedError

    def get_atoms(self, atoms=None, include_hetatms=False, exclude_atoms=None, include_atoms=None):
        """Enumerate over all atoms with options to filter

        Parameters:
        -----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        data = self.data if atoms is None else atoms
        if include_atoms is not None:
            if not isinstance(include_atoms, (list, tuple)):
                include_atoms = [include_atoms]
            data = data[np.in1d(data["serial_number"], include_atoms)]
        elif exclude_atoms is not None:
            if not isinstance(exclude_atoms, (list, tuple)):
                exclude_atoms = [exclude_atoms]
            data = data[~np.in1d(data["serial_number"], exclude_atoms)]

        for a in data:
            yield a
    
    def get_residues(self):
        residues = np.unique(np.stack([a["residue_id"] for a in self.data]))
        for r in residues:
            yield self.data[self.data["residue_id"]==r]

    def unfold_entities(self, entity_list, target_level="A"):
        """Map lower level such as atoms (single row) into higher entites such as 
        residues (multiple rows). Only works for atoms and chainsAdapted from BioPython

        Parameters:
        -----------
        entity_list : list
            List of entites to unfold
        target_level : str 
            level to map to, eg: 'A' for atom, 'R' for residue

        Yields:
        -------
        Either single row atoms or multiple rows for residues
        """
        if target_level in ["C", "M", "S"]:
            #We only allow single chain, single model
            return self.features
        elif target_level not in ["R", "A"]:
            raise RuntimeError(f"{target_level}: Not an entity level.")

        if not isinstance(entity_list, (list, tuple)):
            entity_list = [entity_list]

        if target_level=="A":
            for e in entity_list:
                yield e

        else:
            residues = np.unique(np.stack([e["residue_id"] for e in entity_list]))
            for r in residues:
                yield self.data[self.data["residue_id"]==r]

    def save_pdb(self, path=None, header=None, file_like=False, rewind=True):
        """Write PDB to file

        Parameters:
        -----------
        path : None or str
            Path to save PDB file. If None, file_like needs to be True.
        header : str or list of strs
            Header string to write to the beginning of each PDB file
        file_like : boolean
            Return a StringIO object of the PDB file, do not write to disk. Default False.
        rewind : boolean
            If returning a file-like object, rewind the beginning of the file

        Returns:
        --------
        None or file-like object of PDB file data
        """
        writer = PDBIO()
        lines = not path and not file_like
        if path is None:
            path = StringIO()
        elif isinstance(path, str):
            path = open(path, "w")
        elif not hasattr(path, "write"):
            raise RuntimeError("path must be a filename, file-like object, or None (interpreted as StringIO)")

        if header is not None:
            new_header = ""
            for line in header.splitlines():
                if not line.startswith("REMARK"):
                    line = "REMARK {}".format(line)
                print(line.rstrip(), file=path)

        for atom in self.data:
            res_id = atom["residue_id"].astype(str)
            if res_id[-1].isalpha():
                resseq, icode = int(res_id[:-1]), res_id[-1]
            else:
                resseq, icode = int(res_id), " "

            s = writer._get_atom_line(
                Atom(name=atom["atom_name"].decode("utf-8").strip(), coord=atom[["X", "Y", "Z"]],
                     bfactor=atom["bfactor"], occupancy=1.0, altloc=" ",
                     fullname=atom["atom_name"].decode("utf-8"), serial_number=atom["serial_number"]),
                " ", #hetfield empty
                " ", #segid empty
                atom["serial_number"],
                atom["residue_name"].decode("utf-8") if "residue_name" in self.data.dtype.names else "UNK",
                resseq,
                icode,
                atom["chain"].decode("utf-8"),
            )
            path.write(s)

        if file_like:
            if rewind:
                path.seek(0)
            output = path
        elif lines:
            path.seek(0)
            output = path.read()
            path.close()
        else:
            output = path.name
            path.close()

        return output

    def get_bfactors(self):
        """Get bfactors for all atoms
        """
        return self.data["bfactor"]

    def write_features(self, path=None, key=None, features=None, coarse_grained=False, name=None, work_dir=None, force=None, multiple=False):
        """Write features to an hdf file

        Parameters:
        -----------
        path : str
            Path to sve HDF file
        key : str
            Key to save dataset inside HDF file
        features : str or list of strs
            Features to write
        course_grained : boolean
            Include features only at the residue level. Default False
        name : str
            File name to write file
        work_dir : None or str
            Directory to save file
        force : bool
            Not used
        multiple bool
            Not used
        """
        if path is None:
            if work_dir is None:
                work_dir = os.getcwd()

            if name is None:
                name = self.cath_domain_dataset

            path = os.path.join(work_dir, f"{name}_features.h5")

        if key is None:
            key = self.full_key if multiple else self.cath_domain_dataset_key
        
        pd.DataFrame.from_records(self.features).to_hdf(path, key)
            

        # if not path.startswith("/"):
        #     path = os.path.join(os.path.basename(self.path), path)

        # if not path.endswith(".h5"):
        #     name += ".h5"

        # if distributed and force == 2:
        #     with h5pyd.Folder(os.path.basename(self.path)) as folder:
        #         if os.path.basename(path) in folder:
        #             del folder[os.path.basename(path)]

        # if key is None:
        #     key = self.full_key

        # key = os.path.join(key, feature_key)

        # if features is None:
        #     data = self.features
        # else:
        #     assert set(features).issubset(set(self.features.dtype.names))
        #     data = self.features[features]

        # if not distributed:
        #     import h5py as h5
        # else:
        #     h5 = h5pyd

        # with h5.File(name, "a") as f:
        #     if force == 0 and key in f:
        #         raise RuntimeError(f"{key} already exists in {name}, not overwriting. Set force=True to overwrite")

        #     group = f.require_group(os.path.dirname(key))

        #     ds1 = f.create_table(self.key, data=data, chunks=True,
        #         compression="gzip", compression_opts=9)

    def add_features(self, coarse_grained=False, **features):
        """Add a feature column to dataset
        """
        assert [len(f)==len(self.features) for f in features.values()], "Features must contain same number of atoms (or residues is coarse grained)"

        str_type = np.dtype('O', metadata={'vlen': str})
        na = {float:np.nan, str_type:"", int:9999999}

        new_dt = [(name, "<f8") for name in features]
        new_dt = np.dtype(self.features.dtype.descr + new_dt)

        new_df = np.empty(self.features.shape, dtype=new_dt)

        for c in self.features.dtype.names:
            new_df[c] = self.features[c]

        for col, value in features.items():
            new_df[col] = value

        self.feature_names += list(features.keys())
        self.features = new_df #self.data[self.feature_names]

    def _to_unstructured(self, x):
        """Convert a rec array for atom(s) toa  regualr numpy array"""
        return np.lib.recfunctions.structured_to_unstructured(x)

    def get_coords(self, include_hetatms=False, exclude_atoms=None):
        if self.coords is None:
            self.coords = self._to_unstructured(
                self.pdb_info[["X", "Y", "Z"]]).round(decimals=4)
        return self.coords
    #
    # def orient_to_pai(self, random_flip=False, flip_axis=(0.2, 0.2, 0.2)):
    #     self.shift_coords_to_origin()
    #
    #     coords = PCA(n_components = 3).fit_transform(self.get_coords())
    #     if random_flip:
    #         coords = flip_around_axis(coords, axis=flip_axis)
    #
    #     self.update_coords(coords)
    #
    # def rotate(self, rvs=None, num=1, return_to=None):
    #     """Rotate structure in randomly in place"""
    #     for r in range(num):
    #         if rvs is None:
    #             M, theta, phi, z = special_ortho_group.rvs(3) #rotation_matrix(random=True)
    #         else:
    #             M=rvs
    #         self.shift_coords_to_origin()
    #         old_coords = self.get_coords()
    #         coords = np.dot(self.coords, M).round(decimals=4)
    #         self.update_coords(coords)
    #
    #         self.shift_coords(self.get_mean_coord() if return_to is None else return_to)
    #
    #         yield r, M
    #
    # def update_coords(self, coords):
    #     self.coords = coords
    #     self.mean_coord = None
    #     self.mean_coord_updated = False

    def update_bfactors(self, b_factors):
        self.data["bfactor"] = b_factors

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

    def remove_loops(self, verbose=False):
        ss_groups = self.get_secondary_structures_groups(verbose=verbose)
        ss_groups, leading_trailing_residues = ss_groups[0], ss_groups[-2]
        leading_trailing_residues = list(leading_trailing_residues.values())
        start = [np.concatenate(leading_trailing_residues[0])] if len(leading_trailing_residues)>0 else []
        assert len(ss_groups)>0 and max([len(g) for g in ss_groups])>0, f"Error with {self.name} {self.get_secondary_structures_groups(verbose=True)}"
        ss = np.concatenate([r for ss in ss_groups for r in ss])
        end = [np.concatenate(leading_trailing_residues[1])] if len(leading_trailing_residues)>1 else []
        self.data = np.concatenate((*start, ss, *end), dtype=self.data.dtype)
        self.pdb_info = self.data[atom_columns]
        self.features = self.data[self.feature_names]
