import os
import traceback
from functools import partial
from typing import Union, Any

import h5pyd
from Prop3D.common.featurizer import ProteinFeaturizer

from Prop3D.util import safe_remove
from Prop3D.generate_data.data_stores import data_stores

from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

class CalculateFeaturesError(RuntimeError):
    """A new Error to handle errors during protein featurization

    Parameters
    ----------
    job : toil Job
        Current running job
    cath_domain : str
        Name of cath domain
    stage : str
        function name where error occured
    message : str
        Error message
    errors : list or None
        Orignal error objects
    *args : any
        any xtra args
    **kwds : any 
        ANY XTRA KWDS
    """
    def __init__(self, job: Job, cath_domain: str, stage: str, message:str, errors: Union[list, None] = None, *args: Any, **kwds: Any) -> None:
        super().__init__(*args, **kwds)
        self.cath_domain = cath_domain
        self.stage = stage
        self.message = message
        self.errors = errors if isinstance(errors, list) else []
        self.jobStoreName = os.path.basename(job.fileStore.jobStore.config.jobStore.split(":")[-1])

    def __str__(self):
        """Convert errors into string
        """
        return "Error during {}: {}\nErrors:\n".format(self.stage, self.message,
            "\n".join(map(str, self.errors)))

    def save(self, store=None):
        """Save errors to file in an IOStore

        Parameters
        ----------
        store : IOStore
        """
        if store is None:
            store = data_stores(job).cath_features
        fail_file = "{}.{}".format(self.cath_domain, self.stage)
        with open(fail_file, "w") as f:
            print(self.message, file=f)
            print(self.errors, file=f)

        store.write_output_file(fail_file,
            f"errors/{self.jobStoreName}/{os.path.basename(fail_file)}")
        safe_remove(fail_file)

def calculate_features(job: Job, cath_full_h5: str, cath_domain: str, cathcode: str, update_features: Union[list[str], None] = None, 
                       domain_file: Union[str, None] = None, work_dir: Union[str, None] = None, edge_features: bool = True) -> None:
    """Featurize a protein at the atom, residue, and graph level saving all data into the h5 file on HSDS endpoint

    Parameters
    ----------
    job : toil Job
        Currently running job
    cath_full_h5 : str
        Path to h5 on hsds endpoint
    cath_domain : str
        CATH domain (7-letter code) PDB ID, CHAIN, Domain ID, eg. 1zyzA00
    cathcode : str
        Superfamily cath domain belongs to (Use / instead of .)
    update_features : list of str or None
        Select which features update (either indidual feature names or whole group names). If None, all features will be calculated.
        Default is None.
    domain_file : str or None
        Path to pdb file. If None, it will be downloaded from the raw IOStore (see data_stores)
    work_dir : str
        Where to save temp files
    edge_features: bool
        Include edge feature or not
    """
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    to_remove = []

    if cathcode is not None:
        cath_key = f"/{cathcode}/domains/{cath_domain}"
        s3_cath_key = "{}/{}".format(cathcode, cath_domain)
    elif os.path.isfile(cath_domain):
        cath_key = os.path.splitext(os.path.basename(cath_domain))[0]
        s3_cath_key = None
    else:
        cath_key = cath_domain
        s3_cath_key = None


    if update_features is not None:
        #Save features for cath domain in new seperate h5 files to be red in by the Featurizer
        store = h5pyd.File(cath_full_h5, mode="r", use_cache=False)
        try:
            feat_files = list(store[cath_key].keys())
            if len(feat_files) == 3 and "atom" in feat_files and \
                "residue" in feat_files and feat_files and "edges":
                for feature_type, index_col in (("atoms", "serial_number"), ("residues", "residue_id")):
                    df = pd.DataFrame(store[f"{cath_key}/{feature_type}"]).set_index(index_col)
                    feature_file = os.path.join(work_dir, f"{cath_domain}_{feature_type:-1]}.h5")
                    df.to_hdf(feature_file, "table")
                    del df
                    to_remove.append(feature_file)
            else:
                update_features = None
        except KeyError:
            feats_exist = False
        finally:
            store.close()

    if s3_cath_key is not None:
        domain_file = os.path.join(work_dir, "{}.pdb".format(cath_domain))
        try:
            data_stores(job).prepared_cath_structures.read_input_file(
                s3_cath_key+".pdb", domain_file)
        except Exception as e:
            RealtimeLogger.info("Failed to download prepared cath file {}".format(
                cath_key+".pdb"))
            raise
        output_name = cath_domain
    else:
        domain_file = domain_file if domain_file is not None else cath_domain
        cath_domain = None
        output_name = cath_key

    try:
        structure = ProteinFeaturizer(
            domain_file, cath_domain, job, work_dir,
            force_feature_calculation=update_features is None,
            update_features=update_features)
    except:
        import traceback as tb
        RealtimeLogger.info(f"{tb.format_exc()}")
        raise

    for ext, calculate in (("atom", structure.calculate_flat_features),
                           ("residue", structure.calculate_flat_residue_features),
                           ("edges", partial(structure.calculate_graph, edgelist=True))):
        try:
            out, _ = calculate(write=False)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            tb = traceback.format_exc()
            raise
            CalculateFeaturesError(job, cath_domain, ext.split(".",1)[0], tb).save()
            return

        if ext=="edges":
            df = out
            special_col_types = {"src":"<S8", "dst":"<S8"}
            df["src"] = df["src"].apply(lambda s: "".join(map(str,s[1:])).strip())
            df["dst"] = df["dst"].apply(lambda s: "".join(map(str,s[1:])).strip())
        else:
            del out
            df = structure.get_pdb_dataframe(include_features=True, coarse_grained = ext=="residue")
            special_col_types = {"serial_number":"<i8", "atom_name":"<S5",
                "residue_id":"<S8", "residue_name":"<S8", "chain":"<S2"}

        column_dtypes = {col:special_col_types.get(col, '<f8') for col in df.columns}
        rec_arr = df.to_records(index=False, column_dtypes=column_dtypes)

        with h5pyd.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
            if f"{cath_key}/{ext}" in store.keys():
                try:
                    del store[f"{cath_key}/{ext}"]
                except OSError:
                    pass

            if f"{cath_key}/{ext}" in store.keys():
                try:
                    del store[f"{cath_key}/{ext}"]
                except:
                    pass
            try:
                ds1 = store.create_table(f"{cath_key}/{ext}", data=rec_arr, dtype=list(column_dtypes.items()),
                    chunks=True, compression="gzip", compression_opts=9)
            except OSError as e:
                if "Request Entity Too Large" in str(e):
                    #Dataset too lareg to pass over http PUT
                    span = 500 #atoms in structure int(len(rec_arr)/4)
                    for i, start in enumerate(range(0, len(rec_arr), span)):
                        small_data = rec_arr[start:start+span]
                        if i==0:
                            RealtimeLogger.info(f"Create small data: with: {len(small_data)}")
                            store.create_table(f"{cath_key}/{ext}", data=small_data, dtype=list(column_dtypes.items()),
                                chunks=True, compression="gzip", compression_opts=9)
                        else:
                            RealtimeLogger.info(f"Add small data: with: {len(small_data)}")
                            store[f"{cath_key}/{ext}"].resize((store[f"{cath_key}/{ext}"].shape[0] + small_data.shape[0]), axis=0)
                            store[f"{cath_key}/{ext}"][-small_data.shape[0]:] = small_data

        RealtimeLogger.info("Finished {} features for: {} {}".format(ext, cathcode, output_name))

    RealtimeLogger.info("Finished features for: {} {}".format(cathcode, output_name))

    safe_remove(domain_file)

    if update_features:
        for f in to_remove:
            safe_remove(f)
