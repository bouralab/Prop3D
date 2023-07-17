import h5py
import h5pyd
from time import sleep
from pathlib import Path
from Prop3D.util import safe_remove
from Prop3D.util.toil import map_job

from toil.job import Job
from toil.common import Toil
from toil.realtimeLogger import RealtimeLogger

def load_h5(job: Job, h5_file: str, hsds_file: str, prefix: str = "", hard_link_keys: list[str] = ["validation", "train", "test"], 
            hard_link_map: dict[str, str] = {"data_splits":"domains", "representatives":"domains"}, hardlink_items: bool = False):
    """A parallelized version of hsload
    """
    #RealtimeLogger.info(f"Adding h5_file ({h5_file_id}) to hsds_file ({hsds_file})")
                        
    if prefix == "":
        #Start of hierarchy
        try:
            with h5pyd.File(hsds_file, use_cache=False, retries=100):
                file_mode = "a"
        except IOError:
            file_mode = "w"
    else:
        file_mode = "a"

    RealtimeLogger.info(f"got file mode {file_mode}")


    #Save link to h5 file in file store
    #h5_file = job.fileStore.readGlobalFile(h5_file_id, cache=False, symlink=True)

    RealtimeLogger.info(f"downlaoded h5_file {h5_file}")

    if prefix.startswith('/'):
        prefix = prefix[1:]

    with h5py.File(h5_file, 'r') as f:
        for key in f[f'/{prefix}'].keys():
            full_key = f"/{prefix}/{key}"
            if all([k not in full_key for k in hard_link_map.keys()]):
                #if "domain" not in full_key or "data_splits" not in full_key or "representatives" not in full_key:
                RealtimeLogger.info(f"Adding {prefix}/{key}")
            elif hardlink_items:
                #old_kwd, new_kwd = next(hard_link_map.items())
                for old_kwd, new_kwd in hard_link_map.items():
                    if old_kwd in full_key:
                        store[full_key] = store[f"{prefix.split(old_kwd, 1)[0]}/{new_kwd}/{key}"]
                        RealtimeLogger.info(f"Hardlinked {full_key}")
                        break
                else:
                    raise RuntimeError(f"hard_link_map keys ('{list(old_kwd.keys())}')")
                
                continue

                
            h5_object = f[full_key]

            if isinstance(h5_object, h5py.Group):
                with h5pyd.File(hsds_file, mode=file_mode, use_cache=False, retries=100) as store:
                    store_object = store.require_group(full_key)

                    for attr, attr_value in h5_object.attrs.items():
                        store_object.attrs[attr] = attr_value

                if key not in hard_link_keys:
                    job.addChildJobFn(load_h5, h5_file, hsds_file, prefix=f"{prefix}/{key}",
                                    hard_link_keys=hard_link_keys, hard_link_map=hard_link_map, hardlink_items=False)
                else:
                    #Add hard links last
                    job.addFollowOnJobFn(load_h5, h5_file, hsds_file, prefix=f"{prefix}/{key}", 
                                    hard_link_keys=hard_link_keys, hard_link_map=hard_link_map, hardlink_items=True)
                    
            elif isinstance(h5_object, h5py.Dataset):
                dset = h5_object[:]
                shape = h5_object.shape
                dtype = h5_object.dtype
                dset_parameters = {k:getattr(h5_object, k) for k in ('compression',
                    'compression_opts', 'scaleoffset', 'shuffle', 'fletcher32', 'fillvalue')}
                dset_parameters['chunks'] = True
                
                retry = False
                should_remove = False
                with h5pyd.File(hsds_file, mode=file_mode, use_cache=False, retries=100) as store:
                    try:
                        store.require_dataset(full_key, shape, dtype, data=dset, **dset_parameters)
                    except (OSError, TypeError) as e:
                        if isinstance(e, OSError) and "Request Entity Too Large" in str(e):
                            #First time loading dataset and it is too large
                            retry = True
                            should_remove = True

                        elif isinstance(e, TypeError) or (isinstance(e, OSError) and "Dataset is not extensible" in str(e)):
                            #Previous error and datasets aren't the same size. Remove everything and try again
                            retry = True
                            should_remove = True

                if should_remove:
                    with h5pyd.File(hsds_file, mode=file_mode, use_cache=False, retries=100) as store:
                        try:
                            del store[full_key]
                        except IOError:
                            pass
                    
                    for _ in range(20):
                        with h5pyd.File(hsds_file, mode=file_mode, use_cache=False, retries=100) as store:
                            try:
                                store[full_key]
                                sleep(1)
                                continue
                            except KeyError:
                                break
                    else:
                        continue
                        raise RuntimeError(f"Unable to remove {full_key} to start over")

                if retry:
                    #Dataset too lareg to pass over http PUT
                    
                    #dset_parameters["chunks"] = (chunk_size,)
                    dset_parameters["maxshape"] = None

                    with h5pyd.File(hsds_file, mode=file_mode, use_cache=False, retries=100) as store:
                        new_dset = store.create_dataset(full_key, dset.shape[0], dtype, **dset_parameters)
                        chunk_size = min(new_dset.chunks[0], 500) #atoms in structure int(len(rec_arr)/4)

                    for start in range(0, len(dset), chunk_size):
                        small_data = dset[start:start+chunk_size]
                        with h5pyd.File(hsds_file, mode=file_mode, use_cache=False, retries=100) as store:
                            new_dset = store[full_key]
                            new_dset.resize(dset.shape[0] + small_data.shape[0], axis=0)
                            new_dset[-1:] = small_data

def save_h5(hsds_file: str, h5_file: str, prefix: str = "", hard_link_keys: list[str] = ["validation", "train", "test"], 
            hard_link_map: dict[str,str] = {"data_splits":"domains", "representatives":"domains"}, hardlink_items: bool = False) -> None:
    """A custom version of hsget to save a H5 file from a an h5 file on an hsds endpoint. Parallel writing does not work, so Toil is not used.

    Parameters
    ----------
    hsds_file : str
        Path to h5 file on HSDS endpoint
    h5_file : str
        New path on local server to save the remote h5 file to
    prefix : str
        The prefix to start saving from. Defualt is ''
    hard_link_keys : list
        If a group name is in the following keys, create a soft link to it so no extra data is saved.
    hard_link_map : dict
        Create a soft link by replacing to the path to new one. The path is split on the key.
    hardlink_items : bool
    
    """
    if prefix.startswith('/'):
        prefix = prefix[1:]

    add_keys = []
    hard_linked_keys = []
    with h5pyd.File(hsds_file, "r", use_cache=False, retries=100) as store, h5py.File(h5_file, 'a') as f:
        for key in f[f'/{prefix}'].keys():
            full_key = f"/{prefix}/{key}"
            if hardlink_items:
                #Convert to softlink
                old_kwd, new_kwd = next(hard_link_map.items())
                store[full_key] = h5py.SoftLink(f"{prefix.split(old_kwd, 1)[0]}/{new_kwd}/{key}")

                for old_kwd, new_kwd in hard_link_map.items():
                    if old_kwd in full_key:
                        f[full_key] = h5py.HardLink(f"{prefix.split(old_kwd, 1)[0]}/{new_kwd}/{key}")
                        RealtimeLogger.info(f"Hardlinked {full_key}")
                        break
                else:
                    raise RuntimeError(f"hard_link_map keys ('{list(old_kwd.keys())}')")

                continue
                
            hsds_object = store[full_key]

            if isinstance(hsds_object, h5pyd.Group):
                h5_object = f.require_group(full_key)
                for attr, attr_value in hsds_object.attrs.items():
                    h5_object.attrs[attr] = attr_value

                if key not in hard_link_keys:
                    add_keys.append(full_key)
                else:
                    hard_linked_keys.append(full_key)

            elif isinstance(h5_object, h5py.Dataset):
                dset = hsds_object[:]
                dset_parameters = {k:getattr(h5_object, k) for k in ('shape', 'dtype', 'chunks', 'compression',
                  'compression_opts', 'scaleoffset', 'shuffle', 'fletcher32', 'fillvalue')}
                f.create_dataset(full_key, data=dset, **dset_parameters)
    
    for key in add_keys:
        save_h5(hsds_file, h5_file, prefix=f"{prefix}/{key}",
                hard_link_keys=hard_link_keys, hard_link_map=hard_link_map, hardlink_items=False)
    
    #Add hard links last
    for hard_linked_key in hard_linked_keys:
        save_h5(hsds_file, h5_file, prefix=f"{prefix}/{hard_linked_key}",      
                hard_link_keys=hard_link_keys, hard_link_map=hard_link_map, hardlink_items=True)

if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--load", action="store_true", default=False)
    action.add_argument("--save", action="store_true", default=False)
    parser.add_argument("infile")
    parser.add_argument("outfile")
    options = parser.parse_args()

    if options.save:
        save_h5(options.infile, options.outfile)
    else:
        with Toil(options) as workflow:
            if not workflow.options.restart:
                #inputH5FileID = workflow.importFile(f"file://{str(Path(options.infile).absolute())}")
                job = Job.wrapJobFn(load_h5, str(Path(options.infile).absolute()), options.outfile)
                workflow.start(job)
            else:
                workflow.restart()