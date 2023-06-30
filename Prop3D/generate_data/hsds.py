import h5py
import h5pyd
from pathlib import Path
from Prop3D.util import safe_remove
from Prop3D.util.toil import map_job
from toil.realtimeLogger import RealtimeLogger

def load_h5(job, h5_file, hsds_file, prefix="", hard_link_keys=["validation", "training", "test"], hard_link_map={"data_splits":"domains"}, hardlink_items=False):
    """A parallelized version of hsload
    """
    #RealtimeLogger.info(f"Adding h5_file ({h5_file_id}) to hsds_file ({hsds_file})")
                        
    if prefix == "":
        #Start of hierarchy
        try:
            with h5pyd.File(hsds_file, use_cache=False):
                file_mode = "a"
        except IOError:
            file_mode = "w"
    else:
        file_mode = "a"

    RealtimeLogger.info(f"got file mode {file_mode}")


    #Save link to h5 file in file store
    #h5_file = job.fileStore.readGlobalFile(h5_file_id, cache=False, symlink=True)

    RealtimeLogger.info(f"downlaoded h5_file {h5_file}")

    with h5py.File(h5_file, 'r') as f:
        with h5pyd.File(hsds_file, mode=file_mode, use_cache=False) as store:
            for key in f[f'/{prefix}'].keys():
                if "domain" not in prefix or "data_splits" not in prefix or "representatives" not in prefix:
                    RealtimeLogger.info(f"Adding {prefix}/{key}")

                if hardlink_items:
                    old_kwd, new_kwd = next(hard_link_map.items())
                    store[f"{prefix}/{key}"] = store[f"{prefix.split(old_kwd, 1)[0]}/{new_kwd}/{key}"]
                    continue
                    
                h5_object = f[f'/{prefix}/{key}']

                if isinstance(h5_object, h5py.Group):
                    store_object = store.require_group(f'/{prefix}/{key}')

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
                    dset_parameters = {k:getattr(h5_object, k) for k in ('chunks', 'compression',
                        'compression_opts', 'scaleoffset', 'shuffle', 'fletcher32', 'fillvalue')}
                    
                    try:
                        dset = store.require_dataset(f'/{prefix}/{key}', shape, dtype, data=dset, **dset_parameters)
                    except (OSError, TypeError) as e:
                        retry = False
                        if isinstance(e, OSError) and "Request Entity Too Large" in str(e):
                            #First time loading dataset and it is too large
                            retry = True
                        elif isinstance(e, TypeError):
                            #Previous error and datasets aren't the same size. Remove everything and try again
                            retry = True
                            try:
                                del store[f'/{prefix}/{key}']
                            except IOError:
                                pass

                        if retry:
                            #Dataset too lareg to pass over http PUT
                            chunk_size = 500 #atoms in structure int(len(rec_arr)/4)
                            for i, start in enumerate(range(0, len(dset), chunk_size)):
                                small_data = dset[start:start+chunk_size]
                                if i==0:
                                    store.create_table(f'/{prefix}/{key}', dtype=dtype, shape=small_data.shape, data=small_data, **dset_parameters)
                                else:
                                    store[key].resize((store[f'/{prefix}/{key}'].shape[0] + small_data.shape[0]), axis=0)
                                    store[key][-small_data.shape[0]:] = small_data

def save_h5(hsds_file, h5_file, prefix="", hard_link_keys=["validation", "training", "test"], hard_link_map={"data_splits":"domains"}, hardlink_items=False):
    """A custom version of hsget to save a H5 file from a an h5 file on an hsds endpoint. Parallel writing does not work, so Toil is not used.

    Parameters:
    -----------
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
    add_keys = []
    hard_linked_keys = []
    with h5pyd.File(hsds_file, "r", use_cache=False) as store, h5py.File(h5_file, 'a') as f:
        for key in f[f'/{prefix}'].keys():
            if hardlink_items:
                #Convert to softlink
                old_kwd, new_kwd = next(hard_link_map.items())
                store[f"{prefix}/{key}"] = h5py.SoftLink(f"{prefix.split(old_kwd, 1)[0]}/{new_kwd}/{key}")
                continue
                
            hsds_object = store[f'/{prefix}/{key}']

            if isinstance(hsds_object, h5pyd.Group):
                h5_object = f.require_group(f'/{prefix}/{key}')
                for attr, attr_value in hsds_object.attrs.items():
                    h5_object.attrs[attr] = attr_value

                if key not in hard_link_keys:
                    add_keys.append(f'/{prefix}/{key}')
                else:
                    hard_linked_keys.append(f'/{prefix}/{key}')

            elif isinstance(h5_object, h5py.Dataset):
                dset = hsds_object[:]
                dset_parameters = {k:getattr(h5_object, k) for k in ('shape', 'dtype', 'chunks', 'compression',
                  'compression_opts', 'scaleoffset', 'shuffle', 'fletcher32', 'fillvalue')}
                f.create_dataset(f'/{prefix}/{key}', data=dset, **dset_parameters)
    
        for key in add_keys:
            save_h5(hsds_file, h5_file, prefix=f"{prefix}/{key}",
                    hard_link_keys=hard_link_keys, hard_link_map=hard_link_map, hardlink_items=False)
        
        #Add hard links last
        for hard_linked_key in hard_linked_keys:
            save_h5(hsds_file, h5_file, prefix=f"{prefix}/{hard_linked_key}",      
                    hard_link_keys=hard_link_keys, hard_link_map=hard_link_map, hardlink_items=True)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

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