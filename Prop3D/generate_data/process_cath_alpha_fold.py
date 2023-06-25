import os
from pathlib import Path
import h5pyd

from Prop3D.generate_data.data_stores import data_stores
from Prop3D.parsers.foldseek import EasyCluster
from Prop3D.util.pdb import get_first_chain

from Prop3D.util.toil import map_job

def run_custom_hierarchy(job, prefix, func, full_h5_path, *args, **kwds):
    with h5py.File(full_h5_path, mode="a", use_cache=False, retries=100) as store:
        keys = [f"{prefix}/{k}" for k in store[prefix].keys()]
        if "domains" in store[keys[0]]:
            #It's a leaf node with domains
            map_job(job, func, keys, full_h5_path, *args, **kwds)
        else:
            #Continue down hierarchy
            map_job(job, run_custom_hierarchy, keys, func, full_h5_path, *args, **kwds)

def split_superfamily_at_level(job, cath_full_h5, superfamily, input_store_name, seq_id=.35, work_dir=None, 
  split_size={"train":0.8, "validation":0.1, "test":0.1}, work_dir=None):
    from Prop3D.generate_data.set_cath_h5_toil import split_superfamily_at_level as split_superfamily_at_level_main
    
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    structure_dir = Path(work_dir) / "structures"
    structure_dir.mkdirs(exist_ok=True)


    input_store = data_stores(job).custom_input(input_store_name)
    input_store.download_input_directory(superfamily, work_dir) 

    clusterer = EasyCluster(job=job, wokr_dir=work_dir)
    clusters_df = clusterer.cluster(str(structure_dir.relative_to(work_dir)), min_seq_id=seq_id)
    clusters_df = clusters_df.rename(columns={"structure":"cath_domain"})
    
    if seq_id==.35:
        cluster_representatives = clusters_df["cluster_representatives"].drop_duplicates()
        key = f"{superfamily}/representatives"
        missing_domains = []
        with h5pyd.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
            group = store.require_group(key)

            for domain in cluster_representatives:
                try:
                    group[domain] = store[f"{superfamily}/domains/{domain}'"]
                except KeyError:
                    missing_domains.append(domain)

            if len(missing_domains) > 0:
                store[key].attrs["missing_domains"] = missing_domains
            store[key].attrs["total_domains"] = len(representatives)


    return split_superfamily_at_level_main(job, cath_full_h5, superfamily, clusters_df, "cluster_representatives", 
        f"S{seq_id*100}", split_size)

def create_splits_for_superfamily_levels(job, superfamily, cath_full_h5, input_store_name):
    RealtimeLogger.info(f"Start splits for {superfamily}")

    with h5pyd.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group(f"{superfamily}/data_splits")

    for seq_id in [.35, .60, .95, 1]:
        job.addChildJobFn(split_superfamily_at_level, cath_full_h5, superfamily, seq_id=seq_id)

def setup_custom_cath_file_for_sfam(job, group_data, cath_full_h5, input_store_name):
    prefix, group_df = group_data
    prefix = "/".join(prefix)

    with h5py.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
        group = store.require_group(prefix)

        if len(prefix) == 4:
            store.require_group(f"{prefix}/domains")
            for _, row in group_df.itertuples():
                store.require_group(f"{row.prefix}/domains/{row.cath_domain}")
            
            #Create splits
            job.addFollowOnJobFn(create_splits_for_superfamily_levels, prefix, cath_full_h5, input_store_name)
    
    if len(prefix) < 4:
        group_names = ["C", "A", "T", "H"]
        groups = group_df.groupby(group_names[:len(prefix)+1])
        map_job(job, setup_custom_cath_file_for_sfam, groups, cath_full_h5, input_store_name)

def create_custom_hierarchy(job, prefix, cath_full_h5, input_store_name):
    input_store = data_stores(job).custom_input(input_store_name)
    children = list(input_store.list_input_directory(prefix))
    
    assert len(children) > 0, "Error somewhere"

    if all([input_store.get_number_of_items(c)==0 for c in children]):
        #Is leaf directory, e.g. one superfamily

        children = [Path(f) for f in children]
        if prefix=="" and all([c.stem.count(".")==3 for c in children]):
            #Are all CATH directories
            names = pd.DataFrame([(*f.parent.name.split(), f.stem) for f in child_files], columns=["C", "A", "T", "H", "cath_domain"])
            return map_job(job, setup_custom_cath_file_for_sfam, names.groupby("C"), cath_full_h5, input_store_name)

        else:
            with h5py.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
                store.require_group(f"{prefix}/domains")
                for c in children:
                    store.require_group(f"{prefix}/domains/{c}")
            return job.addFollowOnJobFn(create_splits_for_superfamily_levels, prefix, cath_full_h5, input_store_name)


    children = [f"{prefix}/{c}" for c in children]

    with h5py.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
        for c in children:
            store.require_group(c)
    
    map_job(job, create_hierarchy, children, input_store_name)

def preprocess_files(jobStore, input_data):
    #Run in home direcotry before toil starts

    def create_new_store(files):
        input_store_name = jobStore.split(":")[-1]+"--input"
        input_store = data_stores(jobStore).custom_input(input_store_name)
        for input_file in files:
            input_store.write_output_file(input_file, Path(input_file).name)
        return input_store_name

    if isinstance(input_data, (list, tuple)):
        #List of files, create new input bucket
        input_store_name = create_new_store(input_data)
    elif isinstance(input_data, str):
        input_data = Path(input_data)
        if input_data.is_dir():
            input_store = data_stores(jobStore).custom_input(str(input_data))
            input_store_name = str(input_data)
            if input_store.use_s3:
                input_store.write_output_file(str(input_path), "/")
        elif input_data.is_file():
            with input_data.open() as f:
                line = next(f)
                
            if Path(line).is_file():
                #A file with list of paths on enw lines
                with input_data.open() as f:
                    input_store_name = create_new_store([l.rstrip() for l in f])
            elif get_first_chain(str(input_data)) is not None:
                #Single PDB file
                input_store_name = create_new_store([input_data])
            else:
                raise RuntimeError("Invalid pdb input. Must be a list of pdb files, a direcotry of pdb files (nested is OK), or a single with a list of pdb files one each lines")
        else:
            #Check if file stoe with name already exists
            try:
                data_stores(jobStore).custom_input(str(input_data), create=False)
                input_store_name = str(input_data)
            except RuntimeError:
                try:
                    data_stores(jobStore).custom_input(str(input_data)+"--input", create=False)
                    input_store_name = str(input_data)+"--input"
                except RuntimeError:
                    raise RuntimeError("Invalid pdb input. Must be a list of pdb files, a direcotry of pdb files (nested is OK), or a single with a list of pdb files one each lines")

    else:
        raise RuntimeError("Invalid pdb input. Must be a list of pdb files, a direcotry of pdb files (nested is OK), or a single with a list of pdb files one each lines")

    return input_store_name