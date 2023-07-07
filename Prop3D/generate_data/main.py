import os
from datetime import datetime
from collections import defaultdict

import h5pyd
import pandas as pd
from toil.realtimeLogger import RealtimeLogger

from Prop3D.util import safe_remove, str2boolorlist
from Prop3D.util.iostore import IOStore
from Prop3D.util.cath import run_cath_hierarchy, run_cath_hierarchy_h5
from Prop3D.util.hdf import get_file, filter_hdf_chunks
from Prop3D.util.toil import map_job, map_job_follow_ons
from Prop3D.util.pdb import get_atom_lines

from Prop3D.generate_data.prepare_protein import process_domain
from Prop3D.generate_data.calculate_features import calculate_features
from Prop3D.generate_data.calculate_features_hsds import calculate_features as calculate_features_hsds
from Prop3D.generate_data.set_cath_h5_toil import create_h5_hierarchy

from Prop3D.generate_data.data_stores import data_stores



import logging
logging.getLogger('boto3').setLevel(logging.WARNING)
logging.getLogger('botocore').setLevel(logging.WARNING)
logging.getLogger('s3transfer').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

def get_domain_structure_and_features(job, cath_domain, superfamily,
  cathFileStoreID, update_features=None, further_parallelize=False, force=False, use_hsds=True):
    """1) Run Features depends on prepared strucutre Structure"""
    RealtimeLogger.info("get_domain_structure_and_features Process domain "+cath_domain)

    if use_hsds:
        calc_features_func = calculate_features_hsds
    else:
        calc_features_func = calculate_features

    if superfamily is not None:
        key = "{}/{}.pdb".format(superfamily, cath_domain)
        should_prepare_structure = not data_stores(job).prepared_cath_structures.exists(key) or \
            data_stores(job).prepared_cath_structures.get_size(key)==0
        h5_key = f"/{superfamily}/domains/{cath_domain}"
    else:
        should_prepare_structure = True
        further_parallelize = False
        force=True
        h5_key = os.path.basename(os.path.basename(cath_domain))[0]

    if force or (isinstance(force, int) and force==2) or should_prepare_structure:
        #Get Processed domain file first
        RealtimeLogger.info(f"Processing domain {cath_domain}")
        if further_parallelize:
            job.addChildJobFn(process_domain, cath_domain, superfamily,
                cathFileStoreID=cathFileStoreID)
        else:
            try:
                prepared_file, _, _, local_file = process_domain(job, cath_domain, superfamily, cathFileStoreID=cathFileStoreID)
                local_domain_file = prepared_file if local_file else None
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                #Failed, do not proceed in calculating features
                raise
    # else:
    #     assert 0, (should_prepare_structure, force)

    if not force and not use_hsds:
        #Check if any feature files exist
        feat_files = ["{}/{}_{}".format(superfamily, cath_domain, ext) for ext in \
            ('atom.h5', 'residue.h5', 'edges.h5')]
        feats_exist = all([data_stores(job).cath_features.exists(f) for f in feat_files])
    elif not force and use_hsds:
        with h5pyd.File(cathFileStoreID, mode="r", use_cache=False, retries=100) as store:
            try:
                feat_files = list(store[f"/{superfamily}/domains/{cath_domain}"].keys())
                if len(feat_files) == 3 and "atom" in feat_files and \
                  "residue" in feat_files and feat_files and "edges":
                    feats_exist = True
                else:
                    feats_exist = False
            except KeyError:
                feats_exist = False
    else:
        feats_exist = False

    if not should_prepare_structure and not feats_exist:
        local_domain_file = None #Will download in features function

    RealtimeLogger.info(f"get_domain_structure_and_features FEATURES {force} {update_features}, {feats_exist}")

    if force or update_features is not None or not feats_exist:
        #Calculate features Processed domain file
        if further_parallelize:
            job.addFollowOnJobFn(calc_features_func, cathFileStoreID, cath_domain,
                superfamily, update_features=update_features, domain_file=local_domain_file)
        else:
            RealtimeLogger.info("get_domain_structure_and_features calculate_features")
            try:
                calc_features_func(job, cathFileStoreID, cath_domain, superfamily,
                    update_features=update_features, domain_file=local_domain_file)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                #Failed
                raise

    if update_features is not None:
        jobStoreName = os.path.basename(job.fileStore.jobStore.config.jobStore.split(":")[-1])
        done_file = job.fileStore.getLocalTempFile()
        data_stores(job).data_eppic_cath_features.write_output_file(done_file,
            f"updates/{jobStoreName}/{superfamily}/{cath_domain}")
        safe_remove(done_file)

def process_superfamily(job, superfamily, cathFileStoreID, update_features=None,
  force=False, use_hsds=True, further_parallize=True):

    cathcode = superfamily.replace("/", ".")
    if not use_hsds:
        cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
        cath_domains = filter_hdf_chunks(
            cath_file,
            "table",
            columns=["cath_domain"],
            drop_duplicates=True,
            cathcode=cathcode)["cath_domain"].tolist()
    else:
        with h5pyd.File(cathFileStoreID, mode="a", use_cache=False) as store:
            cath_domains = list(store[f"{superfamily}/domains"].keys())
        cath_domains = cath_domains

    if not force and update_features is None:
        if use_hsds:
            try:
                RealtimeLogger.info("Finding completed domains...")
                # with h5pyd.File(cathFileStoreID, mode="a", use_cache=False, retries=100) as store:
                #     processed_domains = defaultdict(int)
                #     def count_feature_types(name):
                #         if name.count("/") == 6:
                #             #Name looks like 1/2/3/4/domains/xx/atom
                #             domain = name.split("/")[-2]
                #             processed_domains[domain] += 1
                #     store[superfamily].visit(count_feature_types)
                #     done_domains = [domain for domain, count in processed_domains.items() if count==3]
                done_domains = []
            except KeyError as e:
                raise RuntimeError(f"Must create hsds file first. Key not found {e}")
        else:
            assert 0
        cath_domains = list(set(cath_domains)-set(done_domains))

    if False and not force and update_features is None:
        #Get domians that have edge features uploaded (last feature file to be uploaded so we know its done)
        done_domains = [os.path.basename(domain).split("_")[0] for domain in \
            data_stores(job).cath_features.list_input_directory(superfamily) \
            if domain.endswith("edges.txt.gz")]
        n_domains = len(cath_domains)
        cath_domains = list(set(cath_domains)-set(done_domains))
        RealtimeLogger.info("Running {}/{} domains from {}".format(len(cath_domains), n_domains, cathcode))
    elif False and update_features is not None:
        #oldest = datetime.date(2021, 2, 9)
        jobStoreName = os.path.basename(job.fileStore.jobStore.config.jobStore.split(":")[-1])
        updated_domains = [cath_domain for cath_domain, last_modified in \
            data_stores(job).data_eppic_cath_features.list_input_directory(
                f"updates/{jobStoreName}/{superfamily}/", with_times=True)] # \
                #if datetime.strptime(last_modified, '%Y-%m-%dT%H:%M:%S')<oldest]

        cath_domains = list(set(cath_domains)-set(updated_domains))
    else:
        RealtimeLogger.info("Running {} domains from {}".format(len(cath_domains), cathcode))

    if further_parallize:
        map_job(job, get_domain_structure_and_features, cath_domains,
            superfamily, cathFileStoreID, update_features=update_features,
            further_parallelize=False, use_hsds=use_hsds, force=force)
    else:
        RealtimeLogger.info("Looping over domain")
        for domain in cath_domains:
            try:
                RealtimeLogger.info("Processing domain "+domain)
                get_domain_structure_and_features(job, domain, superfamily,
                    cathFileStoreID, update_features=update_features,
                    further_parallelize=False, use_hsds=use_hsds, force=force)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                pass

def start_domain_and_features(job, cathFileStoreID, cathcode=None, skip_cathcode=None,
  update_features=None, use_hsds=True, work_dir=None, pdbs=None, force=False):
    if not use_hsds:
        cath_hierarchy_runner = run_cath_hierarchy
    else:
        cath_hierarchy_runner = run_cath_hierarchy_h5

    done_domains = None
    if cathcode is not None:
        if not isinstance(cathcode, (list, tuple)):
            cathcode = [cathcode]
        RealtimeLogger.info("Running sfams: {}".format(cathcode))

        fixed_sfams = ["/".join(sfam) if isinstance(sfam, (list,tuple)) else sfam.replace(".", "/") \
            for sfam in cathcode]
        num_sfams_to_run = len(fixed_sfams)

        done_domains = None
        if update_features is None and not force:
            if not use_hsds:
                cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)

                cath_domains = filter_hdf_chunks(
                    cath_file,
                    "table",
                    columns=["cath_domain", "cathcode"],
                    drop_duplicates=True)

                done_domains = [os.path.basename(domain).split("_")[0] for sfam in cathcode \
                    for domain in data_stores(job).cath_features.list_input_directory("/".join(sfam) if \
                        isinstance(sfam, (list,tuple)) else sfam.replace(".", "/")) \
                        if domain.endswith("edges.txt.gz")]

            else:
                with h5pyd.File(cathFileStoreID, mode="r", use_cache=False) as store:
                    try:
                        cath_domains = pd.DataFrame([(cath_domain, sfam) for sfam in fixed_sfams \
                            for cath_domain in store[f"{sfam}/domains"].keys()], columns=["cath_domain", "cathcode"])
                    except KeyError:
                        raise RuntimeError(f"Must create hsds file first. Key not found {sfam+'/domain' for sfam in fixed_sfams}")

                    # try:
                    #     check = [f"{sfam}/domains" for sfam in fixed_sfams]
                    #     RealtimeLogger.info(f"Checking {check}")
                    #     done_domains = [domain for sfam in fixed_sfams for domain in \
                    #         store[f"{sfam}/domains"].keys() if \
                    #         len(store[f"{sfam}/domains/{domain}"].keys())==3]
                    # except KeyError as e:
                    #     raise RuntimeError(f"Must create hsds file first. Key not found {e}")

            num_domains_to_run = 500000
            #domains_to_run = cath_domains[~cath_domains["cath_domain"].isin(done_domains)]
            #num_domains_to_run = len(domains_to_run)
            num_sfams_to_run = 7000 #len(cath_domains["cathcode"].drop_duplicates())
        else:
            #No update and forced
            RealtimeLogger.info(f"Counting # of domains from {fixed_sfams}")
            with h5pyd.File(cathFileStoreID, mode="a", use_cache=False) as store:
                try:
                    cath_domains = sum(1 for sfam in fixed_sfams for cath_domain in \
                        store[f"{sfam}/domains"].keys())
                except KeyError:
                    raise RuntiemError("Must create hsds file first")
    elif pdbs is not None:
        assert use_hsds
        with h5pyd.File(cathFileStoreID, mode="r", use_cache=False) as store:
            all_domains = list(store["domains"].keys()) 
            try:
                done_domains = [k for k in all_domains if len(store[f"domains/{k}"].keys())==3]
            except KeyError:
                raise RuntimeError(f"Must create hsds file first.")
            
        if not force: # and (isinstance(pdbs, bool) and pdbs):
            domains_to_run = list(set(all_domains)-set(done_domains))
        else:
            domains_to_run = all_domains
        # elif isinstance(pdbs, (list, tuple)):

        #     domains_to_run = [f for f in pdbs if os.path.splitext(os.path.basename(f))[0] \
        #     not in done_domains]
        RealtimeLogger.info(f"Running domains: {domains_to_run}")
        map_job(job, get_domain_structure_and_features, domains_to_run, None, cathFileStoreID,
            update_features=update_features, force=force, further_parallelize=True, use_hsds=use_hsds)
        return
    else:
        #No cath code proved, use all
        if update_features is None and not force:
            num_sfams_to_run = "all"
            num_domains_to_run = 500000
        else:
            #FIX: find all completed domains?
            num_sfams_to_run = "all"
            num_domains_to_run = 500000

    if update_features is None and done_domains is not None:
        domains_to_run = cath_domains[~cath_domains["cath_domain"].isin(done_domains)]
        num_domains_to_run = len(domains_to_run)
        num_sfams_to_run = len(cath_domains["cathcode"].drop_duplicates())
    else:
        num_domains_to_run = 500000

    RealtimeLogger.info("Domains to from {} superfamilies to run: {}".format(
        num_sfams_to_run, num_domains_to_run))

    if num_domains_to_run > 500:
        #Start CATH hiearchy
        RealtimeLogger.info("Starting CATH Hierachy")
        cath_hierarchy_runner(job, cathcode, process_superfamily, cathFileStoreID,
            skip_cathcode=skip_cathcode, update_features=update_features, use_hsds=use_hsds, force=force)
    else:
        superfamilies = domains_to_run["cathcode"].drop_duplicates().str.replace(".", "/")
        if skip_cathcode is not None and len(skip_cathcode) > 0:
            superfamilies = superfamilies[~superfamilies.isin(skip_cathcode)]
        RealtimeLogger.info("Superfamilies to run: {}".format(len(superfamilies)))
        map_job(job, process_superfamily, superfamilies, cathFileStoreID,
            update_features=update_features, force=force)

    #Build Interactome
    #job.addChildJobFn()

def start_toil(job, cathFileStoreID, cathcode=None, skip_cathcode=None, pdbs=None, update_features=None, use_hsds=True, work_dir=None, force=False):
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    if pdbs is not None and (cathcode is not None or skip_cathcode is not None):
        raise RuntimeError("Cannot use --pdbs with --cathcode or --skip_cathcode")

    if use_hsds:
        job.addChildJobFn(create_h5_hierarchy, cathFileStoreID, cathcode=cathcode,
            skip_cathcode=skip_cathcode, pdbs=pdbs, work_dir=work_dir, force=force)

    job.addFollowOnJobFn(start_domain_and_features, cathFileStoreID, cathcode=cathcode,
        skip_cathcode=skip_cathcode, update_features=update_features, use_hsds=use_hsds,
        pdbs=pdbs, force=force)

def str2boolorval(v):
    def is_num(a):
        try:
            int(a)
            return True
        except ValueError:
            return False

    if isinstance(v, bool):
       return v
    elif is_num(v):
        return int(v)
    elif v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        return False
    
def start_workflow(options):
    with Toil(options) as workflow:
        if not workflow.options.restart:
            if options.no_hsds:
                cathFileStoreID = workflow.importFile("file://" + os.path.abspath(sfam_file))
            else:
                cathFileStoreID = options.hsds_file
            job = Job.wrapJobFn(start_toil, cathFileStoreID, cathcode=options.cathcode,
                skip_cathcode=options.skip_cathcode, pdbs=options.pdb, update_features=options.features,
                use_hsds=not options.no_hsds, work_dir=options.work_dir, force=options.force)
            workflow.start(job)
        else:
            workflow.restart()

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job
    from toil.jobStores.abstractJobStore import NoSuchJobStoreException

    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument(
        "-c", "--cathcode",
        nargs='+',
        default=None)
    parser.add_argument(
        "-s", "--skip_cathcode",
        nargs='+',
        default=None)
    #parser.add_argument("--pdb", default=None, nargs="+")
    parser.add_argument(
        '--pdb', 
        type=str2boolorlist, 
        nargs='*', 
        default=None)
    parser.add_argument(
        "--features",
        nargs="+",
        default=None
    )
    parser.add_argument(
        "--force",
        type=str2boolorval,
        nargs='?',
        const=True,
        default=False,
        help="Which files should be overwrriten. (0, False)=No files overwritten; (1, True)=Only Features, 2=Structure and Features, 3=Entire database. Default is False.")
    parser.add_argument(
        "--hsds_file",
        default=None,
        help="Path to h5 file on the HSDS endpoint")
    parser.add_argument(
        "--no_hsds",
        action="store_true",
        default=False)
    parser.add_argument(
        "--work_dir",
        default=os.getcwd())
    parser.add_argument(
        "--restartable",
        action="store_true",
        default=False,
        help="Create a new workflow or restart if it already exists")
    
    parser.add_argument("--hs_username", default=None, help="HSDS username. If not provided it will use hsinfo")
    parser.add_argument("--hs_password", default=None, help="HSDS username. If not provided it will use hsinfo")
    parser.add_argument("--hs_endpoint", default=None, help="HSDS username. If not provided it will use hsinfo")

    options = parser.parse_args()

    options.logLevel = "DEBUG"
    #options.clean = "always"
    options.targetTime = 1

    if options.pdb is not None and (options.cathcode is not None or options.skip_cathcode is not None):
        raise RuntimeError("Cannot use --pdbs with --cathcode or --skip_cathcode")

    if options.cathcode is not None:
        options.cathcode = [c.split(".") for c in options.cathcode]

    if options.skip_cathcode is not None:
        options.skip_cathcode = [c.replace(".", "/") for c in options.skip_cathcode if c.count(".")==3]

    if options.pdb is not None:
        if isinstance(options.pdb, str):
            options.pdb = [options.pdb]
        if isinstance(options.pdb, (list, tuple)) and len(options.pdb)==0:
            #Empty -> use whole pdb
            options.pdb = True
        if isinstance(options.pdb, (list, tuple)) and len(options.pdb) == 1:
            if os.path.isfile(options.pdb[0]):
                for _ in get_atom_lines(options.pdb[0]):
                    break
                else:
                    #Not a pdb, list of files or IDs
                    with open(options.pdb[0]) as f:
                        options.pdb = [pdb.rstrip() for pdb in f]
            elif os.path.isdir(options.pdb[0]):
                options.pdb = [pdb for pdb in next(os.walk(options.pdb[0]))[2] if pdb.endswith(".pdb")]
            else:
                raise RuntimeError("Invalid option for --pdb. Must be paths to pdbs files as arguments, a single file with path on each line, or a directory with pdb files")

    if options.no_hsds:
        sfam_file = os.path.abspath("cath.h5")
        if not os.path.isfile(sfam_file):
            store = IOStore.get("aws:us-east-1:Prop3D-cath")
            store.read_input_file("cath-domain-description-file-small.h5", sfam_file)
    elif options.hsds_file is None:
        raise RuntimeError("Must specify hsds file")
    else:
        #Copy HSDS info so other workers can access it
        from h5pyd import Config
        config = Config()
        if options.hs_username is None:
            options.hs_username = config["hs_username"]
        if options.hs_password is None:
            options.hs_password = config["hs_password"]
        if options.hs_endpoint is None:
            options.hs_endpoint = config["hs_endpoint"]

        #Set as
    # elif "HS_USERNAME" not in os.environ:
    #     raise RuntimeError("Must specify HSDS username env variable: HS_USERNAME")
    # elif "HS_PASSWORD" not in os.environ:
    #     raise RuntimeError("Must specify HSDS username env variable: HS_PASSWORD")
    # elif "HS_ENDPOINT" not in os.environ or not os.environ["HS_ENDPOINT"].startswith("http"):
    #     raise RuntimeError("Must specify HSDS endpoint env variable: HS_ENDPOINT and it must begin with http")

    os.environ["TOIL_START_DIR"] = options.work_dir

    if options.restartable:
        options.restart = True
        try:
            start_workflow(options)
        except NoSuchJobStoreException:
            options.restart = False
            start_workflow(options)

    else:
        start_workflow(options)

    
