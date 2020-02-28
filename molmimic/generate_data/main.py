import os

from molmimic.util.iostore import IOStore
from molmimic.util.cath import run_cath_hierarchy
from molmimic.util.hdf import get_file, filter_hdf_chunks
from molmimic.util.toil import map_job

from molmimic.generate_data.prepare_protein import process_domain
from molmimic.generate_data.calculate_features import calculate_features

from molmimic.generate_data import data_stores

from toil.realtimeLogger import RealtimeLogger

import logging
logging.getLogger('boto3').setLevel(logging.WARNING)
logging.getLogger('botocore').setLevel(logging.WARNING)
logging.getLogger('s3transfer').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

def get_domain_structure_and_features(job, cath_domain, superfamily,
  cathFileStoreID, update_features=None, force=False):
    """1) Run Feagtures depends on prepared strucutre Structure"""
    RealtimeLogger.info("get_domain_structure_and_features Process domain "+cath_domain)

    if force or not data_stores.prepared_cath_structures.exists(
      "{}/{}.pdb".format(superfamily, cath_domain)):
        #Get Processed domain file first
        job.addChildJobFn(process_domain, cath_domain, superfamily,
            cathFileStoreID=cathFileStoreID)

    #Check if any feature files exist
    feat_files = ["{}/{}_{}".format(superfamily, cath_domain, ext) for ext in \
        ('atom.h5', 'residue.h5', 'edges.h5')]
    feats_exist = [data_stores.cath_features.exists(f) for f in feat_files]

    if force or update_features is not None or not all(feats_exist):
        #Calculate features Processed domain file
        job.addFollowOnJobFn(calculate_features, cath_domain, superfamily,
            update_features=update_features)

def process_superfamily(job, superfamily, cathFileStoreID, update_features=None,
  force=False, further_parallize=True):
    cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
    cathcode = superfamily.replace("/", ".")

    cath_domains = filter_hdf_chunks(
        cath_file,
        "table",
        columns=["cath_domain"],
        drop_duplicates=True,
        cathcode=cathcode)["cath_domain"].iloc[:2]

    if further_parallize:
        map_job(job, get_domain_structure_and_features, cath_domains,
            superfamily, cathFileStoreID, update_features=update_features,
            force=force)
    else:
        RealtimeLogger.info("Looping over domain")
        for domain in cath_domains:
            try:
                RealtimeLogger.info("Processing domain "+domain)
                get_domain_structure_and_features(job, domain, superfamily,
                    cathFileStoreID, update_features=update_features,
                    force=force)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                raise

def start_toil(job, cathFileStoreID, cathcode=None, update_features=None, force=False):
    RealtimeLogger.info("started")

    #Start CATH hiearchy
    run_cath_hierarchy(job, cathcode, process_superfamily, cathFileStoreID,
        update_features=update_features, force=force)

    #Build Interactome
    #job.addChildJobFn()

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument(
        "-c", "--cathcode",
        nargs='+',
        default=None)
    parser.add_argument(
        "--features",
        nargs="+",
        default=None
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False)
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    #options.clean = "always"
    options.targetTime = 1

    if options.cathcode is not None:
        options.cathcode = [c.split(".") for c in options.cathcode]

    sfam_file = os.path.abspath("cath.h5")
    if not os.path.isfile(sfam_file):
        store = IOStore.get("aws:us-east-1:molmimic-cath")
        store.read_input_file("cath-domain-description-file-small.h5", sfam_file)

    with Toil(options) as workflow:
        if not workflow.options.restart:
            cathFileStoreID = workflow.importFile("file://" + os.path.abspath(sfam_file))
            job = Job.wrapJobFn(start_toil, cathFileStoreID, cathcode=options.cathcode,
                update_features=options.features, force=options.force)
            workflow.start(job)
        else:
            workflow.restart()
