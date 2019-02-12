import os
import pandas as pd
from itertools import chain
import datetime
from toil.realtimeLogger import RealtimeLogger
from molmimic.generate_data.calculate_ddi import download_pdb
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.job_utils import map_job

import logging
logging.getLogger('boto').setLevel(logging.CRITICAL)
logging.getLogger('boto3').setLevel(logging.CRITICAL)
logging.getLogger('botocore').setLevel(logging.CRITICAL)
logging.getLogger('iotore').setLevel(logging.CRITICAL)
logging.getLogger().setLevel(logging.CRITICAL)

def compare_sfam(job, sfam, useExisting=False, observed=True):
    work_dir = job.fileStore.getLocalTempDir()
    out_store = IOStore.get("aws:us-east-1:molmimic-missing-structures")
    inf_store = IOStore.get("aws:us-east-1:molmimic-interfaces")
    struc_store = IOStore.get("aws:us-east-1:molmimic-full-structures")

    all_missing = "missing_{}.h5".format("observed" if observed else "inferred")
    all_missing_f = os.path.join(work_dir, all_missing)

    if not useExisting or not out_store.exists(all_missing):
        obs_key = "{sfam}/{sfam}.{type}_interactome".format(sfam=sfam, type="observed" if observed else "inferred")
        obs_f = os.path.join(work_dir, os.path.basename(obs_key))

        try:
            inf_store.read_input_file(obs_key, obs_f)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            job.log("Unable to open {} ({})".format(obs_key, e))
            return []

        st = pd.HDFStore(obs_f)
        df = st.get("/table")
        st.close()

        get_key = lambda f, p, c, s, d: "{}/{}/{}_{}_sdi{}_d{}.pdb".format(int(f),
            p[1:3].lower(), p.upper(), c, s, d)

        mol_ints = df[["mol_pdb", "mol_chain", "mol_sdi_id", "mol_domNo", "mol_superfam_id"]]
        mol_ints = mol_ints.rename(columns={"mol_pdb":"pdb", "mol_chain":"chain",
            "mol_sdi_id":"sdi", "mol_domNo":"domNo", "mol_superfam_id":"sfam_id"})
        int_ints = df[["int_pdb", "int_chain", "int_sdi_id", "int_domNo", "int_superfam_id"]]
        int_ints = int_ints.rename(columns={"int_pdb":"pdb", "int_chain":"chain",
            "int_sdi_id":"sdi", "int_domNo":"domNo", "int_superfam_id":"sfam_id"})
        pdbs = pd.concat((mol_ints, int_ints)).drop_duplicates()

    else:
        out_store.read_input_file(all_missing, all_missing_f)
        sfams = pd.read_hdf(all_missing_f, "table")
        ibis_store = IOStore.get("aws:us-east-1:molmimic-ibis")
        pdbf = ibis_store.read_input_file("PDB.h5", os.path.join(work_dir, "PDB.h5"))
        pdb = pd.read_hdf(os.path.join(work_dir, "PDB.h5"), "merged", columns=["sdi", "sfam_id"]).drop_duplicates()
        pdbs = pd.merge(sfams, pdb, on="sdi")

    missing = []

    for i, row in pdbs.iterrows():
        try:
            if not struc_store.exists(get_key(row.sfam_id, row.pdb, row.chain, row.sdi, row.domNo)):
                raise
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            missing.append(row)
            continue
        RealtimeLogger.info("    Found - {} {} {} {} {}".format(row.sfam_id, row.pdb, row.chain, row.sdi, row.domNo))

    if len(missing)>0 and not useExisting:
        RealtimeLogger.info("{} Missing {} entries".format(sfam, len(missing)))

        missing = pd.DataFrame(missing)

        key = "{}_{}.h5".format("observed" if observed else "inferred", int(sfam))
        path = os.path.join(work_dir, key)
        missing.to_hdf(path, "table")
        out_store.write_output_file(path, key)
    elif len(missing)>0 and useExisting:
        missing = pd.DataFrame(missing)
        file = "missing_{}_{}.h5".format("observed" if observed else "inferred",
            str(datetime.datetime.now()).replace(" ", "-"))
        outfile = os.path.join(work_dir, file)
        missing.to_hdf(outfile, "table")
        out_store.write_output_file(outfile, file)


def get_missing(job, observed=True):
    import datetime
    work_dir = job.fileStore.getLocalTempDir()

    file = "missing_{}_{}.h5".format("observed" if observed else "inferred",
        str(datetime.datetime.now()).replace(" ", "-"))
    outfile = os.path.join(work_dir, file)

    store = IOStore.get("aws:us-east-1:molmimic-missing-structures")

    to_remove = []
    for k in store.list_input_directory("observed" if observed else "inferred"):
        path = os.path.join(work_dir, k)
        store.read_input_file(k, path)
        df = pd.read_hdf(path, "table", mode="r")
        df.to_hdf(outfile, "table", table=True, format="table", append=True,
            mode="a", complib="bzip2", complevel=9)
        to_remove.append(k)
        try:
            os.remove(path)
        except OSError:
            pass


    if not os.path.isfile(outfile):
        df = pd.DataFrame()
        df.to_hdf(outfile, "table", table=True, format="table", append=True,
            mode="a", complib="bzip2", complevel=9)

    store.write_output_file(outfile, file)

    try:
        os.remove(outfile)
    except OSError:
        pass

    for f in to_remove:
        try:
            store.remove(k)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            pass

def compare_sfams(job, useExisting=False, observed=True):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-missing-structures")
    all_missing = "missing_{}.h5".format("observed" if observed else "inferred")
    all_missing_f = os.path.join(work_dir, all_missing)

    if not useExisting or not store.exists(all_missing):
        inf_store = IOStore.get("aws:us-east-1:molmimic-interfaces")
        ending = ".{}_interactome".format("observed" if observed else "inferred")
        sfams = [k.split("/",1)[0] for k in inf_store.list_input_directory() \
            if k.endswith(ending)]
    else:
        store.read_input_file(all_missing, all_missing_f)
        sfams = pd.read_hdf(all_missing_f, "table", columns=["sfam"])["sfams"].drop_duplicates()

    map_job(job, compare_sfam, sfams)

    job.addFollowOnJobFn(get_missing, observed=observed)

def start_toil(job, useExisting=False):
    if not useExisting:
        job.addChildJobFn(compare_sfams, useExisting=useExisting)
        #job.addFollowOnJobFn(compare_sfams, useExisting=useExisting, observed=False)
    else:
        job.addChildJobFn(compare_sfam, None, useExisting=True)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--useExisting", default=False, action="store_true")
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"

    job = Job.wrapJobFn(start_toil, useExisting=options.useExisting)
    with Toil(options) as toil:
        toil.start(job)
