import os, sys
import glob
from multiprocessing.pool import ThreadPool
import shutil
import itertools as it

import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS

def process_inferred_interaction(job, inf_int_id, mol_sfam, table, tableInfStoreID, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    tableInfFile = get_file(job, "IBIS_inferred_{}.h5".format(table), tableInfStoreID)

    inf_int_id_file = os.path.basename(job.fileStore.getLocalTempFileName())

    try:
        inf_int_id_file = "_".join(map(str, inf_int_id))
        inf_int_id, nbr_obs_int_id, nbr_sdi, mol_sdi = inf_int_id
        inferred_interfaces = filter_hdf(tableInfFile, "Intrac{}".format(table),
            nbr_obs_int_id = nbr_obs_int_id,
            nbr_sdi_id = nbr_sdi,
            mol_sdi_id = mol_sdi)
        print inf_int_id_file, inferred_interfaces
        job.log("{} {}".format(inf_int_id_file, inferred_interfaces))
        inferred_interfaces = inferred_interfaces.iloc[inf_int_id].to_frame().T
        job.log("FRAME: {}".format(inferred_interfaces))
        print "FRAME: {}".format(inferred_interfaces)
        inferred_interfaces["mol_sdi_id"] = inferred_interfaces["mol_sdi_id"].astype(int)
        #inferred_interfaces = pd.read_hdf([unicode(tableInfFile)], "/table", where=[
        #    "nbr_obs_int_id=={}".format(nbr_obs_int_id),
        #    "nbr_sdi=={}".format(nbr_sdi),
        #    "mol_sdi=={}".format(mol_sdi),
        #    "int_sdi=={}".format(int_sdi),
        #    ])

        if inferred_interfaces.shape[0] == 0:
            return

        pdb_file = get_file(job, "PDB.h5", pdbFileStoreID)
        struct_domains = filter_hdf(pdb_file, "StructuralDomains", "sdi", mol_sdi)

        #Add PDB, chain, and sdi information. Inner join to only allow sdi's that are in databses, not obsolete
        inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="inner", left_on="mol_sdi_id", right_on="sdi")

        obspath = job.fileStore.getLocalTempFileName()
        out_store.read_input_file("{0}/{0}.observed_interactome".format(mol_sfam), obspath)

        try:
            observed_interactome = filter_hdf(obspath, "table", "obs_int_id", nbr_obs_int_id)
        except TypeError:
            job.log("Failed reading {}".format(obspath))
            raise

        try:
            os.remove(obspath)
        except OSError:
            pass

        observed_interactome = observed_interactome.rename(columns={
            'mol_domNo': 'nbr_domNo',
            #'mol_gi':'nbr_gi',
            'mol_pdb': 'nbr_pdb',
            'mol_chain': 'nbr_chain'})

        #Add in neghibor information from observed interactome
        inferred_interfaces["nbr_obs_int_id"] = inferred_interfaces["nbr_obs_int_id"].astype(int)
        inferred_interfaces = pd.merge(inferred_interfaces, observed_interactome,
            how="left", left_on="nbr_obs_int_id", right_on="obs_int_id",
            suffixes=["_inf", "_obs"])
        del observed_interactome

        #Select relevant columns
        if "int_superfam_id_inf" in inferred_interfaces.columns:
            int_superfam_id_col = "int_superfam_id_inf"
        elif "int_superfam_id_x" in inferred_interfaces.columns:
            #Suffix not working?
            int_superfam_id_col = "int_superfam_id_x"
        else:
            raise RuntimeError("Merge faled for obs and inf")
        try:
            inferred_interfaces = inferred_interfaces[["sdi", "nbr_sdi_id", "nbr_score", "nbr_taxid",
                "nbr_superfam_id", int_superfam_id_col, "nbr_obs_int_id", "resn", "resi",
                "domNo", "pdbId", "chnLett", "from", "to", "int_sdi_id", "int_taxid", "int_res",
                "int_domNo", "int_pdb", "int_chain", "int_sdi_from", "int_sdi_to",
                "nbr_domNo", "nbr_pdb", "nbr_chain"]]
        except KeyError as e:
            job.log("Unable to filter df. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
            raise

        #Rename columns
        inferred_interfaces = inferred_interfaces.rename(columns={
            int_superfam_id_col:"int_superfam_id",
            "int_sdi_id":"int_sdi",
            "resn":"mol_resn",
            "sdi":"mol_sdi",
            "resi":"mol_resi",
            "domNo":"mol_domNo",
            "pdbId":"mol_pdb",
            "chnLett":"mol_chain",
            "from":"mol_sdi_from",
            "to":"mol_sdi_to"})

        try:
            inferred_interfaces = inferred_interfaces[ \
                ~inferred_interfaces["mol_sdi"].isnull() & \
                ~inferred_interfaces["int_sdi"].isnull()]
        except KeyError as e:
            job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
            raise

        try:
            resi = [decode_residues(job, row.pdbId, row.chnLett, row.resi, row) \
                for row in inferred_interfaces.itertuples()]
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            return None

        inferred_interfaces = inferred_interfaces.assign(resi=resi)

        df_file = job.fileStore.getLocalTempFileName()
        inferred_interfaces.to_hdf(unicode(df_file), "table", format="table",
            table=True, complevel=9, complib="bzip2", min_itemsize=1024,
            data_coumns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"],
            )
        job.log("Wrote "+df_file)

        #Add ibis info into out store
        out_store.write_output_file(df_file, "{}/_infrows/Intrac{}/{}.h5".format(int(mol_sfam), table, inf_int_id_file))
    
        for f in (df_file, pdb_file, tableInfFile, obspath):
            try:
                os.remove(f)
            except OSError:
                pass

    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        import traceback
        tb = traceback.format_exc()
        job.log("FAILED {} {}".format(mol_sfam, tb))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write(str(inf_int_id))
            f.write("\n")
            f.write(tb)
        out_store.write_output_file(fail_file, "{}/_infrows/Intrac{}/{}.failed".format(int(mol_sfam), table, inf_int_id_file))
        try:    
            os.remove(fail_file)
        except OSError:
             pass    

def merge_table_sfam(job, sfam_id, table):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    print "Start merge", sfam_id

    status = "inferred (table {}) {}".format(table, sfam_id)
    resi_prefix = "Intrac{}_{}.inferred_interactome".format(table, int(sfam_id))
    new_cols = ["mol_res"]
    data_cols = ["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"]

    resi_path = os.path.join(work_dir, resi_prefix)

    #Combine residues into dataframe
    possible_error = []
    for row_prefix in out_store.list_input_directory("{}/_infrows/Intrac{}".format(int(sfam_id), table)):
        job.log("Running  {} {}".format(sfam_id, row_prefix))

        row_file = os.path.join(work_dir, os.path.basename(row_prefix))
        out_store.read_input_file(row_prefix, row_file)

        df = pd.read_hdf(row_file, "table")
        try:
            df.to_hdf(unicode(resi_path), "table", table=True, format="table", mode="a",
                data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)
            job.fileStore.deleteGlobalFile(conv_store)
        except (SystemExit, KeybaordInterrupt):
            raise
        except:
            import traceback
            tb = traceback.format_exc()
            job.log("Failed writing {}: {} {}".format(sfam_id, resi_path, tb))
            possible_errors.append(tb)

        try:
            os.remove(row_file)
        except OSError:
            pass

        out_store.remove_file(row_prefix)

    if os.path.isfile(resi_path):
        #Upload to S3
        out_store.write_output_file(resi_path, "{}/_inftables/{}".format(int(sfam_id), resi_prefix))

        #Cleanup
        os.remove(resi_path)
        print "End merge", sfam_id
    else:
        job.log("Failed merging: {}".format(resi_path))
        print "Failed merging: {}".format(resi_path)
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write("no rows?\n")
            for e in possible_errors:
                f.write(e)
                f.write("\n")
        out_store.write_output_file(fail_file, "{}/_inftables/{}.failed".format(int(sfam_id), resi_prefix))


def get_table_sfams(job, mol_sfam_id, table, tableInfStoreID, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    tableInfPath = get_file(job, "IBIS_inferred_{}.h5".format(table), tableInfStoreID)

    #inf_job_ids = (row for pd.read_hdf([unicode(tableInfPath)], "/table", mode="r", \
    #    chunksize=100, columns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"], \
    #    where=["nbr_obs_int_id=={}".format(nbr_obs_int_id), \
    #           "nbr_sdi=={}".format(nbr_sdi), \
    #           "mol_sdi=={}".format(mol_sdi), \
    #           "int_sdi=={}".format(int_sdi), \
    #        ]) for row in df.itertuples(index=False))
    
    #def inf_job_ids():
    #    for df in pd.read_hdf(unicode(tableInfPath), "/table", mode="r", chunksize=100, \
    #      columns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"], \
    #      where=["nbr_obs_int_id=={}".format(nbr_obs_int_id), "nbr_sdi=={}".format(nbr_sdi), \
    #             "mol_sdi=={}".format(mol_sdi), "int_sdi=={}".format(int_sdi)]):
    #        for row in df.itertuples(index=False):
    #            yield row
    inf_int_ids = filter_hdf_chunks(tableInfPath, "Intrac{}".format(table), chunksize=100,
        columns=["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi_id"],
        nbr_superfam_id=mol_sfam_id).iloc[:1]
    inf_int_ids = [list(row) for row in inf_int_ids.itertuples()]
    print inf_int_ids
    map_job(job, process_inferred_interaction, inf_int_ids, mol_sfam_id, table, tableInfStoreID, pdbFileStoreID)

    job.addFollowOnJobFn(merge_table_sfam, mol_sfam_id, table)

    try:
        os.remove(tableInfPath)
    except OSError:
        pass

def get_inferred_structural_interactome_by_table(job, table, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    #Read in H5 for entire table
    tableInfPath = get_file(job, "IBIS_inferred_{}.h5".format(table), in_store)
    tableInfPathFileStoreID = job.fileStore.writeGlobalFile(tableInfPath)
    sfams = filter_hdf_chunks(tableInfPath, "Intrac{}".format(table), columns=["nbr_superfam_id"]).drop_duplicates().dropna().iloc[:1]
    #sfams = set(sfam for df in pd.read_hdf([unicode(tableInfPath)], "/table", chunksize=1000, mode="r", \
    #    columns=["mol_superfam_id"]) for sfam in df["mol_superfam_id"].drop_duplicates.dropna())
    map_job(job, get_table_sfams, sfams["nbr_superfam_id"].tolist(), table, tableInfPathFileStoreID, pdbFileStoreID)

    try:
        os.remove(tableInfPath)
    except OSError:
        pass

def merge_inferred_interactome_sfam(job, sfam_id):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    iostore = IOStore.get("{}:molmimic-ibis".format(prefix))

    sfam_prefix = "{}/_inftables".format(int(sfam_id), resi_prefix)

    sfam_file = "{s}/{s}.inferred_interactome".format(s=int(sfam_id))

    for table_prefix in iostore.list_input_directory(sfam_prefix):
        if table_prefix.endswith("failed"): continue
        table_file = os.path.join(work_dir, table_prefix)
        iostore.read_input_file(table_prefix, table_file)
        try:
            for df in pd.read_hdf(unicode(table_file), "table", chunksize=1000):
                df[cols] = df[cols].astype(float)
                df.to_hdf(unicode(sfam_id), "table", mode="a", format="table",
                    table=True, complevel=9, complib="bzip2", min_itemsize=1024)
        except (IOError, ValueError):
            job.log("Failed to read {}".format(table_file))

        try:
            os.remove(table_file)
        except OSError:
            pass

        iostore.remove_file(table_prefix)

    out_store.write_output_file(resi_path, sfam_file)

def merge_inferred_interactome(job):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))
    pdb_path = job.fileStore.readGlobalFile(pdbFileStoreID)
    all_sfam = [os.path.basname(f).split(".") for f in out_store.list_input_directory() if not f.endswith("failed")]
    skip_sfam = [f[0] for f in all_sfam if f[1] == "inferred_interactome"]
    sfam_to_run = [f[0] for f in all_sfam if f[1] == "observed_interactome" \
        and f[0] not in skip_sfam]

    map_job(job, merge_inferred_interactome_sfam, sfam_to_run)

def start_toil(job):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    #Download PDB info
    pdb_file = os.path.join(work_dir, "PDB.h5")
    in_store.read_input_file("PDB.h5", pdb_file)

    #Add pdb info into local job store
    pdbFileStoreID = job.fileStore.writeGlobalFile(pdb_file)

    tables = [1] #xrange(1, 87)
    map_job(job, get_inferred_structural_interactome_by_table, tables, pdbFileStoreID)
    job.addFollowOnJobFn(merge_inferred_interactome)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    print "Running"

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
