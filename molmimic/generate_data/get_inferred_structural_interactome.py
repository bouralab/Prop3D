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

from toil.realtimeLogger import RealtimeLogger

def process_inferred_interaction(job, inf_int_id, nbr_sfam, table, tableInfStoreID, pdbFileStoreID, taxFileStoreID, isrow=False, work_dir=None):
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    inf_int_id_file = os.path.basename(job.fileStore.getLocalTempFileName())
    try:
        if not isrow:
            inf_int_id_file = "_".join(map(str, inf_int_id))
            inf_int_id, nbr_obs_int_id, nbr_sdi, mol_sdi = inf_int_id
            tableInfFile = get_file(job, "IBIS_inferred_{}.h5".format(table), tableInfStoreID)
            inferred_interfaces = filter_hdf_chunks(tableInfFile, "Intrac{}".format(table),
                nbr_obs_int_id = nbr_obs_int_id,
                nbr_sdi_id = nbr_sdi,
                mol_sdi_id = mol_sdi)
            inferred_interfaces = inferred_interfaces.loc[inf_int_id].to_frame().T
        else:
            inferred_interfaces = inf_int_id.to_frame().T
            #inf_int_id = inferred_interfaces[""]
            nbr_obs_int_id = inf_int_id["nbr_obs_int_id"]
            nbr_sdi = inf_int_id["nbr_sdi_id"]
            mol_sdi = inf_int_id["mol_sdi_id"]

        inferred_interfaces["mol_sdi_id"] = inferred_interfaces["mol_sdi_id"].astype(float)

        if inferred_interfaces.shape[0] == 0:
            return

        pdb_file = get_file(job, "PDB.h5", pdbFileStoreID) if not isrow else pdbFileStoreID
        tax_file = get_file(job, "pdb_chain_taxonomy.h5", taxFileStoreID) if not isrow else taxFileStoreID

        try:
            struct_domains = filter_hdf(
                pdb_file, "merged",
                columns = ["sdi", "domNo", "gi", "pdbId", "chnLett", "from", "to", "sfam_id"],
                sdi = mol_sdi).drop_duplicates()
        except RuntimeError:
            job.log("SDI {} is obsolete".format(mol_sdi))
            return

        #Add PDB, chain, and sdi information. Inner join to only allow sdi's that are in databses, not obsolete
        inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="inner", left_on="mol_sdi_id", right_on="sdi")

        if isrow:
            obspath = tableInfStoreID
        else:
            obspath = job.fileStore.getLocalTempFileName()
            out_store.read_input_file("{0}/{0}.observed_interactome".format(nbr_sfam), obspath)

        try:
            observed_interactome = filter_hdf_chunks(obspath, "table", "obs_int_id", nbr_obs_int_id)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            try:
                observed_interactome = filter_hdf_chunks(obspath, "table", "obs_int_id", float(nbr_obs_int_id))
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                job.log("Failed reading {}".format(obspath))
                raise

        try:
            os.remove(obspath)
        except OSError:
            pass

        observed_interactome = observed_interactome.rename(columns={
            'mol_pdb': 'nbr_pdb',
            'mol_chain': 'nbr_chain',
            'mol_domNo': 'nbr_domNo',
            'mol_res': 'nbr_res',
            'mol_gi':'nbr_gi',
            })

        #Add in neghibor information from observed interactome
        inferred_interfaces["nbr_obs_int_id"] = inferred_interfaces["nbr_obs_int_id"].astype(int)
        inferred_interfaces = pd.merge(inferred_interfaces, observed_interactome,
            how="left", left_on="nbr_obs_int_id", right_on="obs_int_id",
            suffixes=["_inf", "_obs"])
        del observed_interactome

        #Select relevant columns
        if "int_superfam_id" in inferred_interfaces.columns:
            int_superfam_id_col = "int_superfam_id"
        elif "int_superfam_id_inf" in inferred_interfaces.columns:
            int_superfam_id_col = "int_superfam_id_inf"
        elif "int_superfam_id_x" in inferred_interfaces.columns:
            #Suffix not working?
            int_superfam_id_col = "int_superfam_id_x"
        else:
            raise RuntimeError("Merge faled for obs and inf: {}".format(inferred_interfaces.columns))
        try:
            inferred_interfaces = inferred_interfaces[[
                "sdi", "pdbId", "chnLett", "domNo", "from", "to", "resn", "resi",
                "gi", "sfam_id",
                "nbr_sdi_id", "nbr_pdb", "nbr_chain", "nbr_domNo", "nbr_res",
                "nbr_taxid", "nbr_gi", "nbr_superfam_id", "nbr_obs_int_id", "nbr_score",
                "int_sdi_id", "int_pdb", "int_chain", "int_domNo", "int_res",
                "int_sdi_from", "int_sdi_to", "int_taxid", "int_gi_x", "int_gi_y",
                int_superfam_id_col]]
        except KeyError as e:
            job.log("Unable to filter df. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
            raise

        #Rename columns
        inferred_interfaces = inferred_interfaces.rename(columns={
            "sdi":"mol_sdi_id",
            "pdbId":"mol_pdb",
            "chnLett":"mol_chain",
            "domNo":"mol_domNo",
            "from":"mol_sdi_from",
            "to":"mol_sdi_to",
            "resn":"mol_resn",
            "resi":"mol_resi",
            "gi":"mol_gi",
            "sfam_id":"mol_superfam_id",
            int_superfam_id_col:"int_superfam_id",
        })

        try:
            inferred_interfaces = inferred_interfaces[ \
                ~inferred_interfaces["mol_sdi"].isnull() & \
                ~inferred_interfaces["int_sdi"].isnull()]
        except KeyError as e:
            job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
            raise

        taxa = []
        for i, row in inferred_interfaces.itertuples():
            #Should only be one row, but maybe not
            try:
                taxas = filter_hdf(tax_file, "table",
                    pdbId = inferred_interfaces.iloc[0]["mol_pdb"],
                    chain = inferred_interfaces.iloc[0]["mol_chain"])
                taxa.append(float(taxas.iloc[0]["taxId"]))
            except (KeyboardInterrupt, SystemExit):
                raise
            else:
                taxa.append(np.NaN)
        taxa = pd.Series(taxa, index=inferred_interfaces.index)
        inferred_interfaces = inferred_interfaces.assign(mol_taxid=taxa)

        try:
            resi = pd.Series([decode_residues(job, row.mol_pdb, row.mol_chain, row.mol_resi, row) \
                for row in inferred_interfaces.itertuples()], index=inferred_interfaces.index)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            import traceback
            tb = traceback.format_exc()
            job.log("FAILED", tb)
            print tb
            raise

        inferred_interfaces = inferred_interfaces.assign(mol_res=resi)
        del resi
        del inferred_interfaces["mol_resi"]

        str_cols = ["mol_res", "int_res", "mol_resn", "int_resn", "mol_pdb",
                    "int_pdb", "mol_chain", "int_chain", "nbr_pdb", "nbr_chain"]

        for col in inferred_interfaces.columns:
            inferred_interfaces[col] = inferred_interfaces[col].astype(str if col in str_cols else float)

        df_file = job.fileStore.getLocalTempFileName()
        inferred_interfaces.to_hdf(unicode(df_file), "table", format="table",
            table=True, complevel=9, complib="bzip2", min_itemsize=1024,
            data_coumns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"])
        job.log("Wrote "+df_file)

        mol_sfams = inferred_interfaces["mol_sfam"].drop_duplicates()
        if len(mol_sfams) == 1:
            #Write to HDF file
            df_file = job.fileStore.getLocalTempFileName()
            inferred_interfaces.to_hdf(unicode(df_file), "table", format="table",
                table=True, complevel=9, complib="bzip2", min_itemsize=1024,
                data_coumns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"])
            job.log("Wrote "+df_file)
            df_files = [df_file]

            #Add ibis info into out store
            out_store.write_output_file(df_file, "{}/_infrows/Intrac{}/{}.h5".format(int(mol_sfams[0]), table, inf_int_id_file))
        else:
            df_files = []
            for inf_row in inferred_interfaces.iterrows():
                mol_sfam = inf_row["mol_sfam"]
                inf_row = inf_row.to_frame().T

                #Write to HDF file
                df_file = job.fileStore.getLocalTempFileName()
                inf_row.to_hdf(unicode(df_file), "table", format="table",
                    table=True, complevel=9, complib="bzip2", min_itemsize=1024,
                    data_coumns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"])
                job.log("Wrote "+df_file)
                df_files.append(df_file)

                #Add ibis info into out store
                out_store.write_output_file(df_file, "{}/_infrows/Intrac{}/{}.h5".format(int(mol_sfam), table, inf_int_id_file))


    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        import traceback
        tb = traceback.format_exc()
        job.log("FAILED {} {}".format(nbr_sfam, tb))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write(str(inf_int_id))
            f.write("\n")
            f.write(tb)
        out_store.write_output_file(fail_file, "{}/_infrows/Intrac{}/{}.failed".format(int(nbr_sfam), table, inf_int_id_file))
        try:
            os.remove(fail_file)
        except OSError:
             pass


        try:
            for f in (df_file, pdb_file, tableInfFile, obspath):
                try:
                    os.remove(f)
                except OSError:
                    pass
        except:
            pass
    finally:
        try:
            files = df_files if isrow else df_files+[pdb_file, tableInfFile, obspath]
            for f in files:
                try:
                    os.remove(f)
                except OSError:
                    pass
        except:
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
    possible_errors = []
    to_delete = []
    for row_prefix in out_store.list_input_directory("{}/_infrows/Intrac{}".format(int(sfam_id), table)):
        if row_prefix.endswith("failed"): continue
        job.log("Running  {} {}".format(sfam_id, row_prefix))

        row_file = os.path.join(work_dir, os.path.basename(row_prefix))

        try:
            out_store.read_input_file(row_prefix, row_file)
            df = pd.read_hdf(row_file, "table")
            for col, _ in df.dtypes[df.dtypes == 'int64'].iteritems():
                df[col] = df[col].astype(float)
            df.to_hdf(unicode(resi_path), "table", table=True, format="table", append=True, mode="a",
                data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)
            job.fileStore.deleteGlobalFile(row_prefix)
        except (SystemExit, KeyboardInterrupt):
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
        to_delete.append(row_prefix)


    if os.path.isfile(resi_path):
        #Upload to S3
        out_store.write_output_file(resi_path, "{}/_inftables/{}".format(int(sfam_id), resi_prefix))

        #Cleanup
        os.remove(resi_path)
        for key in to_delete:
            out_store.remove_file(key)

        print "End merge", sfam_id, table
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


def get_table_sfams(job, mol_sfam_id, table, tableInfStoreID, pdbFileStoreID, taxFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    pdbFilePath = get_file(job, "PDB.h5", pdbFileStoreID)
    taxFilePath = get_file(job, "pdb_chain_taxonomy.h5", taxFileStoreID)
    # obsFile = get_file(job, "IBIS_observed.h5", in_store)
    #
    # try:
    #     observed_interactome = filter_hdf_chunks("IBIS_observed.h5", "table", "obs_int_id", mol_sfam_id)
    # except (SystemExit, KeyboardInterrupt):
    #     raise
    # except:
    #     try:
    #         observed_interactome = filter_hdf_chunks("IBIS_observed.h5", "table", "obs_int_id", float(mol_sfam_id))
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except:
    #         job.log("Failed reading IBIS_observed.h5")
    #         return
    obsPath = os.path.join(work_dir, "{0}.observed_interactome".format(int(mol_sfam_id)))
    out_store.read_input_file("{0}/{0}.observed_interactome".format(int(mol_sfam_id)), obsPath)

    tableInfPath = get_file(job, "IBIS_inferred_{}.h5".format(table), tableInfStoreID)
    skip_int = set([tuple(map(int, os.path.basename(f)[:-3].split("_"))) for f in out_store.list_input_directory(
        "{}/_infrows/Intrac{}".format(int(mol_sfam_id),  table)) if f.endswith(".h5")])
    try:
        inf_int_ids = filter_hdf_chunks(tableInfPath, "Intrac{}".format(table), chunksize=100,
            columns=["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi_id"],
            nbr_superfam_id=mol_sfam_id)
    except (RuntimeError, TypeError):
        job.log("Unable to find sfam {} in table {}, Skipping".format(mol_sfam_id, table))
        return

    #inf_int_ids = set([tuple(row) for row in inf_int_ids.itertuples()])
    #inf_int_ids -= skip_int
    #print "Starting table sfam", mol_sfam_id, inf_int_ids

    #would this be better to just ran as a loop?
    #map_job(job, process_inferred_interaction, list(inf_int_ids), mol_sfam_id, table, tableInfStoreID, pdbFileStoreID, taxFileStoreID)
    try:
        for i, row in inf_int_ids.iterrows():
            #if tuple(row) in skip_int: continue
            RealtimeLogger.info("Running {}".format(row))
            process_inferred_interaction(job, row, mol_sfam_id, table, obsPath,
                pdbFilePath, taxFilePath, isrow=True, work_dir=work_dir)

        #job.addFollowOnJobFn(merge_table_sfam, mol_sfam_id, table)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise
    finally:
        for f in [tableInfPath, obsPath, pdbFilePath]:
            try:
                os.remove(f)
            except OSError:
                pass

def get_inferred_structural_interactome_by_table(job, table, pdbFileStoreID, taxFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    #Read in H5 for entire table
    tableInfPath = get_file(job, "IBIS_inferred_{}.h5".format(table), in_store)
    tableInfPathFileStoreID = job.fileStore.writeGlobalFile(tableInfPath)

    sfams = filter_hdf_chunks(tableInfPath, "Intrac{}".format(table), columns=["nbr_superfam_id"]).drop_duplicates().dropna()
    skip_sfam = set([int(f.split("/", 1)[0]) for f in out_store.list_input_directory() \
       if f.endswith(".inferred_interactome")])
    sfams = sfams[~sfams["nbr_superfam_id"].isin(skip_sfam)]["nbr_superfam_id"].drop_duplicates().dropna().astype(int).tolist()

    partial_sfams = set(int(k.split("/")[0]) for sfam in sfams for k in \
        out_store.list_input_directory(
            "{sfam}/_inftables/Intrac{table}_{sfam}.inferred_interactome".format( \
            sfam=sfam, table=table)) if not k.endswith("failed"))

    sfams = list(set(sfams)-partial_sfams)
    sfams = [305194]

    if len(sfams) > 0:
        map_job(job, get_table_sfams, sfams, table, tableInfPathFileStoreID, pdbFileStoreID, taxFileStoreID)

    try:
        os.remove(tableInfPath)
    except OSError:
        pass

def merge_inferred_interactome_sfam(job, sfam_id):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    iostore = IOStore.get("{}:molmimic-interfaces".format(prefix))

    sfam_prefix = "{}/_inftables".format(int(sfam_id))
    sfam_file = "{s}/{s}.inferred_interactome".format(s=int(sfam_id))
    merged_file = job.fileStore.getLocalTempFileName()

    to_delete = []
    for table_prefix in iostore.list_input_directory(sfam_prefix):
        if table_prefix.endswith("failed"): continue

        print "Running table sfam", table_prefix

        try:
            merge_table_sfam(job, mol_sfam_id, table)
            table_file = os.path.join(work_dir, os.path.basename(table_prefix))
            iostore.read_input_file(table_prefix, table_file)
            for df in pd.read_hdf(unicode(table_file), "table", chunksize=1000):
                df.to_hdf(unicode(merged_file), "table", mode="a", append=True, format="table",
                    table=True, complevel=9, complib="bzip2", min_itemsize=1024)
        except (IOError, ValueError):
            job.log("Failed to read {}".format(table_file))

        try:
            os.remove(table_file)
        except OSError:
            pass

        to_delete.append(table_prefix)


    if os.path.isfile(merged_file):
        #Write output file
        iostore.write_output_file(merged_file, sfam_file)

        #Cleanup
        try:
            os.remove(merged_file)
        except OSError:
            pass
        for key in to_delete:
            iostore.remove_file(key)

    else:
        failed_file = os.path.join(work_dir, "failed_file")
        with open(failed_file, "w") as f:
            print >>f, "No merged_file present"
        iostore.write_output_file(failed_file, sfam_file+".failed")
        try:
            os.remove(failed_file)
        except OSError:
            pass

def cleanup(job, sfam_ids):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    for sfam_id in sfam_ids:
        infrows = store.list_input_directory("{}/_infrows".format(int(sfam_id)))
        inftables = store.list_input_directory("{}/_inftables".format(int(sfam_id)))
        finished_tables = [k.split("/")[-1].split("_")[0] for k in inftables]
        for k in infrows:
            if k.split("/")[2] in finished_tables:
                store.remove_file(k)

def merge_inferred_interactome(job):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))
    all_sfam = [os.path.basename(f).split(".") for f in out_store.list_input_directory() if not f.endswith("failed")]
    skip_sfam = [f[0] for f in all_sfam if f[1] == "inferred_interactome"]
    sfam_to_run = [f[0] for f in all_sfam if f[1] == "observed_interactome" \
        and f[0] not in skip_sfam]
    #sfam_to_run = [305194]
    map_job(job, merge_inferred_interactome_sfam, sfam_to_run)
    job.addFollowOnJobFn(cleanup, sfam_to_run)

def start_toil(job):
    print "Starting job"
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    #Download PDB info
    pdb_file = os.path.join(work_dir, "PDB.h5")
    in_store.read_input_file("PDB.h5", pdb_file)

    #Add pdb info into local job store
    pdbFileStoreID = job.fileStore.writeGlobalFile(pdb_file)

    #Download PDB Taxonomy information
    tax_file = os.path.join(work_dir, "pdb_chain_taxonomy.h5")
    in_store.read_input_file("pdb_chain_taxonomy.h5", tax_file)

    #Add tax info into local job store
    taxFileStoreID = job.fileStore.writeGlobalFile(tax_file)

    tables = set(range(1,87))-set([51])

    job.log("Running tables: {}".format(tables))
    map_job(job, get_inferred_structural_interactome_by_table, tables,
        pdbFileStoreID, taxFileStoreID)
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
