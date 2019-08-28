,import os, sys
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
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

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
            inferred_interfaces = inf_int_id[1].to_frame().T
            #inf_int_id = inferred_interfaces[""]
            nbr_obs_int_id = inf_int_id[1]["nbr_obs_int_id"]
            nbr_sdi = inf_int_id[1]["nbr_sdi_id"]
            mol_sdi = inf_int_id[1]["mol_sdi_id"]
            inf_int_id_file = "{}_{}_{}_{}".format(inf_int_id[0], nbr_obs_int_id, nbr_sdi, mol_sdi)

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

        pdbs = struct_domains["pdbId"]
        taxa = pd.read_hdf(tax_file, "table", where="pdbId in pdbs")
        struct_domains = pd.merge(struct_domains, taxa, on=["pdbId", "chnLett"], how="left")

        #Add PDB, chain, and sdi information. Inner join to only allow sdi's that are in databses, not obsolete
        inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="inner", left_on="mol_sdi_id", right_on="sdi")
        #RealtimeLogger.info("INTERFACE: {}".format(inferred_interfaces.iloc[0]))

        if isrow:
            #RealtimeLogger.info("tableInfStoreID is {} {}".format(type(tableInfStoreID), nbr_obs_int_id))
            #RealtimeLogger.info("{}".format(tableInfStoreID.iloc[0]))
            observed_interactome = tableInfStoreID[tableInfStoreID["obs_int_id"]==nbr_obs_int_id]
            #RealtimeLogger.info("Got obs file from row")

        else:
            obspath = job.fileStore.getLocalTempFileName()
            out_store.read_input_file("{0}/{0}.observed_interactome".format(nbr_sfam), obspath)
            raise RuntimeError("Not a row!")

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

        observed_interactome = observed_interactome.rename(columns={
            'mol_pdb': 'nbr_pdb',
            'mol_chain': 'nbr_chain',
            'mol_domNo': 'nbr_domNo',
            'mol_res': 'nbr_res',
            'mol_gi_x': 'nbr_gi',
            #'mol_taxid': 'nbr_taxid',
            #'mol_superfam_id': 'nbr_superfam_id'
            })

        #Add in neghibor information from observed interactome
        inferred_interfaces["nbr_obs_int_id"] = inferred_interfaces["nbr_obs_int_id"].astype(int)
        inferred_interfaces = pd.merge(inferred_interfaces, observed_interactome,
            how="left", left_on="nbr_obs_int_id", right_on="obs_int_id",
            suffixes=["_inf", "_obs"])
        del observed_interactome
        #RealtimeLogger.info("INTERFACE: {}".format(inferred_interfaces.iloc[0]))

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
                "taxid", "gi", "sfam_id",
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
            "taxid":"mol_taxid",
            "gi":"mol_gi",
            "sfam_id":"mol_superfam_id",
            int_superfam_id_col:"int_superfam_id",
        })

        try:
            inferred_interfaces = inferred_interfaces[ \
                ~inferred_interfaces["mol_sdi_id"].isnull() & \
                ~inferred_interfaces["int_sdi_id"].isnull()]
        except KeyError as e:
            job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
            raise

        # taxa = []
        #
        # for row in inferred_interfaces.itertuples():
        #     #Should only be one row, but maybe not
        #     try:
        #         taxas = filter_hdf(tax_file, "table",
        #             pdbId = row.mol_pdb,
        #             chnLett = row.mol_chain)
        #         taxa.append(float(taxas.iloc[0]["taxId"]))
        #     except (KeyboardInterrupt, SystemExit):
        #         raise
        #     else:
        #         taxa.append(np.NaN)
        # taxa = pd.Series(taxa, index=inferred_interfaces.index)
        # inferred_interfaces = inferred_interfaces.assign(mol_taxid=taxa)

        try:
            resi = pd.Series([decode_residues(job, row.mol_pdb, row.mol_chain, row.mol_resi, row) \
                for row in inferred_interfaces.itertuples()], index=inferred_interfaces.index)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            import traceback
            tb = traceback.format_exc()
            RealtimeLogger.info("FAILED {}".format(tb))
            raise
            resi = inferred_interfaces["mol_resi"].copy()


        inferred_interfaces = inferred_interfaces.assign(mol_res=resi)
        del resi
        del inferred_interfaces["mol_resi"]
        del inferred_interfaces["mol_resn"]

        str_cols = ["mol_pdb", "mol_chain", "mol_res",
                    "int_pdb",  "int_chain", "int_res",
                    "nbr_pdb", "nbr_chain", "nbr_res"]

        #RealtimeLogger.info("INTERFACE: {}".format(inferred_interfaces.iloc[0]))

        for col in inferred_interfaces.columns:
            inferred_interfaces[col] = inferred_interfaces[col].astype(
                str if col in str_cols else float)

        mol_sfams = inferred_interfaces["mol_superfam_id"].drop_duplicates()
        if len(mol_sfams) == 0:
            return
        elif len(mol_sfams) == 1:
            #Write to HDF file
            df_file = job.fileStore.getLocalTempFileName()
            inferred_interfaces.to_hdf(str(df_file), "table", format="table",
                table=True, complevel=9, complib="bzip2", min_itemsize=1024,
                data_coumns=["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi_id", "int_sdi_id"])
            job.log("Wrote "+df_file)
            df_files = [df_file]

            #Add ibis info into out store
            out_store.write_output_file(df_file, "{}/_infrows/Intrac{}/{}.inf.h5".format(int(mol_sfams[0]), table, inf_int_id_file))
        else:
            df_files = []
            for i, inf_row in inferred_interfaces.iterrows():
                mol_sfam = inf_row["mol_superfam_id"]
                inf_row = inf_row.to_frame().T
                for col in inf_row.columns:
                    inf_row[col] = inf_row[col].astype(str if col in str_cols else float)

                #Write to HDF file
                df_file = job.fileStore.getLocalTempFileName()
                inf_row.to_hdf(str(df_file), "table", format="table",
                    table=True, complevel=9, complib="bzip2", min_itemsize=1024,
                    data_coumns=["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi_id", "int_sdi_id"])
                job.log("Wrote "+df_file)
                df_files.append(df_file)

                #Add ibis info into out store
                out_store.write_output_file(df_file, "{}/_infrows/Intrac{}/{}.inf.h5".format(int(mol_sfam), table, inf_int_id_file))


    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        import traceback
        tb = traceback.format_exc()
        job.log("FAILED {} {}".format(nbr_sfam, tb))

        if not isrow:
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
        else:
            return inf_int_id, tb


        try:
            for f in (df_file, pdb_file, tableInfFile):
                try:
                    os.remove(f)
                except OSError:
                    pass
        except:
            pass
    finally:
        try:
            files = df_files if isrow else df_files+[pdb_file, tableInfFile]
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
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

    status = "inferred (table {}) {}".format(table, sfam_id)
    resi_prefix = "Intrac{}_{}.inf.h5".format(table, int(sfam_id))
    new_cols = ["mol_res"]
    data_cols = ["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"]

    resi_path = os.path.join(work_dir, resi_prefix)

    #Combine residues into dataframe
    possible_errors = []
    to_delete = []
    row_key = "{}/_infrows/Intrac{}".format(int(sfam_id), table)
    RealtimeLogger.info("Merginf rows from {}".format(row_key))
    for row_prefix in out_store.list_input_directory(row_key):
        if not row_prefix.endswith(".inf.h5"): continue
        RealtimeLogger.info("Running  {} {}".format(sfam_id, row_prefix))

        row_file = os.path.join(work_dir, os.path.basename(row_prefix))

        try:
            out_store.read_input_file(row_prefix, row_file)
            df = pd.read_hdf(row_file, "table")
            for col, _ in list(df.dtypes[df.dtypes == 'int64'].items()):
                df[col] = df[col].astype(float)
            df.to_hdf(str(resi_path), "table", table=True, format="table", append=True, mode="a",
                data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)

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
        #to_delete.append(row_prefix)


    if os.path.isfile(resi_path):
        #Upload to S3
        outfile = "{}/_inftables/{}".format(int(sfam_id), resi_prefix)
        outprefix = "{}/_inftables/{}".format(int(sfam_id), resi_prefix)
        out_store.write_output_file(resi_path, outprefix)

        #Cleanup
        os.remove(resi_path)
        for key in to_delete:
            out_store.remove_file(key)

        return outprefix
    else:
        RealtimeLogger.info("Failed merging: {}".format(resi_path))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write("no rows?\n")
            for e in possible_errors:
                f.write(e)
                f.write("\n")
        out_store.write_output_file(fail_file, "{}/_inftables/{}.failed".format(int(sfam_id), resi_prefix))


def get_table_sfams(job, mol_sfam_id, table, tableInfStoreID, pdbFileStoreID, taxFileStoreID, sfamFileStoreIDs):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

    RealtimeLogger.info("Running table {} sfam {}".format(table, mol_sfam_id))


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
    sfamFileStoreID = sfamFileStoreIDs[mol_sfam_id]
    obsFilePath = get_file(job, "{}_obs.h5".format(int(mol_sfam_id)),
        sfamFileStoreID, work_dir=work_dir)

    observed_interactome = pd.read_hdf(obsFilePath, "table")
    RealtimeLogger.info("Obs has {} rows".format(observed_interactome.shape))

    # obsFilePath = os.path.join(work_dir, "{0}.observed_interactome".format(int(mol_sfam_id)))
    # out_store.read_input_file("{0}/{0}.observed_interactome".format(int(mol_sfam_id)), obsPath)

    tableInfPath = get_file(job, "IBIS_inferred_{}.h5".format(table), tableInfStoreID)
    # skip_int = set([tuple(map(int, os.path.basename(f)[:-3].split("_"))) for f in out_store.list_input_directory(
    #     "{}/_infrows/Intrac{}".format(int(mol_sfam_id),  table)) if f.endswith(".h5")])
    try:
        inf_int_ids = filter_hdf_chunks(tableInfPath, "Intrac{}".format(table), chunksize=100,
            nbr_superfam_id=mol_sfam_id)
    except (RuntimeError, TypeError):
        job.log("Unable to find sfam {} in table {}, Skipping".format(mol_sfam_id, table))
        return

    #inf_int_ids = set([tuple(row) for row in inf_int_ids.itertuples()])
    #inf_int_ids -= skip_int

    #would this be better to just ran as a loop?
    #map_job(job, process_inferred_interaction, list(inf_int_ids), mol_sfam_id, table, tableInfStoreID, pdbFileStoreID, taxFileStoreID)
    try:
        fail_file = os.path.join(work_dir, "fail_file")
        for row in inf_int_ids.iterrows():
            #if tuple(row) in skip_int: continue
            nbr_obs_int_id = row[1]["nbr_obs_int_id"]
            nbr_sdi = row[1]["nbr_sdi_id"]
            mol_sdi = row[1]["mol_sdi_id"]
            inf_int_id_file = "{}_{}_{}_{}".format(row[0], nbr_obs_int_id, nbr_sdi, mol_sdi)

            if out_store.exists("{}/_infrows/Intrac{}/{}.inf.h5".format(int(mol_sfam_id), table, inf_int_id_file)):
                continue

            RealtimeLogger.info("Running {}".format(row))
            out = process_inferred_interaction(
                job,
                row,
                mol_sfam_id,
                table,
                observed_interactome,
                pdbFilePath,
                taxFilePath,
                isrow=True, work_dir=work_dir)
            if out is not None:
                inf_int_id, tb = out
                with open(fail_file, "a") as f:
                    f.write(str(inf_int_id))
                    f.write("\n")
                    f.write(tb)
        if os.path.isfile(fail_file):
            out_store.write_output_file(fail_file, "{}/_infrows/Intrac{}/failed".format(int(mol_sfam_id), table))
            try:
                os.remove(fail_file)
            except OSError:
                 pass
        #job.addFollowOnJobFn(merge_table_sfam, mol_sfam_id, table)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise
    finally:
        for f in [tableInfPath, obsFilePath, pdbFilePath, taxFilePath]:
            try:
                os.remove(f)
            except OSError:
                pass

def get_inferred_structural_interactome_by_table(job, table, pdbFileStoreID, taxFileStoreID, sfamFileStoreIDs):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

    RealtimeLogger.info("Running table {}".format(table))

    #Read in H5 for entire table
    tableInfPath = get_file(job, "IBIS_inferred_{}.h5".format(table), in_store)
    tableInfPathFileStoreID = job.fileStore.writeGlobalFile(tableInfPath)

    sfams = filter_hdf_chunks(tableInfPath, "Intrac{}".format(table),
        columns=["nbr_superfam_id"]).drop_duplicates().dropna()
    skip_sfam = set([s for s in sfams["nbr_superfam_id"] if \
        out_store.exists("{0}/{0}.inferred_interactome".format(int(s))) or \
        not out_store.exists("{0}/{0}.observed_interactome".format(int(s)))])

    # skip_sfam = set([int(f.split("/", 1)[0]) for f in out_store.list_input_directory() \
    #    if f.endswith(".inferred_interactome")])

    sfams = sfams[~sfams["nbr_superfam_id"].isin(skip_sfam)]
    sfams = sfams["nbr_superfam_id"].drop_duplicates().dropna().astype(int).tolist()

    # partial_sfams = set(int(k.split("/")[0]) for sfam in sfams for k in \
    #     out_store.list_input_directory(
    #         "{sfam}/_inftables/Intrac{table}_{sfam}.inferred_interactome".format( \
    #         sfam=sfam, table=table)) if not k.endswith("failed"))

    #sfams = list(set(sfams)-partial_sfams)

    if len(sfams) > 0:
        map_job(job, get_table_sfams, sfams, table, tableInfPathFileStoreID,
            pdbFileStoreID, taxFileStoreID, sfamFileStoreIDs)

    try:
        os.remove(tableInfPath)
    except OSError:
        pass

def merge_inferred_interactome_sfam(job, sfam_id):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    iostore = IOStore.get("aws:us-east-1:molmimic-interfaces")

    sfam_prefix = "{}/_infrows".format(int(sfam_id))
    sfam_file = "{s}/{s}.inferred_interactome".format(s=int(sfam_id))
    merged_file = job.fileStore.getLocalTempFileName()

    to_delete = []
    done_tables = []
    for table_prefix in iostore.list_input_directory(sfam_prefix):
        if table_prefix.endswith("failed"): continue

        table = int(os.path.basename(os.path.dirname(table_prefix)).replace("Intrac", ""))

        if table in done_tables:
            continue

        RealtimeLogger.info("Running table sfam {}".format(table_prefix))

        try:
            RealtimeLogger.info("Merge {} {}".format(sfam_id, table))
            table_sfam_prefix = merge_table_sfam(job, sfam_id, table)
            if table_sfam_prefix is None:
                RealtimeLogger.info("Merging failed for {} {}".format(sfam_id, table))
            table_file = os.path.join(work_dir, os.path.basename(table_sfam_prefix))
            iostore.read_input_file(table_sfam_prefix, table_file)
            for df in pd.read_hdf(str(table_file), "table", chunksize=1000):
                df.to_hdf(str(merged_file), "table", mode="a", append=True, format="table",
                    table=True, complevel=9, complib="bzip2", min_itemsize=1024)
            to_delete.append(table_sfam_prefix)
        except (IOError, ValueError) as e:
            raise
            job.log("Failed to read {} bc {}".format(table_prefix, e))

        try:
            os.remove(table_file)
        except OSError:
            pass

        done_tables.append(table)


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
            rows_prefix = "{}/_infrows".format(int(sfam_id))
            if iostore.exists(rows_prefix):
                for f in iostore.list_input_directory(rows_prefix):
                    iostore.remove(f)

    else:
        failed_file = os.path.join(work_dir, "failed_file")
        with open(failed_file, "w") as f:
            print("No merged_file present", file=f)
        iostore.write_output_file(failed_file, sfam_file+".failed")
        try:
            os.remove(failed_file)
        except OSError:
            pass

def cleanup(job, sfam_ids):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    store = IOStore.get("aws:us-east-1:molmimic-interfaces")

    for sfam_id in sfam_ids:
        infrows = store.list_input_directory("{}/_infrows".format(int(sfam_id)))
        inftables = store.list_input_directory("{}/_inftables".format(int(sfam_id)))
        finished_tables = [k.split("/")[-1].split("_")[0] for k in inftables]
        for k in infrows:
            if k.split("/")[2] in finished_tables:
                store.remove_file(k)

def merge_inferred_interactome(job, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

    pdb_file = get_file(job, "PDB.h5", pdbFileStoreID)
    sfams = pd.read_hdf(pdb_file, "Superfamilies", columns=
        ["sfam_id"]).drop_duplicates()["sfam_id"]
    os.remove(pdb_file)

    skip_sfam = [s for s in sfams if out_store.exists(
        "{0}/{0}.inferred_interactome".format(s))]

    sfam_to_run = [s for s in sfams if out_store.exists(
        "{0}/{0}.observed_interactome".format(s)) and s not in skip_sfam]

    # all_sfam = [os.path.basename(f).split(".") for f in out_store.list_input_directory() if not f.endswith("failed")]
    # skip_sfam = [f[0] for f in all_sfam if f[1] == "inferred_interactome"]
    # sfam_to_run = [f[0] for f in all_sfam if f[1] == "observed_interactome" \
    #     and f[0] not in skip_sfam]
    map_job(job, merge_inferred_interactome_sfam, sfam_to_run)
    job.addFollowOnJobFn(cleanup, sfam_to_run)

def start_toil(job):
    work_dir = job.fileStore.getLocalTempDir()
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    int_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

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

    sfams = pd.read_hdf(pdb_file, "Superfamilies", columns=
        ["sfam_id"]).drop_duplicates().dropna()["sfam_id"].sort_values()
    #RealtimeLogger.info("SFAMS: {}".format(sfams.shape[0]))
    sfamFileStoreIDs = {}
    for s in sfams:
        k = "{0}/{0}.observed_interactome".format(int(s))
        if int_store.exists(k):
            RealtimeLogger.info("Loading {}".format(s))
            f = job.fileStore.getLocalTempFileName()
            int_store.read_input_file(k, f)
            sfamFileStoreIDs[int(s)] = job.fileStore.writeGlobalFile(f)
            os.remove(f)
        else:
            RealtimeLogger.info("FAILED Loading {}".format(s))

    assert len(sfamFileStoreIDs) > 0

    os.remove(tax_file)
    os.remove(pdb_file)

    job.log("Running tables: {}".format(tables))
    j = job
    for table in table:
        j.addFollowOnJobFn(get_inferred_structural_interactome_by_table, table,
            pdbFileStoreID, taxFileStoreID, sfamFileStoreIDs)
    # map_job(job, get_inferred_structural_interactome_by_table, tables,
    #     pdbFileStoreID, taxFileStoreID, sfamFileStoreIDs)
    j.addFollowOnJobFn(merge_inferred_interactome, pdbFileStoreID)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
