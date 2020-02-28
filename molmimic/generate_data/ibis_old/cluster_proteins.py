from __future__ import print_function
import os, sys
import subprocess
import itertools as it
import glob
import re
from collections import defaultdict

from toil.realtimeLogger import RealtimeLogger

try:
    import pandas as pd
    from molmimic.generate_data.iostore import IOStore
    from molmimic.generate_data.job_utils import map_job
    from molmimic.generate_data.util import get_file, PDB_TOOLS, izip_missing, get_pdb_residues, natural_keys, remove_ter_lines
    from molmimic.parsers.USEARCH import run_usearch
    from molmimic.parsers import superpose
except ImportError as e:
    RealtimeLogger.info(e)

def merge_all_sfam_clusters(job, sfam_id, interaction_type, id):
    assert interaction_type in ("observed", "inferred")
    work_dir = job.fileStore.getLocalTempDir()
    cluster_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")

    cluster_interfaces = None

    to_delete = []
    for cluster_key in cluster_store.list_input_directory("{}/interface_clusters_{}_{}".format(sfam_id, id, interaction_type)):
        cluster_file = os.path.join(work_dir, os.path.basename(cluster_key))

        try:
            cluster_store.read_input_file(cluster_key, cluster_file)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            RealtimeLogger.info("Error getting cluster for : {}".format(cluster_key))
            continue

        try:
            row = pd.read_csv(cluster_file, index_col=False)
            cluster_interfaces = pd.concat((cluster_interfaces, row), axis=0) if \
                cluster_interfaces is not None else row.copy()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            RealtimeLogger.info("Error reading cluster for : {}".format(cluster_key))
            continue

        del row
        try:
            os.remove(cluster_file)
        except OSError:
            continue

        to_delete.append(cluster_key)

    if cluster_interfaces is None:
        RealtimeLogger.info("Error reading clusters for : {} {} {}".format(sfam_id, interaction_type, id))
        return

    interface_file = os.path.join(work_dir, "{}_{}.csv".format(sfam_id, id))
    interface_key = "{}/interface_clusters_{}.h5".format(sfam_id, id)
    if interaction_type == "inferred":
        interface_key += ".inferred"

    cluster_interfaces.to_csv(interface_file, sep="\t")
    cluster_store.write_output_file(interface_file, interface_key)

    try:
        os.remove(interface_file)
    except OSError:
        pass

    # for k in to_delete:
    #     cluster_store.remove_file(k)


def merge_sfam_cluster(job, cluster_num, sfam_id, interaction_type, id, interfaceStoreId, structure_map=True):
    assert interaction_type in ("observed", "inferred")
    from Bio import SeqIO
    work_dir = job.fileStore.getLocalTempDir()
    cluster_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")
    structure_store = IOStore.get("aws:us-east-1:molmimic-full-structures")

    if interaction_type == "inferred":
        centroid_key = "{}/interface_clusters_{}_{}/{}.csv".format(sfam_id, id, interaction_type, cluster_num)
    else:
        centroid_key = "{}/interface_clusters_{}/{}.csv".format(sfam_id, id, cluster_num)

    if cluster_store.exists(centroid_key):
        return

    resolution_re = re.compile("\[resolution: (\d\.)+\]")

    interface_file = get_file(job, "{}.h5".format(sfam_id), interfaceStoreId)
    interfaces = pd.read_hdf(interface_file, "table")

    cluster_key = "{}/clusters_{}/{}.fasta".format(sfam_id, id, cluster_num)
    cluster_file = os.path.join(work_dir, os.path.basename(cluster_key))

    try:
        cluster_store.read_input_file(cluster_key, cluster_file)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        RealtimeLogger.info("Error getting cluster for : {}".format(cluster_key))
        return

    interfaces["mol_pdb"] = interfaces["mol_pdb"].astype(str)
    interfaces["mol_chain"] = interfaces["mol_chain"].astype(str)


    msa = list(SeqIO.parse(cluster_file, "fasta"))
    cluster = {}
    centroid_res = None
    centroid_id = None
    for i, record in enumerate(msa):
        pdb, chain, sdi, domNo = record.id[:-4].split("_")
        sdi, domNo = int(sdi[3:]), int(domNo[1:])

        m = resolution_re.match(record.description)
        resolution = float(m.groups(1)) if m else -1.

        if centroid_id is None or (resolution != -1. and resolution<centroid_res):
            centroid_id = (pdb, chain, sdi, domNo)
            centroid_res = resolution

        if interaction_type=="inferred":
            interface = interfaces[ \
                (interfaces.mol_pdb==pdb)&(interfaces.mol_chain==chain)&\
                (interfaces.mol_sdi==sdi)&(interfaces.mol_domNo==domNo)]
        else:
            interface = interfaces[ \
                (interfaces.mol_pdb==pdb)&(interfaces.mol_chain==chain)&\
                (interfaces.mol_sdi_id==sdi)&(interfaces.mol_domNo==domNo)]

        seq_ungapped_to_gapped = [i for i, aa in enumerate(str(record.seq)) if aa != "-"]

        structure_key = "{}/{}/{}_{}_sdi{}_d{}.pdb".format(sfam_id, pdb[1:3].lower(), pdb, chain, sdi, domNo)
        structure_file = os.path.join(work_dir, os.path.basename(structure_key))
        try:
            structure_store.read_input_file(structure_key, structure_file)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            RealtimeLogger.info("Error getting pdb for : {}".format(structure_key))
            continue

        seq_resi = list(get_pdb_residues(structure_file))

        try:
            os.remove(structure_file)
        except OSError:
            pass

        cluster[(pdb, chain, sdi, domNo)] = (interface, seq_ungapped_to_gapped, seq_resi, resolution)

    centroid = cluster[centroid_id]
    centroid_ungapped_to_gapped = centroid[1]
    centroid_gapped_to_ungapped = {j:i for i, j in enumerate(centroid_ungapped_to_gapped)}

    centroid_ungapped_to_resi = centroid[2]

    all_mol_res = set()
    merged_int_ids = defaultdict(list)
    int_sfams = defaultdict(list)
    try:
        centroid_row = centroid[0].iloc[0]
        for i, row in centroid[0].iterrows():
            try:
                mol_res = str(row.mol_res).split(",")
            except AttributeError:
                continue
            for r in mol_res:
                if interaction_type == "observed":
                    merged_int_ids[r].append(str(int(row.obs_int_id)))
                else:
                    merged_int_ids[r].append("{}_{}_{}".format(
                        row.nbr_obs_int_id, row.nbr_sdi_id, row.mol_sdi))

                try:
                    int_sfams[r].append(str(int(row.int_superfam_id)))
                except ValueError:
                    int_sfams[r].append("None")
            all_mol_res |= set(mol_res)
    except IndexError:
        #Centroid has no known interfaces: So start with 2nd best in cluster and
        #replace values with centroid_id
        sorted_clusters = sorted(cluster.iteritems(), key=lambda x: x[1][-1])

        #Look for xray structures
        for c in sorted_clusters:
            if c[1][-1] >= 0 or c[1][0].shape[0] > 0:
                centroid_row = c[1][0].iloc[0]
                break
        else:
            #Look for NMR structures
            for c in sorted_clusters:
                if c[1][0].shape[0] > 0:
                    centroid_row = c[1][0].iloc[0]
                    break
            else:
                RealtimeLogger.info("Centroid ({}) has no interfaces".format(centroid_id))
                return

        centroid_row["mol_pdb"] = centroid_id[0]
        centroid_row["mol_chain"] = centroid_id[1]
        if interaction_type == "inferred":
            centroid_row["mol_sdi"] = centroid_id[2]
        else:
            centroid_row["mol_sdi_id"] = centroid_id[2]
        centroid_row["mol_domNo"] = centroid_id[3]

    for int_id, (rows, u2g, resi, resolu) in cluster.iteritems():
        if int_id == centroid_id: continue

        if rows is not None:
            if structure_map:
                _structure_map = sequence_structure_align(mol_sfam, int_id, centroid_id)

            for _, row in rows.iterrows():
                if isinstance(row["mol_res"], str):
                    _mol_res = row["mol_res"].split(",")
                elif isinstance(row["mol_res"], (float, int)):
                    _mol_res = [str(row["mol_res"])]
                else:
                    RealtimeLogger.info("mol_res is in a non supported format (str, float, int): {} {}".format(
                        row["mol_res"], type(row["mol_res"])))
                    continue
                for mol_res in _mol_res:
                    mol_res = natural_keys(mol_res, use_int=True)
                    try:
                        ungapped = resi.index(mol_res)
                        if structure_map:
                            try:
                                centroid_ungapped = _structure_map[ungapped]
                            except KeyError:
                                continue
                        else:
                            gapped = u2g[ungapped]
                            centroid_ungapped = centroid_gapped_to_ungapped[gapped]
                        centroid_resi = centroid_ungapped_to_resi[centroid_ungapped]

                    except (ValueError, KeyError, IndexError) as e:
                        RealtimeLogger.info("Failed converting residue {}".format(mol_res))
                        continue

                    all_mol_res.add(centroid_resi)

                    if interaction_type == "observed":
                        merged_int_ids[centroid_resi].append(str(int(row.obs_int_id)))
                    else:
                        merged_int_ids[centroid_resi].append("{}_{}_{}".format(
                            row.nbr_obs_int_id, row.nbr_sdi_id, row.mol_sdi))

                    try:
                        int_sfams[centroid_resi].append(str(int(row.int_superfam_id)))
                    except ValueError:
                        int_sfams[centroid_resi].append("None")



    centroid_row["mol_res"] = ",".join(sorted(["".join(map(str, x)).strip() for x in all_mol_res], key=natural_keys))

    col_name = "obs_int_ids" if interaction_type=="observed" else "inf_int_ids"
    new_cols = pd.Series({
        col_name: ",".join(["{}:{}".format("".join(map(str, r)).strip(), \
            ",".join(int_ids)) for r, int_ids in merged_int_ids.iteritems()]),
        "int_sfam_ids": ",".join(["{}:{}".format("".join(map(str, r)).strip(), \
            ",".join(sfam_ids)) for r, sfam_ids in int_sfams.iteritems()])})

    centroid_row = pd.concat((centroid_row, new_cols))

    centroid_file = os.path.join(work_dir, os.path.basename(centroid_key))

    if len(centroid_row) > 0:
        centroid_row.to_frame().T.to_csv(centroid_file, index=False)
        cluster_store.write_output_file(centroid_file, centroid_key)

    for f in [interface_file, centroid_file]:
        try:
            os.remove(f)
        except OSError:
            pass

def merge_sfam_clusters(job, sfam_id, interaction_type, id=0.9):
    assert interaction_type in ("observed", "inferred")
    work_dir = job.fileStore.getLocalTempDir()
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")
    cluster_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")

    interface_key = "{0}/{0}.{1}_interactome".format(int(sfam_id), interaction_type)
    interface_file = os.path.join(work_dir, "{0}.{1}_interactome".format(int(sfam_id), interaction_type))

    try:
        interface_store.read_input_file(interface_key, interface_file)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        RealtimeLogger.info("Error getting interfaces for : {}".format(interface_key))
        return

    interfaceStoreId = job.fileStore.writeGlobalFile(interface_file, cleanup=True)

    clusters = [int(os.path.splitext(os.path.basename(key))[0]) for key in \
        cluster_store.list_input_directory("{}/clusters_{}".format(int(sfam_id), id))]
    clusters = [c for c in clusters if not cluster_store.exists(
        "{}/interface_clusters_{}/{}.csv".format(int(sfam_id), id, c))]

    #clusters = [484]
    RealtimeLogger.info("CLUSTERS TO RUN: {}".format(clusters))

    map_job(job, merge_sfam_cluster, clusters, str(int(sfam_id)), interaction_type, id, interfaceStoreId)

    job.addFollowOnJobFn(merge_all_sfam_clusters, sfam_id, interaction_type, id)

    try:
        os.remove(interface_file)
    except OSError:
        pass

def pdb2seq(pdb_fname, sfam_id, f, fasta):
    try:
        seq = subprocess.check_output([sys.executable, os.path.join(PDB_TOOLS,
            "pdb_toseq.py")], stdin=f)
        fasta.write(">{}\n{}\n".format(pdb_fname, "\n".join(seq.splitlines()[1:])))
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        RealtimeLogger.info("Error getting fasta for : {} {} {}".format(sfam_id, pdb_fname, e))
        pass

def cluster(job, sfam_id, id, work_dir=None):
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()
    out_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")

    uclust_file_key = "{0}/{0}_clusters_{1}.h5".format(int(sfam_id), id)
    if not out_store.exists(uclust_file_key):
        domain_fasta_key = "{0}/{0}.fasta".format(int(sfam_id))
        domain_fasta = os.path.join(work_dir, "{0}.fasta".format(int(sfam_id)))

        try:
            out_store.read_input_file(domain_fasta_key, domain_fasta)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            RealtimeLogger.info("Error getting fasta for : {}".format(sfam_id))
            return

        clusters_file, uclust_file, msa_file = run_usearch(["-cluster_fast",
            "{i}"+domain_fasta, "-id", str(id),
            "-centroids", "{{o}}{}_clusters.centroids".format(int(sfam_id)),
            "-uc", "{{o}}{}_clusters.uc".format(int(sfam_id)),
            "-msaout", "{{o}}{}_msa.".format(int(sfam_id))], work_dir=work_dir, job=job)

        for f in glob.glob(msa_file+"*"):
            cluster = f.rsplit(".", 1)[-1]
            out_store.write_output_file(f, "{}/clusters_{}/{}.fasta".format(int(sfam_id), id, cluster))

        #Convert uclust to h5
        uclust = pd.read_table(str(uclust_file), comment="#", header=None, names=[
            "record_type",
            "cluster",
            "length",
            "pctId",
            "strand",
            "unk1",
            "unk2",
            "alignment",
            "label_query",
            "label_target"
        ])
        del uclust["unk1"]
        del uclust["unk2"]
        hdf_base = "{}_clusters_{}.h5".format(int(sfam_id), id)
        hdf_file = os.path.join(work_dir, hdf_base)
        uclust.to_hdf(str(hdf_file), "table", complevel=9, complib="bzip2")
        uclust_file_key = "{0}/{0}_clusters_{1}.h5".format(int(sfam_id), id)
        out_store.write_output_file(hdf_file, uclust_file_key)
        os.remove(uclust_file)
    else:
        uclust_file = os.path.join(work_dir, os.path.basename(uclust_file_key))
        out_store.read_input_file(uclust_file_key, uclust_file)
        uclust = pd.read_hdf(str(uclust_file), "table")

    interface_key = "{}/interface_clusters_{}.h5".format(sfam_id, id)
    # if not out_store.exists(interface_key):
    #     job.addFollowOnJobFn(merge_sfam_cluster s, sfam_id, "observed", id=id)

    structure_alignments = []
    for cluster, cluster_group in uclust[uclust["record_type"]=="H"].groupby("cluster"):
        centroid_file = cluster_group.iloc[0]["label_target"].split(" ", 1)[0]
        centroid_chain = centroid_file.split("_")[1]
        centroid_file = os.path.join(work_dir, centroid_file)

        struct_align_file = os.path.join(work_dir, "{}.struct_align".format(cluster))
        with open(struct_align_file, "w") as struct_aln_fh:
            for hit in cluster_group["label_query"]:
                hit_file = hit.split(" ", 1)[0]
                hit_chain = hit_file.split("_")[1]
                hit_file = os.path.join(work_dir, hit_file)
                _, _, _, _, moving_ungapped_to_gapped, target_ungapped_to_gapped = superpose.align(
                    centroid_file, centroid_chain, hit_file, hit_chain,
                    extract=False, parse_postions=True, work_dir=work_dir, job=job)
                print(">{} {}\n{}\{}".format(
                    hit_file, centroid_file,
                    ",".join(map(str, moving_ungapped_to_gapped)),
                    ",".join(map(str, target_ungapped_to_gapped)),
                ), file=struct_aln_fh)

        out_store.write_output_file(
            struct_align_file,
            "{}/clusters_{}/{}.struct_aln".format(int(sfam_id), id, cluster))


    if not out_store.exists(interface_key+".inferred"):
        job.addFollowOnJobFn(merge_sfam_clusters, sfam_id, "inferred", id=id)

def cluster_structures(job, mol_sfam, pdb_complexes, work_dir=None):
    from molmimic.parsers.MaxCluster import run_maxcluster, get_centroid, get_aligned_residues
    from molmimic.generate_data.util import PDB_TOOLS, SubprocessChain

    cluster_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")

    if work_dir is None:
        work_dir = os.getcwd()

    cluster_dir = os.path.join(work_dir, "{}_cluster".format(int(mol_sfam)))
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    file_list = os.path.join(cluster_dir, "{}_file_list.txt".format(int(mol_sfam)))
    log_file = os.path.join(cluster_dir, "{}.max_cluster_logs".format(int(mol_sfam)))
    distance_file = os.path.join(cluster_dir, "{}_distances.txt".format(int(mol_sfam)))
    transformations_file = os.path.join(cluster_dir, "{}_transformations.txt".format(int(mol_sfam)))

    file_key = "{}/{}".format(int(mol_sfam), os.path.basename(file_list))
    log_key = "{}/{}".format(int(mol_sfam), os.path.basename(log_file))
    distance_key = "{}/{}".format(int(mol_sfam), os.path.basename(distance_file))
    transformation_key = "{}/{}".format(int(mol_sfam), os.path.basename(transformations_file))

    with open(file_list, "w") as fh:
        for f in pdb_complexes:
            if os.path.isfile(f):
                print(f, file=fh)
                outfile, _, _, _, moving_ungapped_to_gapped, target_ungapped_to_gapped = superpose.align(
                    centroid_file, centorid_chain, moving_file, moving_chain,
                    extract=False, parse_postions=True, work_dir=work_dir, job=job)

    if cluster_store.exists(log_key):
        try:
            cluster_store.read_input_file(log_key, log_file)
            cluster_store.read_input_file(distance_key, distance_file)
        except ClientError:
            return None, None, None
        logs = log_file
    else:
        if len(pdb_complexes)>1:
            logs = run_maxcluster("rmsd", file_list=file_list, log=True, R=distance_file,
                M=transformations_file, work_dir=cluster_dir, C=1, P=10, job=job)
        elif len(pdb_complexes):
            return 0, pdb_complexes[0], None
        else:
            RealtimeLogger.info("No files exits!!")
            return None, None, None

        if logs is None:
            return None, None, None

        cluster_store.write_output_file(file_list, file_key)
        cluster_store.write_output_file(logs, log_key)
        cluster_store.write_output_file(distance_file, distance_key)
        cluster_store.write_output_file(transformation_file, transformation_key)

    centroid = get_centroid(logs, job=job)
    for centroid_index, p in successful_pdbs:
        if p and p.endswith(centroid):
            break
    else:
        centroid_index = None

    aligned_residues = get_aligned_residues(file_list, transformations_file)

    return None, centroid

def sequence_structure_align(mol_sfam, source, target, work_dir=None):
    RealtimeLogger.info("Starting struct align")
    if work_dir is None:
        work_dir = os.getcwd()

    cluster_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")

    file_list = os.path.join(work_dir, "{}_file_list.txt".format(int(mol_sfam)))
    transformations_file = os.path.join(work_dir, "{}_{}_transformations.txt".format(int(mol_sfam)))

    file_key = "{}/{}".format(int(mol_sfam), os.path.basename(file_list))
    transformation_key = "{}/{}".format(int(mol_sfam), os.path.basename(transformations_file))

    cluster_store.read_input_file(file_key, file_list)
    cluster_store.read_input_file(transformation_key, transformation_file)

    aligned_residues = get_aligned_residues(file_list, transformations_file)

    for key, positions in aligned_residues.iteritems():
        try:
            source_index = key.index(source)
            target_index = 1-source_index
            positions = {p[source_index]:p[target_index] for p in pos}
            break
        except (IndexError, ValueError):
            continue
    else:
        return None

    for f in (file_list, transformation_key):
        try:
            os.remove(f)
        except OSError:
            pass

    return positions

def setup_clustering(job, sfam_id, pdbFileStoreID, resoluFileStoreID, jobStoreIDs=None, structure_map=True):
    from Bio import SeqIO
    work_dir = job.fileStore.getLocalTempDir()
    in_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
    out_store = IOStore.get("aws:us-east-1:molmimic-clustered-structures")

    # sdoms_file = get_file(job, "PDB.h5", pdbFileStoreID)
    #
    # sdoms = pd.read_hdf(str(sdoms_file), "merged")
    # sdoms = sdoms[sdoms["sfam_id"]==sfam_id]
    # sdoms = sdoms[["pdbId", "chnLett", "sdi", "domNo"]].drop_duplicates().dropna()

    sfam_seq_key = "{0}/{0}.fasta".format(int(sfam_id))
    sfam_struct_key = "{0}/{0}.max_cluster_logs".format(int(sfam_id))

    domain_fasta = os.path.join(work_dir, "{}.fasta".format(int(sfam_id)))
    if out_store.exists(sfam_seq_key) and out_store.exists(sfam_struct_key):
        out_store.read_input_file("{0}/{0}.fasta".format(int(sfam_id)), domain_fasta)
        with open(domain_fasta) as fh:
            domain_ids = {l[1:]:l[1:] for l in fh if l.startswith(">")}
    else:
        #Save all domains to fasta
        domain_ids = {}
        with open(domain_fasta, "w") as fasta:
            if jobStoreIDs is not None:
                for pdb_fname, jobStoreID in jobStoreIDs:
                    pdb_file = os.path.join(work_dir, pdb_fname)
                    job.fileStore.readGlobalFile(jobStoreID, userPath=pdb_file+".tmp")
                    remove_ter_lines(pdb_file+".tmp", pdb_file)
                    with open(pdb_file) as f:
                        pdb2seq(pdb_fname, sfam_id, f, fasta)
                        domain_ids[pdb_fname] = pdb_file
            else:
                for i, key in enumerate(in_store.list_input_directory(str(int(sfam_id)))):
                    if not key.endswith(".pdb"): continue
                    if i%10==0:
                        RealtimeLogger.info("{} {}".format(i, key))
                    fname = os.path.basename(key)
                    try:
                        in_store.read_input_file(key, fname)
                    except (KeyboardInterrupt, SystemExit):
                        raise
                    except Exception as e:
                        continue

                    with open(fname) as f:
                        pdb2seq(fname, sfam_id, f, fasta)
                        domain_ids[fname] = fname

                    try:
                        os.remove(fname)
                    except OSError:
                        pass


        #Order domains by resolution so the ones with the highest resolutions are centroids
        resolu_file = get_file(job, "resolu.txt", resoluFileStoreID)
        resolutions = pd.read_csv(str(resolu_file))

        try:
            pdbs, ids, sequences = list(zip(*[(s.id.split("_", 1)[0].upper(), s.id, str(s.seq)) \
                for s in SeqIO.parse(domain_fasta, "fasta")]))
        except ValueError:
            RealtimeLogger.info("Unable to cluster {}. No PDBs passed the protonation/minization steps.".format(sfam_id))
            out_store.write_output_file(domain_fasta, "{}/all.fasta".format(int(sfam_id)))
            return

        domains = pd.DataFrame({"pdbId":pdbs, "domainId":ids, "sequence":sequences})
        domains = pd.merge(domains, resolutions, how="left", on="pdbId")
        domains["resolution"] = domains["resolution"].fillna(-1.)

        xray = domains[domains["resolution"] >= 0.].sort_values("resolution")
        nmr = domains[(domains["resolution"] < 0.) | (domains["resolution"].isna())]

        with open(domain_fasta, "w") as f:
            for row in it.chain(xray.itertuples(index=False), nmr.itertuples(index=False)):
                print(">{} [resolution={}]\n{}".format(row.domainId, row.resolution, row.sequence), file=f)

        out_store.write_output_file(domain_fasta, sfam_seq_key)

    # if structure_map:
    #     cluster_structures(job, sfam_id, domain_ids.values(), work_dir=work_dir)

    for pctId in (0.8, 0.85, 0.9, 0.95, 0.99):
        # if out_store.exists("{0}/{0}_clusters_{1}.h5".format(int(sfam_id), id)) and \
        #   out_store.exists("{0}/{0}_clusters_{1}.h5.inferred".format(int(sfam_id), id)):
        #     continue
        #job.addChildJobFn(cluster, sfam_id, pctId)
        cluster(job, sfam_id, pctId, work_dir=work_dir)

def start_toil(job, pdbs=None, structure_map=True):
    """Start the workflow to process PDB files"""
    work_dir = job.fileStore.getLocalTempDir()
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")

    #Download PDB info
    sdoms_file = os.path.join(work_dir, "PDB.h5")
    in_store.read_input_file("PDB.h5", sdoms_file)

    #Add pdb info into local job store
    pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

    #Downlaod Resolution
    resolu_file = os.path.join(work_dir, "resolu.csv")
    in_store.read_input_file("resolu.csv", resolu_file)

    #Add resolution info into local job store
    resoluFileStoreID = job.fileStore.writeGlobalFile(resolu_file, cleanup=True)

    #Get all unique superfamilies
    #sdoms = pd.read_hdf(str(sdoms_file), "merged")

    #sfams = sdoms["sfam_id"].drop_duplicates().dropna().tolist()

    #map_job(job, setup_clustering, sfams, pdbFileStoreID, resoluFileStoreID)

    # #Add jobs for each sdi
    setup_clustering(job, "299845", pdbFileStoreID, resoluFileStoreID, pdbs)

    #del sdoms
    os.remove(sdoms_file)
    os.remove(resolu_file)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.clean = "always"
    options.targetTime = 1

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as workflow:
        pdbs = [(os.path.basename(f), workflow.importFile('file://' + f)) for f in glob.glob("/root/ig/*/*.pdb")]
        workflow.start(Job.wrapJobFn(start_toil, pdbs))
