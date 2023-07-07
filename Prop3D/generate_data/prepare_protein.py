import os
import re
import sys
import glob
import gzip
import time
import shutil
import logging
import requests
import traceback
import subprocess
import itertools as it
from collections import defaultdict
from multiprocessing.pool import ThreadPool


import pandas as pd
from Bio import SeqIO
from toil.job import JobFunctionWrappingJob

from Prop3D.parsers.SCWRL import SCWRL
from Prop3D.parsers.MODELLER import MODELLER
from Prop3D.parsers.Electrostatics import Pdb2pqr

from Prop3D.util.toil import map_job
from Prop3D.util.hdf import get_file
from Prop3D.util.iostore import IOStore
from Prop3D.util.cath import download_cath_domain
from Prop3D.util import SubprocessChain, safe_remove
from Prop3D.generate_data.data_stores import data_stores
from Prop3D.util.toil import RealtimeLogger_ as RealtimeLogger
from Prop3D.util.pdb import get_first_chain, get_all_chains, PDB_TOOLS, s3_download_pdb


class PrepareProteinError(RuntimeError):
    def __init__(self, cath_domain, stage, message, job,  errors=None, *args, **kwds):
        super().__init__(*args, **kwds)
        self.cath_domain = cath_domain
        self.stage = stage
        self.message = message
        self.job = job
        self.errors = errors if isinstance(errors, list) else []

    def __str__(self):
        return "Error during {}: {}\nErrors:\n".format(self.stage, self.message,
            "\n".join(map(str, self.errors)))

    def save(self, store=None):
        if store is None:
            store = data_stores(self.job).prepared_cath_structures
        fail_file = "{}.{}".format(self.cath_domain, self.stage)
        with open(fail_file, "w") as f:
            print(self.message, file=f)
            print(self.errors, file=f)
        store.write_output_file(fail_file, "errors/"+os.path.basename(fail_file))
        safe_remove(fail_file)

def extract_domain(pdb_file, cath_domain, sfam_id, chain=None, rename_chain=None,
  striphet=True, rslices=None, work_dir=None):
    """Extract a domain from a protein structure and cleans the output to make
    it in standard PDB format. No information in changed or added

    Prepare a protein structure for use in Prop3D. This will
    0) Unzip if gzipped
    1) Cleanup PDB file and add TER lines in between gaps
    2) Remove HETATMS
    3) Split Chains into separate files
    4) Each chain is then protonated and minimized. See 'run_single_chain'
    for more info

    Parameters
    ----------
    pdb_file : str
        Path to PDB file
    chain : str or None
        Chain ID to split out, protonate, and mininize. If None (default),
        all chains will be split out, protonated, and mininized.

    Returns
    -------
    If no chain was specified a list of all chains is returned with the
    path to the prepared file, the pdb code, and chain. If a chain is
    specified, 3 values are returned: the path to the prepared file,
    the pdb code, and chain. See 'run_single_chain' for info.
    """
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    if work_dir is None:
        work_dir = os.getcwd()

    #Name out output domain file
    domain_file = os.path.join(work_dir, "{}.pdb".format(cath_domain))

    files_to_remove = []

    if pdb_file.endswith(".gz"):
        #Open gzip file
        input = domain_file+".full"
        with gzip.open(pdb_file, 'rt') as zipf, open(input, "w") as pdbf:
            pdbf.write(zipf.read())

        #Remove ungzipped file at end
        files_to_remove.append(input)
    else:
        input = pdb_file

    #Make sure input is not empty
    with open(input) as f:
        if f.read() == "":
            raise RuntimeError("Error processing PDB: {}".format(input))

    all_chains = list(get_all_chains(input))
    if chain is None:
        if sfam_id is not None:
            chain = cath_domain[4]
        elif len(all_chains)>0:
            chain = all_chains[0]
        else:
            raise RuntimeError("Error processing chains in PDB: {}".format(input))


    commands = [
        #Pick first model
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", input]
    ]
    prep_steps = ["pdb_selmodel.py -1 {}".format(input)]

    if chain in all_chains:
        #Select desired chain
        commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain)])
        prep_steps.append("pdb_selchain.py -{}".format(chain))
    elif len(all_chains) == 1 and all_chains[0] == " ":
        #No chain specified, rechain
        commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-"+chain[:1]])
        prep_steps.append( "pdb_chain.py -{}".format(chain[:1]))
    else:
        #Get Auth chain instead of PDB chain
        try:
            with requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{cath_domain[:4]}/{chain}") as r:
                result = r.json()
        except:
            raise RuntimeError("Invalid PDB, chain specified ({}) not in chains ({})".format(chain, all_chains))
        try:
            auth_chain = result["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"]
        except KeyError:
            raise RuntimeError("Invalid PDB, chain specified ({}) not in chains ({})".format(chain, all_chains))
        
        if auth_chain in all_chains:
            #Replace auth chain with pdb chain
            commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(auth_chain)])
            commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-"+chain[:1]])
            prep_steps.append("pdb_selchain.py -{}".format(chain))
            prep_steps.append("pdb_chain.py -{}".format(chain[:1]))
        else:
            raise RuntimeError("Invalid PDB, chain specified ({}) not in chains ({})".format(chain, all_chains))

    commands += [
        #Remove altLocs
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
        #Remove HETATMS
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_striphet.py")],
    ]

    prep_steps += ["pdb_delocc.py", "pdb_striphet.py"]

    if rslices is not None:
        #Slice up chain with given ranges
        commands += [[sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+rslices]
        prep_steps += ["pdb_rslice.py {}".format(" ".join(rslices))]

    #Make it tidy
    commands += [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]]
    prep_steps += ["pdb_tidy.py"]

    if rename_chain is not None and (isinstance(rename_chain, str) or \
      (isinstance(rename_chain, bool) and rename_chain)):
        #Rename chain to given name if set in arguments

        if isinstance(rename_chain, bool) and rename_chain:
            new_chain = "-1"
        else:
            new_chain = "-{}".format(rename_chain)

        commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"),
            new_chain])
        prep_steps += ["pdb_chain.py {}".format(new_chain)]

    #Run all commands
    print("Running", prep_steps)
    with open(domain_file, "w") as output:
        SubprocessChain(commands, output)

    #Make sure output domain_file is not empty
    with open(domain_file) as f:
        content = f.read().rstrip()
        if content == "":
            raise RuntimeError("Error processing PDB: {}".format(domain_file))

    safe_remove(files_to_remove)

    return domain_file, prep_steps

def prepare_domain(pdb_file, chain, cath_domain, sfam_id=None,
  perform_cns_min=False, work_dir=None, job=None, cleanup=True):
    """Prepare a single domain for use in Prop3D. This method modifies a PDB
    file by adding hydrogens with PDB2PQR (ff=parse, ph=propka) and minimizing
    using rosetta (lbfgs_armijo_nonmonotone with tolerance 0.001). Finally,
    the domain PDB file is cleaned so that can be parsed by simple PDB parsers.

    Method called during 'run_protein'

    Note: AssertError raised if pdb_file is gzipped or contains 2+ chains

    Parameters
    ----------
    pdb_file : str
        Path to PDB file of chain to prepare. Must not be gzipped or contain
        more than one chain.

    Returns
    -------
    cleaned_minimized_file : str
        Path to prepared PDB file
    pdb : str
        PDB ID
    chain : str
        Chain ID
    domain : 2-tuple
        Structural domain ID and domain Number.
    """
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    if pdb_file.endswith(".gz"):
        raise RuntimeError("Cannot be a gzip archive, try 'extract_domain' instead")

    if work_dir is None:
        work_dir = os.getcwd()

    all_chains = list(get_all_chains(pdb_file))
    num_chains = len(all_chains)
    if not len(get_all_chains(pdb_file)) == 1:
        raise RuntimeError("Must contain only one chain. There are {} chains in {}".format(
            num_chains, pdb_file))
    if chain is None:
        chain = all_chains[0]

    prefix = pdb_file.split(".", 1)[0]

    pdb2pqr = Pdb2pqr(work_dir=work_dir, job=job)
    scwrl = SCWRL(work_dir=work_dir, job=job).fix_rotamers
    modeller = MODELLER(work_dir=work_dir, job=job).remodel_structure

    errors = []
    files_to_remove = []
    for attempt, fixer in enumerate((None, scwrl, modeller)):
        RealtimeLogger.info("Running Fixer: {} {}".format(fixer, attempt))
        try:
            #If failed 1st time, add correct sidechain rotamers (SCWRL)
            #If failed 2nd time, turn CA models into full atom models (MODELLER)
            fixed_pdb = fixer(pdb_file) if fixer is not None else pdb_file
            if fixer is not None:
                #Remove "fixed" pdb file at end
                files_to_remove.append(fixed_pdb)
        except (SystemExit, KeyError):
            raise
        except Exception as e:
            tb = traceback.format_exc()
            errors.append(tb)
            RealtimeLogger.info("Fixer {} failed: {}".format(fixer.__class__.__name__, tb))
            continue

        try:
            #Protonate PDB, Minimize structure, and assign partial charges
            protonated_pdb = pdb2pqr.debump_add_h(fixed_pdb, ff="parse")
            break
        except Exception as error:
            #Failed, try again with different fixer
            tb = traceback.format_exc()
            errors.append(tb)
            RealtimeLogger.info("Pdb2pqr failed: {}".format(tb))


    else:
        #Completely failed, raise error
        raise PrepareProteinError(cath_domain, "prepare", "Unable to protonate" + \
            "{} using pdb2pqr. Please check pdb2pqr error logs.".format(pdb_file) + \
            "Most likeley reason for failing is that the structure is missing " + \
            "too many heavy atoms.", job, errors)

    if perform_cns_min:
        #Remove pdb2pqr file at end
        files_to_remove.append(protonated_pdb)

        #Perform CNS minimization
        minimize = CNSMinimize(work_dir=work_dir, job=job)
        protonated_pdb = minimize(protonated_pdb)

    #Remove pdb2pqr file (or CNS file) at end
    files_to_remove.append(protonated_pdb)

    #Clean
    commands = [
        #Remove header
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_stripheader.py"), protonated_pdb],
        #Rename chain to given chain
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(chain)],
        #Make it tidy
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
    ]

    cleaned_file = prefix+".pdb"
    with open(cleaned_file, "w") as cleaned:
        SubprocessChain(commands, cleaned)

    if cleanup:
        safe_remove(files_to_remove)

    prep_steps = [fixer.__class__.__name__.rsplit(".")[-1]] if fixer is not None \
        else []
    prep_steps += ["pdb2pqr", "pdb_stripheader.py", "pdb_chain.py -{}".format(chain),
        "pdb_tidy.py"]

    return cleaned_file, prep_steps

def _process_domain(job, cath_domain, cathcode, cathFileStoreID=None, force_chain=None,
  force_rslices=None, force=False, work_dir=None, get_from_pdb=False, cleanup=True,
  memory="12G", preemptable=True):
    """Process domains by extracting them, protonating them, and then minimizing.
    PrepareProteinError is raised on error"""
    if work_dir is None:
        if job and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    files_to_remove = []

    cath_key = ""

    if os.path.isfile(cath_domain):
        domain_file = cath_domain
        local_file = True
        chain = None
        seq_range = None
    else:
        local_file = False

        if cathcode is None and "/" in cath_domain:
            cath_domain, seq_range = cath_domain.rsplit("/", 1)
            seq_range = [r.replace("-", ":") for r in seq_range.split("_")]
        else:
            seq_range = None

        if seq_range is not None:
            seq_range = seq_range.split("_")

        try:
            assert len(cath_domain)==7
            #Get cath domain file
            cath_key = "{}/{}.pdb".format(cathcode.replace(".", "/"), cath_domain)
            if cathcode is not None and data_stores(job).prepared_cath_structures.exists(cath_key) and \
              data_stores(job).prepared_cath_structures.get_size(cath_key)>0:
                return None, None, None, False

            #Download cath domain from s3 bucket or cath api
            domain_file = download_cath_domain(job, cath_domain, cathcode, work_dir=work_dir)
            chain = cath_domain[4]
            if cathcode is not None:
                files_to_remove.append(domain_file)
            assert os.path.isfile(domain_file), "Domain file not found: {}".format(domain_file)
        except (SystemExit, KeyboardInterrupt):
            raise
        except (KeyError, AssertionError):
            if cathcode is None:
                #If not using CATH, process pdb files from pdb or alphafold

                try:
                    domain_file, chain, file_type = s3_download_pdb(cath_domain,
                        work_dir=work_dir, job=job)
                    local_file = True
                except (SystemExit, KeyboardInterrupt):
                    raise
                except:
                    #Cannot download
                    tb = traceback.format_exc()
                    raise PrepareProteinError(cath_domain, "s3download", tb, job)

            else:
                #Cannot download
                tb = traceback.format_exc()
                raise PrepareProteinError(cath_domain, "download", tb, job)
        except Exception as e:
            tb = traceback.format_exc()
            raise PrepareProteinError(cath_domain, "unk_download", tb, job)

    if local_file or not data_stores(job).prepared_cath_structures.exists(cath_key+".raw") or \
      data_stores(job).prepared_cath_structures.get_size(cath_key+".raw")==0:
        #Extract domain; cleaned but atomic coordinates not added or changed
        try:
            domain_file, prep_steps = extract_domain(domain_file, cath_domain, cathcode,
                chain=chain, rslices=seq_range, work_dir=work_dir)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            tb = traceback.format_exc()
            raise PrepareProteinError(cath_domain, "extract", tb, job)

        if not local_file:
            #Write raw domain file to store
            data_stores(job).prepared_cath_structures.write_output_file(domain_file, cath_key+".raw")

        #Write preperation steps
        prep_steps_file = os.path.join(work_dir, "{}.raw.prep".format(cath_domain))
        with open(prep_steps_file, "w") as fh:
            for step in prep_steps:
                print(step, file=fh)

        if not local_file:
            data_stores(job).prepared_cath_structures.write_output_file(prep_steps_file,
                cath_key+".raw.prep")

            files_to_remove.append(domain_file)
            files_to_remove.append(prep_steps_file)

        RealtimeLogger.info("Finished extracting domain: {}".format(domain_file))

    chain = cath_domain[4] if not local_file else None

    #Protonate and minimize, raises PrepareProteinError on error
    prepared_file, prep_steps = prepare_domain(domain_file, chain, cath_domain,
        work_dir=work_dir, job=job)

    if os.path.getsize(prepared_file) == 0:
        raise PrepareProteinError(cath_domain, "empty_file", "", job)

    if not local_file:
        #Write prepared domain file to store
        data_stores(job).prepared_cath_structures.write_output_file(prepared_file, cath_key)
        files_to_remove.append(prepared_file)
    RealtimeLogger.info("Finished preparing domain: {}".format(domain_file))

    #Write preperation steps
    prep_steps_file = os.path.join(work_dir, "{}.prep".format(cath_domain))
    with open(prep_steps_file, "w") as fh:
        for step in prep_steps:
            print(step, file=fh)

    if not local_file:
        data_stores(job).prepared_cath_structures.write_output_file(prep_steps_file,
            cath_key+".prep")
        files_to_remove.append(prep_steps_file)

    if cleanup:
        safe_remove(files_to_remove)

    print("RETURN?", local_file)

    return prepared_file, prep_steps, domain_file, local_file

def process_domain(job, cath_domain, cathcode, cathFileStoreID=None, force_chain=None,
  force_rslices=None, force=False, work_dir=None, get_from_pdb=False, cleanup=True,
  memory="12G", preemptable=True):
    try:
        try:
            return _process_domain(job, cath_domain, cathcode, cathFileStoreID=cathFileStoreID,
                force_chain=force_chain, force_rslices=force_rslices, force=force,
                work_dir=work_dir, get_from_pdb=get_from_pdb, cleanup=cleanup,
                memory=memory, preemptable=preemptable)
        except:
            raise
        # except (SystemExit, KeyboardInterrupt):
        #     raise
        # except PrepareProteinError as e:
        #     e.save()
        # except:
        #     tb = traceback.format_exc()
        #     raise PrepareProteinError(cath_domain, "unk_error", tb, job)
    except (SystemExit, KeyboardInterrupt):
        raise
    except PrepareProteinError as e:
        raise
        e.save()

    print("HERE")

def convert_pdb_to_mmtf(job, sfam_id, jobStoreIDs=None, clustered=True, preemptable=True):
    raise NotImplementedError()

    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    clustered = "clustered" if clustered else "full"

    pdb_path = os.path.join(work_dir, "pdb")
    if not os.path.isdir(pdb_path):
        os.makedirs(pdb_path)

    #Download all with same sfam
    if jobStoreIDs is None:
        in_store = IOStore.get("{}:Prop3D-{}-structures".format(prefix), clustered)
        for f in in_store.list_input_directory(sfam_id):
            if f.endswith(".pdb"):
                in_store.read_input_file(f, os.path.join(work_dir, f))
    else:
        for jobStoreID in jobStoreIDs:
            job.fileStore.readGlobalFile(fileStoreID, userPath=pdb_path)

    PdbToMmtfFull(pdb_path, mmtf_path, work_dir=work_dir, job=job)

    out_store = IOStore.get("{}:Prop3D-{}-mmtf".format(prefix, clustered))
    out_store.write_output_directory(mmtf_path, sfam_id)

def create_data_loader(job, sfam_id, preemptable=True):
    """Create H5 for Molmimic3dCNN to read

    Note: move this somewhere else
    """
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]

    pdb_path = os.path.join(work_dir, "pdb")
    if not os.path.isdir(pdb_path):
        os.makedirs(pdb_path)

    id_format = re.compile("^([A-Z0-9]{4})_([A-Za-z0-9]+)_sdi([0-9]+)_d([0-9]+)$")

    #Get all with keys same sfam, but do not download

    in_store = IOStore.get("{}:Prop3D-clustered-structures".format(prefix))
    keys = [id_format.match(f).groups() for f in in_store.list_input_directory(sfam_id) \
        if f.endswith(".pdb") and id_format.match(f)]

    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)))
    clusters_file = os.path.join(pdb_path, "{}_nr.fasta".format(int(sfam_id)))


    try:
        pdb, chain, sdi, domain = list(zip(*[id_format.match(seq.id[:-2]).groups() \
            for s in SeqIO.parse(clusters_file, "fasta")]))
    except ValueError:
        RealtimeLogger.info("Unable to create data loading file for {}.".format(sfam_id))
        return

    domains = pd.DataFrame({"pdb":pdb, "chain":chain, "domNo":domain, "sdi":sdi})

    data_loader = os.path.join(pdb_path, "{}.h5".format(int(sfam_id)))
    domains.to_hdf(str(data_loader), "table", complevel=9, complib="bzip2")

def copy_pdb_h5(job, path_or_pdbFileStoreID):
    if os.path.isfile(path_or_pdbFileStoreID):
        return path_or_pdbFileStoreID
    else:
        work_dir = job.fileStore.getLocalTempDir()
        sdoms_file = os.path.join(work_dir, "PDB.h5")

        with job.fileStore.readGlobalFileStream(path_or_pdbFileStoreID) as fs_sdoms, open(sdoms_file, "w") as f_sdoms:
            for line in fs_sdoms:
                f_sdoms.write(line)

        return sdoms_file

def process_sfam(job, sfam_id, pdbFileStoreID, cath, clusterFileStoreID, further_parallize=False, cores=1):
    work_dir = job.fileStore.getLocalTempDir()

    if cath:
        clusters_file = get_file(job, "S100.txt", clusterFileStoreID, work_dir=work_dir)
        all_sdoms = pd.read_csv(clusters_file, delim_whitespace=True, header=None,
            usecols=[0, 1, 2, 3, 4], names=["cath_domain", "C", "A", "T", "H"])
        RealtimeLogger.info("ALL SDOMS {}".format(all_sdoms.head()))
        groups = all_sdoms.groupby(["C", "A", "T", "H"])
        RealtimeLogger.info("Running sfam {}".format(sfam_id))
        sdoms = groups.get_group(tuple(sfam_id))["cath_domain"].drop_duplicates().dropna()
    else:
        sdoms_file = copy_pdb_h5(job, pdbFileStoreID)
        sdoms = pd.read_hdf(str(sdoms_file), "merged")
        sdoms = sdoms[sdoms["sfam_id"]==float(sfam_id)]["sdi"].drop_duplicates().dropna()
    #sdoms = sdoms[:1]
    if further_parallize:
        if False: #cores > 2:
            #Only makes sense for slurm or other bare-matal clsuters
            # setup_dask(cores)
            # d_sdoms = dd.from_pandas(sdoms, npartitions=cores)
            # RealtimeLogger.info("Running sfam dask {}".format(sdoms))
            # processed_domains = d_sdoms.apply(lambda row: process_domain(job, row.sdi,
            #     sdoms_file), axis=1).compute()
            pass
        else:
            map_job(job, process_domain, sdoms, pdbFileStoreID)
    else:
        for sdom in sdoms:
            process_domain(job, sdom, pdbFileStoreID)


def post_process_sfam(job, sfam_id, jobStoreIDs):
    cluster(job, sfam_id, jobStoreIDs)

    if False:
        convert_pdb_to_mmtf(job, sfam_id, jobStoreIDs)
        create_data_loader(job, sfam_id, jobStoreIDs)

def start_toil(job, cath=True):
    """Start the workflow to process PDB files"""
    work_dir = job.fileStore.getLocalTempDir()


    if cath:
        in_store = IOStore.get("aws:us-east-1:Prop3D-cath")
        sdoms_file = os.path.join(work_dir, "cath-domain-description-file-small.h5")
        in_store.read_input_file("cath-domain-description-file-small.h5", sdoms_file)

        #Add pdb info into local job store
        pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

        clusters_file = os.path.join(work_dir, "cath-domain-list-S100.txt")
        in_store.read_input_file("cath-domain-list-S35.txt", clusters_file)

        with open(clusters_file) as f:
            sfams = list(set([tuple(map(int, l.split()[1:5])) for l in f if l and not l.startswith("#")]))

        clusterFileStoreID = job.fileStore.writeGlobalFile(clusters_file)

    else:
        in_store = IOStore.get("aws:us-east-1:Prop3D-ibis")

        #Download PDB info
        sdoms_file = os.path.join(work_dir, "PDB.h5")
        in_store.read_input_file("PDB.h5", sdoms_file)

        #Add pdb info into local job store
        pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

        #Get all unique superfamilies
        sdoms = pd.read_hdf(str(sdoms_file), "merged")

        # skip_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "keep.csv")
        # if os.path.isfile(skip_file):
        #     skip = pd.read_csv(skip_file)
        #     sdoms = sdoms[sdoms["sdi"].isin(skip["sdi"])]
        #     RealtimeLogger.info("SKIPPING {} sdis; RUNIING {} sdis".format(skip.shape[0], sdoms.shape[0]))
        #
        sfams = sdoms["sfam_id"].drop_duplicates().dropna()

        clusterFileStoreID = None
    #sfams = sfams[:1]
    #sfams = ["653504"]

    # max_cores = job.fileStore.jobStore.config.maxCores if \
    #     job.fileStore.jobStore.config.maxCores > 2 else \
    #     job.fileStore.jobStore.config.defaultCores

    max_cores = job.fileStore.jobStore.config.defaultCores
    #Add jobs for each sdi
    job.addChildJobFn(map_job, process_sfam, sfams, pdbFileStoreID, cath, clusterFileStoreID,
        cores=max_cores)

    #Add jobs for to post process each sfam
    #job.addFollowOnJobFn(map_job, post_process_sfam, sfams, pdbFileStoreID,
    #    cores=max_cores)

    del sfams
    safe_remove(sdoms_file)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    detector = logging.StreamHandler()
    logging.getLogger().addHandler(detector)

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
