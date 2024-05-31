import os
import sys
import requests
import traceback
from typing import Union, Any, Optional, Sequence

from toil.job import Job

from Prop3D.parsers.SCWRL import SCWRL
from Prop3D.parsers.MODELLER import MODELLER
from Prop3D.parsers.Electrostatics import Pdb2pqr

from Prop3D.util.iostore import IOStore
from Prop3D.util.cath import download_cath_domain
from Prop3D.util import SubprocessChain, safe_remove
from Prop3D.generate_data.data_stores import data_stores
from Prop3D.util.toil import RealtimeLogger_ as RealtimeLogger
from Prop3D.util.pdb import get_all_chains, PDB_TOOLS, s3_download_pdb

class PrepareProteinError(RuntimeError):
    """A new Error to handle errors during protein preparation

    Parameters
    ----------
    cath_domain : str
        Name of cath domain
    stage : str
        function name where error occured
    message : str
        Error message
    job : toil Job
        Current running job
    errors : list or None
        Orignal error objects
    *args : any
        any xtra args
    **kwds : any 
        ANY XTRA KWDS
    """
    def __init__(self, cath_domain: str, stage: str, message: str, job: Job,  
                 errors: Union[list, None] = None, *args: Any, **kwds: Any) -> None:
        super().__init__(*args, **kwds)
        self.cath_domain = cath_domain
        self.stage = stage
        self.message = message
        self.job = job
        self.errors = errors if isinstance(errors, list) else []

    def __str__(self) -> str:
        """Convert errors into string
        """
        return "Error during {}: {}\nErrors:\n".format(self.stage, self.message,
            "\n".join(map(str, self.errors)))

    def save(self, store: Union[IOStore, None] = None) -> None:
        """Save errors to file in an IOStore

        Parameters
        ----------
        store : IOStore
        """
        if store is None:
            store = data_stores(self.job).prepared_cath_structures
        fail_file = "{}.{}".format(self.cath_domain, self.stage)
        with open(fail_file, "w") as f:
            print(self.message, file=f)
            print(self.errors, file=f)
        store.write_output_file(fail_file, "errors/"+os.path.basename(fail_file))
        safe_remove(fail_file)

def extract_domain(pdb_file: str, cath_domain: str, sfam_id: Optional[str] = None, chain: Optional[str] = None, 
                   rename_chain: Optional[str] = None, striphet: bool = True, rslices: Optional[Sequence[str]] = None, 
                   work_dir: Optional[str] = None) -> tuple[str, Sequence[str]]:
    """Extract a domain from a protein structure and cleans the output to make
    it in standard PDB format. No information in changed or added.

    Prepare a protein structure for use in Prop3D. This will
    0) Unzip if gzipped
    1) Select first model
    2) select chain
    3) Remove all alternate conformations (altlocs)
    4) Remove HETATMS
    5) Clip structure into regions if given
    6) Cleanup PDB file and add TER lines in between gaps

    Parameters
    ----------
    pdb_file : str
        Path to PDB file
    cath_domain : str 
        Name of cath domain (e.g. 1xyzA00)
    sfam_id : str
        Name of CATH sfam domains belong to. Not used, can be None.
    chain : str or None
        Chain ID to use
    rename_chain : str or None
        Chain to renmae old one to
    striphet : bool
        Remove HETATMS, default is True
    rslices : list of str
        Residue regions to slice input structure. Format follows pdb_tool.rslice
    work_dir : str
        Where to save files

    Returns
    -------
    domain_file : str
        Path to newly cleaned file
    prep_steps : list of str
        All commands used to clean the structure
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

def prepare_domain(pdb_file: str, chain: str, cath_domain: str, sfam_id: Union[str,None] = None,
                   perform_cns_min: bool = False, work_dir: Union[str, None] = None, 
                   job: Union[Job, None] = None, cleanup: bool = True) ->  tuple[str, list[str]]:
    """Prepare a single domain for use in Prop3D. This method modifies a PDB
    file by adding missing atoms (SCWRL), residues (MODELLER) and adding hydrogens 
    and energy minimizr (simple debump) with PDB2PQR (ff=parse, ph=propka).

    Note: AssertError raised if pdb_file is gzipped or contains 2+ chains

    Parameters
    ----------
    pdb_file : str
        Path to PDB file of chain to prepare. Must not be gzipped or contain
        more than one chain.
    chain : str or None
        Chain ID to use
    cath_domain : str 
        Name of cath domain (e.g. 1xyzA00)
    sfam_id : str
        Name of CATH sfam domains belong to. Not used, can be None.
    perform_cns_min : bool
        Apply more rigorous minimimization with CNS that with just Pdb2Pqr. Default False.
    work_dir : str
        Where to save files
    job : toil.Job
        The toil job that is used to call this function. Not needed, can be None
    cleanup: bool
        Cleanup/remove all temp files
    
    Returns
    -------
    cleaned_file : str
        Path to newly cleaned/modified structure
    prep_steps : list of str
        All commands used to clean/modify the structure
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
        try:
            protonated_pdb = pdb2pqr.debump_add_h(fixed_pdb, ff="parse", noopt=True)
        except Exception:
            #Completely failed, raise error
            raise PrepareProteinError(cath_domain, "prepare", "Unable to protonate" + \
                "{} using pdb2pqr. Please check pdb2pqr error logs.".format(pdb_file) + \
                " Most likeley reason for failing is that the structure is missing " + \
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

def _process_domain(job: Job, cath_domain: str, cathcode: Union[str, None], cathFileStoreID: Union[str, None] = None, 
                    force_chain: Union[str, None] = None, force_rslices: Union[list[str], None] = None, 
                    force: Union[bool, int] = False, work_dir: Union[str, None] = None, get_from_pdb: bool = False, 
                    cleanup: bool = True, memory="12G", preemptable=True) -> tuple[str, list[str], str, str]:
    """Process domains by extracting them, protonating them, and then minimizing. 
    PrepareProteinError is raised on error
    
    Parameters
    ----------
    job : toil job
        Currently running job
    cath_domain : str
        Name of cath domain
    cathcode : str or None
        Name of superfamily domain belongs to. Can be None.
    cathFileStoreID : str or None
        Path to h5 file on hsds endpoint
    force_chain : str or None
        DEPRECATED. not used. Orginally used to modify chain name
    force_rslices : list of str of None
        DEPRECATED, not used. Residue regions to slice input structure. Format follows pdb_tool.rslice
    force : bool or int
        DEPRECATED, not used. Force processing domain by re-ecxtacting and preparing.
    work_dir : str ot None
        Path to save temp files
    get_from_pdb : bool 
        DEPRECATED. not used
    cleanup : bool
        Cleanup temp failes. Default True
    memory : str
        Human readable memory string, only used if running jobon AWS. Default is '12G'
    preemptable : bool
        Job can be restarted immediately if on spot instance.
    """
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

    return prepared_file, prep_steps, domain_file, local_file

def process_domain(job: Job, cath_domain: str, cathcode: Union[str, None], cathFileStoreID: Union[str, None] = None, 
                   force_chain: Union[str, None] = None, force_rslices: Union[list[str], None] = None, 
                   force: bool = False, work_dir: Union[str, None] = None, get_from_pdb: bool = False, 
                   cleanup: bool = True, memory: str ="12G", preemptable: bool = True) -> tuple[str, list[str], str, str]:
    """'Safe' method to process domains. Calls _precess_domain, but all errors are saved so toil jobs can continue.

    Not recommened anymore, and all errors are actually raised. Uncommon if truly needed.

    Parameters
    ----------
    job : toil job
        Currently running job
    cath_domain : str
        Name of cath domain
    cathcode : str or None
        Name of superfamily domain belongs to. Can be None.
    cathFileStoreID : str or None
        Path to h5 file on hsds endpoint
    force_chain : str or None
        DEPRECATED. not used. Orginally used to modify chain name
    force_rslices : list of str of None
        DEPRECATED, not used. Residue regions to slice input structure. Format follows pdb_tool.rslice
    force : bool or int
        DEPRECATED, not used. Force processing domain by re-ecxtacting and preparing.
    work_dir : str ot None
        Path to save temp files
    get_from_pdb : bool 
        DEPRECATED. not used
    cleanup : bool
        Cleanup temp failes. Default True
    memory : str
        Human readable memory string, only used if running jobon AWS. Default is '12G'
    preemptable : bool
        Job can be restarted immediately if on spot instance.
    """
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
