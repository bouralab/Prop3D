try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess
    haddock_path = os.path.dirname(os.path.dirname(subprocess.check_output(["which", "RunHaddock.py"])))

start_file = """<html>
<head>
<title>HADDOCK - start</title>
</head>
<body bgcolor=#ffffff>
<h2>Parameters for the start:</h2>
<BR>
<h4><!-- HADDOCK -->
AMBIG_TBL={TBL_FILE}<BR>
HADDOCK_DIR={HADDOCK_DIR}<BR>
N_COMP=2<BR>
PDB_FILE1={PDB1}<BR>
PDB_FILE2={PDB2}<BR>
PROJECT_DIR={WORKDIR}<BR>
PROT_SEGID_1={CHAIN1}<BR>
PROT_SEGID_2={CHAIN2}<BR>
RUN_NUMBER=1<BR>
submit_save=Save updated parameters<BR>
</h4><!-- HADDOCK -->
</body>
</html>
"""

def run_modeller(work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    assert any(os.path.isfile(os.path.join(work_dir, f)) for f in ("new.html", "run.cns"))

    if docker and apiDockerCall is not None and job is not None:
        try:
            apiDockerCall(job,
                          image='edraizen/haddock:latest',
                          working_dir=work_dir,
                          parameters=[])
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
            #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
            #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)

    else:
        try:
            subprocess.check_output([sys.executable, "RunHaddock.py"])
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

def dock(int_id, pdb1, chain1, sites1, pdb2, chain2, sites2, work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    #Write AIR TBL_FILE
    tbl_file = os.path.join(work_dir, "{}.tbl".format(int_id))
    with open(tbl_file, "w") as tbl:
        generate_air_restraints(chain1, sites1, chain2, sites2, out)

    #Write settings file
    settings_file = os.path.join(work_dir, "new.html")
    with open(settings_file, "w") as settings:
        start_file.format(
            TBL_FILE=tbl_file,
            HADDOCK_DIR="/opt/haddock" if docker else haddock_path,
            PDB1=os.path.join("/data", os.path.basename(pdb1)) if docker else pdb1,
            PDB2=os.path.join("/data", os.path.basename(pdb2)) if docker else pdb2,
            CHAIN1=chain1, CHAIN2=chain2
            )

    #Run haddock first to create run.cns
    try:
        run_haddock()
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        if docker:
            #Retry with subprocess
            return dock(int_id, pdb1, chain1, sites1, pdb2, chain2, sites2,
                work_dir=work_dir, docker=False, job=job)

    #Run haddock again to dock
    run_haddock()

def generate_air_restraints(chain1, binding_site1, chain2, binding_site2, out):
    sites = [(chain1, binding_site1), (chain2, binding_site2)]

    for i, (molchain, molsites) in enumerate(sites):
        intchain, intsites = sites[1-i]
        for r1 in molsites:
            print >> out, "assign ( resid {:3d} and segid {})".format(r1, molchain)
            print >> out, "       ("
            for j, r2 in enumerate(intsites):
                print >> out, "        ( resid {:3d} and segid {})".format(r2, intchain)
                if j<len(binding_site2)-1:
                    print >> out, "     or"
