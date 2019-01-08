import os

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
STRUCTURES_0={structures_0}<BR>
STRUCTURES_1={structures_1}<BR>
ANASTRUC_1={anastruc_1}<BR>
cmrest={cmrest}<BR>
surfrest={surfrest}<BR
structures_0={structures_0}<BR>
structures_1={structures_1}<BR>
anastruc_1={anastruc_1}<BR>
rotate180_0={rotate180_0}<BR>
crossdock={crossdock}<BR>
rigidmini={rigidmini}<BR>
randorien={randorien}<BR>
rigidtrans={rigidtrans}<BR>
initiosteps={initiosteps}<BR>
cool1_steps={cool1_steps}<BR>
cool2_steps={cool2_steps}<BR>
cool3_steps={cool3_steps}<BR>
submit_save=Save updated parameters<BR>
</h4><!-- HADDOCK -->
</body>
</html>
"""

def run_haddock(dock_name, work_dir=None, docker=True, toil=False, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    if job:
        job.log(str(os.listdir(work_dir)))
    assert any(os.path.isfile(os.path.join(work_dir, f)) for f in ("new.html", "run.cns"))

    if toil:
        image = 'edraizen/haddock-toil:latest'
        parameters = [
          "--provisioner", "aws",
          "--nodeTypes", "t2.small,t2.small:0.0069",
          "--defaultCores", "1",
          "--maxCores", "1",
          "--maxNodes", "1,99",
          "--maxLocalJobs", "100",
          "--targetTime", "1",
          "--batchSystem", "mesos",
          "--defaultMemory", "997Mi",
          "--defaultDisk", "42121Mi",
          "--logFile", "{}.log".format(dock_name),
          "aws:us-east-1:haddock-{}".format(dock_name),
        ]
    else:
        image = 'edraizen/haddock:latest'
        parameters = []

    if docker and apiDockerCall is not None and job is not None:
        oldcwd = os.getcwd()
        os.chdir(work_dir)
        try:
            out = apiDockerCall(job,
                          image,
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters
                          )
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
            #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
            #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)
        os.chdir(oldcwd)
    else:
        try:
            out = subprocess.check_output([sys.executable, "RunHaddock.py"])
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))
    return out

def dock(int_id, pdb1, chain1, sites1, pdb2, chain2, sites2, structures0=1000,
  structures1=200, anastruc1=200, refine=False, small_refine=False,
  tbl_file=None, settings_file=None, work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    #Write AIR TBL_FILE
    print "TBL FILE", tbl_file
    if tbl_file is None:
        tbl_file = os.path.join(work_dir, "{}.tbl".format(int_id))
        with open(tbl_file, "w") as tbl:
            generate_air_restraints(chain1, sites1, chain2, sites2, tbl)

    parameters = {
        "TBL_FILE": tbl_file,
        "HADDOCK_DIR": "/opt/haddock2.2" if docker else haddock_path,
        "PDB1": os.path.join("/data", os.path.basename(pdb1)) if docker else pdb1,
        "PDB2": os.path.join("/data", os.path.basename(pdb2)) if docker else pdb2,
        "CHAIN1": chain1,
        "CHAIN2": chain2
    }

    if small_refine:
        refine=True

    if refine:
        parameters.update({
            "cmrest": "true",
            "surfrest":  "true",
            "structures_0": 5 if small_refine else structures0,
            "structures_1": 5 if small_refine else structures_1,
            "anastruc_1": 1 if small_refine else anastruc_1,
            "rotate180_0": "false",
            "crossdock": "false",
            "rigidmini": "false",
            "randorien": "false",
            "rigidtrans": "false",
            "initiosteps": 0,
            "cool1_steps": 0,
            "cool2_steps": 0,
            "cool3_steps": 0
        })
    else:
        parameters.update({
            "cmrest": "false",
            "surfrest": "false",
            "rotate180_0": "true",
            "crossdock": "true",
            "rigidmini": "true",
            "randorien": "true",
            "rigidtrans": "true",
            "initiosteps": 500,
            "cool1_steps": 500,
            "cool2_steps": 1000,
            "cool3_steps": 1000
        })


    #Write settings file
    if settings_file is None:
        settings_file = os.path.join(work_dir, "new.html")
        with open(settings_file, "w") as settings:
            start_file.format(**parameters)

    #Run haddock first to create run.cns
    try:
        print run_haddock(int_id, work_dir=work_dir, job=job)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        raise
        if docker:
            #Retry with subprocess
            return dock(int_id, pdb1, chain1, sites1, pdb2, chain2, sites2,
                tbl_file=tbl_file, settings_file=tbl_file, work_dir=work_dir, docker=False, job=job)

    run_dir = os.path.join(work_dir, "run1")
    print os.listdir(work_dir)
    assert os.path.isdir(run_dir)

    #Run haddock again in the run direcotry to dock
    run_haddock(int_id, work_dir=run_dir, job=job)

def generate_air_restraints(chain1, binding_site1, chain2, binding_site2, out):
    sites = [(chain1, binding_site1), (chain2, binding_site2)]

    for i, (molchain, molsites) in enumerate(sites):
        intchain, intsites = sites[1-i]
        for j, r1 in enumerate(molsites):
            print >> out, "assign ( resid {0} and segid {1})".format(r1, molchain)
            print >> out, "       ("
            for k, r2 in enumerate(intsites):
                print >> out, "        ( resid {0} and segid {1})".format(r2, intchain)
                if k<len(intsites)-1:
                    print >> out, "     or"
            print >> out, "       )  2.0 2.0 0.0"
            if j<len(molsites)-1:
                print >> out, "!"
