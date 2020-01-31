from __future__ import print_function
import os
import sys
import tarfile
import shutil
import pandas as pd
import tempfile

from toil.realtimeLogger import RealtimeLogger

from molmimic.util.pdb import PDB_TOOLS
from molmimic.util import SubprocessChain

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
waterrefine={waterrefine}<BR>
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

def run_haddock(dock_name, setup=False, work_dir=None, docker=True, toil=False, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

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
                          parameters=parameters,
                          detach=True
                          )
            # if not setup:
            for line in out.logs(stream=True):
                # stream=True makes this loop blocking; we will loop until
                # the container stops and there is no more output.
                RealtimeLogger.info(line)
            # else:
            #     job.log(out)

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
    #job.log(out)
    return out

def analyze_haddock(analysis_dir, docker=True, job=None):
    if docker and apiDockerCall is not None and job is not None:
        oldcwd = os.getcwd()
        os.chdir(analysis_dir)
        try:
            out = apiDockerCall(job,
                          'edraizen/haddock:latest',
                          entrypoint="csh",
                          working_dir="/data",
                          volumes={analysis_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=["/opt/haddock2.2/tools/ana_structures.csh"],
                          )
            job.log(out)

        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            if "RMSD: Undefined variable" not in str(e):
                raise
        os.chdir(oldcwd)
    else:
        try:
            oldcwd = os.getcwd()
            os.chdir(analysis_dir)
            out = subprocess.check_output(["csh", "/opt/haddock2.2/tools/ana_structures.csh"])
            os.chdir(oldcwd)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            if "RMSD: Undefined variable" not in str(e):
                raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    job.log("ANALYSIS_DIR: {}".format(os.listdir(analysis_dir)))
    results = pd.read_table(
        os.path.join(analysis_dir, "structures_haddock-sorted.stat"),
        nrows=1,
        delim_whitespace=True
    )
    results.columns = ["haddock_"+c for c in results.columns]
    results = results.iloc[0]

    return results

def dock(int_id, pdb1, chain1, sites1, pdb2, chain2, sites2, structures0=1000,
  structures1=200, waterrefine=200, anastruc1=200, refine=False, small_refine=False,
  tbl_file=None, settings_file=None, work_dir=None, clean_docked_file=True,
  docker=True, cores=2, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    #Write AIR TBL_FILE
    if tbl_file is None:
        tbl_file = os.path.join(work_dir, "{}.tbl".format(int_id))
        with open(tbl_file, "w") as tbl:
            generate_air_restraints(chain1, sites1, chain2, sites2, tbl)

    job.log("Wrote tbl file")
    with open(tbl_file) as f:
        job.log(f.read())

    parameters = {
        "TBL_FILE": os.path.join("/data", os.path.basename(tbl_file)) if docker else tbl_file,
        "HADDOCK_DIR": "/opt/haddock2.2" if docker else haddock_path,
        "WORKDIR":"/data",
        "PDB1": os.path.join("/data", os.path.basename(pdb1)) if docker else pdb1,
        "PDB2": os.path.join("/data", os.path.basename(pdb2)) if docker else pdb2,
        "CHAIN1": chain1,
        "CHAIN2": chain2,
        "cpunumber_1": cores
    }

    if refine or small_refine:
        job.log("REFINING INTERFACE")
        parameters.update({
            "cmrest": "true",
            "surfrest":  "true",
            "structures_0": 2 if small_refine else structures0,
            "structures_1": 2 if small_refine else structures1,
            "waterrefine": 2 if small_refine else waterrefine,
            "anastruc_1": 2 if small_refine else anastruc1,
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
            "structures_0": structures0,
            "structures_1": structures1,
            "waterrefine": waterrefine,
            "anastruc_1": anastruc1,
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

    if job is not None:
        job.log(str(parameters))

    #Write settings file
    if settings_file is None:
        settings_file = os.path.join(work_dir, "new.html")
        with open(settings_file, "w") as settings:
            settings.write(start_file.format(**parameters))

    with open(settings_file) as f:
        job.log(f.read())

    #Run haddock first to create run.cns
    try:
        run_haddock(int_id, setup=True, work_dir=work_dir, job=job)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        raise
        if docker:
            #Retry with subprocess
            return dock(int_id, pdb1, chain1, sites1, pdb2, chain2, sites2,
                tbl_file=tbl_file, settings_file=tbl_file, work_dir=work_dir, docker=False, job=job)

    run_dir = os.path.join(work_dir, "run1")
    job.log("DIR AFTER SETUP: {}".format(os.listdir(work_dir)))
    assert os.path.exists(run_dir)
    job.log("DONE")

    #Run haddock again in the run direcotry to dock
    run_haddock(int_id, work_dir=run_dir, job=job)

    #complex_file = os.path.join(run_dir, "structures", "water", "analysis", "molmimicfit_1.pdb")
    #assert os.path.isfile(complex_file)

    water_dir = os.path.join(run_dir, "structures", "it1", "water")
    results = analyze_haddock(water_dir, job=job)

    _complex_file = os.path.join(water_dir, results["haddock_#struc"])
    complex_file = os.path.join(work_dir, "{}.min.pdb".format(int_id))

    _orig_file = os.path.join(run_dir, "begin", "molmimic_1.pdb")
    orig_file = os.path.join(work_dir, "{}.min.pdb".format(int_id))

    for pre_f, post_f in [(_complex_file, complex_file), (_orig_file, orig_file)]:
        assert os.path.isfile(pre_f)

        if clean_docked_file:
            cmds = [
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_segxchain.py"), pre_f],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
            ]
        else:
            cmds = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_segxchain.py"), pre_f]]
        with open(post_f, "w") as f:
            SubprocessChain(cmds, f)

        assert os.path.isfile(post_f)

    results_zip = os.path.join(work_dir, "{}.tar.gz".format(int_id))
    tar = tarfile.open(results_zip, "w:gz")
    tar.add(os.path.join(run_dir, "structures"), arcname=str(int_id))
    tar.close()

    shutil.rmtree(run_dir)

    return orig_file, complex_file, results, results_zip

def generate_air_restraints(chain1, binding_site1, chain2, binding_site2, out):
    sites = [(chain1, binding_site1), (chain2, binding_site2)]

    for i, (molchain, molsites) in enumerate(sites):
        intchain, intsites = sites[1-i]
        for j, r1 in enumerate(molsites):
            print("assign ( resid {0} and segid {1})".format(r1, molchain), file=out)
            print("       (", file=out)
            for k, r2 in enumerate(intsites):
                print("        ( resid {0} and segid {1})".format(r2, intchain), file=out)
                if k<len(intsites)-1:
                    print("     or", file=out)
            print("       )  2.0 2.0 0.0", file=out)
            if j<len(molsites)-1:
                print("!", file=out)

def score_complex(pdb_file, chain, iteration=None, work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    if docker and apiDockerCall is not None and job is not None:
        if not os.path.abspath(os.path.dirname(pdb_file)) == os.path.abspath(work_dir):
            shutil.copy(pdb_file, work_dir)
        try:
            parameters = ["score", os.path.basename(pdb_file), "--chain"]+list(chain)
            if isinstance(iteration, int) and iteration in range(3):
                parameters += ["--iteration", str(iteration)]
            score = apiDockerCall(job,
                          "edraizen/haddock:latest",
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise

    results = pd.read_table(
        os.path.join(work_dir, pdb_file+".haddock-sorted.stat"),
        nrows=1,
        delim_whitespace=True
    )
    results.columns = ["haddock_"+c for c in results.columns]
    results = results.iloc[0]

    return results
