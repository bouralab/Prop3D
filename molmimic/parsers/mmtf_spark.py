import os

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

def PdbToMmtfFull(pdb_path, mmtf_path, work_dir=None, job=None):
    if apiDockerCall is not None and work_dir is not None and job is not None:
        if not os.path.abspath(os.path.dirname(pdb_path)) == os.path.abspath(work_dir):
            raise RuntimeError("Error pdb_path must be inside work_dir")

        pdb_path = os.path.join("/data", os.path.basename(pdb_path))
        mmtf_path = os.path.join("/data", os.path.basename(mmtf_path))

        parameters = ["edu.sdsc.mmtf.spark.io.demos.PdbToMmtfFull", pdb_path,
            mmtf_path]

        apiDockerCall(job,
                      image='edraizen/mmtf-spark:latest',
                      working_dir=work_dir,
                      parameters=parameters)

    else:
        with open(full_pdb_path) as f:
            cx_f = subprocess.check_output("cx", stdin=f)
        cx_f = iter(cx_f.splitlines())
        delCx = None

        chemcomp = os.path.join(work_dir, "pdb", "chemcomp")
        if not os.path.isdir(chemcomp):
            os.makedirs(chemcomp)

        spark_env = os.environ.copy()
        spark_env["PDB_DIR"] = chemcomp
        spark_env["PDB_CACHE_DIR"] = chemcomp

        if not os.path.isdir(mmtf_path):
            os.makedirs(mmtf_path)

        #Convert PDBs to MMTF-Hadoop Sequence file directory
        try:
            subprocess.call(["spark-submit",
                "--class", "edu.sdsc.mmtf.spark.io.demos.PdbToMmtfFull",
                "{}/target/mmtf-spark-0.2.0-SNAPSHOT.jar".format(os.environ["MMTFSPARK_HOME"]),
                pdb_path, mmtf_path],
                env=spark_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            job.log("Error converting to MMTF: {}".format(e))
            pass
