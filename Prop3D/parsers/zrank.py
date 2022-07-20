from __future__ import print_function
import os
import tempfile

from Prop3D.parsers.Electrostatics import Pdb2pqr
from Prop3D.parsers.container import Container
from Prop3D.util import safe_remove

from toil.realtimeLogger import RealtimeLogger

class ZRank(Container):
    IMAGE = "docker://edraizen/zrank"
    PARAMETERS = [
        (":refinement", "store_true", ["-R"]), 
        ("list_file", "path:in")]
    LOCAL = ["zrank"]

    def rank(self, complex_path, refinement=False, retry_with_protonatation=True):
        if not isinstance(complex_path, (list, tuple)):
            complex_path = [complex_path]

        listfile = self.make_listfile(complex_path)

        try:
            self(list_file=list_file, refinement=refinement)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            pdb2pqr = Pdb2pqr(work_dir=self.work_dir, job=self.job)
            if retry_with_protonatation:
                protonated_complex_file = pdb2pqr.debump_add_h(complex_path)
                listfile = self.make_listfile(protonated_complex_file)
                self(list_file=list_file, refinement=refinement)

        assert os.path.isfile(listfile+".zr.out"), "No output for zrank"

        with open(listfile+".zr.out") as f:
            scores = dict(line.rstrip().split() for line in f)

        if len(complex_path) == 1:
            scores = list(scores.values())[0]

        safe_remove((listfile, listfile.name+".zr.out"))

        return scores

    def make_listfile(self, complex_path):
        listfile = tempfile.NamedTemporaryFile(dir=work_dir, prefix="listfile", suffix=".txt", delete=False)
        for pdb in complex_path:
            print(os.path.basename(pdb), file=listfile)
        listfile.close()

        return listfile.name

# def run_zrank(complex_path, refinement=False, retry_with_protonatation=True,
#   work_dir=None, docker=True, job=None):
#     if work_dir is None:
#         work_dir = os.getcwd()
#
#     _parameters = ["-R"] if refinement else []
#
#     if not isinstance(complex_path, (list, tuple)):
#         complex_path = [complex_path]
#
#     listfile = tempfile.NamedTemporaryFile(dir=work_dir, prefix="listfile", suffix=".txt", delete=False)
#     for pdb in complex_path:
#         print(os.path.basename(pdb), file=listfile)
#     listfile.close()
#
#     needs_retry = None
#     if docker and apiDockerCall is not None and job is not None:
#         parameters = _parameters + [os.path.basename(listfile.name)]
#         try:
#             out = apiDockerCall(job,
#                           image='edraizen/zrank:latest',
#                           working_dir="/data",
#                           volumes={work_dir:{"bind":"/data", "mode":"rw"}},
#                           parameters=parameters)
#             RealtimeLogger.info(out)
#         except (SystemExit, KeyboardInterrupt):
#             raise
#         except Exception as e:
#             needs_retry = e
#     else:
#         cmd = [zrank_path] + _parameters + [os.path.join("/data", listfile.name)]
#         try:
#             subprocess.call(cmd)
#         except (SystemExit, KeyboardInterrupt):
#             raise
#         except Exception as e:
#             needs_retry = e
#
#     if needs_retry is not None:
#         if retry_with_protonatation:
#             protonated_complex_file = run_pdb2pqr(complex_path, whitespace=False,
#                 ff="amber", work_dir=work_dir, docker=docker, job=job)
#             return run_zrank(protonated_complex_file, refinement=refinement,
#                 work_dir=work_dir, retry_with_protonatation=False, docker=docker, job=job)
#         else:
#             raise RuntimeError("Cannot run zrank for {}. Error: {}".format(complex_path, needs_retry))
#
#     assert os.path.isfile(listfile.name+".zr.out"), "No output for zrank"
#
#     with open(listfile.name+".zr.out") as f:
#         scores = dict(line.rstrip().split() for line in f)
#
#     if len(complex_path) == 1:
#         scores = list(scores.values())[0]
#
#     for f in (listfile.name, listfile.name+".zr.out"):
#         try:
#             os.remove(f)
#         except OSError:
#             pass
#
#     return scores
