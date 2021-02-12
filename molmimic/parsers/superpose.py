import os
import sys
import glob
import re

from molmimic.parsers.container import Container
from molmimic.util.pdb import PDB_TOOLS, extract_chains, get_all_chains, \
    get_atom_lines, update_xyz
from molmimic.util import SubprocessChain

class Superpose(object):
    def get_chain(self, pdb_file, chain, new_chain):
        chain_file = self.tempfile()
        extract_chains(pdb_file, chain, rename=new_chain*len(chain),
            new_file=chain_file)
        return chain_file

    def align(self, fixed_pdb, fixed_chain, moving_pdb, moving_chain, use_aln=None):
        self.fixed_file = self.get_chain(fixed_pdb, fixed_chain, "B")
        self.moving_file = self.get_chain(moving_pdb, moving_chain, "A")
        self.out_prefix = self.tempfile()
        self.matrix_file = self.out_prefix+".matrix.txt" if use_aln is None else use_aln

        #Perform alignment
        if use_aln is None:
            self.stdout = self(moving_pdb_file=self.moving_file,
                fixed_pdb_file=self.fixed_file, out_file=self.out_prefix,
                matrix_file=self.matrix_file)
        else:
            self.stdout = self(moving_pdb_file=self.moving_file,
                fixed_pdb_file=self.fixed_file, out_file=self.out_prefix,
                input_alignment=self.matrix_file)

        self.raw_aligned_pdb = self.out_prefix+self._get_ending()

        return self.get_alignment_stats()

    def extract(self):
        outfile = os.path.join(self.work_dir, "{}.aligned.pdb".format(
            os.path.splitext(os.path.basename(self.moving_file))[0]))

        #Chain A had the the moving_pdb rottrans
        outfile = extract_chains(outfile, "A")

        #Copy the updated XYZ coords into moving_pdb file to ensure chains are correct
        update_xyz(moving_file, self.raw_aligned_pdb, updated_pdb=outfile)

    def get_full_aligned_file(self):
        outfile = os.path.join(self.work_dir, "{}__{}.aligned.pdb".format(
             os.path.splitext(os.path.basename(self.fixed_file))[0],
             os.path.splitext(os.path.basename(self.moving_file))[0]))

        #Remove extraneous lines
        cmd = [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
        with open(self.raw_aligned_pdb, "w") as out:
            subprocess.call(cmd, stdout=out)

    def clean(self):
        for f in glob.glob(self.out_prefix+"*"):
            try:
                os.remove(f)
            except OSError:
                pass

class TMAlign(Container, Superpose):
    IMAGE = "edraizen/mmalign:latest"
    PARAMETERS = [
        ("moving_pdb_file", "path:in", ""),
        ("fixed_pdb_file", "path:in", ""),
        ("out_file", "path:out", "o"),
        (":input_alignment", "path:in", "I"),
        (":matrix_file", "path:out", "m")]
    ARG_START="-"

    def get_alignment_stats(self, stdout=None):
        if stdout is None:
            stdout = self.stdout
        return self.get_rmsd_tm(stdout)

    def get_rmsd_tm(self, stdout=None):
        if stdout is None:
            stdout = self.stdout

        rmsd_re = re.compile("^Aligned length=.+, RMSD=(.+), Seq_ID=")
        tmscore_re = re.compile("^TM-score=(.+) \(if normalized by length of Chain_2\)")
        rmsd_tm_re = re.compile("^Aligned length=.+, RMSD=(.+), TM-score=(.+), ID=")

        rmsd = tm_score = -1.
        lines = iter(stdout.splitlines())
        for line in lines:
            job.log(line.rstrip())
            m = rmsd_re.match(line)
            if m:
                rmsd = float(m.group(1).strip())
                job.log("RMSD is {}".format(rmsd))
                continue
            m = rmsd_tm_re.match(line)
            if m:
                rmsd = float(m.group(1).strip())
                tm_score = float(m.group(2).strip())
                job.log("RMSD is {}".format(rmsd))
                job.log("TM-score is {}".format(tm_score))
            m = tmscore_re.match(line)
            if m:
                tm_score = float(m.group(1).strip())
                job.log("TM-score is {}".format(tm_score))

        return rmsd, tm_score

    def get_aligned_positions(self, stdout=None):
        if stdout is None:
            stdout = self.stdout
        moving_ungapped_to_gapped = []
        target_ungapped_to_gapped = []
        lines = iter(stdout.splitlines())
        for line in lines:
            if line.startswith("(\":\" denotes aligned residue"):
                moving_ungapped_to_gapped = [i for i, aa in enumerate(next(lines)) if aa != "-"]
                next(lines)
                target_ungapped_to_gapped = [i for i, aa in enumerate(next(lines)) if aa != "-"]
                break
        return  moving_ungapped_to_gapped, target_ungapped_to_gapped

    def _get_ending(self):
        return ".sup_all_atm_lig"

class MMAlign(TMAlign, Superpose):
    IMAGE = "edraizen/mmalign:latest"
    PARAMETERS = [
        ("moving_pdb_file", "path:in", ""),
        ("fixed_pdb_file", "path:in", ""),
        ("out_file", "path:out", "o")]
    ARG_START="-"

    def align(self, *args, **kwds):
        output = super(self, MMAlign).align(*args, **kwds)
        self.parse_matrix(self.matrix_file)
        return output

    def parse_matrix(self, matrix_file, stdout=None):
        if stdout is None:
            stdout = self.stdout
        lines = iter(stdout.splitlines())
        for line in lines:
            if "rotation matrix" in line:
                with open(matrix_file, "w") as matfile:
                    for i, mat_line in enumerate(lines):
                        if i>4: break
                        matfile.write(mat_line)
                break
        return matrix_file

    def _get_ending(self):
        return ".sup_all"

class CEAlign(Container, Superpose):
    IMAGE = "edraizen/ce:latest"
    PARAMETERS = [
        ("fixed_pdb_file", "path:in", "file1"),
        ("moving_pdb_file", "path:in", "file2"),
        "-outputPDB",
        ("out_file", "path:out", ["-outFile", "{}"])]
        ARG_START="--"

    def align(self, *args, **kwds):
        output = super(self, CEAlign).align(*args, **kwds)
        self.matrix_file = None
        return output

    def get_alignment_stats(self, stdout=None):
        if stdout is None:
            stdout = self.stdout
        return 1.0

    def _get_ending(self):
        return ""

class Align(TMAlign):
    pass


# def align(fixed_file, fixed_chain, moving_file, moving_chain, method="tmalign", force_alignment=None, extract=True, parse_postions=False, docker=True, work_dir=None, job=None):
#     # if work_dir is None:
    #     work_dir = os.getcwd()

    # import shutil
    # numalign = sum(1 for f in os.listdir("/root") if f.startswith("align"))/3.+1
    # shutil.copy(fixed_file, os.path.join("/root", "align{}_fixed_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))
    # shutil.copy(moving_file, os.path.join("/root", "align{}_moving_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    # _chain1 = tempfile.NamedTemporaryFile(suffix=".pdb", dir=work_dir, delete=False)
    # extract_chains(fixed_file, fixed_chain, rename="B"*len(fixed_chain), new_file=_chain1.name)
    # _chain1.close()
    #
    # job.log("MOVE CHAINS: {} {}".format(get_all_chains(moving_file), moving_chain))
    # job.log("FIX CHAINS: {} {}".format(get_all_chains(fixed_file), fixed_chain))
    #
    # _chain2 = tempfile.NamedTemporaryFile(suffix=".pdb", dir=work_dir, delete=False)
    # extract_chains(moving_file, moving_chain, rename="A"*len(moving_chain), new_file=_chain2.name)
    # _chain2.close()
    #
    # job.log("MOVE CHAIN A: {}".format(next(get_atom_lines(_chain2.name))))
    # job.log("MOVE CHAIN A: {}".format(get_all_chains(_chain2.name)))
    # job.log("MOVE CHAIN B: {}".format(next(get_atom_lines(_chain1.name))))
    # job.log("MOVE CHAIN B: {}".format(get_all_chains(_chain1.name)))
    #
    # shutil.copy(_chain1.name, os.path.join("/root", "align{}_seg_fixed_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))
    # shutil.copy(_chain2.name, os.path.join("/root", "align{}_seg_moving_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    # _outf = tempfile.NamedTemporaryFile(dir=work_dir, delete=False)
    # _outf.close()
    # _outf = _outf.name

    # if method in ["tmalign", "mmalign"]:
    #     image = "edraizen/{}:latest".format(method)
    #     parameters = [os.path.basename(_chain2.name), os.path.basename(_chain1.name),
    #         "-o", os.path.basename(_outf+".sup")]
    #     if method == "tmalign":
    #         if force_alignment is not None:
    #             parameters += ["-I", os.path.basename(force_alignment)]
    #         else:
    #             parameters += ["-m", os.path.basename(_outf+".matrix.txt")]
    # elif method == "ce":
    #     image = "edraizen/ce:latest"
    #     parameters = lambda f, m, o: ["--file1", os.path.basename(_chain1.name),
    #         "--file2", os.path.basename(_chain2.name), "-outputPDB", "-outFile",
    #         os.path.basename(_outf+".sup")]
    #
    # if docker and apiDockerCall is not None and job is not None:
    #     try:
    #         stdout = apiDockerCall(job,
    #                       image,
    #                       working_dir="/data",
    #                       volumes={work_dir:{"bind":"/data", "mode":"rw"}},
    #                       parameters=parameters
    #         )
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except Exception as e:
    #         print(dir(e))
    #         raise Exception(str(e).decode('utf-8').encode("ascii", "ignore"))
    #
    # else:
    #     raise RuntimeError("Only docker works at the moment")
    # job.log("OUTPUT: "+stdout)
    #
    # rmsd_re = re.compile("^Aligned length=.+, RMSD=(.+), Seq_ID=")
    # tmscore_re = re.compile("^TM-score=(.+) \(if normalized by length of Chain_2\)")
    # rmsd_tm_re = re.compile("^Aligned length=.+, RMSD=(.+), TM-score=(.+), ID=")
    #
    # rmsd = tm_score = -1.
    # lines = iter(stdout.splitlines())
    # for line in lines:
    #     job.log(line.rstrip())
    #     m = rmsd_re.match(line)
    #     if m:
    #         rmsd = float(m.group(1).strip())
    #         job.log("RMSD is {}".format(rmsd))
    #         continue
    #     m = rmsd_tm_re.match(line)
    #     if m:
    #         rmsd = float(m.group(1).strip())
    #         tm_score = float(m.group(2).strip())
    #         job.log("RMSD is {}".format(rmsd))
    #         job.log("TM-score is {}".format(tm_score))
    #     m = tmscore_re.match(line)
    #     if m:
    #         tm_score = float(m.group(1).strip())
    #         job.log("TM-score is {}".format(tm_score))
    #     if method=="mmalign" and "rotation matrix" in line:
    #         with open(_outf+".matrix.txt", "w") as matfile:
    #             for i, mat_line in enumerate(lines):
    #                 if i>4: break
    #                 matfile.write(mat_line)
    #     if parse_postions and line.startswith("(\":\" denotes aligned residue"):
    #         moving_ungapped_to_gapped = [i for i, aa in enumerate(next(lines)) if aa != "-"]
    #         next(lines)
    #         target_ungapped_to_gapped = [i for i, aa in enumerate(next(lines)) if aa != "-"]


    # ending = ".sup_all_atm_lig" if method == "tmalign" else ".sup_all"
    # _outfile = _outf+ending

    # job.log("ALL chains: {}".format(get_all_chains(_outfile)))
    #
    # job.log("FIRST LINE: {}".format(next(get_atom_lines(_outfile))))
    #
    #
    # job.log("ALL alinged files: {}".format(os.listdir(work_dir)))
    #
    # shutil.copy(_outfile, os.path.join("/root",
    #     "align{}_raw_aligned_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    # if method == "tmalign":
    #     shutil.copy(_outf+".matrix.txt", os.path.join("/root",
    #         "align{}_matrix_aligned_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    # if extract:
    #     outfile = os.path.join(work_dir, "{}.aligned.pdb".format(
    #         os.path.splitext(os.path.basename(moving_file))[0]))
    # else:
    #     outfile = os.path.join(work_dir, "{}__{}.aligned.pdb".format(
    #          os.path.splitext(os.path.basename(fixed_file))[0],
    #          os.path.splitext(os.path.basename(moving_file))[0]))
    #
    # if extract:
    #     #Chain A had the the moving_pdb rottrans
    #     _outfile = extract_chains(_outfile, "A")
    #
    #     #Copy the updated XYZ coords into moving_pdb file to ensure chains are correct
    #     update_xyz(moving_file, _outfile, updated_pdb=outfile)
    # else:
    #     #Remove extraneous lines
    #     with open(outfile, "w") as out:
    #         subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")], stdout=out)

    # if force_alignment is None:
    #     matrix_file = outfile+".matrix"
    #     assert os.path.isfile(_outf+".matrix.txt")
    #     shutil.move(_outf+".matrix.txt", matrix_file)
    #     assert os.path.isfile(matrix_file)
    # else:
    #     matrix_file = force_alignment
    #
    # job.log("NEW chains: {}".format(get_all_chains(outfile)))

    # for f in glob.glob(os.path.join(work_dir, _outf+"*")):
    #     try:
    #         os.remove(f)
    #     except OSError:
    #         pass

    # assert os.path.isfile(outfile)
    # if not parse_postions:
    #     return outfile, rmsd, tm_score, matrix_file
    # else:
    #     return outfile, rmsd, tm_score, matrix_file, moving_ungapped_to_gapped, target_ungapped_to_gapped
