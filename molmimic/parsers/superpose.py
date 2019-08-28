import os
import sys
import subprocess
import tempfile
import glob
import re

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    haddock_path = os.path.dirname(os.path.dirname(subprocess.check_output(["which", "RunHaddock.py"])))

from molmimic.generate_data.util import PDB_TOOLS, extract_chains, SubprocessChain, get_all_chains, get_atom_lines, update_xyz

def align(fixed_file, fixed_chain, moving_file, moving_chain, method="tmalign", force_alignment=None, extract=True, parse_postions=False, docker=True, work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    import shutil
    numalign = sum(1 for f in os.listdir("/root") if f.startswith("align"))/3.+1
    shutil.copy(fixed_file, os.path.join("/root", "align{}_fixed_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))
    shutil.copy(moving_file, os.path.join("/root", "align{}_moving_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    _chain1 = tempfile.NamedTemporaryFile(suffix=".pdb", dir=work_dir, delete=False)
    extract_chains(fixed_file, fixed_chain, rename="B"*len(fixed_chain), new_file=_chain1.name)
    _chain1.close()

    job.log("MOVE CHAINS: {} {}".format(get_all_chains(moving_file), moving_chain))
    job.log("FIX CHAINS: {} {}".format(get_all_chains(fixed_file), fixed_chain))

    _chain2 = tempfile.NamedTemporaryFile(suffix=".pdb", dir=work_dir, delete=False)
    extract_chains(moving_file, moving_chain, rename="A"*len(moving_chain), new_file=_chain2.name)
    _chain2.close()

    job.log("MOVE CHAIN A: {}".format(next(get_atom_lines(_chain2.name))))
    job.log("MOVE CHAIN A: {}".format(get_all_chains(_chain2.name)))
    job.log("MOVE CHAIN B: {}".format(next(get_atom_lines(_chain1.name))))
    job.log("MOVE CHAIN B: {}".format(get_all_chains(_chain1.name)))

    shutil.copy(_chain1.name, os.path.join("/root", "align{}_seg_fixed_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))
    shutil.copy(_chain2.name, os.path.join("/root", "align{}_seg_moving_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    _outf = tempfile.NamedTemporaryFile(dir=work_dir, delete=False)
    _outf.close()
    _outf = _outf.name

    if method in ["tmalign", "mmalign"]:
        image = "edraizen/{}:latest".format(method)
        parameters = [os.path.basename(_chain2.name), os.path.basename(_chain1.name),
            "-o", os.path.basename(_outf+".sup")]
        if method == "tmalign":
            if force_alignment is not None:
                parameters += ["-I", os.path.basename(force_alignment)]
            else:
                parameters += ["-m", os.path.basename(_outf+".matrix.txt")]
    elif method == "ce":
        image = "edraizen/ce:latest"
        parameters = lambda f, m, o: ["--file1", os.path.basename(_chain1.name),
            "--file2", os.path.basename(_chain2.name), "-outputPDB", "-outFile",
            os.path.basename(_outf+".sup")]

    if docker and apiDockerCall is not None and job is not None:
        try:
            stdout = apiDockerCall(job,
                          image,
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters
            )
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            print(dir(e))
            raise Exception(str(e).decode('utf-8').encode("ascii", "ignore"))

    else:
        raise RuntimeError("Only docker works at the moment")
    job.log("OUTPUT: "+stdout)

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
        if method=="mmalign" and "rotation matrix" in line:
            with open(_outf+".matrix.txt", "w") as matfile:
                for i, mat_line in enumerate(lines):
                    if i>4: break
                    matfile.write(mat_line)
        if parse_postions and line.startswith("(\":\" denotes aligned residue"):
            moving_ungapped_to_gapped = [i for i, aa in enumerate(next(lines)) if aa != "-"]
            next(lines)
            target_ungapped_to_gapped = [i for i, aa in enumerate(next(lines)) if aa != "-"]


    ending = ".sup_all_atm_lig" if method == "tmalign" else ".sup_all"
    _outfile = _outf+ending

    job.log("ALL chains: {}".format(get_all_chains(_outfile)))

    job.log("FIRST LINE: {}".format(next(get_atom_lines(_outfile))))


    job.log("ALL alinged files: {}".format(os.listdir(work_dir)))

    shutil.copy(_outfile, os.path.join("/root",
        "align{}_raw_aligned_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    if method == "tmalign":
        shutil.copy(_outf+".matrix.txt", os.path.join("/root",
            "align{}_matrix_aligned_{}_{}.pdb".format(numalign, fixed_chain, moving_chain)))

    if extract:
        outfile = os.path.join(work_dir, "{}.aligned.pdb".format(
            os.path.splitext(os.path.basename(moving_file))[0]))
    else:
        outfile = os.path.join(work_dir, "{}__{}.aligned.pdb".format(
             os.path.splitext(os.path.basename(fixed_file))[0],
             os.path.splitext(os.path.basename(moving_file))[0]))

    if extract:
        #Chain A had the the moving_pdb rottrans
        _outfile = extract_chains(_outfile, "A")

        #Copy the updated XYZ coords into moving_pdb file to ensure chains are correct
        update_xyz(moving_file, _outfile, updated_pdb=outfile)
    else:
        #Remove extraneous lines
        with open(outfile, "w") as out:
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")], stdout=out)

    if force_alignment is None:
        matrix_file = outfile+".matrix"
        assert os.path.isfile(_outf+".matrix.txt")
        shutil.move(_outf+".matrix.txt", matrix_file)
        assert os.path.isfile(matrix_file)
    else:
        matrix_file = force_alignment

    job.log("NEW chains: {}".format(get_all_chains(outfile)))

    for f in glob.glob(os.path.join(work_dir, _outf+"*")):
        try:
            os.remove(f)
        except OSError:
            pass

    assert os.path.isfile(outfile)
    if not parse_postions:
        return outfile, rmsd, tm_score, matrix_file
    else:
        return outfile, rmsd, tm_score, matrix_file, moving_ungapped_to_gapped, target_ungapped_to_gapped
