import os
import re
from molmimic.parsers.superpose import Superpose

from toil.realtimeLogger import RealtimeLogger

class TMAlign(Superpose):
    IMAGE = "docker://edraizen/tmalign-cpp"
    ARG_START="-"
    ENTRYPOINT="/opt/TMalign/TMalign"
    RULES = {"make_file_list": "make_file_list"}
    DETACH=True
    PARAMETERS = [
        (":dir", "make_file_list", ["-dir", "{}", "{}"]),
        (":moving_pdb_file", "path:in", ""),
        (":dir1", "make_file_list", ["-dir1", "{}", "{}"]),
        (":fixed_pdb_file", "path:in", ""),
        (":dir2", "make_file_list", ["-dir2", "{}", "{}"]),
        (":normalize_length", "str", "u"),
        (":average_length", "store_true", ["-a", "T"]),
        (":start_alignment", "path:in", "i"),
        (":input_alignment", "path:in", "I"),
        (":matrix_file", "path:out", ["-m", "{}"]),
        (":d", "str"),
        (":out_file", "path:out:ignore", ["-o", "{}"]),
        (":fast", "store_true"),
        (":circular_permuation", "store_true", "cp"),
        (":suffix", "str"),
        (":atom", "str"),
        (":ter", "str"),
        (":split", "str"),
        (":outfmt", "str"),
        (":byresi", "str"),
        (":TMcut", "str"),
        (":mirror", "store_true", ["-mirror","1"]),
        (":het", "store_true", ["-het","1"]),
        (":infmt1", "str"),
        (":infmt2", "str")
    ]

    def make_file_list(self, key, file_list, format_in_path=True):
        file_list = super().make_file_list(key, file_list, format_in_path=format_in_path, basename=True)
        chain_folder = self.work_dir if self.is_local else "/data/"
        if not chain_folder.endswith("/"):
            chain_folder += "/"
        return [chain_folder, file_list]

    def all_vs_all(self, pdb_list, table_out_file=None, distributed=False, **kwds):
        if int(distributed)>0:
            return self._distribute_all_vs_all(pdb_list, table_out_file=table_out_file, n_jobs=int(distributed), **kwds)
        return self(table_out_file=table_out_file, total=len(pdb_list)**2, dir=pdb_list, **kwds)

    def one_vs_all(self, experiment, pdb_list, table_out_file=None, in_all_vs_all=False, **kwds):
        if in_all_vs_all and "out_file" in kwds:
            kwds["out_file"] += f"_{os.path.splitext(os.path.basename(experiment))[0]}"
        return self(table_out_file=table_out_file, total=len(pdb_list)**2, moving_pdb_file=experiment, dir2=pdb_list, **kwds)

    def all_vs_one(self, pdb_list, experiment, table_out_file=None, **kwds):
        return self(table_out_file=table_out_file, total=len(pdb_list)**2, dir1=pdb_list, fixed_pdb_file=experiment, **kwds)

    def one_vs_one(self, moving_pdb_file, fixed_pdb_file, table_out_file=None, **kwds):
        return self(table_out_file=table_out_file, moving_pdb_file=moving_pdb_file, fixed_pdb_file=fixed_pdb_file, **kwds)

    def __call__(self, table_out_file=None, total=None, include_fixed=True, **kwds):
        output = super().__call__(**kwds)
        if table_out_file is None:
            table_out_file = self.tmpdir()

        print("WORK DIR", self.work_dir)
        with open(table_out_file, "w") as f:
            print("chain1,chain2,rmsd,moving_tm_score,fixed_tm_score,chain1_aln,chain2_aln,aln_info", file=f)
            n_results = 0
            for i, line in enumerate(output):
                try:
                    results = self.get_rmsd_tm(stdout=output, include_fixed=include_fixed)
                except AssertionError:
                    continue
                print(",".join(map(str, results)), file=f)
                if n_results%50==0:
                    print(f"Finished {n_results}/{total}")
                n_results += 1
        return table_out_file

    rmsd_re = re.compile("^Aligned length=.+, RMSD=(.+), Seq_ID=")
    moving_tmscore_re = re.compile("^TM-score=(.+) \(if normalized by length of Chain_1")
    fixed_tmscore_re = re.compile("^TM-score=(.+) \(if normalized by length of Chain_2")
    rmsd_tm_re = re.compile("^Aligned length=.+, RMSD=(.+), TM-score=(.+), ID=")
    ending_re = re.compile("^")
    def get_rmsd_tm(self, include_fixed=False, stdout=None):
        if stdout is None:
            stdout = self.stdout

        chain1 = chain2 = None
        rmsd = tm_score = moving_tm_score = fixed_tm_score = None
        chain1_aln = chain2_aln = aln_info = None

        if isinstance(stdout, str):
            lines = iter(stdout.splitlines())
        else:
            lines = iter(stdout)

        test = ""

        for line in lines:
            test += line
            print("SEARCH LINE:", line)
            if line.startswith("Name of Chain_1"):
                chain1 = line.split(" ")[3].rstrip()
                ch2_line = next(lines)
                test+=ch2_line
                chain2 = ch2_line.split(" ")[3].rstrip()
                continue
            m = self.rmsd_re.match(line)
            if m:
                rmsd = float(m.group(1).strip())
                continue
            m = self.rmsd_tm_re.match(line)
            if m:
                rmsd = float(m.group(1).strip())
                tm_score = float(m.group(2).strip())
                continue
            m = self.moving_tmscore_re.match(line)
            if m:
                moving_tm_score = float(m.group(1).strip())
                continue
            m = self.fixed_tmscore_re.match(line)
            if m:
                fixed_tm_score = float(m.group(1).strip())
                continue
            if line.startswith("(\":\" denotes residue pairs"):
                chain1_aln = next(lines).rstrip()
                aln_info = f"\"{next(lines).rstrip()}\""
                chain2_aln = next(lines).rstrip()
                next(lines)
                break

        assert [chain1, chain2, rmsd, moving_tm_score, fixed_tm_score, chain1_aln,
            chain2_aln, aln_info].count(None) == 0, test

        if include_fixed:
            return chain1, chain2, rmsd, moving_tm_score, fixed_tm_score, chain1_aln, \
                chain2_aln, aln_info

        return chain1, chain2, rmsd, moving_tm_score
