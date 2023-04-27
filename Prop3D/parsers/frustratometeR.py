import os
import shutil
from pathlib import Path

from Prop3D.parsers.container import Container
from Prop3D.util import safe_remove
from Prop3D.util.pdb import reres_pdb, get_pdb_residues, get_all_chains

from toil.realtimeLogger import RealtimeLogger

class FrustratometeR(Container):
    IMAGE = "docker://edraizen/frustratometer:latest"
    PARAMETERS = [
        ("pdb_dir", "str", ["{}"]),
        (":mode:singleresidue", "str", ["{}"]),
        (":modeller_key:", "str", ["{}"]), 
        ]
    CONTAINER_FILE_PREFIX = "/pdb"

    def run(self, pdb_file, mode="singleresidue", modeller_key="", bioPDB=None, parse=True, clean=False):
        assert mode in ["singleresidue", "configurational", "mutations"]
        assert len(get_all_chains(pdb_file)) == 1, get_all_chains(pdb_file)

        original_dir = Path(pdb_file).parent

        #frustartometeR removes all insertion codes, so rename them
        inscode_mapping = {str(i+1):r for i, r in enumerate(get_pdb_residues(str(pdb_file), use_raw_resi=True))}

        if modeller_key == "":
            modeller_key = os.environ.get("KEY_MODELLER", "")

        run_dir = Path(self.work_dir) / Path(self.tempfile(delete=True)).stem
        run_dir.mkdir(mode=0o777)
        run_dir.chmod(0o777)
        print(run_dir)
        assert run_dir.is_dir()
        # print(str(run_dir/f"{Path(pdb_file).stem}.reres.pdb"))
        #run_dir = Path(self.work_dir)

        pdb_file = reres_pdb(str(pdb_file), updated_pdb=str(run_dir/f"{Path(pdb_file).stem}.pdb"))

        original_work_dir = self.work_dir
        self.work_dir = run_dir

        pdb_file = self.format_in_path(None, pdb_file, move_files_to_work_dir=True, absolute_path=False)
        pdb_file = Path(pdb_file)
        pdb_stem = pdb_file.stem.split('.')[0]
        outfile = run_dir / pdb_file.with_suffix(".done").name / "FrustrationData" / f"{pdb_file.name}_{mode}"

        pdb_dir = str(pdb_file.parent) if not str(pdb_file).startswith("/pdb") else "/pdb"

        try:
            self(pdb_dir="/pdb", mode=mode, modeller_key=modeller_key)
        except RuntimeError:
            #Check if files are present below
            pass

        assert run_dir.is_dir(), run_dir
        assert (run_dir / pdb_file.with_suffix(".done").name).is_dir(), run_dir / pdb_file.with_suffix(".done").name
        assert (run_dir / pdb_file.with_suffix(".done").name / "FrustrationData").is_dir(), run_dir / pdb_file.with_suffix(".done").name / "FrustrationData"
        assert (run_dir / pdb_file.with_suffix(".done").name / "FrustrationData" / f"{pdb_file.name}_{mode}").is_file(), run_dir / pdb_file.with_suffix(".done").name / "FrustrationData" / f"{pdb_file.name}_{mode}"

        self.work_dir = original_work_dir

        if not parse:
            return outfile

        output = self.parse(outfile, inscode_mapping=inscode_mapping)

        if clean:
            shutil.rmtree(str(run_dir))

        return output

    @classmethod
    def parse(cls, outfile, inscode_mapping=None):
        import pandas as pd
        data = pd.read_csv(outfile, sep=" ")
        if inscode_mapping is not None:
            try:
                data["Res"] = data["Res"].astype(str).map(inscode_mapping)
            except KeyError:
                #Configurational
                data["Res1"] = data["Res1"].astype(str).map(inscode_mapping)
                data["Res2"] = data["Res2"].astype(str).map(inscode_mapping)
        return data