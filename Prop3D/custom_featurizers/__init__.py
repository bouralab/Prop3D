from pathlib import Path
modules = Path(__file__).parent.glob("*.py")
__all__ = [f.stem for f in modules if f.is_file() and not f.stem.startswith("_")]

from typing import Union, TypeVar
import pandas as pd
from toil.job import Job
from Bio.PDB import Structure

StructureType = TypeVar('StructureType', bound='Structure')

class CustomFeaturizer(object):
    def __init__(self, job: Union[Job, None] = None, work_dir: Union[str, None] = None):
        self.job = job
        self.work_dir = os.getcwd() if work_dir is None else work_dir    

    def calculate_prop3D(self, path: str, structure: StructureType) -> tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
        raise NotImplementedError