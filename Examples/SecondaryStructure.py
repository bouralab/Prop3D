from Bio.PDB.Structure import Structure
from Bio.PDB.DSSP import DSSP
import h5pyd

from Prop3D.parsers.container import containCall

from Prop3D.util import safe_remove
from Prop3D.util.cath import run_cath_hierarchy_h5


def calculate_secondary_structure(job, cathcode, domain_pdb):
    raise NotImplementedError
    #Not working 

    work_dir = job.fileStore.getLocalTempDir()

    #Call DSSP
    pdb_name = os.path.splitext(os.path.basename(domain_pdb))[0]
    dssp_file = f"{domain_pdb}.dssp"
    containerCall(job,
        image="docker://edraizen/dssp"
        working_dir="/data",
        volumes={work_dir:{"bind":"/data", "mode":"rw"}},
        parameters=["-i", domain_pdb, "-o", dssp_file],
    )

    #Load DSSP results
    structure = Structure(pdb_name, domain_pdb)
    model = structure[0]
    dssp = DSSP(model, dssp_file, file_type='DSSP')
    keys = list(dssp.keys())

    with h5pyd.File("/home/ed4bu/cath-dssp.h5", "a", use_cache=False) as f:
        group = f.require_group(cathcode)
        group.create_dataset()

if __name__ == "__main__":
    raise NotImplementedError

    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()

    #Local
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.maxLocalJobs = 20

    with Toil(options) as workflow:
        #Traverse CATH Hierarchy creating child jobs for each subsidiary level,
        #bottoming out at each domain, running specified function (python or containerized)
        root = Job.wrapJobFn(cath_hierarchy_runner,
            cath_endpoint="hdf://uvaarc01.virginia.edu/bournelab/cath.h5",
            result_endpoint="/home/ed4bu/cath-ss.h5",
            level="domain",
            run=calculate_secondary_structure,
            result_type="pdb", #Can be 'pdb' for path to pdb file,
                               # 'DistributedStructure' to use HSDS structure,
                               # or 'Structure' to use pdb+BioPython structure
        )
