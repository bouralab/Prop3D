import os
from toil.job import Job
from Prop3D.generate_data.prepare_protein import process_domain
from Prop3D.generate_data.calculate_features import calculate_features

def create_input_files(*pdbs, input_file=None):
    job = Job()
    work_dir = os.getcwd()
    for pdb in pdbs:
        processed_domain, a, b = process_domain(job, pdb, None, work_dir=work_dir, force=True)
        print("CALC feats", processed_domain)
        atom, residue, edges = calculate_features(job, processed_domain, None, work_dir=work_dir)
