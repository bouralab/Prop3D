
def start_toil(option, use_data=False):
    from molmimic.generate_data import download_data
    from molmimic.generate_data import convert_mmdb_to_pdb
    from molmimic.generate_data import get_structural_interactome
    from molmimic.generate_data import calculate_bsa
    from molmimic.generate_data import prepare_protein
    from molmimic.generate_data import calculate_features
    from molmimic.generate_data import filter_dataset

    if use_data:
        data = Job.wrapJobFn(download_data.start_toil).encapsulate()
        mmdb2pdb = data.addFollowOnJobFn(convert_mmdb_to_pdb.start_toil).encapsulate()
    else:
        mmdb2pdb = Job.wrapJobFn(convert_mmdb_to_pdb.start_toil).encapsulate()

    interactome = mmdb2pdb.addChildJobFn(get_structural_interactome.start_toil).encapsulate()
    bsa = interactome.addFollowOnJobFn(calculate_bsa.start_toil).encapsulate()

    prep_protein = mmdb2pdb.addChildJobFn(prepare_protein.start_toil).encapsulate()
    features = mmdb2pdb.addFollowOnJobFn(calculate_features.start_toil).encapsulate()

    filter = mmdb2pdb.addFollowOnJobFn(filter_dataset.start_toil).encapsulate()

    with Toil(options) as toil:
        toil.start(mmdb2pdb if not use_data else data)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"

    start_toil(options)
