from molmimic.util.iostore import IOStore
from itertools import groupby

def fix_s3():
    eppic_store = IOStore.get("aws:us-east-1:molmimic-eppic-service")
    all_files = list(eppic.list_input_directory("sequences"))
    for orig_file, similar_files in groupby(sorted(all_files), key=lambda k: k.split(".json")[0]):
        similar_files = list(similar_files)
        if len(similar_files) == 1 and similar_files[0].endswith(".json"):
            f = os.path.basename(orig_file)
            eppic_store.read_input_file(similar_files[0], f)
            eppic_store.write_output_file(f, orig_file)
            eppic_store.remove_file(similar_files[0])
        elif len(similar_files) == 2:
            for k in similar_files:
                if k.endswith(".json"):
                    pass
