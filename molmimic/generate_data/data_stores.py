from molmimic.util.iostore import IOStore

prepared_cath_structures = IOStore.get("aws:us-east-1:prepared-cath-structures")
cath_api_service = IOStore.get("aws:us-east-1:cath-api-service")
cath_features = IOStore.get("aws:us-east-1:cath-features")

uniref = IOStore.get("aws:us-east-1:molmimic-uniref")
