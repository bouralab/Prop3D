from molmimic.util.iostore import IOStore

prepared_cath_structures = IOStore.get("aws:us-east-1:prepared-cath-structures")
cath_api_service = IOStore.get("aws:us-east-1:cath-api-service")
cath_features = IOStore.get("aws:us-east-1:cath-features")
data_eppic_cath_features = IOStore.get("aws:us-east-1:data-eppic-cath-features")

eppic_interfaces = IOStore.get("aws:us-east-1:eppic-interfaces")
eppic_store = IOStore.get("aws:us-east-1:molmimic-eppic-service")
pdbe_store = IOStore.get("aws:us-east-1:molmimic-pdbe-service")
eppic_local_store = IOStore.get("aws:us-east-1:molmimic-eppic-local")
raw_pdb_store = IOStore.get("aws:us-east-1:molmimic-raw-pdb")

def eppic_interfaces_sync(output_dir=None):
    if output_dir is not None:
        return IOStore.get("file-aws:us-east-1:eppic-interfaces:{}".format(output_dir))
    return eppic_interfaces

uniref = IOStore.get("aws:us-east-1:molmimic-uniref")
