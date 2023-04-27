import os
from pathlib import Path
from toil.job import Job
from Prop3D.util.iostore import IOStore

class data_stores(object):
    def __init__(self, prefix):
        if isinstance(prefix, Job):
            prefix = prefix._fileStore.jobStore.config.jobStore

        if prefix.startswith("file:") or prefix.startswith("aws:"):
            if prefix.startswith("file:"):
                if "TOIL_S3_HOST" in os.environ:
                    #Save data to S3 bucket even if using file jobStore if you are using a custom S3 server suchas as MinIO
                    prefix = f"aws:us-east-1:{os.environ['USER']}-"
                else:
                    localPath = Path(prefix.split(":")[1]).parent
                    prefix = f"file:{localPath}/"
            else:
                prefix = f"{prefix.rsplit(1)[0]}:{os.environ['USER']}-"
        
        self.prefix = prefix

    @property
    def prepared_cath_structures(self):
        return IOStore.get(f"{self.prefix}prepared-cath-structures")
    
    @property
    def cath_api_service(self):
        return IOStore.get(f"{self.prefix}cath-api-service")

    @property
    def cath_features(self):
        return IOStore.get(f"{self.prefix}cath-features")

    @property
    def data_eppic_cath_features(self): 
        return IOStore.get(f"{self.prefix}data-eppic-cath-features")

    @property
    def eppic_interfaces(self): 
        return IOStore.get(f"{self.prefix}eppic-interfaces")

    @property
    def eppic_store(self): 
        return IOStore.get(f"{self.prefix}Prop3D-eppic-service")

    @property
    def pdbe_store(self): 
        return IOStore.get(f"{self.prefix}Prop3D-pdbe-service")

    @property
    def eppic_local_store(self): 
        return IOStore.get(f"{self.prefix}Prop3D-eppic-local")

    @property
    def raw_pdb_store(self): 
        return IOStore.get(f"{self.prefix}Prop3D-raw-pdb")

    @property
    def eppic_interfaces_sync(self, output_dir=None):
        assert 0
        if output_dir is not None and "S3_ENDPOINT" in os.environ:
            return IOStore.get("file-aws:us-east-1:eppic-interfaces:{}".format(output_dir))
        return eppic_interfaces

    @property
    def uniref(self):
        return IOStore.get(f"{self.prefix}Prop3D-uniref")
