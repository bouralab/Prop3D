import os
import shutil
import yaml

from molmimic.parsers.json import JSONApi

class CATHApi(JSONApi):
    def __init__(self, cath_store, work_dir=None, download=True, max_attempts=2):
        super().__init__("https://www.cathdb.info/version/v4_2_0/",
            cath_store, work_dir=work_dir, download=download, clean=False,
            max_attempts=max_attempts)

    def parse(self, file_path, key):
        if ".pdb" in key:
            #Make sure file can open and is PDB file
            with open(file_path) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        break
                else:
                    raise ValueError("{} not a PDB file".format(file_path))

            #Make sure raw files aren't removed
            if key not in self.files:
                self.files[key] = (file_path, False)

            return file_path
        elif ".png" in key or "stockholm" in key:
            if os.path.getsize(file_path) == 0:
                raise ValueError("{} is empty".format(file_path))

            #Make sure raw files aren't removed
            if key not in self.files:
                self.files[key] = (file_path, False)

            #Don't read in image or alingment
            return file_path

        elif "from_cath_id_to_depth" in key or "cluster_summary_data" in data:
            with open(file_path) as fh:
                return yaml.safe_load(fh)
        else:
            #Read json
            return super(self, CATHApi).parse(file_path, key)

    def parse_pdb(self, domain_file):
        try:
            with open(domain_file) as f, open(domain_file+".atom", "w") as f2:
                for line in f:
                    if line.startswith("ATOM"):
                        print(line.rstrip(), file=f2)
        except IOError:
            #Make sure download restarts
            raise ValueError

        try:
            os.remove(fname)
        except (OSError, FileNotFoundError):
            pass

        shutil.move(domain_file+".atom", domain_file)

        return domain_file

    def check_line(self, key, line, attempts):
        if ".pdb" in key:
            return not line.startswith("ATOM")
        else:
            return False


    def extension(self, key):
        if ".pdb" in key or ".png" in key:
            return ""
        elif "stockholm" in key:
            return ".sto"
        elif "from_cath_id_to_depth" in key or "cluster_summary_data" in data:
            return ".yaml"
        else:
            return ".json"

    def get(self, key, no_api=False):
        if no_api:
            return self.get(key)
        else:
            return self.get("api/rest/"+key)

    def get_domain_pdb_file(self, cath_domain):
        return self.get("id/{}.pdb".format(cath_domain))

    def get_domain_summary(self, cath_domain):
        return self.get("id/{}".format(cath_domain))

    def get_domain_static_image(self, cath_domain, size="L"):
        assert size in "SML", "invalid size"
        return self.get("id/{}.png?size={}".format(cath_domain, size))

    def list_cath_domains(self, superfamily):
        return self.get("superfamily/{}".format(superfamily))

    def get_superfamily_clusters(self, superfamily):
        return self.get("superfamily/{}/cluster_summary_data".format(superfamily), no_api=True)

    #Functional Families (FunFams)
    def list_funfams_in_superfamily(self, superfamily):
        return self.get("superfamily/{}/funfam".format(superfamily))

    def get_funfam_info(self, superfamily, funfam):
        return self.get("superfamily/{}/funfam/{}".format(superfamily, funfam))

    def get_funfam_alignment(self, superfamily, funfam):
        return self.get("superfamily/{}/funfam/{}/files/stockholm".format(superfamily, funfam))

    #Classification
    def list_children_in_heirarchy(self, cathcode, depth=9):
        assert 1<=depth<=9, "There are 9 depths in total that correspond to " + \
            "clusters with increasing levels of similarity: C, A, T, H, S, O, L, I, D" + \
            "(see CATH documentation for more info)."
        return self.get("cathtree/from_cath_id_to_depth/{}/{}".format(cathcode, depth))

    #Function
    def list_ec_terms(self, superfamily):
        return self.get("superfamily/{}/ec".format(superfamily))
