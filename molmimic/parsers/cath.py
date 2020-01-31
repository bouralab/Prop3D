import shutil
import yaml

from molmimic.parsers.json import JSONApi

class CATHApi(JSONApi):
    def __init__(self, cath_store, work_dir=None, download=True, clean=True, max_attempts=2):
        super(self, CATHApi).__init__("https://www.cathdb.info/version/v4_2_0/api/rest/",
            cath_store, work_dir=work_dir, download=download, clean=clean,
            max_attempts=max_attempts)

    def parse(self, file_path, key):
        if ".pdb" in key:
            return self.parse_pdb(file_path)
        elif ".png" in key or "stockholm" in key:
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
            with open(domain_file) as f, open(domain_file+".atom") as f2:
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
        if ".pdb" in key:
            return return ".pdb"
        elif ".png" in key:
            return ".png"
        elif "stockholm" in key:
            return ".sto"
        elif "from_cath_id_to_depth" in key or "cluster_summary_data" in data:
            return ".yaml"
        else:
            return ".json"

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
        return self.get("superfamily/{}/cluster_summary_data".format(superfamily))

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
