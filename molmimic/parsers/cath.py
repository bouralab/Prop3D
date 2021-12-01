import os
import shutil
import yaml
import pandas as pd

from molmimic.generate_data import data_stores
from molmimic.parsers.json import JSONApi, WebService
from molmimic.parsers.container import Container

class CATHApi(JSONApi):
    def __init__(self, cath_store=None, work_dir=None, download=True, max_attempts=2):
        if cath_store is None:
            cath_store = data_stores.cath_api_service
        super().__init__("https://www.cathdb.info/version/v4_3_0/",
            cath_store, work_dir=work_dir, download=download, clean=False,
            max_attempts=max_attempts)

        self.ftp = CATHFTP(cath_store=cath_store, work_dir=work_dir,
            download=download, max_attempts=max_attempts)

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

        elif "from_cath_id_to_depth" in key or "cluster_summary_data" in key:
            with open(file_path) as fh:
                result = yaml.safe_load(fh)
                if "500 Internal Server Error" in result:
                    raise ValueError("{} is invalid".format(file_path))
                return result
        else:
            #Read json
            return super().parse(file_path, key)

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
        if "500 Internal Server Error" in line:
            return True
        if ".pdb" in key:
            return not line.startswith("ATOM")
        else:
            return False


    def extension(self, key):
        if ".pdb" in key or ".png" in key:
            return ""
        elif "stockholm" in key:
            return ".sto"
        elif "from_cath_id_to_depth" in key or "cluster_summary_data" in key:
            return ".yaml"
        else:
            return ".json"

    def fix_superfamily(self, sfam):
        if isinstance(sfam, (tuple, list)):
            return ".".join(map(str, sfam))
        return sfam

    def get_domain_pdb_file(self, cath_domain):
        return self.get("api/rest/id/{}.pdb".format(cath_domain))

    def get_domain_summary(self, cath_domain):
        data = self.get("api/rest/domain_summary/{}".format(cath_domain))
        return pd.Series(data["data"])

    def get_domain_static_image(self, cath_domain, size="L"):
        assert size in "SML", "invalid size"
        return self.get("api/rest/id/{}.png?size={}".format(cath_domain, size))

    def list_cath_domains(self, superfamily):
        superfamily = self.fix_superfamily(superfamily)
        return self.get("api/rest/superfamily/{}".format(superfamily))

    def get_superfamily_sequences(self, superfamily):
        superfamily = self.fix_superfamily(superfamily)
        return self.ftp.get("sequence-data/sequence-by-superfamily/" + \
            "cath-superfamily-seqs-{}.fa".format(superfamily))

    def get_superfamily_clusters(self, superfamily):
        superfamily = self.fix_superfamily(superfamily)
        return self.get("superfamily/{}/cluster_summary_data".format(superfamily))

    #Functional Families (FunFams)
    def list_funfams_in_superfamily(self, superfamily):
        superfamily = self.fix_superfamily(superfamily)
        return self.get("api/rest/superfamily/{}/funfam".format(superfamily))

    def get_funfam_info(self, superfamily, funfam):
        superfamily = self.fix_superfamily(superfamily)
        return self.get("api/rest/superfamily/{}/funfam/{}".format(superfamily, funfam))

    def get_funfam_alignment(self, superfamily, funfam):
        superfamily = self.fix_superfamily(superfamily)
        key = "api/rest/superfamily/{}/funfam/{}/files/stockholm".format(superfamily, funfam)
        return self.get(key, no_api=True)

    #Classification
    def list_children_in_heirarchy(self, cathcode, depth=9):
        cathcode = self.fix_superfamily(cathcode)
        assert 1<=depth<=9, "There are 9 depths in total that correspond to " + \
            "clusters with increasing levels of similarity: C, A, T, H, S, O, L, I, D" + \
            "(see CATH documentation for more info)."
        return self.get("api/rest/cathtree/from_cath_id_to_depth/{}/{}".format(cathcode, depth))

    #Function
    def list_ec_terms(self, superfamily):
        superfamily = self.fix_superfamily(superfamily)
        return self.get("api/rest/superfamily/{}/ec".format(superfamily))

class CATHFTP(WebService):
    def __init__(self, cath_store=None, work_dir=None, download=True, max_attempts=2):
        if cath_store is None:
            cath_store = data_stores.cath_api_service
        super().__init__("ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/",
            cath_store, work_dir=work_dir, download=download, clean=False,
            max_attempts=max_attempts)

    def parse(self, file_path, key):
        if os.path.getsize(file_path) == 0:
            raise ValueError("{} is empty".format(file_path))

        #Make sure raw files aren't removed
        if key not in self.files:
            self.files[key] = (file_path, False)

        #Don't read in image or alingment
        return file_path

    def extension(self, key):
        if ".fa" in key or ".pdb" in key or ".png" in key:
            return ""
        elif "stockholm" in key:
            return ".sto"
        elif "from_cath_id_to_depth" in key or "cluster_summary_data" in key:
            return ".yaml"
        else:
            return ".json"

class FunFamScan(Container):
    IMAGE = "docker://edraizen/cath-tools-genomescan"
    PARAMETERS = [
        "-i", ("input_file", "path:in"),
        "-l", ("hmm_library", "path:in"),
        "-o", ("output_dir", "path:out")]

    def scan(self, input_file, hmm_library):
        results_dir  = os.path.join(self.working_dir, "results_{}".format(
            os.path.splitext(os.path.basname(input_file))
        ))
        results_file = os.path.join(results_dir, "{}.crh".format(
            os.path.basname(input_file).rsplit(".", 1)[0]))

        self(input_file=input_file, hmm_library=hmm_library, output_dir=results_dir)

        return FunFamScan.parse_chr(results_file)

    @staticmethod
    def parse_crh(crh_file):
        results = pd.read_csv(crh_file, delim_whitespace=True, comment="#", header=None,
            names=["query-id", "match-id", "score", "boundaries", "resolved",
                "cond-evalue", "indp-evalue"],
            usecols=["query-id", "match-id"])
        results = results.rename(columns={"query-id":"cath_domain", "match-id": "funfam"})
        results["cath_domain"] = results["cath_domain"].str.split("|").str[-1].str.split("/").str[0]
        results["funfam"] = results["funfam"].str.split("/").str[-1]
        return results
