import os

from molmimic.util import safe_remove
from molmimic.generate_data import data_stores
from molmimic.parsers.json import WebService
from molmimic.parsers.container import Container

class Uniref(WebService):
    def __init__(self, uniprot_store, work_dir=None, download=True,
      max_attempts=2, job=None):
        self.job = job
        super(WebService, self).__init__("ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/",
            uniprot_store, work_dir=work_dir, download=download, clean=False,
            max_attempts=max_attempts)

    def download(self, key, fname):
        if key.startswith("blast/"):
            #Fake error to ensure it fails and can start running blast
            raise KeyError(key)

        return super().download(key, fname)

    def parse(self, file_path, key):
        return file_path

    def get_uniref100_db(self):
        make_blast_db = MakeBlastDB(work_dir=self.work_dir, job=self.job)
        name = "uiref100"
        try:
            db_files = [self.get(name+"."+ext) for ext in make_blast_db.blast_db_ext]
        except KeyError:
            #Need to run
            uniref100_file_gz = self.get("uniref100/uniref100.fasta.gz")
            db_name, db_files =  make_blast_db(in_file=uniref100_file_gz,
                out_file="uniref100", title="uniref100")
            safe_remove(uniref100_file_gz)
            for db_file in db_files:
                store.write_output_file(db_file, "blast/{}".format(os.path.basename(v))

        return db_files

class BLAST(Container):
    IMAGE = "docker://edraizen/blast"

class MakeBlastDB(BLAST):
    PARAMETERS = ["-in", ("in_file", "path:in"),
                  (":dbtype:prot", "str", "-dbtype {}"),
                  "-parse_seqids",
                  "-out", ("out_file", "path:out")
                  (":title", "str", "-title {}")]
    ENTRYPOINT = "makeblastdb"
    RETURN_FILES = True

    blast_db_ext = ["phr", "pin", "pog", "psd", "psi", "psq"]

    def __call__(self, *args, **kwds):
        name = super().__call__(*args, **kwds)
        files = [os.path.isfile(out+"."+ext) for ext in self.blast_db_ext]
        if not all(files):
            raise RuntimeError("Unable to make blast db, not all files are present: {}".format(files))
        return name, files
