import os
from molmimic.util import safe_remove
from molmimic.parsers.container import Container

PDB_SECTIONS = [
    #Title
    "HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND", "SOURCE", "KEYWDS",
    "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARKS", "REMARK",

    #Primary Structure
    "DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES", "MODRES",

    #Heterogen
    "HET", "FORMUL", "HETNAM", "HETSYN",

    #Secondary Structure
    "HELIX", "SHEET",

    #Connectivity Annotation
    "SSBOND", "LINK", "CISPEP",

    #Misc Features
    "SITE",

    #Crystallographic and Coordinate Transformation
    "CRYST1", "MTRIXn", "ORIGXn", "SCALEn",

    #Coordinate Section
    "MODEL", "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL",

    #Connectivity
    "CONECT",

    #Bookkeeping
    "MASTER", "END"
]

class Reduce(Container):
    IMAGE = 'docker://edraizen/reduce:latest'
    LOCAL = ["reduce"]
    ENTRYPOINT = "/opt/reduce/reduce"
    PARAMETERS = [
        (":trim", "store_true", ["-Trim"]),
        (":his", "store_true", ["-HIS"]),
        ("in_file", "path:in", ["{}"]),
        ]

    def run(self, in_file, out_file=None, **kwds):
        """Run reduce and save output to file"""

        output = self(in_file, **kwds)

        if out_file is None:
            out_file = os.path.join(self.work_dir, os.path.basename(in_file)+".reduced.pdb")

        with open(out_file, "w") as f:
            for line in output.splitlines():
                if line[:6].strip() in PDB_SECTIONS:
                    print(line.rstrip(), file=f)

        return out_file

    def protonate(self, in_file, out_file=None,):
        return self.run(in_file=in_file, out_file=out_file, his=True)

    def deprotonate(self, in_file, out_file=None):
        return self.run(in_file=in_file, out_file=out_file, trim=True)

    def reprotonate(self, in_file, out_file=None):
        deprotonated_file = self.deprotonate(in_file)
        protonated_file = self.protonate(deprotonated_file, out_file=out_file)
        safe_remove(deprotonated_file)
        return protonated_file
