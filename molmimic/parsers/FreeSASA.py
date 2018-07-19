from tempfile import mkdtemp

try:
    import freesasa
except ImportError:
    freesasa = None

from joblib import Memory

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

from molmimic.util import silence_stdout, silence_stderr

@memory.cache
def run_freesasa(struct, modified=False):
    if freesasa is None:
        print("SASA not installed! SASA will be 0")
        return None, None

    with silence_stdout(), silence_stderr():
        if modified:
            #Note: need to remove hydrogens
            #self.sasa_struct = freesasa.structureFromBioPDB(self.structure)
            sasa_struct = Structure()
            classifier = Classifier()

            atoms = struct.structure.get_atoms()

            for a in atoms:
                if a.element == "H":
                    #Ignore Hydrogens
                    continue

                r = a.get_parent()
                hetflag, resseq, icode = r.get_id()

                c = r.get_parent()
                v = a.get_vector()

                structure.addAtom(a.get_fullname(), r.get_resname(), resseq, c.get_id(),
                                  v[0], v[1], v[2])

            structure.setRadiiWithClassifier(classifier)
        else:
            #Automatically removes hydrogens
            sasa_struct = freesasa.Structure(struct.path)
        sasa = freesasa.calc(sasa_struct)

    return sasa, sasa_struct
