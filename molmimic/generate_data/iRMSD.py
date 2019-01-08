from collections import defaultdict
from prodigy.predict_IC import calculate_ic

from molmimic.parsers.superpose import align
from molmimic.generate_data.util import read_pdb, replace_chains, extract_chains

def rmsd(p1, p2):
    return np.sqrt((p1-p2).square().sum().mean())

def get_cm_5(pdb):
    cm = np.mean(read_pdb(pdb), axis=0)
    points = np.tile(cm, 7)
    for i in xrange(3):
        for k, j in enumerate((-1,1)):
            points[2*i+k] += j*5
    return points

class Complex(Select):
    parser = Bio.PDB.PDBParser(QUIET=1)
    writer = Bio.PDB.PDBIO()

    def __init__(self, pdb, chain1="M", chain2="I", temp=25.0, rename_chains=False):
        selection = "{} {}".format(chain1, chain2)
        self.pdb = pdb
        self.chain1 = chain1
        self.chain2 = chain2
        self.s = Complex.parser.get_structure("ref", pdb)
        self.prodigy = Prodigy(pdb, selection, temp)
        self.prodigy.predict(distance_cutoff=5.5, acc_threshold=0.05)
        self.interface = set([r.id for rs in calculate_ic(s, distance_cutoff=10) for r in rs])
        self.neighbors = {c:defaultdict(list) for c in s.get_chains()}

        for resA, resB in self.interface:
            cA, cB = resA.parent.id, resB.parent.id
            ## Checking for 1st residue (if his chain exist, then if
            ## it is referenced and finally if his partner is already present)
            if resA not in self.neighbors[cA]:
                self.neighbors[cA][resA].append(resB)
            elif resB not in self.neighbors[cA][resA]:
                self.neighbors[cA][resA].append(resB)
            ## Checking for 2nd residue
            if resB not in self.neighbors[cB]:
                self.neighbors[cB][resB].append(resB)
            elif resA not in self.neighbors[cB][resB]:
                self.neighbors[cB][resB].append(resA)

        self.neighbors_id={}

        for key in self.neighbors[chain]:
            self.neighbors_id[key.id]=[res.id[1] for res in self.neighbors[chain][key]]
            self.neighbors_id[key.id].sort()

    @staticmethod
    def compare_pdbs(pdb_ref, pdb_moving):
        complex_ref = Complex(pdb_ref)
        complex_moving = Complex(pdb_moving)
        return complex_ref.compare(complex_moving)

    def get_stats(self, ref=None, ref2=None):
        results = self.prodigy.as_dict()
        if ref is not None:
            comp1 = ref.compare(self)
            if ref2 is not None:
                comp1 = {"obs_"+k:v for k, v in comp1.iteritems()}
            result.update(comp1)
        if ref2 is not None:
            comp2 = {"inf_"+k:v for k, v in ref2.compare(self).iteritems()}
        return pd.Series(results)

    def compare(self, moving, re_rename_chains=True):
        c1, c2 = self.chain1, self.chain2
        self.set_chains("1", "2")
        irmsd_A, irmsd_B, irmsd_avg, irsmd_best = self.iRMSD(moving)
        results = {
            "iRMSD_A": irmsd_A,
            "iRMSD_B": irmsd_B,
            "iRMSD_avg": irmsd_avg,
            "iRMSD_best": irmsd_best,
            "I_RMS": self.I_RMS(moving),
            "L_RMS": self.L_RMS(moving),
            "fcc": self.fcc(moving)
        }
        if re_rename_chains:
            self.set_chains(c1, c2, reset=True)
        return results

    def save_interface(self):
        interface = self.pdb+".interface"
        if not os.path.isfile(interface):
            Complex.writer.set_structure(self.s)
            Complex.writer.save(interface, self)
        return interface

    def set_chains(self, chain1, chain2, reset=False):
        if not reset:
            replace_chains(self.pdb, self.pdb+".newchain",
                **{self.chain1:chain1, self.chain2:chain2})
            self.pdb += ".newchain"
        else:
            self.pdb = self.pdb[:-9]
        self.chain1 = chain1
        self.chain2 = chain2

    def accept_residue(self, residue):
        return int(residue.id in self.interface)

    def L_RMS(self, other):
        f, rmsd = align(self, [self.chain1, self.chain2],
            other, [moving.chain1, moving.chain2])
        os.remove(f)
        return rmsd

    def I_RMS(self, moving):
        interface1 = self.save_interface()
        interface2 = moving.save_interface()
        f, rmsd = align(self.interface, [self.chain1, self.chain2],
            moving.interface, [moving.chain1, moving.chain2])
        os.remove(f)
        os.remove(interface1)
        os.remove(interface2)
        return rmsd

    def fcc(self, moving):
        "Defined the fraction of native contacts between 2 interfaces"

        #Creates neighbors dictionary for each interface residue
        if not self.neighbors:
            self.set_neighbors()
        if not moving.neighbors:
            moving.set_neighbors()

        #Calculation will be done only one one chain thanks to symetry
        chain = next(s.get_chains())

        #Gets the number of contacts for the reference interface
        total = sum(len(l) for l in self.neighbors_id.itervalues())

        #Finds each common pairs for the 2 interfaces
        common = sum(len(set(self.neighbors_id[res_id]).intersection(set(moving.neighbors_id[res_id]))) \
            for res_id in self.neighbors_id if res_id in moving.neighbors_id)

        fcc = float(common/total)

        return fcc

    def iRMSD(self, moving):
        #Extract domains
        c1_1f = extract_chains(self.pdb, self.chain1)
        c1_2f = extract_chains(self.pdb, self.chain2)
        c2_1f = extract_chains(moving.pdb, moving.chain1)
        c2_2f = extract_chains(moving.pdb, moving.chain2)

        #Superimpose A' to A and B' to B
        best1 = align(self.pdb, self.chain1, c2_1f, moving.chain1)
        best2 = align(self.pdb, self.chain2, c2_2f, moving.chain2)

        #Rename chains so they can be aligned and extreacted with complex with same chains
        rename_chains(best1, best1[:-4], **{moving.chain1:moving.chain1*2})
        rename_chains(best2, best2[:-4], **{moving.chain2:moving.chain2*2})
        best1, best2 = best1[:-4], best2[:-4]

        #Re-align complex2 to best oriented domains
        best11, best12 = align(best1, moving.chain1*2, moving.pdb, "".join([moving.chain1, moving.chain2]))
        best21, best22 = align(best2, moving.chain2*2, moving.pdb, "".join([moving.chain1, moving.chain2]))

        #Center of Mass of ref
        cm_ref = np.stack(get_cm_5(c1_1f), get_cm_5(c1_2f))

        #Center of Mass of A' -> A
        cm_A = np.stack(get_cm_5(best11), get_cm_5(best12))

        #Center of Mass of B' -> B
        cm_B = np.stack(get_cm_5(best21), get_cm_5(best22))

        cm_best = np.stack(cm_A[:7,], cm_B[7:,])

        rmsdA = rmsd(cm_ref, cm_A)
        rmsdB = rmsd(cm_ref, cm_B)
        irmsd_avg = np.mean(rmsdA, rmsdB)
        irmsd_best = rmsd(cm_rf, cm_best)
        return rmsdA, rmsdB, irmsd_avg, irmsd_best

    def match_residues(self):
