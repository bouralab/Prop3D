import os
from collections import defaultdict

from prodigy.predict_IC import calculate_ic, Prodigy
import Bio.PDB
from Bio.PDB import Select
import pandas as pd
import freesasa
import numpy as np

from molmimic.parsers.superpose import align
from molmimic.parsers.zrank import run_zrank
from molmimic.generate_data.util import read_pdb, replace_chains, extract_chains, rottrans, get_all_chains

from toil.realtimeLogger import RealtimeLogger

def rmsd(p1, p2):
    return np.sqrt(np.square(p1-p2).sum().mean())

def get_cm(coords):
    return np.mean(coords, axis=0)

def get_cm_5(pdb):
    coords = read_pdb(pdb)
    cm = get_cm(coords)
    points = np.tile(cm, (7,1))
    for i in range(3):
        for k, j in enumerate((-1,1)):
            points[2*i+k] += j*5
    return points

def get_coords(residues):
    return np.array([a.get_coord() for r in residues for a in r])

def radius_of_gyration(coords):
    cm = get_cm(coords)
    p2 = np.tile(cm, (coords.shape[0],1))
    return rmsd(coords, p2)

class Complex(Select):
    parser = Bio.PDB.PDBParser(QUIET=1)
    writer = Bio.PDB.PDBIO()

    def __init__(self, pdb, chain1="M", chain2="I", face1=None, face2=None, temp=25.0, method=None, work_dir=None, job=None):
        selection = "{} {}".format(chain1, chain2)
        self.work_dir = work_dir or os.getcwd()
        self.job = job
        self.pdb = pdb
        self.chain1 = chain1
        self.chain2 = chain2
        self.temp = temp
        self.method = method
        self.s = Complex.parser.get_structure("ref", pdb)
        try:
            self.face1 = [r.id[1] for r in self.s[0][chain1] if r.id[1] in face1]
            self.face2 = [r.id[1] for r in self.s[0][chain2] if r.id[1] in face2]
        except KeyError:
            job.log("CHAINS: {} {}".format([c.id for c in self.s.get_chains()], get_all_chains(pdb)))
            raise
        self.prodigy = Prodigy(self.s, selection, temp)
        self.prodigy.predict(distance_cutoff=5.5, acc_threshold=0.05)
        self.interface = set([r.id for rs in calculate_ic(self.s,
            d_cutoff=10.) for r in rs])
        self.neighbors = {c.id:defaultdict(list) for c in self.s.get_chains()}

        for resA, resB in self.prodigy.ic_network:
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

        self.neighbors_id = {}

        for key in self.neighbors[chain1]:
            self.neighbors_id[key.id[1]] = sorted([res.id[1] for res in \
                self.neighbors[chain1][key]])

        self.results = self._get_stats()
        self["method"] = method

    def __getitem__(self, item):
        return self.results[item]

    def __setitem__(self, idx, value):
        self.results[idx] = value

    def compare_to_reference(self, **references):
        for ref_name, ref_complex in list(references.items()):
            rename = lambda l: "{}_{}".format(ref_name, l)
            #self.results = self.results.append(ref_complex.results.rename(rename))
            self.results = self.results.append(ref_complex.compare(self).rename(rename))
        return self.results

    def _get_stats(self):
        results = pd.concat((
            pd.Series(self.prodigy.as_dict()),
            self.BSA(),
            self.radius_of_gyration(),
            self.zrank(),
            self.cns_energy(),
            self.haddock_score()))
        return results

    def set_prefix(self):
        return self.results.rename(lambda l: "{}_{}".format(self.method, l))

    def BSA(self, param_f=None, atom=False):
        cutoffs = {
            "weak transient": (0, 1500),
            "transient": (1500, 2500),
            "permanent": (2500, float("inf"))
        }

        try:
            # Run freesasa
            classifier = freesasa.Classifier(param_f)
            full_structure = freesasa.Structure(self.pdb, classifier)
            chains = freesasa.structureArray(self.pdb, {'separate-chains': True}, classifier)
            chain1 = freesasa.calc(chains[0])
            chain2 = freesasa.calc(chains[1])
            full = freesasa.calc(full_structure)
        except (KeyboardInterrupt, SystemExit):
            #freesasa will throw a regular exception
            raise
        except (AssertionError, IOError, Exception) as e:
            raise Exception(e)

        if not atom:
            c1_asa = chain1.totalArea()
            c2_asa = chain2.totalArea()
            complex_asa = full.totalArea()
            complex_bsa = c1_asa+c2_asa-complex_asa
            for ppi_type, (low_cut, high_cut) in list(cutoffs.items()):
                if low_cut <= complex_bsa < high_cut:
                    break
            else:
                ppi_type = "unknown"
        else:
            sep_atoms = chain[0].nAtoms()
            parse_chain2 = False
            bsa = {}
            for i in range(result.nAtoms()):
                at_id = (structure.chainLabel(i), structure.residueName(i),
                    structure.residueNumber(i), structure.atomName(i))
                complex_asa = result.atomArea(i)
                monomer_asa = chain1.atomArea(i if i<sep_atoms else i-sep_atoms+1)
                at_bsa = monomer_asa-complex_asa
                bsa[at_id] = at_bsa
            ppi_type = None

        resi1 = "+".join(map(str, self.face1))
        sel1 = "face1, chain {} and resi {}".format(self.chain1, resi1)
        sel1 = freesasa.selectArea((sel1,), chains[0], chain1)
        face1_asa = sel1["face1"]

        resi2 = "+".join(map(str, self.face2))
        sel2 = "face2, chain {} and resi {}".format(self.chain2, resi2)
        sel2 = freesasa.selectArea((sel2,), chains[1], chain2)
        face2_asa = sel2["face2"]

        return pd.Series({
            "c1_asa": c1_asa,
            "c2_asa": c2_asa,
            "complex_asa": complex_asa,
            "face1_asa": face1_asa,
            "face2_asa": face2_asa,
            "complex_bsa": complex_bsa,
            "ppi_type": ppi_type,
        })

    def radius_of_gyration(self):
        face1 = get_coords(list(self.neighbors[self.chain1].keys()))
        face2 = get_coords(list(self.neighbors[self.chain2].keys()))
        r1 = radius_of_gyration(face1)
        r2 = radius_of_gyration(face2)
        rI = radius_of_gyration(np.concatenate((face1, face2), axis=0))
        return pd.Series({"Rg_1":r1, "Rg_2":r2, "Rg_interface":rI})

    def zrank(self):
        self.job.log("ZRANK PDB CHAINS: {}".format(get_all_chains(self.pdb)))
        zrank_initial_score = run_zrank(self.pdb, work_dir=self.work_dir, job=self.job)
        zrank_refinement_score = run_zrank(self.pdb, refinement=True,
            work_dir=self.work_dir, job=self.job)
        return pd.Series({
            "zrank_initial_score": zrank_initial_score,
            "zrank_refinement_score": zrank_refinement_score
        })

    def cns_energy(self):
        from molmimic.parsers.CNS import Minimize, calculate_energy
        return calculate_energy(self.pdb, work_dir=self.work_dir, job=self.job)

    def haddock_score(self):
        from molmimic.parsers.haddock import score_complex
        return score_complex(self.pdb, [self.chain1, self.chain2],
            work_dir=self.work_dir, job=self.job)

    def compare(self, moving):
        c1, c2 = self.chain1, self.chain2
        self_chain = self.set_chains("X", "Y")
        irmsd, irmsd_A, irmsd_B, irmsd_avg, irmsd_best = self_chain.iRMSD(moving)
        i_rms, i_tm = self_chain.I_RMS(moving)
        l_rms, l_tm = self_chain.L_RMS(moving)
        mm_rmsd, mm_tm_score = self_chain.MM_TM_score(moving)
        results = {
            "iRMSD": irmsd,
            "iRMSD_A": irmsd_A,
            "iRMSD_B": irmsd_B,
            "iRMSD_avg": irmsd_avg,
            "iRMSD_best": irmsd_best,
            "I_RMS": i_rms,
            "I_TM": i_tm,
            "L_RMS": l_rms,
            "L_TM": l_tm,
            "mm_tm-score": mm_tm_score,
            "mm_rmsd": mm_rmsd,
            "fcc": self_chain.fcc(moving)
        }
        return pd.Series(results)

    def save_interface(self):
        interface = self.pdb+".interface"
        if not os.path.isfile(interface):
            Complex.writer.set_structure(self.s)
            Complex.writer.save(interface, self)
        return interface

    def set_chains(self, chain1, chain2, reset=False):
        if not reset:
            new_pdb = replace_chains(self.pdb, self.pdb+".newchain",
                **{self.chain1:chain1, self.chain2:chain2})
        else:
            new_pdb = self.pdb[:-9]
        return Complex(
            new_pdb,
            chain1=chain1,
            chain2=chain2,
            face1=self.face1,
            face2=self.face2,
            temp=self.temp,
            work_dir=self.work_dir,
            job=self.job
        )

    def accept_residue(self, residue):
        return int(residue.id in self.interface)

    def L_RMS(self, other):
        f, rmsd, tm_score, _ = align(
            self.pdb, self.chain1+self.chain2,
            other.pdb, other.chain1+other.chain2,
            work_dir=self.work_dir,
            job=self.job)
        os.remove(f)
        return rmsd, tm_score

    def I_RMS(self, moving):
        interface1 = self.save_interface()
        interface2 = moving.save_interface()
        f, rmsd, tm_score, _ = align(
            interface1, self.chain1+self.chain2,
            interface2, moving.chain1+moving.chain2,
            work_dir=self.work_dir,
            job=self.job)
        os.remove(f)
        os.remove(interface1)
        os.remove(interface2)
        return rmsd, tm_score

    def MM_TM_score(self, other):
        f, mm_rmsd, mm_tm_score, _ = align(
            self.pdb, self.chain1+self.chain2,
            other.pdb, other.chain1+other.chain2,
            method="mmalign",
            extract="false",
            work_dir=self.work_dir,
            job=self.job)
        os.remove(f)
        return mm_rmsd, mm_tm_score

    def fcc(self, moving):
        "Defined the fraction of native contacts between 2 interfaces"

        #Gets the number of contacts for the reference interface
        total = sum(len(l) for l in list(self.neighbors_id.values()))

        #Assume chain 1 and chain 2 match in both complexes, but they might have different IDs
        chain_map = {self.chain1:moving.chain1, self.chain2:moving.chain2}

        RealtimeLogger(chain_map)
        RealtimeLogger(list(self.neighbors_id.keys()))
        RealtimeLogger(list(moving.neighbors_id.keys()))

        #Finds each common pairs for the 2 interfaces
        common = sum(len(set([r for r in self.neighbors_id[res_id]]).intersection(\
            set([r for r in moving.neighbors_id[res_id]]))) \
            for res_id in self.neighbors_id if res_id in moving.neighbors_id)

        fcc = float(common/total)

        return fcc

    def iRMSD(self, moving):
        #Extract domains
        c1_1f = extract_chains(self.pdb, self.chain1)
        c1_2f = extract_chains(self.pdb, self.chain2)
        c2_1f = extract_chains(moving.pdb, moving.chain1, rename="1")
        c2_2f = extract_chains(moving.pdb, moving.chain2, rename="2")

        import shutil
        shutil.copy(c2_1f, os.path.join("/root", "c2_1.pdb"))

        shutil.copy(c2_2f, os.path.join("/root", "c2_2.pdb"))

        #Superimpose A' to A and B' to B
        best2_1, _, _, matrix2_1 = align(
            c1_1f, self.chain1,
            c2_1f, "1",
            work_dir=self.work_dir,
            job=self.job)
        best2_2, _, _, matrix2_2 = align(
            c1_2f, self.chain2,
            c2_2f, "2",
            work_dir=self.work_dir,
            job=self.job)

        #Rotate A' to A and B' to B using the rotaion matrixes above
        best2_1_2 = rottrans(c2_2f, matrix2_1)
        best2_2_1 = rottrans(c2_1f, matrix2_2)

        #Center of Mass of ref
        cm_ref = np.vstack((get_cm_5(c1_1f), get_cm_5(c1_2f)))

        #Center of Mass of A' -> A
        cm_A = np.vstack((get_cm_5(best2_1), get_cm_5(best2_1_2)))

        #Center of Mass of B' -> B
        cm_B = np.vstack((get_cm_5(best2_2_1), get_cm_5(best2_2)))

        cm_best = np.vstack((cm_A[:7,], cm_B[7:,]))

        rmsdA = rmsd(cm_ref, cm_A)
        rmsdB = rmsd(cm_ref, cm_B)
        irmsd_avg = np.mean((rmsdA, rmsdB))
        irmsd_best = rmsd(cm_ref, cm_best)
        irmsd = min(rmsdA, rmsdB)
        return irmsd, rmsdA, rmsdB, irmsd_avg, irmsd_best

    def match_residues(self, other, r=5.5):
        from sklearn.metrics import jaccard_similarity_score
        face1, face2 = list(zip(*calculate_ic(s, distance_cutoff=r)))
        face1_pred = set([r.id[1] for r in face1])
        face2_pred = set([r.id[1] for r in face2])
        face1_score = jaccard_similarity_score(face1, face1_pred)
        face2_score = jaccard_similarity_score(face2, face2_pred)
        return face1_score, face2_score
