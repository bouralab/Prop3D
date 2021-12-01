import os
import sys
import random
from itertools import combinations, permutations, tee

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB import Selection
from Bio.PDB.Polypeptide import three_to_one

from molmimic.common.Structure import Structure
from molmimic.parsers.MODELLER import MODELLER
from molmimic.generate_data.create_input_files import create_input_files

"""MLP: Multiple Loop Permutations

Modified by Eli Draizen

Dai L, Zhou Y. Characterizing the existing and potential structural space of
proteins by large-scale multiple loop permutations. J Mol Biol. 2011 May 6;
408(3):585-95. doi: 10.1016/j.jmb.2011.02.056. Epub 2011 Mar 2.
PMID: 21376059; PMCID: PMC3075335.
"""

DATA_DIR = os.environ.get("DATA_DIR", "/home/bournelab/data-eppic-cath-features")

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def rearrange_ss(cath_domain, cathcode, cutoff=15):
    pdb_file = os.path.join(DATA_DIR, "prepared-cath-structures", *cathcode.split("."), f"{cath_domain}.pdb")
    features_path = os.path.join(DATA_DIR, "cath_features", *cathcode.split("."))

    structure = Structure(pdb_file, cath_domain, features_path=features_path)

    #Step 1: Get Secondary Structure
    ss_groups, loop_for_ss, original_order = get_secondary_structures(structure)

    #Step 2: Find allowable combinations
    allowed_ss_combos, prohibited_ss_combos = get_allowable_loop_combinations(ss_groups, cutoff=cutoff)

    #Step 3: Permute strands
    input_file = None
    for perm_pdb_name, loop_data, sequence in get_possible_loop_permutations(cath_domain, ss_groups, allowed_ss_combos, prohibited_ss_combos, loop_for_ss, original_order):
        #Step 4: Model loops correctly with MODELLER
        files = os.listdir(os.getcwd())
        for f in files:
            if os.path.splitext(os.path.basename(perm_pdb_name))[0] in f and ".BL" in f and f.endswith(".pdb"):
                permuted_pdb = f
                break
        else:
            permuted_pdb = model_loops(perm_pdb_name, loop_data, sequence)
        input_file = create_input_files(permuted_pdb, input_file=input_file)

    return input_file

def get_secondary_structures(structure):
    """1. Secondary structure for each domain was assigned by the program DSSP.
    Short helical and strand segments (<4 residues) were treated as coils to
    decrease the number of loops for a given protein by reducing the number of
    secondary structure segments (SSSs).
    """
    ss_type = structure.atom_features[["is_helix", "is_sheet", "Unk_SS"]]

    ss_type = ss_type.rename(columns={"is_helix":"H", "is_sheet":"E", "Unk_SS":"X"})
    ss_type = ss_type.idxmax(axis=1)

    ss_groups = ss_type.groupby([(ss_type != ss_type.shift()).cumsum()-1])

    #Merge group shorter than 4 residues
    for i, ss_group in ss_groups:
        if 0<i<ss_groups.ngroups-1:
            this_group = ss_group.iloc[0]
            prev_group = ss_groups.get_group(i-1).iloc[0]
            next_group = ss_groups.get_group(i+1).iloc[0]
            if len(ss_group)<25 and prev_group == next_group:
                ss_type.loc[ss_group.index] = prev_group
            elif len(ss_group)<50 and next_group == "X" and this_group != prev_group:
                ss_type.loc[ss_group.index] = "X"

    #Regroup with correct SS
    ss_groups = ss_type.groupby([(ss_type != ss_type.shift()).cumsum()-1])

    ss_atom_groups = []
    loop_for_ss = {}
    original_order = {}
    for i, ss_group in ss_groups:
        #Get all atoms from SS and loops
        ss_atoms = tuple(structure.get_atoms(include_atoms=ss_group.index))

        if ss_group.iloc[0] != "X":
            ss_atom_groups.append(ss_atoms)
            original_order[ss_atoms] = len(ss_atom_groups)
        elif len(ss_atom_groups)>0 and ss_group.iloc[0] == "X":
            loop_for_ss[ss_atom_groups[-1]] = ss_atoms

    first_group = ss_groups.get_group(0)
    if first_group.iloc[0] == "X":
        loop_for_ss[1] = tuple(structure.get_atoms(include_atoms=first_group.index))

    last_group = ss_groups.get_group(ss_groups.ngroups-1)
    if last_group.iloc[0] == "X":
        loop_for_ss[ss_groups.ngroups] = tuple(structure.get_atoms(include_atoms=last_group.index))

    return ss_atom_groups, loop_for_ss, original_order

def get_allowable_loop_combinations(ss_groups, cutoff=15):
    """2. The distances between the N-terminus of one SSS and the C-terminus of
    another SSS were calculated for all SSS pairs. The N and C termini of
    two SSSs were allowed to connect by building a new loop between them if
    their distance is less than a cutoff distance (15 angstroms initially).
    The connection between two N (or C) termini of two SSSs was not allowed
    in order to maintain the original N to C direction. The original loops
    longer than 15 angstroms were unchanged."""

    allowed_ss_combos = []
    prohibited_ss_combos = []

    for ss1, ss2 in permutations(ss_groups, 2):
        atom1 = ss1[0].get_parent()["CA"]
        atom2 = ss2[-1].get_parent()["CA"]
        dist = atom1-atom2

        if dist < cutoff:
            allowed_ss_combos.append((ss1, ss2))
        else:
            prohibited_ss_combos.append((ss1, ss2))

    return allowed_ss_combos, prohibited_ss_combos

writer = PDBIO()
def write_atom(fp, atom, atom_number, resseq):
    hetfield, _, _ = atom.get_parent().get_id()
    resname = atom.get_parent().get_resname()
    s = writer._get_atom_line(
                    atom,
                    hetfield,
                    " ",
                    atom_number,
                    resname,
                    resseq,
                    " ",
                    atom.get_parent().get_parent().get_id(),
                )
    fp.write(s)

class Alanine(Residue):
    def __init__(self):
        Residue.__init__(self, (" ", 1, " "), "ALA", 0)
        for atom in (" N  ", " H  ", " CA ", " HA ", " CB ", " HB1", " HB2", " HB3"):
            self.add(Atom(atom.strip(), np.random.rand(3), 20., 1., " ", atom, 1))

def get_possible_loop_permutations(pdb_name, ss_groups, allowed_ss_combos, prohibited_ss_combos, loop_for_ss, original_order, short_loops=True, random_loops=True):
    """3. A combinatorial search was made for all possible loop permutations
    allowed. If two SSSs change from sequence neighbor to non-neighbor after
    rearrangement, their connection loop will be removed. Meanwhile, new loops
    will be built to connect two SSSs that become sequence neighbors after
    rearrangement. For example, a protein with 6 SSSs is arranged in a native
    structure as 1-2-3-4-5-6. One possible rearrangement of this sequence is
    6-5-2-3-4-1. This rearrangement requires retaining two native loops for
    unchanged neighboring SSSs between 2-3 and 3-4, removing three native loops
    (1-2, 4-5 and 5-6) because they are no longer sequence neighbors (5-6 is not
    same as 6-5 because of the N to C direction), and building three new loops
    between 6-5, 5-2, and 4-1. In this study, we limited ourselves to generate
    100 MLP structures and a maximum of five permutated loops per proteins. If
    the number of permutations is greater than 100, we decreased the cutoff
    distance with a step size of 0.5Å to reduce the number of loops allowed to
    permute until the number of permutations is less than or equal to 100."""

    for perm in permutations(ss_groups, len(ss_groups)):
        allowed_ss_combos
        #[((ss1, atom1), (ss2, atom2)), ...]
        order = list(perm) #[ss for p in perm for ss in p]

        """

        ((1,2),(2,3),(3,4),(4,5),(5,6),
        """

        perm_name = [original_order[ss] for ss in order]
        if perm_name == list(range(1,len(ss_groups)+1)):
            #Same as orignal PDB
            print("Skipped -- same as original")
            continue

        skip = False
        for ss1, ss2 in pairwise(order):
            if (ss1, ss2) not in allowed_ss_combos:
                skip = True
                print(f"Skipped {perm_name} -- {original_order[ss1]} and {original_order[ss2]} cannot be joined")
                break

        if skip:
            continue

        perm_name_str = "-".join(map(str, perm_name))
        print(perm_name_str)
        loop_data = []
        sequence = ""
        atom_number = 1
        resseq = 1
        prev_res = None
        perm_pdb_name = f"{pdb_name}_{perm_name_str}_no_loops.pdb"
        with open(perm_pdb_name, "w") as perm_pdb:
            if perm_name[0] in [1, len(ss_groups)]:
                for atom in loop_for_ss[perm_name[0]]:
                    o_res = atom.get_parent()
                    if prev_res is not None and o_res != prev_res:
                        resseq += 1

                    write_atom(perm_pdb, atom, atom_number, resseq)

                    atom_number += 1
                    prev_res = o_res

            for i, ss in enumerate(order):
                for atom in ss:
                    o_res = atom.get_parent()
                    if prev_res is not None and o_res != prev_res:
                        resseq += 1

                    write_atom(perm_pdb, atom, atom_number, resseq)

                    atom_number += 1
                    prev_res = o_res

                loop = loop_for_ss[ss]
                if loop is not None:
                    loop_aa = Selection.unfold_entities(loop, "R")
                    if short_loops and i<len(order)-1:
                        #Create shorter loops
                        next_ss = order[i+1]
                        next_start = next_ss[0].get_parent()["CA"]
                        this_end = ss[-1].get_parent()["CA"]
                        dist = this_end-next_start
                        n_res = min(int(np.ceil(dist/2.5)), 6)
                        if random_loops:
                            loop_aa = random.choices(loop_aa, k=n_res) #, replace=False)
                        else:
                            loop_aa = [Alanine() for i in range(n_res)]

                #else:
                    loop_start = (resseq, loop_aa[0].get_parent().get_id())
                    for o_res in loop_aa:
                        if prev_res is not None and o_res != prev_res:
                            resseq += 1

                        for atom in o_res:
                            write_atom(perm_pdb, atom, atom_number, resseq)
                            atom_number += 1

                        prev_res = o_res

                    loop_end = (resseq, o_res.get_parent().get_id())
                    loop_data.append((loop_start, loop_end))

        yield perm_pdb_name, loop_data, sequence

modeller = None
def model_loops(perm_pdb_name, loop_data, sequence):
    """4. All new loops were built by the program Modloop47. We estimated the
    number of residues for a new loop by dividing the end-to-end distance with
    2.5Å. This approximate formula was obtained from a statistical analysis of
    the end-to-end distances of short loops. Because the maximum end-to-end
    distance for a loop to be permutated is 15Å, the maximum number of residues
    for a rebuilt loop is 6. That is, we have avoided building potentially
    unrealistic long loops (>6) 47. All loops were built with alanine residues
    for computational efficiency."""
    global modeller
    if modeller is None:
        modeller = MODELLER()
    permuted_pdb = modeller.model_loops(perm_pdb_name, loop_data, sequence)
    return permuted_pdb

if __name__ == "__main__":
    assert len(sys.argv)==3, f"Error. Must run '{sys.argv[0]} [cath_domain] [cathcode]'. E.g. '{sys.argv[0]} 4unuA00 2.60.40.10'"

    cath_domain = sys.argv[1]
    cathcode = sys.argv[2]

    rearrange_ss(cath_domain, cathcode)
