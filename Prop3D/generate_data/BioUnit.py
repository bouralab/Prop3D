import os

import numpy as np
import Bio.PDB
from prodigy.predict_IC import calculate_ic

from Prop3D.util.iostore import IOStore
from Prop3D.util.pdb import rottrans_from_matrix, tidy

from toil.realtimeLogger import RealtimeLogger

parser = Bio.PDB.PDBParser(QUIET=1, PERMISSIVE=True)

def check_contacts(original_complex, mol_file, mol_res, int_file, int_res, return_vars=False):
    try:
        s = parser.get_structure("ref", original_complex)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        return False

    chains = [c.id for c in s.get_chains()]
    assert len(chains)==2

    try:
        ic = calculate_ic(s, d_cutoff=8.)
    except ValueError:
        return False

    if chains[0]<chains[1]:
        face1, face2 = zip(*ic)
    else:
        face2, face1 = zip(*ic)
    face1, face2 = set([r.id[1] for r in face1]), set([r.id[1] for r in face2])
    RealtimeLogger.info("FACE1 {} {}".format(mol_res, face1))
    RealtimeLogger.info("FACE2 {} {}".format(int_res, face2))

    full_face1 = len(face1.intersection(mol_res))/float(len(mol_res))
    full_face2 = len(face2.intersection(int_res))/float(len(int_res))
    if not return_vars:
        return len(ic) > 0 and full_face1 >= 0.5 and full_face2 >= 0.5
    else:
        return ic, face1, face2, full_face1, full_face2

def build_biounit(moving_mol, mol_pdb, mol_chain, moving_int, int_pdb, int_chain, work_dir):
    #Original files may be correct
    yield moving_int

    #Get BioUnit info for interactint domain (since inferred are Superimposed into it)
    pdbId = int_pdb
    in_store = IOStore.get("aws:us-east-1:Prop3D-pdb")
    pdb_file_base = os.path.join(pdbId[1:3].lower(), "pdb{}.ent.gz".format(pdbId.lower()))
    pdb_file = os.path.join(work_dir, "pdb{}.ent.gz".format(pdbId.lower()))
    in_store.read_input_file(pdb_file_base, pdb_file)

    remarks = pdbremarks(pdb_file)
    if 350 not in remarks:
        raise RuntimeError("No REMARK 350")
    quat = quat350(remarks[350])

    # for i, (Mi, ti) in enumerate(quat):
    #     if not (Mi == np.eye(3)).all() and not (ti == np.zeros(3)).all():
    #         rt_mol = "{}.rottrans.{}.pdb".format(os.path.splitext(moving_mol)[0], i)
    #         tidy(rottrans_from_matrix(moving_mol, Mi, ti, rt_mol), replace=True)
    #     else:
    #         rt_mol = moving_mol

    for j, (Mj, tj) in enumerate(quat):
        #if not (Mj == np.eye(3)).all(): # and not (tj == np.zeros(3)).all():
        rt_int = "{}.rottrans.{}.pdb".format(os.path.splitext(moving_int)[0], j)
        tidy(rottrans_from_matrix(moving_int, Mj, tj, rt_int), replace=True)
        # else:
        #     continue
        #     rt_int = moving_int

        yield rt_int

    del remarks
    del quat

def build_sym_transforms(moving_mol, mol_pdb, mol_chain, moving_int, int_pdb, int_chain,
  sym_op=None, return_matrices=False, work_dir=None):

    if work_dir is None:
        work_dir = os.getcwd()

    #Get BioUnit info for interactint domain (since inferred are Superimposed into it)
    pdbId = int_pdb

    if sym_op == "X,Y,Z":
        return moving_int
    elif sym_op is not None:
        _sym_op = [dim.split("+") for dim in sym_op.split(",")]
        _sym_op_rot, sym_op_trans = zip(*_sym_op)
        sym_op_trans = np.array(sym_op_trans, dtype=float)

        sym_op_rot = np.eye(3)
        for i, (r,d) in enumerate(zip(_sym_op_rot, ("X","Y","Z"))):
            sym_op_rot[i,i] = float(r.replace(d, "1"))

    pdb_file, format = s3_download_pdb(pdbId)

    assert format == "pdb"

    remarks = pdbremarks(pdb_file)
    if 290 not in remarks:
        raise RuntimeError("No REMARK 350")
    sym_ops = get_sym_ops(remarks[350])

    transformed_pdbs = []
    for i, transform in enumerate(quat):
        if sym_op is not None and sym_op_rot==transform["M_au"]:
            for j in range(3):
                for pm in (1, -1):
                    t = np.array([t+pm if j==jj else t for r in transform["M_au"]])
                    if sym_op_trans == t:
                        if return_matrices:
                            return transform
                        else:
                            rt_int = "{}.rottrans.{}.pdb".format(os.path.splitext(moving_int)[0], sym_op)
                            return tidy(rottrans_from_matrix(moving_int,
                                transform["M"], transform["t"], rt_int), replace=True)

        else:
            if return_matrices:
                return transform
            else:
                rt_int = "{}.rottrans.{}.pdb".format(os.path.splitext(moving_int)[0], transform["op"])
                transformed_pdbs.append(tidy(rottrans_from_matrix(moving_int,
                    transform["M"], transform["t"], rt_int), replace=True))

    return transformed_pdbs

def pdbremarks(filename):
    '''
    Read REMARK lines from PDB file. Return dictionary with remarkNum as key
    and list of lines as value.
    '''
    remarks = dict()
    if filename[-3:] == '.gz':
        import gzip
        f = gzip.open(filename)
    else:
        f = open(filename)
    for line in f:
        recname = line[0:6]
        if recname == 'REMARK':
            num = int(line[7:10])
            lstring = line[11:]
            remarks.setdefault(num, []).append(lstring)
    f.close()
    return remarks

def quat350(rem350):
    '''
    Get transformation matrices for biomolecule 1 from REMARK 350.
    '''
    biomt = []
    chains = []
    seenbiomolecule = False

    it = iter(rem350)

    for line in it:
        if line.startswith('BIOMOLECULE:'):
            if seenbiomolecule:
                break
            seenbiomolecule = True
        elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
            chains =[chain.strip() for chain in line[30:].split(',')]
        elif line.startswith('                   AND CHAINS:'):
            chains += [chain.strip() for chain in line[30:].split(',')]
        elif line.startswith('  BIOMT'):
            current_M = np.eye(3)
            current_t = np.zeros(3)

            for i in range(3):
                l = next(it) if i > 0 else line
                RealtimeLogger("LINE IS {}".format(l))
                row = int(l[7])
                num = int(l[8:12])
                vec = l[12:].split()
                vec = map(float, vec)
                current_M[i, :] = vec[:-1]
                current_t[i] = vec[-1]
            biomt.append((current_M.T, current_t))

    return biomt

def get_sym_ops(rem290):
    symops = {}
    seen_symop = False

    it = iter(rem290)

    for line in it:
        if line.startswith('CRYSTALLOGRAPHIC SYMMETRY:'):
            if seenbiomolecule:
                break
            seenbiomolecule = True
            for _ in range(4):
                next(it)

        elif seen_symop:
            if line.strip() == "":
                seen_symop = False

            sym_num, sym_op = line.strip().split()
            symops[sym_num[0]] = {"op": sym_op}

            sym_op = [dim.split("+") for dim in sym_op.split(",")]
            rot, trans = zip(*sym_op)
            trans = np.array(trans, dtype=float)

            rot_mat = np.eye(3)
            for i, (r,d) in enumerate(zip(rot, ("X","Y","Z"))):
                rot_mat[i,i] = float(r.replace(d, "1"))

            symops[sym_num[0]]["M_au"] = rot_mat
            symops[sym_num[0]]["t_au"] = trans


        elif line.strip().startswith("SMTRY"):
            current_M = np.eye(3)
            current_t = np.zeros(3)

            for i in range(3):
                l = next(it)
                row = int(l[7])
                sym_num = int(l[8:12])
                vec = l[12:].split()
                vec = map(float, vec)
                current_M[i, :] = vec[:-1]
                current_t[i] = vec[-1]

            symops[sym_num]["M"] = current_M
            symops[sym_num]["t"] = current_t

    return symops
