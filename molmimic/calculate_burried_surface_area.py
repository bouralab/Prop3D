import os, sys
sys.path.append("/data/draizene/pdb-tools")
sys.path.append("/data/draizene/molmimic")

import subprocess
from itertools import izip
import gzip
import tempfile
import re

import pandas as pd
import numpy as np

from Bio.PDB import PDBList

from molmimic.calculate_features import SwarmJob
from molmimic.calculate_ibis_dataset import parse_ibis_from_pdb
from molmimic.util import get_interfaces_path

cutoffs = {
    "weak transient": (0, 1500),
    "transient": (1500, 2500),
    "permanent": (2500, float("inf"))
}

def get_pdb(pdb, chain1, chain2):
    tmp_pdb_path = "freesasa_{}_{}_{}.pdb".format(pdb, chain1, chain2)
    pdb_path = '/pdb/pdb/{}/pdb{}.ent.gz'.format(pdb[1:3].lower(), pdb.lower())
    if os.path.isfile(pdb_path):
        with open(tmp_pdb_path, 'w') as tmp, gzip.open(pdb_path, 'rb') as pdbf:
            for line in pdbf:
                tmp.write(line)
    else:
        PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=os.getcwd(), file_format="pdb")
        tmp_pdb_path = "pdb{}.ent".format(pdb.lower())
        if os.path.getsize(tmp_pdb_path) == 0:
            raise IOError

    return tmp_pdb_path

def calculate_buried_surface_area(pdb_file, chain1, chain2, residues=None):
    chains = [chain1+chain2, chain1, chain2]
    command = ["freesasa", "--chain-groups={}".format("+".join(chains))]

    if residues is not None:
        residues = residues.replace(",", "+")
        selection = "binding-site, chain {} and resi {}".format(chain1, residues)
        command.append("--select={}".format(selection))

    command.append(pdb_file)

    FNULL = open(os.devnull, 'w')
    freesasa = subprocess.check_output(command, stderr=FNULL)

    freesasa = iter(freesasa.splitlines())
    chain_results = {}
    for line in freesasa:
        try:
            chain_result = process_chain_section(freesasa)
        except StopIteration:
            break
        chain_results[chain_result["chain_names"]] = chain_result

    try:
        bsa = chain_results[frozenset([chain1])]["total"]+chain_results[frozenset([chain2])]["total"]-chain_results[frozenset(chain1+chain2)]["total"]
    except KeyError:
        print pdb, chain1, chain2
        print subprocess.check_output(command)
        raise

    for ppi_type, (low_cut, high_cut) in cutoffs.iteritems():
        if low_cut <= bsa < high_cut:
            return bsa, chain_results, ppi_type
    else:
        return bsa, chain_results, "unknown"

def calculate_surface_area_chain(pdb_file, chain=None, residues=None):
    command = ["freesasa", "--separate-chains"]

    if chain is not None and residues is not None:
        residues = residues.replace(",", "+")
        selection = "binding-site, chain {} and resi {}".format(chain, residues)
        command.append("--select={}".format(selection))

    command.append(pdb_file)

    FNULL = open(os.devnull, 'w')
    freesasa = subprocess.check_output(command, stderr=FNULL)

    freesasa = iter(freesasa.splitlines())
    chain_results = {}
    for line in freesasa:
        try:
            chain_result = process_chain_section(freesasa)
        except StopIteration:
            print chain_results

        chain_results[chain_result["chain_names"]] = chain_result

    if chain is not None:
        return chain_results[frozenset([chain])]
    else:
        return chain_results


def process_chain_section(freesasa_chain):
    line = next(freesasa_chain)
    while not line.startswith("INPUT"):
        line = next(freesasa_chain)

    chain = {}
    chain["source"] = next(freesasa_chain).split(":")[1].strip()
    chain["chain_names"] = frozenset(list(next(freesasa_chain).split(":")[1].strip()))
    chain["model"] = int(next(freesasa_chain).split(":")[1].strip())
    chain["atoms"] = int(next(freesasa_chain).split(":")[1].strip())

    line = next(freesasa_chain)
    while not line.startswith("RESULTS"):
        line = next(freesasa_chain)

    chain["total"] = float(next(freesasa_chain).split(":")[1].strip())
    chain["Apolar"] = float(next(freesasa_chain).split(":")[1].strip())
    chain["polar"] = float(next(freesasa_chain).split(":")[1].strip())
    chain["chains"] = {}
    for c in chain["chain_names"]:
        line = next(freesasa_chain)
        chain["chains"][c] = float(line.split(":")[1].strip())

    process_selections = False
    while True:
        try:
            line = next(freesasa_chain)
        except StopIteration:
            break

        if line.startswith("SELECTIONS"):
            process_selections = True
            chain["selections"] = {}
            continue
        elif line.startswith("#"):
            break

        if process_selections:
            try:
                name, value = line.replace(" ", "").split(":")
            except ValueError:
                continue
            chain["selections"][name] = float(value)


    return chain

coord_re = re.compile('^(ATOM|HETATM)')
def get_binding_site_from_resi(pdb_file, chain, resi):
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            line = line.strip()
            if coord_re.match(line) and line[21] == chain and line[22:27].strip() in resi:
                atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms

def get_partner_binding_site(pdb_file, binding_site1, chain2, dist=5.0):
    resi = []
    with open(pdb_file) as f:
        for line in f:
            line = line.strip()
            if coord_re.match(line) and line[21] == chain2:
                atom = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                for a in binding_site1:
                    if _calculate_atom_distance(atom, a) <= dist:
                        _resi = line[22:27].strip()
                        if len(resi) == 0 or resi[-1] != _resi:
                            resi.append(_resi)
    if len(resi) > 0:
        return ",".join(resi)
    else:
        return None


def _calculate_atom_distance(i, j):
    """Euclidean distance between two 3d points"""
    return np.sqrt((i[0] - j[0])*(i[0] - j[0]) +
                     (i[1] - j[1])*(i[1] - j[1]) +
                     (i[2] - j[2])*(i[2] - j[2]))

def find_partner(dataset_name, pdb_file, pdb, chain1, sdi1, cdd1, resi1, chain2, cdd2):
    interfaces = get_interfaces_path(dataset_name)
    bs1 = get_binding_site_from_resi(pdb_file, chain1, resi1.split(","))
    #print dataset_name, pdb_file, pdb, chain1, sdi1, cdd1, resi1, chain2, cdd2

    possible_matches = {}
    ipath = os.path.join(interfaces, "{}.extended_sdi".format(cdd2))
    with open(ipath) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            #print fields
            if fields[0] == pdb and fields[1] == chain2 and fields[7] == cdd2 and fields[8] == cdd1 and fields[9] in ("{0}{1}_{0}{2}".format(pdb.upper(),chain1,chain2), "{0}{2}_{0}{1}".format(pdb.upper(),chain1,chain2)):
                resi2 = fields[4].split(",")
                bs2 = get_binding_site_from_resi(pdb_file, chain2, resi2)
                for i in bs1:
                    for j in bs2:
                        dist = _calculate_atom_distance(i, j)
                        try:
                            if dist < possible_matches[tuple(fields)]:
                                possible_matches[tuple(fields)] = dist
                        except KeyError:
                            possible_matches[tuple(fields)] = dist
    if len(possible_matches) > 0:
        return min(possible_matches.iteritems(), key=lambda x: x[1])[0]
    else:
        #print "==> Binding site not in IBIS, calculating based on distance",
        resi2 = get_partner_binding_site(pdb_file, bs1, chain2)
        if resi2 is not None:
            print
            return pdb, chain2, None, None, resi2, None, cdd1==cdd2, cdd2, cdd1, "{0}{1}_{0}{2}".format(pdb, chain2, chain1), "1"
        else:
            #print "but failed"
            return None

def find_homolog(dataset_name, query_cdd, query_resi, target_pdb, target_chain):
    from Bio import pairwise2
    from Bio.SubsMat.MatrixInfo import blosum62
    interfaces = get_interfaces_path(dataset_name)

    print query_cdd, query_resi, target_pdb, target_chain

    target = None
    target_score = 0.

    ipath = os.path.join(interfaces, "{}.extended_sdi".format(query_cdd))
    with open(ipath) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            if fields[-1] == "0": continue
            if fields[0] == target_pdb: print fields
            if fields[0] == target_pdb and fields[1] == target_chain:
                #print "posible match", fields
                _target_resi = fields[5]
                try:
                    _target_score = pairwise2.align.globaldx(query_resi, _target_resi, blosum62, score_only=1)
                except StopIteration:
                    continue

                if target is None or _target_score > _target_score:
                    target_score = _target_score
                    target = fields

    return target

def add_bsa_to_cdd(dataset_name, cdd, cdd_file):
    prefix = os.path.splitext(cdd_file)[0]

    open_files = {
        "bsa_file": open("{}.bsa".format(prefix), "w"),
        "weak_transient_bsa_file": open("{}.weak_transient_bsa".format(prefix), "w"),
        "transient_bsa_file": open("{}.transient_bsa".format(prefix), "w"),
        "permanent_bsa_file": open("{}.permmanent_bsa".format(prefix), "w"),
        "obs_weak_transient_bsa_file": open("{}.obs_weak_transient_bsa".format(prefix), "w"),
        "obs_transient_bsa_file": open("{}.obs_transient_bsa".format(prefix), "w"),
        "obs_permanent_bsa_file": open("{}.obs_permmanent_bsa".format(prefix), "w"),
        "inf_weak_transient_bsa_file": open("{}.inf_weak_transient_bsa".format(prefix), "w"),
        "inf_transient_bsa_file": open("{}.inf_transient_bsa".format(prefix), "w"),
        "inf_permanent_bsa_file": open("{}.inf_permmanent_bsa".format(prefix), "w")
    }

    with open(cdd_file) as cdd:
        for f in open_files.values():
            print >> f, "pdb\tchain\tsdi\tdomNum\tresidues\tis_multimer\tcdd\tpartner_cdd\tdb_evidence\tobserved\tsdi2\tdomNum2\tresidues2\tbsa\tppi_type\tbs_vs_bsa\tinf_vs_obs_monomer"

        write_to = ["bsa_file"]
        next(cdd)
        for line in cdd:

            try:
                pdb, chain, sdi1, domNo1, resi1, resn1, is_multimer, cdd, partner_cdd, pdb_evidence, observed = line.rstrip().split("\t")
            except ValueError:
                continue

            if pdb == "3WBD": continue

            print pdb, chain, sdi1, domNo1, "==>",

            if partner_cdd == "No Domain":
                #FIXME How to find domains with no cdd annotation?
                print "Skipped -- No Domain partner"
                continue

            chain1, chain2 = pdb_evidence.split("_")
            evidence_pdb, evidence_chain1, evidence_chain2 = chain1[:-1], chain1[-1], chain2[-1]


            this_pdb_file = get_pdb(pdb, chain, None)

            if chain1 == chain2:
                print "Skippeed -- Same chain not allowed"
                continue

            if not bool(int(observed)):
                #FIXME: get obersrved working first
                #continue
                print pdb_evidence, evidence_chain2, cdd, partner_cdd,
                evidence_pdb_file = get_pdb(evidence_pdb, evidence_chain1, evidence_chain2)
                homolog = find_homolog(dataset_name, cdd, resn1, evidence_pdb, evidence_chain1)

                if homolog is None:
                    print "Skipped -- Cannot find homolog"
                    if sdi1 != "635784":
                        print sdi1
                        assert 0
                    continue

                evidence_sdi1, evidence_domNo1, evidence_resi1 = homolog[2:5]
            else:
                evidence_pdb_file = this_pdb_file
                evidence_sdi1, evidence_domNo1, evidence_resi1 = sdi1, domNo1, resi1

            

            partner = find_partner(dataset_name, evidence_pdb_file, evidence_pdb, evidence_chain1, evidence_sdi1, cdd, evidence_resi1, evidence_chain2, partner_cdd)
            if partner is None:
                print "Skipped -- Cannot find partner"
                continue

            sdi2, domNo2, resi2 = partner[2:5]

            try:
                bsa, all_results, ppi_type = calculate_buried_surface_area(evidence_pdb_file, evidence_chain1, evidence_chain2, residues=evidence_resi1 if observed == "1" else None)
            except IOError:
                print "Skipped -- Cannot open files for BSA calc"
                continue

            bound_source_chain = all_results[frozenset([evidence_chain1, evidence_chain2])]["chains"][evidence_chain1]
            bound_target_chain = all_results[frozenset([evidence_chain1, evidence_chain2])]["chains"][evidence_chain2]
            bsa_monomer = all_results[frozenset([evidence_chain1])]["total"]-bound_source_chain-bound_target_chain

            if ppi_type == "weak transient":
                write_to.append("weak_transient_bsa_file")
            elif ppi_type == "transient":
                write_to.append("transient_bsa_file")
            elif ppi_type == "permanent":
                write_to.append("permanent_bsa_file")

            if observed == "0":
                #Inferred interface
                try:
                    chain_results = calculate_surface_area_chain(this_pdb_file, chain, resi1)
                except IOError:
                    continue

                bs_vs_bsa = chain_results["selections"]["binding-site"]/bsa_monomer
                infMonomer_vs_obsMonomer = chain_results["total"]/all_results[frozenset([evidence_chain1])]["total"]
                if ppi_type == "weak transient":
                    write_to.append("inf_weak_transient_bsa_file")
                elif ppi_type == "transient":
                    write_to.append("inf_transient_bsa_file")
                elif ppi_type == "permanent":
                    write_to.append("inf_permanent_bsa_file")
            else:
                bs_vs_bsa = all_results[frozenset([evidence_chain1])]["selections"]["binding-site"]/bsa_monomer
                infMonomer_vs_obsMonomer = 1.
                if ppi_type == "weak transient":
                    write_to.append("obs_weak_transient_bsa_file")
                elif ppi_type == "transient":
                    write_to.append("obs_transient_bsa_file")
                elif ppi_type == "permanent":
                    write_to.append("obs_permanent_bsa_file")

            for f in write_to:
                print >> open_files[f], "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line.rstrip(), sdi2, domNo2, resi2, bsa, ppi_type, bs_vs_bsa, infMonomer_vs_obsMonomer)

            print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line.rstrip(), sdi2, domNo2, resi2, bsa, ppi_type, bs_vs_bsa, infMonomer_vs_obsMonomer)

            os.remove(this_pdb_file)
            try:
                os.remove(evidence_pdb_file)
            except OSError:
                pass

    for f in open_files.values():
        f.close()

def submit_ibis_cdd(dataset_name, job_name="calc_bsa", dependency=None):
    ibis_data = get_interfaces_path(dataset_name)
    CDD = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "MMDB", "StructDomSfam.csv"), usecols=["label"]).drop_duplicates().dropna()
    CDD = sorted(CDD["label"].apply(lambda cdd: cdd.replace("/", "").replace("'", "\'")).tolist())

    job = SwarmJob(job_name+"_full", walltime="18:00:00")
    for cdd in CDD:
        cdd_f = os.path.join(ibis_data, "{}.extended_sdi".format(cdd.replace("/", "")))
        job += "/data/draizene/3dcnn-torch-py2 python {} run {} \"{}\" {}\n".format(__file__, dataset_name, cdd.replace("/", ""), cdd_f)

    jid = job.run(dependency=dependency)
    print jid

if __name__ == "__main__":
    if len(sys.argv) == 2:
        submit_ibis_cdd(sys.argv[1])
    elif len(sys.argv) == 5 and sys.argv[1] == "run":
        add_bsa_to_cdd(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print len(sys.argv), sys.argv
