import os
import pandas as pd
from molmimic.generate_data.iostore import IOStore

def show_predicted_interface(sfam_id, prediction_file):
    structure_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")

    interface_key = "{0}/{0}.observed_interactome".format(int(sfam_id))
    interface_file = "{}.observed_interactome".format(int(sfam_id))
    try:
        interface_store.read_input_file(interface_key, interface_file)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise

    interfaces = pd.read_hdf(interface_file, "table")
    print interfaces[["mol_pdb", "mol_chain", "mol_sdi_id", "mol_domNo"]].head()

    scores = {}

    with open(prediction_file)  as f:
        for line in f:
            fields = line.rstrip().split(",")
            id, pred_resi = fields[0], fields[1:]
            print id, pred_resi
            pdb, chain, sdi, domNo = id.split("_")
            sdi, domNo = sdi[3:], domNo[1:]

            pdb_key = "{}/{}/{}.pdb".format(int(sfam_id), id[1:3].lower(), id)
            pdb_file = os.path.abspath("{}.pdb".format(id))

            interface = interfaces[
                (interfaces["mol_pdb"]==pdb)&(interfaces["mol_chain"]==chain)&\
                (interfaces["mol_sdi_id"]==float(sdi))&(interfaces["mol_domNo"]==int(domNo))\
                ]

            if interface.shape[0] == 0:
                print "Skipped", id
                continue

            true_resi = []
            for resi in interface["mol_res"]:
                true_resi += resi.split(",")
            true_resi = set(true_resi)
            pred_resi = set(pred_resi)

            print "   ", true_resi
            print "   ", pred_resi

            tp_resi = true_resi.intersection(pred_resi)
            fp_resi = pred_resi-true_resi
            fn_resi = true_resi-pred_resi

            f1 = (2*len(tp_resi))/float((2*len(tp_resi)+len(fn_resi)+len(fp_resi)))

            if len(tp_resi) == 0:
                print "    Failed"
                
            else:
                print id, "f1:", f1

            if not os.path.isfile(pdb_file):
                try:
                    structure_store.read_input_file(pdb_key, pdb_file)
                except (SystemExit, KeyboardInterrupt):
                    raise
                except:
                    raise

            pymol_file = """#f1 score: {f1}
load {pdb_file}
hide everything
show surface, {id}
color gray90, {id}

select {id}_tp, resi {tp_resi}
select {id}_fp, resi {fp_resi}
select {id}_fn, resi {fn_resi}
color marine, {id}_tp
color firebrick, {id}_fp
color tv_orange, {id}_fn
""".format(f1=f1, pdb_file=pdb_file, id=id, tp_resi="+".join(tp_resi),
    fp_resi="+".join(fp_resi), fn_resi="+".join(fn_resi))

            with open("{}.pymol".format(id), "w") as pymol:
                pymol.write(pymol_file)

            scores[id] = f1
    return scores
