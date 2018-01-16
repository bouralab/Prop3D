import sys
sys.path.append("/data/draizene/flask")
sys.path.append("/data/draizene/molmimic")

from flask import Flask
app = Flask(__name__)

from flask import render_template, request, jsonify

from itertools import izip

import numpy as np

@app.route('/')
def voxel_browser():
    return render_template("voxel_browser.html")

ibis_data = None
def get_ibis_data():
    global ibis_data
    if ibis_data is None:
        from molmimic.torch_model.torch_loader import IBISDataset
        ibis_data = IBISDataset("/data/draizene/molmimic/molmimic/ibis_luca.tab", input_shape=(512,512,512)).data

    return ibis_data

@app.route("/search_results")
def search():
    ibis_data = get_ibis_data()

    fields = request.args.get("search", "").split(".")
    pdb = fields[0]

    if len(fields) > 1:
        result = ibis_data.loc[(ibis_data["pdb"]==pdb)&(ibis_data["chain"]==fields[1])]
    else:
        result = ibis_data.loc[ibis_data["pdb"]==pdb]

    return result.to_json()

@app.route("/browse")
def browse():
    return render_template("ibis_browser.html")

current_protein = None
current_features = None
current_truth = None
@app.route("/get_protein")
def get_protein():
    from molmimic.biopdbtools import Structure, InvalidPDB
    features_name = Structure.get_feature_names()

    use_current = bool(int(request.args.get("use_current", 0)))
    binding_site = request.args.get("binding_site", None)
    input_shape = request.args.get("shape", None)
    if input_shape is None:
        input_shape = (512, 512, 512)
    else:
        if input_shape.count(",")==2:
            input_shape = map(int, input_shape.split(","))
        else:
            return jsonify({"error":"Input shape must 3-comma separted integers"})

    full_protein = bool(int(request.args.get("full_protein", 1)))
    only_atom = True
    only_aa = False #bool(int(request.args.get("only_aa", 1)))
    rotate = bool(int(request.args.get("rotate", 0)))
    orient_pai = bool(int(request.args.get("orient_pai", 1)))
    voxel_size = float(request.args.get("voxel_size", 0.5))
    expand_atom = bool(request.args.get("expand_atom", "true"))

    color_by = request.args.get("color_by", "atom")
    if color_by == "truth" and binding_site is not None:
        color_type = "intensity"
        get_color = None
    if color_by == "atom" or only_atom:
        color_type = "category"
        if only_atom:
            get_color = lambda d: np.argmax(d)
        else:
            get_color = lambda d: np.argmax(d[13:18])
    elif color_by == "residue" or only_aa:
        color_type = "category"
        if only_aa:
            get_color = lambda d: np.argmax(d)
        else:
            get_color = lambda d: np.argmax(d[35:56])
    elif color_by in features_name:
        color_type = "intensity"
        index = features_name.index(color_by)
        get_color = lambda d: d[index]
    else:
        color_by = "residue"
        color_type = "intensity"
        if only_aa or only_atom:
            get_color = lambda d: np.argmax(d)
        else:
            get_color = lambda d: np.argmax(d[35:56])

    if current_protein is None:
        use_current = False

    if not use_current:
        pdb = request.args.get("pdb", None)
        if pdb is None:
            return jsonify({"error":"Must specify pdb"})

        chain = request.args.get("chain", None)
        if chain is None:
            #Get first chain
            chain = get_ibis_data().loc[get_ibis_data()["pdb"]==pdb].iloc[0]["chain"]

        global current_protein
        current_protein = Structure(pdb, chain)
    else:
        pdb = current_protein.pdb
        chain = current_protein.chain

    if orient_pai:
        current_protein.orient_to_pai()

    if rotate:
        current_protein.rotate(1).next()

    if voxel_size != current_protein.voxel_size:
        current_protein.set_voxel_size(voxel_size)

    if binding_site is not None:
        if "," not in binding_site:
            binding_site = get_ibis_data().loc[get_ibis_data()["unique_obs_int"]==binding_site].iloc[0]["resi"]

        binding_site = s.align_seq_to_struc(resi, return_residue=True)

    try:
        global current_features, current_truth
        features = current_protein.get_features(
            input_shape=input_shape,
            residue_list=binding_site,
            only_atom=only_atom,
            only_aa=only_aa,
            expand_atom=expand_atom,
            include_full_protein=full_protein,
            return_truth=binding_site is not None)
    except (KeyboardInterrupt, SystemExit):
        raise
    except InvalidPDB:
        return jsonify({"error":"Invalid pdb"})

    colors = None
    if binding_site is not None:
        (data_idx, current_features), (_, current_truth) = features
        if color_by == "truth":
            colors = current_truth
    else:
        data_idx, current_features = features

    if colors is None:
        colors = [get_color(d) for d in current_features]

    return jsonify({"pdb":pdb, "chain":chain, "indices":data_idx.tolist(), "colors":colors, "color_type":color_type})

def rotate_current_protein():
    if current_protein is None:
        return jsonify({"error":"No current_protein"})

@app.route("/get_atom_features")
def get_atom_features():
    global current_protein, current_features, current_truth

    from molmimic.biopdbtools import Structure
    feature_names = Structure.get_feature_names(only_atom=True)

    atom_serial = int(request.args.get("atom_serial", 0))

    if current_features is not None:
        pdb = current_protein.pdb
        chain = current_protein.chain
        results = [{"name":name, "value":value} for name, value in \
            izip(feature_names, current_features[atom_serial])]
        #results = {"total":len(results), "pdb":pdb, "chain":chain, "atom":atom_serial, "rows":results}
        return jsonify(results)
    else:
        return jsonify({"error":"No current_protein"})

@app.route("/get_atoms_in_voxel")
def get_atoms_in_voxel():
    if current_protein is None:
        return jsonify({"error":"No current_protein"})

    i = request.args.get("i", None)
    j = request.args.get("j", None)
    k = request.args.get("k", None)

    try:
        coord = map(int, (i,j,k))
    except TypeError:
        return jsonify({"error":"Must specify i, j, and k as ints"})


    atoms = [(atom.get_serial_number(), ", ".join(map(str, atom.get_coord().round(decimals=3).tolist()))) for atom in current_protein.get_atoms_from_grids([coord])]
    return jsonify(atoms)
