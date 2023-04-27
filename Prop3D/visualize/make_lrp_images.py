import os
import sys
sys.path.append("/usr/local/Cellar/pymol/2.3.0/libexec/lib/python3.7/site-packages")

from pymol import cmd

import numpy as np
from sklearn.preprocessing import MinMaxScaler
from collections import defaultdict

def make_lrp_image(pdb_file, hue="blue_white_red", map_to_ca=False):

    pdb = os.path.basename(pdb_file)
    cmd.load(pdb_file, pdb)
    cmd.color("gray70", pdb)



    #cmd.copy("spheres", pdb)

    myspace = {'bfactors': []}
    cmd.iterate(pdb, 'bfactors.append(b)', space=myspace)

    low = np.quantile(myspace["bfactors"], 0.5)
    #mid = np.quantile(myspace["bfactors"], 0.75)
    high = np.quantile(myspace["bfactors"], 0.8)

    print("low", low)
    print("high", high)

    if map_to_ca:
        bfact = {'bfactors': defaultdict(dict)}
        cmd.iterate("structure", 'bfactors[resi][name] = b', space=bfact)

        for resi, bfactors in bfact["bfactors"].items():
            best_b = max(bfactors.values())
            cmd.alter(f"structure and resi {resi}", f"b={best_b}")

    cmd.spectrum("b", hue, pdb, low, high) #magenta_white_cyan, 0, 100

    cmd.hide("all", pdb)
    cmd.show("cartoon", pdb)
    cmd.refresh()

    # cmd.alter("spheres", "vdw=1")

    # scaler = MinMaxScaler()
    # sphere_radii = scaler.fit_transform(np.array(myspace["bfactors"]).reshape(-1,1)).flatten()
    # for i, (sphere_radius, b) in enumerate(zip(sphere_radii, myspace["bfactors"])):
    #     if b<mid:
    #         cmd.set("sphere_scale", 0, f"ID {i}")
    #         cmd.set("sphere_transparency", 1.0, f"ID {i}")
    #         cmd.set_bond("stick_transparency", 0.8, f"ID {i}")
    #     else:
    #         cmd.set("sphere_scale", sphere_radius/2, f"ID {i}")
    #         cmd.set("sphere_transparency", 1-sphere_radius, f"ID {i}")
    #         cmd.set_bond("stick_transparency", 1-sphere_radius, f"ID {i}")

    
    # cmd.show("spheres", "spheres")
    # cmd.show("sticks", "spheres")
    # cmd.set_bond("stick_radius", 0.1, "spheres")
    #cmd.set("stick_transparency", 0.50, "spheres")

    cmd.refresh()
    #cmd.bg_color("opaque_background")
    #cmd.set("ray_opaque_background", "off")

    # if sfam is not None:
    #     if sfam == "ig":
    #         #Based on 1TEN orientation in doi:10.1006/jmbi.1994.1582
    #         cmd.load(os.path.join(prefix, "1TEN.pdb"), "target")
    #     elif sfam == "ob":
    #         #Based on 1c4q.A orientation in https://doi.org/10.1002/pro.3742
    #         cmd.load(os.path.join(prefix, "1c4qA.pdb"), "target")
    #     elif sfam == "sh3":
    #         #Based on 1kq2.A orientation in https://doi.org/10.1002/pro.3742
    #         cmd.load(os.path.join(prefix, "1kq2A.pdb"), "target")
    #
    #     #cmd.align(pdb, "target")
    #     cmd.cealign("target", pdb)
    #     cmd.delete("target")



    cmd.center(pdb)
    cmd.zoom(pdb, buffer=8.0)
    cmd.refresh()

    cmd.ray(2400, 2400)
    cmd.png(os.path.join(prefix, "lrp-images/{}-lrp-total_relevance-cartoon.png".format(pdb)))

    # cmd.show("surface")
    # cmd.center(pdb)
    # cmd.zoom(pdb, buffer=8.0)
    # cmd.set("transparency", 0.75)
    # cmd.refresh()

    # cmd.ray(2400, 2400)
    # cmd.png(os.path.join(prefix, "lrp-images/{}-lrp-total_relevance-surface.png".format(pdb)))

    cmd.delete(pdb)

    cmd.save()


def make_lrp_images(scores_csv):
    with open(scores_csv) as f:
        header = next(f)
        for line in f:
            fields = line.rstrip().split(",")
            sfam, pdb = fields[1].split("/")
            pdb = pdb[:-4]
            make_lrp_image(pdb, sfam, os.path.join(os.path.dirname(scores_csv), "lrp"))


    os.chdir(ocwd)

if __name__ == "__main__":
    make_lrp_images(sys.argv[1])
