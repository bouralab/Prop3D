import os
import sys
sys.path.append("/usr/local/Cellar/pymol/2.3.0/libexec/lib/python3.7/site-packages")

from pymol import cmd

def make_lrp_image(pdb, sfam, prefix):
    print("Running", sfam, pdb)

    pdb_file = os.path.join(prefix, "{}-lrp-total_relevance.pdb".format(pdb))
    img_dir = os.path.join(prefix, "lrp-images")
    if not os.path.isdir(img_dir):
        os.makedirs(img_dir)

    cmd.load(pdb_file, pdb)
    cmd.spectrum("b", "magenta_white_cyan", "all", 0, 100)
    cmd.png("test.png")

    cmd.hide("all")
    cmd.show("cartoon")
    cmd.refresh()
    #cmd.bg_color("opaque_background")
    #cmd.set("ray_opaque_background", "off")

    if sfam == "ig":
        #Based on 1TEN orientation in doi:10.1006/jmbi.1994.1582
        cmd.load(os.path.join(prefix, "1TEN.pdb"), "target")
    elif sfam == "ob":
        #Based on 1c4q.A orientation in https://doi.org/10.1002/pro.3742
        cmd.load(os.path.join(prefix, "1c4qA.pdb"), "target")
    elif sfam == "sh3":
        #Based on 1kq2.A orientation in https://doi.org/10.1002/pro.3742
        cmd.load(os.path.join(prefix, "1kq2A.pdb"), "target")

    #cmd.align(pdb, "target")
    cmd.cealign("target", pdb)
    cmd.delete("target")
    cmd.center(pdb)
    cmd.zoom(pdb, buffer=8.0)
    cmd.refresh()

    cmd.ray(2400, 2400)
    cmd.png(os.path.join(prefix, "lrp-images/{}-lrp-total_relevance-cartoon.png".format(pdb)))

    cmd.show("surface")
    cmd.center(pdb)
    cmd.zoom(pdb, buffer=8.0)
    cmd.set("transparency", 0.75)
    cmd.refresh()

    cmd.ray(2400, 2400)
    cmd.png(os.path.join(prefix, "lrp-images/{}-lrp-total_relevance-surface.png".format(pdb)))

    cmd.delete(pdb)


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
