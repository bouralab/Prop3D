# -*- coding: utf-8 -*-
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
import itertools
from random import randint
import pymol

from chempy import cpv



import numpy as np
from collections import defaultdict, Counter

from scipy import spatial

import sys
molpath = os.path.realpath(os.path.join(os.path.dirname(pymol.__script__), "..", ".."))
sys.path.append(molpath)


from molmimic.torch_model.torch_loader import IBISDataset


#############################################################################
#
# drawgridbox.py -- Draw a grid box around a selection
#
# AUTHOR: Cunliang Geng, Troels Schwarz-Linnet
# DATE  : 19/5/2016
#
# Acknowledgement:
# This scirpt was written based on the DrawBoundingBox by Jason Vertrees
#
#############################################################################

vdw_radii = {
    "H" : 1.2,
    "Li" : 1.82,
    "Na" : 2.27,
    "K" : 2.75,
    "C" : 1.7,
    "N" : 1.55,
    "O" : 1.52,
    "F" : 1.47,
    "P" : 1.80,
    "S" : 1.80,
    "Cl" : 1.75,
    "Br" : 1.85,
    "Se" : 1.90,
    "Zn" : 1.39,
    "Cu" : 1.4,
    "Ni" : 1.63,
}

def start_cgo():
    return [
        LINEWIDTH, float(lw),
        BEGIN, LINES,
        COLOR, 0.0, 0.0, 1.0
        ]

class CGO(object):
    def __init__(self, name=None, ):
        self.name = name or "gridbox_" + str(randint(0,10000))
        while self.name in cmd.get_names():
            self.name += "_{}".format(randint(0,10000))
        self.cgo = []

    def add_box(self, i, j, k, voxel=1.0, lw=2.0, r=1.0, g=0.0, b=0.0, style=LINES):
        self.cgo += [
            LINEWIDTH, float(lw),
            BEGIN, style,
            COLOR, r, g, b,
            VERTEX, i, j, k,
            VERTEX, i, j, k+voxel,
            VERTEX, i, j+voxel, k,
            VERTEX, i, j+voxel, k+voxel,
            VERTEX, i+voxel, j, k,
            VERTEX, i+voxel, j, k+voxel,
            VERTEX, i+voxel, j+voxel, k,
            VERTEX, i+voxel, j+voxel, k+voxel,
            VERTEX, i, j, k,
            VERTEX, i+voxel, j, k,
            VERTEX, i, j+voxel, k,
            VERTEX, i+voxel, j+voxel, k,
            VERTEX, i, j+voxel, k+voxel,
            VERTEX, i+voxel, j+voxel, k+voxel,
            VERTEX, i, j, k+voxel,
            VERTEX, i+voxel, j, k+voxel,
            VERTEX, i, j, k,
            VERTEX, i, j+voxel, k,
            VERTEX, i+voxel, j, k,
            VERTEX, i+voxel, j+voxel, k,
            VERTEX, i, j, k+voxel,
            VERTEX, i, j+voxel, k+voxel,
            VERTEX, i+voxel, j, k+voxel,
            VERTEX, i+voxel, j+voxel, k+voxel,
            END
        ]

    def add_cube(self, i, j, k, r=1.0, g=0.0, b=0.0, a=0.4, voxel=1.0):
        # r = voxel/2.
        # r *= 3 ** -.5
        # x = i+r
        # y = j+r
        # z = k+r
        x, y,z = i, j, k

        settings = {'COLOR':[r,g,b], 'INVERT':True}

        # YZ Plane
        self.cgo += make_plane_points(name='p1', l1=[x, y+voxel, z], l2=[x, y, z], l3=[x, y, z+voxel], center=False, makepseudo=False, settings=settings)

        # YZ Plane, shifted in X,
        self.cgo += make_plane_points(name='p6', l1=[x+voxel, y+voxel, z], l2=[x+voxel, y, z], l3=[x+voxel, y, z+voxel], center=False, makepseudo=False, settings=settings)

        # XZ Plane
        self.cgo += make_plane_points(name='p2', l1=[x+voxel, y, z], l2=[x, y, z], l3=[x, y, z+voxel], center=False, makepseudo=False, settings=settings)

        # XZ Plane, shifted in Y, #orange
        self.cgo += make_plane_points(name='p5', l1=[x+voxel, y+voxel, z], l2=[x, y+voxel, z], l3=[x, y+voxel, z+voxel], center=False, makepseudo=False, settings=settings)

        # XY Plane, forest
        self.cgo += make_plane_points(name='p4', l1=[x+voxel, y, z], l2=[x, y, z], l3=[x, y+voxel, z], center=False, makepseudo=False, settings=settings)

        # XY Plane, shifted in Z, red
        self.cgo += make_plane_points(name='p3', l1=[x+voxel, y, z+voxel], l2=[x, y, z+voxel], l3=[x, y+voxel, z+voxel], center=False, makepseudo=False, settings=settings)

    def load(self):
        az = cmd.get('auto_zoom', quiet=1)
        cmd.set('auto_zoom', 0, quiet=1)
        cmd.load_cgo(self.cgo, self.name)
        cmd.set('auto_zoom', az, quiet=1)

def drawBoundingBox(selection="(all)", minX=None, minY=None, minZ=None, maxX=None, maxY=None, maxZ=None, padding=0.0, linewidth=2.0, r=1.0, g=1.0, b=1.0):
    """
    DESCRIPTION
            Given selection, draw the bounding box around it.

    USAGE:
            drawBoundingBox [selection, [padding, [linewidth, [r, [g, b]]]]]

    PARAMETERS:
            selection,              the selection to enboxen.  :-)
                                    defaults to (all)

            padding,                defaults to 0

            linewidth,              width of box lines
                                    defaults to 2.0

            r,                      red color component, valid range is [0.0, 1.0]
                                    defaults to 1.0

            g,                      green color component, valid range is [0.0, 1.0]
                                    defaults to 1.0

            b,                      blue color component, valid range is [0.0, 1.0]
                                    defaults to 1.0

    RETURNS
            string, the name of the CGO box

    NOTES
            * This function creates a randomly named CGO box that minimally spans the protein. The
            user can specify the width of the lines, the padding and also the color.
    """
    if None in [minX, minY, minZ, maxX, maxY, maxZ]:
        ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

    print("Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))

    minX = minX - float(padding)
    minY = minY - float(padding)
    minZ = minZ - float(padding)
    maxX = maxX + float(padding)
    maxY = maxY + float(padding)
    maxZ = maxZ + float(padding)

    if padding != 0:
        print("Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))

    boundingBox = [
            LINEWIDTH, float(linewidth),

            BEGIN, LINES,
            COLOR, float(r), float(g), float(b),

            VERTEX, minX, minY, minZ,       #1
            VERTEX, minX, minY, maxZ,       #2

            VERTEX, minX, maxY, minZ,       #3
            VERTEX, minX, maxY, maxZ,       #4

            VERTEX, maxX, minY, minZ,       #5
            VERTEX, maxX, minY, maxZ,       #6

            VERTEX, maxX, maxY, minZ,       #7
            VERTEX, maxX, maxY, maxZ,       #8


            VERTEX, minX, minY, minZ,       #1
            VERTEX, maxX, minY, minZ,       #5

            VERTEX, minX, maxY, minZ,       #3
            VERTEX, maxX, maxY, minZ,       #7

            VERTEX, minX, maxY, maxZ,       #4
            VERTEX, maxX, maxY, maxZ,       #8

            VERTEX, minX, minY, maxZ,       #2
            VERTEX, maxX, minY, maxZ,       #6


            VERTEX, minX, minY, minZ,       #1
            VERTEX, minX, maxY, minZ,       #3

            VERTEX, maxX, minY, minZ,       #5
            VERTEX, maxX, maxY, minZ,       #7

            VERTEX, minX, minY, maxZ,       #2
            VERTEX, minX, maxY, maxZ,       #4

            VERTEX, maxX, minY, maxZ,       #6
            VERTEX, maxX, maxY, maxZ,       #8

            END
    ]

    boxName = "box_" + str(randint(0,10000))
    while boxName in cmd.get_names():
            boxName = "box_" + str(randint(0,10000))

    cmd.load_cgo(boundingBox,boxName)
    return boxName

def update_selection(old, new):
    if old == "(all)":
        return new
    elif old[0] == "(" and old[-1] == ")":
        return "{} and {})".format(old[:-1], new)
    else:
        return "{} and {}".format(old, new)

def drawgridbox(selection="(all)", volume=256.0, voxel_size=1.0, lw=2.0, r=1.0, g=1.0, b=1.0, show_binding_sites=True, oversample=False, undersample=False):
    """
    DESCRIPTION
        Given selection, draw a grid box around it.

    USAGE:
        drawgridbox [selection, [nx, [ny, [nz, [padding, [lw, [r, [g, b]]]]]]]]

    PARAMETERS:
        selection,    the selection to enboxen
                      defaults to (all)

        nx,           number of grids on axis X
                      defaults to 10

        ny,           number of grids on axis Y
                      defaults to 10

        nz,           number of grids on axis Z
                      defaults to 10

        padding,      defaults to 0

        lw,           line width
                      defaults to 2.0

        r,            red color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        g,            green color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        b,            blue color component, valid range is [0.0, 1.0]
                      defaults to 1.0

    RETURNS
        string, the name of the CGO box

    NOTES
        * This function creates a randomly named CGO grid box. The user can
        specify the number of grids on X/Y/Z axis, the width of the lines,
        the padding and also the color.
    """

    all_coords = Counter()

    extent_min, extent_max = cmd.get_extent(selection)
    extent_x = np.arange(np.floor(extent_min[0])-5, np.ceil(extent_max[0])+5, voxel_size)
    extent_y = np.arange(np.floor(extent_min[1])-5, np.ceil(extent_max[1])+5, voxel_size)
    extent_z = np.arange(np.floor(extent_min[2])-5, np.ceil(extent_max[2])+5, voxel_size)

    com = cmd.centerofmass(selection)
    r = 256/2.
    min_coord = com-np.floor(r)
    max_coord = com+np.ceil(r)
    drawBoundingBox(minX=min_coord[0], minY=min_coord[1], minZ=min_coord[2],
                    maxX=max_coord[0], maxY=max_coord[1], maxZ=max_coord[2])

    model = cmd.get_model(selection)

    coords = np.array([a.coord for a in model.atom])
    mean_coord = np.mean(coords, axis=0)
    volume_center = np.array((extent_x.shape[0]/2., extent_y.shape[0]/2., extent_z.shape[0]/2.))
    shift_by = volume_center-mean_coord
    print(volume_center, shift_by)

    if show_binding_sites:
        volume_center = np.array((128.,128.,128.))
        shift_by = volume_center-mean_coord
        pdb = cmd.get_names("objects", selection=selection)[0].split(".",1)[0].upper()
        for a in model.atom:
            chain = a.chain
            break

        print(pdb, chain)

        ibis_dataset = IBISDataset(os.path.join(molpath, "molmimic", "ibis_luca.tab"), transform=False, input_shape=(256,256,256), balance_classes=bool(undersample))

        print(ibis_dataset.data.iloc[0])

        index = ibis_dataset.data.loc[(ibis_dataset.data['pdb']==pdb)&(ibis_dataset.data["chain"]==chain)].index[0]
        print(index)
        grids = ibis_dataset[0]


        voxels = CGO("{}_VOXELS".format(pdb))

        for grid, truth in zip(grids["indices"], grids["truth"]):
            grid -= shift_by
            g = float(not truth)
            b = float(truth)
            if oversample and not truth:
                continue
            voxels.add_cube(grid[0], grid[1], grid[2], r=0, g=g, b=b)
        voxels.load()

    else:
        xs = np.arange(0, extent_x.shape[0]+1)
        ys = np.arange(0, extent_y.shape[0]+1)
        zs = np.arange(0, extent_z.shape[0]+1)

        mx, my, mz = np.meshgrid(xs, ys, zs)
        voxel_tree = spatial.cKDTree(list(zip(mx.ravel(), my.ravel(), mz.ravel())))

        for a in model.atom:
            vdw = vdw_radii.get(a.name.strip()[0].title(), 2.0)

            neighbors = voxel_tree.query_ball_point(np.array(a.coord)+shift_by, r=vdw)

            for idx in neighbors:
                voxel = voxel_tree.data[idx].astype(int)
                voxel = voxel-shift_by
                all_coords[tuple(voxel.tolist())] += 1

        bonded_boxes = CGO(name="bonded")
        nonbonded_boxes = CGO(name="nonbonded", r=0.0, b=1.0)

        for voxel, count in list(all_coords.items()):
            if count > 1:
                bonded_boxes.add_box(voxel[0], voxel[1], voxel[2], voxel_size)
            else:
                nonbonded_boxes.add_box(voxel[0], voxel[1], voxel[2], voxel_size)

        nonbonded_boxes.load()
        bonded_boxes.load()

    cmd.rotate()

def point(p):
    x, y, z = p
    return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]


def line(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]


def plane(corner1, corner2, corner3, corner4, normal, settings):
    planeObj = []
    # planeObj.extend(point(corner1))
    # planeObj.extend(point(corner2))
    # planeObj.extend(point(corner3))
    # planeObj.extend(point(corner4))
    # planeObj.extend(line(corner1, corner2))
    # planeObj.extend(line(corner2, corner3))
    # planeObj.extend(line(corner3, corner4))
    # planeObj.extend(line(corner4, corner1))

    # Make settings
    if 'ALPHA' in settings:
        planeObj.extend([ALPHA, settings['ALPHA']])

    if 'COLOR' in settings:
        planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
    else:
        planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish

    planeObj.extend([BEGIN, TRIANGLE_STRIP])
    planeObj.append(NORMAL)

    if 'INVERT' in settings:
        if settings['INVERT']==True:
            planeObj.extend(cpv.negate(normal))
        else:
            planeObj.extend(normal)
    else:
        planeObj.extend(normal)


    for corner in [corner1, corner2, corner3, corner4, corner1]:
        planeObj.append(VERTEX)
        planeObj.extend(corner)
    planeObj.append(END)
    return planeObj


def planeFromPoints(p1, p2, p3, vm1=None, vm2=None, center=True, settings={}):
    v1 = cpv.sub(p1, p2)
    v2 = cpv.sub(p3, p2)
    normal = cpv.cross_product(v1, v2)

    if 'translate' in settings:
        vtran = cpv.scale(cpv.normalize(normal), settings['translate'])
        p1_t = cpv.sub(p1, vtran)
        p2_t = cpv.sub(p2, vtran)
        p3_t = cpv.sub(p3, vtran)
        print("New coordinates are:")
        print_info("New", p1_t, p2_t, p3_t)
        print("New coordinates are for normalized plane:")
        v1_t = cpv.normalize(cpv.sub(p1_t, p2_t))
        v2_t = cpv.normalize(cpv.sub(p3_t, p2_t))
        normal_t = cpv.normalize(cpv.cross_product(v1_t, v2_t))
        v2_t = cpv.normalize(cpv.cross_product(normal_t, v1_t))
        p1_t2 = cpv.add(v1_t, p2_t)
        p3_t2 = cpv.add(v2_t, p2_t)
        print_info("Newnormal", p1_t2, p2_t, p3_t2)

    if vm1!=None:
        v1 = cpv.scale(cpv.normalize(v1), vm1)
    if vm2!=None:
        v2 = cpv.scale(cpv.normalize(v2), vm2)

    centrum = p2
    if center:
        corner1 = cpv.add(cpv.add(centrum, v1), v2)
        corner2 = cpv.sub(cpv.add(centrum, v1), v2)
        corner3 = cpv.sub(cpv.sub(centrum, v1), v2)
        corner4 = cpv.add(cpv.sub(centrum, v1), v2)
    else:
        corner1 = cpv.add(cpv.add(centrum, v1), v2)
        corner2 = cpv.add(centrum, v1)
        corner3 = centrum
        corner4 = cpv.add(centrum, v2)

    return plane(corner1, corner2, corner3, corner4, normal, settings)


def print_info(name, coor1, coor2, coor3):
    cs1 = (list(map(float, [ '%.2f' % elem for elem in coor1 ])) )
    cs2 = (list(map(float, [ '%.2f' % elem for elem in coor2 ])) )
    cs3 = (list(map(float, [ '%.2f' % elem for elem in coor3 ])) )
    print("You can also use the function calls with these coordinates")
    print("plane.make_plane_points(name='%s', l1=%s, l2=%s, l3=%s)"%(name, cs1, cs2, cs3))


def make_plane(name,a1='(pk1)',a2='(pk2)',a3='(pk3)', vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, a1, a2, a3

    where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
    """
    # get coordinates from atom selections
    coor1 = cmd.get_model(a1).get_coord_list()[0]
    coor2 = cmd.get_model(a2).get_coord_list()[0]
    coor3 = cmd.get_model(a3).get_coord_list()[0]

    # Help with alternative
    print_info(name, coor1, coor2, coor3)

    # Get the plane
    plane = planeFromPoints(p1=coor1, p2=coor2, p3=coor3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    #makePrimitive(plane, name)
    #cmd.show("cgo", "plane*")
    return plane

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=coor1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=coor2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=coor3)

# Extend function to be called inside pymol
cmd.extend("make_plane", make_plane)

def make_plane_points(name,l1=None,l2=None,l3=None, vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, l1, l2, l3

    where each xys is a list with floats of x,y,z coordinates
    """
    if l1==None or l2==None or l3==None:
        print("Please provide a list of xyz floats for each 3 positions")
        return
    if type(l1) is not list or type(l2) is not list or type(l3) is not list:
        print(type(l1),type(l2),type(l3))
        print("Please provide 3 list of xyz floats for each 3 positions")
        return

    plane = planeFromPoints(p1=l1, p2=l2, p3=l3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    #makePrimitive(plane, name)
    return plane

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=l1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=l2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=l3)


#@jit
def get_centers(xs, ys, zs, voxel_size):
    centers = []
    half_voxel = voxel_size/2.
    _x = xs + half_voxel
    _y = ys + half_voxel
    _z = zs + half_voxel
    for x in _x:
        for y in _y:
            for z in _z:
                centers.append((x,y,z))
    return np.array(centers)

#@jit
def distance(a,b):
    d = (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    d = np.sqrt(d)
    return d

#@jit
def grid_for_atom(coord, vdw_radii, centers):
    """Modified from DeepSite paper
    """
    best_center_index = None
    best_center = None
    best_occupancy = None

    for i, grid_center in enumerate(centers):
        dist_to_center = distance(coord, grid_center)
        x = float(vdw_radii)/dist_to_center
        n = 1-np.exp(-x**12)
        print(n)
        if best_occupancy is None or n>best_occupancy:
            best_occupancy = n
            best_center = grid_center
            best_center_index = i
    return best_center_index, best_center

cmd.extend ("drawgridbox", drawgridbox)
