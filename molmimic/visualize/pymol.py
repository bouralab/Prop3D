# -*- coding: utf-8 -*-
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
import itertools
from random import randint
from numba import jit

import sys
print sys.executable

import numpy as np
from collections import defaultdict, Counter

from scipy import spatial

#############################################################################
#
# drawgridbox.py -- Draw a grid box around a selection
#
# AUTHOR: Cunliang Geng
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
    def __init__(self, name=None, lw=2.0, r=1.0, g=0.0, b=0.0, style=LINES):
        self.name = name or "gridbox_" + str(randint(0,10000))
        while self.name in cmd.get_names():
            self.name += "_{}".format(randint(0,10000))
        self.cgo = [
            LINEWIDTH, float(lw),
            BEGIN, style,
            COLOR, r, g, b
        ]

    def add_box(self, i, j, k, voxel=1.0):
        self.cgo += [
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
        ]

    def load(self):
        cmd.load_cgo(self.cgo+[END], self.name)

def update_selection(old, new):
    if old == "(all)":
        return new 
    elif old[0] == "(" and old[-1] == ")":
        return "{} and {})".format(old[:-1], new)
    else:
        return "{} and {}".format(old, new)

def drawgridbox(selection="(all)", volume=256.0, voxel_size=1.0, lw=2.0, r=1.0, g=1.0, b=1.0):
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

    xs = np.arange(0, extent_x.shape[0]+1)
    ys = np.arange(0, extent_y.shape[0]+1)
    zs = np.arange(0, extent_z.shape[0]+1)

    mx, my, mz = np.meshgrid(xs, ys, zs)
    voxel_tree = spatial.cKDTree(zip(mx.ravel(), my.ravel(), mz.ravel()))
    
    model = cmd.get_model(selection)

    coords = np.array([a.coord for a in model.atom])
    mean_coord = np.mean(coords, axis=0)
    volume_center = np.array((extent_x.shape[0]/2., extent_y.shape[0]/2., extent_z.shape[0]/2.))
    shift_by = volume_center-mean_coord
    print volume_center, shift_by

    for a in model.atom:
        vdw = vdw_radii.get(a.name.strip()[0].title(), 2.0)

        neighbors = voxel_tree.query_ball_point(np.array(a.coord)+shift_by, r=vdw)

        for idx in neighbors:
            voxel = voxel_tree.data[idx].astype(int)
            voxel = voxel-shift_by
            all_coords[tuple(voxel.tolist())] += 1

    bonded_boxes = CGO(name="bonded")
    nonbonded_boxes = CGO(name="nonbonded", r=0.0, b=1.0)

    for voxel, count in all_coords.iteritems():
        if count > 1:
            bonded_boxes.add_box(voxel[0], voxel[1], voxel[2], voxel_size)
        else:
            nonbonded_boxes.add_box(voxel[0], voxel[1], voxel[2], voxel_size)
            
    nonbonded_boxes.load()
    bonded_boxes.load()

    return bonded_boxes.name

@jit
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

@jit
def distance(a,b):
    d = (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    d = np.sqrt(d)
    return d

@jit
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
        print n
        if best_occupancy is None or n>best_occupancy:
            best_occupancy = n
            best_center = grid_center
            best_center_index = i
    return best_center_index, best_center

cmd.extend ("drawgridbox", drawgridbox)