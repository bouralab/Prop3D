import numpy as np
import mayavi.mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Modified from Fatiando a Terra fatiando.vis.myv (http://www.fatiando.org/v0.2/index.html)
"""

def plot_v(volume):
    xx, yy, zz = np.where(volume == True)
    plot_voxels(xx,yy,zz)

def plot_voxels(xx,yy,zz, colors=None, grayscale=True):
    nodes = mayavi.mlab.points3d(xx, yy, zz, mode="cube", color=(0, 1, 0), scale_factor=1)

    if colors:
        nodes.mlab_source.dataset.point_data.scalars = colors

    mayavi.mlab.show()

def plot_volume(volume, binding_site=None):
    fig = mayavi.mlab.figure(bgcolor=(1,1,1), size=(500, 500))

    tx, ty, tz, _ = np.where(volume == 1)
    mayavi.mlab.points3d(tx, ty, tz, mode="cube", scale_factor=1, color=(0,1,0), opacity=1)

    if binding_site is not None:
        bx, by, bz, _ = np.where(binding_site == 1)
        mayavi.mlab.points3d(bx, by, bz, mode="cube", scale_factor=1, color=(0,0,1), opacity=1)

    wall_west((0,96,0,96,0,96))
    wall_south((0,96,0,96,0,96))
    wall_top((0,96,0,96,0,96))

def wall_north(bounds, color=(0,0,0), opacity=0.1, scale=(1, 1, 1)):
    """
    Draw a 3D wall in Mayavi2 on the North side.

    .. note:: Remember that x->North, y->East and z->Down

    Parameters:

    * bounds : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        The extent of the region where the wall is placed
    * color : tuple = (r, g, b)
        RGB of the color of the wall
    * opacity : float
        Decimal percentage of opacity
    * scale : (slon, slat, sz)
        Scale factors used to exaggerate on a particular direction, e.g., if
        scale = (1, 1, 2), the vertical dimension will be 2x larger than the
        others

    """
    s, n, w, e, t, b = bounds
    _wall([n, n, w, e, b, t], color, opacity, scale)

def wall_south(bounds, color=(0,0,0), opacity=0.1, scale=(1, 1, 1)):
    """
    Draw a 3D wall in Mayavi2 on the South side.

    .. note:: Remember that x->North, y->East and z->Down

    Parameters:

    * bounds : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        The extent of the region where the wall is placed
    * color : tuple = (r, g, b)
        RGB of the color of the wall
    * opacity : float
        Decimal percentage of opacity
    * scale : (slon, slat, sz)
        Scale factors used to exaggerate on a particular direction, e.g., if
        scale = (1, 1, 2), the vertical dimension will be 2x larger than the
        others

    """
    s, n, w, e, t, b = bounds
    _wall([s, s, w, e, b, t], color, opacity, scale)

def wall_east(bounds, color=(0,0,0), opacity=0.1, scale=(1, 1, 1)):
    """
    Draw a 3D wall in Mayavi2 on the East side.

    .. note:: Remember that x->North, y->East and z->Down

    Parameters:

    * bounds : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        The extent of the region where the wall is placed
    * color : tuple = (r, g, b)
        RGB of the color of the wall
    * opacity : float
        Decimal percentage of opacity
    * scale : (slon, slat, sz)
        Scale factors used to exaggerate on a particular direction, e.g., if
        scale = (1, 1, 2), the vertical dimension will be 2x larger than the
        others

    """
    s, n, w, e, t, b = bounds
    _wall([s, n, e, e, b, t], color, opacity, scale)

def wall_west(bounds, color=(0,0,0), opacity=0.1, scale=(1, 1, 1)):
    """
    Draw a 3D wall in Mayavi2 on the West side.

    .. note:: Remember that x->North, y->East and z->Down

    Parameters:

    * bounds : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        The extent of the region where the wall is placed
    * color : tuple = (r, g, b)
        RGB of the color of the wall
    * opacity : float
        Decimal percentage of opacity
    * scale : (slon, slat, sz)
        Scale factors used to exaggerate on a particular direction, e.g., if
        scale = (1, 1, 2), the vertical dimension will be 2x larger than the
        others

    """
    s, n, w, e, t, b = bounds
    _wall([s, n, w, w, b, t], color, opacity, scale)

def wall_top(bounds, color=(0,0,0), opacity=0.1, scale=(1, 1, 1)):
    """
    Draw a 3D wall in Mayavi2 on the Top side.

    .. note:: Remember that x->North, y->East and z->Down

    Parameters:

    * bounds : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        The extent of the region where the wall is placed
    * color : tuple = (r, g, b)
        RGB of the color of the wall
    * opacity : float
        Decimal percentage of opacity
    * scale : (slon, slat, sz)
        Scale factors used to exaggerate on a particular direction, e.g., if
        scale = (1, 1, 2), the vertical dimension will be 2x larger than the
        others

    """
    s, n, w, e, t, b = bounds
    _wall([s, n, w, e, t, t], color, opacity, scale)

def wall_bottom(bounds, color=(0,0,0), opacity=0.1, scale=(1, 1, 1)):
    """
    Draw a 3D wall in Mayavi2 on the Bottom side.

    .. note:: Remember that x->North, y->East and z->Down

    Parameters:

    * bounds : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        The extent of the region where the wall is placed
    * color : tuple = (r, g, b)
        RGB of the color of the wall
    * opacity : float
        Decimal percentage of opacity
    * scale : (slon, slat, sz)
        Scale factors used to exaggerate on a particular direction, e.g., if
        scale = (1, 1, 2), the vertical dimension will be 2x larger than the
        others

    """
    s, n, w, e, t, b = bounds
    _wall([s, n, w, e, b, b], color, opacity, scale)

def _wall(bounds, color, opacity, scale):
    """Generate a 3D wall in Mayavi"""
    p = mayavi.mlab.pipeline.builtin_surface()
    p.source = 'outline'
    p.data_source.bounds = bounds
    p.data_source.generate_faces = 1
    su = mayavi.mlab.pipeline.surface(p)
    su.actor.property.color = color
    su.actor.property.opacity = opacity
    su.actor.actor.scale = scale

def axes(fig, nlabels=5, extent=None, ranges=None, color=(0,0,0),
             width=2, fmt="%-#.2f"):
    """
    Add an Axes module to a Mayavi2 plot or dataset.

    Parameters:

    * plot
        Either the plot (as returned by one of the plotting functions of this
        module) or a TVTK dataset.
    * nlabels : int
        Number of labels on the axes
    * extent : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        Default if the objects extent.
    * ranges : list = [xmin, xmax, ymin, ymax, zmin, zmax]
        What will be display in the axes labels. Default is *extent*
    * color : tuple = (r, g, b)
        RGB of the color of the axes and text
    * width : float
        Line width
    * fmt : str
        Label number format

    Returns:

    * axes : Mayavi axes instace
        The axes object in the pipeline

    """
    a = mayavi.mlab.axes(fig, nb_labels=nlabels, color=color)
    a.label_text_property.color = color
    a.title_text_property.color = color
    if extent is not None:
        a.axes.bounds = extent
    if ranges is not None:
        a.axes.ranges = ranges
        a.axes.use_ranges = True
    a.property.line_width = width
    a.axes.label_format = fmt
    a.axes.x_label, a.axes.y_label, a.axes.z_label = "X", "Y", "Z"
    return a

def plot_full_volume(volume, truth, grayscale=True):
    xx, yy, zz = np.where(truth[0,...,0] != 1)
    nodes = mayavi.mlab.points3d(xx, yy, zz, mode="cube", scale_factor=1, opacity=0.01)
    nodes.glyph.scale_mode = 'scale_by_vector'
    nodes.mlab_source.dataset.point_data.scalars = volume[0,xx,yy,zz,0].flatten()

    xx, yy, zz = np.where(truth[0,...,0] == 1)
    nodes = mayavi.mlab.points3d(xx, yy, zz, mode="cube", scale_factor=1)
    nodes.glyph.scale_mode = 'scale_by_vector'
    nodes.mlab_source.dataset.point_data.scalars = volume[0,xx,yy,zz,0].flatten()

    nodes = mayavi.mlab.points3d(volumes_negpts[1][0], volumes_negpts[1][1], volumes_negpts[1][2], mode="cube", scale_factor=1, color=(0, 1, 0), opacity=0.01)

def plot_volume_matplotlib(volume):
    """Plots volume in 3D, interpreting the coordinates as voxels
    From: EnzyNet
    """
    # Initialization
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection = '3d')
    ax.set_aspect('equal')

    # Parameters
    len_vol = volume.shape[0]

    # Set position of the view
    ax.view_init(elev = 20, azim = 135)

    # Hide tick labels
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    # Plot
    plot_matrix(ax, volume)

    # Tick at every unit
    ax.set_xticks(np.arange(volume.shape[0]))
    ax.set_yticks(np.arange(volume.shape[1]))
    ax.set_zticks(np.arange(volume.shape[2]))

    # Min and max that can be seen
    ax.set_xlim(0, volume.shape[0]-1)
    ax.set_ylim(0, volume.shape[1]-1)
    ax.set_zlim(0, volume.shape[2]-1)

    # Clear grid
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Change thickness of grid
    ax.xaxis._axinfo["grid"]['linewidth'] = 0.1
    ax.yaxis._axinfo["grid"]['linewidth'] = 0.1
    ax.zaxis._axinfo["grid"]['linewidth'] = 0.1

    # Change thickness of ticks
    ax.xaxis._axinfo["tick"]['linewidth'] = 0.1
    ax.yaxis._axinfo["tick"]['linewidth'] = 0.1
    ax.zaxis._axinfo["tick"]['linewidth'] = 0.1

    # Change tick placement
    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.2
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.2
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.2
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.2

    plt.show()

def plot_cube_at(pos = (0,0,0), ax = None, color=1):
    """Plots a cube element at position pos
    From: EnzyNet
    """
    assert 0 <= color <= 1
    if ax != None:
        X, Y, Z = cuboid_data(pos)
        ax.plot_surface(X, Y, Z, color=(0,color,0), rstride=1, cstride=1, alpha=1)

def plot_matrix(ax, matrix):
    'Plots cubes from a volumic matrix'
    for i in xrange(matrix.shape[0]):
        for j in xrange(matrix.shape[1]):
            for k in xrange(matrix.shape[2]):
                if matrix[i,j,k] == 1:
                    #print "Plotting voxel at", i, j, k
                    plot_cube_at(pos = (i-0.5,j-0.5,k-0.5), ax = ax)

def cuboid_data(pos, size = (1,1,1)):
    """Gets coordinates of cuboid
    From: EnzyNet
    """
    # Gets the (left, outside, bottom) point
    o = [(a - b) / 2. for a, b in zip(pos, size)]

    # Get the length, width, and height
    l, w, h = size
    x = np.array([[o[0], o[0] + l, o[0] + l, o[0], o[0]] for i in range(4)])
    y = np.array([[o[1], o[1], o[1] + w, o[1] + w, o[1]],
         [o[1], o[1], o[1] + w, o[1] + w, o[1]],
         [o[1], o[1], o[1], o[1], o[1]],
         [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]])
    z = np.array([[o[2], o[2], o[2], o[2], o[2]],
         [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],
         [o[2], o[2], o[2] + h, o[2] + h, o[2]],
         [o[2], o[2], o[2] + h, o[2] + h, o[2]]])

    return x, y, z
