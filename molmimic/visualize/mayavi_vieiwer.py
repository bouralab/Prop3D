import numpy as np
try:
    import mayavi.mlab
except ImportError:
    pass

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

from scipy.ndimage.interpolation import rotate
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

from Bio.SVDSuperimposer import SVDSuperimposer

import seaborn as sns

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

def plot_full_volume(volume, truth, grayscale=True, walls=True):
    xx, yy, zz = np.where(truth[0,...,0] != 1)
    nodes = mayavi.mlab.points3d(xx, yy, zz, mode="cube", scale_factor=1, opacity=0.01)
    nodes.glyph.scale_mode = 'scale_by_vector'
    nodes.mlab_source.dataset.point_data.scalars = volume[0,xx,yy,zz,0].flatten()

    xx, yy, zz = np.where(truth[0,...,0] == 1)
    nodes = mayavi.mlab.points3d(xx, yy, zz, mode="cube", scale_factor=1)
    nodes.glyph.scale_mode = 'scale_by_vector'
    nodes.mlab_source.dataset.point_data.scalars = volume[0,xx,yy,zz,0].flatten()

    nodes = mayavi.mlab.points3d(volumes_negpts[1][0], volumes_negpts[1][1], volumes_negpts[1][2], mode="cube", scale_factor=1, color=(0, 1, 0), opacity=0.01)

def create_figure(n_samples, size=(96,96,96), walls=True, elev=20, azim=135, no_prediction=False):
    # Initialization
    plt.clf()
    plt.cla()
    plt.close('all')


    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    if no_prediction:
        fig = plt.figure(figsize=(n_samples*10,10))
    else:
        fig = plt.figure(figsize=(12,12))
    axes = []
    n_plots = n_samples if no_prediction else 2*n_samples
    for iax in xrange(n_plots):
        ax = fig.add_subplot(1 if no_prediction else 2, n_samples, iax+1, projection='3d')
        ax.set_aspect('equal')

        # Set position of the view
        ax.view_init(elev = elev, azim = azim)

        # Hide tick labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])

        ax.set_xlim(0, size[0])
        ax.set_ylim(0, size[1])
        ax.set_zlim(0, size[2])

        # Clear grid
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        if walls:
            # Tick at every unit
            ax.set_xticks(np.arange(size[0]))
            ax.set_yticks(np.arange(size[1]))
            ax.set_zticks(np.arange(size[2]))

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
        axes.append(ax)

    return fig, axes

def save_fig(fig, fname_prefix):
    #pp = PdfPages()
    plt.savefig('{}.png'.format(fname_prefix), format='png')
    plt.close(fig)
    plt.clf()
    plt.cla()
    plt.close('all')
    del fig

def plot_volume_matplotlib(ax, volume, truth=None, rot_z180=None, rot_x45=None, colors=None):
    """Plots volume in 3D, interpreting the coordinates as voxels
    From: EnzyNet
    """
    if truth is not None:
       truth, volume, rot_z180, rot_x45 = move_to_camera_center(truth, volume)

    plot_matrix(ax, volume, colors=colors)

    if truth is not None:
        return rot_z180, rot_x45


def plot_cube_at(pos = (0,0,0), ax = None, color=(0,1,0), alpha=0.4):
    """Plots a cube element at position pos
    From: EnzyNet
    """
    if ax != None:
        X, Y, Z = cuboid_data(pos)
        ax.plot_surface(X, Y, Z, color=color, rstride=1, cstride=1, alpha=alpha, shade=False)

set2 = sns.color_palette("Set2", 8)
def plot_matrix(ax, matrix, truth=False, colors=False):
    'Plots cubes from a volumic matrix'
    if len(matrix.shape) >= 3:
        if len(matrix.shape) == 4 and matrix.shape[3]==3:
            use_raw_color = True
        else:
            use_raw_color = False

        half_k = matrix.shape[2]/2.
        for i in xrange(matrix.shape[0]):
            for j in xrange(matrix.shape[1]):
                for k in xrange(matrix.shape[2]):
                    #if matrix[i,j,k] == 1:
                    #print "Plotting voxel at", i, j, k
                    if truth:
                        color = (0,0,1)
                    elif use_raw_color:
                        color = matrix[i,j,k]
                    elif colors is not None:
                        color = colors[n]
                        if len(color) > 3:
                            color = set2[np.argmax(color)]
                    else:
                         color = (0,1,0)

                    alpha = abs(half_k-k)/half_k
                    plot_cube_at(pos=(i,j,k), ax=ax, color=color, alpha=alpha)
    elif len(matrix.shape) == 2 and matrix.shape[1] == 3:
        for n, point in enumerate(matrix):
            i,j,k = point
            if truth:
                color = (0,0,1)
            elif colors is not None:
                color = colors[n]
                if len(color) > 3:
                    color = set2[np.argmax(color)]
            else:
                 color = (0,1,0)
            plot_cube_at(pos = (i,j,k), ax = ax, color=color)


def cuboid_data(o, size = (1,1,1)):
    """Gets coordinates of cuboid
    From: EnzyNet
    """
    # Gets the (left, outside, bottom) point
    #o = [(a - b) / 2. for a, b in zip(pos, size)]

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

def to_spherical(xyz, elev_first=False):
    ptsnew = np.zeros(xyz.shape)
    xy = xyz[0]**2 + xyz[1]**2
    ptsnew[0] = np.sqrt(xy + xyz[2]**2)
    if elev_first:
        ptsnew[1] = np.arctan2(np.sqrt(xy), xyz[2]) # for elevation angle defined from Z-axis down
        #ptsnew[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
        ptsnew[2] = np.arctan2(xyz[1], xyz[0])
    else:
        ptsnew[1] = np.arctan2(xyz[1], xyz[0])
        ptsnew[2] = np.arctan2(np.sqrt(xy), xyz[2]) # for elevation angle defined from Z-axis down
        #ptsnew[:,2] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    return ptsnew

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

def get_best_rotation(truth_pts, all_pts):
    mean = np.mean(all_pts)
    r, theta, phi = np.degrees(np.mean(to_spherical(truth_pts-mean), axis=0))

    center_pt = np.degrees(sph2cart(np.pi/4., np.pi/4., r))


    M = rotation_matrix(theta=-phi, phi=-theta, z=1.0)
    new_truth = np.dot(truth_pts, M)
    new_pts = np.dot(all_pts, M)
    return new_truth, all_pts

def rotation_matrix(random = False, theta = 0, phi = 0, z = 0):
    'Creates a rotation matrix'
    # Adapted from: EnzyNet and http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
    # Initialization
    if random == True:
        randnums = np.random.uniform(size=(3,))
        theta, phi, z = randnums
    theta = theta * 2.0*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0  # For magnitude of pole deflection.
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    st = np.sin(theta)
    ct = np.cos(theta)
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    M = (np.outer(V, V) - np.eye(3)).dot(R)

    return M

def move_to_camera_center(truth_pts, all_pts, rot_z180=None, rot_x45=None, return_matrix=False):
    """http://nghiaho.com/?page_id=671"""
    #clusters = max([KMeans(n_clusters=k).fit(truth_pts) for k in xrange(1, 4)], \
    #    key=lambda km: silhouette_score(truth_pts, km.labels_))
    if rot_z180 is None and rot_x45 is None:
        cluster = KMeans(n_clusters=2).fit(truth_pts)

        #Truth
        centroidA = cluster.cluster_centers_[0]
        mean = np.array((48., 48., 48.)) #np.mean(all_pts)

        r, theta, phi = to_spherical(centroidA-mean)

        if theta < 0:
            theta += 2*np.pi
        if theta > np.pi:
            rot_z180 = z_axis_rotation(np.pi)

        if phi < 0:
            phi += 2*np.pi
        elif phi > np.pi/4.:
            rot_x45 = x_axis_rotation(-np.pi/4.)

    changed = True
    if rot_z180 is not None:
        truth_pts = np.dot(truth_pts-mean, rot_z180)+mean
        all_pts = np.dot(all_pts-mean, rot_z180)+mean
        changed = True

    if rot_x45 is not None:
        truth_pts = np.dot(truth_pts-mean, rot_x45)+mean
        all_pts = np.dot(all_pts-mean, rot_x45)+mean
        changed = True

    return truth_pts, all_pts, rot_z180, rot_x45

def rotation_matrix_axis(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def z_axis_rotation(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta),  0],
                     [0,             0,              1]])

def y_axis_rotation(theta):
    return np.array([[np.cos(theta),  0, np.sin(theta)],
                     [0,              1, 0            ],
                     [-np.sin(theta), 0, np.cos(theta)]])

def x_axis_rotation(theta):
    return np.array([[1, 0,             0             ],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])
