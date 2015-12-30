"""
grid_utils
----------

Utilities for gridding and ungridding data efficiently in Numpy.

"""

import numpy as np
from scipy.interpolate import griddata

def meshgrid_pairs(x, y):
    """ Create meshed grid, but as pairs of coordinates

    Parameters
    ----------
    x: np.array
        1D array representing coordinates of a grid (i.e. axis values)
    y: np.array
        1D array representing coordinates of a grid (i.e. axis values)

    Returns
    -------
    gridded_pairs: np.array
        Returns a 3D array of gridded points, with shape (len(x), len(y), 2)
        where each index pair (i, j) returns the corresponding coordinates
        pair.

    Notes
    -----
    'x' is the inner dimension, i.e. image dimensions are (n_y, n_x). This is
    counterintuitive (to me at least) but in line with numpy definitions.
    """
    xg, yg = np.meshgrid(x, y)
    gridded_pairs = np.dstack((xg, yg))
    return gridded_pairs


def ungrid(x, y, d):
    """ Unravel a 2D mesh D with 1D coord axes X and Y into a list of length X*Y

    Use this to get back to a coordinate list [(x1, y1, z1), (x2, y2, z2), ...].

    Parameters
    ----------
    x: np.array
        1D array representing coordinates of a grid (i.e. axis values)
    y: np.array
        1D array representing coordinates of a grid (i.e. axis values)
    d: np.array
        2D gridded data, where x and y are the axes of the data.

    Returns
    -------
    xyz: np.array
        Coordinate list of (x, y, z) data. Shape is (len(x)*len(y), 3),
        where each entry is an (x, y, z) coordinate.

    Notes
    -----
    'x' is the inner dimension, i.e. image dimensions are (n_y, n_x). This is
    counterintuitive (to me at least) but in line with numpy definitions.
    """

    gp = meshgrid_pairs(x, y)
    gp = np.dstack((gp, d))
    lx, ly, lz = gp.shape[0], gp.shape[1], gp.shape[2]
    xyz = gp.reshape((lx*ly, lz))
    #xyz = xyz[:, (1, 0, 2)]                             # Need YXZ -> XYZ?
    return xyz


def grid_xyz(xyz, n_x, n_y, **kwargs):
    """ Grid data as a list of X,Y,Z coords into a 2D array

    Parameters
    ----------
    xyz: np.array
        Numpy array of X,Y,Z values, with shape (n_points, 3)
    n_x: int
        Number of points in x direction (fastest varying!)
    n_y: int
        Number of points in y direction

    Returns
    -------
    gridded_data: np.array
        2D array of gridded data, with shape (n_x, n_y)

    Notes
    -----
    'x' is the inner dimension, i.e. image dimensions are (n_y, n_x). This is
    counterintuitive (to me at least) but in line with numpy definitions.
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    x_ax = np.linspace(np.min(x), np.max(x), n_x)
    y_ax = np.linspace(np.min(y), np.max(y), n_y)

    xg, yg = np.meshgrid(x_ax, y_ax)

    data = griddata(xyz[:, :2], z, (xg, yg), **kwargs)
    return data    


def scatter(idx, vals, target):
    """ Efficient scatter operation using Numpy
    
    target[idx] += vals, but allowing for repeats in idx
    
    Parameters
    ----------
    idx: np.array
        Array of index values into which to scatter
    vals: np.array
        Array of values which shall be scattered into
    target: np.array
        Target array to scatter values into
    """
    np.add.at(target, idx, vals)
    

def gather(idx, vals):
    """ Efficient gather operation using Numpy
    
    Parameters
    ----------
    idx: np.array
        Array of index values from which to gather
    vals: np.array
        Array of values that are being gathered
    """
    return np.take(vals, idx, axis=0)