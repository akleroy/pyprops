import numpy as np
import math
from scipy.ndimage import label

# Utilities for working with shapes in a data cube.

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DEFINE STRUCTURING ELEMENTS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def ellipse_connectivity(
    major=None,
    minor=None,
    posang=None,
    ):
    """
    ...
    """

    # ------------------------------------------------------------
    # Error Checking on Inputs
    # ------------------------------------------------------------

    if major==None:
        print "Requires a major axis."
        return

    if minor > major:
        print "Minor axis must be <= major axis."
        return

    if minor==None:
        minor = minor

    if posang==None:
        posang=0.0

    # ------------------------------------------------------------
    # Build the ellipse
    # ------------------------------------------------------------
    
    npix = 2*np.ceil(major/2.0)+1
    y,x = np.indices((npix, npix))
    y -= np.mean(y)*1.
    x -= np.mean(x)*1.
    dtor = math.pi/180.
    xp = x*np.cos(posang*dtor) - y*np.sin(posang*dtor)
    yp = x*np.sin(posang*dtor) + y*np.cos(posang*dtor)
    
    return (((xp/(major/2.0))**2 + (yp/(minor/2.0))**2) <= 1.0)


def rectangle_connectivity(
    major=None,
    minor=None,
    posang=None,
    ):
    """
    ...
    """

    # ------------------------------------------------------------
    # Error Checking on Inputs
    # ------------------------------------------------------------

    if major==None:
        print "Requires a major axis."
        return

    if minor > major:
        print "Minor axis must be <= major axis."
        return

    if minor==None:
        minor = minor

    if posang==None:
        posang=0.0

    # ------------------------------------------------------------
    # Build the rectangle
    # ------------------------------------------------------------
    
    npix = 2*np.ceil(major/2.0)+1
    y,x = np.indices((npix, npix))
    y -= np.mean(y)*1.
    x -= np.mean(x)*1.
    dtor = math.pi/180.
    xp = x*np.cos(posang*dtor) - y*np.sin(posang*dtor)
    yp = x*np.sin(posang*dtor) + y*np.cos(posang*dtor)
    
    return ((np.abs(xp) <= (major/2.0))*(np.abs(yp) <= (minor/2.0)))


def simple_connectivity(
    ndim=3,
    skip_axes=None, 
    corners=False):
    """
    Return simple connectivity either with or without corners for n
    dimensions (default 3). Can suppress connectivity along one or
    more axes.
    """

    # ------------------------------------------------------------
    # Error Checking on Inputs
    # ------------------------------------------------------------

    if (type(ndim) != type(0)):
        print "Requires an integer number of axes."
        return
    
    # ... boolean data type
    if (type(corners) != type(True)):
        print "Requires boolean data type for corners flag."
        return
    
    # ------------------------------------------------------------
    # Generate the connectivity
    # ------------------------------------------------------------

    if corners == True:
        # ... connect along diagonals
        connect = np.ones(np.ones(ndim)*3)
    else:
        # ... suppress diagonals
        connect = np.zeros(np.ones(ndim)*3)
        center = np.ones(ndim,dtype=np.dtype('int'))
        connect[tuple(center)] = 1
        for axis in np.arange(ndim):
            pixel = np.ones(ndim)
            pixel[axis] = 0
            connect[tuple(pixel)] = 1
            pixel[axis] = 2
            connect[tuple(pixel)] = 1

    # ------------------------------------------------------------
    # Suppress connect along specified axes
    # ------------------------------------------------------------

    if skip_axes != None:
    
        # ... catch the case of an integer
        if (type(skip_axes) == type(1)):
            skip_axes = [skip_axes]

        ind = np.indices(connect.shape)
        for axis in skip_axes:
            blank = (ind[axis] == 0) + (ind[axis] == 2)
            connect[blank] = 0                        

    return connect

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# LABEL CONTIGUOUS REGIONS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#
# ... where "contiguous" is defined by the structuring elements.
#

def blob_color(mask,
               corners=False,
               connect=None):
    """
    Basic blob-coloring wrapper. Feed it a binary mask and it will
    return a colored array. Also accepts a custom connectivity
    structure. If this is not used, it generates a default
    connectivity of +/- one step in each dimension direction. The flag
    corners tells it whether to also connect along diagonals.
    """

    # ------------------------------------------------------------
    # Error Checking on Inputs
    # ------------------------------------------------------------
    
    # ... existence
    try:
        mask
    except NameError:
        print "Requires data."
        return
    
    # ... numpy array
    if (type(mask) != type(np.arange(0))):
        print "Requires a numpy array."
        return

    # ... boolean data type
    if (mask.dtype != np.dtype('bool')):
        print "Requires boolean data type."
        return
    
    # ------------------------------------------------------------
    # Define connectivitiy
    # ------------------------------------------------------------
    #
    # Can be supplied by the user, else calculate it here

    if connect == None:
        connect = connectivity(mask.ndim, corners=corners)

    # Wrap the SciPy implementation
    color, ncolors = label(mask,structure=connect)

    return color

# ------------------------------------------------------------
# BLOB EXTRACTION
# ------------------------------------------------------------

# Needs work.

# These focus on remapping blobs to vectors or dictionaries for easy
# manipulation. Definitely an area that could be either sped up or
# merged into moment calculation. Overlap with find_objects but not
# totally clear that those make it so much faster.

def blob_to_vec(
    data,
    mask=None,
    origin=None
    ):
    """
    Return coordinates and data value given a mask. Needs a rewrite
    using find_object fast slicer.
    """

    # ... a default mask
    if mask == None:
        mask = np.isfinite(data)

    # ... note dimensions
    ndim = data.ndim

    # ... one-d case
    if ndim == 1:
        x = (mask).nonzero()
        if origin != None:
            x += origin
        val = data[x]
        return x, val

    # ... two-d case
    if ndim == 2:
        x, y = (mask).nonzero()
        if origin != None:
            x += origin[0]
            y += origin[1]
        val = data[x,y]
        return x, y, val

    # ... three-d case
    if ndim == 3:
        x, y, z = (mask).nonzero()
        if origin != None:
            x += origin[0]
            y += origin[1]
            z += origin[2]
        val = data[x,y,z]
        return x, y, z, val
    
    return

def vec_to_image(
    x = None,
    y = None,
    z = None,
    val = None,
    pad = 1,
    ):
    """
    Convert a 1, 2, or 3d vector to a numpy image.
    """

    # Note dimensionality
    if y == None:
        ndim = 1
    elif z == None:
        ndim = 2
    else:
        ndim = 3

    # Get data type of input
    our_dtype = val.dtype

    # Get the extent of the vector
    minx = np.min(x)
    maxx = np.max(x)
    x0 = minx - pad
    dx = (maxx - minx + 1) + 2*pad

    if ndim >= 2:
        miny = np.min(y)
        maxy = np.max(y)
        y0 = miny - pad
        dy = (maxy - miny + 1) + 2*pad

    if ndim >= 3:
        minz = np.min(z)
        maxz = np.max(z)
        z0 = minz - pad
        dz = (maxz - minz + 1) + 2*pad

    # Work out the size of the smaller cube
    if ndim == 1:
        size_vec = [dx]
    elif ndim == 2:
        size_vec = [dx, dy]
    else:
        size_vec = [dx, dy, dz]

    # Make the smaller data array
    data = np.array(size_vec, our_dtype)

    # Fill the data structure
    if ndim == 1:
        data[x-x0] = val
        origin = (x0)
    elif ndim == 2:
        data[x-x0,y-y0] = val
        origin = (x0, y0)
    elif ndim == 3:
        data[x-x0,y-y0,z-z0] = val
        origin = (x0, y0, z0)

    # Attach extra information? / TBD

    return data, origin

# ------------------------------------------------------------
# BLOB CHARACTERIZATION
# ------------------------------------------------------------

def get_all_blob_shapes(
    color_image,
    save_coords=False):
    """
    Return statistics dictionaries for all blobs in a color image.
    """

    # note dimensions
    ndim = color_image.ndim

    # vectorize the color image
    if ndim == 1:
        x, color = blob_to_vec(color_image, mask=(color_image > 0))

    if ndim == 2:
        x, y, color = blob_to_vec(color_image, mask=(color_image > 0))

    if ndim == 3:
        x, y, z, color = blob_to_vec(color_image, mask=(color_image > 0))
    
    # note the set of unique colors
    unique_colors = np.unique(color)
        
    # initialize output
    all_shapes = {}

    # loop over colors
    for this_color in unique_colors:

        ind = (color == this_color).nonzero()
        
        if ndim==1:
            this_shape = mask_shape(x[ind])

        if ndim==2:
            this_shape = mask_shape(x[ind], y[ind])

        if ndim==3:
            this_shape = mask_shape(x[ind], y[ind], z[ind])

        # note the color in the dictionary
        this_shape["color"] = this_color

        # if requested save the coordinates in the data
        if save_coords:
            if ndim == 1:
                this_shape["x"] = x
                
            if ndim == 2:
                this_shape["x"] = x
                this_shape["y"] = y
                    
            if ndim == 3:
                this_shape["x"] = x
                this_shape["y"] = y
                this_shape["z"] = z

        # place the dictionary in the parent structure
        all_shapes[this_color] = this_shape

    # return
    return all_shapes
