import numpy as np
from scipy.ndimage import label

# ------------------------------------------------------------
# WRAPPER TO THE BLOB COLORER
# ------------------------------------------------------------

# The guts of the blob-handling process. Wraps scipy's ndimage
# function "label." Extra utility can be added by playing with the
# connectivity but the default usage is to map simply connected
# regions.

def blob_color(mask,
               corners=False,
               connect=None):
    """
    Basic blob-coloring wrapper.
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
# DEFINE CONNECTIVITY
# ------------------------------------------------------------

def connectivity(ndim=3, 
                 skip_axes=None, 
                 corners=False):
    """
    Return simple connectivity either with or without corners for n
    dimensions (default 3). Can suppress connectivity along one or
    more axes.
    """
    if corners == True:
        connect = np.ones(np.ones(ndim)*3)
    else:
        connect = np.zeros(np.ones(ndim)*3)
        center = np.ones(ndim,dtype=np.dtype('int'))
        connect[tuple(center)] = 1
        for axis in np.arange(ndim):
            pixel = np.ones(ndim)
            pixel[axis] = 0
            connect[tuple(pixel)] = 1
            pixel[axis] = 2
            connect[tuple(pixel)] = 1
    return connect

# ------------------------------------------------------------
# BLOB EXTRACTION
# ------------------------------------------------------------

# These focus on remapping blobs to vectors or dictionaries for easy
# manipulation. Definitely an area that could be either sped up or
# merged into moment calculation.

def blob_to_vec(
    data,
    mask=None
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
        val = data[x]
        return x, val

    # ... two-d case
    if ndim == 2:
        x, y = (mask).nonzero()
        val = data[x,y]
        return x, y, val

    # ... three-d case
    if ndim == 3:
        x, y, z = (mask).nonzero()
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

def stat_blob(
    x, 
    y = None, 
    z = None
    ):
    """
    """

    # Initialize the output
    blob_stats = {}

    # -=-=-=-=-=
    # 1-D extent
    # -=-=-=-=-=

    # ... X
    blob_stats["minx"] = np.min(x)
    blob_stats["maxx"] = np.max(x)
    blob_stats["deltax"] = blob_stats["maxx"]-blob_stats["minx"]+1

    # ... Y
    if y != None:
        blob_stats["miny"] = np.min(y)
        blob_stats["maxy"] = np.max(y)
        blob_stats["deltay"] = blob_stats["maxy"]-blob_stats["miny"]+1

    # ... Z
    if z != None:
        blob_stats["minz"] = np.min(z)
        blob_stats["maxz"] = np.max(z)
        blob_stats["deltaz"] = blob_stats["maxz"]-blob_stats["minz"]+1

    # -=-=-=-=-=
    # 2-D extent
    # -=-=-=-=-=    

    # Create a number that maps to a unique 2-d position out of pairs
    # of coordinates and then note the area by counting the unique
    # set of such numbers.

    # ... XY
    if y != None:
        blob_stats["areaxy"] = len(np.unique(y*(blob_stats["maxx"]+1) + x))

    if z != None:
        # ... XZ
        blob_stats["areaxz"] = len(np.unique(z*(blob_stats["maxx"]+1) + x))

        # ... YZ
        blob_stats["areayz"] = len(np.unique(z*(blob_stats["maxy"]+1) + y))

    # -=-=-=-=-=
    # 3-D extent
    # -=-=-=-=-=

    # Equate the volume to the set of pixels
    
    blob_stats["volume"] = len(x)

    # Return the dictionary
    return blob_stats

def stat_all_blobs(
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
    all_stats = {}

    # loop over colors
    for this_color in unique_colors:

        ind = (color == this_color).nonzero()
        
        if ndim==1:
            this_stats = stat_blob(x[ind])

        if ndim==2:
            this_stats = stat_blob(x[ind], y[ind])

        if ndim==3:
            this_stats = stat_blob(x[ind], y[ind], z[ind])

        # note the color in the dictionary
        this_stats["color"] = this_color

        # if requested save the coordinates in the data
        if ndim == 1:
            this_stats["x"] = x

        if ndim == 2:
            this_stats["x"] = x
            this_stats["y"] = y

        if ndim == 3:
            this_stats["x"] = x
            this_stats["y"] = y
            this_stats["z"] = z

        # place the dictionary in the parent structure
        all_stats[this_color] = this_stats

    # return
    return all_stats
