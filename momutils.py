import numpy as np
from scipy.ndimage import label

# ------------------------------------------------------------
# SHAPE OF A BINARY OBJECT
# ------------------------------------------------------------

def mask_shape(
    x, 
    y = None, 
    z = None
    ):
    """
    """

    # Initialize the output
    obj_shape = {}

    # -=-=-=-=-=
    # 1-D extent
    # -=-=-=-=-=

    # ... X
    obj_shape["minx"] = np.min(x)
    obj_shape["maxx"] = np.max(x)
    obj_shape["deltax"] = obj_shape["maxx"]-obj_shape["minx"]+1

    # ... Y
    if y != None:
        obj_shape["miny"] = np.min(y)
        obj_shape["maxy"] = np.max(y)
        obj_shape["deltay"] = obj_shape["maxy"]-obj_shape["miny"]+1

    # ... Z
    if z != None:
        obj_shape["minz"] = np.min(z)
        obj_shape["maxz"] = np.max(z)
        obj_shape["deltaz"] = obj_shape["maxz"]-obj_shape["minz"]+1

    # -=-=-=-=-=
    # 2-D extent
    # -=-=-=-=-=    

    # Create a number that maps to a unique 2-d position out of pairs
    # of coordinates and then note the area by counting the unique
    # set of such numbers.

    # ... XY
    if y != None:
        obj_shape["areaxy"] = len(np.unique(y*(obj_shape["maxx"]+1) + x))

    if z != None:
        # ... XZ
        obj_shape["areaxz"] = len(np.unique(z*(obj_shape["maxx"]+1) + x))

        # ... YZ
        obj_shape["areayz"] = len(np.unique(z*(obj_shape["maxy"]+1) + y))

    # -=-=-=-=-=
    # 3-D extent
    # -=-=-=-=-=

    # Equate the volume to the set of pixels
    
    obj_shape["volume"] = len(x)

    # -=-=-=-=-=-=-=
    # Richer measures
    # -=-=-=-=-=-=-=

    # ... perimeter by plane
    # ... longest chord by plane
    # ... best fit ellipse
    # ... roundness
    # ... would cost more, too

    # Return the dictionary
    return obj_shape

# ------------------------------------------------------------
# STRUCTURE OF AN OBJECT WITH DATA VALUES
# ------------------------------------------------------------

def stat_0d():
    pass

def stat_1d():
    pass

def stat_2d():
    pass

def stat_3d():
    pass
