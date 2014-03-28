import numpy as np
import math
from scipy.ndimage import label

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Object to hold structuring elements
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Struct:
    """ 
    ...
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # The element
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    struct = None
    ndim = None
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Build the element
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(
        self,       
        whatami="simple",
        ndim=3,
        **kwargs
        ):
        """
        Construct a new structuring object.
        """

        self.ndim = ndim

        if whatami=="simple":
            if kwargs.has_key("corners"):                
                self.struct = simple(ndim,corners=kwargs["corners"])
            else:
                self.struct = simple(ndim,corners=False)
    
        if whatami=="cylinder":
            if kwargs.has_key("major") == False:
                print "Need a major axis."
                return
            self.struct = ellipse(major=kwargs["major"])
            if ndim > 2:

                if kwargs.has_key("zaxis") == False:
                    print "Assuming z axis is axis 0."
                    zaxis=0
                else:
                    zaxis = kwargs["zaxis"]

                if kwargs.has_key("depth") == False:
                    print "Need a depth."
                    return None

                self.extend_along_axis(
                    axis=zaxis, 
                    depth=kwargs["depth"])

        if whatami=="rectangle":
            if kwargs.has_key("major") == False:
                print "Need a major axis."
                return
            self.struct = rectangle(major=kwargs["major"])
            if ndim > 2:

                if kwargs.has_key("zaxis") == False:
                    print "Assuming z axis is axis 0."
                    zaxis=0
                else:
                    zaxis = kwargs["zaxis"]

                if kwargs.has_key("depth") == False:
                    print "Need a depth."
                    return None

                self.extend_along_axis(
                    axis=zaxis, 
                    depth=kwargs["depth"])

    # ------------------------------------------------------------
    # Suppress connect along specified axis
    # ------------------------------------------------------------

    def suppress_axis(
        self,
        axis=0
        ):
        """
        """
    
        rolled_view = self.struct.swapaxes(axis,0)

        len_first_axis = rolled_view.shape[0]
        
        if len_first_axis % 2 != 1:
            print "WARNING! Expected odd structuring element. Got even."

        midplane = np.floor(len_first_axis / 2)

        for plane in np.arange(len_first_axis):
            if plane == midplane:
                continue
            rolled_view[plane,...] *= 0

    # ------------------------------------------------------------
    # Extend connection along specified axis
    # ------------------------------------------------------------

    def extend_along_axis(
        self,
        axis=0,
        depth=0
        ):
        """
        """
        
        # ... work out a new shape
        old_shape = self.struct.shape
        new_shape = [depth]
        for dim in old_shape:
            new_shape.append(dim)
        new_shape = tuple(new_shape)

        # ... resize to that shape, broadcasting will copy
        self.struct = np.resize(self.struct, new_shape)
        
        # ... move the new axis 0 to the desired position
        new_order = []
        for this_axis in range(self.ndim):
            if this_axis != 0:
                new_order.append(this_axis)
            if axis == this_axis:
                new_order.append(0)
        
        # ... transpose to the desired order
        self.struct.transpose(new_order)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DEFINE STRUCTURING ELEMENTS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# .............................
# Simple one-pixel connectivity
# .............................

def simple(
    ndim=3,
    corners=False
    ):
    """
    Return simple connectivity either with or without corners.
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

    return connect

# .............................
# Ball
# .............................

# .............................
# Ellipse connectivity
# .............................

def ellipse(
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
        minor = major

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

# .............................
# Rectangle connectivity
# .............................

def rectangle(
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
        minor = major

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
        connect = simple(mask.ndim, corners=corners)

    # Wrap the SciPy implementation
    color, ncolors = label(mask,structure=connect)

    return color

