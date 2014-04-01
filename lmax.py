# Local maximum class. Indended to identify, reject, and compare local
# maxima based on an associated data set and mask.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
import copy
import numpy as np

from scipy.ndimage import maximum_filter
from scipy.ndimage import uniform_filter

import matplotlib.pyplot as plt

from pyprops import cube, mask, noise
from struct import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# LOCAL MAXIMA LIST OBJECT
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# The maxima are treated as a set held inside a single object

class Lmax():
    """
    ...
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Attributes
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    num = None
    pix = None
    merger = None
    
    data = None
    mask = None

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(
        self,
        data = None,
        mask = None
        ):
        """
        Construct a new local maximum object.
        """
        if data != None:
            self.link_to_data(data)
        if mask != None:
            self.link_to_mask(mask)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Links to data cube
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def link_to_data(
        self,
        val=None
        ):
        """
        Link the lmax object to a data object.
        """
        if val != None:
            self.data = val

    def link_to_mask(
        self,
        val=None
        ):
        """
        Link the lmax object to a mask object.
        """
        if val != None:
            self.mask = val

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read/write
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Find local maxima
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # There are several ways to do this: rank filters, rolling the
    # cube. The quickest numpy way seems to be a maximum filter
    # followed by a uniform filter. The downside of this is that it
    # only allows square search kernels. The rank filter way to do it
    # is included below (commented out) but is very slow.

    def all_local_max(
        self,
        sky_halfbox = 3,
        spec_halfbox = 1,
        search_kernel = None,
        timer=False,
        ):
        """
        Extract local maxima using a maximum filter.
        """
        

        # Start the timer if requested
        if timer:
            start=time.time()
            full_start=time.time()

        # Copy the data to suppress emission outside the mask by
        # setting it to a low value (less than the original minimum of
        # the data).
        data = copy.deepcopy(self.data.data)
        low_value = np.min(data)-1.
        data[self.data.valid == False] = low_value
        if self.mask != None:
            data[mask == False] = low_value

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Generate the filter and handle dimensions
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if self.data.spec_axis == 0:
            uniform_size = (spec_halfbox*2+1,
                            sky_halfbox*2+1,
                            sky_halfbox*2+1)

        if self.data.spec_axis == 1:
            uniform_size = (sky_halfbox*2+1,
                            spec_halfbox*2+1,
                            sky_halfbox*2+1)

        if self.data.spec_axis == 2:
            uniform_size = (sky_halfbox*2+1,
                            sky_halfbox*2+1,
                            spec_halfbox*2+1)

        footprint = np.ones(uniform_size)
        uniform_total = (sky_halfbox*2+1)* \
            (sky_halfbox*2+1)* \
            (spec_halfbox*2+1)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Apply the filter
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        # Make an image that's equal to the maximum
        max_image = maximum_filter(
            data, 
            footprint=footprint,
            mode="constant", 
            cval=np.min(data))
        lmax_cube = (data == max_image)

        # Fast... but forces a square search kernel
        max_count = uniform_filter(
            lmax_cube*1.,
            size=uniform_size,
            mode="constant", cval=0.)* \
            uniform_total
        lmax_cube *= (max_count==1)

        # This works and is clean and arbitrary-shaped, but it's
        # incredibly slow:

        # max_image = rank_filter(data, 
        #                        rank=-1,
        #                        footprint=footprint,
        #                        mode="constant", 
        #                        cval=np.min(data))
        
        # rank2_image = rank_filter(data,
        #                          rank=-2,
        #                          footprint=footprint,
        #                          mode="constant", 
        #                          cval=np.min(data))
        # lmax_cube = (data == max_image)*(max_image != rank2_image)
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Extract the maxima
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
        
        self.pix = (lmax_cube).nonzero()
        
        self.num = np.arange(len(self.pix[0]))

        if timer:        
            stop=time.time()
            print "Generating local max candidates took ", stop-start


