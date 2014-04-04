# Local maximum class. Indended to identify, reject, and compare local
# maxima based on an associated data set and mask.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time, copy, sys
import numpy as np

from scipy.ndimage import maximum_filter
from scipy.ndimage import uniform_filter
from scipy.ndimage import label, find_objects

import matplotlib.pyplot as plt

from pyprops import cube, mask, noise
from struct import *
from levutils import *

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
    merger_matrix = None
    
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

        # Copy the data to suppress emission outside the mask by
        # setting it to a low value (less than the original minimum of
        # the data).

        data = copy.deepcopy(self.data.data)
        low_value = np.min(data[self.data.valid])-1.
        data[self.data.valid == False] = low_value
        if self.mask != None:
            data[self.mask.data == False] = low_value

        # Generate the filter and handle dimensions

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
            (spec_halfbox*2+1)*1.

        # Apply the filter

        # Start the timer if requested
        if timer:
            start_max=time.time()

        # Make an image that's equal to the maximum
        max_image = maximum_filter(
            data, 
            footprint=footprint,
            mode="constant", 
            cval=low_value)
        lmax_cube = (data == max_image)*(data != low_value)

        if timer:
            stop_max=time.time()
            start_count=time.time()
            print "Maximum filter took ", stop_max-start_max

        # Fast... but forces a square search kernel
        max_count = uniform_filter(
            lmax_cube*uniform_total,
            size=uniform_size,
            mode="constant", cval=0.)
        lmax_cube *= (max_count == 1)

        if timer:
            stop_count=time.time()
            print "Maximum count took ", stop_count-start_count

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
        
        # Extract the maxima

        self.pix = (lmax_cube).nonzero()
        
        self.num = np.arange(len(self.pix[0]))

        if timer:        
            stop=time.time()
            print "Generating local max candidates took ", stop-start

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Calculations Involving Merger Levels
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def calc_merger(
        self,
        levels=None,
        corners=False,
        verbose = False,
        timer = False
        ):
        
        # Initialization
        if timer:
            start=time.time()
            full_start=time.time()

        if self.num == None or self.pix == None:
            print "Find candidate local maxima before calculating mergers."
            return

        # Copy the data to suppress emission outside the mask by
        # setting it to a low value (less than the original minimum of
        # the data).

        data = copy.deepcopy(self.data.data)
        if self.mask != None:
            use = self.mask.data*self.data.valid
        else:
            use = self.data.valid
        min_use = np.min(self.data.data[use])
        max_use = np.max(self.data.data[use])
        low_value = min_use-1.
        data[(use==False)] = low_value

        # Initialize default levels - if we have a noise value, treat
        # the RMS as a reasonable spacing (arguing that we cannot
        # distinguish much more than this finely). Else space 100
        # levels between the minimum and maximum value.

        if levels == None:
            if self.data.noise != None:
                levels = contour_values(
                    linspace = True,
                    maxval = max_use,
                    minval = min_use,                
                    spacing = 1.0*self.data.noise.scale
                    )
            else:
                levels = contour_values(
                    linspace = True,
                    maxval = max_use,
                    minval = min_use,                
                    nlev = 100
                    )

        # Build the connectivity structure

        structure = (Struct(
                "simple", 
                ndim=data.ndim,                
                corners=corners)).struct

        # Initialize the output

        self.merger_matrix = np.zeros((len(self.num), len(self.num)))*np.nan

        if timer:
            stop=time.time()
            print "Prep took ", stop-start
            start=time.time()

        # Loop over levels

        count = 0
        nlev = len(levels)
        for level in levels:
            
            perc = count*1./nlev
            sys.stdout.write('\r')            
            sys.stdout.write("Calculating merger for level %d out of %d" % (count, nlev))
            sys.stdout.flush()
            count += 1

            # Label this level
            thresh = (data >= level)

            labels, ncolors = label(
                thresh,
                structure=structure)

            # Get the assignments for the seeds
            seed_labels = labels[self.pix]

            # Get the number of discrete assignments
            max_label = np.max(seed_labels)

            if max_label == 0:
                continue

            # Histogram and look for merged cases
            bins = np.arange(0.5,max_label+0.5,1)
            hist_label = (np.histogram(seed_labels,bins=bins))[0]
            digi_label = np.digitize(seed_labels,bins=bins)

            multi_ind = ((hist_label > 1).nonzero())[0]

            # Loop over merged cases
            for ind in multi_ind:

                shared_seeds = (((digi_label-1) == ind).nonzero())[0]

                # Note the threshold in the matrix                        
                for seed in shared_seeds:

                    # this step (weirdly) can take a lot of time
                    self.merger_matrix[seed, shared_seeds] = level
                    self.merger_matrix[shared_seeds, seed] = level

        # Clean up diagonal
        for i in range(len(self.num)):
            self.merger_matrix[i,i] = np.nan

        # Finish
        if timer:
            stop=time.time()
            print "Merger calculations took ", stop-start
            start=time.time()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Reject Local Maxima
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def keep_lmax_subset(
        self,
        keep=None
        ):
        """
        ...
        """

        if keep == None:
            return

        if len(pix) == 3:
            new_pix = (
                self.pix[0][keep],
                self.pix[1][keep],
                self.pix[2][keep],            
                )
            if self.merger_matrix == None:                
                new_merger = self.merger_matrix[keep]
                new_merger = new_merger[:,keep]
        if len(pix) == 2:
            new_pix = (
                self.pix[0][keep],
                self.pix[1][keep],
                )
            if self.merger_matrix == None:                
                new_merger = self.merger_matrix[keep]
                new_merger = new_merger[:,keep]
        self.pix = new_pix
        self.merger = new_merger
        self.num = np.arange(len(self.pix[0]))

    def reject_on_value(
        self,
        thresh=3.0,
        snr=True,
        verbose=False
        ):
        """
        Reject maxima on absolute value.
        """

        if snr:
            data = self.data.snr()
        else:
            data = self.data.data
        
        keep = (self.num == self.num)
        i = 0
        for seed in self.num:
            coords = (
                self.pix[0][seed],
                self.pix[1][seed],
                self.pix[2][seed])
                
            if data[coords] < thresh:
                keep[i] = False
            
            i += 1
        
        self.keep_lmax_subset(keep)

    def reject_on_delta(
        self,
        thresh=3.0,
        snr=True,
        verbose=False
        ):
        """
        """

        if self.merger_matrix == None:
            print "You need to calculate the merger matrix."

        if snr:
            data = self.data.snr()
        else:
            data = self.data.data

        keep = (self.num == self.num)
        i = 0
        for seed in self.num:
            
            # TBD

            #coords = (
            #    self.pix[0][seed],
            #    self.pix[1][seed],
            #    self.pix[2][seed])
                
            # TBD

            #if data[coords] < thresh:
            #    keep[i] = False
            
            i += 1

        self.keep_lmax_subset(keep)

