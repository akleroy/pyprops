# Assignment class, intended to hold object assignments. Extends the
# cube class.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
import copy
import sys
import numpy as np

import scipy.ndimage as ndimage
import matplotlib.pyplot as plt

from pyprops import cube, noise, mask, lmax
from struct import *
from levutils import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# ASSIGNMENT OBJECT
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Assign(cube.Cube):
    """
    ...
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Attributes (in addition to those in Cube)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    linked_data = None
    linked_mask = None    
    linked_lmax = None

    nclouds = None
    levels = None
    slices = None

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(
        self,
        *args,
        **kwargs
        ):
        """
        Construct a new assignment object.
        """
        self.set_linked_data(kwargs.pop("linked_data", None))
        self.set_linked_lmax(kwargs.pop("linked_lmax", None))
        self.set_linked_mask(kwargs.pop("linked_mask", None))
        cube.Cube.__init__(self, *args, **kwargs)
        self.valid = None
        self.data = None

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Copy from another cube
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Modify the lower level call to link to the data

    def init_from_mask(
        self, 
        prev):
        """
        Initialize a new cube from another cube. Copy the data.
        """
        cube.Cube.init_from_cube(self, prev)
        self.linked_mask = prev

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Links to data cube
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def set_linked_data(
        self,
        val=None
        ):
        """
        Link the mask object to a data cube object.
        """
        if val != None:
            self.linked_data = val

    def set_linked_mask(
        self,
        val=None
        ):
        """
        Link the mask object to a data cube object.
        """
        if val != None:
            self.linked_mask = val

    def set_linked_lmax(
        self,
        val=None
        ):
        """
        Link the mask object to a data cube object.
        """
        if val != None:
            self.linked_lmax = val

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Backup/undo
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def step_back(
        self):
        """
        Restore the backup mask, setting it to be the new mask.
        """
        if self.backup != None:
            self.data = self.backup

    def save_backup(
        self):
        """
        Restore the backup mask, setting it to be the new mask.
        """
        self.backup = self.data

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read/write
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # ... inherits cube, write is fine

    # ... read needs to get the number of clouds and calculate the slices

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manipulate Assignment
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # ... merge two clouds

    # ... delete
    
    # ... renumber to maximum compactness

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Island Assignment
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def islands(
        self,
        corners=False
        ):
        """
        Generate an assignment mask using simple connectedness.
        """

        if self.linked_mask == None:
            print "Island assignment requires a mask."
            return
                
        structure = (Struct(
                "simple", 
                ndim=self.linked_mask.data.ndim,                
                corners=corners)).struct

        labels, nlabels = ndimage.label(self.linked_mask.data,
                                        structure=structure)
        
        self.data = labels
        self.nclouds = nlabels
        self.slices = ndimage.find_objects(self.data)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Clumpfind Assignment
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def clumpfind(
        self,
        levels=None,
        corners=False,
        seeded=False,
        allow_new_peaks=True,
        timer=True
        ):
        """
        Generate a clumpfind assignment mask.
        """

        # ...................................................
        # Check user options
        # ...................................................

        if self.linked_data == None:
            print "Clumpfind assignment requires data."
            return
        
        if seeded == True:
            if self.linked_lmax == None:
                print "Seeded clumpfind assignment requires local maxima."
                return
        
        if seeded == False and allow_new_peaks == False:
            print "Cannot run an unseeded (classic) clumpfind without being able to add seeds."
            return

        # ...................................................
        # Get data to use
        # ...................................................            

        # Get the data and set the values we will not use to a low
        # number that will be ignored by the algorithm.

        data = copy.deepcopy(self.linked_data.data)
        if self.linked_mask != None:
            use = self.linked_mask.data*self.linked_data.valid
        else:
            use = self.linked_data.valid
        min_use = np.min(self.linked_data.data[use])
        max_use = np.max(self.linked_data.data[use])
        low_value = min_use-1.
        data[(use==False)] = low_value

        # ...................................................
        # Calculate contour levels
        # ...................................................

        if levels == None:
            if self.linked_data.noise != None:
                print "Defaulting to 2 sigma spacing."
                levels = contour_values(
                    linspace = True,
                    maxval = max_use,
                    minval = min_use,                
                    spacing = 2.0*self.linked_data.noise.scale
                    )
            else:
                print "Need a noise estimate."
                return

        self.levels = levels

        # ...................................................
        # Build the structuring element
        # ...................................................

        structure = (Struct(
                "simple", 
                ndim=self.linked_data.data.ndim,                
                corners=corners)).struct

        # ...................................................
        # Initialize the output
        # ...................................................

        # ... data
        self.data = np.zeros_like(data, dtype=np.int)

        # ... local maxima
        if seeded == False:
            print "Initializing a new set of local maxima"
            self.linked_lmax = \
                lmax.Lmax(self.linked_data, self.linked_mask)

        # ...................................................
        # Loop over levels (from high to low)
        # ...................................................

        nlev = len(levels)
        count = 0

        for level in levels:            

            # ........................
            # Print a counter
            # ........................

            perc = count*1./nlev
            sys.stdout.write('\r')            
            sys.stdout.write("Clumpfind level %d out of %d" % (count, nlev))
            sys.stdout.flush()
            count += 1

            # ............................
            # Label regions for this level
            # ............................

            thresh = (data >= level)
            labels, ncolors = ndimage.label(
                thresh,
                structure=structure)
            
            # ...........................
            # Vectorize the labeled data
            # ...........................

            # This gives a big speedup for sparse data.

            ind = np.where(thresh)
            val = self.linked_data.data[ind]
            ind_arr = cube.xyztup_to_array(ind, coordaxis=1)
            label_vec = labels[ind]

            # Get the assignments for the current seeds
            if self.linked_lmax.num > 0:
                seed_labels = labels[self.linked_lmax.as_tuple()]
            
            # ........................................
            # Loop over discrete regions at this level
            # ........................................

            for label in range(1,ncolors+1):
                
                # ........................................
                # Get the indices for this region
                # ........................................

                this_color = np.where(label_vec == label)
                this_val = val[this_color]
                this_ind_arr = ind_arr[this_color[0],:]
                this_ind = cube.xyzarr_to_tuple(this_ind_arr,coordaxis=1)

                # ........................................
                # Check if we should add a new peak
                # ........................................

                # If there are no peaks or if there are no peaks in
                # this region, we want to add a new one --- but only
                # if that's allowed! 

                # A future extension is to add additional criteria
                # that must be met to add a peak (volume, area, etc.)

                if self.linked_lmax.num == 0:
                    if allow_new_peaks:
                        add_a_new_peak = True
                    else:
                        continue
                elif np.sum(seed_labels == label) == 0:
                    if allow_new_peaks:
                        add_a_new_peak = True
                    else:
                        continue
                else:
                    add_a_new_peak = False
                    
                # ........................................
                # Add a new peak
                # ........................................

                if add_a_new_peak:

                    # Find the location of the maximum value
                    maxind = np.argmax(this_val)

                    # Get the corresponding coordinates
                    peak_index = this_ind_arr[maxind,:]

                    # Add a local maximum
                    new_name = self.linked_lmax.add_local_max(peak_index)

                    # Label these data in the assignment cube
                    self.data[this_ind] = new_name

                    continue

                # ........................................
                # Deal with the case of a signle seed
                # ........................................

                if np.sum(seed_labels == label) == 1:
                    
                    maxind = np.where((seed_labels == label))

                    self.data[this_ind] = self.linked_lmax.name[maxind]

                    continue

                # ........................................
                # Deal with the case of competing seeds
                # ........................................

                # Several matching labels
                if np.sum(seed_labels == label) > 1:

                    # Initialize an assignment vector
                    this_assign = np.zeros_like(this_val)
                    best_dist = np.zeros_like(this_val)

                    # Identify the competing seeds
                    maxind = np.where((seed_labels == label))

                    n_max = len(maxind[0])

                    for i in range(n_max):
                        
                        this_max_name = self.linked_lmax.name[maxind[0][i]]

                        this_max_coord = self.linked_lmax.indices[this_max_name-1]

                        dist_to_this_max = \
                            np.sum((this_ind_arr - this_max_coord)**2,axis=1)
                        
                        if i == 0:
                            # ... all true for the first test
                            is_closest = (dist_to_this_max == dist_to_this_max)
                        else:
                            is_closest = (dist_to_this_max < best_dist)

                        this_assign[is_closest] = this_max_name
                        best_dist[is_closest] = dist_to_this_max[is_closest]


                    self.data[this_ind] = this_assign

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # CPROPS Assignment
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

