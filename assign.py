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

    # Seeded

    def clumpfind_orig(
        self,
        levels=None,
        corners=False,
        timer=True
        ):
        """
        Generate a clumpfind assignment mask.
        """

        if self.linked_data == None:
            print "Clumpfind assignment requires data."
            return
        
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

        # Calculate contour levels to use

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

        print levels

        # Build the structuring element

        structure = (Struct(
                "simple", 
                ndim=self.linked_data.data.ndim,                
                corners=corners)).struct

        # Initialize the output

        # ... data
        self.data = np.zeros_like(data, dtype=np.int)
        # ... local maxima
        self.linked_lmax = \
            lmax.Lmax(self.linked_data, self.linked_mask)

        # Loop over levels (from high to low)
        next_cloud = 1
        nlev = len(levels)
        count = 0
        for level in levels:            

            # Counter
            perc = count*1./nlev
            sys.stdout.write('\r')            
            sys.stdout.write("Clumpfind level %d out of %d" % (count, nlev))
            sys.stdout.flush()
            count += 1

            # Label this level
            thresh = (data >= level)
            labels, ncolors = ndimage.label(
                thresh,
                structure=structure)
            
            # Vectorize
            ind_vec = np.where(thresh)
            val_vec = self.linked_data.data[ind_vec]
            ind_vec_arr = np.vstack(ind_vec).transpose()
            label_vec = labels[ind_vec]
            
            # Get the assignments for the current seeds
            if self.linked_lmax.num > 0:
                seed_labels = labels[self.linked_lmax.as_tuple()]
            else:
                noseeds = True
            
            # Slow, improve later
            for label in range(1,ncolors+1):
                
                this_ind = np.where(label_vec == label)
                sub_val_vec = val_vec[this_ind]                
                sub_ind_vec_arr = ind_vec_arr[this_ind[0],:]
                sub_ind = tuple((ind_vec[i])[this_ind] for i in range(ind_vec_arr.shape[1]))

                # No matching labels
                if noseeds == True:
                    add_new_peak = True
                elif np.sum(seed_labels == label) == 0:
                    add_new_peak = True
                else:
                    add_new_peak = False
                    
                if add_new_peak:
                    # find the 
                    maxind = np.where(sub_val_vec == np.max(sub_val_vec))
                    # ... keep only first one
                    maxind = maxind[0]
                    if len(maxind) > 1:
                        maxind = maxind[0]
                    maxind = int(maxind)
                    peak_index = sub_ind_vec_arr[maxind,:]
                    self.linked_lmax.add_local_max(peak_index)
                    self.data[sub_ind] = \
                        np.max(self.linked_lmax.name)
                    continue

                # One matching labels, label
                if np.sum(seed_labels == label) == 1:
                    self.data[sub_ind] = \
                        self.linked_lmax.name[np.where((seed_labels == label))]
                    continue

                # Several matching labels
                if np.sum(seed_labels == label) > 1:
                    pass

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # CPROPS Assignment
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

