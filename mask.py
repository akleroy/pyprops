# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
import copy
import numpy as np
from scipy.ndimage import histogram
from scipy.ndimage import binary_dilation
from scipy.ndimage import binary_erosion
from scipy.ndimage import label, find_objects
import matplotlib.pyplot as plt

from struct import *
from cube import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MASK OBJECT
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Mask:
    """
    ...
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Contents
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    data = None
    mask = None
    backup = None
    spec_axis = None
    deg_axis = None
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize and infrastructure
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(
        self,
        data=None,
        spec_axis=None
        ):
        """
        Construct a new mask object.
        """
        self.data = data
        if data != None:
            try:
                self.spec_axis = data.spec_axis
            except NameError:
                self.spec_axis = spec_axis
        else:
            self.spec_axis = spec_axis

    def set_data(
        self,
        val=None
        ):
        """
        Link the mask object to a data object.
        """
        if val != None:
            self.data = val

    def set_spectral_axis(
        self,
        val=None
        ):
        """
        Set the spectral axis.
        """
        if val != None:
            self.spec_axis = val

    def step_back(
        self):
        """
        Restore the backup mask, setting it to be the new mask.
        """
        if self.backup != None:
            self.mask = self.backup

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read/write
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def from_casa_image(
        self,
        infile=None,
        append=False):
        """
        Read a mask from a CASA image file. Any value above 0.5 will
        be treated as True in the resulting mask.
        """
        
        # ... use the method inside the cube module. Keep only the
        # first thing returned
        new_mask, new_valid, new_cs = \
            (from_casa_image(
                infile=infile,
                transpose=False,
                dropdeg=True)
             
        # ... read the data and cast as boolean
        if append:
            self.mask *= (new_mask > 0.5)
        else:
            self.mask = (new_mask > 0.5)
 
        self.cs = new_cs
        ia.close()

        return

    def to_casa_image(
        self,
        outfile=None,
        overwrite=False,        
        template=None):
        """
        Write a mask to a CASA image file.
        """

        if template == None:
            try:
                template = self.filename
            except AttributeError:
                try:
                    template = self.data.filename
                except AttributeError:
                    print "Need a template in order to write."
                    return
        
        to_casa_image(
            data = self.mask*1.,
            template = template,
            outfile = outfile,
            overwrite = overwrite
            )

    def from_fits_file(self,
                       infile=None,
                       append=False):
        """
        Read a mask from a FITS file.
        """

        if astropy_ok == False:
            print "Cannot read FITS files without astropy."
            return

        # ... open the file
        hdulist = fits.open(infile)

        # ... read the data and cast as boolean
        if append:
            self.mask *= (hdulist[0].data > 0.5)
        else:
            self.mask = hdulist[0].data > 0.5

        # ... save the header
        self.hdr = hdulist[0].header

        return

    def to_fits_file(self,
                     outfile=None,
                     overwrite=False):
        """
        Write the cube to a FITS file.
        """


        if astropy_ok == False:
            print "Cannot write FITS files without astropy."
            return

        if self.data.filemode != "astropy":
            print "Can only/read write FITS in astropy mode."
            return

        try:
            hdr_copy = self.hdr.copy()
        except NameError:
            hdr_copy = self.data.hdr.copy()

        hdr_copy["BUNIT"] = "Mask"

        fits.writeto(
            outfile, 
            self.mask.astype(np.int16), 
            hdr_copy, 
            clobber=overwrite)
     
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Expose the mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def twod(
        self,
        axis=None):
        """
        """
        if axis==None:
            if self.spec_axis==None:
                if self.data.spec_axis!=None:
                    axis=self.data.spec_axis
            else:
                axis=self.spec_axis

        if self.mask.ndim == 2:
            return self.mask

        if axis == None:
            return None
        
        return (np.sum(self.mask, axis=axis) >= 1)
        

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manipulate the mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def pare_on_volume(
        self,
        thresh=None,
        corners=False,
        timer=False,
        backup=False
        ):
        """
        Remove discrete regions that do not meet a pixel-volume
        threshold. Requires a pixel threshold be set.
        """

        # .............................................................
        # Error checking
        # .............................................................

        if thresh==None:
            print "Need a threshold."
            return

        # .............................................................
        # Back up the mask first if requested
        # .............................................................

        if backup==True:
            self.backup=self.mask

        # .............................................................
        # Label the mask
        # .............................................................
        
        structure = (Struct(
                "simple", 
                ndim=self.mask.ndim,                
                corners=corners)).struct

        labels, nlabels = label(self.mask,
                                structure=structure)

        # .............................................................
        # Histogram the labels
        # .............................................................

        hist = histogram(
            labels, 0.5, nlabels+0.5, nlabels)
        
        # .............................................................
        # Identify the low-volume regions
        # .............................................................
        
        if np.sum(hist < thresh) == 0:
            return
        
        loc = find_objects(labels)

        for reg in np.arange(1,nlabels):
            if hist[reg-1] > thresh:
                continue
            self.mask[loc[reg-1]] *= (labels[loc[reg-1]] != reg)

#        for reg in below_thresh:
#            self.mask[(labels == reg)] = False

            
    def erode_small_regions(
        self,        
        major=3,
        depth=2,
        timer=False,
        backup=False
        ):
        """
        Use 'morphological opening' (erosion followed by dilation) to
        remove small regions from the mask.
        """

        # .............................................................
        # Back up the mask first if requested
        # .............................................................

        if backup==True:
            self.backup=self.mask

        # .............................................................
        # Time the operation if requested.
        # .............................................................

        if timer:
            start=time.time()
            full_start=time.time()

        # .............................................................
        # Construction of structuring element
        # .............................................................        

        structure = Struct(
            "rectangle", 
            major=major, 
            zaxis=self.spec_axis, 
            depth=depth)
        
        # .............................................................
        # Erosion
        # .............................................................

        self.mask = binary_erosion(
            self.mask, 
            structure=structure.struct,
            iterations=1
            )

        # .............................................................
        # Dilation
        # .............................................................

        self.mask = binary_erosion(
            self.mask, 
            structure=structure.struct,
            iterations=1
            )

        # .............................................................
        # Finish timing
        # .............................................................

        if timer:
            full_stop=time.time()
            print "Small region suppression took ", full_stop-full_start

        return

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Grow the mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def grow(
        self,
        iters=-1,
        xy_only=False,
        z_only=False,
        corners=False,
        constraint=None,
        timer=False,
        verbose=False,
        backup=False
        ):
        """
        Manipulate an existing mask. Mostly wraps binary dilation operator
        in scipy with easy flags to create structuring elements.
        """

        # .............................................................
        # Back up the mask first if requested
        # .............................................................

        if backup==True:
            self.backup=self.mask

        # .............................................................
        # Time the operation if requested.
        # .............................................................

        if timer:
            start=time.time()
            full_start=time.time()

        # .............................................................
        # Construct the dilation structure (calls "connectivity" in the
        # blobutils).
        # .............................................................
   
        skip_axes = []

        # ... if 2d-only then blank the spectral axis (if there is
        # one) in the connectivity definition.

        if xy_only == True:
            if self.spec_axis != None:
                skip_axes.append(self.spec_axis)

        # ... if 1d-only then blank the position axes (if they exist)
        # in the connectivity definition.

        if z_only == True:
            axes = range(self.mask.ndim)
            for axis in axes:
                if axis != self.spec_axis:
                    skip_axes.append(axis)

        # ... build the sturcturing element

        structure = Struct(
            "simple", 
            ndim=self.mask.ndim,                
            corners=corners)
        for skip in skip_axes:
            structure.suppress_axis(skip)

        # .............................................................
        # Apply the binary dilation with the constructed parameters
        # .............................................................

        self.mask = binary_dilation(
            self.mask, 
            structure=structure.struct,
            iterations=iters,
            mask=constraint,
            )

        # .............................................................
        # Finish timing
        # .............................................................

        if timer:
            full_stop=time.time()
            print "Mask expansion took ", full_stop-full_start

        return

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Generate a new mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def joint_threshold(
        self,
        usesnr=True,
        scale=1.0,
        thresh=4.0,
        nchan=2,
        outof=None,
        append=False,
        backup=False,
        timer=False
        ):
        """
        Masking using joint velocity channel conditions.
        
        usesnr (default True) : use a signal-to-noise cube if one can be 
        derived from the data (requires associated data and noise).

        scale (default 1.0) : a factor used to scale the threshold. Set by
        default to 1.0, appropriate for a threshold in S/N units being
        applied to mask a signal-to-noise cube. Set it to

        threshold (default 4.0) : the threshold (times the scale) that
        must be exceeded for a pixel to be included in the mask.
                
        nchan (default 2) : number of channels that must be above the
        specified threshold in order for a region to be included in the
        mask.
        
        out_of (default to nchan) : relaxes the requirement for the number
        of channels searched at once for emission. This can be used to
        allow for noisy / spiky spectra. E.g., require 3 channels out of 5
        to match the threshold with out_of=5, nchan=3. Note that then all
        five channels will included in the final mask.
                
        Defaults to requiring two out of two channels (nchan=2,
        out_of=None) above 4.0 (thresh=4.0), assuming that it has been fed
        a signal-to-noise mask (scale=1.0).
        """

        # .............................................................
        # Back up the mask first if requested
        # .............................................................

        if backup==True:
            self.backup=self.mask

        # .............................................................
        # Time the operation if requested.
        # .............................................................

        if timer:
            start=time.time()
            full_start=time.time()

        # .............................................................
        # Set defaults and catch errors
        # .............................................................

        # default out_of to nchan
        if outof==None:
            outof = nchan

        # catch error case
        if outof < nchan:
            outof = nchan 

        # .............................................................
        # Get the data that we will work with
        # .............................................................

        if usesnr:
            working_data = self.data.snr()
        else:
            working_data = self.data.data

        # .............................................................
        # Build the mask
        # .............................................................

        # initial mask set by threshold
        new_mask = np.int_(working_data >= thresh*scale)

        # If we have a spectral axis apply the joint conditions
        if self.spec_axis != None and self.data.data.ndim > 2:

            # roll the cube "out_of" times along the spectral axis and keep a
            # running tally of the number of points above the threshold by
            # summing mask.
            for i in np.arange(outof):
                new_mask += np.roll(new_mask,i,self.spec_axis)
    
            # keep only points in the mask which meet the "nchan" criteria
            new_mask = (new_mask >= nchan)

            # roll the mask in the other direction to ensure that all points
            # that contributed to the valid point are included in the final
            # mask
            for i in np.arange(outof):
                new_mask += np.roll(new_mask,-i,self.spec_axis)

            # calculate the final mask, adding a finite check, now a bool
            new_mask = (new_mask >= 1)

        # .............................................................
        # Append or replace
        # .............................................................

        if append:
            self.mask *= new_mask
        else:
            self.mask = new_mask

        # .............................................................
        # Finish timing
        # .............................................................

        if timer:
            full_stop=time.time()
            print "Joint thresholding took ", full_stop-full_start

        # return
        return

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # The CPROPS Recipe
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def cprops_mask(
        self,
        hithresh=4.0,
        lothresh=2.0,
        nchan=2,
        usesnr=True,
        scale=1.0,
        corners=False,
        append=False,
        backup=False,
        timer=False        
        ):
        pass
        
        # .............................................................
        # Back up the mask first if requested
        # .............................................................

        if backup==True:
            self.backup=self.mask

        # .............................................................
        # Time the operation if requested.
        # .............................................................

        if timer:
            start=time.time()
            full_start=time.time()

        # .............................................................
        # Build the mask
        # .............................................................

        inner_mask = Mask(self.data)
        outer_mask = Mask(self.data)
        inner_mask.joint_threshold(            
            usesnr=usesnr,
            scale=scale,
            thresh=hithresh,
            nchan=nchan,
            append=append,
            timer=timer
            )

        inner_mask.erode_small_regions(
            major=3,
            depth=2,
            timer=timer)
            
        outer_mask.joint_threshold(            
            usesnr=usesnr,
            scale=scale,
            thresh=lothresh,
            nchan=nchan,
            append=append,
            timer=timer
            )

        inner_mask.grow(
            corners=corners,
            constraint=outer_mask.mask,
            timer=timer
            )
            
        # .............................................................
        # Append or replace
        # .............................................................

        if append:
            self.mask *= inner_mask.mask
        else:
            self.mask = inner_mask.mask

        # .............................................................
        # Finish timing
        # .............................................................

        if timer:
            full_stop=time.time()
            print "CPROPS-style masking took ", full_stop-full_start

        # return
        return

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Visualize mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def contour_on_peak(
        self,
        scale=10
        ):

        plt.figure()

        map = self.data.peak_map()

        vmax = np.max(map[np.isfinite(map)])
        vmin = 0.0
        if scale != None:
            if self.data.noise != None:
                if self.data.noise.scale != None:
                    vmax = scale*self.data.noise.scale
                    vmin = 0.0

        plt.imshow(
            map,
            vmin=vmin,
            vmax=vmax,
            origin='lower')
        plt.contour(
            self.twod(),
            linewidths=0.5,
            colors='white'
            )
        plt.show()
        
