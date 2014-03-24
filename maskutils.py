# from astropy.io import fits
import numpy as np
from scipy.ndimage import binary_dilation
from scipy.ndimage import label
from blobutils import *
import time

# Routines to generate and manipulate masks. Heavily leverages
# scipy.ndimage, some via the related blobutils module.

# ------------------------------------------------------------
# MAIN EXPOSED FUNCTION
# ------------------------------------------------------------

def make_cprops_mask(
    cube,
    prior=None,
    hi_thresh=4.0,
    lo_thresh=2.0,
    scale=1.0,
    n_chan=2,
    out_of=None,
    spec_axis=2,  
    min_vol=None,
    min_area=None,
    min_vchan=None,
    timer=False,
    verbose=False
    ):
    """
    General purpose masking routine. Identify regions above a given
    threshold across several channels. Reject regions based on volume,
    area, and velociy extent. Returns a mask matchedin size to the
    input image.

    Rather than wrap a read/writer or signal-to-noise cube generator
    in here, these are held externally.

    The scale parameter gives the option to invert the cube by using
    -1 or to scale, e.g., by the RMS.
    """

    # ------------------------------------------------------------
    # Generate the candidate mask
    # ------------------------------------------------------------

    if timer:
        start=time.time()
        full_start=time.time()

    # High threshold
    hi_mask = make_spec_mask(
        cube,
        thresh=hi_thresh,
        scale=scale,
        n_chan=n_chan,
        out_of=out_of,
        spec_axis=spec_axis)

    if timer:
        stop=time.time()
        print "High threshold mask took ", stop-start
        start=time.time()

    # Low threshold
    lo_mask = make_spec_mask(
        cube,
        thresh=lo_thresh,
        scale=scale,
        n_chan=n_chan,
        out_of=out_of,
        spec_axis=spec_axis)

    if timer:
        stop=time.time()
        print "Low threshold mask took ", stop-start
        start=time.time()
                   
    # Expand regions that meet the high threshold to includ all
    # connected low threshold regions
    mask = binary_dilation(hi_mask, mask=lo_mask, iterations=-1)

    if timer:
        stop=time.time()
        print "Mask dilation took ", stop-start
        start=time.time()

    # ------------------------------------------------------------
    # Pare small regions
    # ------------------------------------------------------------

    # Note that this step could be accomplished to some degree by
    # using a (binary) opening operator. It's not perfect, it will
    # clip your edges but it might be much faster and better defined.

    if min_area != None or min_vol != None or min_vchan != None:
        reject_regions(mask,
                       min_area=min_area,
                       min_vol=min_vol,
                       min_vchan=min_vchan,
                       spec_axis=spec_axis,
                       inplace=True,
                       timer=timer)

    if timer:
        stop=time.time()
        print "Small region rejection took ", stop-start

    # ------------------------------------------------------------
    # Apply any prior mask supplied
    # ------------------------------------------------------------

    if prior != None:
        mask *= prior

    # ------------------------------------------------------------
    # Return final mask
    # ------------------------------------------------------------

    if timer:
        full_stop=time.time()
        print "Whole masking took ", full_stop-full_start

    return mask

# ------------------------------------------------------------
# MASK-MAKING FUNCTIONS
# ------------------------------------------------------------

def make_spec_mask(
    cube,
    thresh=4.0,
    scale=1.0,              
    n_chan=2,
    out_of=None,
    spec_axis=0):
    """
    Masking using joint velocity channel conditions.

    threshold (default 4.0) : the threshold (times the scale) that
    must be exceeded for a pixel to be included in the mask.

    scale (default 1.0) : a factor used to scale the threshold. Set by
    default to 1.0, appropriate for a threshold in S/N units being
    applied to mask a signal-to-noise cube. Set it to

    n_chan (default 2) : number of channels that must be above the
    specified threshold in order for a region to be included in the
    mask.

    out_of (default to n_chan) : relaxes the requirement for the number
    of channels searched at once for emission. This can be used to
    allow for noisy / spiky spectra. E.g., require 3 channels out of 5
    to match the threshold with out_of=5, n_chan=3. Note that then all
    five channels will included in the final mask.

    spec_axis (default to 0) : the axis to roll along

    Defaults to requiring two out of two channels (n_chan=2,
    out_of=None) above 4.0 (thresh=4.0), assuming that it has been fed
    a signal-to-noise mask (scale=1.0). Assumes that the spectral axis
    is the second axis (spec_axis=1).
    """

    # default out_of to n_chan
    if out_of==None:
        out_of = n_chan

    # catch error case
    if out_of < n_chan:
        out_of = n_chan 

    # initial mask set by threshold
    mask = np.int_(cube >= thresh*scale)

    # roll the cube "out_of" times along the spectral axis and keep a
    # running tally of the number of points above the threshold by
    # summing mask.
    for i in np.arange(out_of):
        mask += np.roll(mask,i,spec_axis)
    
    # keep only points in the mask which meet the "n_chan" criteria
    mask = (mask >= n_chan)

    # roll the mask in the other direction to ensure that all points
    # that contributed to the valid point are included in the final
    # mask
    for i in np.arange(out_of):
        mask += np.roll(mask,-i,spec_axis)

    # calculate the final mask, adding a finite check, now a bool
    mask = (mask >= 1)*np.isfinite(cube)

    # return
    return mask

# ------------------------------------------------------------
# MASK MANIPULATION
# ------------------------------------------------------------

def grow_mask(mask_in,
              iters=-1,
              xy_only=False,
              z_only=False,
              corners=False,
              spec_axis=0,
              constraint=None,
              timer=False,
              verbose=False
              ):
    """
    Manipulate an existing mask. Mostly wraps binary dilation operator
    in scipy with easy flags to create structuring elements.
    """
    if timer:
        start=time.time()
        full_start=time.time()

    # Construct the dilation structure (calls "connectivity" in the
    # blobutils).
   
    skip_axes = []

    if xy_only == True:
        skip_axes.append(spec_axis)

    if z_only == True:
        axes = range(mask_in.ndim)
        for axis in axes:
            if axis != spec_axis:
                skip_axes.appen(axis)

    structure = connectivity(ndim=mask_in.ndim,
                             skip_axes=skip_axes,
                             corners=corners) 

    print structure

    # Apply the dilation

    mask = binary_dilation(mask_in, 
                           structure=structure,
                           iterations=iters,
                           mask=constraint,
                           )

    # Return
    
    if timer:
        full_stop=time.time()
        print "Mask expansion took ", full_stop-full_start

    return mask

def mask2d_to_3d():
    pass

def mask3d_to_2d():
    pass

def reject_regions(mask,
                   assign=None,
                   min_area=None,
                   min_vol=None,
                   min_vchan=None,
                   spec_axis=2,
                   inplace=True,
                   timer=False):
    """
    Reject small regions from a mask. Some of this functionality could
    be replicated at increased speed by using binary erosion (with a
    kernel customized according to the axes).
    """

    if timer:
        start=time.time()

    ndim = mask.ndim

    # Identify discrete regions in the mask if not already provided
    if assign == None:
        assign = blob_color(mask)

    if timer:
        stop=time.time()
        print "... blob coloring took ", stop-start
        start=time.time()

    # Derive statistics for the regions
    shapes = get_all_blob_shapes(assign, save_coords=True)

    if timer:
        stop=time.time()
        print "... blob statistics took ", stop-start
        start=time.time()

    # Initialize a new mask
    if inplace == False:
        new_mask = mask

    # Loop over all regions
    for this_color in shapes.keys():  

        blob_shape = shapes[this_color]

        reject = False

        # Reject on volume
        if min_vol != None and reject == False:
            if blob_shape["volume"] < min_vol:
                reject = True
        
        # Reject on velocity extent
        if min_vchan != None and reject == False:

            # ... switch on spectral axis
            if spec_axis == 0:
                if blob_shape["deltax"] < min_vchan:
                    reject = True                
            elif spec_axis == 1:
                if blob_shape["deltay"] < min_vchan:
                    reject = True                
            elif spec_axis == 2:
                if blob_shape["deltaz"] < min_vchan:
                    reject = True

        # Reject on areal extent
        if min_area != None and reject == False:

            # ... switch on spectral axis
            if spec_axis == 2:
                if blob_shape["areaxy"] < min_area:
                    reject = True                    
            elif spec_axis == 0:
                if blob_shape["areayz"] < min_area:
                    reject = True                
            elif spec_axis == 1:
                if blob_shape["areaxz"] < min_area:
                    reject = True
                    
        # Pare a rejected region from the mask
        if reject == True:
            if inplace == True:
#                mask *= (assign!=this_blob)
                if ndim == 1:
                    mask[blob_shape["x"]] = False
                if ndim == 2:
                    mask[blob_shape["x"],blob_shape["y"]] = False
                if ndim == 3:
                    mask[blob_shape["x"],blob_shape["y"],blob_shape["z"]] = False
            else:
                if ndim == 1:
                    new_mask[blob_shape["x"]] = False
                if ndim == 2:
                    new_mask[blob_shape["x"],blob_shape["y"]] = False
                if ndim == 3:
                    new_mask[blob_shape["x"],blob_shape["y"],blob_shape["z"]] = False
#                new_mask *= (assign!=this_blob)
#                new_mask[assign==this_blob] = False

    if timer:
        stop=time.time()
        print "... blob parsing took ", stop-start

    # Return the new mask
    if inplace == False:
        return new_mask
    else:
        return

