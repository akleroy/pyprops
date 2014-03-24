# from astropy.io import fits
import numpy as np
from scipy.ndimage import binary_dilation
from scipy.ndimage import label
from scipy.ndimage import maximum_filter
from scipy.ndimage import rank_filter
from scipy.ndimage import uniform_filter
from scipy.ndimage import convolve
from blobutils import *
import time

# Routines to identify and select local maxima

# ------------------------------------------------------------
# HIGH LEVEL ROUTINES
# ------------------------------------------------------------

def all_local_max(
    data_in,
    mask=None,
    sky_halfbox = 3,
    spec_halfbox = 1,
    search_kernel = None,
    spec_axis = 0,
    roll=False,
    timer=False,
    ):
    """
    """

    # Define a default mask (finite values)
    if mask == None:
        mask = np.isfinite(data_in)
    else:
        mask *= np.isfinite(data_in)

    # Copy the data to suppress emission outside the mask
    data = data_in
    blank = (mask).nonzero()
    data[blank] = np.min(data)-1.

    # Call the local max finder
    if roll == False:
        candidates = \
            all_local_max_dilation(
            data,
            sky_halfbox = sky_halfbox,
            spec_halfbox = spec_halfbox,
            search_kernel = search_kernel,
            spec_axis = spec_axis,
            timer=timer)
    else:
        candidates = \
            all_local_max_roll(
            data,
            sky_halfbox = sky_halfbox,
            spec_halfbox = spec_halfbox,
            search_kernel = search_kernel,
            spec_axis = spec_axis,
            timer=timer)

    # Return the local maxima
    return candidates

def reject_local_max():
    pass

def read_local_max():
    pass

def write_local_max():
    pass

# ------------------------------------------------------------
# FIND LOCAL MAXIMA BY COMPARING TO A DILATING FILTER
# ------------------------------------------------------------

def all_local_max_dilation(
    data,
    sky_halfbox = 3,
    spec_halfbox = 1,
    search_kernel = None,
    spec_axis = 0,
    timer=False,
    flatokay=False,
    ):
    """
    Extract local maxima using a maximum filter. Fast, preferred option.
    """

    # Start the timer if requested
    if timer:
        start=time.time()
        full_start=time.time()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Generate the filter and handle dimensions
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if spec_axis == 0:
        uniform_size = (spec_halfbox*2+1,
                        sky_halfbox*2+1,
                        sky_halfbox*2+1)
    elif spec_axis == 1:
        uniform_size = (sky_halfbox*2+1,
                        spec_halfbox*2+1,
                        sky_halfbox*2+1)
    elif spec_axis == 2:
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
    max_image = maximum_filter(data, 
                               footprint=footprint,
                               mode="constant", 
                               cval=np.min(data))
    lmax_cube = (data == max_image)

    # Fast... but forces a square search kernel
    max_count = uniform_filter(lmax_cube*1.,
                               size=uniform_size,
                               mode="constant", cval=0.)* \
                               uniform_total
    lmax_cube *= (max_count==1)

    # This works and is clean, but it's incredibly slow:
    #
    #max_image = rank_filter(data, 
    #                        rank=-1,
    #                        footprint=footprint,
    #                        mode="constant", 
    #                        cval=np.min(data))
    #
    #rank2_image = rank_filter(data,
    #                          rank=-2,
    #                          footprint=footprint,
    #                          mode="constant", 
    #                          cval=np.min(data))
    #lmax_cube = (data == max_image)*(max_image != rank2_image)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Extract the maxima
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    lmax_coords = (lmax_cube).nonzero()

    if timer:        
        stop=time.time()
        print "Generating local max candidates took ", stop-start

    return lmax_coords

# ------------------------------------------------------------
# FIND LOCAL MAXIMA BY ROLLING THE CUBE
# ------------------------------------------------------------

def all_local_max_roll(
    data,
    sky_halfbox = 3,
    spec_halfbox = 1,
    search_kernel = None,
    spec_axis = 0,
    timer=False,
    flat_rolls=False
    ):
    """
    """

    # Start the timer if requested
    if timer:
        start=time.time()
        full_start=time.time()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle dimensions and orientation
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if spec_axis == 0:
        sky_1 = 1
        sky_2 = 2

    if spec_axis == 1:
        sky_1 = 0
        sky_2 = 2        

    if spec_axis == 2:
        sky_1 = 0
        sky_2 = 1

    sky_steps = np.arange(-sky_halfbox,sky_halfbox+1,1)
    spec_steps = np.arange(-spec_halfbox,spec_halfbox+1,1)
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Roll the data
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Initialize the output
    lmax_cube = (data == data)

    # Initialize counter
    total_rolls = len(sky_steps)*len(sky_steps)*len(spec_steps)
    roll_count = 0

    # Roll in each dimension
    for roll_sky_1 in sky_steps:
        rolled_sky_1 = np.roll(data,roll_sky_1,axis=sky_1)
        for roll_sky_2 in sky_steps:                       
            rolled_sky_2 = np.roll(rolled_sky_1,roll_sky_2,axis=sky_2)
            for roll_spec in spec_steps:

                # track progress
                roll_count += 1
                print roll_count, " out of ", total_rolls

                # don't compare to self
                if roll_sky_1 == roll_sky_2 == roll_spec == 0:
                    continue

                # check the search kernel if that's our approach
                
                # try rolling two ways - doesn't seem to be a time difference to me

                if flat_rolls == True:
                    if spec_axis == 0:
                        total_roll = roll_sky_2 + \
                            roll_sky_1*data.shape[1] + \
                            roll_spec*data.shape[1]*data.shape[2]                                                
                    elif spec_axis == 1:
                        total_roll = roll_sky_2 + \
                            roll_spec*data.shape[1] + \
                            roll_sky_1*data.shape[1]*data.shape[2]                                                
                    else:
                        total_roll = roll_spec + \
                            roll_sky_2*data.shape[1] + \
                            roll_sky_1*data.shape[1]*data.shape[2]
                    rolled = np.roll(data, total_roll)
                else:                    
                    rolled = np.roll(rolled_sky_2,roll_spec,axis=spec_axis)
                    
                lmax_cube *= (data > rolled)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Extract the maxima
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    lmax_coords = (lmax_cube).nonzero()

    if timer:        
        stop=time.time()
        print "Generating local max candidates took ", stop-start

    return lmax_coords

# ------------------------------------------------------------
# REJECT LOCAL MAXIMA ON A VARIETY OF CONDITIONS
# ------------------------------------------------------------

