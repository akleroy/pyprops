import numpy as np
from scipy.ndimage import label
from scipy.ndimage import histogram
from blobutils import connectivity
from blobutils import get_all_blob_shapes
import time

# Utilities to deal with contour levels in an image.

# ------------------------------------------------------------
# Build a set of contour levels
# ------------------------------------------------------------

def contour_values(
    data_in = None,
    mask = None,
    allvals = False,
    histeq = False,
    linspace = True,
    logspace = False,
    nlev = 100,
    spacing = None,
    minval = None,
    maxval = None,
    ):
    """
    Generate contour levels from data and a series of switches.
    """

    if data_in != None:
        data = data_in[mask*np.isfinite(data_in)]
    
    if minval == None:
        minval = np.min(data[mask])

    if maxval == None:
        maxval = np.min(data[mask])

    if nlev == None:
        nlev = 100

    if not logspace and not linspace and not histeq and not allvals:
        linspace = True

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Return sorted values as contour levels
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if allvals or (histeq and nlev > len(data)): 
        levels = np.sort(data)
        keep = ((levels > minval) < maxval)
        levels = levels[keep]
        levels = np.unique(levels)
        levels = levels[::-1]
        return levels

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Return histogram-equalized levels
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if histeq:
        levels = np.sort(data)
        keep = ((levels > minval) < maxval)
        levels = levels[keep]
        # ... unclear whether to use unique here - when it matters it
        # could break histogram equalitzation but otherwise it has the
        # potential to break the number of levels. For now put the
        # unique here as the code is cleaner.
        levels = np.unique(levels)
        level_index = np.arange(nlev)/(nlev-1)*levels        
        levels = levels[level_index]
        levels = levels[::-1]
        return levels

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Space levels linearly across the data range
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if linspace:
        if spacing == None:
            spacing = (maxval - minval)/(nlev-1)
        levels = minval + \
            spacing*np.arange(np.floor((maxval-minval)/spacing)+1)
        levels = levels[::-1]
        return levels

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Space levels logarithmically across the data range
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if logspace:            
        logminval = np.log10(minval)
        if not np.isfinite(logminval):
            logminval = np.nanmin(np.log10(data))
        logmaxval = np.log10(maxval)
        # ... really shouldn't need the trap for the max
        
        # read the spacing as a factor
        if spacing == None:
            logspacing = (logmaxval-logminval)/(nlevels-1)
        else:
            logspacing = np.log10(spacing)

        loglevels = logminval + \
            logspacing* \
            np.arange(np.floor((logmaxval-logminval)/logspacing)+1)
        levels = 10**loglevels
        levels = levels[::-1]
        return levels

# ------------------------------------------------------------
# IDENTIFY MERGER MATRIX FROM A SET OF LOCAL MAXIMA
# ------------------------------------------------------------

def merge_given_levels(
    data,
    seeds,
    levels,
    mask=None,
    corners=False,
    verbose = False,
    timer = False
    ):
    
    # Initialization
    
    if timer:
        start=time.time()
        full_start=time.time()
        
    if mask == None:
        mask = np.isfinite(data)

    connect = connectivity(data.ndim, corners=corners)

    merger_matrix = np.zeros((len(seeds[0]), len(seeds[0])))*np.nan

    # Loop over levels

    for level in levels:

        print level
        # Label this level

        labels, ncolors = label(mask*(data >= level) , structure=connect)

        # Get the assignments for the seeds

        seed_labels = labels[seeds]

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
#            print len(shared_seeds)

            # Note the threshold in the matrix                        

            for seed in shared_seeds:

                # this step (weirdly) can take a lot of time

                merger_matrix[seed, shared_seeds] = level
                merger_matrix[shared_seeds, seed] = level

    # Clean up diagonal
    for i in range(len(seeds[0])):
        merger_matrix[i,i] = np.nan

    # Finish

    if timer:
        stop=time.time()
        print "Merger calculations took ", stop-start
        start=time.time()

    return merger_matrix
