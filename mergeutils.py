import numpy as np
from scipy.ndimage import label
from scipy.ndimage import histogram
import time

# ------------------------------------------------------------
# IDENTIFY MERGER MATRIX FROM A SET OF LOCAL MAXIMA
# ------------------------------------------------------------

    
#   , levels = levels_in $
#   , spacing = spacing $
#   , minlevel = minlevel $
#   , maxlevel = maxlevel $
#   , nlevels = nlevels $
#   , nmin = nmin $

def merge_given_levels(
    data,
    seeds,
    levels,
    mask=None,
    corners=None,
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

    merger_matrix = np.array((len(seeds), len(seeds)))

    # Loop over levels

    for level in levels:

        # Label this level

        labels = label(mask*(data >= level) , connectivity=connect)

        # Get the assignments for the seeds

        seed_labels = labels[seeds]

        # Get the number of discrete assignments

        max_label = np.max(seed_labels)

        if max_label == 0:
            continue

        # Histogram and look for merged cases

        bins = np.arange(0.5,max_label+0.5,1)
        hist_label = np.histogram(seed_labels,bins=bins)
        digi_label = np.digitize(seed_labels,bins=nins)

        multi_ind = (hist_label > 1).nonzero()

        # Loop over merged cases

        for ind in multi_ind:

            shared_seeds = (digi_label == ind).nonzero()

            # Note the threshold in the matrix

            for seed_1 in shared_seeds:
                for seed_2 in shared_seeds:
                    if seed_1 == seed_2:
                        continue
                    merger_matrix[seed_1, seed_2] = level

    # Finish

    if timer:
        stop=time.time()
        print "Merger calculations took ", stop-start
        start=time.time()

    return merger_matrix
