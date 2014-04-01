import numpy as np
from scipy.ndimage import label
from scipy.ndimage import histogram
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

