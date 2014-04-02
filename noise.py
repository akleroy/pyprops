# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Imports
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from cube import *
import time

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# The Noise Object
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Noise:
    """ 
    ...
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Contents
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    data = None
    signal = None
    backup = None
    spec_axis = None
    
    # Noise
    scale = None
    spec = None    
    map = None
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize and infrastructure
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(
        self,
        data
        ):
        """
        Construct a new mask object.
        """
        self.data = data        
        self.spec_axis = data.spec_axis

    def set_data(
        self,
        val=None
        ):
        """
        Link the noise object to a data object.
        """
        if val != None:
            self.data = val

    def set_signal_mask(
        self,
        val=None
        ):
        """
        Link the noise object to a mask object.
        """
        if val != None:
            self.signal = val

    def set_spectral_axis(
        self,
        val=None
        ):
        """
        Set the spectral axis.
        """
        if val != None:
            self.spec_axis = val

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Generate a noise estimate
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    def calc_1d(
        self, 
        method="ROBUST",
        timer=False,
        verbose=False,
        show=False,
        showbins=100,
        showsig=5.,
        showlog=True):
        """
        Calculate a single noise estimate for a data set.
        """

        # .............................................................
        # Time the operation if requested.
        # .............................................................

        if timer:
            start=time.time()

        # .............................................................
        # Fit the 1-d noise distribution
        # .............................................................

        # Identify which values to fit - the data are valid and there
        # is no signal associated with them.

        use = self.data.valid
        if self.signal != None:
            use *= (self.signal.data == False)
        if self.data.signal != None:
            use *= (self.data.signal.data == False)
        
        # Call the external noise fitter
        self.scale = est_noise_1d(
            self.data.data[use],
            method=method)

        # .............................................................
        # Report and/or plot
        # .............................................................

        if verbose:
            print "Fit 1d noise distribution value: ", self.scale

        if show:
            # Plot a histogram from -5 to +5 sigma
            med = np.median(self.data.data[use])
            low = med - showsig*self.scale
            high = med + showsig*self.scale
            bin_vals = np.arange(showbins)*(high-low)/(showbins-1.) + low
            fig = plt.figure()
            ax = fig.add_subplot(111)
            n, bins, patches = ax.hist(
                self.data.data[use], 
                bin_vals,
                facecolor='blue',
                log=showlog)
            x_fid = np.arange(10.*showbins)*(high-low)/(showbins*10.) + low
            y_fid = np.exp(-1.*(x_fid-med)**2/(2.*self.scale**2))* \
                np.max(n)
            ax.plot(x_fid,y_fid,'red',linewidth=4,label='fit RMS + median')
            ax.set_xlabel('Data Value')
            ax.set_ylabel('Counts')
            ax.grid=True
            plt.show()

        # .............................................................
        # Finish timing
        # .............................................................

        if timer:
            stop=time.time()
            print "Fitting the noise (1d) took ", stop-start

        return

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Noise Routines
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# These may be of general use and so are not part of the noise
# class. Instead they can be called piecemeal.

# ------------------------------------------------------------
# NOISE ESTIMATORS
# ------------------------------------------------------------

def est_noise_1d(
    data,
    method="ROBUST"):
    """
    Calculate a single noise estimate for a vector.
    """        

    # Use a robust estimator (sigma clipping + M.A.D.)
    if method == "ROBUST":
        sig_thresh = \
            sig_n_outliers(len(data), 
                           n_out=1.0, 
                           pos_only=True)        
        return sigma_rob(data, 
                         iterations=1, 
                         thresh=sig_thresh)

    # Use the M.A.D. only.
    if method == "MAD":
        return mad(data)

    # Default to the standard deviation
    return numpy.std(data)

# ------------------------------------------------------------
# STASTICS HELPER PROCEDURES
# ------------------------------------------------------------

def mad(data, sigma=True):
    """
    Return the median absolute deviation.
    """
    med = np.median(data)
    mad = np.median(np.abs(data - med))
    if sigma==False:
        return mad
    else:
        return mad*1.4826

def sigma_rob(data, iterations=1, thresh=3.0):
    """
    Iterative m.a.d. based sigma with positive outlier rejection.
    """
    noise = mad(data)
    for i in range(iterations):
        ind = (data <= thresh*noise).nonzero()
        noise = mad(data[ind])
    return noise

def sig_n_outliers(n_data, n_out=1.0, pos_only=True):
    """
    Return the sigma needed to expect n (default 1) outliers given
    n_data points.
    """
    perc = float(n_out)/float(n_data)
    if pos_only == False:
        perc *= 2.0
    return abs(scipy.stats.norm.ppf(perc))

# ------------------------------------------------------------
# Commentary
# ------------------------------------------------------------

# In theory the masked array class inside of numpy should expedite
# handling of blanked data (similarly the scipy.stats.nanmedian or
# nanstd functions). However, the masked array median operator seems
# to be either broken or infeasibly slow. This forces us into loops,
# which (shockingly) work out to be the fastest of the ways I have
# tried, but are still far from good.
    
