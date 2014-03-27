import math
import numpy as np
import os
import scipy.stats

# ------------------------------------------------------------
# COMBINE DATA AND NOISE ESTIMATE INTO A SIGNAL-TO-NOISE CUBE
# ------------------------------------------------------------

def snr_cube(data, 
             noise,
             specaxis=0):
    """
    Return a signal-to-noise ratio cube from data + noise
    estimates. Accepts single-value noise, noise spectrum, noise map
    or noise cube. This is mostly here to handle dimensions. It will
    treat 1-d noise as a spectrum, 2-d noise as a spatial map, and
    otherwise directly divide the two.
    """
    
    # We have a single noise estimate
    if len(noise) == 1:
        snr = data / noise
        return snr

    # The noise and data match dimensions
    if data.ndim == noise.ndim:
        snr = data / noise
        return snr

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Noise has lower dimensionality than data - consider two cases.
    # If the noise is 1d assume that we have a spectrum and divide
    # along the spectral axis.  If the noise is 2d assume that we have
    # a map and divide perpendicular to the spectral axis.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Case of a noise spectrum
    if noise.ndim == 1:
        snr = data.copy()
        for i in np.arange(data.shape[specaxis]):
            if specaxis == 0:
                if snr.ndim == 2:
                    snr[i,:] = data[i,:]/noise[i]
                elif snr.ndim == 3:
                    snr[i,:,:] = data[i,:,:]/noise[i]
            elif specaxis == 1:
                if snr.ndim == 1:
                    snr[:,i] = data[:,i]/noise[i]
                elif snr.ndim == 3:
                    snr[:,i,:] = data[:,i,:]/noise[i]
            elif specaxis == 2:
                # ... only possible with three dimensions
                snr[:,:,i] = data[:,i,:]/noise[i]                    

    # Case of a noise map
    if noise.ndim == 2:
        pass

    # We should not be here
    return []

# ------------------------------------------------------------
# EXTRACT A SINGLE NOISE VALUE
# ------------------------------------------------------------

def noise(
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
        return sigma_rob(data[valid], 
                         iterations=1, 
                         thresh=sig_thresh)

    # Use the M.A.D. only.
    if method == "MAD":
        return mad(data[valid])

    # Default to the standard deviation
    return numpy.std(data[not_noise == False])

# ------------------------------------------------------------
# EXTRACT A NOISE MAP
# ------------------------------------------------------------

def noise_map(data, signal=None, blank=None, method="ROBUST"):
    """
    Calculate a two-dimensional noise estimate for the cube, giving
    options for region to avoid and use of robust statistics.
    """

    # TBD - could add a default estimate and require # of channels for
    # estimation.

    # Initialize a default blanking mask
    if blank==None:
        blank = np.isfinite(data) == False
        
    # Identify valid data to estimate the noise
    if signal==None:
        valid = (blank == False)
    else:
        valid = ((blank == False) * (signal == False))

    # initialize the output
    noise_map = np.zeros((data.shape[0],data.shape[1]))

    # Use a robust estimator (sigma clipping + M.A.D.)
    if method == "ROBUST":

        # work out expected sigma threshold for one outlier
        n_spec = data.shape[2]
        sig_thresh = sig_n_outliers(n_spec, 
                                    n_out=1.0, 
                                    pos_only=True)        
        for i in np.arange(data.shape[0]):
            for j in np.arange(data.shape[1]):        
                noise_map[i,j] = sigma_rob((data[i,j])[valid[i,j]], 
                                           iterations=1, 
                                           thresh=sig_thresh)
        return noise_map

    # Use the median absolute deviation
    if method == "MAD":        
        for i in np.arange(data.shape[0]):
            for j in np.arange(data.shape[1]):
                noise_map[i,j] = \
                    mad((data[i,j])[valid[i,j]])   
        return noise_map

    # Default to standard deviation
    for i in np.arange(data.shape[0]):
        for j in np.arange(data.shape[1]):
            noise_map[i,j] = np.std((data[i,j])[valid[i,j]])
        
    return noise_map

# ------------------------------------------------------------
# DERIVE A NOISE SPECTRUM
# ------------------------------------------------------------
    
def noise_spec(data, signal=None, blank=None, method="ROBUST"):
    """
    Work out the noise along the third dimension, giving options for
    region to avoid and use of robust statistics.
    """

    # TBD - could add a default estimate and require # of channels for
    # estimation.

    # Initialize a default blanking mask
    if blank==None:
        blank = np.isfinite(data) == False
        
    # Identify valid data to estimate the noise
    if signal==None:
        valid = (blank == False)
    else:
        valid = ((blank == False) * (signal == False))

    # Initialize output
    noise_spec = np.zeros((data.shape[2]))

    # Estimate via outlier rejection
    if method=="ROBUST":
        # work out expected sigma threshold for one outlier
        n_spec = data.shape[0]*data.shape[1]
        sig_thresh = sig_n_outliers(n_spec, 
                                    n_out=1.0, 
                                    pos_only=True)
        
        for i in range(data.shape[2]):
            noise_spec[i] = sigma_rob((data[:,:,i])[valid[:,:,i]], 
                                      iterations=1, 
                                      thresh=sig_thresh)
        return noise_spec

    # Estimate via median absolute deviation
    if method=="MAD":
        for i in range(data.shape[2]):
            noise_spec[i] = \
                mad((data[:,:,i])[valid[:,:,i]])
        return noise_spec

    # Default to standard deviation
    for i in range(data.shape[2]):
        noise_spec[i] = np.std((data[:,:,i])[valid[:,:,i]])

    return noise_spec

# ------------------------------------------------------------
# DERIVE A NOISE CUBE
# ------------------------------------------------------------

def noise_cube(data, signal=None, blank=None, method="ROBUST",
               smooth_map=None, smooth_spec=None):
    """
    Work out the pixel-by-pixel noise estimates, treating the sky and
    spectral dimensions as separable.
    """

    from scipy.ndimage import median_filter
    
    # First work out and take out the 2-d structure
    nmap = noise_map(data, signal=signal, blank=blank, method=method)
    if smooth_map != None:
        nmap = median_filter(nmap, size=(smooth_map,smooth_map))
    ndata = snr_cube(data, nmap)

    nspec = noise_spec(ndata, signal=signal, blank=blank, method=method)
    if smooth_spec != None:
        nspec = median_filter(nspec, size=smooth_spec)
    
    noise_cube = data.copy()
    for i in range(noise_cube.shape[2]):
        noise_cube[:,:,i,0] = nmap * nspec[i]
        
    return noise_cube

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
# TBD
# ------------------------------------------------------------

# Commentary: in theory the masked array class inside of numpy should
# expedite handling of blanked data (similarly the
# scipy.stats.nanmedian or nanstd functions). However, the masked
# array median operator seems to be either broken or infeasibly
# slow. This forces us into loops, which (shockingly) work out to be
# the fastest of the ways I have tried, but are still far from good.
