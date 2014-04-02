# Data class, intended to hold data (as opposed to masks or
# assignment) links to masks, noise estimation, etc. Extends the cube
# class.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import numpy as np
import mask
import noise
from pyprops import cube, mask, noise

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CUBE CLASS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Data(cube.Cube):
    """ 
    
    Class to hold data and supporting information for a data cube
    (images and spectra should also work). Initialize it with a
    filename or a data set.
    
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Attributes (beyond Cube)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # ... associated noise
    noise = None

    # ... a boolean array indicating regions with signal
    signal = None

    # ... the unit
    unit = None
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialization
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(self,
                 *args,
                 **kwargs
                 ):
        cube.Cube.__init__(self, *args, **kwargs)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read/write/copy Overwrites
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    def init_from_cube(
        self,
        prev):
        cube.Cube.copy_from_cube(self, prev)
        self.signal = copy.deepcopy(prev.signal)
        self.unit = copy.deepcopy(prev.unit)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle units
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # ... Jy/beam -> K

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Expose data
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    
    def peak_map(
        self,
        axis=None):
        """
        Return a peak value map along the specified axis.
        """
        if axis==None:
            axis=self.spec_axis

        if axis==None:
            return self.data
        else:            
            data = self.data.copy()
            data[self.valid==False] = np.min(self.data[self.valid])
            return np.max(data, axis=axis)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Interaction with other classes
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # Noise
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    def snr(
        self):
        """
        Calculate a signal to noise cube from the data cube and the
        associated noise object.
        """

        # We need an associated noise object
        if self.noise == None:
            print "Need associated noise object to get S/N."
            return None
        
        if self.noise.scale == None:
            print "Noise object has no scale. Estimate the noise first."
            return

        # Scale the cube by the single-value noise
        snr = self.data / self.noise.scale

        # If there is a spectral component to the noise, scale by that
        if self.noise.spec != None:
            if self.spec_axis != None:
                # ... a view with spectral axis last
                snr_speclast = snr.swapaxes(self.spec_axis,-1)
                # ... this will alter "snr" because it's a view
                snr_speclast /= self.noise.spec

        # If there is a spatial component to the noise, also scale by that.
        if self.noise.map != None:
            if self.data.ndim == 2:
                snr /= self.noise.map
            if self.data.ndim == 3:
                # ... a view with spectral axis first
                snr_specfirst = snr.swapaxes(self.spec_axis,0)
                # ... this will alter SNR
                snr_specfirst /= self.noise.map
        
        # Return the signal to noise cube
        return snr
        
    def estimate_noise(
        self,   
        method="ROBUST",
        show=False,
        timer=False,
        verbose=False):
        """
        Estimate the noise (currently in one dimension). Saves the
        results to a noise object linked to the data.
        """
        self.noise = noise.Noise(self)
        self.noise.calc_1d(
            method=method,
            show=show,
            timer=timer,
            verbose=verbose)

    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # Masks
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    def set_signal(
        self,
        val = None):
        """
        Set the signal mask associated with the data. The signal mask
        is understood to be the location where significant
        astronomical signal is likely present.
        """
        if val != None:
            self.signal = val

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTERNAL FUNCTIONS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Janskies to Kelvin converted (and vice versa)

