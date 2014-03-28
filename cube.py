# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# .............................
# Try to import astropy modules
# .............................

try:
    from astropy.io import fits
    from astropy.wcs import wcs
except ImportError:
    astropy_ok = False
    print "WARNING! astropy import failed. astropy mode is disabled."
else:
    astropy_ok = True
    

# .............................
# Try to import CASA modules
# .............................

try:
    from tasks import *
    from taskinit import *
    import casac
except ImportError:
    casa_ok = False
    print "WARNING! CASA import failed. CASA mode is disabled."
else:
    casa_ok = True

# ...............................
# Imports that should always work
# ...............................

import numpy as np
import mask
import noise

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CUBE CLASS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Cube:
    """ 
    
    Class to hold data and supporting information for a data cube
    (images and spectra should also work). Initialize it with a
    filename or a data set.
    
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Core data
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    data = None
    valid = None
    signal = None
    spec_axis = None
    unit = None
    
    # ... filename (if any)
    filename = None
    filemode = None

    # ... astropy-specific
    wcs = None
    hdr = None

    # ... casa-specific
    cs = None

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize and load data
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def __init__(self,
                 data_in=None,
                 spec_axis=None,
                 ):
        if type(data_in) == type("file.fits"):
            # ... we got a file name
            if (data_in[-4:]).upper() == "FITS":
                # ... it's a fits file
                self.from_fits_file(data_in)
            else:
                # ... treat it as a CASA image?
                self.from_casa_image(data_in)
        elif type(data_in) == type(np.array([1,1])):
            # ... we got an array
            self.data = data_in
            self.spec_axis = spec_axis

    def from_fits_file(self,
                       filename=None
                       ):
        """
        Read in a cube from a FITS file using astropy.
        """

        if astropy_ok == False:
            print "Cannot read FITS files without astropy."
            return

        # ... record the filename
        self.filename=filename
        self.filemode="astropy"

        # ... open the file
        hdulist = fits.open(filename)

        # ... read the data
        self.data = hdulist[0].data

        # ... note where data is valid
        self.valid = np.isfinite(self.data)

        # ... note the header and WCS information
        self.wcs = wcs.WCS(hdulist[0].header)
        self.hdr = hdulist[0].header

        # ... figure out the spectral axis
        self.spec_axis = \
            self.find_spec_axis_fits()

    def to_fits_file(self,
                     outfile=None,
                     overwrite=False):
        """
        Write the cube to a FITS file.
        """

        if astropy_ok == False:
            print "Cannot write FITS files without astropy."
            return

        if self.filemode != "astropy":
            print "Can only/read write FITS in astropy mode."
            return
        if ((outfile == self.filename) or 
            (outfile == None)) and \
            (overwrite == False):
            print "Will not overwrite with overwite=False."
            return        
        fits.writeto(outfile, self.data, self.hdr, clobber=overwrite)

    def from_casa_image(self,
                        infile=None
                        ):
        
        if casa_ok == False:
            print "Cannot read CASA files without CASA."
            return

        # ... record the filename
        self.filename=infile
        self.filemode="casa"

        ia.open(infile)
        self.data = ia.getchunk(dropdeg=True)
        self.valid = ia.getchunk(getmask=True, dropdeg=True)
        self.cs = ia.coordsys()
        ia.close()

        # ... figure out the spectral axis
        self.spec_axis = \
            get_casa_axis(self.cs, "Spectral")        

    def to_casa_image(self,
                      outfile=None,
                      overwrite=False
                      ):
        
        to_casa_image(
            data = self.data,
            template = self.filename,
            outfile = outfile,
            overwrite = overwrite
            )

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle units
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # ... Jy/beam -> K (calculator external?)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Get beam information
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # ... set beam components, round out beam

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle coordinates
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Figure out which axis is spectral
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    def set_spec_axis(
        self,
        val=None
        ):
        if val != None:
            self.spec_axis = val

    def find_spec_axis_fits(
        self
        ):
        if self.filemode == "astropy":
            axis_types = self.wcs.get_axis_types()
            count = 0
            for axis in axis_types:            
                if axis['coordinate_type'] == "spectral":
                    return self.data.ndim - count - 1
                count += 1
        return

    def find_spec_axis_casa(
        self
        ):
        if self.filemode == "casa":
            axis_types = self.cs.axiscoordinatetypes()
            count = 0
            for axis in axis_types: 
                # ... ignore (degenerate) stokes axis
                if axis == "Stokes":
                    continue
                if axis == "Spectral":
                    return count
                count += 1
        return

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manipulate the cube
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def set_invalid_to(
        self,
        val=np.nan):
        """
        Set all invalid data to a specified number (default nan).
        """
        self.data[self.valid == False] = val

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Expose data
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    
    def peak_map(
        self,
        axis=None):
        """
        """
        if axis==None:
            axis=self.spec_axis

        if axis==None:
            return self.data
        else:
            return np.max(self.data, axis=axis)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Access other classes
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # Noise
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    def snr(
        self):
        """
        """
        if self.noise == None:
            return None
        
        if self.noise.scale == None:
            return

        snr = self.data / self.noise.scale

        if self.noise.spec != None:
            if self.spec_axis != None:
                # ... a view with spectral axis last
                snr_speclast = snr.swapaxes(self.spec_axis,-1)
                # ... this will alter SNR
                snr_speclast /= self.noise.spec

        if self.noise.map != None:
            if self.data.ndim == 2:
                snr /= self.noise.map
            if self.data.ndim == 3:
                # ... a view with spectral axis first
                snr_specfirst = snr.swapaxes(self.spec_axis,0)
                # ... this will alter SNR
                snr_specfirst /= self.noise.map
        
        return snr
        
    def estimate_noise(
        self,   
        method="ROBUST",
        show=False,
        timer=False,
        verbose=False):
        """
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
        if val != None:
            self.signal = val

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTERNAL HELPER FUNCTIONS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def get_casa_axis(
    coordsys,    
    desired_type="Spectral",
    skipdeg=True,
    ):
    """
    Get the position the first axis of a given type. Returns None if
    there is not an axis of this type.
    """
    if casa_ok == False:
        print "Needs CASA."
        return

    axis_types = coordsys.axiscoordinatetypes()
    count = 0
    for axis in axis_types: 
        if skipdeg:
            if axis == "Stokes":
                continue
        if axis == desired_type:
            return count
        count += 1
    return None

def to_casa_image(
    data=None,
    outfile=None,
    template=None,
    overwrite=False):
    """
    ...
    """
    
    if casa_ok == False:
        print "Cannot write CASA files without CASA."
        return
    
    if template == None:
        print "Need a valid template."
        return

    if data == None:
        print "Need data."
        return
        
    myimage = ia.newimagefromimage(
        infile=template,
        outfile=outfile,
        overwrite=overwrite)

    mycs = myimage.coordsys()
    stokes = get_casa_axis(mycs,"Stokes",skipdeg=False)
    if stokes == None:
        myimage.putchunk(data)
    else:
        myimage.putchunk(
            np.expand_dims(data, stokes))
        
    myimage.close()
