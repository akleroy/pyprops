# Cube class, parent for data, mask, and assign.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# ...............................
# Imports that should always work
# ...............................

import numpy as np

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

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CUBE CLASS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

class Cube:
    """ 
    
    Parent class to hold data and supporting information for a data
    cube.
    
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Attributes
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    # Data value. Either None or a numpy array.
    data = None

    # Valid data flag. Either None or a boolean array matched in size
    # to data. Indicates where data has a valid value. False is
    # equivalent to not-a-number. Useful for numpy operations.
    valid = None

    # Flag indicating that the cube is in "python order", meaning
    # zyx. Convention is to hold all data in this order, though this
    # requires some contortion working with CASA.
    pyorder = None

    # Spectral axis index. None indicates no known spectral axis.
    spec_axis = None

    # Filename that the cube came from.
    filename = None

    # Header for the cube (if it was a FITS file)
    hdr = None

    # The astropy-generated wcs object for the cube
    astropy_wcs = None

    # The CASA coordinate system for the cube.
    casa_cs = None

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialization
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    

    def __init__(self,
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

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # FITS input/output
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    def from_fits_file(
        self,
        filename=None,
        skipdata=False,
        skipvalid=False,
        skiphdr=False
        ):
        """
        Read in a cube from a FITS file using astropy.
        """
        
        # Require astropy to be loaded.
        if astropy_ok == False:
            print "Cannot read FITS files without astropy."
            return

        # ... record the filename
        self.filename=filename

        # ... open the file
        hdulist = fits.open(filename)

        # ... read the data - assume first extension
        if skipdata == False:
            self.data = hdulist[0].data

        # ... note where data is valid
        if skipvalid == False:
            self.valid = np.isfinite(self.data)

        # ... note the header and WCS information
        if skiphdr == False:
            self.astropy_wcs = wcs.WCS(hdulist[0].header)
            self.hdr = hdulist[0].header
            # ... figure out the spectral axis
            self.find_spec_axis()

        # ... astropy always reads in pyorder
        self.pyorder = True

    def to_fits_file(self,
                     outfile=None,
                     overwrite=False,
                     data=None,
                     hdr=None,
                     ):
        """
        Write the cube to a FITS file.
        """

        # Require astropy to be loaded.
        if astropy_ok == False:
            print "Cannot write FITS files without astropy."
            return

        # Allow user supplied header but default to stored header
        if hdr == None:
            hdr = self.hdr

        # Write to disk
        if self.pyorder:
            if data == None:
                fits.writeto(
                    outfile, 
                    self.data, 
                    hdr, 
                    clobber=overwrite)
            else:
                fits.writeto(
                    outfile, 
                    data, 
                    hdr, 
                    clobber=overwrite)                
        else:
            # If not in python order for some reason, transpose            
            if data == None:
                fits.writeto(
                    outfile, 
                    np.transpose(self.data), 
                    hdr, 
                    clobber=overwrite)
            else:
                fits.writeto(
                    outfile, 
                    np.transpose(self.data), 
                    hdr, 
                    clobber=overwrite)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # CASA input/output
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Read and write from a CASA image. This has a few
    # complications. First, by default CASA does not return the
    # "python order" and so we either have to transpose the cube on
    # read or have dueling conventions. Second, CASA often has
    # degenerate stokes axes present in unpredictable places (3rd or
    # 4th without a clear expectation). We need to replicate these
    # when writing but don't want them in memory. By default, try to
    # yield the same array in memory that we would get from astropy.

    def from_casa_image(
        self,
        infile=None,
        transpose=True,
        dropdeg=True,
        skipdata=False,
        skipvalid=False,
        skipcs=False
        ):

        """
        Load a cube into memory from a CASA image. By default it will
        tranpose the cube into a 'python' order and drop degenerate
        axes. These options can be suppressed. The object holds the
        coordsys object from the image in memory.
        """

        # We need to be in CASA to use this functionality
        if casa_ok == False:
            print "Cannot read CASA files without CASA loaded."
            return

        # ... record the filename
        self.filename=infile

        # ... use the ia tool to get the file contents
        ia.open(infile)

        # ... read in the data
        if skipdata == False:
            self.data = ia.getchunk(dropdeg=dropdeg)

        # ... CASA stores validity of data as a mask
        if skipvalid == False:
            self.valid = ia.getchunk(
                getmask=True, 
                dropdeg=dropdeg)

        # ... transpose if desired
        if transpose:
            self.data = np.transpose(self.data)
            self.valid = np.tranpose(self.valid)
            self.pyorder = True
        else:
            self.pyorder = False
            
        # ... read in coordinate system object
        if skipcs == False:
            self.casa_cs = ia.coordsys()
            self.find_spec_axis()
            
        # ... close the ia tool
        ia.close()
        
    def to_casa_image(self,
                      outfile=None,
                      overwrite=False,
                      data=None,
                      template=None,
                      paddegaxes=True
                      ):
        """
        Write the data to a casa image file using the original file as
        a template (required for now).
        """

        # ... check that CASA is loaded
        if casa_ok == False:
            print "Cannot write CASA files without CASA loaded."
            return
    
        # ... set the template        
        if template != None:
            havetemplate = True

        if havetemplate == False and self.filename != None:
            template = self.filename
            havetemplate = True

        if havetempalte == False:
            print "Need a valid template."
            return
        
        # ... use CASA to create a new image based on the template

        # In theory there are other ways to create images using arrays or
        # other ia objects. Right now, require that a new image be based
        # off an existing CASA file.

        myimage = ia.newimagefromimage(
            infile=template,
            outfile=outfile,
            overwrite=overwrite)

        # ... place the supplied data into the new image

        # Read the coordsys of the new image to see if we have degenerate
        # axes. If so and if the "paddegaxes" keyword is true then add
        # these degenerate axes before writing. If requested, also
        # transpose the data before writing.

        # ... get the CASA coord sys tool for this image.
        image_cs = myimage.coordsys()

        # ... identify any stokes axis
        stokes = get_casa_axis(mycs,"Stokes",skipdeg=False)

        # ... branch on:
        # (1) whether we are in "python" order
        # (2) whether we have a stokes axis to pad
        # (3) whether the user supplied data
        if stokes == None:
            if self.pyorder == True:
                if data == None:
                    myimage.putchunk(np.tranpose(self.data))
                else:
                    myimage.putchunk(np.tranpose(data))
            else:
                if data == None:
                    myimage.putchunk(self.data)
                else:
                    myimage.putchunk(data)
        else:
            if self.pyorder == True:
                if data == None:
                    myimage.putchunk(
                        np.expand_dims(np.tranpose(self.data), stokes))
                else:
                    myimage.putchunk(
                        np.expand_dims(np.tranpose(data), stokes))
            else:
                if data == None:
                    myimage.putchunk(
                        np.expand_dims(self.data, stokes))
                else:
                    myimage.putchunk(
                        np.expand_dims(self.data, stokes))

        # ... close the CASA ia tool.
        myimage.close()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle coordinates
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # ... generate arrays (sparser than imhead) of coordinates

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Figure out which axis is spectral
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # Allow the user to manually set the spectral axis or look at the
    # CASA coordsys or Astropy FITS information to figure out which
    # axis is spectral.

    def set_spec_axis(
        self,
        val=None
        ):
        if val != None:
            self.spec_axis = val

    def find_spec_axis(
        self,
        mode=None,
        ):
        """
        Identify the spectral axis.
        """

        # Figure out which kind of data we have.
        if mode == None:
            if self.wcs != None:
                mode=="FITS"
            elif self.casa_cs != None:
                mode=="CASA"
            else:
                return

        if mode.upper() == "FITS":
            axis_types = self.wcs.get_axis_types()
            count = 0
            for axis in axis_types:            
                if axis['coordinate_type'] == "spectral":
                    if pyorder == True:
                        self.spec_axis = self.data.ndim - count - 1
                    else:
                        self.spec_axis = count
                count += 1
            self.spec_axis = None
            return

        if mode.upper() == "CASA":
            # ... figure out the spectral axis
            spec_in_casa = get_casa_axis(self.cs, "Spectral")

            if spec_in_casa == None:
                self.spec_axis = None
                return

            if self.pyorder == True:
                # ... in this case we are transposed from the original order            
                self.spec_axis = self.data.ndim - spec_in_casa - 1
                return
            else:
                # ... here we are in the original order
                self.spec_axis = spec_in_casa
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

    # Potentially set blank values (user or header) to invalid or clip
    # below some threshold (often useful).


# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTERNAL CASA FUNCTIONS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Given a coordinate system object, find the first occurence of a
# specified type of axis. This might be generall useful and coordys
# objects should be small, so put it out here.

def get_casa_axis(
    coordsys,    
    wanttype="Spectral",
    skipdeg=True,
    ):
    """
    Get the position the first axis of a given type from a CASA
    coordsys tool. Returns None if there is no axis of this type. Is
    not case sensitive. By default it returns the index SKIPPING any
    Stokes axis. Set skipdeg to false to change this.
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Check inputs
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if casa_ok == False:
        print "Needs CASA."
        return

    valid = ["SPECTRAL", "STOKES", "DIRECTION"]
    if valid.count(wanttype.upper()) != 1:
        return "Need a valid axis type."

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Get coordinate types and loop
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    axis_types = coordsys.axiscoordinatetypes()
    count = 0
    for axis in axis_types: 
        if skipdeg:
            if axis.upper() == "STOKES":
                continue
        if axis.upper() == wanttype.upper():
            return count
        count += 1

    # ... return none if no match is found.
    return None
