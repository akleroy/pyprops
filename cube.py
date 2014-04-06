# Cube class, parent for data, mask, and assign.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# ...............................
# Imports that should always work
# ...............................

import numpy as np
import copy

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

class Cube(object):
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
                 data = None,                 
                 spec_axis = None,                 
                 pyorder = True,
                 hdr = None,
                 wcs = None,
                 casa_cs = None,
                 ):

        # If we have another cube object, then copy it
        if isinstance(data, Cube):
            self.init_from_cube(data)

        # If we have a string then read a file
        if type(data) == type("file.fits"):

            if (data[-4:]).upper() == "FITS":
                # ... if it ends in FITS call astropy
                self.from_fits_file(data)
            else:
                # ... else call CASA
                self.from_casa_image(
                    data,
                    transpose=pyorder)

        # If we have an array, save it
        if type(data) == type(np.array([1,1])):
            self.data = data

        # Allow the user to force/supply attributes but note that
        # these may already be set by the file reader. If so, respect
        # that.

        if spec_axis != None and self.spec_axis == None:
            self.spec_axis = spec_axis    

        if hdr != None and self.hdr == None:
            self.hdr == hdr

        if wcs != None and self.astropy_wcs == None:
            self.astropy_wcs == wcs

        if casa_cs != None and self.casa_cs == None:
            self.casa_cs == wcs
        
        if self.pyorder == None:
            self.pyorder = pyorder

        # Some filling-out logic
            
        if self.astropy_wcs == None and self.hdr != None:
            if astropy_ok:
                self.astropy_wcs = wcs.WCS(self.hdr)

        if self.spec_axis == None:
            self.find_spec_axis()

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

        # ... astropy always reads in pyorder
        self.pyorder = True

        # ... note the header and WCS information
        if skiphdr == False:
            self.astropy_wcs = wcs.WCS(hdulist[0].header)
            self.hdr = hdulist[0].header
            # ... figure out the spectral axis
            self.find_spec_axis()

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
        transpose the cube into a 'python' order and drop degenerate
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
            self.valid = np.transpose(self.valid)
            self.pyorder = True
        else:
            self.pyorder = False
            
        # ... read in coordinate system object
        if skipcs == False:
            if dropdeg:
                temp_cs = ia.coordsys()

                stokes = get_casa_axis(
                    temp_cs,
                    wanttype="Stokes",
                    skipdeg=False,
                    )

                if stokes == None:
                    order = np.arange(self.data.ndim)
                else:
                    order = []
                    for ax in np.arange(self.data.ndim+1):
                        if ax == stokes:
                            continue
                        order.append(ax)

                self.casa_cs = ia.coordsys(order)
                
                # This should work, but coordsys.reorder() has a bug
                # on the error checking. JIRA filed. Until then the
                # axes will be reversed from the original.

                #if transpose == True:
                #    new_order = np.arange(self.data.ndim)
                #    new_order = new_order[-1*np.arange(self.data.ndim)-1]
                #    print new_order
                #    self.casa_cs.reorder(new_order)
            
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
        havetemplate = False
        if template != None:
            havetemplate = True

        if havetemplate == False and self.filename != None:
            template = self.filename
            havetemplate = True

        if havetemplate == False:
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
        stokes = get_casa_axis(image_cs,"Stokes",skipdeg=False)

        # ... branch on:
        # (1) whether we are in "python" order
        # (2) whether we have a stokes axis to pad
        # (3) whether the user supplied data
        if stokes == None:
            if self.pyorder == True:
                if data == None:
                    myimage.putchunk(np.transpose(self.data))
                else:
                    myimage.putchunk(np.transpose(data))
            else:
                if data == None:
                    myimage.putchunk(self.data)
                else:
                    myimage.putchunk(data)
        else:
            if self.pyorder == True:
                if data == None:
                    myimage.putchunk(
                        np.expand_dims(np.transpose(self.data), stokes))
                else:
                    myimage.putchunk(
                        np.expand_dims(np.transpose(data), stokes))
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
    # Copy from another cube
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Could use deepcopy() but I don't totally get how to do this and
    # make it jive with the subclasses. We want to be able, e.g., to
    # initialize a mask off of a data object. This is blunt and
    # inelegant but works.

    def init_from_cube(
        self, 
        prev):
        """
        Initialize a new cube from another cube. Copy the data.
        """
        self.data = copy.deepcopy(prev.data)
        self.valid = copy.deepcopy(prev.valid)
        self.pyorder = copy.deepcopy(prev.pyorder)
        self.spec_axis = copy.deepcopy(prev.spec_axis)
        self.filename = copy.deepcopy(prev.filename)
        self.hdr = copy.deepcopy(prev.hdr)
        self.astropy_wcs = copy.deepcopy(prev.astropy_wcs)

        # CASA objects don't play well with deepcopy, use their own
        # copy method:
        if casa_ok:            
            self.casa_cs = prev.casa_cs.copy()
        else:
            self.casa_cs = None

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle coordinates
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    def make_axes(
        self,
        mode=None,
        sky_only=False,
        freq_only=False,
        restfreq=None):
        """
        Make coordinate axes from th WCS or coordsys information.
        """

        # Figure out which kind of data we have.
        if mode == None:
            if self.astropy_wcs != None:
                mode="FITS"
            elif self.casa_cs != None:
                mode="CASA"
            else:
                print "Cannot determine coordinate system mode."
                return

        # Work out the axis mapping between image on disk and array
        # given an arbitrary spectral axis and the possibility of not
        # being in "pyorder."

        shape = self.data.shape        
        ndim = self.data.ndim        
        if ndim == 2:
            x_im = 0
            y_im = 1

        if ndim == 3 and self.spec_axis != None:
            space_axes = []
            for ax in np.arange(3):
                if self.spec_axis == ax:
                    continue
                space_axes.append(ax)
            if self.pyorder:
                x_ar = space_axes[-1]
                y_ar = space_axes[-2]
            else:
                x_ar = space_axes[0]
                y_ar = space_axes[1]

        if self.pyorder:
            x_im = ndim - x_ar - 1
            y_im = ndim - y_ar - 1
            z_ar = self.spec_axis
            z_im = ndim - z_ar - 1
        else:
            x_im = x_ar
            y_im = y_ar
            z_ar = self.spec_axis
            z_im = z_ar        

        # This can use some generalization. right now it makes sense
        # for images aligned with a cardinal direction and still
        # expects a lot out of the user. Some cleanup would not be
        # very hard; e.g., use the radec flag on the pix-to-world
        # conversion to end up with RA and DEC images. The spectral
        # axis needs more help (CASA has the reverse situation), ask
        # around about whether there's good code to handle that
        # elsewhere.

        if mode.upper() == "FITS":

            # Build the spectral axis
            if self.spec_axis != None and sky_only == False:

                # Make an array of indices along the velocity axis
                ind = np.indices((shape[z_ar],))        

                # Make an array of the same length with zeros for each coord.
                pix = np.zeros((shape[z_ar],ndim))

                # Put the spectral coord into the right place
                pix[:,z_im] = ind.flatten()
                
                # The origin refers to the array, right? That's
                # 0-indexed because we are in numpy... but I'm only
                # like 70% on this. When would you use 1?
                world = self.astropy_wcs.wcs_pix2world(pix,0)

                # Save the axis - for now assume it's velocity
                self.vaxis = world[:,z_im].flatten()
                            
            # Build the spatial axes
            if freq_only == False:

                # ... MAKE IMAGES OF INDICES IN THE (X,Y)
                if self.pyorder:
                    ind = np.indices((shape[y_ar],shape[x_ar]))
                    pix = np.zeros((shape[y_ar]*shape[x_ar], ndim))
                    pix[:, x_im] = (ind[1,:,:]).flatten()
                    pix[:, y_im] = (ind[0,:,:]).flatten()
                else:
                    ind = np.indices((shape[x_ar],shape[y_ar]))
                    pix = np.zeros((shape[x_ar]*shape[y_ar],ndim))
                    pix[:, x_im] = (ind[0,:,:]).flatten()
                    pix[:, y_im] = (ind[1,:,:]).flatten()

                world = self.astropy_wcs.wcs_pix2world(pix,1)

                # ... REFORM THE OUTPUT INTO X AND Y IMAGES
                if self.pyorder:
                    self.ximg = world[:, x_im].reshape(shape[y_ar], shape[x_ar])
                    self.yimg = world[:, y_im].reshape(shape[y_ar], shape[x_ar])
                else:
                    self.ximg = world[:, x_im].reshape(shape[x_ar], shape[y_ar])
                    self.yimg = world[:, y_im].reshape(shape[x_ar], shape[y_ar])

        # This needs testing with degenerate axes
        
        if mode.upper() == "CASA":
            
            # Copy the coordinate system (we will edit it)
            cs_copy = self.casa_cs.copy()

            # Set the units to degrees and Hz
            units = []
            for ctype in cs_copy.axiscoordinatetypes():
                if ctype == "Direction":
                    units.append('deg')
                if ctype == "Spectral":
                    units.append('Hz')
                if ctype == "Stokes":
                    units.append('')
            cs_copy.setunits(units)

            if restfreq != None:
                cs_copy.setrestfreq(restfreq)
                self.restfreq=restfreq
            else:
                # Else store the rest frequency in Hz.
                self.restfreq=qa.convert(mycube.casa_cs.restfrequency(), "Hz")["value"]

            # Build the spectral axis
            if self.spec_axis != None and sky_only == False:

                # Make an array of indices along the velocity axis
                ind = np.indices((shape[z_ar],))        

                # Make an array of the same length with zeros for each coord.
                pix = np.zeros((ndim,shape[z_ar]))
                
                # Place the spectral index into the array
                pix[z_im,:] = ind.flatten()

                # Extract the frequency axis
                world = cs_copy.toworldmany(pix)['numeric']
                self.faxis = world[z_im,:].flatten()
                
                # Convert to velocity
                self.vaxis = np.array(cs_copy.frequencytovelocity(self.faxis))

            # Build the spatial axes
            if freq_only == False:

                if ndim == 1:
                    print "Can't yet to one-d spatial coordinates."
                    
                # ... MAKE IMAGES OF INDICES IN THE (X,Y)
                if self.pyorder:
                    ind = np.indices((shape[y_ar],shape[x_ar]))
                    pix = np.zeros((ndim,shape[y_ar]*shape[x_ar]))                        
                    pix[x_im,:] = (ind[1,:,:]).flatten()
                    pix[y_im,:] = (ind[0,:,:]).flatten()
                else:
                    ind = np.indices((shape[x_ar],shape[y_ar]))
                    pix = np.zeros((ndim,shape[x_ar]*shape[y_ar]))                        
                    pix[x_im,:] = (ind[0,:,:]).flatten()
                    pix[y_im,:] = (ind[1,:,:]).flatten()

                world = cs_copy.toworldmany(pix)['numeric']

                # ... REFORM THE OUTPUT INTO X AND Y IMAGES
                if self.pyorder:
                    self.ximg = world[x_im,:].reshape(shape[y_ar], shape[x_ar])
                    self.yimg = world[y_im,:].reshape(shape[y_ar], shape[x_ar])
                else:
                    self.ximg = world[x_im,:].reshape(shape[x_ar], shape[y_ar])
                    self.yimg = world[y_im,:].reshape(shape[x_ar], shape[y_ar])

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
            if self.astropy_wcs != None:
                mode="FITS"
            elif self.casa_cs != None:
                mode="CASA"
            else:
                print "Cannot determine file mode."
                return

        if mode.upper() == "FITS":
            axis_types = self.astropy_wcs.get_axis_types()
            count = 0
            for axis in axis_types:            
                if axis['coordinate_type'] == "spectral":
                    if self.pyorder == True:
                        self.spec_axis = self.data.ndim - count - 1
                        return
                    else:
                        self.spec_axis = count
                        return
                count += 1
            self.spec_axis = None
            return

        if mode.upper() == "CASA":
            # ... figure out the spectral axis
            spec_in_casa = get_casa_axis(self.casa_cs, "Spectral")

            if spec_in_casa == None:
                self.spec_axis = None
                return

            self.spec_axis = spec_in_casa
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
# NUMPY HELPERS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# It's often useful to move between arrays of coordinates and tuples
# (for indexing). I'm sure there's a better way to do this, but for
# now these two helper functions will make things simpler.

def xyzarr_to_tuple(
    xyzarr=None,
    coordaxis=0):
    """
    Convert a coordinate array to a tuple of arrays.
    """
    
    if xyzarr == None:
        return
    
    if coordaxis == 0:
        return tuple(xyzarr[i,:] for i in range(xyzarr.shape[coordaxis]))
    else:
        return tuple(xyzarr[:,i] for i in range(xyzarr.shape[coordaxis]))

def xyztup_to_array(
    xyztup=None,
    coordaxis=0):
    """
    Convert a coordinate tuple to an array.
    """
    
    if xyztup == None:
        return
    
    if coordaxis == 0:
        return np.vstack(xyztup)
    else:
        return np.vstack(xyztup).transpose()
    
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

def axis_num(
    axis="x",
    system="CASA",
    ndim=3,
    spec_axis=0):
    
    pass
        
