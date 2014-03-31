# Example applications.

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ON A FITS FILE
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

from pyprops import data, mask, noise

# Read a new cube object from FITS
mycube = data.Data("ant_COline_5kms_evenbeam.fits")

# Estimate the noise
mycube.estimate_noise(
    method="MAD", show=True, verbose=True, timer=True)

# Make a mask object out of this data.
mymask = mask.Mask(mycube)

# Use the CPROPS thresholding reciped to make a new mask.
mymask.cprops_mask(hithresh=4.,lothresh=2.)

# Plot the mask over the cube
mymask.contour_on_peak(scale=10)

# Try paring on volume
mymask.pare_on_volume(100)

# Plot the mask again
mymask.contour_on_peak(scale=10)

# Link the mask back in to the cube
mycube.set_signal(mymask)

# Re-estimate the noise (will now avoid known signal)
mycube.estimate_noise(show=True, verbose=True, timer=True)

# Write the mask to disk
mymask.to_fits_file("ant_mask.fits",overwrite=True)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ON A CASA FILE
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Now for CASA

from pyprops import data, mask, noise

# Load data from a CASA image
mycube = data.Data("protostar_13co.image")

# Estimate the noise
mycube.estimate_noise(verbose=True, show=True)

# Make a new mask object
mymask = mask.Mask(mycube)

# Generate a CPROPS mask
mymask.cprops_mask(hithresh=3.,lothresh=2)

# Plot the mask
mymask.contour_on_peak(scale=10)

# Write to disk
mymask.to_casa_image(outfile="protostar_13co.mask", overwrite=True)

# Open the viewer. Go ahead and load the image then contour the mask
viewer()
