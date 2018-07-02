# GOGREEN_catalogs_code
preimage_coord_transform_SpARCS0035.py is an obsolete file.  It has been replaced by:
preimage_coord_transform_sky_SpARCS0035.py, which does the transformations in sky coordinates explicitly, with a goal of transforming catalogs and not images.

preimage_coord_transform_im_<clust>.py transforms the images and catalogs in pixel coordinates.  It requires an initial guess for the transformation.

preimage_coord_transform_sky_<clust>.py transforms the catalogs in WCS coordinates.
  
The southern clusters, i.e. SpARCS0035, SpARCS0219, SpARCS0335, SPT2106, SPT0546, SPT0205 are transformed using reference NIR images

The norther clusters, i.e. SpARCS1033, SpARCS1034, SpARCS1051, SpARCS1616, SpARCS1634, SpARCS1638 are transformed using reference coordinates from GAIA.

If run in ipython, it is important to run this without the --pylab command.  For some reason the code crashes if matplotlib.pyplot is imported before geomap is run for the first time.  
