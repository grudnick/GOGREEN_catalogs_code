# GOGREEN_catalogs_code
preimage_coord_transform_SpARCS0035.py is an obsolete file.  It has been replaced by:
preimage_coord_transform_sky_SpARCS0035.py, which does the transformations in sky coordinates explicitly, with a goal of transforming catalogs and not images.

If run in ipython, it is important to run this without the --pylab command.  For some reason the code crashes if matplotlib.pyplot is imported before geomap is run for the first time.  
