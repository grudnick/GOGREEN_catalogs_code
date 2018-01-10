import numpy as np
from astropy.io import ascii
from astropy.io import fits

#routine to read catalogs
#rouine to match catalogs and plot result
#output geomap input file
#routine to run geomap and transform images
#routine to match transformed images
#routine to plot result
#routine to write output

def cat_match_SpARCS0035():

    '''Written by Gregory Rudnick 10 January 2018

    PURPOSE:

    Read in two catalogs and put them in the 
    them within some tolerance.

    Plot the differences in each coordinate.

    INPUT PARAMETERS:

    raref, decref: the reference coordinates in degrees.  Numpy arrays.

    rain, decin: the input coordinates in degrees.  If performing
    coordinate transforms, these would be the ones to be transformed.
    Numpy arrays

    septol: the maximum separation allowed for a match in the native
    units of the catalog

    '''

    (gg_dat, ref_dat) = cat_read_SpARCS0035();
    print(gg_dat)

def cat_read_SpARCS0035():

    #read in catalogs

    catgogreen = '/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS0035/SPARCS0035_phot_v2.0_USE.fits'

    gg_hdul = fits.open(catgogreen)
    gg_dat = gg_hdul[1].data

    refcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS0035/SpARCS0035_J1.v0.sexcat'

    ref_dat = ascii.read(refcat)


    return gg_dat, ref_dat;
