import numpy as np

def cat_match_SpARCS1034(catgogreen, catgclass):

    '''Written by Gregory Rudnick 9 January 2018

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
