import numpy as np
from astropy.io import ascii
from astropy.io import fits
#import matplotlib.pyplot as plt
import catalog_match_sky_trans as cmts
from pyraf import iraf
from astropy.io.votable import parse   
import os


'''Run with
import preimage_coord_transform_sky_SpARCS1034 as pcts

#this matches the input and reference catalog and runs geomap and
#geoxytran on these catalogspct.cat_match_SpARCS1034(septol, fullcat_trans = 1)
pcts.cat_match_sky_SpARCS1034(3.0, fullcat_trans = 1)


'''

#main program#
def cat_match_sky_SpARCS1034(septol, **kwargs):

    '''Written by Gregory Rudnick 10 January 2018

    PURPOSE:

    Read in two catalogs and match them within some tolerance.  Write
    the set of matched coordinates out in a format suitable for geomap input

    Plot the differences in each coordinate.

    INPUT PARAMETERS:

    septol: the maximum separation allowed for a match in arcseconds

    OPTIONAL KEYWORD PARAMETERS

    fullcat_trans:  if set ==1 then transform full catalog

    speczcat_trans:  if set ==1 then transform the spectroscopic catalog

    '''

    #read in the z-band GOGREEN catalog and the reference catalog.
    clustname = 'SpARCS1034'
    (gg_dat, ref_dat, gg_catname, ref_catname, speczcatname) = cat_read_SpARCS1034()
    #print(gg_dat)

    #rename inputs to make code more readable
    raref = np.array(ref_dat['RA_ICRS'])
    decref = np.array(ref_dat['DE_ICRS'])
    rain = gg_dat['RA']
    decin = gg_dat['DEC']

    #match catalogs against each other
    #return matched values
    #catpath = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS1034'
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords_sky.' + clustname + '.in'
    (rarefm,decrefm,rainm,decinm, lims) = cmts.cat_sky_match(raref, decref, rain, decin, septol, matchfile = geomap_infile)

    print("limits of coordinates are",lims)

    #run geomap to compute the transformation and geoxytran to
    #transform the input coordinates using that solution.  This also
    #plots the transformed coordinates.
    (dbfile, geomap_infile) =  georun_sky(clustname, lims)

    #I'm putting the import statement here as it seems I need to
    #import pyraf and run geomap before I run pyplot routines,
    #otherwise the code crashes.
    #import matplotlib.pyplot as plt

    
    #make a plot of the residuals
    pretrans_plotfile = catpath + '/pretrans.' + clustname + '_coordiff_sky.pdf'
    cmts.match_diff_sky_plot(rarefm,decrefm,rainm,decinm, plotfile = pretrans_plotfile)

    if 'fullcat_trans' in kwargs.keys():
        if kwargs['fullcat_trans'] == 1:
            cat_trans_sky(gg_catname, dbfile, geomap_infile, ref_catname, clustname, septol, \
                      ramin = lims['ramin'], ramax = lims['ramax'], decmin = lims['decmin'], \
                      decmax = lims['decmax'])

    if 'speczcat_trans' in kwargs.keys():
        if kwargs['speczcat_trans'] == 1:
            speczcat_trans_sky(speczcatname, dbfile, geomap_infile, ref_catname, clustname, septol, \
                      ramin = lims['ramin'], ramax = lims['ramax'], decmin = lims['decmin'], \
                      decmax = lims['decmax'])
            
def georun_sky(clustname, lims):

    #run geomap and geoxytran.  I have to run this within this module
    #as geomap crashes when run from an external module.

    #in this version 
    
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords_sky.' + clustname + '.in'
    dbfile = catpath + '/' + clustname + '_geomap_sky.db'
    iraf.geomap.xmin=lims['ramin']
    iraf.geomap.xmax=lims['ramax']
    iraf.geomap.ymin=lims['decmin']
    iraf.geomap.ymax=lims['decmax']
    iraf.geomap.maxiter=3
    iraf.geomap.reject=3.0
    iraf.geomap.fitgeometry='rxyscale'
       
    iraf.geomap(geomap_infile, dbfile)
    
    geotran_outfile = catpath + '/geomap_coords_sky.' + clustname + '.out'

    #remove geotran output file if it already exists
    if os.path.isfile(geotran_outfile) is True:
        cmdstr = 'rm ' + geotran_outfile
        os.system(cmdstr)

    #transform the input coordinates in the geomap input file.  
    iraf.geoxytran(geomap_infile, geotran_outfile, dbfile, geomap_infile, direction="backward",xcolumn=3,ycolumn=4)
    
    #read in the geomap input catalog that was transformed using geoxytran.
    trans_ggdat = ascii.read(geotran_outfile)
    ratrans = trans_ggdat['rain']
    dectrans = trans_ggdat['decin']
    raref = trans_ggdat['raref']
    decref = trans_ggdat['decref']

    #make a plot of the new residuals
    posttrans_plotfile = catpath + '/posttrans.' + clustname + '_coordiff_sky.pdf'
    cmts.match_diff_sky_plot(raref,decref,ratrans,dectrans, plotfile = posttrans_plotfile)

    return dbfile, geomap_infile

#def cat_trans_sky_run():
    
#    cat_trans_sky('/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS1034/SPARCS0219_phot_v2.0_USE.fits', 'SpARCS1034_geomap.db', './geomap_coords.SpARCS1034.in','./SpARCS1034_J1.v0.cat', 'SpARCS1034', 3.0)

def cat_trans_sky(incat, dbfile, geomap_infile, refcat, clustname, septol, **kwargs):

    '''This routine reads in a FITS catalog, writes a temporary output
    ASCII catalog, then transforms this catalog using geoxytran.  It
    then appends these new coordinates as new columns to the original catalog.

    It is usually run by hand after the geomap solution has been derived.

    INPUT

    incat: the input fits catalog

    dbfile: the geomap database file

    geomap_infile: used as the database record

    refcat: the original astrometric reference catalog

    clustname: the name of the cluster

    septol: the maximum separation allowed for a match in arcseconds

    OPTIONAL KEYWORDS

    ramin, ramax, decmin, decmax.  These are the limits over which the
    transform was originally computed.  If these are given then it
    uses those limits to color the points in the ra and decdiff plots.
    If one is given, all must be given.

    '''

    #read in GMOS photometry catalog
    cat_hdul = fits.open(incat)
    cat_dat = cat_hdul[1].data

    #read in original astrometric catalog
    #assume that this is a GAIA VOTable
    votable = parse(refcat)
    table = votable.get_first_table()
    ref_dat = table.array
    #convert masked arrays to normal arrays
    raref = np.array(ref_dat['RA_ICRS'])
    decref = np.array(ref_dat['DE_ICRS'])
    
    #ref_dat = ascii.read(refcat)
    #raref = np.array(ref_dat['RA'])
    #decref = np.array(ref_dat['DEC'])

    #tmp coordinate files for geoxytran
    tmpin = 'tmp_geoxytran_sky_in'
    tmpout = 'tmp_geoxytran_sky_out'

    #remove the existing transformed file 
    newcat = incat.replace('.fits','.trans.fits')

    if os.path.isfile(tmpout) is True:
        cmdstr = 'rm ' + tmpout
        os.system(cmdstr)
        cmdstr = 'rm ' + newcat
        os.system(cmdstr)

    #output temporary ASCII file with coordinates
    fo = open(tmpin, "w")
    fo.write("# ra dec\n")
    for i,val in enumerate(cat_dat['RA']):
        fo.write('{} {}\n'.format(cat_dat['RA'][i],cat_dat['DEC'][i]))
    fo.close()

    iraf.geoxytran(tmpin, tmpout, dbfile, geomap_infile, direction="backward",\
                   xcolumn=1, ycolumn = 2)

    #read in ascii output file
    trans_dat = ascii.read(tmpout)
    ratrans = np.array(trans_dat['ra'])
    dectrans = np.array(trans_dat['dec'])
    
    #replace the old RAs and DECs with new RAs and DECs in place    
    tcat_hdul = cat_hdul
    tcat_hdul[1].data['RA'] = ratrans
    tcat_hdul[1].data['DEC'] = dectrans
    
    #write the new fits file
    tcat_hdul.writeto(newcat)

    #match new catalog against original catalog
    mfile = "allcat_sky_match.txt"
    (rarefm,decrefm,ratransm,dectransm, translims) = cmts.cat_sky_match(raref, decref, \
                                                                       ratrans, dectrans, septol, \
                                                                       matchfile = mfile)

    #make a plot of the residuals
    allcattrans_plotfile = 'allcat_trans.' + clustname + '_coordiff_sky.pdf'
    #passes ra and dec limits if they are defined to find source
    #outside of ra and dec lims.  Assumes that if one keyword is given
    #that all are given
    if 'ramin' in kwargs.keys():
        cmts.match_diff_sky_plot(rarefm,decrefm,ratransm,dectransm, plotfile = allcattrans_plotfile, \
                            ramin = kwargs['ramin'], ramax = kwargs['ramax'], \
                            decmin = kwargs['decmin'], decmax = kwargs['decmax'])
    else:
        cmts.match_diff_sky_plot(rarefm,decrefm,ratransm,dectransm, plotfile = allcattrans_plotfile)

def speczcat_trans_sky(incat, dbfile, geomap_infile, refcat, clustname, septol, **kwargs):

    '''This routine reads in a FITS catalog, writes a temporary output
    ASCII catalog, then transforms this catalog using geoxytran.  It
    then writes a new catalog with  old coordinates replaced by the new ones.

    INPUT

    incat: the input fits catalog

    dbfile: the geomap database file

    geomap_infile: used as the database record

    refcat: the original astrometric reference catalog

    clustname: the name of the cluster

    septol: the maximum separation allowed for a match in arcseconds

    OPTIONAL KEYWORDS

    ramin, ramax, decmin, decmax.  These are the limits over which the
    transform was originally computed.  If these are given then it
    uses those limits to color the points in the ra and decdiff plots.
    If one is given, all must be given.

    '''

    #read in the spectroscopic redshift catalog
    cat_hdul = fits.open(incat)
    cat_dat = cat_hdul["MDF"].data

    #read in original astrometric catalog
    #assume that this is a GAIA VOTable
    votable = parse(refcat)
    table = votable.get_first_table()
    ref_dat = table.array
    #convert masked arrays to normal arrays
    raref = np.array(ref_dat['RA_ICRS'])
    decref = np.array(ref_dat['DE_ICRS'])

    #ref_dat = ascii.read(refcat)
    #raref = np.array(ref_dat['ALPHA_SKY'])
    #decref = np.array(ref_dat['DELTA_SKY'])

    #tmp coordinate files for geoxytran
    tmpin = 'tmp_geoxytran_sky_in'
    tmpout = 'tmp_geoxytran_sky_out'

    #remove the existing transformed file 
    newcat = incat.replace('.fits','.trans.fits')

    if os.path.isfile(tmpout) is True:
        cmdstr = 'rm ' + tmpout
        os.system(cmdstr)
        cmdstr = 'rm ' + newcat
        os.system(cmdstr)

    
    #output temporary ASCII file with coordinates
    fo = open(tmpin, "w")
    fo.write("# ra dec\n")
    for i,val in enumerate(cat_dat['RA']):
        #convert speczcat RA to degrees when writing out
        fo.write('{} {}\n'.format(cat_dat['RA'][i] * 15.0 ,cat_dat['DEC'][i]))
    fo.close()

    iraf.geoxytran(tmpin, tmpout, dbfile, geomap_infile, direction="backward",\
                   xcolumn=1, ycolumn = 2)

    #read in ascii output file
    trans_dat = ascii.read(tmpout)
    ratrans = np.array(trans_dat['ra'])
    dectrans = np.array(trans_dat['dec'])
    
    #replace the old RAs and DECs with new RAs and DECs in place    
    tcat_hdul = cat_hdul
    #convert RA back to hours for output
    tcat_hdul["MDF"].data['RA'] = ratrans / 15.0
    tcat_hdul["MDF"].data['DEC'] = dectrans 
    
    #write the new fits file
    tcat_hdul.writeto(newcat)

    #I commented this out for the GAIA catalog as the spetroscopic
    #catalog will have no matches with GAIA
    
    # #match new catalog against original catalog
    # mfile = "speczcat_sky_match.txt"
    # (rarefm,decrefm,ratransm,dectransm, translims) = cmts.cat_sky_match(raref, decref, \
    #                                                                    ratrans, dectrans, septol, \
    #                                                                    matchfile = mfile)

    # #make a plot of the residuals
    # allcattrans_plotfile = 'speczcat_trans.' + clustname + '_coordiff_sky.pdf'
    # #passes ra and dec limits if they are defined to find source
    # #outside of ra and dec lims.  Assumes that if one keyword is given
    # #that all are given
    # if 'ramin' in kwargs.keys():
    #     cmts.match_diff_sky_plot(rarefm,decrefm,ratransm,dectransm, plotfile = allcattrans_plotfile, \
    #                         ramin = kwargs['ramin'], ramax = kwargs['ramax'], \
    #                         decmin = kwargs['decmin'], decmax = kwargs['decmax'])
    # else:
    #     cmts.match_diff_sky_plot(rarefm,decrefm,ratransm,dectransm, plotfile = allcattrans_plotfile)

    
def cat_read_SpARCS1034():

    #read in catalogs

    #official GOGREEN imaging catalog
    catgogreen = '/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS1034/SPARCS1034_zband_IRAC_v1_handcheck_cat.fits'

    gg_hdul = fits.open(catgogreen)
    gg_dat = gg_hdul[1].data
    
    #select the subset of data with a z-band detection
    izdet = np.where(gg_dat['zmag'] < 90.)
    gg_dat = gg_dat[izdet]

    #refcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS1034/SpARCS1034_J.v0.sexcat' 
    #ref_dat = ascii.read(refcat)

    #this is a GAIA reference catalog
    refcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/GOGREEN_GAIA/sparcs1034_gaia_votable.xml'
    votable = parse(refcat)
    table = votable.get_first_table()
    ref_dat = table.array

    #read the catalog spectroscopic redshifts
    speczcat = "/Users/grudnick/Work/GOGREEN/Data/Spectroscopy/v0.3/SpARCS1034_final.fits"

    return gg_dat, ref_dat, catgogreen, refcat, speczcat;

