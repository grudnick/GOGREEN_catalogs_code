import numpy as np
from astropy.io import ascii
from astropy.io import fits
import catalog_match_trans as cmt
from pyraf import iraf
import os

#routine to read catalogs - done
#rouine to match catalogs and plot result - done
#output geomap input file - done
#routine to run geomap and transform catalogs  - done
      
#routine to match transformed catalogs - done
#routine to plot result - done
#routine to write transform full photometric catalog. - done

#for some reason I can't show the plots interactively and have to just
#save them.

#match new catalog against original catalog - done
#plot residuals - done

#Modify the geoxytran calls to indicate the column number of the input file
#figure out whether I should be transforming "forward" or "backward"
#make x-y plot file for new transformed catalog to see where the
#poorly transformed galaxies are

#**************

#transform z-band image
#move generic routines to separate module and keep cluster-specific
#routines as a file with cluster name

'''Run with
import preimage_coord_transform_SpARCS0035 as pct

#this matches the input and reference catalog and runs geomap and
#geoxytran on these catalogs
pct.cat_match_SpARCS0035(septol)

#this applies the transformation to the whole catalog
pct.cat_trans_run()
or
pct.cat_trans('/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS0035/SPARCS0035_phot_v2.0_USE.fits', 'SpARCS0035_geomap.db', './geomap_coords.SpARCS0035.in','./SpARCS0035_J1.v0.cat', 'SpARCS0035', 3.0)

'''

def cat_match_SpARCS0035(septol, **kwargs):

    '''Written by Gregory Rudnick 10 January 2018

    PURPOSE:

    Read in two catalogs and match them within some tolerance.  Write
    the set of matched coordinates out in a format suitable for geomap input

    Plot the differences in each coordinate.

    INPUT PARAMETERS:

    septol: the maximum separation allowed for a match in arcseconds

    OPTIONAL KEYWORD PARAMETERS

    fullcat_trans:  if set ==1 then transform full catalog


    '''

    #read in the z-band GOGREEN catalog and the reference catalog.
    clustname = 'SpARCS0035'
    (gg_dat, ref_dat, gg_catname, ref_catname) = cat_read_SpARCS0035()
    #print(gg_dat)

    #rename inputs to make it more readable
    raref = np.array(ref_dat['ALPHA_SKY'])
    decref = np.array(ref_dat['DELTA_SKY'])
    rain = gg_dat['RA']
    decin = gg_dat['DEC']

    #match catalogs against each other
    #return matched values
    #catpath = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS0035'
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords.' + clustname + '.in'
    (rarefm,decrefm,rainm,decinm, lims) = cmt.cat_sky_match(raref, decref, rain, decin, septol, matchfile = geomap_infile)

    print("limits of coordinates are",lims)

    #make a plot of the residuals
    pretrans_plotfile = catpath + '/pretrans.' + clustname + '_coordiff.pdf'
    cmt.match_diff_plot(rarefm,decrefm,rainm,decinm, plotfile = pretrans_plotfile)

    #run geomap to compute the transformation and geoxytran to
    #transform the input coordinates using that solution.  This also
    #plots the transformed coordinates.
    (dbfile, geomap_infile) =  georun(clustname, lims)
    
    if 'fullcat_trans' in kwargs.keys():
        if kwargs['fullcat_trans'] == 1:
            cat_trans(gg_catname, dbfile, geomap_infile, ref_catname, clustname, septol, \
                      ramin = lims['ramin'], ramax = lims['ramax'], decmin = lims['decmin'], \
                      decmax = lims['decmax'])

    
def georun(clustname, lims):

    #run geomap and geoxytran.  I have to run this within this module
    #as geomap crashes when run from an external module.

    #in this version 
    
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords.' + clustname + '.in'
    dbfile = catpath + '/' + clustname + '_geomap.db'
    iraf.geomap.xmin=lims['ramin']
    iraf.geomap.xmax=lims['ramax']
    iraf.geomap.ymin=lims['decmin']
    iraf.geomap.ymax=lims['decmax']
    iraf.geomap.maxiter=3
    iraf.geomap.reject=3.0
    iraf.geomap.fitgeometry='rxyscale'
       
    iraf.geomap(geomap_infile, dbfile)
    
    geotran_outfile = catpath + '/geomap_coords.' + clustname + '.out'

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
    posttrans_plotfile = catpath + '/posttrans.' + clustname + '_coordiff.pdf'
    cmt.match_diff_plot(raref,decref,ratrans,dectrans, plotfile = posttrans_plotfile)

    return dbfile, geomap_infile

def cat_trans_run():
    
    cat_trans('/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS0035/SPARCS0035_phot_v2.0_USE.fits', 'SpARCS0035_geomap.db', './geomap_coords.SpARCS0035.in','./SpARCS0035_J1.v0.cat', 'SpARCS0035', 3.0)

def cat_trans(incat, dbfile, geomap_infile, refcat, clustname, septol, **kwargs):

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
    ref_dat = ascii.read(refcat)
    raref = np.array(ref_dat['ALPHA_SKY'])
    decref = np.array(ref_dat['DELTA_SKY'])

    #tmp coordinate files for geoxytran
    tmpin = 'tmp_geoxytran_in'
    tmpout = 'tmp_geoxytran_out'

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
    mfile = "allcat_match.txt"
    (rarefm,decrefm,ratransm,dectransm, translims) = cmt.cat_sky_match(raref, decref, \
                                                                       ratrans, dectrans, septol, \
                                                                       matchfile = mfile)

    #make a plot of the residuals
    allcattrans_plotfile = 'allcat_trans.' + clustname + '_coordiff.pdf'
    #passes ra and dec limits if they are defined to find source
    #outside of ra and dec lims.  Assumes that if one keyword is given
    #that all are given
    if 'ramin' in kwargs.keys():
        cmt.match_diff_plot(rarefm,decrefm,ratransm,dectransm, plotfile = allcattrans_plotfile, \
                            ramin = kwargs['ramin'], ramax = kwargs['ramax'], \
                            decmin = kwargs['decmin'], decmax = kwargs['decmax'])
    else:
        cmt.match_diff_plot(rarefm,decrefm,ratransm,dectransm, plotfile = allcattrans_plotfile)

    
def cat_read_SpARCS0035():

    #read in catalogs

    catgogreen = '/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS0035/SPARCS0035_phot_v2.0_USE.fits'

    gg_hdul = fits.open(catgogreen)
    gg_dat = gg_hdul[1].data

    #select the subset of data with a z-band detection
    izdet = np.where(gg_dat['zmag'] < 90.)
    gg_dat = gg_dat[izdet]

    refcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS0035/SpARCS0035_J1.v0.sexcat'

    ref_dat = ascii.read(refcat)


    return gg_dat, ref_dat, catgogreen, refcat;

    
