import numpy as np
from astropy.io import ascii
from astropy.io import fits
import catalog_match_im_trans as cmti
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

#make x and y coordinate catalog for z-band image 

#rename all functions in preimage_coord_transform_<clust> to have "sky" in the name

#rename catalog_match_trans to catalog_match_sky_trans and rename all
#routines to have sky in name

#go through new preimage_coord_transform_sky_<clust> and change cmt
#routines to have "sky" in name

#check that code works

#make new modules that transform image called
#preimage_coord_transform_im_<clust>.py and do everything in x and y
#coordinates
#finish modifying cat_trans to write out transformed coordinates to
#SEXtractor file and everything after in module

#make version of catalog_match_sky_trans to catalog_match_im_trans and rename all
#routines to have im in name.  Make sure they work on x and y-coordinates.

#write cat_im_match (version of cat_sky_match) with lims in x and
#y coordinate.  Part of this matches against existing photometric
#catalog in z-band to get an identical set of coordinates.

# pass file name with reference and input coordinates of three
# reference points to cmti.cat_im_match.  

#In that code generate files that can be read in by xyxymatch

#run xyxymatch to generate geomap input files

#run geomap on xyxymatch output files.  Change georun_im to take input
#of xyxyoutput file


#figure out why xyxymatch has so few matches by using region files

#test on x and y catalogs 

#transform z-band image

#**************

#transform spectroscopic catalog

#move generic routines to separate module and keep cluster-specific
#routines as a file with cluster name

'''Run with
import preimage_coord_transform_im_SpARCS0219 as pcti

#this matches the input and reference catalog and runs geomap and
#geoxytran on these catalogs
pcti.cat_match_SpARCS0219(septol, fullcat_trans = 1, im_trans = 1)


'''

def cat_match_im_SpARCS0219(septol, **kwargs):

    '''Written by Gregory Rudnick 10 January 2018

    PURPOSE:

    Read in two catalogs and match them within some tolerance.  Write
    the set of matched coordinates out in a format suitable for geomap input

    Plot the differences in each coordinate.

    INPUT PARAMETERS:

    septol: the maximum separation allowed for a match in pixels (0.16asec/pix for GMOS)

    OPTIONAL KEYWORD PARAMETERS

    fullcat_trans:  if set ==1 then transform full catalog

    im_trans: if set==1 then transform image using geotran

    '''

    #read in the z-band GOGREEN catalog and the reference catalog.
    clustname = 'SpARCS0219'
    (gg_dat, ref_dat, gg_catname, ref_catname, zcat_dat, zcatname, initcoordfile) \
        = cat_read_SpARCS0219()
    #print(gg_dat)

    #rename inputs to make code more readable
    xref = np.array(ref_dat['X_IMAGE'])
    yref = np.array(ref_dat['Y_IMAGE'])
    xin = np.array(zcat_dat['X_IMAGE'])
    yin = np.array(zcat_dat['Y_IMAGE'])

    #match catalogs against each other
    #return matched values
    #catpath = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS0219'
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords_im.' + clustname + '.in'
    (xrefm,yrefm,xinm,yinm, lims) \
        = cmti.cat_im_match(xref, yref, xin, yin, septol, \
                            icfile = initcoordfile, matchfile = geomap_infile)

    print("limits of coordinates are",lims)

    #make a plot of the residuals
    pretrans_plotfile = catpath + '/pretrans.' + clustname + '_coordiff_im.pdf'
    cmti.match_diff_im_plot(xrefm,yrefm,xinm,yinm, plotfile = pretrans_plotfile)

    #run geomap to compute the transformation and geoxytran to
    #transform the input coordinates using that solution.  This also
    #plots the transformed coordinates.
    (dbfile, geomap_infile) =  georun_im(clustname, lims)
    
    if 'fullcat_trans' in kwargs.keys():
        if kwargs['fullcat_trans'] == 1:
            cat_trans_im(zcatname, dbfile, geomap_infile, ref_catname, clustname, septol, \
                         xmin = lims['xmin'], xmax = lims['xmax'], ymin = lims['ymin'], \
                         ymax = lims['ymax'])
            
    #transform image using geotran
    if 'im_trans' in kwargs.keys():
        if kwargs['im_trans'] == 1:
            (impath, refimpath) = preimage_read(clustname)
            imtrans(impath, refimpath, dbfile,geomap_infile,lims)
            
def imtrans(impath, refimpath, dbfile, geomap_infile, lims):

    '''Run geotran on an image

    INPUT:

    impath : the full path to the image to be transformed

    refimpath : the full path to the reference image

    dbfile: the geomap database file
    
    geomap_infile: the name of the input coordinate file for the
    geomap transformation, which is used as the record identifier for
    the transfor

    lims: a dictionary that defines the minimum and maximum 

    '''

    #find the x and y limits of the reference image
    refim = fits.open(refimpath)
    xmax = refim[0].header['NAXIS1']
    ymax = refim[0].header['NAXIS2']
    xmin = 1
    ymin = 1
    refim.close()

    outpath = impath.replace('.fits[1]','.trans.fits')

    iraf.geotran(impath, outpath, dbfile, geomap_infile, xmin = xmin, xmax = xmax, \
                 ymin = ymin, ymax = ymax, xscale = 1.0, yscale = 1.0)
            
def georun_im(clustname, lims):

    #run geomap and geoxytran.  I have to run this within this module
    #as geomap crashes when run from an external module.

    #in this version 
    
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords_im.' + clustname + '.in'
    dbfile = catpath + '/' + clustname + '_geomap_im.db'
    iraf.geomap.xmin=lims['xmin']
    iraf.geomap.xmax=lims['xmax']
    iraf.geomap.ymin=lims['ymin']
    iraf.geomap.ymax=lims['ymax']
    iraf.geomap.maxiter=3
    iraf.geomap.reject=3.0
    iraf.geomap.fitgeometry='rxyscale'
       
    iraf.geomap(geomap_infile, dbfile)
    
    geotran_outfile = catpath + '/geomap_coords_im.' + clustname + '.out'

    #remove geotran output file if it already exists
    if os.path.isfile(geotran_outfile) is True:
        cmdstr = 'rm ' + geotran_outfile
        os.system(cmdstr)

    #transform the input coordinates in the geomap input file.  
    iraf.geoxytran(geomap_infile, geotran_outfile, dbfile, geomap_infile, direction="backward",xcolumn=3,ycolumn=4)
    
    #read in the geomap input catalog that was transformed using geoxytran.
    trans_ggdat = ascii.read(geotran_outfile)
    xtrans = trans_ggdat['xin']
    ytrans = trans_ggdat['yin']
    xref = trans_ggdat['xref']
    yref = trans_ggdat['yref']

    #make a plot of the new residuals
    posttrans_plotfile = catpath + '/posttrans.' + clustname + '_coordiff_im.pdf'
    cmti.match_diff_im_plot(xref,yref,xtrans,ytrans, plotfile = posttrans_plotfile)

    return dbfile, geomap_infile


def cat_trans_im(incat, dbfile, geomap_infile, refcat, clustname, septol,  **kwargs):

    '''This routine reads in a FITS catalog, writes a temporary output
    ASCII catalog, then transforms this catalog using geoxytran.  It
    then appends these new coordinates as new columns to the original catalog.

    It is usually run by hand after the geomap solution has been derived.

    INPUT

    incat: the input SExtractor catalog

    dbfile: the geomap database file

    geomap_infile: used as the database record

    refcat: the original astrometric reference catalog

    clustname: the name of the cluster

    septol: the maximum separation allowed for a match in pixels

    OPTIONAL KEYWORDS

    xmin, xmax, ymin, ymax.  These are the limits over which the
    transform was originally computed.  If these are given then it
    uses those limits to color the points in the ra and decdiff plots.
    If one is given, all must be given.

    '''

    #read in GMOS sextractor photometry catalog with (x,y) coordinates
    cat_dat = ascii.read(incat)

    #read in original astrometric catalog
    ref_dat = ascii.read(refcat)
    xref = np.array(ref_dat['X_IMAGE'])
    yref = np.array(ref_dat['Y_IMAGE'])

    #tmp coordinate files for geoxytran
    tmpin = 'tmp_geoxytran_im_in'
    tmpout = 'tmp_geoxytran_im_out'

    #make a new name for the transformed file
    newcat = incat.replace('.sexcat','.trans.cat')

    if os.path.isfile(tmpout) is True:
        cmdstr = 'rm ' + tmpout
        os.system(cmdstr)
        cmdstr = 'rm ' + newcat
        os.system(cmdstr)

    #output temporary ASCII file with coordinates
    fo = open(tmpin, "w")
    fo.write("# x y\n")
    for i,val in enumerate(cat_dat['X_IMAGE']):
        fo.write('{} {}\n'.format(cat_dat['X_IMAGE'][i],cat_dat['Y_IMAGE'][i]))
    fo.close()

    iraf.geoxytran(tmpin, tmpout, dbfile, geomap_infile, direction="backward",\
                   xcolumn=1, ycolumn = 2)

    #read in ascii output file
    trans_dat = ascii.read(tmpout)
    xtrans = np.array(trans_dat['x'])
    ytrans = np.array(trans_dat['y'])

    #replace the old RAs and DECs with new RAs and DECs in place    
    tcat_dat = cat_dat
    tcat_dat['X_IMAGE'] = xtrans
    tcat_dat['Y_IMAGE'] = ytrans

    #write the new catalog file with the transformed coordinates
    ascii.write(tcat_dat, newcat, format='commented_header')
    
    #match new catalog against original catalog
    mfile = "allcat_im_match.txt"
    (xrefm,yrefm,xtransm,ytransm, translims) = cmti.cat_im_match(xref, yref, \
                                                                 xtrans, ytrans, septol, \
                                                                 matchfile = mfile)

    #make a plot of the residuals
    allcattrans_plotfile = 'allcat_trans.' + clustname + '_coordiff_im.pdf'
    #passes ra and dec limits if they are defined to find source
    #outside of ra and dec lims.  Assumes that if one keyword is givenSpARCS0219_GMOS_z.v0.sexcat'
    #that all are given
    
    if 'xmin' in kwargs.keys():
        cmti.match_diff_im_plot(xrefm,yrefm,xtransm,ytransm, plotfile = allcattrans_plotfile, \
                            xmin = kwargs['xmin'], xmax = kwargs['xmax'], \
                            ymin = kwargs['ymin'], ymax = kwargs['ymax'])
    else:
        cmti.match_diff_im_plot(xrefm,yrefm,xtransm,ytransm, plotfile = allcattrans_plotfile)

    
def cat_read_SpARCS0219():

    #read in catalogs

    #official GOGREEN imaging catalog
    catgogreen = '/Users/grudnick/Work/GOGREEN/Catalogs/Preimaging/SpARCS0219/SPARCS0219_phot_v2.1_USE.fits'

    gg_hdul = fits.open(catgogreen)
    gg_dat = gg_hdul[1].data

    #read in z-band image with x and y coordinates
    zcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS0219/SpARCS0219_GMOS_z.v0.sexcat'
    zcat_dat = ascii.read(zcat)
    
    #select the subset of data with a z-band detection
    izdet = np.where(gg_dat['zmag'] < 90.)
    gg_dat = gg_dat[izdet]

    refcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SpARCS0219/SpARCS0219_J.v0.sexcat'

    ref_dat = ascii.read(refcat)

    #file with coordinates of three objects in the input and reference image
    initcoordfile = 'SpARCS0219_initref.coord.txt'


    return gg_dat, ref_dat, catgogreen, refcat, zcat_dat, zcat, initcoordfile;

def preimage_read(clustname):

    '''
    Read in the correct preimage using the GOGREEN.fits file.

    INPUT PARAMETERS:

    clustname: the name of the cluster that is used in the GOGREEN.fits file

    OUTPUT

    the name and path of the fits image.
    
    '''
    configpath = '/Users/grudnick/Work/GOGREEN/Repo/gogreen/config/GOGREEN.fits'
    
    info_hdul = fits.open(configpath)
    info_dat = info_hdul[1].data

    i = np.where(info_dat['cluster'] == clustname)
    imdir = '/Users/grudnick/Work/GOGREEN/Data/Preimaging/' + clustname + '/GMOS/Z/'
    imroot = info_dat['z_image'][i[0]]
    imroot = imroot[0]
    imname = imroot + '.fits[1]'
    impath = imdir + imname

    #also get the path of the FOURSTAR reference image
    refimpath = '/Users/grudnick/Work/GOGREEN/Data/Imaging/Fourstar/Reduced/Sep2016/SpARCS0219_20160910_J_v02_ipe.fits'
    
    return impath,refimpath
