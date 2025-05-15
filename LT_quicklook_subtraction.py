import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from astropy.io import fits #FITS files handling
import os,sys,requests,glob,scipy,urllib,datetime,warnings
from scipy.signal import convolve as scipy_convolve
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import wcs
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats,sigma_clip,SigmaClip
from astropy.nddata import Cutout2D
from photutils import IRAFStarFinder,CircularAperture,CircularAnnulus,aperture_photometry,SExtractorBackground,Background2D
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
#from BeautifulSoup import BeautifulSoup python2.7
from bs4 import BeautifulSoup
from image_registration import chi2_shift
from functions import *
#from credentials import *
from credentials import *
from plot_formating import *
#from pyzapspec.py import *
import astroscrappy
#from lacosmic import lacosmic
import multiprocessing
import time,resource
from termcolor import colored
import astroquery
from astroquery import sdss
import datetime
from matplotlib.colors import LogNorm

from datetime import date,timedelta
if not sys.warnoptions:
    warnings.simplefilter("ignore")
norm = ImageNormalize(stretch=SqrtStretch())

def writeswarpdefaultconfigfile():
    configs='''

# Default configuration file for SWarp 2.38.0
# EB 2015-04-22
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME        swarpout.fits      # Output filename
WEIGHTOUT_NAME       coadd.weight.fits # Output weight-map filename
 
HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers
 
#------------------------------- Input Weights --------------------------------
 
WEIGHT_TYPE            NONE      # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
RESCALE_WEIGHTS        Y               # Rescale input weights/variances (Y/N)?
WEIGHT_SUFFIX          .resamp.weight.fits             # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)
WEIGHT_THRESH                         # Bad pixel weight-threshold
 
#------------------------------- Co-addition ----------------------------------
 
COMBINE                N               # Combine resampled images (Y/N)?
COMBINE_TYPE           AVERAGE          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR
CLIP_AMPFRAC           0.3             # Fraction of flux variation allowed
                                       # with clipping
CLIP_SIGMA             4.0             # RMS error multiple variation allowed
                                       # with clipping
CLIP_WRITELOG          N               # Write output file with coordinates of
                                       # clipped pixels (Y/N) 
CLIP_LOGNAME           clipped.log     # Name of output file with coordinates
                                       # of clipped pixels
BLANK_BADPIXELS        N               # Set to 0 pixels having a weight of 0
 
#-------------------------------- Astrometry ----------------------------------
 
CELESTIAL_TYPE         NATIVE           # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL             # MANUAL, ALL or MOST
CENTER                00:00:00.0, 00:00:00.0 # Coordinates of the image center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.3037             # Pixel scale
IMAGE_SIZE             0               # Image size (0 = AUTOMATIC)
 
#-------------------------------- Resampling ----------------------------------
 
RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images
 
RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)
 
FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header
 
GAIN_KEYWORD           ARAWGAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           4             # Default gain if no FITS keyword found
                                       # 0 = infinity (all or for each image)
SATLEV_KEYWORD         SATURATE        # FITS keyword for saturation level (ADU)
SATLEV_DEFAULT         40000.0         # Default saturation if no FITS keyword
 
#--------------------------- Background subtraction ---------------------------
 
SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?
                                       # (all or for each image)
 
BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE              128             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)
BACK_FILTTHRESH        0.0             # Threshold above which the background-
                                       # map filter operates
 
#------------------------------ Memory management -----------------------------
 
VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               1047            # Maximum amount of virtual memory (MB)
MEM_MAX                1256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        1256             # RAM dedicated to co-addition(MB)
 
#------------------------------ Miscellaneous ---------------------------------
 
DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              N               # Write XML file (Y/N)?
XML_NAME               swarp_join.xml       # Filename for XML output
XSL_URL                    file:///lustre/projects/p025_swin/iandreoni/swarp/swarp-2.38.0/share/swarp/swarp.xsl
                                       # Filename for XSL style-sheet
VERBOSE_TYPE           QUIET          # QUIET,LOG,NORMAL, or FULL
NNODES                 1               # Number of nodes (for clusters)
NODE_INDEX             0               # Node index (for clusters)
 
NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic
NOPENFILES_MAX         512             # Maximum number of files opened by SWarp
'''
    if not os.path.exists(path+'config_files/config.swarp'):
     pf = open(path+'config_files/config.swarp','w')
     pf.write(configs)
     pf.close()


def writeswarpconfigfile():
    configs='''
# Default configuration file for SWarp 2.38.0
# EB 2015-04-22
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME        swarpout.fits      # Output filename
WEIGHTOUT_NAME       coadd.weight.fits # Output weight-map filename
 
HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers
 
#------------------------------- Input Weights --------------------------------
 
WEIGHT_TYPE            NONE      # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
RESCALE_WEIGHTS        Y               # Rescale input weights/variances (Y/N)?
WEIGHT_SUFFIX          .resamp.weight.fits             # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)
WEIGHT_THRESH                         # Bad pixel weight-threshold
 
#------------------------------- Co-addition ----------------------------------
 
COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           AVERAGE          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR
CLIP_AMPFRAC           0.3             # Fraction of flux variation allowed
                                       # with clipping
CLIP_SIGMA             4.0             # RMS error multiple variation allowed
                                       # with clipping
CLIP_WRITELOG          N               # Write output file with coordinates of
                                       # clipped pixels (Y/N) 
CLIP_LOGNAME           clipped.log     # Name of output file with coordinates
                                       # of clipped pixels
BLANK_BADPIXELS        N               # Set to 0 pixels having a weight of 0
 
#-------------------------------- Astrometry ----------------------------------
 
CELESTIAL_TYPE         NATIVE           # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL             # MANUAL, ALL or MOST
CENTER                00:00:00.0, 00:00:00.0 # Coordinates of the image center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.3037              # Pixel scale
IMAGE_SIZE             0               # Image size (0 = AUTOMATIC)
 
#-------------------------------- Resampling ----------------------------------
 
RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images
 
RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)
 
FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header
 
GAIN_KEYWORD           ARAWGAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           4             # Default gain if no FITS keyword found
                                       # 0 = infinity (all or for each image)
SATLEV_KEYWORD         SATURATE        # FITS keyword for saturation level (ADU)
SATLEV_DEFAULT         40000.0         # Default saturation if no FITS keyword
 
#--------------------------- Background subtraction ---------------------------
 
SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?
                                       # (all or for each image)
 
BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE              128             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)
BACK_FILTTHRESH        0.0             # Threshold above which the background-
                                       # map filter operates
 
#------------------------------ Memory management -----------------------------
 
VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               1047            # Maximum amount of virtual memory (MB)
MEM_MAX                1256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        1256             # RAM dedicated to co-addition(MB)
 
#------------------------------ Miscellaneous ---------------------------------
 
DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              N               # Write XML file (Y/N)?
XML_NAME               '''+path+'''config_files/swarp_out.xml       # Filename for XML output
XSL_URL                file:///lustre/projects/p025_swin/iandreoni/swarp/swarp-2.38.0/share/swarp/swarp.xsl
                                       # Filename for XSL style-sheet
VERBOSE_TYPE           QUIET          # QUIET,LOG,NORMAL, or FULL
NNODES                 1               # Number of nodes (for clusters)
NODE_INDEX             0               # Node index (for clusters)
 
NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic
NOPENFILES_MAX         512             # Maximum number of files opened by SWarp
''' 
    if not os.path.exists(path+'config_files/config_comb.swarp'):
     pf = open(path+'config_files/config_comb.swarp','w')
     pf.write(configs)
     pf.close()

def writesdssswarpconfigfile(path):
    configs='''
IMAGEOUT_NAME         sdss_coadd.fits      # Output filename CHANGEME
WEIGHTOUT_NAME         coadd.weight.fits # Output weight-map filename
HEADER_ONLY            N               # Only a header as an output file (Y/N)?
WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
WEIGHT_SUFFIX          weight.fits     # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CHI2
COMBINE_BUFSIZE        1024            # Buffer size for combine (MB)
CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
CENTER_TYPE            MANUAL          # MANUAL, ALL or MOST
RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images
RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
OVERSAMPLING           0               # Oversampling in each dimension
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header
GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found
SUBTRACT_BACK          N               # Subtraction sky background (Y/N)?
BACK_TYPE              AUTO            # AUTO or MANUAL
BACK_DEFAULT           0.0             # Default background value in MANUAL
BACK_SIZE              128             # Background mesh size (pixels)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
VMEM_DIR               .               # Directory path for swarp files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                2048            # Maximum amount of usable RAM (MB)
DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
WRITE_FILEINFO         Y               # Write information about each input
WRITE_XML              N               # Write XML file (Y/N)?
XML_NAME               '''+path+'''config_files/swarp_out_sdss.xml       # Filename for XML output
VERBOSE_TYPE           NORMAL           # QUIET,NORMAL or FULL
NTHREADS               8               # Number of simultaneous threads for
IMAGE_SIZE             4000,4000     # scale = 0.396127 arcsec/pixel CHANGEME
PIXELSCALE_TYPE        manual          # OPTIONAL
PIXEL_SCALE            0.396127         # OPTIONAL
''' 
    if not os.path.exists(path+'config_files/swarp_sdss.conf'):
     pf = open(path+'config_files/swarp_sdss.conf','w')
     pf.write(configs)
     pf.close()

def writepsfexparfile():
    params = '''X_IMAGE
VIGNET(60,60)
X_IMAGE
Y_IMAGE
FLUX_APER
FLUXERR_APER
FLUX_RADIUS
ELONGATION
FLAGS
SNR_WIN
'''
    pf = open('config_files/default.param','w')
    pf.write(params)
    pf.close()

def prepsexfile(path):
    params = '''# Simple configuration file for SExtractor prior to PSFEx use
# only non-default parameters are present.
# EB 2007-08-01
#

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME     '''+path+'''config_files/prepsfex.cat   # Catalog filename
CATALOG_TYPE     FITS_LDAC      # FITS_LDAC format
PARAMETERS_NAME  '''+path+'''config_files/default.param # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_MINAREA   5              # minimum number of pixels above threshold
DETECT_THRESH    3              # a fairly conservative threshold
ANALYSIS_THRESH  3              # idem

FILTER           Y              # apply filter for detection ("Y" or "N")?
FILTER_NAME      '''+path+'''config_files/sex.conv   # name of the file containing the filter

#-------------------------------- WEIGHTing ----------------------------------
#-------------------------------- FLAGging -----------------------------------
#------------------------------ Photometry -----------------------------------

PHOT_APERTURES   30             # <- put the referrence aperture diameter here
SATUR_LEVEL      30000.0        # <- put the right saturation threshold here
GAIN             1.6            # <- put the detector gain in e-/ADU here

#------------------------- Star/Galaxy Separation ----------------------------
#------------------------------ Background -----------------------------------
#------------------------------ Check Image ----------------------------------
#--------------------- Memory (change with caution!) -------------------------
#------------------------------- ASSOCiation ---------------------------------
#----------------------------- Miscellaneous ---------------------------------'''
    pf = open(path+'config_files/prepsfex.sex','w')
    pf.write(params)
    pf.close()

def psfexfile(path):
    params = '''
# Default configuration file for PSFEx 3.9.0
# EB 2010-10-10
#

#-------------------------------- PSF model ----------------------------------

BASIS_TYPE      NONE            # NONE, PIXEL, GAUSS-LAGUERRE or FILE
BASIS_NUMBER    11              # Basis number or parameter
PSF_SAMPLING    1.0             # Sampling step in pixel units (0.0 = auto)
PSF_ACCURACY    0.001            # Accuracy to expect from PSF "pixel" values
PSF_SIZE        59,59           # Image size of the PSF model
CENTER_KEYS     X_IMAGE,Y_IMAGE # Catalogue parameters for source pre-centering
PHOTFLUX_KEY    FLUX_APER(1)    # Catalogue parameter for photometric norm.
PHOTFLUXERR_KEY FLUXERR_APER(1) # Catalogue parameter for photometric error

PSF_RECENTER    Y               # Allow recentering of PSF-candidates Y/N ?


#----------------------------- PSF variability -----------------------------

PSFVAR_KEYS     X_IMAGE,Y_IMAGE # Catalogue or FITS (preceded by :) params
PSFVAR_GROUPS   1,1             # Group tag for each context key
PSFVAR_DEGREES  1               # Polynom degree for each group
PSFVAR_NSNAP    1               # Number of PSF snapshots per axis


#----------------------------- Sample selection ------------------------------

SAMPLE_AUTOSELECT  Y            # Automatically select the FWHM (Y/N) ?
SAMPLEVAR_TYPE     SEEING       # File-to-file PSF variability: NONE or SEEING
SAMPLE_FWHMRANGE   2.0,50.0     # Allowed FWHM range
SAMPLE_VARIABILITY 0.4          # Allowed FWHM variability (1.0 = 100%)
SAMPLE_MINSN       2           # Minimum S/N for a source to be used
SAMPLE_MAXELLIP    0.4          # Maximum (A-B)/(A+B) for a source to be used


#----------------------- PSF homogeneisation kernel --------------------------

HOMOBASIS_TYPE     GAUSS-LAGUERRE # NONE or GAUSS-LAGUERRE
HOMOBASIS_NUMBER   10           # Kernel basis number or parameter
HOMOBASIS_SCALE    1.0          # GAUSS-LAGUERRE beta parameter
HOMOPSF_PARAMS     2.0, 3.0     # Moffat parameters of the idealised PSF
HOMOKERNEL_DIR                  # Where to write kernels (empty=same as input)
HOMOKERNEL_SUFFIX  .homo.fits   # Filename extension for homogenisation kernels

#----------------------------- Output catalogs -------------------------------

OUTCAT_TYPE        ASCII_HEAD        # NONE, ASCII_HEAD, ASCII, FITS_LDAC
OUTCAT_NAME        '''+path+'''out/psfex_out.cat  # Output catalog filename




#------------------------------- Check-plots ----------------------------------

CHECKPLOT_DEV       NULL         # NULL, XWIN, TK, PS, PSC, XFIG, PNG,
                                # JPEG, AQT, PDF or SVG
CHECKPLOT_TYPE      FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS # or NONE
CHECKPLOT_NAME      fwhm, ellipticity, counts, countfrac, chi2, resi

#------------------------------ Check-Images ---------------------------------

CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL
                                # Check-image types
CHECKIMAGE_NAME '''+path+'''out/chi.fits,'''+path+'''out/proto.fits,'''+path+'''out/samp.fits,'''+path+'''out/resi.fits,'''+path+'''out/snap.fits,'''+path+'''out/moffat.fits,'''+path+'''out/submoffat.fits,'''+path+'''out/subsym.fits
                                # Check-image filenames
CHECKIMAGE_CUBE Y

# CHI (square-root of) chi^2 maps for all input vignettes

#----------------------------- Miscellaneous ---------------------------------

PSF_DIR         '''+path+'''out/ # Where to write PSFs (empty=same as input)
PSF_SUFFIX      .psf            # Filename extension for output PSF filename
VERBOSE_TYPE    QUIET          # can be QUIET,NORMAL,LOG or FULL
WRITE_XML       Y               # Write XML file (Y/N)?
XML_NAME        '''+path+'''out/psfex.xml       # Filename for XML output
NTHREADS        0               # Number of simultaneous threads for
                                # the SMP version of PSFEx
                                # 0 = automatic'''
    pf = open(path+'config_files/psfex_conf.psfex','w')
    pf.write(params)
    pf.close()



t_start = time.time() #Time the subtraction code
#date_start=datetime.datetime.fromtimestamp(t_start).strftime("%A_%B_%d,_%Y_%I:%M:%S")
'''
# record any pipeline failures
if not os.path.exists(path+'night_log'):
    print('Making')
    os.makedirs(path+'night_log')
errors=open(path+'night_log/errors_'+date_start+'.txt','a')
errors.write('Errors on: '+date_start+'\n')

if not os.path.exists(path+'out'):
    os.makedirs(path+'out')
'''
# config files
if not os.path.exists(path+'config_files'):
    os.makedirs(path+'config_files')
if not os.path.exists(path+'config_files/config.swarp'):
    writeswarpdefaultconfigfile()
if not os.path.exists(path+'config_files/config_comb.swarp'):
    writeswarpconfigfile()
if not os.path.exists(path+'config_files/swarp_sdss.conf'):
    writesdssswarpconfigfile(path)
if not os.path.exists(path+'config_files/default.param'):
    writepsfexparfile()
if not os.path.exists(path+'config_files/prepsfex.sex'):
    prepsexfile(path)
if not os.path.exists(path+'config_files/psfex_conf.psfex'):
    psfexfile(path)

#################################
# load in LT images and combine
#################################
header_Arr=[]
if len(sys.argv)==1:
    print('Add image names without fits extension to run eg.:  python LT_subtraction.py filename1 filename2')
    sys.exit(1)

if len(sys.argv)>2:
    name=[path+s + '.fits' for s in sys.argv[1:len(sys.argv)]]
    name_joined=' '.join(name)
    sci_img_hdu=fits.open(name[1])[0]
    sci_filt=sci_img_hdu.header['FILTER1'][5].lower()
    sci_ra=sci_img_hdu.header['CAT-RA']
    sci_dec=sci_img_hdu.header['CAT-DEC']
    sci_obj=sci_img_hdu.header['OBJECT']
    sci_inst = sci_img_hdu.header['INSTRUME']
    sci_utstart = sci_img_hdu.header['UTSTART']
    sci_prop = sci_img_hdu.header['PROPID']
    sci_obj, sep, tail = sci_obj.partition('_')
    sci_airmass = sci_img_hdu.header['AIRMASS']
    #tail should be the request id from marshal triggering if there is one
    
    sci_seeing=sci_img_hdu.header['L1SEESEC']
    sci_est_seeing=sci_img_hdu.header['SCHEDSEE']
    sci_exp_time=sci_img_hdu.header['EXPTIME']*(len(sys.argv)-1)
    sci_c= SkyCoord(sci_ra+' '+sci_dec, unit=(u.hourangle, u.deg))
    sci_ps=sci_img_hdu.header['CCDSCALE']
    sci_jd=np.mean([fits.open(s)[0].header['mjd']+2400000.5 for s in name])

    
    sci_img_name=sci_obj+'_'+sci_filt+sci_img_hdu.header['DATE-OBS'][:-13]+'_'+str(datetime.timedelta(hours=int(sci_img_hdu.header['DATE-OBS'][11:13]), minutes=int(sci_img_hdu.header['DATE-OBS'][14:16]),seconds=float(sci_img_hdu.header['DATE-OBS'][17:21])).seconds)+'.fits'
    
    swarp_command=swarp_path+" "+name_joined+" -c "+path+"config_files/config_comb.swarp -COPY_KEYWORDS DATE-OBS -CENTER '"+sci_ra+" "+sci_dec+"' -SUBTRACT_BACK Y -VERBOSE_TYPE QUIET -IMAGEOUT_NAME "+path+"combined_imgs"+'/'+sci_img_name+" -RESAMPLE Y -RESAMPLE_DIR '"+path+"' -COMBINE Y -IMAGE_SIZE '"+str(2000)+","+str(2000)+"'"

    if not os.path.exists(path+'combined_imgs/'+sci_img_name):
     print('Combining',len(sys.argv)-1,sci_filt,'images...')
     os.system(swarp_command)
     status=os.system(swarp_command)
     print('Combined',len(sys.argv)-1,' ',sci_filt,' images!')
    if os.path.exists(path+'combined_imgs/'+sci_img_name):
     print('Images:',sci_filt,'already combined')
    sci_img_hdu=fits.open(path+'combined_imgs/'+sci_img_name)[0]
    sci_img=sci_img_hdu.data

if len(sys.argv)==2:
    
    name=sys.argv[1]+'.fits'
    sci_img_hdu=fits.open(path+name)[0]
    sci_filt=sci_img_hdu.header['FILTER1'][5].lower()
    print(colored('Single filter:', 'green'),sci_filt)
    sci_ra=sci_img_hdu.header['CAT-RA']
    sci_dec=sci_img_hdu.header['CAT-DEC']
    sci_obj=sci_img_hdu.header['OBJECT']
    sci_inst = sci_img_hdu.header['INSTRUME']
    sci_utstart = sci_img_hdu.header['UTSTART']
    sci_prop = sci_img_hdu.header['PROPID']
    sci_obj, sep, tail = sci_obj.partition('_')
    sci_airmass = sci_img_hdu.header['AIRMASS']

    sci_obj, sep, tail = sci_obj.partition('_') #tail should be the request id from marshal triggering if there is one
    
    sci_seeing=sci_img_hdu.header['L1SEESEC']
    sci_est_seeing=sci_img_hdu.header['SCHEDSEE']
    sci_img_name=sci_obj+'_'+sci_filt+'comb.fits'
    sci_exp_time=sci_img_hdu.header['EXPTIME']
    sci_c= SkyCoord(sci_ra+' '+sci_dec, unit=(u.hourangle, u.deg))
    sci_ps=sci_img_hdu.header['CCDSCALE']
    sci_jd=sci_img_hdu.header['mjd']+2400000.5
    
    sci_img_name=sci_obj+'_'+sci_filt+sci_img_hdu.header['DATE-OBS'][:-13]+'_'+str(datetime.timedelta(hours=int(sci_img_hdu.header['DATE-OBS'][11:13]), minutes=int(sci_img_hdu.header['DATE-OBS'][14:16]), seconds=float(sci_img_hdu.header['DATE-OBS'][17:21])).seconds)+'.fits'
    sci_img=sci_img_hdu.data


print('Seeing (observed,estimated):',"%.2f %.2f" %(sci_seeing,sci_est_seeing))
if sci_seeing>=sci_est_seeing*20 and sci_seeing>5.0:
  print('Poor seeing')
  #errors.write('Seeing error: '+sci_obj+','+sci_filt+' band, Actual seeing:'+str(sci_seeing)+', Estimated seeing:'+str(sci_est_seeing)+' (ie. guiding/focus error or seeing got much worse): '+'\n')
  #errors.close()
  sys.exit(1)

ra_string="%.6f" % round(sci_c.ra.deg, 6)
dec_string="%.6f" % round(sci_c.dec.deg, 6)
ref_width=sci_img_hdu.header['NAXIS2']*sci_ps*1.1/60.
print(ra_string)
#####################################################
# firstly is it in SDSS footprint? exit if not since we cannot calibrate/download reference image
# If it is in the footprint:
#    Download a grid of SDSS images
#    Stitch them together
#

def sdss_query_image(ra_string,dec_string,filt):

 sdss_url='https://dr12.sdss.org'
 url=sdss_url+'/fields/raDec?ra='+str(ra_string)+'&dec='+str(dec_string)
 #html_page = urllib2.urlopen(url) python2.7 version
 html_page = urllib.request.urlopen(url)
 image_link=[]
 soup = BeautifulSoup(html_page)
 for link in soup.findAll('a'):
    if 'frame-'+str(filt) in link.get('href'):
        image_link.append(sdss_url+link.get('href'))
 
 try:
    image_link=image_link[0]
    image_name=image_link.rsplit('/', 1)[-1]
 except IndexError:
    print(Exception('Exiting... Not in the SDSS footprint, no u-band!'))
    #errors.write(sci_obj+','+sci_filt+': Not in the SDSS footprint '+'\n')
    #errors.close()
    sys.exit(1)
 try:
        r=requests.get(image_link)
        r.raise_for_status()
        if os.path.exists(image_name[:-4]):
            print('SDSS image already downloaded')
        if not os.path.exists(image_name[:-4]):
            zname=image_name
            zfile = open(path+'ref_imgs/'+image_name, 'wb')
            zfile.write(r.content)
            zfile.close()
            os.system('bzip2 -d '+path+'ref_imgs/'+image_name)
            print('downloading new SDSS ',str(filt),'-band..',image_name[:-4])
            ref_path=path+'ref_imgs/'+image_name[:-4]
            os.system('rm '+path+'ref_imgs/'+image_name)
 except requests.exceptions.HTTPError as err:
        print('Not in SDSS footprint! Exiting..')
        #errors.write(sci_obj+','+sci_filt+': Not in SDSS footprint '+'\n')
        #errors.close()
        sys.exit(1)
 print('---------------------------')
 return(ref_path)


if sci_filt=='u':
#if sci_filt=='g' or sci_filt=='r' or sci_filt=='i' or sci_filt=='z':
 spacing=0.1
 u_images=[]
 grid_size=3
 image_name=sci_obj+'_ref.fits'
 ref_path=path+'ref_imgs/'+image_name
 if not os.path.exists(ref_path):
    #download a grid of images
    for x in range(-1,2):
     for y in range(-1,2):
      sci_c.ra.deg+(x*spacing),sci_c.dec.deg+(y*spacing)
      sdss_query_image(sci_c.ra.deg+(x*spacing),sci_c.dec.deg+(y*spacing),sci_filt)
      u_images.append(sdss_query_image(sci_c.ra.deg+(x*spacing),sci_c.dec.deg+(y*spacing),sci_filt))

#if sci_filt=='u':
 image_name=sci_obj+'_ref.fits'
 ref_path=path+'ref_imgs/'+image_name
 if not os.path.exists(ref_path):
  u_images = list(dict.fromkeys(u_images))
  with open(path+'ref_imgs/sdss_list.txt', 'w') as f:
    for item in u_images:
        f.write("%s[0]\n" %(item))
        print("%s[0]\n" %(item))
  f.close()

 image_name=sci_obj+'_ref.fits'
 ref_path=path+'ref_imgs/'+image_name
 if not os.path.exists(ref_path):
  swarp_u_command=swarp_path+" @"+path+'ref_imgs/sdss_list.txt'+" -c "+path+"config_files/swarp_sdss.conf -CENTER '"+sci_ra+" "+sci_dec+"'"+" -IMAGE_SIZE '"+'3000'+","+'3000'+"'"+" -IMAGEOUT_NAME "+path+'ref_imgs/'+sci_obj+'_ref.fits'
  print(swarp_u_command)
  os.system(swarp_u_command)



##############################################
#  Background subtraction
#################################
if not os.path.exists(path+'bkg_subtracted_science'):
    os.makedirs(path+'bkg_subtracted_science')

t_bkg_sub_start = time.time()

sig_clip = SigmaClip(sigma=3.)
bkg_estimator = SExtractorBackground(sig_clip)
bkg = Background2D(sci_img, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
sci_bkgsb=sci_img-bkg.background

bkg_update = Background2D(sci_bkgsb, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
median_bkg=bkg_update.background_median

t_bkg_sub_end = time.time()
bkg_sub_removal_time = t_bkg_sub_end - t_bkg_sub_start

print(colored('Background subtraction time:', 'green'),"%.1f" % round(bkg_sub_removal_time, 0),'seconds')
#print(f'SCI IMAGE SIZE after subtraction {np.shape(sci_bkgsb)}')
##############################################
#  Cosmic ray removal
##############################################

t_cosmic_start = time.time()
#based on how long the exposure times are
if sci_exp_time>60.0:
 sci_cr_mask,sci_bkgsb=astroscrappy.detect_cosmics(sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, gain=1.62, readnoise=8.0, satlevel=65536.0, niter=4, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=True)
 
if sci_exp_time<=60.0:
 sci_cr_mask,sci_bkgsb=astroscrappy.detect_cosmics(sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, gain=1.62, readnoise=8.0, satlevel=65536.0, niter=2, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=True)
 
sci_img_name=sci_img_name[:-5]+"bkgsub.fits"
fits.writeto(path+'bkg_subtracted_science/'+sci_img_name, sci_bkgsb, sci_img_hdu.header,overwrite=True)

t_cosmic_end = time.time()
cosmic_removal_time = t_cosmic_end - t_cosmic_start
print(colored('Cosmic Removal time:', 'green'),"%.1f" % round(cosmic_removal_time, 0),'seconds')

#print(len(sci_cr_mask),np.shape(sci_bkgsb))

#################################
# Get panstamps cutout
#################################

if sci_filt=='u':
 image_name=sci_obj+'_ref.fits'
 print(image_name[:-5])
 print(image_name)
 ref_path=path+'ref_imgs/'+image_name

if sci_filt=='g' or sci_filt=='r' or sci_filt=='i' or sci_filt=='z':
 ref_path=path+'ref_imgs/stack_'+str(sci_filt)+'_ra'+ra_string+'_dec'+dec_string+'_arcsec'+str(int(ref_width*60))+'*'
 print(ref_path)
 if len(glob.glob(ref_path))==0:
   print('Downloading reference image from PS...',sci_ra,sci_dec)
   os.system(panstamps_path+' -f --width='+str(ref_width)+' --filters='+sci_filt+' --downloadFolder='+path+'ref_imgs stack '+sci_ra+' '+sci_dec)
   print(panstamps_path+' -f --width='+str(ref_width)+' --filters='+sci_filt+' --downloadFolder='+path+'ref_imgs stack '+sci_ra+' '+sci_dec)
   #os.system('panstamps -f --width='+str(ref_width)+' --filters='+sci_filt+' --downloadFolder='+path+'ref_imgs stack '+sci_ra+' '+sci_dec)
 if len(glob.glob(ref_path))>0:
    print('PS1 reference image already exists')

ref_img_name=glob.glob(ref_path)[0]
ref_img_hdu=fits.open(ref_img_name)[0]
ref_img=ref_img_hdu.data[0:ref_img_hdu.header['NAXIS2'],0:ref_img_hdu.header['NAXIS1']]


cutouts_dict = {}
cutouts_dict['ref'] = fits.open(ref_img_name)[0].data
cutouts_dict['sci'] = fits.open(path+name)[0].data



#############################
# Align images and make them the same size
############################

if not os.path.exists(path+'aligned_images'):
    os.makedirs(path+'aligned_images')
if os.path.exists(path+'aligned_images'):
    os.system('rm '+path+'aligned_images/*.fits')

print(colored('Aligning science with reference image..', 'green'))

wcs_command='python '+path+'align_quick.py'+" "+path+'bkg_subtracted_science/'+sci_img_name+" "+ref_img_name+"  -m relative -r 100"
print(wcs_command)

try:
 os.system(wcs_command)
except Exception as e:
    print("Error with alignment...")

swarp_command=swarp_path+" "+path+'bkg_subtracted_science/'+sci_img_name+" "+ref_img_name+ " -c "+path+"config_files/config.swarp -CENTER '"+sci_ra+" "+sci_dec+"' -SUBTRACT_BACK Y -VERBOSE_TYPE QUIET -RESAMPLE Y -RESAMPLE_DIR '"+path+"aligned_images/' -COMBINE N -IMAGE_SIZE '"+str(image_size)+","+str(image_size)+"'"
status=os.system(swarp_command)

sci_ali_name=glob.glob(path+'aligned_images/'+sci_img_name[:-5]+'.resamp.fits')[0]

if sci_filt=='u':
    print(path+'aligned_images/'+image_name[:-5]+'.resamp.fits')
    print(glob.glob(path+'aligned_images/'+image_name[:-5]+'.resamp.fits'))
    ref_ali_name=glob.glob(path+'aligned_images/'+image_name[:-5]+'.resamp.fits')[0]
    print(ref_ali_name)

if sci_filt=='g' or sci_filt=='r' or sci_filt=='i' or sci_filt=='z':
    ref_ali_name=glob.glob(path+'aligned_images/stack_'+str(sci_filt)+'_ra'+ra_string+'_dec'+dec_string+'_arcsec'+str(int(ref_width*60))+'_skycell*.resamp.fits')[0]




#if len(ref_ali_name)==0:
    #errors.write(sci_obj+','+sci_filt+': Problem finding the aligned reference image '+'\n')
    #errors.close()
    #sys.exit(1)

#################################################
# Convolve images
#################################################
print('Convolve images..')

if not os.path.exists(path+'convolved_sci'):
    os.makedirs(path+'convolved_sci')
    print('Created folder for convolved_sci images')

if not os.path.exists(path+'convolved_ref'):
    os.makedirs(path+'convolved_ref')
    print('Created folder for convolved_ref images')

if not os.path.exists(path+'out'):
    os.makedirs(path+'out')
    print('Created folder for output PSF models')

sci_conv_name=path+"convolved_sci/"+sci_obj+'_'+sci_filt+'_'+sci_img_hdu.header['DATE-OBS'][:-13]+'_'+str(datetime.timedelta(hours=int(sci_img_hdu.header['DATE-OBS'][11:13]), minutes=int(sci_img_hdu.header['DATE-OBS'][14:16]), seconds=float(sci_img_hdu.header['DATE-OBS'][17:21])).seconds)+'sci_convolved.fits'
ref_conv_name=path+"convolved_ref/"+sci_obj+'_'+sci_filt+'_'+sci_img_hdu.header['DATE-OBS'][:-13]+'_'+str(datetime.timedelta(hours=int(sci_img_hdu.header['DATE-OBS'][11:13]), minutes=int(sci_img_hdu.header['DATE-OBS'][14:16]), seconds=float(sci_img_hdu.header['DATE-OBS'][17:21])).seconds)+'ref_convolved.fits'

#################################################
# CONVOLVE REFERENCE WITH PSF OF SCIENCE IMAGE

if os.path.exists(path+"config_files/prepsfex.cat"):
 os.system("rm "+path+"config_files/prepsfex.cat")

print('Convolving the reference with the PSF of the science image..')

# SExtractor command for the science image
sextractor_command=sex_path+" "+sci_ali_name+" -c "+path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+path+"config_files/prepsfex.cat -MAG_ZEROPOINT 25.0"
#print(sextractor_command)
os.system(sextractor_command)

if os.path.exists(path+'out/proto_prepsfex.fits'):
   print('Deleting prexisting PSF model in out folder')
   os.system('rm '+path+'out/proto_prepsfex.fits')

os.system(psfex_path+" "+path+"config_files/prepsfex.cat -c "+path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")

psf_sci_image_name=path+'out/proto_prepsfex.fits'
psf_sci_image = fits.open(psf_sci_image_name)

hdu_psf_model_sci= fits.open(path+'out/prepsfex.psf')
chi_sq_psf=hdu_psf_model_sci[1].header['CHI2']
print(colored('Reduced Chi^2 of science image PSF fit:','green'),"%.1f" % chi_sq_psf)
if chi_sq_psf>3:
 print(colored('Warning: PSF model may not be accurate','green'))

#Convolve the reference image with a Gaussian kernel
kernel_sci = psf_sci_image[0].data[0]

# Read the REFERENCE image and convolve it with the  kernel science
ref_image_aligned=fits.open(ref_ali_name)
ref_conv = scipy_convolve(ref_image_aligned[0].data, kernel_sci, mode='same', method='fft')
ref_conv=np.nan_to_num(ref_conv)

fits.writeto(ref_conv_name, data=ref_conv, header=ref_image_aligned[0].header,overwrite=True)
#################################################
#CONVOLVE SCIENCE WITH PSF OF REFERENCE IMAGE

os.system("rm "+path+"config_files/prepsfex.cat")

print('Convolving the science with the PSF of the reference image..')

#SExtractor command for the sci image
sextractor_command=sex_path+" "+ref_ali_name+" -c "+path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+path+"config_files/prepsfex.cat -MAG_ZEROPOINT 25.0"
os.system(sextractor_command)

#Select point sources
#Read the SExtractor catalog.
#Note that astropy.io has the option to select explicitely SExtractor format for reading files
if os.path.exists(path+'out/proto_prepsfex.fits'):
   print('psf exists.. deleting')
   os.system('rm '+path+'out/proto_prepsfex.fits')

os.system(psfex_path+" "+path+"config_files/prepsfex.cat -c "+path+"config_files/psfex_conf.psfex")

psf_ref_image_name=path+'out/proto_prepsfex.fits'
psf_ref_image = fits.open(psf_ref_image_name)

hdu_psf_model_ref= fits.open(path+'out/prepsfex.psf')
chi_sq_psf_ref=hdu_psf_model_ref[1].header['CHI2']
print(colored('Reduced Chi^2 of science image PSF fit:','green'),"%.1f" % chi_sq_psf_ref)
if chi_sq_psf_ref>3:
 print(colored('Warning: PSF ref model may not be accurate','green'))



#Convolve the reference image with a Gaussian kernel
kernel_ref = psf_ref_image[0].data[0]

# Read the REFERENCE image and convolve it with the Gaussian kernel
sci_image_aligned=fits.open(sci_ali_name)
sci_conv = scipy_convolve(sci_image_aligned[0].data, kernel_ref, mode='same', method='fft')
sci_conv=np.nan_to_num(sci_conv)

#print(f'SCI_IMG before convolution {np.shape(sci_image_aligned[0].data)}')
#print(f'SCI_CONV SHAPE after convolution {np.shape(sci_conv)}')

fits.writeto(sci_conv_name, data=sci_conv, header=sci_image_aligned[0].header,overwrite=True)

sci_conv_hdu,ref_conv_hdu,ref_ali_hdu,sci_ali_hdu=fits.open(sci_conv_name)[0],fits.open(ref_conv_name)[0],fits.open(ref_ali_name)[0],fits.open(sci_ali_name)[0]
sci_conv,ref_conv,ref_ali,sci_ali=sci_conv_hdu.data,ref_conv_hdu.data,ref_ali_hdu.data,sci_ali_hdu.data

############################################################################
ref_width=sci_img_hdu.header['NAXIS2']*sci_ps/60

if sci_filt=='g' or sci_filt=='r' or sci_filt=='i' or sci_filt=='z':
 ref_cat=panstarrs_query(ra_deg=round(sci_c.ra.deg,6),dec_deg=round(sci_c.dec.deg,6), rad_deg=round(ref_width/60.,6))
 stars=np.where((ref_cat['iMeanPSFMag']-ref_cat['iMeanKronMag']<0.05) & (ref_cat[str(sci_filt)+'MeanPSFMagErr']<0.06) & (ref_cat[str(sci_filt)+'MeanPSFMag']<23.5)  & (ref_cat[str(sci_filt)+'MeanPSFMag']!=-999.))[0]
 #print(ref_cat)
 ref_cat=np.array(ref_cat[stars])
 ref_coords_wcs_sky = SkyCoord(ra=ref_cat['raMean']*u.deg, dec=ref_cat['decMean']*u.deg,frame='fk5')
 ref_coords_wcs=np.column_stack((ref_cat['raMean'],ref_cat['decMean']))
 ref_coords_pix=wcs_to_pixels(sci_ali_name,ref_coords_wcs)

if sci_filt=='u':
 
 ref_cat=sdss_query(ra_deg=round(sci_c.ra.deg,6),dec_deg=round(sci_c.dec.deg,6), rad_deg=round((ref_width)*2/60.,6))
 #print(ref_cat)
 #stars=np.where(ref_cat['type']==6)[0] sdss_query already has criteria type==6 if star so not needed here
 #ref_cat=np.array(ref_cat[stars])
 stars=ref_cat
 ref_coords_wcs_sky = SkyCoord(ra=ref_cat['ra']*u.deg, dec=ref_cat['dec']*u.deg,frame='fk5')
 ref_coords_wcs=np.column_stack((ref_cat['ra'],ref_cat['dec']))
 ref_coords_pix=wcs_to_pixels(sci_ali_name,ref_coords_wcs)
 
 
  


print('catalog stars in PS1/SDSS found=',len(stars),',',round(3600*(ref_width/60),6), 'arcsec search radius')
mean, median, std = sigma_clipped_stats(sci_conv, sigma=3.0)

#######################################################
# iraf star finder
#threshold: The absolute image value above which to select sources.
#fwhm: The full-width half-maximum (FWHM) of the 2D circular Gaussian kernel in units of pixels.
#sigma_radius=1.5,minsep_fwhm=2.5,sharplo=0.5,sharphi=2.0,roundlo=0.0,roundhi=0.2, how round to find sources, default is 0.2 ... roundhi=1-(b/a)
iraffind= IRAFStarFinder(threshold=starscale*std,fwhm=3.0,roundhi=0.3)
print('threshold for detecting stars:',starscale*std)

try:
 sources = iraffind(sci_conv[60:len(sci_conv)-60,60:len(sci_conv)-60] - median)
 fwhm=sources['fwhm']

 star_coords_pix=np.column_stack((sources['xcentroid']+60.,sources['ycentroid']+60.))
 star_coords_wcs=load_wcs_from_file(filename=sci_ali_name,coord=star_coords_pix)
 star_coords_wcs_sky = SkyCoord(ra=star_coords_wcs[:,0]*u.deg, dec=star_coords_wcs[:,1]*u.deg,frame='fk5')

 #find crossover between panstarrs ad reference images
 indx, d2d, d3d =star_coords_wcs_sky.match_to_catalog_sky(ref_coords_wcs_sky)
 #where stars match by search rad in arcseconds!
 upd_indx=np.where(d2d<search_rad/3600.*u.deg)[0]

 matched_ref_coords_pix=ref_coords_pix[indx[upd_indx]]
 matched_star_coords_pix=star_coords_pix[upd_indx]
 matched_fwhm=fwhm[upd_indx]
 matched_catalog=ref_cat[indx[upd_indx]]
 print()
 print(f'ref cat',len(ref_cat))
 print(f'matched cat',len(indx),len(upd_indx))
 print()

 if sci_filt=='g' or sci_filt=='r' or sci_filt=='i' or sci_filt=='z':
  matched_catalog_mag=matched_catalog[str(sci_filt)+'MeanPSFMag']
  #print(matched_catalog)

 if sci_filt=='u':
  print(sci_filt,' matching!',len(ref_cat))
  string_band='psfMag_'+str(sci_filt)
  matched_catalog=ref_cat
  matched_catalog_mag=np.asarray(matched_catalog['mag'])
  print(matched_catalog)
  matched_star_coords_pix = wcs_to_pixels(sci_ali_name,matched_catalog[['ra','dec']])
  print(matched_star_coords_pix)
except Exception as e:
  print(e)
  #errors.write(sci_obj+','+sci_filt+': No stars found in LT image '+'\n')
  #errors.close()
  #sys.exit(1)

print("stars detected in the image=",len(star_coords_pix))
print("length of matched catalog=",len(matched_catalog_mag))

if len(matched_catalog_mag)<10:
    print("WARNING: few stars ("+str(len(matched_catalog_mag))+") matched!")

if not os.path.exists(path+'calb_stars'):
    os.makedirs(path+'calb_stars')
#fig, ax = plt.subplots(figsize=(10,10))

#stars found iraf star finder in convolved scince image
CircularAperture(star_coords_pix, r=12.).plot(color='yellow', lw=4, alpha=0.4)
#stars in referecne in coorinates of the science
CircularAperture(ref_coords_pix, r=9.).plot(color='red', lw=4, alpha=0.4)
# matched stars
CircularAperture(matched_star_coords_pix, r=6.).plot(color='blue', lw=4, alpha=0.4)

'''
plt.text(0.05, 0.94, 'matched - blue, science stars - yellow,reference catalog - red', horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
plt.imshow(sci_conv, cmap='Greys',vmin=0,vmax=50, origin='lower', norm=norm)
plt.savefig(path+'calb_stars/'+str(len(matched_star_coords_pix))+'stars_'+sci_obj+'_'+sci_filt+'_'+sci_img_hdu.header['DATE-OBS'][0:13]+'matched_stars.pdf')
'''
if len(matched_catalog_mag)<3:
    #errors.write(sci_obj+','+sci_filt+': Less than 5 matched calibration stars '+'\n')
    #errors.close()
    sys.exit(1)
    print("Exiting: not enough stars to calibrate!")
print('finished the star matching process')

###################################################################
#  Combine science PSF with reference PSF, flatten PSF to 1D and gaussian fit

if not os.path.exists(path+'convolved_psf'):
    os.makedirs(path+'convolved_psf')
    print('created folder for convolved_psf')
if os.path.exists(path+'convolved_psf'):
    os.system('rm '+path+'convolved_psf/*.fits')

comb_psf = scipy_convolve(kernel_sci, kernel_ref, mode='same', method='fft')
comb_psf=comb_psf/np.sum(comb_psf)
#print('Check: psf normalised to 1, psf sum=',str(np.sum(comb_psf)))
hdu_comb_psf= fits.PrimaryHDU(comb_psf)
hdu_comb_psf.writeto(path+"convolved_psf/"+sci_obj+'_'+sci_filt+'_'+sci_img_hdu.header['DATE-OBS'][0:13]+"comb_psf.fits",overwrite=True)




def cutout_psf(data,psf_array,xpos,ypos):
	xcutout,ycutout=np.shape(psf_array)[0],np.shape(psf_array)[1]
	position=(xpos,ypos)
	size = (xcutout, ycutout)
	cutout = Cutout2D(data, position, size)
	return(cutout.data)

def psf_fit(data_cutout,psf_array):
#fit PSF with x or y-shift
# shift cutout psf to be aligned with psf model
# restrict it to move <5 pixels in each direction so the fit doesn't go wild and fit a nearby bright star etc.
 xoff, yoff, exoff, eyoff = chi2_shift(psf_array,data_cutout, 10,return_error=True, upsample_factor='auto')
 if xoff>5.0 or yoff>5.0:
  xoff,yoff=0.0,0.0
 data_cutout_shift=scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
 resize_sci=np.reshape(data_cutout_shift,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1])
 #resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1],1)
 resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1])

 slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)
 return(slope, intercept, r_value*r_value,xoff, yoff)

def psf_fit_noshift(data_cutout,psf_array):
	# forced photometry - fit with no x or y-shift for limits calculations
	resize_sci=np.reshape(data_cutout,np.shape(data_cutout)[0]*np.shape(data_cutout)[1])
	resize_psf=np.reshape(psf_array,np.shape(data_cutout)[0]*np.shape(data_cutout)[1])
	slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)
	return(slope, intercept, r_value*r_value)



###########################################################
# Zero points!
zp_sci,zp_ref,rsq_sci,rsq_ref=[],[],[],[]

for i in range(0,len(matched_star_coords_pix)):
 try:
  cutout_sci=cutout_psf(data=sci_conv,psf_array=comb_psf,xpos=matched_star_coords_pix[i][0],ypos=matched_star_coords_pix[i][1])
  zp_sci.append(2.5*np.log10(psf_fit(data_cutout=cutout_sci,psf_array=comb_psf)[0])+matched_catalog_mag[i])
  rsq_sci.append(psf_fit(data_cutout=cutout_sci,psf_array=comb_psf)[2])

  cutout_ref=cutout_psf(data=ref_conv,psf_array=comb_psf,xpos=matched_star_coords_pix[i][0],ypos=matched_star_coords_pix[i][1])
  zp_ref.append(2.5*np.log10(psf_fit(data_cutout=cutout_ref,psf_array=comb_psf)[0])+matched_catalog_mag[i])  
  rsq_ref.append(psf_fit(data_cutout=cutout_ref,psf_array=comb_psf)[2])



 except (IndexError,NameError) as error:
  print('error in zp')
  zp_sci.append(0)
  zp_ref.append(0)
  rsq_sci.append(0)
  rsq_ref.append(0)

zp_sci,zp_ref=np.array(zp_sci),np.array(zp_ref)
rsq_ref,rsq_sci=np.array(rsq_ref),np.array(rsq_sci)


#print(np.count_nonzero(~np.isnan(zp_ref)),np.count_nonzero(~np.isnan(zp_sci)))

if np.count_nonzero(~np.isnan(zp_ref))<3:
    if np.count_nonzero(~np.isnan(zp_sci))<3:
     print(colored('few stars in the field','green'))
     #errors.write('Few (< 3) field stars to calibrate in science and reference image. \n')
    if np.count_nonzero(~np.isnan(zp_sci))>3:
     print(colored('Something wrong with reference since no zero point stars detected, but detected in science','green'))
     print(colored('Exiting..','red'))
     #errors.write('Few (< 3) field stars to calibrate in reference image, but  (> 3) in scienc image \n')
     #sys.exit()
#print('rsq_ref',rsq_ref)
#print()
#print('rsq_sci',rsq_sci)
#print()
#print('zp ref',zp_ref)
#print()
#print('zp sci', zp_sci)

if sci_filt=='u':
  thresh=0.1
if sci_filt!='u':
  thresh=0.7
arr=(rsq_ref >=thresh) & (rsq_sci >=thresh) & (zp_sci >=0) & (zp_ref >=0) #filter out stars poorly fitted with PSF fit & weird zp measurement
'''
print('Reference')
print(rsq_ref,zp_ref)
print('Science')
print(rsq_sci,zp_sci)
'''
#print('arr',arr)

if all(val ==False for val in arr)==True:
    weird_img = 'science and reference'
    if all(val_sci_rsq <0.7 for val_sci_rsq in rsq_sci)==True:
        weird_img = 'science'
    else:
        weird_img = 'reference'
    print(f'All stars in the {weird_img} images fit poorly or returned anomalous zeropoint measurements')
    sys.exit(1)

rsq_ref,rsq_sci = rsq_ref[arr],rsq_sci[arr]

zp_sci,zp_ref = zp_sci[arr],zp_ref[arr]
zp_sci,zp_ref=sigma_clip(zp_sci,sigma=3,maxiters=4),sigma_clip(zp_ref,sigma=3,maxiters=4)

if not os.path.exists(path+'scaled_subtracted_imgs'):
    os.makedirs(path+'scaled_subtracted_imgs')
    print('Created folder for scaled and subtracted images')

#print(zp_sci)
zp_sci=zp_sci[~zp_sci.mask]
zp_ref=zp_ref[~zp_ref.mask]

print('zp sci=%.3f std=%.3f num_stars=%i'%(np.nanmedian(zp_sci),np.nanstd(zp_sci),len(zp_sci)))
print('zp ref=%.3f std=%.3f num_stars=%i'%(np.nanmedian(zp_ref),np.nanstd(zp_ref),len(zp_ref)))

###########################################################
# Scale and subtract

#this is another check on the std of the zero point of the image! np.nanstd(zp_ref)/np.sqrt(len(zp_ref))
scale_factor=np.power(10,(np.nanmedian(zp_sci)-np.nanmedian(zp_ref))/2.5)
print('Scale factor: %.3f'%scale_factor)
print('Image dimensions: Science %s, Reference: %s'%(np.shape(sci_conv),np.shape(ref_conv)))

sub=(sci_conv)-(scale_factor*ref_conv)
hdu_fits_sub= fits.PrimaryHDU(sub)
hdu_fits_sub.writeto(path+'scaled_subtracted_imgs/'+sci_img_name[:-5]+'_scaled_subtraction.fits',overwrite=True)

#scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
#sub=(sci_conv)-(scale_factor*scipy.ndimage.shift(ref_conv, [-0.5, -0.5], order=3, mode='reflect', cval=0.0, prefilter=True))

coords_sn=wcs_to_pixels(sci_conv_name,np.column_stack((sci_c.ra.deg,sci_c.dec.deg)))


cutouts_dict['sci_sub'] = fits.open(path+'scaled_subtracted_imgs/'+sci_img_name[:-5]+'_scaled_subtraction.fits')[0].data

#[print(os.path.isfile(cutouts_dict[i])) for i in cutouts_dict.keys()]

'''
print(coords_sn)

sn_x=coords_sn[0][0]
sn_y=coords_sn[0][1]

sci_img_mean, sci_img_median, std = sigma_clipped_stats(cutouts_dict['sci'])
vmin=sci_img_median-2*std
vmax=sci_img_median+2*std

plt.close()
plt.figure()
plt.imshow(cutouts_dict['sci'],cmap='gray',vmax=vmax,vmin=vmin,)
plt.ylim(sn_y-30,sn_y+30)
plt.xlim(sn_x-30,sn_x+30)
plt.show()
'''
def mag_err_function(data,psf,zp_sci,num,sn_x,sn_y):
    
    print(colored('Calculating magnitudes', 'green'))
    magerr,maglim,flux_bkg_list,flux_new_sn_list=[],[],[],[]
    
    # cutout a part of the image the same size as the measured PSF
    sn_cutout=cutout_psf(data=data,psf_array=psf,xpos=sn_x,ypos=sn_y)
    # calculate the magnitude of the object
    sn_flux=psf_fit(sn_cutout,psf_array=psf)[0]
    sn_mag=-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)
    # print the offset of the shifted fit
    print('offset of shifted PSF fit (xpos,ypos)=%.3f %.3f'%(psf_fit(sn_cutout,psf_array=psf)[3],psf_fit(sn_cutout,psf_array=psf)[4]))
    # set up a grid for the values
    psf_size=np.shape(psf)[0]+1
    x=np.linspace(-num*psf_size,num*psf_size,(num*2)+1)
    y=x
    x0, y0, radius = 0.0, 0.0, psf_size/2
    x, y = np.meshgrid(x, y)
    r = np.sqrt((x - x0)**2 + (y - y0)**2)
    outside = r > radius
    x=x[outside]
    y=y[outside]
    # plot of the grid used for uncertainties..
    '''
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set(xlabel=r'X', ylabel=r'Y')
    ax.scatter(x+psf_size*(num+1), y+psf_size*(num+1),s=11,color='orange',edgecolors='black')
    ax.scatter(x+psf_size*(num+1), y+psf_size*(num+1),s=1150,color='orange',edgecolors='black',marker='s',alpha=0.2)

    plt.imshow(data[int(sn_x)-(psf_size*(num+1)):int(sn_x)+(psf_size*(num+1)),int(sn_y)-(psf_size*(num+1)):int(sn_y)+(psf_size*(num+1))], cmap='Greys',vmin=0,vmax=50, origin='lower', norm=norm)
    plt.bar(left=0+(psf_size*(num+1)),height=psf_size, width=psf_size,bottom=-(psf_size/2)+(psf_size*(num+1)),alpha=0.2,edgecolor='black',linewidth=2)
    plt.savefig(path+'calb_stars/grid_for_uncertainties.pdf',dpi=400,rasterized=True,bbox_inches='tight')
    '''
    #
    for i in range(0,len(x)):
       bkg_cutout=cutout_psf(data=data,psf_array=psf,xpos=sn_x+x[i],ypos=sn_y+y[i])
       flux_bkg=psf_fit_noshift(data_cutout=bkg_cutout,psf_array=psf)[0]
       #flux_bkg_list is the PSF fitted to the sky
       flux_bkg_list.append(flux_bkg)
       #
       new_sn=bkg_cutout+((psf)*psf_fit(sn_cutout,psf_array=psf)[0])+psf_fit(sn_cutout,psf_array=psf)[1]
       flux_new_sn=psf_fit(data_cutout=new_sn,psf_array=psf)[0]
       #flux_new_sn_list is the PSF fitted to the artificial supernova
       flux_new_sn_list.append(flux_new_sn)

    flux_bkg_list=sigma_clip(flux_bkg_list,sigma=3,maxiters=3)
    flux_new_sn_list=sigma_clip(flux_new_sn_list,sigma=3,maxiters=3)
    
    #print('AFTER SIGMA CLIPPING...2 iterations, 5 sigma')
    print('S/N (std background)= %.3f'%(sn_flux/np.std(flux_bkg_list))) 
    print('S/N (std artifical sn)= %.3f'%(sn_flux/np.std(flux_new_sn_list)))
    #print('S/N (std art sn)=',sn_flux/np.std(flux_new_sn_list))
    
    if sn_flux/np.std(flux_bkg_list)<=2:
        #not detected
        minimal_mag=sn_mag
        mag=99.0
        magstd=99.0
        magerr=99.0
        print('less than 2.5 sigma')
        print('mag=%.3f+/-%.3f ' %(sn_mag,-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.std(flux_new_sn_list))))
        print('flux corresponding to the (marginal) detection + 2*sigma =%.3f'%(-2.5*np.log10(np.median(sn_flux)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci)))
        print('5-sig limit =%.3f'%(-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci)))
        print('3-sig limit =%.3f'%(-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci)))
        # pick the most conservative out of:
        #flux corresponding to the (marginal) detection + 2*sigma')
        #'3-sig limit
        limits=[-2.5*np.log10(np.median(sn_flux)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci),-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci)]
        maglim=np.min(limits)
        print('l')
    print('y',sn_flux/np.std(flux_bkg_list))
    if sn_flux/np.std(flux_bkg_list)>=2:
        minimal_mag=sn_mag
        print('mag=%.3f+/-%.3f ' %(sn_mag,-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.std(flux_new_sn_list))))
        print('flux corresponding to the (marginal) detection + 2*sigma',-2.5*np.log10(np.median(sn_flux)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci))
        print('5-sig limit =',-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci))
        print('3-sig limit =',-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci))
        mag=sn_mag
        magstd=-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.std(flux_new_sn_list))
        magerr=-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.std(flux_new_sn_list))
        maglim=-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci)

    #print 'background: median flux, std flux=',np.median(flux_bkg_list),np.std(flux_bkg_list)
    #print 'background 1-sig =',-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*1))+np.nanmedian(zp_sci)
    #print 'background 3-sig =',-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci)
    #print 'background 5-sig =',-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci)
    
    # * If >3 sigma detection: report to Marshal with 1-sigma errors.
    # * If <3 sigma detection:r eport whichever of these two is more conservative:
    # * The flux corresponding to the (marginal) detection + 2*sigma.
    # * The flux corresponding to 3*sigma.
 
    #print 'artifical sn (sig-clip): median flux, std flux=',np.median(flux_new_sn_list),np.std(flux_new_sn_list)
    #print 'artifical sn median mag=',-2.5*np.log10(np.median(flux_new_sn_list))+np.nanmedian(zp_sci)
    #print 'medium_artifical_sn_mag - (medium_artifical_sn_mag+std) =',-2.5*np.log10(np.median(flux_new_sn_list))+2.5*np.log10(np.median(flux_new_sn_list)+np.std(flux_new_sn_list))
    #print 'sn_mag - (sn_mag+std) =',-2.5*np.log10(psf_fit(sn_cutout,psf_array=psf)[0])+2.5*np.log10(psf_fit(sn_cutout,psf_array=psf)[0]+np.std(flux_new_sn_list))
 
    return(mag,magstd,magerr,maglim,sn_flux/np.std(flux_bkg_list),sn_flux/np.std(flux_new_sn_list),-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*1))+np.nanmedian(zp_sci),-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci),-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci),minimal_mag)


t_photometry_start = time.time()

mag=mag_err_function(data=sub,psf=comb_psf,zp_sci=zp_sci,num=4,sn_x=coords_sn[0][0],sn_y=coords_sn[0][1])
mag_all_err=np.power(mag[1]**2+np.power(np.nanstd(zp_ref),2)+np.power(np.nanstd(zp_sci),2),0.5)
print("(mag_err**2+std_zp**2+std_ref**2)^(0.5) %.3f \n" %(mag_all_err))

print(colored('----------------------------------------', 'blue'))
#print('mag=%.3f+/-%.3f lim=%.3f' %(mag[0],mag[1],mag[3]))
print('mag=%.3f+/-%.3f lim=%.3f' %(mag[0],mag_all_err,mag[3]))
print(colored('----------------------------------------', 'blue'))

t_photometry_end = time.time()

photometry_time = t_photometry_end - t_photometry_start

print(colored('Photometry time:', 'green'),"%.1f" % round(photometry_time, 0),'seconds')

memkb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(colored('Memory usage Mbyte:', 'green')," %5.1f" %(memkb/1024.0/1024.0))
#print(colored('Memory usage kbyte:', 'green'),"%s" %(memkb))
print(f"sdss{sci_filt}")

if not os.path.exists(path+'phot_fits_info'):
    os.mkdir(path+'phot_fits_info')

d = [str(sci_jd-2400000.5),sci_utstart,sci_ra,sci_dec,sci_prop,sci_inst,sci_exp_time,sci_airmass,sci_seeing,sci_est_seeing]
t=open(path+'phot_fits_info/'+sci_img_name[:-11]+"_phot_info.txt","w")
t.write(f'{d[0]},{d[1]},{d[2]},{d[3]},{d[4]},{d[5]},{d[6]},{d[7]},{d[8]},{d[9]}')
t.close()

mjd = sci_jd-2400000.5
#header_df = pd.DataFrame(d,columns=['utcstart','ra','dec','propid','inst','expt','airmass','seeing','est_seeing'])
#print(header_df)
#header_df.to_csv(path+'photo_fits_info/'+sci_img_name[:-11]+"_phot_info.txt")
t_end = time.time()
if not os.path.exists(path+'photometry'):
    os.makedirs(path+'photometry')
    print('created the photometry folder')
if os.path.exists(path+'photometry'):
    text_file = open(path+'photometry/'+sci_img_name[:-11]+"_photometry.txt", "w")
    text_file.write(sci_obj+" "+f"sdss{sci_filt}"+" "+str(mjd)+" %.3f %.3f %.3f %s %s %s \n" % (mag[0],mag_all_err,mag[3],ra_string,dec_string,sci_exp_time))
    text_file.write("S/N %.3f, %.3f sigma \n" % (mag[4],mag[5]))
    text_file.write("1 sigma limit: %.3f\n" % (mag[6]))
    text_file.write("3 sigma limit: %.3f\n" % (mag[7]))
    text_file.write("5 sigma limit: %.3f\n" % (mag[8]))
    text_file.write("mag=%.3f\n" %(mag[9]))
    text_file.write("time photometry: %.1f seconds\n" % (t_end-t_start))
    text_file.write("memory useage: %.1f Mbyte\n" % (memkb/1024.0/1024.0))
    text_file.write("zeropoint reference %.3f sd %.3f \n" % (np.median(zp_ref),np.nanstd(zp_ref)))
    text_file.write("zeropoint science %.3f sd %.3f \n" % (np.median(zp_sci),np.nanstd(zp_sci)))
    text_file.write("(mag_err**2+std_zp**2+std_ref**2)^(0.5) %.3f \n" %(mag_all_err))
    text_file.close()

#setting data to todays data in format YYYYMMDD
t = date.today()
TIME = datetime.datetime.now().strftime("%H:%M:%S")
year,month,dayy = t.strftime("%Y"),t.strftime("%m"),t.strftime("%d")
today = Time(f'{year}-{month}-{dayy} {TIME}')
TODAY = t.strftime("%Y%m%d") #todays date in YYYYMMDD format
apo = Observer.at_site("lapalma")
sun_set_today = apo.sun_set_time(today, which="nearest") #sun set on day of observing
time_suns_today = "{0.iso}".format(sun_set_today)[-12:]
sun_set_tomorrow = apo.sun_set_time(today,which="next")
time_suns_tomorrow = "{0.iso}".format(sun_set_tomorrow)[-12:]

print(t)
if time_suns_today<TIME<'23:59:59':
    date_ = TODAY
    DATE = re.sub("-","",date_)
if '00:00:00'<TIME<time_suns_tomorrow:
    date_ = str(t - datetime.timedelta(days=1))
    DATE=re.sub("-","",date_)

 
#folder is for easy way to display photometry by date of observation
if os.path.exists(path+'photometry_date')==False:
    os.mkdir(path+'photometry_date')
if not os.path.exists(path+f'photometry_date/{DATE}'):
    os.mkdir(path+f'photometry_date/{DATE}')
if not os.path.exists(path+f'photometry_date/{DATE}/cut_outs'):
    os.mkdir(path+f'photometry_date/{DATE}/cut_outs')


if os.path.exists(path+f'photometry_date/{DATE}'):
    text_file = open(path+f'photometry_date/{DATE}/'+sci_img_name[:-11]+"_photometry.txt", "w")
    text_file.write(sci_obj+" "+f"sdss{sci_filt}"+" "+str(mjd)+" %.3f %.3f %.3f %s %s %s \n" % (mag[0],mag_all_err,mag[3],ra_string,dec_string,sci_exp_time))
    text_file.write("S/N %.3f, %.3f sigma \n" % (mag[4],mag[5]))
    text_file.write("1 sigma limit: %.3f\n" % (mag[6]))
    text_file.write("3 sigma limit: %.3f\n" % (mag[7]))
    text_file.write("5 sigma limit: %.3f\n" % (mag[8]))
    text_file.write("mag=%.3f\n" %(mag[9]))
    text_file.write("time photometry: %.1f seconds\n" % (t_end-t_start))
    text_file.write("memory useage: %.1f Mbyte\n" % (memkb/1024.0/1024.0))
    text_file.write("zeropoint reference %.3f sd %.3f \n" % (np.median(zp_ref),np.nanstd(zp_ref)))
    text_file.write("zeropoint science %.3f sd %.3f \n" % (np.median(zp_sci),np.nanstd(zp_sci)))
    text_file.write("(mag_err**2+std_zp**2+std_ref**2)^(0.5) %.3f \n" %(mag_all_err))
    text_file.close()


print(colored('Total time:', 'green'),"%.1f" % round(t_end-t_start, 0),'seconds')

#print(cutouts_dict)

#plt.imshow(cutouts_dict['sci_sub'][0].data,norm=LogNorm(),cmap='gray')
#plt.show()
'''
print(coords_sn)

sn_x=coords_sn[0][0]
sn_y=coords_sn[0][1]

psf_file = os.listdir(path+"convolved_psf")[0]
conv_sci_psf = fits.open(path+f"convolved_psf/{psf_file}")[0].data

fig=plt.figure(figsize=(10,10))

plt.imshow(conv_sci_psf[0].data,)
plt.xlim(sn_x-10,sn_x+10)
plt.ylim(sn_y-10,sn_y+10)
plt.show()
plt.savefig(path+f'photometry_date/{DATE}/cut_outs/'+sci_img_name[:-11]+"_psf.png")
'''
null=None

if mag[0]<90 and sci_filt!='u':
    data = {"filter":f"sdss{sci_filt}","magerr": mag_all_err,"obj_id": sci_obj,"origin": "LT_IOO_PIPE","mag": mag[0],"limiting_mag": mag[3],"mjd": mjd,"instrument_id": 33,"magsys": "ab","group_ids":"all"}
    lt_data().upload_phot(data=data,name=sci_obj)
    print(data)







