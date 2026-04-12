
from astropy.wcs import WCS
import astroscrappy
import requests
from astropy import wcs
from astropy.io import fits
from astropy.io.votable import parse_single_table
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from sedm_credentials import *
import csv
import astroplan
import typing
from astroplan import Observer
from typing import Mapping, Optional
from astropy.time import Time
from bs4 import BeautifulSoup
import re
import pandas as pd
import urllib.parse
import json
from astroquery import sdss, vo_conesearch
from astroquery.sdss import SDSS
from astroquery.vo_conesearch import conesearch
import datetime
import smtplib
from email.mime.text import MIMEText
from tabulate import tabulate
from email.mime.application import MIMEApplication
from os.path import basename
from email.mime.multipart import MIMEMultipart
from termcolor import colored
from astropy.table import Table
import time
from datetime import date,timedelta
from io import StringIO
from colorama import Fore,Style
# import astroalign as aa
from reproject import reproject_exact,reproject_interp
import glob
# term_color = {0:Fore.CYAN,1:Fore.GREEN,2:Fore.YELLOW,3:Fore.MAGENTA,4:Fore.WHITE,5:Fore.BLUE}Style.RESET_ALL
# f"{Fore.RED}[WARNING]  ::{Style.RESET_ALL}"+
# f"{Fore.YELLOW}[WARNING]  ::{Style.RESET_ALL}"+
# f"{Fore.GREEN}[INFO]  ::{Style.RESET_ALL}"+
# f"{Fore.CYAN}[INFO]  ::{Style.RESET_ALL}"+
info_g = f"{Fore.GREEN}[INFO]    ::{Style.RESET_ALL}"
info_b = f"{Fore.CYAN}[INFO]    ::{Style.RESET_ALL}"
warn_r = f"{Fore.RED}[WARNING] ::{Style.RESET_ALL}"
warn_y = f"{Fore.YELLOW}[WARNING] ::{Style.RESET_ALL}"
process_g = f"{Fore.GREEN}[PROCESS] ::{Style.RESET_ALL}"

def estimate_seeing(filename):
    print(info_g+f' Estimating the seeing based on the mode of the FWHM of point sources in the field')
    sextracted = sextract(filename,0, 0, 3, 12, maxellip=0.7, saturation=-1,delete=True)
    see = mode([sex.fwhm for sex in sextracted])
    print(info_g+f' Estimated seeing: ',see )
    return see




def writeswarpdefaultconfigfile(verbose='QUIET',RA='00:00:00.0',DEC='00:00:00.0'):
    verbose=verbose.upper()
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
CENTER                '''+RA+''','''+DEC+''' # Coordinates of the image center
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
VERBOSE_TYPE           '''+verbose+'''          # QUIET,LOG,NORMAL, or FULL
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


def writeswarpconfigfile(verbose='QUIET'):
    verbose=verbose.upper()
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
 
RESAMPLE               N               # Resample input images (Y/N)?
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
VERBOSE_TYPE           '''+verbose+'''          # QUIET,LOG,NORMAL, or FULL
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

def writesdssswarpconfigfile(verbose='quiet'):
    verbose=verbose.upper()
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
VERBOSE_TYPE           '''+verbose+'''           # QUIET,NORMAL or FULL
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
    pf = open(path+'config_files/default.param','w')
    pf.write(params)
    pf.close()

def prepsexfile(verbose='QUIET',gain=1.6):
    verbose=verbose.upper()
    params = '''# Simple configuration file for SExtractor prior to PSFEx use
# only non-default parameters are present.
# EB 2007-08-01
#

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME     '''+data1_path+'''config_files/prepsfex.cat   # Catalog filename
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
GAIN             '''+str(gain)+'''            # <- put the detector gain in e-/ADU here

#------------------------- Star/Galaxy Separation ----------------------------
#------------------------------ Background -----------------------------------
#------------------------------ Check Image ----------------------------------
#--------------------- Memory (change with caution!) -------------------------
#------------------------------- ASSOCiation ---------------------------------
#----------------------------- Miscellaneous ---------------------------------'''
    pf = open(path+'config_files/prepsfex.sex','w')
    pf.write(params)
    pf.close()

def psfexfile(verbose='QUIET'):
    verbose=verbose.upper()
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
OUTCAT_NAME        '''+data1_path+'''out/psfex_out.cat  # Output catalog filename




#------------------------------- Check-plots ----------------------------------

CHECKPLOT_DEV       NULL         # NULL, XWIN, TK, PS, PSC, XFIG, PNG,
                                # JPEG, AQT, PDF or SVG
CHECKPLOT_TYPE      FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS # or NONE
CHECKPLOT_NAME      fwhm, ellipticity, counts, countfrac, chi2, resi

#------------------------------ Check-Images ---------------------------------

CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL
                                # Check-image types
CHECKIMAGE_NAME '''+data1_path+'''out/chi.fits,'''+data1_path+'''out/proto.fits,'''+data1_path+'''out/samp.fits,'''+data1_path+'''out/resi.fits,'''+data1_path+'''out/snap.fits,'''+data1_path+'''out/moffat.fits,'''+data1_path+'''out/submoffat.fits,'''+data1_path+'''out/subsym.fits
                                # Check-image filenames
CHECKIMAGE_CUBE Y

# CHI (square-root of) chi^2 maps for all input vignettes

#----------------------------- Miscellaneous ---------------------------------

PSF_DIR         '''+data1_path+'''out/ # Where to write PSFs (empty=same as input)
PSF_SUFFIX      .psf            # Filename extension for output PSF filename
VERBOSE_TYPE    '''+verbose+'''          # can be QUIET,NORMAL,LOG or FULL
WRITE_XML       Y               # Write XML file (Y/N)?
XML_NAME        '''+data1_path+'''out/psfex.xml       # Filename for XML output
NTHREADS        0               # Number of simultaneous threads for
                                # the SMP version of PSFEx
                                # 0 = automatic'''
    pf = open(path+'config_files/psfex_conf.psfex','w')
    pf.write(params)
    pf.close()




##################################
# PS1 QUERYY

def panstarrs_query(ra_deg, dec_deg, rad_deg, logger=None,
                    mindet=1, 
                    maxsources=10000,
                    server=('https://archive.stsci.edu/panstarrs/search.php')): 
    """
    Query Pan-STARRS DR1 @ MAST
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field 
                                          radius in degrees
                mindet: minimum number of detection (optional)
                maxsources: maximum number of sources
                server: servername
    returns: astropy.table object
    """
    if not os.path.exists(data1_path+'ps_catalogs'):
      os.makedirs(data1_path+'ps_catalogs')
    

    if not os.path.exists(data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml'):  
        r = requests.get(server, 
                params= {'RA': ra_deg, 'DEC': dec_deg, 
                'SR': rad_deg, 'max_records': maxsources, 
                'outputformat': 'VOTable', 
                'ndetections': ('>%d' % mindet)})
                
        outf = open(data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml', 'w') 
        outf.write(r.text) 
        outf.close() 
        if logger!=None:logger.info(info_g+f" PS1 Catalog downloaded to: "+data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml')
        else:print(info_g+f" PS1 Catalog downloaded to: "+data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml')
    else:
        if logger!=None:logger.info(info_g+f" PS1 Catalog already downloaded to: "+data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml')
        else:print(info_g+f" PS1 Catalog already downloaded to: "+data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml')
    # write query data into local file
    # parse local file into astropy.table object 
    data = parse_single_table(data1_path+'ps_catalogs/ps_'+str(ra_deg)+'_'+str(dec_deg)+'_'+str(rad_deg)+'.xml')
    return data.to_table(use_names_over_ids=True) 

##################################
# SDSS QUERYY
def sdss_query(ra_deg, dec_deg, rad_deg):
    
    data = []
    catra=[]
    catdec=[]
    catmag=[]
    # Get catalog from USNO
    queryurl = 'http://skyserver.sdss.org/dr16/en/tools/search/x_results.aspx?searchtool=Radial&TaskName=Skyserver.Search.Radial&whichphotometry=optical&coordtype=equatorial&ra='+str(ra_deg)+'&dec='+str(dec_deg)+'&radius='+str(rad_deg)+'&min_u=0&max_u=30&min_g=0&max_g=30&min_r=0&max_r=30&min_i=0&max_i=30&min_z=0&max_z=30&min_j=0&max_j=30&min_h=0&max_h=30&min_k=0&max_k=30&format=csv&limit=5000'
    # print(queryurl)
    with requests.Session() as s:
        download = s.get(queryurl)
        #  print
        decoded_content = download.content.decode('utf-8')
        cr = csv.reader(decoded_content.splitlines(), delimiter=',')
        fields2= next(cr)
        #print(fields2)
        #cr.readlines()
        #  print(cr.readlines())
        # print(cr)
        for ind,row in enumerate(cr):
            #  print(row)
            #  print(ind)
            if row[6]=='6':# its a star
                catra.append(float(row[7]))
                catdec.append(float(row[8]))
                catmag.append(float(row[9]))
                data.append([float(row[7]),float(row[8]),float(row[9])])

    print(info_g+" Catalog SDSS dr16; length:",len(catra))       
    cat_table = pd.DataFrame(data=data,columns=['ra','dec','mag'])
    # print(cat_table)
    #return catra,catdec,catmag
    return cat_table

# def sdss_query_image(ra_string,dec_string,filt,nx,ny,log=None):


#     print(info_g+f" Querrying SDSS for reference imaging in {filt}-band",(nx,ny))
#     sdss_url='https://dr12.sdss.org'
#     url=sdss_url+'/fields/raDec?ra='+str(ra_string)+'&dec='+str(dec_string)
#     # print(url)
#     # sys.exit(1)
#     #html_page = urllib2.urlopen(url) python2.7 version
#     html_page = urllib.request.urlopen(url)

#     image_link=[]
#     soup = BeautifulSoup(html_page)
#     for link in soup.findAll('a'):
#         if 'frame-'+str(filt) in link.get('href'):
#             image_link.append(sdss_url+link.get('href'))

#     try:
#         image_link=image_link[0]
#         image_name=image_link.rsplit('/', 1)[-1]
#         # sys.exit(1)
#     except IndexError:
#         print(warn_r+f' Exiting... Not in the SDSS footprint, no {filt}-band!')
#         return

#     try:
#         if filt=='i': 
#             image_link = re.sub('irg','i',image_link)
#             image_link = re.sub('.jpg','.fits.bz2',image_link)
#             image_name = image_link.rsplit('/', 1)[-1]
#             # sys.exit(1)
        
#         # print(image_link)
#         # sys.exit(1)

#         r=requests.get(image_link)
#         r.raise_for_status()
#         if os.path.exists(image_name[:-4]):
#             # if self.termoutp!='quiet':
#             print(info_b+' SDSS image already downloaded',image_name[:-4], (nx,ny))
#         if not os.path.exists(image_name[:-4]):
#             zname=image_name
#             zfile = open(data1_path+'ref_imgs/'+image_name, 'wb')
#             zfile.write(r.content)
#             zfile.close()
#             os.system('bzip2 -d '+data1_path+'ref_imgs/'+image_name )
#             # if self.termoutp!='quiet':
#             print(info_g+' Downloading new SDSS ',str(filt),'-band..',image_name[:-4], (nx,ny))
#             ref_path=data1_path+'ref_imgs/'+image_name[:-4]
#             # os.system('rm '+path+'ref_imgs/'+image_name+'.bz2')
#     except requests.exceptions.HTTPError as err:
#         print(warn_r+' Not in SDSS footprint! Exiting..')

#         return
#     # sys.exit(1)
#     return(ref_path)

def submit_sql_query(ra,dec):
    url = 'https://skyserver.sdss.org/dr18/en/tools/search/x_results.aspx'
    # Your SQL query
    query = f'''SELECT TOP 50 
p.run,p.rerun,p.camCol,p.field,p.obj 
FROM ..PhotoObj AS p 
JOIN dbo.fGetNearbyObjEq({str(ra)},{str(dec)},5) AS b ON  b.objID = P.objID 
WHERE  ( p.type = 3 OR p.type = 6) '''
    payload = {
        'cmd': query,
        'format': 'json',
        'searchtool': 'SQL',
        'taskname': 'Skyserver.Search.SQL',
        'syntax': 'NoSyntax',
        'returnquery': 'true'
    }
    
    response = requests.post(url, data=payload)
    
    if response.status_code == 200:
        return response.json()
    else:
        return f"Error: {response.status_code}, {response.text}"

def response_to_dataframe(response,filt='u'):
    if isinstance(response, list) and len(response) > 0:
        first_item = response[0]
        if isinstance(first_item, dict) and 'Rows' in first_item:
            df = pd.DataFrame(first_item['Rows'])
            if len(df) == 0:return None
            # print(df)
            df['frame_name'] = 'http://dr16.sdss.org/sas/dr16/eboss/photoObj/frames/'+df['rerun'].astype(str)+'/'+df['run'].astype(str)+'/'+df['camCol'].astype(str)+'/frame-'+filt+\
                '-'+df['run'].astype(str).str.zfill(6)+'-'+df['camCol'].astype(str)+'-'+df['field'].astype(str).str.zfill(4)+'.fits.bz2'
            # df = df.loc[requests.head(df['frame_name'].values).status_code==200]
            drop_list = []
            for i in range(len(df)):
                # print(df['frame_name'].values[i],requests.head(df['frame_name'].values[i]).status_code)
                if requests.head(df['frame_name'].values[i]).status_code not in [200,301]:
                    # df.drop(i,inplace=True)
                    drop_list.append(i)
            df.drop(drop_list,inplace=True)
            return df,df['frame_name'].unique()
    
    print(warn_y+" Error: Unable to convert response to DataFrame")
    # print("Response structure:")
    # print(json.dumps(response, indent=2))
    return None

def sdss_query_image(ra_string,dec_string,filt,nx,ny,log=None,lnks_done=[]): 


    log.info(info_g+f" Querrying SDSS for reference imaging in {filt}-band ("+str(nx)+str(',')+str(ny)+f') ({ra_string:.2f},{dec_string:.2f})')

    image_links=[]

    links_sql = submit_sql_query(ra_string,dec_string)
    df,links = response_to_dataframe(links_sql,filt=filt)
    # print(links)
    if len(links)==0:
        log.warning(warn_r+f' Exiting... Not in the SDSS footprint, no {filt}-band!')
        return
    if all([lnk in lnks_done for lnk in links]):
        log.info(info_b+f' SDSS image already downloaded'+str((nx,ny)))
        return
    for link in links:
        image_links.append(link)
    # try:
    for image_link in image_links:
        if filt=='i': 
            image_link = re.sub('irg','i',image_link)
            image_link = re.sub('.jpg','.fits.bz2',image_link)
            image_name = image_link.rsplit('/', 1)[-1]
            # sys.exit(1)
        
        print(image_link)
        # sys.exit(1)
        image_name=image_link.rsplit('/', 1)[-1]
        r=requests.get(image_link)
        r.raise_for_status()
        if os.path.exists(data1_path+'ref_imgs/'+image_name[:-4]):
            # if self.termoutp!='quiet':
            log.info(info_b+f' SDSS image already downloaded'+image_name[:-4]+str((nx,ny)))
        if not os.path.exists(image_name[:-4]):
            zname=image_name
            zfile = open(data1_path+'ref_imgs/'+image_name, 'wb')
            zfile.write(r.content)
            zfile.close()
            os.system('bzip2 -d '+data1_path+'ref_imgs/'+image_name )
            # if self.termoutp!='quiet':
            log.info(info_g+' Downloading new SDSS '+str(filt)+'-band..'+image_name[:-4]+str((nx,ny)))
            ref_path=data1_path+'ref_imgs/'+image_name[:-4]
            # os.system('rm '+path+'ref_imgs/'+image_name+'.bz2')
    # except requests.exceptions.HTTPError as err:
    #     log.warning(warn_r+' Migth not be in SDSS footprint! Exiting..')

    #     return
    # sys.exit(1)
    return ref_path, links





def load_wcs_from_file(filename,coord):
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename)    
    w = wcs.WCS(hdulist[0].header)
    pixcrd = np.array(coord, np.float_)
    world = w.wcs_pix2world(pixcrd, 1)
    return (world)

def wcs_to_pixels(filename,coord):
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename)    
    w = wcs.WCS(hdulist[0].header)
    world= np.array(coord, np.float_)
    pixcrd = w.wcs_world2pix(world, 1)
    return (pixcrd)



 

ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"



def api(method,endpoint,data=None):
    headers={"Authorization":f"token {token}"}
    response=requests.request(method.upper(),urllib.parse.urljoin("https://fritz.science", endpoint),json=data,headers=headers)

    return response

def save_to_group(obj_id,group_id):
    url = "https://fritz.science/api/source_groups"

    payload = {
        "objId": obj_id,
        "inviteGroupIds": [group_id],
            }
    headers = {
            "Content-Type": "application/json",
            "Authorization": f"token {token}"
        }

    response = requests.post(url, json=payload, headers=headers)
    if response.status_code==200:
        print(info_g+f" Saved to group {group_id}:",response.json())
    else:
        print(warn_r+f" Failed to save to group {group_id}:",response.json())

    return

def SN_data_phot(name):
    '''Retrieves individual events, form needs to be an array or a list, puts data into a dataframe'''
    data_array = []

    if 'ZTF' in name:
        while True:
            response = api("get", f"api/sources/{name}/photometry")
            if response.status_code==200:
                break
        data = response.json().get("data", None)
        if len(data) ==0:
            df = pd.DataFrame()
        else:
            df = pd.DataFrame(data).sort_values(by=['mjd'], ascending=True)

        return df

class subphot_data():
    def batch_update_astrometry(self,file_names):
        
        orig_file = fits.open(file_names[0])
        orig_header, orig_data = orig_file[0].header, orig_file[0].data
        output_files = [file_names[i].replace('.fits','_anet.fits').split('/')[-1] for i in range(len(file_names))]
        orig_height,orig_width = orig_header['NAXIS1'],orig_header['NAXIS2']
        try:
            ra,dec = orig_header['CAT-RA'],orig_header['CAT-DEC']
        except:
            try:
                ra,dec = orig_header['OBJCTRA'],orig_header['OBJCTDEC']
            except:
                ra,dec = orig_header['OBJRA'],orig_header['OBJDEC']
        sf_args = f" --new-fits {str(None)} --index-xyls {str(None)} --axy {str(None)} --scamp {str(None)} --corr {str(None)} --rdl {str(None)} --match {str(None)} --solved {str(None)} --no-plots --no-verify "


        solve_field_command = solve_field_path+' '+' '.join(file_names)+' '+'--ra '+str(ra)+' '+'--dec '+str(dec)+' '+'--overwrite '+'-o '+' '.join(output_files)+' '+'--backend-config '+solve_field_config_path+' '+'--scale-low 0.25 --scale-units degwidth '#+sf_args
        print(solve_field_command)
        os.system(solve_field_command)
        sys.exit()

        # new_file = 
        new_file = fits.open(file_name.replace('.fits','_anet.new'))
        new_header = new_file[0].header




        keywords = ['AP_0_0', 'AP_0_1', 'AP_0_2', 'AP_1_0', 'AP_1_1', 'AP_2_0',
       'AP_ORDER', 'A_0_0', 'A_0_1', 'A_0_2', 'A_1_0', 'A_1_1', 'A_2_0',
       'A_ORDER', 'BP_0_0', 'BP_0_1', 'BP_0_2', 'BP_1_0', 'BP_1_1',
       'BP_2_0', 'BP_ORDER', 'B_0_0', 'B_0_1', 'B_0_2', 'B_1_0', 'B_1_1',
       'B_2_0', 'B_ORDER', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1',
       'CDELT2', 'CRPIX1', 'CRPIX1', 'CRPIX2', 'CRPIX2', 'CRVAL1',
       'CRVAL1', 'CRVAL2', 'CRVAL2', 'CTYPE1', 'CTYPE1', 'CTYPE2',
       'CTYPE2', 'CUNIT1', 'CUNIT1', 'CUNIT2', 'CUNIT2', 'EQUINOX',
       'LATPOLE', 'LATPOLE', 'LONPOLE', 'LONPOLE', 'MJDREF', 'PC1_1',
       'PC1_2', 'PC2_1', 'PC2_2', 'PROJP1', 'PROJP3', 'PV1_1', 'PV1_2',
       'PV2_1', 'PV2_2', 'RADECSYS', 'RADESYS', 'WAT0_001', 'WAT1_001',
       'WAT2_001', 'WCSAXES', 'WCSAXES', 'WCSDIM']

        for i in keywords:
            try:
                orig_header[i] = ((wcs_fits[i]),'WCS by APT')
            except:
                continue

        # sys.exit()
        # print(output_file.split('.')[0]+'*')
        temp_files = [i for i in glob.glob(file_name.split('.')[0]+'*') if not i.endswith('.fits') and not i.endswith('.new')]
        # print(temp_files)
        for file in temp_files:
            if file==file_name:
                continue
            os.remove(file)
        
        new_file = new_file
        # new_file[0].header = new_header
        new_upwcs_fits = file_name.replace('.fits','_anet.fits')
        os.remove(new_upwcs_fits.replace('.fits','.new'))
        new_file[0].header = orig_header
        # print('new_upwcs_fits',new_file[0].header)
        # sys.exit()
        new_file.writeto(new_upwcs_fits,overwrite=True)
        return new_upwcs_fits,new_header

    def update_astrometry(self,file_name,ra,dec):
        
    # def update_astrometry_net(self,file_name,ra,dec):
        #this functions uses command line to update the astrometry of a fits file and returns the new fits object with th eupdate header but same data

        orig_file = fits.open(file_name)
        orig_header, orig_data = orig_file[0].header, orig_file[0].data
        output_file = file_name.replace('.fits','_anet.fits').split('/')[-1]
        orig_height,orig_width = orig_header['NAXIS1'],orig_header['NAXIS2']
        sf_args = f" --new-fits {str(None)} --index-xyls {str(None)} --axy {str(None)} --scamp {str(None)} --corr {str(None)} --rdl {str(None)} --match {str(None)} --solved {str(None)} --no-plots --no-verify "


        solve_field_command = solve_field_path+' '+file_name+' '+'--ra '+str(ra)+' '+'--dec '+str(dec)+' '+'--overwrite '+'-o '+output_file+' '+'--backend-config '+solve_field_config_path+' '+'--scale-low 0.25 --scale-units degwidth '#+sf_args
        print(solve_field_command)
        os.system(solve_field_command)


        # new_file = 
        new_file = fits.open(file_name.replace('.fits','_anet.new'))
        new_header = new_file[0].header

        # keywords = ['CDELT1','CDELT2',]

        # for i in keywords:
        #     try:
        #         new_header[i] = ((orig_header[i]),'WCS by APT')
        #     except Exception as e:
        #         print(e)
        #         continue

        # print(new_header)
        # return
        # new_file[0].header = new_header
        # new_file.writeto(file_name.replace('.fits','_anet.fits'),overwrite=True)


        # return file_name.replace('.fits','_anet.fits'),new_header


        keywords = ['AP_0_0', 'AP_0_1', 'AP_0_2', 'AP_1_0', 'AP_1_1', 'AP_2_0',
       'AP_ORDER', 'A_0_0', 'A_0_1', 'A_0_2', 'A_1_0', 'A_1_1', 'A_2_0',
       'A_ORDER', 'BP_0_0', 'BP_0_1', 'BP_0_2', 'BP_1_0', 'BP_1_1',
       'BP_2_0', 'BP_ORDER', 'B_0_0', 'B_0_1', 'B_0_2', 'B_1_0', 'B_1_1',
       'B_2_0', 'B_ORDER', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1',
       'CDELT2', 'CRPIX1', 'CRPIX1', 'CRPIX2', 'CRPIX2', 'CRVAL1',
       'CRVAL1', 'CRVAL2', 'CRVAL2', 'CTYPE1', 'CTYPE1', 'CTYPE2',
       'CTYPE2', 'CUNIT1', 'CUNIT1', 'CUNIT2', 'CUNIT2', 'EQUINOX',
       'LATPOLE', 'LATPOLE', 'LONPOLE', 'LONPOLE', 'MJDREF', 'PC1_1',
       'PC1_2', 'PC2_1', 'PC2_2', 'PROJP1', 'PROJP3', 'PV1_1', 'PV1_2',
       'PV2_1', 'PV2_2', 'RADECSYS', 'RADESYS', 'WAT0_001', 'WAT1_001',
       'WAT2_001', 'WCSAXES', 'WCSAXES', 'WCSDIM']

        for i in keywords:
            try:
                orig_header[i] = ((wcs_fits[i]),'WCS by APT')
            except:
                continue

        # sys.exit()
        # print(output_file.split('.')[0]+'*')
        temp_files = [i for i in glob.glob(file_name.split('.')[0]+'*') if not i.endswith('.fits') and not i.endswith('.new')]
        # print(temp_files)
        for file in temp_files:
            if file==file_name:
                continue
            os.remove(file)
        
        new_file = new_file
        # new_file[0].header = new_header
        new_upwcs_fits = file_name.replace('.fits','_anet.fits')
        os.remove(new_upwcs_fits.replace('.fits','.new'))
        new_file[0].header = orig_header
        # print('new_upwcs_fits',new_file[0].header)
        # sys.exit()
        new_file.writeto(new_upwcs_fits,overwrite=True)
        return new_upwcs_fits,new_header

    def stack_images(self,images,save_name='',clean=True,s_type='average'):
        # print(images)
        # print(save_name)

        images_ = []
        for k in range(len(images)):
            if not images[k].startswith(data1_path) and not images[k].startswith('/'):
                images_.append(data1_path+images[k])
            else:
                images_.append(images[k])

        images = images_
        # print(images)

        # Load the first image as reference
        primary_image = fits.open(images[0])
        primary_header,primary_data = primary_image[0].header,primary_image[0].data

        if clean:
            if primary_header['EXPTIME'] >60.0:
                cr_sigclip,cr_sigfrac,cr_objlim,cr_gain,cr_readnoise,cr_satlevel = 5.0,0.3,5.0,primary_header['GAIN'],8.0,65536.0
                if isinstance(cr_gain,float):cr_gain=cr_gain
                else: cr_gain=2.2
                cr_niter,cr_sepmed,cr_cleantype,cr_fsmode,cr_psfmodel,cr_psffwhm,cr_psfsize,cr_psfk,cr_psfbeta,cr_verbose=4,True,'meanmask','median','gauss',2.5,7,None,4.765,False
            else:
                cr_sigclip,cr_sigfrac,cr_objlim,cr_gain,cr_readnoise,cr_satlevel = 5.0,0.3,5.0,primary_header['GAIN'],8.0,65536.0
                if isinstance(cr_gain,float):cr_gain=cr_gain
                else: cr_gain=2.2
                cr_niter,cr_sepmed,cr_cleantype,cr_fsmode,cr_psfmodel,cr_psffwhm,cr_psfsize,cr_psfk,cr_psfbeta,cr_verbose=2,True,'meanmask','median','gauss',2.5,7,None,4.765,False

            cr_primary_mask,cr_primary_data=astroscrappy.detect_cosmics(primary_data,sigclip=cr_sigclip, sigfrac=cr_sigfrac, objlim=cr_objlim, gain=cr_gain, readnoise=cr_readnoise, satlevel=cr_satlevel, niter=cr_niter, 
                                            sepmed=cr_sepmed, cleantype=cr_cleantype, fsmode=cr_fsmode, psfmodel=cr_psfmodel, psffwhm=cr_psffwhm, psfsize=cr_psfsize, psfk=cr_psfk, psfbeta=cr_psfbeta, verbose=cr_verbose)
        
            primary_data = cr_primary_data

        # Initialize an empty array to store aligned images
        aligned_images = [primary_data]
        shifted_data = {}
        shifted_data['primary'] = primary_data


        # Align the remaining images to the reference image
        for im_path in images[1:]:
            image = fits.getdata(im_path)
            try:
                aligned_image,fp = aa.register(np.array(image,dtype="<f4"), np.array(primary_data,dtype="<f4"), fill_value=np.nan)
            except:
                aligned_image,fp = reproject_interp(fits.open(im_path)[0], primary_header)
                aligned_image[np.isnan(aligned_image)] = np.nanmedian(aligned_image)

            shifted_data[path] = aligned_image
            aligned_images.append(aligned_image)

        # Perform the average stacking
        if s_type == 'median':
            stacked_image = np.median(aligned_images, axis=0)
        elif s_type == 'mean':
            stacked_image = np.mean(aligned_images, axis=0)
        elif s_type == 'stack':
            stacked_image = np.stack(aligned_images, axis=0)
            stacked_image = np.nanmedian(stacked_image,axis=0)
        elif s_type == 'average':
            stacked_image = np.average(aligned_images, axis=0)

        if save_name != '':
            fits.writeto(save_name,stacked_image,header=primary_header,overwrite=True)


        return stacked_image,shifted_data

    def down_archive(self,names='',RA='',DEC='', proposals='', dates=['All'],instrument='IO:O',downl=True):  
        #dictionary for archive proposal indexing, format proposal:[password, archive data index, archive cat.txt index]


        final_dict = {}  #dictionary containing all fits asscoiated with given names
        if proposals == '' and RA == '' and DEC == '' and names != '':  #if no proposals specified then will search through all proposals in all_proposals
            proposals = proposals_arc.keys()
            print(info_b+f' No proposal specified, searching all proposals for {str(names)}')
            print()

        if proposals != '' and RA != '' and DEC != '' and names == '':
            proposals = proposals
            print(info_g+f' Searching {proposals} for objects at RA={RA} DEC={DEC}')
            print()

        if proposals == '' and RA != '' and DEC != '' and names == '':
            proposals = proposals_arc.keys()
            print(info_b+f' No proposal specified, searching all proposals for objects at RA={RA} DEC={DEC}')
            print()
            rah,ram,ras = RA[0:2],RA[3:5],RA[6:8]
            dech,decm,decs = DEC[1:3],DEC[4:6],DEC[7:9]

            #searching within +-2 as coordinates differ slightly
            RA0=RA[:-3]
            RA1 = f'{rah}:{ram}:{float(ras)-3}'[:-2]
            RA2 = f'{rah}:{ram}:{float(ras)-2}'[:-2]
            RA3 = f'{rah}:{ram}:{float(ras)-1}'[:-2]
            RA4 = f'{rah}:{ram}:{float(ras)+1}'[:-2]
            RA5 = f'{rah}:{ram}:{float(ras)+2}'[:-2]
            RA6 = f'{rah}:{ram}:{float(ras)+3}'[:-2]
            right_ascension = [RA0,RA1,RA2,RA3,RA4,RA5,RA6]

            dech,decm,decs = DEC[1:3],DEC[4:6],DEC[7:9]
            DEC0=DEC[:-2]
            DEC1 = f'{dech}:{decm}:{float(ras)-3}'[:-2]
            DEC2 = f'{dech}:{decm}:{float(ras)-2}'[:-2]
            DEC3 = f'{dech}:{decm}:{float(ras)-1}'[:-2]
            DEC4 = f'{dech}:{decm}:{float(ras)+1}'[:-2]
            DEC5 = f'{dech}:{decm}:{float(ras)+2}'[:-2]
            DEC6 = f'{dech}:{decm}:{float(ras)+3}'[:-2]
            declination = [DEC0,DEC1,DEC2,DEC3,DEC4,DEC5,DEC6]


        for proposal in proposals:

            #opening specific proposal in archive
            url1 = f'https://telescope.livjm.ac.uk/data/archive/scratch/{proposals_arc[proposal][1]}/{proposal}/Priv/'
            response1 = requests.get(url1,auth=(proposal,proposals_arc[proposal][0]))

            #opening proposal catalog
            url_cat = f'https://telescope.livjm.ac.uk/data/archive/scratch/{proposals_arc[proposal][1]}/cat{proposals_arc[proposal][1]}.txt'
            response_cat = requests.get(url_cat,auth=(proposal,proposals_arc[proposal][0]))
            # print(response_cat.status_code)

            #saving catalog for later reference
            soup_cat = BeautifulSoup(response_cat.content, 'html.parser')
            np.savetxt(f'{path}{proposal}_cat.txt', soup_cat, delimiter=',', fmt='%s')

            if names!='':
                for n in range(len(names)):
                    name = names[n]
                    #opening and reading catalog
                    proposal_catalog = open(f'{path}{proposal}_cat.txt')
                    catalog_arr = []
                    for line in proposal_catalog.readlines()[1::]: #extracting fits name and appending to array

                        if dates[0] == 'All':
                            if name in line and instrument in line:    #if correct name and correct instrument
                                if line[21]==' ':
                                    catalog_arr.append(line[0:21])
                                if line[21]!=' ':
                                    catalog_arr.append(line[0:22])
        
                        if dates[0] != 'All':
                            for d in range(len(dates)):
                                date = dates[d]        #if correct name and correct instrument and correct date/s
                                if name in line and instrument in line and date in line: 
                                    catalog_arr.append(line[0:21])

                    #creating array of all urls for later use
                    all_urls = []
                    for fits in catalog_arr:
                        fit_url = re.sub(' ','',f'{url1}{fits}.fits')
                        all_urls.append(fit_url)

                final_dict[f'{proposal}: {name}'] = all_urls    #specifying proposal in case same object in multiple proposals


            elif names == '':
                #opening and reading catalog
                proposal_catalog = open(f'{path}{proposal}_cat.txt')
                catalog_arr = []
                for line in proposal_catalog.readlines()[1::]: #extracting fits name and appending to array
                    if RA != '' and DEC != '':
                        if dates[0] == 'All':
                            if any(val1 in line for val1 in right_ascension)==True or any(val2 in line for val2 in declination)==True and instrument in line:
                                if line[21]==' ':
                                    catalog_arr.append(line[0:21])
                                if line[21]!=' ':
                                    catalog_arr.append(line[0:22])
                                name = re.sub(' ','',line[-15:])
                                name = re.sub('\n','',name)
                                while name[0] == ' ':     #removing spaces
                                    name = re.sub(' ','',name)
                                    name = re.sub('\n','',name)

                        if dates[0] != 'All':
                            for d in range(len(dates)):
                                date = dates[d]        #if correct name and correct instrument and correct date/s
                                if name in line and instrument in line and date in line: 
                                    if any(val1 in line for val1 in right_ascension)==True or any(val2 in line for val2 in declination)==True and instrument in line:
                                        if line[21]==' ':
                                            catalog_arr.append(line[0:21])
                                        if line[21]!=' ':
                                            catalog_arr.append(line[0:22])
                                        while name[0] == ' ':     #removing spaces
                                            name = re.sub(' ','',name)
                                            name = re.sub('\n','',name)

                if len(catalog_arr)!=0:
                    all_urls = []
                    for fits in catalog_arr:
                        fit_url = re.sub(' ','',f'{url1}{fits}.fits')
                        all_urls.append(fit_url)

                    final_dict[f'{proposal}: {name}'] = all_urls    #specifying proposal in case same object in multiple proposals



        for prop_name in final_dict.keys():
            new_folder = re.sub(' ','',prop_name[9:])
            if os.path.exists(f'{path}data/Archive/{new_folder}') == False:
                os.mkdir(f'{path}data/Archive/{new_folder}')

            proposal = re.sub(':','',prop_name[:9])
            proposal = re.sub(' ','',proposal)
            for url in final_dict[prop_name]:
                response_fits = requests.get(f'{url}', auth=(proposal,proposals_arc[proposal][0]),allow_redirects=False)
                if str(response_fits.status_code)=='200':
                    fits_name = re.sub(f'https://telescope.livjm.ac.uk/data/archive/scratch/{proposals_arc[proposal][1]}/{proposal}/Priv/','',url)
             
                    open(f'{path}data/Archive/{new_folder}/{fits_name}','wb').write(response_fits.content)
                elif str(response_fits.status_code)=='404':
                    response_fits = requests.get(f'{url}.gz', auth=(proposal,proposals_arc[proposal][0]),allow_redirects=False)
                    fits_name = re.sub(f'https://telescope.livjm.ac.uk/data/archive/scratch/{proposals_arc[proposal][1]}/{proposal}/Priv/','',url)
                    open(f'{path}data/Archive/{new_folder}/{fits_name}.gz','wb').write(response_fits.content)
                    gunzip_command = f"gunzip {path}data/Archive/{new_folder}/{fits_name}.gz"
                    os.system(gunzip_command)
                    
                    

            if len(final_dict[prop_name]) != 0:
                print(info_g+f' Found {len(final_dict[prop_name])} observations in {new_folder}, downloaded to data/Archive/{new_folder}')
                print()
            else:
                print(warn_y+f' Found no observations of {new_folder} in {prop_name[:7]} with {instrument}')
                print()

        print(info_g+f' Completed download for {name}')
        return 



    def down_recentdata(self,day=''):
        if not os.path.exists(f'{data1_path}RecentData'):
            os.mkdir(f'{data1_path}RecentData')
        
        #setting data to todays data in format YYYYMMDD
        t = date.today()
        TIME = datetime.datetime.now().strftime("%H:%M:%S")
        TODAY = t.strftime("%Y%m%d") #todays date in YYYYMMDD format
        
        #setting correct date of observation
        year,month,dayy = t.strftime("%Y"),t.strftime("%m"),t.strftime("%d")
        today = Time(f'{year}-{month}-{dayy} {TIME}')
        apo = Observer.at_site("lapalma")
        sun_set_today = apo.sun_set_time(today, which="nearest") #sun set on day of observing
        time_suns_today = "{0.iso}".format(sun_set_today)[-12:]
        sun_set_tomorrow = apo.sun_set_time(today,which="next")
        time_suns_tomorrow = "{0.iso}".format(sun_set_tomorrow)[-12:]

        
        if time_suns_today<TIME<'23:59:59' and day == '':
            DAY = TODAY
        if '00:00:00'<TIME<time_suns_tomorrow and day == '':
            DAY = lt_data.yesterday(TODAY)
  
        if day !='':
            DAY = str(day)

        for proposal in proposals_arc.keys():
            try:
                if os.path.exists(f'{data1_path}RecentData/{DAY}') == False:
                    os.mkdir(f'{data1_path}RecentData/{DAY}')
                url1 = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}'
                response1 = requests.get(url1,auth=(proposal,proposals_arc[proposal][0]))
                avail = {} #dictionary containing all observations available on quicklook

                for i in range(len(response1.text)):
                    if year in response1.text[i:i+4]:
                        url = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}/{response1.text[i:i+8]}'
                        avail[f'{response1.text[i:i+8]}'] = url

                ordered_avail = {} #ordered dictionary containing all observations available on quicklook
                for m,n in sorted(avail.items(),reverse=True):   #reversing the order to get newest first
                    ordered_avail[f'{m}'] = n

                url2 = f'{ordered_avail[DAY]}/{DAY}_1.tgz'
                print(info_g+f' Data has been taken for {proposal} and uploaded to LT Recent Data')
                print(url2)
                response2 = requests.get(url2,auth=(proposal,proposals_arc[proposal][0]))
                open(f'{data1_path}RecentData/{proposal}_{DAY}_1.tgz','wb').write(response2.content)
                print(info_g+f' Downloaded data as tarball at {datetime.datetime.now().strftime("%H:%M:%S")}')
                print()
            except:
                print(warn_y+f' {proposal} observations not found for last night ({DAY}) on Recent Data')
                print()
                pass
        return





    def down_quicklook(self,day=''):
        if not os.path.exists(f'{data1_path}Quicklook'):
            os.mkdir(f'{data1_path}Quicklook')

        #setting data to todays data in format YYYYMMDD
        t = date.today()
        TIME = datetime.datetime.now().strftime("%H:%M:%S")
        TODAY = t.strftime("%Y%m%d") #todays date in YYYYMMDD format
        
        #setting correct date of observation
        year,month,dayy = t.strftime("%Y"),t.strftime("%m"),t.strftime("%d")
        today = Time(f'{year}-{month}-{dayy} {TIME}')
        apo = Observer.at_site("lapalma")
        sun_set_today = apo.sun_set_time(today, which="nearest") #sun set on day of observing
        time_suns_today = "{0.iso}".format(sun_set_today)[-12:]
        sun_set_tomorrow = apo.sun_set_time(today,which="next")
        time_suns_tomorrow = "{0.iso}".format(sun_set_tomorrow)[-12:]

        
        if time_suns_today<TIME<'23:59:59' and day == '':
            DAY = TODAY
        if '00:00:00'<TIME<time_suns_tomorrow and day == '':
            DAY = re.sub('-','',str(t-timedelta(days=1)))
            

        if day !='':
            if str(type(day))=="<class 'set'>":
                day=list(day)
                DAY = day[0]
            else:
                DAY=day

        
        print(info_g+f' Last nights observation date is {DAY}, looking for data from the following proposals: {", ".join(proposals_arc.keys())}')
        print()
  
        if os.path.exists(f"{data1_path}photometry_date/{DAY}")==False:
            os.mkdir(f"{data1_path}photometry_date/{DAY}")

        
                
        for proposal in proposals_arc.keys():
            new_fits_only=[]
            try:
                if os.path.exists(f'{data1_path}Quicklook/{DAY}') == False:
                    os.mkdir(f'{data1_path}Quicklook/{DAY}')
                    os.mkdir(f'{data1_path}Quicklook/{DAY}/spec')
                #opening quicklook homepage
                url1 = f'https://telescope.livjm.ac.uk/DataProd/quicklook/{proposal}'
                response1 = requests.get(url1,auth=(proposal,proposals_arc[proposal][0]))
                # print(response1.status_code)
                avail = {} #dictionary containing all observations available on quicklook

                #searching for string '2022' in the response and appending the full YYYYMMDD if on quicklook
                for i in range(len(response1.text)):
                    if year in response1.text[i:i+4]:
                        url = f'https://telescope.livjm.ac.uk/DataProd/quicklook/{proposal}/{response1.text[i:i+8]}'
                        avail[f'{response1.text[i:i+8]}'] = url

                ordered_avail = {} #ordered dictionary containing all observations available on quicklook
                for m,n in sorted(avail.items(),reverse=True):   #reversing the order to get newest first
                    ordered_avail[f'{m}'] = n

                #opening quicklook homepage for specififc proposal
                url2 = ordered_avail[DAY]
                response2 = requests.get(url2,auth=(proposal,proposals_arc[proposal][0]))
                print(colored('---------------------------------------------------','green'))
                print(info_g+f' {proposal} observations found on Quickoook for ({DAY})')

                #collecting gzipped names of fits images
                gzip_fits = []
                for j in range(len(response2.text)):
                    if 'fits.gz' in response2.text[j:j+7]:
                        fits_name = re.sub('"','',response2.text[j-22:j+7])
                        if fits_name[0]!='v' or fits_name[0]!='h':
                            fits_name = re.sub('"','',response2.text[j-23:j+7])
                        fits_name = re.sub('>','',fits_name)
                        fits_name = re.sub('	','',fits_name)
                        fits_name = re.sub('=','',fits_name)
                        if fits_name not in gzip_fits:
                            gzip_fits.append(fits_name)
                        #print(fits_name)

                #creating a txt file with all fits.gz files and appending to it everytime it checks
                if os.path.exists(f'{data1_path}Quicklook/{DAY}/{proposal}_{DAY}_fits_log.txt')==False:
                    print(info_g+f' Fits log for {DAY} does not yet exist, creating file and writing')
                    fits_name_txt = open(f'{data1_path}Quicklook/{DAY}/{proposal}_{DAY}_fits_log.txt', 'w+')
                    for name in gzip_fits:
                        name=re.sub(' ','',name)
                        name=re.sub('=','',name)

                        if name not in os.listdir(f'{data1_path}Quicklook/{DAY}/'):
                            #gzip_fits = re.sub('=','',gzip_fits)
                            fits_name_txt.write(f'{name}\n')
                            new_fits_only.append(name)
                    fits_name_txt.close()
                    print(info_g+f' Written {len(new_fits_only)} new fits into {proposal}_{DAY}_fits_log.txt')

                elif os.path.exists(f'{data1_path}Quicklook/{DAY}/{proposal}_{DAY}_fits_log.txt')==True:
                    print(info_g+f' Fits log for {DAY} exists, opening and writing new fits')
                    fits_name_txt = open(f'{data1_path}Quicklook/{DAY}/{proposal}_{DAY}_fits_log.txt', 'a+')
                    fits_in = open(f'{data1_path}Quicklook/{DAY}/{proposal}_{DAY}_fits_log.txt').readlines()
                    for name in gzip_fits:
                        if all(name not in f for f in fits_in):
                            name=re.sub('=','',name)
                            fits_name_txt.write(f'{name}\n')
                            new_fits_only.append(name)
                    fits_name_txt.close()
                    print(info_g+f' Written {len(new_fits_only)} new fits into {proposal}_{DAY}_fits_log.txt')

                #reading in and downloading fits files to data/{TODAY} directory
                down_fits = open(f'{data1_path}Quicklook/{DAY}/{proposal}_{DAY}_fits_log.txt').readlines()

                if len(new_fits_only)>0:
                    for fits in new_fits_only:
                        #fits = fits[:-1]
    
                        if os.path.exists(f'{data1_path}Quicklook/{DAY}/{fits}')==False and os.path.exists(f'{data1_path}Quicklook/{DAY}/{fits}.gz')==False and os.path.exists(f"{data1_path}Quicklook/{DAY}/{re.sub('.gz','',fits)}")==False:# and os.path.exists(f'data/Quicklook/{DAY}/{fits}')==False:
                            url3 = f'https://telescope.livjm.ac.uk/DataProd/quicklook/{proposal}/{DAY}/{fits}'
                            url3 = re.sub('	','',url3)
                            print(url3)

                            response3 = requests.get(url3,auth=(proposal,proposals_arc[proposal][0]))
                            if fits[0]=="h":
                                open(f'{data1_path}Quicklook/{DAY}/{fits}', 'wb').write(response3.content)
                                os.system(f"gunzip {data1_path}Quicklook/{DAY}/{fits}")
                            elif (fits[0:3]=="v_e")==True:
                                open(f'{data1_path}Quicklook/{DAY}/spec/{fits}', 'wb').write(response3.content)
                                os.system(f"gunzip {data1_path}Quicklook/{DAY}/spec/{fits}")
         
                print()
            except Exception as e:
                print(warn_y+f' {proposal} observations not found for last night ({DAY}) on quicklook')
                print()
                pass
        return DAY

    def SN_data_phot(self,name):
        data_array = []

        response = api("get", f"api/sources/{name}/photometry")
        data = response.json().get("data", None)
        if len(data) ==0:
            df = pd.DataFrame(columns=['id','mjd','mag','magerr','limitting_mag','filter','instrument_name','origin','altdata','ra','dec','UTC','instrument_id'])
        else:
            df = pd.DataFrame(data).sort_values(by=['mjd'], ascending=True)

        return df



    









import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter, label, center_of_mass
from skimage import measure
from shapely.geometry import Polygon, MultiPoint
from astropy.io import fits

# ----------------------------
# Helper functions
# ----------------------------
def extract_curved_border(mask):
    """Extract the longest contour from a binary mask (curved border)."""
    contours = measure.find_contours(mask.astype(float), 0.5)
    if len(contours) == 0:
        raise RuntimeError("No valid border found")
    return max(contours, key=len)[:, ::-1]  # convert (row,col) -> (x,y)

def detect_stars_with_size(image, threshold_sigma=5.0, fwhm=3.0, border_mask=None, n_max=None, avoid_radius=None):
    """
    Detect stars and measure approximate size. Avoid overlapping stars.

    Parameters
    ----------
    image : 2D array
        Science image
    threshold_sigma : float
        Detection threshold in sigma above background
    fwhm : float
        Approximate FWHM of stars in pixels (used for max filter)
    border_mask : 2D bool array, optional
        Mask of valid pixels
    n_max : int, optional
        Max stars to keep by flux
    avoid_radius : float, optional
        Pixels to avoid around each detected star

    Returns
    -------
    xy : (N,2) array of centroids
    flux : (N,) array of fluxes
    size : (N,) array of star sizes (approximate sigma)
    """
    # Background statistics
    mean = np.mean(image)
    median = np.median(image)
    std = np.std(image)
    threshold = median + threshold_sigma * std

    # Initialize mask of pixels allowed for detection
    if border_mask is None:
        valid_mask = np.ones_like(image, dtype=bool)
    else:
        valid_mask = border_mask.copy()

    xy_list = []
    flux_list = []
    size_list = []

    while True:
        # Apply maximum filter to find local maxima
        size_filter = int(np.ceil(fwhm))
        local_max = (maximum_filter(image*valid_mask, size=size_filter) == image*valid_mask)
        detect = local_max & (image > threshold) & valid_mask

        labeled, n_labels = label(detect)
        if n_labels == 0:
            break

        # Pick the brightest remaining peak
        peak_flux = 0
        peak_idx = -1
        for region_id in range(1, n_labels+1):
            mask = labeled == region_id
            total_flux = image[mask].sum()
            if total_flux > peak_flux:
                peak_flux = total_flux
                peak_idx = region_id

        if peak_idx == -1:
            break

        mask = labeled == peak_idx
        y_c, x_c = center_of_mass(image, labels=mask, index=1)

        # Compute approximate size: sigma_x, sigma_y
        ys, xs = np.nonzero(mask)
        fluxes = image[ys, xs]
        x_mean = np.sum(xs * fluxes) / np.sum(fluxes)
        y_mean = np.sum(ys * fluxes) / np.sum(fluxes)
        sigma_x = np.sqrt(np.sum(fluxes * (xs - x_mean)**2) / np.sum(fluxes))
        sigma_y = np.sqrt(np.sum(fluxes * (ys - y_mean)**2) / np.sum(fluxes))
        star_size = np.mean([sigma_x, sigma_y])

        xy_list.append([x_c, y_c])
        flux_list.append(np.sum(image[mask]))
        size_list.append(star_size)

        # Mask out pixels within avoid_radius or star_size
        if avoid_radius is not None:
            y_grid, x_grid = np.indices(image.shape)
            dist2 = (x_grid - x_c)**2 + (y_grid - y_c)**2
            mask_radius = dist2 <= (avoid_radius + star_size)**2
            valid_mask[mask_radius] = False
        else:
            # If no avoid_radius, mask only the detected pixels
            valid_mask[mask] = False

        # Stop if reached n_max
        if n_max is not None and len(xy_list) >= n_max:
            break

    return np.array(xy_list), np.array(flux_list), np.array(size_list)

# ----------------------------
# Example usage
# ----------------------------


def get_stars(sci_file,ref_file):
    sci_data, sci_hdr = fits.getdata(sci_file, header=True)
    sci_data = sci_data.astype(float)

    # print(f"10 pixels = {pix_10_arcsec:.2f} arcsec")
    ref_data, ref_hdr = fits.getdata(ref_file, header=True)
    ref_data = ref_data.astype(float)

    # Define valid science pixels
    sci_mask = np.abs(sci_data) > 1e-6
    sci_border_pix = extract_curved_border(sci_mask)

    sci_wcs = WCS(sci_hdr)
    ref_wcs = WCS(ref_hdr)

    # Detect stars, avoid overlapping within 5 pixels
    # print("Detecting stars in science image...")
    xy_sci, flux_sci, size_sci = detect_stars_with_size(
        sci_data, threshold_sigma=3.0, fwhm=3.0, border_mask=sci_mask, n_max=50, avoid_radius=30
    )
    print(f"Detected {len(xy_sci)} stars in science image.")

    # print("Detecting stars in reference image...")
    xy_ref, flux_ref, size_ref = detect_stars_with_size(
        ref_data, threshold_sigma=2.0, fwhm=3.0, border_mask=sci_mask, n_max=50, avoid_radius=20
    )
    print(f"Detected {len(xy_ref)} stars in reference image.")
    #keep only the stars that are inside border -10 pixels
    BL=15
    sci_keep = []
    # print("Filtering stars near border...")
    for i in range(len(xy_sci)):
        xl_border = sci_border_pix[:,0][sci_border_pix[:,1]==xy_sci[i,1]].min()+BL
        xr_border = sci_border_pix[:,0][sci_border_pix[:,1]==xy_sci[i,1]].max()-BL
        yl_border = sci_border_pix[:,1][sci_border_pix[:,0]==xy_sci[i,0]].min()+BL
        yr_border = sci_border_pix[:,1][sci_border_pix[:,0]==xy_sci[i,0]].max()-BL
        # print(xy_sci[i])
        # print(xl_border,xr_border,xl_border<xy_sci[i,0]<xr_border)
        # print(yl_border,yr_border,yl_border<xy_sci[i,1]<yr_border)

            
        if (xl_border<xy_sci[i,0]<xr_border) and (yl_border<xy_sci[i,1]<yr_border):
            sci_keep.append(True)
        else:
            sci_keep.append(False)
    if sum(sci_keep)>4:
        sci_keep = np.array(sci_keep)
        xy_sci = xy_sci[sci_keep]
        flux_sci = flux_sci[sci_keep]
        size_sci = size_sci[sci_keep]
    ref_keep = []
    # print("Filtering reference stars near border...")
    for i in range(len(xy_ref)):
        xl_border = sci_border_pix[:,0][sci_border_pix[:,1]==xy_ref[i,1]].min()+BL
        xr_border = sci_border_pix[:,0][sci_border_pix[:,1]==xy_ref[i,1]].max()-BL
        yl_border = sci_border_pix[:,1][sci_border_pix[:,0]==xy_ref[i,0]].min()+BL
        yr_border = sci_border_pix[:,1][sci_border_pix[:,0]==xy_ref[i,0]].max()-BL
        if (xl_border<xy_ref[i,0]<xr_border) and (yl_border<xy_ref[i,1]<yr_border):
            ref_keep.append(True)
        else:
            ref_keep.append(False)
    
    if sum(ref_keep)>4:
        ref_keep = np.array(ref_keep)
        xy_ref = xy_ref[ref_keep]
        flux_ref = flux_ref[ref_keep]
        size_ref = size_ref[ref_keep]
    radec_sci = sci_wcs.all_pix2world(xy_sci[:,0], xy_sci[:,1], 0)
    radec_ref = ref_wcs.all_pix2world(xy_ref[:,0], xy_ref[:,1], 0)


    return (xy_sci, flux_sci, size_sci, sci_border_pix, sci_data, radec_sci), (xy_ref, flux_ref, size_ref, sci_border_pix, ref_data, radec_ref)

def plot_stars(xy, flux, size, border_pix, data,ax):
    # ----------------------------
    # Plot detected stars
    # ----------------------------
    vmin, vmax = np.percentile(data, [5, 99])
    # fig, ax = plt.subplots(figsize=(10,10))
    ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    ax.plot(border_pix[:,0], border_pix[:,1], 'r-', lw=2, label='Science border')
    # ax.axvline(78.5-25,c='pink')
    # ax.axvline(892.5+25,c='pink')

    for i, (x, y) in enumerate(xy):
        circle = plt.Circle((x, y), 6, color='yellow', fill=False, lw=2.5)
        ax.add_artist(circle)
        ax.text(x+3, y+3, f"{i+1}", color='red', fontsize=12)

    # ax.set_title("Detected Stars with Sizes and Avoidance Radius")
    # ax.legend()
    # plt.show()

    # Print
    # for i in range(len(xy)):
    #     print(f"Star {i+1}: xy=({xy[i,0]:.1f},{xy[i,1]:.1f}), flux={flux[i]:.1f}")

    return ax

def get_cover(sci_border_pix, sci_xy_matched):
    sci_border_polygon = Polygon(sci_border_pix)
    border_poly = Polygon(sci_border_pix)
    if not border_poly.is_valid:
        border_poly = border_poly.buffer(0)

    # Convex hull of matched stars
    hull = MultiPoint(sci_xy_matched).convex_hull

    # Intersection area
    intersection = hull.intersection(border_poly)

    border_area = border_poly.area
    covered_area = intersection.area

    coverage_fraction = covered_area / border_area if border_area > 0 else 0
    return coverage_fraction,hull,intersection
    
def match_stars(radec_sci, radec_ref,xy_sci,xy_ref, border_pix,max_sep_arcsec=2.0):
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    sci_coords = SkyCoord(ra=radec_sci[0]*u.deg, dec=radec_sci[1]*u.deg, frame='fk5')
    ref_coords = SkyCoord(ra=radec_ref[0]*u.deg, dec=radec_ref[1]*u.deg, frame='fk5')

    idx, sep2d, _ = sci_coords.match_to_catalog_sky(ref_coords)

    matched_sci = []
    matched_ref = []
    for i in range(len(sci_coords)):
        if sep2d[i].arcsec < max_sep_arcsec:
            matched_sci.append((sci_coords[i].ra.deg, sci_coords[i].dec.deg))
            matched_ref.append((ref_coords[idx[i]].ra.deg, ref_coords[idx[i]].dec.deg))


    sci_xy_matched = []
    ref_xy_matched = []
    # print(len(sci_radec,len(matched_sci)))
    for sc, rc in zip(matched_sci, matched_ref):
        sci_idx = np.where((radec_sci[0]==sc[0]) & (radec_sci[1]==sc[1]))[0][0]
        ref_idx = np.where((radec_ref[0]==rc[0]) & (radec_ref[1]==rc[1]))[0][0]
        sci_xy_matched.append(xy_sci[sci_idx])
        ref_xy_matched.append(xy_ref[ref_idx])
    sci_xy_matched = np.array(sci_xy_matched)
    ref_xy_matched = np.array(ref_xy_matched)
    return np.array(matched_sci), np.array(matched_ref), sci_xy_matched, ref_xy_matched

def plot_matches(sci_xy, sci_data,ref_xy,ref_data,border_pix,plot=False):
    # ax.plot(sci_xy[:,0], sci_xy[:,1], 'ro', label='Science Stars')
    # ax.plot(ref_xy[:,0], ref_xy[:,1], 'bx', label='Reference Stars')
    axs=None
    if plot:
        cols = plt.cm.tab20.colors
        sci_vmin, sci_vmax = np.percentile(sci_data, [5, 99])
        ref_vmin, ref_vmax = np.percentile(ref_data, [5, 99])
        fig, axs= plt.subplots(1,2,figsize=(14,7))
        axs[0].imshow(sci_data, origin='lower', cmap='gray', vmin=sci_vmin, vmax=sci_vmax)
        axs[0].plot(border_pix[:,0], border_pix[:,1], 'r-', lw=2, label='Science border')
        for i, (x, y) in enumerate(sci_xy):
            # print(x, y)
            circle = plt.Circle((x, y), 6, color=cols[i%len(cols)], fill=False, lw=2.5)
            axs[0].add_artist(circle)
            axs[0].text(x+3, y+3, f"{i+1}", color=cols[i%len(cols)], fontsize=12)

        axs[1].imshow(ref_data, origin='lower', cmap='gray', vmin=ref_vmin, vmax=ref_vmax)
        axs[1].plot(border_pix[:,0], border_pix[:,1], 'r-', lw=2, label='Science border')
        for i, (x, y) in enumerate(ref_xy):
            circle = plt.Circle((x, y), 6, color=cols[i%len(cols)], fill=False, lw=2.5)
            axs[1].add_artist(circle)
            axs[1].text(x+3, y+3, f"{i+1}", color=cols[i%len(cols)], fontsize=12)

        # axs[0].set_title("Matched Stars Between Science and Reference")
        # axs[1].set_title("Matched Stars Between Science and Reference")

        hx, hy = hull.exterior.xy
        axs[0].plot(hx, hy, 'b--', lw=2, label='Matched-star hull')

        # Intersection (filled)
        if not intersection.is_empty:
            ix, iy = intersection.exterior.xy
            axs[0].fill(ix, iy, color='cyan', alpha=0.3, label='Covered area')

        # axs.legend()
        # plt.show()
        plt.close()
    cov,hull,intersection = get_cover(border_pix, sci_xy)
    return axs,cov



def transform_stars(
    sci_xy,
    ref_xy,
    sci_data,
    ref_data,
    order=1,
    transform_type="poly",
    print_errors=True
):
    """
    Transform science image into reference frame using matched stars.

    Parameters
    ----------
    sci_xy : (N,2) array
        Science star pixel positions
    ref_xy : (N,2) array
        Reference star pixel positions
    sci_data : 2D array
        Science image
    ref_data : 2D array
        Reference image
    order : int
        Polynomial order (poly only)
    transform_type : str
        'poly', 'affine', 'similarity', 'projective'
    print_errors : bool
        Print per-star residuals

    Returns
    -------
    corrected_sci : 2D array
        Science image warped into reference frame
    trans : transform object
    """

    from skimage.transform import (
        warp,
        PolynomialTransform,
        AffineTransform,
        SimilarityTransform,
        ProjectiveTransform,
    )

    sci_xy = np.asarray(sci_xy)
    ref_xy = np.asarray(ref_xy)

    # ----------------------------
    # Select transform
    # ----------------------------
    if transform_type == "poly":
        trans = PolynomialTransform()
        success = trans.estimate(sci_xy,ref_xy, order=order)

    elif transform_type == "affine":
        trans = AffineTransform()
        success = trans.estimate(sci_xy, ref_xy)

    elif transform_type == "similarity":
        trans = SimilarityTransform()
        success = trans.estimate(sci_xy, ref_xy)

    elif transform_type == "projective":
        trans = ProjectiveTransform()
        success = trans.estimate(sci_xy, ref_xy)

    else:
        raise ValueError(f"Unknown transform_type: {transform_type}")

    if not success:
        raise RuntimeError(f"{transform_type} transform estimation failed")

    # ----------------------------
    # Warp science → reference
    # ----------------------------

    inverse = trans.inverse if transform_type != "poly" else trans
    
    corrected_sci = warp(
        sci_data,
        inverse_map=inverse,
        output_shape=ref_data.shape,
        preserve_range=True,
    )

    # ----------------------------
    # Print diagnostics
    # ----------------------------
    if print_errors:
        print(info_g+f"=== {transform_type.upper()} TRANSFORM ===")
        print(info_g+" idx | sci(x,y) → trans(x,y) | ref(x,y) | residual (pix)")
        # print(INFO_"-" * 72)

        total_err = []

        for i, (sci_pt, ref_pt) in enumerate(zip(sci_xy, ref_xy)):
            x_s, y_s = sci_pt
            x_r, y_r = ref_pt

            x_t, y_t = trans(np.array([[x_s, y_s]]))[0]

            err = np.hypot(x_t - x_r, y_t - y_r)
            total_err.append(err)

            print(
                f"{i:3d} | "
                f"({x_s:7.2f},{y_s:7.2f}) → "
                f"({x_t:7.2f},{y_t:7.2f}) | "
                f"({x_r:7.2f},{y_r:7.2f}) | "
                f"{err:6.3f}"
            )


        total_err = np.array(total_err)
        # print("-" * 72)
        print(info_g+
            f"Mean residual: {total_err.mean():.3f} px | "
            f"Median: {np.median(total_err):.3f} px | "
            f"Max: {total_err.max():.3f} px"
        )

    return corrected_sci, trans, total_err

def plot_correct(corrected_sci,ref_data,sci_border_pix,sci_xy_matched,ref_xy_matched,trans):
    # ----------------------------
    # Plot corrected science image
    # ----------------------------
    fig,axs = plt.subplots(1,2,figsize=(14,7))
    sci_vmin, sci_vmax = np.percentile(corrected_sci, [5, 99])
    ref_vmin, ref_vmax = np.percentile(ref_data, [5, 99])
    axs[0].imshow(corrected_sci, origin='lower', cmap='gray', vmin=sci_vmin, vmax=sci_vmax)
    cols = plt.cm.tab20.colors
    axs[0].plot(sci_border_pix[:,0], sci_border_pix[:,1], 'r-', lw=2, label='Science border')
    axs[1].plot(sci_border_pix[:,0], sci_border_pix[:,1], 'r-', lw=2, label='Science border')
    for i, (x, y) in enumerate(ref_xy_matched):
        circle = plt.Circle((x, y), 6, color=cols[i%len(cols)], fill=False, lw=2.5)
        axs[1].add_artist(circle)
        axs[1].text(x+3, y+3, f"{i+1}", color=cols[i%len(cols)], fontsize=12)
    for i, (x,y) in enumerate(sci_xy_matched):
        # circle = plt.Circle((x, y), 6, color=cols[i%len(cols)], fill=False, lw=2.5)
        # axs[0].add_artist(circle)
        #plot an arrow showing movement from sci to transformed position
        x_t,y_t = trans(np.array([[x, y]]))[0]
        axs[0].arrow(x, y, x_t - x, y_t - y, color=cols[i%len(cols)], head_width=2, head_length=2)
        circle = plt.Circle((x_t, y_t), 6, color=cols[i%len(cols)], fill=False, lw=2.5, ls='--')
        axs[0].add_artist(circle)
        axs[0].text(x_t+3, y_t+3, f"{i+1}", color=cols[i%len(cols)], fontsize=12)

    axs[1].imshow(ref_data, origin='lower', cmap='gray', vmin=ref_vmin, vmax=ref_vmax)
    # axs[1].set_title('Transformed Science Image')
    plt.show()

    return




def correct_distrotion(science_fits,reference_fits,):
    print('Finding stars in science and reference images...')
    sci_stars,ref_stars = get_stars(science_fits,reference_fits)
    (sci_xy, sci_flux, sci_size, 
        sci_border_pix, sci_data, sci_radec) = sci_stars

    (ref_xy, ref_flux, ref_size, 
        ref_border_pix, ref_data, ref_radec) = ref_stars 

    # figs,axs = plt.subplots(1,2,figsize=(14,7))
    # axs[0]=plot_stars(sci_xy, sci_flux, sci_size, sci_border_pix, sci_data,axs[0])
    # axs[1]=plot_stars(ref_xy, ref_flux, ref_size, ref_border_pix, ref_data,axs[1])
    # plt.show()

    match_tries = 0
    print('Matching stars between science and reference images...')
    [print(sci_radec[o]) for o in range(len(sci_radec))]
    [print(ref_radec[o]) for o in range(len(ref_radec))]

    while match_tries<=7:
        match_sci, match_ref, sci_xy_matched, ref_xy_matched = match_stars(sci_radec, ref_radec, sci_xy, ref_xy, sci_border_pix, max_sep_arcsec=3+match_tries)
        # match_sci, match_ref = match_stars(sci_radec, ref_radec, max_sep_arcsec=3+match_tries)
        cov = get_cover(sci_border_pix, sci_xy_matched)[0]
        print(f"Sep: {3+match_tries}, Matched stars: {len(match_sci)}, Coverage: {cov:.3f}")
        if cov>=0.5:
            break
        match_tries+=1


    axs,cov=plot_matches(sci_xy_matched, sci_data, ref_xy_matched, ref_data, sci_border_pix)

    # print(sci_xy_matched,ref_xy_matched
    # get_cover(sci_border_pix, sci_xy_matched)


    errs = {}
    for mode in ['poly']:#'affine', 'similarity', 'projective', ]:
        for order in [1,2]:
            if mode !='poly' and order!=1:
                continue
            print(f"\n\n=== Testing {mode.upper()} Transform ===")
            corrected_sci,trans,err = transform_stars(sci_xy_matched, ref_xy_matched, sci_data, ref_data, order=order, transform_type=mode)
            errs[(mode, order)] = err,corrected_sci,trans

            # plot_correct(corrected_sci,ref_data,sci_border_pix,sci_xy_matched,ref_xy_matched,trans)


    #choose best transform
    mean_errs = {k: np.mean(v[0]) for k,v in errs.items()}
    best_transform = min(mean_errs, key=mean_errs.get)
    print(f"Best transform: {best_transform[0].upper()} order {best_transform[1]} with mean residual {mean_errs[best_transform]:.3f} px")

    return errs[best_transform][1], errs[best_transform][2]

# science_fits = "/Users/kryanhinds/sedm_phot/aligned_images/ZTF25acchxhv_g2025-12-28T13_49728bkgsub_padded.fits"
# reference_fits = "/Users/kryanhinds/sedm_phot/aligned_images/stack_g_ra247.269893_dec2.509080_arcsec366_skycell1384.066.resamp.fits"
# reference_catalog_file ='/Users/kryanhinds/sedm_phot/ps_catalogs/ps_247.269893_2.50908_0.092553.xml'

# raw = '/Users/kryanhinds/sedm_phot/./data/ZTF25acchxhv/rc20251228_13_51_10_f_b_ZTF25acchxhv_r_r.fits'
# bkg = '/Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF25acchxhv_r2025-12-28T13_49870bkgsub.fits'
# science_fits = "/Users/kryanhinds/sedm_phot/aligned_images/ZTF25acchxhv_r2025-12-28T13_49870bkgsub_padded.fits"
# reference_fits = "/Users/kryanhinds/sedm_phot/aligned_images/stack_r_ra247.269893_dec2.509080_arcsec419_skycell1384.066.resamp.fits"

def plot_img(img):
    img=fits.getdata(img) if isinstance(img,str) else img
    plt.imshow(img, origin='lower', cmap='gray', vmin=np.percentile(img,5), vmax=np.percentile(img,95))
    plt.show()

def find_border_high(sci_hdu):
    border_xl,border_xr = [],[]
    border_yl,border_yt = [],[]

    pos = np.arange(250,850,50)

    for Y in pos:
        middle_med_x = np.nanmedian(sci_hdu.data[Y][250:750])
        break_left, break_right = False, False
        for i in range(len(sci_hdu.data[Y])-5,0,-1):
            j = len(sci_hdu.data[Y])-i
            if all(x>middle_med_x*1.2 for x in sci_hdu.data[Y][i:i+5]) and break_left ==False and i< len(sci_hdu.data[Y])/2:
                border_xl.append(i)
                break_left = True

            if all(x>middle_med_x*1.2 for x in sci_hdu.data[Y][j:j+5]) and break_right ==False and j>len(sci_hdu.data[Y])/2:
                border_xr.append(j)
                break_right = True
            
            if break_left and break_right:
                break
                
    if len(border_xl)==0:border_xl = 0
    else:border_xl = int(np.median(border_xl))+10

    if len(border_xr)==0:border_xr = np.shape(sci_hdu.data)[1]
    else:border_xr = int(np.median(border_xr))-10
    if border_xr<850: border_xr = 950
    if border_xl>150: border_xl = 50

    print(f"Left X border at {border_xl}")
    print(f'Right X border at {border_xr}')


    for X in pos:
        break_top, break_bottom = False, False
        middle_med_y = np.nanmedian(sci_hdu.data[250:750,X])
        for i in range(len(sci_hdu.data[:,X])-5,0,-1):
            j = len(sci_hdu.data[:,X])-i
            if all(x>middle_med_y*1.2 for x in sci_hdu.data[i:i+5,X]) and break_bottom ==False and i< len(sci_hdu.data[:,X])/2:
                border_yl.append(i)
                break_bottom = True

            if all(x>middle_med_y*1.2 for x in sci_hdu.data[j:j+5,X]) and break_top ==False and j>len(sci_hdu.data[:,X])/2:
                border_yt.append(j)
                break_top = True
                
            if break_bottom and break_top:
                break
    if len(border_yl)==0:border_yl = 0
    else:border_yl = int(np.median(border_yl))+10

    if len(border_yt)==0:border_yt = np.shape(sci_hdu.data)[0]
    else:border_yt = int(np.median(border_yt))-10

    if border_yt<850: border_yt = 850
    if border_yl>150: border_yl = 150
    print(f"Bottom Y border at {border_yl}")
    print(f'Top Y border at {border_yt}')

    return int(border_xl),int(border_xr),int(border_yl),int(border_yt)

def find_border_0(sci_hdu):
    print(info_g+"Finding borders using gradient method...")
    inds_every_10_x = np.arange(10, np.shape(sci_hdu.data)[1], 10)
    inds_every_10_y = np.arange(10, np.shape(sci_hdu.data)[0], 10)
    grads = {'xl':{},'xr':{},'yl':{}, 'yt':{}}
    pos = np.arange(250,850,50)
    for Y in pos:
        gl,gr = [],[]
        for i in inds_every_10_x:
            gl.append([int(i)-5, (sci_hdu.data[int(Y),int(i)]-sci_hdu.data[int(Y),int(i)-10])/10])
            # print(np.shape(sci_hdu.data)[1] - i+5)
            try:gr.append([np.shape(sci_hdu.data)[1] - i+5, (sci_hdu.data[int(Y),int(i)]-sci_hdu.data[int(Y),int(i)+10])/10])
            except:pass
        gl,gr = np.array(gl),np.array(gr)
        grads['xl'][Y] = gl[np.where(gl[:,1]!=0)[0][0],0]
        grads['xr'][Y] = gr[np.where(gr[:,1]!=0)[0][0],0]
        # grads['xr'][Y] = 
        # grads = np.array(grads)

    for X in pos:
        gl,gt = [],[]
        for i in inds_every_10_y:
            gl.append([int(i)-5, (sci_hdu.data[int(i),int(X)]-sci_hdu.data[int(i)-10,int(X)])/10])
            try:gt.append([np.shape(sci_hdu.data)[0]-int(i)+5, (sci_hdu.data[int(i),int(X)]-sci_hdu.data[int(i)+10,int(X)])/10])
            except:pass
        gl,gt = np.array(gl),np.array(gt)

        grads['yl'][X] = gl[np.where(gl[:,1]!=0)[0][0],0]
        grads['yt'][X] = gt[np.where(gt[:,1]!=0)[0][0],0]

    xl,xr = np.median(list(grads['xl'].values()))+5, np.median(list(grads['xr'].values()))-5
    print(f"X left border at {xl}, X right border at {xr}")
    yl,yt = np.median(list(grads['yl'].values()))+5, np.median(list(grads['yt'].values()))-5
    print(f"Y bottom border at {yl}, Y top border at {yt}")

    return xl,xr,yl,yt

def get_borders(sci_hdu):
    print(info_g+"Finding borders using two methods...")
    xl_1,xr_1,yl_1,yt_1 = find_border_0(sci_hdu)
    print(info_g+"Method 1 done ")
    xl_2,xr_2,yl_2,yt_2 = find_border_high(sci_hdu)
    print(info_g+"Method 2 done ")

    XL = np.max([xl_1, xl_2])
    XR = np.min([xr_1, xr_2])
    YL = np.max([yl_1, yl_2])
    YT = np.min([yt_1, yt_2])

    return int(XL), int(XR), int(YL), int(YT)

def new_cutout(sci_hdu):
    print(info_g+"Cutting out the image using the borders...")
    XL, XR, YL, YT = get_borders(sci_hdu)
    data = sci_hdu.data.astype(float)
    header = sci_hdu.header.copy()
    print(f"Cutting out image with borders: XL={XL}, XR={XR}, YL={YL}, YT={YT}")
    blank = np.zeros_like(data)
    blank[int(YL):int(YT), int(XL):int(XR)] = data[int(YL):int(YT), int(XL):int(XR)]

    trimmed_data = np.nan_to_num(blank)



    return trimmed_data
