#!/home/arikhind/miniconda3/envs/ltsub/bin python3
'''#!/Users/kryanhinds/opt/miniconda/envs/ltsubphot/bin python3'''

import numpy as np
#import pyfits
import astropy.io.fits as pyfits
import os, sys
import glob, math
from math import sin, cos, tan, asin, sqrt, pi
import time, datetime
import string
import urllib
from subphot_functions import *
from subphot_credentials import *
from subphot_plot_formating import *
import random 
from astropy.io import ascii
rand_nums_string0 = ''.join(random.choice(string.digits) for i in range(5))
rand_nums_string1 = ''.join(random.choice(string.digits) for i in range(5))
files_to_remove = []
def writeparfile(path):
    params = '''X_IMAGE
Y_IMAGE
ALPHA_J2000
DELTA_J2000
MAG_AUTO
MAGERR_AUTO
ELLIPTICITY
FWHM_IMAGE
FLAGS'''
    pf = open(path+'config_files/align_temp.param','w')
    pf.write(params)
    pf.close()



def writeconfigfile(satlevel=55000.):
    configs='''
#-------------------------------- Catalog ------------------------------------
 
CATALOG_NAME     '''+path+'''config_files/align_temp.cat       # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
PARAMETERS_NAME  '''+path+'''config_files/align_temp.param     # name of the file containing catalog contents
 
#------------------------------- Extraction ----------------------------------
 
DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   5              # minimum number of pixels above threshold
DETECT_THRESH    3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
 
FILTER           Y              # apply filter for detection (Y or N)?
FILTER_NAME      '''+path+'''config_files/align_sex.conv       # name of the file containing the filter
 
DEBLEND_NTHRESH  16             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.02           # Minimum contrast parameter for deblending
 
CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency
 
MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                # NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
 
PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                # <min_radius>
 
 
 
MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             0.0            # detector gain in e-/ADU
PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
 
SEEING_FWHM      1.2            # stellar FWHM in arcsec
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename
 
#------------------------------ Background -----------------------------------
 
BACK_SIZE        64             # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>
 
BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL
 
#------------------------------ Check Image ----------------------------------
 
CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  check.fits     # Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
 
MEMORY_OBJSTACK  3000           # number of objects in stack
MEMORY_PIXSTACK  300000         # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
 
VERBOSE_TYPE     QUIET          # can be QUIET, NORMAL or FULL
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output
'''
    #SATUR_LEVEL      '''+str(satlevel)+'''        # level (in ADUs) at which arises saturation
    pf = open(path+'config_files/align_sex.config','w')
    pf.write(configs)
    pf.close()

    convol='''CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
'''

    cf = open(path+'config_files/align_sex.conv','w')
    cf.write(convol)
    cf.close()



class Obj:
    ra = 0.0
    dec = 0.0
    mag = 0.0
    
    ra_rad = 0.0
    dec_rad = 0.0
    
    def __init__(self, inra, indec, inmag):
        self.ra = inra
        self.dec = indec
        self.ra_rad = inra * math.pi/180
        self.dec_rad = indec * math.pi/180
        self.mag = inmag
    
    def rotate(self, dpa_deg, ra0, dec0):
        dpa_rad = dpa_deg * math.pi/180
        sindpa = sin(dpa_rad)
        cosdpa = cos(dpa_rad)
        rascale = cos(dec0*math.pi/180)
        
        #this is only valid for small fields away from the pole.
        x = (self.ra  - ra0 ) * rascale
        y = (self.dec - dec0)
        
        xrot = cosdpa * x - sindpa * y
        yrot = sindpa * x + cosdpa * y
        
        self.ra   = (xrot / rascale) + ra0
        self.dec  =  yrot + dec0 
        self.ra_rad  = self.ra  * math.pi/180
        self.dec_rad =  self.dec * math.pi/180

class SexObj(Obj):
    x = 0.
    y = 0.
    mag = 0.0
    magerr = 0.0
    ellip = 0.0
    fwhm = 0.0
    flag = 0
    
    def __init__(self, inline):
        inlinearg = inline.split()

        self.x = float(inlinearg[0])
        self.y = float(inlinearg[1])
        self.ra = float(inlinearg[2])
        self.dec = float(inlinearg[3])
        self.mag = float(inlinearg[4])
        self.magerr = float(inlinearg[5])
        self.ellip = float(inlinearg[6])
        self.fwhm = float(inlinearg[7])
        self.flag = int(inlinearg[8])

        self.ra_rad = self.ra * math.pi/180
        self.dec_rad = self.dec * math.pi/180
        
        self.row = [self.x, self.y, self.ra, self.dec, self.mag, self.magerr, self.ellip, self.fwhm, self.flag, self.ra_rad, self.dec_rad]



#Pixel distance
def imdistance(obj1, obj2):
    return ((obj1.x - obj2.x)**2 + (obj1.y - obj2.y)**2)**0.5

#Great circle distance between two points.
def distance(obj1, obj2):
    # both must be Obj's.
    
    ddec = obj2.dec_rad - obj1.dec_rad
    dra  = obj2.ra_rad - obj1.ra_rad
    dist_rad = 2 * asin(sqrt( (sin(ddec/2.))**2 + cos(obj1.dec_rad) * cos(obj2.dec_rad) * (sin(dra/2.))**2))

    dist_deg = dist_rad * 180. / math.pi
    dist_sec = dist_deg * 3600.
    return dist_sec

#Non-great-circle distance is much faster
def quickdistance(obj1, obj2, cosdec):
    ddec = obj2.dec - obj1.dec
    dra  = obj2.ra  - obj1.ra
    if dra > 180: dra = 360 - dra
    return 3600 * sqrt(ddec**2 + (cosdec*dra)**2)


# #Compare objects using magnitude.
# def magcomp(obj1, obj2): #useful for sorting
#     return (obj1.mag > obj2.mag) - (obj1.mag < obj2.mag)

#Compare objects using magnitude.
def magcomp(obj): #useful for sorting; Altered by KD for compatibility with python 3
    return obj.mag
    #return (obj1.mag > obj2.mag) - (obj1.mag < obj2.mag)

#Check if two values are the same to within a fraction specified.
def fuzzyequal(v1, v2, tolerance):
   return abs(v1/v2 - 1) < tolerance
        
def median(l):
   a = np.array(l)
   return np.median(a)

def stdev(l):
   a = np.array(l)
   return np.std(a)

def mode(l):
   if len(l) == 0: return
   s = np.array(sorted(l))
   d = s[1:] - s[:-1]
   nd = len(d)
   if   nd >= 32: g = nd/16
   elif nd >= 6: g = 2
   else:         g = 1
   #g = max(nd / 16,1)  #sensitive to clusters up to a little less than 1/16 of the data set
   minmean = d.sum()
   imean = nd / 2
   for i in range(nd):
       r = [max(i-g,0),min(i+g,nd)]
       m = d[int(r[0]):int(r[1])].mean()
       if m < minmean: 
          minmean = m
          imean = i
   mode = s[imean] #+ s[imean+1])/2
   return mode

def rasex2deg(rastr):
    rastr = str(rastr).strip()
    ra=rastr.split(':')
    if len(ra) == 1: return float(rastr)
    return 15*(float(ra[0])+float(ra[1])/60.0+float(ra[2])/3600.0)
    
def decsex2deg(decstr):
    decstr = str(decstr).strip()
    dec=decstr.split(':')
    if len(dec) == 1: return float(decstr)
    sign=1
    if (decstr[0] == '-'): sign=-1
    return sign*(abs(float(dec[0]))+float(dec[1])/60.0+float(dec[2])/3600.0)
    
def ffunc(u, v, h):
  (A20, A02, A11, A21, A12, A30, A03) = (h.get('A_2_0',0),h.get('A_0_2',0),h.get('A_1_1',0),h.get('A_2_1',0),h.get('A_1_2',0),h.get('A_3_0',0),h.get('A_0_3',0))
  return A20*u**2  + A02*v**2  + A11*u*v  + A21*u**2*v  + A12*u*v**2  + A30*u**3  + A03*v**3

def gfunc(u, v, h):
  (B20, B02, B11, B21, B12, B30, B03) = (h.get('B_2_0',0),h.get('B_0_2',0),h.get('B_1_1',0),h.get('B_2_1',0),h.get('B_1_2',0),h.get('B_3_0',0),h.get('B_0_3',0))
  return B20*u**2  + B02*v**2  + B11*u*v  + B21*u**2*v  + B12*u*v**2  + B30*u**3  + B03*v**3

def finvfunc(u, v, h):
  (AP20, AP02, AP11, AP21, AP12, AP30, AP03) = (h.get('AP_2_0',0),h.get('AP_0_2',0),h.get('AP_1_1',0),h.get('AP_2_1',0),h.get('AP_1_2',0),h.get('AP_3_0',0),h.get('AP_0_3',0))
  return AP20*u**2 + AP02*v**2 + AP11*u*v + AP21*u**2*v + AP12*u*v**2 + AP30*u**3 + AP03*v**3

def ginvfunc(u, v, h):
  (BP20, BP02, BP11, BP21, BP12, BP30, BP03) = (h.get('BP_2_0',0),h.get('BP_0_2',0),h.get('BP_1_1',0),h.get('BP_2_1',0),h.get('BP_1_2',0),h.get('BP_3_0',0),h.get('BP_0_3',0))
  return BP20*u**2 + BP02*v**2 + BP11*u*v + BP21*u**2*v + BP12*u*v**2 + BP30*u**3 + BP03*v**3

def pix2wcs_sip2(x, y, h):
  u = x - h['CRPIX1']
  v = y - h['CRPIX2']
  
  up = u + ffunc(u,v,h)
  vp = v + gfunc(u,v,h)
  
  drax = h['CD1_1'] * up + h['CD1_2'] * vp
  ddec = h['CD2_1'] * up + h['CD2_2'] * vp

  cosdec = math.cos(math.pi * h['CRVAL2'] / 180.)
  ra  = h['CRVAL1'] + drax/cosdec
  dec = h['CRVAL2'] + ddec

  return (ra, dec)

#Calculate the (spherical) position angle between two objects.
def posangle(obj1, obj2):
    
    dra  = obj2.ra_rad - obj1.ra_rad
    pa_rad = numpy.arctan2(cos(obj1.dec_rad)*tan(obj2.dec_rad)-sin(obj1.dec_rad)*cos(dra), sin(dra));
    pa_deg = pa_rad * 180./math.pi;
    pa_deg = 90. - pa_deg  #defined as degrees east of north
    while pa_deg > 200: pa_deg -= 360.   # make single-valued
    while pa_deg < -160: pa_deg += 360.  # note there is a crossing point at PA=200, images at this exact PA
    return pa_deg                        # will have the number of matches cut by half at each comparison level



def sextract(sexfilename, nxpix, nypix, border=3, corner=12, minfwhm=1.5, maxfwhm=25, maxellip=0.5, saturation=-1,delete=False):
    if 'ref_imgs' in sexfilename:rand_nums_string = rand_nums_string0
    else:rand_nums_string = rand_nums_string1
    if maxellip == -1: maxellip = 0.5
    if saturation > 0: 
       sexsaturation = saturation
    else:
       sexsaturation = 64000

    if nxpix <= 0: nxpix = 10000
    if nypix <= 0: nypix = 10000

    if os.path.isfile(path+'config_files/temp.cat'):
       os.system('rm -f '+path+'config_files/temp.cat')

    im = pyfits.open(sexfilename) # get the header to check for SIP distortions
    try:h = im[0].header # could be in the extension, of course
    except:im[1].header
    siporder = max(h.get('A_ORDER',0), h.get('B_ORDER',0))
    try:bkg=np.nanpercentile(im[0].data,50)
    except:bkg=np.nanpercentile(im[1].data,50)
    im.close()
    
    

    try:
       # Sextract the image !
       if not os.path.exists(path+'config_files/align_temp.param'):
           writeparfile(path)
       if not os.path.exists(path+'config_files/align_sex.config'):
           writeconfigfile(saturation)
       
       os.system("sex " + sexfilename + " -c "+path+"config_files/align_sex.config -SATUR_LEVEL "+str(sexsaturation)+' -BACK_TYPE MANUAL -BACK_VALUE '+str(bkg) +f" -CATALOG_NAME "+data1_path+f"config_files/prepsfex_{rand_nums_string}.cat ")
       files_to_remove.append(data1_path+f"config_files/prepsfex_{rand_nums_string}.cat")
       #print("sex " + sexfilename + " -c "+path+"config_files/align_sex.config -SATUR_LEVEL "+str(sexsaturation)+' -BACK_TYPE MANUAL -BACK_VALUE '+str(bkg) +f" -CATALOG_NAME "+path+f"config_files/prepsfex_{rand_nums_string}.cat ")
    except Exception as e:
       print(warn_r+' Error: Problem running sextractor',e)
       print(warn_r+' Check that program is installed and runs at command line using ' + sex_path)
       sys.exit(1)

    # Read in the sextractor catalog
    try:
       cat = open(data1_path+f"config_files/prepsfex_{rand_nums_string}.cat",'r')
       #cat = open(path+"config_files/temp.cat",'r')
       catlines = cat.readlines()
    #    print(len(catlines), 'objects detected in image', sexfilename)
    #    print(catlines)


       cat.close()
    except Exception as e:
       print(warn_r+' Cannot load sextractor output file!',e)
       sys.exit(1)

    if len(catlines) == 0:
       print(warn_r+'Sextractor catalog is empty: try a different catalog?')
       sys.exit(1)

    minx = border
    miny = border
    maxx = nxpix - border    # This should be generalized
    maxy = nypix - border
    
    
    l = -1
    nsexinit = 0
    nsexpass = 0
    xlist = []
    ylist = []
    sexlist = []
    fwhmlist = []
    elliplist = []
    flaglist = []
    while l < len(catlines)-1:
        l += 1
        if (len(catlines[l]) <= 1 or catlines[l][0] == '#'):
            continue
        
        iobj = SexObj(catlines[l]) #process the line into an object
        if siporder == 2: (iobj.ra, iobj.dec) = pix2wcs_sip2(iobj.x, iobj.y, h)
        nsexinit += 1
        
        #Initial filtering
        if iobj.ellip > maxellip : continue
        if iobj.fwhm < minfwhm: continue
        if iobj.fwhm > maxfwhm: continue
        if iobj.x < minx: continue
        if iobj.y < miny: continue
        if iobj.x > maxx: continue
        if iobj.y > maxy: continue
        if iobj.x + iobj.y < corner: continue
        if iobj.x + (nypix-iobj.y) < corner: continue
        if (nxpix-iobj.x) < corner: continue
        if (nxpix-iobj.x) + (nypix-iobj.y) < corner: continue
        if saturation > 0:
           if iobj.flag > 0: continue  # this will likely overdo it for very deep fields.
        
        sexlist.append(iobj)
        xlist.append(iobj.x)
        ylist.append(iobj.y)
        fwhmlist.append(iobj.fwhm)
        elliplist.append(iobj.ellip)
        flaglist.append(iobj.flag)
        nsexpass += 1

    #print nsexinit, 'raw sextractor detections'
    #print nsexpass, 'pass initial critiera'

     # Remove detections along bad columns

    threshprob = 0.0001
    ctbadcol = 0
    for i in range(5):
        txp = 1.0
        xthresh = 1
        while txp > threshprob: 
          txp *= min((len(sexlist)*1.0/nxpix),0.8) # some strange way of estimating the threshold.
          xthresh += 1                          #what I really want is a general analytic expression for
        removelist = []                         #the 99.99% prob. threshold for value of n for >=n out 
        modex = mode(xlist)                     #of N total sources to land in the same bin (of NX total bins)
        for j in range(len(sexlist)):
           if (sexlist[j].x > modex-1) and (sexlist[j].x < modex+1):
             removelist.append(j)
        removelist.reverse()
        if len(removelist) > xthresh:
        #  print(removelist)
         for k in removelist:
           del xlist[k]
           del ylist[k]
           del sexlist[k]
           del fwhmlist[k]
           del elliplist[k]
           del flaglist[k]
           ctbadcol += 1

        typ = 1.0
        ythresh = 1
        while typ > threshprob: 
          typ *= min((len(sexlist)*1.0/nypix),0.8)
          ythresh += 1
        removelist = []
        modey = mode(ylist)
        for j in range(len(sexlist)):
           if (sexlist[j].y > modey-1) and (sexlist[j].y < modey+1):
             removelist.append(j)
        removelist.reverse()
        if len(removelist) > ythresh:
        #  print(removelist)
         for k in removelist:
           del xlist[k]
           del ylist[k]
           del sexlist[k]
           del fwhmlist[k]
           del elliplist[k]
           del flaglist[k]
           ctbadcol += 1
    if ctbadcol > 0: print(info_g+' Removed ', ctbadcol, ' detections along bad columns.')

    
    # Remove galaxies and cosmic rays

    if len(fwhmlist) > 5:
       fwhmlist.sort()
       fwhm20 = np.percentile(fwhmlist,20)
       fwhm25 = np.percentile(fwhmlist,25)
       fwhm50 = np.percentile(fwhmlist,50)     #percentile values
       fwhm75 = np.percentile(fwhmlist,75)
       fwhmmode = mode(fwhmlist)
    else:
       fwhmmode = minfwhm
       fwhm20 = minfwhm

    refinedminfwhm = median([0.75*fwhmmode,0.9*fwhm20,minfwhm]) # if CR's are bigger and more common than stars, this is dangerous...
    #print 'Refined min FWHM:', refinedminfwhm, 'pix'

    ngood = 0
    goodsexlist = []
    for sex in sexlist:
       #print(sex.x, sex.y, sex.ra, sex.dec, sex.mag, sex.ellip, sex.fwhm, sex.flag)
       if sex.fwhm > refinedminfwhm: # and sex.ellip < ellipmode +elliptol:
          goodsexlist.append(sex)
          #print ' o',
          ngood += 1
       #print    
    
    #Sort by magnitude
    #goodsexlist.sort(key=magcomp)
    #print('length of science list ',len(goodsexlist))
    goodsexlist.sort(key=lambda x: x.mag)

    #i = 0
    #for sex in goodsexlist:
    #      if i < 1000: print i, sex.mag, sex.fwhm, sex.ellip
    #      i += 1
    # print(len(sexlist), 'objects detected in image ('+ str(len(sexlist)-len(goodsexlist)) +' discarded)',sexfilename)

    if delete==True: os.remove(data1_path+f"config_files/prepsfex_{rand_nums_string}.cat")
    return goodsexlist 

def getcatalog(catalog, ra, dec, boxsize, minmag=8.0, maxmag=-1, maxpm=60., tryc=False):
    # Get catalog from USNO

    if maxmag == -1:
        maxmag = 999 #default (custom catalog)
        if catalog == 'ub2': maxmag = 21.0#19.5
        if catalog == 'sdss': maxmag = 22.0
        if catalog == 'tmc': maxmag = 20.0

    sepchar = ''    
    if (catalog =='ub2' or catalog=='sdss' or catalog=='tmc'):
        usercat = 0
        racolumn = 1
        deccolumn = 2
        magcolumn = 6
        if catalog=='tmc': magcolumn=3
        if catalog=='sdss': 
           racolumn=7
           deccolumn=8
           magcolumn=11
           sepchar=','
        pmracolumn = 10
        pmdeccolumn = 11    
        if catalog!='sdss':
           queryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + catalog +  "&ra=" + str(ra) + "&dec=" + str(dec) + "&system=J2000&rad=" + str(-boxsize) + "&sort=mag&epoch=2000.00000&nstar=6400"
        else:
           # The WCSTOOLS SDSS lookup doesn't function reliably anymore.
           queryurl = 'http://skyserver.sdss.org/dr12/en/tools/search/x_radial.aspx?whichway=equitorial&ra='+str(ra)+'&dec='+str(dec)+'&radius='+str(boxsize/60.)+'&min_u=0&max_u=30&min_g=0&max_g=30&min_r=0&max_r='+str(maxmag)+'&min_i=0&max_i=30&min_z=0&max_z=30&format=csv&limit=5000' #; DR12
        #print queryurl
        cat = urllib.request.urlopen(queryurl)
        catlines = cat.readlines()
        cat.close()
        if tryc==True:
           if len(catlines) < 15: return []  # give up right away and try another catalog

        #if len(catlines) > 6400-20:
        #   print 'WARNING: Reached maximum catalog query size.'
        #   print '         Gaps may be present in the catalog, leading to a poor solution or no solution.'
        #   print '         Decrease the search radius.'
    else:
        usercat = 1
        try:
           cat = open(catalog,'r')
           #print('Reading user catalog ', catalog)
        except Exception as e:
           print(warn_r+' Failed to open user catalog ', catalog)
           print(warn_r+' File not found or invalid online catalog.  Specify ub2, sdss, or tmc.')
           print(warn_r+" "+ e)
           return []
        racolumn = 0
        deccolumn = 1   # defaults
        magcolumn = -1    #  (to override, specify in first line using format #:0,1,2)  
        catlines = cat.readlines()
        cat.close()
    #print '                ', maxmag, maxpm
    l = -1
    catlist = []
    fwhmlist = []

    while l < len(catlines)-1:
        l += 1
        inline = catlines[l].strip()
        if len(inline) <= 2: continue
        if inline[0:2] == '#:':
             inlinearg = inline[2:].split(',')
             racolumn = int(inlinearg[0])-1
             deccolumn = int(inlinearg[1])-1
             if len(inlinearg) > 2: magcolumn = int(inlinearg[2])-1
             continue
        if (inline[0] < '0' or inline[0] > '9') and inline[0]!='.': continue #this may be too overzealous about
        if (inline[1] < '0' or inline[1] > '9') and inline[1]!='.': continue # removing comments...

        if sepchar == '': inlinearg = inline.split()
        if sepchar != '': inlinearg = inline.split(sepchar)
        narg = len(inlinearg)
    
        if inlinearg[racolumn].find(':') == -1:
           ra = float(inlinearg[racolumn])
        else:
           ra = rasex2deg(inlinearg[racolumn])
        if inlinearg[deccolumn].find(':') == -1:
           dec = float(inlinearg[deccolumn])
        else:
           dec = decsex2deg(inlinearg[deccolumn])
        if magcolumn >= 0 and narg > magcolumn: 
            try:
               mag = float(inlinearg[magcolumn])
            except:
               mag = float(inlinearg[magcolumn][0:-2])
        else:
            mag = maxmag
        if usercat == 0 and narg > pmracolumn and narg > pmdeccolumn:
            pmra = float(inlinearg[pmracolumn])
            pmdec = float(inlinearg[pmdeccolumn])
        else:
            pmra = pmdec = 0
        #print
        #print ra, dec, mag,
        #print pmra, pmdec,
        if mag > maxmag: continue #don't believe anything this faint
        if mag < minmag: continue #ignore anything this bright
        if abs(pmra) > maxpm or abs(pmdec) > maxpm: continue
        #print ' OK',
        iobj = Obj(ra, dec, mag) #process the line into an object
        catlist.append(iobj)
        
    #catlist.sort(magcomp)   
    catlist.sort(key=lambda x: x.mag)
    #print('length of catalog list ',len(catlist))
    #print
    return catlist




def matchlists(list1, list2, maxdist=1.0, logfile='', quiet=1):
    
    n1 = len(list1)
    n2 = len(list2)

    if logfile != '':  f = open(logfile, 'w')
 
    totmatch = 0
    radiff = []
    decdiff = []
    #print n1, n2
    for i in range(n1):
       jmaxdist = maxdist
       matchj = -1
       nmatch = 0
       cosdec = cos(list1[i].dec_rad)
       for j in range(n2):
          dist = quickdistance(list1[i], list2[j], cosdec)
          if dist > jmaxdist: continue
          jmaxdist = dist
          matchj = j
          nmatch += 1
       if nmatch == 1:
          dra = list1[i].ra - list2[matchj].ra
          ddec = list1[i].dec - list2[matchj].dec
          dra_arcsec = dra * cosdec * 3600.
          ddec_arcsec = ddec * 3600.
          #if logfile == '':
             #if quiet==0: print('%3i = %3i:  %6.2f %6.2f' % (i, matchj, dra_arcsec, ddec_arcsec))
          #else:
             #f.write('%3i [ %8.3f %8.3f ] = %3i [ %12.8f %12.8f ] : %7.3f %7.3f\n' % (i, list1[i].x, list1[i].y, matchj, list2[matchj].ra, list2[matchj].dec, dra_arcsec, ddec_arcsec))
          radiff.append(dra)
          decdiff.append(ddec)
          totmatch += 1

    if totmatch==0:
       print(warn_y+' No matches.')
       return (0., 0., totmatch)

    avradiff = median(radiff)
    avdecdiff = median(decdiff)

    stdradiff = stdev(radiff)
    stddecdiff = stdev(decdiff)
    rms = np.std(np.sqrt(np.array(radiff)**2 + np.array(decdiff)**2))

    #print('  Average:  %6.2f %6.2f' % (avradiff*cosdec*3600., avdecdiff*3600.))
    #print(' StdDev: %6.2f %6.2f' % (stdradiff*cosdec*3600., stddecdiff*3600.))
    #print(' (%i matches)' % totmatch)

    if logfile != '':  f.close()

    return (avradiff, avdecdiff, totmatch)

def writeregionfile(filename, objlist, color="green",sys=''):
    if sys == '': sys = 'wcs'
    out = open(filename,'w')
    i = -1
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    if sys == 'wcs': 
      out.write('fk5\n')
      for ob in objlist:
        i += 1
        out.write("point(%.7f,%.7f) # point=boxcircle text={%i}\n" % (ob.ra, ob.dec, i))
    if sys == 'img': 
      out.write('image\n')
      for ob in objlist:
        i += 1
        out.write("point(%.3f,%.3f) # point=boxcircle text={%i}\n" % (ob.x, ob.y, i))
    out.close()


def nudge(image, raoffset, decoffset, rms=-1, acc=-1, n=-1, source='', outfile=''):

    im = pyfits.open(image, mode='update')
    im.verify('silentfix')
    next = len(im)
    ext = 0
    if next > 1: ext = 1
        #print 'nudging, ', next
    im[ext].header['crval1'] = im[ext].header['crval1'] + raoffset/cos((3.1415926/180.)*decoffset)
    im[ext].header['crval2'] = im[ext].header['crval2'] + decoffset

    im[ext].header.set('RANUDGE', raoffset*3600.,  comment='RA shift in arcsec')
    im[ext].header.set('DECNUDGE', decoffset*3600., comment='Dec shift in arcsec')
    # print(info_g+f" Nudging X by:")
    if acc > 0:
          im[ext].header.set('RMSNUDGE', rms, comment='Per-star scatter in shift')
    if rms > 0:
          im[ext].header.set('ACCNUDGE', acc, comment='Uncertainty in shift')
    if n > 0:
          im[ext].header.set('NNUDGE', n, comment='Number of stars used in shift')
    if len(source) > 0:
          im[ext].header.set('NUDGESRC', source, comment='Image or catalog aligned against')
 
    if outfile=='' or outfile==image:
          im.flush()
    else:
          im.writeto(outfile, output_verify='silentfix')
        #   print(np.shape(im[0].data))

def main():
    
    files=[]

    if (len(sys.argv)==1):
        usage()
        sys.exit(1)

    i=1
    method = 'hybrid' # align first image to catalog, then subsequent images to first image
    aligncenter = [-1., -1.]
    relsearchrad = 60. # for aligning images against each other
    abssearchrad = 60. # for aligning images against the catalog
    nosave = 0
    quiet = 0
    catalog = ''
    logfile = ''
    prefix = ''
    maxdist = 1.5 #for matching
    saturation = 40000.
    reqmatch = 2

    while (i<len(sys.argv)):
        arg=sys.argv[i]
        isarg=0
        if (arg.find("-r") == 0):
            searchrad = float(sys.argv[i+1])
            relsearchrad = searchrad
            abssearchrad = searchrad
            i+=1
            isarg = 1
        if (arg.find("-x") == 0):
            maxdist = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-a") == 0):
            incra = sys.argv[i+1]
            incdec = sys.argv[i+2]
            if incra.find(':') > 0: aligncenter[0] = rasex2deg(incra)
            if incdec.find(':') > 0: aligncenter[1] = decsex2deg(incdec)
            if incra.find(':') == -1: aligncenter[0] = float(incra)
            if incdec.find(':') == -1: aligncenter[1] = float(incdec)
            i+=2
            isarg = 1
        if arg.find("-m") == 0:
            newmethod = sys.argv[i+1]
            method = ''
            if newmethod.find('rel') >= 0: method = 'relative'
            if newmethod.find('abs') >= 0: method = 'absolute'
            if newmethod.find('hyb') >= 0: method = 'hybrid'
            # "chain" method might be useful sometimes: align each image versus the previous one (e.g. for a depth or wavelength series)
            i += 1
            isarg = 1
            if method == '':
               print(warn_y+" Cannot recognize method ", newmethod)
               print(warn_y+" Choose relative, absolute, or hybrid")
               return
        if (arg.find("-n") == 0):
            nosave = 1
            isarg = 1
        if (arg.find("-q") == 0):
            quiet = 1
            isarg = 1
        if (arg.find("-c") == 0):
            catalog = sys.argv[i+1]
            i+=1
            isarg = 1
        if (arg.find("-l") == 0):
            logfile = sys.argv[i+1]
            i+=1
            isarg = 1
        if (arg.find("-p") == 0):
            prefix = sys.argv[i+1]
            i+=1
            isarg = 1
        if (arg.find("-e") == 0):
            reqmatch = int(sys.argv[i+1])
            i+=1
            isarg = 1
        if (not isarg):
            if (len(files) > 0):
                files+=",%s" % arg
            else:
                files=arg
        i+=1

    if len(files) == 0:
        print(warn_r+' No files selected!')
        return

    filenames=files.split(",")


    if aligncenter[0] <= -1.:
        if quiet==0:
            print(warn_r+' No alignment center specified!')
    im = pyfits.open(filenames[0])
    next = len(im)
    ext = 0
    for ext in range(next):
        try:
         aligncenter[0] = im[ext].header['CRVAL1']
         aligncenter[1] = im[ext].header['CRVAL2']
        # aligncenter[0] = im[ext].data[]
         break
        except:
         pass
         #if quiet==0: 
            #print('Using '+filenames[0]+' rotation center CRVAL1='+str(aligncenter[0])+', CRVAL2='+str(aligncenter[1]))


    

    center = Obj(aligncenter[0], aligncenter[1], 0.0)
    cosdec = cos(center.dec_rad)
    nfiles = len(filenames)


    if method == 'hybrid' or method == 'absolute':
        # Download catalog from the web


        # Load in reference star catalog
        if catalog != '':
            catlist = getcatalog(catalog, aligncenter[0], aligncenter[1], abssearchrad+30, tryc=False)
        # If no catalog specified, try SDSS, and if nothing returns try USNO
        else:
            trycats = ['sdss', 'ub2', 'tmc']  # need to update this.  tdc is really slow.
            for trycat in trycats:
                catlist = getcatalog(trycat, aligncenter[0], aligncenter[1], abssearchrad+30, tryc=True)
                if len(catlist) >= 3:
                    mindist = 999999.
                    cosdec = cos(aligncenter[1]*3.14159/180.)
                    for c in catlist:
                        dist = abs(c.ra - aligncenter[0])*cosdec/3600. + abs(c.dec - aligncenter[1])/3600.
                        if dist < mindist: mindist = dist
                        if mindist > 90.: continue  # sometimes just off an edge of SDSS
                        #print('Using catalog', trycat)
                        catalog = trycat
                        break
            if (catalog == ''):
                #print('No catalog is available.  Check your internet connection.')
                return -1

        writeregionfile('icat.reg',catlist,sys='wcs')


        # Match image(s) to catalog
        #   print(method+' alignment method')
        if method == 'hybrid': frange = [0]
        if method == 'absolute': frange = range(nfiles)
        for f in frange:
            inlist = sextract(filenames[f], 0, 0, 3, 12, maxellip=0.7, saturation=-1)
            for i in range(len(inlist)-1,0,-1):cdist = quickdistance(center, inlist[i], cosdec)
            if cdist > abssearchrad: inlist.remove(inlist[i])
            writeregionfile('initpos.reg',inlist,sys='img')

            # Match and align image against the catalog
            #print(catalog+' catalog', 'vs.', filenames[f])
            (raoff, decoff, nmatch) = matchlists(inlist, catlist, logfile=logfile, maxdist=maxdist)
            # print('  ', nmatch, 'matches')
            # sys.exit(1)
            if nmatch >= reqmatch:
                try:
                    if nosave==0: nudge(filenames[f], -raoff, -decoff, source=catalog, outfile=prefix+filenames[f])
                except:
                    if quiet==0: 
                        print(warn_r+ ' XXX')
                        print(warn_r+ ' FAILURE: Could not change image coordinates - probably could not find the offset')
                        print(warn_r+ ' XXX')
                    else:
                        if nmatch > 0: 
                            # print(poo)
                            print(warn_r+ ' Not enough matches to nudge ('+str(nmatch)+'  found, '+str(reqmatch)+'  required)')

    # Match images to reference (first) image
    if method == 'hybrid' or method == 'relative':
    
        match_attempt,nmatch = 0,0
        while nmatch <= 2:
            print(info_g+f' Attempting to match science image to reference image {match_attempt}')
            if match_attempt ==0: print(info_g+' Using search radius of '+str(relsearchrad)+' pixels')
            if match_attempt > 0: 
                relsearchrad+=30
                print(info_g+' Expanding search radius to '+str(relsearchrad)+' pixels')
            match_attempt += 1
            if match_attempt > 3:
                print(warn_r+' Could not find enough matches to nudge')
                return
            for f in range(nfiles):
            
                inlist = sextract(filenames[f], 0, 0, 3, 12, maxellip=0.7, saturation=-1)  # sextracting again for image 0 slightly inefficient

                if 'ref_img' in  filenames[f]:print(info_g+f' {len(inlist)} objects detected in reference image', filenames[f])
                else:print(info_g+f' {len(inlist)} objects detected in science image', filenames[f])

                #print len(inlist)
                for i in range(len(inlist)-1,-1,-1):
                    cdist = quickdistance(center, inlist[i], cosdec)
                    if cdist > relsearchrad: inlist.remove(inlist[i])

                if f == 0: 
                    reflist = inlist  # don't match image 0, but transfer inlist to reference list
                else:
                    (raoff, decoff, nmatch) = matchlists(inlist, reflist)
                    print(info_g+f' {nmatch} matches found between science image {f} and reference image 0')
                    if nmatch >= reqmatch:
                        try:
                            if nosave==0: nudge(filenames[f], -raoff, -decoff, source=filenames[0], outfile=prefix+filenames[f])
                        except:
                            if quiet==0: print(warn_r+ ' XXX')
                            print(warn_r+ ' Could not change image coordinates - probably could not find the offset')
                            if quiet==0: print(warn_r+ ' XXX')
                    else:
                        if match_attempt<=3:print(warn_y+' Not enough matches to nudge ('+str(nmatch)+'found, '+str(reqmatch)+'required)')
                        if match_attempt>3:
                            print(warn_r+' Not enough matches to nudge ('+str(nmatch)+'found, '+str(reqmatch)+'required), stopping x-y nudge')
                            # return
                            break

    print(info_g+" X-Y nudge complete")
    print(info_g+" Removing temporary files")
    for remove_file in np.unique(files_to_remove):
        if os.path.exists(remove_file):os.remove(remove_file)
    return

######################################################################
# Running as executable
if __name__=='__main__':
    main()

######################################################################


