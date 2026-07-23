
#!/Users/kryanhinds/miniconda/envs/sedm/bin python3
# from drizzle import drizzle,doblot

from skimage import transform
import logging
from sedm_align_quick import sextract,mode
import matplotlib
import matplotlib.dates as mdates
import pandas as pd
import statistics as stats
import numpy as np
import astropy, photutils
from astropy.io import fits, ascii #FITS files handling
import os,sys,requests,glob,scipy,urllib,datetime,warnings
from scipy.signal import convolve as scipy_convolve
import matplotlib.pyplot as plt
from astropy import units as u,wcs,visualization
from astropy.table import Table,Column,vstack
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, fit_wcs_from_points
from astropy.coordinates import SkyCoord,FK5
from astropy.stats import sigma_clipped_stats,sigma_clip,SigmaClip
from astropy.nddata import Cutout2D
from photutils import IRAFStarFinder,CircularAperture,CircularAnnulus,aperture_photometry,SExtractorBackground,Background2D
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch,simple_norm
from sedm_telescopes import *
#from BeautifulSoup import BeautifulSoup python2.7
from bs4 import BeautifulSoup
from image_registration import chi2_shift
from sedm_functions import *
from sedm_credentials import *
from sedm_plot_formating import *
#from pyzapspec.py import *
import astroscrappy
#from lacosmic import lacosmic
# from multiprocessing import Pool
import concurrent.futures
import time,resource
from termcolor import colored
import astroquery
from astroquery import sdss
from matplotlib.colors import LogNorm
import shutil,re,random,string
from datetime import date,timedelta
from photutils.utils import calc_total_error
import astroalign as aa
from reproject import reproject_exact,reproject_interp
# import progressbar
from time import sleep
from tabulate import tabulate
if not sys.warnoptions:
    warnings.simplefilter("ignore")

import matplotlib as mpl
mpl.rc('font', family='serif')



def get_image(args):
    nx, ny,sci_c,sci_filt,log,lnks_done = args
    return sdss_query_image(sci_c.ra.deg+(nx*0.1), 
                                sci_c.dec.deg+(ny*0.1),
                                sci_filt, nx=nx, ny=ny,log=log,lnks_done=lnks_done)
                                

def estimate_seeing(filename):
    sp_logger.info(info_g+f' Estimating the seeing based on the mode of the FWHM of point sources in the field')
    sextracted = sextract(filename,0, 0, 3, 12, maxellip=0.7, saturation=-1,delete=True)
    see = mode([sex.fwhm for sex in sextracted])
    sp_logger.info(info_g+f' Estimated seeing: '+str(see ))
    return see

def check_seeing(ims,s=5,sp_logger=None):
    #function checks the seeing on all images in ims (list) and rejects the image is the seeing >=5,
    #returns a list of images with seeing <5 in the same order and format as ims
    good_seeing=[]
    im_stack_path = '/'.join(ims[0].split('/')[:-1])
    # try:
    # except:pass 
    for i in range(len(ims)):
        if ims[i].endswith('.fits')==False:
            ims[i] = ims[i]+'.fits'
        if ims[i].startswith(im_stack_path)==False:
            ims[i] = im_stack_path+'/'+ims[i]
        if ims[i].startswith(path)==False:
            ims[i] = path+ims[i]

        if i==0:
            if fits.open(ims[0])[0].header['INSTRUME'] in ['HiPERCAM','OSIRIS']:s=7
        # print(im_stack_path,ims[i])
        im = fits.open(ims[i])
        try:
            seeing = im[0].header['L1SEESEC']
        except:
            try:
                seeing = im[0].header['FWHM']
            except:
                try:
                    seeing = im[0].header['ESO TEL AMBI FWHM END']
                except:
                    sp_logger.warning(warn_r+f" Seeing keywords not found in "+str(ims[i])+", estimating based of FWHM mode")
                    seeing=estimate_seeing(ims[i])
        if seeing < s:
            good_seeing.append(ims[i][:-5])
        else:
            sp_logger.warning(warn_y+f' {ims[i]} has seeing of {seeing}, rejecting')
    # sys.exit(1)
    return good_seeing

def quick_app_phot(data,stars,gain=1.62):
    data_sig_clip = SigmaClip(sigma=3.)
    data_bkg_estimator = SExtractorBackground(data_sig_clip)        
    data_bkg = Background2D(data, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=data_bkg_estimator)
    data_bkg_error = data_bkg.background_rms

    data_err = calc_total_error(data,data_bkg_error,gain)
    app_pos = [(stars[i][0],stars[i][1]) for i in range(len(stars))]
    data_photaps = photutils.CircularAperture(app_pos, r=10)
    data_photTab = photutils.aperture_photometry(data, data_photaps, error=data_err)
    data_photTab['SNR'] = data_photTab['aperture_sum']/data_photTab['aperture_sum_err']

    return data_photTab

class subtracted_phot(subphot_data):
    def __init__(self,ims,args):
        self.bright_cat_mag_lim = 16
        self.zp_sci_lim = 30
        self.path=path
        self.path=path
        self.path=path
        self.ims_path = '/'.join(ims[0].split('/')[:-1])
        self.args=args
        self.ims,self.cutout_tf,self.unsubtract = ims,self.args.cutout,self.args.unsubtract
        self.out_dir,self.termoutp,self.stack_tf,self.upfritz,self.cleandirs = self.args.output, self.args.termoutp,self.args.stack,self.args.upfritz,self.args.cleandirs
        self.morning_round_up,self.zp_only,self.auto_ref,self.auto_cat = self.args.mroundup,self.args.zp_only,self.args.ref_img,self.args.ref_cat
        self.relative_flux,self.special_case = self.args.relative_flux,self.args.special_case
        self.use_sdss = self.args.use_sdss
        self.user_zp_sci,self.user_zp_ref = self.args.user_zp_sci,self.args.user_zp_ref
        # print(self.args)
        self.max_psf_offset = 0.5
        self.upfritz_f = self.args.upfritz_f
        self.ra_dec_pos = self.args.position
        self.redo_astrometry = self.args.redo_astrometry
        self.forced_phot = self.args.forced_phot
        self.sci_bkg_med=None
        # print(self.redo_astrometry)
        if len(self.ra_dec_pos)==1:
            self.ra_dec_pos = self.args.position[0].split(' ')


        # print(self.ra_dec_pos[0],self.ra_dec_pos[1])
        # if any(self.ra_dec_pos[0].startswith(J) for J in ['"',"'"]):
        #     self.ra_dec_pos = self.ra_dec_pos[1:-1]
        #     self.ra_dec_pos = self.ra_dec_pos.split(' ')
        #     print(self.ra_dec_pos)

        # sys.exit()
        self.sys_exit=False
        self.matched_star_coords_pix = []

        

        if not os.path.exists(self.path+'phot_fits_info'):os.mkdir(self.path+'phot_fits_info')
        if not os.path.exists(self.path+'combined_imgs'):os.mkdir(self.path+'combined_imgs')
        if not os.path.exists(self.path+'out'):os.mkdir(self.path+'out')
        if not os.path.exists(self.path+'temp_config_files'):os.mkdir(self.path+'temp_config_files')
        if not os.path.exists(self.path+'aligned_images'):os.mkdir(self.path+'aligned_images')
        if not os.path.exists(self.path+'ref_imgs'):os.mkdir(self.path+'ref_imgs')
        if not os.path.exists(self.path+'bkg_subtracted_science'):os.mkdir(self.path+'bkg_subtracted_science')
        if not os.path.exists(self.path+'scaled_subtracted_imgs'):os.mkdir(self.path+'scaled_subtracted_imgs')
        if not os.path.exists(self.path+'convolved_sci'):os.mkdir(self.path+'convolved_sci')
        if not os.path.exists(self.path+'convolved_ref'):os.mkdir(self.path+'convolved_ref')
        if not os.path.exists(self.path+'convolved_psf'):os.mkdir(self.path+'convolved_psf')
        if not os.path.exists(self.path+'trimmed_sci_imgs'):os.mkdir(self.path+'trimmed_sci_imgs')
        if not os.path.exists(self.path+'photometry_data'):os.mkdir(self.path+'photometry_data')
        if not os.path.exists(self.path+'photometry'):os.mkdir(self.path+'photometry')
        if not os.path.exists(self.path+'ps_catalogs'):os.mkdir(self.path+'ps_catalogs')




        

        if self.args.show_plots!=True:
            matplotlib.use('Agg')        
            # print('created the photometry folder')

       
        #setting data to todays data in format YYYYMMDD
        self.t = date.today()
        self.TIME = datetime.datetime.now().strftime("%H:%M:%S")
        self.year,self.month,self.dayy = self.t.strftime("%Y"),self.t.strftime("%m"),self.t.strftime("%d")
        self.today = Time(f'{self.year}-{self.month}-{self.dayy} {self.TIME}')
        self.TODAY = self.t.strftime("%Y%m%d") #todays date in YYYYMMDD format
        self.apo = Observer.at_site("lapalma")
        self.sun_set_today = self.apo.sun_set_time(self.today, which="nearest") #sun set on day of observing
        self.time_suns_today = "{0.iso}".format(self.sun_set_today)[-12:]
        self.sun_set_tomorrow = self.apo.sun_set_time(self.today,which="next")
        self.time_suns_tomorrow = "{0.iso}".format(self.sun_set_tomorrow)[-12:]


        if self.time_suns_today<self.TIME<'23:59:59':
            self.date_ = self.TODAY
            self.DATE = re.sub("-","",self.date_)
        if '00:00:00'<self.TIME<self.time_suns_tomorrow:
            self.date_ = str(self.t - datetime.timedelta(days=1))
            self.DATE=re.sub("-","",self.date_)


        if self.args.sp_logger!=False:self.sp_logger = self.args.sp_logger
        else:
            sp_logger = logging
            # if os.path.exists(log_name):os.remove(log_name)
            print(info_g+f' Logging to {log_name}')
            sp_logger.basicConfig(level=logging.INFO,encoding='utf-8',handlers=[
            logging.StreamHandler(),logging.FileHandler(log_name,encoding='utf-8')],format='%(message)s')
            
            self.sp_logger=sp_logger
        self.t_start = time.time()
        self.to_subtract=True

        if self.unsubtract==True:
            self.to_subtract=False

            
        self.if_stacked=False
        self.MJD,self.files_to_clean = [],[]
        self.rand_nums_string = ''.join(random.choice(string.digits) for i in range(5))
        


        if len(self.ims)<1:
            self.sp_logger.warning(warn_r+' Add image names without fits extension to run eg.:  python sedm_subtract.py -i filename1 filename2')
            ##sys.exit(1)
            self.sys_exit=True

        if len(self.ims)==1:
            try:
                if self.ims[0].startswith('h_')==False:    #checking whether name is of a single objet or the name starts with full path
                    # print(self.ims[1].split('/'))
                    self.arg_split = self.ims[0].split('/')
                    for s in range(len(self.arg_split)):
                        if any(x in self.arg_split[s] for x in ['h_e_','h_s_']):
                            self.folder = re.sub('_1','',self.arg_split[s].split('_')[2])
                            self.name_split = self.arg_split[s].split('_')
                            self.exp_identity = self.name_split[4]
                            self.single_data_path = '/'.join(self.arg_split[:s])

                    self.img_type=""
                    

                    #checking whether this is part of a stack and then checking whether the stacks are available, appending to ims to proceed with stacking else will continue with the one
                    if self.exp_identity!='1' and self.stack_tf==True:
                        attempts = [1,2,3,4,5,6]
                        for j in attempts:
                            fits_name = '_'.join(self.name_split[0:4])+'_'+f'{j}'+'_'+'_'.join(self.name_split[5:])
                            try_name = self.path+f"{self.single_data_path}/"+'_'.join(self.name_split[0:4])+'_'+f'{j}'+'_'+'_'.join(self.name_split[5:])+".fits"
                            if os.path.exists(try_name)==True and re.sub('.fits','',try_name) not in self.ims:
                                self.ims.append(fits_name)
            except:
                pass        
    
            if len(self.ims)!=1:

                self.sp_logger.info(info_g+f' Detected {self.ims[0]} is part of a stack. Found {len(self.ims)-2} other images in the stack, proceeding with the stack of {len(self.ims)-1} images total')
                pass

            else:
                self.img_type=""
                if self.ims[0].endswith('.fits')==False and self.ims[0].endswith('.fits.gz')==False:
                    self.name=self.ims[0]+'.fits'
                else:
                    self.name=self.ims[0]
                
                # if self.name.startswith(
                # self.sp_logger.info(info_g+' Single image to subtract:',self.name)

                if self.args.telescope_facility not in SEDM:
                    if '_' in self.name: 
                        self.folder = re.sub('_1','',self.name.split('/')[-1].split('_')[2])
                    else:
                        self.folder = self.DATE

                if self.path not in self.name and self.name.startswith('/')==False:
                    self.name=self.path+self.name
                
                
                self.sci_path = self.name
                print(self.name)
                self.sci_img_hdu=fits.open(self.name.replace('.fits.fits','.fits'))[0]
                self.telescope = self.sci_img_hdu.header['TELESCOP']

                if self.telescope in SEDM:self.telescope,self.sci_int='SEDM-P60','P60'

                elif self.telescope=='GTC':
                    if self.sci_img_hdu.header['INSTRUME']=='HIPERCAM':self.telescope='GTC-HIPERCAM'
                    else: self.telescope='GTC-OSIRIS'
                    


                # self.sp_logger.info(header_kw)


                if any(t in header_kw.keys() for t in [self.telescope,self.telescope.upper()]):
                    self.sp_logger.info(info_g+' Telescope: \033[1m'+self.telescope+'\033[0m')
                    self.tele_kw = header_kw[self.telescope]
                    
                    self.FILT_kw,self.OBJ_kw,self.RA_kw = self.tele_kw['filter'],self.tele_kw['object'],self.tele_kw['ra']
                    self.DEC_kw,self.AIRM_kw,self.UTS_kw = self.tele_kw['dec'],self.tele_kw['airmass'],self.tele_kw['utstart']
                    self.INST_kw,self.EXPT_kw,self.DATE_kw,self.MJD_kw = self.tele_kw['instrument'],self.tele_kw['exptime'],self.tele_kw['date'],self.tele_kw['mjd']
                    self.PIX_kw,self.GAIN_kw = self.tele_kw['pixscale'],self.tele_kw['gain']
                    self.SEE_kw = self.tele_kw['seeing']
                    self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = None,None,None
                    self.ESEE_kw = None
                    self.sci_prop = None

                    
                    if self.telescope=='Liverpool Telescope':
                        self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = self.tele_kw['propid'],self.tele_kw['bkg_med'],self.tele_kw['bkg_mean']
                        self.ESEE_kw = self.tele_kw['est_seeing']
                        self.sci_filt=self.sci_img_hdu.header['FILTER1'][5].lower()
                        self.sci_inst = 'IO:O'
                        self.sci_prop = self.sci_img_hdu.header[self.PRO_kw]

                    elif self.telescope in ['eso-ntt','ntt-efosc','ntt','efosc','ESO-NTT','NTT-EFOSC','NTT','EFOSC']:
                        self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = '-','-','-'
                        self.ESEE_kw = self.tele_kw['seeing']
                        self.sci_filt=self.sci_img_hdu.header['ESO INS FILT1 NAME'][0]
                        self.sci_inst = 'EFOSC'
                        self.sci_prop = None
                        
                        
                    
                    elif self.telescope=='HCT':
                        self.sci_filt=self.sci_img_hdu.header['FILTER'] #Bessel filters
                        self.sci_inst = 'HFOSC'
                        self.sci_prop = None
                        self.folder = str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))
                        # self.sp_logger.info(self.folder)

                    elif self.telescope in SEDM:
                        self.sp_logger.info(info_g+f" Trying to convert SIP to TPV for {self.sci_path}")
                        try:
                            from sip_tpv import sip_to_pv
                            for cd,pc in zip(['CD1_1','CD1_2','CD2_1','CD2_2'],['PC1_1','PC1_2','PC2_1','PC2_2']):
                                self.sci_img_hdu.header[cd]=self.sci_img_hdu.header[pc]
                                del self.sci_img_hdu.header[pc]

                            sip_to_pv(self.sci_img_hdu.header)
                            self.sci_img_hdu.writeto(self.sci_path.replace('.fits','_tpv.fits'),overwrite=True)
                            self.files_to_clean.append(self.sci_path.replace('.fits','_tpv.fits'))
                            self.sci_path = self.sci_path.replace('.fits','_tpv.fits')
                            self.sci_img_hdu=fits.open(self.sci_path)[0]
                            self.sp_logger.info(info_g+f" Successfully converted SIP to TPV for {self.sci_path}")
                        except Exception as e:
                            self.sp_logger.warning(warn_r+f" Error converting SIP to TPV for {self.sci_path}, error: {e}")


                        self.redo_astrometry = True
                        self.sci_filt=self.sci_img_hdu.header['FILTER']
                        self.sci_inst='SEDM-P60'
                        if self.name.split('/')[-1].startswith('rc'):
                            self.sci_prop=None
                            self.folder=self.sci_path.split('rc')[1][2:10]
                        else:
                            self.sci_prop=None
                            self.folder=self.sci_path.split('/')[-1][0:8]


                        if self.redo_astrometry==True:
                            self.sci_img_hdu=fits.open(self.sci_path)[0]
                            # if 'ACQ' in self.sci_path.split('/')[-1]:
                            self.sp_logger.info(info_g+f" This is an ACQ image, trimming the image to remove the edges")
                            #trim the left and bottom of the image (60 pixels)
                            self.orig_sci_wcs = WCS(self.sci_img_hdu.header)
                            vmin,vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                            height,width = self.sci_img_hdu.data.shape
                            trim_size=0
                            # self.sci_trimmed_img = np.zeros((height,width))*np.nan
                            self.sp_logger.info(info_g+f" Original image size: "+str(self.sci_img_hdu.data.shape))
                            self.sci_trimmed_img= new_cutout(self.sci_img_hdu)
                            # self.sci_trimmed_img.data = np.nan_to_num(self.sci_trimmed_img.data)
                            # if not os.path.exists(self.path+'trimmed_sci_imgs'):os.makedirs(self.path+'trimmed_sci_imgs')
                            self.sci_path = self.path+'trimmed_sci_imgs/'+self.sci_path.split('/')[-1]
                            self.sci_path_o = self.sci_path
                            fits.writeto(self.sci_path.replace('.fits','_trimmed.fits'),self.sci_trimmed_img,overwrite=True,header=self.sci_img_hdu.header)
                            # fits.writeto(self.sci_path.replace('.fits','_trimmed.fits'),self.sci_trimmed_img.data,overwrite=True,header=self.sci_img_hdu.header)
                            self.sci_path = self.sci_path.replace('.fits','_trimmed.fits')
                            self.files_to_clean.append(self.sci_path)
                            # sys.exit()

                            self.sci_img_hdu.data = self.sci_trimmed_img.data
                            # self.sp_logger.info(info_g+' Updating astrometry for',self.sci_path,' using local version of astrometry.net')
                            self.im_size1,self.im_size2 = self.sci_img_hdu.header['NAXIS1'],self.sci_img_hdu.header['NAXIS2']
                            self.sci_ra,self.sci_dec = self.sci_img_hdu.header[self.RA_kw],self.sci_img_hdu.header[self.DEC_kw]
                            self.sci_c=SkyCoord(self.sci_ra,self.sci_dec,unit=(u.hourangle, u.deg),frame='fk5')
                            self.sci_ra_d,self.sci_dec_d = self.sci_c.ra.deg,self.sci_c.dec.deg
                            self.ali_center_ra,self.ali_center_dec = self.sci_ra_d,self.sci_dec_d
                            self.sci_img_center = self.sci_img_hdu.data.shape[0]/2,self.sci_img_hdu.data.shape[1]/2
                            self.sci_wcs = WCS(self.sci_img_hdu.header)
                            self.sci_cnt = self.sci_wcs.all_pix2world(self.sci_img_center[0],self.sci_img_center[1],1)

                            self.sci_img_hdu=fits.open(self.sci_path)[0]
                            # self.sp_logger.info(self.sci_path)
                            self.redo_astrometry=False
            
                        if self.sci_img_hdu.header['ONTARGET']!=True:
                            # self.sp_logger.info(self.sci_img_hdu.header['ONTARGET'])
                            self.sp_logger.warning(warn_r+' This image is not on target, please check images, passing on this image')
                            self.sys_exit=True
                            return
                        self.sci_inst='SEDM-P60'
                        self.sci_prop=None

                        if 'rc' in self.sci_path.split("/")[-1]:
                            self.folder=self.sci_path.split('rc')[1][2:10]
                        elif self.sci_path.split('/')[-1].startswith('20') and any(self.sci_path.split('/')[-1].endswith(strn) for strn in ['p.fits','p_trimmed.fits']):
                            self.folder=self.sci_path.split('/')[-1][:8]
                        else:
                            self.folder=self.sci_path.split('/')[-1].split('_')[2]

                        # self.sp_logger.info(self.folder)

                    elif self.telescope=='SLT':
                        self.sci_filt=self.sci_img_hdu.header['FILTER'][0]
                        self.sci_inst = 'SLT-Andor'
                        self.sci_prop = None
                        self.folder = str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))

                    elif 'GTC' in self.telescope:
                        if self.sci_img_hdu.header['INSTRUME']=='OSIRIS':
                            self.sci_filt=self.sci_img_hdu.header['FILTER2'][-1]
                            self.sci_inst='OSIRIS'
                            self.folder=str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))

                        else:
                            self.sci_filt=input('Filter?: ')
                            self.sci_inst='HiPERCAM'
                            self.sci_prop=None
                            self.folder=str(re.sub('-','',self.sci_img_hdu.header['DATE'].split('T')[0]))
                    
                    elif self.telescope=='TJO':
                        self.sci_filt=self.sci_img_hdu.header['FILTER'][-1]
                        self.sci_inst='MEIA3'
                        self.folder=str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS']))


    
                else:
                    self.sp_logger.warning(warn_r+' This script is not designed for '+self.telescope)
                    self.sp_logger.warning(warn_r+' Please add the header keywords for '+self.telescope+' to the header_kw dictionary in the script')
                    self.sys_exit=True




                # sys.exit()
                #check if ONTARGET in header is True
                if self.sci_inst in SEDM and self.sci_img_hdu.header['ONTARGET']!=True:
                    self.sp_logger.warning(warn_r+' This image is not on target, please check images, passing on this image')
                    self.sys_exit=True
                    return
                self.sci_inst='SEDM-P60'
                self.sci_prop=None
                self.sci_obj=self.sci_img_hdu.header[self.OBJ_kw]
                if self.sci_obj=='2024afav':
                    self.sci_img_hdu.header[self.RA_kw] = "12:49:12.10"
                    self.sci_img_hdu.header[self.DEC_kw] = "-18:06:13.20"

                if 'rc' in self.sci_path.split("/")[-1]:
                    self.folder=self.sci_path.split('rc')[1][2:10]
                elif self.sci_path.split('/')[-1].startswith('20') and any(self.sci_path.split('/')[-1].endswith(strn) for strn in ['p.fits','p_trimmed.fits']):
                    self.folder=self.sci_path.split('/')[-1][:8]
                else:
                    self.folder=self.sci_path.split('/')[-1].split('_')[2]

                
                self.sp_logger.info(info_g+f' Single filter: '+'\033[1m'+self.sci_filt+'\033[0m')

                try:
                    self.crval1,self.crval2,self.crpix1,self.crpix2 = self.sci_img_hdu.header['CRVAL1'],self.sci_img_hdu.header['CRVAL2'],self.sci_img_hdu.header['CRPIX1'],self.sci_img_hdu.header['CRPIX2']
                    self.crdelt1,self.crdelt2,self.naxis1,self.naxis2 = self.sci_img_hdu.header['CDELT1'],self.sci_img_hdu.header['CDELT2'],self.sci_img_hdu.header['NAXIS1'],self.sci_img_hdu.header['NAXIS2']
                except:
                    pass
                if self.ra_dec_pos=='header':
                    self.sci_ra=self.sci_img_hdu.header[self.RA_kw]
                    self.sci_dec=self.sci_img_hdu.header[self.DEC_kw]
                    # self.sp_logger.info(self.sci_ra,self.sci_dec)
                else:
                    self.sp_logger.info(info_g+f' RA & Dec specified by user (J2000): '+'\033[1m'+self.ra_dec_pos[0]+'\033[0m'+'\033[1m'+self.ra_dec_pos[1]+'\033[0m')
                    self.sp_logger.info(info_g+f' Catlog RA & Dec from header (J2000): '+'\033[1m'+self.sci_img_hdu.header[self.RA_kw]+'\033[0m'+'\033[1m'+self.sci_img_hdu.header[self.DEC_kw]+'\033[0m')
                    self.sci_ra,self.sci_dec = self.ra_dec_pos[0],self.ra_dec_pos[1]
                    self.sci_c= SkyCoord(self.sci_ra,self.sci_dec, unit=(u.hourangle, u.deg),frame='fk5')
                    self.sci_ra_d,self.sci_dec_d = self.sci_c.ra.deg,self.sci_c.dec.deg
                    # self.sp_logger.info(info_g+' RA & Dec (deg):','\033[1m'+str(round(self.sci_ra_d,2))+'\033[0m','\033[1m'+str(round(self.sci_dec_d,2))+'\033[0m')
                    self.sp_logger.info(info_g+f' Checking if RA & Dec are within the image')
                    self.X_pix_co = (self.sci_ra_d - self.crval1) / self.crdelt1 + self.crpix1
                    self.X_in = 0<self.X_pix_co<self.naxis1

                    self.Y_pix_co = (self.sci_dec_d - self.crval2) / self.crdelt2 + self.crpix2
                    self.Y_in = 0<self.Y_pix_co<self.naxis2

                    if any([self.X_in,self.Y_in])==False:
                        self.sp_logger.info(warn_r+f' RA or Dec are outside the image')
                        if self.X_in==False:
                            self.sp_logger.warning(warn_r+f' RA is outside the image')
                        if self.Y_in==False:
                            self.sp_logger.warning(warn_r+f' Dec is outside the image')
                        self.sp_logger.warning(warn_r+f" Exiting ")
                        self.sys_exit=True
                        return
                    else:
                        self.sp_logger.info(info_g+f' RA & Dec are within the image')


                self.sci_date_obs = self.sci_img_hdu.header[self.DATE_kw]

                if self.UTS_kw == '-':
                    self.sci_utstart = self.sci_img_hdu.header[self.DATE_kw].split('T')[1]
                else:
                    self.sci_utstart = self.sci_img_hdu.header[self.UTS_kw]
                self.sci_obj, self.sep, self.tail = self.sci_obj.partition('_')
                
                if self.AIRM_kw=='-':
                    self.sci_airmass = ''
                else:
                    self.sci_airmass = self.sci_img_hdu.header[self.AIRM_kw]


                try:self.sci_seeing = self.sci_img_hdu.header[self.SEE_kw]
                except Exception as e:self.sci_seeing = estimate_seeing(self.sci_path)

                # self.sp_logger.info(type(self.GAIN_kw)==float)

                if self.GAIN_kw=='-': 
                    self.sci_gain=1.5
                    self.sp_logger.warning(info_g+f' No GAIN keyword found, using GAIN=1.5 as default')
                elif type(self.GAIN_kw)==float:
                    self.sci_gain = self.GAIN_kw
                else:self.sci_gain = float(self.sci_img_hdu.header[self.GAIN_kw])

                if self.ESEE_kw != None:
                    self.sci_est_seeing = self.sci_img_hdu.header[self.ESEE_kw]
                    try:self.sci_bkg_med,self.sci_bkg_mean = self.sci_img_hdu.header[self.BKGMED_kw],self.sci_img_hdu.header[self.BKGMEA_kw]
                    except:self.sci_bkg_med,self.sci_bkg_mean = None,None
                    try:self.sci_prop = self.sci_img_hdu.header[self.PRO_kw]
                    except:self.sci_prop = None
                else:
                    self.sci_est_seeing = None
                    self.sci_bkg_med,self.sci_bkg_mean = None,None
                    self.sci_prop = None

            
                self.sci_obj, self.sep, self.tail = self.sci_obj.partition('_') #tail should be the request id from marshal triggering if there is one
                if self.sci_obj=='ZTFaaqousn':self.sci_obj='ZTF25aaqousn'
                
                if ' ' in self.sci_obj:
                    self.sci_obj=self.sci_obj.split(' ')[0]
                if self.sci_obj=='2023vyl':self.sci_obj='ZTF23abnprwj'
                if 'ACQ-' in self.sci_obj:self.sci_obj=self.sci_obj.split('ACQ-')[1]
                
                self.sci_img_name=self.sci_obj+'_'+self.sci_filt+'comb.fits'
                self.sci_exp_time=self.sci_img_hdu.header[self.EXPT_kw]

                if ':' in str(self.sci_ra):
                    # self.sp_logger.info(info_g+' RA and DEC not in decimal degrees')
                    #RA & DEC already noy decimal degrees

                    if all(s not in str(self.sci_dec) for s in ['+','-']):
                        # self.sp_logger.info(warn_y+' No sign on declination, assuming to be positive')
                        self.sci_dec = '+'+self.sci_dec
                    self.sci_c= SkyCoord(self.sci_ra,self.sci_dec, unit=(u.hourangle, u.deg),frame='fk5')
                else:
                    # print(self.sci_ra,self.sci_dec)
                    self.sci_c= SkyCoord(str(self.sci_ra)+' '+str(self.sci_dec), unit=(u.deg, u.deg),frame='fk5')

                if self.PIX_kw == 'in_comments':
                    self.scale_comments = [s for s in self.sci_img_hdu.header["COMMENT"] if 'scale:' in s]
                    self.scale_comment = float(self.scale_comments[0].split(' ')[1])
                    self.sci_ps = self.scale_comment
                elif self.PIX_kw =='-':
                    try:
                        self.sci_ps = self.estimate_ps()
                        self.sp_logger.info(info_g+f' Estimating the pixel scale to be" {self.sci_ps}')
                    except Exception as e:
                        self.sp_logger.warning(e)
                else:
                    if type(self.PIX_kw)!=str: self.sci_ps=float(self.PIX_kw)
                    else:self.sci_ps=self.sci_img_hdu.header[self.PIX_kw]
                self.sci_mjd = self.sci_img_hdu.header[self.MJD_kw]

                if self.sci_mjd>400000:
                    self.sci_jd=self.sci_mjd
                    self.sci_mjd=self.sci_jd-2400000.5
                else:
                    self.sci_mjd=self.sci_mjd
                    self.sci_jd=self.sci_mjd+2400000.5

                self.MJD = [self.sci_mjd]
                
                self.sci_img_name=self.sci_obj+'_'+self.sci_filt+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'.fits'
                self.sci_img=self.sci_img_hdu.data

        if len(self.ims)>=2:
            self.sp_logger.info(info_g+' Multiple images to stack')
            if self.ims[0].startswith('h_')==False:    #checking whether name is of a single objet or the name starts with full path
                self.sp_logger.info(info_g+' Assuming all images to be stacked are in the same directory')
                self.arg_split = self.ims[0].split('/')
                for s in range(len(self.arg_split)):
                    if 'h_' in self.arg_split[s]:
                        self.folder = self.arg_split[s].split('_')[2]
                        # self.stack_data_path = '/'.join(self.arg_split[:s]) #prefix to the first fits file showing the directory that contains all the fits to be stacked
                    # else:
            self.stack_data_path = '/'.join(self.ims[0].split('/')[:-1])
            # self.sp_logger.info(self.stack_data_path)

            self.img_type="_stacked"
            self.name = [self.ims[0]+".fits"]
            for s in range(1,len(self.ims)):
                if not self.ims[s].startswith(self.stack_data_path): self.ims[s]=self.stack_data_path+self.ims[s]
                self.name.append(str(self.ims[s])+".fits")

            for k in range(len(self.name)):
                if not self.name[k].startswith(self.path): self.name[k]=self.path+self.name[k]
            self.sci_img_hdu = fits.open(self.name[0])[0]
            self.telescope = self.sci_img_hdu.header['TELESCOP']

            if self.telescope=='GTC':
                if self.sci_img_hdu.header['INSTRUME']=='HIPERCAM':self.telescope='GTC-HIPERCAM'
                else: self.telescope='GTC-OSIRIS'

            if self.telescope in header_kw.keys():
                self.sp_logger.info(info_g+' Telescope:'+self.telescope)
                self.tele_kw = header_kw[self.telescope]
                
                self.FILT_kw,self.OBJ_kw,self.RA_kw = self.tele_kw['filter'],self.tele_kw['object'],self.tele_kw['ra']
                self.DEC_kw,self.AIRM_kw,self.UTS_kw = self.tele_kw['dec'],self.tele_kw['airmass'],self.tele_kw['utstart']
                self.INST_kw,self.EXPT_kw,self.DATE_kw,self.MJD_kw = self.tele_kw['instrument'],self.tele_kw['exptime'],self.tele_kw['date'],self.tele_kw['mjd']
                self.PIX_kw,self.GAIN_kw = self.tele_kw['pixscale'],self.tele_kw['gain']
                self.SEE_kw = self.tele_kw['seeing']
                self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = None,None,None
                self.ESEE_kw = None
                self.sci_prop = None

                self.sci_date_obs = self.sci_img_hdu.header[self.DATE_kw]
                if self.UTS_kw == '-':
                    self.sci_utstart = self.sci_img_hdu.header[self.DATE_kw].split('T')[1]
                else:
                    self.sci_utstart = self.sci_img_hdu.header[self.UTS_kw]
                
                # self.sp_logger.info(self.telescope),sys.exit()
                if self.telescope=='Liverpool Telescope':
                    self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = self.tele_kw['propid'],self.tele_kw['bkg_med'],self.tele_kw['bkg_mean']
                    self.ESEE_kw = self.tele_kw['est_seeing']
                    self.sci_filt=self.sci_img_hdu.header['FILTER1'][5].lower()
                    self.sci_inst = 'IO:O'

                elif self.telescope in ['eso-ntt','ntt-efosc','ntt','efosc','ESO-NTT','NTT-EFOSC','NTT','EFOSC']:
                    self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = '-','-','-'
                    self.ESEE_kw = self.tele_kw['seeing']
                    self.sci_filt=self.sci_img_hdu.header['ESO INS FILT1 NAME'][0]
                    self.sci_inst = 'EFOSC'
                    self.sci_prop = None

                
                elif self.telescope=='HCT':
                    self.sci_filt=self.sci_img_hdu.header['FILTER'].lower()
                    self.sci_inst = 'HFOSC'

                elif self.telescope=='SLT':
                    self.sci_filt=self.sci_img_hdu.header['FILTER'][0]
                    self.sci_inst = 'SLT-Andor'
                    self.sci_prop = None
                    self.folder = str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))

                elif 'GTC' in self.telescope:
                    if self.sci_img_hdu.header['INSTRUME']=='OSIRIS':
                        self.sci_filt=self.sci_img_hdu.header['FILTER2'][-1]
                        self.sci_inst='OSIRIS'
                        self.folder=str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))

                    else:
                        self.sci_filt=input('Filter?')
                        self.sci_inst='HiPERCAM'
                        self.sci_prop=None
                        self.folder=str(re.sub('-','',self.sci_img_hdu.header['DATE'].split('T')[0]))

                elif self.telescope=='TJO':
                    self.sci_filt=self.sci_img_hdu.header['FILTER'][-1]
                    self.sci_inst='MEIA3'
                    self.folder=str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS']))


            else:
                self.sp_logger.warning(warn_r+' This script is not designed for '+self.telescope)
                self.sp_logger.warning(warn_r+' Please add the header keywords for '+self.telescope+' to the header_kw dictionary in the script')
                self.sys_exit=True

            self.MJD = [np.round(fits.open(self.name[n])[0].header[self.MJD_kw],6) for n in range(len(self.name))]

            #this code will convert the self.MJD to self.MJD if it is in JD format
            if self.MJD[0]>400000:
                self.MJD = [self.MJD[n]-2400000.5 for n in range(len(self.MJD))]
            

            try:self.see_all = [fits.open(self.name[n])[0].header[self.SEE_kw] for n in range(len(self.name))]
            except Exception as e:self.see_all = [estimate_seeing(self.name[n]) for n in range(len(self.name))]

            self.see_std,self.see_med = np.std(self.see_all),np.median(self.see_all)
            # print(self.see_all,self.see_std,self.see_med)
            self.names_final = []
            if self.see_std>0:
                for n in range(len(self.see_all)):
                    if self.see_all[n]<self.see_med+5 * self.see_std:
                        print(info_g+f' Image {self.name[n]} has seeing of {self.see_all[n]} arcsec, adding to stack')
                        self.names_final.append(self.name[n])
            
                self.name = self.names_final
            # print(self.name)

            if len(self.name)>0:
                self.no_stacked = len(self.name)
                self.if_stacked=True
                self.name_joined=' '.join(self.name)
                self.sci_img_hdu=fits.open(self.name[0])[0]
                self.sci_obj=self.sci_img_hdu.header[self.OBJ_kw]
                if self.sci_obj=='2023vyl':self.sci_obj='ZTF23abnprwj'
                if self.sci_obj=='ZTFaaqousn':self.sci_obj='ZTF25aaqousn'


                if self.sci_obj=='2024afav':
                    self.sci_img_hdu.header[self.RA_kw] = "12:49:12.10"
                    self.sci_img_hdu.header[self.DEC_kw] = "-18:06:13.20"


                if self.ra_dec_pos=='header':
                    self.sp_logger.info(info_g+' RA & Dec specified by header')
                    self.sci_ra=self.sci_img_hdu.header[self.RA_kw]
                    self.sci_dec=self.sci_img_hdu.header[self.DEC_kw]
                else:
                    self.sp_logger.info(info_g+' RA & Dec specified by user: '+'\033[1m'+self.ra_dec_pos[0]+'\033[0m'+'\033[1m'+self.ra_dec_pos[1]+'\033[0m')
                    self.sp_logger.info(info_g+' Catlog RA & Dec from header (J2000): '+'\033[1m'+self.sci_img_hdu.header[self.RA_kw]+'\033[0m'+'\033[1m'+self.sci_img_hdu.header[self.DEC_kw]+'\033[0m')
                    self.sci_ra,self.sci_dec = self.ra_dec_pos[0],self.ra_dec_pos[1]
                    self.sci_c= SkyCoord(self.sci_ra,self.sci_dec, unit=(u.hourangle, u.deg),frame='fk5')
                    self.sci_ra_d,self.sci_dec_d = self.sci_c.ra.deg,self.sci_c.dec.deg
                    # self.sp_logger.info(info_g+' RA & Dec (deg):','\033[1m'+str(self.sci_ra_d)+'\033[0m','\033[1m'+str(self.sci_dec_d)+'\033[0m')
                    self.sp_logger.info(info_g+' Checking if RA & Dec are within the image')
                    self.X_pix_co = (self.sci_ra_d - self.crval1) / self.crdelt1 + self.crpix1
                    self.X_in = 0<self.X_pix_co<self.naxis1

                    self.Y_pix_co = (self.sci_dec_d - self.crval2) / self.crdelt2 + self.crpix2
                    self.Y_in = 0<self.Y_pix_co<self.naxis2

                    if any([self.X_in,self.Y_in])==False:
                        self.sp_logger.warning(warn_r+' RA or Dec are outside the image')
                        if self.X_in==False:
                            self.sp_logger.warning(warn_r+' RA is outside the image')
                        if self.Y_in==False:
                            self.sp_logger.warning(warn_r+' Dec is outside the image')
                        self.sp_logger.warning(warn_r+" Exiting ")
                        self.sys_exit=True
                        return

            
                try:self.sci_utstart = self.sci_img_hdu.header[self.UTS_kw]
                except:self.sci_utstart = self.sci_img_hdu.header[self.DATE_kw].split('T')[1]
                self.sci_obj, self.sep, self.tail = self.sci_obj.partition('_')
                if self.GAIN_kw=='-': 
                    self.sci_gain=1.5
                    self.sp_logger.warning(info_g+f' No GAIN keyword found, using GAIN=1.5 as default')
                elif type(self.GAIN_kw)==float:
                    self.sci_gain = self.GAIN_kw
                else:self.sci_gain = float(self.sci_img_hdu.header[self.GAIN_kw])
                

                if self.ESEE_kw != None:
                    self.sci_est_seeing = self.sci_img_hdu.header[self.ESEE_kw]
                    try:self.sci_bkg_med,self.sci_bkg_mean = self.sci_img_hdu.header[self.BKGMED_kw],self.sci_img_hdu.header[self.BKGMEA_kw]
                    except:self.sci_bkg_med,self.sci_bkg_mean = None,None
                    try:self.sci_prop = self.sci_img_hdu.header[self.PRO_kw]
                    except:self.sci_prop = None
                    self.sci_airmass = self.sci_img_hdu.header[self.AIRM_kw]

                else:
                    self.sci_est_seeing = None
                    self.sci_bkg_med,self.sci_bkg_mean = None,None
                    self.sci_prop = None
                    self.sci_airmass = None
                
                try:self.sci_seeing=self.sci_img_hdu.header[self.SEE_kw]
                except:self.sci_seeing=np.median(self.see_all)
                self.sci_exp_time=self.sci_img_hdu.header[self.EXPT_kw]*(len(self.ims)-1)
                if ':' in str(self.sci_ra):
                    if all(s not in str(self.sci_dec) for s in ['+','-']):
                        self.sp_logger.info(warn_y+' No sign on declination, assuming to be positive')
                        self.sci_dec = re.sub(' ','','+'+self.sci_dec)
                    self.sci_c= SkyCoord(ra=self.sci_ra,dec=self.sci_dec, frame='fk5', unit=(u.hourangle, u.deg))
                else:
                    self.sci_c= SkyCoord(str(self.sci_ra)+' '+str(self.sci_dec), unit=(u.deg, u.deg),frame=FK5)

                if self.PIX_kw == 'in_comments':
                    self.scale_comments = [s for s in self.sci_img_hdu.header[self.PIX_kw] if 'scale:' in s]
                    self.scale_comment = float(self.scale_comments[0].split(' ')[1])
                else:
                    if type(self.PIX_kw)!=str: self.sci_ps=float(self.PIX_kw)
                    else:self.sci_ps=self.sci_img_hdu.header[self.PIX_kw]

                self.sci_mjd = self.sci_img_hdu.header[self.MJD_kw]

                if self.sci_mjd>400000:
                    self.sci_jd = self.sci_mjd
                    self.sci_jd=np.mean([fits.open(s)[0].header[self.MJD_kw] for s in self.name])
                    self.sci_mjd-=2400000.5
                else:
                    self.sci_jd = self.sci_mjd+2400000.5
                    self.sci_jd=np.mean([fits.open(s)[0].header[self.MJD_kw]+2400000.5 for s in self.name])
                
                if self.sci_obj=='ZTFaaqousn':self.sci_obj='ZTF25aaqousn'
                self.sci_img_name=self.sci_obj+'_'+self.sci_filt+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]),seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'.fits'
                
                # plt.figure(figsize=(10,10))
                # vmin,vmax = visualization.ZScaleInterval().get_limits(stacked_data_new)
                # plt.imshow(stacked_data_new,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)
                # plt.show()
                # sys.exit()

            
                if not os.path.exists(self.path+'combined_imgs/'+self.sci_img_name):

                    self.sp_logger.info(info_g+' Combining '+str(len(self.ims))+' '+self.sci_filt+' images...')
                    self.sp_logger.info(info_g+f" Stacking {self.sci_obj} {self.sci_filt} images with Swarp")

                    self.ra_string="%.6f" % round(self.sci_c.ra.deg, 6)
                    # self.sp_logger.info(str(self.sci_c.dec)[0]=='-')
                    if self.sci_c.dec.deg>0:
                        self.dec_string="+"+"%.6f" % round(self.sci_c.dec.deg, 6)
                    elif str(self.sci_c.dec)[0]=='-':
                        self.dec_string="%.6f" % round(self.sci_c.dec.deg, 6)
                    else:
                        self.dec_string='-'+"%.6f" % round(self.sci_c.dec.deg, 6)
                    # self.sp_logger.info(self.ra_string,self.dec_string)
                    if self.sci_img_name.startswith(' '):self.sci_img_name = self.sci_img_name[1:]
                    self.img_size1,self.img_size2 = self.sci_img_hdu.header['NAXIS1']*0.9999,self.sci_img_hdu.header['NAXIS2']*0.9999
                    self.comb_centre = str(self.sci_img_hdu.header['CRVAL1'])+" "+str(self.sci_img_hdu.header['CRVAL2'])

                    swarp_command=swarp_path+" "+self.name_joined+" -c "+self.path+"config_files/config_comb.swarp -COPY_KEYWORDS DATE-OBS -CENTER '"+self.comb_centre+"' -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -IMAGEOUT_NAME "+self.path+"combined_imgs"+'/'+self.sci_img_name+" -RESAMPLE Y -RESAMPLE_DIR '"+self.path+"' -COMBINE Y -IMAGE_SIZE '"+str(self.img_size1)+","+str(self.img_size2)+"'"
                    # self.sp_logger.info(swarp_command)
                    # sys.exit()
                    # self.sp_logger.info(self.name)
                    # self.stacked_path, _ = self.stack_images(self.name,save_name=path+"combined_imgs"+'/'+self.sci_img_name)
                    os.system(swarp_command)
                    # self.sp_logger.info('stacked')
                    # sys.exit()
                    
                    # self.status=os.system(swarp_command)

                    self.sp_logger.info(info_g+' Combined '+str(len(self.ims))+' '+self.sci_filt+' images!')
                    self.sci_img_hdu=fits.open(self.path+'combined_imgs/'+self.sci_img_name)[0]
                    self.sci_img=self.sci_img_hdu.data
                    self.sci_path = self.path+'combined_imgs/'+self.sci_img_name
                    self.files_to_clean.append(self.path+'combined_imgs/'+self.sci_img_name)
                    self.name=self.sci_path
        
                elif os.path.exists(self.path+'combined_imgs/'+self.sci_img_name):

                    self.sp_logger.info(info_b+' Images: '+self.sci_filt+' already combined')
                    self.sci_img_hdu=fits.open(self.path+'combined_imgs/'+self.sci_img_name)[0]
                    self.sci_img=self.sci_img_hdu.data
                    self.sci_path = self.path+'combined_imgs/'+self.sci_img_name
                    self.name=self.sci_path

            else:
                self.sys_exit=True
                pass
    

        
        
        if self.sci_bkg_med!=None and self.sci_bkg_med<self.sci_exp_time/100 and self.sci_bkg_mean<self.sci_exp_time/100 and self.sci_filt!='u' and self.sys_exit!=True:  
            self.sp_logger.warning(warn_r+f' Median of background counts is {self.sci_bkg_med}, possible shutter issue')
            self.sp_logger.warning(warn_r+' Exiting!!')
            self.sys_exit=True


        if self.sys_exit!=True:
            if all(x!=None for x in [self.sci_seeing,self.sci_est_seeing]):
                self.sp_logger.info(info_g+' Seeing (observed,estimated): '+'\033[1m'+"%.2f %.2f" %(self.sci_seeing,self.sci_est_seeing)+'\033[0m')
            elif self.sci_seeing!=None and self.sci_est_seeing==None:
                self.sp_logger.info(info_g+' Seeing (observed,estimated): '+'\033[1m'+"%.2f %.2f" %(self.sci_seeing,0)+'\033[0m')

        if self.sys_exit!=True:
            # self.sp_logger.info(self.sci_seeing>5)
            if self.sci_seeing!=None and self.sci_est_seeing!=None:
                # self.sp_logger.info((self.sci_seeing>=self.sci_est_seeing*20 and self.sci_seeing>5.0),self.sci_seeing>7.5)
                if (self.sci_seeing>=self.sci_est_seeing*20 and self.sci_seeing>5.0)==True or self.sci_seeing>5:
                    self.sp_logger.warning(warn_r+' Poor seeing, exiting')
                
                    self.sp_logger.warning(warn_r+f" Seeing error: "+self.sci_obj+","+self.sci_filt+" band, Actual seeing:"+str(self.sci_seeing)+", Estimated seeing:"+str(self.sci_est_seeing)+" (ie. guiding/focus error or seeing got much worse)\n")

                    self.sys_exit=True
                    # return

        if self.out_dir=='by_name' and self.morning_round_up==False:
            self.folder = self.sci_obj
            self.out_dir='photometry/'
            if not os.path.exists(self.path+self.out_dir):
                self.sp_logger.info(info_g+' Creating photometry directory: '+self.path+self.out_dir)
                os.mkdir(self.path+self.out_dir)
            
            #folder is for easy way to display photometry by date of observation
            if os.path.exists(self.path+'photometry_date')==False:
                self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+'photometry_date')
                os.mkdir(self.path+'photometry_date')
            if not os.path.exists(self.path+f'photometry_date/{self.folder}'):
                self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+f'photometry_date/{self.folder}')
                os.mkdir(self.path+f'photometry_date/{self.folder}')  

            if self.morning_round_up!=False:
                if not os.path.exists(self.path+f'photometry_date/{self.folder}/morning_rup'):
                    self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+f'photometry_date/{self.folder}/morning_rup')
                    os.mkdir(self.path+f'photometry_date/{self.folder}/morning_rup')

            if self.cutout_tf!=False:
                if not os.path.exists(self.path+f'photometry_date/{self.folder}/cut_outs'):
                    self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+f'photometry_date/{self.folder}/cut_outs')
                    os.mkdir(self.path+f'photometry_date/{self.folder}/cut_outs')

        if self.out_dir=='by_obs_date' or self.morning_round_up!=False:
            self.out_dir='photometry/'
            # self.folder = self.sci_obj
            if not os.path.exists(self.path+self.out_dir):
                self.sp_logger.info(info_g+' Creating photometry directory: '+self.path+self.out_dir)
                os.mkdir(self.path+self.out_dir)
            
            #folder is for easy way to display photometry by date of observation
            if os.path.exists(self.path+'photometry_date')==False:
                self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+'photometry_date')
                os.mkdir(self.path+'photometry_date')
            if not os.path.exists(self.path+f'photometry_date/{self.folder}'):
                self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+f'photometry_date/{self.folder}')
                os.mkdir(self.path+f'photometry_date/{self.folder}')  

            if self.morning_round_up!=False:
                if not os.path.exists(self.path+f'photometry_date/{self.folder}/morning_rup'):
                    self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+f'photometry_date/{self.folder}/morning_rup')
                    os.mkdir(self.path+f'photometry_date/{self.folder}/morning_rup')

            if self.cutout_tf!=False:
                if not os.path.exists(self.path+f'photometry_date/{self.folder}/cut_outs'):
                    self.sp_logger.info(info_g+' Creating photometry_date directory: '+self.path+f'photometry_date/{self.folder}/cut_outs')
                    os.mkdir(self.path+f'photometry_date/{self.folder}/cut_outs')
        else:
            self.out_dir=str(self.out_dir)+'/'
            if not os.path.exists(self.path+self.out_dir):
                self.sp_logger.info(info_g+' Creating photometry directory: '+self.path+self.out_dir)
                os.mkdir(self.path+self.out_dir)

            if self.morning_round_up!=False:
                if not os.path.exists(self.path+self.out_dir+'morning_rup'):
                    self.sp_logger.info(info_g+' Creating photometry directory: '+self.path+self.out_dir+'morning_rup')
                    os.mkdir(self.path+self.out_dir+'morning_rup')
            
            if self.cutout_tf!=False:
                if not os.path.exists(self.path+f'{self.out_dir}cut_outs'):
                    self.sp_logger.info(info_g+' Creating photometry directory: '+self.path+f'{self.out_dir}cut_outs')
                    os.mkdir(self.path+f'{self.out_dir}cut_outs')


        if ':' in str(self.sci_ra):self.sci_c=SkyCoord(self.sci_ra,self.sci_dec,unit=(u.hourangle, u.deg),frame='fk5')
        else:self.sci_c=SkyCoord(self.sci_ra,self.sci_dec,unit=(u.deg, u.deg),frame='fk5')
        self.ra_string="%.6f" % round(self.sci_c.ra.deg, 6)
        self.dec_string="%.6f" % round(self.sci_c.dec.deg, 6)
        self.ref_width=self.sci_img_hdu.header['NAXIS2']*self.sci_ps*1.1/60.
        # self.ref_width*=4
        # self.sp_logger.info(self.ref_width)
        # self.sp_logger.info(image_size*self.sci_ps*1.1/60.)
        # self.ref_width = image_size*self.sci_ps*1.1/60.
        self.sci_mean, self.sci_median, self.sci_std = sigma_clipped_stats(self.sci_img_hdu.data)
        self.coords_sn_sci=wcs_to_pixels(self.sci_path,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
        self.coords_sn_sci_x,self.coords_sn_sci_y = self.coords_sn_sci
        self.sp_logger.info(info_g+' Object name:'+' \033[1m'+self.sci_obj+'\033[0m')
        if self.sci_prop!=None: self.sp_logger.info(info_g+' Object proposal ID:'+' \033[1m'+self.sci_prop+'\033[0m')
        self.sp_logger.info(info_g+' Object position (RA,Dec) J2000:'+'\033[1m'+' ('+str(self.sci_c.to_string("hmsdms"))+')\033[0m')
        self.sp_logger.info(info_g+' Object position (RA,Dec) deg:'+'\033[1m'+' ('+self.ra_string+','+self.dec_string+')'+ '\033[0m')
        self.sp_logger.info(info_g+' Object position (x,y)'+'\033[1m'+' ('+str(int(self.coords_sn_sci_x))+','+str(int(self.coords_sn_sci_y))+')'+ '\033[0m')
        self.sp_logger.info(info_g+' Image dimensions (x,y):'+'\033[1m'+' ('+str(self.sci_img_hdu.header['NAXIS1'])+'x'+str(self.sci_img_hdu.header['NAXIS2'])+')'+ '\033[0m')
        self.sp_logger.info(info_g+f' {self.ref_width}'+' arcmin reference width')
        self.sp_logger.info(info_g+' Observation date:'+' \033[1m'+str(round(self.sci_mjd,3))+'\033[0m')
        self.sp_logger.info(info_g+' Observation time:'+' \033[1m'+str(self.sci_utstart)+'\033[0m')
        self.sp_logger.info(info_g+' Exposure time:'+' \033[1m'+str(self.sci_exp_time)+'\033[0m')
        self.name=self.sci_path

        if self.forced_phot in ['cat','fits','head',]:
            self.forced_phot = [self.sci_ra,self.sci_dec]
        if self.sci_obj in ['SN2023ixf','ZTF23aaklqou','2023ixf','23aaklqou'] or self.special_case!=None:
            self.special_case='SN2023ixf'
            self.sp_logger.info(info_g+' Special case: '+self.special_case)
        else:
            self.special_case=self.special_case

        


        if self.if_stacked==False:
            self.phot_path = self.path+'phot_fits_info/'+self.sci_img_name[:-5]+"_phot_info.txt"
            self.fits_info = [str(self.sci_mjd),self.sci_utstart,self.sci_ra, self.sci_dec,self.sci_prop,self.sci_inst,self.sci_exp_time,self.sci_airmass, self.sci_seeing, self.sci_est_seeing,0]
            self.fits_info_txt =open(self.phot_path,"w")
            self.fits_info_txt.write(f'{self.fits_info[0]},{self.fits_info[1]},{self.fits_info[2]},{self.fits_info[3]},{self.fits_info[4]},{self.fits_info[5]},{self.fits_info[6]},{self.fits_info[7]},{self.fits_info[8]},{self.fits_info[9]},1')
            self.fits_info_txt.close()

        else:
            self.phot_path = self.path+'phot_fits_info/'+self.sci_img_name[:-5]+"_stacked_phot_info.txt"
            self.fits_info = [str(self.sci_mjd),self.sci_utstart,self.sci_ra,self.sci_dec,self.sci_prop,self.sci_inst,self.sci_exp_time, self.sci_airmass,self.sci_seeing,self.sci_est_seeing,self.no_stacked]
            self.fits_info_txt=open(self.phot_path,"w")
            self.fits_info_txt.write(f'{self.fits_info[0]},{self.fits_info[1]},{self.fits_info[2]},{self.fits_info[3]},{self.fits_info[4]},{self.fits_info[5]},{self.fits_info[6]},{self.fits_info[7]},{self.fits_info[8]},{self.fits_info[9]},{self.no_stacked}')
        


    def estimate_ps(self):
        CD1_1,CD1_2 = self.sci_img_hdu.header['CD1_1'],self.sci_img_hdu.header['CD1_2']
        CD2_1,CD2_2 = self.sci_img_hdu.header['CD2_1'],self.sci_img_hdu.header['CD2_2']
        # print(np.sqrt(CD1_1*2 + CD1_2*2 +CD2_1*2 + CD2_2*2 )),sys.exit()
        return 


    def make_sdss_ref(self,sdss_filt='u'):
        self.sp_logger.info(info_g+f" Making SDSS reference image")
        #if self.sci_filt=='u':
        self.spacing=0.1
        self.sdss_images=[]
        self.grid_size=3
        self.image_name=self.sci_obj+f'_ref_{sdss_filt}.fits'
        self.ref_path=self.path+'ref_imgs/'+self.image_name
        if os.path.exists(self.ref_path):
            self.sp_logger.info(info_b+' SDSS reference image already exists')
            return self.ref_path
    
        self.nx_range = range(-2,2)
        self.ny_range = range(-2,2)
        # sdss_args = [(nx, ny) for nx in self.nx_range for ny in self.ny_range]

        sdss_args = [(nx, ny,self.sci_c,self.sci_filt,self.sp_logger) for nx in self.nx_range for ny in self.ny_range]

        # for i in range(len(self.nx_range)):
        #     sdss_args.append((self.nx_range[i],self.ny_range[i],self.sci_c,self.sci_filt))

        for i in range(len(sdss_args)):
            if i==0:
                lnks_done=[]

            image = get_image(sdss_args[i] + (lnks_done,))
            image = image[0]
            lnks_done.append(image[1])

            self.sdss_images.append(image)
        # with concurrent.futures.ThreadPoolExecutor() as executor:
        #     futures = [executor.submit(get_image, sdss_args[i])#[0], sdss_args[i][1],sdss_args[i][2],sdss_args[i][3]) )
        #             # for nx in range(-3, 2) 
        #             # for ny in range(-3, 2)
        #             for i in range(len(sdss_args))]
            
        #     self.sdss_images = [f.result() for f in futures]

        # self.sp_logger.info(self.sdss_images)
        sdss_unique=[]
        [sdss_unique.append(s) for s in self.sdss_images if s not in sdss_unique and s is not None]
        self.sdss_images=np.array(sdss_unique)
        # sys.exit(1)
        # self.image_name=self.sci_obj+f'_ref.fits'
        # self.ref_path=self.path+'ref_imgs/'+self.image_name
        if not os.path.exists(self.ref_path):
            self.sdss_images = list(dict.fromkeys(self.sdss_images))
            with open(self.path+f'ref_imgs/{self.sci_obj}_sdss_list_{sdss_filt}.txt', 'w') as f:
                for item in self.sdss_images:
                    f.write("%s[0]\n" %(item))
                    self.sp_logger.info(info_b+" %s\n" %(item))
            f.close()
    
        # self.image_name=self.sci_obj+f'_ref_{sdss_filt}.fits'
        # self.ref_path=self.path+'ref_imgs/'+self.image_name
        if not os.path.exists(self.ref_path):
            self.sp_logger.info(info_g+' Combining SDSS images with swarp centering on the object at '+self.ra_string+' '+self.dec_string)
            self.swarp_sdss_command=swarp_path+" @"+self.path+f'ref_imgs/{self.sci_obj}_sdss_list_{sdss_filt}.txt'+" -c "+self.path+"config_files/swarp_sdss.conf -CENTER '"+self.ra_string+" "+self.dec_string+"'"+" -IMAGE_SIZE '"+'4000'+","+'4000'+"'"+" -IMAGEOUT_NAME "+self.path+'ref_imgs/'+self.sci_obj+f'_ref_{sdss_filt}.fits -VERBOSE_TYPE QUIET'
            self.sp_logger.info(self.swarp_sdss_command)
            # sys.exit(1)
            os.system(self.swarp_sdss_command)

        # self.sp_logger.info(self.ref_path)
        # sys.exit(1)
        return self.ref_path




    def bkg_subtract(self,sigma=3.):
        self.sp_logger.info(info_g+f" Subtracting background")
        self.t_bkg_sub_start = time.time()
        ##############################################
            #  Background subtraction
        #################################
        if not os.path.exists(self.path+'bkg_subtracted_science'):
            os.makedirs(self.path+'bkg_subtracted_science')

        self.sig_clip = SigmaClip(sigma=sigma)
        self.bkg_estimator = SExtractorBackground(self.sig_clip)
        self.sci_img = np.nan_to_num(self.sci_img)
        self.bkg = Background2D(
            self.sci_img, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=self.bkg_estimator)

        self.bkg_median = []
        for row in self.bkg.background:
            [self.bkg_median.append(i) for i in row]
            # row_median = stats.median(row)

        self.bkg_median = stats.median(self.bkg_median)
        if self.sci_bkg_med!=None and self.bkg_median<5 and self.sci_bkg_med<5 and self.sci_bkg_mean<5 and self.sci_filt!='u':  
            self.sp_logger.warning(warn_r+f' Median of background counts is {self.bkg_median}, possible shutter issue')
            self.sp_logger.warning(warn_r+f' {self.sci_obj} {self.sci_filt} band, Median background counts: {self.sci_bkg_med}, Mean background counts: {self.sci_bkg_mean}')
            self.sp_logger.warning(warn_r+' Exiting!!')
            self.sys_exit=True

            
            return

        
        if self.to_subtract==True: #subtracting the background
            # self.sp_logger.info('subtractiog background',self.bkg.background)
            # self.sp_logger.info(self.bkg_median)
            # self.sci_bkgsb=self.sci_img-np.median(self.sci_img)#
            self.sci_bkgsb = self.sci_img - self.bkg.background
            self.bkg_update = Background2D(self.sci_bkgsb, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=self.bkg_estimator)
            self.median_bkg=self.bkg_update.background_median

        if self.to_subtract==False: #not subtracting the background image
            self.sci_bkgsb=self.sci_img
            self.bkg_update = Background2D(self.sci_bkgsb, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=self.bkg_estimator)
            self.median_bkg=self.bkg_update.background_median


        self.t_bkg_sub_end = time.time()
        self.bkg_sub_removal_time = self.t_bkg_sub_end - self.t_bkg_sub_start


        self.sp_logger.info(info_g+f' Background subtraction time: {round(self.bkg_sub_removal_time, 0)} seconds')

       
        self.bkg_error = self.bkg.background_rms
        self.bkg_err = calc_total_error(self.sci_bkgsb,self.bkg_error,self.sci_gain)

        # print(self.median_bkg)
        # sys.exit()
        return self.sig_clip,self.bkg,self.sci_bkgsb,self.median_bkg


    def remove_cosmic(self):
        if self.termoutp!='quiet':
            self.sp_logger.info(info_g+f" Removing cosmic rays")
        self.t_cosmic_start = time.time()
        #based on how long the exposure times are
        if self.telescope not in SEDM:
            if 60<self.sci_exp_time<120.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=4, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)
    
            elif 30<self.sci_exp_time<=60.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=3, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)

            elif self.sci_exp_time>=120:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=6, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)

            elif self.sci_exp_time<=30:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=2, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)

        if self.telescope in SEDM:
            if 60<self.sci_exp_time<120.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=4, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)
              
            elif 30<self.sci_exp_time<=60.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=3, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)

            elif self.sci_exp_time>=120:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=6, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)
                
            elif self.sci_exp_time<=30:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, 
                gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=2, sepmed=True, cleantype='meanmask', 
                fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)
                

        self.sci_img_name=self.sci_img_name[:-5]+"bkgsub.fits"
        fits.writeto(self.path+'bkg_subtracted_science/'+self.sci_img_name, self.sci_bkgsb, self.sci_img_hdu.header,overwrite=True)

        self.files_to_clean.append(self.path+'bkg_subtracted_science/'+self.sci_img_name)

        self.t_cosmic_end = time.time()
        self.cosmic_removal_time = self.t_cosmic_end - self.t_cosmic_start

        self.sp_logger.info(info_g+' Cosmic Removal time: '+"%.1f" % round(self.cosmic_removal_time, 0)+' seconds')
        # sys.exit(1)
        
        return self.sci_img_name



    def find_dist_to_border(self,sci_data,ref_data,sci_match_coords,sci_phot_tab,ref_match_coords,ref_phot_tab,X=15,D=60):
        # X - distance from actual image edges
        # D - distance from border between images and 0 padding from SWarp of in swarp_ref_align
        # print(len(sci_match_coords),len(ref_match_coords))
        sci_keep_pix,sci_keep_sky = [],[]
        ref_keep_pix,ref_keep_sky = [],[]

        x_sci_len,y_sci_len = sci_data.shape
        x_ref_len,y_ref_len = ref_data.shape

        low_sci,high_sci = 0,y_sci_len
        left_sci,right_sci= 0,x_sci_len

        low_ref,high_ref = 0,y_ref_len
        left_ref,right_ref= 0,x_ref_len

        for h in range(len(sci_match_coords)):
            found_sci_low,found_sci_high,found_sci_left,found_sci_right = False,False,False,False
            found_ref_low,found_ref_high,found_ref_left,found_ref_right = False,False,False,False

            #check if they are close to the boundary where all pixels left, right, up or down are 0 or at a bound where all pixels left, right, up or down are <0
            sci_x,sci_y = int(sci_match_coords[h][0]),int(sci_match_coords[h][1])
            sci_ra,sci_dec = sci_phot_tab['ra'][h],sci_phot_tab['dec'][h]
            ref_x,ref_y = int(ref_match_coords[h][0]),int(ref_match_coords[h][1])
            ref_ra,ref_dec = ref_phot_tab['ra'][h],ref_phot_tab['dec'][h]

            for i in range(X,x_sci_len-X):
                j=x_sci_len-i-1
                if all(sci_data[:i,sci_y]==0) and any(sci_data[i:i+X,sci_y]==0)==False and found_sci_high==False:low_sci,found_sci_high = i,True
                if all(sci_data[j:,sci_y]==0) and any(sci_data[j-X:j,sci_y]==0)==False and found_sci_low==False:high_sci,found_sci_low = j,True
                if all(ref_data[:i,ref_y]==0) and any(ref_data[i:i+X,ref_y]==0)==False and found_ref_high==False:low_ref,found_ref_high = i,True
                if all(ref_data[j:,ref_y]==0) and any(ref_data[j-X:j,ref_y]==0)==False and found_ref_low==False:high_ref,found_ref_low = j,True
            
            for i in range(X,y_sci_len-X):
                j=y_sci_len-i-1
                if all(sci_data[sci_x,:i]==0) and any(sci_data[sci_x,i:i+X]==0)==False and found_sci_left==False:left_sci,found_sci_left = i,True
                if all(sci_data[sci_x,j:]==0) and any(sci_data[sci_x,j-X:j]==0)==False and found_sci_right==False:right_sci,found_sci_right = j,True
                if all(ref_data[ref_x,:i]==0) and any(ref_data[ref_x,i:i+X]==0)==False and found_ref_left==False:left_ref,found_ref_left = i,True
                if all(ref_data[ref_x,j:]==0) and any(ref_data[ref_x,j-X:j]==0)==False and found_ref_right==False:right_ref,found_ref_right = j,True


            # print('for the sci star',sci_x,sci_y)
            self.sp_logger.info(info_g+f" For star in scicence #{h} at {sci_x}, {sci_y}:"),self.sp_logger.info(info_g+f" - Borders left & right: {left_sci}, {right_sci}"),self.sp_logger.info(info_g+f" - Borders up & down: {low_sci}, {high_sci} ")
            sci_left_dist,sci_right_dist,sci_up_dist,sci_down_dist = abs(sci_x-left_sci),abs(sci_x-right_sci),abs(sci_y-low_sci),abs(sci_y-high_sci)
            self.sp_logger.info(info_g+f" - Distance left & right: {sci_left_dist}, {sci_right_dist}"),self.sp_logger.info(info_g+f" - Distance up, down: {sci_up_dist}, {sci_down_dist}")
            sci_far_from_border = (sci_left_dist>=D and sci_right_dist>=D and sci_up_dist>=D and sci_down_dist>=D)
            self.sp_logger.info(info_g+f" Star in sci {X} pix from all border? {sci_far_from_border}")

            self.sp_logger.info(info_g+f" For star in reference at {ref_x}, {ref_y}:"),self.sp_logger.info(info_g+f" - Borders left & right: {left_ref}, {right_ref}"),self.sp_logger.info(info_g+f" - Borders up & down: {low_ref}, {high_ref} ")
            ref_left_dist,ref_right_dist,ref_up_dist,ref_down_dist = abs(ref_x-left_ref),abs(ref_x-right_ref),abs(ref_y-low_ref),abs(ref_y-high_ref)
            self.sp_logger.info(info_g+f" - Distance left & right: {ref_left_dist}, {ref_right_dist}"),self.sp_logger.info(info_g+f" - Distance up, down: {ref_up_dist}, {ref_down_dist}")
            ref_far_from_border = (ref_left_dist>=D and ref_right_dist>=D and ref_up_dist>=D and ref_down_dist>=D)
            self.sp_logger.info(info_g+f" Star in ref {X} pix from all border? {ref_far_from_border}")

            if sci_far_from_border==True and ref_far_from_border==True:
                self.sp_logger.info(info_g+f" Keeping star sci: {sci_x}, {sci_y}, ref: {ref_x}, {ref_y}")
                sci_keep_pix.append([sci_x,sci_y]),ref_keep_pix.append([ref_x,ref_y])
                sci_keep_sky.append([sci_ra.value,sci_dec.value]),ref_keep_sky.append([ref_ra.value,ref_dec.value])
            else:
                self.sp_logger.info(warn_y+f" Removing star sci: {sci_x}, {sci_y}, ref: {ref_x}, {ref_y}")

            self.sp_logger.info(info_g+f" ")
        
        return {'sci_keep_pix':sci_keep_pix,'ref_keep_pix':ref_keep_pix,'sci_keep_sky':sci_keep_sky,'ref_keep_sky':ref_keep_sky}

    def distortion_correction(self):
        # return
        self.sp_logger.info(info_g+' Attempting to solve distortion in science image')
        self.ref_width=self.sci_img_hdu.header['NAXIS2']*self.sci_ps/60
            

        if self.auto_cat=="auto":
            if self.use_sdss==False and (self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z'):
                self.ref_cat=panstarrs_query(ra_deg=round(self.sci_c.ra.deg,6),dec_deg=round(self.sci_c.dec.deg,6), rad_deg=round(self.ref_width/60.,6))
                # self.ref_cat['iMeanPSFMag']-self.ref_cat['iMeanKronMag']<0.05 is to remove galaxies
                self.stars=np.where(np.abs((self.ref_cat['iMeanPSFMag']-self.ref_cat['iMeanKronMag'])<0.05) & (self.ref_cat[str(self.sci_filt)+'MeanPSFMagErr']<0.06) & (self.ref_cat[str(self.sci_filt)+'MeanPSFMag']<23.5)  & (self.ref_cat[str(self.sci_filt)+'MeanPSFMag']!=-999.))[0]
                # print(self.ref_cat[['objInfoFlag','qualityFlag','gFlags','rFlags','iFlags','zFlags',]][0])
                # sys.exit()
                self.ref_cat=np.array(self.ref_cat[self.stars])
                if len(self.stars)==0:
                    self.sp_logger.warning(warn_r+f' {self.sci_obj} {self.sci_mjd} {self.sci_filt}: Length of PS1 reference catalog is 0')
                    self.sys_exit=True
                    return
                self.ref_coords_wcs_sky = SkyCoord(ra=self.ref_cat['raMean']*u.deg, dec=self.ref_cat['decMean']*u.deg,frame='fk5')
                self.ref_coords_wcs=np.column_stack((self.ref_cat['raMean'],self.ref_cat['decMean']))
                self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)
                # self.sp_logger.info(self.ref_coords_wcs)

            if self.sci_filt=='u' or self.use_sdss==True:
                self.ref_cat=sdss_query(ra_deg=round(self.sci_c.ra.deg,6),dec_deg=round(self.sci_c.dec.deg,6), rad_deg=round((self.ref_width)*3,6))
                self.stars=self.ref_cat

                if len(self.stars)==0:
                    self.sp_logger.warning(warn_r+f' {self.sci_obj} {self.sci_mjd} {self.sci_filt}: Length of SDSS reference catalog is 0')
                    self.sys_exit=True
                    return

                self.ref_coords_wcs_sky = SkyCoord(ra=np.array(self.ref_cat['ra'])*u.deg, dec=np.array(self.ref_cat['dec'])*u.deg,frame='fk5')
                self.ref_coords_wcs=np.column_stack((self.ref_cat['ra'],self.ref_cat['dec']))
                self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)
        else:
            #If a reference catalogue is passed in with it's full path as self.auto_cat
            if self.auto_cat.startswith(self.path):
                self.auto_cat = self.path+self.auto_cat

            if not os.path.exists(self.auto_cat):
                self.sp_logger.info(warn_r+f' Reference catalogue {self.auto_ref} not found')
                self.sys_exit=True
                return

            # self.sp_logger.info(pd.read_csv(self.auto_cat,skiprows=1,header=None))
            self.cat_arr=  np.array(ascii.read(self.auto_cat,data_start=1,names=["filt","mag","magerr","xpos","ypos","ra","dec"]))
            self.ref_cat = pd.DataFrame(self.cat_arr,columns=["filt","mag","magerr","xpos","ypos","ra","dec"]) #mag,magerr,(pixel)xpos,(pixel)ypos,RA,DEC
            self.stars=self.ref_cat
            # self.ref_coords_wcs_sky = np.column_stack((self.ref_cat["RA"],self.ref_cat["DEC"]))
            
            self.ref_coords_wcs_sky = SkyCoord(ra=self.ref_cat["ra"]*u.deg,dec=self.ref_cat["dec"]*u.deg,frame='fk5')
            self.ref_coords_wcs = self.ref_coords_wcs_sky
            self.ref_coords_pix = np.column_stack((self.ref_cat["xpos"],self.ref_cat["ypos"]))
            
        self.ref_ali_wcs = WCS(self.ref_ali_name)  
        self.sp_logger.info(info_g+' Catalog stars in PS1/SDSS found='+str(len(self.stars))+','+str(round(3600*(self.ref_width/60),6))+ 'arcsec search radius')
        self.mean,self.median, self.std = sigma_clipped_stats(self.sci_ali_img_hdu.data, sigma=0.5,)
        # self.sp_logger.info(info_g+' Mean, Median, Std of sci_ali: '+str(round(self.mean,6))+' '+str(round(self.median,6))+' '+str(round(self.std,6)))
        # self.sp_logger.info(info_g+' Detecting stars with SExtractor')

        self.sp_logger.info(info_g+' Detecting stars with IRAFStarFinder')
        OUTEDGE = 0
        self.ref_coords_pix = np.column_stack((self.ref_coords_pix[:,0]-OUTEDGE,self.ref_coords_pix[:,1]-OUTEDGE))
        cond = (0<self.ref_coords_pix[:,0]) & (self.ref_coords_pix[:,0]<self.sci_ali_img_hdu.data.shape[0]-OUTEDGE) & (0<self.ref_coords_pix[:,1]) & (self.ref_coords_pix[:,1]<self.sci_ali_img_hdu.data.shape[1]-OUTEDGE)
        self.ref_coords_wcs,self.ref_coords_wcs_sky = self.ref_coords_wcs[cond],self.ref_coords_wcs_sky[cond]
        self.ref_coords_pix = self.ref_coords_pix[cond]
        # print(self.ref_coords_pix)
        self.iraffind1= IRAFStarFinder(threshold=0,fwhm=3.0,roundhi=0.3,min_separation=20.0,xycoords=self.ref_coords_pix)
        self.iraffind2= IRAFStarFinder(threshold=0,fwhm=3.0,roundhi=0.3,min_separation=20.0)
        # self.sp_logger.info(info_g+' Threshold for detecting stars: '+str(int(starscale*self.std)))
        self.sources1 = self.iraffind1(self.sci_ali_img_hdu.data[OUTEDGE:len(self.sci_ali_img_hdu.data)-OUTEDGE,OUTEDGE:len(self.sci_ali_img_hdu.data)-OUTEDGE] - self.median)
        self.sources2 = self.iraffind2(self.sci_ali_img_hdu.data[OUTEDGE:len(self.sci_ali_img_hdu.data)-OUTEDGE,OUTEDGE:len(self.sci_ali_img_hdu.data)-OUTEDGE] - self.median)

        self.sp_logger.info(info_g+f" Running SExtractor to detect stars in science image")
        print(sex_path + " " + self.sci_ali_name + " -c "+path+"config_files/align_sex.config -SATUR_LEVEL 50000 -BACK_TYPE MANUAL -BACK_VALUE "+str(self.median_bkg) +f" -CATALOG_NAME "+path+f"config_files/sci_dist_{self.rand_nums_string}.cat")
        os.system(sex_path + " " + self.sci_ali_name + " -c "+path+"config_files/align_sex.config -SATUR_LEVEL 50000 -BACK_TYPE MANUAL -BACK_VALUE "+str(self.median_bkg) +f" -CATALOG_NAME "+path+f"config_files/sci_dist_{self.rand_nums_string}.cat")
        self.sources3 = ascii.read(path+f"config_files/sci_dist_{self.rand_nums_string}.cat")#,
        
        self.sources=Table()
        self.sources['xcentroid'] = np.concatenate([self.sources3['X_IMAGE'].data,self.sources1['xcentroid'].data,self.sources2['xcentroid'].data],dtype=np.float64)
        self.sources['ycentroid'] = np.concatenate([self.sources3['Y_IMAGE'].data,self.sources1['ycentroid'].data,self.sources2['ycentroid'].data],dtype=np.float64)
        #set the dtype of the sources to ints

        #returns id, xcentroid, ycentroid, fwhm, sharpness, roundness, pa, npix, sky, peak, flux, mag
        self.sp_logger.info(info_g+' Found '+str(len(self.sources))+' stars in science image')

        
        self.star_coords_pix=np.column_stack((self.sources['xcentroid']+OUTEDGE,self.sources['ycentroid']+OUTEDGE))
        self.ref_coords_pix = np.column_stack((self.ref_coords_pix[:,0]+OUTEDGE,self.ref_coords_pix[:,1]+OUTEDGE))

        self.star_coords_wcs=load_wcs_from_file(filename=self.sci_ali_name,coord=self.star_coords_pix)

        self.star_coords_wcs_sky = SkyCoord(ra=self.star_coords_wcs[:,0]*u.deg, dec=self.star_coords_wcs[:,1]*u.deg,frame='fk5')
        self.scp = self.sources
        self.scp['xcentroid']=self.scp['xcentroid']+OUTEDGE
        self.scp['ycentroid']=self.scp['ycentroid']+OUTEDGE

        self.matched_catalog_mag=[]

        self.sci_ali_photTab = quick_app_phot(self.sci_ali_img_hdu.data,self.star_coords_pix,gain=self.sci_gain)
        self.sci_ali_photTab['ra']=self.star_coords_wcs_sky.ra
        self.sci_ali_photTab['dec']=self.star_coords_wcs_sky.dec
        # self.ref_ali_photTab['ra']=self.ref_coords_wcs_sky.ra
        # self.ref_ali_photTab['dec']=self.ref_coords_wcs_sky.dec

        cond = (np.isnan(self.sci_ali_photTab['SNR'])==False)&(self.sci_ali_photTab['SNR']>10)
        self.len_before_SNR=len(self.sci_ali_photTab)
        self.sci_ali_photTab=self.sci_ali_photTab[cond]
        self.sources=self.sources[cond]
        self.star_coords_pix=self.star_coords_pix[cond]
        self.star_coords_wcs=self.star_coords_wcs[cond]
        self.len_after_SNR=len(self.sci_ali_photTab)
        self.sp_logger.info(info_g+f' Keeping {self.len_after_SNR}/{self.len_before_SNR} stars with SNR>10 ({round(100*self.len_after_SNR/self.len_before_SNR,2)}%)')
    
        
        self.star_coords_wcs=load_wcs_from_file(filename=self.sci_ali_name,coord=self.star_coords_pix)
        self.star_coords_wcs_sky = SkyCoord(ra=self.star_coords_wcs[:,0]*u.deg, dec=self.star_coords_wcs[:,1]*u.deg,frame='fk5')

        #find crossover between panstarrs ad reference images
        self.sp_logger.info(info_g+' Searching for stars in reference catalog')
        self.indx, self.d2d, self.d3d =self.star_coords_wcs_sky.match_to_catalog_sky(self.ref_coords_wcs_sky)
        self.upd_indx=np.where(self.d2d<=1/3600.*u.deg)[0]
        d2d_ = self.d2d[self.upd_indx]

        if self.sci_filt=='g':m=1
        elif self.sci_filt=='r':m=0.5
        elif self.sci_filt=='i':m=0.1
        else:m=3
        self.upd_indx=np.where(self.d2d<=m/3600.*u.deg)[0]
        d2d_ = self.d2d[self.upd_indx]
        d2d_ = d2d_.to(u.arcsec)
        self.sp_logger.info(info_g+' Catalog stars in PS1/SDSS found = '+str(len(self.upd_indx))+', in a '+str(round(3600*(m/60),6))+ ' arcsec search radius')


    
        self.matched_ref_coords_pix=self.ref_coords_pix[self.indx[self.upd_indx]]
        self.matched_star_coords_pix=self.star_coords_pix[self.upd_indx]

        self.ref_coords_wcs_sky=self.ref_coords_wcs_sky[self.indx[self.upd_indx]]
        self.star_coords_wcs_sky=self.star_coords_wcs_sky[self.upd_indx]

        # print(len(self.matched_star_coords_pix))
        # print(len(self.matched_ref_coords_pix))
        # sys.exit()
        # save_to_reg=False
        # if save_to_reg:
        #     if not os.path.exists(self.path+'region_files'):os.makedirs(self.path+'region_files')
        #     self.files_to_clean.append(self.path+'region_files/'+self.sci_obj+'_sci_all_sky.reg')
        #     self.files_to_clean.append(self.path+'region_files/'+self.sci_obj+'_ref_all_sky.reg')
        #     with open(self.path+'region_files/'+self.sci_obj+'_sci_all_sky.reg','w') as f:
        #         f.write('global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
        #         f.write('fk5'+'\n')
        #         for i in range(len(self.matched_star_coords_pix)):
        #             f.write(f'circle({self.star_coords_wcs_sky[i].ra.deg},{self.star_coords_wcs_sky[i].dec.deg},3")'+' # text={'+f'xy={int(self.matched_star_coords_pix[i][0]),int(self.matched_star_coords_pix[i][1])}'+'}'+'\n')
        #             # f.write(f'physical;circle({self.matched_star_coords_pix[i][0]},{self.matched_star_coords_pix[i][1]},10)'+' # text={'+f'xy={int(self.matched_star_coords_pix[i][0]),int(self.matched_star_coords_pix[i][1])}'+'}'+'\n')
        #         f.close()
        #     with open(self.path+'region_files/'+self.sci_obj+'_ref_all_sky.reg','w') as f:
        #         f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
        #         f.write('fk5'+'\n')
        #         for i in range(len(self.matched_ref_coords_pix)):
        #             f.write(f'circle({self.ref_coords_wcs_sky[i].ra.deg},{self.ref_coords_wcs_sky[i].dec.deg},3")'+' # text={'+f'xy={int(self.matched_ref_coords_pix[i][0]),int(self.matched_ref_coords_pix[i][1])}'+'}'+'\n')
        #             # f.write(f'physical;circle({self.matched_ref_coords_pix[i][0]},{self.matched_ref_coords_pix[i][1]},10)'+' # text={'+f'xy={int(self.matched_ref_coords_pix[i][0]),int(self.matched_ref_coords_pix[i][1])}'+'}'+'\n')

        #         f.close()

        # sys.exit()

        # [print(f'circle({self.ref_coords_wcs_sky[i].ra.deg},{self.ref_coords_wcs_sky[i].dec.deg},3")'+' # text={'+f'xy={self.matched_ref_coords_pix[i]}'+'}') for i in range(len(self.ref_coords_wcs_sky))]
        # print()
        # [print(f'circle({self.star_coords_wcs_sky[i].ra.deg},{self.star_coords_wcs_sky[i].dec.deg},3")'+' # text={'+f'xy={self.matched_star_coords_pix[i]}'+'}') for i in range(len(self.ref_coords_wcs_sky))]


        self.sci_ali_photTab = quick_app_phot(self.sci_ali_img_hdu.data,self.matched_star_coords_pix,gain=self.sci_gain)
        self.ref_ali_photTab = quick_app_phot(self.ref_ali_img_hdu.data,self.matched_ref_coords_pix,4)

        self.sci_ali_photTab['ra']=self.star_coords_wcs_sky.ra
        self.sci_ali_photTab['dec']=self.star_coords_wcs_sky.dec
        self.ref_ali_photTab['ra']=self.ref_coords_wcs_sky.ra
        self.ref_ali_photTab['dec']=self.ref_coords_wcs_sky.dec



        self.uniq_inds = np.unique(self.matched_ref_coords_pix, axis=0, return_index=True)[1]
        self.matched_star_coords_pix = self.matched_star_coords_pix[self.uniq_inds]
        self.matched_ref_coords_pix = self.matched_ref_coords_pix[self.uniq_inds]
        self.star_coords_wcs_sky = self.star_coords_wcs_sky[self.uniq_inds]
        self.ref_coords_wcs_sky = self.ref_coords_wcs_sky[self.uniq_inds]

        bord_dists = self.find_dist_to_border(sci_data=self.sci_ali_img_hdu.data,
                                                ref_data=self.ref_ali_img_hdu.data,
                                                sci_match_coords=self.matched_star_coords_pix,
                                                ref_match_coords=self.matched_ref_coords_pix,
                                                sci_phot_tab=self.sci_ali_photTab,
                                                ref_phot_tab=self.ref_ali_photTab,
                                                X=20,D=30) 
        #        return {'sci_keep_pix':sci_keep_pix,'ref_keep_pix':ref_keep_pix,'sci_keep_sky':sci_keep_sky,'ref_keep_sky':ref_keep_sky}
        self.orig_len_match = len(self.matched_star_coords_pix)

        if len(self.matched_star_coords_pix)>5:
            self.matched_star_coords_pix = bord_dists['sci_keep_pix']
            self.matched_ref_coords_pix = bord_dists['ref_keep_pix']

            self.star_coords_wcs_sky = bord_dists['sci_keep_sky']
            self.ref_coords_wcs_sky = bord_dists['ref_keep_sky']
            # print(self.star_coords_wcs_sky)
            self.new_len_match = len(self.matched_star_coords_pix)
            self.sp_logger.info(info_g+f' Kept {self.new_len_match} stars out of {self.orig_len_match} ({self.new_len_match/self.orig_len_match*100:.2f}%)')

        else:
            self.sp_logger.info(info_g+f" Not rejecting any stars to allow correction for distortion")





        if self.auto_cat == 'auto' and self.use_sdss==False and (self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z'):
            self.matched_catalog=self.ref_cat[self.indx[self.upd_indx]]
            self.matched_catalog_mag=self.matched_catalog[str(self.sci_filt)+'MeanPSFMag']

        elif self.auto_cat == 'auto' and (self.sci_filt=='u' or self.use_sdss==True):
            self.string_band='psfMag_'+str(self.sci_filt)
            self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
            self.matched_catalog_mag=np.asarray(self.matched_catalog['mag'])
            self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

        elif self.auto_cat !='auto':
            # self.sp_logger.info(self.ref_cat.loc[self.indx[self.upd_indx]])
            self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
            self.matched_catalog_mag=np.asarray([float(m) for m in self.matched_catalog['mag']])
            self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

    

        save_to_reg = False
        if save_to_reg: 
            if not os.path.exists(self.path+'region_files'):os.makedirs(self.path+'region_files')
            self.files_to_clean.append(self.path+'region_files/'+self.sci_obj+'_sci.reg')
            self.files_to_clean.append(self.path+'region_files/'+self.sci_obj+'_ref.reg')
            with open(self.path+'region_files/'+self.sci_obj+'_sci.reg','w') as f:
                f.write('global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
                f.write('fk5'+'\n')
                for i in range(len(self.matched_star_coords_pix)):
                    # f.write(f'circle({self.sci_keep_sky[i].ra.deg},{self.sci_keep_sky[i].dec.deg},3")'+' # text={'+f'xy={self.sci_keep_pix[i]}'+'}'+'\n')
                    f.write(f'physical;circle({self.matched_star_coords_pix[i][0]},{self.matched_star_coords_pix[i][1]},20)'+' # text={'+f'xy={self.matched_star_coords_pix[i]}'+'}'+'\n')
                f.close()

            with open(self.path+'region_files/'+self.sci_obj+'_ref.reg','w') as f:
                f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
                f.write('fk5'+'\n')
                for i in range(len(self.matched_ref_coords_pix)):
                    # f.write(f'circle({self.ref_keep_sky[i].ra.deg},{self.ref_keep_sky[i].dec.deg},3")'+' # text={'+f'xy={self.ref_keep_pix[i]}'+'}'+'\n')
                    f.write(f'physical;circle({self.matched_ref_coords_pix[i][0]},{self.matched_ref_coords_pix[i][1]},20)'+' # text={'+f'xy={self.matched_ref_coords_pix[i]}'+'}'+'\n')
                f.close()

        
        # sys.exit()
        # for i in range(len(self.matched_star_coords_pix)):
            # print(self.matched_star_coords_pix[i],self.matched_ref_coords_pix[i])
        # from skimage.metrics import mean_squared_error
        # print(len(self.matched_star_coords_pix))
        # for i in range(len(self.matched_star_coords_pix)):
        #     print(self.matched_star_coords_pix[i],self.matched_ref_coords_pix[i])

        # print()
        # for i in range(len(self.uniq_inds)):
        #     print(self.matched_star_coords_pix[self.uniq_inds[i]],self.matched_ref_coords_pix[self.uniq_inds[i]])
        
        if len(self.matched_star_coords_pix)>=3:
            self.sp_logger.info(info_g+f" Finding the transformation between science and reference to correct for distortion")
            self.tform_ref_sci,self.tform_sci_ref = transform.PolynomialTransform(),transform.PolynomialTransform()

            # print(len(self.matched_star_coords_pix),len(self.matched_ref_coords_pix))
            
            # sys.exit()
            self.tform_ref_sci.estimate(self.matched_ref_coords_pix,self.matched_star_coords_pix, order=1)
            self.tform_ref_sci = transform.warp(self.sci_ali_img_hdu.data, self.tform_ref_sci,
                                                output_shape=(self.sci_ali_img_hdu.data.shape[1],self.sci_ali_img_hdu.data.shape[0]),
                                                order=3,
                                                mode='constant', cval=0.0)

            # self.tform_sci_ref.estimate(self.matched_star_coords_pix,self.matched_ref_coords_pix, order=1)
            # self.tform_sci_ref = transform.warp(self.sci_ali_img_hdu.data, self.tform_sci_ref,
            #                                     output_shape=(self.sci_ali_img_hdu.data.shape[1],self.sci_ali_img_hdu.data.shape[0]),
            #                                     order=3,)
            self.align_success = True
            # new_sci_ali = fits.open(self.sci_ali_name)
            # new_sci_ali[0].data = self.tform_ref_sci
            # new_sci_ali.writeto(self.sci_ali_name.replace('.fits',f'_ref_sci_warped_order{1}.fits'), 
            #                     overwrite=True)
            # self.sp_logger(info_g+f" Warped science using polynomial order 1")


    
            # self.sp_logger(info_g+f" Warped science using polynomial order 2")
            self.align_success = True
            new_sci_ali = fits.open(self.sci_ali_name)
            new_sci_ali[0].data = self.tform_ref_sci
            new_sci_ali[0].header = self.ref_ali_img_hdu.header
            new_sci_ali.writeto(self.sci_ali_name.replace('.fits',f'_ref_sci_warped_order{1}_snum{len(self.matched_star_coords_pix)}.fits'), overwrite=True)

            self.sci_ali_name = self.sci_ali_name.replace('.fits',f'_ref_sci_warped_order{1}_snum{len(self.matched_star_coords_pix)}.fits')
            self.files_to_clean.append(self.sci_ali_name)
            self.sp_logger.info(info_g+f" Warped science saved to {self.sci_ali_name}")
            
            def add_circles(ax, coords, color='yellow'):
                for x, y in coords:
                    circle = plt.Circle((x, y), 10, color=color, fill=False)
                    ax.add_artist(circle)

            self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
            self.rvmin,self.rvmax = visualization.ZScaleInterval().get_limits(self.ref_ali_img_hdu.data)
            self.rvmin,self.rvmax = np.nanpercentile(self.ref_ali_img_hdu.data,[5,95])
            # fig_poly_coll, axes = plt.subplots(2, 2, figsize=(12, 12),)
            
            import matplotlib.gridspec as gridspec
            fig_poly_coll = plt.figure(figsize=(12, 12))
            gs = gridspec.GridSpec(2, 4)
            gs.update(wspace=0.5)
            a1 = plt.subplot(gs[0, :2], )
            a2 = plt.subplot(gs[0, 2:])
            a3 = plt.subplot(gs[1, 1:3])
            # plt.show()
            

            self.sci_ali_sn_coords = wcs_to_pixels(self.sci_ali_name,np.column_stack((self.sci_c.ra.value,self.sci_c.dec.value)))[0]
            self.ref_ali_sn_coords = wcs_to_pixels(self.ref_ali_name,np.column_stack((self.sci_c.ra.value,self.sci_c.dec.value)))[0]

            a1.imshow(self.sci_ali_img_hdu.data,cmap='gray',vmin=self.vmin,vmax=self.vmax)
            #plot the sn position
            a1.add_artist(plt.Circle((self.sci_ali_sn_coords[0],self.sci_ali_sn_coords[1]),20,color='red',fill=False))
            a1.set_title('Original Science')
            add_circles(a1,self.matched_ref_coords_pix)
            
            a2.imshow(self.ref_ali_img_hdu.data,cmap='gray',vmin=self.rvmin,vmax=self.rvmax)
            a2.add_artist(plt.Circle((self.ref_ali_sn_coords[0],self.ref_ali_sn_coords[1]),20,color='red',fill=False))
            a2.set_title('Original Reference')
            add_circles(a2,self.matched_ref_coords_pix)

            a3.imshow(self.tform_ref_sci,cmap='gray',vmin=self.vmin,vmax=self.vmax)
            a3.add_artist(plt.Circle((self.ref_ali_sn_coords[0],self.ref_ali_sn_coords[1]),20,color='red',fill=False))
            a3.set_title('Warped Science')
            add_circles(a3,self.matched_ref_coords_pix)
            
            #remove axis ticks
            a1.set_xticks([]),a1.set_yticks([])
            a2.set_xticks([]),a2.set_yticks([])
            a3.set_xticks([]),a3.set_yticks([])
            # a4.imshow(self.tform_sci_ref,cmap='gray',vmin=self.vmin,vmax=self.vmax)
            # a4.set_title('Sci Ref Trans ')
            # add_circles(a4,self.matched_star_coords_pix)



            fig_sci_ali,axes_sci = plt.subplots(figsize=(12,12))
            axes_sci.imshow(self.sci_ali_img_hdu.data,cmap='gray',vmin=self.vmin,vmax=self.vmax)
            add_circles(axes_sci,self.matched_ref_coords_pix)

            fig_ref_ali,axes_ref = plt.subplots(figsize=(12,12))
            axes_ref.imshow(self.ref_ali_img_hdu.data,cmap='gray',vmin=self.rvmin,vmax=self.rvmax)
            add_circles(axes_ref,self.matched_ref_coords_pix)



                
                

            
            


            if not os.path.exists(self.path+'poly_comps/'):os.mkdir(self.path+'poly_comps/')
            fig_poly_coll.savefig(self.path+'poly_comps/'+self.sci_img_name[:-11]+'_poly_collage.pdf',bbox_inches='tight')
            self.sp_logger.info(info_g+f" Saved image comparing polynomial orders 1 and 2 as "+self.path+'poly_comps/'+self.sci_img_name[:-11]+'_poly_collage.pdf') 
            # fig_sci_ali.savefig(self.path+'sedm_comps2/'+self.sci_img_name[:-11]+'_sci.pdf')
            # fig_ref_ali.savefig(self.path+'sedm_comps2/'+self.sci_img_name[:-11]+'_ref.pdf')


        else:
            self.sp_logger.info(warn_y+f" Not enough stars matched for distorion correction ")
        # sys.exit(1)
        return



    def swarp_ref_align(self,image_size=image_size):

        # prepsexfile(gain=self.sci_gain)

        if self.telescope in SEDM:
            self.image_size=1000
        else:
            self.image_size=1500

        if self.sci_obj=='ZTF25aasjeza':# and self.sci_filt=='z':
            self.image_size=900
        if self.sci_obj=='ZTF24abdiwwv': #02:22:10.96 -20:23:21.01
            # self.image_size=1500
            self.image_size=1650
            # self.image_size=1300

        # [override] --force_image_size N : use the user-supplied size as the
        # target SWarp frame and pad both sci & ref to it in the mismatch
        # branch below, instead of shrinking to whichever resamp came out
        # smaller (which can crop the transient out when the target is far
        # from the science image centre).
        try:
            _fis = getattr(self.args, 'force_image_size', None)
        except Exception:
            _fis = None
        if _fis is not None and int(_fis) > 0:
            self.image_size = int(_fis)
            self.sp_logger.info(info_g + f' [force] --force_image_size override: '
                                         f'using IMAGE_SIZE = {self.image_size} px')
        self._force_image_size_active = (_fis is not None and int(_fis) > 0)
        self._initial_image_size = self.image_size

        self.sp_logger.info(info_g+f' Aligning science with reference image')

        if self.auto_ref=='auto':
            if self.sci_filt=='u' or self.use_sdss==True: #if filt ==u, need to call make sdss_ref before this or can call it here
                self.make_sdss_ref(sdss_filt=self.sci_filt)
                if self.sys_exit==True:
                    return

                self.image_name=self.sci_obj+f'_ref_{self.sci_filt}.fits'
                self.ref_path=self.path+'ref_imgs/'+self.image_name
                self.ref_img_name=self.ref_path
                self.sp_logger.info(info_g+f' Using SDSS reference image for alignment: '+self.ref_img_name)
                # self.sp_logger.info(self.ref_path)
                # sys.exit(1)
                panstamps_status=100 #no need to worry


            elif self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z':
                self.ref_folder = self.path+'ref_imgs/stack_'+str(self.sci_filt)
                self.ref_path = self.path+'ref_imgs/stack_'+str(self.sci_filt)+'_ra'+self.ra_string+'_dec'+self.dec_string+'_arcsec*'
                # sys.exit()
                if len(glob.glob(self.ref_path))>0:
                    self.sp_logger.info(info_b+' PS1 reference image already exists: '+glob.glob(self.ref_path)[0])
                    panstamps_status=0
                if len(glob.glob(self.ref_path))==0:
                    self.sp_logger.info(info_g+f' Running panstamps: '+str(panstamps_path+' -f --width='+str(self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg)))
                    self.sp_logger.info(info_g+' Downloading reference image from PS...')
                    try:
                        panstamps_status = os.system(panstamps_path+' -f --width='+str(self.ref_width*3)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))


                        self.sp_logger.info(info_g+' Reference image download status '+str(panstamps_status))
                        self.sp_logger.info(info_g+' Reference image name: '+glob.glob(self.ref_path)[0])

                    except Exception as e:
                        self.sp_logger.warning('error'+str(e))
                        self.sp_logger.warning(warn_r+f" {self.sci_obj}, {self.sci_mjd}, {self.sci_filt}: Unable to download reference image from panstarrs")

                        self.sp_logger.warning(warn_r+f" Error downloading PS reference image, error encountered :{e}")
                        self.sp_logger.warning(warn_r+f" Exiting...")
                        ##sys.exit(1)
                        self.sys_exit=True
                        return
                
            if self.sys_exit==True:
                return

            # self.sp_logger.info(self.ref_path)
            if self.use_sdss==False and panstamps_status==0:
                c=0
                while len(glob.glob(self.ref_path))==0:
                    #the image has downloaded but the name is not the same as the one we want so we need to find it by shortening the name we search for
                    self.ref_path=self.ref_path[:-2]+'*'
                    # self.sp_logger.info(self.ref_path)
                    # self.sp_logger.info(glob.glob(self.ref_path))
                    c+=1
                    if c>10:
                        break
            # self.sp_logger.info(glob.glob(self.ref_path))
            # self.sp_logger.info(self.ref_path)
            self.ref_img_name=glob.glob(self.ref_path)[0]
            self.ref_img_hdu=fits.open(self.ref_img_name)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
        
        else:
            #If a reference image is passed in with it's full path as self.auto_ref
            if not self.auto_ref.startswith(self.path):
                self.auto_ref = self.path+self.auto_ref

            if not os.path.exists(self.auto_ref):
                self.sp_logger.warning(warn_r+f' Reference image {self.auto_ref} not found')
                self.sp_logger.warning(warn_r+f" {self.sci_obj}, {self.sci_mjd}, {self.sci_filt}: Unable to find reference image from path specified {self.auto_ref}")

                self.sys_exit=True
                return

            self.ref_img_name=glob.glob(self.auto_ref)[0]
            # try:
            self.ref_img_hdu=fits.open(self.auto_ref)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
            # except:
            #     # print(self.auto_ref)
            #     try:
            #         self.ref_img_hdu=fits.open(self.auto_ref)[1]
            #         print(self.ref_img_hdu.header['NAXIS1'],self.ref_img_hdu.header['NAXIS2'])
            #         self.ref_img=self.ref_img_hdu.data#[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
            #         print(np.shape(self.ref_img))
            #     except:
            #         self.ref_img_hdu=fits.open(self.auto_ref)[1]
            #         print(self.ref_img_hdu.header['NAXIS1'],self.ref_img_hdu.header['NAXIS2'])
            #         self.ref_img=fits.open(self.auto_ref)[0].data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
            #         print(np.shape(self.ref_img))

                # sys.exit(1)
                # self.sp_logger.info(self.ref_img==self.sci_bkgsb)

        self.coords_sn_ref=wcs_to_pixels(self.ref_img_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
        self.coords_sn_ref_x,self.coords_sn_ref_y = self.coords_sn_ref




        if not os.path.exists(self.path+'aligned_images'):
            os.makedirs(self.path+'aligned_images')

        
        wcs_command='python3 '+self.path+'sedm_align_quick.py'+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+"  -m relative -r 100"
        # wcs_command='python3 '+self.path+'sedm_align_quick.py'+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+"  -m relative -r 30 -a "+self.ra_string+" "+self.dec_string
        # self.sp_logger.info(wcs_command)

        self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
        self.ref_img_hdu=fits.open(self.ref_img_name)[0]
        self.align_success=False
        self.align_fail_count=0
        self.valid_mask=None      # True where science aligned image has real data (not padding)
        self.ref_valid_mask=None  # True where reference aligned image has real data

        
        os.system(wcs_command)
        self.sp_logger.info(info_g+" Original science size: "+str(np.shape(self.sci_img_hdu.data)))
        self.sp_logger.info(info_g+" Original reference size: "+str(np.shape(self.ref_img_hdu.data)))
        # sys.exit(1)
        # print('hi')
        # sys.exit(1)

        while self.align_success==False:
            self.align_fail_count+=1



            self.sp_logger.info(info_g+" Aligning images with SWarp making a resampled image of size "+str(self.image_size)+"x"+str(self.image_size)+" pixels")

            # self.ali_center_ra,self.ali_center_dec = self.ra_string,self.dec_string
            #convert the center of the image to degrees
            # self.sci_shape= np.shape(self.sci_img_hdu.data)
            # self.ra_center,self.dec_center = self.sci_shape[0]/2,self.sci_shape[1]/2
            # if float(self.dec_string)>0:self.ali_center_dec='+'+self.dec_string
            # elif self.dec_string[0]=='-': self.ali_center_dec=self.dec_string
            # else:self.ali_center_dec='-'+self.dec_string

            self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
            self.ref_img_hdu=fits.open(self.ref_img_name)[0]

            #find centre of original reference image
            self.ref_img_size = np.shape(self.ref_img_hdu.data)
            self.ref_img=self.ref_img_hdu.data
            self.ref_img_center = np.shape(self.ref_img)[0]/2,np.shape(self.ref_img)[1]/2
            self.ref_wcs = WCS(self.ref_img_hdu.header)
            self.ref_cnt = self.ref_wcs.all_pix2world(self.ref_img_center[0],self.ref_img_center[1],1)
            self.ref_img_center_ra,self.ref_img_center_dec = self.ref_cnt[0],self.ref_cnt[1]

            self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)
            self.sci_img_size = np.shape(self.sci_img_hdu[0].data)
            self.sci_img=self.sci_img_hdu[0].data
            self.sci_img_center = np.shape(self.sci_img)[0]/2,np.shape(self.sci_img)[1]/2
            self.sci_wcs = WCS(self.sci_img_hdu[0].header)
            self.sci_cnt = self.sci_wcs.all_pix2world(self.sci_img_center[0],self.sci_img_center[1],1)
            self.sci_img_center_ra,self.sci_img_center_dec = self.sci_cnt[0],self.sci_cnt[1]
            self.ali_center_ra = str(self.sci_img_center_ra)
            self.ali_center_dec = str(self.sci_img_center_dec)
            if float(self.ali_center_dec) > 0: self.ali_center_dec = '+' + self.ali_center_dec
            self.sp_logger.info(info_g+" Aligning images to a centre of "+self.ali_center_ra+","+self.ali_center_dec)


            self.image_size1,self.image_size2=np.shape(self.sci_img_hdu[0].data)[0],np.shape(self.sci_img_hdu[0].data)[1]
            # self.image_size1,self.image_size2=np.shape(self.ref_img_hdu.data)[0],np.shape(self.ref_img_hdu.data)[1]

            swarp_command=swarp_path+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+" -c "+self.path+"config_files/config.swarp -CENTER '"+str(self.ali_center_ra)+" "+str(self.ali_center_dec)+"' -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -RESAMPLE Y -RESAMPLE_DIR '"+self.path+"aligned_images/' -COMBINE N -IMAGE_SIZE '"+str(self.image_size)+","+str(self.image_size)+"'"
            status=os.system(swarp_command)
            self.sp_logger.info(swarp_command)
            # self.sp_logger.info(info_g+" SWarp took %2f seconds" % (time.time() - swarp_time)    )
            self.sp_logger.info(info_g+" SWarp status: "+str(status))

            



            
            self.sci_ali_name=glob.glob(self.path+'aligned_images/'+self.sci_img_name[:-5]+'.resamp.fits')[0]
            

            if self.auto_ref=="auto" and self.sci_filt in ['g','r','i','z'] and self.use_sdss==False:
                self.ref_ali_name=self.path+'aligned_images/stack_'+str(self.sci_filt)+'_ra'+self.ra_string+'_dec'+self.dec_string+'_arcsec*.resamp.fits'
                try:
                    self.ref_ali_name=glob.glob(self.ref_ali_name)[0]
                except:
                    self.ref_ali_name=self.ref_img_name.replace('.fits','.resamp.fits')
                    self.ref_ali_name=self.ref_ali_name.replace('ref_imgs','aligned_images')

            elif self.auto_ref=="auto" and (self.sci_filt=='u' or self.use_sdss==True):
                # self.ref_ali_name=glob.glob(self.path+'aligned_images/'+self.sci_img_name[:-5]+'.resamp.fits')[0]
                self.ref_ali_name=self.ref_img_name.replace('.fits','.resamp.fits')
                self.ref_ali_name=self.ref_ali_name.replace('ref_imgs','aligned_images')

            elif self.auto_ref!="auto":
                self.ref_ali_name = re.sub(self.path,'',self.auto_ref[:-5])
                self.ref_ali_name = self.path+'aligned_images/'+re.sub('ref_imgs/','',self.ref_ali_name)+'.resamp.fits'

            if self.align_fail_count==1:
                self.sp_logger.info(info_g+' Alignment attempt 1')
            else: 
                self.sp_logger.info(warn_y+' Alignment attempt '+str(self.align_fail_count))

            self.ref_ali_img_hdu = fits.open(self.ref_ali_name)[0]
            self.sci_ali_img_hdu = fits.open(self.sci_ali_name)[0]
            self.sp_logger.info(info_g+" Aligned science image: "+self.sci_ali_name+' '+str(np.shape(self.sci_ali_img_hdu.data)))
            self.sp_logger.info(info_g+" Aligned reference image: "+self.ref_ali_name+' '+str(np.shape(self.ref_ali_img_hdu.data)))
            


            self.sp_logger.info(info_g+' Science & Reference image sizes : '+str(np.shape(fits.open(self.sci_ali_name)[0].data))+' '+str(np.shape(fits.open(self.ref_ali_name)[0].data)))
            self.sp_logger.info(info_b+' Target size : '+str(float(self.image_size))+'x'+str(float(self.image_size)))
            self.sp_logger.info(info_g+' Aligned Science image size : '+str(np.shape(fits.open(self.sci_ali_name)[0].data)))
            self.sp_logger.info(info_g+' Aligned Reference image size : '+str(np.shape(fits.open(self.ref_ali_name)[0].data)))

            # sys.exit()
            
            try_aa_instead = False
            if self.align_fail_count>15:
                self.sp_logger.warning(warn_y+' Alignment failed 15 times')
                break

            if np.shape(fits.open(self.sci_ali_name)[0].data)!=np.shape(fits.open(self.ref_ali_name)[0].data) and ((self.image_size,self.image_size)!=np.shape(fits.open(self.sci_ali_name)[0].data) or (self.image_size,self.image_size)!=np.shape(fits.open(self.ref_ali_name)[0].data)):
                self.sp_logger.warning(warn_y+' Science image or reference image and aligned image have different sizes. This may be due to the reference image being too far away from the science image')

                # [force_image_size] branch: pad BOTH sci and ref to the
                # original requested IMAGE_SIZE using their respective COMIN
                # offsets, instead of shrinking to the smaller resamp shape
                # (which crops the transient out when the target is offset
                # from the science image centre).
                if getattr(self, '_force_image_size_active', False):
                    try:
                        _target = int(self._initial_image_size)
                        _sci_hdu = fits.open(self.sci_ali_name)[0]
                        _ref_hdu = fits.open(self.ref_ali_name)[0]
                        _scx = int(_sci_hdu.header.get('COMIN1', 1))
                        _scy = int(_sci_hdu.header.get('COMIN2', 1))
                        _rcx = int(_ref_hdu.header.get('COMIN1', 1))
                        _rcy = int(_ref_hdu.header.get('COMIN2', 1))

                        def _pad_to(data, cx, cy, target):
                            out = np.zeros((target, target), dtype=data.dtype)
                            y0, x0 = cy - 1, cx - 1
                            y1, x1 = y0 + data.shape[0], x0 + data.shape[1]
                            sy0, sy1 = max(0, -y0), min(data.shape[0], target - y0)
                            sx0, sx1 = max(0, -x0), min(data.shape[1], target - x0)
                            dy0, dy1 = max(0, y0), min(target, y1)
                            dx0, dx1 = max(0, x0), min(target, x1)
                            if sy1 > sy0 and sx1 > sx0:
                                out[dy0:dy1, dx0:dx1] = data[sy0:sy1, sx0:sx1]
                            return out

                        _padded_sci = _pad_to(_sci_hdu.data, _scx, _scy, _target)
                        _padded_ref = _pad_to(_ref_hdu.data, _rcx, _rcy, _target)

                        # New WCS = SWarp combined-frame WCS. Each .resamp header
                        # carries the projection of the COMBINEd output but its
                        # CRPIX is relative to the cropped resamp. Add COMIN-1
                        # to recover the combined-frame CRPIX.
                        _new_hdr = _sci_hdu.header.copy()
                        _new_hdr['CRPIX1'] = float(_new_hdr['CRPIX1']) + (_scx - 1)
                        _new_hdr['CRPIX2'] = float(_new_hdr['CRPIX2']) + (_scy - 1)
                        _new_hdr['NAXIS1'] = _target
                        _new_hdr['NAXIS2'] = _target

                        # Write padded files
                        _sci_padded_name = self.path + 'aligned_images/' + \
                            self.sci_img_name.split('/')[-1].replace('.fits', '_padded.fits')
                        _ref_padded_name = self.path + 'aligned_images/' + \
                            os.path.basename(self.ref_ali_name).replace('.resamp.fits', '_padded.fits')
                        fits.writeto(_sci_padded_name, _padded_sci, header=_new_hdr, overwrite=True)
                        fits.writeto(_ref_padded_name, _padded_ref, header=_new_hdr, overwrite=True)

                        # Wire them in as the official aligned outputs
                        self.sci_ali_img_hdu = fits.open(_sci_padded_name)[0]
                        self.ref_ali_img_hdu = fits.open(_ref_padded_name)[0]
                        self.sci_ali_name = _sci_padded_name
                        self.ref_ali_name = _ref_padded_name

                        # Valid masks (track where each input actually had data)
                        self.valid_mask = np.zeros((_target, _target), dtype=bool)
                        _y0s, _y1s = max(0, _scy-1), min(_target, _scy-1 + _sci_hdu.data.shape[0])
                        _x0s, _x1s = max(0, _scx-1), min(_target, _scx-1 + _sci_hdu.data.shape[1])
                        self.valid_mask[_y0s:_y1s, _x0s:_x1s] = True
                        self.ref_valid_mask = np.zeros((_target, _target), dtype=bool)
                        _y0r, _y1r = max(0, _rcy-1), min(_target, _rcy-1 + _ref_hdu.data.shape[0])
                        _x0r, _x1r = max(0, _rcx-1), min(_target, _rcx-1 + _ref_hdu.data.shape[1])
                        self.ref_valid_mask[_y0r:_y1r, _x0r:_x1r] = True

                        self.sp_logger.info(info_g + f' [force] Padded sci+ref to '
                                                     f'common {_target}x{_target} frame: '
                                                     f'sci@COMIN=({_scx},{_scy}) '
                                                     f'ref@COMIN=({_rcx},{_rcy})')
                        self.sp_logger.info(info_g + f' [force] Sci valid: '
                                                     f'{100*self.valid_mask.mean():.1f}% | '
                                                     f'Ref valid: {100*self.ref_valid_mask.mean():.1f}%')
                        self.image_size = _target
                        self.align_success = True
                        continue
                    except Exception as _e_force:
                        self.sp_logger.warning(warn_y + f' [force] pad-to-target failed '
                                                       f'({_e_force}); falling back to shrink-and-retry')

                self.min_align_size = np.min([np.shape(fits.open(self.sci_ali_name)[0].data),np.shape(fits.open(self.ref_ali_name)[0].data)])
                self.sp_logger.warning(warn_y+' Trying image size of '+str(self.min_align_size)+'x'+str(self.min_align_size)+' pixels')
                self.image_size=self.min_align_size
                self.align_success=False

                try:
                    # if self.telescope in SEDM:
                    #find centre of original reference image
                    self.ref_img_hdu=fits.open(self.ref_img_name)[0]
                    self.ref_img_size = np.shape(self.ref_img_hdu.data)
                    self.ref_img=self.ref_img_hdu.data
                    # self.ref_img_center = np.shape(self.ref_img)[0]/2,np.shape(self.ref_img)[1]/2
                    self.ref_wcs = WCS(self.ref_img_hdu.header)
                    # self.ref_cnt = self.ref_wcs.all_pix2world(self.ref_img_center[0],self.ref_img_center[1],1)
                    # self.ref_img_center_ra,self.ref_img_center_dec = self.ref_cnt[0],self.ref_cnt[1]

                    self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
                    self.sci_img_size = np.shape(self.sci_img_hdu.data)
                    self.sci_img=self.sci_img_hdu.data
                    # self.sci_img_center = np.shape(self.sci_img)[0]/2,np.shape(self.sci_img)[1]/2
                    # self.sci_wcs = WCS(self.sci_img_hdu.header)
                    # self.sci_cnt = self.sci_wcs.all_pix2world(self.sci_img_center[0],self.sci_img_center[1],1)
                    # self.sci_img_center_ra,self.sci_img_center_dec = self.sci_cnt[0],self.sci_cnt[1]

                    self.sci_ali_img_hdu = fits.open(self.sci_ali_name)[0]
                    self.sci_ali_img_size = np.shape(self.sci_ali_img_hdu.data)
                    # self.sci_path,_ = self.update_astrometry(self.sci_path,ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])

                    # sys.exit()
                    requested_size = np.shape(self.ref_ali_img_hdu.data)
                    # requested_size = (self.image_size1,self.image_size2)
                    # self.sp_logger.info('requested_size',requested_size)
                    self.sp_logger.info(info_g+' Padding science image to match reference image size of '+str(requested_size[0])+'x'+str(requested_size[1])+' pixels')
                    cominx,cominy = int(self.sci_ali_img_hdu.header['COMIN1']),int(self.sci_ali_img_hdu.header['COMIN2'])
                    try:
                        self.sp_logger.info(info_g+' Science COMIN1,COMIN2 : '+str(cominx)+' '+str(cominy))
                        pad_sci = np.zeros((requested_size[0],requested_size[1]))#*np.nan
                        pad_sci[cominy:cominy+self.sci_ali_img_size[0],cominx:cominx+self.sci_ali_img_size[1]] = self.sci_ali_img_hdu.data
                        # pad_ref[cominy:cominy+self.ref_ali_img_size[0],cominx:cominx+self.ref_ali_img_size[1]] = self.ref_ali_img_hdu.data
                        new_sci_img = pad_sci
                        self.sp_logger.warning(info_g+' Successfully padded science image')
                    except Exception as e:
                        try:
                            self.sp_logger.warning(warn_y+f" Error padding science image, error encountered, changing cominy+=1 {e}")
                            pad_sci = np.zeros((requested_size[0],requested_size[1]))#*np.nan
                            cominy-=1
                            cominx-=1
                            pad_sci[cominy:cominy+self.sci_ali_img_size[0],cominx:cominx+self.sci_ali_img_size[1]] = self.sci_ali_img_hdu.data
                            new_sci_img = pad_sci
                            self.sp_logger.info(info_g+' Successfully padded science image')
                            # sys.exit()
                        except Exception as e:
                            self.sp_logger.warning(warn_y+f" Failed padding, trying to Swarp again {e}")
                            # sys.exit()
                            self.align_success=False
                            continue

                    # new_ref_img = pad_ref
                    self.sci_ali_img_hdu.data = pad_sci
                    # self.ref_ali_img_hdu.data = pad_ref
                    #give the science image the same wcs as the reference image
                    self.ref_wcs = WCS(self.ref_ali_img_hdu.header)
                    self.sci_ali_img_hdu.header = self.ref_ali_img_hdu.header
                    # except Exception as e:
                

                    for keyw in ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2']:
                        self.sci_ali_img_hdu.header[keyw] = self.ref_ali_img_hdu.header[keyw]


                    # self.sp_logger.info(self.sci_ali_img_hdu.data)
                    self.sci_ali_dupe = self.sci_ali_img_hdu.data
                    # Track which pixels are real data vs padding before converting NaN/zeros
                    self.valid_mask = np.zeros(pad_sci.shape, dtype=bool)
                    self.valid_mask[cominy:cominy+self.sci_ali_img_size[0], cominx:cominx+self.sci_ali_img_size[1]] = True
                    self.ref_valid_mask = np.ones(self.ref_ali_img_hdu.data.shape, dtype=bool)
                    self.sp_logger.info(info_g+f' Valid pixel mask created: {self.valid_mask.sum()} / {self.valid_mask.size} pixels are real data ({100*self.valid_mask.mean():.1f}%)')
                    self.sci_ali_img_hdu.data = np.nan_to_num(self.sci_ali_img_hdu.data)#np.nanmedian(self.sci_ali_img_hdu.data))
                    fits.writeto(self.path+'aligned_images/'+self.sci_img_name.split('/')[-1].replace('.fits','_padded.fits'),self.sci_ali_img_hdu.data,overwrite=True,header=self.sci_ali_img_hdu.header)
                    # self.sci_ali_name,_ = self.update_astrometry(self.path+'aligned_images/'+self.sci_img_name.split('/')[-1].replace('.fits','_padded.fits'),ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])
                    # self.sp_logger.info('new sci ali name',self.sci_ali_name)
                    # sys.exit()

                    self.sci_ali_name = self.path+'aligned_images/'+self.sci_img_name.split('/')[-1].replace('.fits','_padded.fits')
                    self.align_success=True
                    # print(np.shape(self.sci_ali_img_hdu.data))
                    # print(np.shape(self.ref_ali_img_hdu.data))
                    # self.sp_
                    sci_ali_x,sci_ali_y = np.shape(self.sci_ali_img_hdu.data)
                    ref_ali_x,ref_ali_y = np.shape(self.ref_ali_img_hdu.data)
                    
                    self.sp_logger.info(info_g+' Alignment with padding successful')
                    self.sp_logger.info(info_g+' Aligned Science image: '+self.sci_ali_name+' '+str(sci_ali_x)+'x'+str(sci_ali_y))
                    self.sp_logger.info(info_g+' Aligned Reference image: '+self.ref_ali_name+' '+str(ref_ali_x)+'x'+str(ref_ali_y))
                    # sys.exit()
                    # self.align_fail_count=+1
                    if all(s<500 for s in self.sci_ali_img_hdu.data.shape):
                        self.sp_logger.warning(warn_r+' Aligned science image size too small after padding, alignment failed')
                        try_aa_instead=True
                except Exception as e:
                    self.sp_logger.warning(warn_r+' Alignment failed '+e)
                    self.align_success=False
                    self.align_fail_count=+1
                    # sys.exit()
                    # self.align_fail_count=+1

            if try_aa_instead:
                self.sp_logger.info(info_g+' Trying alignment with reproject and astroalign instead')
                ref_mask = np.where(np.isnan(self.ref_img_hdu.data),False,True) #True where there is data, False where there is not

            
                self.ref_img_hdu.data = self.ref_img_hdu.data*ref_mask
                fig = plt.figure(figsize=(10,10))

                self.sp_logger.info(info_g+f" Reprojecting the reference onto the science image")
                self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
                self.ref_resampled_o, self.footprint= reproject_interp(self.ref_img_hdu, self.sci_img_hdu.header,order=1)
                self.ref_resampled = self.ref_resampled_o
                self.footprint= self.footprint.astype(bool)==False

                self.ref_resampled[np.isnan(self.ref_resampled)] = np.nanmedian(self.ref_resampled)

                try:
                    # self.registered, self.footprint= aa.register(self.ref_resampled*self.footprint, self.sci_img_hdu.data)
                    self.sp_logger.info(info_g+" Aligning images with astroalign")
                    self.registered, self.footprint= aa.register(self.ref_resampled, np.nan_to_num(self.sci_img_hdu.data))
                except Exception as e:
                    self.sp_logger.info(warn_y+" Extra alignment with astroalign failed, using reprojected image instead")
                    self.sp_logger.info(e)
                self.registered = self.ref_resampled


                self.ref_masked = np.ma.masked_array(self.registered, self.footprint, fill_value=np.nanmedian(self.ref_resampled)).filled()
                self.ref_masked[np.isnan(self.ref_masked)] = np.nanmedian(self.ref_resampled)


                self.ref_ali_hdu = fits.PrimaryHDU(self.ref_masked,header=self.sci_img_hdu.header)

                self.sp_logger.info(info_g+" Saving aligned images to "+self.path+'aligned_images/')
                self.sci_ali_name=self.path+'aligned_images/'+self.sci_img_name[:-5]+'.aa.fits'
                self.ref_ali_name=self.path+'aligned_images/'+self.ref_img_name[:-5].split('/')[-1]+'.aa.fits'

                self.ref_ali_hdu.writeto(self.ref_ali_name,overwrite=True)
                self.sci_img_hdu.writeto(self.sci_ali_name,overwrite=True)

                # Science is written unchanged; reference was filled with nanmedian in gap regions
                self.valid_mask = np.ones(self.sci_img_hdu.data.shape, dtype=bool)
                self.ref_valid_mask = np.ones(self.sci_img_hdu.data.shape, dtype=bool)
                self.sp_logger.info(info_g+' Valid pixel masks set for aa/reproject fallback path (science unpadded, reference filled with sky)')

                # print(self.sci_ali_name,self.ref_ali_name)
                # sys.exit()

            # else:
            # sys.exit()
            self.align_success=True
            self.sp_logger.info(info_g+' Alignment successful')

            # For normal SWarp path (no padding, no aa), derive valid masks from the actual outputs
            if self.valid_mask is None:
                _sci_data = fits.open(self.sci_ali_name)[0].data
                self.valid_mask = np.isfinite(_sci_data) & (_sci_data != 0)
                self.sp_logger.info(info_g+f' Normal SWarp path: valid mask covers {100*self.valid_mask.mean():.1f}% of science image')
            if self.ref_valid_mask is None:
                _ref_data = fits.open(self.ref_ali_name)[0].data
                self.ref_valid_mask = np.isfinite(_ref_data) & (_ref_data != 0)

            self.sp_logger.info(info_g+' Attempting distortion correction if needed')
            # try:
            #     x,_= correct_distrotion(self.sci_ali_name,self.ref_ali_name)
            #     self.sci_ali_img_hdu.data = x
            #     fits.writeto(self.sci_ali_name.replace('.fits','_distortion_corrected.fits'),self.sci_ali_img_hdu.data,overwrite=True,header=self.sci_ali_img_hdu.header)
            #     self.sci_ali_name = self.sci_ali_name.replace('.fits','_distortion_corrected.fits')
            #     self.sp_logger.info(info_g+f" Distortion correction applied, saved to {self.sci_ali_name}")

            # except Exception as e:
            #     self.sp_logger.warning(warn_y+f" Distortion correction failed, error encountered: {e}")
            #     pass
            # try:
            #     if self.telescope in SEDM:
            #         self.distortion_correction()
            # except Exception as e:
            #     pass
            # sys.exit()
            # self.align_fail_count=+1
                
            if self.align_success and self.args.show_plots==True:
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                ali_fig = plt.figure(figsize=(10,10))
                plt.subplot(2,2,1)
                plt.imshow(self.sci_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Original Science image')
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.ref_img_hdu.data)
                plt.subplot(2,2,2)
                plt.imshow(self.ref_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Original Reference image')
                plt.subplot(2,2,3)
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_ali_img_hdu.data)
                plt.imshow(self.sci_ali_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Aligned Science image')
                plt.subplot(2,2,4)
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.ref_ali_img_hdu.data)
                plt.imshow(self.ref_ali_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Aligned Reference image')
                plt.show()
                # sys.exit()
                plt.close(ali_fig)

        if self.cutout_tf!=False:
            self.fig,self.ax = plt.subplots(nrows=1, ncols=3)
            self.sci_img_ali_hdu = fits.open(self.sci_ali_name)[0]
            self.percentile = np.percentile(self.sci_img_ali_hdu.data,[5,10,20,30,40,50,60,70,80,90]) 
            self.cmap = 'gray'
            self.coords_sn_sci_ali=wcs_to_pixels(self.sci_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
            self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_y = self.coords_sn_sci_ali
            # self.vmin,self.vmax,self.cmap = -((1.5*self.percentile[1])-self.percentile[3]),(3*self.percentile[5]) - self.percentile[7],'gray'
            self.vmin_sci,self.vmax_sci = visualization.ZScaleInterval().get_limits(self.sci_img_ali_hdu.data)
            self.coords_sn_sub=wcs_to_pixels(self.sci_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))
            self.coords_sn_sub_x,self.coords_sn_sub_y = self.coords_sn_sub[0]
            self.ax[0].imshow(self.sci_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_sci, vmax=self.vmax_sci)
            # self.sp_logger.info(self.coords_sn_sub_x,self.coords_sn_sub_y)
            self.ax[0].plot([self.coords_sn_sci_ali_x-30,self.coords_sn_sci_ali_x-10],[self.coords_sn_sci_ali_y,self.coords_sn_sci_ali_y],color='lime',lw=2.5),self.ax[0].plot([self.coords_sn_sci_ali_x+10,self.coords_sn_sci_ali_x+30],[self.coords_sn_sci_ali_y,self.coords_sn_sci_ali_y],color='lime',lw=2.5)
            self.ax[0].plot([self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_x],[self.coords_sn_sci_ali_y-30,self.coords_sn_sci_ali_y-10],color='lime',lw=2.5),self.ax[0].plot([self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_x],[self.coords_sn_sci_ali_y+10,self.coords_sn_sci_ali_y+30],color='lime',lw=2.5)
            self.ax[0].set_xlim(self.coords_sn_sci_ali_x-50,self.coords_sn_sci_ali_x+50)
            self.ax[0].set_ylim(self.coords_sn_sci_ali_y-50,self.coords_sn_sci_ali_y+50)

            self.ax[0].axis('off') 
            self.ax[1].axis('off') 
            self.ax[2].axis('off') 
            self.ax[0].set_title('New')
            if self.out_dir=="photometry/":       
                self.cutout_name = self.path+f'photometry_date/{self.folder}/cut_outs/'+self.sci_img_name[:-11]+"_cutout_panel"+self.img_type+".png"
            else:
                self.cutout_name = self.path+f'{self.out_dir}cut_outs/'+self.sci_img_name[:-11]+"_cutout_panel"+self.img_type+".png"

            self.fig.savefig(self.cutout_name)

        if self.cutout_tf!=False:
            self.ref_img_ali_hdu = fits.open(self.ref_ali_name)[0]
            self.coords_sn_ref_ali=wcs_to_pixels(self.ref_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
            self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_y = self.coords_sn_ref_ali

            self.vmin_ref,self.vmax_ref = visualization.ZScaleInterval().get_limits(self.ref_img_ali_hdu.data)
            self.ax[1].imshow(self.ref_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_ref, vmax=self.vmax_ref)
            self.ax[1].plot([self.coords_sn_ref_ali_x-30,self.coords_sn_ref_ali_x-10],[self.coords_sn_ref_ali_y,self.coords_sn_ref_ali_y],color='lime',lw=4),self.ax[1].plot([self.coords_sn_ref_ali_x+10,self.coords_sn_ref_ali_x+30],[self.coords_sn_ref_ali_y,self.coords_sn_ref_ali_y],color='lime',lw=2.5)
            self.ax[1].plot([self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_x],[self.coords_sn_ref_ali_y-30,self.coords_sn_ref_ali_y-10],color='lime',lw=2.5),self.ax[1].plot([self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_x],[self.coords_sn_ref_ali_y+10,self.coords_sn_ref_ali_y+30],color='lime',lw=2.5)
            self.ax[1].set_xlim(self.coords_sn_ref_ali_x-50,self.coords_sn_ref_ali_x+50)
            self.ax[1].set_ylim(self.coords_sn_ref_ali_y-50,self.coords_sn_ref_ali_y+50)
            self.ax[1].set_title('Ref')
            self.ax[1].axis('off')        

            self.fig.savefig(self.cutout_name)




        return self.ref_ali_name,self.ref_img_hdu,self.ref_img, self.sci_ali_name
    



    
    # def py_ref_aa_align(self,image_size=image_size):
    def py_ref_align(self,image_size=image_size):
    

        if image_size==1500:
            self.image_size=image_size
        else:
            self.image_size=image_size

        self.sp_logger.info(info_g+f' Aligning science with reference image')

        if self.auto_ref=='auto':
            if self.sci_filt=='u' or self.use_sdss==True: #if filt ==u, need to call make sdss_ref before this or can call it here
                self.make_sdss_ref()
                if self.sys_exit==True:
                    return
                self.image_name=self.sci_obj+'_ref.fits'
                self.ref_path=self.path+'ref_imgs/'+self.image_name


            if self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z':
                self.ref_folder = self.path+'ref_imgs/stack_'+str(self.sci_filt)
                self.ref_path = self.path+'ref_imgs/stack_'+str(self.sci_filt)+'_ra'+self.ra_string+'_dec'+self.dec_string+'_arcsec*'
                if len(glob.glob(self.ref_path))>0:
                    self.sp_logger.info(info_b+' PS1 reference image already exists')
                if len(glob.glob(self.ref_path))==0:
                    self.sp_logger.info(info_g+' Downloading reference image from PS...')
                    try:
                        # self.sp_logger.info(panstamps_path+' -f --width='+str(self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))
                        os.system(panstamps_path+' -f --width='+str(self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))
                        
                        self.sp_logger.info(info_g+' Reference image downloaded')
                    except Exception as e:
                        self.sp_logger.warning(warn_r+f" {self.sci_obj}, {self.sci_mjd}, {self.sci_filt}: Unable to download reference image from panstarrs, error encountered :{e}")
                        self.sp_logger.info(warn_r+f" Error downloading PS reference image, error encountered :{e}")
                        self.sp_logger.info(warn_r+f" Exiting...")
                        ##sys.exit(1)
                        self.sys_exit=True
                        return
                
            if self.sys_exit==True:
                return

            self.ref_path = glob.glob(self.ref_path)[0]
            # self.sp_logger.info(glob.glob(self.ref_path))
            # self.sp_logger.info(self.ref_path)
            self.ref_img_name=glob.glob(self.ref_path)[0]
            self.ref_img_hdu=fits.open(self.ref_img_name)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
        
        else:
            #If a reference image is passed in with it's full path as self.auto_ref
            if self.auto_ref.startswith(self.path):
                self.auto_ref = self.path+self.auto_ref

            # self.sp_logger.info(self.auto_ref)
            # self.sp_logger.info(self.sci_path)
            if not os.path.exists(self.auto_ref):
                self.sp_logger.info(warn_r+f' Reference image {self.auto_ref} not found')
                self.sys_exit=True
                return

            self.ref_img_name=glob.glob(self.auto_ref)[0]
            self.ref_img_hdu=fits.open(self.auto_ref)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
            # self.sp_logger.info(self.ref_img==self.sci_bkgsb)

        self.coords_sn_ref=wcs_to_pixels(self.ref_img_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
        self.coords_sn_ref_x,self.coords_sn_ref_y = self.coords_sn_ref




        if not os.path.exists(self.path+'aligned_images'):
            os.makedirs(self.path+'aligned_images')

        wcs_command='python3 '+self.path+'sedm_align.py'+" -sci "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" -ref "+self.ref_img_name+"  -m relative -r 100"


        self.align_success=False
        self.align_fail_count=0

        
        os.system(wcs_command)
        # sys.exit(1)


        while self.align_success==False:
            self.align_fail_count+=1
            self.sp_logger.info(info_g+" Aligning images with reprojection and astroalgin making a resampled image of size "+str(self.image_size)+"x"+str(self.image_size)+" pixels")

            

            self.sci_img_hdu = fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
            self.sci_img_header = self.sci_img_hdu.header

            #detect the edges of the reference image
            self.ref_img_edge_x = [i for i in range(self.ref_img_hdu.header['NAXIS1']) if np.isnan(self.ref_img_hdu.data[:,i]).all()]
            self.ref_img_edge_y = [i for i in range(self.ref_img_hdu.header['NAXIS2']) if np.isnan(self.ref_img_hdu.data[i,:]).all()]


            if len(self.ref_img_edge_y)==0:
                self.ref_img_edge_y = [0, self.ref_img_hdu.header['NAXIS2']-1]
            elif len(self.ref_img_edge_y)>0:
                if self.ref_img_edge_y[0] == 0: self.ref_img_edge_y = [self.ref_img_edge_y[-1],self.ref_img_hdu.header['NAXIS2']-1]
                else: self.ref_img_edge_y = [0,self.ref_img_edge_y[0]]
            if len(self.ref_img_edge_x)==0:
                self.ref_img_edge_x = [0, self.ref_img_hdu.header['NAXIS1']-1]
            elif len(self.ref_img_edge_x)>0:
                if self.ref_img_edge_x[0] == 0: self.ref_img_edge_x = [self.ref_img_edge_x[-1],self.ref_img_hdu.header['NAXIS1']-1]
                else: self.ref_img_edge_x = [0,self.ref_img_edge_x[0]]

            self.ref_img_hdu_0 = self.ref_img_hdu.copy()
            ref_mask = np.where(np.isnan(self.ref_img_hdu.data),False,True) #True where there is data, False where there is not

            
            # fig = plt.figure(figsize=(10,10))
            # plt.title('Reference image')
            # vmin,vmax = visualization.ZScaleInterval().get_limits(self.ref_img_hdu.data)
            # plt.imshow(self.ref_img_hdu.data,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)

    
            self.ref_img_hdu.data = self.ref_img_hdu.data*ref_mask
            fig = plt.figure(figsize=(10,10))

            # plt.title('Reference image masked')
            # plt.imshow(self.ref_img_hdu.data,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)
            # plt.show()

            self.sp_logger.info(info_g+f" Reprojecting the reference onto the science image")
            self.ref_resampled_o, self.footprint= reproject_interp(self.ref_img_hdu, self.sci_img_header,order=1)
            self.ref_resampled = self.ref_resampled_o
            self.footprint= self.footself.sp_logger.info.astype(bool)==False

            self.ref_resampled[np.isnan(self.ref_resampled)] = np.nanmedian(self.ref_resampled)

            try:
                # self.registered, self.footprint= aa.register(self.ref_resampled*self.footself.sp_logger.info, self.sci_img.data)
                self.sp_logger.info(info_g+" Aligning images with astroalign")
                self.registered, self.footprint= aa.register(self.ref_resampled, self.sci_img.data)
            except Exception as e:
                self.sp_logger.info(warn_y+" Extra alignment with astroalign failed, using reprojected image instead")
                self.sp_logger.info(e)
            self.registered = self.ref_resampled


            self.ref_masked = np.ma.masked_array(self.registered, self.footself.sp_logger.info, fill_value=np.nanmedian(self.ref_resampled)).filled()
            self.ref_masked[np.isnan(self.ref_masked)] = np.nanmedian(self.ref_resampled)


            self.ref_ali_hdu = fits.PrimaryHDU(self.ref_masked,header=self.sci_img_header)

            self.sp_logger.info(info_g+" Saving aligned images to "+self.path+'aligned_images/')
            self.sci_ali_name=self.path+'aligned_images/'+self.sci_img_name[:-5]+'.resamp.fits'
            self.ref_ali_name=self.path+'aligned_images/'+self.ref_img_name[:-5].split('/')[-1]+'.resamp.fits'

            self.ref_ali_hdu.writeto(self.ref_ali_name,overwrite=True)
            self.sci_img_hdu.writeto(self.sci_ali_name,overwrite=True)

        
            if self.align_fail_count==1:
                self.sp_logger.info(info_g+' Alignment attempt 1')
            else: 
                self.sp_logger.info(warn_y+' Alignment attempt '+str(self.align_fail_count))




            self.ref_ali_img_hdu = fits.open(self.ref_ali_name)[0]
            self.sci_ali_img_hdu = fits.open(self.sci_ali_name)[0]

            
            


            self.sp_logger.info(info_g+' Science & Reference image sizes : ',np.shape(fits.open(self.sci_ali_name)[0].data),np.shape(fits.open(self.ref_ali_name)[0].data))
            self.sp_logger.info(info_b+' Target size : ',self.image_size,'x',self.image_size)
            self.sp_logger.info(info_g+' Aligned Science image size : ',np.shape(fits.open(self.sci_ali_name)[0].data))
            self.sp_logger.info(info_g+' Aligned Reference image size : ',np.shape(fits.open(self.ref_ali_name)[0].data))

            if self.args.show_plots==True:
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                ali_fig = plt.figure(figsize=(10,10))
                plt.subplot(2,2,1)
                plt.imshow(self.sci_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Original Science image')
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.ref_img_hdu.data)
                plt.subplot(2,2,2)
                plt.imshow(self.ref_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Original Reference image')
                plt.subplot(2,2,3)
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                plt.imshow(self.sci_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Aligned Science image')
                plt.subplot(2,2,4)
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.ref_ali_img_hdu.data)
                plt.imshow(self.ref_ali_img_hdu.data,vmin=self.vmin,vmax=self.vmax,cmap='gray')
                plt.title('Aligned Reference image')
                plt.show()
                # sys.exit()
                plt.close(ali_fig)

            if self.align_fail_count>15:
                self.sp_logger.info(warn_y+' Alignment failed 15 times')
                break

            self.align_success=True
            self.sp_logger.info(info_g+' Alignment successful')

            # if self.telescope in SEDM:
            #     self.distortion_correction()

                
        
            if self.cutout_tf!=False:
                self.fig,self.ax = plt.subplots(nrows=1, ncols=3)
                self.sci_img_ali_hdu = fits.open(self.sci_ali_name)[0]
                self.percentile = np.percentile(self.sci_img_ali_hdu.data,[5,10,20,30,40,50,60,70,80,90]) 
                self.cmap = 'gray'
                self.coords_sn_sci_ali=wcs_to_pixels(self.sci_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
                self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_y = self.coords_sn_sci_ali
                self.vmin_sci,self.vmax_sci = visualization.ZScaleInterval().get_limits(self.sci_img_ali_hdu.data)

                self.ax[0].imshow(self.sci_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_sci, vmax=self.vmax_sci)
                self.ax[0].plot([self.coords_sn_sci_ali_x-30,self.coords_sn_sci_ali_x-10],[self.coords_sn_sci_ali_y,self.coords_sn_sci_ali_y],color='lime',lw=2.5),self.ax[0].plot([self.coords_sn_sci_ali_x+10,self.coords_sn_sci_ali_x+30],[self.coords_sn_sci_ali_y,self.coords_sn_sci_ali_y],color='lime',lw=2.5)
                self.ax[0].plot([self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_x],[self.coords_sn_sci_ali_y-30,self.coords_sn_sci_ali_y-10],color='lime',lw=2.5),self.ax[0].plot([self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_x],[self.coords_sn_sci_ali_y+10,self.coords_sn_sci_ali_y+30],color='lime',lw=2.5)
                self.ax[0].set_xlim(self.coords_sn_sci_ali_x-115,self.coords_sn_sci_ali_x+115)
                self.ax[0].set_ylim(self.coords_sn_sci_ali_y-115,self.coords_sn_sci_ali_y+115)
                self.ax[0].axis('off') 
                self.ax[1].axis('off') 
                self.ax[2].axis('off') 
                self.ax[0].set_title('New')
                if self.out_dir=="photometry/":       
                    self.cutout_name = self.path+f'photometry_date/{self.folder}/cut_outs/'+self.sci_img_name[:-11]+"_cutout_panel"+self.img_type+".png"
                else:
                    self.cutout_name = self.path+f'{self.out_dir}cut_outs/'+self.sci_img_name[:-11]+"_cutout_panel"+self.img_type+".png"

                self.fig.savefig(self.cutout_name)

            if self.cutout_tf!=False:
                self.ref_img_ali_hdu = fits.open(self.ref_ali_name)[0]
                self.coords_sn_ref_ali=wcs_to_pixels(self.ref_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
                self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_y = self.coords_sn_ref_ali

                self.vmin_ref,self.vmax_ref = visualization.ZScaleInterval().get_limits(self.ref_img_ali_hdu.data)
                self.ax[1].imshow(self.ref_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_ref, vmax=self.vmax_ref)
                self.ax[1].plot([self.coords_sn_sci_ref_x-30,self.coords_sn_sci_ref_x-10],[self.coords_sn_sci_ref_y,self.coords_sn_sci_ref_y],color='lime',lw=2.5),self.ax[1].plot([self.coords_sn_sci_ref_x+10,self.coords_sn_sci_ref_x+30],[self.coords_sn_sci_ref_y,self.coords_sn_sci_ref_y],color='lime',lw=2.5)
                self.ax[1].plot([self.coords_sn_sci_ref_x,self.coords_sn_sci_ref_x],[self.coords_sn_sci_ref_y-30,self.coords_sn_sci_ref_y-10],color='lime',lw=2.5),self.ax[1].plot([self.coords_sn_sci_ref_x,self.coords_sn_sci_ref_x],[self.coords_sn_sci_ref_y+10,self.coords_sn_sci_ref_y+30],color='lime',lw=2.5)
                self.ax[1].set_xlim(self.coords_sn_ref_ali_x-50,self.coords_sn_ref_ali_x+50)
                self.ax[1].set_ylim(self.coords_sn_ref_ali_y-50,self.coords_sn_ref_ali_y+50)
                self.ax[1].set_title('Ref')
                self.ax[1].axis('off')        

                self.fig.savefig(self.cutout_name)


            return self.ref_ali_name,self.ref_img_hdu,self.ref_img, self.sci_ali_name

            # except Exception as e:
            #     self.sp_logger.info(warn_y+f" Alignment failed with error: {e}")
            #     # self.align_fail_count+=1
            #     self.sp_logger.info(warn_y+f" Trying again with different parameters. Attempt {self.align_fail_count} of {20}")
            #     if self.align_fail_count==20:
            #         self.sp_logger.info(warn_y+f" Alignment failed 20 times. Aborting")
            #         self.align_success=False
            #         return None,None,None,None


    def psfex_convolve_images(self):
        self.sp_logger.info(info_g+' Beginning to convolve images with SeXtractor and PSFEx')

        if not os.path.exists(self.path+'convolved_sci'):
            os.makedirs(self.path+'convolved_sci')

        if not os.path.exists(self.path+'convolved_ref'):
            os.makedirs(self.path+'convolved_ref')

        if not os.path.exists(self.path+'out'):
            os.makedirs(self.path+'out')

        try:
            self.sci_conv_name=self.path+"convolved_sci/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'sci_convolved.fits'
        except:
            self.sci_img_hdu=self.sci_img_hdu[0]
            self.sci_conv_name=self.path+"convolved_sci/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'sci_convolved.fits'
        self.ref_conv_name=self.path+"convolved_ref/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'ref_convolved.fits'

        
        self.files_to_clean.append(self.sci_conv_name)
        self.files_to_clean.append(self.ref_conv_name)
         #################################################
        # CONVOLVE REFERENCE WITH PSF OF SCIENCE IMAGE

        if os.path.exists(self.path+f"config_files/sci_prepsfex_{self.rand_nums_string}.cat"):
            os.system("rm "+self.path+f"config_files/sci_prepsfex_{self.rand_nums_string}.cat")


        self.sp_logger.info(info_g+' Convolving the reference with the PSF of the science image')

        prepsexfile()
        psfexfile()

        # SExtractor command for the science image
        # print(self.sci_ali_name)
        # sys.exit(1)
        MAGZP = 25.0
        sextractor_command=sex_path+" "+self.sci_ali_name+" -c "+self.path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT"+" "+str(MAGZP)
        self.sp_logger.info(info_g+f' Running SExtractor: {sextractor_command}')
        self.sp_logger.info(info_g+' Creating PSFex catalog with SExtractor')

        sex_status = os.system(sextractor_command)
        self.sp_logger.info(info_g+' SExtractor status: '+str(sex_status))

        

        if os.path.exists(self.path+f'out/sci_proto_prepsfex_{self.rand_nums_string}.fits'):
            os.system('rm '+self.path+f'out/sci_proto_prepsfex_{self.rand_nums_string}.fits')
        self.files_to_clean.append(self.path+f'out/sci_proto_prepsfex_{self.rand_nums_string}.fits')
        self.files_to_clean.append(self.path+f'out/sci_prepsfex_{self.rand_nums_string}.psf')
        self.files_to_clean.append(self.path+f'out/sci_resi_{self.rand_nums_string}.fits')
        self.files_to_clean.append(self.path+f'out/sci_subsym_{self.rand_nums_string}.fits')
        self.files_to_clean.append(self.path+f'out/sci_moffat_{self.rand_nums_string}.fits')

        self.sp_logger.info(info_g+' Running PSFex with SExtractor catalog: '+self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}.cat")
        #check the catalog output by sextractor and change the stars matched to have a flag of 0


        self.sp_logger.info(info_g+' Running PSFex with SExtractor catalog:'+" "+self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")

        psfex_status = os.system(psfex_path+" "+self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")
        self.sp_logger.info(info_g+' PSFex status: '+str(psfex_status))

        self.files_to_clean.append(self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}.cat")
        # self.sp_logger.info(psfex.PSFEx(path+f'config_files/prepsfex_{self.rand_nums_string}.cat'))
        # sys.exit(1)


        self.psf_sci_image_name=self.path+f'out/proto_sci_prepsfex_{self.rand_nums_string}.fits'
        self.sp_logger.info(info_g+ ' PSFEx science image created: '+self.psf_sci_image_name)
        self.files_to_clean.append(self.path+f'out/proto_sci_prepsfex_{self.rand_nums_string}.fits')
        self.psf_sci_image = fits.open(self.psf_sci_image_name)

        # self.sp_logger.info(self.psf_sci_image[0].data)
        # fig = plt.figure(figsize=(12,8))
        # plt.imshow(self.psf_sci_image[0].data[0])
        # fig.savefig('ZTF21aceqrju_g_psf.png')
        # plt.show()

        self.hdu_psf_model_sci= fits.open(self.path+f'out/sci_prepsfex_{self.rand_nums_string}.psf')
        self.files_to_clean.append(self.path+f'out/sci_prepsfex_{self.rand_nums_string}.psf')
        self.chi_sq_psf=self.hdu_psf_model_sci[1].header['CHI2']


        self.sp_logger.info(info_g+' Reduced Chi^2 of science image PSF fit: '+"%.1f" % self.chi_sq_psf)
        # ── ePSF fallback when PSFEx chi² is poor ────────────────────────────
        # chi² > 3 indicates PSFEx struggled (no bright isolated stars, bad
        # seeing, crowded field, etc.).  In that case we replace the PSFEx
        # kernel with an EPSFBuilder kernel from build_psf.py before it gets
        # assigned to self.kernel_sci below.
        self._epsf_sci_kernel = None
        if self.chi_sq_psf>3:
            # print(colored('Warning: PSF model may not be accurate','green'))
            self.sp_logger.warning(warn_y+' Warning: PSF model may not be accurate')
            self.sp_logger.info(info_g+' Attempting ePSF fallback (build_psf) for science image …')
            try:
                from build_psf import build_psf_from_fits as _build_epsf
                _fwhm_guess = getattr(self, 'sci_seeing', 2.0) / getattr(self, 'sci_ps', 0.37)
                _epsf_res = _build_epsf(
                    self.sci_ali_name,
                    fwhm_guess  = max(1.5, _fwhm_guess),
                    threshold_sigma = 4.0,
                    min_snr     = 10,
                    max_stars   = 40,
                    min_stars   = 3,
                    epsf_iters  = 3,
                )
                self._epsf_sci_kernel = _epsf_res['kernel']
                self.sp_logger.info(
                    info_g +
                    f' ePSF science fallback OK: FWHM={_epsf_res["fwhm"]:.2f} px'
                    f'  elong={_epsf_res["elongation"]:.2f}'
                    f'  n_stars={_epsf_res["n_stars"]}'
                    f'  scatter={_epsf_res["fwhm_scatter"]:.3f} px'
                )
            except Exception as _epsf_e:
                self.sp_logger.warning(warn_y + f' ePSF science fallback failed: {_epsf_e}')
        if self.chi_sq_psf<1e-10:
            self.sp_logger.warning(warn_r+' Warning: PSF model is a perfect fit, may be overfitting')
            self.sp_logger.warning(warn_y+' Trying again with only the highest signal-to-noise stars')
            # try:
            self.sci_prepsfex_cat = fits.open(self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}.cat")
            self.sci_prepsfex_cat_tab = Table(self.sci_prepsfex_cat[2].data)
            self.sci_prepsfex_cat_snr_min = 37#np.percentile(self.sci_prepsfex_cat_tab['SNR_WIN'],10)
            # print(self.sci_prepsfex_cat_snr_min)
            # print(np.min(self.sci_prepsfex_cat_tab['SNR_WIN']))
            self.sci_prepsfex_cat[2].data = self.sci_prepsfex_cat[2].data[(self.sci_prepsfex_cat_tab['SNR_WIN']>self.sci_prepsfex_cat_snr_min)]#&(self.sci_prepsfex_cat_tab['ELONGATION']<1.2)]
            # print(Table(self.sci_prepsfex_cat[2].data))
            self.sci_prepsfex_cat.writeto(self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}_highSNR.cat",overwrite=True)
            self.sp_logger.info(info_g+' Writing high SNR stars to '+self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}_highSNR.cat")
            self.sp_logger.info(info_g+' Running PSFex with only the highest signal-to-noise stars')
            psfex_status = os.system(psfex_path+" "+self.path+f"temp_config_files/sci_prepsfex_{self.rand_nums_string}_highSNR.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET -SAMPLE_MINSN 10")
            self.sp_logger.info(info_g+' High SNR PSFex status: '+str(psfex_status))
            
            self.hdu_psf_model_sci= fits.open(self.path+f'out/sci_prepsfex_{self.rand_nums_string}_highSNR.psf')
            self.files_to_clean.append(self.path+f'out/sci_prepsfex_{self.rand_nums_string}_highSNR.psf')
            self.chi_sq_psf=self.hdu_psf_model_sci[1].header['CHI2']
            self.sp_logger.info(info_g+' Reduced Chi^2 of science image PSF fit with high SNR stars: '+"%.1f" % self.chi_sq_psf)
            # sys.exit(1)
            self.sys_exit=True
            return
        # sys.exit(1)
        # Use ePSF fallback kernel if PSFEx chi² was bad (set above), otherwise
        # use the standard PSFEx proto kernel.
        if getattr(self, '_epsf_sci_kernel', None) is not None:
            self.kernel_sci = self._epsf_sci_kernel
            self.sp_logger.info(info_g+' Using ePSF kernel for science convolution (PSFEx chi² too high)')
        else:
            self.kernel_sci = self.psf_sci_image[0].data[0]

        if self.to_subtract!=False:
        # Read the REFERENCE image and convolve it with the  kernel science
            self.ref_image_aligned=fits.open(self.ref_ali_name)
            self.files_to_clean.append(self.ref_ali_name)
            self.ref_conv = scipy_convolve(self.ref_image_aligned[0].data, self.kernel_sci, mode='same', method='fft')
            self.ref_conv=np.nan_to_num(self.ref_conv)
            self.sp_logger.info(info_g+' Saving convolved reference image to '+self.ref_conv_name)
            fits.writeto(self.ref_conv_name, data=self.ref_conv, header=self.ref_image_aligned[0].header,overwrite=True)

            self.sp_logger.info(info_g+' Convolving the science with the PSF of the reference image')

            sextractor_command=sex_path+" "+self.ref_ali_name+" -c "+self.path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"temp_config_files/ref_prepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT 25.0"
            os.system(sextractor_command)

            if self.sci_obj in ['2023ixf','SN2023ixf']:
                #open the sextractor catalog and keep only the point sources that are in calstars4e.cat
                self.good_ps1_ixf_stars = ascii.read(self.path+'calstars4e.cat',names=['ra','dec','u','g','r','i','z','PS1g','PS1r','PS1i','PS1z'],format='no_header',delimiter=' ')
                self.sp_logger.info(info_g+' Using calstars4e.cat to select point sources for PSFEx for SN2023ixf')
                self.good_ps1_ixf_stars['ra'].unit = u.deg
                self.good_ps1_ixf_stars['dec'].unit = u.deg

                self.sex_ref_orig = fits.open(self.path+f"temp_config_files/ref_prepsfex_{self.rand_nums_string}.cat" ) #read in the original sextractor catalog
                self.sex_ref = Table(self.sex_ref_orig[2].data)
                self.ref_ali_wcs = WCS(self.ref_ali_name)   
                self.sex_ref_wcs = self.ref_ali_wcs.all_pix2world(np.column_stack((self.sex_ref['X_IMAGE'],self.sex_ref['Y_IMAGE'])),1) #find the RA and DEC coordinates of the sextractor catalog in wcs
                self.sex_ref_wcs = SkyCoord(ra=self.sex_ref_wcs[:,0]*u.deg, dec=self.sex_ref_wcs[:,1]*u.deg,frame='fk5')
                self.sex_ref['RA_DEG'],self.sex_ref['DEC_DEG'] = self.sex_ref_wcs.ra,self.sex_ref_wcs.dec
                self.sex_ref_wcs = SkyCoord(ra=np.array(self.sex_ref['RA_DEG'])*u.deg, dec=np.array(self.sex_ref['DEC_DEG'])*u.deg,frame='fk5')
                self.good_ps1_ixf_stars_wcs = SkyCoord(ra=np.array(self.good_ps1_ixf_stars['ra'])*u.deg, dec=np.array(self.good_ps1_ixf_stars['dec'])*u.deg,frame='fk5')

                
                self.indx, self.d2d, self.d3d = self.sex_ref_wcs.match_to_catalog_sky(self.good_ps1_ixf_stars_wcs)
                self.upd_indx=np.where(self.d2d<(search_rad)/3600.*u.deg)[0]

                self.sp_logger.info(info_g+' Bright stars found='+str(len(self.upd_indx))+', '+str(round(3600*(self.ref_width/60),6))+'arcsec search radius')
                self.bright_stars_sc = Table()
                self.bright_stars_sc['ra'] = self.good_ps1_ixf_stars_wcs.ra.deg
                self.bright_stars_sc['dec'] = self.good_ps1_ixf_stars_wcs.dec.deg
                self.bright_stars_sc_pix = self.ref_ali_wcs.all_world2pix(np.column_stack((self.bright_stars_sc['ra'],self.bright_stars_sc['dec'])),1)
                self.bright_stars_sc_pix_x,self.bright_stars_sc_pix_y = self.bright_stars_sc_pix[:,0],self.bright_stars_sc_pix[:,1]
                self.bright_stars_sc['xcentroid'],self.bright_stars_sc['ycentroid'] = self.bright_stars_sc_pix_x,self.bright_stars_sc_pix_y

                for col in ['u','g','r','i','z','PS1g','PS1r','PS1i','PS1z']:
                    self.bright_stars_sc[col] = self.good_ps1_ixf_stars[col]



                self.crossmatch_sex_ref = self.sex_ref[self.upd_indx]

                self.sex_ref_orig[2].data = self.sex_ref_orig[2].data[self.upd_indx]
                # print(self.crossmatch_sex_ref['RA_DEG','DEC_DEG','X_IMAGE','Y_IMAGE'])
                # sys.exit(1)
                self.sex_ref_orig[2].data['FLAGS'] = 0
                # self.sp_logger.info(self.sex_ref_orig[2].data['X_IMAGE'])

                # for col in self.sex_ref_orig[2].columns.names:
                #     for ind in self.upd_indx:
                #         self.sex_ref_orig[2].data[col] = self.sex_ref_orig[2].data[col][self.upd_indx] 

                #delete rows if not in upd_indx
                # for k in range(len(self.sex_ref_orig[2].data)):
                #     if k not in self.upd_indx:
                        

                self.sex_ref_orig[2].header['NAXIS2'] = len(self.upd_indx)





                self.sex_ref_orig.writeto(self.path+f"temp_config_files/ref_prepsfex_{self.rand_nums_string}.cat",overwrite=True)


            if os.path.exists(self.path+f'out/proto_ref_prepsfex_{self.rand_nums_string}.fits'):
                os.system('rm '+self.path+f'out/proto_ref_prepsfex_{self.rand_nums_string}.fits')

            

            # os.system(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex")


            # sextractor_command=sex_path+" "+self.ref_ali_name+" -c "+self.path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"config_files/pooprepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT 25.0"
            # self.sp_logger.info(sextractor_command)

            # self.sp_logger.info(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")
            # sys.exit(1)
            # self.sp_logger.info(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE FULL")
            # print(psfex_path+" "+self.path+f"config_files/ref_prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")
            os.system(psfex_path+" "+self.path+f"temp_config_files/ref_prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")
            # sys.exit(1)
            self.files_to_clean.append(self.path+f"temp_config_files/ref_prepsfex_{self.rand_nums_string}.cat")

            self.psf_ref_image_name=self.path+f'out/proto_ref_prepsfex_{self.rand_nums_string}.fits'
            self.sp_logger.info(info_g+' PSFEx reference image created '+self.psf_ref_image_name)
            self.files_to_clean.append(self.path+f'out/proto_ref_prepsfex_{self.rand_nums_string}.fits')
            self.psf_ref_image = fits.open(self.psf_ref_image_name)

            self.hdu_psf_model_ref= fits.open(self.path+f'out/ref_prepsfex_{self.rand_nums_string}.psf')
            self.files_to_clean.append(self.path+f'out/ref_prepsfex_{self.rand_nums_string}.psf')
            self.chi_sq_psf_ref=self.hdu_psf_model_ref[1].header['CHI2']


            # self.sp_logger.info(colored('Reduced Chi^2 of science image PSF fit:','green'),"%.1f" % self.chi_sq_psf_ref)
            self.sp_logger.info(info_g+' Reduced Chi^2 of reference image PSF fit: '+str(self.chi_sq_psf_ref))
            # ── ePSF fallback for reference image ─────────────────────────────
            self._epsf_ref_kernel = None
            if self.chi_sq_psf_ref>3:
                # self.sp_logger.warning(colored('Warning: PSF ref model may not be accurate','green'))
                self.sp_logger.warning(warn_y+' Warning: PSF ref model may not be accurate')
                self.sp_logger.info(info_g+' Attempting ePSF fallback (build_psf) for reference image …')
                try:
                    from build_psf import build_psf_from_fits as _build_epsf
                    _fwhm_guess = getattr(self, 'sci_seeing', 2.0) / getattr(self, 'sci_ps', 0.37)
                    _epsf_ref_res = _build_epsf(
                        self.ref_ali_name,
                        fwhm_guess  = max(1.5, _fwhm_guess),
                        threshold_sigma = 4.0,
                        min_snr     = 10,
                        max_stars   = 40,
                        min_stars   = 3,
                        epsf_iters  = 3,
                    )
                    self._epsf_ref_kernel = _epsf_ref_res['kernel']
                    self.sp_logger.info(
                        info_g +
                        f' ePSF reference fallback OK: FWHM={_epsf_ref_res["fwhm"]:.2f} px'
                        f'  elong={_epsf_ref_res["elongation"]:.2f}'
                        f'  n_stars={_epsf_ref_res["n_stars"]}'
                        f'  scatter={_epsf_ref_res["fwhm_scatter"]:.3f} px'
                    )
                except Exception as _epsf_ref_e:
                    self.sp_logger.warning(warn_y + f' ePSF reference fallback failed: {_epsf_ref_e}')

            # Use ePSF fallback kernel if PSFEx chi² was bad, otherwise use PSFEx proto kernel.
            if getattr(self, '_epsf_ref_kernel', None) is not None:
                self.kernel_ref = self._epsf_ref_kernel
                self.sp_logger.info(info_g+' Using ePSF kernel for reference convolution (PSFEx chi² too high)')
            else:
                self.kernel_ref = self.psf_ref_image[0].data[0]

            # Read the Science image and convolve it with the Gaussian kernel
            # print(self.sci_ali_name)
            self.sci_image_aligned=fits.open(self.sci_ali_name)
            self.files_to_clean.append(self.sci_ali_name)
            self.sci_conv = scipy_convolve(self.sci_image_aligned[0].data,self.kernel_ref, mode='same', method='fft')
            self.sci_conv=np.nan_to_num(self.sci_conv)

            self.sp_logger.info(info_g+' Saving convolved science image to '+self.sci_conv_name)
            fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)
            # sys.exit(1)

            # self.sp_logger.info(self.sci_conv_name)
            # self.sp_logger.info(self.ref_conv_name)
            self.sci_conv_hdu,self.ref_conv_hdu,self.ref_ali_hdu,self.sci_ali_hdu=fits.open(self.sci_conv_name)[0],fits.open(self.ref_conv_name)[0],fits.open(self.ref_ali_name)[0],fits.open(self.sci_ali_name)[0]
            self.sci_conv,self.ref_conv,self.ref_ali,self.sci_ali=self.sci_conv_hdu.data,self.ref_conv_hdu.data,self.ref_ali_hdu.data,self.sci_ali_hdu.data
            self.files_to_clean.append(self.sci_conv_name),self.files_to_clean.append(self.ref_conv_name)

            # self.sp_logger.info(self.sci_conv==self.ref_conv)
            # sys.exit(1)
            if self.args.show_plots:

                psf_fig,psf_ax = plt.subplots(figsize=(10,10),nrows=1,ncols=2)
                self.sci_norm = simple_norm(self.kernel_sci, 'log', percent=99.)
                self.ref_norm = simple_norm(self.kernel_ref, 'log', percent=99.)
                psf_ax[0].imshow(self.kernel_sci, norm=self.sci_norm, origin='lower', cmap='viridis')
                psf_ax[0].set_title('Science PSF')
                psf_ax[1].imshow(self.kernel_ref, norm=self.ref_norm, origin='lower', cmap='viridis')
                psf_ax[1].set_title('Reference PSF')
                plt.show()
                plt.close(psf_fig)



            return self.sci_conv,self.sci_ali,self.ref_conv,self.ref_ali
        
        else:
            # Read the REFERENCE image and convolve it with the Gaussian kernel
            self.sci_image_aligned=fits.open(self.sci_ali_name)
            self.sci_conv = scipy_convolve(self.sci_image_aligned[0].data,self.kernel_sci, mode='same', method='fft')
            self.sci_conv=np.nan_to_num(self.sci_conv)

            fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)

            self.sci_conv_hdu,self.sci_ali_hdu=fits.open(self.sci_conv_name)[0],fits.open(self.sci_ali_name)[0]
            self.sci_conv,self.sci_ali=self.sci_conv_hdu.data,self.sci_ali_hdu.data


            return self.sci_conv,self.sci_ali



    def py_convolve_images(self,plot_psfs=False,progress_bar=False):
        # bar.start()
        self.sp_logger.info(info_g+f" Beginning cross convolutions of PSFs and images")

        if not os.path.exists(self.path+'convolved_sci'):
            os.makedirs(self.path+'convolved_sci')
        if not os.path.exists(self.path+'convolved_ref'):
            os.makedirs(self.path+'convolved_ref')

        try:
            self.sci_conv_name=self.path+"convolved_sci/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'sci_convolved.fits'
        except:
            self.sci_img_hdu=self.sci_img_hdu[0]
            self.sci_conv_name=self.path+"convolved_sci/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'sci_convolved.fits'
        self.ref_conv_name=self.path+"convolved_ref/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'ref_convolved.fits'
        
        #########################################
        #BUILDING SCIENCE PSF
        self.sci_ali_hdu = fits.open(self.sci_ali_name)[0]
        self.sp_logger.info(info_g+' Measuring PSF of science image...')
        self.sci_ali_psf_built = False
        try:
            from psf_measure import measure_psf as _measure_psf
            _fwhm_guess = getattr(self, 'sci_seeing', 2.0) / getattr(self, 'sci_ps', 0.37)
            _sci_data = self.sci_ali_hdu.data.astype(float)
            if self.valid_mask is not None and not np.all(self.valid_mask):
                _sci_data = _sci_data.copy()
                _sci_data[~self.valid_mask] = np.nan
            _sci_res = _measure_psf(_sci_data, fwhm_guess=_fwhm_guess, min_stars=1, threshold_sigma=4.0)
            self.sci_ali_psf = _sci_res['psf']
            self.sp_logger.info(info_g+f" Science PSF: FWHM={_sci_res['fwhm']:.2f}px, n_stars={_sci_res['n_stars']}")
            self.sci_ali_psf_built = True
        except Exception as _e:
            self.sp_logger.warning(warn_y+f" measure_psf failed ({_e})")
        
        if self.sci_ali_psf_built == True:
            self.sp_logger.info(info_g+f" Science PSF built successfully")
            if plot_psfs!=False or self.args.show_plots==True:
                sci_psf_fig = plt.figure(figsize=(12,8))
                norm = simple_norm(self.sci_ali_psf, 'log', percent=99.)
                plt.imshow(self.sci_ali_psf, norm=norm, origin='lower', cmap='viridis')
                plt.title(f"Science PSF")
                plt.colorbar()
                plt.show()
                plt.close(sci_psf_fig)

        if self.sci_ali_psf_built==False:
            self.sp_logger.warning(warn_y+f" Unable to build empirical science PSF — falling back to Gaussian PSF from seeing estimate")
            try:
                from astropy.convolution import Gaussian2DKernel
                _fwhm_pix = self.sci_seeing / self.sci_ps
                _sigma    = _fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
                _gauss    = Gaussian2DKernel(_sigma, x_size=31, y_size=31)
                self.sci_ali_psf = _gauss.array / _gauss.array.sum()
                self.sci_ali_psf_built = True
                self.sp_logger.warning(warn_y+f" Gaussian science PSF: sigma={_sigma:.2f} pix (FWHM={_fwhm_pix:.2f} pix)")
            except Exception as _e:
                self.sp_logger.info(warn_r+f" Gaussian PSF fallback also failed: {_e} — exiting")
                self.sys_exit=True
                return
        else:
            self.ref_image_aligned=fits.open(self.ref_ali_name)
            self.ref_ali_img_hdu = self.ref_image_aligned[0]
            self.ref_conv = scipy_convolve(self.ref_ali_img_hdu.data, self.sci_ali_psf, mode='same', method='fft')
            fits.writeto(self.ref_conv_name, data=self.ref_conv, header=self.ref_image_aligned[0].header,overwrite=True)

            self.sp_logger.info(info_g+f" Convolved science PSF with reference image with size {self.ref_conv.shape}")
            self.sp_logger.info(info_g+f" Convolved reference image saved to {self.ref_conv_name}")

            
                
            

            if self.to_subtract!=False:
                #########################################
                #BUILDING REFERENCE PSF
                self.ref_ali_hdu = fits.open(self.ref_ali_name)[0]
                self.sp_logger.info(info_g+' Measuring PSF of reference image...')
                try:
                    from psf_measure import measure_psf as _measure_psf
                    _fwhm_ref = getattr(self, 'sci_seeing', 2.0) / getattr(self, 'sci_ps', 0.37)
                    _ref_data = self.ref_ali_hdu.data.astype(float)
                    if self.ref_valid_mask is not None and not np.all(self.ref_valid_mask):
                        _ref_data = _ref_data.copy()
                        _ref_data[~self.ref_valid_mask] = np.nan
                    _ref_res = _measure_psf(_ref_data, fwhm_guess=_fwhm_ref, min_stars=1, threshold_sigma=4.0)
                    self.ref_ali_psf = _ref_res['psf']
                    self.sp_logger.info(info_g+f" Reference PSF: FWHM={_ref_res['fwhm']:.2f}px, n_stars={_ref_res['n_stars']}")
                except Exception as _e:
                    self.sp_logger.warning(warn_y+f" measure_psf failed for reference ({_e})")

                # If the empirical reference PSF failed, fall back to a Gaussian using science seeing
                _ref_psf_built = hasattr(self,'ref_ali_psf') and self.ref_ali_psf is not None
                if not _ref_psf_built:
                    self.sp_logger.warning(warn_y+f" Unable to build empirical reference PSF — falling back to Gaussian PSF")
                    try:
                        from astropy.convolution import Gaussian2DKernel
                        _fwhm_pix = self.sci_seeing / self.sci_ps
                        _sigma    = _fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
                        _gauss    = Gaussian2DKernel(_sigma, x_size=31, y_size=31)
                        self.ref_ali_psf = _gauss.array / _gauss.array.sum()
                        self.sp_logger.warning(warn_y+f" Gaussian reference PSF: sigma={_sigma:.2f} pix (FWHM={_fwhm_pix:.2f} pix)")
                    except Exception as _e:
                        self.sp_logger.info(warn_r+f" Gaussian reference PSF fallback failed: {_e} — exiting")
                        self.sys_exit=True
                        return

                self.sp_logger.info(info_g+f" Reference PSF built succesfully")
                plot_psfs=True
                if plot_psfs!=False  or self.args.show_plots==True:
                    ref_psf_fig = plt.figure(figsize=(12,8))
                    norm = simple_norm(self.ref_ali_psf, 'log', percent=99.)
                    plt.imshow(self.ref_ali_psf, norm=norm, origin='lower', cmap='viridis')
                    plt.title(f"Reference PSF")
                    plt.colorbar()
                    plt.show()
                    plt.close(ref_psf_fig)
                # bar.update(1)

                self.sci_image_aligned=fits.open(self.sci_ali_name)
                self.sci_conv = scipy_convolve(self.sci_image_aligned[0].data, self.ref_ali_psf, mode='same', method='fft')
                self.sci_conv=np.nan_to_num(self.sci_conv)
                fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)
                self.sci_conv_hdu = fits.open(self.sci_conv_name)[0]
                self.ref_conv_hdu = fits.open(self.ref_conv_name)[0]

                self.sp_logger.info(info_g+f" Convolved reference PSF with science image with size {self.sci_conv.shape}")
                self.sp_logger.info(info_g+f" Convolved science image saved to {self.sci_conv_name}")

                

                self.kernel_sci,self.kernel_ref=self.sci_ali_psf,self.ref_ali_psf
                # self.sp_logger.info(self.ref_conv)
                # bar.finish()
                self.files_to_clean.append(self.ref_conv_name)
                self.files_to_clean.append(self.sci_conv_name)
                return self.sci_conv,self.kernel_sci,self.ref_conv,self.kernel_ref
            else:
                self.sci_image_aligned=fits.open(self.sci_ali_name)
                self.sci_conv = self.sci_image_aligned
                self.sci_conv=np.nan_to_num(self.sci_conv)

                fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)

                self.sci_conv_hdu,self.sci_ali_hdu=fits.open(self.sci_conv_name)[0],fits.open(self.sci_ali_name)[0]
                self.sci_conv,self.sci_ali=self.sci_conv_hdu.data,self.sci_ali_hdu.data

                self.kernel_sci,self.kernel_ref=self.sci_ali_psf,self.sci_ali_psf

                self.files_to_clean.append(self.sci_conv_name)
                return self.sci_conv,self.kernel_sci



    def gen_ref_cat(self):
        #show the fer_conv and sci_conv images
        self.sp_logger.info(info_g+f" Reference convolved image size: {self.ref_conv.shape}")
        self.sp_logger.info(info_g+f" Science convolved image size: {self.sci_conv.shape}")


        self.sp_logger.info(info_g+f" Generating reference catalogs")
        self.ref_width=self.sci_img_hdu.header['NAXIS2']*self.sci_ps/60

        if self.auto_cat=="auto":
            if self.use_sdss==False and (self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z'):
                self.ref_cat=panstarrs_query(ra_deg=round(self.sci_c.ra.deg,6),dec_deg=round(self.sci_c.dec.deg,6), rad_deg=round(self.ref_width/60.,6))

                self.stars=np.where((self.ref_cat['iMeanPSFMag']-self.ref_cat['iMeanKronMag']<0.05) & (self.ref_cat[str(self.sci_filt)+'MeanPSFMagErr']<0.06) & (self.ref_cat[str(self.sci_filt)+'MeanPSFMag']<23.5)  & (self.ref_cat[str(self.sci_filt)+'MeanPSFMag']!=-999.))[0]
                self.ref_cat=np.array(self.ref_cat[self.stars])
                if len(self.stars)==0:
                    self.sp_logger.warning(warn_r+f' {self.sci_obj} {self.sci_mjd} {self.sci_filt}: Length of PS1 reference catalog is 0')
                    self.sys_exit=True
                    return
                self.ref_coords_wcs_sky = SkyCoord(ra=self.ref_cat['raMean']*u.deg, dec=self.ref_cat['decMean']*u.deg,frame='fk5')
                self.ref_coords_wcs=np.column_stack((self.ref_cat['raMean'],self.ref_cat['decMean']))
                self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)
                # self.sp_logger.info(self.ref_coords_wcs)

            if self.sci_filt=='u' or self.use_sdss==True or self.args.sdsscat==True:
                # if not os.path.exists(
                self.ref_cat=sdss_query(ra_deg=round(self.sci_c.ra.deg,6),dec_deg=round(self.sci_c.dec.deg,6), rad_deg=round((self.ref_width),6))
                # self.sp_logger.info(self.ref_cat)
                self.stars=self.ref_cat

                if len(self.stars)==0:
                    self.sp_logger.warning(warn_r+f' {self.sci_obj} {self.sci_mjd} {self.sci_filt}: Length of SDSS reference catalog is 0')
                    self.sys_exit=True
                    return
                
                # self.sp_logger.info(len(self.stars),round(3600*(self.ref_width/60),6), 'arcsec search radius')
                self.ref_coords_wcs_sky = SkyCoord(ra=np.array(self.ref_cat['ra'])*u.deg, dec=np.array(self.ref_cat['dec'])*u.deg,frame='fk5')

                # self.sp_logger.info(len(self.ref_coords_wcs_sky))
                self.ref_coords_wcs=np.column_stack((self.ref_cat['ra'],self.ref_cat['dec']))
                # self.sp_logger.info(len(self.ref_coords_wcs))
                self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)
                # self.sp_logger.info(self.ref_coords_pix)
        else:
            #If a reference catalogue is passed in with it's full path as self.auto_cat
            if not self.auto_cat.startswith(self.path):
                self.auto_cat = self.path+self.auto_cat

            if not os.path.exists(self.auto_cat):
                self.sp_logger.info(warn_r+f' Reference catalogue {self.auto_ref} not found')
                
                self.sys_exit=True
                return

            # self.sp_logger.info(pd.read_csv(self.auto_cat,skiprows=1,header=None))
            self.sp_logger.info(info_g+f' Guessing format of reference catalog')
            #only need ra dec and mag
            # print(pd.read_csv(self.auto_cat))
            self.cat_tab = ascii.read(self.auto_cat)
            #convert the table to a pandas dataframe
            # try:
            try:
                self.ref_cat = self.cat_tab['ra','dec',self.sci_filt].to_pandas()
                self.ref_cat['mag'] = self.ref_cat[self.sci_filt]
            except:
                try:self.ref_cat = self.cat_tab['ra','dec','mag'].to_pandas()
                except:
                    self.sp_logger.info(warn_r+f' Reference catalogue {self.auto_cat} does not have the required columns')
                    self.sys_exit=True
                    return

            # print(type(self.ref_cat['ra'].iloc[0]))
            self.ref_cat = self.ref_cat.loc[self.ref_cat['mag']!='-'].reset_index(drop=True)
            self.ref_cat['ra'] = self.ref_cat['ra'].astype(float)*u.deg
            self.ref_cat['dec'] = self.ref_cat['dec'].astype(float)*u.deg


            self.stars=self.ref_cat
            # self.ref_coords_wcs_sky = np.column_stack((self.ref_cat["RA"],self.ref_cat["DEC"]))
            
            self.ref_coords_wcs_sky = SkyCoord(ra=self.ref_cat["ra"],dec=self.ref_cat["dec"],frame='fk5',unit='deg')
            self.ref_coords_wcs = self.ref_coords_wcs_sky
            self.ref_coords_wcs=np.column_stack((self.ref_cat['ra'],self.ref_cat['dec']))


            self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)

            # self.ref_coords_pix = np.column_stack((self.ref_cat["xpos"],self.ref_cat["ypos"]))
            # 
        # print(self.ref_coords_wcs)
        # print(self.sci_conv_name)
        # print(self.ref_conv_name)
        # sys.exit()

        self.sp_logger.info(info_g+' Catalog stars in PS1/SDSS found='+str(len(self.stars))+','+str(round(3600*(self.ref_width/60),6))+ 'arcsec search radius')
        # self.sp_logger.info(self.ref_coords_pix)
        # self.sp_logger.info(self.sci_conv)
        self.mean,self.median, self.std = sigma_clipped_stats(self.sci_conv, sigma=3.0)
        self.sp_logger.info(info_g+' Mean, Median, Std of sci_conv: '+str(round(self.mean,6))+' '+str(round(self.median,6))+' '+str(round(self.std,6)))
        self.sp_logger.info(info_g+' Detecting stars with SExtractor')

        self.sp_logger.info(info_g+' Detecting stars with IRAFStarFinder')
        self.iraffind= IRAFStarFinder(threshold=abs(starscale*self.std),fwhm=3.0,roundhi=0.3)

        self.sp_logger.info(info_g+' Threshold for detecting stars: '+str(int(starscale*self.std)))
        
        
        if np.shape(self.sci_conv)[0]>1100:
            self.sources = self.iraffind(self.sci_conv[60:len(self.sci_conv)-60,60:len(self.sci_conv)-60] - self.median)
            self.star_coords_pix=np.column_stack((self.sources['xcentroid']+60.,self.sources['ycentroid']+60.))
        else:
            self.sources = self.iraffind(self.sci_conv - self.median)
            self.star_coords_pix=np.column_stack((self.sources['xcentroid'],self.sources['ycentroid']))

        self.matched_catalog_mag=[]


    
        # print(self.sci_ali_name)
        self.star_coords_wcs=load_wcs_from_file(filename=self.sci_conv_name,coord=self.star_coords_pix)
        # self.sp_logger.info(self.star_coords_wcs)
        self.star_coords_wcs_sky = SkyCoord(ra=self.star_coords_wcs[:,0]*u.deg, dec=self.star_coords_wcs[:,1]*u.deg,frame='fk5')

        if self.special_case!=None:
            if any(phrase == self.special_case for phrase in ['SN2023ixf','bright','brightstars.cat','sn2023ixf','2023ixf','brightstars']):
                #find crossover between self.star_coords_wcs_sky and bright_stars_sc
                self.bright_stars_sc_wcs = SkyCoord(ra=self.bright_stars_sc['ra']*u.deg, dec=self.bright_stars_sc['dec']*u.deg,frame='fk5')

                self.indx, self.d2d, self.d3d =self.star_coords_wcs_sky.match_to_catalog_sky(self.bright_stars_sc_wcs)
                #where stars match by search rad in arcseconds!
                self.upd_indx=np.where(self.d2d<search_rad/3600.*u.deg)[0]
                # self.sp_logger.info(self.indx[self.upd_indx])
                self.sp_logger.info(info_g+' Bright stars found= '+str(len(self.upd_indx))+','+str(round(3600*(self.ref_width/60),6))+ 'arcsec search radius')
                self.star_coords_pix=self.star_coords_pix[self.upd_indx]
                self.star_coords_wcs=self.star_coords_wcs[self.upd_indx]
                self.star_coords_wcs_sky=self.star_coords_wcs_sky[self.upd_indx]

                self.matched_ref_coords_pix=self.bright_stars_sc['xcentroid','ycentroid'][self.indx[self.upd_indx]]
                self.matched_star_coords_pix=self.star_coords_pix

                self.matched_catalog_mag=self.bright_stars_sc[self.sci_filt][self.indx[self.upd_indx]]
                self.matched_catalog = self.bright_stars_sc[self.indx[self.upd_indx]]

                # for i in range(len(self.matched_ref_coords_pix)):
                #     self.sp_logger.info('sci',self.star_coords_wcs[i])
                #     self.sp_logger.info('ref',self.matched_catalog['ra','dec'][i])
                #     self.sp_logger.info(self.matched_star_coords_pix[i])
                #     self.sp_logger.info(self.matched_ref_coords_pix[i])
                #     self.sp_logger.info()
                # # self.sp_logger.info(self.bright_stars_sc[self.indx[self.upd_indx]])
                # sys.exit()

        else:
            #find crossover between panstarrs ad reference images
            self.sp_logger.info(info_g+' Searching for stars in reference catalog')
            self.indx, self.d2d, self.d3d =self.star_coords_wcs_sky.match_to_catalog_sky(self.ref_coords_wcs_sky)
            #where stars match by search rad in arcseconds!
            # self.upd_indx=np.where(self.d2d<=7.5/3600.*u.deg)[0]
            self.upd_indx=np.where(self.d2d<=1/3600.*u.deg)[0]
            d2d_ = self.d2d[self.upd_indx]
            if len(self.upd_indx)<=8:
                self.sp_logger.info(warn_y+f' {len(self.upd_indx)}(<=7) stars found in reference catalog, increasing search radius: 1->2')
                self.upd_indx=np.where(self.d2d<=2/3600.*u.deg)[0]
                d2d_ = self.d2d[self.upd_indx]
                if len(self.upd_indx)<=8:
                    self.sp_logger.info(warn_y+f' {len(self.upd_indx)}(<=7) stars found in reference catalog, increasing search radius again: 2->5')
                    self.upd_indx=np.where(self.d2d<=5/3600.*u.deg)[0]
                    d2d_ = self.d2d[self.upd_indx]
                    if len(self.upd_indx)<=9:
                        self.sp_logger.info(warn_y+f' {len(self.upd_indx)}(<=7) stars found in reference catalog, increasing search radius again: 5->7')
                        self.upd_indx=np.where(self.d2d<=7./3600.*u.deg)[0]
                        
                        d2d_ = self.d2d[self.upd_indx]
            # self.sp_logger.info(self.indx[self.upd_indx])
            self.sp_logger.info(info_g+' Catalog stars in PS1/SDSS found = '+str(len(self.upd_indx))+', in a '+str(round(3600*(7.5/60),6))+ 'arcsec search radius')

        
            # print(self.upd_indx)
            # print(self.indx)
            self.matched_ref_coords_pix=self.ref_coords_pix[self.indx[self.upd_indx]]
            self.matched_star_coords_pix=self.star_coords_pix[self.upd_indx]
            # print(self.upd_indx)
            # print(self.indx)

            # [print(f'physical;circle({self.matched_ref_coords_pix[i][0]},{self.matched_ref_coords_pix[i][1]},20)'+' # text={'+f'xy={int(self.matched_ref_coords_pix[i][0]),int(self.matched_ref_coords_pix[i][1])}'+'}') for i in range(len(self.matched_ref_coords_pix))]
            # print()
            # [print(f'physical;circle({self.matched_star_coords_pix[i][0]},{self.matched_star_coords_pix[i][1]},20)'+' # text={'+f'xy={int(self.matched_star_coords_pix[i][0]),int(self.matched_star_coords_pix[i][1])}'+'}') for i in range(len(self.matched_star_coords_pix))]
            # sys.exit()
            # self.matched_fwhm=self.fwhm[self.upd_indx]

            
            # self.sp_logger.info(self.sci_filt=='r')
            # self.sp_logger.info(self.matched_ref_coords_pix,self.matched_star_coords_pix)
            if self.auto_cat == 'auto' and self.use_sdss==False and self.args.sdsscat==False and (self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z'):
                self.matched_catalog=self.ref_cat[self.indx[self.upd_indx]]
                # self.sp_logger.info(self.matched_catalog)
                # sys.exit()
                print(1)
                self.matched_catalog_mag=self.matched_catalog[str(self.sci_filt)+'MeanPSFMag']

                # for i in range(len(self.matched_catalog_mag)):
                #     self.sp_logger.info(self.matched_catalog[i])
                #     self.sp_logger.info(self.bright_stars_sc[i])
                #     self.sp_logger.info()


            elif self.auto_cat == 'auto' and (self.sci_filt=='u' or self.use_sdss==True or self.args.sdsscat==True):
                print(2)
                self.string_band='psfMag_'+str(self.sci_filt)
                self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
                self.matched_catalog_mag=np.asarray(self.matched_catalog['mag'])
                self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

            elif self.sci_filt=='u' or self.use_sdss==True or self.args.sdsscat==True:
                print(0)
                self.string_band='psfMag_'+str(self.sci_filt)
                self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
                self.matched_catalog_mag=np.asarray(self.matched_catalog['mag'])
                self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

            elif self.auto_cat !='auto':
                print(3)
                # self.sp_logger.info(self.ref_cat.loc[self.indx[self.upd_indx]])
                self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
                self.matched_catalog_mag=np.asarray([float(m) for m in self.matched_catalog['mag']])
                self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

        

        #remove_duplicates
        # self.matched_ind_uniq = np.unique(self.matched_star_coords_pix, axis=1, return_index=True)[1]

        # print(self.matched_ind_uniq,len(self.matched_ind_uniq),len(self.matched_star_coords_pix))

        # self.matched_star_coords_pix = self.matched_star_coords_pix[self.matched_ind_uniq]
        # self.matched_ref_coords_pix = self.matched_ref_coords_pix[self.matched_ind_uniq]
        # self.matched_catalog_mag = self.matched_catalog_mag[self.matched_ind_uniq]
        # print(self.matched_star_coords_pix)
        # print(self.matched_ref_coords_pix)

    
    

        self.sp_logger.info(info_g+" Stars detected in the image= "+str(len(self.star_coords_pix)))
        self.sp_logger.info(info_g+" Length of matched catalog= "+str(len(self.matched_catalog_mag)))



        if len(self.matched_catalog_mag)<5:
            self.sp_logger.info(warn_y+" Few stars ("+str(len(self.matched_catalog_mag))+") matched!")

    

        if (len(self.matched_catalog_mag)<1 and self.special_case==None) or (self.special_case!=None and len(self.matched_catalog_mag)<1):
            # self.sp_logger.info(self.special_case, len(self.matched_catalog_mag))
            self.sp_logger.warning(warn_r+f" {self.sci_obj} {self.sci_mjd} {self.sci_filt}: Less than 2 matched calibration stars ")

            self.sp_logger.warning(warn_r+" Exiting: not enough stars to calibrate!")
            self.sys_exit=True
            return
            

        self.sp_logger.info(info_g+' Finished the star matching process')

        keep_indx=[]
        new_matched_catalog_mag,new_matched_catalog,new_matched_star_coords_pix=[],[],[]
        # self.sp_logger.info(self.matched_catalog_mag)
        # [self.sp_logger.info(self.matched_catalog_mag[k]) for k in range(len(self.matched_catalog_mag)) if self.matched_catalog_mag[k]!='-' and float(self.matched_catalog_mag[k])>13.5]
        [keep_indx.append(k) for k in range(len(self.matched_catalog_mag)) if self.matched_catalog_mag[k]!='-' and float(self.matched_catalog_mag[k])>13.5]
        # self.sp_logger.info(keep_indx)
        for i in keep_indx:
            new_matched_catalog_mag.append(float(self.matched_catalog_mag[i]))
            try:new_matched_catalog.append(self.matched_catalog[i]) 
            except: new_matched_catalog.append(self.matched_catalog.iloc[i])
            new_matched_star_coords_pix.append(self.matched_star_coords_pix[i])
        self.sp_logger.info(info_g+' Number of stars kept after magnitude cut: '+str(len(new_matched_catalog_mag)))
        # self.sp_logger.info(new_matched_catalog)
        # self.matched_catalog_mag,self.matched_catalog,self.matched_star_coords_pix=np.array(new_matched_catalog_mag),np.array(new_matched_catalog),np.array(new_matched_star_coords_pix)
        self.matched_catalog_mag,self.matched_star_coords_pix=np.array(new_matched_catalog_mag),np.array(new_matched_star_coords_pix)
        
        #remove stars that fall within an ellipse centred at the centre of the image with 

        
        return self.ref_cat,self.matched_catalog,self.matched_catalog_mag,self.matched_star_coords_pix




    def combine_psf(self):

        self.sp_logger.info(info_g+f" Combining PSFs")
        ###################################################################
        #  Combine science PSF with reference PSF, flatten PSF to 1D and gaussian fit

        if not os.path.exists(self.path+'convolved_psf'):
            os.makedirs(self.path+'convolved_psf')
        
        # try:
        #     if os.path.exists(self.path+'convolved_psf'):
        #         os.system('rm '+self.path+'convolved_psf/*.fits')
        # except:
        #     pass

        self.comb_psf = scipy_convolve(self.kernel_sci, self.kernel_ref, mode='same', method='fft')
        self.comb_psf=self.comb_psf/np.sum(self.comb_psf)
        self.hdu_comb_psf_name =self.path+"convolved_psf/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][0:13]+f"comb_psf_{self.rand_nums_string}.fits"
        self.hdu_comb_psf= fits.PrimaryHDU(self.comb_psf)
        self.sp_logger.info(info_g+f" Saving combined PSF to {self.hdu_comb_psf_name}")
        self.hdu_comb_psf.writeto(self.hdu_comb_psf_name,overwrite=True)




    def cutout_psf(self,data,psf_array,xpos,ypos):
        all_cutouts=[]
        for d in range(len(xpos)):
            # self.sp_logger.info('x',xpos,'y',ypos)
            xcutout,ycutout=np.shape(psf_array)[0],np.shape(psf_array)[1]
            position=(xpos[d],ypos[d])
            size = (xcutout, ycutout)
            # self.sp_logger.info('size',size,'position',position,'data shape',np.shape(data))

            cutout = Cutout2D(data, position, size)#, mode='partial')
            all_cutouts.append(cutout.data)
        return np.array(all_cutouts)

    def psf_fit(self,data_cutout,psf_array):
        params = []
        for d in range(len(data_cutout)):
            #fit PSF with x or y-shift
            # shift cutout psf to be aligned with psf model
            # restrict it to move <5 pixels in each direction so the fit doesn't go wild and fit a nearby bright star etc.
            xoff, yoff, exoff, eyoff = chi2_shift(psf_array,data_cutout[d], 10,return_error=True, upsample_factor='auto')
            #convert to arcsec
            xoff_arc=abs(xoff)*self.sci_ps
            yoff_arc=abs(yoff)*self.sci_ps
            if (xoff>5.0 or yoff>5.0) and self.telescope not in SEDM:
                xoff,yoff=0.0,0.0
            data_cutout_shift=scipy.ndimage.shift(data_cutout[d], [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
            resize_sci=np.reshape(data_cutout_shift,np.shape(data_cutout_shift[d])[0]*np.shape(data_cutout_shift)[1])
            #resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1],1)
            resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1]) 
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)
            params.append([slope, intercept, r_value*r_value,xoff, yoff, abs(xoff_arc), abs(yoff_arc)])

        return np.array(params)

    # import numpy as np
    # from astropy.stats import SigmaClip

    def psf_fit_bg(self,data, psf_array,sn_x,sn_y,shift=True):
        from photutils import Background2D, MedianBackground,BkgZoomInterpolator
        from astropy.modeling import models, fitting
        params = []
        for d in range(len(data)):
            sigma_clip = SigmaClip(sigma=3.)
            bkg_estimator = SExtractorBackground(sigma_clip=sigma_clip)
            box = int(np.ceil(np.shape(psf_array)[0]/4))
            if float(box)%2==0:box+=1
            # else:
            filter_size = (5*int(box), 5*int(box))

            app = CircularAperture([sn_x,sn_y], r=1.5*box)
            masks = app.to_mask(method='center') 
            mask_arr = masks.to_image(shape=data[d].shape).astype(bool)
            data_cutout = self.cutout_psf(data=data[d],psf_array=self.comb_psf,xpos=[sn_x],ypos=[sn_y])[0]
            
            self.sp_logger.info(info_g+' Estimating spatially varying background')
            self.sp_logger.info(info_g+' Box size: '+str(box))
            self.sp_logger.info(info_g+' Filter size: '+str(filter_size))
            self.sp_logger.info(info_g+' Mask size: '+str(np.shape(mask_arr)))


            # try:
            bkg1 = Background2D(data[d], box_size=box, filter_size=filter_size,
                                sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,interpolator=BkgZoomInterpolator(order=1))
                
                
            surface_function_init = models.Polynomial2D(degree=2)
            fitter = fitting.LevMarLSQFitter()
            fit_x,fit_y = np.arange(0,np.shape(data[d])[1]),np.arange(0,np.shape(data[d])[0])
            fit_xx,fit_yy = np.meshgrid(fit_x,fit_y)
            surface_fit = fitter(surface_function_init, fit_xx, fit_yy, data[d])
            surface= surface_fit(fit_xx, fit_yy)

            # import properimage.single_image as si

            # si_bkg_poly = si.SingleImage(self.sci_scale_sub_name).background
            import sep

            bkg = sep.Background(np.array(data[d]).astype(np.float64)).back()
            bkg = self.cutout_psf(data=bkg,psf_array=self.comb_psf,xpos=[sn_x],ypos=[sn_y])[0]

            data_cutout_mbkg = data_cutout - bkg
            #save the background subtracted image
            fits.writeto(self.path+'scaled_subtracted_imgs/'+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][0:13]+f"complex_bkg_subtracted_{self.rand_nums_string}.fits", data=data_cutout_mbkg, header=self.sci_img_hdu.header,overwrite=True)
            
            #save the surface fit
            data_cutout_msurf = data_cutout - self.cutout_psf(data=surface,psf_array=self.comb_psf,xpos=[sn_x],ypos=[sn_y])[0]
            fits.writeto(self.path+'scaled_subtracted_imgs/'+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][0:13]+f"surface_fit_{self.rand_nums_string}.fits", data=data_cutout_msurf, header=self.sci_img_hdu.header,overwrite=True)

            data_cutout = data_cutout - bkg
            xoff, yoff, exoff, eyoff = chi2_shift(psf_array, data_cutout, 10, 
                                                    return_error=True, upsample_factor='auto')

            xoff_arc=abs(xoff)*self.sci_ps
            yoff_arc=abs(yoff)*self.sci_ps

            if shift!=False:
                if (xoff>5.0 or yoff>5.0) and self.telescope not in SEDM:
                    xoff,yoff=0.0,0.0
                data_cutout_shift=scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
                resize_sci=np.reshape(data_cutout_shift,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1])
                resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1]) 
            else:
                resize_sci=np.reshape(data_cutout,np.shape(data_cutout)[0]*np.shape(data_cutout)[1])
                resize_psf=np.reshape(psf_array,np.shape(data_cutout)[0]*np.shape(data_cutout)[1])

            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)
            params.append([slope, intercept, r_value*r_value,xoff, yoff, abs(xoff_arc), abs(yoff_arc)])

        return np.array(params)

    def psf_fit_noshift(self,data_cutout,psf_array):
        params = []
        for d in range(len(data_cutout)):
            # forced photometry - fit with no x or y-shift for limits calculations
            # print(np.shape(data_cutout[d]))
            resize_sci=np.reshape(data_cutout[d],np.shape(data_cutout[d])[0]*np.shape(data_cutout[d])[1])
            resize_psf=np.reshape(psf_array,np.shape(data_cutout[d])[0]*np.shape(data_cutout[d])[1])

            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)
            params.append([slope, intercept, r_value*r_value])
        return np.array(params)



    def get_zeropts(self):
        if self.termoutp!='quiet':
            self.sp_logger.info(info_g+f' Calculating zeropoints')
            # self.sp_logger.info(info_g+' Size of sci_conv: '+str(np.shape(self.sci_conv)))
            # self.sp_logger.info(info_g+' Size of ref_conv: '+str(np.shape(self.ref_conv)))
            # self.sp_logger.info(info_g+' Size of comb_psf: '+str(np.shape(self.comb_psf)))
        ###########################################################
        # Zero points!

        self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref=[],[],[],[]
        # self.sp_logger.info(self.matched_star_coords_pix)

        if len(self.matched_star_coords_pix)==0:
            #check if args.relative_flux is set to True
            if self.relative_flux:
                #if so, then we set a default zeropoint of 1
                self.sp_logger.warning(info_b+f' No matched stars found in {self.sci_obj} image, but relative flux is set to True, so calculating relative fluxes and ignoring zeropoint')
                self.zp_sci.append(1)
                self.rsq_sci.append(1)
                return
            else:
                self.sp_logger.warning(warn_r+f' No matched stars found in {self.sci_obj} image, exiting')
                self.sys_exit=True
                return

        if len(self.matched_star_coords_pix)>0:
            #check the countsto reject saturated stars in the science image
            self.counts = []

            for i in range(0,len(self.matched_star_coords_pix)):
                self.cutout_sci=self.cutout_psf(data=self.sci_conv,psf_array=self.comb_psf,xpos=[self.matched_star_coords_pix[i][0]],ypos=[self.matched_star_coords_pix[i][1]])[0]
                # self.sp_logger.info(info_g+' Cutout shape',np.shape(self.cutout_sci))
                if np.shape(self.cutout_sci)!=np.shape(self.comb_psf):
                    #pad with 0s to make it the same size as the PSF
                    self.cutout_sci = np.pad(self.cutout_sci,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_sci)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_sci)[1])),'constant',constant_values=0)
                self.sci_psf_fit = self.psf_fit_noshift(data_cutout=[self.cutout_sci],psf_array=self.comb_psf)[0]
                self.counts.append(self.sci_psf_fit[0])
            self.counts=np.array(self.counts)
            self.counts_o = self.counts
            self.counts=self.counts[~np.isnan(self.counts)]
            self.counts=self.counts[~np.isinf(self.counts)]
            self.counts=self.counts[self.counts>0]
            # self.sp_logger.info(self.counts)
            # self.sp_logger.info(self.counts[np.where(self.counts<=50000)])

            # print(self.special_case==None)
            if self.telescope not in SEDM and self.special_case==None and len(self.counts)>3:
                self.sp_logger.info(info_g+' Number of stars before saturation check: '+str(len(self.matched_star_coords_pix)))
                c=0
                if self.sci_filt in ['g','r','z']:sat_mag,count_lim = 13.5,47000
                elif self.sci_filt in ['i']:sat_mag,count_lim = 15.0,40000
                else:sat_mag,count_lim = 16.5,45000

                # count_lim=40000
                # sat_mag=12.5

                if len(self.matched_star_coords_pix)>=20: keep_min=20
                elif len(self.matched_star_coords_pix)>=10: keep_min=10
                elif len(self.matched_star_coords_pix)>=5: keep_min=5
                else: keep_min=1
                while True: #iterate over this and reduce u-band sat mag by 0.5 2 times if len matched catalog>=20
                    self.sp_logger.info(info_g+f' #### Iteration {c} ####')
                    self.sp_logger.info(info_g+' Saturation magnitude:'+str(sat_mag))
                    self.sp_logger.info(info_g+' Count limit:'+str(count_lim))
                    # self.sp_logger.info(info_g+' Keep min:',keep_min)
                    self.matched_new_pix,self.matched_new_mag = [],[]
                    # print(self.counts_o)
                    for k in range(0,len(self.counts_o)):
                        # print(self.counts_o[k])
                        # self.sp_logger.info(info_g+f' Star {k} kept : mag:{self.matched_catalog_mag[k]:.2f} x: {int(self.matched_star_coords_pix[k][0])} y: {int(self.matched_star_coords_pix[k][1])}')
                        if not np.isnan(self.counts_o[k]) and not np.isinf(self.counts_o[k]) and self.counts_o[k]>0 and self.counts_o[k]<=count_lim:
                            if self.matched_catalog_mag[k]>sat_mag:
                                self.matched_new_pix.append(self.matched_star_coords_pix[k])
                                self.matched_new_mag.append(self.matched_catalog_mag[k])
                                # self.sp_logger.info(info_g+f' Star {k} kept : mag:{self.matched_catalog_mag[k]:.2f} x: {int(self.matched_star_coords_pix[k][0])} y: {int(self.matched_star_coords_pix[k][1])}')

                    if len(self.matched_new_pix)>=keep_min or c==3 or count_lim>=47000:
                        break
                    else:
                        if keep_min!=1:
                            keep_min-=1
                        if self.sci_filt in ['u'] and len(self.matched_star_coords_pix)>3:
                            sat_mag-=0.5
                            # count_lim+=5000
                        else:
                            sat_mag-=0.5
                            # count_lim+=5000

                    c+=1
                self.matched_star_coords_pix = np.array(self.matched_new_pix)
                self.matched_catalog_mag = np.array(self.matched_new_mag)

                self.sp_logger.info(info_g+' Number of stars after saturation check: '+str(len(self.matched_star_coords_pix)))


            # self.sp_logger.info(self.special_case)
            if self.special_case!=None and any(phrase in self.special_case for phrase in ['SN2023ixf','bright','brightstars.cat','sn2023ixf','2023ixf','brightstars']):
                self.sp_logger.info(info_g+f' Checking for stars brighter than: {self.bright_cat_mag_lim}')

            #remove duplicates
            
            self.uniq_inds = np.unique(self.matched_star_coords_pix,axis=0,return_index=True)[1]
            self.matched_star_coords_pix = self.matched_star_coords_pix[self.uniq_inds]
            self.matched_catalog_mag = self.matched_catalog_mag[self.uniq_inds]
            self.sp_logger.info(info_g+' Number of stars after removing duplicates: '+str(len(self.matched_star_coords_pix)))

            self.zp_size_check=1
            self._zp_loop_attempts=0
            while self.zp_size_check!=0 and self._zp_loop_attempts<3:
                for i in range(0,len(self.matched_star_coords_pix)):

                    if self.special_case!=None:
                        if any(phrase in self.special_case for phrase in ['SN2023ixf','bright','brightstars.cat','sn2023ixf','2023ixf','brightstars']):
                            if self.matched_catalog_mag[i]>self.bright_cat_mag_lim+0.35:continue
                    
                    self.cutout_sci=self.cutout_psf(data=self.sci_conv,psf_array=self.comb_psf,xpos=[self.matched_star_coords_pix[i][0]],ypos=[self.matched_star_coords_pix[i][1]])[0]
                    self.cutout_ref=self.cutout_psf(data=self.ref_conv,psf_array=self.comb_psf,xpos=[self.matched_star_coords_pix[i][0]],ypos=[self.matched_star_coords_pix[i][1]])[0]
                    if np.shape(self.cutout_sci)!=np.shape(self.comb_psf) or np.shape(self.cutout_ref)!=np.shape(self.comb_psf):continue
                        # self.sp_logger.info(info_g+f' Cutout shape mismatch for star {i}, padding')
                        # if np.shape(self.cutout_sci)!=np.shape(self.comb_psf):
                        #     self.sp_logger.info(info_g+f' Cutout shape mismatch for science image')
                        #     # print(np.shape(self.cxutout_sci),np.shape(self.comb_psf))
                        #     # if self.zp_size_check==2:
                        #     self.cutout_sci = np.pad(self.cutout_sci,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_sci)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_sci)[1])),'constant',constant_values=0)
                        #     # print(np.shape(self.cutout_sci),np.shape(self.comb_psf))
                        # if np.shape(self.cutout_ref)!=np.shape(self.comb_psf):
                        #     self.sp_logger.info(info_g+f' Cutout shape mismatch for reference image')
                        #     # if self.zp_size_check==2:
                        #     self.cutout_ref = np.pad(self.cutout_ref,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_ref)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_ref)[1])),'constant',constant_values=0)

                    
                    # self.cutout_ref = np.pad(self.cutout_ref,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_ref)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_ref)[1])),'constant',constant_values=0)
                    # print(np.shape(self.cutout_sci),np.shape(self.cutout_ref),np.shape(self.comb_psf))
                    self.sci_psf_fit = self.psf_fit(data_cutout=[self.cutout_sci],psf_array=self.comb_psf)[0]
                    self.zp_sci.append(2.5*np.log10(self.sci_psf_fit[0])+self.matched_catalog_mag[i])
                    self.rsq_sci.append(self.sci_psf_fit[2])


                    self.ref_psf_fit = self.psf_fit(data_cutout=[self.cutout_ref],psf_array=self.comb_psf)[0]
                    self.zp_ref.append(2.5*np.log10(self.ref_psf_fit[0])+self.matched_catalog_mag[i]) 
                    self.rsq_ref.append(self.ref_psf_fit[2])
            
                if len(self.zp_sci)!=0:
                    self.zp_size_check=0
                    break
                else:
                    self.zp_size_check=2
                    self._zp_loop_attempts+=1
                

        self.zp_sci,self.zp_ref=np.array(self.zp_sci),np.array(self.zp_ref)
        self.rsq_ref,self.rsq_sci=np.array(self.rsq_ref),np.array(self.rsq_sci)
        # print(self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref)

        if self.sci_filt=='u':
            self.thresh=0.05
        if self.sci_filt!='u':
            self.thresh=0.75
        if self.telescope in SEDM:
            self.thresh=0.1
            if self.sci_filt=='u':
                self.thresh=0.5

        print(self.zp_sci)
        print(self.rsq_sci)
        print(self.zp_ref)
        print(self.rsq_ref)
        # sys.exit()
        # if len(self.zp_sci) and len(self.zp_ref)>=10:
        #     #remove the smallest and largest values
        #     self.zp_sci,self.zp_ref=sigma_clip(self.zp_sci,sigma=5,maxiters=4),sigma_clip(self.zp_ref,sigma=5,maxiters=4)

        if np.count_nonzero(~np.isnan(self.zp_ref))<2 and self.special_case==None:
            if np.count_nonzero(~np.isnan(self.zp_sci))<3 and self.telescope not in SEDM:
                self.sp_logger.info(warn_y+' Few stars in the field')
            if np.count_nonzero(~np.isnan(self.zp_sci))<3 and self.telescope in SEDM:
                self.sp_logger.info(warn_y+' Few stars in the field, continuing if 1 or 2 stars in science')
                


             #errors.write('Few (< 3) field stars to calibrate in science and reference image. \n')
            if np.count_nonzero(~np.isnan(self.zp_sci))>3:
                self.sp_logger.info(warn_r+' Something wrong with reference since no zero point stars detected, but detected in science')
                self.sp_logger.info(warn_r+' Exiting..')

        self.zp_df = pd.DataFrame({'zp_sci':self.zp_sci,'zp_ref':self.zp_ref,'rsq_sci':self.rsq_sci,'rsq_ref':self.rsq_ref}) #show in reverse order
        self.sp_logger.info(tabulate(self.zp_df.loc[self.zp_df['rsq_sci']>0.5].sort_values(by=['rsq_sci','rsq_ref'],ascending=False), headers='keys', tablefmt='psql',showindex=False))
        # print(tabulate(self.zp_df.sort_values(by=['rsq_sci','rsq_ref']), headers='keys', tablefmt='psql'))
        # for i in range(len(self.zp_sci)):
        #     print(self.zp_sci[i],self.zp_ref[i],self.rsq_sci[i],self.rsq_ref[i])

        self.rsq_cont=False
        #check how many stars have rsq values above the threshold and lower by 0.05 until there are at least 2 stars
        # self.thresh=0.4
        if self.sci_filt=='u' and self.telescope not in SEDM:
            self.zp_sci_lim=27.5
        while self.rsq_cont==False:
            self.arr=(self.rsq_ref >=self.thresh) & (self.rsq_sci >=self.thresh) & (self.zp_sci >=0) & (self.zp_ref >=0) & (self.zp_sci<=self.zp_sci_lim)#filter out stars poorly fitted with PSF fit & weird zp measurement
            # self.sp_logger.info(self.rsq_ref[self.arr],sum(self.arr)) 
            if (len([i for i in self.rsq_ref if i>=self.thresh])<1 or len([i for i in self.rsq_sci if i>=self.thresh])<1) and len(self.zp_sci)>=5:
                self.sp_logger.info(warn_y+f' Lowering threshold for rsq values to {self.thresh-0.05:.2f}')
                self.thresh-=0.05
            elif (len([i for i in self.rsq_ref if i>=self.thresh])<1 or len([i for i in self.rsq_sci if i>=self.thresh])<1) and len(self.zp_sci)<5:
                self.sp_logger.info(warn_y+f' No stars with rsq values above threshold of {self.thresh:.2f}. Lowering threshold for rsq values to {self.thresh-0.05:.2f}')
                self.thresh-=0.05
            else:
                self.rsq_cont=True
                break

            if self.thresh<0.45 and self.sci_filt!='u':
                self.sp_logger.warning(warn_r+' No stars with rsq values above threshold')
                self.sp_logger.warning(warn_r+' Exiting..')
                self.sys_exit=True
                self.rsq_cont=True
                break
            if self.thresh<0.1:
                self.sys_exit=True
                self.sp_logger.warning(warn_r+' No stars with rsq values above threshold, exiting..')
                return

        self.arr=(self.rsq_ref >=self.thresh) & (self.rsq_sci >=self.thresh) & (self.zp_sci >=0) & (self.zp_ref >=0) & (self.zp_ref<=34)#filter out stars poorly fitted with PSF fit & weird zp measurement
        # print(self.rsq_ref>=self.thresh)
        # print(self.rsq_sci>=self.thresh)
        # print(self.zp_sci>=0)
        # print(self.zp_ref>=0)

        if all(val ==False for val in self.arr)==True:
            weird_img = 'science and reference'
            if all(val_sci_rsq <self.thresh for val_sci_rsq in self.rsq_sci)==True:
                self.weird_img = 'science'
            else:
                self.weird_img = 'reference'


            self.sp_logger.warning(warn_r+f' All stars in the {self.weird_img} images fit poorly or returned anomalous zeropoint measurements!')
            self.sp_logger.warning(warn_r+f' {self.sci_obj} {self.sci_filt} {self.sci_mjd} {self.weird_img} images fit poorly or returned anomalous zeropoint measurements!')

            self.sp_logger.warning(warn_r+" Exiting...")
            self.sys_exit=True
            return

        self.rsq_ref,self.rsq_sci = self.rsq_ref[self.arr],self.rsq_sci[self.arr]
        self.zp_sci,self.zp_ref = self.zp_sci[self.arr],self.zp_ref[self.arr]
        self.sp_logger.info(info_g+' Number of stars after rsq cut: '+str(len(self.zp_sci)))

        self.zp_sci_new,self.zp_ref_new=sigma_clip(self.zp_sci,sigma=3,maxiters=4),sigma_clip(self.zp_ref,sigma=3,maxiters=4)
        self.zp_sci_new=self.zp_sci_new[~self.zp_sci_new.mask]
        self.zp_ref_new=self.zp_ref_new[~self.zp_ref_new.mask]
        if len(self.zp_sci_new)==len(self.zp_sci) or len(self.zp_ref_new)==len(self.zp_ref):
            self.sp_logger.info(warn_y+' No stars removed after sigma clipping')
            if len(self.zp_sci)>=5 and len(self.zp_ref)>=5 and self.sci_filt!='u':
                self.sp_logger.info(info_g+f' Removing stars with zp values outside 5th and 95th percentiles')
                self.zp_sci_new,self.zp_ref_new = [i for i in self.zp_sci if i>np.percentile(self.zp_sci,5) and i<np.percentile(self.zp_sci,95)],[i for i in self.zp_ref if i>np.percentile(self.zp_ref,5) and i<np.percentile(self.zp_ref,95)]
            
        self.zp_sci,self.zp_ref = np.array(self.zp_sci_new),np.array(self.zp_ref_new)
        self.sp_logger.info(info_g+' Number of stars after sigma clipping: '+str(len(self.zp_sci)))
        # sys.exit()

        # self.sp_logger.info(self.zp_sci,np.nanmedian(self.zp_sci),np.nanstd(self.zp_sci))
        # self.sp_logger.info(self.zp_ref,np.nanmedian(self.zp_ref),np.nanstd(self.zp_ref))
        # self.sp_logger.info()
        if (self.sci_filt=='u' or self.use_sdss) and len(self.zp_sci)>=5:N=1.75
        else:N=3
        
        # self.sp_logger.info(N)
        # for i in range(len(self.zp_sci)-1):
        #     self.sp_logger.info((self.zp_sci[i],self.zp_ref[i],self.rsq_sci[i],self.rsq_ref[i],np.nanmedian(self.zp_sci),np.nanstd(self.zp_sci),abs(self.zp_sci[i]-np.nanmedian(self.zp_sci))<N*np.nanstd(self.zp_sci),abs(self.zp_ref[i]-np.nanmedian(self.zp_ref))<N*np.nanstd(self.zp_ref)))

        # sys.exit()
        while True:
            if len(self.zp_sci)==len(self.zp_ref) and np.nanstd(self.zp_sci)>1:
                self.new_zp_sci,self.new_zp_ref,self.new_rsq_sci,self.new_rsq_ref=[],[],[],[]
                for k in range(len(self.zp_sci)):
                    if abs(self.zp_sci[k]-np.nanmedian(self.zp_sci))<N*np.nanstd(self.zp_sci) or abs(self.zp_ref[k]-np.nanmedian(self.zp_ref))<N*np.nanstd(self.zp_ref):
                        self.new_zp_sci.append(self.zp_sci[k]),self.new_zp_ref.append(self.zp_ref[k]),self.new_rsq_sci.append(self.rsq_sci[k]),self.new_rsq_ref.append(self.rsq_ref[k])

                try:
                    if len(self.new_zp_sci)==len(self.zp_sci) and all(zp>10 for zp in [len(self.new_zp_sci),len(self.zp_sci)]):
                        N-=0.25
                except Exception as e:
                    if len(self.new_zp_sci)==len(self.zp_sci) and len(self.zp_sci)>10:
                        N-=0.25
                    pass
                if N<=1:break
                self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref=np.array(self.new_zp_sci),np.array(self.new_zp_ref),np.array(self.new_rsq_sci),np.array(self.new_rsq_ref)
                break
            else:
                if len(self.zp_sci[np.where(self.zp_sci<=np.nanmedian(self.zp_sci)+N*np.nanstd(self.zp_sci))])>=2:self.new_zp_sci = self.zp_sci[np.where(self.zp_sci<=np.nanmedian(self.zp_sci)+N*np.nanstd(self.zp_sci))]
                if len(self.zp_ref[np.where(self.zp_ref<=np.nanmedian(self.zp_ref)+N*np.nanstd(self.zp_ref))])>=2:self.new_zp_ref = self.zp_ref[np.where(self.zp_ref<=np.nanmedian(self.zp_ref)+N*np.nanstd(self.zp_ref))]
                
                break


        self.sp_logger.info(info_g+' Number of stars used for zeropoint calculation: '+str(len(self.zp_sci)))
        # print(self.zp_sci)
        # print(self.rsq_sci)
        # print(self.zp_ref)
        # print(self.rsq_ref)
        # print(np.nanmedian(self.zp_sci),np.nanstd(self.zp_sci),np.nanmedian(self.zp_ref),np.nanstd(self.zp_ref))
        if len(self.zp_sci)==1:self.zp_sys_err = ((1-self.rsq_sci[0])**0.5+(1-self.rsq_ref[0])**0.5)**2
        else:self.zp_sys_err = 0
            # else:self.zp_sys_err = 0.15

        # sys.exit()
        # print(self.zp_sci)
        # print(self.rsq_sci)
        self.sp_logger.info(info_g+f' Science Zeropoints from ({len(self.zp_sci)} stars)')
        self.sp_logger.info(info_g+f'   - ZP Mean: {np.nanmean(self.zp_sci):.3f}')
        self.sp_logger.info(info_g+f'   - ZP Std: {np.nanstd(self.zp_sci):.3f}')
        self.sp_logger.info(info_g+f'   - ZP Range: [{np.nanmin(self.zp_sci):.3f},{np.nanmax(self.zp_sci):.3f}]')
        self.sp_logger.info(info_g+f'   - ZP Mode: {stats.mode(self.zp_sci):.3f}')
        self.sp_logger.info(info_g+f'   - ZP Median: {np.nanmedian(self.zp_sci):.3f}')
        self.sp_logger.info(info_g+f'   - ZP 16th & 84th percentiles: {np.nanpercentile(self.zp_sci,16):.3f}, {np.nanpercentile(self.zp_sci,84):.3f}')
        self.sp_logger.info(info_g+f'   - # above threshold ({self.thresh}): {len(self.rsq_sci[self.rsq_sci>self.thresh])}')
        self.sp_logger.info(info_g+f' Science RSQ')
        self.sp_logger.info(info_g+f'   - RSQ Mean: {np.nanmean(self.rsq_sci):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ Std: {np.nanstd(self.rsq_sci):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ Range: [{np.nanmin(self.rsq_sci):.3f},{np.nanmax(self.rsq_sci):.3f}]')
        self.sp_logger.info(info_g+f'   - RSQ Mode: {stats.mode(self.rsq_sci):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ Median: {np.nanmedian(self.rsq_sci):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ 16th & 84th percentiles: {np.nanpercentile(self.rsq_sci,16):.3f}, {np.nanpercentile(self.rsq_sci,84):.3f}')
        self.sp_logger.info(info_g+f' Reference Zeropoints from ({len(self.zp_ref)} stars)')
        self.sp_logger.info(info_g+f'   - ZP Mean: {np.nanmean(self.zp_ref):.3f}')
        self.sp_logger.info(info_g+f'   - ZP Std: {np.nanstd(self.zp_ref):.3f}')
        self.sp_logger.info(info_g+f'   - ZP Range: [{np.nanmin(self.zp_ref):.3f},{np.nanmax(self.zp_ref):.3f}]')
        self.sp_logger.info(info_g+f'   - ZP Mode: {stats.mode(self.zp_ref):.3f}')
        self.sp_logger.info(info_g+f'   - ZP Median: {np.nanmedian(self.zp_ref):.3f}')
        self.sp_logger.info(info_g+f'   - ZP 16th & 84th percentiles: {np.nanpercentile(self.zp_ref,16):.3f}, {np.nanpercentile(self.zp_ref,84):.3f}')
        self.sp_logger.info(info_g+f'   - # above threshold ({self.thresh}): {len(self.rsq_ref[self.rsq_ref>self.thresh])}')
        self.sp_logger.info(info_g+f' Reference RSQ')
        self.sp_logger.info(info_g+f'   - RSQ Mean: {np.nanmean(self.rsq_ref):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ Std: {np.nanstd(self.rsq_ref):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ Range: [{np.nanmin(self.rsq_ref):.3f},{np.nanmax(self.rsq_ref):.3f}]')
        self.sp_logger.info(info_g+f'   - RSQ Mode: {stats.mode(self.rsq_ref):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ Median: {np.nanmedian(self.rsq_ref):.3f}')
        self.sp_logger.info(info_g+f'   - RSQ 16th & 84th percentiles: {np.nanpercentile(self.rsq_ref,16):.3f}, {np.nanpercentile(self.rsq_ref,84):.3f}')

        # sys.exit()
        if len(self.zp_sci)==1:
            self.sp_logger.info(info_g+' ZP Sci=%.3f std=%.3f No. stars=%i'%(np.nanmedian(self.zp_sci),self.zp_sys_err,len(self.zp_sci)))
            self.sp_logger.info(info_g+' ZP Ref=%.3f std=%.3f No. stars=%i'%(np.nanmedian(self.zp_ref),self.zp_sys_err,len(self.zp_ref)))
        else:
            self.sp_logger.info(info_g+' ZP Sci=%.3f std=%.3f No. stars=%i'%(np.nanmedian(self.zp_sci),np.nanstd(self.zp_sci),len(self.zp_sci)))
            self.sp_logger.info(info_g+' ZP Ref=%.3f std=%.3f No. stars=%i'%(np.nanmedian(self.zp_ref),np.nanstd(self.zp_ref),len(self.zp_ref)))

        # sys.exit()
        if self.zp_only!=False:
            if str(np.nanmedian(self.zp_sci))!='nan':
                self.zp_med,self.zp_std,self.zp_len = np.nanmedian(self.zp_sci),np.std(self.zp_sci),len(self.zp_sci)
                if len(self.name)>0:
                    self.zp_name = re.sub(f'data/IOO_Stands/{self.folder}_1','',self.name[0])
                    self.zp_file_name = self.path+"zeropoints/"+self.folder+f"/{re.sub('.fits','',self.zp_name)}zpts_{self.sci_filt}.txt"
                    
                else:
                    if not os.path.exists(self.path+'zeropoints/'+self.folder):
                        os.mkdir(self.path+'zeropoints/'+self.folder)
                    self.zp_name = re.sub(f'data/IOO_Stands/{self.folder}_1/','',self.name)
                    self.zp_file_name = self.path+"zeropoints/"+self.folder+f"/{re.sub('.fits','',self.zp_name)}zpts_{self.sci_filt}.txt"
                self.zp_file = open(self.zp_file_name,"w")
                self.zp_file.write(f"{self.zp_med},{self.zp_std},{self.zp_len}")
                self.zp_file.close()

            self.sys_exit=True
            return self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref
        if self.user_zp_sci!=None:
            self.sp_logger.info(info_g+' Using user defined zeropoint for science image')
            self.zp_sci,self.rsq_sci = [float(self.user_zp_sci)],[0.999999999]
        if self.user_zp_ref!=None:
            self.sp_logger.info(info_g+' Using user defined zeropoint for reference image')
            self.zp_ref,self.rsq_ref = [float(self.user_zp_ref)],[0.999999999]

        return self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref

    def check_nearby_resid(self,data,psf,sn_x,sn_y,sn_phot):
        #thid function checks if there are any nearby residuals that could be affecting the fit within 5 arcsec of sn pos
        #if there are, it will return a warning and exit the program
        #if there are not, it will return the residuals
        #it does this by measuring the photometry at 20 random locations outside the psf fit
        #if the photometry is significantly different from the photometry at the sn position, it will return a warning and exit the program
        #if the photometry is not significantly different, it will return the residuals

        if self.termoutp!='quiet':
            self.sp_logger.info(info_g+' Checking for nearby residuals')

        #create a list of 20 random positions outside the psf fit but within 5 arcsec of the sn position
        rand_x,rand_y=[],[]
        #calculate 5 arcsec in pixels
        r = 5/self.sci_ps
        for i in range(50):
            rand_x.append(random.uniform(sn_x-5,sn_x+r/2))
            rand_y.append(random.uniform(sn_y-5,sn_y+r/2))

        #calculate the photometry at the random positions
        rand_phot=[]
        # for i in range(len(rand_x)): #_fit_noshift(data_cutout=bkg_cutout,psf_array=psf)[0]
        cutouts=self.cutout_psf(data=data,psf_array=psf,xpos=rand_x[i],ypos=rand_y[i])
        rand_phot=self.psf_fit_noshift(cutouts,psf_array=psf)[:,0]

        #check the distribution of the photometry at the random positions
        rand_phot=np.array(rand_phot)
        # self.sp_logger.info(rand_phot)
        # self.sp_logger.info(np.nanmedian(rand_phot),np.nanstd(rand_phot))
        # self.sp_logger.info(np.max(rand_phot),np.min(rand_phot))
        # self.sp_logger.info()
        # self.sp_logger.info('sn_phot',sn_phot)
        # self.sp_logger.info()
        #find the difference between the photometry at the sn position and the photometry at the random positions
        diff_phot=sn_phot-rand_phot
        for j in range(len(diff_phot)):
            #find thr median difference excluding the jth value
            temp_diff_phot=np.delete(diff_phot,j)
            med_diff_phot=np.nanmedian(temp_diff_phot)
            std_diff_phot=np.nanstd(temp_diff_phot)
            mean_diff_phot=np.nanmean(temp_diff_phot)

            #if the jth value is more than 3 sigma away from the median diff, self.sp_logger.info a warning and exit the program
            if diff_phot[j]>med_diff_phot+3*std_diff_phot or diff_phot[j]<med_diff_phot-3*std_diff_phot:
                #see if the diff_phot[j] is more negative than the surrounding background
                if diff_phot[j]<np.nanmedian(rand_phot)-np.nanstd(rand_phot):

                    if self.termoutp!='quiet':
                        self.sp_logger.info('Residuals detected within 5 arcsec of sn position')
                        self.sp_logger.info('')
                




    def mag_err_function(self,data,psf,zp_sci,num,sn_x,sn_y):
        mag=99.0
        magstd=99.0
        magerr=99.0
        minimal_mag=99.0

        self.sp_logger.info(info_g+' Calculating magnitudes')
        magerr,maglim,flux_bkg_list,flux_new_sn_list=[],[],[],[]
        self.ra_off,self.dec_off = None,None
    
        # cutout a part of the image the same size as the measured PSF
        self.sn_cutout=self.cutout_psf(data=data,psf_array=psf,xpos=[sn_x],ypos=[sn_y])[0]
        # calculate the magnitude of the object
        if self.forced_phot==False: 
            # self.main_sn_psf_fit = self.psf_fit_bg([data],psf_array=psf,sn_x=sn_x,sn_y=sn_y)[0]
            self.main_sn_psf_fit = self.psf_fit([self.sn_cutout],psf_array=psf)[0]

        else: 
            if len(self.forced_phot)==1:self.forced_phot=self.forced_phot[0].split(',')
            self.sp_logger.info(info_g+f" Performing forced photometry at {self.forced_phot[0]} {self.forced_phot[1]}")
            sn_x,sn_y = self.forced_phot[0],self.forced_phot[1]
            new_sci_c = SkyCoord(ra=sn_x,dec=sn_y,unit=(u.hourangle, u.deg),frame='fk5')
            sn_x,sn_y = wcs_to_pixels(self.sci_conv_name,np.column_stack((new_sci_c.ra.deg,new_sci_c.dec.deg)))[0]

            self.sn_cutout = self.cutout_psf(data=data,psf_array=psf,xpos=[sn_x],ypos=[sn_y])[0]
            # print(len(self.forced_phot))
            self.main_sn_psf_fit = self.psf_fit_noshift([self.sn_cutout],psf_array=psf)[0]

        sn_flux= self.main_sn_psf_fit[0]
        sn_mag=-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)
        # self.sp_logger.info the offset of the shifted fit
        self.sp_logger.info(info_g+f' Preliminary:')
        self.sp_logger.info(info_g+f' MJD: {self.sci_mjd}')
        self.sp_logger.info(info_g+f' BACKGROUND: {self.main_sn_psf_fit[1]:.3f} counts')
        self.sp_logger.info(info_g+f' SN FLUX {sn_flux:.3f} counts')
        self.sp_logger.info(info_g+f' SN MAG {(-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)):.3f} mag')
        self.sp_logger.info(info_g+f' SN MAG - BACKGROUND {(-2.5*np.log10(sn_flux-np.abs(self.main_sn_psf_fit[1]))+np.nanmedian(zp_sci)):.3f} mag')
        # self.sp_logger.info('-'*20)
        if self.forced_phot==False:
            self.sp_logger.info(info_g+' Offset of shifted PSF fit (xpos,ypos)=%.3f %.3f'%(self.main_sn_psf_fit[3],self.main_sn_psf_fit[4]))
            # if the xoff_arc and yoff_arc are greater than 1.5 arcsec, then the fit is shifted too much, 
            # defaulting to the original position so will use psf_fit_noshift
            self.sp_logger.info(info_g+' xoff_arc=%.3f arcsec, yoff_arc=%.3f arcsec'%(self.main_sn_psf_fit[5],self.main_sn_psf_fit[6]))
            self.ra_off,self.dec_off = self.main_sn_psf_fit[5],self.main_sn_psf_fit[6]
            self.max_psf_offset=0.6
            if (sn_mag>15 and self.forced_phot==False) or (self.forced_phot!=False) or np.isnan(sn_mag):
                if  self.telescope not in SEDM and any(abs(offset)>self.max_psf_offset for offset in [self.main_sn_psf_fit[5], self.main_sn_psf_fit[6]]):
                    self.sp_logger.warning(warn_y+f" PSF fit is shifted too much, defaulting to original position")
                    self.sp_logger.warning(warn_y+" Magnitude before shifting to original position = %.3f"%sn_mag)

                    # self.main_sn_psf_fit = self.psf_fit_bg([data],psf_array=psf,sn_x=sn_x,sn_y=sn_y,shift=False)[0]
                    self.main_sn_psf_fit = self.psf_fit_noshift([self.sn_cutout],psf_array=psf)[0]

                    sn_flux= self.main_sn_psf_fit[0]
                    sn_mag=-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)
                    self.sp_logger.info(warn_y+" New magnitude = %.3f"%sn_mag)
            if self.telescope in SEDM and any(abs(offset)>1 for offset in [self.main_sn_psf_fit[5], self.main_sn_psf_fit[6]]):
                self.sp_logger.warning(warn_y+f" Mag before shifting to original position = %.3f"%sn_mag)
                self.sp_logger.warning(warn_y+f" PSF fit is shifted too much, defaulting to original position")
                self.sp_logger.warning(warn_y+" xoff =%.3f arsec, yoff_arc=%.3f arcsec"%(self.main_sn_psf_fit[5],self.main_sn_psf_fit[6]))
                self.main_sn_psf_fit = self.psf_fit_noshift([self.sn_cutout],psf_array=psf)[0]
                sn_flux= self.main_sn_psf_fit[0]
                sn_mag=-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)
                self.sp_logger.info(warn_y+" New magnitude = %.3f"%sn_mag)

        psf_size=np.shape(psf)[0]+1
        #chose num number of coordinates to calculate the magnitude error within the image size but outside the psf
        x=np.linspace(-num*psf_size,num*psf_size,(num*2)+1)
        len_x1 = len(x)
        # self.sp_logger.info('x min max',x.min(),x.max())
        # self.sp_logger.info('data shape',np.shape(data))

        x0, y0, radius = 0.0, 0.0, psf_size/2
        x, y = np.meshgrid(x, x)
        r = np.sqrt((x - x0)**2 + (y - y0)**2)
        outside = (r > radius)
        x=x[outside].flatten()
        y=y[outside].flatten()

        # print(len(x),len_x1)


        _inj_mask = self.valid_mask if (hasattr(self,'valid_mask') and self.valid_mask is not None) else None
        _has_padding = _inj_mask is not None and not np.all(_inj_mask)

        # Pre-build a pool of valid injection positions when padding is present.
        # This avoids the old approach of rejection-sampling (30 random tries per position)
        # which could still land on bad pixels near padding boundaries.
        _half_psf = int(np.shape(psf)[0] / 2) + 1
        _h, _w = np.shape(data)
        if _has_padding:
            from scipy.ndimage import binary_erosion
            # Erode the valid mask by the PSF half-width so every sampled position
            # guarantees a fully-valid PSF-sized cutout.
            _struct = np.ones((2*_half_psf+1, 2*_half_psf+1), dtype=bool)
            _eroded = binary_erosion(_inj_mask, structure=_struct)
            _vy, _vx = np.where(_eroded)
            # Express positions relative to SN position (matching the x/y offset grid)
            _valid_offsets = list(zip((_vx - sn_x).astype(int), (_vy - sn_y).astype(int)))
            np.random.shuffle(_valid_offsets)
            self.sp_logger.info(info_g+f' Injection pool: {len(_valid_offsets)} valid positions in valid-data region')
        else:
            _valid_offsets = None

        _pool_idx = 0
        all_cutouts = []
        for i in range(0,len(x)):
            if _has_padding and _valid_offsets is not None:
                if _pool_idx >= len(_valid_offsets):
                    break  # exhausted the pool
                x__, y__ = _valid_offsets[_pool_idx]
                _pool_idx += 1
                # Guard: skip if within PSF exclusion radius of SN
                if np.sqrt(x__**2 + y__**2) <= radius:
                    continue
            else:
                x__,y__=x[i],y[i]
                # Guard: skip out-of-bounds positions
                cx, cy = sn_x+x__, sn_y+y__
                if (cx-_half_psf < 0 or cx+_half_psf >= _w or
                        cy-_half_psf < 0 or cy+_half_psf >= _h):
                    continue

            bkg_cutout=self.cutout_psf(data=data,psf_array=psf,xpos=[sn_x+x__],ypos=[sn_y+y__])[0]
            # all_cutouts.append(bkg_cutout)

            flux_bkg=self.psf_fit_noshift(data_cutout=[bkg_cutout],psf_array=psf)[0][0]
            #flux_bkg_list is the PSF fitted to the sky
            flux_bkg_list.append(flux_bkg)

            # Cache the source PSF fit result (psf_fit is called once, not twice).
            _sn_fit = self.psf_fit([self.sn_cutout], psf_array=psf)[0]
            new_sn = bkg_cutout + (psf * _sn_fit[0]) + _sn_fit[1]
            # Use psf_fit_noshift for artificial source recovery, matching the
            # background measurement method. The injection position is exactly known,
            # so no centroid shift is needed. Using psf_fit (shift-allowed) lets
            # chi2_shift latch onto subtraction artefacts at each injection site,
            # inflating std(flux_new_sn_list) and producing unrealistically low S/N.
            flux_new_sn=self.psf_fit_noshift(data_cutout=[new_sn],psf_array=psf)[0][0]
            #flux_new_sn_list is the PSF fitted to the artificial supernova
            flux_new_sn_list.append(flux_new_sn)

        # Remove masked/non-finite values that sigma_clip wraps as '--' or NaN before statistics
        flux_bkg_list    = [f for f in flux_bkg_list    if np.isfinite(float(f)) if f != '--']
        flux_new_sn_list = [f for f in flux_new_sn_list if np.isfinite(float(f)) if f != '--']

        if len(flux_bkg_list) < 3:
            self.sp_logger.warning(warn_r+' Too few valid background injection positions (<3) — limiting magnitude unreliable')
            flux_bkg_list = [0.0]  # prevent downstream crash; limits will be flagged as 99
        if len(flux_new_sn_list) < 3:
            self.sp_logger.warning(warn_r+' Too few valid artificial SN positions (<3) — magnitude error unreliable')
            flux_new_sn_list = [0.0]

        self.sp_logger.info(info_g+f' Background injection positions used: {len(flux_bkg_list)}')
        self.sp_logger.info(info_g+' Sigma clipping the background flux')
        flux_bkg_list=sigma_clip(flux_bkg_list,sigma=2.5,maxiters=5)
        self.sp_logger.info(info_g+' Sigma clipping the artificial supernova flux')
        flux_new_sn_list=sigma_clip(flux_new_sn_list,sigma=2.5,maxiters=5)

        # Scatter-based limits: use MAD for robustness against residual-inflated noise
        # MAD is much less sensitive to outliers from subtraction artifacts, cosmic rays,
        # and bright star residuals that survive sigma clipping.
        _bkg_std  = np.nanstd(flux_bkg_list)
        _bkg_mad  = 1.4826 * np.nanmedian(np.abs(np.asarray(flux_bkg_list) - np.nanmedian(flux_bkg_list)))

        _bkg_std_original = _bkg_std
        # If std is significantly higher than MAD (>2.5x), residuals are inflating the noise
        if _bkg_mad > 0 and _bkg_std > 2.5 * _bkg_mad:
            self.sp_logger.info(warn_y+f' Background noise inflated by residuals: std={_bkg_std_original:.1f} vs MAD={_bkg_mad:.1f}')
            _bkg_std = _bkg_mad  # Use the robust MAD estimate
            self.sp_logger.info(info_g+f' Using robust MAD-based noise estimate: {_bkg_std:.1f} counts')
        else:
            _bkg_std = _bkg_std if _bkg_std > 0 else 1e-30  # guard against zero std
        _flux_1sig = 1.0 * _bkg_std
        _flux_3sig = 3.0 * _bkg_std
        _flux_5sig = 5.0 * _bkg_std
        _zp = np.nanmedian(zp_sci)
        _lim_1sig = -2.5*np.log10(_flux_1sig) + _zp
        _lim_3sig = -2.5*np.log10(_flux_3sig) + _zp
        _lim_5sig = -2.5*np.log10(_flux_5sig) + _zp

        SNR = sn_flux / _bkg_std
        self.SNR = SNR

        self.sp_logger.info(info_g+' S/N (std background)= %.3f'%SNR)
        self.sp_logger.info(info_g+' S/N (std artifical sn)= %.3f'%(sn_flux/max(np.nanstd(flux_new_sn_list),1e-30)))

        if SNR <= 2:
            #not detected
            minimal_mag=sn_mag
            mag=99.0
            magstd=99.0
            magerr=99.0

            self.sp_logger.info(info_g+' Less than 2 sigma — reporting limit')
            # Most conservative of: flux of marginal detection + 2*sigma, or 3-sigma scatter limit
            _lim_marginal = -2.5*np.log10(max(abs(sn_flux) + 2*_bkg_std, 1e-30)) + _zp
            maglim = np.nanmin([_lim_marginal, _lim_3sig])
            self.sp_logger.info(info_g+' 5-sig limit =%.3f'%_lim_5sig)
            self.sp_logger.info(info_g+' 3-sig limit =%.3f'%_lim_3sig)

        if SNR >= 2:
            minimal_mag=sn_mag
            mag=sn_mag
            magstd=-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.nanstd(flux_new_sn_list))
            magerr=magstd
            maglim=_lim_3sig
            self.sp_logger.info(info_g+f' Limiting magnitude (3-sigma scatter): {maglim:.3f} mag')
            self.sp_logger.info(info_g+' Mag = %.3f+/-%.3f ' %(sn_mag,magerr))
            self.sp_logger.info(info_g+f' 5-sig limit = {_lim_5sig:.3f}')
            self.sp_logger.info(info_g+f' 3-sig limit = {_lim_3sig:.3f}')

        # * If >3 sigma detection: report to Marshal with 1-sigma errors.
        # * If <3 sigma detection: report whichever of these two is more conservative:
        # * The flux corresponding to the (marginal) detection + 2*sigma.
        # * The flux corresponding to 3*sigma.

        sn_flux=sn_flux/3631
        sn_flux_err = abs(np.nanstd(flux_new_sn_list))/3631

        try:self.check_nearby_resid(data,psf,sn_x,sn_y,sn_mag)
        except:pass

        self.find_align_err()
        return(mag,magstd,magerr,maglim,SNR,sn_flux/max(np.nanstd(flux_new_sn_list)/3631,1e-30),
        _lim_1sig,_lim_3sig,_lim_5sig,
        minimal_mag,SNR,sn_flux,sn_flux_err)
    
    def find_align_err(self):
        self.find_astrometric_error()

        # self.sci_wcs_rms = [i for i in self.sci_img_hdu.header['COMMENT'] if 'code error' in i][0]

        self.xshifts,self.yshifts = np.random.normal(0,self.ra_error,10), np.random.normal(0,self.dec_error,10)
        # [print(self.xshifts[i],self.yshifts[i]) for i in range(100)]
        # sys.exit()
        return

        self.shifted_ref = []
        for k in range(len(self.xshifts)):
            self.shifted_ref.apend(scipy.ndimage.shift(self.ref_conv, [self.xshifts[k], self.yshifts[k]], order=3, mode='reflect', cval=0.0, prefilter=True))
        
        self.shifted_sub=(self.sci_conv)-(self.scale_factor*np.array(self.shifted_ref))
        self.shifted_cutouts = self.cutout_psf(data=self.shifted_sub,psf_array=self.comb_psf,xpos=[self.coords_sn[0][0]],ypos=[self.coords_sn[0][1]])

        if self.forced_phot:self.shifted_sn_psf_fit = self.psf_fit_noshift(self.shifted_cutouts,psf_array=self.comb_psf)
        else:self.shifted_sn_psf_fit = self.psf_fit(self.shifted_cutouts,psf_array=self.comb_psf)


        self.shifted_sn_fluxs = np.array(self.shifted_sn_psf_fit)[:,0]


    def find_astrometric_error(self):
        # This function finds the astrometric error between the science and reference images
        # by detecting stars in the reference and finding the difference in their positions in the science image
        # using sextractor

        self.sp_logger.info(info_g+ f' Finding astrometric error')
        # self.sci_sex_command = 
        os.system(sex_path + " " + self.sci_ali_name + " -c "+path+"config_files/align_sex.config -SATUR_LEVEL 50000 -BACK_TYPE MANUAL"+f" -CATALOG_NAME "+path+f"config_files/sci_asterr_{self.rand_nums_string}.cat")
        self.sci_sources = ascii.read(path+f"config_files/sci_asterr_{self.rand_nums_string}.cat")
        os.system(sex_path + " " + self.ref_ali_name + " -c "+path+"config_files/align_sex.config -SATUR_LEVEL 50000 -BACK_TYPE MANUAL"+f" -CATALOG_NAME "+path+f"config_files/ref_asterr_{self.rand_nums_string}.cat")
        self.ref_sources = ascii.read(path+f"config_files/ref_asterr_{self.rand_nums_string}.cat")
        
        self.files_to_clean.append(path+f"config_files/sci_asterr_{self.rand_nums_string}.cat"),self.files_to_clean.append(path+f"config_files/ref_asterr_{self.rand_nums_string}.cat")
        self.sci_sources,self.ref_sources = self.sci_sources[self.sci_sources['FLAGS']==0],self.ref_sources[self.ref_sources['FLAGS']==0]

        # If either image yields zero clean SExtractor detections, the
        # match_to_catalog_sky() call below crashes with
        # "catalog cannot be a scalar or length-0".  Skip the diagnostic
        # gracefully rather than aborting the whole reduction.
        if len(self.sci_sources) == 0 or len(self.ref_sources) == 0:
            self.sp_logger.warning(warn_y + f' Astrometric-error skipped: '
                                            f'sci sources={len(self.sci_sources)}, '
                                            f'ref sources={len(self.ref_sources)}')
            self.ra_error = float('nan')
            self.dec_error = float('nan')
            self.matched_sources = None
            return

        # print(self.sci_sources.columns)
        self.sci_sc = SkyCoord(ra=self.sci_sources['ALPHA_J2000'],dec=self.sci_sources['DELTA_J2000'],unit=(u.deg,u.deg))
        self.ref_sc = SkyCoord(ra=self.ref_sources['ALPHA_J2000'],dec=self.ref_sources['DELTA_J2000'],unit=(u.deg,u.deg))

        self.indx, self.d2d, self.d3d = self.sci_sc.match_to_catalog_sky(self.ref_sc)
        self.matches = np.where(self.d2d<=2./3600*u.deg)[0]


        self.matched_sources = self.sci_sources[self.matches]
        self.matched_ref_sources = self.ref_sources[self.indx[self.matches]]

        self.matched_sources['diff_RA'] = self.matched_sources['ALPHA_J2000'] - self.matched_ref_sources['ALPHA_J2000']
        self.matched_sources['diff_DEC'] = self.matched_sources['DELTA_J2000'] - self.matched_ref_sources['DELTA_J2000']
        self.matched_sources['diff_RA_arcsec'] = self.matched_sources['diff_RA']*3600
        self.matched_sources['diff_DEC_arcsec'] = self.matched_sources['diff_DEC']*3600

        self.matched_sources['diff_X'] = self.matched_sources['X_IMAGE'] - self.matched_ref_sources['X_IMAGE']
        self.matched_sources['diff_Y'] = self.matched_sources['Y_IMAGE'] - self.matched_ref_sources['Y_IMAGE']
        self.matched_sources['diff_X_arcsec'] = self.matched_sources['diff_X']*self.sci_ps
        self.matched_sources['diff_Y_arcsec'] = self.matched_sources['diff_Y']*self.sci_ps

        # print(self.matched_sources)

        self.ra_error,self.dec_error = np.sqrt(np.median(self.matched_sources['diff_X_arcsec']**2)), np.sqrt(np.median(self.matched_sources['diff_Y_arcsec']**2))
        # print(self.astrometric_error)
        self.sp_logger.info(info_g+f' RA error: {self.ra_error:.3f} arcsec ({self.ra_error/self.sci_ps:.3f} pixels) ')
        self.sp_logger.info(info_g+f' DEC error: {self.dec_error:.3f} arcsec ({self.dec_error/self.sci_ps:.3f} pixels) ')

        # self.sys_exit=True
        return


    def scaled_subtract(self):

        self.sp_logger.info(info_g+f" Subtracting scaled & convolved refrence image")


        #this is another check on the std of the zero point of the image! np.nanstd(zp_ref)/np.sqrt(len(zp_ref))
        if self.relative_flux==True:
            self.scale_factor=1
        else:
            self.scale_factor=np.power(10,(np.nanmedian(self.zp_sci)-np.nanmedian(self.zp_ref))/2.5)

        self.sp_logger.info(info_g+' Scale factor: %.3f'%self.scale_factor)

        if self.to_subtract==True:
            self.sci_scale_sub_name = self.path+'scaled_subtracted_imgs/'+self.sci_img_name[:-5]+'_scaled_subtraction.fits'


            self.sub=(self.sci_conv)-(self.scale_factor*self.ref_conv)
            self.hdu_fits_sub= fits.PrimaryHDU(self.sub)
            self.hdu_fits_sub.header=self.sci_conv_hdu.header
            self.hdu_fits_sub.writeto(self.sci_scale_sub_name,overwrite=True)

            #scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
            #sub=(sci_conv)-(scale_factor*scipy.ndimage.shift(ref_conv, [-0.5, -0.5], order=3, mode='reflect', cval=0.0, prefilter=True))

            self.coords_sn_sub=wcs_to_pixels(self.sci_conv_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))
            self.coords_sn_sub_x,self.coords_sn_sub_y = self.coords_sn_sub[0]


        if self.to_subtract==False:
            self.sci_scale_sub_name = self.path+'no_scaled_subtracted_imgs/'+self.sci_img_name[:-5]+'no_scaled_subtraction.fits'
            self.sub=(self.sci_conv)
            self.hdu_fits_sub= fits.PrimaryHDU(self.sub)
            #add a header
            self.hdu_fits_sub.header=self.sci_conv_hdu.header
            self.hdu_fits_sub.writeto(self.sci_scale_sub_name,overwrite=True)

            #scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
            #sub=(sci_conv)-(scale_factor*scipy.ndimage.shift(ref_conv, [-0.5, -0.5], order=3, mode='reflect', cval=0.0, prefilter=True))

            self.coords_sn_sub=wcs_to_pixels(self.sci_conv_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))
            self.coords_sn_sub_x,self.coords_sn_sub_y = self.coords_sn_sub[0]

        
        self.sp_logger.info(info_g+f' Saved scaled subtracted science image to: {self.sci_scale_sub_name}')
        self.sub_mean, self.sub_median, self.sub_std = sigma_clipped_stats(self.hdu_fits_sub.data)
        self.coords_sn=self.coords_sn_sub

        self.sub_img = fits.open(self.sci_scale_sub_name)[0]

        
        # self.files_to_clean.append(self.sci_scale_sub_name)
        return self.sub,self.sci_scale_sub_name,self.coords_sn_sub



    def get_photometry(self):

        self.sp_logger.info(info_g+f" Extracting photometry")
        self.t_photometry_start = time.time()

        self.mag=self.mag_err_function(data=self.sub,psf=self.comb_psf,zp_sci=self.zp_sci,num=10,sn_x=self.coords_sn[0][0],sn_y=self.coords_sn[0][1])
        # self.sp_logger.info(self.mag[2])
        # self.sp_logger.info(np.nanstd(self.zp_ref),len(self.zp_ref),np.nanstd(self.zp_ref)/len(self.zp_ref)**0.5)
        # self.sp_logger.info(np.nanstd(self.zp_sci),len(self.zp_sci),np.nanstd(self.zp_sci)/len(self.zp_sci)**0.5)
        # self.sp_logger.info()

        if self.mag[0]<25:
            self.mag_all_err=self.mag[2]**2
        else:
            self.mag_all_err=99

        # if self.cutout_tf==True:
        #   fig = plt.figure(figsize=(10,10))
        #   self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sub_img.data)
        #   plt.imshow(self.sub_img.data,cmap=self.cmap,vmin=self.vmin, vmax=self.vmax)
        #   plt.xlim(self.coords_sn_sub_x-115,self.coords_sn_sub_x+115)
        #   plt.ylim(self.coords_sn_sub_y-115,self.coords_sn_sub_y+115)
        #   #plot the position of the SN
        #   plt.scatter(self.coords_sn_sub_x,self.coords_sn_sub_y,marker='x',color='red',s=100)
        #   #plot the position of the SN + offset in pixels
        #   plt.scatter(self.coords_sn_sub_x+self.main_sn_psf_fit[3],self.coords_sn_sub_y+self.main_sn_psf_fit[4],marker='x',color='blue',s=100)
        #   plt.axis('off')
        #   plt.show()



        if self.cutout_tf!=False:
            self.fig,self.ax = plt.subplots(nrows=1, ncols=3)
            self.sci_img_ali_hdu = fits.open(self.sci_ali_name)[0]
            self.cmap = 'gray'
            self.coords_sn_sci_ali=wcs_to_pixels(self.sci_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
            self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_y = self.coords_sn_sci_ali
            self.vmin_sci,self.vmax_sci = visualization.ZScaleInterval().get_limits(self.sci_img_ali_hdu.data)

            self.ax[0].imshow(self.sci_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_sci, vmax=self.vmax_sci)
            self.ax[0].plot([self.coords_sn_sci_ali_x-30,self.coords_sn_sci_ali_x-10],[self.coords_sn_sci_ali_y,self.coords_sn_sci_ali_y],color='lime',lw=2.5),self.ax[0].plot([self.coords_sn_sci_ali_x+10,self.coords_sn_sci_ali_x+30],[self.coords_sn_sci_ali_y,self.coords_sn_sci_ali_y],color='lime',lw=2.5)
            self.ax[0].plot([self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_x],[self.coords_sn_sci_ali_y-30,self.coords_sn_sci_ali_y-10],color='lime',lw=2.5),self.ax[0].plot([self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_x],[self.coords_sn_sci_ali_y+10,self.coords_sn_sci_ali_y+30],color='lime',lw=2.5)
            self.ax[0].set_xlim(self.coords_sn_sci_ali_x-50,self.coords_sn_sci_ali_x+50)
            self.ax[0].set_ylim(self.coords_sn_sci_ali_y-50,self.coords_sn_sci_ali_y+50)
            self.ax[0].axis('off') 
            self.ax[1].axis('off') 
            self.ax[2].axis('off') 
            self.ax[0].set_title('New')

            if self.out_dir=="photometry/":       
                self.cutout_name = self.path+f'photometry_date/{self.folder}/cut_outs/'+self.sci_img_name[:-11]+"_cutout_panel"+self.img_type+".png"
            else:
                self.cutout_name = self.path+f'{self.out_dir}cut_outs/'+self.sci_img_name[:-11]+"_cutout_panel"+self.img_type+".png"

            self.sp_logger.info(info_g+f" Saving cutout panel to {self.cutout_name}")
            self.fig.savefig(self.cutout_name)

            self.ref_img_ali_hdu = fits.open(self.ref_ali_name)[0]
            self.coords_sn_ref_ali=wcs_to_pixels(self.ref_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
            self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_y = self.coords_sn_ref_ali

            self.vmin_ref,self.vmax_ref = visualization.ZScaleInterval().get_limits(self.ref_img_ali_hdu.data)
            self.ax[1].imshow(self.ref_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_ref, vmax=self.vmax_ref)
            self.ax[1].plot([self.coords_sn_ref_ali_x-30,self.coords_sn_ref_ali_x-10],[self.coords_sn_ref_ali_y,self.coords_sn_ref_ali_y],color='lime',lw=2.5),self.ax[1].plot([self.coords_sn_ref_ali_x+10,self.coords_sn_ref_ali_x+30],[self.coords_sn_ref_ali_y,self.coords_sn_ref_ali_y],color='lime',lw=2.5)
            self.ax[1].plot([self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_x],[self.coords_sn_ref_ali_y-30,self.coords_sn_ref_ali_y-10],color='lime',lw=2.5),self.ax[1].plot([self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_x],[self.coords_sn_ref_ali_y+10,self.coords_sn_ref_ali_y+30],color='lime',lw=2.5)
            self.ax[1].set_xlim(self.coords_sn_ref_ali_x-50,self.coords_sn_ref_ali_x+50)
            self.ax[1].set_ylim(self.coords_sn_ref_ali_y-50,self.coords_sn_ref_ali_y+50)
            self.ax[1].set_title('Ref')
            self.ax[1].axis('off')        

            self.fig.savefig(self.cutout_name)

            try:
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sub_img.data)
                self.ax[2].imshow(self.sub_img.data,cmap=self.cmap,vmin=self.vmin, vmax=self.vmax)
                #add crosshairs to the image
                self.ax[2].plot([self.coords_sn_sub_x-30,self.coords_sn_sub_x-10],[self.coords_sn_sub_y,self.coords_sn_sub_y],color='lime',lw=2.5),self.ax[2].plot([self.coords_sn_sub_x+10,self.coords_sn_sub_x+30],[self.coords_sn_sub_y,self.coords_sn_sub_y],color='lime',lw=2.5)
                self.ax[2].plot([self.coords_sn_sub_x,self.coords_sn_sub_x],[self.coords_sn_sub_y-30,self.coords_sn_sub_y-10],color='lime',lw=2.5),self.ax[2].plot([self.coords_sn_sub_x,self.coords_sn_sub_x],[self.coords_sn_sub_y+10,self.coords_sn_sub_y+30],color='lime',lw=2.5)
                self.ax[2].set_xlim(self.coords_sn_sub_x-50,self.coords_sn_sub_x+50)
                self.ax[2].set_ylim(self.coords_sn_sub_y-50,self.coords_sn_sub_y+50)
                self.ax[2].axis('off') 
                self.ax[2].set_title('Sub')      
                if self.out_dir=="photometry/":
                    self.sub_cutout_name = self.path+f'photometry_date/{self.folder}/cut_outs/'+self.sci_img_name[:-11]+"_cutout_sub"+self.img_type+".png"
                else:
                    self.sub_cutout_name = self.path+f'{self.out_dir}cut_outs/'+self.sci_img_name[:-11]+"_cutout_sub"+self.img_type+".png" 
                self.fig.savefig(self.cutout_name)
            except:
                pass
            plt.show()
            plt.close()

        # self.sp_logger.info(self.mag_all_err**0.5)
        # self.sp_logger.info(np.nanstd(self.zp_ref),len(self.zp_ref),(np.nanstd(self.zp_ref)/(len(self.zp_ref)**0.5)))
        self.mag_all_err+=(np.nanstd(self.zp_ref))**2
        self.sp_logger.info(info_g+f' Systematic zeropoint error: {np.nanstd(self.zp_sys_err):.3f}')
        self.mag_all_err+=self.zp_sys_err**2
        # self.sp_logger.info(self.zp_ref)
        self.sp_logger.info(info_g+f' Reference zeropoint std: {np.nanstd(self.zp_ref):.3f}')
        # self.sp_logger.info(self.mag_all_err**0.5)
        # self.sp_logger.info(np.nanstd(self.zp_sci),len(self.zp_sci),(np.nanstd(self.zp_sci)/(len(self.zp_sci)**0.5)))
        self.mag_all_err+=(np.nanstd(self.zp_sci))**2
        self.sp_logger.info(info_g+f' Science zeropoint std: {np.nanstd(self.zp_sci):.3f}')
        # self.sp_logger.info(self.zp_sci)
        # self.sp_logger.info(self.mag_all_err**0.5)
        self.mag_all_err=self.mag_all_err**0.5
        # self.sp_logger.info(self.mag_all_err)

        # self.sp_logger.info(info_g+f" Magnitude: {self.mag[0]} +/- {self.mag_all_err}")
        self.SNR = self.mag[10]
        if self.mag[0]>40:
            self.sp_logger.warning(warn_y+f" Magnitude of 99 measured, SNR= {self.SNR}")

        self.t_photometry_end = time.time()
        self.photometry_time = self.t_photometry_end - self.t_photometry_start
        

            # self.sp_logger.info("(mag_err**2+std_zp**2+std_ref**2)^(0.5) %.3f \n" %(self.mag_all_err))

        if self.relative_flux!=True:
            self.sp_logger.info(colored('--------------------------------------------------------------------------------', 'blue'))
            # self.sp_logger.info(self.mag[0],self.mag_all_err,self.mag[3])
            self.sp_logger.info(' Mag = %.3f+/-%.3f lim=%.3f MJD=%.3f' %(self.mag[0],self.mag_all_err,self.mag[3],self.sci_mjd))
            self.sp_logger.info(colored('--------------------------------------------------------------------------------', 'blue'))

        if self.relative_flux==True or self.mag[0]>40:
            self.sp_logger.info(colored('--------------------------------------------------------------------------------', 'blue'))
            self.sp_logger.info('Flux = %.3f+/-%.3f' %(self.mag[11],self.mag[12]))
            self.sp_logger.info(colored('--------------------------------------------------------------------------------', 'blue'))

            self.sp_logger.info(info_g+' Photometry time: '+"%.1f" % round(self.photometry_time, 0)+'seconds')

        self.memkb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        self.sp_logger.info(info_g+' Memory usage Mbyte: '+" %5.1f" %(self.memkb/1024.0/1024.0))


        self.t_end = time.time()



        self.img_type=""
        if self.if_stacked==True:
            self.img_type = "_stacked"
        if os.path.exists(self.path+self.out_dir):
            if self.morning_round_up==False:
                if self.to_subtract==False:
                    if self.out_dir=="photometry/":
                        self.phot_file_name = self.path+'photometry/'+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                        self.phot_date_name = self.path+f'photometry_date/{self.folder}/'+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                    elif self.out_dir!="photometry/":
                        self.phot_file_name = self.path+self.out_dir+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                        
                elif self.to_subtract!=False:
                    if self.out_dir=="photometry/":
                        self.phot_file_name = self.path+'photometry/'+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"
                        self.phot_date_name = self.path+f'photometry_date/{self.folder}/'+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"
                    elif self.out_dir!="photometry/":
                        self.phot_file_name = self.path+self.out_dir+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"

            elif self.morning_round_up!=False:
                if self.to_subtract==False:
                    if self.out_dir=="photometry/":
                        self.phot_file_name = self.path+'photometry/'+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                        self.phot_date_name = self.path+f'photometry_date/{self.folder}/morning_rup/'+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                    elif self.out_dir!="photometry/":
                        self.phot_file_name = self.path+self.out_dir+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                        self.phot_date_name = self.path+self.out_dir+f'morning_rup/'+self.sci_img_name[:-11]+self.img_type+"_no_sb_photometry.txt"
                        
                elif self.to_subtract!=False:
                    if self.out_dir=="photometry/":
                        self.phot_file_name = self.path+'photometry/'+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"
                        self.phot_date_name = self.path+f'photometry_date/{self.folder}/morning_rup/'+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"
                    elif self.out_dir!="photometry/":
                        self.phot_file_name = self.path+self.out_dir+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"
                        self.phot_date_name = self.path+self.out_dir+f'morning_rup/'+self.sci_img_name[:-11]+self.img_type+"_photometry.txt"

            self.text_file = open(self.phot_file_name, "w")
            self.text_file.write(self.sci_obj+" "+f"sdss{self.sci_filt}"+" "+str(self.sci_mjd)+" %.3f %.3f %.3f %s %s %s  %.3f %.3f \n" % (self.mag[0],self.mag_all_err,self.mag[3],self.ra_string,self.dec_string,self.sci_exp_time,self.mag[11],self.mag[12]))
            self.text_file.write(str(f'TELESCOPE: {self.telescope}')+'\n')
            self.text_file.write(f"PSF SHIFT: {self.ra_off},{self.dec_off} \n")
            self.text_file.write("S/N %.3f, %.3f sigma \n" % (self.mag[4],self.mag[5]))
            self.text_file.write("1 sigma limit: %.3f\n" % (self.mag[6]))
            self.text_file.write("3 sigma limit: %.3f\n" % (self.mag[7]))
            self.text_file.write("5 sigma limit: %.3f\n" % (self.mag[8]))
            self.text_file.write("mag=%.3f\n" %(self.mag[9]))
            self.text_file.write("time photometry: %.1f seconds\n" % (self.t_end-self.t_start))
            self.text_file.write("memory useage: %.1f Mbyte\n" % (self.memkb/1024.0/1024.0))
            self.text_file.write("zeropoint reference %.3f sd %.3f \n" % (np.median(self.zp_ref),np.nanstd(self.zp_ref)))
            self.text_file.write("zeropoint science %.3f sd %.3f \n" % (np.median(self.zp_sci),np.nanstd(self.zp_sci)))
            self.text_file.write("(mag_err**2+std_zp**2+std_ref**2)^(0.5) %.3f \n" %(self.mag_all_err))
            self.text_file.write("#object,filter,mjd,mag,mag_err,mag_lim,ra,dec,exp_t,flux,flux_err \n")
            self.text_file.close()
    
            self.sp_logger.info(info_g+f' Photometry results written to {self.phot_file_name}')
            # shutil.copy2(self.phot_file_name, self.phot_file_name.replace)
        if self.out_dir=="photometry/":
            shutil.copy2(self.phot_file_name, self.phot_date_name)
            self.sp_logger.info(info_g+f' Photometry results copied to {self.phot_date_name}')
            # sys.exit()


        self.sp_logger.info(info_g+' Total time: '+"%.1f" % round(self.t_end-self.t_start, 0)+' seconds')

        return {"obj":self.sci_obj,"filt":f"sdss{self.sci_filt}",'mjd':str(self.sci_mjd),"mag":self.mag[0],"mag_err":self.mag_all_err,"lim_mag":self.mag[3],
                "ra":self.ra_string,"dec":self.dec_string,"exp_t":self.sci_exp_time,"flux":self.mag[11],"flux_err":self.mag[12],'seeing':self.sci_seeing,'SNR':self.SNR}


    def delete_photometry(self,phot_id_data):
        self.response = api("delete", f"api/photometry/{phot_id_data['phot_id']}")
        if str(self.response.status_code)=='200':
            self.sp_logger.info(info_g+f" Successfully deleted photometry for phot ID: {phot_id_data['phot_id']}, mjd: {phot_id_data['mjd']}")
        else:
            self.sp_logger.info(warn_r+f" Unable to delete photometry for phot ID: {phot_id_data['phot_id']},mjd: {phot_id_data['mjd']} STATUS CODE={self.response.status_code}")


    def save_to_group(obj_id):
        save_group_id = '1754'
        save_url = "http://localhost:5000/api/source_groups"

        save_payload = {
            "objId": obj_id,
            "inviteGroupIds": [save_group_id],
            "unsaveGroupIds": []
        }

        save_headers = {"Content-Type": "application/json"}

        # save_response = requests.post(save_url, json=save_payload, headers=save_headers)
        save_response = requests.post(url=f"https://fritz.science/api/source_groups",headers = {'Authorization': f'token {token}'} ,json=save_payload)

        if save_response.status_code==200:self.sp_logger.info(info_g+f" Saved {obj_id} to group {save_group_id}")
        else:self.sp_logger.info(warn_r+f" Unable to save {obj_id} to group {save_group_id}, STATUS CODE={save_response.status_code}")
        return

    def sedm_uplaod_phot(self):
        #essentially the same as uplaod phot but specific to SEDM/P60
        return

    def upload_phot(self):
        try:save_to_group(self.sci_obj,1754)
        except:pass
        print(info_g+f" Attempting to upload photometry")

        # print(self.mag[0]+self.mag_all_err,self.mag[3],self.mag[0]+self.mag_all_err>self.mag[3])
        if self.SNR<=3 or self.mag[0]+self.mag_all_err>self.mag[3]:
            print(warn_r+f' S/N of {self.SNR:.2f} or mag+magerr of {self.mag[0]+self.mag_all_err:.2f} > maglim of {self.mag[3]:.2f}, uploading limit')
            self.mag=[99.0,90.0,99.0,self.mag[3]]
        

        
        # if self.sci_obj.startswith('AT') or self.sci_obj.startswith('SN'):self.sci_obj = self.sci_obj[2:]
        self.data = {"filter":f"sdss{self.sci_filt}","magerr": self.mag_all_err,"obj_id": self.sci_obj,
                    "mag":self.mag[0],"limiting_mag": self.mag[3],"mjd": self.sci_mjd,"magsys": "ab","group_ids":'all'}
        # return 
        self.data['instrument_id'],self.data['origin']='2','SEDM_SUBPHOT_KPIPE'
        # elif self.telescope in ['SLT','Lulin','slt','lulin']:self.data['instrument_id'],self.data['origin']='1104','SLT_SUBPHOT_KPIPE'
        # elif self.telescope in ['TJO','tjo','MEIA3','meia3']:self.data['instrument_id'],self.data['origin']='1103','TJO_SUBPHOT_KPIPE'
        # elif self.telescope in ['OSIRIS','osiris','GTC-OSIRIS','gtc-osiris']:self.data['instrument_id'],self.data['origin']='37','OSIRIS_SUBPHOT_KPIPE'
        # elif self.telescope in ['HIPERCAME','hipercame','GTC-HIPERCAM','gtc-hipercam']:self.data['instrument_id'],self.data['origin']='1105','HIPERCAM_SUBPHOT_KPIPE'
        self.upload_new=True
        self.data['altdata'] = {}
        if self.if_stacked==True:self.data['altdata']['stacked'],self.data['altdata']['no_in_stack']= True,self.no_stacked

        self.data['altdata']['Reducer'] = 'K-Ryan Hinds'
        # self.data['altdata']['Proposal PI'] = 'Dan Perley'
        # self.data['altdata']['exptime'] = self.sci_exp_time
        # self.data['altdata']['seeing'] = self.sci_seeing
        # if self.sci_prop!=None:
        #     self.data['altdata']['Proposal ID']=self.sci_prop
        #     if self.sci_prop in ['JL24A04','JL24B14']:self.data['altdata']['Proposal PI'] = 'K-Ryan Hinds'
        #     elif self.sci_prop in ['JL24B15','JL25A01']:self.data['altdata']['Proposal PI'] = 'Jacob Wise'
        #     elif self.sci_prop in ['JL24B10']:self.data['altdata']['Proposal PI'] = 'Chris Copperwheat'



        self.name=self.sci_obj

        if any(X in self.name for X in ['-ugriz','-griz','-gri']):
            self.name = re.sub(r'-ugriz|-griz|-gri', '', self.name)

        if self.name.startswith('SN') or self.name.startswith('AT'):
            self.name = self.name[2:]

        self.data['obj_id'] = self.name

        [print(key+': '+str(item)) for key,item in self.data.items()]
        print(info_g+f' Fritz page: https://fritz.science/source/{self.name}')

        if self.mag[0]>40:
            print(warn_r+f' Magnitude of {self.mag[0]} recorded, not uploading magnitude')
            if self.mag[3]<30:
                print(info_g+f' Limiting magnitude of {np.round(self.mag[3],3)} recorded, attempting to upload, checking if already uploaded')
                # self.upload_new=False
                # if self.upload_new==False:
                #     print(info_g+f' Limiting magnitude of {np.round(self.mag[3],3)} already uploaded, not uploading')
                #     [print(key,item) for key,item in self.data.items()]
                #     return
                
                # return 
                self.all_fritz_data = self.SN_data_phot(self.name)
                self.fritz_ind = [ind for ind,val in enumerate(self.all_fritz_data['instrument_id']) if val==2 and self.all_fritz_data['filter'].iloc[ind]==self.data['filter'] and self.all_fritz_data['origin'].iloc[ind]=='SEDM_SUBPHOT_KPIPE']

                self.on_fritz_mjd = [np.round(self.all_fritz_data['mjd'].iloc[ind],6) for ind in self.fritz_ind]
                self.on_fritz_id = [self.all_fritz_data['id'].iloc[ind] for ind in self.fritz_ind]
                self.on_fritz_mag = [np.round(self.all_fritz_data['limiting_mag'].iloc[ind],3) for ind in self.fritz_ind]
                self.upload_new=True  #identifier to update photometry in the case of stacking

                if self.sci_mjd in self.on_fritz_mjd:
                    print(info_g+f' Limiting magnitude of {np.round(self.mag[3],3)} already uploaded, checking if it is the same')
                    self.upload_new=False
                    if any(abs(self.on_fritz_mag[ind]-self.mag[3])<0.1 for ind in range(len(self.on_fritz_mag))) or any(self.on_fritz_mag[ind]==np.round(self.mag[3],3) for ind in range(len(self.on_fritz_mag))):
                        print(info_g+f' Limiting magnitude of {np.round(self.mag[3],3)} already uploaded, not uploading')
                        self.sys_exit=True 
                        return
                    else:   
                        print(info_g+f' Limiting magnitude of {np.round(self.mag[3],3)} not uploaded, uploading')
                        self.upload_new=True

                if self.upload_new==True:
                    self.data = {'filter':f"sdss{self.sci_filt}",'mag':None,'magerr':None,'mjd':self.sci_mjd,'limiting_mag_nsigma':5,'obj_id':self.sci_obj,'origin':'SEDM_SUBPHOT_KPIPE',
                                        'magsys':'ab','group_ids':'all','limiting_mag':self.mag[3],'instrument_id':'2',}
                    self.data['altdata'] = {}


                    self.data['altdata']['Reducer'] = 'K-Ryan Hinds'
                    # self.data['altdata']['Proposal PI'] = 'Dan Perley'
                    self.data['altdata']['exptime'] = self.sci_exp_time
                    self.data['altdata']['seeing'] = self.sci_seeing
                    # if self.sci_prop!=None:
                        # self.data['altdata']['Proposal ID']=self.sci_prop
                        # if self.sci_prop in ['JL24A04','JL24B14']:self.data['altdata']['Proposal PI'] = 'K-Ryan Hinds'
                        # elif self.sci_prop in ['JL24B15']:self.data['altdata']['Proposal PI'] = 'Jacob Wise'

                    # if self.if_stacked==True:self.data['altdata']['stacked'],self.data['altdata']['no_in_stack']= True,self.no_stacked
                    

                    self.start_upl_time = time.time()
                    attempts=0
                    while True:
                        self.response = requests.post(url=f"https://fritz.science/api/photometry",headers = {'Authorization': f'token {token}'} ,json=self.data)

                        if self.response.status_code==200:
                            break
                        if time.time()-self.start_upl_time>60 and attempts<=1:
                            print(warn_y+f' Upload timed out, waiting 30s and trying again')
                            time.sleep(30)
                            attempts+=1
                            self.start_upl_time = time.time()
                        elif time.time()-self.start_upl_time>60 and attempts>1:
                            print(warn_r+f' Upload timed out, giving up')
                            break
                    
            return self.data
      
        self.t = date.today()
        self.now = datetime.datetime.now()
        self.current_time = self.now.strftime("%H:%M:%S")
        self.data_mjd = np.round(self.data['mjd'],6)
        self.smallest_mjd=np.min(self.MJD)
        if len(self.MJD)==1:
            self.MJD = [self.data_mjd+0.000001,self.data_mjd,self.data_mjd-0.000001]
        else:
            for m in range(len(self.MJD)):
                self.MJD.append(self.MJD[m]+0.000001)
                self.MJD.append(self.MJD[m]-0.000001)

        self.MJD = np.sort(self.MJD)
        
        try:
            self.all_fritz_data = self.SN_data_phot(self.name)
            # print(self.all_fritz_data)
            # if len(self.all_fritz_data)==0:
            #     self.all_fritz_data = self.SN_data_phot(self.name[2:])
            self.fritz_ind = [ind for ind,val in enumerate(self.all_fritz_data['instrument_id']) if val==2 and self.all_fritz_data['filter'].iloc[ind]==self.data['filter'] and self.all_fritz_data['origin'].iloc[ind]=='SEDM_SUBPHOT_KPIPE']

            self.on_fritz_mjd = [np.round(self.all_fritz_data['mjd'].iloc[ind],6) for ind in self.fritz_ind]
            self.on_fritz_id = [self.all_fritz_data['id'].iloc[ind] for ind in self.fritz_ind]
            self.on_fritz_mag = [np.round(self.all_fritz_data['mag'].iloc[ind],3) for ind in self.fritz_ind]
            self.upload_new=True  #identifier to update photometry in the case of stacking

            if self.if_stacked==True or self.upfritz_f==True:
                try:
                    self.del_ind = np.abs(self.on_fritz_mjd-self.smallest_mjd).argmin()
                    self.del_id,self.del_mjd,self.del_mag = self.on_fritz_id[self.del_ind], self.on_fritz_mjd[self.del_ind], self.on_fritz_mag[self.del_ind]
            
                    #if there is photometry on fritz that was uploaded with an mjd thats within 0.005 days of this current measurement, we want to delete and replace with the stacked version
                    #this uses the delete photometry function and then uploads this new photometric point 
                    if abs(np.round(self.sci_mjd,6)-self.del_mjd)<0.005:  
                        deletion_dict = {'phot_id':self.del_id,'mjd':self.del_mjd}
                        print(info_g+f" Updating photometric point with phot ID: {self.del_id}")
                        self.upload_new=True
                        self.delete_photometry(deletion_dict)

                        #update the photometry downloadeds
                        self.all_fritz_data = self.SN_data_phot(self.name)
                        self.fritz_ind = [ind for ind,val in enumerate(self.all_fritz_data['instrument_id']) if val==2 and self.all_fritz_data['filter'].iloc[ind]==self.data['filter']  and self.all_fritz_data['origin'].iloc[ind]=='SEDM_SUBPHOT_KPIPE']

                        self.on_fritz_mjd = [np.round(self.all_fritz_data['mjd'].iloc[ind],6) for ind in self.fritz_ind]
                        self.on_fritz_id = [self.all_fritz_data['id'].iloc[ind] for ind in self.fritz_ind]
                        self.on_fritz_mag = [np.round(self.all_fritz_data['mag'].iloc[ind],3) for ind in self.fritz_ind]



                except:
                    print(info_g+f" No photometry to delete for {self.sci_obj} in {self.data['filter']} at {self.data['mjd']}, proceeding with upload")

            if any(x in self.on_fritz_mjd for x in self.MJD)==True and self.upload_new!=True:
                print(info_g+f" Data taken for {self.sci_obj} in {self.data['filter']} at {self.data['mjd']} already uploaded, not uploading")
                # print(0,'-------------------')
                self.upload_new=False


            elif (all(x not in self.on_fritz_mjd for x in self.MJD)==True and self.upload_new!=False)==True or (self.if_stacked==True and self.upload_new!=False and np.round(self.sci_mjd,6) not in self.on_fritz_mjd)==True:
                print(info_g+f" Data taken for {self.sci_obj} in {self.data['filter']} at {self.data['mjd']} not uploaded, proceeding with upload")
                
                #uploading data
                self.start_upl_time = time.time()
                self.response = requests.post(url=f"https://fritz.science/api/photometry",headers = {'Authorization': f'token {token}'} ,json=self.data)
                print(self.response)
                # while True:
                #     #sleep for 5 seconds and try again if the upload fails
                #     time.sleep(5)
                #     self.response = requests.post(url=f"https://fritz.science/api/photometry",headers = {'Authorization': f'token {token}'} ,json=self.data)
                #     if self.response.status_code==200:
                #         break
                #     if time.time()-self.start_upl_time>60:
                #         print(warn_r+f' Upload timed out, giving up')
                #         break
                    

                if self.response.status_code == 200:
                    print(info_g+f" Successfully uploaded photometry for {self.name} in {self.data['filter']} taken at {self.data['mjd']}, uploaded at {self.current_time}")
                else:
                    print(warn_r+f" Failed to upload photometry for {self.name} in {self.data['filter']}, status code ={self.response.status_code} ")
            
                return self.data

            else:
                print(info_b+f" Data taken for {self.sci_obj} in {self.data['filter']} at {self.data['mjd']} already uploaded, not uploading")

                return self.data  

        except Exception as e:
            print(warn_y+f" Unable to find Fritz photometry for {self.sci_obj} "+e)
            print(warn_r+f" Fritz error for {self.sci_obj} in {self.sci_filt} band, (ie. event not on Fritz or error retrieving photometry from Fritz, check token and Fritz status)")
        
    def clean_directory(self):
        print(info_g+f" Cleaning directories")
        for self.out_file in os.listdir(self.path+"out"):
            if self.rand_nums_string in self.out_file:
                try:
                    if os.path.exists(self.path+'out/'+self.out_file):os.system('rm '+self.path+'out/'+self.out_file)
                except:
                    pass

        for self.config_file in os.listdir(self.path+"temp_config_files"):
            if self.rand_nums_string in self.config_file:
                try:
                    if os.path.exists(self.path+'temp_config_files/'+self.config_file):os.system('rm '+self.path+'temp_config_files/'+self.config_file)
                except:
                    pass
        
        try:
            for self.temp_file in self.files_to_clean:
                if os.path.exists(self.temp_file):os.system('rm '+self.temp_file)  
        except:
            pass 
        

    