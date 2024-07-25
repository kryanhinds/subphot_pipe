
#!/Users/kryanhinds/miniconda/envs/subphot/bin python3
# from drizzle import drizzle,doblot
from subphot_align_quick import sextract
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
from subphot_telescopes import *
#from BeautifulSoup import BeautifulSoup python2.7
from bs4 import BeautifulSoup
from image_registration import chi2_shift
from subphot_functions import *
from subphot_credentials import *
from subphot_plot_formating import *
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
from photutils.psf import IntegratedGaussianPRF
from photutils.utils import calc_total_error
from skimage.transform import AffineTransform, warp, warp_coords
import astroalign as aa
from reproject import reproject_exact,reproject_interp
# import progressbar
from time import sleep
# bar = progressbar.ProgressBar(maxval=20, \
#     widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])


if not sys.warnoptions:
    warnings.simplefilter("ignore")

norm = ImageNormalize(stretch=SqrtStretch())


def sdss_query_image(ra_string,dec_string,filt,nx,ny,log=None):


    print(info_g+f" Querrying SDSS for reference imaging in {filt}-band",(nx,ny))
    sdss_url='https://dr12.sdss.org'
    url=sdss_url+'/fields/raDec?ra='+str(ra_string)+'&dec='+str(dec_string)
    # print(url)
    # sys.exit(1)
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
        # sys.exit(1)
    except IndexError:
        print(Exception(warn_r+f' Exiting... Not in the SDSS footprint, no {filt}-band!'))
        # self.sys_exit=True
        # if self.output_to_log!=False:
        #     self.log_error_2 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Not in the SDSS footprint '+'\n'
        #     self.night_log = open(self.log_name,"r+")
        #     if self.log_error_2 in self.night_log.readlines():
        #         pass
        #     else:
        #         self.night_log.write(self.log_error_2)
        #         self.night_log.write('\n')
        #         self.night_log.close()
        # ##sys.exit(1)
        # self.sys_exit=True
        return

    try:
        if filt=='i': 
            image_link = re.sub('irg','i',image_link)
            image_link = re.sub('.jpg','.fits.bz2',image_link)
            image_name = image_link.rsplit('/', 1)[-1]
            # sys.exit(1)
        
        # print(image_link)
        # sys.exit(1)

        r=requests.get(image_link)
        r.raise_for_status()
        if os.path.exists(image_name[:-4]):
            # if self.termoutp!='quiet':
            print(info_b+' SDSS image already downloaded',image_name[:-4], (nx,ny))
        if not os.path.exists(image_name[:-4]):
            zname=image_name
            zfile = open(path+'ref_imgs/'+image_name, 'wb')
            zfile.write(r.content)
            zfile.close()
            os.system('bzip2 -d '+path+'ref_imgs/'+image_name )
            # if self.termoutp!='quiet':
            print(info_g+' Downloading new SDSS ',str(filt),'-band..',image_name[:-4], (nx,ny))
            ref_path=path+'ref_imgs/'+image_name[:-4]
            # os.system('rm '+path+'ref_imgs/'+image_name+'.bz2')
    except requests.exceptions.HTTPError as err:
        print(warn_r+' Not in SDSS footprint! Exiting..')
        # if self.output_to_log!=False:
        #     self.night_log = open(self.log_name,"r+")
        #     self.log_error_2 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Not in the SDSS footprint '+'\n'
        #     if self.log_error_2 in self.night_log.readlines():
        #         pass
        #     else:
        #         self.night_log.write(self.log_error_2)
        #         self.night_log.write('\n')
        #         self.night_log.close()
        # if self.termoutp!='quiet':
        #     print(warn_r+f" Error encountered :{err}")
        # self.sys_exit=True
        return
    # sys.exit(1)
    return(ref_path)

def get_image(args):
    nx, ny,sci_c,sci_filt = args
    return sdss_query_image(sci_c.ra.deg+(nx*0.1), 
                                sci_c.dec.deg+(ny*0.1),
                                sci_filt, nx=nx, ny=ny)
                                
# def sdss_query_image(self,ra_string,dec_string,filt,nx,ny):
#     self.sys_exit=False
#     if self.termoutp!='quiet':
#         print(info_g+f" Querrying SDSS for reference imaging in {filt}-band",(nx,ny))
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
#         print(Exception(warn_r+f' Exiting... Not in the SDSS footprint, no {filt}-band!'))
#         self.sys_exit=True
#         if self.output_to_log!=False:
#             self.log_error_2 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Not in the SDSS footprint '+'\n'
#             self.night_log = open(self.log_name,"r+")
#             if self.log_error_2 in self.night_log.readlines():
#                 pass
#             else:
#                 self.night_log.write(self.log_error_2)
#                 self.night_log.write('\n')
#                 self.night_log.close()
#         ##sys.exit(1)
#         self.sys_exit=True
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
#             if self.termoutp!='quiet':
#                 print(info_b+' SDSS image already downloaded',image_name[:-4], (nx,ny))
#         if not os.path.exists(image_name[:-4]):
#             zname=image_name
#             zfile = open(path+'ref_imgs/'+image_name, 'wb')
#             zfile.write(r.content)
#             zfile.close()
#             os.system('bzip2 -d '+path+'ref_imgs/'+image_name )
#             if self.termoutp!='quiet':
#                 print(info_g+' Downloading new SDSS ',str(filt),'-band..',image_name[:-4], (nx,ny))
#             ref_path=path+'ref_imgs/'+image_name[:-4]
#             # os.system('rm '+path+'ref_imgs/'+image_name+'.bz2')
#     except requests.exceptions.HTTPError as err:
#         print(warn_r+' Not in SDSS footprint! Exiting..')
#         if self.output_to_log!=False:
#             self.night_log = open(self.log_name,"r+")
#             self.log_error_2 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Not in the SDSS footprint '+'\n'
#             if self.log_error_2 in self.night_log.readlines():
#                 pass
#             else:
#                 self.night_log.write(self.log_error_2)
#                 self.night_log.write('\n')
#                 self.night_log.close()
#         if self.termoutp!='quiet':
#             print(warn_r+f" Error encountered :{err}")
#         self.sys_exit=True
#         return
#     # sys.exit(1)
#     return(ref_path)

# def get_image(self,args):
#     nx, ny = args
#     return self.sdss_query_image(self.sci_c.ra.deg+(nx*self.spacing), 
#                                 self.sci_c.dec.deg+(ny*self.spacing),
#                                 self.sci_filt, nx=nx, ny=ny)

class subtracted_phot(subphot_data):
    def __init__(self,ims,args):
        self.bright_cat_mag_lim = 16
        self.path=path
        self.data1_path=data1_path
        self.ims_path = '/'.join(ims[0].split('/')[:-1])
        self.args=args
        self.ims,self.cutout_tf,self.unsubtract = ims,self.args.cutout,self.args.unsubtract
        self.out_dir,self.termoutp,self.stack_tf,self.upfritz,self.cleandirs = self.args.output, self.args.termoutp,self.args.stack,self.args.upfritz,self.args.cleandirs
        self.morning_round_up,self.zp_only,self.auto_ref,self.auto_cat = self.args.mroundup,self.args.zp_only,self.args.ref_img,self.args.ref_cat
        self.relative_flux,self.special_case = self.args.relative_flux,self.args.special_case
        self.use_sdss = self.args.use_sdss
        # print(self.args)
        self.max_psf_offset = 0.9
        self.ra_dec_pos = self.args.position
        self.redo_astrometry = self.args.redo_astrometry
        self.forced_phot = self.args.forced_phot
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

        

        if not os.path.exists(self.path+'phot_fits_info'):
            os.mkdir(path+'phot_fits_info')

        if not os.path.exists(self.path+'combined_imgs'):
            os.mkdir(path+'combined_imgs')

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



        self.t_start = time.time()
        self.to_subtract=True

        if self.unsubtract==True:
            self.to_subtract=False

            
        self.if_stacked=False
        self.MJD,self.files_to_clean = [],[]
        self.rand_nums_string = ''.join(random.choice(string.digits) for i in range(5))
        


        if len(self.ims)<1:
            print(warn_r+' Add image names without fits extension to run eg.:  python LT_subtraction.py filename1 filename2')
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
                            try_name = self.data1_path+f"{self.single_data_path}/"+'_'.join(self.name_split[0:4])+'_'+f'{j}'+'_'+'_'.join(self.name_split[5:])+".fits"
                            if os.path.exists(try_name)==True and re.sub('.fits','',try_name) not in self.ims:
                                self.ims.append(fits_name)
            except:
                pass        
    
            if len(self.ims)!=1:
                if self.termoutp !='quiet':

                    print(info_g+f' Detected {self.ims[0]} is part of a stack. Found {len(self.ims)-2} other images in the stack, proceeding with the stack of {len(self.ims)-1} images total')
                pass

            else:
                self.img_type=""
                self.name=self.ims[0]+'.fits'
                # print(info_g+' Single image to subtract:',self.name)

                if '_' in self.name: 
                    self.folder = re.sub('_1','',self.name.split('/')[-1].split('_')[2])
                else:
                    self.folder = self.DATE

                if self.data1_path not in self.name:
                    self.name=self.data1_path+self.name
                
                
                self.sci_path = self.name
                self.sci_img_hdu=fits.open(self.name)[0]
                self.telescope = self.sci_img_hdu.header['TELESCOP']

                if self.telescope in SEDM:
                    self.telescope='SEDM-P60'

                # print(header_kw)
                if self.telescope in header_kw.keys():
                    if self.termoutp!='quiet':
                        # print(info_g+' Telescope:',self.telescope) 
                        #print telescpe in bold
                        print(info_g+' Telescope: \033[1m'+self.telescope+'\033[0m')
                    self.tele_kw = header_kw[self.telescope]
                    
                    self.FILT_kw,self.OBJ_kw,self.RA_kw = self.tele_kw['filter'],self.tele_kw['object'],self.tele_kw['ra']
                    self.DEC_kw,self.AIRM_kw,self.UTS_kw = self.tele_kw['dec'],self.tele_kw['airmass'],self.tele_kw['utstart']
                    self.INST_kw,self.EXPT_kw,self.DATE_kw,self.MJD_kw = self.tele_kw['instrument'],self.tele_kw['exptime'],self.tele_kw['date'],self.tele_kw['mjd']
                    self.PIX_kw,self.GAIN_kw = self.tele_kw['pixscale'],self.tele_kw['gain']
                    self.SEE_kw = self.tele_kw['seeing']

                    
                    if self.telescope=='Liverpool Telescope':
                        self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = self.tele_kw['propid'],self.tele_kw['bkg_med'],self.tele_kw['bkg_mean']
                        self.ESEE_kw = self.tele_kw['est_seeing']
                        self.sci_filt=self.sci_img_hdu.header['FILTER1'][5].lower()
                        self.sci_inst = 'IO:O'
                        self.sci_prop = self.sci_img_hdu.header[self.PRO_kw]
                    else:
                        self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = None,None,None
                        self.ESEE_kw = None
                        self.sci_prop = None
                    
                    if self.telescope=='HCT':
                        self.sci_filt=self.sci_img_hdu.header['FILTER'] #Bessel filters
                        self.sci_inst = 'HFOSC'
                        self.sci_prop = None
                        self.folder = str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))
                        # print(self.folder)

                    if self.telescope in SEDM:
                        self.redo_astrometry = True
                        self.sci_filt=self.sci_img_hdu.header['FILTER']
                        if self.name.split('/')[-1].startswith('rc'):
                            self.sci_inst='SEDM-P60'
                            self.sci_prop=None
                            self.folder=self.sci_path.split('rc')[1][2:10]
                        else:
                            self.sci_inst='SEDM-P60'
                            self.sci_prop=None
                            self.folder=self.sci_path.split('/')[-1][0:8]


                        if self.redo_astrometry==True:
                            self.sci_img_hdu=fits.open(self.sci_path)[0]
                            # if 'ACQ' in self.sci_path.split('/')[-1]:
                            print(info_g+f" This is an ACQ image, trimming the image to remove the edges")
                            #trim the left and bottom of the image (60 pixels)
                            self.orig_sci_wcs = WCS(self.sci_img_hdu.header)
                            vmin,vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                            height,width = self.sci_img_hdu.data.shape
                            trim_size=100
                            self.sci_trimmed_img = np.zeros((height,width))*np.nan
                            self.sci_trimmed_img[trim_size:,trim_size:] = self.sci_img_hdu.data[trim_size:,trim_size:]
                            # self.sci_img_hdu.header = self.sci_trimmed_img.wcs.to_header()
                            # self.sci_img_hdu.header['NAXIS1'],self.sci_img_hdu.header['NAXIS2'], = self.sci_img_hdu.data.shape[0],self.sci_img_hdu.data.shape[1]
                            # self.sci_img_hdu.header['CRPIX1'],self.sci_img_hdu.header['CRPIX2'] = self.sci_img_hdu.header['CRPIX1']-width+60,self.sci_img_hdu.header['CRPIX2']-height+60
                            # self.sci_img_hdu.header['CRVAL1'],self.sci_img_hdu.header['CRVAL2'] = self.sci_img_hdu.header['CRVAL1'],self.sci_img_hdu.header['CRVAL2']
                            # self.sci_img_hdu.header['CD_1'],self.sci_img_hdu.header['CD_2'] = self.sci_img_hdu.header['CD_1'],self.sci_img_hdu.header['CD_2']
                            # self.sci_img_hdu.header['CD2_1'],self.sci_img_hdu.header['CD2_2'] = self.sci_img_hdu.header['CD2_1'],self.sci_img_hdu.header['CD2_2']
                            # for kw in ['NAXIS1','NAXIS2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2']:
                                # self.sci_img_hdu.header[kw] = self.sci_img_hdu.header[kw
                            #show the original image next to the trimmed image
                            # fig,ax = plt.subplots(1,2)
                            # ax[0].imshow(self.sci_img_hdu.data,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)
                            # ax[1].imshow(self.sci_trimmed_img.data,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)
                            # plt.show()
                            # save the image as a new file with name *_trimmmed.fits
                            fits.writeto(self.sci_path.replace('.fits','_trimmed.fits'),self.sci_trimmed_img.data,overwrite=True,header=self.sci_img_hdu.header)
                            self.sci_path = self.sci_path.replace('.fits','_trimmed.fits')
                            # sys.exit()

                            self.sci_img_hdu.data = self.sci_trimmed_img.data
                            print(info_g+' Updating astrometry for',self.sci_path,' using local version of astrometry.net')
                            self.im_size1,self.im_size2 = self.sci_img_hdu.header['NAXIS1'],self.sci_img_hdu.header['NAXIS2']
                            self.sci_ra,self.sci_dec = self.sci_img_hdu.header[self.RA_kw],self.sci_img_hdu.header[self.DEC_kw]
                            self.sci_c=SkyCoord(self.sci_ra,self.sci_dec,unit=(u.hourangle, u.deg),frame='fk5')
                            self.sci_ra_d,self.sci_dec_d = self.sci_c.ra.deg,self.sci_c.dec.deg
                            self.ali_center_ra,self.ali_center_dec = self.sci_ra_d,self.sci_dec_d
                            self.sci_img_center = self.sci_img_hdu.data.shape[0]/2,self.sci_img_hdu.data.shape[1]/2
                            self.sci_wcs = WCS(self.sci_img_hdu.header)
                            self.sci_cnt = self.sci_wcs.all_pix2world(self.sci_img_center[0],self.sci_img_center[1],1)
                            self.sci_img_center_ra,self.sci_img_center_dec = self.sci_cnt[0],self.sci_cnt[1]
                            try:
                                self.sci_path,_ = self.update_astrometry(self.sci_path,ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])
                                print(info_g+' Updated astrometry for',self.sci_path)
                            except Exception as e:
                                # print(e)
                                print(warn_r+' Astrometry.net failed, please check the image and that the image is on target')
                            # swarp_command=swarp_path+" "+self.path+self.sci_path+" -c "+self.path+"config_files/config_resize.swarp -CENTER '"+str(self.sci_img_center_ra)+" "+str(self.sci_img_center_dec)+"' -SUBTRACT_BACK Y -VERBOSE_TYPE QUIET -RESAMPLE Y -IMAGEOUT_NAME "+self.path+self.sci_path.replace('.fits','_large.fits')+" -COMBINE N "#-IMAGE_SIZE '"+str(500)+","+str(500)+"'" #original
                            # print(swarp_command)
                            # print(self.sci_path.replace('.fits','_large.fits'))
                            # print(self.sci_path.split('/')[-1])
                            # os.system(swarp_command)
                            # sys.exit()
                            # # new_sci_img = np.pad
                            # #pad around the science image to create a 2000,2000 image, fill with zeros, place the science image in the center
                            # pad_size1 = 1500
                            # pad_size2 = 1500

                            # # Calculate pad widths
                            # pad_width1 = (pad_size1 - self.im_size1) // 2
                            # pad_width2 = (pad_size2 - self.im_size2) // 2

                            # # Pad the image  
                            # new_sci_img = np.pad(self.sci_img_hdu.data, 
                            #                     ((pad_width1, pad_size1 - self.im_size1 - pad_width1),
                            #                     (pad_width2, pad_size2 - self.im_size2 - pad_width2)),
                            #                     mode='constant',
                            #                     constant_values=1)
                            # # new_sci_img = np.pad(self.sci_img_hdu.data,((0,2000-self.im_size1),(0,2000-self.im_size2)),mode='constant',constant_values=0)
                            # print(new_sci_img.shape)
                            # #save the new science image 
                            # #adjust the header to match the new image
                            # # self.sci_img_hdu.header['NAXIS1'] = new_sci_img.shape[0]
                            # # self.sci_img_hdu.header['NAXIS2'] = new_sci_img.shape[1]
                            # # self.sci_img_hdu.header['CRPIX1'] = new_sci_img.shape[0]/2
                            # # self.sci_img_hdu.header['CRPIX2'] = new_sci_img.shape[1]/2
                            # # self.sci_img_hdu.header['CRVAL1'] = self.sci_img_center_ra
                            # # self.sci_img_hdu.header['CRVAL2'] = self.sci_img_center_dec


                            # new_sci_hdu = fits.PrimaryHDU(new_sci_img,header=self.sci_img_hdu.header)
                            
                            # new_sci_hdu.writeto(self.sci_path.replace('.fits','_large.fits'),overwrite=True)

                            # vmin,vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                            # plt.imshow(new_sci_img,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)
                            # plt.show()
                            # self.sci_path,_ = self.update_astrometry(self.sci_path.replace('.fits','_large.fits'),ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])

                            # self.sci_path = self.sci_path
                            # self.sci_path = self.name
                            self.sci_img_hdu=fits.open(self.sci_path)[0]
                            # print(self.sci_path)
                            self.redo_astrometry=False
                            # sys.exit()
                    
                    





                        # sys.exit()
                        #check if ONTARGET in header is True
                        if self.sci_img_hdu.header['ONTARGET']!=True:
                            print(self.sci_img_hdu.header['ONTARGET'])
                            print(warn_r+' This image is not on target, please check images, passing on this image')
                            self.sys_exit=True
                            return
                        self.sci_inst='SEDM-P60'
                        self.sci_prop=None
                        if 'rc' in self.sci_path.split("/")[-1]:
                            self.folder=self.sci_path.split('rc')[1][2:10]
                        else:
                            self.folder=self.sci_path.split('/')[-1].split('_')[2]

                        # print(self.folder)

                    if self.telescope=='SLT':
                        self.sci_filt=self.sci_img_hdu.header['FILTER'][0]
                        self.sci_inst = 'SLT-Andor'
                        self.sci_prop = None
                        self.folder = str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))


    
                else:
                    print(warn_r+' This script is not designed for',self.telescope)
                    print(warn_r+' Please add the header keywords for',self.telescope,'to the header_kw dictionary in the script')
                    self.sys_exit=True


                if self.redo_astrometry==True:
                    print(info_g+' Updating astrometry for',self.sci_path,' using local version of astrometry.net')
                    self.sci_path,_ = self.update_astrometry(self.sci_path,ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])
                    # self.sci_path = self.sci_path
                    # self.sci_path = self.name
                    self.sci_img_hdu=fits.open(self.sci_path)[0]
                    print(info_g+' Updated astrometry for',self.sci_path)





                # sys.exit()
                #check if ONTARGET in header is True
                if self.sci_inst in SEDM and self.sci_img_hdu.header['ONTARGET']!=True:
                    print(warn_r+' This image is not on target, please check images, passing on this image')
                    self.sys_exit=True
                    return
                self.sci_inst='SEDM-P60'
                self.sci_prop=None
                self.sci_obj=self.sci_img_hdu.header[self.OBJ_kw]

                if 'rc' in self.sci_path.split("/")[-1]:
                    self.folder=self.sci_path.split('rc')[1].split('_')[0]
                else:
                    self.folder=self.sci_path.split('/')[-1].split('_')[2]

                
                if self.termoutp!='quiet':
                    print(info_g+f' Single filter:','\033[1m'+self.sci_filt+'\033[0m')

                try:
                    self.crval1,self.crval2,self.crpix1,self.crpix2 = self.sci_img_hdu.header['CRVAL1'],self.sci_img_hdu.header['CRVAL2'],self.sci_img_hdu.header['CRPIX1'],self.sci_img_hdu.header['CRPIX2']
                    self.crdelt1,self.crdelt2,self.naxis1,self.naxis2 = self.sci_img_hdu.header['CDELT1'],self.sci_img_hdu.header['CDELT2'],self.sci_img_hdu.header['NAXIS1'],self.sci_img_hdu.header['NAXIS2']
                except:
                    pass
                if self.ra_dec_pos=='header':
                    self.sci_ra=self.sci_img_hdu.header[self.RA_kw]
                    self.sci_dec=self.sci_img_hdu.header[self.DEC_kw]
                    # print(self.sci_ra,self.sci_dec)
                else:
                    print(info_g+f' RA & Dec specified by user (J2000):','\033[1m'+self.ra_dec_pos[0]+'\033[0m','\033[1m'+self.ra_dec_pos[1]+'\033[0m')
                    print(info_g+f' Catlog RA & Dec from header (J2000):','\033[1m'+self.sci_img_hdu.header[self.RA_kw]+'\033[0m','\033[1m'+self.sci_img_hdu.header[self.DEC_kw]+'\033[0m')
                    self.sci_ra,self.sci_dec = self.ra_dec_pos[0],self.ra_dec_pos[1]
                    self.sci_c= SkyCoord(self.sci_ra,self.sci_dec, unit=(u.hourangle, u.deg),frame='fk5')
                    self.sci_ra_d,self.sci_dec_d = self.sci_c.ra.deg,self.sci_c.dec.deg
                    # print(info_g+' RA & Dec (deg):','\033[1m'+str(round(self.sci_ra_d,2))+'\033[0m','\033[1m'+str(round(self.sci_dec_d,2))+'\033[0m')
                    print(info_g+f' Checking if RA & Dec are within the image')
                    self.X_pix_co = (self.sci_ra_d - self.crval1) / self.crdelt1 + self.crpix1
                    self.X_in = 0<self.X_pix_co<self.naxis1

                    self.Y_pix_co = (self.sci_dec_d - self.crval2) / self.crdelt2 + self.crpix2
                    self.Y_in = 0<self.Y_pix_co<self.naxis2

                    if any([self.X_in,self.Y_in])==False:
                        print(warn_r+f' RA or Dec are outside the image')
                        if self.X_in==False:
                            print(warn_r+f' RA is outside the image')
                        if self.Y_in==False:
                            print(warn_r+f' Dec is outside the image')
                        print(warn_r+f" Exiting ")
                        self.sys_exit=True
                        return
                    else:
                        print(info_g+f' RA & Dec are within the image')


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


                self.sci_seeing = self.sci_img_hdu.header[self.SEE_kw]
                if self.GAIN_kw=='-': self.sci_gain=1.088
                else:self.sci_gain = float(self.sci_img_hdu.header[self.GAIN_kw])

                if self.ESEE_kw != None:
                    self.sci_est_seeing = self.sci_img_hdu.header[self.ESEE_kw]
                    self.sci_bkg_med,self.sci_bkg_mean = self.sci_img_hdu.header[self.BKGMED_kw],self.sci_img_hdu.header[self.BKGMEA_kw]
                    self.sci_prop = self.sci_img_hdu.header[self.PRO_kw]
                else:
                    self.sci_est_seeing = None
                    self.sci_bkg_med,self.sci_bkg_mean = None,None
                    self.sci_prop = None

            
                self.sci_obj, self.sep, self.tail = self.sci_obj.partition('_') #tail should be the request id from marshal triggering if there is one
                if ' ' in self.sci_obj:
                    self.sci_obj=self.sci_obj.split(' ')[0]
                if self.sci_obj=='2023vyl':self.sci_obj='ZTF23abnprwj'
                
                self.sci_img_name=self.sci_obj+'_'+self.sci_filt+'comb.fits'
                self.sci_exp_time=self.sci_img_hdu.header[self.EXPT_kw]

                if ':' in str(self.sci_ra):
                    # print(info_g+' RA and DEC not in decimal degrees')
                    #RA & DEC already noy decimal degrees

                    if all(s not in str(self.sci_dec) for s in ['+','-']):
                        # print(warn_y+' No sign on declination, assuming to be positive')
                        self.sci_dec = '+'+self.sci_dec
                    self.sci_c= SkyCoord(self.sci_ra,self.sci_dec, unit=(u.hourangle, u.deg),frame='fk5')
                else:
                    self.sci_c= SkyCoord(self.sci_ra+' '+self.sci_dec, unit=(u.deg, u.deg),frame='fk5')

                if self.PIX_kw == 'in_comments':
                    self.scale_comments = [s for s in self.sci_img_hdu.header["COMMENT"] if 'scale:' in s]
                    self.scale_comment = float(self.scale_comments[0].split(' ')[1])
                    self.sci_ps = self.scale_comment
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
                
                # print(self.sci_obj)
                # print(self.sci_filt)
                # print(self.sci_img_hdu.header[self.DATE_kw])
                # sys.exit(1)
                self.sci_img_name=self.sci_obj+'_'+self.sci_filt+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'.fits'
                self.sci_img=self.sci_img_hdu.data

        if len(self.ims)>=2:
            print(info_g+' Multiple images to stack')
            if self.ims[0].startswith('h_')==False:    #checking whether name is of a single objet or the name starts with full path
                print(info_g+' Assuming all images to be stacked are in the same directory')
                self.arg_split = self.ims[0].split('/')
                for s in range(len(self.arg_split)):
                    if 'h_' in self.arg_split[s]:
                        self.folder = self.arg_split[s].split('_')[2]
                        # self.stack_data_path = '/'.join(self.arg_split[:s]) #prefix to the first fits file showing the directory that contains all the fits to be stacked
                    # else:
            self.stack_data_path = '/'.join(self.ims[0].split('/')[:-1])
            # print(self.stack_data_path)

            self.img_type="_stacked"
            self.name = [self.ims[0]+".fits"]
            for s in range(1,len(self.ims)):
                if not self.ims[s].startswith(self.stack_data_path): self.ims[s]=self.stack_data_path+self.ims[s]
                self.name.append(str(self.ims[s])+".fits")

            for k in range(len(self.name)):
                if not self.name[k].startswith(self.data1_path): self.name[k]=self.data1_path+self.name[k]
            self.sci_img_hdu = fits.open(self.name[0])[0]
            self.telescope = self.sci_img_hdu.header['TELESCOP']

            
            if self.telescope in header_kw.keys():
                if self.termoutp!='quiet':
                    print(info_g+' Telescope:',self.telescope)
                self.tele_kw = header_kw[self.telescope]
                
                self.FILT_kw,self.OBJ_kw,self.RA_kw = self.tele_kw['filter'],self.tele_kw['object'],self.tele_kw['ra']
                self.DEC_kw,self.AIRM_kw,self.UTS_kw = self.tele_kw['dec'],self.tele_kw['airmass'],self.tele_kw['utstart']
                self.INST_kw,self.EXPT_kw,self.DATE_kw,self.MJD_kw = self.tele_kw['instrument'],self.tele_kw['exptime'],self.tele_kw['date'],self.tele_kw['mjd']
                self.PIX_kw,self.GAIN_kw = self.tele_kw['pixscale'],self.tele_kw['gain']
                self.SEE_kw = self.tele_kw['seeing']

                self.sci_date_obs = self.sci_img_hdu.header[self.DATE_kw]
                if self.UTS_kw == '-':
                    self.sci_utstart = self.sci_img_hdu.header[self.DATE_kw].split('T')[1]
                else:
                    self.sci_utstart = self.sci_img_hdu.header[self.UTS_kw]
                
                if self.telescope=='Liverpool Telescope':
                    self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = self.tele_kw['propid'],self.tele_kw['bkg_med'],self.tele_kw['bkg_mean']
                    self.ESEE_kw = self.tele_kw['est_seeing']
                    self.sci_filt=self.sci_img_hdu.header['FILTER1'][5].lower()
                    self.sci_inst = 'IO:O'
                else:
                    self.PRO_kw,self.BKGMEA_kw,self.BKGMED_kw = None,None,None
                    self.ESEE_kw = None
                
                if self.telescope=='HCT':
                    self.sci_filt=self.sci_img_hdu.header['FILTER'].lower()
                    self.sci_inst = 'HFOSC'

                if self.telescope=='SLT':
                    self.sci_filt=self.sci_img_hdu.header['FILTER'][0]
                    self.sci_inst = 'SLT-Andor'
                    self.sci_prop = None
                    self.folder = str(re.sub('-','',self.sci_img_hdu.header['DATE-OBS'].split('T')[0]))


            else:
                print(warn_r+' This script is not designed for',self.telescope)
                print(warn_r+' Please add the header keywords for',self.telescope,'to the header_kw dictionary in the script')
                self.sys_exit=True

            self.MJD = [np.round(fits.open(self.name[n])[0].header[self.MJD_kw],6) for n in range(len(self.name))]

            #this code will convert the self.MJD to self.MJD if it is in JD format
            if self.MJD[0]>400000:
                self.MJD = [self.MJD[n]-2400000.5 for n in range(len(self.MJD))]

            self.see_all = [fits.open(self.name[n])[0].header[self.SEE_kw] for n in range(len(self.name))]

            self.see_std,self.see_med = np.std(self.see_all),np.median(self.see_all)

            self.names_final = []
            for n in range(len(self.see_all)):
                if self.see_all[n]<self.see_med+5 * self.see_std:
                    self.names_final.append(self.name[n])
            
            self.name = self.names_final


            if len(self.name)>0:
                self.no_stacked = len(self.name)
                self.if_stacked=True
                self.name_joined=' '.join(self.name)
                self.sci_img_hdu=fits.open(self.name[0])[0]
                self.sci_obj=self.sci_img_hdu.header[self.OBJ_kw]
                if self.sci_obj=='2023vyl':self.sci_obj='ZTF23abnprwj'


                if self.ra_dec_pos=='header':
                    print(info_g+' RA & Dec specified by header')
                    self.sci_ra=self.sci_img_hdu.header[self.RA_kw]
                    self.sci_dec=self.sci_img_hdu.header[self.DEC_kw]
                else:
                    print(info_g+' RA & Dec specified by user:','\033[1m'+self.ra_dec_pos[0]+'\033[0m','\033[1m'+self.ra_dec_pos[1]+'\033[0m')
                    print(info_g+' Catlog RA & Dec from header (J2000):','\033[1m'+self.sci_img_hdu.header[self.RA_kw]+'\033[0m','\033[1m'+self.sci_img_hdu.header[self.DEC_kw]+'\033[0m')
                    self.sci_ra,self.sci_dec = self.ra_dec_pos[0],self.ra_dec_pos[1]
                    self.sci_c= SkyCoord(self.sci_ra,self.sci_dec, unit=(u.hourangle, u.deg),frame='fk5')
                    self.sci_ra_d,self.sci_dec_d = self.sci_c.ra.deg,self.sci_c.dec.deg
                    # print(info_g+' RA & Dec (deg):','\033[1m'+str(self.sci_ra_d)+'\033[0m','\033[1m'+str(self.sci_dec_d)+'\033[0m')
                    print(info_g+' Checking if RA & Dec are within the image')
                    self.X_pix_co = (self.sci_ra_d - self.crval1) / self.crdelt1 + self.crpix1
                    self.X_in = 0<self.X_pix_co<self.naxis1

                    self.Y_pix_co = (self.sci_dec_d - self.crval2) / self.crdelt2 + self.crpix2
                    self.Y_in = 0<self.Y_pix_co<self.naxis2

                    if any([self.X_in,self.Y_in])==False:
                        print(warn_r+' RA or Dec are outside the image')
                        if self.X_in==False:
                            print(warn_r+' RA is outside the image')
                        if self.Y_in==False:
                            print(warn_r+' Dec is outside the image')
                        print(warn_r+" Exiting ")
                        self.sys_exit=True
                        return

            
                self.sci_utstart = self.sci_img_hdu.header[self.UTS_kw]
                self.sci_obj, self.sep, self.tail = self.sci_obj.partition('_')
                if self.GAIN_kw=='-': self.sci_gain=1.088
                else:self.sci_gain = float(self.sci_img_hdu.header[self.GAIN_kw])
                

                if self.ESEE_kw != None:
                    self.sci_est_seeing = self.sci_img_hdu.header[self.ESEE_kw]
                    self.sci_bkg_med,self.sci_bkg_mean = self.sci_img_hdu.header[self.BKGMED_kw],self.sci_img_hdu.header[self.BKGMEA_kw]
                    self.sci_prop = self.sci_img_hdu.header[self.PRO_kw]
                    self.sci_airmass = self.sci_img_hdu.header[self.AIRM_kw]

                else:
                    self.sci_est_seeing = None
                    self.sci_bkg_med,self.sci_bkg_mean = None,None
                    self.sci_prop = None
                    self.sci_airmass = None
                
                self.sci_seeing=self.sci_img_hdu.header[self.SEE_kw]
                self.sci_exp_time=self.sci_img_hdu.header[self.EXPT_kw]*(len(self.ims)-1)
                if ':' in str(self.sci_ra):
                    if all(s not in str(self.sci_dec) for s in ['+','-']):
                        print(warn_y+' No sign on declination, assuming to be positive')
                        self.sci_dec = re.sub(' ','','+'+self.sci_dec)
                    self.sci_c= SkyCoord(ra=self.sci_ra,dec=self.sci_dec, frame='fk5', unit=(u.hourangle, u.deg))
                else:
                    self.sci_c= SkyCoord(self.sci_ra+' '+self.sci_dec, unit=(u.deg, u.deg),frame=FK5)

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

                self.sci_img_name=self.sci_obj+'_'+self.sci_filt+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]),seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'.fits'
                
                # plt.figure(figsize=(10,10))
                # vmin,vmax = visualization.ZScaleInterval().get_limits(stacked_data_new)
                # plt.imshow(stacked_data_new,origin='lower',cmap='gray',vmin=vmin,vmax=vmax)
                # plt.show()
                # sys.exit()

            
                if not os.path.exists(path+'combined_imgs/'+self.sci_img_name):
                    if self.termoutp!='quiet':
                        print(info_g+' Combining',len(self.ims),self.sci_filt,'images...')
                        print(info_g+f" Stacking {self.sci_obj} {self.sci_filt} images with Swarp")

                    self.ra_string="%.6f" % round(self.sci_c.ra.deg, 6)
                    if self.sci_c.dec.deg>0:
                        self.dec_string="+"+"%.6f" % round(self.sci_c.dec.deg, 6)
                    else:
                        self.dec_string='-'+"%.6f" % round(self.sci_c.dec.deg, 6)

                    # swarp_command=swarp_path+" "+self.name_joined+" -c "+self.path+"config_files/config_comb.swarp -COPY_KEYWORDS DATE-OBS -CENTER '"+self.ra_string+" "+self.dec_string+"' -SUBTRACT_BACK N -VERBOSE_TYPE NORMAL -IMAGEOUT_NAME "+path+"combined_imgs"+'/'+self.sci_img_name+" -RESAMPLE Y -RESAMPLE_DIR '"+path+"' -COMBINE Y -IMAGE_SIZE '"+str(4000)+","+str(4000)+"'"
                    # print(swarp_command)
                    # print(self.name)
                    self.stacked_path, _ = self.stack_images(self.name,save_name="combined_imgs"+'/'+self.sci_img_name)
                    # os.system(swarp_command)
                    # print('stacked')
                    # sys.exit()
                    
                    # self.status=os.system(swarp_command)
                    if self.termoutp!='quiet':
                        print(info_g+' Combined',len(self.ims),' ',self.sci_filt,' images!')
                    self.sci_img_hdu=fits.open(self.path+'combined_imgs/'+self.sci_img_name)[0]
                    self.sci_img=self.sci_img_hdu.data
                    self.sci_path = self.path+'combined_imgs/'+self.sci_img_name
                    self.name=self.sci_path
        
                elif os.path.exists(self.path+'combined_imgs/'+self.sci_img_name):
                    if self.termoutp!='quiet':
                        print(info_b+' Images:',self.sci_filt,'already combined')
                    self.sci_img_hdu=fits.open(self.path+'combined_imgs/'+self.sci_img_name)[0]
                    self.sci_img=self.sci_img_hdu.data
                    self.sci_path = self.path+'combined_imgs/'+self.sci_img_name
                    self.name=self.sci_path

            else:
                self.sys_exit=True
                pass
    
        if self.args.make_log!=False:
            self.log_dest=self.args.make_log
            if not os.path.exists(self.path+self.log_dest):
                os.mkdir(self.path+self.log_dest)

            self.output_to_log = True
            if self.log_dest=='night_log':
                self.log_name = f"{self.path}night_log/{self.DATE}_night_log.txt"
            else:
                self.log_name = f"{self.path}{self.log_dest}/{self.out_dir}_{self.sci_filt}_log.txt"

            if not os.path.exists(self.log_name):
                self.night_log = open(self.log_name,"a")
                self.night_log.close()
            else:
                self.night_log = open(self.log_name,"r+")
                self.night_log.close()

        else:
            self.output_to_log = False
        
        
        if self.sci_bkg_med!=None and self.sci_bkg_med<self.sci_exp_time/100 and self.sci_bkg_mean<self.sci_exp_time/100 and self.sci_filt!='u' and self.sys_exit!=True:  
            print(warn_r+f' Median of background counts is {self.sci_bkg_med}, possible shutter issue')
            print(warn_r+' Exiting!!')
            self.sys_exit=True

            if self.output_to_log!=False:
                self.night_log = open(self.log_name,"r+")
                self.log_error_0 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+f': Background counts <threshold, median background {self.sci_bkg_med}  '+'\n'
                if self.log_error_0 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_0)
                    self.night_log.close()

        if self.termoutp!='quiet' and self.sys_exit!=True:
            if all(x!=None for x in [self.sci_seeing,self.sci_est_seeing]):
                print(info_g+' Seeing (observed,estimated):','\033[1m'+"%.2f %.2f" %(self.sci_seeing,self.sci_est_seeing),'\033[0m')
            elif self.sci_seeing!=None and self.sci_est_seeing==None:
                print(info_g+' Seeing (observed,estimated):','\033[1m'+"%.2f %.2f" %(self.sci_seeing,0),'\033[0m')

        if self.sys_exit!=True:
            # print(self.sci_seeing,self.sci_est_seeing)
            if self.sci_seeing!=None and self.sci_est_seeing!=None:
                # print((self.sci_seeing>=self.sci_est_seeing*20 and self.sci_seeing>5.0),self.sci_seeing>7.5)
                if (self.sci_seeing>=self.sci_est_seeing*20 and self.sci_seeing>5.0)==True or self.sci_seeing>7.5:
                    if self.termoutp!='quiet':
                        print(warn_r+' Poor seeing, exiting')
                
                    print(warn_r+f" Seeing error: "+self.sci_obj+","+self.sci_filt+" band, Actual seeing:"+str(self.sci_seeing)+", Estimated seeing:"+str(self.sci_est_seeing)+" (ie. guiding/focus error or seeing got much worse)\n")
                    if self.output_to_log!=False:
                        self.log_error_1 = f"Seeing error: "+self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+" band, Actual seeing:"+str(self.sci_seeing)+", Estimated seeing:"+str(self.sci_est_seeing)+" (ie. guiding/focus error or seeing got much worse)\n"
                        self.night_log = open(self.log_name,"r+")
                        if self.log_error_1 in self.night_log.readlines():
                            pass
                        else: 
                            self.night_log.write(self.log_error_1)
                            self.night_log.close()
                    ##sys.exit(1)
                    self.sys_exit=True
                    # return


        if self.out_dir=='by_obs_date':
            self.out_dir='photometry/'
            if not os.path.exists(self.path+self.out_dir):
                print(info_g+' Creating photometry directory: ',self.path+self.out_dir)
                os.mkdir(self.path+self.out_dir)
            
            #folder is for easy way to display photometry by date of observation
            if os.path.exists(self.path+'photometry_date')==False:
                print(info_g+' Creating photometry_date directory: ',self.path+'photometry_date')
                os.mkdir(self.path+'photometry_date')
            if not os.path.exists(self.path+f'photometry_date/{self.folder}'):
                print(info_g+' Creating photometry_date directory: ',self.path+f'photometry_date/{self.folder}')
                os.mkdir(self.path+f'photometry_date/{self.folder}')  

            if self.morning_round_up!=False:
                if not os.path.exists(self.path+f'photometry_date/{self.folder}/morning_rup'):
                    print(info_g+' Creating photometry_date directory: ',self.path+f'photometry_date/{self.folder}/morning_rup')
                    os.mkdir(self.path+f'photometry_date/{self.folder}/morning_rup')

            if self.cutout_tf!=False:
                if not os.path.exists(self.path+f'photometry_date/{self.folder}/cut_outs'):
                    print(info_g+' Creating photometry_date directory: ',self.path+f'photometry_date/{self.folder}/cut_outs')
                    os.mkdir(self.path+f'photometry_date/{self.folder}/cut_outs')
        else:
            self.out_dir=str(self.out_dir)+'/'
            if not os.path.exists(self.path+self.out_dir):
                print(info_g+' Creating photometry directory: ',self.path+self.out_dir)
                os.mkdir(self.path+self.out_dir)

            if self.morning_round_up!=False:
                if not os.path.exists(self.path+self.out_dir+'morning_rup'):
                    print(info_g+' Creating photometry directory: ',self.path+self.out_dir+'morning_rup')
                    os.mkdir(self.path+self.out_dir+'morning_rup')
            
            if self.cutout_tf!=False:
                if not os.path.exists(self.path+f'{self.out_dir}cut_outs'):
                    print(info_g+' Creating photometry directory: ',self.path+f'{self.out_dir}cut_outs')
                    os.mkdir(self.path+f'{self.out_dir}cut_outs')

        
        self.sci_c=SkyCoord(self.sci_ra,self.sci_dec,unit=(u.hourangle, u.deg),frame='fk5')
        self.ra_string="%.6f" % round(self.sci_c.ra.deg, 6)
        self.dec_string="%.6f" % round(self.sci_c.dec.deg, 6)
        self.ref_width=self.sci_img_hdu.header['NAXIS2']*self.sci_ps*1.1/60.
        # self.ref_width*=4
        # print(self.ref_width)
        # print(image_size*self.sci_ps*1.1/60.)
        # self.ref_width = image_size*self.sci_ps*1.1/60.
        self.sci_mean, self.sci_median, self.sci_std = sigma_clipped_stats(self.sci_img_hdu.data)
        self.coords_sn_sci=wcs_to_pixels(self.sci_path,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
        self.coords_sn_sci_x,self.coords_sn_sci_y = self.coords_sn_sci
        print(info_g+' Object name:','\033[1m'+self.sci_obj+'\033[0m')
        if self.sci_prop!=None: print(info_g+' Object proposal ID:','\033[1m'+self.sci_prop+'\033[0m')
        print(info_g+' Object position (RA,Dec) J2000:','\033[1m'+' (',self.sci_ra,',',self.sci_dec,')'+ '\033[0m')
        print(info_g+' Object position (RA,Dec) deg:','\033[1m'+' (',self.ra_string,',',self.dec_string,')'+ '\033[0m')
        print(info_g+' Object position (x,y)','\033[1m'+' (',int(self.coords_sn_sci_x),',',int(self.coords_sn_sci_y),')'+ '\033[0m')
        print(info_g+' Image dimensions (x,y):','\033[1m'+' (',self.sci_img_hdu.header['NAXIS1'],'x',self.sci_img_hdu.header['NAXIS2'],')'+ '\033[0m')
        print(info_g+f' {self.ref_width}'+' arcmin reference width')
        print(info_g+' Observation date:','\033[1m'+str(round(self.sci_mjd,3))+'\033[0m')
        print(info_g+' Observation time:','\033[1m'+self.sci_utstart+'\033[0m')
        print(info_g+' Exposure time:','\033[1m'+str(self.sci_exp_time)+'\033[0m')
        self.name=self.sci_path

        
        if self.sci_obj in ['SN2023ixf','ZTF23aaklqou','2023ixf','23aaklqou'] or self.special_case!=None:
            self.special_case='SN2023ixf'
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
        





    def make_sdss_ref(self,sdss_filt='u'):
        if self.termoutp!='quiet':
            print(info_g+f" Making SDSS reference image")
        #if self.sci_filt=='u':
        self.spacing=0.1
        self.sdss_images=[]
        self.grid_size=3
        self.image_name=self.sci_obj+f'_ref_{sdss_filt}.fits'
        self.ref_path=self.path+'ref_imgs/'+self.image_name
        if os.path.exists(self.ref_path):
            print(info_b+' SDSS reference image already exists')
            return self.ref_path
        # if not os.path.exists(self.ref_path):
        #     #download a grid of images
        #     for nx in range(-3,2):
        #         for ny in range(-3,2):
        #             # print(nx,ny)
        #             self.sci_c.ra.deg+(nx*self.spacing),self.sci_c.dec.deg+(ny*self.spacing)
        #             # self.sdss_query_image(self.sci_c.ra.deg+(nx*self.spacing),self.sci_c.dec.deg+(ny*self.spacing),self.sci_filt,nx=nx,ny=ny)
        #             self.sdss_images.append(self.sdss_query_image(self.sci_c.ra.deg+(nx*self.spacing),self.sci_c.dec.deg+(ny*self.spacing),self.sci_filt,nx=nx,ny=ny))

        # if __name__ == '__main__':
        self.nx_range = range(-5,5)
        self.ny_range = range(-5,5)
        # sdss_args = [(nx, ny) for nx in self.nx_range for ny in self.ny_range]
# nx, ny,sci_c,sci_filt = args
        sdss_args = [(nx, ny,self.sci_c,self.sci_filt) for nx in self.nx_range for ny in self.ny_range]

        # for i in range(len(self.nx_range)):
        #     sdss_args.append((self.nx_range[i],self.ny_range[i],self.sci_c,self.sci_filt))


        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(get_image, sdss_args[i])#[0], sdss_args[i][1],sdss_args[i][2],sdss_args[i][3]) )
                    # for nx in range(-3, 2) 
                    # for ny in range(-3, 2)
                    for i in range(len(sdss_args))]
            
            self.sdss_images = [f.result() for f in futures]

        # print(self.sdss_images)
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
                    if self.termoutp!='quiet':
                        print(info_b+" %s\n" %(item))
            f.close()
    
        # self.image_name=self.sci_obj+f'_ref_{sdss_filt}.fits'
        # self.ref_path=self.path+'ref_imgs/'+self.image_name
        if not os.path.exists(self.ref_path):
            print(info_g+' Combining SDSS images with swarp centering on the object at '+self.ra_string+' '+self.dec_string)
            self.swarp_sdss_command=swarp_path+" @"+self.path+f'ref_imgs/{self.sci_obj}_sdss_list_{sdss_filt}.txt'+" -c "+self.path+"config_files/swarp_sdss.conf -CENTER '"+self.ra_string+" "+self.dec_string+"'"+" -IMAGE_SIZE '"+'4000'+","+'4000'+"'"+" -IMAGEOUT_NAME "+self.path+'ref_imgs/'+self.sci_obj+f'_ref_{sdss_filt}.fits -VERBOSE QUIET'
            print(self.swarp_sdss_command)
            # sys.exit(1)
            os.system(self.swarp_sdss_command)

        # print(self.ref_path)
        # sys.exit(1)
        return self.ref_path




    def bkg_subtract(self,sigma=3.):
        if self.termoutp!='quiet':
            print(info_g+f" Subtracting background")
        self.t_bkg_sub_start = time.time()
        ##############################################
            #  Background subtraction
        #################################
        if not os.path.exists(path+'bkg_subtracted_science'):
            os.makedirs(path+'bkg_subtracted_science')

        self.sig_clip = SigmaClip(sigma=sigma)
        self.bkg_estimator = SExtractorBackground(self.sig_clip)
        
        self.bkg = Background2D(self.sci_img, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=self.bkg_estimator)

        self.bkg_median = []
        for row in self.bkg.background:
            [self.bkg_median.append(i) for i in row]
            # row_median = stats.median(row)

        self.bkg_median = stats.median(self.bkg_median)
        if self.sci_bkg_med!=None and self.bkg_median<5 and self.sci_bkg_med<5 and self.sci_bkg_mean<5 and self.sci_filt!='u':  
            print(warn_r+f' Median of background counts is {self.bkg_median}, possible shutter issue')
            print(warn_r+' Exiting!!')
            self.sys_exit=True

            if self.output_to_log!=False:
                self.night_log = open(self.log_name,"r+")
                self.log_error_2_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+f': Background counts <threshold, median background {self.bkg_median}  '+'\n'
                if self.log_error_2_5 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_2_5)
                    self.night_log.write('\n')
                    self.night_log.close()
            
            return

        
        if self.to_subtract==True: #subtracting the background
            # print('subtractiog background',self.bkg.background)
            # print(self.bkg_median)
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

        if self.termoutp!='quiet':
            print(info_g+' Background subtraction time:', "%.1f" % round(self.bkg_sub_removal_time, 0),'seconds')

       
        self.bkg_error = self.bkg.background_rms
        self.bkg_err = calc_total_error(self.sci_bkgsb,self.bkg_error,self.sci_gain)
        return self.sig_clip,self.bkg,self.sci_bkgsb,self.median_bkg


    def remove_cosmic(self):
        if self.termoutp!='quiet':
            print(info_g+f" Removing cosmic rays")
        self.t_cosmic_start = time.time()
        if self.termoutp!='quiet':
            self.verbose = False
        else:
            self.verbose=True
        #based on how long the exposure times are
        if self.telescope not in SEDM:
            if 60<self.sci_exp_time<120.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=4, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)
    
            elif 30<self.sci_exp_time<=60.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=3, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)

            elif self.sci_exp_time>=120:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=6, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)

            elif self.sci_exp_time<=30:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=3.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=16.0, satlevel=45000.0, niter=2, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)

        if self.telescope in SEDM:
            if 60<self.sci_exp_time<120.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=4, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)
              
            elif 30<self.sci_exp_time<=60.0:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=3, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)

            elif self.sci_exp_time>=120:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=6, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)
                
            elif self.sci_exp_time<=30:
                self.sci_cr_mask,self.sci_bkgsb=astroscrappy.detect_cosmics(self.sci_bkgsb,sigclip=5.0, sigfrac=0.3, objlim=5.0, gain=self.sci_gain, readnoise=20.0, satlevel=45000.0, niter=2, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=self.verbose)
                

        self.sci_img_name=self.sci_img_name[:-5]+"bkgsub.fits"
        fits.writeto(self.path+'bkg_subtracted_science/'+self.sci_img_name, self.sci_bkgsb, self.sci_img_hdu.header,overwrite=True)

        self.files_to_clean.append(self.path+'bkg_subtracted_science/'+self.sci_img_name)

        self.t_cosmic_end = time.time()
        self.cosmic_removal_time = self.t_cosmic_end - self.t_cosmic_start
        if self.termoutp!='quiet':
            print(info_g+' Cosmic Removal time:',"%.1f" % round(self.cosmic_removal_time, 0),'seconds')
        # sys.exit(1)
        
        return self.sci_img_name

    def swarp_ref_align(self,image_size=image_size):
        

        if self.telescope in SEDM:
            self.image_size=3000
        else:
            self.image_size=image_size

        self.image_size=2500
        if self.termoutp!='quiet':
            print(info_g+f' Aligning science with reference image')

        if self.auto_ref=='auto':
            if self.sci_filt=='u' or self.use_sdss==True: #if filt ==u, need to call make sdss_ref before this or can call it here
                self.make_sdss_ref(sdss_filt=self.sci_filt)
                if self.sys_exit==True:
                    return
                self.image_name=self.sci_obj+f'_ref_{self.sci_filt}.fits'
                self.ref_path=self.path+'ref_imgs/'+self.image_name
                self.ref_img_name=self.ref_path
                print(info_g+f' Using SDSS reference image for alignment: ',self.ref_img_name)
                # print(self.ref_path)
                # sys.exit(1)
                panstamps_status=100 #no need to worry


            elif self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z':
                self.ref_folder = self.path+'ref_imgs/stack_'+str(self.sci_filt)
                self.ref_path = self.path+'ref_imgs/stack_'+str(self.sci_filt)+'_ra'+self.ra_string+'_dec'+self.dec_string+'_arcsec*'
                if len(glob.glob(self.ref_path))>0:
                    if self.termoutp!='quiet':
                        print(info_b+' PS1 reference image already exists :',glob.glob(self.ref_path)[0])
                        panstamps_status=0
                if len(glob.glob(self.ref_path))==0:
                    if self.termoutp!='quiet':
                        print(info_g+' Downloading reference image from PS...')
                    try:
                        # print(panstamps_path+' -f --width='+str(self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))
                        panstamps_status = os.system(panstamps_path+' -f --width='+str(2*self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))
                        # elif self.telescope in SEDM:  panstamps_status = os.system(panstamps_path+' -f --width='+str(2*self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))



                        print(info_g+' Reference image downloaded',panstamps_status)
                    except Exception as e:
                        print('error', e)
                        if self.output_to_log!=False:
                            self.log_error_3 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Unable to download reference image from panstarrs '+'\n'
                            self.night_log = open(self.log_name,"r+")
                            if self.log_error_3 in self.night_log.readlines():
                                pass
                            else:
                                self.night_log.write(self.log_error_3)
                                self.night_log.write('\n')
                                self.night_log.close()
                        if self.termoutp!='quiet':
                            print(warn_r+f" Error downloading PS reference image, error encountered :{e}")
                            print(warn_r+f" Exiting...")
                        ##sys.exit(1)
                        self.sys_exit=True
                        return
                
            if self.sys_exit==True:
                return

            # print(self.ref_path)
            if self.use_sdss==False and panstamps_status==0:
                c=0
                while len(glob.glob(self.ref_path))==0:
                    #the image has downloaded but the name is not the same as the one we want so we need to find it by shortening the name we search for
                    self.ref_path=self.ref_path[:-2]+'*'
                    # print(self.ref_path)
                    # print(glob.glob(self.ref_path))
                    c+=1
                    if c>10:
                        break
            # print(glob.glob(self.ref_path))
            # print(self.ref_path)
            self.ref_img_name=glob.glob(self.ref_path)[0]
            self.ref_img_hdu=fits.open(self.ref_img_name)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
        
        else:
            #If a reference image is passed in with it's full path as self.auto_ref
            if self.auto_ref.startswith(self.path):
                self.auto_ref = self.path+self.auto_ref

            # print(self.auto_ref)
            # print(self.sci_path)
            if not os.path.exists(self.auto_ref):
                print(warn_r+f' Reference image {self.auto_ref} not found')
                self.log_error_3 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+f': Unable to find reference image from path specified {self.auto_ref} '+'\n'
                self.night_log = open(self.log_name,"r+")
                if self.log_error_3 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_3)
                    self.night_log.write('\n')
                    self.night_log.close()
                
                self.sys_exit=True
                return

            self.ref_img_name=glob.glob(self.auto_ref)[0]
            self.ref_img_hdu=fits.open(self.auto_ref)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
            # print(self.ref_img==self.sci_bkgsb)

        self.coords_sn_ref=wcs_to_pixels(self.ref_img_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
        self.coords_sn_ref_x,self.coords_sn_ref_y = self.coords_sn_ref




        if not os.path.exists(path+'aligned_images'):
            os.makedirs(path+'aligned_images')

        
        wcs_command='python3 '+self.path+'subphot_align_quick.py'+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+"  -m relative -r 100"
        # wcs_command='python3 '+self.path+'subphot_align.py'+" -sci "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" -ref "+self.ref_img_name+"  -m relative -r 100"
        # print(wcs_command)

        self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
        self.ref_img_hdu=fits.open(self.ref_img_name)[0]
        self.align_success=False
        self.align_fail_count=0

        
        os.system(wcs_command)


        while self.align_success==False:
            self.align_fail_count+=1

            # if self.align_fail_count>1:
            #     #save the aligned images as the original names so Swarp can find them and align them
            #     fits.writeto(self.path+'aligned_images/'+self.sci_img_name.split('/')[-1],self.sci_ali_img_hdu.data,overwrite=True)
            #     fits.writeto(self.path+'aligned_images/'+self.ref_img_name.split('/')[-1],self.ref_ali_img_hdu.data,overwrite=True)
            # try:
            # except Exception as e:
            #     if self.termoutp!='quiet':
            #         print(warn_y+" Error with alignment...")
                

            print(info_g+" Aligning images with SWarp making a resampled image of size "+str(self.image_size)+"x"+str(self.image_size)+" pixels")
            #image size is from credentials file
            #time swarp
            # swarp_time = time.time()
            # if self.telescope=='SEDM-P60' or self.telescope=='SEDM':
            #     #find the center of the science image and convert to degree
            #     self.sci_y_center,self.sci_x_center = np.shape(self.sci_img_hdu.data)[0]/2,np.shape(self.sci_img_hdu.data)[1]/2
            #     print(info_g+" Science image center is at pixel coordinates "+str(self.sci_x_center)+","+str(self.sci_y_center))
            #     #convert to ra and dec
            #     self.sci_wcs = WCS(self.sci_img_hdu.header)
            #     self.sci_c = self.sci_wcs.pixel_to_world(self.sci_x_center,self.sci_y_center)
            #     # ax = fig.add_subplot(111)

            #     # wcs.utils.pixel_to_skycoord(self.sci_x_center,self.sci_y_center,self.sci_img_hdu.header)
            #     self.ali_center_ra,self.ali_center_dec = str(self.sci_c.ra.deg),str(self.sci_c.dec.deg)

            #     # print(self.ali_center_ra,self.ali_center_dec)
            #     print(info_g+" Science image center is at ra,dec coordinates "+self.ali_center_ra+","+self.ali_center_dec)
            #     # sys.exit()

            #     # self.ali_center_ra,self.ali_center_dec = self.sci_c.ra.deg,self.sci_c.dec.deg
            # else:
            self.ali_center_ra,self.ali_center_dec = self.ra_string,self.dec_string
            if float(self.dec_string)>0:self.ali_center_dec='+'+self.dec_string
            elif self.dec_string[0]=='-': self.ali_center_dec=self.dec_string 
            else:self.ali_center_dec='-'+self.dec_string
            
            # print(self.ali_center_ra,self.ali_center_dec)
            print(info_g+" Aligning images to a centre of "+self.ali_center_ra+","+self.ali_center_dec)

            self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
            self.ref_img_hdu=fits.open(self.ref_img_name)[0]

            #find centre of original reference image
            self.ref_img_size = np.shape(self.ref_img_hdu.data)
            self.ref_img=self.ref_img_hdu.data
            self.ref_img_center = np.shape(self.ref_img)[0]/2,np.shape(self.ref_img)[1]/2
            self.ref_wcs = WCS(self.ref_img_hdu.header)
            self.ref_cnt = self.ref_wcs.all_pix2world(self.ref_img_center[0],self.ref_img_center[1],1)
            self.ref_img_center_ra,self.ref_img_center_dec = self.ref_cnt[0],self.ref_cnt[1]

            self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
            self.sci_img_size = np.shape(self.sci_img_hdu.data)
            self.sci_img=self.sci_img_hdu.data
            self.sci_img_center = np.shape(self.sci_img)[0]/2,np.shape(self.sci_img)[1]/2
            self.sci_wcs = WCS(self.sci_img_hdu.header)
            self.sci_cnt = self.sci_wcs.all_pix2world(self.sci_img_center[0],self.sci_img_center[1],1)
            self.sci_img_center_ra,self.sci_img_center_dec = self.sci_cnt[0],self.sci_cnt[1]

            self.image_size1,self.image_size2=np.shape(self.sci_img_hdu.data)[0],np.shape(self.sci_img_hdu.data)[1]
            # self.image_size1,self.image_size2=np.shape(self.ref_img_hdu.data)[0],np.shape(self.ref_img_hdu.data)[1]

            # print(np.shape(self.ref_img_hdu.data),np.shape(self.sci_img_hdu.data))
            # swarp_command=swarp_path+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+" -c "+self.path+"config_files/config.swarp -CENTER '"+str(self.ali_center_ra)+" "+str(self.ali_center_dec)+"' -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -RESAMPLE Y -RESAMPLE_DIR '"+self.path+"aligned_images/' -COMBINE N -IMAGE_SIZE '"+str(self.image_size)+","+str(self.image_size)+"'" #original
            # if self.telescope in SEDM:
                # if 'ACQ' not in self.sci_img_name:
            
            swarp_command=swarp_path+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+" -c "+self.path+"config_files/config.swarp -CENTER '"+str(self.ali_center_ra)+" "+str(self.ali_center_dec)+"' -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -RESAMPLE Y -RESAMPLE_DIR '"+self.path+"aligned_images/' -COMBINE N -IMAGE_SIZE '"+str(self.image_size)+","+str(self.image_size)+"'" #ref centre
                # elif 'ACQ' in self.sci_img_name:
                    #move the image to the centre of the science image
                    # swarp_command=swarp_path+" "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" "+self.ref_img_name+" -c "+self.path+"config_files/config.swarp -CENTER '"+str(self.sci_img_center_ra)+" "+str(self.sci_img_center_dec)+"' -SUBTRACT_BACK Y -VERBOSE_TYPE QUIET -RESAMPLE Y -RESAMPLE_DIR '"+self.path+"aligned_images/' -COMBINE N -IMAGE_SIZE '"+str(self.image_size1*2)+","+str(self.image_size2*2)+"'" #sci centre



            # fits.writeto(self.path+'aligned_images/stack_g_ra8.690419_dec-8.386699_arcsec366_skycell1056.093.resamp.fits',self.ref_img_hdu.data,self.ref_img_hdu.header,overwrite=True)
            status=os.system(swarp_command)
            print(swarp_command)
            # sys.exit()
            # print(info_g+" SWarp took %2f seconds" % (time.time() - swarp_time)    )
            print(info_g+" SWarp status: "+str(status))

            



            
            self.sci_ali_name=glob.glob(self.path+'aligned_images/'+self.sci_img_name[:-5]+'.resamp.fits')[0]
            

            if self.auto_ref=="auto" and self.sci_filt in ['g','r','i','z'] and self.use_sdss==False:
                # print(self.path+'aligned_images/stack_'+str(self.sci_filt)+'_ra'+self.ra_string+'_dec'+self.dec_string+'_arcsec*.resamp.fits')
                # print(glob.glob(self.path+'aligned_images/stack_'+str(self.sci_filt)+'_ra'+self.ra_string+'_dec'+self.dec_string+'_arcsec*.resamp.fits'))
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
                print(info_g+' Alignment attempt 1')
            else: 
                print(warn_y+' Alignment attempt '+str(self.align_fail_count))


            # print(info_g+" Aligned science image: "+self.sci_ali_name,np.shape(fits.open(self.sci_ali_name)[0].data))
            # print(info_g+" Aligned reference image: "+self.ref_ali_name,np.shape(fits.open(self.ref_ali_name)[0].data))

            self.ref_ali_img_hdu = fits.open(self.ref_ali_name)[0]
            self.sci_ali_img_hdu = fits.open(self.sci_ali_name)[0]
            print(info_g+" Aligned science image: "+self.sci_ali_name,np.shape(self.sci_ali_img_hdu.data))
            print(info_g+" Aligned reference image: "+self.ref_ali_name,np.shape(self.ref_ali_img_hdu.data))
            


            print(info_g+' Science & Reference image sizes : ',np.shape(fits.open(self.sci_ali_name)[0].data),np.shape(fits.open(self.ref_ali_name)[0].data))
            print(info_b+' Target size : ',self.image_size,'x',self.image_size)
            print(info_g+' Aligned Science image size : ',np.shape(fits.open(self.sci_ali_name)[0].data))
            print(info_g+' Aligned Reference image size : ',np.shape(fits.open(self.ref_ali_name)[0].data))

            

            if self.align_fail_count>15:
                print(warn_y+' Alignment failed 15 times')
                break

            if np.shape(fits.open(self.sci_ali_name)[0].data)!=np.shape(fits.open(self.ref_ali_name)[0].data) and ((self.image_size,self.image_size)!=np.shape(fits.open(self.sci_ali_name)[0].data) or (self.image_size,self.image_size)!=np.shape(fits.open(self.ref_ali_name)[0].data)):
                print(warn_y+' Science image or reference image and aligned image have different sizes. This may be due to the reference image being too far away from the science image')
                self.min_align_size = np.min([np.shape(fits.open(self.sci_ali_name)[0].data),np.shape(fits.open(self.ref_ali_name)[0].data)])
                print(warn_y+' Trying image size of '+str(self.min_align_size)+'x'+str(self.min_align_size)+' pixels')
                self.image_size=self.min_align_size
                self.align_success=False

                try:
                    # if self.telescope in SEDM:
                    #find centre of original reference image
                    self.ref_img_hdu=fits.open(self.ref_img_name)[0]
                    self.ref_img_size = np.shape(self.ref_img_hdu.data)
                    self.ref_img=self.ref_img_hdu.data
                    self.ref_img_center = np.shape(self.ref_img)[0]/2,np.shape(self.ref_img)[1]/2
                    self.ref_wcs = WCS(self.ref_img_hdu.header)
                    self.ref_cnt = self.ref_wcs.all_pix2world(self.ref_img_center[0],self.ref_img_center[1],1)
                    self.ref_img_center_ra,self.ref_img_center_dec = self.ref_cnt[0],self.ref_cnt[1]

                    self.sci_img_hdu=fits.open(self.path+'bkg_subtracted_science/'+self.sci_img_name)[0]
                    self.sci_img_size = np.shape(self.sci_img_hdu.data)
                    self.sci_img=self.sci_img_hdu.data
                    self.sci_img_center = np.shape(self.sci_img)[0]/2,np.shape(self.sci_img)[1]/2
                    self.sci_wcs = WCS(self.sci_img_hdu.header)
                    self.sci_cnt = self.sci_wcs.all_pix2world(self.sci_img_center[0],self.sci_img_center[1],1)
                    self.sci_img_center_ra,self.sci_img_center_dec = self.sci_cnt[0],self.sci_cnt[1]

                    self.sci_ali_img_hdu = fits.open(self.sci_ali_name)[0]
                    self.sci_ali_img_size = np.shape(self.sci_ali_img_hdu.data)
                    # self.sci_path,_ = self.update_astrometry(self.sci_path,ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])

                    # sys.exit()
                    requested_size = np.shape(self.ref_ali_img_hdu.data)
                    # requested_size = (self.image_size1,self.image_size2)
                    # print('requested_size',requested_size)
                    print(info_g+' Padding science image to match reference image size of '+str(requested_size[0])+'x'+str(requested_size[1])+' pixels')
                    cominx,cominy = int(self.sci_ali_img_hdu.header['COMIN1']),int(self.sci_ali_img_hdu.header['COMIN2'])
                    try:
                        print(info_g+' Science COMIN1,COMIN2 : ',cominx,cominy)
                        pad_sci = np.zeros((requested_size[0],requested_size[1]))*np.nan
                        pad_sci[cominy:cominy+self.sci_ali_img_size[0],cominx:cominx+self.sci_ali_img_size[1]] = self.sci_ali_img_hdu.data
                        # pad_ref[cominy:cominy+self.ref_ali_img_size[0],cominx:cominx+self.ref_ali_img_size[1]] = self.ref_ali_img_hdu.data
                        new_sci_img = pad_sci
                        print(info_g+' Successfully padded science image')
                    except Exception as e:
                        try:
                            print(warn_y+f" Error padding science image, error encountered, changing cominy+=1",e)
                            pad_sci = np.zeros((requested_size[0],requested_size[1]))*np.nan
                            cominy-=1
                            cominx-=1
                            pad_sci[cominy:cominy+self.sci_ali_img_size[0],cominx:cominx+self.sci_ali_img_size[1]] = self.sci_ali_img_hdu.data
                            new_sci_img = pad_sci
                            print(info_g+' Successfully padded science image')
                            # sys.exit()
                        except Exception as e:
                            print(warn_y+f" Failed padding, trying to Swarp again",e)
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
                    #     try:
                    #         print(warn_y+f" Error padding science image, error encountered, changing cominy-=1",e)
                    #         # cominy-=1
                    #         pad_sci = np.zeros((requested_size[0],requested_size[1]))*np.nan

                    #         pad_sci[cominy:cominy+self.sci_ali_img_size[0]-1,cominx:cominx+self.sci_ali_img_size[1]-1] = self.sci_ali_img_hdu.data
                    #         new_sci_img = pad_sci
                    #         self.sci_ali_img_hdu.data = new_sci_img
                    #     except Exception as e:
                    #         print(warn_y+f" Failed changing cominy-=1, changing cominx-=1",e)
                    #         try:
                    #             cominx-=1
                    #             cominy+=1
                    #             pad_sci[cominy:cominy+self.sci_ali_img_size[0],cominx:cominx+self.sci_ali_img_size[1]] = self.sci_ali_img_hdu.data
                    #             new_sci_img = pad_sci
                    #             self.sci_ali_img_hdu.data = new_sci_img
                    #         except Exception as e:
                    #             print(warn_y+f" Failed changing cominx-=1, changing cominx+=1",e)

                    # print('new_sci_img shape',np.shape(new_sci_img))
                    # fig = plt.figure(figsize=(10,10))
                    # plt.subplot(2,2,1)
                    # vmin,vmax = visualization.ZScaleInterval().get_limits(self.sci_img_hdu.data)
                    # plt.imshow(self.sci_img_hdu.data,cmap='gray',vmin=vmin,vmax=vmax)
                    # plt.title('Original Science image')
                    # plt.scatter(self.sci_img_center[0],self.sci_img_center[1],marker='x',color='r')
                    # plt.subplot(2,2,2)
                    # vmin,vmax = visualization.ZScaleInterval().get_limits(self.ref_img_hdu.data)
                    # plt.imshow(self.ref_img_hdu.data,cmap='gray',vmin=vmin,vmax=vmax)
                    # plt.title('Original Reference image')
                    # plt.scatter(self.ref_img_center[0],self.ref_img_center[1],marker='x',color='r')
                    # plt.subplot(2,2,3)
                    # vmin,vmax = visualization.ZScaleInterval().get_limits(self.sci_ali_img_hdu.data)
                    # plt.imshow(self.sci_ali_img_hdu.data,cmap='gray',vmin=vmin,vmax=vmax)
                    # plt.title('Aligned Science image')
                    # plt.scatter(self.sci_ali_img_hdu.data.shape[0]/2,self.sci_ali_img_hdu.data.shape[1]/2,marker='x',color='r')
                    # plt.subplot(2,2,4)
                    # vmin,vmax = visualization.ZScaleInterval().get_limits(self.ref_ali_img_hdu.data)
                    # plt.imshow(self.ref_ali_img_hdu.data,cmap='gray',vmin=vmin,vmax=vmax)
                    # plt.title('Aligned Reference image')
                    # plt.scatter(self.ref_ali_img_hdu.data.shape[0]/2,self.ref_ali_img_hdu.data.shape[1]/2,marker='x',color='r')
                    # # plt.savefig(self.path+'aligned_images/'+self.sci_img_name[:-5]+'.png')
                    # # plt.savefig(self.path+'aligned_images/'+self.sci_img_name[:-5]+'.png')
                    # # plt.show()
                    # plt.close()
                    #save the aligned images as the original names so Swarp can find them and align them
                    #change the sci header to match new wcs
                    for keyw in ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2']:
                        self.sci_ali_img_hdu.header[keyw] = self.ref_ali_img_hdu.header[keyw]


                    # print(self.sci_ali_img_hdu.data)
                    self.sci_ali_img_hdu.data = np.nan_to_num(self.sci_ali_img_hdu.data)#np.nanmedian(self.sci_ali_img_hdu.data))
                    fits.writeto(self.path+'aligned_images/'+self.sci_img_name.split('/')[-1].replace('.fits','_padded.fits'),self.sci_ali_img_hdu.data,overwrite=True,header=self.sci_ali_img_hdu.header)
                    # self.sci_ali_name,_ = self.update_astrometry(self.path+'aligned_images/'+self.sci_img_name.split('/')[-1].replace('.fits','_padded.fits'),ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw])
                    # print('new sci ali name',self.sci_ali_name)
                    # sys.exit()

                    self.sci_ali_name = self.path+'aligned_images/'+self.sci_img_name.split('/')[-1].replace('.fits','_padded.fits')
                    self.align_success=True
                    print(info_g+' Alignment successful')
                    # self.align_fail_count=+1
                except Exception as e:
                    print(warn_r+' Alignment failed',e)
                    self.align_success=False
                    self.align_fail_count=+1
                    # sys.exit()
                    # self.align_fail_count=+1


            # else:
            self.align_success=True
            print(info_g+' Alignment successful')
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

            self.ax[0].imshow(self.sci_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_sci, vmax=self.vmax_sci)

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

        if self.termoutp!='quiet':
            print(info_g+f' Aligning science with reference image')

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
                    if self.termoutp!='quiet':
                        print(info_b+' PS1 reference image already exists')
                if len(glob.glob(self.ref_path))==0:
                    if self.termoutp!='quiet':
                        print(info_g+' Downloading reference image from PS...')
                    try:
                        # print(panstamps_path+' -f --width='+str(self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))
                        os.system(panstamps_path+' -f --width='+str(self.ref_width)+' --filters='+self.sci_filt+' --downloadFolder='+self.path+'ref_imgs stack '+str(self.sci_c.ra.deg)+' '+str(self.sci_c.dec.deg))

                        print(info_g+' Reference image downloaded')
                    except Exception as e:
                        print('error', e)
                        if self.output_to_log!=False:
                            self.log_error_3 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Unable to download reference image from panstarrs '+'\n'
                            self.night_log = open(self.log_name,"r+")
                            if self.log_error_3 in self.night_log.readlines():
                                pass
                            else:
                                self.night_log.write(self.log_error_3)
                                self.night_log.write('\n')
                                self.night_log.close()
                        if self.termoutp!='quiet':
                            print(warn_r+f" Error downloading PS reference image, error encountered :{e}")
                            print(warn_r+f" Exiting...")
                        ##sys.exit(1)
                        self.sys_exit=True
                        return
                
            if self.sys_exit==True:
                return

            self.ref_path = glob.glob(self.ref_path)[0]
            # print(glob.glob(self.ref_path))
            # print(self.ref_path)
            self.ref_img_name=glob.glob(self.ref_path)[0]
            self.ref_img_hdu=fits.open(self.ref_img_name)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
        
        else:
            #If a reference image is passed in with it's full path as self.auto_ref
            if self.auto_ref.startswith(self.path):
                self.auto_ref = self.path+self.auto_ref

            # print(self.auto_ref)
            # print(self.sci_path)
            if not os.path.exists(self.auto_ref):
                print(warn_r+f' Reference image {self.auto_ref} not found')
                self.log_error_3 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+f': Unable to find reference image from path specified {self.auto_ref} '+'\n'
                self.night_log = open(self.log_name,"r+")
                if self.log_error_3 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_3)
                    self.night_log.write('\n')
                    self.night_log.close()
                
                self.sys_exit=True
                return

            self.ref_img_name=glob.glob(self.auto_ref)[0]
            self.ref_img_hdu=fits.open(self.auto_ref)[0]
            self.ref_img=self.ref_img_hdu.data[0:self.ref_img_hdu.header['NAXIS2'],0:self.ref_img_hdu.header['NAXIS1']] 
            # print(self.ref_img==self.sci_bkgsb)

        self.coords_sn_ref=wcs_to_pixels(self.ref_img_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
        self.coords_sn_ref_x,self.coords_sn_ref_y = self.coords_sn_ref




        if not os.path.exists(path+'aligned_images'):
            os.makedirs(path+'aligned_images')

        wcs_command='python3 '+self.path+'subphot_align.py'+" -sci "+self.path+'bkg_subtracted_science/'+self.sci_img_name+" -ref "+self.ref_img_name+"  -m relative -r 100"


        self.align_success=False
        self.align_fail_count=0

        
        os.system(wcs_command)
        # sys.exit(1)


        while self.align_success==False:
            self.align_fail_count+=1
            print(info_g+" Aligning images with reprojection and astroalgin making a resampled image of size "+str(self.image_size)+"x"+str(self.image_size)+" pixels")

            

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

            print(info_g+f" Reprojecting the reference onto the science image")
            self.ref_resampled_o, self.footprint = reproject_interp(self.ref_img_hdu, self.sci_img_header,order=1)
            self.ref_resampled = self.ref_resampled_o
            self.footprint = self.footprint.astype(bool)==False

            self.ref_resampled[np.isnan(self.ref_resampled)] = np.nanmedian(self.ref_resampled)

            try:
                # self.registered, self.footprint = aa.register(self.ref_resampled*self.footprint, self.sci_img.data)
                print(info_g+" Aligning images with astroalign")
                self.registered, self.footprint = aa.register(self.ref_resampled, self.sci_img.data)
            except Exception as e:
                print(warn_y+" Extra alignment with astroalign failed, using reprojected image instead")
                print(e)
            self.registered = self.ref_resampled


            self.ref_masked = np.ma.masked_array(self.registered, self.footprint, fill_value=np.nanmedian(self.ref_resampled)).filled()
            self.ref_masked[np.isnan(self.ref_masked)] = np.nanmedian(self.ref_resampled)


            self.ref_ali_hdu = fits.PrimaryHDU(self.ref_masked,header=self.sci_img_header)

            print(info_g+" Saving aligned images to "+self.path+'aligned_images/')
            self.sci_ali_name=self.path+'aligned_images/'+self.sci_img_name[:-5]+'.resamp.fits'
            self.ref_ali_name=self.path+'aligned_images/'+self.ref_img_name[:-5].split('/')[-1]+'.resamp.fits'

            self.ref_ali_hdu.writeto(self.ref_ali_name,overwrite=True)
            self.sci_img_hdu.writeto(self.sci_ali_name,overwrite=True)

        
            if self.align_fail_count==1:
                print(info_g+' Alignment attempt 1')
            else: 
                print(warn_y+' Alignment attempt '+str(self.align_fail_count))




            self.ref_ali_img_hdu = fits.open(self.ref_ali_name)[0]
            self.sci_ali_img_hdu = fits.open(self.sci_ali_name)[0]

            
            


            print(info_g+' Science & Reference image sizes : ',np.shape(fits.open(self.sci_ali_name)[0].data),np.shape(fits.open(self.ref_ali_name)[0].data))
            print(info_b+' Target size : ',self.image_size,'x',self.image_size)
            print(info_g+' Aligned Science image size : ',np.shape(fits.open(self.sci_ali_name)[0].data))
            print(info_g+' Aligned Reference image size : ',np.shape(fits.open(self.ref_ali_name)[0].data))

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
                # plt.show()
                # sys.exit()
                plt.close(ali_fig)

            if self.align_fail_count>15:
                print(warn_y+' Alignment failed 15 times')
                break

            self.align_success=True
            print(info_g+' Alignment successful')

                
        
            if self.cutout_tf!=False:
                self.fig,self.ax = plt.subplots(nrows=1, ncols=3)
                self.sci_img_ali_hdu = fits.open(self.sci_ali_name)[0]
                self.percentile = np.percentile(self.sci_img_ali_hdu.data,[5,10,20,30,40,50,60,70,80,90]) 
                self.cmap = 'gray'
                self.coords_sn_sci_ali=wcs_to_pixels(self.sci_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
                self.coords_sn_sci_ali_x,self.coords_sn_sci_ali_y = self.coords_sn_sci_ali
                self.vmin_sci,self.vmax_sci = visualization.ZScaleInterval().get_limits(self.sci_img_ali_hdu.data)

                self.ax[0].imshow(self.sci_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_sci, vmax=self.vmax_sci)

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
                self.ax[1].set_xlim(self.coords_sn_ref_ali_x-115,self.coords_sn_ref_ali_x+115)
                self.ax[1].set_ylim(self.coords_sn_ref_ali_y-115,self.coords_sn_ref_ali_y+115)
                self.ax[1].set_title('Ref')
                self.ax[1].axis('off')        

                self.fig.savefig(self.cutout_name)


            return self.ref_ali_name,self.ref_img_hdu,self.ref_img, self.sci_ali_name

            # except Exception as e:
            #     print(warn_y+f" Alignment failed with error: {e}")
            #     # self.align_fail_count+=1
            #     print(warn_y+f" Trying again with different parameters. Attempt {self.align_fail_count} of {20}")
            #     if self.align_fail_count==20:
            #         print(warn_y+f" Alignment failed 20 times. Aborting")
            #         self.align_success=False
            #         return None,None,None,None


    def psfex_convolve_images(self):
        if self.termoutp!='quiet':
            print(info_g+' Beginning to convolve images with SeXtractor and PSFEx')

        if not os.path.exists(path+'convolved_sci'):
            os.makedirs(path+'convolved_sci')

        if not os.path.exists(path+'convolved_ref'):
            os.makedirs(path+'convolved_ref')

        if not os.path.exists(path+'out'):
            os.makedirs(path+'out')

        self.sci_conv_name=self.path+"convolved_sci/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'sci_convolved.fits'
        self.ref_conv_name=self.path+"convolved_ref/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'ref_convolved.fits'
        
        self.files_to_clean.append(self.sci_conv_name)
        self.files_to_clean.append(self.ref_conv_name)
         #################################################
        # CONVOLVE REFERENCE WITH PSF OF SCIENCE IMAGE

        if os.path.exists(path+f"config_files/prepsfex_{self.rand_nums_string}.cat"):
            os.system("rm "+path+f"config_files/prepsfex_{self.rand_nums_string}.cat")

        if self.termoutp!='quiet':
            print(info_g+' Convolving the reference with the PSF of the science image')

        # SExtractor command for the science image
        sextractor_command=sex_path+" "+self.sci_ali_name+" -c "+self.path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT 25.0"

        os.system(sextractor_command)
        

        if os.path.exists(path+f'out/proto_prepsfex_{self.rand_nums_string}.fits'):
            os.system('rm '+self.path+f'out/proto_prepsfex_{self.rand_nums_string}.fits')

        os.system(psfex_path+" "+path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")
        self.files_to_clean.append(path+f"config_files/prepsfex_{self.rand_nums_string}.cat")
        # print(psfex.PSFEx(path+f'config_files/prepsfex_{self.rand_nums_string}.cat'))
        # sys.exit(1)


        self.psf_sci_image_name=self.path+f'out/proto_prepsfex_{self.rand_nums_string}.fits'
        self.files_to_clean.append(path+f'out/proto_prepsfex_{self.rand_nums_string}.fits')
        self.psf_sci_image = fits.open(self.psf_sci_image_name)

        # print(self.psf_sci_image[0].data)
        # fig = plt.figure(figsize=(12,8))
        # plt.imshow(self.psf_sci_image[0].data[0])
        # fig.savefig('ZTF21aceqrju_g_psf.png')
        # plt.show()

        self.hdu_psf_model_sci= fits.open(self.path+f'out/prepsfex_{self.rand_nums_string}.psf')
        self.files_to_clean.append(self.path+f'out/prepsfex_{self.rand_nums_string}.psf')
        self.chi_sq_psf=self.hdu_psf_model_sci[1].header['CHI2']

        if self.termoutp!='quiet':
            print(colored(info_g+' Reduced Chi^2 of science image PSF fit:','green'),"%.1f" % self.chi_sq_psf)
            if self.chi_sq_psf>3:
                # print(colored('Warning: PSF model may not be accurate','green'))
                print(warn_y+' Warning: PSF model may not be accurate')

        #Convolve the reference image with a Gaussian kernel
        self.kernel_sci = self.psf_sci_image[0].data[0]

        if self.to_subtract!=False:
        # Read the REFERENCE image and convolve it with the  kernel science
            self.ref_image_aligned=fits.open(self.ref_ali_name)
            self.ref_conv = scipy_convolve(self.ref_image_aligned[0].data, self.kernel_sci, mode='same', method='fft')
            self.ref_conv=np.nan_to_num(self.ref_conv)

            fits.writeto(self.ref_conv_name, data=self.ref_conv, header=self.ref_image_aligned[0].header,overwrite=True)


            #################################################
            #CONVOLVE SCIENCE WITH PSF OF REFERENCE IMAGE

            os.system("rm "+path+f"config_files/prepsfex_{self.rand_nums_string}.cat")

            if self.termoutp!='quiet':
                print(info_g+' Convolving the science with the PSF of the reference image')

            sextractor_command=sex_path+" "+self.ref_ali_name+" -c "+self.path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT 25.0"
            os.system(sextractor_command)

            if self.sci_obj in ['2023ixf','SN2023ixf']:
                #open the sextractor catalog and keep only the point sources that are in calstars4e.cat
                self.good_ps1_ixf_stars = pd.read_csv(path+'calstars4e.cat',names=['ra','dec','u','g','r','i','z','PS1g','PS1r','PS1i','PS1z'],header=None,sep='\s+')
                self.good_ps1_ixf_stars = Table.from_pandas(self.good_ps1_ixf_stars)
                print(info_g+' Using calstars4e.cat to select point sources for PSFEx for SN2023ixf')
                self.good_ps1_ixf_stars = Table.from_pandas(self.good_ps1_ixf_stars)
                self.good_ps1_ixf_stars['ra'].unit = u.deg
                self.good_ps1_ixf_stars['dec'].unit = u.deg

                # self.sextracted_ref = sextract(self.ref_ali_name, 0, 0, 3, 12, maxellip=0.7, saturation=-1)
                # self.sex_ref = pd.DataFrame([i.row for i in self.sextracted_ref],columns=['x','y','ra', 'dec', 'mag', 'magerr', 'ellip', 'fwhm', 'flag', 'ra_rad', 'dec_rad'])

                self.sex_ref_orig = fits.open(path+f"config_files/prepsfex_{self.rand_nums_string}.cat" )
                # self.sex_ref = Table.read(path+f"config_files/prepsfex_{self.rand_nums_string}.cat", format="fits",hdu='LDAC_OBJECTS')
                self.sex_ref = Table(self.sex_ref_orig[2].data)
                # print(self.sex_ref)
                self.ref_ali_wcs = WCS(self.ref_ali_name)   
                #convert X_image and Y_image to ra and dec
                self.sex_ref_wcs = self.ref_ali_wcs.all_pix2world(np.column_stack((self.sex_ref['X_IMAGE'],self.sex_ref['Y_IMAGE'])),1)
                # self.ref_ali_photTab_coords = SkyCoord(ra=self.ref_ali_photTab_coords[:,0]*u.deg, dec=self.ref_ali_photTab_coords[:,1]*u.deg,frame='fk5')
                        # self.ref_ali_photTab['ra'],self.ref_ali_photTab['dec'] = self.ref_ali_photTab_coords.ra,self.ref_ali_photTab_coords.dec
                self.sex_ref_wcs = SkyCoord(ra=self.sex_ref_wcs[:,0]*u.deg, dec=self.sex_ref_wcs[:,1]*u.deg,frame='fk5')
                self.sex_ref['RA_DEG'],self.sex_ref['DEC_DEG'] = self.sex_ref_wcs.ra,self.sex_ref_wcs.dec
                # print(self.sex_ref)
                # print(np.array(self.sex_ref[['ra','dec']]))
                self.sex_ref_wcs = SkyCoord(ra=np.array(self.sex_ref['RA_DEG'])*u.deg, dec=np.array(self.sex_ref['DEC_DEG'])*u.deg,frame='fk5')
                self.good_ps1_ixf_stars_wcs = SkyCoord(ra=np.array(self.good_ps1_ixf_stars['ra'])*u.deg, dec=np.array(self.good_ps1_ixf_stars['dec'])*u.deg,frame='fk5')

                self.indx, self.d2d, self.d3d = self.sex_ref_wcs.match_to_catalog_sky(self.good_ps1_ixf_stars_wcs)
                #where stars match by search rad in arcseconds!
                self.upd_indx=np.where(self.d2d<search_rad/3600.*u.deg)[0]
                # print(self.upd_indx)
                print(info_g+' Bright stars found=',len(self.upd_indx),',',round(3600*(self.ref_width/60),6), 'arcsec search radius')
                self.bright_stars_sc = Table()
                self.bright_stars_sc['ra'] = self.good_ps1_ixf_stars_wcs.ra.deg
                self.bright_stars_sc['dec'] = self.good_ps1_ixf_stars_wcs.dec.deg
                self.bright_stars_sc_pix = self.ref_ali_wcs.all_world2pix(np.column_stack((self.bright_stars_sc['ra'],self.bright_stars_sc['dec'])),1)
                self.bright_stars_sc_pix_x,self.bright_stars_sc_pix_y = self.bright_stars_sc_pix[:,0],self.bright_stars_sc_pix[:,1]
                self.bright_stars_sc['xcentroid'],self.bright_stars_sc['ycentroid'] = self.bright_stars_sc_pix_x,self.bright_stars_sc_pix_y

                for col in ['u','g','r','i','z','PS1g','PS1r','PS1i','PS1z']:
                    self.bright_stars_sc[col] = self.good_ps1_ixf_stars[col]
                    # set columns to float
                    # self.bright_stars_sc[col] = self.bright_stars_sc[col].astype(float)

                # for i in self.upd_indx:
                #     print(self.sex_ref_wcs[i].ra.deg,self.sex_ref_wcs[i].dec.deg,self.sex_ref['X_IMAGE','Y_IMAGE'][i].values)
                #     print(self.good_ps1_ixf_stars_wcs[self.indx[i]])
                    # print( )
                # sys.exit(1)
                # self.star_coords_pix=self.star_coords_pix[self.upd_indx]



                self.crossmatch_sex_ref = self.sex_ref[self.upd_indx]
                print(self.crossmatch_sex_ref['X_IMAGE','Y_IMAGE'])
                # self.crossmatch_sex_ref['']


                # self.sex_ref_orig[2].data = self.sex_ref_orig[2].data[self.upd_indx]

                # for col in ['X_IMAGE','Y_IMAGE','VIGNET','FLUX_APER','FLUXERR_APER','FLUX_RADIUS','ELONGATION','FLAGS','SNR_WIN']:
                #     col_ = self.crossmatch_sex_ref[col]
                #     self.sex_ref_orig[2].data[col] = fits.column(array=col_.value, name=col_.name, unit=col_.unit,)

                self.sex_ref_orig[2].data = self.sex_ref_orig[2].data[self.upd_indx]
                self.sex_ref_orig[2].data['FLAGS'] = 0
                # print(self.sex_ref_orig[2].data['X_IMAGE'])

                # for col in self.sex_ref_orig[2].columns.names:
                #     for ind in self.upd_indx:
                #         self.sex_ref_orig[2].data[col] = self.sex_ref_orig[2].data[col][self.upd_indx] 

                #delete rows if not in upd_indx
                # for k in range(len(self.sex_ref_orig[2].data)):
                #     if k not in self.upd_indx:
                        

                self.sex_ref_orig[2].header['NAXIS2'] = len(self.upd_indx)

                print()



                self.sex_ref_orig.writeto(path+f"config_files/prepsfex_{self.rand_nums_string}.cat",overwrite=True)
            # os.system(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE FULL -OUTCAT_NAME "+self.path+f"out/proto_prepsfex_{self.rand_nums_string}.cat")
            # os.system(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}_new.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE FULL -OUTCAT_NAME "+self.path+f"out/proto_prepsfex_{self.rand_nums_string}_new.cat")


                #write the crossmatched stars to a file that is a fits_ldac catalog
                # import astromatic_wrapper as aw
                # aw.utils.ldac.save_table_as_ldac(self.crossmatch_sex_ref,path+f"config_files/prepsfex_{self.rand_nums_string}.cat",overwrite=True)#,format='fits_ldac')





            # sys.exit(1)
            if os.path.exists(path+f'out/proto_prepsfex_{self.rand_nums_string}.fits'):
                os.system('rm '+path+f'out/proto_prepsfex_{self.rand_nums_string}.fits')

            

            # os.system(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex")


            # sextractor_command=sex_path+" "+self.ref_ali_name+" -c "+self.path+"config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"config_files/pooprepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT 25.0"
            # print(sextractor_command)

            # print(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET")
            # sys.exit(1)
            # print(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex -VERBOSE_TYPE FULL")
            os.system(psfex_path+" "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -c "+self.path+"config_files/psfex_conf.psfex")# -VERBOSE_TYPE")# FULL")
            # sys.exit(1)
            self.files_to_clean.append(self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat")

            self.psf_ref_image_name=self.path+f'out/proto_prepsfex_{self.rand_nums_string}.fits'
            self.files_to_clean.append(self.path+f'out/proto_prepsfex_{self.rand_nums_string}.fits')
            self.psf_ref_image = fits.open(self.psf_ref_image_name)

            self.hdu_psf_model_ref= fits.open(self.path+f'out/prepsfex_{self.rand_nums_string}.psf')
            self.files_to_clean.append(self.path+f'out/prepsfex_{self.rand_nums_string}.psf')
            self.chi_sq_psf_ref=self.hdu_psf_model_ref[1].header['CHI2']

            if self.termoutp!='quiet':
                # print(colored('Reduced Chi^2 of science image PSF fit:','green'),"%.1f" % self.chi_sq_psf_ref)
                print(info_g+' Reduced Chi^2 of science image PSF fit:',self.chi_sq_psf_ref)
                if self.chi_sq_psf_ref>3:
                    print(colored('Warning: PSF ref model may not be accurate','green'))
                    print(warn_y+' Warning: PSF ref model may not be accurate')

            #Convolve the reference image with a Gaussian kernel
            self.kernel_ref = self.psf_ref_image[0].data[0]

            # Read the Science image and convolve it with the Gaussian kernel
            self.sci_image_aligned=fits.open(self.sci_ali_name)
            self.sci_conv = scipy_convolve(self.sci_image_aligned[0].data,self.kernel_ref, mode='same', method='fft')
            self.sci_conv=np.nan_to_num(self.sci_conv)

            fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)

            # print(self.sci_conv_name)
            # print(self.ref_conv_name)
            self.sci_conv_hdu,self.ref_conv_hdu,self.ref_ali_hdu,self.sci_ali_hdu=fits.open(self.sci_conv_name)[0],fits.open(self.ref_conv_name)[0],fits.open(self.ref_ali_name)[0],fits.open(self.sci_ali_name)[0]
            self.sci_conv,self.ref_conv,self.ref_ali,self.sci_ali=self.sci_conv_hdu.data,self.ref_conv_hdu.data,self.ref_ali_hdu.data,self.sci_ali_hdu.data

            # print(self.sci_conv==self.ref_conv)
            # sys.exit(1)
            if self.args.show_plots:
                #show the psf of the science image and the reference image
                # norm = simple_norm(self.ref_ali_epsf.data, 'log', percent=99.)
                #     plt.imshow(self.ref_ali_epsf.data, norm=norm, origin='lower', cmap='viridis')
                #     plt.title(f"Reference PSF")
                #     plt.colorbar()

                # psf_fig = plt.figure(figsize=(10,10))
                # plt.subplot(1,2,1)
                # self.sci_norm = simple_norm(self.kernel_sci, 'log', percent=99.)
                # plt.imshow(self.kernel_sci, norm=self.sci_norm, origin='lower', cmap='viridis')
                # plt.title('Science PSF')
                # plt.subplot(1,2,2)
                # self.ref_norm = simple_norm(self.kernel_ref, 'log', percent=99.)
                # plt.imshow(self.kernel_ref, norm=self.ref_norm, origin='lower', cmap='viridis')
                # plt.title('Reference PSF')
                # plt.colorbar()
                # plt.show()
                # plt.close(psf_fig)

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
        if self.termoutp!='quiet':
            print(info_g+f" Beginning cross convolutions of PSFs and images")
            
        self.sci_conv_name=self.path+"convolved_sci/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'sci_convolved.fits'
        self.ref_conv_name=self.path+"convolved_ref/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][:-13]+'_'+str(datetime.timedelta(hours=int(self.sci_img_hdu.header[self.DATE_kw][11:13]), minutes=int(self.sci_img_hdu.header[self.DATE_kw][14:16]), seconds=float(self.sci_img_hdu.header[self.DATE_kw][17:21])).seconds)+'ref_convolved.fits'
        
        #########################################
        #BUILDING SCIENCE PSF
        self.sci_ali_hdu = fits.open(self.sci_ali_name)[0]
        self.sci_ali_mean,self.sci_ali_median, self.sci_ali_std = sigma_clipped_stats(self.sci_ali_hdu.data, sigma=5.0)
        self.sci_ali_iraffind= IRAFStarFinder(threshold=abs(starscale*self.sci_ali_std)*3,fwhm=3,roundhi=0.3,minsep_fwhm=1.5,peakmax=45000)
        self.sci_ali_co = self.sci_ali_iraffind(self.sci_ali_hdu.data)
        if self.termoutp!='quiet':
            print(info_g+' Measuring PSF of science image...')
            print(info_g+f" {len(self.sci_ali_co)} stars detected in aligned science image")
        self.sci_ali_cox,self.sci_ali_coy = self.sci_ali_co['xcentroid'],self.sci_ali_co['ycentroid']


        self.sci_ali_sig_clip = SigmaClip(sigma=3.)
        self.sci_ali_bkg_estimator = SExtractorBackground(self.sci_ali_sig_clip)        
        self.sci_ali_bkg = Background2D(self.sci_ali_hdu.data, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=self.sci_ali_bkg_estimator)
        self.sci_ali_bkg_error = self.sci_ali_bkg.background_rms

        self.sci_ali_err = calc_total_error(self.sci_ali_hdu.data,self.sci_ali_bkg_error,self.sci_gain)
        self.ali_pos = [(self.sci_ali_cox[i],self.sci_ali_coy[i]) for i in range(len(self.sci_ali_coy))]
        self.sci_ali_photaps = photutils.CircularAperture(self.ali_pos, r=20)
        self.sci_ali_photTab = photutils.aperture_photometry(self.sci_ali_hdu.data, self.sci_ali_photaps, error=self.sci_ali_err)
        self.sci_ali_photTab['SNR'] = self.sci_ali_photTab['aperture_sum']/self.sci_ali_photTab['aperture_sum_err']
        #remove any where aperture_sum is nan,negative, or 0 and where SNR is negative or 0
        # self.ref_ali_photTab = self.ref_ali_photTab[(self.ref_ali_photTab['SNR']<1000)&(self.ref_ali_photTab['SNR']>0)&(self.ref_ali_photTab['aperture_sum']>0)&(self.ref_ali_photTab['aperture_sum_err']>0)]
        # self.ref_ali_cox,self.ref_ali_coy = self.ref_ali_photTab['xcenter'],self.ref_ali_photTab['ycenter']
        self.sci_ali_photTab = self.sci_ali_photTab[(self.sci_ali_photTab['aperture_sum']>0)&(self.sci_ali_photTab['aperture_sum_err']>0)&(self.sci_ali_photTab['SNR']>0)&(self.sci_ali_photTab['SNR']<4000)&(self.sci_ali_photTab['aperture_sum']<45000)]
        self.sci_ali_cox,self.sci_ali_coy = self.sci_ali_photTab['xcenter'],self.sci_ali_photTab['ycenter']
        self.sci_ali_goodStars = (np.isnan(self.sci_ali_photTab['aperture_sum'])==False) & (self.sci_ali_photTab['aperture_sum']>0) & (self.sci_ali_photTab['aperture_sum_err']>0)
        self.sci_ali_goodStars = (self.sci_ali_photTab['aperture_sum']>80*self.sci_ali_photTab['aperture_sum_err'])
        self.sci_ali_seq_SNR = 40 #intial SNR threshold
        
        #lowering SNR if intital threshold leaves less than 10 stars
        # print(self.sci_ali_photTab['aperture_sum','aperture_sum_err','xcenter','ycenter'][self.sci_ali_goodStars])
        if len(self.sci_ali_goodStars[self.sci_ali_goodStars]) < self.sci_ali_seq_SNR:
            print(info_g+' Less than 20 stars with SNR >=',self.sci_ali_seq_SNR,'Lowering SNR threshold to',self.sci_ali_seq_SNR/2)
            self.sci_ali_goodStars = (self.sci_ali_photTab['aperture_sum']>self.sci_ali_seq_SNR*self.sci_ali_photTab['aperture_sum_err'])
            self.sci_ali_seq_SNR/=2
            if len(self.sci_ali_goodStars[self.sci_ali_goodStars]) < self.sci_ali_seq_SNR:
                print(info_g+' Less than 20 stars with SNR >=',self.sci_ali_seq_SNR,'Lowering SNR threshold to',self.sci_ali_seq_SNR/2)
                self.sci_ali_goodStars = (self.sci_ali_photTab['aperture_sum']>self.sci_ali_seq_SNR*self.sci_ali_photTab['aperture_sum_err'])
                self.sci_ali_seq_SNR/=2


        print(info_g+f" Found {len(self.sci_ali_goodStars)} stars in science with SNR>{self.sci_ali_seq_SNR}")
        self.sci_ali_psf_built = False
        self.sci_ali_nddata = astropy.nddata.NDData(data=self.sci_ali_hdu.data)
        self.sci_ali_psfinput = astropy.table.Table()
        self.sci_ali_psfinput['x'],self.sci_ali_psfinput['y'] = self.sci_ali_cox,self.sci_ali_coy

        self.sci_ali_psfthresh=5 #number of stars needed for fit
        if self.sci_filt=='u':
            self.sci_ali_psfthresh=2

        # psf quality tests
        self.sci_ali_sumpsf,self.sci_ali_minpsf,self.sci_ali_x_peak,self.sci_ali_y_peak,self.sci_ali_psf_iter,self.sci_ali_stamprad = -1,-1,0,0,0,15
        
        while (self.sci_ali_sumpsf < 0) or (self.sci_ali_minpsf < -0.01) or ((self.sci_ali_x_peak/len(self.sci_ali_psf) < 0.4) or (self.sci_ali_x_peak/len(self.sci_ali_psf) > 0.6)) or ((self.sci_ali_y_peak/len(self.sci_ali_psf) < 0.4) or (self.sci_ali_y_peak/len(self.sci_ali_psf) > 0.6)) and (self.sci_ali_psf_iter < 5):
            if self.sci_ali_psf_iter > 0:
                print(warn_y+f' #{self.sci_ali_psf_iter}: PSF failed quality checked, randomly varying parameters and trying again')
                self.sci_ali_stamprad += np.random.randint(11)-5
                self.sci_ali_stamprad = max([self.sci_ali_stamprad,10])
                self.sci_ali_psfthresh += np.random.randint(10)

            #extract stars from image
            self.sci_ali_psfstars = photutils.psf.extract_stars(self.sci_ali_nddata, self.sci_ali_psfinput[self.sci_ali_photTab['aperture_sum']>self.sci_ali_psfthresh*self.sci_ali_photTab['aperture_sum_err']], size=2*self.sci_ali_stamprad+5)

            while(len(self.sci_ali_psfstars))<5 and self.sci_ali_psfthresh>0:
                print(warn_y+f' #{self.sci_ali_psf_iter}: Warning: too few PSF stars with threshold '+str(self.sci_ali_psfthresh)+' sigma, trying lower sigma')
                self.sci_ali_psfthresh -= 1
                self.sci_ali_psfstars = photutils.psf.extract_stars(self.sci_ali_nddata, self.sci_ali_psfinput[self.sci_ali_photTab['aperture_sum']>self.sci_ali_psfthresh*self.sci_ali_photTab['aperture_sum_err']], size=2*self.sci_ali_stamprad+5)
            if len(self.sci_ali_psfstars)<5:
                self.sci_ali_psfthresh = 5
                print(warn_y+f' #{self.sci_ali_psf_iter}: Could not find 5 PSF stars, trying for 3...')
                while(len(self.sci_ali_psfstars))<3 and self.sci_ali_psfthresh>0:
                    self.sci_ali_psfthresh -= 1
                    self.sci_ali_psfstars = photutils.psf.extract_stars(self.sci_ali_nddata, self.sci_ali_psfinput[self.sci_ali_photTab['aperture_sum']>self.sci_ali_psfthresh*self.sci_ali_photTab['aperture_sum_err']], size=2*self.sci_ali_stamprad+5)
                if self.sci_ali_psfthresh < 3:
                    self.sci_ali_empirical = False
                    break


            try:
                self.sci_ali_epsf_builder = photutils.EPSFBuilder(maxiters=10,recentering_maxiters=10,oversampling=1,smoothing_kernel='quadratic',progress_bar=progress_bar,shape=2*self.sci_ali_stamprad-1)
                self.sci_ali_psf_iter+=1
                self.sci_ali_epsf, self.sci_ali_fitted_stars = self.sci_ali_epsf_builder(self.sci_ali_psfstars) #science PSF
                self.sci_ali_psf = self.sci_ali_epsf.data
                self.sci_ali_x_peak,self.sci_ali_y_peak = np.where(self.sci_ali_psf==self.sci_ali_psf.max())[1][0],np.where(self.sci_ali_psf==self.sci_ali_psf.max())[0][0]
                self.sci_ali_minpsf,self.sci_ali_sumpsf = np.min(self.sci_ali_psf),np.sum(self.sci_ali_psf)
                self.sci_ali_psf_built=True

                

                if self.sci_ali_psf_iter>20:
                    print(warn_r+f" #{self.sci_ali_psf_iter}: Error building PSF from science image, error encountered iteration{self.sci_ali_psf_iter}>=max iterations:{e}")
                    print(warn_r+f" Exiting...")
                    self.sys_exit=True
                    return


            except Exception as e:
                print(warn_y+f" Building PSF failed, error encountered {e}")
            #     if self.output_to_log!=False:
            #         self.log_error_4_1 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Unable to build science PSF from science image '+'\n'
            #         self.night_log = open(self.log_name,"r+")
            #         if self.log_error_4_1 in self.night_log.readlines():
            #             pass
            #         else:
            #             self.night_log.write(self.log_error_4_1)
            #             self.night_log.write('\n')
            #             self.night_log.close()
            #     if self.termoutp!='quiet':
            #         print(warn_r+f" #{self.sci_ali_psf_iter}: Error building PSF from science image, error encountered :{e}")
            #         # print(f"Exiting...")

        
        if self.sci_ali_psf_built == True:
            print(info_g+f" Science PSF built successfully")
            if plot_psfs!=False or self.args.show_plots==True:
                sci_psf_fig = plt.figure(figsize=(12,8))
                norm = simple_norm(self.sci_ali_epsf.data, 'log', percent=99.)
                plt.imshow(self.sci_ali_epsf.data, norm=norm, origin='lower', cmap='viridis')
                plt.title(f"Science PSF")
                plt.colorbar()
                plt.show()
                plt.close(sci_psf_fig)

        if self.sci_ali_psf_built==False:
            self.log_error_4_1 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Unable to build science PSF from science image, <3 bright stars to build PSF '+'\n'
            self.night_log = open(self.log_name,"r+")
            if self.log_error_4_1 in self.night_log.readlines():
                pass
            else:
                self.night_log.write(self.log_error_4_1)
                self.night_log.write('\n')
                self.night_log.close()
            print(warn_r+f" Unable to build science PSF from science image, <3 bright stars to build PSF")
            print(warn_r+f" Exiting...")
            self.sys_exit=True
            return
        else:
            self.ref_image_aligned=fits.open(self.ref_ali_name)
            self.ref_ali_img_hdu = self.ref_image_aligned[0]
            self.ref_conv = scipy_convolve(self.ref_ali_img_hdu.data, self.sci_ali_psf, mode='same', method='fft')
            fits.writeto(self.ref_conv_name, data=self.ref_conv, header=self.ref_image_aligned[0].header,overwrite=True)

            if self.termoutp!='quiet':
                print(info_g+f" Convolved science PSF with reference image with size {self.ref_conv.shape}")
                print(info_g+f" Convolved reference image saved to {self.ref_conv_name}")

            
                
            

            if self.to_subtract!=False:
                #########################################
                #BUILDING REFERENCE PSF
                self.ref_ali_hdu = fits.open(self.ref_ali_name)[0]
                self.ref_ali_mean,self.ref_ali_median, self.ref_ali_std = sigma_clipped_stats(self.ref_ali_hdu.data, sigma=5.0)
                self.ref_ali_iraffind= IRAFStarFinder(threshold=abs(starscale*self.ref_ali_std)*3,fwhm=3,roundhi=0.3,minsep_fwhm=1.5,peakmax=50000)
                self.ref_ali_co = self.ref_ali_iraffind(self.ref_ali_hdu.data)
                if self.termoutp!='quiet':
                    print(info_g+f" {len(self.ref_ali_co)} stars detected in aligned refrence image")
                    



                self.ref_ali_cox,self.ref_ali_coy = self.ref_ali_co['xcentroid'],self.ref_ali_co['ycentroid']


                self.ref_ali_sig_clip = SigmaClip(sigma=3.)
                self.ref_ali_bkg_estimator = SExtractorBackground(self.ref_ali_sig_clip)        
                self.ref_ali_bkg = Background2D(self.ref_ali_hdu.data, (150, 150), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=self.ref_ali_bkg_estimator)
                self.ref_ali_bkg_error = self.ref_ali_bkg.background_rms

                self.ref_ali_err = calc_total_error(self.ref_ali_hdu.data,self.ref_ali_bkg_error,self.sci_gain)
                self.ali_pos = [(self.ref_ali_cox[i],self.ref_ali_coy[i]) for i in range(len(self.ref_ali_coy))]
                self.ref_ali_photaps = photutils.CircularAperture(self.ali_pos, r=20)
                self.ref_ali_photTab = photutils.aperture_photometry(self.ref_ali_hdu.data, self.ref_ali_photaps, error=self.ref_ali_err)
                self.ref_ali_photTab['SNR'] = self.ref_ali_photTab['aperture_sum']/self.ref_ali_photTab['aperture_sum_err']
                #remove saturated stars, i.e. stars with SNR>1000, SNR<0, aperture_sum<0, aperture_sum_err<0
                self.ref_ali_photTab = self.ref_ali_photTab[(self.ref_ali_photTab['SNR']<1000)&(self.ref_ali_photTab['SNR']>0)&(self.ref_ali_photTab['aperture_sum']>0)&(self.ref_ali_photTab['aperture_sum_err']>0)&(self.ref_ali_photTab['aperture_sum']<50000)]
                # print(np.sort(self.ref_ali_photTab)[0:15])
                # sys.exit()
                self.ref_ali_cox,self.ref_ali_coy = self.ref_ali_photTab['xcenter'],self.ref_ali_photTab['ycenter']
                #check in the aperture of the stars for any NaN values or negative values
                # for k in range(len(self.ref_ali_photTab)):
                #     if any(np.isnan(x) for x in self.ref_ali_hdu.data[int(self.ref_ali_photTab['ycenter'][k])-20:int(self.ref_ali_photTab['ycenter'][k])+20,
                #                                                         int(self.ref_ali_photTab['xcenter'][k])-20:int(self.ref_ali_photTab['xcenter'][k])+20]) == True:
                #         print('saturated star')
                # print(self.ref_ali_photTab['aperture_sum','aperture_sum_err','SNR','xcenter','ycenter'])
                self.ref_ali_goodStars = (self.ref_ali_photTab['aperture_sum']>80*self.ref_ali_photTab['aperture_sum_err'])
                if self.special_case!=None:
                    if any(phrase == self.special_case for phrase in ['SN2023ixf','bright','brightstars.cat','sn2023ixf','2023ixf','brightstars']):
                        self.bright_stars_sc = ascii.read(self.path+'calstars4e.cat')#,skip_lines=1,names=['ra','dec','B','g','V','r','R','i','I','z','y'])
                        print(info_g+f" Using bright stars catalog")
                        self.bright_stars_sc_wcs = SkyCoord(ra=self.bright_stars_sc['ra']*u.deg, dec=self.bright_stars_sc['dec']*u.deg,frame='fk5')
                        print(self.ref_ali_name)
                        self.bs_coords_pix=wcs_to_pixels(self.ref_ali_name,np.column_stack((self.bright_stars_sc['ra'],self.bright_stars_sc['dec'])))



                        # plt.show()

                        self.bright_stars_sc['xcentroid']=self.bs_coords_pix[:,0]
                        self.bright_stars_sc['ycentroid']=self.bs_coords_pix[:,1]

                        #make an empty table with the same columns as the catalog

                        self.bright_stars_sc_new = Table(names=self.bright_stars_sc.colnames,dtype=[float for i in range(len(self.bright_stars_sc.colnames))])
                        #append the data from the catalog to the empty table

                        for i in range(len(self.bright_stars_sc)-1):
                            row = []
                            for j in range(len(self.bright_stars_sc[i])):
                                if self.bright_stars_sc[i][j]!='-':
                                    row.append(float(self.bright_stars_sc[self.bright_stars_sc.colnames[j]].data[i]))
                                else:
                                    row.append(None)
                            self.bright_stars_sc_new.add_row(row)


                        self.bright_stars_sc=self.bright_stars_sc_new
                        print(self.bright_stars_sc)
                        self.bright_stars_sc = self.bright_stars_sc[self.bright_stars_sc[self.sci_filt]<self.bright_cat_mag_lim]
                        self.ref_aligned_wcs = WCS(self.ref_ali_hdu.header)
                        self.ref_ali_photTab_coords = self.ref_aligned_wcs.all_pix2world(np.column_stack((self.ref_ali_photTab['xcenter'],self.ref_ali_photTab['ycenter'])),1)
                        self.ref_ali_photTab_coords = SkyCoord(ra=self.ref_ali_photTab_coords[:,0]*u.deg, dec=self.ref_ali_photTab_coords[:,1]*u.deg,frame='fk5')
                        self.ref_ali_photTab['ra'],self.ref_ali_photTab['dec'] = self.ref_ali_photTab_coords.ra,self.ref_ali_photTab_coords.dec
                        # self.bright_stars_sc['ra'],self.bright_stars_sc['dec'] = self.bright_stars_sc_coords.ra,self.bright_stars_sc_coords.dec 
                        # print(self.bright_stars_sc)
                        # print(self.ref_ali_photTab)
                        # sys.exit()


                        self.indx, self.d2d, self.d3d =self.ref_ali_photTab_coords.match_to_catalog_sky(self.bright_stars_sc_wcs)
                        #where stars match by search rad in arcseconds!
                        self.upd_indx=np.where(self.d2d<search_rad/3600.*u.deg)[0]
                        # print(self.indx[self.upd_indx])
                        print(info_g+' Bright stars found=',len(self.upd_indx),',',round(3600*(self.ref_width/60),6), 'arcsec search radius')
                        self.ref_ali_photTab=self.ref_ali_photTab[self.upd_indx]
                        self.ref_ali_cox,self.ref_ali_coy = self.ref_ali_photTab['xcenter'],self.ref_ali_photTab['ycenter']



                        #find the common stars between the catalog and the detected stars
                        # print(int(self.ref_ali_cox.value),type(self.ref_ali_cox))
                        # print(int(self.ref_ali_coy.value),type(self.ref_ali_coy))
                        #convert self.ref_ali_cox and self.ref_ali_coy to arrays of integers
                        self.ref_ali_cox = np.array([int(i.value) for i in self.ref_ali_cox])
                        self.ref_ali_coy = np.array([int(i.value) for i in self.ref_ali_coy])

                        self.bright_stars_sc['xcentroid']=np.array([int(i) for i in self.bright_stars_sc['xcentroid']])
                        self.bright_stars_sc['ycentroid']=np.array([int(i) for i in self.bright_stars_sc['ycentroid']])

                        print(self.ref_ali_cox,self.ref_ali_coy)
                        # print(self.bright_stars_sc['xcentroid'],self.bright_stars_sc['ycentroid'])
                        # print(len(self.ref_ali_cox),len(self.ref_ali_coy))

                        # print(int(self.bright_stars_sc['xcentroid']),type(self.bright_stars_sc['xcentroid']))
                        # print(int(self.bright_stars_sc['ycentroid']),type(self.bright_stars_sc['ycentroid']))
                        # print
                        # self.ref_ali_goodStars = np.array([i for i in range(len(self.ref_ali_cox)-1) if any(np.sqrt((self.ref_ali_cox[i]-self.bright_stars_sc['xcentroid'])**2+(self.ref_ali_coy[i]-self.bright_stars_sc['ycentroid'])**2)<12)])
                        # print([(np.sqrt((self.ref_ali_cox[i]-self.bright_stars_sc['xcentroid'])**2+(self.ref_ali_coy[i]-self.bright_stars_sc['ycentroid'])**2),self.ref_ali_photTab['xcenter','ycenter'][i]) for i in range(len(self.ref_ali_cox)) if any(np.sqrt((self.ref_ali_cox[i]-self.bright_stars_sc['xcentroid'])**2+(self.ref_ali_coy[i]-self.bright_stars_sc['ycentroid'])**2)<100)])
                        self.ref_ali_seq_SNR = 20 #intial SNR threshold
                        print(self.ref_ali_goodStars)
                        sys.exit()


                # sys.exit()
                # if len(self.ref_ali_goodStars[self.ref_ali_goodStars]) >1e3: #choose the top 500 stars with the highest SNR to speed upo the psf building process
                    # self.ref_ali_goodStars = np.argsort(self.ref_ali_photTab['aperture_sum']>10*self.ref_ali_photTab['aperture_sum_err'])[::-1][:500]
                
                # lowering SNR if intital threshold leaves less than 10 stars
                if self.special_case==None:
                    self.ref_ali_seq_SNR = 80 #intial SNR threshold
                    if len(self.ref_ali_goodStars[self.ref_ali_goodStars]) < self.ref_ali_seq_SNR:
                        print(info_g+f" Lowering SNR threshold to {self.ref_ali_seq_SNR/2}")
                        self.ref_ali_goodStars = (self.ref_ali_photTab['aperture_sum']>5*self.ref_ali_photTab['aperture_sum_err'])
                        self.ref_ali_seq_SNR/=2
    
                # print(self.ref_ali_photTab)
                if self.termoutp!='quiet':
                    print(info_g+f" Found {len(self.ref_ali_goodStars)} stars in refrence with SNR>{self.ref_ali_seq_SNR}")
                # if len(self.ref_ali_photTab)>250:
                #     self.ref_ali_photTab_sort = self.ref_ali_photTab.sort('SNR')

                # self.ref_ali_photTab = vstack([row for ind,row in enumerate(self.ref_ali_photTab_sort)])
                # print(self.ref_ali_photTab)
                # sys.exit()
                self.ref_ali_nddata = astropy.nddata.NDData(data=self.ref_ali_hdu.data)
                self.ref_ali_psfinput = astropy.table.Table()
                self.ref_ali_psfinput['x'],self.ref_ali_psfinput['y'] = self.ref_ali_cox,self.ref_ali_coy
                # print(len(self.ref_ali_goodStars))
                # sys.exit()

                self.ref_ali_psfthresh=5 #number of stars needed for fit
                self.ref_psfthresh_tried = []
                # if self.sci_filt=='u':
                #     self.ref_ali_psfthresh=2

                # psf quality tests
                self.ref_ali_sumpsf,self.ref_ali_minpsf,self.ref_ali_x_peak,self.ref_ali_y_peak,self.ref_ali_psf_iter,self.ref_ali_stamprad = -1,-1,0,0,0,15
                self.ref_stamprad_tried = []
                while (self.ref_ali_sumpsf < 0) or (self.ref_ali_minpsf < -0.01) or ((self.ref_ali_x_peak/len(self.ref_ali_psf) < 0.4) or (self.ref_ali_x_peak/len(self.ref_ali_psf) > 0.6)) or ((self.ref_ali_y_peak/len(self.ref_ali_psf) < 0.4) or (self.ref_ali_y_peak/len(self.ref_ali_psf) > 0.6)) and (self.ref_ali_psf_iter < 5):
                    if self.ref_ali_psf_iter > 0:
                        print(warn_y+f' #{self.ref_ali_psf_iter}: PSF failed quality checked, randomly varying parameters and trying again')
                        # while self.ref_ali_stamprad in self.ref_stamprad_tried:
                        self.ref_ali_stamprad += np.random.randint(11)-5
                        self.ref_ali_stamprad = max([self.ref_ali_stamprad,10])
                        # while self.ref_ali_psfthresh in self.ref_psfthresh_tried:
                        self.ref_ali_psfthresh += np.random.randint(10)

                    #extract stars from image
                    self.ref_ali_psfstars = photutils.psf.extract_stars(self.ref_ali_nddata, self.ref_ali_psfinput[self.ref_ali_photTab['aperture_sum']>self.ref_ali_psfthresh*self.ref_ali_photTab['aperture_sum_err']], size=2*self.ref_ali_stamprad+5)
                    self.ref_psfthresh_tried.append(self.ref_ali_psfthresh),self.ref_stamprad_tried.append(self.ref_ali_stamprad)
                    while(len(self.ref_ali_psfstars))<5 and self.ref_ali_psfthresh>0:
                        print(warn_y+f' #{self.ref_ali_psf_iter}: Warning: too few PSF stars with threshold '+str(self.ref_ali_psfthresh)+' sigma, trying lower sigma')
                        self.ref_ali_psfthresh -= 1
                        self.ref_ali_psfstars = photutils.psf.extract_stars(self.ref_ali_nddata, self.ref_ali_psfinput[self.ref_ali_photTab['aperture_sum']>self.ref_ali_psfthresh*self.ref_ali_photTab['aperture_sum_err']], size=2*self.ref_ali_stamprad+5)
                        # self.ref_psf_thresh_tried.append(self.ref_ali_psfthresh)
                    if len(self.ref_ali_psfstars)<5:
                        self.ref_ali_psfthresh = 5
                        print(warn_y+f' #{self.ref_ali_psf_iter}: Could not find 5 PSF stars, trying for 3...')
                        while(len(self.ref_ali_psfstars))<3 and self.ref_ali_psfthresh>0:
                            self.ref_ali_psfthresh -= 1
                            self.ref_ali_psfstars = photutils.psf.extract_stars(self.ref_ali_nddata, self.ref_ali_psfinput[self.ref_ali_photTab['aperture_sum']>self.ref_ali_psfthresh*self.ref_ali_photTab['aperture_sum_err']], size=2*self.ref_ali_stamprad+5)
                        if self.ref_ali_psfthresh < 3:
                            self.ref_ali_empirical = False
                            break


                    try:
                        self.ref_ali_epsf_builder = photutils.EPSFBuilder(maxiters=15,recentering_maxiters=15,oversampling=1,smoothing_kernel='quadratic',progress_bar=progress_bar,shape=2*self.ref_ali_stamprad-1)
                        self.ref_ali_psf_iter+=1
                        self.ref_ali_epsf, self.ref_ali_fitted_stars = self.ref_ali_epsf_builder(self.ref_ali_psfstars) #science PSF
                        self.ref_ali_psf = self.ref_ali_epsf.data
                        self.ref_ali_x_peak,self.ref_ali_y_peak = np.where(self.ref_ali_psf==self.ref_ali_psf.max())[1][0],np.where(self.ref_ali_psf==self.ref_ali_psf.max())[0][0]
                        self.ref_ali_minpsf,self.ref_ali_sumpsf = np.min(self.ref_ali_psf),np.sum(self.ref_ali_psf)

                        # print(self.ref_ali_sumpsf,self.ref_ali_minpsf,self.ref_ali_x_peak,self.ref_ali_y_peak,self.ref_ali_psf_iter,self.ref_ali_stamprad)
                        # print('sumpsf:',self.ref_ali_sumpsf)
                        # print('minpsf:',self.ref_ali_minpsf)
                        # print('x_peak:',self.ref_ali_x_peak)
                        # print('y_peak:',self.ref_ali_y_peak)
                        # print('psf_iter:',self.ref_ali_psf_iter)
                        # print('stamprad:',self.ref_ali_stamprad)
                        # print(self.ref_ali_psf)
                        

                        


                    except Exception as e:
                        print(warn_r+f" #{self.ref_ali_psf_iter}: Building reference PSF failed, error encountered {e}")
                        if self.output_to_log!=False:
                            self.log_error_4_2 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Unable to build reference PSF from reference image '+'\n'
                            self.night_log = open(self.log_name,"r+")
                            if self.log_error_4_2 in self.night_log.readlines():
                                pass
                            else:
                                self.night_log.write(self.log_error_4_2)
                                self.night_log.write('\n')
                                self.night_log.close()
                        if self.termoutp!='quiet':
                            print(warn_r+f" #{self.ref_ali_psf_iter}: Error building PSF from reference image, error encountered iteration{self.ref_ali_psf_iter}:{e}")
                            # print(f"Exiting...")
                        
                        if self.ref_ali_psf_iter>40:
                            print(warn_r+f" #{self.ref_ali_psf_iter}: Error building PSF from reference image, error encountered iteration{self.ref_ali_psf_iter}>=max iterations:{e}")
                            print(warn_r+f" Exiting...")
                            self.sys_exit=True
                            return

                print(info_g+f" Reference PSF built succesfully")
                plot_psfs=True
                if plot_psfs!=False  or self.args.show_plots==True:
                    ref_psf_fig = plt.figure(figsize=(12,8))
                    norm = simple_norm(self.ref_ali_epsf.data, 'log', percent=99.)
                    plt.imshow(self.ref_ali_epsf.data, norm=norm, origin='lower', cmap='viridis')
                    plt.title(f"Reference PSF")
                    plt.colorbar()
                    plt.show()
                    plt.close(ref_psf_fig)
                # bar.update(1)

                self.sci_image_aligned=fits.open(self.sci_ali_name)
                self.sci_conv = scipy_convolve(self.sci_image_aligned[0].data, self.ref_ali_psf, mode='same', method='fft')
                self.sci_conv=np.nan_to_num(self.sci_conv)
                fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)

                if self.termoutp!='quiet':
                    print(info_g+f" Convolved reference PSF with science image with size {self.sci_conv.shape}")
                    print(info_g+f" Convolved science image saved to {self.sci_conv_name}")

                

                self.kernel_sci,self.kernel_ref=self.sci_ali_psf,self.ref_ali_psf
                # print(self.ref_conv)
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
        print(info_g+f" Reference convolved image size: {self.ref_conv.shape}")
        print(info_g+f" Science convolved image size: {self.sci_conv.shape}")
        # if self.args.show_plots==True:
            # sci_conv_fig = plt.figure(figsize=(12,8))
            # vmin,vmax=visualization.ZScaleInterval().get_limits(self.sci_conv[~np.isnan(self.sci_conv)])
            # plt.imshow(self.sci_conv, origin='lower', cmap='viridis',vmin=vmin,vmax=vmax)
            # plt.title(f"Convolved science image")
            # plt.colorbar()
            # # plt.show()

            # ref_conv_fig = plt.figure(figsize=(12,8))
            # vmin,vmax=visualization.ZScaleInterval().get_limits(self.ref_conv)
            # plt.imshow(self.ref_conv, origin='lower', cmap='viridis',vmin=vmin,vmax=vmax)
            # plt.title(f"Convolved reference image")
            # plt.colorbar()
            # plt.show()
            # plt.close(sci_conv_fig),plt.close(ref_conv_fig)

        if self.termoutp!='quiet':
            print(info_g+f" Generating reference catalogs")
        self.ref_width=self.sci_img_hdu.header['NAXIS2']*self.sci_ps/60

        if self.auto_cat=="auto":
            if self.use_sdss==False and (self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z'):
                self.ref_cat=panstarrs_query(ra_deg=round(self.sci_c.ra.deg,6),dec_deg=round(self.sci_c.dec.deg,6), rad_deg=round(self.ref_width/60.,6))

                self.stars=np.where((self.ref_cat['iMeanPSFMag']-self.ref_cat['iMeanKronMag']<0.05) & (self.ref_cat[str(self.sci_filt)+'MeanPSFMagErr']<0.06) & (self.ref_cat[str(self.sci_filt)+'MeanPSFMag']<23.5)  & (self.ref_cat[str(self.sci_filt)+'MeanPSFMag']!=-999.))[0]
                self.ref_cat=np.array(self.ref_cat[self.stars])
                if len(self.stars)==0:
                    if self.output_to_log!=False:
                        self.log_error_3_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Length of PS1 reference catalog is 0 '+'\n'
                        self.night_log = open(self.log_name,"r+")
                        self.night_log.write('\n')
                        if self.log_error_3_5 in self.night_log.readlines():
                            pass
                        else:
                            self.night_log.write(self.log_error_3_5)
                            self.night_log.close()
                    
                    self.sys_exit=True
                    return
                self.ref_coords_wcs_sky = SkyCoord(ra=self.ref_cat['raMean']*u.deg, dec=self.ref_cat['decMean']*u.deg,frame='fk5')
                self.ref_coords_wcs=np.column_stack((self.ref_cat['raMean'],self.ref_cat['decMean']))
                self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)
                # print(self.ref_coords_wcs)

            if self.sci_filt=='u' or self.use_sdss==True:
                self.ref_cat=sdss_query(ra_deg=round(self.sci_c.ra.deg,6),dec_deg=round(self.sci_c.dec.deg,6), rad_deg=round((self.ref_width)*3,6))
                # print(self.ref_cat)
                self.stars=self.ref_cat

                if len(self.stars)==0:
                    if self.output_to_log!=False:
                        self.log_error_3_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Length of SDSS reference catalog is 0 '+'\n'
                        self.night_log = open(self.log_name,"r+")
                        if self.log_error_3_5 in self.night_log.readlines():
                            pass
                        else:
                            self.night_log.write(self.log_error_3_5)
                            self.night_log.write('\n')

                            self.night_log.close()
                    self.sys_exit=True
                    return
                
                # print(len(self.stars),round(3600*(self.ref_width/60),6), 'arcsec search radius')
                self.ref_coords_wcs_sky = SkyCoord(ra=np.array(self.ref_cat['ra'])*u.deg, dec=np.array(self.ref_cat['dec'])*u.deg,frame='fk5')

                # print(len(self.ref_coords_wcs_sky))
                self.ref_coords_wcs=np.column_stack((self.ref_cat['ra'],self.ref_cat['dec']))
                # print(len(self.ref_coords_wcs))
                self.ref_coords_pix=wcs_to_pixels(self.ref_ali_name,self.ref_coords_wcs)
                # print(self.ref_coords_pix)
        else:
            #If a reference catalogue is passed in with it's full path as self.auto_cat
            if self.auto_cat.startswith(self.path):
                self.auto_cat = self.path+self.auto_cat

            if not os.path.exists(self.auto_cat):
                print(warn_r+f' Reference catalogue {self.auto_ref} not found')
                self.log_error_3_5_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+f': Unable to find reference catalogue from path specified {self.auto_cat} '+'\n'
                self.night_log = open(self.log_name,"r+")
                if self.log_error_3_5_5 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_3_5_5)
                    self.night_log.write('\n')
                    self.night_log.close()
                
                self.sys_exit=True
                return

            # print(pd.read_csv(self.auto_cat,skiprows=1,header=None))
            self.cat_arr=  np.array(ascii.read(self.auto_cat,data_start=1,names=["filt","mag","magerr","xpos","ypos","ra","dec"]))
            self.ref_cat = pd.DataFrame(self.cat_arr,columns=["filt","mag","magerr","xpos","ypos","ra","dec"]) #mag,magerr,(pixel)xpos,(pixel)ypos,RA,DEC
            self.stars=self.ref_cat
            # self.ref_coords_wcs_sky = np.column_stack((self.ref_cat["RA"],self.ref_cat["DEC"]))
            
            self.ref_coords_wcs_sky = SkyCoord(ra=self.ref_cat["ra"]*u.deg,dec=self.ref_cat["dec"]*u.deg,frame='fk5')
            self.ref_coords_wcs = self.ref_coords_wcs_sky
            self.ref_coords_pix = np.column_stack((self.ref_cat["xpos"],self.ref_cat["ypos"]))
            
 

        if self.termoutp!='quiet':
            print(info_g+' Catalog stars in PS1/SDSS found=',len(self.stars),',',round(3600*(self.ref_width/60),6), 'arcsec search radius')
        # print(self.ref_coords_pix)
        # print(self.sci_conv)
        self.mean,self.median, self.std = sigma_clipped_stats(self.sci_conv, sigma=3.0)
        print(info_g+' Mean, Median, Std of sci_conv:',round(self.mean,6),round(self.median,6),round(self.std,6))
        # print(self.sci_conv_name)
        # sys.exit()
        # self.mean,self.median, self.std = abs(self.mean),abs(self.median),abs(self.std)
        # if all([self.mean,self.median, self.std]) == 0:
        # self.mean,self.median, self.std = np.nanmean(self.sci_conv),np.nanmedian(self.sci_conv),np.nanstd(self.sci_conv)
        # print(self.sci_conv)
        #save self.sci_conv
        # fits.writeto(self.sci_conv_name, data=self.sci_conv, header=self.sci_image_aligned[0].header,overwrite=True)
        #######################################################
        # iraf star finder
        #threshold: The absolute image value above which to select sources.
        #fwhm: The full-width half-maximum (FWHM) of the 2D circular Gaussian kernel in units of pixels.
        #sigma_radius=1.5,minsep_fwhm=2.5,sharplo=0.5,sharphi=2.0,roundlo=0.0,roundhi=0.2, how round to find sources, default is 0.2 ... roundhi=1-(b/a)
        # if self.std*starscale<=3:
        print(info_g+' Detecting stars with SExtractor')
        # self.iraffind= IRAFStarFinder(threshold=(starscale*5)*self.std,fwhm=3.0)#,roundhi=0.15)
        # sex_command2 = sex_path+" "+self.sci_conv_name+" -c "+self.path+"config_files/align_sex.config -VERBOSE_TYPE QUIET -CATALOG_NAME "+self.path+f"config_files/prepsfex_{self.rand_nums_string}.cat -MAG_ZEROPOINT 25.0"
        # os.system(sex_command2)
        # cat = open(path+f"config_files/prepsfex_{self.rand_nums_string}.cat",'r')
        # catlines = cat.readlines()
        # cat.close()

        # print(catlines)
        self.source_objs = sextract(self.sci_conv_name, 0, 0, 3, 12, maxellip=0.7, saturation=-1)
        self.sources = np.array([(self.source_objs[i].x,self.source_objs[i].y) for i in range(len(self.source_objs))])
        if self.termoutp!='quiet':
            print(info_g+' Threshold for detecting stars:',(starscale*5)*self.std)
        self.star_coords_pix=np.column_stack((self.sources[:,0],self.sources[:,1]))
        # else:
        #     print(info_g+' Detecting stars with IRAFStarFinder')
        #     self.iraffind= IRAFStarFinder(threshold=abs(starscale*self.std),fwhm=3.0,roundhi=0.3)
        #     if self.termoutp!='quiet':
        #         print(info_g+' Threshold for detecting stars:',int(starscale*self.std))
        #     self.sources = self.iraffind(self.sci_conv[60:len(self.sci_conv)-60,60:len(self.sci_conv)-60] - self.median)
        
        
        #     self.star_coords_pix=np.column_stack((self.sources['xcentroid']+60.,self.sources['ycentroid']+60.))

        self.matched_catalog_mag=[]
        # try:
        # print(info_g+'Size of sci_conv: '+str(np.shape(self.sci_conv)))
        # if self.telescope in SEDM:
        #     self.iraffind= IRAFStarFinder(threshold=(starscale*5)*self.std,fwhm=3.0,roundhi=0.15)

        #     self.sources = self.iraffind(self.sci_conv)# - self.median)#[60:len(self.sci_conv)-60,60:len(self.sci_conv)-60])# - self.median)
        #     search_rad=1
        #     # self.sources = self.iraffind(self.sci_conv)# - self.median)
        # else:
        #     self.sources = self.iraffind(self.sci_conv[60:len(self.sci_conv)-60,60:len(self.sci_conv)-60] - self.median)



    
        
        self.star_coords_wcs=load_wcs_from_file(filename=self.sci_ali_name,coord=self.star_coords_pix)
        # print(self.star_coords_wcs)
        self.star_coords_wcs_sky = SkyCoord(ra=self.star_coords_wcs[:,0]*u.deg, dec=self.star_coords_wcs[:,1]*u.deg,frame='fk5')

        if self.special_case!=None:
            if any(phrase == self.special_case for phrase in ['SN2023ixf','bright','brightstars.cat','sn2023ixf','2023ixf','brightstars']):
                #find crossover between self.star_coords_wcs_sky and bright_stars_sc
                self.bright_stars_sc_wcs = SkyCoord(ra=self.bright_stars_sc['ra']*u.deg, dec=self.bright_stars_sc['dec']*u.deg,frame='fk5')

                self.indx, self.d2d, self.d3d =self.star_coords_wcs_sky.match_to_catalog_sky(self.bright_stars_sc_wcs)
                #where stars match by search rad in arcseconds!
                self.upd_indx=np.where(self.d2d<search_rad/3600.*u.deg)[0]
                # print(self.indx[self.upd_indx])
                print(info_g+' Bright stars found=',len(self.upd_indx),',',round(3600*(self.ref_width/60),6), 'arcsec search radius')
                self.star_coords_pix=self.star_coords_pix[self.upd_indx]
                self.star_coords_wcs=self.star_coords_wcs[self.upd_indx]
                self.star_coords_wcs_sky=self.star_coords_wcs_sky[self.upd_indx]

                self.matched_ref_coords_pix=self.bright_stars_sc['xcentroid','ycentroid'][self.indx[self.upd_indx]]
                self.matched_star_coords_pix=self.star_coords_pix

                self.matched_catalog_mag=self.bright_stars_sc[self.sci_filt][self.indx[self.upd_indx]]
                self.matched_catalog = self.bright_stars_sc[self.indx[self.upd_indx]]

                # for i in range(len(self.matched_ref_coords_pix)):
                #     print('sci',self.star_coords_wcs[i])
                #     print('ref',self.matched_catalog['ra','dec'][i])
                #     print(self.matched_star_coords_pix[i])
                #     print(self.matched_ref_coords_pix[i])
                #     print()
                # # print(self.bright_stars_sc[self.indx[self.upd_indx]])
                # sys.exit()

        else:
            #find crossover between panstarrs ad reference images
            print(info_g+' Searching for stars in reference catalog')
            # print(info_g+' Using '+self.auto_ref+' catalog')
            # print(self.star_coords_wcs_sky,self.ref_coords_wcs_sky)
            self.indx, self.d2d, self.d3d =self.star_coords_wcs_sky.match_to_catalog_sky(self.ref_coords_wcs_sky)
            #where stars match by search rad in arcseconds!
            self.upd_indx=np.where(self.d2d<7.5/3600.*u.deg)[0]
            # print(self.indx[self.upd_indx])
            print(info_g+' Catalog stars in PS1/SDSS found=',len(self.upd_indx),',',round(3600*(self.ref_width/60),6), 'arcsec search radius')

        

            self.matched_ref_coords_pix=self.ref_coords_pix[self.indx[self.upd_indx]]
            self.matched_star_coords_pix=self.star_coords_pix[self.upd_indx]

            # self.matched_fwhm=self.fwhm[self.upd_indx]

            
            # print(self.sci_filt=='r')
            # print(self.matched_ref_coords_pix,self.matched_star_coords_pix)
            if self.auto_cat == 'auto' and self.use_sdss==False and (self.sci_filt=='g' or self.sci_filt=='r' or self.sci_filt=='i' or self.sci_filt=='z'):
                self.matched_catalog=self.ref_cat[self.indx[self.upd_indx]]
                # print(self.matched_catalog)
                # sys.exit()
                self.matched_catalog_mag=self.matched_catalog[str(self.sci_filt)+'MeanPSFMag']

                # for i in range(len(self.matched_catalog_mag)):
                #     print(self.matched_catalog[i])
                #     print(self.bright_stars_sc[i])
                #     print()


            elif self.auto_cat == 'auto' and (self.sci_filt=='u' or self.use_sdss==True):
                self.string_band='psfMag_'+str(self.sci_filt)
                self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
                self.matched_catalog_mag=np.asarray(self.matched_catalog['mag'])
                self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

            elif self.auto_cat !='auto':
                # print(self.ref_cat.loc[self.indx[self.upd_indx]])
                self.matched_catalog= self.ref_cat.loc[self.indx[self.upd_indx]] 
                self.matched_catalog_mag=np.asarray([float(m) for m in self.matched_catalog['mag']])
                self.matched_star_coords_pix = wcs_to_pixels(self.sci_ali_name,self.matched_catalog[['ra','dec']])

            
    
        if self.termoutp!='quiet':
            print(info_g+" Stars detected in the image=",len(self.star_coords_pix))


            print(info_g+" Length of matched catalog=",len(self.matched_catalog_mag))



        if len(self.matched_catalog_mag)<5:
            if self.termoutp!='quiet':
                print(warn_y+" Few stars ("+str(len(self.matched_catalog_mag))+") matched!")

        # sys.exit()

        # if not os.path.exists(path+'calb_stars'):
        #     os.makedirs(self.path+'calb_stars')


        #stars found iraf star finder in convolved scince image
        # CircularAperture(self.star_coords_pix, r=12.).plot(color='yellow', lw=4, alpha=0.4)
        #stars in referecne in coorinates of the science
        # CircularAperture(self.ref_coords_pix, r=9.).plot(color='red', lw=4, alpha=0.4)
        # matched stars
        # CircularAperture(self.matched_star_coords_pix, r=6.).plot(color='blue', lw=4, alpha=0.4)


        if (len(self.matched_catalog_mag)<=1 and self.special_case==None) or (self.special_case!=None and len(self.matched_catalog_mag)<1):
            # print(self.special_case, len(self.matched_catalog_mag))
            if self.output_to_log!=False:
                self.log_error_4 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Less than 3 matched calibration stars '+'\n'
                self.night_log = open(self.log_name,"r+")
                if self.log_error_4 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_4)
                    self.night_log.write('\n')
                    self.night_log.close()
            print(warn_r+" Exiting: not enough stars to calibrate!")
            self.sys_exit=True
            return
            
        if self.termoutp!='quiet':
            print(info_g+' Finished the star matching process')

        keep_indx=[]
        new_matched_catalog_mag,new_matched_catalog,new_matched_star_coords_pix=[],[],[]
        # print(self.matched_catalog_mag)
        # [print(self.matched_catalog_mag[k]) for k in range(len(self.matched_catalog_mag)) if self.matched_catalog_mag[k]!='-' and float(self.matched_catalog_mag[k])>13.5]
        [keep_indx.append(k) for k in range(len(self.matched_catalog_mag)) if self.matched_catalog_mag[k]!='-' and float(self.matched_catalog_mag[k])>13.5]
        # print(keep_indx)
        for i in keep_indx:
            new_matched_catalog_mag.append(float(self.matched_catalog_mag[i]))
            try:new_matched_catalog.append(self.matched_catalog[i]) 
            except: new_matched_catalog.append(self.matched_catalog.iloc[i])
            new_matched_star_coords_pix.append(self.matched_star_coords_pix[i])
        print(info_g+' Number of stars kept after magnitude cut:',len(new_matched_catalog_mag))
        # print(new_matched_catalog)
        # self.matched_catalog_mag,self.matched_catalog,self.matched_star_coords_pix=np.array(new_matched_catalog_mag),np.array(new_matched_catalog),np.array(new_matched_star_coords_pix)
        self.matched_catalog_mag,self.matched_star_coords_pix=np.array(new_matched_catalog_mag),np.array(new_matched_star_coords_pix)
        
        #remove stars that fall within an ellipse centred at the centre of the image with 

        
        return self.ref_cat,self.matched_catalog,self.matched_catalog_mag,self.matched_star_coords_pix




    def combine_psf(self):
        if self.termoutp!='quiet':
            print(info_g+f" Combining PSFs")
        ###################################################################
        #  Combine science PSF with reference PSF, flatten PSF to 1D and gaussian fit

        if not os.path.exists(self.path+'convolved_psf'):
            os.makedirs(self.path+'convolved_psf')
        
        try:
            if os.path.exists(self.path+'convolved_psf'):
                os.system('rm '+self.path+'convolved_psf/*.fits')
        except:
            pass

        self.comb_psf = scipy_convolve(self.kernel_sci, self.kernel_ref, mode='same', method='fft')
        self.comb_psf=self.comb_psf/np.sum(self.comb_psf)
        self.hdu_comb_psf_name =self.path+"convolved_psf/"+self.sci_obj+'_'+self.sci_filt+'_'+self.sci_img_hdu.header[self.DATE_kw][0:13]+f"comb_psf_{self.rand_nums_string}.fits"
        self.hdu_comb_psf= fits.PrimaryHDU(self.comb_psf)
        self.hdu_comb_psf.writeto(self.hdu_comb_psf_name,overwrite=True)




    def cutout_psf(self,data,psf_array,xpos,ypos):
        # print('x',xpos,'y',ypos)
        xcutout,ycutout=np.shape(psf_array)[0],np.shape(psf_array)[1]
        position=(xpos,ypos)
        size = (xcutout, ycutout)
        # print('size',size,'position',position,'data shape',np.shape(data))

        cutout = Cutout2D(data, position, size)#, mode='partial')
        # print('cutout shape',np.shape(cutout.data))
        return(cutout.data)

    def psf_fit(self,data_cutout,psf_array):
        #fit PSF with x or y-shift
        # shift cutout psf to be aligned with psf model
        # restrict it to move <5 pixels in each direction so the fit doesn't go wild and fit a nearby bright star etc.
        xoff, yoff, exoff, eyoff = chi2_shift(psf_array,data_cutout, 10,return_error=True, upsample_factor='auto')
        #convert to arcsec
        xoff_arc=abs(xoff)*self.sci_ps
        yoff_arc=abs(yoff)*self.sci_ps
        if (xoff>5.0 or yoff>5.0) and self.telescope not in SEDM:
            xoff,yoff=0.0,0.0
        data_cutout_shift=scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
        resize_sci=np.reshape(data_cutout_shift,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1])
        #resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1],1)
        resize_psf=np.reshape(psf_array,np.shape(data_cutout_shift)[0]*np.shape(data_cutout_shift)[1]) 
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)

        return(slope, intercept, r_value*r_value,xoff, yoff, abs(xoff_arc), abs(yoff_arc))

    def psf_fit_noshift(self,data_cutout,psf_array):
        # forced photometry - fit with no x or y-shift for limits calculations
        resize_sci=np.reshape(data_cutout,np.shape(data_cutout)[0]*np.shape(data_cutout)[1])
        # print('data_cutout',np.shape(data_cutout))
        # print('resize_sci',np.shape(resize_sci))
        # print('resize_psf',np.shape(psf_array))
        resize_psf=np.reshape(psf_array,np.shape(data_cutout)[0]*np.shape(data_cutout)[1])

        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(resize_psf, resize_sci)
        return(slope, intercept, r_value*r_value)



    def get_zeropts(self):
        if self.termoutp!='quiet':
            print(info_g+f' Calculating zeropoints')
            print(info_g+' Size of sci_conv: '+str(np.shape(self.sci_conv)))
            print(info_g+' Size of ref_conv: '+str(np.shape(self.ref_conv)))
            print(info_g+' Size of comb_psf: '+str(np.shape(self.comb_psf)))
        ###########################################################
        # Zero points!

        self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref=[],[],[],[]
        # print(self.matched_star_coords_pix)

        if len(self.matched_star_coords_pix)==0:
            #check if args.relative_flux is set to True
            if self.relative_flux:
                #if so, then we set a default zeropoint of 1
                print(info_b+f' No matched stars found in {self.sci_obj} image, but relative flux is set to True, so calculating relative fluxes and ignoring zeropoint')
                self.zp_sci.append(1)
                self.rsq_sci.append(1)
                return
            else:
                print(info_r+f' No matched stars found in {self.sci_obj} image, exiting')
                self.sys_exit=True
                return

        if len(self.matched_star_coords_pix)>0:
            #check the countsto reject saturated stars in the science image
            self.counts = []

            for i in range(0,len(self.matched_star_coords_pix)):
                self.cutout_sci=self.cutout_psf(data=self.sci_conv,psf_array=self.comb_psf,xpos=self.matched_star_coords_pix[i][0],ypos=self.matched_star_coords_pix[i][1])
                # print(info_g+' Cutout shape',np.shape(self.cutout_sci))
                if np.shape(self.cutout_sci)!=np.shape(self.comb_psf):
                    #pad with 0s to make it the same size as the PSF
                    self.cutout_sci = np.pad(self.cutout_sci,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_sci)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_sci)[1])),'constant',constant_values=0)
                self.sci_psf_fit = self.psf_fit_noshift(data_cutout=self.cutout_sci,psf_array=self.comb_psf)
                self.counts.append(self.sci_psf_fit[0])
            self.counts=np.array(self.counts)
            self.counts_o = self.counts
            self.counts=self.counts[~np.isnan(self.counts)]
            self.counts=self.counts[~np.isinf(self.counts)]
            self.counts=self.counts[self.counts>0]
            # print(self.counts)
            # print(self.counts[np.where(self.counts<=50000)])

            if self.telescope not in SEDM:
                print(info_g+' Number of stars before saturation check:',len(self.matched_star_coords_pix))
                self.matched_new_pix,self.matched_new_mag = [],[]
                for k in range(0,len(self.counts_o)):
                    # print(self.counts_o[k])
                    if not np.isnan(self.counts_o[k]) and not np.isinf(self.counts_o[k]) and self.counts_o[k]>0 and self.counts_o[k]<=450000:
                        if self.sci_filt in ['g','r','i','z']:
                            sat_mag = 14.5
                        else:
                            sat_mag = 15
                        if self.matched_catalog_mag[k]>sat_mag:
                            self.matched_new_pix.append(self.matched_star_coords_pix[k])
                            self.matched_new_mag.append(self.matched_catalog_mag[k])

                self.matched_star_coords_pix = np.array(self.matched_new_pix)
                self.matched_catalog_mag = np.array(self.matched_new_mag)
                print(info_g+' Number of stars after saturation check:',len(self.matched_star_coords_pix))


            
            for i in range(0,len(self.matched_star_coords_pix)):

                if self.special_case!=None:
                    if any(phrase in self.special_case for phrase in ['SN2023ixf','bright','brightstars.cat','sn2023ixf','2023ixf','brightstars']):
                        if self.matched_catalog_mag[i]>self.bright_cat_mag_lim:
                            continue
                self.cutout_sci=self.cutout_psf(data=self.sci_conv,psf_array=self.comb_psf,xpos=self.matched_star_coords_pix[i][0],ypos=self.matched_star_coords_pix[i][1])
                if np.shape(self.cutout_sci)!=np.shape(self.comb_psf):
                    self.cutout_sci = np.pad(self.cutout_sci,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_sci)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_sci)[1])),'constant',constant_values=0)
                    print(info_g+' New shape',np.shape(self.cutout_sci))
                self.sci_psf_fit = self.psf_fit(data_cutout=self.cutout_sci,psf_array=self.comb_psf)
                self.zp_sci.append(2.5*np.log10(self.sci_psf_fit[0])+self.matched_catalog_mag[i])
                self.rsq_sci.append(self.sci_psf_fit[2])


                self.cutout_ref=self.cutout_psf(data=self.ref_conv,psf_array=self.comb_psf,xpos=self.matched_star_coords_pix[i][0],ypos=self.matched_star_coords_pix[i][1])
                if np.shape(self.cutout_ref)!=np.shape(self.comb_psf):
                    self.cutout_ref = np.pad(self.cutout_ref,((0,np.shape(self.comb_psf)[0]-np.shape(self.cutout_ref)[0]),(0,np.shape(self.comb_psf)[1]-np.shape(self.cutout_ref)[1])),'constant',constant_values=0)
                    print(info_g+' New shape',np.shape(self.cutout_ref))
                self.ref_psf_fit = self.psf_fit(data_cutout=self.cutout_ref,psf_array=self.comb_psf)
                self.zp_ref.append(2.5*np.log10(self.ref_psf_fit[0])+self.matched_catalog_mag[i]) 
                self.rsq_ref.append(self.ref_psf_fit[2])



        self.zp_sci,self.zp_ref=np.array(self.zp_sci),np.array(self.zp_ref)
        self.rsq_ref,self.rsq_sci=np.array(self.rsq_ref),np.array(self.rsq_sci)

        if np.count_nonzero(~np.isnan(self.zp_ref))<2 and self.special_case==None:
            if np.count_nonzero(~np.isnan(self.zp_sci))<3:
                if self.termoutp!='quiet':
                    print(warn_y+' Few stars in the field')
                    self.log_error_4_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': <3 stars found in the field '+'\n'
                    self.night_log = open(self.log_name,"r+")
                    if self.log_error_4_5 in self.night_log.readlines():
                        pass
                    else:
                        self.night_log.write(self.log_error_4_5)
                        self.night_log.write('\n')
                        self.night_log.close()

             #errors.write('Few (< 3) field stars to calibrate in science and reference image. \n')
            if np.count_nonzero(~np.isnan(self.zp_sci))>3:
                print(warn_r+' Something wrong with reference since no zero point stars detected, but detected in science')
                print(warn_r+' Exiting..')



        if self.sci_filt=='u':
            self.thresh=0.1
        if self.sci_filt!='u':
            self.thresh=0.8
        self.rsq_cont=False
        #check how many stars have rsq values above the threshold and lower by 0.05 until there are at least 2 stars
        while self.rsq_cont==False:
            self.arr=(self.rsq_ref >=self.thresh) & (self.rsq_sci >=self.thresh) & (self.zp_sci >=0) & (self.zp_ref >=0) #filter out stars poorly fitted with PSF fit & weird zp measurement
            # print(self.rsq_ref[self.arr],sum(self.arr))
            if len([i for i in self.rsq_ref if i>=self.thresh])<2 or len([i for i in self.rsq_sci if i>=self.thresh])<2:
                print(warn_y+' Lowering threshold for rsq values')
                self.thresh-=0.05
            else:
                self.rsq_cont=True
                break

            if self.thresh<0.45:
                print(warn_r+' No stars with rsq values above threshold')
                print(warn_r+' Exiting..')
                self.sys_exit=True
                self.rsq_cont=True
                break

        self.arr=(self.rsq_ref >=self.thresh) & (self.rsq_sci >=self.thresh) & (self.zp_sci >=0) & (self.zp_ref >=0) #filter out stars poorly fitted with PSF fit & weird zp measurement

        if all(val ==False for val in self.arr)==True:
            weird_img = 'science and reference'
            if all(val_sci_rsq <self.thresh for val_sci_rsq in self.rsq_sci)==True:
                self.weird_img = 'science'
            else:
                self.weird_img = 'reference'


            print(warn_r+f' All stars in the {self.weird_img} images fit poorly or returned anomalous zeropoint measurements!')
            if self.output_to_log!=False:
                self.log_error_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+f': All stars in the {self.weird_img} images fit poorly or returned anomalous zeropoint measurements '+'\n'
                self.night_log = open(self.log_name,"r+")
                if self.log_error_5 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_5)
                    self.night_log.write('\n')
                    self.night_log.close()
            print(warn_r+" Exiting...")
            self.sys_exit=True
            return

        self.rsq_ref,self.rsq_sci = self.rsq_ref[self.arr],self.rsq_sci[self.arr]
        self.zp_sci,self.zp_ref = self.zp_sci[self.arr],self.zp_ref[self.arr]

        self.zp_sci,self.zp_ref=sigma_clip(self.zp_sci,sigma=5,maxiters=4),sigma_clip(self.zp_ref,sigma=5,maxiters=4)
        self.zp_sci=self.zp_sci[~self.zp_sci.mask]
        self.zp_ref=self.zp_ref[~self.zp_ref.mask]

        if self.termoutp!='quiet':
            print(info_g+' zp sci=%.3f std=%.3f num_stars=%i'%(np.nanmedian(self.zp_sci),np.nanstd(self.zp_sci),len(self.zp_sci)))
            print(info_g+' zp ref=%.3f std=%.3f num_stars=%i'%(np.nanmedian(self.zp_ref),np.nanstd(self.zp_ref),len(self.zp_ref)))

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


        return self.zp_sci,self.zp_ref,self.rsq_sci,self.rsq_ref

    def check_nearby_resid(self,data,psf,sn_x,sn_y,sn_phot):
        #thid function checks if there are any nearby residuals that could be affecting the fit within 5 arcsec of sn pos
        #if there are, it will return a warning and exit the program
        #if there are not, it will return the residuals
        #it does this by measuring the photometry at 20 random locations outside the psf fit
        #if the photometry is significantly different from the photometry at the sn position, it will return a warning and exit the program
        #if the photometry is not significantly different, it will return the residuals

        if self.termoutp!='quiet':
            print(info_g+' Checking for nearby residuals')

        #create a list of 20 random positions outside the psf fit but within 5 arcsec of the sn position
        rand_x,rand_y=[],[]
        #calculate 5 arcsec in pixels
        r = 5/self.sci_ps
        for i in range(50):
            rand_x.append(random.uniform(sn_x-5,sn_x+r/2))
            rand_y.append(random.uniform(sn_y-5,sn_y+r/2))

        #calculate the photometry at the random positions
        rand_phot=[]
        for i in range(len(rand_x)): #_fit_noshift(data_cutout=bkg_cutout,psf_array=psf)[0]
            cutout=self.cutout_psf(data=data,psf_array=psf,xpos=rand_x[i],ypos=rand_y[i])
            rand_phot.append(self.psf_fit_noshift(cutout,psf_array=psf)[0])

        #check the distribution of the photometry at the random positions
        rand_phot=np.array(rand_phot)
        # print(rand_phot)
        print(np.nanmedian(rand_phot),np.nanstd(rand_phot))
        print(np.max(rand_phot),np.min(rand_phot))
        print()
        print('sn_phot',sn_phot)
        print()
        #find the difference between the photometry at the sn position and the photometry at the random positions
        diff_phot=sn_phot-rand_phot
        for j in range(len(diff_phot)):
            #find thr median difference excluding the jth value
            temp_diff_phot=np.delete(diff_phot,j)
            med_diff_phot=np.nanmedian(temp_diff_phot)
            std_diff_phot=np.nanstd(temp_diff_phot)
            mean_diff_phot=np.nanmean(temp_diff_phot)

            #if the jth value is more than 3 sigma away from the median diff, print a warning and exit the program
            if diff_phot[j]>med_diff_phot+3*std_diff_phot or diff_phot[j]<med_diff_phot-3*std_diff_phot:
                #see if the diff_phot[j] is more negative than the surrounding background
                if diff_phot[j]<np.nanmedian(rand_phot)-np.nanstd(rand_phot):

                    if self.termoutp!='quiet':
                        print('residuals detected within 5 arcsec of sn position')
                        print()
                




    def mag_err_function(self,data,psf,zp_sci,num,sn_x,sn_y):
        mag=99.0
        magstd=99.0
        magerr=99.0
        minimal_mag=99.0
        if self.termoutp!='quiet':
            print(info_g+' Calculating magnitudes')
        magerr,maglim,flux_bkg_list,flux_new_sn_list=[],[],[],[]
        self.ra_off,self.dec_off = None,None
    
        # cutout a part of the image the same size as the measured PSF
        self.sn_cutout=self.cutout_psf(data=data,psf_array=psf,xpos=sn_x,ypos=sn_y)
        # calculate the magnitude of the object
        if self.forced_phot==False: 
            self.main_sn_psf_fit = self.psf_fit(self.sn_cutout,psf_array=psf)
        else: 
            print(info_g+f" Performing forced photometry at {self.forced_phot[0]} {self.forced_phot[1]}")
            sn_x,sn_y = self.forced_phot[0],self.forced_phot[1]
            new_sci_c = SkyCoord(ra=sn_x,dec=sn_y,unit=(u.hourangle, u.deg),frame='fk5')
            sn_x,sn_y = wcs_to_pixels(self.sci_conv_name,np.column_stack((new_sci_c.ra.deg,new_sci_c.dec.deg)))[0]

            # self.coords_sn_sub=wcs_to_pixels(self.sci_conv_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))
            # self.coords_sn_sub_x,self.coords_sn_sub_y = self.coords_sn_sub[0]

            self.sn_cutout = self.cutout_psf(data=data,psf_array=psf,xpos=sn_x,ypos=sn_y)
            self.main_sn_psf_fit = self.psf_fit_noshift(self.sn_cutout,psf_array=psf)

        sn_flux= self.main_sn_psf_fit[0]
        sn_mag=-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)
        # print the offset of the shifted fit
        if self.termoutp!='quiet' and self.forced_phot==False:
            print(info_b+' offset of shifted PSF fit (xpos,ypos)=%.3f %.3f'%(self.main_sn_psf_fit[3],self.main_sn_psf_fit[4]))
            # print(main_sn_psf_fit[3],main_sn_psf_fit[4])
            # print(main_sn_psf_fit[5],main_sn_psf_fit[6])

            # if the xoff_arc and yoff_arc are greater than 1.5 arcsec, then the fit is shifted too much, 
            # defaulting to the original position so will use psf_fit_noshift
            # print(self.main_sn_psf_fit[5],self.main_sn_psf_fit[6])
            print(info_g+' xoff_arc=%.3f arcsec, yoff_arc=%.3f arcsec'%(self.main_sn_psf_fit[5],self.main_sn_psf_fit[6]))
            self.ra_off,self.dec_off = self.main_sn_psf_fit[5],self.main_sn_psf_fit[6]
            # self.main_sn_psf_fit = main_sn_psf_fit
            # if np.abs(main_sn_psf_fit[5])>1.5 or np.abs(main_sn_psf_fit[6])>1.5 or any(offset>2 for offset in [main_sn_psf_fit[5], main_sn_psf_fit[6]]):
            if (sn_mag>17.5 and self.forced_phot==False and self.telescope not in SEDM) or (self.telescope in SEDM and self.forced_phot!=False):
                if any(offset>self.max_psf_offset for offset in [self.main_sn_psf_fit[5], self.main_sn_psf_fit[6]]):
                    print(warn_y+f" PSF fit is shifted too much, defaulting to original position")
                    print(warn_y+" xoff =%.3f arsec, yoff_arc=%.3f arcsec"%(self.main_sn_psf_fit[5],self.main_sn_psf_fit[6]))
                    self.main_sn_psf_fit = self.psf_fit_noshift(self.sn_cutout,psf_array=psf)
                    sn_flux= self.main_sn_psf_fit[0]
                    sn_mag=-2.5*np.log10(sn_flux)+np.nanmedian(zp_sci)
                    print(warn_y+" New magnitude = %.3f"%sn_mag)
            

        # set up a grid for the values
        psf_size=np.shape(psf)[0]+1
        # x=np.linspace(-num*psf_size,num*psf_size,(num*2)+1)
        # y=x
        # x0, y0, radius = 0.0, 0.0, psf_size/2
        # x, y = np.meshgrid(x, y)
        # r = np.sqrt((x - x0)**2 + (y - y0)**2)
        # outside = (r > radius)
        # x=x[outside]
        # y=y[outside]

        np.shape(psf)[0]+1
        #chose num number of coordinates to calculate the magnitude error within the image size but outside the psf
        x=np.linspace(-num*psf_size,num*psf_size,(num*2)+1)
        # print('x min max',x.min(),x.max())
        # print('data shape',np.shape(data))
        y=x
        x0, y0, radius = 0.0, 0.0, psf_size/2
        x, y = np.meshgrid(x, y)
        r = np.sqrt((x - x0)**2 + (y - y0)**2)
        outside = (r > radius)
        x=x[outside]
        y=y[outside]


   
        #measure the background flux and the flux of the artificial supernova
        # print(info_g+' Calculating magnitude error for %i of %i'%(i+1,len(x)))
        for i in range(0,len(x)):
            x__,y__=x[i],y[i]
            if sn_x+x__<0 or sn_x+x__>np.shape(data)[1] or sn_y+y__<0 or sn_y+y__>np.shape(data)[0] or (sn_x+x__-np.shape(psf)[0]/2)<0 or (sn_y+y__-np.shape(psf)[1]/2)<0 or (sn_x+x__+np.shape(psf)[0]/2)>np.shape(data)[1] or (sn_y+y__+np.shape(psf)[1]/2)>np.shape(data)[0]:
                # print(warn_y+' Artificial supernova is outside the image')
                # print(x__,y__)
                # continue
                # print(warn_y+' Randomly selecting a new position')
                # continue
                x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))

                if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                    c=0
                    while sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                        x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))
                        c+=1
                        if c>20:
                            break
                        if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                            # print(warn_y+f' Artificial supernova is outside the image, skipping')
                            continue
                    # x__,y__=abs(x__),abs(y__)
                    # print(warn_y+' New position is x=%.3f, y=%.3f'%(sn_x+x__,sn_y+y__),'#',1)
                    # if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                    #     x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))
                    #     # x__,y__=abs(x__),abs(y__)
                    #     print(warn_y+' New position is x=%.3f, y=%.3f'%(sn_x+x__,sn_y+y__), '#',2)
                    #     if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                    #         x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))
                    #         # x__,y__=abs(x__),abs(y__)
                    #         print(warn_y+' New position is x=%.3f, y=%.3f'%(sn_x+x__,sn_y+y__), '#',3)
                    #         if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                    #             x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))
                    #             # x__,y__=abs(x__),abs(y__)
                    #             print(warn_y+' New position is x=%.3f, y=%.3f'%(sn_x+x__,sn_y+y__), '#',4)
                    #             if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                    #                 x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))
                    #                 # x__,y__=abs(x__),abs(y__)
                    #                 print(warn_y+' New position is x=%.3f, y=%.3f'%(sn_x+x__,sn_y+y__), '#',5)
                    #                 if sn_x+x__-radius<0 or sn_x+x__+radius>np.shape(data)[1] or sn_y+y__-radius<0 or sn_y+y__+radius>np.shape(data)[0]:
                    #                     x__,y__=np.random.randint(-num*psf_size,num*psf_size,(2))
                    #                     # x__,y__=abs(x__),abs(y__)
                    #                     print(warn_y+' New position is x=%.3f, y=%.3f'%(sn_x+x__,sn_y+y__), '#',6)
            # print('data',np.shape(data),'psf',np.shape(psf))
            bkg_cutout=self.cutout_psf(data=data,psf_array=psf,xpos=sn_x+x__,ypos=sn_y+y__)
            # print('bkg_cutout',np.shape(bkg_cutout))
            flux_bkg=self.psf_fit_noshift(data_cutout=bkg_cutout,psf_array=psf)[0]
            #flux_bkg_list is the PSF fitted to the sky
            flux_bkg_list.append(flux_bkg)

            new_sn=bkg_cutout+((psf)*self.psf_fit(self.sn_cutout,psf_array=psf)[0])+self.psf_fit(self.sn_cutout,psf_array=psf)[1]
            flux_new_sn=self.psf_fit(data_cutout=new_sn,psf_array=psf)[0]
            #flux_new_sn_list is the PSF fitted to the artificial supernova
            flux_new_sn_list.append(flux_new_sn)

        flux_bkg_list=sigma_clip(flux_bkg_list,sigma=3,maxiters=6)

        flux_new_sn_list=sigma_clip(flux_new_sn_list,sigma=3,maxiters=6)

        SNR = sn_flux/np.std(flux_bkg_list)
        if self.termoutp!='quiet':
            print(info_g+' S/N (std background)= %.3f'%(sn_flux/np.std(flux_bkg_list))) 
            print(info_g+' S/N (std artifical sn)= %.3f'%(sn_flux/np.std(flux_new_sn_list)))
        
        
        if sn_flux/np.std(flux_bkg_list)<=2:
            #not detected
            minimal_mag=sn_mag
            mag=99.0
            magstd=99.0
            magerr=99.0
            if self.termoutp!='quiet':
                print(info_g+' Less than 2.5 sigma')
                print(info_g+' Mag=%.3f+/-%.3f ' %(sn_mag,-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.std(flux_new_sn_list))))
                print(info_g+' Flux corresponding to the (marginal) detection + 2*sigma =%.3f'%(-2.5*np.log10(np.median(sn_flux)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci)))
                print(info_g+' 5-sig limit =%.3f'%(-2.5*np.log10(abs(np.median(flux_bkg_list))+abs(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci)))
                print(info_g+' 3-sig limit =%.3f'%(-2.5*np.log10(abs(np.median(flux_bkg_list))+abs(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci)))
            # pick the most conservative out of:
            #flux corresponding to the (marginal) detection + 2*sigma')
            #'3-sig limit
            limits=[-2.5*np.log10(np.median(sn_flux)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci),-2.5*np.log10(np.median(flux_bkg_list)+abs(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci)]
            maglim=np.min(limits)

        
        if sn_flux/np.std(flux_bkg_list)>=2:
            minimal_mag=sn_mag
            mag=sn_mag
            
            magstd=-2.5*np.log10(sn_flux)+2.5*np.log10(sn_flux+np.std(flux_new_sn_list))
            # print('magstd',magstd)
            if self.sci_filt=='g' or self.sci_filt=='u':
                self.sys_fact = 1.07
            else:
                self.sys_fact = 1.03
            
            # print('fact',sn_mag,magstd,magstd*self.sys_fact)
            magerr=magstd*self.sys_fact
            # print('magerr',magerr)
            # print('magerr',magerr,magstd/(len(flux_new_sn_list)**0.5),len(flux_new_sn_list),len(flux_new_sn_list)**0.5)
            maglim=-2.5*np.log10(abs(np.median(flux_bkg_list))+abs(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci)
            # print('maglim',maglim)

            if self.termoutp!='quiet':
                flux_bkg_list = [i for i in flux_bkg_list if i !='--']
                flux_new_sn_list = [i for i in flux_new_sn_list if i !='--']
                # print(final_flux_bkg_list,np.std(final_flux_bkg_list),np.median(final_flux_bkg_list))
                # print(final_flux_new_sn_list,np.std(final_flux_new_sn_list),np.median(final_flux_new_sn_list))

                print(info_g+' Mag=%.3f+/-%.3f ' %(sn_mag,magerr))
                print(info_g+' Flux corresponding to the (marginal) detection + 2*sigma',-2.5*np.log10(np.median(sn_flux)+(np.std(flux_bkg_list)*2))+np.nanmedian(zp_sci))
                print(info_g+' 5-sig limit =',-2.5*np.log10(abs(np.median(flux_bkg_list))+abs(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci))
                print(info_g+' 3-sig limit =',-2.5*np.log10(abs(np.median(flux_bkg_list))+abs(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci))
                # print(flux_bkg_list)
                # print('lims',np.std(flux_bkg_list),np.median(flux_bkg_list),np.median(flux_bkg_list)+np.std(flux_bkg_list)*2)

        # * If >3 sigma detection: report to Marshal with 1-sigma errors.
        # * If <3 sigma detection:r eport whichever of these two is more conservative:
        # * The flux corresponding to the (marginal) detection + 2*sigma.
        # * The flux corresponding to 3*sigma.

        sn_flux=sn_flux/3631
        sn_flux_err = abs(np.std(flux_new_sn_list))/3631
        # print(magerr,magstd)

        self.check_nearby_resid(data,psf,sn_x,sn_y,sn_mag)
        return(mag,magstd,magerr,maglim,sn_flux/np.std(flux_bkg_list),sn_flux/np.std(flux_new_sn_list),-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*1))+np.nanmedian(zp_sci),
        -2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*3))+np.nanmedian(zp_sci),-2.5*np.log10(np.median(flux_bkg_list)+(np.std(flux_bkg_list)*5))+np.nanmedian(zp_sci),
        minimal_mag,SNR,sn_flux,sn_flux_err)




    def scaled_subtract(self):
        if self.termoutp!='quiet':
            print(info_g+f" Subtracting scaled & convolved refrence image")
        ###########################################################
        # Scale and subtract

        if not os.path.exists(self.path+'scaled_subtracted_imgs'):
            os.makedirs(self.path+'scaled_subtracted_imgs')

        if not os.path.exists(self.path+'no_scaled_subtracted_imgs'):
            os.makedirs(self.path+'no_scaled_subtracted_imgs')

        #this is another check on the std of the zero point of the image! np.nanstd(zp_ref)/np.sqrt(len(zp_ref))
        if self.relative_flux==True:
            self.scale_factor=1
        else:
            self.scale_factor=np.power(10,(np.nanmedian(self.zp_sci)-np.nanmedian(self.zp_ref))/2.5)
        if self.termoutp!='quiet':
            print(info_g+' Scale factor: %.3f'%self.scale_factor)
            # print(info_g+' Image dimensions: Science %s, Reference: %s'%(np.shape(self.sci_conv),np.shape(self.ref_conv)))

        if self.to_subtract==True:
            self.sci_scale_sub_name = self.path+'scaled_subtracted_imgs/'+self.sci_img_name[:-5]+'_scaled_subtraction.fits'
            self.sub=(self.sci_conv)-(self.scale_factor*self.ref_conv)
            self.hdu_fits_sub= fits.PrimaryHDU(self.sub)
            self.hdu_fits_sub.writeto(self.sci_scale_sub_name,overwrite=True)

            #scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
            #sub=(sci_conv)-(scale_factor*scipy.ndimage.shift(ref_conv, [-0.5, -0.5], order=3, mode='reflect', cval=0.0, prefilter=True))

            self.coords_sn_sub=wcs_to_pixels(self.sci_conv_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))
            self.coords_sn_sub_x,self.coords_sn_sub_y = self.coords_sn_sub[0]


        if self.to_subtract==False:
            self.sci_scale_sub_name = self.path+'no_scaled_subtracted_imgs/'+self.sci_img_name[:-5]+'no_scaled_subtraction.fits'
            self.sub=(self.sci_conv)
            self.hdu_fits_sub= fits.PrimaryHDU(self.sub)
            self.hdu_fits_sub.writeto(self.sci_scale_sub_name,overwrite=True)

            #scipy.ndimage.shift(data_cutout, [-yoff, -xoff], order=3, mode='reflect', cval=0.0, prefilter=True)
            #sub=(sci_conv)-(scale_factor*scipy.ndimage.shift(ref_conv, [-0.5, -0.5], order=3, mode='reflect', cval=0.0, prefilter=True))

            self.coords_sn_sub=wcs_to_pixels(self.sci_conv_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))
            self.coords_sn_sub_x,self.coords_sn_sub_y = self.coords_sn_sub[0]

        
        self.sub_mean, self.sub_median, self.sub_std = sigma_clipped_stats(self.hdu_fits_sub.data)
        self.coords_sn=self.coords_sn_sub

        self.sub_img = fits.open(self.sci_scale_sub_name)[0]

        
        # self.files_to_clean.append(self.sci_scale_sub_name)
        return self.sub,self.sci_scale_sub_name,self.coords_sn_sub



    def get_photometry(self):
        if self.termoutp!='quiet':
            print(info_g+f" Extracting photometry")
        self.t_photometry_start = time.time()

        self.mag=self.mag_err_function(data=self.sub,psf=self.comb_psf,zp_sci=self.zp_sci,num=10,sn_x=self.coords_sn[0][0],sn_y=self.coords_sn[0][1])
        # print(self.mag[2])
        # print(np.nanstd(self.zp_ref),len(self.zp_ref),np.nanstd(self.zp_ref)/len(self.zp_ref)**0.5)
        # print(np.nanstd(self.zp_sci),len(self.zp_sci),np.nanstd(self.zp_sci)/len(self.zp_sci)**0.5)
        # print()

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

            self.ref_img_ali_hdu = fits.open(self.ref_ali_name)[0]
            self.coords_sn_ref_ali=wcs_to_pixels(self.ref_ali_name,np.column_stack((self.sci_c.ra.deg,self.sci_c.dec.deg)))[0]
            self.coords_sn_ref_ali_x,self.coords_sn_ref_ali_y = self.coords_sn_ref_ali

            self.vmin_ref,self.vmax_ref = visualization.ZScaleInterval().get_limits(self.ref_img_ali_hdu.data)
            self.ax[1].imshow(self.ref_img_ali_hdu.data,cmap=self.cmap,vmin=self.vmin_ref, vmax=self.vmax_ref)
            self.ax[1].set_xlim(self.coords_sn_ref_ali_x-50,self.coords_sn_ref_ali_x+50)
            self.ax[1].set_ylim(self.coords_sn_ref_ali_y-50,self.coords_sn_ref_ali_y+50)
            self.ax[1].set_title('Ref')
            self.ax[1].axis('off')        

            self.fig.savefig(self.cutout_name)

            try:
                self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sub_img.data)
                self.ax[2].imshow(self.sub_img.data,cmap=self.cmap,vmin=self.vmin, vmax=self.vmax)
                self.ax[2].set_xlim(self.coords_sn_sub_x-115,self.coords_sn_sub_x+115)
                self.ax[2].set_ylim(self.coords_sn_sub_y-115,self.coords_sn_sub_y+115)
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

        # print(self.mag_all_err**0.5)
        # print(np.nanstd(self.zp_ref),len(self.zp_ref),(np.nanstd(self.zp_ref)/(len(self.zp_ref)**0.5)))
        self.mag_all_err+=(np.nanstd(self.zp_ref))**2
        print(info_g+' Reference zeropoint std: '+str(np.nanstd(self.zp_ref)))
        # print(self.mag_all_err**0.5)
        # print(np.nanstd(self.zp_sci),len(self.zp_sci),(np.nanstd(self.zp_sci)/(len(self.zp_sci)**0.5)))
        self.mag_all_err+=(np.nanstd(self.zp_sci))**2
        print(info_g+' Science zeropoint std: '+str(np.nanstd(self.zp_sci)))
        # print(self.mag_all_err**0.5)
        self.mag_all_err=self.mag_all_err**0.5
        # print(self.mag_all_err)

        # print(info_g+f" Magnitude: {self.mag[0]} +/- {self.mag_all_err}")
        self.SNR = self.mag[10]
        if self.mag[0]>40:
            if self.output_to_log!=False:
                self.log_error_5_5 = self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+': Magnitude of 99 measured, SNR= '+str(self.SNR)+'\n'
                self.night_log = open(self.log_name,"r+")
                if self.log_error_5_5 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_5_5)
                    self.night_log.write('\n')
                    self.night_log.close()
        

        self.t_photometry_end = time.time()
        self.photometry_time = self.t_photometry_end - self.t_photometry_start
        
        if self.termoutp!='quiet':
            # print("(mag_err**2+std_zp**2+std_ref**2)^(0.5) %.3f \n" %(self.mag_all_err))

            if self.relative_flux!=True:
                print(colored('----------------------------------------', 'blue'))
                # print(self.mag[0],self.mag_all_err,self.mag[3])
                print(' Mag = %.3f+/-%.3f lim=%.3f' %(self.mag[0],self.mag_all_err,self.mag[3]))
                print(colored('----------------------------------------', 'blue'))

            if self.relative_flux==True or self.mag[0]>40:
                print(colored('----------------------------------------', 'blue'))
                print('Flux = %.3f+/-%.3f' %(self.mag[11],self.mag[12]))
                print(colored('----------------------------------------', 'blue'))

            print(info_g+' Photometry time:',"%.1f" % round(self.photometry_time, 0),'seconds')

            self.memkb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print(info_g+' Memory usage Mbyte:'," %5.1f" %(self.memkb/1024.0/1024.0))


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
    
            print(info_g+f' Photometry results written to {self.phot_file_name}')
        if self.out_dir=="photometry/":
            shutil.copy2(self.phot_file_name, self.phot_date_name)


        print(info_g+' Total time:',"%.1f" % round(self.t_end-self.t_start, 0),'seconds')

        return {"obj":self.sci_obj,"filt":f"sdss{self.sci_filt}",'mjd':str(self.sci_mjd),"mag":self.mag[0],"mag_err":self.mag_all_err,"lim_mag":self.mag[3],
                "ra":self.ra_string,"dec":self.dec_string,"exp_t":self.sci_exp_time,"flux":self.mag[11],"flux_err":self.mag[12],'seeing':self.sci_seeing}


    def delete_photometry(self,phot_id_data):
        self.response = api("delete", f"api/photometry/{phot_id_data['phot_id']}")
        if str(self.response.status_code)=='200':
            print(info_g+f" Successfully deleted photometry for phot ID: {phot_id_data['phot_id']}, mjd: {phot_id_data['mjd']}")
        else:
            print(warn_r+f" Unable to delete photometry for phot ID: {phot_id_data['phot_id']},mjd: {phot_id_data['mjd']} STATUS CODE={self.response.status_code}")


    def upload_phot(self):
        if self.sci_obj=='SN2023ixf' or self.sci_obj=='2023ixf' or self.sci_obj=='ZTF23aaklqou':
            self.sci_obj='2023ixf'
            return

        if self.termoutp!='quiet':
            print(info_g+f" Attempting to upload photometry")
        self.data = {"filter":f"sdss{self.sci_filt}","magerr": self.mag_all_err,"obj_id": self.sci_obj,"origin": "LT_IOO_PIPE",
                    "mag":self.mag[0],"limiting_mag": self.mag[3],"mjd": self.sci_mjd,"instrument_id": 33,"magsys": "ab","group_ids":"all"}
        self.upload_new=True
        if self.if_stacked==True:
            self.data['altdata'] = {'stacked':True,'no_in_stack':self.no_stacked}
        self.name=self.sci_obj

        if any(X in self.name for X in ['-ugriz','-griz','-gri']):
            self.name = re.sub(r'-ugriz|-griz|-gri', '', self.name)
        [print(key,item) for key,item in self.data.items()]

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
                self.fritz_ind = [ind for ind,val in enumerate(self.all_fritz_data['instrument_id']) if val==33 and self.all_fritz_data['filter'].iloc[ind]==self.data['filter']]

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
                    self.data = {'filter':f"sdss{self.sci_filt}",'mag':None,'magerr':None,'mjd':self.sci_mjd,'limiting_mag_nsigma':5,'obj_id':self.sci_obj,'origin':'LT_IOO_PIPE',
                                        'magsys':'ab','group_ids':'all','limiting_mag':self.mag[3],'instrument_id':33}

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
            self.fritz_ind = [ind for ind,val in enumerate(self.all_fritz_data['instrument_id']) if val==33 and self.all_fritz_data['filter'].iloc[ind]==self.data['filter']]

            self.on_fritz_mjd = [np.round(self.all_fritz_data['mjd'].iloc[ind],6) for ind in self.fritz_ind]
            self.on_fritz_id = [self.all_fritz_data['id'].iloc[ind] for ind in self.fritz_ind]
            self.on_fritz_mag = [np.round(self.all_fritz_data['mag'].iloc[ind],3) for ind in self.fritz_ind]
            self.upload_new=True  #identifier to update photometry in the case of stacking

            if self.if_stacked==True:
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
                        self.fritz_ind = [ind for ind,val in enumerate(self.all_fritz_data['instrument_id']) if val==33 and self.all_fritz_data['filter'].iloc[ind]==self.data['filter']]

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
                while True:
                    self.response = requests.post(url=f"https://fritz.science/api/photometry",headers = {'Authorization': f'token {token}'} ,json=self.data)
                    if self.response.status_code==200:
                        break

                if self.response.status_code == 200:
                    print(info_g+f" Successfully uploaded photometry for {self.name} in {self.data['filter']} taken at {self.data['mjd']}, uploaded at {self.current_time}")
                else:
                    print(warn_r+f" Failed to upload photometry for {self.name} in {self.data['filter']}, status code ={self.response.status_code} ")
            
                return self.data

            else:
                print(info_b+f" Data taken for {self.sci_obj} in {self.data['filter']} at {self.data['mjd']} already uploaded, not uploading")

                return self.data  

        except Exception as e:
            print(warn_y+f" Unable to find Fritz photometry for {self.sci_obj}",e)
            if self.output_to_log!=False:
                self.log_error_6 = f"Fritz error: "+self.sci_obj+","+str(self.sci_mjd)+','+self.sci_filt+" band, (ie. event not on Fritz or error retrieving photometry from Fritz, check token and Fritz status)\n"
                self.night_log = open(self.log_name,"r+")
                if self.log_error_6 in self.night_log.readlines():
                    pass
                else:
                    self.night_log.write(self.log_error_6)
                    self.night_log.write('\n')
                    self.night_log.close()
        
    def clean_directory(self):
        if self.termoutp!='quiet':
            print(info_g+f" Cleaning directories")
        for self.out_file in os.listdir(self.path+"out"):
            if self.rand_nums_string in self.out_file:
                try:
                    if os.path.exists(path+'out/'+self.out_file):os.system('rm '+path+'out/'+self.out_file)
                except:
                    pass

        for self.config_file in os.listdir(self.path+"config_files"):
            if self.rand_nums_string in self.config_file:
                try:
                    if os.path.exists(self.path+'config_files/'+self.config_file):os.system('rm '+path+'config_files/'+self.config_file)
                except:
                    pass
        
        try:
            for self.temp_file in self.files_to_clean:
                if os.path.exists(self.temp_file):os.system('rm '+self.temp_file)  
        except:
            pass 
        

                  



    def calib_app_phot(self,catalog=None,r=None,calib_stand=False,calib_sci=False,sep=10,star_mask=10,show_gc=False,zp_fit_params={},ref_only=False):
        if self.termoutp!='quiet':
            print(info_g+f" Performing aperture photometry")

        if calib_stand!=False:
            self.mean,self.median, self.std = sigma_clipped_stats(self.sci_bkgsb, sigma=5.0)
            self.iraffind= IRAFStarFinder(threshold=abs(starscale*self.std)*3,fwhm=3,roundhi=0.3,minsep_fwhm=1.5)
            self.sources = self.iraffind(self.sci_bkgsb)
            self.all_zerop,self.all_zerop_err = [],[]
            self.companions = {'2MASS J02581300+0110409':["02:58:13.0144912656","+01:10:41.039608692"],
                            'Companion 1':["22:21:27.0924","+01:10:50.992"],
                            'Companion 2':["22:21:25.6474","+01:10:33.148"],
                            'SDSS J224146.86+011204.0':["22:41:46.863","+01:12:04.01"],
                            "[HMT98] AS 38-1": ["22:41:46.4","+01:11:53"],
                            'USNO-A2.0 0900-20357327':["22:41:46.4004867624","+01:11:52.195938360"],
                            '[HMT98] AS 38-2':["22:41:50.2","+01:12:44"],
                            'USNO-A2.0 0900-20357673':["22:41:50.2359466680","+01:12:43.248414384"],
                            'SA 114-' :["22:41:38.8","+01:12:14"],
                            'SDSS J224151.80+011229.0' :["22:41:51.8002660032","+01:12:29.034676656"],
                            '2MASS J22414003+0111094':["22:41:40.0293893472","+01:11:09.486159900"],
                            'Companion 3':["22:41:44.660","+01:12:15.45"],
                            'Companion 4':["22:41:26.622","+01:09:52.42"],
                            'Companion 5':["22:41:27.095","+01:10:52.25"],
                            'Companion 6':["22:41:25.597","+01:10:32.05"],
                            'Companion 7':["22:41:46.275","+01:12:29.98"],
                            'Companion 8':["22:41:44.78","+01:12:12.13"]}

            
            comp_ra,comp_dec=[],[]
            self.init_ref_cat = []
            for c in self.companions.keys():
                companion_ra,companion_dec=self.companions[c]
                comp_ra.append(companion_ra),comp_dec.append(companion_dec)

            self.companions_wcs = SkyCoord(ra=comp_ra,dec=comp_dec,unit=(u.hourangle, u.deg),frame='fk5')
            self.star_coords_pix=np.column_stack((self.sources['xcentroid'],self.sources['ycentroid']))
            self.star_coords_wcs=load_wcs_from_file(filename=self.name,coord=self.star_coords_pix)
            self.star_coords_wcs_sky = SkyCoord(ra=self.star_coords_wcs[:,0]*u.deg, dec=self.star_coords_wcs[:,1]*u.deg,frame='fk5')

            # print(catalog[['ra','dec']])
            # print(list(catalog['ra']))
            # self.cat_ra,self.cat_dec = Table(list(catalog['ra']),names=('ra'),),Table(list(catalog['dec']),names=('dec'))
            self.cat_ra,self.cat_dec = [float(ra_str) for ra_str in list(catalog['ra'])],[float(dec_str) for dec_str in list(catalog['dec'])]
            self.catalog_coords_wcs = SkyCoord(ra=self.cat_ra*u.deg,dec=self.cat_dec*u.deg,frame='fk5')
            self.final_standards,counter ={},1

            #find crossover between panstarrs and reference images
            self.max_sep=sep*u.arcsec
            for i in range(len(self.catalog_coords_wcs)):
                self.d2d=self.catalog_coords_wcs[i].separation(self.star_coords_wcs_sky)
                self.index= self.d2d < self.max_sep

                if any(self.d2d[k].arcsec<=5 for k in np.where(self.index==True)[0]):
                    pos=[self.star_coords_pix[j] for j in range(len(self.index)) if self.index[j]!=False]
                    if 'name' in catalog.columns:
                        print(info_b+" "+catalog['name'].iloc[i], f'standard found at {int(pos[0][0]),int(pos[0][1])}')
                        standard_name = catalog['name'].iloc[i]
                    else:
                        print(info_b+f' standard found at {int(pos[0][0]),int(pos[0][1])}')
                        standard_name = f'standard {counter}'
                    catalog_row = catalog.iloc[i]
                    standard_pos = pos[0]
                    self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_bkgsb)                

                    new_ra,new_dec = [],[]
                    for k in range(len(self.catalog_coords_wcs)):
                        if k!=i:
                            new_ra.append(float(catalog['ra'].iloc[k])),new_dec.append(float(catalog['dec'].iloc[k]))
                        else:
                            new_ra.append(0),new_dec.append(0)

                    new_catalog_coords_wcs = SkyCoord(ra=new_ra*u.deg,dec=new_dec*u.deg,frame='fk5')
                    next_closest_d2d = self.catalog_coords_wcs[i].separation(self.companions_wcs)
                    next_ind = next_closest_d2d <= 60*u.arcsec

                    #masking stars that are within 15" of the standrad (positions held in self.companions and need to be updated)
                    for m_ind in range(len(next_ind)):
                        if next_ind[m_ind]==True:
                            mask_ra_wcs,mask_dec_wcs = self.companions_wcs[m_ind].ra,self.companions_wcs.dec[m_ind]
                            mask_ra_pix,mask_dec_pix = wcs_to_pixels(self.path+'bkg_subtracted_science/'+self.sci_img_name,np.column_stack((mask_ra_wcs,mask_dec_wcs)))[0]
                            mask_ra_pix,mask_dec_pix = int(mask_ra_pix),int(mask_dec_pix)

                            x,y = np.arange(0, np.shape(self.sci_bkgsb)[1]),np.arange(0, np.shape(self.sci_bkgsb)[0])
                            mask = (y[np.newaxis,:]-mask_ra_pix)**2 + (x[:,np.newaxis]-mask_dec_pix)**2 <= star_mask**2

                            if np.max(x)==2047:
                                mask_row,mask_col=np.where(mask==True)
                            else:
                                mask_col,mask_row=np.where(mask==True)

                            for l in range(len(mask_row)):
                                self.sci_bkgsb[mask_row[l],mask_col[l]] = abs(self.median)

                            print(info_b+f" {str(catalog['name'].iloc[i])} has a companion within 1 arcmin, masking star at {mask_ra_pix,mask_dec_pix} with circle of radius {star_mask} pixels")

                    fig = plt.figure(figsize=(20,20))
                    plt.title(f"{standard_name} at {standard_pos} Bkg Subtracted With Companions Masked")
                    plt.imshow(self.sci_bkgsb,cmap='gray',vmin=self.vmin, vmax=self.vmax,origin='lower')

                    for s in range(len(pos)):
                        self.radii,self.apertures = range(1,150,1), []
                        x_pos,y_pos= pos[s]

                        for r in self.radii:
                            self.apertures.append(CircularAperture((int(x_pos),int(y_pos)),r))
                        self.phot_table = aperture_photometry(self.sci_bkgsb,self.apertures)
                        self.apertures[30].plot(color='red',lw=1.5)

                        if show_gc!=False:
                            self.app_counts = []
                            for j in range(len(self.radii)):
                                self.app_counts.append(self.phot_table[f'aperture_sum_{j}'])

                            fig,axs = plt.subplots(figsize=(20,20))
                            y_=np.sort(self.app_counts)
                            x_=self.radii
                            locator = mdates.HourLocator(interval=1)
                            locator.MAXTICKS = 100000
                            axs.xaxis.set_minor_locator(locator)
                            axs.set_title(f"Growth Curve for {standard_name} at {x_pos,y_pos}")
                            axs.set_xlabel('Pixel Radius')
                            axs.set_ylabel('Counts')
                            axs.plot(x_,y_)
                            # axs.legend()
                            plt.show()
                            plt.close()
                        plt.close()

                        

                        if self.sci_filt=='u':#choosing aperture r~20 pixels for u and 30 pixels otherwise
                            final_r=25
                        else:
                            final_r=25

                        self.final_phot_table = aperture_photometry(self.sci_bkgsb,self.apertures[final_r],error=self.bkg_err)
                        self.final_counts,self.final_counts_err = float(self.final_phot_table[f'aperture_sum']),float(self.final_phot_table[f'aperture_sum_err']) 
                        self.zp_err = (self.final_counts_err/self.sci_exp_time)/(self.final_counts/self.sci_exp_time)*1.086

                        # print(catalog_row)
                        if catalog_row[f"{self.sci_filt}'"]=='-':
                            print(warn_r+f" {standard_name} has no {self.sci_filt} magnitude in the catalog, skipping...")
                            counter+=1
                            continue
                        self.filt_mag,self.filt_mag_err = float(catalog_row[f"{self.sci_filt}'"]),float(catalog_row[f"sig {self.sci_filt}'"])

                        if ref_only!=True:
                            self.field_mag = -2.5*np.log10((self.final_counts/self.sci_exp_time))
                            self.zerop_err = np.sqrt(self.zp_err**2 + self.filt_mag_err**2)
                            self.all_zerop.append(self.filt_mag - self.field_mag)
                            self.all_zerop_err.append(self.zerop_err)
                            print(info_g+f" Zeropoint of {standard_name}: {np.round(self.filt_mag - self.field_mag,3)} ± {np.round(self.zerop_err,3)}")
                            print(info_g+f" Airmass: {np.round(self.sci_airmass,3)}, SNR: {int(self.final_counts/self.final_counts_err)}")

                            self.final_standards[float(self.rand_nums_string)+counter] = {'zp':self.filt_mag - self.field_mag,'zp_err':self.zerop_err,'airmass':self.sci_airmass,
                                                                                        'cat_row':catalog_row,'cat_name':standard_name,'filt':self.sci_filt,
                                                                                        'SNR': int(self.final_counts/self.final_counts_err),
                                                                                        'sci_field_mag':self.field_mag}

                        else:
                            self.init_ref_cat.append([self.sci_filt,self.filt_mag,self.filt_mag_err,x_pos,y_pos,catalog_row['ra'],catalog_row['dec']])
                        counter+=1

            if len(self.init_ref_cat)>0:
                self.ref_cat = open(self.path+self.args.output+self.sci_obj+'_'+self.sci_filt+'_LT_ref_cat.txt','w+')
                self.ref_cat.write("#filt,mag,magerr,xpos,ypos,RA,DEC"+"\n")
                for elem in self.init_ref_cat:
                    # print(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]},{elem[5]},{elem[6]}"+"\n")
                    self.ref_cat.write(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]},{elem[5]},{elem[6]}"+"\n")
                self.ref_cat.close()

                print(info_g+ f" Written reference catalog for {self.sci_obj} to ",self.path+self.args.output+self.sci_obj+'_'+self.sci_filt+'_LT_ref_cat.txt')

                return self.init_ref_cat

            if len(self.all_zerop)>0:
                self.mean_zerop = np.average(self.all_zerop,weights=1/np.array(self.all_zerop_err))#np.mean(self.all_zerop)
                self.mean_zerop_err = np.mean(self.all_zerop_err)
                self.mean_SNR = np.mean([self.final_standards[i]['SNR'] for i in self.final_standards.keys()])

                # self.final_standards[self.rand_nums_string+"_avg"] = {'zp':self.mean_zerop,'zp_err':self.mean_zerop_err,'airmass':self.sci_airmass,'filt':self.sci_filt,'SNR':self.mean_SNR}
                brightest = [self.final_standards[i] for i in self.final_standards.keys() if all(self.final_standards[i]['zp']>= self.final_standards[j]['zp'] for j in self.final_standards.keys())]
                self.final_standards[self.rand_nums_string+"_avg"] = brightest[0]
                # print(brightest)
            return self.final_standards        
        
        elif calib_sci!=False:
            filt_low = {'g':'SDSS-G','r':'SDSS-R','i':'SDSS-I','z':'SDSS-Z','u':'SDSS-U'}
            self.mean,self.median, self.std = sigma_clipped_stats(self.sci_bkgsb, sigma=5.0)
            self.iraffind= IRAFStarFinder(threshold=abs(starscale*self.std)*30,fwhm=3,roundhi=0.3,minsep_fwhm=1.5)
            self.sources = self.iraffind(self.sci_bkgsb)

            self.star_coords_pix=np.column_stack((self.sources['xcentroid'],self.sources['ycentroid']))
            self.star_coords_wcs=load_wcs_from_file(filename=self.name,coord=self.star_coords_pix)
            self.star_coords_wcs_sky = SkyCoord(ra=self.star_coords_wcs[:,0]*u.deg, dec=self.star_coords_wcs[:,1]*u.deg,frame='fk5')

            if self.sci_filt=='u' and r==None:#choosing aperture r~20 pixels for u and 30 pixels otherwise
                r=30
            elif self.sci_filt!='u' and r==None:
                r=40
            else:
                r=r

            fig = plt.figure(figsize=(20,20))
            self.vmin,self.vmax = visualization.ZScaleInterval().get_limits(self.sci_bkgsb)
            plt.imshow(self.sci_bkgsb,cmap='gray',vmin=self.vmin, vmax=self.vmax,origin='lower')

            self.sci_phot_tables,j = [],0  
            for xy_pos in self.star_coords_pix:
                if all(0+r<int(star_xy)<np.shape(self.sci_bkgsb)[0]-r for star_xy in xy_pos)==True: #excluding stars within r of edge of image
                    pos_aperture = CircularAperture((int(xy_pos[0]),int(xy_pos[1])),r)
                    phot_table = aperture_photometry(self.sci_bkgsb,pos_aperture,error=self.bkg_err)
                    cal_app_plot = pos_aperture.plot(color='blue',lw=1.5)
                    cal_ra, cal_dec = self.star_coords_wcs[j]
                    # print(self.star_coords_wcs[j])
                    self.sci_phot_tables.append([float(phot_table['aperture_sum'])/self.sci_exp_time,float(phot_table['aperture_sum_err'])/self.sci_exp_time,xy_pos[0],xy_pos[1],cal_ra,cal_dec]) #flux, flux_err, x, y, Ra, DEC
                    j+=1

            self.sn_aperture = CircularAperture((int(self.coords_sn_sci_x),int(self.coords_sn_sci_y)),r)
            sn_app_plot = self.sn_aperture.plot(color='yellow',lw=1.8)

            plt.legend([cal_app_plot,sn_app_plot],['Calibration Stars','SN Position'])
            # plt.show()    
            plt.close()

            print(info_g+f" {len(self.sources)} sources found in the reference image")

            self.cal_mags_arr = []
            for i in range(len(self.sci_phot_tables)):
                if not np.isnan(self.sci_phot_tables[i][0]) and self.sci_phot_tables[i][0]>0:
                    self.cal_inst_mag,self.cal_inst_mag_err = -2.5*np.log10(self.sci_phot_tables[i][0]), (self.sci_phot_tables[i][1]/self.sci_phot_tables[i][0])*1.086

                    #calculating ZP from zp_fit_params y=mx+c y=ZP & x=airmass
                    if all(key_ in zp_fit_params.keys() for key_ in ['slope','intercept']):
                        self.slope, self.intercept = zp_fit_params['slope'],zp_fit_params['intercept']
                        self.field_zp = self.slope*self.sci_airmass + self.intercept
                        self.cal_inst_mag+=self.field_zp
                        # self.cal_inst_mag_err
                    else:
                        print(warn_y+f" Best fit parameters for {self.sci_obj} noy passed in correctly, returning instrumental magnitudes")
                else:
                    self.cal_inst_mag,self.cal_inst_mag_err = '-','-'
                
                self.cal_mags_arr.append([self.sci_filt,self.cal_inst_mag,self.cal_inst_mag_err,self.sci_phot_tables[i][2],self.sci_phot_tables[i][3],self.sci_phot_tables[i][4],self.sci_phot_tables[i][5]])

            self.ref_cat = open(self.path+self.args.output+self.sci_obj+'_'+self.sci_filt+'_LT_ref_cat.txt','w+')
            self.ref_cat.write("#filt,mag,magerr,xpos,ypos,RA,DEC"+"\n")
            for elem in self.cal_mags_arr:
                # print(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]},{elem[5]},{elem[6]}"+"\n")
                self.ref_cat.write(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]},{elem[5]},{elem[6]}"+"\n")
            self.ref_cat.close()

            print(info_g+ f" Written reference catalog for {self.sci_obj} to ",self.path+self.args.output+self.sci_obj+'_'+self.sci_filt+'_LT_ref_cat.txt')

            return self.cal_mags_arr  #magnitude, magnitude error, xpos pixel, ypos pixel, RA, DEC



                    
        

