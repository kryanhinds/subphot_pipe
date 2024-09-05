
'''#!/home/arikhind/miniconda3/envs/ltsub/bin python3'''

from io import StringIO
from tabulate import tabulate
import logging
import concurrent.futures
import sys
for k in sys.path:
    if 'homebrew' in k:
        sys.path.remove(k)
for k in sys.path:
    if 'homebrew' in k:
        sys.path.remove(k)

# sys.path.append('Users/kryanhinds/miniconda/)
# print(sys.path)
import numpy as np
from astropy.io import fits #FITS files handling
import os  #Call commands from outside Python
import re
# import termcolor
from termcolor import colored
from subphot_credentials import *
import datetime
from datetime import date,timedelta
import time
import pandas as pd
from subphot_functions import *
import argparse
from subphot_quicklook_pipe import *
from subphot_telescopes import header_kw,SEDM
from astropy.time import Time
# from mpi4py import MPI

parser = argparse.ArgumentParser()

parser.add_argument('--calib_folders','-cf',default=[], nargs='*',
                    help='List of folders to calibrate (accepts wildcards or '
                    'space-delimited list)')

parser.add_argument('--stands_cat','-stc',default=None,
                    help='Standards catalogue ')

parser.add_argument('--survey','-s',default='SDSS',
                    help='Catalogue origins (i.e. SDSS or PS1)')

parser.add_argument('--folder','-f',default='', nargs='*',
                    help='List of folders to reduce (accepts wildcards or '
                    'space-delimited list)')

parser.add_argument('--email','-e',action='store_true',default=False,
                    help="Send email with photometry results, requires correct credentials in credentials.py (email username and password)")

parser.add_argument('--multipro','-mp',default='inactive',
                    help="Speed up subtraction by utilising multiprocessing and using 5 nodes")

parser.add_argument('--bands','-fb',default=['All'],nargs='+',
                    help="Filters to process, default is All e.g. g,r,i,z,u for SDSS-G, SDSS-R, SDSS-I, SDSS-Z & SDSS-U")

parser.add_argument('--sci_names','-sn',default=['All'],nargs='+',
                    help="Science Names - names of science object to process, default is all found")

parser.add_argument('--new_only','-new',action='store_true',default=False,
                    help="Processing new images only")

parser.add_argument('--make_log','-log',default='0',
                    help="Outputting log, default is off")

parser.add_argument('--output','-o',default='by_name',#'by_obs_date',
                    help="Dest. for output photometry files, default creates a folder for photometry by observation date e.g. for h_e_20220206_***_1_1.fits output photometry dest will be photometry_data/20220206")

parser.add_argument('--ref_img','-refimg',default='auto',
                    help="Name of reference image, default will be to automatically downloading from PS1/SDSS")

parser.add_argument('--ref_cat','-refcat',default='auto',
                    help="Name of reference catalog, default will be to automatically downloading from PS1/SDSS")

parser.add_argument('--down_date_ql','-qdl',default=[], nargs='+',
                    help='Download data from quicklook, for closest night of observations parse "current_obs" else parse 20220202 for observations on the night of 2022/02/02 ')#(e.g. for 2022/02/02 sunset at La Palma-11:59pm --> observations on the night of 2022/02/02, 12am-sunset at La Palma the next day --> observations on the night of 2022/02/02). Specify another date - note Quicklook only shows the last 7 days with quick reductions')
 
parser.add_argument('--down_date_re','-rdl',default=[], nargs='+',
                    help='Download data from recent data, for closest night of observations parse "current_obs" else parse 20220202 for observations on the night of 2022/02/02 ')# (e.g. for 2022/02/02 sunset at La Palma-11:59pm --> observations on the night of 2022/02/02, 12am-sunset at La Palma the next day --> observations on the night of 2022/02/02). Specify another date - note Recent Data only shows the last 30 days with full reductions')

parser.add_argument('--ims','-i',default='', nargs='+',
                    help='Files to reduce, if stacking then name all in stack but with the first image containing the full path')

parser.add_argument('--cutout','-cut',action='store_true',
                    help='Produce cut outs and deposite in cut_outs, default is False ')

parser.add_argument('--ign_999','-ign999',action='store_true',
                    help='Ignore seeing values ==999 or 966 ')

parser.add_argument('--unsubtract','-un',action='store_true',default=False,
                    help='Unsubtratced option, default is False ')

parser.add_argument('--termoutp','-t',default='normal',
                    help="Type of output to terminal, full/normal/quiet, default is normal")

parser.add_argument('--stack','-stk',action='store_true',default=False,
                    help="Stack where possible, default is False")

parser.add_argument('--upfritz','-up',action='store_true',default=False,
                    help="Upload photometry to Fritz (SkyPortal), default is False. If True, will require correct credentials in credntials.py (Fritz token enabled for uploading)")

parser.add_argument('--upfritz_f','-upf',action='store_true',default=False,
                    help="Replace photometry on Fritz (SkyPortal), default is False. If True, will require correct credentials in credntials.py (Fritz token enabled for uploading)")

parser.add_argument('--cleandirs','-cln',action='store_false',default=True,
                    help="Clean directories by removing all intermediate products (cut outs are not included in this) ")

parser.add_argument('--mroundup','-mrup',action='store_true',default=False,
                    help="Morning round up, multiprocessing subtraction on most recent night of observation looking in data/Quicklook/DATE, default is False ")

parser.add_argument('--zp_only','-zp_only',action='store_true',default=False,
                    help="Calculate the zeropoint only and save to data/zeropoints/by_obs_date")

parser.add_argument('--user_zp_sci','-zp_sci',default=None,
                    help="User supplied zeropoint, default is None")

parser.add_argument('--user_zp_ref','-zp_ref',default=None,
                    help="User supplied zeropoint, default is None")

parser.add_argument('--web_page','-web',action='store_true',default=False,
                    help="Create webpage, default is off")

parser.add_argument('--show_plots','-showp',action='store_true',default=False,
                    help="Make plots pop up")

parser.add_argument('--cal_mode','-cal_mode',default=None,
                    help="Keyword for calibrating science/reference fields")

parser.add_argument('--plc','-lc',default=[None], nargs='*',
                    help="Plots LCs for the given objects, default is off. First arg is folder name, second arg is object name")

parser.add_argument('--plc_space','-lcs',default='mag',
                    help="Space to plot LC's in, default is mag")

parser.add_argument('-relative_flux','-rel',action='store_true',default=False
                    ,help="Plot relative fluxes, default is False")

parser.add_argument('--special_case','-sc',default=None,
                    help="Special case, default is None (e.g. for SN2023ixf where it is bright in a nearby galaxy, we want to use the best bright stars\
                        and they are held in brightstars.cat, so we can use --special_cases brightstars.cat)")

parser.add_argument('--use_swarp','-swarp',action='store_true',default=False,
                    help="Use swarp to stack images and align science image with reference image, default is False")

parser.add_argument('--position','-pos',default='header',nargs='+',
                    help="Position of object in image, default is header, if not in header then parse RA DEC, input as 'HH:MM:SS ±HH:MM:SS")

parser.add_argument('--use_psfex','-psfex',action='store_true',default=False,
                    help="Use psfex to create psf, default is False")

parser.add_argument('--telescope_facility','-tel',default='LT',
                    help="Telescope facility, default is LT. SEDM and HCT are also option")

parser.add_argument('--redo_astrometry','-reastrom',default=False,action='store_true',
                    help="Redo astrometry, default is True, set to False if you want to use the astrometry from the header")

parser.add_argument('--list_fits','-ls',default=None,nargs='+',
                    help="List fits files in directory, and present the name, mjd, seeing, airmass, exptime, filter, RA & DEC")

parser.add_argument('--forced_phot','-fp',default=False,nargs='+',
                    help="Perform forced photometry on the images, default is False. If True, forced photometry performed at CAT-RA and CAT-DEC positions")

parser.add_argument('--use_sdss','-sdss',default=False,action='store_true',
                    help="Use SDSS for reference image, default is False")

parser.add_argument('--redo_batch_astrometry','-rebatch',default=None,nargs='+',
                    help="Redo astrometry for batch, default is False")

parser.add_argument('--pros_job_id','-pid',default=None,
                    help="Prospero job ID, default is None")
args = parser.parse_args()

# print(args.ims)
# # print(args.forced_phot)

# print(args.redo_batch_astrometry)
# print(args.folder)
if args.redo_batch_astrometry!=None:
    redo_list = []

    for n in range(len(args.redo_batch_astrometry)):
        
        [redo_list.append(i) for i in glob.glob(data1_path+args.redo_batch_astrometry[n]+'/*') if i.endswith('.fits')]
    print(redo_list)

    for i in range(len(redo_list)):
        print(info_g+' Redoing astrometry for',redo_list[i],' &', redo_list[i+1])

        subphot_data().batch_update_astrometry([redo_list[i],redo_list[i+1]])#,ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw]) 

# sys.exit(1)
# print(args.sci_names)
if args.list_fits!=None:
    if args.list_fits=='':
        print('Please specify a folder to list fits files in')
    else:
        if args.list_fits[0]=='.':
            # print(info_g+' Listing fits files in current data list_fits listed in subphot_credentials.py')
            args.list_fits[0]=''
        else:
            # print(info_g+' Listing fits files in '+args.list_fits)
            args.list_fits[0]+='/'
        print(info_g+' Listing fits files in '+data1_path+' '.join(args.list_fits))
        for fold_ in args.list_fits:
            arr=[]
            for file_ in os.listdir(data1_path+fold_):
                print(file_)
                if file_.endswith('.fits') and 'tpv' not in file_:
                    hdul = fits.open(data1_path+fold_+file_)
                    hdr = hdul[0].header
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER'],hdr['MJD-OBS'],hdr['SEEING'],hdr['AIRMASS'],hdr['EXPTIME']])
                        continue#,hdr['RA'],hdr['DEC']])
                    except:pass
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER1'],hdr['MJD'],hdr['L1SEESEC'],hdr['AIRMASS'],hdr['EXPTIME'],])
                        continue#hdr['CAT-RA'],hdr['CAT-DEC']])
                    except:pass
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER'],hdr['MJD_OBS'],hdr['FWHM'],hdr['AIRMASS'],hdr['EXPTIME']])
                        continue#,hdr['OBJRA'],hdr['OBJDEC']])
                    except:pass
            df = pd.DataFrame(arr,columns=['IMG','OBJ','FILT','MJD','SEE','AIRM','EXPTIME'])#,'RA','DEC'])

            print(tabulate(df.sort_values(by=['FILT','MJD']), headers='keys', tablefmt='psql'))

    sys.exit(1)

# print(args.use_psfex)
# print(args.use_swarp)
args.telescope_facility = args.telescope_facility.upper()
# sys.exit(1)
seeing_limit =5
if len(args.down_date_ql)>0:
    if args.down_date_ql[0] in ['current_obs' ,'current','most_recent','last_night']:
        ql_day = subphot_data().down_quicklook()
    else:
        for ql_date in args.down_date_ql:
            # if len(ql_date)
            ql_day = subphot_data().down_quicklook(ql_date)

if len(args.down_date_re)>0:
    if args.down_date_re[0] in ['current_obs' ,'current','most_recent','last_night']:
        re_day = subphot_data().down_recentdata()
    else:
        for re_date in args.down_date_re:
            re_day = subphot_data().down_recentdata(re_date)

# sys.exit(1)
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


if time_suns_today<TIME<'23:59:59':
    date_ = TODAY
    DATE = re.sub("-","",date_)
if '00:00:00'<TIME<time_suns_tomorrow:
    date_ = str(t - datetime.timedelta(days=1))
    DATE=re.sub("-","",date_)


if args.mroundup!=False:
    args.folder=[f"Quicklook/"+DATE]
    args.make_log='night_log'
    

try:
    if store_lc_ims and not os.path.exists(data1_path+'entire_lc_imgs'):os.makedirs(data1_path+'entire_lc_imgs')
except:pass


if args.termoutp!='quiet':
    mk_new = True
# config files
if not os.path.exists(path+'config_files'):
    os.makedirs(path+'config_files')

if not os.path.exists(path+'config_files/config.swarp') or mk_new==True:
    writeswarpdefaultconfigfile(args.termoutp)

if not os.path.exists(path+'config_files/config_comb.swarp')  or mk_new==True:
    writeswarpconfigfile(args.termoutp)

if not os.path.exists(path+'config_files/swarp_sdss.conf')  or mk_new==True:
    writesdssswarpconfigfile(args.termoutp)

if not os.path.exists(path+'config_files/default.param')  or mk_new==True:
    writepsfexparfile()

if not os.path.exists(path+'config_files/prepsfex.sex')  or mk_new==True:
    prepsexfile()

if not os.path.exists(path+'config_files/psfex_conf.psfex')  or mk_new==True:
    psfexfile(args.termoutp)

FILTERS = {'g':'SDSS-G','r':'SDSS-R','i':'SDSS-I','z':'SDSS-Z','u':'SDSS-U'}#,'B':'Bessel-B','V':'Bessel-V','R':'Bessel-R','I':'Bessel-I'}
if args.bands==['All']:
    args.bands_to_process=['g','r','i','z','u']#,'B','V','R','I']
else:
    args.bands_to_process=args.bands

# print(args.make_log==0,args.make_log=='0')
# sys.exit(1)

if args.make_log!=False:

    log_dest=args.make_log
    # print(args.output)
    if args.output in ['by_name']:
        # print('yes')
        if args.folder!='':
            for k in range(len(args.folder)):
                # print(args.folder[k])
                if '/' in args.folder[k]:out_dir=args.folder[k].split('/')[-1]
                else:out_dir=args.folder[k]
                if not os.path.exists(data1_path+out_dir):os.mkdir(data1_path+out_dir)
        else:
            out_dir='photometry'
    else:out_dir=args.output.split('/')[-1]


    if args.output!='by_name' and not os.path.exists(data1_path+out_dir):os.mkdir(data1_path+out_dir)

    if not os.path.exists(data1_path+log_dest):os.mkdir(data1_path+log_dest)
    if log_dest=='night_log':log_name = f"{data1_path}night_log/{DATE}_night_log.log"
    else:log_name = f"{data1_path}{log_dest}/{out_dir}_log.log"
    if log_dest=='0':
        log_dest=out_dir
        log_name = f"{data1_path}{out_dir}/{out_dir}_log.log"
    
    if args.pros_job_id!=None:
        if not os.path.exists(f'/mnt/data1/users/arikhind/phot_data/nightly_routine_logs/{DATE}'):os.mkdir(f'/mnt/data1/users/arikhind/phot_data/nightly_routine_logs/{DATE}')
        log_name = f'/mnt/data1/users/arikhind/phot_data/nightly_routine_logs/{DATE}/SPNR_{args.pros_job_id}_py.log'


    if args.mroundup!=False:
        log_name = f"{data1_path}{log_dest}/{DATE}_night_log.log"

    sp_logger = logging
    sp_logger.basicConfig(level=logging.INFO,encoding='utf-8',handlers=[
      logging.StreamHandler(),logging.FileHandler(log_name,encoding='utf-8')],format='%(message)s')
    
    sp_logger.info(info_g+f' Logging to {log_name}')
    args.sp_logger=sp_logger
else:args.sp_logger=None


# sys.exit()
class multi_subtract():
    def __init__(self,folder,folder_path=None):
        #folder_path is the path from data_1path to where the folders are (so that when cron is called, the same folder path can be used from data1_path instead of having to
        # specify the full path from data1_path each time, it is easier for large jobs)
        if folder_path!=None:
            self.folder_path=folder_path
        else:
            folder_path=''
        self.path=path
        self.data1_path=data1_path
        self.g_fits,self.r_fits,self.i_fits,self.z_fits,self.u_fits=[],[],[],[],[]
        self.bess_V_fits,self.bess_R_fits,self.bess_I_fits,self.bess_B_fits = [],[],[],[]
        self.fits_dict = {'SDSS-G':self.g_fits,'SDSS-R':self.r_fits,'SDSS-I':self.i_fits,'SDSS-Z':self.z_fits,'SDSS-U':self.u_fits}
                        # 'Bessell-V':self.bess_V_fits,'Bessell-R':self.bess_R_fits,'Bessell-I':self.bess_I_fits,'Bessell-B':self.bess_B_fits}
        self.filts = {'SDSS-U':'u','SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','SDSS-Z':'z',
                    'Bessell-V':'V','Bessell-R':'R','Bessell-I':'I','Bessell-B':'B',
                    'up_Astrondon_2018':'u','gp_Astrondon_2018':'g','rp_Astrondon_2018':'r','ip_Astrondon_2018':'i','zp_Astrondon_2018':'z',}


        self.folder = folder #specific folder
        self.FOLDER = str(self.folder_path)+'/'+str(self.folder) #full folder path spliced together based off first argument
        self.folder = re.sub(str(self.folder_path),'',self.folder)
        self.folder = re.sub('//','',self.folder)
        self.FOLDER = re.sub('//','/',self.FOLDER)

    def trunc(self,string_,instrument=args.telescope_facility):
        #truncating fits file to identify fits to be stacked and whether the fits is a quicklook reduction (_1_9.fits) or recent data reduction (_1_1.fits)
        if instrument =='LT' or instrument=='lt':
            # print(string_,'_'.join(string_.split('_')[0:4])+'_','_'+'_'.join(string_.split('_')[5:]))
            return ['_'.join(string_.split('_')[0:4])+'_','_'+'_'.join(string_.split('_')[5:])]
        elif instrument=='HCT' or instrument=='hct':
            return [string_,'']
        elif instrument=='SEDM' or instrument=='sedm':
            return [string_,'']
        elif instrument=='SLT' or instrument=='slt':
            return [string_,'']

    def analyse_folder(self):
        self.fits_files,self.all_fits,self.all_in_dir = [],[],[]

        [self.all_in_dir.append(f) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if f.endswith('.fits')]
        # print(self.all_in_dir)
        if args.telescope_facility=='LT' or args.telescope_facility=='lt':
            [self.all_fits.append(self.trunc(f)) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and f.startswith('h_') and self.trunc(f) not in self.all_fits)
            and not f.startswith('.')]
        elif args.telescope_facility=='HCT' or args.telescope_facility=='hct':
            [self.all_fits.append(self.trunc(f,instrument='HCT')) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and self.trunc(f,instrument='HCT') not in self.all_fits) 
            and not f.startswith('.')]
        elif args.telescope_facility=='SEDM' or args.telescope_facility=='sedm':
            [self.all_fits.append(self.trunc(f,instrument='SEDM')) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and self.trunc(f,instrument='SEDM') not in self.all_fits)
            and not f.startswith('.')]
        elif args.telescope_facility=='SLT' or args.telescope_facility=='slt':
            [self.all_fits.append(self.trunc(f,instrument='SLT')) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and self.trunc(f,instrument='SLT') not in self.all_fits)
            and not f.startswith('.')]

        # print(self.all_fits)
        # # print()
        # print(self.all_in_dir)
        # sys.exit(1)

        if len(self.all_fits)==1 and sum([self.all_fits[0][0] in x for x in self.all_in_dir])==0:
            self.fits_files.append([self.all_fits[0],1])
            # print(info_g+f' {self.all_fits[0][0]} is a single image')
            # sys.exit(1)
        else:
            for self.file_arr in self.all_fits:
                self.file,self.end_string = self.file_arr
                self.end_string=re.sub('.fits','',self.end_string)
                self.number = sum([self.file in x for x in os.listdir(f"{self.data1_path}{self.FOLDER}") if not x.startswith('.')])

                if self.number>1:
                    self.file_name,self.fi = [],[]
                    for k in np.linspace(1,self.number,self.number):
                        k=int(k)
                        if k==1:
                            self.file_name.append(f"{self.file}{k}{self.end_string}") 
                            self.fi.append(f"{self.file}{k}{self.end_string}.fits)")
                        if k!=1:
                            self.file_name.append(f"{self.file}{k}{self.end_string}")
                            self.fi.append(f"{self.file}{k}{self.end_string}.fits)")



                    if all(h in os.listdir(f"{self.data1_path}{self.FOLDER}") for h in self.fi)==False:
                        self.fits_files.append([self.file_name,self.number,f"{self.file}1{self.end_string}.fits"])
                else:
                    if self.file.endswith('.fits')==False:
                        self.fits_files.append([[f"{self.file}1{self.end_string}.fits"],1])
                    else:
                        self.fits_files.append([[self.file],1])
                        

        # print(self.fits_files)
        for i in range(len(self.fits_files)):
            if len(self.fits_files[i])==2:
                try:self.fits_hdu = fits.open(self.data1_path+str(self.FOLDER)+"/"+str(self.fits_files[i][0][0]))[0].header
                except:self.fits_hdu = fits.open(self.data1_path+str(self.FOLDER)+"/"+str(self.fits_files[i][0])).header
            else:
                # print(self.fits_files[i][2])
                self.fits_hdu = fits.open(self.data1_path+str(self.FOLDER)+"/"+str(self.fits_files[i][2]))[0].header
            # print(self.fits_hdu['FILTER'])
            try:
                self.fits_filt,self.fits_obj = self.fits_hdu['FILTER1'],self.fits_hdu['OBJECT']
            except:
                # print(warn_y+' No FILTER1 or OBJECT1 in header, trying FILTER and OBJECT')

                try:
                    self.fits_filt,self.fits_obj = self.fits_hdu['FILTER'],self.fits_hdu['OBJECT']
                    if 'p_Astrodon' in self.fits_filt:
                        self.fits_filt = FILTERS[self.fits_filt.split('p_')[0]]

                    # print(info_g+' Found FILTER and OBJECT in header',self.fits_filt,self.fits_obj)
                except:
                    # print(warn_r+' No FILTER1 or OBJECT in header, skipping file')
                    continue 

            if ' ' and self.fits_hdu['TELESCOP']=='60': #specifically for P60
                # print(self.fits_filt,self.fits_obj)
                if 'ACQ-' in self.fits_obj: 
                    self.fits_obj = self.fits_obj.split('-')[1]
                # print(self.fits_obj.split(' '))
                self.fits_obj = self.fits_obj.split(' ')
                self.fits_obj,self._f  = self.fits_obj[0],self.fits_obj[-1]

                self.fits_filt = {'r':'SDSS-R','g':'SDSS-G','i':'SDSS-I','u':'SDSS-U'}[self.fits_filt]
                # print(self.fits_obj,self.fits_filt)

            if args.bands==['All'] and args.sci_names == ['All']:
                FILT = self.fits_filt
                filter_array = self.fits_dict[FILT]
                filter_array.append(self.fits_files[i])
                self.fits_dict[FILT] = filter_array


            elif args.bands!=['All'] and args.sci_names == ['All']:
                for filt_band in args.bands_to_process:
                    FILT = FILTERS[filt_band]

                    if self.fits_filt== FILT:
                        filter_array = self.fits_dict[FILT]
                        filter_array.append(self.fits_files[i])   
                        self.fits_dict[FILT] = filter_array

            elif args.bands==['All'] and args.sci_names != ['All']:
                FILT = self.fits_filt

                if self.fits_obj in args.sci_names or self.fits_obj in args.sci_names:
                    filter_array = self.fits_dict[FILT]
                    filter_array.append(self.fits_files[i])
                    self.fits_dict[FILT] = filter_array
            

            elif args.bands!=['All'] and args.sci_names != ['All']:
                for filt_band in args.bands_to_process:
                    FILT = FILTERS[filt_band]

                    if self.fits_filt == FILT and self.fits_obj in args.sci_names:
                        filter_array = self.fits_dict[FILT]
                        filter_array.append(self.fits_files[i])
                        self.fits_dict[FILT] = filter_array

        # print(self.fits_dict)
        return self.fits_dict 

filt_kws = ['FILTER1','FILTER']
name_kws = ['OBJECT','TARGET','TCSTGT']
date_obs_kws = ['DATE-OBS','DATE','UTC']
def run_subtraction(data_dict):
    final_phot=[]
    '''Function takes a dictonary input of fits files and performs photometry automatically stacking based on the structure of the fits_files array e.g. length 1==single length >1 == stack'''
    filts = {'SDSS-U':'u','SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','SDSS-Z':'z',
            'Bessel-V':'V','Bessel-R':'R','Bessel-I':'I','Bessel-B':'B'}
    f_time_start = time.time()
    fits_files,filter_,FOLDER,new_only = data_dict['fits'],data_dict['filter'],data_dict['FOLDER'],data_dict['new_only']

    sp_logger.info(colored('---------------------------------------------------------------------------------------------','yellow'))
    sp_logger.info(info_g+f" For {filter_}, there are {len(fits_files)} fits")
    

    if len(fits_files)==0:
        sp_logger.info(warn_y+f' No fits files found in filter: {filter_}')
        pass
    else:
        for file_array in fits_files:
            if file_array[1]=='1' or file_array[1]==1:
                fits_file=file_array[0][0]
                sub_file = [re.sub('.fits','',file_array[0][0])]
        
            else:
                fits_file = file_array[2]      
                sub_file = file_array[0]
            
            if fits_file.endswith('.fits') and 'tpv' not in fits_file:# and fits_file.startswith('h_'):
                fits_name = sub_file 
                sub_file[0] = str(FOLDER)+"/"+sub_file[0]
                if data1_path not in sub_file[0]:
                    sub_file[0] = data1_path+sub_file[0]

                
                if args.termoutp!='quiet':
                    if len(sub_file)==1:
                        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
                        sp_logger.info(info_g+f' Performing image subtraction on {sub_file[0]}')
                    else:
                        sp_logger.info(info_g+f" Performing image subtraction on {', '.join(sub_file)}")
                        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
                

                sp_logger.info(info_g+' Extracting header information from '+sub_file[0])
                sci_img = fits.open(f'{data1_path}{FOLDER}/{fits_file}')
                sci_hdr = sci_img[0].header
                for filt_kw in filt_kws:
                    if filt_kw in sci_hdr.keys():
                        filt = sci_hdr[filt_kw]
                        break

                for name_kw in name_kws:
                    if name_kw in sci_hdr.keys():
                        name = sci_hdr[name_kw]
                        break

                for date_obs_kw in date_obs_kws:
                    if date_obs_kw in sci_hdr.keys():
                        date_obs = sci_hdr[date_obs_kw]
                        break
                
                if ' ' and sci_hdr['TELESCOP']=='60': #specifically for P60
                    if 'ACQ-' in name: 
                        name = name.split('-')[1]
                    name = name.split(' ')
                    name,_f = name[0],name[-1]

                    filt = {'r':'SDSS-R','g':'SDSS-G','i':'SDSS-I','u':'SDSS-U'}[filt]

                if store_lc_ims:
                    if not os.path.exists(f'{data1_path}entire_lc_imgs/{name}'):
                        os.makedirs(f'{data1_path}entire_lc_imgs/{name}')
                        sp_logger.info(info_g+' Making folder to store entire light curve for '+name)
                    
                    sp_logger.info(info_g+' Storing image/s for '+name+' in '+data1_path+'entire_lc_imgs/'+name+' if it does not already exist')
                    try:
                        for fit in sub_file:
                            fit_file_name = fit.split('/')[-1]
                            lc_path = f'{data1_path}entire_lc_imgs/{name}/'
                            fit_ = f'{lc_path}{fit_file_name}'
                            if not fit.endswith('.fits'):fit,fit_ = fit+'.fits',fit_+'.fits'
                            if not os.path.exists(fit_):
                                shutil.copy(fit,data1_path+'entire_lc_imgs/'+name) 
                    except Exception as e:
                        pass
                # print(sub_file)
                # print(fits_file)
                # sys.exit()
                # try:filt,name,date_obs = sci_img[0].header['FILTER1'],sci_img[0].header['OBJECT'],sci_img[0].header['DATE-OBS']
                # except:filt,name,date_obs = sci_img[0].header['FILTER'],sci_img[0].header['OBJECT'],sci_img[0].header['DATE-OBS']
                if 'p_Astrodon' in filt:
                    filt = FILTERS[filt.split('p_')[0]]

                if filt in filts.keys():
                    final_name = name+'_'+filts[filt]+date_obs[:-13]+'_'+str(datetime.timedelta(hours=int(date_obs[11:13]), minutes=int(date_obs[14:16]),seconds=float(date_obs[17:31])).seconds)+'_photometry.txt'
                    final_name_stk = re.sub('photometry','stacked_photometry',final_name)
                    # sp_logger.info(info_g+f' Final name of photometry file is: {final_name}')
                    if len(sub_file)>1:
                        args.stack=True

                    if new_only==False:
                        try:
                        # if True:
                            if len(sub_file)>1:
                                # sp_logger.info(sub_file)
                                sp_logger.info(info_g+f' Checking seeing of all images is below {seeing_limit}')
                                sub_file = check_seeing(sub_file)
                                # sp_logger.info(sub_file)
                                # sys.exit(1)
                                if len(sub_file)==0: sp_logger.warning(warn_r+' No images with seeing < 5, passing on image'); continue
                                if len(sub_file)==1: sp_logger.warning(warn_y+' Only one image with seeing < 5, continuing with single image'); args.stack=False

                            sp_logger.info('---------------------------------------------------------------------------------------------')
                            sp_logger.info(info_g+f' Starting sequence on {sub_file[0]}')
                            sub_obj = subtracted_phot(ims=sub_file,args=args)
                            sys_exit=sub_obj.sys_exit
                            # sp_logger.info(sys_exit)

                            if sys_exit==True:
                                pass
                            else:

                                if args.unsubtract==True:
                                    sub_obj.to_subtract=False

                                sys_exit=sub_obj.sys_exit

                                if sys_exit==True:
                                    pass
                                else:

                                    sub_obj.bkg_subtract()
                                    sys_exit=sub_obj.sys_exit
                                    if sys_exit==True:
                                        pass
                                    else:

                                        sub_obj.remove_cosmic()
                                        sys_exit=sub_obj.sys_exit
                                        if sys_exit==True:
                                            pass
                                        else:
                                            
                                            if args.use_swarp==True: sub_obj.swarp_ref_align()
                                            else: sub_obj.py_ref_align()
                                            sys_exit=sub_obj.sys_exit
                                            if sys_exit==True:
                                                pass
                                            else:

                                                if args.use_psfex==True: sub_obj.psfex_convolve_images()
                                                else: sub_obj.py_convolve_images()
                                                sys_exit=sub_obj.sys_exit
                                                if sys_exit==True:
                                                    pass
                                                else:

                                                    sub_obj.gen_ref_cat()
                                                    sys_exit=sub_obj.sys_exit
                                                    if sys_exit==True:
                                                        pass
                                                    else:

                                                        sub_obj.combine_psf()
                                                        sys_exit=sub_obj.sys_exit
                                                        if sys_exit==True:
                                                            pass
                                                        else:

                                                            sub_obj.get_zeropts()
                                                            sys_exit=sub_obj.sys_exit
                                                            if sys_exit==True:
                                                                pass
                                                            else:

                                                                sub_obj.scaled_subtract()
                                                                sys_exit=sub_obj.sys_exit
                                                                if sys_exit==True:
                                                                    pass
                                                                else:

                                                                    final_phot.append(sub_obj.get_photometry())
                                                                    sys_exit=sub_obj.sys_exit
                                                                    if sys_exit==True:
                                                                        pass
                                                                    else:
                                                                        if args.upfritz==True or args.upfritz_f==True:
                                                                            sub_obj.upload_phot()
                                
                                                                        sys_exit=True
                            
                                    if args.cleandirs!=False:
                                        sub_obj.clean_directory()
                        
                        except Exception as e:
                            sp_logger.warning(e)
                            pass
                            



                    elif new_only==True:
                        if (final_name in os.listdir(f'{data1_path}photometry'))==False and (final_name_stk in os.listdir(f'{data1_path}photometry'))==False:
                            if len(sub_file)>1:
                                args.stack =True
                            # print(final_name)
                            # print(final_name in os.listdir(f'{data1_path}photometry'))
                            # print(final_name_stk in os.listdir(f'{data1_path}photometry'))
                            # sys.exit()
                            try:
                                if len(sub_file)>1:
                                    sub_file = check_seeing(sub_file)
                                    if len(sub_file)==0: sp_logger.info(warn_r+' No images with seeing < 5, passing on image'); continue
                                    if len(sub_file)==1: sp_logger.info(warn_y+' Only one image with seeing < 5, continuing with single image'); args.stack=False
                                sub_obj = subtracted_phot(ims=sub_file,args=args)
                                sys_exit=sub_obj.sys_exit

                                if sys_exit==True:
                                    pass
                                else:

                                    if args.unsubtract==True:
                                        sub_obj.to_subtract=False

                                    sys_exit=sub_obj.sys_exit

                                    if sys_exit==True:
                                        pass
                                    else:

                                        sub_obj.bkg_subtract()
                                        sys_exit=sub_obj.sys_exit
                                        if sys_exit==True:
                                            pass
                                        else:

                                            sub_obj.remove_cosmic()
                                            sys_exit=sub_obj.sys_exit
                                            if sys_exit==True:
                                                pass
                                            else:

                                                if args.use_swarp==True: sub_obj.swarp_ref_align()
                                                else: sub_obj.py_ref_align()
                                                sys_exit=sub_obj.sys_exit
                                                if sys_exit==True:
                                                    pass
                                                else:

                                                    if args.use_psfex==True: sub_obj.psfex_convolve_images()
                                                    else: sub_obj.py_convolve_images()
                                                    sys_exit=sub_obj.sys_exit
                                                    if sys_exit==True:
                                                        pass
                                                    else:

                                                        sub_obj.gen_ref_cat()
                                                        sys_exit=sub_obj.sys_exit
                                                        if sys_exit==True:
                                                            pass
                                                        else:

                                                            sub_obj.combine_psf()
                                                            sys_exit=sub_obj.sys_exit
                                                            if sys_exit==True:
                                                                pass
                                                            else:

                                                                sub_obj.get_zeropts()
                                                                sys_exit=sub_obj.sys_exit
                                                                if sys_exit==True:
                                                                    pass
                                                                else:

                                                                    sub_obj.scaled_subtract()
                                                                    sys_exit=sub_obj.sys_exit
                                                                    if sys_exit==True:
                                                                        pass
                                                                    else:

                                                                        final_phot.append(sub_obj.get_photometry())
                                                                        sys_exit=sub_obj.sys_exit
                                                                        if sys_exit==True:
                                                                            pass
                                                                        else:
                                                                            if args.upfritz==True or args.upfritz_f==True:
                                                                                sub_obj.upload_phot()
                                    
                                                                            sys_exit=True
                                                                            
                                        
                                        if args.cleandirs!=False:
                                            sub_obj.clean_directory()
                            except Exception as e:
                                sp_logger.warning(e)
                                pass
                                
                                    

                        else:
                            sp_logger.info(info_b+f' {name} photometry for data taken at {date_obs} in {filter_} already measured')
    

        f_time_end = time.time()
        f_time_total = f_time_end - f_time_start
        if f_time_total > 60:
            f_time_total = f_time_total/60
            time_unit = 'minutes'
        else:
            time_unit = 'seconds'
        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
        sp_logger.info(info_g+f" Finished subtracted photometry in {filter_}, time taken {np.round(f_time_total,2)} {time_unit}")
        #sp_logger.info a table of the photometry
        if len(final_phot)>0: 
            # sp_logger.info(pd.DataFrame(final_phot,columns=final_phot[0].keys()).sort_values(by=['obj','mjd']))    
            sp_logger.info(tabulate(pd.DataFrame(final_phot,columns=final_phot[0].keys()).sort_values(by=['obj','mjd']), headers='keys', tablefmt='psql'))
    return final_phot
            
FILTS=[]



for filt_band in args.bands_to_process:
    FILTS.append(FILTERS[filt_band]) 

if len(args.ims)>0:
    # sp_logger.info('gfefd',args.ims)
    final_phot=[]
    ims = args.ims
    # sp_logger.info(args.ims)
    ims_path = '/'.join(ims[0].split('/')[:-1])
        
    if args.stack==False:
        if '*' in ims[0]:
            ims_start = ims[0].split('*')[0]
            ims=[]
            sp_logger.info(info_g+' Searching for images in '+data1_path+ims_path+' starting with '+ims_start)
            [ims.append(f) for f in os.listdir(data1_path+ims_path) if f.endswith('.fits')==True and ims_start.split('/')[-1] in f and not f.startswith('.')]
            ims = list(np.sort(ims))
            sp_logger.info(info_g+f" Found {len(ims)} images to process")


        # sys.exit()

        # sp_logger.info(ims)
        for i in range(len(ims)):
            
            image = re.sub('.fits','',ims[i])
            if ims_path not in ims[i]:
                image=ims_path+'/'+re.sub('.fits','',ims[i])

            fits_hdu = fits.open(data1_path+image+'.fits')[0].header
            try:
                fits_filt,fits_obj = fits_hdu['FILTER1'],fits_hdu['OBJECT']
                # sp_logger.info(info_g+' Found FILTER1 and OBJECT in header',fits_filt,fits_obj,fits_obj in args.sci_names)
            except:
                sp_logger.info(warn_y+' No FILTER1 or OBJECT in header, trying FILTER and OBJECT')

                try:
                    fits_filt,fits_obj = fits_hdu['FILTER'],fits_hdu['OBJECT']
                    if 'p_Astrodon' in fits_filt:
                        fits_filt = fits_filt.split('p_')[0]
                    sp_logger.info(info_g+' Found FILTER and OBJECT in header '+fits_filt+' '+fits_obj)
                except:
                    sp_logger.info(warn_r+' No FILTER or OBJECT in header, skipping image')
                    try:fits_filt,fits_obj = fits_hdu['FILTERS'],fits_hdu['OBJECT']
                    except:continue


            if args.bands==['All'] and args.sci_names == ['All']:
                pass

            elif args.bands!=['All'] and args.sci_names == ['All']:
                if fits_filt in FILTS:
                    pass
                else:
                    continue
            
            elif args.bands==['All'] and args.sci_names != ['All']:
                if len(args.sci_names)>1 and (fits_obj in args.sci_names or any(x in fits_obj for x in args.sci_names)):
                    pass
                elif len(args.sci_names)==1 and fits_obj in args.sci_names:
                    pass
                else:
                    continue
            
            elif args.bands!=['All'] and args.sci_names != ['All']:
                # sp_logger.info(fits_obj in args.sci_names,fits_obj,fits_filt,FILTS,fits_filt in FILTS)
                # if (fits_obj in args.sci_names or any(x in fits_obj for x in args.sci_names)) and fits_filt in FILTS:
                if len(args.sci_names)>1 and (fits_obj in args.sci_names or any(x in fits_obj for x in args.sci_names)) and fits_filt in FILTS:
                    pass
                elif len(args.sci_names)==1 and fits_obj in args.sci_names and fits_filt in FILTS:
                    pass
                else:
                    continue

            sp_logger.info(colored(f'------------------------------  {i}  ------------------------------','blue'))
            sp_logger.info(info_g+f' Performing image subtraction on {image}, {fits_filt} filter, {fits_obj}')
            # continue
            sub_obj = subtracted_phot(ims=[image],args=args)
            sys_exit=sub_obj.sys_exit


            if sys_exit==True:
                pass
            else:

                if args.unsubtract==True:
                    sub_obj.to_subtract=False

                sys_exit=sub_obj.sys_exit

                if sys_exit==True:
                    pass
                else:

                    sub_obj.bkg_subtract()
                    sys_exit=sub_obj.sys_exit
                    if sys_exit==True:
                        pass
                    else:

                        sub_obj.remove_cosmic()
                        sys_exit=sub_obj.sys_exit
                        if sys_exit==True:
                            pass
                        else:

                            if args.use_swarp==True: sub_obj.swarp_ref_align()
                            else: sub_obj.py_ref_align()
                            sys_exit=sub_obj.sys_exit
                            if sys_exit==True:
                                pass
                            else:

                                if args.use_psfex==True: sub_obj.psfex_convolve_images()
                                else: sub_obj.py_convolve_images()
                                sys_exit=sub_obj.sys_exit
                                if sys_exit==True:
                                   pass
                                else:

                                    if args.relative_flux!=True:
                                        sub_obj.gen_ref_cat()
                                    else:
                                        pass
                                    sys_exit=sub_obj.sys_exit
                                    if sys_exit==True:
                                        pass
                                    else:

                                        sub_obj.combine_psf()
                                        sys_exit=sub_obj.sys_exit
                                        if sys_exit==True:
                                            pass
                                        else:

                                            sub_obj.get_zeropts()
                                            sys_exit=sub_obj.sys_exit
                                            if sys_exit==True:
                                                pass
                                            else:

                                                sub_obj.scaled_subtract()
                                                sys_exit=sub_obj.sys_exit
                                                if sys_exit==True:
                                                    pass
                                                else:

                                                    final_phot.append(sub_obj.get_photometry())
                                                    sys_exit=sub_obj.sys_exit
                                                    if sys_exit==True:
                                                        pass
                                                    else:

                                                        if args.upfritz==True or args.upfritz_f==True:
                                                            sub_obj.upload_phot()
                                                            sys_exit=True
                                                        
                    if args.cleandirs!=False:
                        sub_obj.clean_directory()
                    

        #sp_logger.info a table of the photometry
        if len(final_phot)>0: sp_logger.info(pd.DataFrame(final_phot,columns=final_phot[0].keys()).sort_values(by=['obj','mjd']))                                          

                    
    elif args.stack==True:
        # sp_logger.info(ims)
        if '*' in ims[0]:
            ims_start = ims[0].split('*')[0]
            # sp_logger.info(ims[0])
            # sp_logger.info(ims_start)
            ims=[]
            sp_logger.info(info_g+' Searching for images in '+data1_path+ims_path+' starting with '+ims_start)

            # for f in os.listdir(data1_path+ims_path):
            #     if f.endswith('.fits')==True and ims_start.split('/')[-1] in f:
            # sp_logger.info(glob.glob(data1_path+ims_path+'/*'))
            [ims.append(f) for f in glob.glob(data1_path+ims_path+'/'+ims_start.split('/')[-1]+'*')]#) if f.endswith('.fits')==True ]#and ims_start.split('/')[-1] in f]
            ims = list(np.sort(ims))
            # sp_logger.info(ims)
            # sys.exit(1)
            # []
            sp_logger.info(info_g+f" Found {len(ims)} images to process")

        sp_logger.info(colored(f'---------------------------------------------------------------------------------------------','blue'))
        sp_logger.info(info_g+f' Performing image subtraction on {", ".join(ims)}')
        sp_logger.info('')

        # sp_logger.info(ims)
        ims = check_seeing(ims)
        if len(ims)==0: sp_logger.info(warn_r+' No images with seeing < 5, continuing with single image'); sys.exit()
        if len(ims)==1: sp_logger.info(warn_y+' Only one image with seeing < 5, exiting'); args.stack=False

        # for i in range(len(ims)):
        #     image = re.sub('.fits','',ims[i])
        #     if ims_path not in ims[i]:
        #         image=ims_path+'/'+re.sub('.fits','',ims[i])

        #     fits_hdu = fits.open(data1_path+image+'.fits')[0].header
        #     try:
        #         fits_filt,fits_obj = fits_hdu['FILTER1'],fits_hdu['OBJECT']
        #         # sp_logger.info(info_g+' Found FILTER1 and OBJECT in header',fits_filt,fits_obj,fits_obj in args.sci_names)
        #     except:
        #         sp_logger.info(warn_y+' No FILTER1 or OBJECT in header, trying FILTER and OBJECT')

        #         try:
        #             fits_filt,fits_obj = fits_hdu['FILTER'],fits_hdu['OBJECT']
        #             if 'p_Astrodon' in fits_filt:
        #                 fits_filt = fits_filt.split('p_')[0]
        #             sp_logger.info(info_g+' Found FILTER and OBJECT in header',fits_filt,fits_obj)
        #         except:
        #             sp_logger.info(warn_r+' No FILTER or OBJECT in header, skipping image')
        #             continue



        sub_obj = subtracted_phot(ims=ims,args=args)
        sys_exit=sub_obj.sys_exit

        if sys_exit==True:
            pass
        else:

            if args.unsubtract==True:
                sub_obj.to_subtract=False

            sys_exit=sub_obj.sys_exit

            if sys_exit==True:
                pass
            else:

                sub_obj.bkg_subtract()
                sys_exit=sub_obj.sys_exit
                if sys_exit==True:
                    pass
                else:

                    sub_obj.remove_cosmic()
                    sys_exit=sub_obj.sys_exit
                    if sys_exit==True:
                        pass
                    else:

                        if args.use_swarp==True: sub_obj.swarp_ref_align()
                        else: sub_obj.py_ref_align()
                        sys_exit=sub_obj.sys_exit
                        if sys_exit==True:
                            pass
                        else:

                            if args.use_psfex==True: sub_obj.psfex_convolve_images()
                            else: sub_obj.py_convolve_images()
                            sys_exit=sub_obj.sys_exit
                            if sys_exit==True:
                                pass
                            else:

                                sub_obj.gen_ref_cat()
                                sys_exit=sub_obj.sys_exit
                                if sys_exit==True:
                                    pass
                                else:

                                    sub_obj.combine_psf()
                                    sys_exit=sub_obj.sys_exit
                                    if sys_exit==True:
                                        pass
                                    else:

                                        sub_obj.get_zeropts()
                                        sys_exit=sub_obj.sys_exit
                                        if sys_exit==True:
                                            pass
                                        else:

                                            sub_obj.scaled_subtract()
                                            sys_exit=sub_obj.sys_exit
                                            if sys_exit==True:
                                                pass
                                            else:

                                                final_phot.append(sub_obj.get_photometry())
                                                sys_exit=sub_obj.sys_exit
                                                if sys_exit==True:
                                                    pass
                                                else:

                                                    if args.upfritz==True or args.upfritz_f==True:
                                                        sub_obj.upload_phot()
                                                        sys_exit=True

                if args.cleandirs!=False:
                    sub_obj.clean_directory()

        #sp_logger.info a table of the photometry
        if len(final_phot)>0: sp_logger.info(pd.DataFrame(final_phot,columns=final_phot[0].keys()).sort_values(by=['obj','mjd']))



elif len(args.folder)>0 and args.plc[0]!=[None]: 
    sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
    folder_path = '/'.join(args.folder[0].split('/')[:-1])
    if args.multipro=='inactive':
        data_folders = args.folder
        if args.mroundup==True:
            sp_logger.info(info_g+f" Performing morning round up on photometry taken last night {DATE}")

        all_phot = []
        for d in range(len(data_folders)):
            t_start = time.time()
            sub_folder = re.sub(str(folder_path+'/'),'',data_folders[d])
            multi_sub_obj = multi_subtract(sub_folder,folder_path)
            folder=multi_sub_obj.folder
            FOLDER=multi_sub_obj.FOLDER
            data_dict = multi_sub_obj.analyse_folder()
            # sp_logger.info(data_dict.keys())
            for filt in data_dict.keys():
                filt_dict =  {'fits':data_dict[filt],'filter':filt,'FOLDER':FOLDER,'new_only':args.new_only}
                all_phot.append(run_subtraction(filt_dict))


            t_end=time.time()

            sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
            sp_logger.info(info_g+f" Completed photometry on {sub_folder} in {','.join(args.bands_to_process)} in {np.round(t_end-t_start,2)} seconds ")
            sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
        
        #combine all photometry into one list
        all_phot = [item for sublist in all_phot for item in sublist]
        #sp_logger.info a table of the photometry
        if len(all_phot)>0: 
            # sp_logger.info(pd.DataFrame(all_phot,columns=all_phot[0].keys()).sort_values(by=['obj','mjd']))
            sp_logger.info(tabulate(pd.DataFrame(all_phot,columns=all_phot[0].keys()).sort_values(by=['obj','mjd']), headers='keys', tablefmt='psql'))







    elif args.multipro=='mpi':
        data_folders = args.folder


        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        


        if args.mroundup==True:
            sp_logger.info(info_g+f" Performing multiprocessing on photometry taken last night {DATE}")


        for d in range(len(data_folders)):
            DATA_DICT = {'fits':[]}
            t_start = time.time()
            sub_folder = re.sub(str(folder_path+'/'),'',data_folders[d])
            mp_sub_obj = multi_subtract(sub_folder,folder_path)
            folder=mp_sub_obj.folder
            FOLDER = mp_sub_obj.FOLDER
            data_dict = mp_sub_obj.analyse_folder()

            for i in range(len(args.bands_to_process)):
                if i==rank:
                    phot_band=FILTERS[args.bands_to_process[i]]
                    # sp_logger.info(i,phot_band)
                    sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
                    sp_logger.info(info_g+f' Node {rank} is measuring photometry in {phot_band}')
                    sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))

                    # DATA_DICT = {'fits':data_dict[phot_band],'filter':phot_band,'FOLDER':FOLDER,'new_only':args.new_only}

                    if len(DATA_DICT['fits'])>0:
                        run_subtraction(DATA_DICT)
                        t_end=time.time()
                        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
                        sp_logger.info(info_g+f" Completed photometry in {np.round(t_end-t_start,2)} seconds ")
                        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
                    else:
                        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
                        sp_logger.info(warn_y+f' No fits files found for {phot_band}, passing')
                        sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
            
    
    elif args.multipro=='pools':
        data_folders = args.folder
        if args.mroundup==True:
            sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
            sp_logger.info(info_g+f"Performing multiprocessing morning round up on photometry taken last night {DATE}")
            sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))


        for d in range(len(data_folders)):
            t_start = time.time()
            sub_folder = re.sub(str(folder_path+'/'),'',data_folders[d])
            mp_sub_obj = multi_subtract(sub_folder,folder_path)
            folder=mp_sub_obj.folder
            FOLDER = mp_sub_obj.FOLDER
            data_dict = mp_sub_obj.analyse_folder()

            g_dict = {'fits':data_dict['SDSS-G'],'filter':'SDSS-G','FOLDER':FOLDER,'new_only':args.new_only}
            r_dict = {'fits':data_dict['SDSS-R'],'filter':'SDSS-R','FOLDER':FOLDER,'new_only':args.new_only}
            i_dict = {'fits':data_dict['SDSS-I'],'filter':'SDSS-I','FOLDER':FOLDER,'new_only':args.new_only}
            z_dict = {'fits':data_dict['SDSS-Z'],'filter':'SDSS-Z','FOLDER':FOLDER,'new_only':args.new_only}
            u_dict = {'fits':data_dict['SDSS-U'],'filter':'SDSS-U','FOLDER':FOLDER,'new_only':args.new_only}

            obj_array = np.array([g_dict,r_dict,i_dict,z_dict,u_dict],dtype=object)

            #multiprocessing bit
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(run_subtraction,obj_array)

            t_end=time.time()

            sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))
            sp_logger.info(info_g+f"Completed photometry in {np.round(t_end-t_start,2)} seconds ")
            sp_logger.info(colored('---------------------------------------------------------------------------------------------','blue'))



    if args.pros_job_id!=None:
        if not os.path.exists(data1_path+'nightly_routine_logs/'+DATE):os.mkdir(data1_path+'nightly_routine_logs/'+DATE)
        os.system(f'cp /mnt/data1/users/arikhind/phot_data/nightly_routine_logs/SPNR.log {data1_path}nightly_routine_logs/{DATE}/SPNR_{args.pros_job_id}_scron.log')
        sp_logger.info(info_g+f' Copied log file for job {args.pros_job_id} to {data1_path}nightly_routine_logs/{DATE}/SPNR_{args.pros_job_id}_scron.log')

        # os.system

def load_photometry(folder,name,filters=['All']):

    if filters[0]!='All':
        filter_list = filters
    else:
        filter_list = ['g','r','i','z','u','R']

    all_phot_files = [data1_path+folder+'/'+f for f in os.listdir(data1_path+folder) if re.search(name,f)]
    # all_phot_files = [path+folder+'/'+f for f in os.listdir(path+folder) if name+'_' in f]


    phot_dict = {}
    all_phot_data = []
    for filt in filter_list:
        filt_arr = [f for f in all_phot_files if re.search('_'+filt+'2',f)] #find all files that contain  "_{filter name}2"
        if len(filt_arr)>0:
            for file_ in filt_arr:
                d1 = open(file_, 'r')
                data1 = np.asarray(d1.readlines()[0].split(' ')[:-1])
                all_phot_data.append(data1)

        phot_dict[filt] = all_phot_data
                # sp_logger.info(data1)

    return phot_dict
        # phot_dict[filt] = 




    




if args.plc[0]!=None:
    if len(args.plc)<2:
        sp_logger.info(warn_y+'Please provide the folder and name for a given object in that order')
    
    else:
        sp_logger.info(info_g+' Creating light curve from photometry found in %s for object %s'%(args.plc[0],args.plc[1:]))
        phot_folder = args.plc[0] #folder containing the photometry
        phot_name = [i for i in args.plc[1:] if i!='save'] #name of the object

        for obj_name in phot_name:
            #find all text files in the photometry folder that contain obj_name
            phot_files = [data1_path+f for f in os.listdir(data1_path+phot_folder) if re.search(obj_name,f)]
            sp_logger.info(info_g+f'Found {len(phot_files)} files for {obj_name}')

            filters = {'g':[],'r':[],'i':[],'z':[],'u':[],'R':[],'B':[]} #dictionary to store the photometry for each filter
            for f in phot_files:
                filter_name = f.split(obj_name+'_')[1][0]
                filters[filter_name].append(f)

            #now we have a dictionary of the photometry for each filter
            #we can now plot the light curve
            for key in filters.keys():
                if len(filters[key])>0:
                    sp_logger.info(info_g+f' Plotting light curve for {obj_name} in {key}')
                    # plot_light_curve(filters[key],obj_name,key)
                else:
                    sp_logger.info(warn_y+f' No photometry found for {obj_name} in {key}')
                    continue

                fig,axs = plt.subplots(figsize=(10,10))
                lc_data = load_photometry(phot_folder,obj_name,filters=[key])
                df = pd.DataFrame(lc_data[key],columns=['name','filt','mjd','mag','mag_err','lim_mag','ra_deg','dec_deg','expt','airmass','flux','flux_err'])
                df = df.sort_values(by=['mjd'])

                x = [float(val) for ind,val in enumerate(df['mjd']) if float(df['mag'].iloc[ind])<30]
                

                if args.plc_space in ['mag','m','MAG','M']:
                    y = [float(val) for ind,val in enumerate(df['mag']) if float(val)<30]
                    y_err = [float(val) for ind,val in enumerate(df['mag_err']) if float(df['mag'].iloc[ind])<30]
                    axs.invert_yaxis()
                    axs.set_ylabel('Magnitude')
                else:
                    y = [float(val) for ind,val in enumerate(df['flux']) if df['mag'][ind]<30]
                    y_err = [float(val) for ind,val in enumerate(df['flux_err']) if df['mag'][ind]<30]
                    axs.set_ylabel('Relative Flux')

                # for i in range(len(x)):
                #     if y[i]/y_err[i]>1:
                #         axs.scatter(x[i],y[i],c='blue')
                #     else:
                #         axs.scatter(x[i],y[i],c='k')
                axs.scatter(x,y,c='k')
                axs.errorbar(x,y,yerr=y_err,fmt='o',c='k',capsize=3)
                axs.set_xlabel('MJD')
                axs.set_title(f'Light curve for {obj_name} in {key}')
                # axs.set_xlim(min(x)-0.015,max(x)+0.015)
                # sp_logger.info(x)
                # if axs.set_ylim(-2,max(y)+2.5)
                # axs.axhline(y=0,c='r',ls='--')
                #show the date on the x axis at 45 degree angle
                xlabels = [np.round(float(val),3) for ind,val in enumerate(df['mjd'])]
                axs.set_xticklabels(xlabels,rotation=45,fontsize=10)
                # axs.set_yticks([-1,0,1,2,3,4,5,6],fontsize=11)
                
                # axs.set_xticks(x,rotation=90)
                axs.grid('x','major')
                axs.grid('y')
                # axs.set_xticks(x,rotation=90)
                plt.show()

                if 'save' in args.plc:
                    sp_logger.info(info_g+f' Saving light curve for {obj_name} in {key}')
                    fig.savefig(data1_path+phot_folder+'/'+obj_name+'_'+key+'_light_curve.png',dpi=300)
                    plt.close(fig)

                    sp_logger.info(info_g+f' Saving photometry for {obj_name} in {key}')
                    df.to_csv(data1_path+phot_folder+'/'+obj_name+'_'+key+'_photometry.csv',index=False)


if args.mroundup==True:
    today_phot_files = [f for f in os.listdir(data1_path+'photometry_date/'+DATE) if f!='cut_outs' and f!='morning_rup']

    if os.path.exists(data1_path+'photometry_date/'+DATE+'/morning_rup'):
        today_mrup_phot_files = [f for f in os.listdir(data1_path+'photometry_date/'+DATE+'/morning_rup') if f!='cut_outs' and f!='morning_rup']

        if len(today_mrup_phot_files)!=len(today_phot_files): #finding the files in today_phot_files that aren't in today_mrup_phot_files
            for i in today_phot_files:
                if i not in today_mrup_phot_files:
                    try:
                        os.system('cp '+data1_path+'photometry_date/'+DATE+'/'+i+' '+data1_path+'photometry_date/'+DATE+'/morning_rup/'+i)
                    except Exception as e:
                        print(f'Failed to copy across {i} to photometry_data/{DATE}/morning_rup',e)

    # LT_proposals = 'JL23A05 JZ21B01 JL23A06 JL23A07 JL23B05 JL23B06 JL24A04 JL24A09'

    # morning_logs = glob.glob(data1_path+'morning_rup_logs/'+DATE+'/*')
    # #find the log with DATE in the name
    # for log in morning_logs:
    #     if DATE in log:
    #         mlog = log
    #         break
    # # print(mlog)

    # try:
    #     os.system(f'python3 {path}subphot_morning_email.py -p '+LT_proposals+' -e K.C.Hinds@2021.ljmu.ac.uk -mlog '+mlog)
    #     # d.a.perley@ljmu.ac.uk J.L.Wise@2022.ljmu.ac.uk A.M.Bochenek@2023.ljmu.ac.uk')
    # except Exception as e:
    #     print('Error with morning email script', e)