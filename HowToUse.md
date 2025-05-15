Set path in subphot_credentials.py to the path of your current working directory
  This is where the code will look for config files and any special catalogs (e.g. 2023ixf catalog)

Set data1_path in subphot_credentials.py to the path of your data folder
  All data should be in the same folder, it can be the same as your current working directory but you still need to set the path
  When you run the code, it will look in here for images if -i is used or -f is used
  Outputs will be saved in the same folder as it assumes this has the largest amount of space

Set token in subphot_credentials.py to your Fritz account token
  This is used to upload results to Fritz as well as retrieving the latest photometry

smtp2go_api_key in subphot_credentials.py
  This is used to send emails when the code is finished running in the mornings
  outlook_account is a dictionary containing the email address and password for the account
  smtp2go_api_key is used as it is used to sending emails from via API, you can replace this with your own email address and password if you want
  but you will need to set up the email account to allow less secure apps, which not all email providers allow now

Set swarp_path, sex_path, sexpath, panstamps_path
  These are the paths to the software used in the code
  You can set these to the path of your own installation if you want, but they should be set to the default installation paths
  sexpath is used somewhere but I haven't found where

image_size in alignment

starscale and search_rad used in detecting stars above a threshold in science frame during catalog matching

store_lc_ims when downloading images will give the option to store all the images in a folder associated with the target so that later analysis is easier

proposals_arc is a dictionary containing the proposal IDs for the targets
  This is used to download images from the Quicklook or RecentData
  It isn't possible to download images from the Archive as there's no way to get the credentials for a specific programme since it changes ever two days
  This words via web scraping, essentially it looks for the proposal ID in the raw html file and if it is on the Quicklook or RecentData page, data has been taken so it'll then scrape through the page to get the names of the images and downloads any new images. It is a very messy system, I have since found better ways to search for images but haven't implemented them yet. The functions to do this are in the subphot_functions called down_quicklook and down_recentdata

Basic Usage

  '-s' : selects which reference catalog to use (either SDSS or PS1) but by default it will use PS1
  '-f' : selects which folder(s) to reduce, assuming that the folder(s) are in the data1_path
  '-i' : selects which images to reduce, assuming that the images are in the data1_path, can supply a list of images or wildcards (e.g. h_e_20220206_*.fits which will select all images with that prefix)
  '-stk' : stacks images if there are multiple images provided, only use if you want to stack all images supplied, when a folder is supplied it will stack all images in that folder that *should be* stacked, (e.g. h_e_20220206_*.fits will be stacked but h_e_20220207_1_1_1.fits will not be stacked in there as it is not apart of the same set)
  '-mp' : uses multiprocessing to speed up the reduction, this is currently set to process a band per node, so if you have 5 nodes it will process 5 bands at a time, this is set to inactive by default as it is not always needed and can cause issues with the reduction --NOTE: not used in a while, will either use pools/threads or can use mpi
  '-e' : sends an email when the reduction is finished, this requires the smtp2go_api_key and outlook_account to be set in the subphot_credentials.py
  '-new' : only processes new images, this is set to False by default
  '-b' : selects which bands to process, this is set to all by default (e.g. -b g r which will process g and r bands)
  'sdsscat' : uses SDSS for the reference catalog, this is set to False by default
  '-sn' : selects which science names to process, this is set to all by default (e.g. -sn SN2022abc which will process only that science name)
  '-log' : outputs a log file, this is set to 0 by default (e.g. -log ZTF000000 which will output a log file with that name in the data1_path)
  '-o' : selects the output folder for the photometry files, this is set to by_obs_date by default (e.g. -o by_name which will create a folder for each science name in the data1_path)
  '-refimg' : selects which reference image to use, this is set to auto by default so will download from SDSS/PS1, to supply your own reference image use the relative path to the image (e.g. -refimg ref_imgs/ref_ZTF000000.fits)
  '-refcat' : selects which reference catalog to use, this is set to auto by default so will download from SDSS/PS1, to supply your own reference catalog use the relative path to the catalog (e.g. -refcat ps_catalogs/ref_ZTF000000.cat)
  




  '--ref_img','-refimg',default='auto',
                    help="Name of reference image, default will be to automatically downloading from PS1/SDSS")

  '--ref_cat','-refcat',default='auto',
                    help="Name of reference catalog, default will be to automatically downloading from PS1/SDSS")

  '--down_date_ql','-qdl',default=[], nargs='+',
                    help='Download data from quicklook, for closest night of observations parse "current_obs" else parse 20220202 for observations on the night of 2022/02/02 ')#(e.g. for 2022/02/02 sunset at La Palma-11:59pm --> observations on the night of 2022/02/02, 12am-sunset at La Palma the next day --> observations on the night of 2022/02/02). Specify another date - note Quicklook only shows the last 7 days with quick reductions')
 
  '--down_date_re','-rdl',default=[], nargs='+',
                    help='Download data from recent data, for closest night of observations parse "current_obs" else parse 20220202 for observations on the night of 2022/02/02 ')# (e.g. for 2022/02/02 sunset at La Palma-11:59pm --> observations on the night of 2022/02/02, 12am-sunset at La Palma the next day --> observations on the night of 2022/02/02). Specify another date - note Recent Data only shows the last 30 days with full reductions')

  '--ims','-i',default='', nargs='+',
                    help='Files to reduce, if stacking then name all in stack but with the first image containing the full path')

  '--cutout','-cut',action='store_true',
                    help='Produce cut outs and deposite in cut_outs, default is False ')

  '--ign_999','-ign999',action='store_true',
                    help='Ignore seeing values ==999 or 966 ')

  '--unsubtract','-un',action='store_true',default=False,
                    help='Unsubtratced option, default is False ')

  '--termoutp','-t',default='normal',
                    help="Type of output to terminal, full/normal/quiet, default is normal")

  '--stack','-stk',action='store_true',default=False,
                    help="Stack where possible, default is False")

  '--upfritz','-up',action='store_true',default=False,
                    help="Upload photometry to Fritz (SkyPortal), default is False. If True, will require correct credentials in credntials.py (Fritz token enabled for uploading)")

  '--upfritz_f','-upf',action='store_true',default=False,
                    help="Replace photometry on Fritz (SkyPortal), default is False. If True, will require correct credentials in credntials.py (Fritz token enabled for uploading)")

  '--cleandirs','-cln',action='store_false',default=True,
                    help="Clean directories by removing all intermediate products (cut outs are not included in this) ")

  '--mroundup','-mrup',action='store_true',default=False,
                    help="Morning round up, multiprocessing subtraction on most recent night of observation looking in data/Quicklook/DATE, default is False ")

  '--zp_only','-zp_only',action='store_true',default=False,
                    help="Calculate the zeropoint only and save to data/zeropoints/by_obs_date")

  '--user_zp_sci','-zp_sci',default=None,
                    help="User supplied zeropoint, default is None")

  '--user_zp_ref','-zp_ref',default=None,
                    help="User supplied zeropoint, default is None")

  '--web_page','-web',action='store_true',default=False,
                    help="Create webpage, default is off")

  '--show_plots','-showp',action='store_true',default=False,
                    help="Make plots pop up")

  '--cal_mode','-cal_mode',default=None,
                    help="Keyword for calibrating science/reference fields")

  '--plc','-lc',default=[None], nargs='*',
                    help="Plots LCs for the given objects, default is off. First arg is folder name, second arg is object name")

  '--plc_space','-lcs',default='mag',
                    help="Space to plot LC's in, default is mag")

  '-relative_flux','-rel',action='store_true',default=False
                    ,help="Plot relative fluxes, default is False")

  '--special_case','-sc',default=None,
                    help="Special case, default is None (e.g. for SN2023ixf where it is bright in a nearby galaxy, we want to use the best bright stars\
                        and they are held in brightstars.cat, so we can use --special_cases brightstars.cat)")

  '--use_swarp','-swarp',action='store_true',default=True,
                    help="Use swarp to stack images and align science image with reference image, default is False")

  '--position','-pos',default='header',nargs='+',
                    help="Position of object in image, default is header, if not in header then parse RA DEC, input as 'HH:MM:SS ±HH:MM:SS")

  '--use_psfex','-psfex',action='store_true',default=True,
                    help="Use psfex to create psf, default is False")

  '--telescope_facility','-tel',default='LT',
                    help="Telescope facility, default is LT. SEDM and HCT are also option")

  '--redo_astrometry','-reastrom',default=False,action='store_true',
                    help="Redo astrometry, default is True, set to False if you want to use the astrometry from the header")

  '--list_fits','-ls',default=None,nargs='+',
                    help="List fits files in directory, and present the name, mjd, seeing, airmass, exptime, filter, RA & DEC")

  '--forced_phot','-fp',default=False,nargs='+',
                    help="Perform forced photometry on the images, default is False. If True, forced photometry performed at CAT-RA and CAT-DEC positions")

  '--use_sdss','-sdss',default=False,action='store_true',
                    help="Use SDSS for reference image, default is False")

  '--redo_batch_astrometry','-rebatch',default=None,nargs='+',
                    help="Redo astrometry for batch, default is False")

  '--pros_job_id','-pid',default=None,
                    help="Prospero job ID, default is None")
