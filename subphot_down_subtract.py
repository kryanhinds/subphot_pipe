#!/home/arikhind/miniconda3/envs/ltsub/bin python3
'''#!/Users/kryanhinds/opt/miniconda/envs/ltsubphot/bin python3'''
# Import the relevant packages
# Import the relevant packages
import sys
import numpy as np
import astropy
from astropy.io import fits #FITS files handling
import os  #Call commands from outside Python
import re
from termcolor import colored
from subphot_credentials import *
#from LT_quicklook_pipe import lt_subtract
import datetime
from datetime import date,timedelta
from astropy.time import Time
from astroplan import Observer
from subphot_functions import info_g, warn_y, warn_r,subphot_data

import time


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

folder="data/Quicklook/"+DATE


new_only=True

folder_path = 'data/Quicklook/'

DAY = subphot_data().down_quicklook()
print(info_g+' Downloaded data for '+DAY+warn_y)

os.system(f"python3 {path}subphot_make_webpage.py {DAY}")

#down_command = f"python3 {path}lt_subtract.py -qdl current_obs"
#os.system(f"python3 {path}make_webpage.py {DATE}")
#os.sytem(down_command)
fits_names = [f for f in os.listdir(data1_path+folder) if f.endswith(".fits")]


if len(fits_names)>0:
    sub_command = f"python3 {path}subphot_subtract.py -f {folder} -up -new -cut -log night_log -swarp -psfex -o by_obs_date"
    os.system(sub_command)
    # print(sub_command)










