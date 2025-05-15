import numpy as np
from subphot_credentials import *

import sys
import astropy
import os  #Call commands from outside Python
import datetime
from datetime import date,timedelta
from astropy.time import Time
from astroplan import Observer
from subphot_functions import info_g, warn_y, warn_r,subphot_data
import re
import time

jobid = sys.argv[1]
time_log = '/mnt/data1/users/arikhind/phot_data/nightly_routine_logs/SPNR'+jobid+'py.log'
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

folder="Quicklook/"+DATE


print(info_g+f" Running nightly routine for {DATE} at {TIME}")
print(info_g+' Storing log in '+time_log)
print(info_g+f' Downloading new data for {DATE} to {data1_path+folder}')
DAY = subphot_data().down_quicklook()
print(info_g+' Downloaded data for '+DATE)

print(info_g+f' Running /users/arikhind/subphot_pipe/subphot_subtract.py for {DATE} at {TIME} as {jobid}')
# print(f'python3 /users/arikhind/subphot_pipe/subphot_subtract.py -f {folder} -o by_obs_date -swarp -psfex -cut -up -new -pid {jobid}') #jobid e.g. 202408151426 from $(date +%Y%m%d%H%M)

os.system(f'python3 /users/arikhind/subphot_pipe/subphot_subtract.py -f {folder} -o by_obs_date -swarp -psfex -cut -up -new -pid {jobid}')

