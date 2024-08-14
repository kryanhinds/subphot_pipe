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
from typing import Mapping, Optional

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

folder="Quicklook/"+DATE


new_only=True

folder_path = 'Quicklook/'

#down_command = f"python3 {path}lt_subtract.py -qdl current_obs" d578adb5-7d1b-49b0-88ad-bdcce71d1050

#os.sytem(down_command)
fits_names = [f for f in os.listdir(data1_path+folder) if f.endswith(".fits")]

try:
    os.system( f"python3 {path}subphot_subtract.py -mrup -up -cut -log night_log -swarp -psfex")

except Exception as e:
    print('Error with subphot_mrup.py script',e)

'''
try:
    os.system( f"python3 {path}subphot_subtract.py -mrup -up -cut -log night_log -swarp")

except Exception as e:
    print('Error with subphot_mrup.py script',e) 
'''
try:
    os.system(f"python3 {path}subphot_make_webpage.py {DATE}")
except Exception as e:
    print('Error with make_webpage.py script',e)



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

LT_proposals = 'JL23A05 JZ21B01 JL23A06 JL23A07 JL23B05 JL23B06 JL24A04 JL24A09'

try:
    os.system(f'python3 {path}subphot_morning_email.py -p '+LT_proposals+' -e K.C.Hinds@2021.ljmu.ac.uk d.a.perley@ljmu.ac.uk J.L.Wise@2022.ljmu.ac.uk A.M.Bochenek@2023.ljmu.ac.uk')
except Exception as e:
    print('Error with morning email script', e)
try:
#     os.system(f'python3 {path}subphot_morning_email.py -p PL24A05 -e K.C.Hinds@2021.ljmu.ac.uk')
# except Exception as e:
#     print('Error with morning email script', e)


# try:
#     os.system(f'python3 {path}morning_email.py -p PL23A11 PL22A17 PL22A12 PL23A13 -e K.C.Hinds@2021.ljmu.ac.uk')
# except Exception as e:
#     print('Error with morning email script', e)




