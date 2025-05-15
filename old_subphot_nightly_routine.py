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

time_log = data1_path+'night_logs/'+DATE+'_'+str(TIME).replace(':','.')+'.log'
job_name = 'SPNR_'+DATE+'_'+str(TIME).replace(':','.')

print(info_g+f" Running nightly routine for {DATE} at {TIME}")
print(info_g+' Storing log in '+time_log)
print(info_g+f' Downloading new data for {DATE} to {data1_path+folder}')
DAY = subphot_data().down_quicklook()
print(info_g+' Downloaded data for '+DATE)


print(info_g+' Creating bash script for subtraction routine for '+DATE)


# print(time_log)
# print(job_name)

sbatch_script=f'''
#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --threads-per-core=1
#SBATCH --partition=transients
#SBATCH --job-name={job_name}
#SBATCH --output={time_log}


flight env activate gridware
bash
module load apps/swarp/2.41.5/gcc-8.5.0 apps/psfex/3.24.1/gcc-8.5.0 apps/sextractor/2.28.0/gcc-8.5.0
conda activate subphot


python3 ~/subphot_pipe/subphot_subtract.py -f {folder} -o by_obs_date -swarp -psfex -cut -up -new
'''

with open(path+f'subphot_shell_scripts/nightly_routines/SPNR_{DATE}.sh','w') as f:
    f.write(sbatch_script)
    f.close()

print(info_g+' Bash script written to '+path+'subphot_shell_scripts/nightly_routines/SPNR_'+DATE+'.sh')
print(info_g+f' Submitting job {job_name}')
# os.system('sbatch '+path+'subphot_shell_scripts/nightly_routines/SPNR_'+DATE+'.sh')
print(info_g+f' Removing script {path}subphot_shell_scripts/nightly_routines/SPNR_{DATE}.sh')

# os.remove(path+f'subphot_shell_scripts/nightly_routines/SPNR_{DATE}.sh')