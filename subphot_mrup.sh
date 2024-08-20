#!/bin/bash -l

#see if the time is 8.01am
if [ $(date +%H%M) -gt 0000 ] && [ $(date +%H%M) -lt 1730 ]; then
    # echo -e $INFO_G Logging to /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log | tee /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log
    echo -e $INFO_G "Time is 8.01am"
    
fi
# # export PYTHON_UNBUFFERRED=True
# flight env activate gridware
# # source ~/.bash_profile
# module load apps/swarp/2.41.5/gcc-8.5.0 apps/psfex/3.24.1/gcc-8.5.0 apps/sextractor/2.28.0/gcc-8.5.0
# conda activate subphot
# INFO_G="\033[0;32m[INFO]  :: \033[0m"
# # echo -e $INFO_G Logging to /mnt/data1/users/arikhind/phot_data/nightly_routine_logs/$(date +%Y%m%d)/nightly_routine_log_$(date +%Y%m%d%H%M)_bash.log | tee /mnt/data1/users/arikhind/phot_data/nightly_routine_logs/$(date +%Y%m%d)/nightly_routine_log_$(date +%Y%m%d%H%M)_bash.log

# echo

# if [ ! -d /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d) ]; then
#     mkdir /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)
#     mkdir /mnt/data1/users/arikhind/phot_data/photometry_date/$(date -d "yesterday" +%Y%m%d)
#     echo -e $INFO_G Created directory /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d) | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log
# fi
# # echo Logging terminal output to /mnt/data1/users/arikhind/phot_data/nightly_routine_logs/$(date +%Y%m%d)/nightly_routine_log_$(date +%Y%m%d%H%M)_bash.log
# echo
# #download data 
# echo -e $INFO_G Downloading any missing data | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log
# python3 -u /users/arikhind/subphot_pipe/subphot_subtract.py -qdl $(date +%Y%m%d%H%M) 2>&1 | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log
# echo -e $INFO_G Running subtraction | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log
# python3 -u /users/arikhind/subphot_pipe/subphot_subtract.py -f Quicklook/$(date +%Y%m%d) -mrup -swarp -psfex -o by_obs_date -cut -upf -pid $(date +%Y%m%d%H%M) 2>&1 | tee /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log

# #run make_LC
# echo -e $INFO_G Making lightcurves | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log

# #run webpage
# echo -e $INFO_G Updating/making webpages | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log

# #run email
# echo -e $INFO_G Sending email | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log
# python3 -u /users/arikhind/subphot_pipe/subphot_morning_email.py -p JL23A05 JZ21B01 JL23A06 JL23A07 JL23B05 JL23B06 JL24A04 JL24A09 JL24B15 -e K.C.Hinds@2021.ljmu.ac.uk -mlog /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log 2>&1 | tee -a /mnt/data1/users/arikhind/phot_data/morning_rup_logs/$(date -d "yesterday" +%Y%m%d)/mrup_$(date +%Y%m%d%H%M)_bash.log

