#!/home/arikhind/miniconda3/envs/ltsub/bin python3
'''#!/Users/kryanhinds/opt/miniconda/envs/ltsubphot/bin python3'''
import os
from subphot_credentials import *
from subphot_functions import *

import sys


DAY = subphot_data().down_quicklook()

os.system(f"python3 {path}subphot_make_webpage.py {DAY}")

#down_command = f"python3 {path}lt_subtract.py -qdl current_obs"
#os.system(f"python3 {path}make_webpage.py {DATE}")
#os.sytem(down_command)



#lt_data().down_archive(names=['GRB190829A'])

#if ret_method == 'archive':
#    lt_data().down_archive(names=[event])


    
        



