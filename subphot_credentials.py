

import requests
from collections import OrderedDict
import os
# directory that the code will be run from
cwd='/users/arikhind/subphot_pipe/'
#cwd= os.getcwd()+'/'
path='/users/arikhind/subphot_pipe/'
data1_path='/mnt/data1/users/arikhind/phot_data/'
token="549e219d-b0eb-4a19-98ee-3a658af2303b"

# for sending reduced data to at the end of the day




# outlook acount for reductions
outlook_account = {"usr":"arikhind@ljmu.ac.uk","pswd":"Avatar99-izzy"}
email_to = ["K.C.Hinds@2021.ljmu.ac.uk","d.a.perley@ljmu.ac.uk","j.l.wise@2022.ljmu.ac.uk"]



#######################################
# also requires python modues: astropy, scipy, photutils, matplotlib, numpy, image_registration, https://github.com/keflavich/image_registration
#astroquery,urllib2,beautifulsoup4,pytest, lacosmic,scikit-image
swarp_path='/opt/apps/alces/swarp/2.41.5/gcc-8.5.0/bin/swarp'
#swarp_path = '/home/arikhind/sn_project1/LT_subtracted_photometry/swarp-2.41.5'
sex_path='/opt/apps/alces/sextractor/2.28.0/gcc-8.5.0/bin/sex'
sexpath='/opt/apps/alces/sextractor/2.28.0/gcc-8.5.0/bin/sex'
psfex_path='/opt/apps/alces/psfex/3.24.1/gcc-8.5.0/bin/psfex'
#panstamps_path='/home/arikhind/.local-u/bin/panstamps'
panstamps_path='/users/arikhind/miniconda3/envs/subphot/bin/panstamps'
scamp_path='/usr/local/bin/scamp'
solve_field_path='/opt/homebrew/Cellar/astrometry-net/0.94_2/bin/solve-field'
solve_field_config_path='/opt/homebrew/Cellar/astrometry-net/0.94_2/etc/astrometry.cfg'

# set the size of the cutout. Occasionally there is bug. The bug means that
# sometimes you get an error
#ValueError: operands could not be broadcast together with shapes (1500,1490) (1500,1500)
# So have to reduce it a little
image_size=1500
# when matching stars to PS1/SDSS, set parameter for finding stars above std of background and search radius to match stars from image to catalog in arcseconds
# set how many times above std of background to find stars
starscale=1.5
# set search radius for PS  matching in arcseconds
search_rad=1

store_lc_ims = True

# PI: Luke Harvey
# Login ID:       PL22A17
# Password:       5735592

# PI: Maxime Deckers
# Login ID:		PL23A12
# Password:       9250074

# PI: Jacco Terwel
# Login ID:		PL23A13
# Password:       3388985

# PI: Georgios Dimitriadis
# Login ID:		PL23A11
# Password:       3411559

proposals_arc={'JL23A05':['8797590'],
               'JL22B10':['2165898'],
               'JL23B05':['5390986'],
               'JL23A06':['2435478'],
               'JZ21B01':['4987478'],
               'JL23A07':['7696732'],
               'JL23B06':['7586251'],
               'JL24A04':['8509523'],
               'PL24A05':['9499033'],# Georgios Dimitriadis
               'JL24A09':['9645313'],
               'JL24B15':['1047928'],
               }
               #'JL21A12':['1092370'],
               #'JL21A14':['4009837'],
               #'JL21B15':['8415098'],

               #'JT20B04a':['1370535'],
               #'JT20B04b':['2849417'],
               #'JT20B05a':['4147333'],
               #'JT20B05b':['4739528'],
               # 'JQ22B02':['15296895'],
                # 'PQ22B01':['4012951'],
            #    'PQ22B02':['3807989']}

