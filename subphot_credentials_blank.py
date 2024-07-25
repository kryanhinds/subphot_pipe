# directory that the code will be run from
path='/users/arikhind/subphot_pipe/'

#path to store photometry data in and temporary files (temporary files include ref_imgs - not deleted - background subtracted images, convolved images etc.)
data1_path=''

#Fritz Sky-Portal Token for uploading, deleting and retrieving data
token="549e219d-b0eb-4a19-98ee-3a658af2303b"

# for sending reduced data to at the end of the day




# outlook acount for reductions
outlook_account = {"usr":" ","pswd":" "}
email_to = ["","",""]



swarp_path='/opt/apps/alces/swarp/2.41.5/gcc-8.5.0/bin/swarp'
sex_path='/opt/apps/alces/sextractor/2.28.0/gcc-8.5.0/bin/sex'
sexpath='/opt/apps/alces/sextractor/2.28.0/gcc-8.5.0/bin/sex'
psfex_path='/opt/apps/alces/psfex/3.24.1/gcc-8.5.0/bin/psfex'
panstamps_path='/users/arikhind/miniconda3/envs/subphot/bin/panstamps'

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


proposals_arc={'':['']}

