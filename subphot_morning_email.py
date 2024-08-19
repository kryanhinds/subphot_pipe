
# '''#!/Users/kryanhinds/opt/miniconda/envs/ltsubphot/bin python3'''
import numpy as np
import glob,os,sys,requests,datetime,smtplib
from email.mime.text import MIMEText
from tabulate import tabulate
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from subphot_credentials import *
import sys
import re
import pandas as pd
from astropy.io import fits
from prettytable import PrettyTable
from termcolor import colored
import argparse
from subphot_functions import info_g,warn_r,warn_y

parser = argparse.ArgumentParser()

parser.add_argument('--date','-d',default='TODAY',help="Date to process")
parser.add_argument('--proposals','-p',default=['all'],help="Proposals to process",nargs='+')
parser.add_argument('--emails','-e',default=['all'],help="Email to send to, default is to look in credentials.py",nargs='+')
parser.add_argument('--mrup_log','-mlog',default=None,help="Log file for morning rup")
args = parser.parse_args()

def query_recent_data(DATE,proposals=['all'],email_to=email_to):
    DATE = '20240815'

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
  
    if proposals==['all']:
        collect_proposals = [key for key in proposals_arc.keys()]
    else:
        collect_proposals = proposals

   
    #getting photometry file names for date specified
    names=[]
    prop_names=[] #the names of objects in the proposals specified
    # print(os.path.exists(f"{path}/photometry_date/{DATE}/morning_rup"))
    # sys.exit()
    if os.path.exists(f"{data1_path}/photometry_date/{DATE}/morning_rup"):
        names = pd.DataFrame([event_name.split("_",1)[0] for event_name in os.listdir(f"{data1_path}photometry_date/{DATE}/morning_rup")],columns=['name'])
        names = np.asarray(names.drop_duplicates())

    data={}
    all_data_phot = pd.DataFrame(columns=['ztfid','filter','mjd','mag','mag_err','prop_id','lim_mag','stack','ra','dec','see',])
    all_data_spec = pd.DataFrame(columns=['ztfid','mjd','ra','dec','propid'])
    measured_phot=False
    if len(names)==0:
        counter=0
        if os.path.exists(data1_path+f'Quicklook/{DATE}'):
            fits_names = [file_ for file_ in os.listdir(data1_path+f'Quicklook/{DATE}') if file_.endswith('.fits')]
            if len(fits_names)>0:
                for k in fits_names:
                    sci_img_hdu = fits.open(data1_path+f'Quicklook/{DATE}/{k}')
                    sci_head = sci_img_hdu[0].header
                    sci_obj,sci_filt,sci_mjd,sci_prop,sci_ra,sci_dec,sci_see = sci_head['OBJECT'], sci_head['FILTER1'],sci_head['MJD'],sci_head['PROPID'],sci_head['CAT-RA'],sci_head['CAT-DEC'],sci_head['L1SEESEC']
      
                    if sci_prop not in collect_proposals:continue

                    # print(all_data_phot.columns)
                    # if counter!=0:print(all_data_phot.iloc[-1])
                    # print(sci_obj,sci_filt,sci_mjd,sci_prop,sci_ra,sci_dec,sci_see)
                    prop_names.append(sci_obj)
                    counter+=1
                    sci_img_name=sci_obj+'_'+sci_filt+sci_head['DATE-OBS'][:-13]+'_'+str(datetime.timedelta(hours=int(sci_head['DATE-OBS'][11:13]), minutes=int(sci_head['DATE-OBS'][14:16]), seconds=float(sci_head['DATE-OBS'][17:21])).seconds)+'.fits'


    if len(names)!=0:
        counter=0
        measured_phot=True
        for n in range(len(names)):
            name=names[n][0]
            today_phot_files = []
            today_phot_files = [x for x in os.listdir(f'{data1_path}photometry_date/{DATE}/morning_rup')]
            today_phot_data = pd.DataFrame(columns=['ztfid','filter','mjd','mag','mag_err','prop_id','stack','ra','dec','see'])

            for f in range(len(today_phot_files)):

                d1 = open(f'{data1_path}photometry_date/{DATE}/morning_rup/{today_phot_files[f]}', 'r')
                data1 = np.asarray(d1.readlines()[0].split(' ')[:-1])
 
             
                today_fits_info = re.sub("photometry.txt","phot_info.txt",today_phot_files[f])
                if os.path.exists(f'{data1_path}phot_fits_info/{today_fits_info}'):pass
                else:today_fits_info = re.sub("photometry.txt","stacked_phot_info.txt",today_phot_files[f])
             
                d2 = open(f'{data1_path}phot_fits_info/{today_fits_info}', 'r')
                data2 = np.asarray(d2.readlines()[0].split(','))


                ztfid,filt,mjd,utcstart,mag,mag_err,lim_mag,ra,dec = data1[0],data1[1],data1[2],data2[1],data1[3],data1[4],data1[5],data2[2],data2[3]
                exp_t,propid,inst,airmass,see,est_see,no_stacked = data2[6],data2[4],data2[5],data2[7],data2[8],data2[9],data2[10]


                if propid not in collect_proposals:continue
                prop_names.append(ztfid)
                today_phot_data.loc[f]=ztfid,filt.replace('-','').lower(),mjd,mag,mag_err,propid,no_stacked,ra,dec,see
                all_data_phot.loc[counter]=ztfid,filt,mjd,mag,mag_err,propid,lim_mag,no_stacked,ra,dec,see
                counter+=1
             

         
            data[f'{name}'] = {'data':today_phot_data,
                            'g_lc':f'{data1_path}light_curves/{name}/{name}_g_LC.png',
                            'r_lc':f'{data1_path}light_curves/{name}/{name}_r_LC.png',
                            'i_lc':f'{data1_path}light_curves/{name}/{name}_i_LC.png',
                            'z_lc':f'{data1_path}light_curves/{name}/{name}_z_LC.png',
                            'u_lc':f'{data1_path}light_curves/{name}/{name}_u_LC.png',}


    
    objects=[]
    spec=[]
    spec_objects=[]
    [objects.append(NAME) for NAME in data.keys()]
    if os.path.exists(f"{data1_path}Quicklook/{DATE}/spec"):
        all_spectra = os.listdir(f"{data1_path}Quicklook/{DATE}/spec")
        spec = [all_spectra[s] for s in range(len(all_spectra)) if all_spectra[s].startswith(f"v_e_{DATE}")]
      

    for t in range(len(spec)):
        spec_fits = fits.open(f"{data1_path}Quicklook/{DATE}/spec/{spec[t]}")
        name_spec,ra_spec,dec_spec,propid_spec,mjd_spec = spec_fits[0].header['OBJECT'], spec_fits[0].header['RA'], spec_fits[0].header['DEC'], spec_fits[0].header['PROPID'], spec_fits[0].header['MJD']
        if  propid_spec not in collect_proposals:
            continue
        all_data_spec.loc[t]=name_spec,mjd_spec,ra_spec,dec_spec,propid_spec

        spec_objects.append(name_spec)


    all_phot_fits = os.listdir(f"{data1_path}Quicklook/{DATE}")
    all_data = pd.DataFrame(columns=['fits','ztfid','mjd','filter','prop_id','ra','dec','see'])
    count=0


    for u in range(len(all_phot_fits)):
        if all_phot_fits[u].endswith(".fits")==True:
            fit=all_phot_fits[u]
            img=fits.open(f"{data1_path}Quicklook/{DATE}/{fit}")
            ztfid,mjd,filt,prop_id,catra,catdec,see = img[0].header['OBJECT'], img[0].header['MJD'], img[0].header['FILTER1'].replace('-','').lower(), img[0].header['PROPID'], img[0].header['CAT-RA'], img[0].header['CAT-DEC'], img[0].header['L1SEESEC']
            if prop_id not in collect_proposals:
                continue
            all_data.loc[count]=fit,ztfid,mjd,filt,prop_id,catra,catdec,see
            count+=1

    if len(spec)>0:
        spec_ztf_urls = [f"https://fritz.science/source/{NAME}" for NAME in data.keys()]
   

    s = smtplib.SMTP('mail.smtp2go.com',587)
    s.set_debuglevel(1)
      
    body1=''

    if len(spec_objects)>0:
        spec_objects= np.unique(spec_objects)
        body1="<p>Objects observed with SPRAT\n</p>"
        for j in range(len(spec_objects)):
            body1=body1+f"<p>{spec_objects[j]}\n</p>"

    body1=body1+"\n"
   
    if len(objects)>0:
        body1=body1+"<p>All objects observed with IO:O:\n</p>"
        body1=body1+"\n"
        all_data = all_data.sort_values(by=['ztfid'])
        all_data = all_data.drop_duplicates()
        tabular_fields1 = ["ZTFID","Filter","RA","DEC", "MJD","Prop ID","Seeing ",]
        tabular_table1 = PrettyTable()
        tabular_table1.field_names = tabular_fields1

        for i in range(len(all_data)):
            tabular_table1.add_row([all_data['ztfid'].iloc[i], 
                                    all_data['filter'].iloc[i],
                                    all_data['ra'].iloc[i], 
                                    all_data['dec'].iloc[i], 
                                    np.round(float(all_data['mjd'].iloc[i]),3), 
                                    all_data['prop_id'].iloc[i], 
                                    all_data['see'].iloc[i]])
        tabular_table1_html = tabular_table1.get_html_string()
    else:
        print(info_g+f" No objects observed</p>")
        tabular_table1_html = ""

    body1=body1+'\n ----------------------------------------------------------------------------------------------- \n'
    body2='\n'

    body2=body2+f"<p>Objects with measured photometry {len(data.keys())}:</p>" 
    body2=body2+"\n"


    if len(data.keys())>0:
        for key in data.keys():

            
            if any(key == nam for nam in prop_names) or any(key in nam for nam in prop_names):
                body2=body2+f"<p>{key}\n</p>"
                body2=body2+f"<p>ZTF url: https://fritz.science/source/{key}\n</p>"
                # body2=body2+f"<p>LT photometry url: https://www.astro.ljmu.ac.uk/~arikhind/lt_subtract/{DATE}/{key}.html\n</p>"
                body2=body2+'\n ----------------------------------------------------------------------------------------------- \n'

    if len(all_data_phot)>0 and measured_phot!=False:
        body2=body2+"<p>All photometry from today\n</P"
        all_data_phot = all_data_phot.sort_values(by=['ztfid'])
        all_data_phot = all_data_phot.drop_duplicates()
        tabular_fields2 = ["ZTFID", "Filter","MJD","Mag","Mag Err", "Lim Mag", "Prop ID", "Stacked"]
        tabular_table2 = PrettyTable()
        tabular_table2.field_names = tabular_fields2
        for ind,val in enumerate(all_data_phot['ztfid']):
            # [print(all_data_phot[i].iloc[ind],type(all_data_phot[i].iloc[ind]),i ) for i in all_data_phot.columns]
            tabular_table2.add_row([val, 
                                    all_data_phot['filter'].iloc[ind], 
                                    np.round(float(all_data_phot['mjd'].iloc[ind]),3), 
                                    np.round(float(all_data_phot['mag'].iloc[ind]),3), 
                                    np.round(float(all_data_phot['mag_err'].iloc[ind]),3), 
                                    np.round(float(all_data_phot['lim_mag'].iloc[ind]),3), 
                                    all_data_phot['prop_id'].iloc[ind],
                                    all_data_phot['stack'].iloc[ind]])

    
        tabular_table2_html = tabular_table2.get_html_string()

    else:
        tabular_table2_html = ""

    body2=body2+'\n \n \n'
    DATE_ = f"{DATE[0:4]}-{DATE[4:6]}-{DATE[6:8]}"
    subj=str(f"{DATE} Photometry and Spectroscopy")
    msg=MIMEMultipart('alternative')

    # if os.path.exists(f"{data1_path}night_log/{DATE}_night_log.log"):
    if os.path.exists(f"{args.mrup_log}"):
        # print
        night_log = open(f"{args.mrup_log}",'rb')
        # night_log = night_log.readlines()
        # print(type(night_log_read[0]))
        # sys.exit()
        # for l in range(len(night_log_read)):
            # if any(ch in night_log_read[l] for ch in ['[32m','[0m','[31m]','[33m','[36m','[1m']):
                # print(night_log_read[l])
                # night_log_read[l] = night_log_read[l].replace('[32m','').replace('[0m','').replace('[31m]','').replace('[33m','').replace('[36m','').replace('[1m','')
                # print(night_log_read[l])
        # night_log_ = [line.replace(ch,'') for line in night_log.readlines()]
        # night_log = night_log_
        log_attachm = MIMEApplication(night_log.read(),Name=f"{args.mrup_log}")
        log_attachm['Content-Disposition'] = 'attachment; filename="%s"' % f"{DATE}_night_log.txt"
        msg.attach(log_attachm)
        print(info_g+f" Attached mrup log to email")

    # sys.exit()
    whole_mail = ""
    whole_mail+='''
        <head>
        <style>
            table, th, td {
                border: 1px solid black;
                border-collapse: collapse;
            }
            th, td {
                padding: 5px;
                text-align: left;    
            }    
        </style>
        </head>'''
    whole_mail+=body1
    whole_mail+="<br><br>"
    whole_mail+=tabular_table1_html
    whole_mail+="<br><br>"
    whole_mail+=body2
    whole_mail+="<br><br>"
    whole_mail+=tabular_table2_html
    whole_mail = MIMEText(whole_mail,'html')
    msg.attach(whole_mail)

    msg['Subject'] = f'LT {DATE_}'
    msg['From'] = outlook_account["usr"]
    s.ehlo()
    s.starttls()
    s.login(outlook_account["usr"],outlook_account["pswd"])

    for n in range(len(email_to)):
        s.sendmail(outlook_account["usr"], email_to[n], msg.as_string()) 
        print(info_g+f' Morning email sent to {email_to[n]}')


    s.quit()


if args.date == 'TODAY':
    now=datetime.datetime.now()
    timedelta=datetime.timedelta(1)
    night_now=(now-timedelta).strftime("%Y%m%d")
else:
    night_now=args.date

if args.emails == ['all']:
    email_to = email_to
else:
    email_to = args.emails

if args.proposals == ['all']:
    proposals = ['all']
else:
    proposals = []
    proposals_input = args.proposals

    for propid in proposals_input: #checking if the proposal exists in proposals_arc dictionary
        if propid not in proposals_arc.keys():
            print(warn_r+f" Proposal {propid} not found in proposals_arc dictionary")
        else:
            proposals.append(propid)

print(info_g+f" Sending morning email for {night_now} to: \033[1m{', '.join(email_to)}\033[0m")
print(info_g+f" Proposals: \033[1m{', '.join(proposals)}\033[0m")
# '\033[1m'+self.sci_obj+'\033[0m'

query_recent_data(DATE=night_now,proposals=proposals,email_to=email_to)

