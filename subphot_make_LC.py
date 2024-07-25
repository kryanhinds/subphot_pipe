# Import the relevant packages
import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import os
from subphot_credentials import *

import pandas as pd

# Plot
import matplotlib.pyplot as plt
import mechanize
from bs4 import BeautifulSoup
# import urllib2 
import cookiejar ## http.cookiejar in python3
import requests
import wget
from typing import Mapping, Optional
import urllib.parse
import json

import time
import astropy
from astropy.time import Time
from termcolor import colored
import re
from subphot_functions import *
#ignore warnings
import warnings
import argparse

parser = argparse.ArgumentParser(description='Make light curves for a given supernova')
parser.add_argument('--name', '-sn', nargs='+', help='Name of supernova', required=True)
args = parser.parse_args()
warnings.filterwarnings('ignore')

def plot_mag_LC(name,axs):
    filters = {}
    ztf_filts = ['ztfg','ztfr','ztfi','ztfz']
    lt_filts = ['sdssg','sdssr','sdssi','sdssz','sdssu']
    sloan_filts=lt_filts
    #not_ztf_cols = ['uvot::v','uvot::uvw2','uvot::uvw1','uvot::uvm2','uvot::u','uvot::b']
    
    graph_cols = {'ztfg':'darkgreen','ztfr':'darkred', 'ztfi':'darkblue', 'ztfz':'orange',
                 'sdssg':'green','sdssr':'red','sdssi':'blue','sdssz':'orange','sdssu':'black',
                 'uvot::v':'blue','uvot::uvw2':'red','uvot::uvw1':'black','uvot::uvm2':'yellow','uvot::u':'green','uvot::b':'purple',
                 'SDSSg':'green','SDSSr':'red','SDSSi':'blue','SDSSz':'orange','SDSSu':'black','atlaso':'orange','atlasc':'lime'}
    
    graph_marks = {'ztfg':'o','ztfr':'o', 'ztfi':'o', 'ztfz':'o',
                 'sdssg':'v','sdssr':'v','sdssi':'v','sdssz':'v','sdssu':'v','other':'x'}
                 
    
    filt_cols = {'ztfg':'g','ztfr':'r', 'ztfi':'i', 'ztfz':'z',
                 'sdssg':'g','sdssr':'r','sdssi':'i','sdssz':'z','sdssu':'u',
                 'uvot::v':'v','uvot::uvw2':'uvw2','uvot::uvw1':'uvw1','uvot::uvm2':'uvm2','uvot::u':'u','uvot::b':'b',
                 'SDSSg':'g','SDSSr':'r','SDSSi':'i','SDSSz':'z','SDSSu':'u'}


    data_dict = {}
    TNS_dict = {'SN2021yja':'ZTF21acaqdee','SN2022ann':'ZTF22aaaihet'}
    if name in TNS_dict.keys():name = TNS_dict[name]
    data = SN_data_phot(f'{name}') #downloading data
    

    x = pd.DataFrame(data=[[val,data['filter'].iloc[ind]] for ind,val in enumerate(data['instrument_name'])],columns=['inst','filter'])
    all_filts = np.asarray(x['filter'].drop_duplicates())


    for c in range(len(all_filts)):
        mag = pd.DataFrame(columns=['mag','flag'])
        T =  pd.DataFrame(columns=['mjd','flag'])
        mag_err = pd.DataFrame(columns=['magerr','flag'])
        instrum = pd.DataFrame(columns=['instrument'])

        data_df = pd.DataFrame(columns=['mag','mjd','magerr','instrument','filter','flag'])
        
        FILT = all_filts[c]
        for x in range(len(data.iloc[:,1])):
            instr = str(data['instrument_name'].iloc[x])
            mag_str = str(data['mag'].iloc[x])
            if mag_str == '' or mag_str == 'None' or mag_str == 'nan':continue #if no magnitude recorded, passing
                
            if data['filter'].iloc[x] == FILT and data['ra'].iloc[x] != 'nan' and str(data['instrument_name'].iloc[x])!='IOO':


                if float(data['magerr'].iloc[x])<float(data['mag'].iloc[x]): #removing duplicate measurements and large errors
                    FLAG='A'
                    data_df.loc[x] = float(data['mag'].iloc[x]),float(data['mjd'].iloc[x]),float(data['magerr'].iloc[x]),instr,FILT,FLAG
            
     
              
        if FILT in sloan_filts:
            FILT=re.sub('sdss','SDSS',FILT)
        data_dict[f'{FILT}'] = data_df

    t_now = Time.now().mjd
    filt_leg ={}
    fig_ind = {'ztfg':0,'sdssg':0,'ztfr':1,'sdssr':1,'ztfi':2,'sdssi':2,'ztfz':3,'sdssz':3,'ztfu':4,'sdssu':4,'other':5}

    for filt in data_dict.keys():
        data_df = data_dict[filt]
        mag,mag_err,T,instrum = data_df['mag'],data_df['magerr'],data_df['mjd'],data_df['instrument']
        #reseting the indexs of the dataframes for ease of use later
        mag.reset_index(drop=True, inplace=True)
        T.reset_index(drop=True, inplace=True)
        mag_err.reset_index(drop=True, inplace=True)
        
        T_fm = []

        if filt in ztf_filts:mark, c, col, label = 'o', fig_ind[filt], graph_cols[filt], filt
        if filt in lt_filts:mark, c, col, label = 'v', fig_ind[filt], graph_cols[filt], filt
        if filt not in ztf_filts and filt not in lt_filts:mark, c, col, label = 'x', fig_ind['other'], graph_cols[filt], filt

        if filt[0]=='S':filt=re.sub('SDSS','sdss',filt)

        if len(mag)>0:
            f = axs[c].errorbar(t_now-np.cfloat(data_df['mjd'].values), np.cfloat(data_df['mag'].values), yerr=np.cfloat(data_df['magerr'].values), color=col,elinewidth=0.7,fmt='o',label=label,marker=mark)
            filters[filt] = f

        
    return data_dict, filters

    



t_now = Time.now().mjd
#path = 'photometry'




columns = ['ztfid','filter','mjd','mag','mag_err','lim_mag','ra_deg','dec_deg','exp_t','flux','flux_err']
all_data_df = pd.DataFrame(columns=columns)
all_data_dict = {}


event_names = args.name
phot_files=os.listdir(f'{path}photometry')
for n in range(len(phot_files)):
    for event_name in event_names:
        if event_name in phot_files[n]:
            with open(f'{path}photometry/{phot_files[n]}', 'r') as photometry:
                phot_data = [i for i in photometry.readlines(0)[0].split(' ') if i not in [' ','',str("\n")]]
                if len(phot_data)!=11:
                    for h in range(len(phot_data),11):
                        phot_data.append(None)
                all_data_df.loc[n] = phot_data
                photometry.close()
            #print(all_data_df)
            #all_data_df['mjd'].iloc[n] = float(all_data_df['mjd'].iloc[n]) - 2400000.5
    
all_data_dict[event_names[0]] = all_data_df

# print(all_data_dict[event_names[0]])
# sys.exit()


# gcount,rcount,icount,zcount,ucount=0,0,0,0,0
keys = [key for key in all_data_dict.keys()]
#print(f'Supernovae with photometry: {keys}')

filters = {}
for event_name in event_names:
    if os.path.exists(f'{path}light_curves/{event_name}')==False:os.mkdir(f'{path}light_curves/{event_name}')
    if os.path.exists(f'/home/arikhind/public_html/light_curves/{event_name}')==False:os.mkdir(f'/home/arikhind/public_html/light_curves/{event_name}')

    all_data_dict[event_name].to_csv(f'{path}phot_tables/{event_name}.txt',index=False)
    all_data_dict[event_name].to_csv(f'{path}light_curves/{event_name}/all_phot_{event_name}.txt',index=False)

    
    #fig = plt.figure(figsize=(10,6))
    fig,axs = plt.subplots(figsize = (40,40), nrows=6)
    fig0,axs[0] = plt.subplots()
    fig1,axs[1] = plt.subplots()
    fig2,axs[2] = plt.subplots()
    fig3,axs[3] = plt.subplots()
    fig4,axs[4] = plt.subplots()
    fig5,axs[5] = plt.subplots()
    #tight layout for all subplots
    # plt.tight_layout()
    # ztf_data = plot_mag_LC(event_name,axs)
    
    try:ztf_data = plot_mag_LC(event_name,axs)
    except Exception as e:print(warn_y+f' Could not retrieve Fritz photometry for {event_name}',e)
    
    # filt_ax = {'ztfg':0,'sdssg':0,'ztfr':1,'sdssr':1,'ztfi':2,'sdssi':2,'ztfz':3,'sdssz':3,'ztfu':4,'sdssu':4,'other':5}
    filt_ax = {'sdssg':[0,'g','s'],
                'sdssr':[1,'r','8'],
                'sdssi':[2,'b','h'],
                'sdssz':[3,'orange','p'],
                'sdssu':[4,'black','*'],
                'other':[5,'purple','x']}

    # filters[]
    for filt in filt_ax.keys():
        filt_mag,filt_mag_err,filt_T = all_data_dict[event_name]['mag'].loc[all_data_dict[event_name]['filter']==filt],all_data_dict[event_name]['mag_err'].loc[all_data_dict[event_name]['filter']==filt],all_data_dict[event_name]['mjd'].loc[all_data_dict[event_name]['filter']==filt]
        
        if len(filt_mag)==0:continue
        
        filt_mag,filt_mag_err,filt_T = np.cfloat(filt_mag),np.cfloat(filt_mag_err),np.cfloat(filt_T)
        filt_mag,filt_mag_err,filt_T = filt_mag[(~np.isnan(filt_mag)) & (filt_mag<=25)],filt_mag_err[(~np.isnan(filt_mag_err)) & (filt_mag<=25)],filt_T[(~np.isnan(filt_T)) & (filt_mag<=25)]
        # print(all_data_dict[event_name].loc[all_data_dict[event_name]['filter']==filt])
        ax_,col,mark = filt_ax[filt]

        if filt in ['sdssg','sdssr','sdssi','sdssz','sdssu']:
            filters[filt] = axs[ax_].errorbar(y=filt_mag,x=t_now-filt_T, yerr=filt_mag_err,color=col,fmt='o',marker=mark,label=filt)
        elif filt=='other':
            instrum = all_data_dict[event_name]['inst'].loc[all_data_dict[event_name]['filter']==filt]
            filters[filt] = axs[ax_].errorbar(y=filt_mag,x=t_now-filt_T, yerr=filt_mag_err,color=col,fmt='o',marker=mark,label=instrum+'/'+filt)



    
    for f,i,fig_ in zip(['g','r','i','z','u','other'],[0,1,2,3,4,5],[fig0,fig1,fig2,fig3,fig4,fig5]):
        try:
            axs[i].invert_yaxis(),axs[i].set_xlabel('Days Ago'),axs[i].invert_xaxis(),axs[i].set_ylabel('Apparent Mag')
            axs[i].legend(loc='upper right', bbox_to_anchor=(1.15, 1))
            axs[i].set_title(f'{event_name} {f}-band Light Curve')
            fig_.savefig(f'{path}light_curves/{event_name}/{event_name}_{f}_LC.png')
            fig_.savefig(f'/home/arikhind/public_html/light_curves/{event_name}/{event_name}_{f}_LC.png')
            print(info_g+f' Light curve for {event_name} {f}-band created and saved to {path}light_curves/{event_name}/{event_name}_{f}_LC.png')
        except:
            pass

    
    print(info_g+f' Light curves for {event_name} created')

