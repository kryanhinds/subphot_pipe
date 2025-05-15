from datetime import date,timedelta,datetime
import numpy as np
import os
import sys
from subphot_credentials import *
import astropy
from astropy.io import fits #FITS files handling
from subphot_functions import *
import re
import traceback
import termcolor
from termcolor import colored
import pandas as pd
from datetime import date
from astropy.time import Time
# from subphot_make_LC import *
import datetime
import shutil

#ignore warnings
import warnings
warnings.filterwarnings('ignore')

        
def write_event_phot_page(name,data,DATE,cutouts):
    '''
    t = date.today()
    TIME = datetime.datetime.now().strftime("%H:%M:%S")
    year,month,dayy = t.strftime("%Y"),t.strftime("%m"),t.strftime("%d")
    today = Time(f'{year}-{month}-{dayy} {TIME}')
    
    apo = Observer.at_site("lapalma")
    sun_set_today = apo.sun_set_time(today, which="nearest") #sun set on day of observing
    time_suns_today = "{0.iso}".format(sun_set_today)[-12:]
    sun_set_tomorrow = apo.sun_set_time(today,which="next")
    time_suns_tomorrow = "{0.iso}".format(sun_set_tomorrow)[-12:]

         
    if time_suns_today<TIME<'23:59:59':
        date_ = t
        DATE = re.sub("-","",date_)
    if '00:00:00'<TIME<time_suns_tomorrow:
        date_ = str(t + datetime.timedelta(days=1))
        DATE = re.sub("-","",date_)
    '''

       

    g_path,r_path,i_path,z_path,u_path = f"/~arikhind/public_html/light_curves/{name}/{name}_g_LC.png",f"/~arikhind/public_html/light_curves/{name}/{name}_r_LC.png",f"/~arikhind/public_html/light_curves/{name}/{name}_i_LC.png",f"/~arikhind/public_html/light_curves/{name}/{name}_z_LC.png",f"/~arikhind/public_html/light_curves/{name}/{name}_u_LC.png"
    g_sci_png = cutouts["g"]["sci"]
    r_sci_png = cutouts["r"]["sci"]
    i_sci_png = cutouts["i"]["sci"]
    z_sci_png = cutouts["z"]["sci"]
    u_sci_png = cutouts["u"]["sci"]

    preamble = '''
    <!DOCTYPE html>
    <html>
    <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <style>
    table {
      font-family: arial, sans-serif;
      border-collapse: collapse;
      width: 100%;
    }

    td, th {
      border: 1px solid #dddddd;
      text-align: left;
      padding: 8px;
    }

    tr:nth-child(even) {
      background-color: #dddddd;
    }
    </style>
    <body>
    '''+f'''
    <h1>{name}</h1>


    <p>Page contains photometry for {name} taken with LT IO:O. Solid squares are data from LT/IO:O, solid circles are data from P48 and x's are other instruments like SEDM/P60. The tables contain the photometry not currently on Fritz.  <p>    <br>
    <a href="https://fritz.science/source/{name}"> {name} Fritz Link </a> <br>
   
    '''

    tab_postamble = '''
    </table>
    '''

    gen_postamble = '''
    </body>
    </html>
    '''
    
    tab_preamble = f'''

    <table>
      <tr>
        <th>mjd</th>
        <th>filter</th>
        <th>mag</th>
        <th>magerr</th>
        <th>limiting_mag</th>
        <th>magsys</th>
     
        
      </tr>
    '''

    trs = '<tr>' #start
    tre = '</tr>' #end 
    ths = '<th>'
    the = '</th>'
    tds = '<td>'
    tde = '</td>'
    

    date_now=datetime.datetime.now()
    lines = []
    t_now = Time.now().mjd

    #data_new=data
    data_new = pd.DataFrame(data.loc[data['new_phot']=='Y'],columns=['ztfid', 'filter', 'mjd', 'mag', 'mag_err', 'lim_mag', 'ra_deg','dec_deg', 'exp_t', 'new_phot', 'magsys'])
    
    if os.path.exists(f'/home/arikhind/public_html/lt_subtract/{DATE}')==False:
        os.mkdir(f'/home/arikhind/public_html/lt_subtract/{DATE}')
    with open(f'/home/arikhind/public_html/lt_subtract/{DATE}/{name}.html', 'w') as f:
        f.write(preamble)
        f.write("<p> Magnitude Light Curves<p><br><br>")

        f.write('''<p> SDSS-G & ZTFG </p>''')     
        f.write(f'''<p>  <img src="{g_sci_png}" alt="New" />  </p><br>''')
        f.write(f'''<p>  <img src="/~arikhind/light_curves/{name}/{name}_g_LC.png" alt="{name} LC g" /> </p><br>''')

        f.write('''<p> SDSS-R & ZTFR </p><br>''')   
        f.write(f'''<p>  <img src="{r_sci_png}" alt="New" />  </p><br>''')
        f.write(f'''<p>  <img src="/~arikhind/light_curves/{name}/{name}_r_LC.png" alt="{name} LC r" /> </p><br>''')

        f.write('''<p> SDSS-I & ZTFI </p><br>''')   
        f.write(f'''<p>  <img src="{i_sci_png}" alt="New" />  </p><br>''')
        f.write(f'''<p>  <img src="/~arikhind/light_curves/{name}/{name}_i_LC.png" alt="{name} LC i" /> </p><br>''')

        f.write('''<p> SDSS-Z & ZTFZ </p><br>''')   
        f.write(f'''<p>  <img src="{z_sci_png}" alt="New" />  </p><br>''')
        f.write(f'''<p>  <img src="/~arikhind/light_curves/{name}/{name}_z_LC.png" alt="{name} LC z" /> </p><br>''')
        try:

            f.write('''<p> SDSS-U </p><br>''')   
            f.write(f'''<p>  <img src="{u_sci_png}" alt="New" />  </p><br>''')
            f.write(f'''<p>  <img src="/~arikhind/light_curves/{name}/{name}_u_LC.png" alt="{name} LC u" />  </p><br>''')
        except:
            pass
        f.write(f'''<p>  <img src="/~arikhind/light_curves/{name}/{name}_other_LC.png" alt="{name} LC Other" /> </p><br>''')
        f.write(f'''<p> As of {date_now} (in mjd {t_now}), the following photonetry was not found on Fritz. </p>''')
        f.write("<p>mjd,filter,mag,magerr,limiting_mag,magsys</p>")
        



        for ind,val in enumerate(data_new['mjd']):
            if float(val)>1000000:
                val = float(val)-2400000.5
                data_new['filter'].iloc[ind] = f"sdss{data_new['filter'].iloc[ind]}"
            mjd,mag,magerr,filt,limiting_mag,magsys = np.round(val,6), np.round(data_new['mag'].iloc[ind],3), np.round(data_new['mag_err'].iloc[ind],3),data_new['filter'].iloc[ind], np.round(data_new['lim_mag'].iloc[ind],3), data_new['magsys'].iloc[ind]
            row=f'''{mjd},{filt},{mag},{magerr},{limiting_mag},{magsys}<br>'''

            lines.append(row)
    


        lines.append(gen_postamble)
        
        for line in lines:
            f.write(line)
        
        f.close()

    print(info_g+f' Webpage for {name} created, displaying {DATE} photometry')
        
    return
    
 
def write_obs_page(DATE):
    observations_path = f"/home/arikhind/public_html/lt_subtract"
    preamble = '''
    <!DOCTYPE html>
    <html>
    <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <style>
    table {
      font-family: arial, sans-serif;
      border-collapse: collapse;
      width: 100%;
    }

    td, th {
      border: 1px solid #dddddd;
      text-align: left;
      padding: 8px;
    }

    tr:nth-child(even) {
      background-color: #dddddd;
    }
    </style>
    <body>
    '''+f'''
    


    <p>Page contains photometry taken on {DATE}  <p>    <br>

    <p> As of {DATE}, the following photonetry was not found on Fritz. </p>
   
    '''

    tab_postamble = '''
    </table>
    '''

    gen_postamble = '''
    </body>
    </html>
    '''
    
    tab_preamble = f'''

    <table>
      <tr>
        <th>mjd</th>
        <th>filter</th>
        <th>mag</th>
        <th>magerr</th>
        <th>limiting_mag</th>
        <th>magsys</th>
     
        
      </tr>
    '''

    trs = '<tr>' #start
    tre = '</tr>' #end 
    ths = '<th>'
    the = '</th>'
    tds = '<td>'
    tde = '</td>'
    


    lines = []
    

    #data_new=data
    events_list = os.listdir(f"{observations_path}/{DATE}")

    
    with open(f'/home/arikhind/public_html/lt_subtract/{DATE}/{DATE}_obs.html', 'w') as f:
        f.write(preamble)
        for event in events_list:
            
            if 'ZTF' in event or 'SN' in event or 'GRB' in event:
            
                NAME = re.sub(".html","",event)
                row=f'''<a href="{event}">{NAME} Observations</a> <br>'''

                lines.append(row)
    


        lines.append(gen_postamble)
        
        for line in lines:
            f.write(line)
        
        f.close()

    print(info_g+f' Webpage for {DATE} observations created')
        
    return

def make_lt_homepage():
    observations_path = f"/home/arikhind/public_html/lt_subtract/"
    preamble = '''
    <!DOCTYPE html>
    <html>
    <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <style>
    table {
      font-family: arial, sans-serif;
      border-collapse: collapse;
      width: 100%;
    }

    td, th {
      border: 1px solid #dddddd;
      text-align: left;
      padding: 8px;
    }

    tr:nth-child(even) {
      background-color: #dddddd;
    }
    </style>
    <body>
    '''+f'''
    <h1>LT Subtracted Photometry Homepage</h1>


    <p>This is the Homepage for Liverpool Telescope image subtraction, containing the latest photometry taken by IO:O. .  <p>    <br>
    
    <p> Available Observations<p>

    '''

    tab_postamble = '''
    </table>
    '''

    gen_postamble = '''
    </body>
    </html>
    '''
    
    tab_preamble = f'''

    <table>
      <tr>
        <th>mjd</th>
        <th>filter</th>
        <th>mag</th>
        <th>magerr</th>
        <th>limiting_mag</th>
        <th>magsys</th>
     
        
      </tr>
    '''

    trs = '<tr>' #start
    tre = '</tr>' #end 
    ths = '<th>'
    the = '</th>'
    tds = '<td>'
    tde = '</td>'
    


    lines = []
    

    obs_list = np.sort(os.listdir(observations_path))[::-1]

    with open(f'/home/arikhind/public_html/lt_subtractions_home.html', 'w') as f:
        f.write(preamble)
        



        for observation in obs_list:
            row=f'''<a href="lt_subtract/{observation}/{observation}_obs.html">{observation} Observations</a> <br>'''
            

            lines.append(row)
    


        lines.append(gen_postamble)
        
        for line in lines:
            f.write(line)
        
        f.close()

    print(info_g+f' Webpage for {name} created, displaying {DATE} photometry')
        
    return
    




date_ = sys.argv[1]
phot_path = 'photometry'

if os.path.exists(path+'photometry_date/'+date_)==True:
    all_phot_data = [f for f in os.listdir(f'{path}photometry_date/{date_}') if f!='morning_rup' and f!='cut_outs']
    if os.path.exists(f'{path}photometry_date/{date_}/morning_rup')==True:
        all_phot_data = [f for f in os.listdir(f'{path}photometry_date/{date_}/morning_rup') if f!='morning_rup' and f!='cut_outs']
else:
    all_phot_data
    print(warn_r+f" No photometry found for {date_}")
    sys.exit()

print(info_g+f" Found photometry for {date_}")
# print(all_phot_data)
# # print(os.listdir(f'{path}photometry_date/{date_}/morning_rup'))
# sys.exit()
DATE = f"{date_[0:4]}-{date_[4:6]}-{date_[6:8]}"
date_now=datetime.datetime.now()

#print('date now',date_now)

timedelta=datetime.timedelta(1)
date_tomo = (date_now+timedelta).strftime("%Y%m%d")
DATE_TOMO = f"{date_tomo[0:4]}-{date_tomo[4:6]}-{date_tomo[6:8]}"

date_yest = (date_now-timedelta).strftime("%Y%m%d")
DATE_YEST = f"{date_yest[0:4]}-{date_yest[4:6]}-{date_yest[6:8]}"
#print('date yesterday',date_yest)   
date_now=date_now.strftime("%Y%m%d")
DATE_NOW = f"{date_now[0:4]}-{date_now[4:6]}-{date_now[6:8]}"
DATES=[DATE_YEST,DATE_NOW,DATE_TOMO]



if os.path.exists(f'/home/arikhind/public_html/lt_subtract/{date}')==False:
    os.mkdir(f'/home/arikhind/public_html/lt_subtract/{date}')

print(all_phot_data)
#getting photometry file names for data specified
names = pd.DataFrame([event_name.split("_",1)[0] for event_name in all_phot_data if "cut" not in event_name],columns=['name'])
names = np.array(names.drop_duplicates())
names=[]
if len(names)==1:
    print(info_g+f" Creating webpage for the following object: {names[0][0]}")
else:
    print(info_g+f' Creating webpages for the following objects: {", ".join([i[0] for i in names])}')
cutouts_dict = {}
# print(len(names))
if len(names)!=0:
    for n in range(len(names)):
        name=names[n][0]
        lc_command = f"python3 {path}subphot_make_LC.py -sn {name}"
        print(info_g+f" Creating light curve for {name} and updating/creating webpage")
        os.system(lc_command)
        data = SN_data_phot(f'{name}')
        
        try:
            data_on_fritz = data[['mjd','mag','magerr','limiting_mag','magsys','instrument_id','filter']]
            lt_on_fritz = pd.DataFrame(data_on_fritz.loc[data_on_fritz['instrument_id']==33])
        except Exception as e:
            data_on_fritz=pd.DataFrame(columns=['mjd','mag','magerr','limiting_mag','magsys','instrument_id','filter'])
            lt_on_fritz = data_on_fritz
    
        try:
            lt_event_data_df = pd.read_csv(f"{path}light_curves/{name}/all_phot_{name}.txt")
            lt_event_data_df['new_phot'] = 'N'
            lt_new_df = pd.DataFrame(columns=lt_event_data_df.columns)
        except Exception as e:
            lt_event_data_df = pd.DataFrame(columns=['ztfid', 'filter', 'mjd', 'mag', 'mag_err', 'lim_mag', 'ra_deg','dec_deg', 'exp_t', 'new_phot', 'magsys'])
            lt_new_df = pd.DataFrame(columns=lt_event_data_df.columns)
            print(warn_r+f" Error reading in photometry for {name} : {e}")
    
        obs_on = np.asarray(lt_on_fritz['mjd'].values)


        for ind,val in enumerate(lt_event_data_df['mjd'].values):

            if np.round(val,6) not in np.round(obs_on,6) and any(obs_on[x]-0.001<val<obs_on[x]+0.001 for x in range(len(obs_on)))==False and lt_event_data_df['mag'].iloc[ind]<90:
                lt_new_df.loc[ind] = lt_event_data_df.loc[ind]
                lt_event_data_df['new_phot'].iloc[ind] = 'Y'
    
        lt_event_data_df['magsys']='ab'
        lt_new_df['magsys']='ab'
        lt_event_data_df = lt_event_data_df.sort_values(by=['mjd'],ascending=False)

        if not os.path.exists(f"/home/arikhind/public_html/light_curves/{name}/cutouts"):
            os.mkdir(f"/home/arikhind/public_html/light_curves/{name}/cutouts")
        if os.path.exists(f"/home/arikhind/public_html/light_curves/{name}/cutouts"):
            os.system("rm "f"/home/arikhind/public_html/light_curves/{name}/cutouts/*.png")
        
        cutouts_dict[name] = {"g":{"sci":"","ref":"","sub":""},"r":{"sci":"","ref":"","sub":""},"i":{"sci":"","ref":"","sub":""},"u":{"sci":"","ref":"","sub":""},"z":{"sci":"","ref":"","sub":""}}
        
        try:
            for fil in cutouts_dict[name].keys():
                new_sci_png,new_ref_png,new_sub_png="","",""
                for co_file in os.listdir(path+f"photometry_date/{date_}/cut_outs"):
                    #print(co_file,str(name+"_"+fil+"20"),str(name+"_"+fil+"20") in co_file)
                    if str(name+"_"+fil+"20") in co_file:
                        if "cutout_panel" in co_file and "stacked" in co_file:
                            sci_png = path+f"photometry_date/{date_}/cut_outs/"+co_file
                            new_sci_png = f"/home/arikhind/public_html/light_curves/{name}/cutouts/"+co_file
                            shutil.copy2(sci_png, new_sci_png)
                            new_sci_png = re.sub("home/","~",new_sci_png)
                            new_sci_png = re.sub("/public_html","",new_sci_png)

                        elif "cutout_panel" in co_file and "stacked" not in co_file:
                            sci_png = path+f"photometry_date/{date_}/cut_outs/"+co_file
                            new_sci_png = f"/home/arikhind/public_html/light_curves/{name}/cutouts/"+co_file
                            shutil.copy2(sci_png, new_sci_png)
                            new_sci_png = re.sub("home/","~",new_sci_png)
                            new_sci_png = re.sub("/public_html","",new_sci_png)
                 

                cutouts_dict[name][fil] = {"sci":new_sci_png}

        except Exception as e:
            print(warn_r+f" Error making webpage for {name} : {e}")
            pass
            



            



        write_event_phot_page(name=name,data=lt_event_data_df,DATE=date_,cutouts=cutouts_dict[name])
 

    write_obs_page(date_)

    make_lt_homepage()

else:
    print(info_b+f" Found no useable photometry for observations on {DATE}")







    
    
