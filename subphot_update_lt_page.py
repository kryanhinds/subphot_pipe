import numpy as np
import os
import sys

import astropy
from astropy.io import fits #FITS files handling
from subphot_functions import *
import re
import traceback
import termcolor
from termcolor import colored
import pandas as pd
from datetime import date,datetime
from astropy.time import Time
from subphot_make_LC import *


def make_lt_homepage():
    observationss_path = f"/home/arikhind/public_html/lt_subtract/"
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
    

    obs_list = os.listdir(observations_path)

    with open(f'/home/arikhind/public_html/lt_subtractions_home.html', 'w') as f:
        f.write(preamble)
        



        for observation in obs_list:
            obs = re.sub(".html","",observation)
            row=f'''<a href="/home/arikhind/public_html/lt_subtract/{obs}">{obs} Observations</a> <br>'''
            

            lines.append(row)
    


        lines.append(gen_postamble)
        
        for line in lines:
            f.write(line)
        
        f.close()

    print(f'Webpage for {name} created, displaying {date} photometry')
        
    return
    



make_lt_homepage()







    
    
 