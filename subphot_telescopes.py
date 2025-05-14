
header_kw = {
    'Liverpool Telescope':{ 'filter':'FILTER1','object':'OBJECT','ra':'CAT-RA','dec':'CAT-DEC','airmass':'AIRMASS','utstart':'UTSTART','propid':'PROPID','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'MJD','seeing':'L1SEESEC','est_seeing':'SCHEDSEE','pixscale':'CCDSCALE',
                            'gain':'GAIN','bkg_med':'L1MEDIAN','bkg_mean':'L1MEAN'},

    'HCT':                { 'filter':'FILTER','object':'OBJECT','ra':'TARRA','dec':'TARDEC','airmass':'-','utstart':'-','propid':'-','instrument':'HFOSC',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'JD','seeing':'FWHM','est_seeing':'-','pixscale':'SCALE','gain':'GAIN','bkg_med':'-','bkg_mean':'-'},

    'SEDM-P60':           { 'filter':'FILTER','object':'OBJECT','ra':'OBJRA','dec':'OBJDEC','airmass':'AIRMASS','utstart':'UTC','propid':'-','instrument':'SEDM-P60',
                            'exptime':'EXPTIME','date':'UTC','mjd':'MJD_OBS','seeing':'FWHM','est_seeing':'-','pixscale':0.370625,'gain':'GAIN','bkg_med':'-','bkg_mean':'-'},

    'SLT':                { 'filter':'FILTER','object':'OBJECT','ra':'OBJCTRA','dec':'OBJCTDEC','airmass':'AIRMASS','utstart':'UT','propid':'-','instrument':'INSTRUMEN',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'jd','seeing':'FWHM','est_seeing':'-','pixscale':0.76,'gain':1.088,'bkg_med':'-','bkg_mean':'-'},
    
    'GTC-OSIRIS':         { 'filter':'FILTER2','object':'OBJECT','ra':'CAT-RA','dec':'CAT-DEC','airmass':'AIRMASS','utstart':'-', 'propid':'-','instrument':'INSTRUME', 
                            'exptime':'EXPTIME','date':'DATE','mjd':'MJD-OBS','seeing':'-','est_seeing':'-','pixscale':0.125,'gain':'GAIN','bkg_med':'-','bkg_mean':'-'},
                            
    'GTC-HIPERCAM':       { 'filter':'user','object':'OBJECT','ra':'RA','dec':'DEC','airmass':'-','utstart':'TIMSTAMP', 'propid':'-','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE','mjd':'MJDUTC','seeing':'-','est_seeing':'-','pixscale':0.35,'gain':1.2,'bkg_med':'-','bkg_mean':'-'  },

    'TJO':                { 'filter':'FILTER','object':'TARGET','ra':'OBJRA','dec':'OBJDEC','airmass':'-','utstart':'DATE','propid':'-','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE','mjd':'MJD-OBS','seeing':'-','est_seeing':'-','pixscale':0.402214,'gain':'GAIN','bkg_med':'-','bkg_mean':'-' },
                        
    'NOT-NOTcam':         { 'filter':'NCFLTNM2','object':'TCSTGT','ra':'OBJRA','dec':'OBJDEC','airmass':'AIRMASS','utstart':'DATE','propoid':'PROPID','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE','mjd':'-', 'seeing':'-','est_seeing':'-','pixscale':0.23,'gain':'GAIN1','bkg_med':'-','bkg_mean':'-'  },
                            
    'NOT-ALFOSC-FASU':    { 'filter':'SEQID', 'object':'TCSTGT','ra':'OBJRA','dec':'OBJDEC','airmass':'AIRMASS','utstart':'DATE','propoid':'PROPID','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE','mjd':'calc', 'seeing':'-','est_seeing':'-','pixscale':0.214028,'gain':'GAIN','bkg_med':'-','bkg_mean':'-' },
                            
    'ESO-NTT':          { 'filter':'ESO INS FILT1 NAME','object':'OBJECT','ra':'RA','dec':'DEC','airmass':'ESO TEL AIRM END','utstart':'DATE-OBS','propid':'-','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'MJD-OBS',
                            'seeing':'ESO TEL AMBI FWHM END','est_seeing':'-','pixscale':0.12,'gain':'ESO DET OUT1 GAIN','bkg_med':'-','bkg_mean':'-','-':None},
                            
                            
                            }

SEDM = ['SEDM-P60','P60','SEDM','sedm-p60','p60','sedm','SEDM-p60','sedm-P60','60']

