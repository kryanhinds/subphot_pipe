
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

    # NOT/ALFOSC imaging (TELESCOP='NOT'); RA/DEC keywords are in decimal degrees.
    # No MJD keyword — the pipeline computes it from DATE-OBS when 'mjd' is '-'.
    'NOT':                { 'filter':'SEQID','object':'OBJECT','ra':'RA','dec':'DEC','airmass':'AIRMASS','utstart':'-','propid':'PROPID','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'-','seeing':'-','est_seeing':'-','pixscale':0.2139,'gain':'GAIN','bkg_med':'-','bkg_mean':'-'},

    # LCOGT 1m/2m network (TELESCOP varies per site e.g. '1m0-01'; normalised to
    # 'LCOGT' via ORIGIN='LCOGT'). BANZAI .fz files are flattened on input.
    'LCOGT':              { 'filter':'FILTER','object':'OBJECT','ra':'CAT-RA','dec':'CAT-DEC','airmass':'AIRMASS','utstart':'UTSTART','propid':'PROPID','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'MJD-OBS','seeing':'L1FWHM','est_seeing':'-','pixscale':'PIXSCALE','gain':'GAIN','bkg_med':'L1MEDIAN','bkg_mean':'L1MEAN'},

    # Lowell Discovery Telescope / LMI (TELESCOP='DCT', normalised to 'LDT');
    # SCALE keyword gives the binned pixel scale.
    'LDT':                { 'filter':'FILTER','object':'OBJECT','ra':'OBJRA','dec':'OBJDEC','airmass':'-','utstart':'-','propid':'-','instrument':'INSTRUME',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'MJD-OBS','seeing':'-','est_seeing':'-','pixscale':'SCALE','gain':'GAIN','bkg_med':'-','bkg_mean':'-'},

    # Lulin One-metre Telescope / SOPHIA (TELESCOP='LOT'); RA/DEC are
    # space-separated sexagesimal (normalised to colons by the pipeline).
    'LOT':                { 'filter':'FILTER','object':'OBJECT','ra':'RA','dec':'DEC','airmass':'AIRMASS','utstart':'UT','propid':'-','instrument':'-',
                            'exptime':'EXPTIME','date':'DATE-OBS','mjd':'MJD-OBS','seeing':'FWHM','est_seeing':'-','pixscale':0.3843,'gain':'GAIN','bkg_med':'-','bkg_mean':'-'},


                            }

SEDM = ['SEDM-P60','P60','SEDM','sedm-p60','p60','sedm','SEDM-p60','sedm-P60','60']

# ---------------------------------------------------------------------------
# Object-name canonicalisation
# ---------------------------------------------------------------------------
# FITS OBJECT values arrive in several shapes:
#     'ZTF26aakjzdt  g'    SEDM: name + filter token
#     'ACQ-ZTF26abgbpyv r' SEDM acquisition frame
#     'AT 2026fgk'         TJO: a name that genuinely contains a space
#     'GRB260310A'         NOT/LT
# The pipeline previously handled these in five different places with three
# different rules — split('-')[1] (truncates names containing a hyphen),
# split('ACQ-')[1], and split(' ')[0] (which turned 'AT 2026fgk' into 'AT'
# and 'GRB 260310A' into 'GRB', so those epochs were written under the wrong
# object entirely).  One function now does it everywhere.
_FILTER_TOKENS = {'u','g','r','i','z','y',
                  'up','gp','rp','ip','zp','ys',
                  'U','B','V','R','I','J','H','K',
                  'sdssu','sdssg','sdssr','sdssi','sdssz'}


def clean_object_name(raw, return_filter=False):
    """Canonical object name from a FITS OBJECT/TARGET value.

    Removes a leading 'ACQ-' acquisition prefix (case-insensitive, once) and a
    trailing filter token, then joins any remaining whitespace so the same
    target is named identically across facilities ('AT 2026fgk' -> 'AT2026fgk',
    matching what LOT already writes).

    return_filter=True also returns the filter token that was stripped (SEDM
    encodes the filter in OBJECT), or None.
    """
    s = str(raw if raw is not None else '').strip()
    if s.upper().startswith('ACQ-'):
        s = s[4:].strip()
    parts = s.split()
    filt = None
    if len(parts) > 1 and parts[-1] in _FILTER_TOKENS:
        filt = parts[-1]
        parts = parts[:-1]
    return (''.join(parts), filt) if return_filter else ''.join(parts)



# ---------------------------------------------------------------------------
# Marshal request-id suffix
# ---------------------------------------------------------------------------
# Triggered observations can arrive as '<name>_<request id>', and the pipeline
# used to drop everything after the first '_'.  EP and GRB pipeline sources are
# named 'EP-YYMMDD_HHMMSS' / 'GRB-YYMMDD_HHMMSS' on Fritz, where the '_HHMMSS'
# part IS the object ID, so 'EP-261003_001444' was being cut to 'EP-261003',
# which does not exist on Fritz.
import re as _re
_TIMESTAMPED_NAME = _re.compile(r'^\s*(?:ACQ-)?(?:EP|GRB)-?\d{6}_\d{6}(?!\d)',
                                _re.IGNORECASE)


def split_request_id(raw):
    """Split a FITS OBJECT value into (name, sep, request_id), like str.partition('_').

    EP-YYMMDD_HHMMSS and GRB-YYMMDD_HHMMSS names keep their '_HHMMSS' suffix;
    only an underscore *after* that suffix is treated as a request id.
    All other names behave exactly as str.partition('_').
    """
    s = str(raw if raw is not None else '')
    m = _TIMESTAMPED_NAME.match(s)
    if m is None:
        return s.partition('_')
    head, rest = s[:m.end()], s[m.end():]
    rest_name, sep, tail = rest.partition('_')
    return head + rest_name, sep, tail
