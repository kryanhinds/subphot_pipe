            # print(tabulate(df.sort_values(by=['FILT','MJD']), headers='keys', tablefmt='psql'))

'''#!/home/arikhind/miniconda3/envs/ltsub/bin python3'''

from io import StringIO
from tabulate import tabulate
import logging
import concurrent.futures
import sys
for k in sys.path:
    if 'homebrew' in k:
        sys.path.remove(k)
for k in sys.path:
    if 'homebrew' in k:
        sys.path.remove(k)

# sys.path.append('Users/kryanhinds/miniconda/)
# print(sys.path)
import numpy as np
from astropy.io import fits #FITS files handling
import os  #Call commands from outside Python
import re
import shutil  # For terminal size detection
# import termcolor
from termcolor import colored
from subphot_credentials import *
import datetime
from datetime import date,timedelta
import time
import pandas as pd
from subphot_functions import *
import argparse
from subphot_quicklook_pipe import *
from subphot_telescopes import header_kw,SEDM,clean_object_name
import subphot_runlog as runlog
from astropy.time import Time

# ---------------------------------------------------------------------------
# [SRV] Path hygiene — must run before ANY directory is created.
#
# Three failure modes seen on the minar cron installation:
#   * an unexpanded '~' in path/data1_path makes os.makedirs() create a
#     LITERAL '~' directory tree under the current working directory;
#   * a missing trailing slash silently concatenates ('.../pipeentire_lc_imgs');
#   * a relative path (or tools that write relative to CWD) drops artefacts
#     wherever cron happened to start, which is '/' for most cron daemons.
# Normalising here fixes all three for every downstream 'path + name' string.
# ---------------------------------------------------------------------------
def _norm_root(p, label):
    try:
        _p = str(p)
    except Exception:
        return p
    _q = os.path.abspath(os.path.expanduser(os.path.expandvars(_p)))
    if not _q.endswith(os.sep):
        _q += os.sep
    if _q != _p:
        print(f'[INFO]    :: [SRV] {label} normalised: {_p!r} -> {_q!r}')
    return _q

path = _norm_root(path, 'path')
# An empty data1_path means "inputs are given as full paths, everything else
# lives in the pipeline root" (minar runs this way).  Normalising '' literally
# turns it into whatever directory the process started in, which is '/' under
# cron and in the preflight smoke test: inputs were then rewritten relative to
# '/', re-joined onto the pipeline root (a doubled path that cannot be opened),
# and entire_lc_imgs/, logs and -o output were all attempted under '/'.
if not str(data1_path if data1_path is not None else '').strip():
    print(f'[INFO]    :: [SRV] data1_path is empty, using the pipeline root {path}')
    data1_path = path
data1_path = _norm_root(data1_path, 'data1_path')
# subphot_functions / subphot_quicklook_pipe imported their own copies from
# subphot_credentials; keep them identical to the normalised values
for _mod_name in ('subphot_functions', 'subphot_quicklook_pipe'):
    _mod = sys.modules.get(_mod_name)
    if _mod is not None:
        _mod.path, _mod.data1_path = path, data1_path

# Anchor the process in the pipeline directory.  All user-supplied image and
# folder arguments are already resolved against data1_path (not CWD), so this
# is behaviour-preserving, and it stops SExtractor/SWarp/PSFEx scratch files
# from landing in '/' or '~' when cron starts us elsewhere.
try:
    if os.path.isdir(path) and os.path.realpath(os.getcwd()) != os.path.realpath(path):
        os.chdir(path)
        print(f'[INFO]    :: [SRV] working directory set to {path}')
except Exception as _cd_e:
    print(f'[WARNING] :: [SRV] could not chdir to {path}: {_cd_e}')


def rel_to_data1(p):
    """Express a user-supplied path RELATIVE to data1_path.

    The pipeline concatenates 'data1_path + <arg>' in many places, so an
    absolute argument (which is what a cron wrapper naturally passes, e.g.
    -i /data/sedmdrp/redux/.../image.fits) produced a doubled path like
    '/data/.../subphot_pipe//data/sedmdrp/...' and a FileNotFoundError.
    Absolute paths inside the data root are rewritten as relative; absolute
    paths outside it are returned unchanged and handled by the callers'
    os.path.exists() checks.
    """
    try:
        _p = str(p)
    except Exception:
        return p
    if not os.path.isabs(_p):
        return _p
    _abs = os.path.abspath(os.path.expanduser(_p))
    _root = os.path.abspath(data1_path)
    if _abs.startswith(_root + os.sep):
        return os.path.relpath(_abs, _root)
    return _p


def safe_makedirs(target, label=''):
    """Create a directory, but never outside the pipeline root.

    Guards against the '~'/'/' pollution above: anything resolving outside
    `path`/`data1_path` is refused with a warning instead of being created.
    Uses exist_ok=True so concurrent cron jobs cannot race each other.
    """
    try:
        _t = os.path.abspath(os.path.expanduser(os.path.expandvars(str(target))))
        _roots = [os.path.realpath(path), os.path.realpath(data1_path)]
        _extra = os.environ.get('SUBPHOT_EXTRA_OUTPUT_ROOTS', '')
        _roots += [os.path.realpath(r) for r in _extra.split(':') if r]
        if not any(os.path.realpath(_t).startswith(r) for r in _roots):
            print(f'[WARNING] :: [SRV] refusing to create {label or "directory"} '
                  f'outside the pipeline root: {_t}')
            return None
        os.makedirs(_t, exist_ok=True)
        return _t
    except Exception as _mk_e:
        print(f'[WARNING] :: [SRV] could not create {label or "directory"} {target}: {_mk_e}')
        return None

# ---------------------------------------------------------------------------
# Optional pinned progress bar (-pb/--progress).  See _PinnedBar below: it
# reserves the terminal's last line with a scroll region, which is the only
# approach that survives output from os.system() subprocesses (SWarp,
# SExtractor, PSFEx) that write straight to the terminal.
# ---------------------------------------------------------------------------
_PBAR = None

class _PinnedBar:
    """Progress bar pinned to the terminal's bottom line.

    Uses a DECSTBM scroll region so that EVERYTHING else — pipeline prints,
    logging records on stderr, and the output of os.system() subprocesses such
    as SWarp/SExtractor/PSFEx, which write straight to the terminal and cannot
    be intercepted by redirecting sys.stdout — scrolls in the region above the
    reserved final line.  The bar therefore stays put without needing to
    capture or proxy any of that output.

    Falls back to plain periodic lines when stdout/stderr is not a terminal
    (e.g. redirected to a log file), so batch runs still record progress.
    """

    def __init__(self, total, desc='reducing'):
        self.total, self.n, self.desc = int(total), 0, desc
        self.t0, self.postfix, self.rows = time.time(), '', None
        try:
            self.tty = open('/dev/tty', 'w')
            self.enabled = self.tty.isatty()
        except Exception:
            self.tty, self.enabled = sys.stderr, False
        if self.enabled:
            self._reserve()
            self.draw()

    # -- terminal plumbing ---------------------------------------------------
    def _reserve(self):
        """Reserve the bottom line: scroll region = rows 1..rows-1."""
        self.rows = shutil.get_terminal_size((80, 24)).lines
        self.tty.write('\n')                       # make room for the bar line
        self.tty.write(f'\x1b[1;{self.rows-1}r')   # DECSTBM scroll region
        self.tty.write(f'\x1b[{self.rows-1};1H')   # put the cursor inside it
        self.tty.flush()

    def _release(self):
        self.tty.write('\x1b[r')                   # restore full-screen scrolling
        self.tty.write(f'\x1b[{self.rows};1H\x1b[2K')  # clear the bar line
        self.tty.flush()

    # -- rendering -----------------------------------------------------------
    @staticmethod
    def _fmt(sec):
        sec = int(max(0, sec)); h, m, s = sec//3600, (sec % 3600)//60, sec % 60
        return f'{h:d}:{m:02d}:{s:02d}' if h else f'{m:02d}:{s:02d}'

    def _text(self, cols):
        frac = (self.n/self.total) if self.total else 0.0
        el = time.time()-self.t0
        eta = (el/self.n)*(self.total-self.n) if self.n else 0.0
        head = f'{self.desc} '
        tail = f' {self.n}/{self.total} [{self._fmt(el)}<{self._fmt(eta)}]'
        if self.postfix:
            tail += f' {self.postfix}'
        width = max(10, cols - len(head) - len(tail) - 3)
        filled = int(round(frac*width))
        bar = '█'*filled + '─'*(width-filled)
        return (head + '|' + bar + '|' + tail)[:cols]

    def draw(self):
        if not self.enabled:
            return
        size = shutil.get_terminal_size((80, 24))
        if size.lines != self.rows:                # terminal was resized
            self._reserve()
        self.tty.write('\x1b7')                    # save cursor
        self.tty.write(f'\x1b[{self.rows};1H\x1b[2K')
        self.tty.write(self._text(size.columns))
        self.tty.write('\x1b8')                    # restore cursor
        self.tty.flush()

    # -- public API ----------------------------------------------------------
    def update(self, n=1, **postfix):
        if postfix:
            self.postfix = ' '.join(f'{k}={v}' for k, v in postfix.items())
        self.n = min(self.n + n, self.total)
        if self.enabled:
            self.draw()
        elif n:                                    # non-tty: one line per item
            print(f'[PROGRESS] {self.n}/{self.total} '
                  f'({100*self.n/max(self.total,1):.0f}%) {self.postfix}', flush=True)

    def set(self, n):
        self.n = max(0, min(int(n), self.total))
        if self.enabled:
            self.draw()

    def close(self):
        if self.enabled:
            self._release()
            try: self.tty.close()
            except Exception: pass


def progress_start(total, desc='reducing'):
    """Create the pinned bar. No-op unless -pb/--progress was given."""
    global _PBAR
    if not getattr(args, 'progress', False) or total <= 0:
        return None
    _PBAR = _PinnedBar(total, desc)
    return _PBAR

def progress_update(n=1, **postfix):
    if _PBAR is not None:
        _PBAR.update(n, **postfix)

def progress_set(n):
    """Set the bar to an absolute position (robust when a loop body may 'continue')."""
    if _PBAR is not None:
        _PBAR.set(n)

def progress_close():
    global _PBAR
    if _PBAR is not None:
        _PBAR.close(); _PBAR = None

# from mpi4py import MPI

parser = argparse.ArgumentParser()

parser.add_argument('--calib_folders','-cf',default=[], nargs='*',
                    help='List of folders to calibrate (accepts wildcards or '
                    'space-delimited list)')

parser.add_argument('--stands_cat','-stc',default=None,
                    help='Standards catalogue ')

parser.add_argument('--survey','-s',default='SDSS',
                    help='Catalogue origins (i.e. SDSS, PS1, or legacy for Legacy Survey/DECaLS)')

parser.add_argument('--folder','-f',default='', nargs='*',
                    help='List of folders to reduce (accepts wildcards or '
                    'space-delimited list)')

parser.add_argument('--email','-e',action='store_true',default=False,
                    help="Send email with photometry results, requires correct credentials in credentials.py (email username and password)")

parser.add_argument('--multipro','-mp',default='inactive',
                    help="Speed up subtraction by utilising multiprocessing and using 5 nodes")

parser.add_argument('--bands','-fb',default=['All'],nargs='+',
                    help="Filters to process, default is All e.g. g,r,i,z,u for SDSS-G, SDSS-R, SDSS-I, SDSS-Z & SDSS-U")

parser.add_argument('--sdsscat','-sdsscat',default=False,action='store_true',
                    help="Use SDSS for reference catalog, default is False. Use if you have the SDSS image but not catalogue")

parser.add_argument('--sci_names','-sn',default=['All'],nargs='+',
                    help="Science Names - names of science object to process, default is all found")

parser.add_argument('--new_only','-new',action='store_true',default=False,
                    help="Processing new images only")

parser.add_argument('--make_log','-log',default='0',
                    help="Outputting log, default is off")

parser.add_argument('--output','-o',default='by_name',#'by_obs_date',
                    help="Dest. for output photometry files, default creates a folder for photometry by observation date e.g. for h_e_20220206_***_1_1.fits output photometry dest will be photometry_data/20220206")

parser.add_argument('--wcs_match_min','-wcsmin',type=float,default=0.15,
                    help="Reject frames from a stack when fewer than this fraction of their detected sources match the reference catalogue (default 0.15; 0 disables). Catches frames whose WCS was never solved, which otherwise dilute the stack.")

parser.add_argument('--trans_max','-transmax',type=float,default=0.75,
                    help="Reject frames from a stack whose transparency is more than this many magnitudes below the best frame in the group (default 0.75; 0 disables). Catches cloud-affected frames, which the seeing check cannot.")

parser.add_argument('--file_exclude','-fexcl',default=None,nargs='+',
                    help="Skip files matching these glob pattern(s) in -f folder mode, e.g. -fexcl '*_astro_*' to drop duplicate products of the same epoch. Applied after -fpat.")

parser.add_argument('--file_pattern','-fpat',default=None,nargs='+',
                    help="Only consider files matching these glob pattern(s) in -f folder mode, e.g. -fpat '*_g_g.fits' '*_r_r.fits'. Useful for SEDM, where a folder holds several products per exposure.")

parser.add_argument('--stack_headers','-stkh',action='store_true',default=False,
                    help="Group images for stacking from FITS HEADERS (object+filter+time) instead of filename patterns. Telescope-agnostic; use for any facility, and for datasets with missing/irregular exposure numbering.")

parser.add_argument('--stack_gap','-stkgap',type=float,default=10.0,
                    help="With -stkh: start a new stack when consecutive exposures are separated by more than this many minutes (default 10, which reproduces per-sequence grouping). Use a large value (e.g. 1440) to stack a whole night together.")

parser.add_argument('--progress','-pb',action='store_true',default=False,
                    help="Show a progress bar pinned to the bottom line; all other output scrolls above it. Requires tqdm.")

parser.add_argument('--ref_img','-refimg',default='auto',
                    help="Name of reference image, default will be to automatically downloading from PS1/SDSS. "
                         "May contain a {filt} placeholder to select a per-filter reference, e.g. "
                         "'ref_imgs/DEEP_PS1{filt}_ZTF26aakjzdt_wcsfix.fits'; if the substituted file is "
                         "missing the pipeline falls back to the automatic PS1/SDSS download for that filter.")

parser.add_argument('--ref_cat','-refcat',default='auto',
                    help="Name of reference catalog, default will be to automatically downloading from PS1/SDSS")

parser.add_argument('--down_date_ql','-qdl',default=[], nargs='+',
                    help='Download data from quicklook, for closest night of observations parse "current_obs" else parse 20220202 for observations on the night of 2022/02/02 ')#(e.g. for 2022/02/02 sunset at La Palma-11:59pm --> observations on the night of 2022/02/02, 12am-sunset at La Palma the next day --> observations on the night of 2022/02/02). Specify another date - note Quicklook only shows the last 7 days with quick reductions')
 
parser.add_argument('--down_date_re','-rdl',default=[], nargs='+',
                    help='Download data from recent data, for closest night of observations parse "current_obs" else parse 20220202 for observations on the night of 2022/02/02 ')# (e.g. for 2022/02/02 sunset at La Palma-11:59pm --> observations on the night of 2022/02/02, 12am-sunset at La Palma the next day --> observations on the night of 2022/02/02). Specify another date - note Recent Data only shows the last 30 days with full reductions')

parser.add_argument('--ims','-i',default='', nargs='+',
                    help='Files to reduce, if stacking then name all in stack but with the first image containing the full path')

parser.add_argument('--cutout','-cut',action='store_true',
                    help='Produce cut outs and deposite in cut_outs, default is False ')

parser.add_argument('--ign_999','-ign999',action='store_true',
                    help='Ignore seeing values ==999 or 966 ')

parser.add_argument('--unsubtract','-un',action='store_true',default=False,
                    help='Unsubtratced option, default is False ')

parser.add_argument('--termoutp','-t',default='normal',
                    help="Type of output to terminal, full/normal/quiet, default is normal")

parser.add_argument('--stack','-stk',action='store_true',default=False,
                    help="Stack where possible, default is False")

parser.add_argument('--stack_all','-stka',action='store_true',default=False,
                    help="Stack ALL files of a band regardless of name matching, default is False")

parser.add_argument('--progress_bar','-progb',action='store_true',default=False,
                    help="Show fixed progress bar at bottom of terminal, default is False")

parser.add_argument('--upfritz','-up',action='store_true',default=False,
                    help="Upload photometry to Fritz (SkyPortal), default is False. If True, will require correct credentials in credntials.py (Fritz token enabled for uploading)")

parser.add_argument('--upfritz_f','-upf',action='store_true',default=False,
                    help="Replace photometry on Fritz (SkyPortal), default is False. If True, will require correct credentials in credntials.py (Fritz token enabled for uploading)")

parser.add_argument('--cleandirs','-cln',action='store_false',default=True,
                    help="Clean directories by removing all intermediate products (cut outs are not included in this) ")

parser.add_argument('--mroundup','-mrup',action='store_true',default=False,
                    help="Morning round up, multiprocessing subtraction on most recent night of observation looking in data/Quicklook/DATE, default is False ")

parser.add_argument('--zp_only','-zp_only',action='store_true',default=False,
                    help="Calculate the zeropoint only and save to data/zeropoints/by_obs_date")

parser.add_argument('--user_zp_sci','-zp_sci',default=None,
                    help="User supplied zeropoint, default is None")

parser.add_argument('--user_zp_ref','-zp_ref',default=None,
                    help="User supplied zeropoint, default is None")

parser.add_argument('--force_image_size','-fis',default=None,type=int,
                    help="Force SWarp output IMAGE_SIZE to this many pixels per side, "
                    "and pad sci+ref into that frame instead of shrinking to the "
                    "smaller resamp shape. Use when the transient is offset from the "
                    "science image centre and the default shrink-and-retry loop "
                    "crops it out. Typical SEDM values: 1000 or larger.")

parser.add_argument('--web_page','-web',action='store_true',default=False,
                    help="Create webpage, default is off")

parser.add_argument('--show_plots','-showp',action='store_true',default=False,
                    help="Make plots pop up")

parser.add_argument('--cal_mode','-cal_mode',default=None,
                    help="Keyword for calibrating science/reference fields")

parser.add_argument('--plc','-lc',default=[None], nargs='*',
                    help="Plots LCs for the given objects, default is off. First arg is folder name, second arg is object name")

parser.add_argument('--plc_space','-lcs',default='mag',
                    help="Space to plot LC's in, default is mag")

parser.add_argument('-relative_flux','-rel',action='store_true',default=False
                    ,help="Plot relative fluxes, default is False")

parser.add_argument('--special_case','-sc',default=None,
                    help="Special case, default is None (e.g. for SN2023ixf where it is bright in a nearby galaxy, we want to use the best bright stars\
                        and they are held in brightstars.cat, so we can use --special_cases brightstars.cat)")

parser.add_argument('--use_swarp','-swarp',action='store_true',default=True,
                    help="Use swarp to stack images and align science image with reference image, default is False")

parser.add_argument('--position','-pos',default='header',nargs='+',
                    help="Position of object in image, default is header, if not in header then parse RA DEC, input as 'HH:MM:SS ±HH:MM:SS")

parser.add_argument('--use_psfex','-psfex',action='store_false',default=True,
                    help="Use PSFEx for PSF modelling (default True). Pass -psfex to disable and use the Python fallback instead.")

parser.add_argument('--psfex_prefer','-psfex_prefer',action='store_true',default=False,
                    help="Legacy behaviour: run SExtractor+PSFEx and prefer the PSFEx stamp when it passes "
                         "quality checks (shape agreement with the subphot_psf measurement, clean stamp, sane chi^2). "
                         "By default the subphot_psf measured kernel is always used and PSFEx is not run.")

parser.add_argument('--telescope_facility','-tel',default='auto',
                    help="Telescope facility. Default 'auto' detects it from the TELESCOP keyword of the first input FITS; pass LT, SEDM, HCT, SLT etc. to override.")

parser.add_argument('--redo_astrometry','-reastrom',default=False,action='store_true',
                    help="Redo astrometry, default is True, set to False if you want to use the astrometry from the header")

parser.add_argument('--list_fits','-ls',default=None,nargs='+',
                    help="List fits files in directory, and present the name, mjd, seeing, airmass, exptime, filter, RA & DEC")

parser.add_argument('--forced_phot','-fp',default=False,nargs='+',
                    help="Perform forced photometry on the images, default is False. If True, forced photometry performed at CAT-RA and CAT-DEC positions")

parser.add_argument('--use_sdss','-sdss',default=False,action='store_true',
                    help="Use SDSS for reference image, default is False")

parser.add_argument('--redo_batch_astrometry','-rebatch',default=None,nargs='+',
                    help="Redo astrometry for batch, default is False")

parser.add_argument('--pros_job_id','-pid',default=None,
                    help="Prospero job ID, default is None")
args = parser.parse_args()

# [SRV] Normalise user-supplied paths so absolute arguments (what a cron
# wrapper naturally passes) work everywhere the code does 'data1_path + arg'.
for _pa in ('ims', 'folder', 'list_fits'):
    _v = getattr(args, _pa, None)
    if isinstance(_v, (list, tuple)):
        setattr(args, _pa, [rel_to_data1(_x) for _x in _v])
    elif isinstance(_v, str) and _v:
        setattr(args, _pa, rel_to_data1(_v))
# -o must stay a bare directory name: 'path + out_dir' is concatenated
# downstream, so an absolute value would escape the pipeline root.
if isinstance(getattr(args, 'output', None), str) and os.path.isabs(args.output):
    _o = args.output.rstrip('/').split('/')[-1]
    print(f'[WARNING] :: [SRV] -o must be a directory NAME, not a path; '
          f'using {_o!r} under {data1_path}')
    args.output = _o

# Enable stacking if stack_all is requested
if hasattr(args, 'stack_all') and args.stack_all:
    args.stack = True

# print(args.ims)
# # print(args.forced_phot)

# print(args.redo_batch_astrometry)
# print(args.folder)
if args.redo_batch_astrometry!=None:
    redo_list = []

    for n in range(len(args.redo_batch_astrometry)):
        
        [redo_list.append(i) for i in glob.glob(data1_path+args.redo_batch_astrometry[n]+'/*') if i.endswith('.fits')]
    print(redo_list)

    for i in range(len(redo_list)):
        print(info_g+' Redoing astrometry for',redo_list[i],' &', redo_list[i+1])

        subphot_data().batch_update_astrometry([redo_list[i],redo_list[i+1]])#,ra=self.sci_img_hdu.header[self.RA_kw],dec=self.sci_img_hdu.header[self.DEC_kw]) 

# sys.exit(1)
# print(args.sci_names)
if args.list_fits!=None:
    if args.list_fits=='':
        print('Please specify a folder to list fits files in')
    else:
        if args.list_fits[0]=='.':
            # print(info_g+' Listing fits files in current data list_fits listed in subphot_credentials.py')
            args.list_fits[0]=''
        else:
            # print(info_g+' Listing fits files in '+args.list_fits)
            args.list_fits[0]+='/'
        print(info_g+' Listing fits files in '+data1_path+' '.join(args.list_fits))
        # print(args.list_fits[1].split('/'))
        for fold_ in args.list_fits:
            arr=[]
            for file_ in os.listdir(data1_path+fold_):
                # print(file_)
                if file_.endswith('.fits') and 'tpv' not in file_ and not file_.startswith('.'):
                    # print(file_)
                    hdul = fits.open(data1_path+fold_+file_)
                    # for key,val in hdul[0].header.items():
                    #     print(key,val)
                    hdr = hdul[0].header
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER'],hdr['MJD-OBS'],hdr['SEEING'],hdr['AIRMASS'],hdr['EXPTIME'],hdr['DATE-OBS']])
                        continue#,hdr['RA'],hdr['DEC']])
                    except:pass
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER1'],hdr['MJD'],hdr['L1SEESEC'],hdr['AIRMASS'],hdr['EXPTIME'],hdr['DATE-OBS']])
                        continue#hdr['CAT-RA'],hdr['CAT-DEC']])
                    except:pass
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER'],hdr['MJD_OBS'],hdr['FWHM'],hdr['AIRMASS'],hdr['EXPTIME'],hdr['UTC']])
                        continue#,hdr['OBJRA'],hdr['OBJDEC']])
                    except:pass
                    try:
                        arr.append([file_,hdr['OBJECT'],hdr['FILTER'],hdr['MJD_OBS'],hdr['FWHM'],hdr['AIRMASS'],hdr['EXPTIME'],hdr['DATE-OBS']])
                        continue#,hdr['OBJRA'],hdr['OBJDEC']])
                    except Exception as e:
                        pass

            df = pd.DataFrame(arr,columns=['IMG','OBJ','FILT','MJD','SEE','AIRM','EXPTIME','DATE-OBS'])

            print(tabulate(df.sort_values(by=['OBJ','FILT']), headers='keys', tablefmt='psql'))
            df = df.sort_values(by=['FILT','MJD'])
            df.to_csv(path+args.list_fits[0].split('/')[1]+'_fits_list.csv',index=False)
            print(info_g+' Saved to '+path+args.list_fits[0].split('/')[1]+'_fits_list.csv')

    sys.exit(1)

# print(args.use_psfex)
# print(args.use_swarp)
args.telescope_facility = args.telescope_facility.upper()

def _detect_telescope_facility():
    """Detect the telescope facility from the TELESCOP keyword of the first
    readable input FITS (from -i or -f). Returns None if nothing is readable."""
    import glob as _glob
    candidates = []
    if len(args.ims) > 0:
        candidates += list(args.ims)
    if args.folder not in ('', None, []):
        for _fold in args.folder:
            for _base in (_fold, path+_fold, data1_path+_fold):
                if os.path.isdir(_base):
                    candidates += sorted(_glob.glob(os.path.join(_base, '*.fits*')))[:5]
                    break
    for _cand in candidates:
        for _p in (_cand, _cand+'.fits', path+_cand, path+_cand+'.fits'):
            if not os.path.exists(_p):
                continue
            # TELESCOP can live in a later HDU (e.g. LCOGT .fz compressed files)
            _tel,_origin = None,''
            try:
                _hdul = fits.open(_p)
                for _h in _hdul:
                    if 'TELESCOP' in _h.header:
                        _tel = _h.header['TELESCOP']
                        _origin = str(_h.header.get('ORIGIN','')).strip()
                        break
                _hdul.close()
            except Exception:
                break
            if _tel is None:
                break
            _tel = str(_tel).strip()
            if _tel in SEDM: return 'SEDM'
            if _tel == 'Liverpool Telescope': return 'LT'
            if _origin=='LCOGT' or _tel.lower().startswith(('0m4','1m0','2m0')): return 'LCOGT'
            if _tel in ('DCT','LDT'): return 'LDT'
            if 'NTT' in _tel.upper() or 'EFOSC' in _tel.upper(): return 'NTT'
            if _tel in ('HCT','SLT','TJO','GTC','NOT','LOT'): return _tel
            return _tel.upper()
    return None

if args.telescope_facility in ('AUTO',''):
    _detected = _detect_telescope_facility()
    if _detected is not None:
        args.telescope_facility = _detected
        print(info_g+f' Telescope facility auto-detected from FITS header: {_detected}')
    else:
        args.telescope_facility = 'SEDM'
        print(warn_y+' Could not auto-detect telescope facility (no readable FITS given); defaulting to SEDM')

# sys.exit(1)
seeing_limit =5

# sys.exit(1)
    #setting data to todays data in format YYYYMMDD
t = date.today()
TIME = datetime.datetime.now().strftime("%H:%M:%S")
year,month,dayy = t.strftime("%Y"),t.strftime("%m"),t.strftime("%d")
today = Time(f'{year}-{month}-{dayy} {TIME}')
TODAY = t.strftime("%Y%m%d") #todays date in YYYYMMDD format
apo = Observer.at_site("palomar" if args.telescope_facility in SEDM else "lapalma")
sun_set_today = apo.sun_set_time(today, which="nearest") #sun set on day of observing
time_suns_today = "{0.iso}".format(sun_set_today)[-12:]
sun_set_tomorrow = apo.sun_set_time(today,which="next")
time_suns_tomorrow = "{0.iso}".format(sun_set_tomorrow)[-12:]


if time_suns_today<TIME<'23:59:59':
    date_ = TODAY
    DATE = re.sub("-","",date_)
if '00:00:00'<TIME<time_suns_tomorrow:
    date_ = str(t - datetime.timedelta(days=1))
    DATE=re.sub("-","",date_)


if args.mroundup!=False:
    args.folder=[f"Quicklook/"+DATE]
    args.make_log='night_log'
    

try:
    if store_lc_ims:
        safe_makedirs(data1_path+'entire_lc_imgs', 'entire_lc_imgs')
except Exception as e:
    print(warn_r+f' Could not create entire_lc_imgs directory at {data1_path}entire_lc_imgs: {e}')


def archive_lc_images(files, name, base_dir=None):
    """[SRV] Copy the raw science frames of an epoch into
    entire_lc_imgs/<object>/ so the full light-curve image set accumulates.

    Previously this lived inline in the folder (-f) branch only, so the
    single-image (-i) mode the server uses for each new source archived
    NOTHING.  Both paths now call this.

    files    : list of paths/basenames for the epoch (a stack group or one file)
    name     : object name (directory under entire_lc_imgs/)
    base_dir : folder to resolve bare basenames against (stack members carry
               no path); defaults to the directory of the first entry.
    """
    if not store_lc_ims or not name:
        return 0
    _clean = re.sub(r'[^\w.+-]', '_', str(name).strip())
    if not _clean:
        return 0
    lc_path = os.path.join(data1_path + 'entire_lc_imgs', _clean) + os.sep
    if safe_makedirs(lc_path, 'entire_lc_imgs/<object>') is None:
        return 0
    if base_dir is None:
        base_dir = os.path.dirname(str(files[0])) if files else ''
    n_done = 0
    for fit in files:
        fit = str(fit)
        if not fit.endswith('.fits'):
            fit += '.fits'
        # resolve in turn against: as-given, the group's folder, data1_path
        if not os.path.exists(fit):
            for _cand in ([os.path.join(base_dir, os.path.basename(fit))] if base_dir else []) + \
                         [data1_path + fit.lstrip('/'),
                          os.path.join(data1_path, base_dir.lstrip('/'), os.path.basename(fit))
                          if base_dir else None]:
                if _cand and os.path.exists(_cand):
                    fit = _cand
                    break
        if not os.path.exists(fit):
            print(warn_y + f' Could not archive {fit} — file not found')
            continue
        dst = os.path.join(lc_path, os.path.basename(fit))
        if os.path.exists(dst):
            n_done += 1
            continue
        try:
            shutil.copy(fit, dst)
            n_done += 1
        except Exception as e:
            print(warn_r + f' Could not archive {fit} to {lc_path}: {e}')
    if n_done:
        print(info_g + f' Archived {n_done} image(s) to {lc_path}')
    return n_done


if args.termoutp!='quiet':
    mk_new = True
# config files
if not os.path.exists(path+'config_files'):
    os.makedirs(path+'config_files')

if not os.path.exists(path+'config_files/config.swarp') or mk_new==True:
    writeswarpdefaultconfigfile(args.termoutp)

if not os.path.exists(path+'config_files/config_comb.swarp')  or mk_new==True:
    writeswarpconfigfile(args.termoutp)

if not os.path.exists(path+'config_files/swarp_sdss.conf')  or mk_new==True:
    writesdssswarpconfigfile(args.termoutp)

if not os.path.exists(path+'config_files/default.param')  or mk_new==True:
    writepsfexparfile()

if not os.path.exists(path+'config_files/prepsfex.sex')  or mk_new==True:
    prepsexfile()

if not os.path.exists(path+'config_files/psfex_conf.psfex')  or mk_new==True:
    psfexfile(args.termoutp)

FILTERS = {'g':'SDSS-G','r':'SDSS-R','i':'SDSS-I','z':'SDSS-Z','u':'SDSS-U','B':'Bessell-B','V':'Bessell-V','R':'Bessell-R','I':'Bessell-I'}

def canonical_filter(header, facility=None):
    """Canonical filter name (e.g. 'SDSS-G') from any telescope's header.

    Reading FILTER1/FILTER only works for a subset of telescopes: NOT stores it
    in SEQID ('g_3x90s'), LOT/SLT use Astrodon names ('gp_Astrodon_2019'), TJO
    writes 'SDSS g'.  Consult the telescope's header_kw entry first, then fall
    back to the usual keywords, and finally reduce whatever we find to its
    ugriz/BVRI letter.
    """
    fac = {'LT':'Liverpool Telescope','SEDM':'SEDM-P60','NTT':'ESO-NTT'}.get(facility, facility)
    keys = []
    kw = header_kw.get(fac, {}).get('filter') if fac else None
    if kw and kw not in ('-', 'user'):
        keys.append(kw)
    keys += ['FILTER1', 'FILTER', 'SEQID', 'FILTERS', 'FILTER2']

    raw = None
    for k in keys:
        if k in header and str(header[k]).strip() not in ('', '-'):
            raw = str(header[k]).strip(); break
    if raw is None:
        return None
    if raw in FILTERS.values():          # already canonical
        return raw

    low = raw.lower()
    if 'astrodon' in low:                # 'gp_Astrodon_2019' -> g
        low = low.split('p_')[0]
    for pre in ('sdss-', 'sdss_', 'sdss ', 'ps1-', 'ps1_'):
        if low.startswith(pre):
            low = low[len(pre):]
    tok = low.replace('-', '_').replace(' ', '_').split('_')[0]
    for cand in (tok, tok[:1]):          # 'rp' -> 'r', 'g' -> 'g'
        if cand in FILTERS:
            return FILTERS[cand]
    for ch in low:                       # last resort: first ugriz letter
        if ch in FILTERS:
            return FILTERS[ch]
    return None
if args.bands==['All']:
    args.bands_to_process=['g','r','i','z','u']#,'B','V','R','I']
else:
    args.bands_to_process=args.bands

# print(args.make_log==0,args.make_log=='0')
# sys.exit(1)

if args.make_log!=False:

    log_dest=args.make_log
    # print(args.output)
    if args.output in ['by_name']:
        # print('yes')
        if args.folder!='':
            for k in range(len(args.folder)):
                # print(args.folder[k])
                if '/' in args.folder[k]:out_dir=args.folder[k].split('/')[-1]
                else:out_dir=args.folder[k]
                safe_makedirs(data1_path+out_dir, 'output directory')
        else:
            out_dir='photometry'
    else:out_dir=args.output.split('/')[-1]


    # 'by_obs_date' is a keyword, not a directory — its outputs/logs live in photometry_date/
    if args.output not in ('by_name','by_obs_date'):
        safe_makedirs(data1_path+out_dir, 'output directory')

    # [SRV] '0' is the DEFAULT value of --make_log (meaning "log beside the
    # output"), not a directory name — mkdir'ing it created a stray '0/'
    # folder in the pipeline root on every single run.
    if log_dest not in ('0', 0, '', None):
        safe_makedirs(data1_path+str(log_dest), 'log directory')
    if log_dest=='night_log':log_name = f"{data1_path}night_log/{DATE}_night_log.log"
    else:log_name = f"{data1_path}{log_dest}/{out_dir}_log.log"
    if log_dest=='0':
        log_dest=out_dir
        if out_dir=='by_obs_date':
            safe_makedirs(data1_path+'photometry_date', 'photometry_date')
            log_name = f"{data1_path}photometry_date/{DATE}_log.log"
        else:
            log_name = f"{data1_path}{out_dir}/{out_dir}_log.log"

    if args.pros_job_id!=None:
        # nightly-routine logs live outside the pipeline root by design
        _nr_dir = f'/mnt/aridata1/users/arikhind/phot_data/nightly_routine_logs/{DATE}'
        try:
            os.makedirs(_nr_dir, exist_ok=True)
            log_name = f'{_nr_dir}/SPNR_{args.pros_job_id}_py.log'
        except Exception as _nr_e:
            print(f'[WARNING] :: [SRV] nightly-routine log dir unavailable ({_nr_e}); '
                  f'falling back to {log_name}')


    if args.mroundup!=False:
        log_name = f"{data1_path}{log_dest}/{DATE}_night_log.log"

    # [SRV] The log file used to receive only the handful of sp_logger.*
    # calls — everything else in the pipeline uses print(), so a run that
    # emitted ~4800 console lines wrote 43 to disk, and those carried raw
    # ANSI colour codes.  Tee stdout/stderr into the log instead, stripping
    # escapes, so the file is a faithful, greppable record of the run.
    safe_makedirs(os.path.dirname(log_name), 'log directory')

    class _Tee:
        _ANSI = re.compile(r'\x1b\[[0-9;]*[A-Za-z]')

        def __init__(self, stream, fh):
            self._stream, self._fh = stream, fh

        def write(self, data):
            try:
                self._stream.write(data)
            except Exception:
                pass
            try:
                self._fh.write(self._ANSI.sub('', data))
                self._fh.flush()
            except Exception:
                pass
            return len(data)

        def flush(self):
            for _o in (self._stream, self._fh):
                try:
                    _o.flush()
                except Exception:
                    pass

        def isatty(self):
            return getattr(self._stream, 'isatty', lambda: False)()

        def fileno(self):
            return self._stream.fileno()

    try:
        _log_fh = open(log_name, 'a', encoding='utf-8', buffering=1)
        _log_fh.write(f'\n{"="*70}\n'
                      f'RUN {datetime.datetime.now().isoformat(timespec="seconds")}  '
                      f'cwd={os.getcwd()}\n'
                      f'argv: {" ".join(sys.argv)}\n{"="*70}\n')
        sys.stdout = _Tee(sys.stdout, _log_fh)
        sys.stderr = _Tee(sys.stderr, _log_fh)
    except Exception as _lg_e:
        print(f'[WARNING] :: [SRV] could not open log file {log_name}: {_lg_e}')

    sp_logger = logging
    # force=True: basicConfig is a no-op if any import already configured the
    # root logger, which silently produced empty logs.
    try:
        sp_logger.basicConfig(level=logging.INFO, encoding='utf-8',
                              handlers=[logging.StreamHandler(sys.stdout)],
                              format='%(message)s', force=True)
    except TypeError:      # python < 3.8 has no force=
        sp_logger.basicConfig(level=logging.INFO, encoding='utf-8',
                              handlers=[logging.StreamHandler(sys.stdout)],
                              format='%(message)s')

    print(info_g+f' Logging to {log_name}')
    args.sp_logger=sp_logger
else:args.sp_logger=None


# sys.exit()
class multi_subtract():
    def __init__(self,folder,folder_path=None):
        #folder_path is the path from data_1path to where the folders are (so that when cron is called, the same folder path can be used from data1_path instead of having to
        # specify the full path from data1_path each time, it is easier for large jobs)
        if folder_path!=None:
            self.folder_path=folder_path
        else:
            folder_path=''
        self.path=path
        self.data1_path=data1_path
        self.g_fits,self.r_fits,self.i_fits,self.z_fits,self.u_fits=[],[],[],[],[]
        self.bess_V_fits,self.bess_R_fits,self.bess_I_fits,self.bess_B_fits = [],[],[],[]
        self.fits_dict = {'SDSS-G':self.g_fits,'SDSS-R':self.r_fits,'SDSS-I':self.i_fits,'SDSS-Z':self.z_fits,'SDSS-U':self.u_fits,
                        'Bessell-V':self.bess_V_fits,'Bessell-R':self.bess_R_fits,'Bessell-I':self.bess_I_fits,'Bessell-B':self.bess_B_fits}
        self.filts = {'SDSS-U':'u','SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','SDSS-Z':'z',
                    'Bessell-V':'V','Bessell-R':'R','Bessell-I':'I','Bessell-B':'B',}
                    # 'up_Astrondon_2018':'u','gp_Astrondon_2018':'g','rp_Astrondon_2018':'r','ip_Astrondon_2018':'i','zp_Astrondon_2018':'z',}


        self.folder = folder #specific folder
        self.FOLDER = str(self.folder_path)+'/'+str(self.folder) #full folder path spliced together based off first argument
        self.folder = re.sub(str(self.folder_path),'',self.folder)
        self.folder = re.sub('//','',self.folder)
        self.FOLDER = re.sub('//','/',self.FOLDER)

    def trunc(self,string_,instrument=args.telescope_facility):
        #truncating fits file to identify fits to be stacked and whether the fits is a quicklook reduction (_1_9.fits) or recent data reduction (_1_1.fits)
        if instrument =='LT' or instrument=='lt':
            # print(string_,'_'.join(string_.split('_')[0:4])+'_','_'+'_'.join(string_.split('_')[5:]))
            return ['_'.join(string_.split('_')[0:4])+'_','_'+'_'.join(string_.split('_')[5:])]
        elif instrument=='HCT' or instrument=='hct':
            return [string_,'']
        elif instrument=='SEDM' or instrument=='sedm':
            return [string_,'']
        elif instrument=='SLT' or instrument=='slt':
            return [string_,'']

    def expand_compressed(self):
        # folder mode only understands plain *.fits: write a single-HDU .fits next
        # to every .fits.gz / .fits.fz / .fz (originals kept, existing .fits reused)
        _dir = f"{self.data1_path}{self.FOLDER}"
        for _f in sorted(os.listdir(_dir)):
            if _f.startswith('.') or not _f.endswith(('.fits.gz','.fits.fz','.fz')):
                continue
            _base = _f[:-len('.fz')] if _f.endswith('.fz') and not _f.endswith('.fits.fz') else _f.rsplit('.',1)[0]
            if not _base.endswith('.fits'): _base += '.fits'
            _out = os.path.join(_dir, _base)
            if os.path.exists(_out):
                continue
            try:
                with fits.open(os.path.join(_dir, _f)) as _hdul:
                    _hdu = next((h for h in _hdul if getattr(h,'is_image',False) and h.data is not None), None)
                    if _hdu is None:
                        print(warn_y+f' No image data in {_f} — skipping')
                        continue
                    fits.PrimaryHDU(data=_hdu.data, header=_hdu.header.copy()).writeto(_out)
                print(info_g+f' Expanded compressed input {_f} -> {_base}')
            except Exception as _e:
                print(warn_y+f' Could not expand {_f}: {_e}')

    def analyse_folder(self):
        self.fits_files,self.all_fits,self.all_in_dir = [],[],[]
        self.expand_compressed()

        [self.all_in_dir.append(f) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if f.endswith('.fits')]
        # print(self.all_in_dir)

        # Stack all mode: collect all FITS files without name matching
        if hasattr(args, 'stack_all') and args.stack_all:
            print(info_g+' Stack ALL mode: stacking all files of matching filter band')
            [self.all_fits.append([f, '']) for f in os.listdir(f"{self.data1_path}{self.FOLDER}")
             if f.endswith('.fits') and not f.startswith('.') and 'tpv' not in f]
        # Standard mode: use file truncation to match related files
        elif args.telescope_facility=='LT' or args.telescope_facility=='lt':
            [self.all_fits.append(self.trunc(f)) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and f.startswith('h_') and self.trunc(f) not in self.all_fits)
            and not f.startswith('.')]
        elif args.telescope_facility=='HCT' or args.telescope_facility=='hct':
            [self.all_fits.append(self.trunc(f,instrument='HCT')) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and self.trunc(f,instrument='HCT') not in self.all_fits)
            and not f.startswith('.')]
        elif args.telescope_facility=='SEDM' or args.telescope_facility=='sedm':
            [self.all_fits.append(self.trunc(f,instrument='SEDM')) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and self.trunc(f,instrument='SEDM') not in self.all_fits)
            and not f.startswith('.') and 'tpv' not in f]
        elif args.telescope_facility=='SLT' or args.telescope_facility=='slt':
            [self.all_fits.append(self.trunc(f,instrument='SLT')) for f in os.listdir(f"{self.data1_path}{self.FOLDER}") if (f.endswith('.fits') and self.trunc(f,instrument='SLT') not in self.all_fits)
            and not f.startswith('.')]

        # print(self.all_fits)
        # # print()
        # print(self.all_in_dir)
        # sys.exit(1)

        # Header-based grouping: build stack groups from OBJECT/FILTER/time in the
        # FITS headers rather than from filename patterns.  Telescope-agnostic and
        # immune to missing or irregular exposure numbering.
        if getattr(args, 'stack_headers', False):
            import fnmatch
            _pats = getattr(args, 'file_pattern', None)
            _excl = getattr(args, 'file_exclude', None)
            _recs, _nseen = [], 0
            for _f in sorted(os.listdir(f"{self.data1_path}{self.FOLDER}")):
                if not _f.endswith('.fits') or _f.startswith('.') or 'tpv' in _f:
                    continue
                _nseen += 1
                # restrict to the requested products (e.g. SEDM writes several
                # files per exposure; only the *_g_g.fits style ones are wanted)
                if _pats and not any(fnmatch.fnmatch(_f, _p) for _p in _pats):
                    continue
                if _excl and any(fnmatch.fnmatch(_f, _p) for _p in _excl):
                    continue
                try:
                    _h = fits.open(f"{self.data1_path}{self.FOLDER}/{_f}")[0].header
                except Exception as _e:
                    print(warn_y+f' Skipping unreadable {_f}: {_e}')
                    continue
                _obj = str(_h.get('OBJECT','?')).split('_')[0].split(' ')[0].strip()
                # canonical, so telescopes that hide the filter in SEQID/Astrodon
                # names still group per band instead of lumping every band together
                _filt = canonical_filter(_h, args.telescope_facility) or '?'
                _mjd = None
                for _kw in ('MJD','MJD-OBS','MJD_OBS'):
                    if _kw in _h:
                        try: _mjd = float(_h[_kw]); break
                        except Exception: pass
                if _mjd is None:
                    try: _mjd = float(Time(str(_h['DATE-OBS']).strip()).mjd)
                    except Exception: _mjd = float('inf')
                if _mjd > 2400000:      # JD supplied instead of MJD
                    _mjd -= 2400000.5
                _recs.append((_mjd, _f, _obj, _filt))

            _groups, _gap_d = {}, float(getattr(args,'stack_gap',10.0))/1440.0
            for _mjd,_f,_obj,_filt in sorted(_recs):
                _key = (_obj,_filt)
                _bucket = _groups.setdefault(_key, [])
                if _bucket and (_mjd - _bucket[-1][-1][0]) <= _gap_d:
                    _bucket[-1].append((_mjd,_f))       # same run
                else:
                    _bucket.append([(_mjd,_f)])          # start a new run
            _ngroups = 0
            for (_obj,_filt),_runs in sorted(_groups.items()):
                for _run in _runs:
                    _names = [re.sub(r'\.fits$','',_x[1]) for _x in _run]
                    if args.stack==False and len(_names)>1:
                        for _n in _names:                # stacking off: one entry per image
                            self.fits_files.append([[_n],1,f"{_n}.fits"]); _ngroups += 1
                    else:
                        self.fits_files.append([_names,len(_names),f"{_names[0]}.fits"]); _ngroups += 1
            if _pats or _excl:
                print(info_g+f' File selection (keep={_pats}, drop={_excl}): '
                      f'kept {len(_recs)} of {_nseen} non-tpv files')
            print(info_g+f' Header grouping: {len(_recs)} images -> {_ngroups} '
                  f'{"stack group(s)" if args.stack else "entry(s)"} '
                  f'(gap {getattr(args,"stack_gap",10.0):.0f} min)')

        # Stack all mode: create single entry with all files in directory
        elif hasattr(args, 'stack_all') and args.stack_all:
            all_files = [f for f in os.listdir(f"{self.data1_path}{self.FOLDER}")
                        if f.endswith('.fits') and not f.startswith('.') and 'tpv' not in f]
            if all_files:
                # Remove .fits extension for consistency with rest of code
                file_names = [re.sub(r'\.fits$', '', f) for f in all_files]
                self.fits_files.append([file_names, len(all_files), all_files[0]])
                print(info_g+f' Stack ALL: collected {len(all_files)} FITS files for processing')
        elif len(self.all_fits)==1 and sum([self.all_fits[0][0] in x for x in self.all_in_dir])==0:
            self.fits_files.append([self.all_fits[0],1])
            # print(info_g+f' {self.all_fits[0][0]} is a single image')
            # sys.exit(1)
        else:
            for self.file_arr in self.all_fits:
                self.file,self.end_string = self.file_arr
                self.end_string = re.sub(r'\.fits$', '', self.end_string)
                self.number = sum([self.file in x for x in os.listdir(f"{self.data1_path}{self.FOLDER}") if not x.startswith('.') and x.endswith('.fits')])

                if self.number>1:
                    self.file_name,self.fi = [],[]
                    for k in np.linspace(1,self.number,self.number):
                        k=int(k)
                        if k==1:
                            self.file_name.append(f"{self.file}{k}{self.end_string}")
                            self.fi.append(f"{self.file}{k}{self.end_string}.fits)")
                        if k!=1:
                            self.file_name.append(f"{self.file}{k}{self.end_string}")
                            self.fi.append(f"{self.file}{k}{self.end_string}.fits)")



                    if all(h in os.listdir(f"{self.data1_path}{self.FOLDER}") for h in self.fi)==False:
                        self.fits_files.append([self.file_name,self.number,f"{self.file}1{self.end_string}.fits"])
                else:
                    if self.file.endswith('.fits')==False:
                        self.fits_files.append([[f"{self.file}1{self.end_string}.fits"],1])
                    else:
                        self.fits_files.append([[self.file],1])
                        

        # print(self.fits_files)
        for i in range(len(self.fits_files)):
            if len(self.fits_files[i])==2:
                try:self.fits_hdu = fits.open(self.data1_path+str(self.FOLDER)+"/"+str(self.fits_files[i][0][0]))[0].header
                except:self.fits_hdu = fits.open(self.data1_path+str(self.FOLDER)+"/"+str(self.fits_files[i][0])).header
            else:
                # print(self.fits_files[i][2])
                self.fits_hdu = fits.open(self.data1_path+str(self.FOLDER)+"/"+str(self.fits_files[i][2]))[0].header
            # print(self.fits_hdu['FILTER'])
            self.fits_filt = canonical_filter(self.fits_hdu, args.telescope_facility)
            self.fits_obj  = str(self.fits_hdu.get('OBJECT', '')).strip()
            if self.fits_filt is None or self.fits_obj == '':
                print(warn_y+f' No usable FILTER/OBJECT in {self.fits_files[i][-1]} — skipping')
                continue 

            if str(self.fits_hdu.get('TELESCOP','')).strip() in SEDM: #specifically for P60
                # [OBJ] canonical: strips 'ACQ-' and the trailing filter token
                # (the old split('-')[1] truncated hyphenated names)
                self.fits_obj, self._f = clean_object_name(self.fits_obj, return_filter=True)
                if self._f is None:
                    self._f = self.fits_filt

                # canonical_filter above may already have produced the
                # canonical name — only remap raw single-letter values
                self.fits_filt = {'r':'SDSS-R','g':'SDSS-G','i':'SDSS-I','u':'SDSS-U',
                                'B':'Bessell-B','V':'Bessell-V','R':'Bessell-R','I':'Bessell-I'}.get(self.fits_filt, self.fits_filt)

                                #   }[self.fits_filt]
                # print(self.fits_obj,self.fits_filt)

            if args.bands==['All'] and args.sci_names == ['All']:
                FILT = self.fits_filt
                filter_array = self.fits_dict[FILT]
                filter_array.append(self.fits_files[i])
                self.fits_dict[FILT] = filter_array


            elif args.bands!=['All'] and args.sci_names == ['All']:
                for filt_band in args.bands_to_process:
                    FILT = FILTERS[filt_band]

                    if self.fits_filt== FILT:
                        filter_array = self.fits_dict[FILT]
                        filter_array.append(self.fits_files[i])   
                        self.fits_dict[FILT] = filter_array

            elif args.bands==['All'] and args.sci_names != ['All']:
                FILT = self.fits_filt

                if self.fits_obj in args.sci_names or self.fits_obj in args.sci_names:
                    filter_array = self.fits_dict[FILT]
                    filter_array.append(self.fits_files[i])
                    self.fits_dict[FILT] = filter_array
            

            elif args.bands!=['All'] and args.sci_names != ['All']:
                for filt_band in args.bands_to_process:
                    FILT = FILTERS[filt_band]

                    if self.fits_filt == FILT and self.fits_obj in args.sci_names:
                        filter_array = self.fits_dict[FILT]
                        filter_array.append(self.fits_files[i])
                        self.fits_dict[FILT] = filter_array

        # Sort files within each filter by date (MJD)
        for filt_key in self.fits_dict.keys():
            if len(self.fits_dict[filt_key]) > 0:
                # Extract MJD for each file and sort
                file_mjd_pairs = []
                for file_entry in self.fits_dict[filt_key]:
                    try:
                        # Handle both single file and multi-file entries
                        if isinstance(file_entry[0], list):
                            fits_path = f"{self.data1_path}{self.FOLDER}/{file_entry[0][0]}"
                        else:
                            fits_path = f"{self.data1_path}{self.FOLDER}/{file_entry[0]}"
                        fits_hdr = fits.open(fits_path)[0].header
                        # Try to extract MJD from various header keywords
                        mjd = None
                        for mjd_kw in ['MJD', 'MJD-OBS', 'MJD_OBS']:
                            if mjd_kw in fits_hdr:
                                mjd = float(fits_hdr[mjd_kw])
                                break
                        if mjd is not None:
                            file_mjd_pairs.append((file_entry, mjd))
                        else:
                            # If no MJD found, append with a placeholder (will sort to end)
                            file_mjd_pairs.append((file_entry, float('inf')))
                    except Exception as e:
                        # If error reading file, append with placeholder
                        file_mjd_pairs.append((file_entry, float('inf')))

                # Sort by MJD
                file_mjd_pairs.sort(key=lambda x: x[1])
                # Update fits_dict with sorted files
                self.fits_dict[filt_key] = [entry[0] for entry in file_mjd_pairs]

        # print(self.fits_dict)
        return self.fits_dict 

filt_kws = ['FILTER1','FILTER']
name_kws = ['OBJECT','TARGET','TCSTGT']
date_obs_kws = ['DATE-OBS','DATE','UTC']
def _run_pipeline(sub_obj, sp_logger):
    """Run the full reduction sequence on a subtracted_phot object.

    Executes each pipeline step in order and returns the photometry dict on
    success, or None if any step sets sys_exit=True.  Using early-returns
    instead of nested if/else keeps the control flow linear and easy to follow.
    """
    if args.unsubtract:
        sub_obj.to_subtract = False

    steps = [
        ('Background subtraction',   sub_obj.bkg_subtract),
        ('Cosmic-ray removal',        sub_obj.remove_cosmic),
        ('Image alignment',           sub_obj.swarp_ref_align if args.use_swarp else sub_obj.py_ref_align),
        ('PSF convolution',           sub_obj.psfex_convolve_images if args.use_psfex else sub_obj.py_convolve_images),
        ('Reference catalog',         sub_obj.gen_ref_cat),
        ('Combined PSF',              sub_obj.combine_psf),
        ('Zeropoint calibration',     sub_obj.get_zeropts),
        ('Scaled subtraction',        sub_obj.scaled_subtract),
    ]
    for step_name, step_fn in steps:
        step_fn()
        if sub_obj.sys_exit:
            print(warn_r+f' Pipeline stopped after: {step_name}')
            return None

    result = sub_obj.get_photometry()
    if sub_obj.sys_exit:
        return None

    if args.upfritz or args.upfritz_f:
        sub_obj.upload_phot()

    return result


def run_subtraction(data_dict):
    """Process all FITS files in one filter band, returning list of phot dicts."""
    filts = {'SDSS-U':'u','SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','SDSS-Z':'z'}
    p60_filt_map = {'r':'SDSS-R','g':'SDSS-G','i':'SDSS-I','u':'SDSS-U','z':'SDSS-Z'}
    f_time_start = time.time()
    fits_files, filter_, FOLDER, new_only = (
        data_dict['fits'], data_dict['filter'], data_dict['FOLDER'], data_dict['new_only'])
    final_phot = []

    print(colored('---------------------------------------------------------------------------------------------','yellow'))
    print(info_g+f" For {filter_}, there are {len(fits_files)} fits")
    print(info_g+f" For {filter_}, there are {len(fits_files)} fits")  # Also log to file

    if len(fits_files) == 0:
        print(warn_y+f' No fits files found in filter: {filter_}')
        return final_phot

    progress_start(len(fits_files), desc=f'{FOLDER}')
    for file_idx, file_array in enumerate(fits_files, 1):
            # Show progress bar if requested (output to stderr to avoid stdout interference)
            if args.progress_bar:
                progress = file_idx / len(fits_files)
                bar_length = 40
                filled = int(bar_length * progress)
                bar = '█' * filled + '░' * (bar_length - filled)
                bar_text = f'\r[{bar}] {file_idx}/{len(fits_files)} ({int(progress*100)}%)'
                sys.stderr.write(bar_text)
                sys.stderr.flush()

            if file_array[1]=='1' or file_array[1]==1:
                # the stored name has no extension, but the guard below requires
                # one — without this every single-frame group was silently skipped
                fits_file = file_array[2] if len(file_array) > 2 else file_array[0][0]
                if not str(fits_file).endswith('.fits'):
                    fits_file = f'{fits_file}.fits'
                sub_file = [re.sub('.fits','',file_array[0][0])]

            else:
                fits_file = file_array[2]      
                sub_file = file_array[0]
            
            if not fits_file.endswith('.fits') or 'tpv' in fits_file:
                continue

            sub_file[0] = str(FOLDER) + '/' + sub_file[0]
            if data1_path not in sub_file[0]:
                sub_file[0] = data1_path + sub_file[0]

            if args.termoutp != 'quiet':
                print(colored('════════════════════════════════════════════════════════════════════════════════════════════════','magenta'))
                print(colored(f'  ⭐ PROCESSING FILE {file_idx} OF {len(fits_files)} [{file_idx}/{len(fits_files)}]', 'cyan'))
                label = sub_file[0] if len(sub_file) == 1 else ', '.join(sub_file)
                print(info_g+f' Performing image subtraction on {label}')
                print(colored('════════════════════════════════════════════════════════════════════════════════════════════════','magenta'))
                # Also log to file
                print(colored(f'  ⭐ PROCESSING FILE {file_idx} OF {len(fits_files)} [{file_idx}/{len(fits_files)}]', 'cyan'))
                print(info_g+f' Performing image subtraction on {label}')

            # Extract filter / object name / date from header
            sci_hdr = fits.open(f'{data1_path}{FOLDER}/{fits_file}')[0].header
            filt = canonical_filter(sci_hdr, args.telescope_facility) \
                   or next((sci_hdr[k] for k in filt_kws if k in sci_hdr), None)
            name = next((sci_hdr[k] for k in name_kws if k in sci_hdr), None)
            date_obs = next((sci_hdr[k] for k in date_obs_kws if k in sci_hdr), None)
            if filt is None or name is None or date_obs is None:
                print(warn_y+f' Could not read filter/name/date from {fits_file} — skipping')
                continue

            # P60/SEDM-specific name and filter normalisation
            if str(sci_hdr.get('TELESCOP','')).strip() in SEDM:
                name = clean_object_name(name) or name      # [OBJ] canonical
                filt = p60_filt_map.get(filt, filt)
            else:
                name = clean_object_name(name) or name

            if 'p_Astrodon' in filt:
                filt = FILTERS[filt.split('p_')[0]]

            if filt not in filts:
                continue

            # Optionally archive raw images for the full LC
            archive_lc_images(sub_file, name, base_dir=os.path.dirname(sub_file[0]))

            if len(sub_file) > 1:
                args.stack = True

            # Build the expected output filename to support new_only mode.
            # Some facilities (e.g. TJO stacks) carry a date-only DATE-OBS —
            # treat a missing time part as midnight instead of crashing.
            _ds = str(date_obs)
            if len(_ds) >= 19 and _ds[10] in 'T ':
                _dt = datetime.timedelta(
                    hours=int(_ds[11:13]),
                    minutes=int(_ds[14:16]),
                    seconds=float(_ds[17:31]))
            else:
                _dt = datetime.timedelta(0)
            final_name = f'{name}_{filts[filt]}{_ds[:10]}_{_dt.seconds}_photometry.txt'
            final_name_stk = final_name.replace('photometry', 'stacked_photometry')

            if new_only:
                phot_dir = f'{data1_path}photometry'
                existing = os.listdir(phot_dir) if os.path.exists(phot_dir) else []
                if final_name in existing or final_name_stk in existing:
                    print(info_b+f' {name} {filter_} {date_obs} already measured — skipping')
                    continue

            # Seeing filter for stacks
            if len(sub_file) > 1:
                sub_file = check_seeing(sub_file, sp_logger=sp_logger)
                if len(sub_file) == 0:
                    print(warn_r+' No images with seeing < 5 — skipping')
                    continue
                if len(sub_file) == 1:
                    print(warn_y+' Only one image with seeing < 5 — proceeding as single')
                    args.stack = False

            # Run the full pipeline via the clean helper
            try:
                print(info_g+f' Starting reduction sequence on {sub_file[0]}')
                runlog.begin(sub_file[0] if len(sub_file)==1 else ', '.join(map(str,sub_file)), mode='folder')
                sub_obj = subtracted_phot(ims=sub_file, args=args)
                if sub_obj.sys_exit:
                    continue
                result = _run_pipeline(sub_obj, sp_logger)
                if result is not None:
                    final_phot.append(result)
                    remaining = len(fits_files) - file_idx
                    print(colored(f'✓ COMPLETED: {file_idx}/{len(fits_files)} | {remaining} remaining', 'green'))
                    try: progress_update(0, last=f"{result['filt']}={result['mag']:.2f}")
                    except Exception: pass
            except Exception as e:
                print(warn_r+f' Unhandled exception on {sub_file[0]}: {e}')
            finally:
                progress_update(1)

    progress_close()

    f_time_end = time.time()
    f_time_total = f_time_end - f_time_start
    time_unit = 'minutes' if f_time_total > 60 else 'seconds'
    if f_time_total > 60:
        f_time_total /= 60
    print(colored('---------------------------------------------------------------------------------------------','blue'))
    print(info_g+f' Finished {filter_} in {np.round(f_time_total, 2)} {time_unit}')
    if final_phot:
        print(tabulate(
            pd.DataFrame(final_phot, columns=final_phot[0].keys()).sort_values(by=['obj','mjd']),
            headers='keys', tablefmt='psql'))
    return final_phot
            
FILTS=[]



for filt_band in args.bands_to_process:
    FILTS.append(FILTERS[filt_band]) 

# per-image run log (run_logs/runlog_YYYYMMDD.jsonl): what was requested, the
# stage each image reached, the photometry and its sanity flags.  Summarise
# with subphot_run_summary.py.
runlog.install(path)   # pipeline root: data1_path may be '' and resolve to cron's cwd '/'
runlog.instrument(subtracted_phot)

if len(args.ims)>0:
    # print('gfefd',args.ims)
    final_phot=[]
    ims = args.ims
    # flatten fpack-compressed / multi-extension inputs (e.g. LCOGT .fz) before any
    # name munging — downstream code assumes simple single-HDU files named *.fits
    for _k,_im in enumerate(ims):
        if str(_im).endswith(('.fz','.fits.gz')):
            _src = _im if os.path.exists(_im) else data1_path+str(_im)
            try:
                _flat = flatten_multiext_fits(_src, path+'trimmed_sci_imgs')
                ims[_k] = os.path.relpath(_flat, data1_path)
                runlog.alias(ims[_k], _im)
                print(info_g+f' Flattened compressed input {_im} -> {ims[_k]}')
            except Exception as _e:
                print(warn_y+f' Could not flatten {_im}: {_e}')
    ims_path = '/'.join(ims[0].split('/')[:-1])
        
    if args.stack==False:
        if '*' in ims[0]:
            ims_start = ims[0].split('*')[0]
            ims=[]
            print(info_g+' Searching for images in '+data1_path+ims_path+' starting with '+ims_start)
            [ims.append(f) for f in os.listdir(data1_path+ims_path) if f.endswith('.fits')==True and ims_start.split('/')[-1] in f and not f.startswith('.') and 'tpv' not in f]
            ims = list(np.sort(ims))
            print(info_g+f" Found {len(ims)} images to process")


        # sys.exit()

        # print(ims)
        progress_start(len(ims), desc='images')
        for i, ims_file in enumerate(ims, 1):
            progress_set(i-1)
            runlog.begin(ims_file, mode='single')

            image = re.sub('.fits','',ims_file)
            if ims_path not in ims_file:
                image=ims_path+'/'+re.sub('.fits','',ims_file)

            fits_hdu = fits.open(data1_path+image+'.fits')[0].header
            # filter keyword differs per telescope — use the header_kw entry when known
            _fac_key = {'LT':'Liverpool Telescope','SEDM':'SEDM-P60','NTT':'ESO-NTT'}.get(args.telescope_facility,args.telescope_facility)
            _filt_kw = header_kw.get(_fac_key,{}).get('filter','FILTER')
            fits_filt = str(fits_hdu.get(_filt_kw,fits_hdu.get('FILTER','?')))
            fits_obj  = str(fits_hdu.get('OBJECT','?'))

            if args.bands==['All'] and args.sci_names == ['All']:
                pass

            elif args.bands!=['All'] and args.sci_names == ['All']:
                if fits_filt in FILTS:
                    pass
                else:
                    continue
            
            elif args.bands==['All'] and args.sci_names != ['All']:
                if len(args.sci_names)>1 and (fits_obj in args.sci_names or any(x in fits_obj for x in args.sci_names)):
                    pass
                elif len(args.sci_names)==1 and fits_obj in args.sci_names:
                    pass
                else:
                    continue
            
            elif args.bands!=['All'] and args.sci_names != ['All']:
                # print(fits_obj in args.sci_names,fits_obj,fits_filt,FILTS,fits_filt in FILTS)
                # if (fits_obj in args.sci_names or any(x in fits_obj for x in args.sci_names)) and fits_filt in FILTS:
                if len(args.sci_names)>1 and (fits_obj in args.sci_names or any(x in fits_obj for x in args.sci_names)) and fits_filt in FILTS:
                    pass
                elif len(args.sci_names)==1 and fits_obj in args.sci_names and fits_filt in FILTS:
                    pass
                else:
                    continue

            print(colored('════════════════════════════════════════════════════════════════════════════════════════════════','magenta'))
            print(colored(f'  ⭐ PROCESSING FILE {i} OF {len(ims)} [{i}/{len(ims)}]', 'cyan'))
            print(info_g+f' Performing image subtraction on {image}, {fits_filt} filter, {fits_obj}')
            print(colored('════════════════════════════════════════════════════════════════════════════════════════════════','magenta'))

            # [SRV] archive the raw frame for the full LC — the single-image
            # path (used per new source on the server) never did this before
            _arch_name = clean_object_name(fits_obj) or fits_obj   # [OBJ] canonical
            archive_lc_images([data1_path+image+'.fits'], _arch_name,
                              base_dir=os.path.dirname(data1_path+image))

            sub_obj = subtracted_phot(ims=[image],args=args)
            sys_exit=sub_obj.sys_exit


            if sys_exit==True:
                pass
            else:

                if args.unsubtract==True:
                    sub_obj.to_subtract=False

                sys_exit=sub_obj.sys_exit

                if sys_exit==True:
                    pass
                else:

                    # SEDM frames are only re-solved when the instrument flagged the WCS
                    # as unsolved (IQWCS=False): then the header holds just the pointing
                    _unsolved_sedm = sub_obj.telescope in SEDM and getattr(sub_obj,'wcs_solved',None) is False
                    if (args.redo_astrometry and sub_obj.telescope not in SEDM) or _unsolved_sedm:
                        sub_obj.resolve_astrometry()
                    sub_obj.bkg_subtract()
                    sys_exit=sub_obj.sys_exit
                    if sys_exit==True:
                        pass
                    else:

                        sub_obj.remove_cosmic()
                        sys_exit=sub_obj.sys_exit
                        if sys_exit==True:
                            pass
                        else:

                            if args.use_swarp==True: sub_obj.swarp_ref_align()
                            else: sub_obj.py_ref_align()
                            sys_exit=sub_obj.sys_exit
                            if sys_exit==True:
                                pass
                            else:

                                if args.use_psfex==True: sub_obj.psfex_convolve_images()
                                else: sub_obj.py_convolve_images()
                                sys_exit=sub_obj.sys_exit
                                if sys_exit==True:
                                   pass
                                else:

                                    if args.relative_flux!=True:
                                        sub_obj.gen_ref_cat()
                                    else:
                                        pass
                                    sys_exit=sub_obj.sys_exit
                                    if sys_exit==True:
                                        pass
                                    else:

                                        sub_obj.combine_psf()
                                        sys_exit=sub_obj.sys_exit
                                        if sys_exit==True:
                                            pass
                                        else:

                                            sub_obj.get_zeropts()
                                            sys_exit=sub_obj.sys_exit
                                            if sys_exit==True:
                                                pass
                                            else:

                                                sub_obj.scaled_subtract()
                                                sys_exit=sub_obj.sys_exit
                                                if sys_exit==True:
                                                    pass
                                                else:

                                                    final_phot.append(sub_obj.get_photometry())
                                                    sys_exit=sub_obj.sys_exit
                                                    if sys_exit==True:
                                                        pass
                                                    else:
                                                        remaining = len(ims) - i
                                                        print(colored(f'✓ COMPLETED: {i}/{len(ims)} | {remaining} remaining', 'green'))
                                                        try: progress_update(0, last=f"{final_phot[-1]['filt']}={final_phot[-1]['mag']:.2f}")
                                                        except Exception: pass

                                                        if args.upfritz==True or args.upfritz_f==True:
                                                            sub_obj.upload_phot()
                                                            sys_exit=True

                    if args.cleandirs!=False:
                        sub_obj.clean_directory()
                    

        #sp_logger.info a table of the photometry
        progress_set(len(ims)); progress_close()
        if len(final_phot)>0: print(pd.DataFrame(final_phot,columns=final_phot[0].keys()).sort_values(by=['obj','mjd']))                                          

                    
    elif args.stack==True:
        # print(ims)
        if '*' in ims[0]:
            ims_start = ims[0].split('*')[0]
            # print(ims[0])
            # print(ims_start)
            ims=[]
            print(info_g+' Searching for images in '+data1_path+ims_path+' starting with '+ims_start)

            # for f in os.listdir(data1_path+ims_path):
            #     if f.endswith('.fits')==True and ims_start.split('/')[-1] in f:
            # print(glob.glob(data1_path+ims_path+'/*'))
            [ims.append(f) for f in glob.glob(data1_path+ims_path+'/'+ims_start.split('/')[-1]+'*')]#) if f.endswith('.fits')==True ]#and ims_start.split('/')[-1] in f]
            ims = list(np.sort(ims))
            # print(ims)
            # sys.exit(1)
            # []
            print(info_g+f" Found {len(ims)} images to process")

        print(colored(f'---------------------------------------------------------------------------------------------','blue'))
        print(info_g+f' Performing image subtraction on {", ".join(ims)}')
        print('')

        # print(ims)
        ims = check_seeing(ims,sp_logger=sp_logger)
        if len(ims)==0: print(warn_r+' No images with seeing < 5, continuing with single image'); sys.exit()
        if len(ims)==1: print(warn_y+' Only one image with seeing < 5, exiting'); args.stack=False

        # for i in range(len(ims)):
        #     image = re.sub('.fits','',ims[i])
        #     if ims_path not in ims[i]:
        #         image=ims_path+'/'+re.sub('.fits','',ims[i])

        #     fits_hdu = fits.open(data1_path+image+'.fits')[0].header
        #     try:
        #         fits_filt,fits_obj = fits_hdu['FILTER1'],fits_hdu['OBJECT']
        #         # print(info_g+' Found FILTER1 and OBJECT in header',fits_filt,fits_obj,fits_obj in args.sci_names)
        #     except:
        #         print(warn_y+' No FILTER1 or OBJECT in header, trying FILTER and OBJECT')

        #         try:
        #             fits_filt,fits_obj = fits_hdu['FILTER'],fits_hdu['OBJECT']
        #             if 'p_Astrodon' in fits_filt:
        #                 fits_filt = fits_filt.split('p_')[0]
        #             print(info_g+' Found FILTER and OBJECT in header',fits_filt,fits_obj)
        #         except:
        #             print(warn_r+' No FILTER or OBJECT in header, skipping image')
        #             continue



        runlog.begin(', '.join(map(str,ims)), mode='stack')
        sub_obj = subtracted_phot(ims=ims,args=args)
        sys_exit=sub_obj.sys_exit

        if sys_exit==True:
            pass
        else:

            if args.unsubtract==True:
                sub_obj.to_subtract=False

            sys_exit=sub_obj.sys_exit

            if sys_exit==True:
                pass
            else:

                # SEDM frames are only re-solved when the instrument flagged the WCS
                # as unsolved (IQWCS=False): then the header holds just the pointing
                _unsolved_sedm = sub_obj.telescope in SEDM and getattr(sub_obj,'wcs_solved',None) is False
                if (args.redo_astrometry and sub_obj.telescope not in SEDM) or _unsolved_sedm:
                    sub_obj.resolve_astrometry()
                sub_obj.bkg_subtract()
                sys_exit=sub_obj.sys_exit
                if sys_exit==True:
                    pass
                else:

                    sub_obj.remove_cosmic()
                    sys_exit=sub_obj.sys_exit
                    if sys_exit==True:
                        pass
                    else:

                        if args.use_swarp==True: sub_obj.swarp_ref_align()
                        else: sub_obj.py_ref_align()
                        sys_exit=sub_obj.sys_exit
                        if sys_exit==True:
                            pass
                        else:

                            if args.use_psfex==True: sub_obj.psfex_convolve_images()
                            else: sub_obj.py_convolve_images()
                            sys_exit=sub_obj.sys_exit
                            if sys_exit==True:
                                pass
                            else:

                                sub_obj.gen_ref_cat()
                                sys_exit=sub_obj.sys_exit
                                if sys_exit==True:
                                    pass
                                else:

                                    sub_obj.combine_psf()
                                    sys_exit=sub_obj.sys_exit
                                    if sys_exit==True:
                                        pass
                                    else:

                                        sub_obj.get_zeropts()
                                        sys_exit=sub_obj.sys_exit
                                        if sys_exit==True:
                                            pass
                                        else:

                                            sub_obj.scaled_subtract()
                                            sys_exit=sub_obj.sys_exit
                                            if sys_exit==True:
                                                pass
                                            else:

                                                final_phot.append(sub_obj.get_photometry())
                                                sys_exit=sub_obj.sys_exit
                                                if sys_exit==True:
                                                    pass
                                                else:

                                                    if args.upfritz==True or args.upfritz_f==True:
                                                        sub_obj.upload_phot()
                                                        sys_exit=True

                if args.cleandirs!=False:
                    sub_obj.clean_directory()

        #sp_logger.info a table of the photometry
        if len(final_phot)>0: print(pd.DataFrame(final_phot,columns=final_phot[0].keys()).sort_values(by=['obj','mjd']))



elif len(args.folder)>0 and args.plc[0]!=[None]: 
    print(colored('---------------------------------------------------------------------------------------------','blue'))
    folder_path = '/'.join(args.folder[0].split('/')[:-1])
    if args.multipro=='inactive':
        data_folders = args.folder
        if args.mroundup==True:
            print(info_g+f" Performing morning round up on photometry taken last night {DATE}")

        all_phot = []
        for d in range(len(data_folders)):
            t_start = time.time()
            sub_folder = re.sub(str(folder_path+'/'),'',data_folders[d])
            # print(sub_folder)
            multi_sub_obj = multi_subtract(sub_folder,folder_path)
            folder=multi_sub_obj.folder
            FOLDER=multi_sub_obj.FOLDER
            data_dict = multi_sub_obj.analyse_folder()
            # print(data_dict)
            # print(data_dict.keys())
            for filt in data_dict.keys():
                filt_dict =  {'fits':data_dict[filt],'filter':filt,'FOLDER':FOLDER,'new_only':args.new_only}
                all_phot.append(run_subtraction(filt_dict))


            t_end=time.time()

            print(colored('---------------------------------------------------------------------------------------------','blue'))
            print(info_g+f" Completed photometry on {sub_folder} in {','.join(args.bands_to_process)} in {np.round(t_end-t_start,2)} seconds ")
            print(colored('---------------------------------------------------------------------------------------------','blue'))
        
        #combine all photometry into one list
        all_phot = [item for sublist in all_phot for item in sublist]
        #sp_logger.info a table of the photometry
        if len(all_phot)>0: 
            # print(pd.DataFrame(all_phot,columns=all_phot[0].keys()).sort_values(by=['obj','mjd']))
            print(tabulate(pd.DataFrame(all_phot,columns=all_phot[0].keys()).sort_values(by=['obj','mjd']), headers='keys', tablefmt='psql'))







    elif args.multipro=='mpi':
        data_folders = args.folder


        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        


        if args.mroundup==True:
            print(info_g+f" Performing multiprocessing on photometry taken last night {DATE}")


        for d in range(len(data_folders)):
            DATA_DICT = {'fits':[]}
            t_start = time.time()
            sub_folder = re.sub(str(folder_path+'/'),'',data_folders[d])
            mp_sub_obj = multi_subtract(sub_folder,folder_path)
            folder=mp_sub_obj.folder
            FOLDER = mp_sub_obj.FOLDER
            data_dict = mp_sub_obj.analyse_folder()

            for i in range(len(args.bands_to_process)):
                if i==rank:
                    phot_band=FILTERS[args.bands_to_process[i]]
                    # print(i,phot_band)
                    print(colored('---------------------------------------------------------------------------------------------','blue'))
                    print(info_g+f' Node {rank} is measuring photometry in {phot_band}')
                    print(colored('---------------------------------------------------------------------------------------------','blue'))

                    # DATA_DICT = {'fits':data_dict[phot_band],'filter':phot_band,'FOLDER':FOLDER,'new_only':args.new_only}

                    if len(DATA_DICT['fits'])>0:
                        run_subtraction(DATA_DICT)
                        t_end=time.time()
                        print(colored('---------------------------------------------------------------------------------------------','blue'))
                        print(info_g+f" Completed photometry in {np.round(t_end-t_start,2)} seconds ")
                        print(colored('---------------------------------------------------------------------------------------------','blue'))
                    else:
                        print(colored('---------------------------------------------------------------------------------------------','blue'))
                        print(warn_y+f' No fits files found for {phot_band}, passing')
                        print(colored('---------------------------------------------------------------------------------------------','blue'))
            
    
    elif args.multipro=='pools':
        data_folders = args.folder
        if args.mroundup==True:
            print(colored('---------------------------------------------------------------------------------------------','blue'))
            print(info_g+f"Performing multiprocessing morning round up on photometry taken last night {DATE}")
            print(colored('---------------------------------------------------------------------------------------------','blue'))


        for d in range(len(data_folders)):
            t_start = time.time()
            sub_folder = re.sub(str(folder_path+'/'),'',data_folders[d])
            mp_sub_obj = multi_subtract(sub_folder,folder_path)
            folder=mp_sub_obj.folder
            FOLDER = mp_sub_obj.FOLDER
            data_dict = mp_sub_obj.analyse_folder()

            g_dict = {'fits':data_dict['SDSS-G'],'filter':'SDSS-G','FOLDER':FOLDER,'new_only':args.new_only}
            r_dict = {'fits':data_dict['SDSS-R'],'filter':'SDSS-R','FOLDER':FOLDER,'new_only':args.new_only}
            i_dict = {'fits':data_dict['SDSS-I'],'filter':'SDSS-I','FOLDER':FOLDER,'new_only':args.new_only}
            z_dict = {'fits':data_dict['SDSS-Z'],'filter':'SDSS-Z','FOLDER':FOLDER,'new_only':args.new_only}
            u_dict = {'fits':data_dict['SDSS-U'],'filter':'SDSS-U','FOLDER':FOLDER,'new_only':args.new_only}

            obj_array = np.array([g_dict,r_dict,i_dict,z_dict,u_dict],dtype=object)

            #multiprocessing bit
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(run_subtraction,obj_array)

            t_end=time.time()

            print(colored('---------------------------------------------------------------------------------------------','blue'))
            print(info_g+f"Completed photometry in {np.round(t_end-t_start,2)} seconds ")
            print(colored('---------------------------------------------------------------------------------------------','blue'))



    if args.pros_job_id!=None:
        if not os.path.exists(data1_path+'nightly_routine_logs/'+DATE):os.mkdir(data1_path+'nightly_routine_logs/'+DATE)
        os.system(f'cp /mnt/aridata1/users/arikhind/phot_data/nightly_routine_logs/SPNR.log {data1_path}nightly_routine_logs/{DATE}/SPNR_{args.pros_job_id}_scron.log')
        print(info_g+f' Copied log file for job {args.pros_job_id} to {data1_path}nightly_routine_logs/{DATE}/SPNR_{args.pros_job_id}_scron.log')

        # os.system

def load_photometry(folder,name,filters=['All']):

    if filters[0]!='All':
        filter_list = filters
    else:
        filter_list = ['g','r','i','z','u','R']

    all_phot_files = [data1_path+folder+'/'+f for f in os.listdir(data1_path+folder) if re.search(name,f)]
    # all_phot_files = [path+folder+'/'+f for f in os.listdir(path+folder) if name+'_' in f]


    phot_dict = {}
    all_phot_data = []
    for filt in filter_list:
        filt_arr = [f for f in all_phot_files if re.search('_'+filt+'2',f)] #find all files that contain  "_{filter name}2"
        if len(filt_arr)>0:
            for file_ in filt_arr:
                d1 = open(file_, 'r')
                data1 = np.asarray(d1.readlines()[0].split(' ')[:-1])
                all_phot_data.append(data1)

        phot_dict[filt] = all_phot_data
                # print(data1)

    return phot_dict
        # phot_dict[filt] = 




    




if args.plc[0]!=None:
    if len(args.plc)<2:
        print(warn_y+'Please provide the folder and name for a given object in that order')
    
    else:
        print(info_g+' Creating light curve from photometry found in %s for object %s'%(args.plc[0],args.plc[1:]))
        phot_folder = args.plc[0] #folder containing the photometry
        phot_name = [i for i in args.plc[1:] if i!='save'] #name of the object

        for obj_name in phot_name:
            #find all text files in the photometry folder that contain obj_name
            phot_files = [data1_path+f for f in os.listdir(data1_path+phot_folder) if re.search(obj_name,f)]
            print(info_g+f'Found {len(phot_files)} files for {obj_name}')

            filters = {'g':[],'r':[],'i':[],'z':[],'u':[],'R':[],'B':[]} #dictionary to store the photometry for each filter
            for f in phot_files:
                filter_name = f.split(obj_name+'_')[1][0]
                filters[filter_name].append(f)

            #now we have a dictionary of the photometry for each filter
            #we can now plot the light curve
            for key in filters.keys():
                if len(filters[key])>0:
                    print(info_g+f' Plotting light curve for {obj_name} in {key}')
                    # plot_light_curve(filters[key],obj_name,key)
                else:
                    print(warn_y+f' No photometry found for {obj_name} in {key}')
                    continue

                fig,axs = plt.subplots(figsize=(10,10))
                lc_data = load_photometry(phot_folder,obj_name,filters=[key])
                df = pd.DataFrame(lc_data[key],columns=['name','filt','mjd','mag','mag_err','lim_mag','ra_deg','dec_deg','expt','airmass','flux','flux_err'])
                df = df.sort_values(by=['mjd'])

                x = [float(val) for ind,val in enumerate(df['mjd']) if float(df['mag'].iloc[ind])<30]
                

                if args.plc_space in ['mag','m','MAG','M']:
                    y = [float(val) for ind,val in enumerate(df['mag']) if float(val)<30]
                    y_err = [float(val) for ind,val in enumerate(df['mag_err']) if float(df['mag'].iloc[ind])<30]
                    axs.invert_yaxis()
                    axs.set_ylabel('Magnitude')
                else:
                    y = [float(val) for ind,val in enumerate(df['flux']) if df['mag'][ind]<30]
                    y_err = [float(val) for ind,val in enumerate(df['flux_err']) if df['mag'][ind]<30]
                    axs.set_ylabel('Relative Flux')

                # for i in range(len(x)):
                #     if y[i]/y_err[i]>1:
                #         axs.scatter(x[i],y[i],c='blue')
                #     else:
                #         axs.scatter(x[i],y[i],c='k')
                axs.scatter(x,y,c='k')
                axs.errorbar(x,y,yerr=y_err,fmt='o',c='k',capsize=3)
                axs.set_xlabel('MJD')
                axs.set_title(f'Light curve for {obj_name} in {key}')
                # axs.set_xlim(min(x)-0.015,max(x)+0.015)
                # print(x)
                # if axs.set_ylim(-2,max(y)+2.5)
                # axs.axhline(y=0,c='r',ls='--')
                #show the date on the x axis at 45 degree angle
                xlabels = [np.round(float(val),3) for ind,val in enumerate(df['mjd'])]
                axs.set_xticklabels(xlabels,rotation=45,fontsize=10)
                # axs.set_yticks([-1,0,1,2,3,4,5,6],fontsize=11)
                
                # axs.set_xticks(x,rotation=90)
                axs.grid('x','major')
                axs.grid('y')
                # axs.set_xticks(x,rotation=90)
                plt.show()

                if 'save' in args.plc:
                    print(info_g+f' Saving light curve for {obj_name} in {key}')
                    fig.savefig(data1_path+phot_folder+'/'+obj_name+'_'+key+'_light_curve.png',dpi=300)
                    plt.close(fig)

                    print(info_g+f' Saving photometry for {obj_name} in {key}')
                    df.to_csv(data1_path+phot_folder+'/'+obj_name+'_'+key+'_photometry.csv',index=False)


if args.mroundup==True:
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

    # LT_proposals = 'JL23A05 JZ21B01 JL23A06 JL23A07 JL23B05 JL23B06 JL24A04 JL24A09'

    # morning_logs = glob.glob(data1_path+'morning_rup_logs/'+DATE+'/*')
    # #find the log with DATE in the name
    # for log in morning_logs:
    #     if DATE in log:
    #         mlog = log
    #         break
    # # print(mlog)

    # try:
    #     os.system(f'python3 {path}subphot_morning_email.py -p '+LT_proposals+' -e K.C.Hinds@2021.ljmu.ac.uk -mlog '+mlog)
    #     # d.a.perley@ljmu.ac.uk J.L.Wise@2022.ljmu.ac.uk A.M.Bochenek@2023.ljmu.ac.uk')
    # except Exception as e:
    #     print('Error with morning email script', e)