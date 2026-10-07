"""Per-image run log for the subphot pipeline.

Every image the pipeline is asked to reduce gets one record, written as JSON
lines to <pipeline root>/run_logs/runlog_YYYYMMDD.jsonl (UTC date the run
started).  A record answers:

  * what was requested           -> 'requested' (file name as given on the CLI)
  * did it produce photometry    -> status 'detection' or 'limit'
  * is the photometry reasonable -> 'reasonable' + 'flags' (sanity checks below)
  * did it fail, and where       -> status 'stopped' (a step rejected the image),
                                    'crashed' (exception), 'killed' (no end
                                    event: the process died), with 'stage' and
                                    'reason'
  * what happened on Fritz       -> a separate 'fritz' event

Each run writes a 'start' event immediately and an 'end' event when it
finishes, so a process killed by a timeout or segfault still shows up (start
without end).  subphot_run_summary.py merges the events and reports counts.

Only the standard library is used, so nothing new has to be installed.
"""
import atexit
import collections
import datetime
import fcntl
import functools
import json
import math
import os
import re
import socket
import sys
import traceback
import uuid

# ---------------------------------------------------------------------------
# Sanity checks: photometry outside these bounds is produced but flagged.
# ---------------------------------------------------------------------------
CHECKS = {
    'mag_min': 8.0,          # brighter than this is almost certainly wrong
    'mag_max': 24.0,         # fainter than this is beyond any facility we use
    'magerr_max': 0.5,       # a "detection" with a larger error is not useful
    'lim_min': 16.0,         # shallower limit than this: cloud or bad frame
    'lim_max': 24.5,         # deeper than this is not credible for our telescopes
    'zp_std_max': 0.15,      # science zeropoint scatter (mag)
    'zp_stars_min': 5,       # calibration stars left after clipping
    'seeing_max': 5.0,       # arcsec
}

# Pipeline methods -> human-readable stage names.  Order is the run order.
STAGES = collections.OrderedDict([
    ('__init__',              'Load & header checks'),
    ('resolve_astrometry',    'Astrometry'),
    ('bkg_subtract',          'Background subtraction'),
    ('remove_cosmic',         'Cosmic-ray removal'),
    ('py_ref_align',          'Image alignment'),
    ('swarp_ref_align',       'Image alignment'),
    ('py_convolve_images',    'PSF convolution'),
    ('psfex_convolve_images', 'PSF convolution'),
    ('gen_ref_cat',           'Reference catalog'),
    ('combine_psf',           'Combined PSF'),
    ('get_zeropts',           'Zeropoint calibration'),
    ('scaled_subtract',       'Scaled subtraction'),
    ('get_photometry',        'Photometry'),
])
PRE_STAGE = 'Input (before pipeline)'

_ANSI = re.compile(r'\x1b\[[0-9;]*[A-Za-z]')
_RED_WARNING = '\x1b[31m[WARNING]'

_state = {
    'dir': None,          # run_logs directory, None until install()
    'open': None,         # the record currently being processed
    'aliases': {},        # flattened/decompressed path -> name requested
    'capture': None,
    'installed': False,
}


# ---------------------------------------------------------------------------
# stdout capture: remembers recent lines so a failure can quote the warning
# the pipeline printed just before it stopped.
# ---------------------------------------------------------------------------
class _LineCapture:
    def __init__(self, stream, maxlen=400):
        self._stream = stream
        self.lines = collections.deque(maxlen=maxlen)
        self.count = 0          # total lines seen, used as a position marker
        self._partial = ''

    def write(self, data):
        try:
            text = self._partial + data
            *done, self._partial = text.split('\n')
            for line in done:
                if line.strip():
                    self.lines.append((self.count, line))
                    self.count += 1
        except Exception:
            pass
        return self._stream.write(data)

    def flush(self):
        return self._stream.flush()

    def isatty(self):
        return getattr(self._stream, 'isatty', lambda: False)()

    def fileno(self):
        return self._stream.fileno()

    def __getattr__(self, name):
        return getattr(self._stream, name)


def _mark():
    cap = _state['capture']
    return cap.count if cap is not None else 0


def _last_warning(since=0):
    """Best explanation printed since position *since*: last red warning,
    else last warning of any colour, else the last line."""
    cap = _state['capture']
    if cap is None:
        return None
    recent = [l for n, l in cap.lines if n >= since]
    for test in (lambda l: _RED_WARNING in l, lambda l: '[WARNING]' in _ANSI.sub('', l)):
        hits = [l for l in recent if test(l)]
        if hits:
            return _clean(hits[-1])
    return _clean(recent[-1]) if recent else None


def _clean(line):
    line = _ANSI.sub('', line).strip()
    line = re.sub(r'^\[(WARNING|INFO|ERROR)\]\s*::\s*', '', line)
    return line[:500]


# ---------------------------------------------------------------------------
# JSON-lines writer (append + flock, safe with several pipeline processes)
# ---------------------------------------------------------------------------
def _now():
    return datetime.datetime.now(datetime.timezone.utc)


def _sanitize(x):
    if isinstance(x, dict):
        return {str(k): _sanitize(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [_sanitize(v) for v in x]
    if hasattr(x, 'item') and not isinstance(x, (str, bytes)):   # numpy scalar
        try:
            x = x.item()
        except Exception:
            return str(x)
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if x is None or isinstance(x, (bool, int, float, str)):
        return x
    return str(x)


def _write(rec, event, **fields):
    if _state['dir'] is None:
        return
    entry = {'event': event, 'run_id': rec['run_id'], 'time': _now().isoformat(timespec='seconds')}
    if event == 'start':
        entry.update({k: rec[k] for k in ('requested', 'mode', 'host', 'pid', 'argv')})
    entry.update(fields)
    line = json.dumps(_sanitize(entry)) + '\n'
    try:
        with open(rec['log_file'], 'a', encoding='utf-8') as fh:
            fcntl.flock(fh, fcntl.LOCK_EX)
            fh.write(line)
            fh.flush()
            fcntl.flock(fh, fcntl.LOCK_UN)
    except Exception as e:
        try:
            sys.__stdout__.write(f'[WARNING] :: run log write failed ({e})\n')
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def install(root):
    """Start logging to <root>/run_logs. Call once, after stdout is set up."""
    if _state['installed']:
        return
    try:
        d = os.path.join(os.path.abspath(os.path.expanduser(root)), 'run_logs')
        os.makedirs(d, exist_ok=True)
        _state['dir'] = d
    except Exception as e:
        print(f'[WARNING] :: run log disabled, cannot create run_logs ({e})')
        return
    _state['capture'] = _LineCapture(sys.stdout)
    sys.stdout = _state['capture']
    _prev_hook = sys.excepthook

    def _hook(exc_type, exc, tb):
        rec = _state['open']
        if rec is not None and not rec['closed']:
            _end(rec, 'crashed', reason=f'{exc_type.__name__}: {exc}', where=_where(tb))
        _prev_hook(exc_type, exc, tb)
    sys.excepthook = _hook
    atexit.register(_at_exit)
    _state['installed'] = True


def alias(actual, requested):
    """Record that *actual* (e.g. a decompressed copy) was requested as *requested*."""
    _state['aliases'][os.path.basename(str(actual))] = str(requested)


def begin(requested, mode='single'):
    """Open a record for one requested image (or stack). Closes any open one."""
    _close_dangling()
    requested = str(requested)
    requested = _state['aliases'].get(os.path.basename(requested), requested)
    start = _now()
    rec = {
        'run_id': uuid.uuid4().hex[:12],
        'requested': requested,
        'mode': mode,
        'host': socket.gethostname().split('.')[0],
        'pid': os.getpid(),
        'argv': ' '.join(sys.argv),
        'start': start,
        'stage': PRE_STAGE,
        'constructed': False,
        'closed': False,
        'log_file': os.path.join(_state['dir'] or '.', f'runlog_{start:%Y%m%d}.jsonl'),
    }
    _state['open'] = rec
    _write(rec, 'start')
    return rec


def instrument(cls):
    """Wrap the pipeline class so each step records its stage and outcome."""
    for meth, label in STAGES.items():
        fn = cls.__dict__.get(meth)
        if fn is not None and not getattr(fn, '_runlogged', False):
            setattr(cls, meth, _wrap_step(fn, label, meth == '__init__'))
    fn = cls.__dict__.get('upload_phot')
    if fn is not None and not getattr(fn, '_runlogged', False):
        setattr(cls, 'upload_phot', _wrap_upload(fn))
    return cls


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))


def _where(tb):
    """Deepest traceback frame in the pipeline's own code (not a library)."""
    try:
        frames = traceback.extract_tb(tb)
        own = [f for f in frames if os.path.abspath(f.filename).startswith(_HERE)
               and os.path.basename(f.filename) != 'subphot_runlog.py']
        f = (own or frames)[-1]
        return f'{os.path.basename(f.filename)}:{f.lineno} in {f.name}'
    except Exception:
        return None


def _rec_for(obj, create=False):
    rec = getattr(obj, '_runlog_rec', None)
    if rec is None and create:
        rec = _state['open']
        if rec is None or rec['closed'] or rec['constructed']:
            ims = getattr(obj, '_runlog_ims', None)
            rec = begin(ims or 'unknown', mode='unknown')
        rec['constructed'] = True
        try:
            object.__setattr__(obj, '_runlog_rec', rec)
        except Exception:
            pass
    return rec


def _header_info(obj):
    out = {}
    for key, attr in (('object', 'sci_obj'), ('filter', 'sci_filt'), ('mjd', 'sci_mjd'),
                      ('telescope', 'telescope'), ('seeing', 'sci_seeing'),
                      ('exptime', 'sci_exp_time')):
        v = getattr(obj, attr, None)
        if v is not None:
            out[key] = v
    return out


def _wrap_step(fn, label, is_init):
    @functools.wraps(fn)
    def step(self, *a, **k):
        if is_init:
            ims = k.get('ims', a[0] if a else None)
            try:
                object.__setattr__(self, '_runlog_ims', ', '.join(map(str, ims)) if isinstance(ims, (list, tuple)) else ims)
            except Exception:
                pass
        rec = _rec_for(self, create=is_init)
        if rec is None or rec['closed'] or _state['dir'] is None:
            return fn(self, *a, **k)
        # a step may call another wrapped step; the inner one must not leave
        # its name behind if the outer step is what later fails
        prev, rec['stage'], mark = rec['stage'], label, _mark()
        nested = rec.get('depth', 0) > 0
        rec['depth'] = rec.get('depth', 0) + 1
        try:
            out = fn(self, *a, **k)
        except BaseException as e:
            rec['depth'] -= 1
            if not rec['closed']:
                _end(rec, 'crashed', obj=self, reason=f'{type(e).__name__}: {e}',
                     where=_where(sys.exc_info()[2]))
            raise
        rec['depth'] -= 1
        if rec['closed']:
            return out
        if getattr(self, 'sys_exit', False):
            _end(rec, 'stopped', obj=self, reason=_last_warning(since=mark))
        elif label == 'Photometry':
            _end_photometry(rec, self, out)
        elif nested:
            rec['stage'] = prev
        return out
    step._runlogged = True
    return step


def _wrap_upload(fn):
    @functools.wraps(fn)
    def upload(self, *a, **k):
        rec = getattr(self, '_runlog_rec', None)
        mark = _mark()
        try:
            out = fn(self, *a, **k)
        except BaseException as e:
            if rec is not None:
                _write(rec, 'fritz', status='crashed', detail=f'{type(e).__name__}: {e}')
            raise
        if rec is not None:
            _write(rec, 'fritz', status=getattr(self, 'fritz_status', 'unknown'),
                   detail=getattr(self, 'fritz_detail', None) or _last_warning(since=mark),
                   fritz_name=getattr(self, 'name', None))
        return out
    upload._runlogged = True
    return upload


def _end(rec, status, obj=None, **fields):
    if rec['closed']:
        return
    rec['closed'] = True
    info = _header_info(obj) if obj is not None else {}
    _write(rec, 'end', status=status, stage=rec['stage'],
           duration_s=round((_now() - rec['start']).total_seconds(), 1),
           **info, **fields)


def _num(x):
    try:
        x = float(x)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def _end_photometry(rec, obj, result):
    mag = _num(getattr(obj, 'mag', [None])[0])
    lim = _num(obj.mag[3]) if len(getattr(obj, 'mag', [])) > 3 else None
    err = _num(getattr(obj, 'mag_all_err', None))
    snr = _num(getattr(obj, 'SNR', None))
    is_limit = (snr is None or snr <= 3 or mag is None or mag > 40
                or err is None or lim is None or mag + err > lim)
    zp = getattr(obj, 'zp_sci', None)
    try:
        zp_n = int(len(zp))
        zp_std = _num(getattr(obj, 'zp_sci_std_eff', None))
        if zp_std is None and zp_n > 1:
            m = sum(zp) / zp_n
            zp_std = (sum((z - m) ** 2 for z in zp) / zp_n) ** 0.5
    except Exception:
        zp_n, zp_std = None, None
    seeing = _num(getattr(obj, 'sci_seeing', None))

    c, flags = CHECKS, []
    if is_limit:
        if lim is None:
            flags.append('no limiting magnitude')
        elif not c['lim_min'] <= lim <= c['lim_max']:
            flags.append(f'limit {lim:.2f} outside {c["lim_min"]}-{c["lim_max"]}')
    else:
        if not c['mag_min'] <= mag <= c['mag_max']:
            flags.append(f'mag {mag:.2f} outside {c["mag_min"]}-{c["mag_max"]}')
        if err is None or err <= 0 or err > c['magerr_max']:
            flags.append(f'mag error {err} not in (0, {c["magerr_max"]}]')
        if lim is not None and not c['lim_min'] <= lim <= c['lim_max']:
            flags.append(f'limit {lim:.2f} outside {c["lim_min"]}-{c["lim_max"]}')
    if zp_std is not None and zp_std > c['zp_std_max']:
        flags.append(f'zeropoint scatter {zp_std:.3f} > {c["zp_std_max"]}')
    if zp_n is not None and zp_n < c['zp_stars_min']:
        flags.append(f'only {zp_n} zeropoint stars (< {c["zp_stars_min"]})')
    if seeing is not None and seeing > c['seeing_max']:
        flags.append(f'seeing {seeing:.2f}" > {c["seeing_max"]}"')

    _end(rec, 'limit' if is_limit else 'detection', obj=obj,
         mag=None if is_limit else mag, magerr=None if is_limit else err,
         lim=lim, snr=snr, zp_std=zp_std, zp_stars=zp_n,
         reasonable=not flags, flags=flags)


def _close_dangling():
    rec = _state['open']
    if rec is None or rec['closed']:
        return
    if not rec['constructed']:
        _end(rec, 'not_processed',
             reason='did not reach the pipeline (skipped by band/object selection or an early exit)')
    else:
        _end(rec, 'stopped', reason=f'run ended during {rec["stage"]} without a result')


def _at_exit():
    try:
        _close_dangling()
    except Exception:
        pass
