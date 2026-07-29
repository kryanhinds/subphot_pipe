#!/usr/bin/env python3
"""
server_preflight.py — verify a subphot_pipe installation before it goes live.

Checks the things that actually broke on the server: credential paths that
create stray '~' / '/' directories, missing localized config files, missing
external binaries, unwritable output directories, logging that captures
nothing, and entire_lc_imgs archiving.

    python server_preflight.py            # checks only
    python server_preflight.py --smoke /path/to/one_image.fits
                                          # + a real single-image reduction

Exit status is 0 when nothing FAILED (warnings are allowed), 1 otherwise.
"""
import os
import re
import shutil
import subprocess
import sys

OK, WARN, FAIL = 'ok', 'warn', 'fail'
_RESULTS = []


def _log(status, title, detail=''):
    _RESULTS.append((status, title, detail))
    mark = {OK: '  ok  ', WARN: ' warn ', FAIL: ' FAIL '}[status]
    colour = {OK: '\033[32m', WARN: '\033[33m', FAIL: '\033[31m'}[status]
    print(f'{colour}[{mark}]\033[0m {title}' + (f'\n         {detail}' if detail else ''))


def check_credentials():
    try:
        import subphot_credentials as cred
    except Exception as e:
        _log(FAIL, 'import subphot_credentials', str(e))
        return None
    _log(OK, 'import subphot_credentials')

    for name in ('path', 'data1_path'):
        raw = getattr(cred, name, None)
        if raw is None:
            _log(FAIL, f'{name} defined', 'missing from subphot_credentials.py')
            continue
        detail = f'{name} = {raw!r}'
        if '~' in str(raw):
            _log(WARN, f'{name} contains "~"',
                 detail + '\n         (auto-expanded at runtime, but set an absolute path)')
        elif not os.path.isabs(str(raw)):
            _log(WARN, f'{name} is not absolute', detail)
        elif not str(raw).endswith('/'):
            _log(WARN, f'{name} has no trailing slash',
                 detail + '\n         (auto-fixed at runtime; set it properly anyway)')
        else:
            _log(OK, f'{name} well-formed', detail)

        resolved = os.path.abspath(os.path.expanduser(str(raw)))
        if not os.path.isdir(resolved):
            _log(FAIL, f'{name} exists', resolved)
        elif not os.access(resolved, os.W_OK):
            _log(FAIL, f'{name} writable', resolved)
        else:
            _log(OK, f'{name} exists and is writable', resolved)
    return cred


def check_binaries(cred):
    for attr, label in (('sex_path', 'SExtractor'), ('swarp_path', 'SWarp'),
                        ('psfex_path', 'PSFEx'), ('panstamps_path', 'panstamps')):
        binp = getattr(cred, attr, None)
        if not binp:
            _log(WARN, f'{label} path set', f'{attr} missing/empty')
            continue
        resolved = os.path.expanduser(str(binp))
        if os.path.isfile(resolved) and os.access(resolved, os.X_OK):
            _log(OK, f'{label} executable', resolved)
        elif shutil.which(resolved):
            _log(OK, f'{label} on PATH', shutil.which(resolved))
        else:
            _log(FAIL, f'{label} executable', f'{attr} = {binp!r} not found/executable')


def check_configs(cred):
    root = os.path.abspath(os.path.expanduser(str(getattr(cred, 'path', '.'))))
    cfg_dir = os.path.join(root, 'config_files')
    if not os.path.isdir(cfg_dir):
        _log(FAIL, 'config_files/ present', cfg_dir)
        return
    required = ['align_sex.config', 'align_sex.conv', 'align_temp.param',
                'config.swarp', 'config_comb.swarp', 'config_resize.swarp',
                'default.param', 'temp.param', 'prepsfex.sex',
                'psfex_conf.psfex', 'sex.config', 'sex.conv']
    missing = [f for f in required if not os.path.isfile(os.path.join(cfg_dir, f))]
    if missing:
        _log(FAIL, 'config_files complete', 'missing: ' + ', '.join(missing))
    else:
        _log(OK, 'config_files complete', f'{len(required)} files in {cfg_dir}')

    # file references inside the configs must resolve on THIS machine.
    # Absolute ones are checked as-is; relative ones against the pipeline root
    # (SExtractor/SWarp are invoked from there) and against config_files/.
    _exts = ('.param', '.conv', '.config', '.sex', '.psfex', '.swarp', '.conf')
    stale = []
    for fname in required:
        fpath = os.path.join(cfg_dir, fname)
        if not os.path.isfile(fpath):
            continue
        try:
            lines = open(fpath, errors='ignore').read().splitlines()
        except Exception:
            continue
        for line in lines:
            line = line.split('#')[0].strip()
            if not line:
                continue
            for tok in line.split()[1:]:
                tok = tok.strip('"\'')
                if not tok.endswith(_exts):
                    continue
                cands = [tok] if os.path.isabs(tok) else [
                    os.path.join(root, tok), os.path.join(cfg_dir, tok)]
                if not any(os.path.exists(c) for c in cands):
                    stale.append(f'{fname}: {tok}')
    if stale:
        _log(FAIL, 'config_files references resolve',
             'referenced files not found:\n         ' + '\n         '.join(sorted(set(stale))[:8]))
    else:
        _log(OK, 'config_files references resolve')


def check_python_deps():
    mods = ['numpy', 'scipy', 'astropy', 'photutils', 'pandas', 'matplotlib',
            'astroscrappy', 'requests', 'termcolor', 'tabulate', 'image_registration']
    missing = []
    for m in mods:
        try:
            __import__(m)
        except Exception:
            missing.append(m)
    if missing:
        _log(FAIL, 'python dependencies', 'missing: ' + ', '.join(missing))
    else:
        _log(OK, 'python dependencies', f'{len(mods)} modules import cleanly')


def check_runtime_dirs(cred):
    """The pipeline creates these on demand; confirm they are creatable."""
    root = os.path.abspath(os.path.expanduser(str(getattr(cred, 'data1_path', '.'))))
    names = ['entire_lc_imgs', 'aligned_images', 'bkg_subtracted_science',
             'convolved_sci', 'convolved_ref', 'convolved_psf', 'out',
             'ps_catalogs', 'ref_imgs', 'temp_config_files', 'photometry_date']
    bad = []
    for n in names:
        d = os.path.join(root, n)
        try:
            os.makedirs(d, exist_ok=True)
            if not os.access(d, os.W_OK):
                bad.append(n + ' (not writable)')
        except Exception as e:
            bad.append(f'{n} ({e})')
    if bad:
        _log(FAIL, 'runtime directories', ', '.join(bad))
    else:
        _log(OK, 'runtime directories', f'{len(names)} present/creatable under {root}')


def check_stray_dirs(cred):
    """The bugs this release fixes: stray '0', '~' and '/'-rooted folders."""
    root = os.path.abspath(os.path.expanduser(str(getattr(cred, 'path', '.'))))
    found = []
    for stray in ('0', '~'):
        p = os.path.join(root, stray)
        if os.path.isdir(p):
            found.append(p)
    if os.path.isdir('/~'):
        found.append('/~')
    home_tilde = os.path.expanduser('~/~')
    if os.path.isdir(home_tilde):
        found.append(home_tilde)
    if found:
        _log(WARN, 'no stray directories',
             'left over from the OLD version — safe to delete:\n         '
             + '\n         '.join(found))
    else:
        _log(OK, 'no stray directories')


def check_fritz(cred):
    tok = getattr(cred, 'token', None)
    if not tok or len(str(tok)) < 10:
        _log(WARN, 'Fritz token present', 'token missing/short — uploads (-up/-upf) will fail')
    else:
        _log(OK, 'Fritz token present', f'{len(str(tok))} chars (value not shown)')


def check_git():
    try:
        br = subprocess.run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                            capture_output=True, text=True, timeout=15)
        st = subprocess.run(['git', 'status', '--porcelain'],
                            capture_output=True, text=True, timeout=15)
    except Exception as e:
        _log(WARN, 'git state', str(e))
        return
    branch = br.stdout.strip()
    dirty = [l for l in st.stdout.splitlines() if l and not l.startswith('??')]
    _log(OK if branch == 'main' else WARN, 'git branch', branch)
    if dirty:
        _log(WARN, 'tracked files unmodified',
             'local edits present — a pull will conflict:\n         '
             + '\n         '.join(dirty[:8]))
    else:
        _log(OK, 'tracked files unmodified')


def smoke_test(image, cred):
    """Run one real reduction and confirm the outputs this release fixes."""
    root = os.path.abspath(os.path.expanduser(str(getattr(cred, 'data1_path', '.'))))
    out_dir = 'preflight_smoke'
    cmd = [sys.executable, os.path.join(root, 'subphot_subtract.py'),
           '-i', image, '-o', out_dir, '-cut']
    print(f'\n--- smoke test: {" ".join(cmd)}\n')
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd='/')
    except Exception as e:
        _log(FAIL, 'smoke test ran', str(e))
        return
    out_path = os.path.join(root, out_dir)
    if r.returncode != 0:
        _log(FAIL, 'smoke test exit code',
             f'rc={r.returncode}\n         ' + '\n         '.join(r.stdout.splitlines()[-6:]))
    else:
        _log(OK, 'smoke test exit code', 'rc=0')

    phot = [f for f in os.listdir(out_path)] if os.path.isdir(out_path) else []
    _log(OK if any(f.endswith('photometry.txt') for f in phot) else FAIL,
         'smoke test wrote photometry', out_path)
    _log(OK if os.path.isdir(os.path.join(out_path, 'cut_outs')) else WARN,
         'smoke test wrote cutouts')
    logf = os.path.join(out_path, f'{out_dir}_log.log')
    if os.path.isfile(logf):
        n = sum(1 for _ in open(logf, errors='ignore'))
        _log(OK if n > 50 else FAIL, 'smoke test log captured output', f'{n} lines in {logf}')
    else:
        _log(FAIL, 'smoke test log captured output', f'no log at {logf}')
    lc = os.path.join(root, 'entire_lc_imgs')
    recent = []
    if os.path.isdir(lc):
        for d in os.listdir(lc):
            dd = os.path.join(lc, d)
            if os.path.isdir(dd) and os.listdir(dd):
                recent.append(d)
    _log(OK if recent else FAIL, 'smoke test archived to entire_lc_imgs',
         ', '.join(recent[:5]) if recent else 'nothing archived')
    check_stray_dirs(cred)


def main():
    image = None
    if '--smoke' in sys.argv:
        i = sys.argv.index('--smoke')
        if i + 1 < len(sys.argv):
            image = sys.argv[i + 1]

    print('=' * 72)
    print('subphot_pipe preflight')
    print(f'  python : {sys.version.split()[0]} ({sys.executable})')
    print(f'  cwd    : {os.getcwd()}')
    print('=' * 72)

    cred = check_credentials()
    if cred is not None:
        check_binaries(cred)
        check_configs(cred)
        check_python_deps()
        check_runtime_dirs(cred)
        check_stray_dirs(cred)
        check_fritz(cred)
        check_git()
        if image:
            smoke_test(image, cred)

    n_fail = sum(1 for s, _, _ in _RESULTS if s == FAIL)
    n_warn = sum(1 for s, _, _ in _RESULTS if s == WARN)
    print('=' * 72)
    print(f'{len(_RESULTS)} checks: {len(_RESULTS)-n_fail-n_warn} ok, {n_warn} warn, {n_fail} FAIL')
    if n_fail:
        print('NOT READY — fix the FAIL items above.')
    elif n_warn:
        print('Ready, with warnings worth reading.')
    else:
        print('Ready.')
    print('=' * 72)
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
