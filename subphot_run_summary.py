#!/usr/bin/env python
"""Summarise the subphot run log (run_logs/runlog_*.jsonl).

Examples
    python subphot_run_summary.py                     # today (UTC)
    python subphot_run_summary.py --date 2026-10-03
    python subphot_run_summary.py --days 7            # last 7 UTC days incl. today
    python subphot_run_summary.py --date 2026-10-03 --csv night.csv
    python subphot_run_summary.py --email             # also email the report

Outcome of each requested image (latest attempt counts if it was rerun):
    PASS     photometry produced (detection or limit) and every sanity check passed
    FLAGGED  photometry produced but at least one sanity check failed
    FAIL     stopped by a pipeline step, crashed, or killed (start but no end)
    SKIPPED  requested but never entered the pipeline (band/object selection);
             excluded from the pass/fail fractions
    RUNNING  started recently on a live process and not finished yet

Email uses the smtp2go HTTP API with smtp2go_api_key and email_to from
subphot_credentials.py (sender: email_from if defined, else email_to).
Only the standard library and requests are needed.
"""
import argparse
import collections
import csv
import datetime
import glob
import json
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
JUMP_MAG = 1.5          # flag a detection this far from the previous one ...
JUMP_DAYS = 10.0        # ... of the same object+filter within this many days
KILLED_AFTER_H = 6.0    # a start with no end older than this is 'killed'


def _root():
    sys.path.insert(0, HERE)
    try:
        import subphot_credentials as cred
        return os.path.abspath(os.path.expanduser(cred.path)), cred
    except Exception:
        return HERE, None


def load_events(log_dir):
    runs = collections.OrderedDict()
    for fn in sorted(glob.glob(os.path.join(log_dir, 'runlog_*.jsonl'))):
        with open(fn, encoding='utf-8') as fh:
            for line in fh:
                try:
                    ev = json.loads(line)
                except ValueError:
                    continue
                r = runs.setdefault(ev['run_id'], {'run_id': ev['run_id']})
                kind = ev.pop('event', None)
                if kind == 'start':
                    r.update(ev)
                    r['start'] = ev['time']
                elif kind == 'end':
                    r.update({k: v for k, v in ev.items() if k != 'time'})
                    r['end'] = ev['time']
                elif kind == 'fritz':
                    r['fritz'] = ev.get('status')
                    r['fritz_detail'] = ev.get('detail')
    return list(runs.values())


def _pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except Exception:
        return False


def classify(r, now):
    st = r.get('status')
    if st in ('detection', 'limit'):
        return 'PASS' if r.get('reasonable', True) and not r.get('flags') else 'FLAGGED'
    if st == 'not_processed':
        return 'SKIPPED'
    if st in ('stopped', 'crashed'):
        return 'FAIL'
    # start without end
    started = _parse(r.get('start'))
    age_h = (now - started).total_seconds() / 3600 if started else 1e9
    if (age_h < KILLED_AFTER_H and r.get('host') == socket.gethostname().split('.')[0]
            and _pid_alive(r.get('pid'))):
        return 'RUNNING'
    r['status'], r['reason'] = 'killed', 'process ended without finishing (killed, timeout or hard crash)'
    return 'FAIL'


def _parse(t):
    try:
        return datetime.datetime.fromisoformat(t)
    except Exception:
        return None


def add_jump_flags(runs):
    """Flag detections far from the previous detection of the same object+filter."""
    last = {}
    for r in sorted(runs, key=lambda r: float(r.get('mjd') or 0)):
        if r.get('status') != 'detection' or r.get('mag') is None or r.get('mjd') is None:
            continue
        key = (str(r.get('object')), str(r.get('filter')))
        prev = last.get(key)
        if prev and 0 < float(r['mjd']) - float(prev['mjd']) <= JUMP_DAYS \
                and abs(r['mag'] - prev['mag']) > JUMP_MAG:
            r.setdefault('flags', []).append(
                f"jump of {r['mag'] - prev['mag']:+.2f} mag vs MJD {float(prev['mjd']):.2f}")
            r['reasonable'] = False
        last[key] = r


def select(runs, d0, d1):
    out = []
    for r in runs:
        t = _parse(r.get('start'))
        if t is not None and d0 <= t.date() <= d1:
            out.append(r)
    return out


def latest_per_request(runs):
    by = collections.OrderedDict()
    for r in sorted(runs, key=lambda r: r.get('start') or ''):
        by.setdefault(r.get('requested'), []).append(r)
    # latest attempt wins, except that a rerun which never entered the pipeline
    # (band/object selection) does not hide an earlier real result
    def pick(v):
        real = [r for r in v if r.get('status') != 'not_processed']
        return (real or v)[-1]
    return {k: pick(v) for k, v in by.items()}, {k: len(v) for k, v in by.items()}


def _pct(n, d):
    return f'{100 * n / d:.0f}%' if d else '-'


def _fmt_phot(r):
    if r.get('status') == 'detection':
        return f"{r.get('mag', float('nan')):.2f}+/-{(r.get('magerr') or float('nan')):.2f}"
    if r.get('status') == 'limit':
        return f">{r.get('lim') or float('nan'):.2f}"
    return ''


def report(runs, label, now):
    latest, attempts = latest_per_request(runs)
    for r in latest.values():
        r['outcome'] = classify(r, now)
    cnt = collections.Counter(r['outcome'] for r in latest.values())
    judged = cnt['PASS'] + cnt['FLAGGED'] + cnt['FAIL']
    phot = [r for r in latest.values() if r['outcome'] in ('PASS', 'FLAGGED')]
    fails = [r for r in latest.values() if r['outcome'] == 'FAIL']
    flagged = [r for r in latest.values() if r['outcome'] == 'FLAGGED']

    L = []
    L.append(f'subphot pipeline summary: {label}')
    L.append('=' * 72)
    L.append(f'Requested images:        {len(latest)}   ({len(runs)} runs incl. reruns)')
    L.append(f'  PASS                   {cnt["PASS"]:4d}  {_pct(cnt["PASS"], judged)}')
    L.append(f'  FLAGGED (check these)  {cnt["FLAGGED"]:4d}  {_pct(cnt["FLAGGED"], judged)}')
    L.append(f'  FAIL                   {cnt["FAIL"]:4d}  {_pct(cnt["FAIL"], judged)}')
    if cnt['SKIPPED']:
        L.append(f'  skipped (not selected) {cnt["SKIPPED"]:4d}  not counted in the fractions')
    if cnt['RUNNING']:
        L.append(f'  still running          {cnt["RUNNING"]:4d}')
    nd = sum(r.get('status') == 'detection' for r in phot)
    L.append(f'Photometry produced:     {len(phot)}  ({nd} detections, {len(phot) - nd} limits)')

    if fails:
        L.append('')
        L.append('Failures by stage:')
        for (stage, st), n in collections.Counter(
                (r.get('stage', '?'), r.get('status')) for r in fails).most_common():
            L.append(f'  {n:4d}  {stage}  ({st})')

    fr = collections.Counter(r.get('fritz') for r in latest.values() if r.get('fritz'))
    if fr:
        L.append('')
        L.append('Fritz: ' + ', '.join(f'{k} {v}' for k, v in fr.most_common()))

    def table(title, rows, cols):
        if not rows:
            return
        L.append('')
        L.append(title)
        L.append('-' * len(title))
        for r in rows:
            L.append('  ' + ' | '.join(str(c(r)) for c in cols))

    _name = lambda r: os.path.basename(str(r.get('requested')))
    _obj = lambda r: f"{r.get('object', '?')} {r.get('filter', '')}".strip()
    table('FAILURES  (file | object filt | stage | reason)', fails,
          [_name, _obj, lambda r: r.get('stage'),
           lambda r: (r.get('reason') or '') + (f" [{r['where']}]" if r.get('where') else '')])
    table('FLAGGED  (file | object filt | phot | flags)', flagged,
          [_name, _obj, _fmt_phot, lambda r: '; '.join(r.get('flags') or [])])
    table('PASSED  (file | object filt | phot | fritz)',
          [r for r in latest.values() if r['outcome'] == 'PASS'],
          [_name, _obj, _fmt_phot, lambda r: r.get('fritz') or '-'])
    reruns = {k: n for k, n in attempts.items() if n > 1}
    if reruns:
        L.append('')
        L.append(f'Rerun images (latest attempt used): {len(reruns)}')
    return '\n'.join(L), cnt, judged, list(latest.values())


def write_csv(rows, path):
    cols = ['outcome', 'requested', 'object', 'filter', 'mjd', 'telescope', 'status', 'stage',
            'reason', 'mag', 'magerr', 'lim', 'snr', 'zp_std', 'zp_stars', 'seeing', 'flags',
            'fritz', 'fritz_detail', 'start', 'end', 'duration_s', 'host', 'mode', 'run_id']
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            r = dict(r)
            r['flags'] = '; '.join(r.get('flags') or [])
            w.writerow(r)


def send_email(cred, subject, body, to=None):
    import requests
    key = getattr(cred, 'smtp2go_api_key', None)
    to = to or getattr(cred, 'email_to', None)
    sender = getattr(cred, 'email_from', None) or to
    if not key or not to:
        raise RuntimeError('smtp2go_api_key and email_to must be set in subphot_credentials.py')
    resp = requests.post('https://api.smtp2go.com/v3/email/send', timeout=30,
                         json={'api_key': key, 'to': [to], 'sender': sender,
                               'subject': subject, 'text_body': body})
    ok = resp.status_code == 200 and resp.json().get('data', {}).get('succeeded', 0) >= 1
    if not ok:
        raise RuntimeError(f'smtp2go HTTP {resp.status_code}: {resp.text[:300]}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--date', help='UTC date YYYY-MM-DD (default today)')
    ap.add_argument('--days', type=int, default=1, help='number of UTC days ending on --date')
    ap.add_argument('--log-dir', help='run_logs directory (default <pipeline path>/run_logs)')
    ap.add_argument('--csv', help='also write one row per requested image to this CSV')
    ap.add_argument('--email', action='store_true', help='email the report (smtp2go)')
    ap.add_argument('--to', help='override the recipient')
    a = ap.parse_args()

    root, cred = _root()
    log_dir = a.log_dir or os.path.join(root, 'run_logs')
    now = datetime.datetime.now(datetime.timezone.utc)
    d1 = datetime.date.fromisoformat(a.date) if a.date else now.date()
    d0 = d1 - datetime.timedelta(days=max(a.days, 1) - 1)
    label = f'{d1} (UTC)' if d0 == d1 else f'{d0} to {d1} (UTC)'

    runs = load_events(log_dir)
    # jump check needs history, so flags are computed over everything loaded
    add_jump_flags(runs)
    sel = select(runs, d0, d1)
    if not sel:
        text, cnt, judged, rows = f'subphot pipeline summary: {label}\nNo runs logged in {log_dir}.', {}, 0, []
    else:
        text, cnt, judged, rows = report(sel, label, now)
    print(text)
    if a.csv:
        write_csv(rows, a.csv)
        print(f'\nCSV written to {os.path.abspath(a.csv)}')
    if a.email:
        subj = (f'subphot {label}: {cnt.get("PASS", 0)} pass, {cnt.get("FLAGGED", 0)} flagged, '
                f'{cnt.get("FAIL", 0)} fail of {judged}') if judged else f'subphot {label}: no runs'
        send_email(cred, subj, text, to=a.to)
        print(f'\nEmailed: {subj}')


if __name__ == '__main__':
    main()
