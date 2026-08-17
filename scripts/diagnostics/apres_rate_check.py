#!/usr/bin/env python3
"""The minimal-processing ApRES tidal check: their rates against the tide rate.

WHY THIS EXISTS. apres_comparison.py chains 277 half-hour displacement
profiles end-to-end and fits the cumulative series against the tide. That
gave +3.79 mm per metre of tide at 100 m - opposite in sign to the radar
(-1.46 +/- 0.55) and the flexure model (-1.2), and it framed weeks of
"why can't multipass see the ApRES signal".

This script asks the same question with nearly no processing of mine: if
strain = A*tide, then the strain RATE of a half-hour pair is A*d(tide)/dt.
Regressing the ApRES project's own vsr_per_year_medfilt against the
CATS2008 tide rate gives A directly - no chaining, no cumulative drift,
no depth-grid handling.

Result (2026-08-17):
    GA04: A = -12.4 +/- 0.4 ue per m of tide (-1.24 mm over 100 m), R2=0.75
    GA01: A =  +9.3 +/- 1.0                  (+0.93 mm),            R2=0.27
GA04 - the only site whose bed pick holds - agrees with the radar and the
model in sign and magnitude. GA01 is the site with the 82 m bed-pick
drift and explains a quarter of its variance; it is listed for contrast,
not used.

CONCLUSION: the +3.79 was an artefact of the cumulative-chaining analysis
in apres_comparison.py (bug not yet localised - candidates: sign handling
in the chain, drift correlated with tide, the shallow grid points the vsr
fit excludes). Until it is found, THIS number is the ApRES comparison,
and all three independent estimates agree: ApRES rates, the four-line
radar network, and thin-plate flexure.

Run from the repo root: python3 scripts/diagnostics/apres_rate_check.py
"""
import csv, math, os
from datetime import datetime, timezone

def dnum(dt):
    epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
    return 719529.0 + (dt - epoch).total_seconds() / 86400.0

HERE = os.path.dirname(os.path.abspath(__file__))
tt, th = [], []
for r in csv.DictReader(open(os.path.join(HERE, 'cats2008_apres_window.csv'))):
    tt.append(float(r['datenum'])); th.append(float(r['tide_m']))

def tide_rate(d):
    import bisect
    i = bisect.bisect_left(tt, d)
    if i < 2 or i > len(tt) - 3:
        return None
    return (th[i+1] - th[i-1]) / (tt[i+1] - tt[i-1])

def main():
    for site in ['GA04', 'GA01']:
        xs, ys = [], []
        f = os.path.expanduser(f'~/projects/EAGER_ApRES/results/{site}/pair_results.csv')
        for r in csv.DictReader(open(f)):
            try:
                dt_days = float(r['dt_days']); vsr = float(r['vsr_per_year_medfilt'])
            except Exception:
                continue
            if not (0.01 < dt_days < 0.1) or not math.isfinite(vsr):
                continue
            t1 = datetime.fromisoformat(r['t1'].replace(' ', 'T'))
            if t1.tzinfo is None:
                t1 = t1.replace(tzinfo=timezone.utc)
            tr = tide_rate(dnum(t1) + dt_days / 2)
            if tr is None:
                continue
            xs.append(tr); ys.append(vsr / 365.25)
        n = len(xs)
        if n < 3:
            print(f'{site}: only {n} usable sample(s) after filtering - cannot fit. '
                  'Check that pair_results.csv has vsr_per_year_medfilt and that '
                  'cats2008_apres_window.csv covers the record.')
            continue
        mx = sum(xs)/n; my = sum(ys)/n
        sxx = sum((x-mx)**2 for x in xs)
        sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
        if sxx == 0:
            print(f'{site}: tide rate is constant over the {n} samples - cannot fit.')
            continue
        A = sxy / sxx
        resid = [(y-my) - A*(x-mx) for x, y in zip(xs, ys)]
        se = math.sqrt(sum(q*q for q in resid) / (n-2) / sxx)
        r2 = 1 - sum(q*q for q in resid) / sum((y-my)**2 for y in ys)
        print(f'{site}: n={n:3d}  A = {1e6*A:+8.1f} +/- {1e6*se:.1f} ue per m of tide '
              f'({1e3*A*100:+.2f} mm over 100 m), R2={r2:.2f}')

if __name__ == '__main__':
    main()
