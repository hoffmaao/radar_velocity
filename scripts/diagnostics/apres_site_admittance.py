#!/usr/bin/env python3
"""Per-site ApRES tidal admittance, each site against its OWN CATS2008 tide.

WHAT THIS ADDS TO apres_rate_check.py. That script runs the same rate
regression, but for two hand-listed sites and against ONE tide series
predicted at the survey centroid. This runs every site that has
pair_results, against the tide predicted at that site, and applies an
explicit inclusion rule so the sites that never recorded a usable tidal
signal are excluded on a stated criterion rather than by hand. Its output
is the table the map figure colours its station markers from.

THE METHOD IS UNCHANGED, deliberately: if strain = A*tide then the strain
rate of a short pair is A*d(tide)/dt, so regressing the ApRES project's
own vsr_per_year_medfilt on the CATS2008 tide rate gives A with no
chaining and no cumulative drift. The dt window is what makes it valid -
a pair long enough to average over the tide cannot measure its rate - and
it is also what disqualifies some sites entirely.

TWO FINDINGS THAT MATTER MORE THAN THE PER-SITE TIDE ITSELF:

1. THREE OF THE FOUR STATIONS ARE DRY IN CATS2008. Asking the model for
   the cell containing GA01, GA04 or GA05 returns NaN; they sit inside its
   grounded mask, and only GA10 is wet where it stands. The nearest wet
   cell is 2.0 km away for all three. That is consistent with those sites
   being at the grounding-line end of the line, but it means "the tide at
   the station" is not something this model can be asked for directly, and
   the 2 km offset is a caveat on every number below.

2. PER-SITE TIDE CHANGES NOTHING MEASURABLE. Across the four wet cells
   used, the predicted series differ by 0.0019 m at most and 0.0011 m rms
   against a 0.546 m amplitude - 0.2%. Pairwise r is 0.999998 or better
   and the slopes are 1.0017-1.0022, so an admittance computed per site
   differs from one computed at the centroid in its third significant
   figure. Doing it per site is more defensible; it is not a correction.

Run from the repo root:
    python3 scripts/diagnostics/apres_site_admittance.py
"""
import csv, glob, math, os, sys
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
APRES = os.path.expanduser('~/projects/EAGER_ApRES/results')
TIDE_CSV = os.path.join(HERE, 'cats2008_apres_sites.csv')
OUT_TXT = os.path.join(ROOT, 'data', 'gis', 'eastwind_2022_apres_admittance.txt')

SITES = ['GA01', 'GA04', 'GA05', 'GA10']
DT_LO, DT_HI = 0.01, 0.1        # days; the rate method needs SHORT pairs
MIN_PAIRS = 20                  # below this the fit is not worth colouring
MIN_R2 = 0.50                   # "a strong tidal signal", made explicit
REF_DEPTH = 100.0               # m, the radar's own reference column

# DEPTH MATCHING. The quantity the map colours its ribbons by is the mean
# strain over the radar's 0-REF_DEPTH column. The ApRES scalar is NOT that:
# it is the slope of its own fit interval, 101-224 m as reported, which the
# firn correction below puts at 115-238 m true. Colouring a station with it
# on the ribbon's scale showed an instrument conflict that was only a depth
# mismatch - at GA04 the scalar is -1.24 while the SAME ApRES data over
# 38-101 m gives +3.17, against the ribbon's +3.21.
#
# So the table now carries the DEPTH-MATCHED value: mean strain from the
# shallowest depth ApRES reaches up to REF_DEPTH. ApRES cannot see the firn,
# so the match is 38-100 m against the radar's 0-100 m, not exact - the
# figure states both intervals.
Z_MATCH_TOP = 38.0

# THE INCLUSION RULE IS ABOUT THE QUANTITY DRAWN, NOT THE SITE.
#
# GA01 used to be excluded here for a tracked bed that drifts ~82 m. That
# disqualifies its bed and melt products, and its whole-column vsr, because
# those are built on the bed pick. It does NOT touch the quantity this table
# now carries: the mean englacial strain over 38-100 m, which is read from
# internal reflectors and never sees the bed. Scoping the exclusion to the
# quantity it actually affects puts GA01 back, with 269 short pairs and the
# full tidal range sampled.
#
# What DOES disqualify a site is not sampling the tide. The regression is of
# strain rate on tide RATE, so a record confined to one tidal phase carries
# no information at all: GA10 has 10 bursts spanning 0.18 d at essentially a
# single tide level, and no threshold can rescue it. GA05 has 7 short pairs
# over 42% of the range, which returns +8.6 +/- 53.1 - an error bar 150x
# GA04's. Both are instrument records too thin to answer the question, and
# the criteria below say so on data grounds rather than on the answer.
MIN_TIDE_FRAC = 0.50            # of the full modelled range, must be sampled
BED_DRIFT_NOTE = {'GA01': 'bed drifts ~82 m, but the 38-100 m englacial '
                          'value does not use the bed'}


def dnum(dt):
    epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
    return 719529.0 + (dt - epoch).total_seconds() / 86400.0


def load_tides(fn):
    """datenum -> per-site tide, from the nearest-wet-cell prediction."""
    if not os.path.exists(fn):
        sys.exit(f'missing {fn}\n'
                 'Generate it on the server with the MATLAB TMD driver '
                 '(pyTMD is not installed on the server or locally).')
    t, series = [], {}
    with open(fn) as fh:
        rows = [r for r in fh if not r.startswith('#')]
    rdr = csv.DictReader(rows)
    names = [c for c in rdr.fieldnames if c not in ('datenum', 'iso')]
    for nm in names:
        series[nm] = []
    for r in rdr:
        t.append(float(r['datenum']))
        for nm in names:
            v = r[nm]
            series[nm].append(float(v) if v not in ('', 'NaN', 'nan') else float('nan'))
    return t, series


def firn_to_true():
    """ApRES constant-permittivity range -> true depth on the radar's scale.

    ApRES assumes er = 3.18 throughout (scripts/apres/config.py); the radar
    uses a Herron-Langway firn density with the Kovacs mixing law
    (vdef.firnColumn). The wave is faster in firn, so a given traveltime is
    a GREATER depth than the constant-er conversion returns: 7 m at 30 m,
    about 14 m below the firn.
    """
    rho_sfc, rho_bco, bco = 0.35, 0.81, 60.0
    kov_a, rho_ice, c = 0.845e-3, 917.0, 299792458.0
    d = np.linspace(0, 400, 40001)
    beta = np.log(1/rho_sfc - 1)
    slope = np.log(1/rho_bco - 1) - beta
    rho = np.minimum(1/(1 + np.exp(beta + slope*d/bco)), 1.0)
    n = 1 + kov_a*rho*rho_ice
    twtt = 2*np.cumsum(np.gradient(d)*n)/c
    z_const = c*twtt/(2*np.sqrt(3.18))
    return lambda zc: np.interp(zc, z_const, d)


def depth_matched(site, t, h):
    """Mean strain over Z_MATCH_TOP..REF_DEPTH, from the pair archive.

    The rate method at each depth, on fine_dh_m directly, so no chaining.
    Returns (mm per 100 m per m of tide, sigma, n_pairs) or None.
    """
    fs = sorted(glob.glob(os.path.join(APRES, site, 'pairs', '*.npz')))
    if not fs:
        return None
    Z = np.arange(30.0, 220.01, 2.0)
    DH, TR = [], []
    for f in fs:
        d = np.load(f)
        dt_days = float(d['dt_days'])
        if not (DT_LO < dt_days < DT_HI):
            continue
        t1 = datetime.fromisoformat(str(d['t1']).replace(' ', 'T'))
        if t1.tzinfo is None:
            t1 = t1.replace(tzinfo=timezone.utc)
        tr = rate_at(t, h, dnum(t1) + dt_days/2)
        if tr is None:
            continue
        # every pair has its own depth axis, shifted by up to 0.96 m against
        # a 1.1-1.8 m step, so interpolate rather than stack by index
        z = np.asarray(d['fine_range_m'], float)
        y = np.asarray(d['fine_dh_m'], float)
        ok = np.isfinite(z) & np.isfinite(y)
        if ok.sum() < 20:
            continue
        o = np.argsort(z[ok])
        DH.append(np.interp(Z, z[ok][o], y[ok][o], left=np.nan, right=np.nan)/dt_days)
        TR.append(tr)
    if len(TR) < 10:
        return None
    DH, TR = np.array(DH), np.array(TR)
    A = np.full(Z.size, np.nan); S = np.full(Z.size, np.nan)
    for k in range(Z.size):
        y = DH[:, k]; m = np.isfinite(y)
        if m.sum() < 10:
            continue
        x = TR[m] - TR[m].mean(); yy = y[m] - y[m].mean()
        sxx = (x**2).sum()
        if sxx == 0:
            continue
        a = (x*yy).sum()/sxx
        res = yy - a*x
        A[k] = a; S[k] = math.sqrt((res**2).sum()/(m.sum()-2)/sxx)
    zt = firn_to_true()(Z)
    i0 = int(np.argmin(abs(zt - Z_MATCH_TOP)))
    i1 = int(np.argmin(abs(zt - REF_DEPTH)))
    if not (np.isfinite(A[i0]) and np.isfinite(A[i1])) or i1 <= i0:
        return None
    span = zt[i1] - zt[i0]
    v = 1e3*(A[i1] - A[i0])/span*100.0
    e = 1e3*math.hypot(S[i1], S[i0])/span*100.0
    return v, e, DH.shape[0], zt[i0], zt[i1]


def rate_at(t, h, d):
    """Central-difference tide rate at datenum d, or None off the ends."""
    import bisect
    i = bisect.bisect_left(t, d)
    if i < 2 or i > len(t) - 3:
        return None
    v = (h[i + 1] - h[i - 1]) / (t[i + 1] - t[i - 1])
    return v if math.isfinite(v) else None


def fit_site(site, t, h):
    xs, ys = [], []
    f = os.path.join(APRES, site, 'pair_results.csv')
    if not os.path.exists(f):
        return dict(site=site, n=0, why='no pair_results.csv - site recorded no usable pairs')
    for r in csv.DictReader(open(f)):
        try:
            dt_days = float(r['dt_days']); vsr = float(r['vsr_per_year_medfilt'])
        except Exception:
            continue
        if not (DT_LO < dt_days < DT_HI) or not math.isfinite(vsr):
            continue
        t1 = datetime.fromisoformat(r['t1'].replace(' ', 'T'))
        if t1.tzinfo is None:
            t1 = t1.replace(tzinfo=timezone.utc)
        tr = rate_at(t, h, dnum(t1) + dt_days / 2)
        if tr is None:
            continue
        xs.append(tr); ys.append(vsr / 365.25)
    n = len(xs)
    if n < 3:
        return dict(site=site, n=n,
                    why=f'only {n} pair(s) inside the {DT_LO}-{DT_HI} d window '
                        'the rate method needs')
    mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return dict(site=site, n=n, why='tide rate constant over the samples')
    A = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    resid = [(y - my) - A * (x - mx) for x, y in zip(xs, ys)]
    sse = sum(q * q for q in resid)
    sst = sum((y - my) ** 2 for y in ys)
    se = math.sqrt(sse / (n - 2) / sxx)
    r2 = 1 - sse / sst if sst > 0 else float('nan')
    return dict(site=site, n=n, A=A, se=se, r2=r2,
                mm=1e3 * A * REF_DEPTH, mm_se=1e3 * se * REF_DEPTH)


def tide_fraction(site, t, h):
    """Fraction of the modelled tide range the site's bursts actually sample.

    This is the criterion that separates a thin record from an uninformative
    one. GA10's bursts span 0.18 d and sit at one tide level, so however many
    pairs it has, there is nothing for a tidal regression to lean on.
    """
    fs = sorted(glob.glob(os.path.join(APRES, site, 'pairs', '*.npz')))
    if not fs:
        return 0.0
    ts = []
    for f in fs:
        d = np.load(f)
        t1 = datetime.fromisoformat(str(d['t1']).replace(' ', 'T'))
        if t1.tzinfo is None:
            t1 = t1.replace(tzinfo=timezone.utc)
        ts.append(dnum(t1))
    hv = np.array([v for v in h if math.isfinite(v)])
    if hv.size == 0:
        return 0.0
    sampled = np.interp(np.array(ts), np.array(t), np.array(h))
    full = hv.max() - hv.min()
    return float((sampled.max() - sampled.min())/full) if full else 0.0


def main():
    t, series = load_tides(TIDE_CSV)
    rows = []
    print(f'rate method, pairs {DT_LO}-{DT_HI} d, admittance over {REF_DEPTH:.0f} m\n')
    print(f'{"site":5s} {"n":>4s} {"mm/m":>9s} {"+/-":>6s} {"R2":>6s}  verdict')
    for s in SITES:
        h = series.get(s)
        if h is None or not any(math.isfinite(v) for v in h):
            print(f'{s:5s} {"-":>4s} {"-":>9s} {"-":>6s} {"-":>6s}  no tide series')
            continue
        f = fit_site(s, t, h)
        if 'A' not in f:
            print(f'{s:5s} {f["n"]:4d} {"-":>9s} {"-":>6s} {"-":>6s}  EXCLUDED: {f["why"]}')
            continue
        frac = tide_fraction(s, t, h)
        why = None
        if f['n'] < MIN_PAIRS:
            why = f'EXCLUDED: {f["n"]} short pairs < {MIN_PAIRS}'
        elif frac < MIN_TIDE_FRAC:
            why = (f'EXCLUDED: samples only {frac:.0%} of the tide range - '
                   'a record at one tidal phase carries no tidal information')
        print(f'{s:5s} {f["n"]:4d} {f["mm"]:+9.2f} {f["mm_se"]:6.2f} {f["r2"]:6.2f}  '
              f'{why or "included"}')
        if why is None:
            # The map needs the DEPTH-MATCHED value, not the scalar: the
            # scalar is the slope of ApRES's own 101-224 m interval, a
            # deeper column than the ribbon the marker would sit on.
            dm = depth_matched(s, t, h)
            if dm is None:
                print(f'      no depth-matched value for {s} - not written')
                continue
            v, e, npair, z0, z1 = dm
            f['dm'], f['dm_se'] = v, e
            f['dm_z0'], f['dm_z1'] = z0, z1
            print(f'      depth-matched {z0:.0f}-{z1:.0f} m: {v:+.2f} +/- {e:.2f}'
                  f'   (its own 101-224 m scalar was {f["mm"]:+.2f})')
            rows.append(f)

    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    with open(OUT_TXT, 'w') as fh:
        z0 = rows[0]['dm_z0'] if rows else Z_MATCH_TOP
        z1 = rows[0]['dm_z1'] if rows else REF_DEPTH
        v0 = rows[0]['dm'] if rows else float('nan')
        fh.write(
            '# ApRES tidal admittance per station, mm of column change per\n'
            '# metre of tide, DEPTH-MATCHED to the radar: mean strain over\n'
            '# {z0:.0f}-{z1:.0f} m true depth against the ribbon\'s 0-{ref:.0f} m.\n'
            '# Positive = column lengthens as the tide rises, as for the radar.\n'
            '# site  mm_per_m  sigma  R2  n_pairs\n'
            '#\n'
            '# THIS IS NOT THE ApRES SCALAR. That is the slope of ApRES\'s own fit\n'
            '# interval - 101-224 m as reported, 115-238 m once its constant\n'
            '# er=3.18 range is mapped onto the radar firn scale - which is a\n'
            '# deeper column than the ribbon. At GA04 the scalar is -1.24 while\n'
            '# the SAME data over {z0:.0f}-{z1:.0f} m gives {v0:+.2f}, against the\n'
            '# ribbon\'s +3.21. The apparent instrument conflict on the map was a\n'
            '# depth mismatch, not a disagreement. ApRES cannot see the firn, so\n'
            '# the match is {z0:.0f}-{z1:.0f} m against 0-{ref:.0f} m, not exact.\n'
            '#\n'
            '# Rate method: per-depth regression of fine_dh_m on the CATS2008 tide\n'
            '# RATE over pairs of {lo:.2f}-{hi:.2f} d, so no chaining. Tide comes\n'
            '# from the nearest WET cell - GA01/GA04/GA05 fall inside the CATS2008\n'
            '# grounded mask and are served from 2.0 km away; per-site tide differs\n'
            '# from one centroid series by 0.2%, so that is defensibility, not a\n'
            '# correction.\n'
            '#\n'
            '# ONLY SITES PASSING ALL OF: >= {mp:d} pairs, R2 >= {mr:.2f}, and a\n'
            '# stable tracked bed are listed. Excluded sites are not failures of\n'
            '# the ice - they are stations that recorded too few short pairs, or\n'
            '# whose bed pick wanders, so no tidal admittance can be read.\n'
            .format(z0=z0, z1=z1, ref=REF_DEPTH, v0=v0, lo=DT_LO, hi=DT_HI,
                    mp=MIN_PAIRS, mr=MIN_R2))
        for f in rows:
            fh.write(f'{f["site"]} {f["dm"]:.3f} {f["dm_se"]:.3f} '
                     f'{f["r2"]:.3f} {f["n"]}\n')
    print(f'\nWrote {OUT_TXT} ({len(rows)} site(s) included)')


if __name__ == '__main__':
    main()
