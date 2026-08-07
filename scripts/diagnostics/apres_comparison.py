#!/usr/bin/env python3
"""Compare the repeat-pass InSAR tidal response against ApRES at the same site.

WHAT THIS DOES, AND WHY IT IS NOT A CONVERSION. The ApRES processing in
~/projects/EAGER_ApRES already fits every pair: `dh(z) = m0 + m1*(z - z_mid)`
by weighted least squares, saved per pair in `results/<site>/pairs/*.npz`
along with the underlying displacement profiles `fine_range_m`/`fine_dh_m`.
An earlier version of this script took only the scalar `vsr_per_year` out of
those fits and converted it with an assumed tidal amplitude. That threw away
the depth information the fits are built on.

This version uses the profiles directly. The pairs chain end to end with no
gaps (t2 of one is t1 of the next), so cumulatively summing `fine_dh_m` gives
displacement against time at EVERY depth on the 24-224 m grid. Each depth is
then fitted with

    dh(z,t) = a + b*t + c*tide(t)

which is exactly `vdef.fitTideAdmittance`, the estimator this project uses -
same model, same reference-invariance, same meaning for c. So the two
instruments are compared with one estimator rather than through a chain of
assumptions, and the result is a PROFILE, not a single number.

The tide is CATS2008 at the survey centroid. The ApRES site coordinates were
never recorded (GPS was off for the whole deployment), but CATS2008 varies by
0.002 m across the 5 km survey array, so any nearby point is equivalent.

WHAT IT SHOWS. ApRES and the InSAR measure the same physical quantity, and
where they overlap in depth the ApRES amplitude is far below the InSAR
systematic floor - which is why the InSAR cannot see it.

CAVEATS
  - Only GA04 is trustworthy. GA01's tracked bed drifts -81.6 m over the
    record and GA05's +158.7 m, i.e. the bed pick jumps between reflectors.
  - The ApRES grid starts at 24 m, so it does NOT constrain the top 24 m,
    where firn compaction is fastest.
  - Cumulative summing propagates any per-pair bias into a drift; the linear
    term b absorbs a constant drift but not a varying one.
  - 5.8 days cannot separate K1 from O1 or M2 from S2. The joint fit against
    the full predicted tide sidesteps constituent splitting, but a phase
    error in CATS2008 would bias c.
"""

import csv
import glob
import math
import os
from datetime import datetime, timezone

import numpy as np

APRES_ROOT = os.path.expanduser("~/projects/EAGER_ApRES")
SITES = ["GA04"]                       # GA01/GA05 bed picks are unreliable
TIDE_CSV = os.path.join(os.path.dirname(__file__), "cats2008_apres_window.csv")

# InSAR side, from scripts/figures/strain_rates.m over the 0-100 m column
INSAR_FLOOR_MM_PER_M = 8.28            # same-leg systematic floor
INSAR_EXPECTED_MM_PER_M = -1.23        # thin-plate flexure estimate


def load_tide(path):
    t, h = [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            t.append(float(r["datenum"]))
            h.append(float(r["tide_m"]))
    return np.array(t), np.array(h)


def to_datenum(iso):
    d = datetime.fromisoformat(iso.replace(" ", "T"))
    if d.tzinfo is None:
        d = d.replace(tzinfo=timezone.utc)
    epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
    return 719529.0 + (d - epoch).total_seconds() / 86400.0


def chain_site(site):
    """Cumulative displacement dh(z, t) from the chained per-pair fits."""
    files = sorted(glob.glob(os.path.join(APRES_ROOT, "results", site, "pairs", "*.npz")))
    recs = []
    for f in files:
        d = np.load(f, allow_pickle=True)
        recs.append(dict(t1=str(d["t1"]), t2=str(d["t2"]),
                         z=np.asarray(d["fine_range_m"], float),
                         dh=np.asarray(d["fine_dh_m"], float)))
    recs.sort(key=lambda r: r["t1"])
    z0 = recs[0]["z"]                      # common grid; grids differ by <0.03 m
    t = [to_datenum(recs[0]["t1"])]
    cum = np.zeros_like(z0)
    D = [cum.copy()]
    for r in recs:
        dh = np.interp(z0, r["z"], r["dh"], left=np.nan, right=np.nan)
        cum = cum + dh
        D.append(cum.copy())
        t.append(to_datenum(r["t2"]))
    return z0, np.array(t), np.array(D)     # D is (ntime, ndepth)


def fit_admittance(t, y, tide):
    """dh = a + b*t + c*tide, returning c and its 1-sigma. Same model as
    vdef.fitTideAdmittance: the intercept absorbs the reference constants so
    b and c are invariant to which epoch is called zero."""
    ok = np.isfinite(y) & np.isfinite(tide) & np.isfinite(t)
    if ok.sum() < 6:
        return np.nan, np.nan
    tc = t[ok] - t[ok].mean()
    hc = tide[ok] - tide[ok].mean()
    if np.std(tc) == 0 or np.std(hc) == 0:
        return np.nan, np.nan
    if abs(np.corrcoef(tc, hc)[0, 1]) > 0.99:
        return np.nan, np.nan
    X = np.column_stack([np.ones(ok.sum()), tc, hc])
    beta, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    resid = y[ok] - X @ beta
    dof = max(1, ok.sum() - 3)
    cov = (resid @ resid / dof) * np.linalg.inv(X.T @ X)
    return beta[2], math.sqrt(abs(cov[2, 2]))


def main():
    tt, th = load_tide(TIDE_CSV)
    print("ApRES tidal response, fitted with the InSAR estimator")
    print(f"  tide: CATS2008, range {th.max()-th.min():.3f} m over the window\n")

    for site in SITES:
        z, t, D = chain_site(site)
        tide = np.interp(t, tt, th)
        adm = np.full(z.size, np.nan)
        err = np.full(z.size, np.nan)
        for j in range(z.size):
            adm[j], err[j] = fit_admittance(t, D[:, j], tide)
        mm = adm * 1e3                      # m per m of tide -> mm per m of tide
        mme = err * 1e3

        print(f"=== {site}: {D.shape[0]} epochs, {z[0]:.0f}-{z[-1]:.0f} m ===")
        print(f"{'depth[m]':>9} {'mm per m of tide':>18} {'1-sigma':>9}")
        for j in range(0, z.size, 10):
            print(f"{z[j]:9.1f} {mm[j]:18.3f} {mme[j]:9.3f}")

        # value at the depth the InSAR figure reports
        j100 = int(np.argmin(np.abs(z - 100.0)))
        print(f"\n  at {z[j100]:.0f} m: {mm[j100]:+.3f} +/- {mme[j100]:.3f} mm per metre of tide")
        print(f"  InSAR expected (flexure model):  {INSAR_EXPECTED_MM_PER_M:+.2f} mm")
        print(f"  InSAR systematic floor:          {INSAR_FLOOR_MM_PER_M:.2f} mm")
        if np.isfinite(mm[j100]) and mm[j100] != 0:
            print(f"  floor / |ApRES|:                 {INSAR_FLOOR_MM_PER_M/abs(mm[j100]):.1f}x")
        print(f"  ApRES 1-sigma / InSAR floor:     {mme[j100]/INSAR_FLOOR_MM_PER_M:.4f}"
              "   (how much finer ApRES resolves it)\n")

        # Export for scripts/figures/strain_rates.m, which draws this as the
        # reference curve instead of a thin-plate model. A measurement from a
        # second instrument beats a model with a guessed flexure wavelength.
        out = os.path.join(os.path.dirname(__file__), f"apres_{site}_tide_profile.csv")
        with open(out, "w") as f:
            f.write("depth_m,mm_per_m_tide,mm_per_m_tide_err\n")
            for j in range(z.size):
                if not np.isfinite(mm[j]):
                    continue
                f.write(f"{z[j]:.3f},{mm[j]:.5f},{mme[j]:.5f}\n")
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
