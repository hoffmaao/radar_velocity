#!/usr/bin/env python3
"""Radar against ApRES at GA04, resolved in depth - and what actually differs.

SUPERSEDES THE SCALAR COMPARISON in convergence.py, which set the ApRES
number against the radar number at the site and reported a ~4 sigma
tension. That framing could not have been right: the two were never the
same quantity. The radar admittance is a 0-100 m column mean; the ApRES
scalar is the slope of its own fit interval, which fit_used_mask shows is
101-224 m on every pair. Comparing them is comparing different depths of a
column whose strain varies with depth - which is the entire signal.

THREE CORRECTIONS THIS FIGURE APPLIES.

 1. DEPTH-RESOLVED, NOT SCALAR. The rate method is run at each depth from
    fine_dh_m directly - dh(z)/dt = A(z) d(tide)/dt - instead of on the
    project's fitted scalar. Every pair contributes independently, so
    there is no chaining and none of the cumulative drift that produced
    the retracted +3.79.

 2. ONE DEPTH AXIS. ApRES converts range with a CONSTANT permittivity of
    3.18 and no firn correction (scripts/apres/config.py); the radar uses
    vdef.firnColumn. The ApRES axis is therefore shallow by 7 m at 30 m
    and about 14 m below the firn, so its "101-224 m" is really 115-238 m
    in the radar's frame. Every depth-matched comparison before this one
    was misaligned by that much. ApRES is mapped onto the radar's scale
    here before anything is compared.

 3. A SMOOTH MODEL WITH AN HONEST ERROR. Interval "strain" was a straight
    line fitted to dh over a hand-chosen window. It has no physics in it -
    bending makes strain LINEAR in depth, so dh is quadratic - and it is
    unstable: the radar gave +2.58 over 30-100 m and -1.44 over 45-100 m,
    a swing of 3.6x its own formal sigma. Its error bar was also
    optimistic, because adjacent depths are correlated by the fast-time
    multilook and the range resolution, exactly as vdef.invertStrainRate
    already allows for with an effective sample count. Here dh(z) is
    fitted with the same LEGENDRE basis, strain is the analytic
    derivative, and the covariance is inflated by the MEASURED residual
    correlation length. Panel (b)'s band is what that buys.

WHAT THE FIGURE SHOWS. The two strain profiles are PARALLEL. Decomposed
in panel (c), the bending curvature - the gradient of strain with depth,
which is what the tide actually drives - agrees between the instruments,
while a depth-uniform offset of about 4 mm per 100 m does not. So the
disagreement is not depth-dependent processing; it is a bulk term. That
is the shape of the mis-registration leak documented in the README:
a shift d leaves d*(dtau'(t) - dtau'(t_ref)), which is linear in depth,
hence a constant strain offset, and tide-proportional, hence in the
admittance.

DATA. The radar profile is committed at
data/EAGER_2022_GA04_radar_depth_profile.csv (built on the server by
pooling the four legs within 300 m of GA04, network inversion and the
reference-invariant fit at each depth). The ApRES side is computed here
from the pair archive, which is local-only, the same arrangement
convergence.py uses.

Run from the repo root:
    python3 scripts/figures/apres_radar_depth.py
"""
import csv, glob, math, os, sys
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grl_style

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RADAR_CSV = os.path.join(ROOT, "data", "EAGER_2022_GA04_radar_depth_profile.csv")
TIDE_CSV = os.path.join(ROOT, "scripts", "diagnostics", "cats2008_apres_sites.csv")
APRES = os.path.expanduser("~/projects/EAGER_ApRES/results")
OUT = os.path.join(ROOT, "figs", "EAGER_2022_apres_radar_depth.png")

SITE = "GA04"
DT_LO, DT_HI = 0.01, 0.1      # d; the rate method needs SHORT pairs
ZLO, ZHI, ZSTEP = 30.0, 220.0, 2.0
FIT_LO, FIT_HI = 45.0, 210.0  # where both instruments have real coverage
ORDER = 2                     # bending: strain linear in depth
Z_REF = 130.0                 # where the depth-uniform level is quoted

RAD, APR, INK, SOFT = "#2A78D6", "#EB6834", "#000000", "#6E6E6E"


# --------------------------------------------------------------------------
def firn_to_true():
    """Map an ApRES constant-permittivity range onto true depth.

    ApRES assumes er = 3.18 everywhere. The radar uses a Herron-Langway
    firn density with the Kovacs mixing law (vdef.firnColumn /
    vdef.defaultParams), where the wave is faster in the firn, so a given
    traveltime corresponds to a GREATER depth than the constant-er
    conversion returns.
    """
    rho_sfc, rho_bco, bco = 0.35, 0.81, 60.0
    kov_a, rho_ice, c = 0.845e-3, 917.0, 299792458.0
    d = np.linspace(0, 400, 40001)
    beta = np.log(1 / rho_sfc - 1)
    slope = np.log(1 / rho_bco - 1) - beta
    rho = np.minimum(1 / (1 + np.exp(beta + slope * d / bco)), 1.0)
    n = 1 + kov_a * rho * rho_ice
    twtt = 2 * np.cumsum(np.gradient(d) * n) / c
    z_const = c * twtt / (2 * np.sqrt(3.18))
    return lambda zc: np.interp(zc, z_const, d)


def dnum(dt):
    ep = datetime(1970, 1, 1, tzinfo=timezone.utc)
    return 719529.0 + (dt - ep).total_seconds() / 86400.0


def apres_profile():
    """ApRES admittance at each depth, by the rate method. No chaining."""
    if not os.path.exists(TIDE_CSV):
        sys.exit(f"missing {TIDE_CSV}")
    t, h = [], []
    rows = [l for l in open(TIDE_CSV) if not l.startswith("#")]
    for r in csv.DictReader(rows):
        t.append(float(r["datenum"])); h.append(float(r[SITE]))
    t, h = np.array(t), np.array(h)

    def tide_rate(x):
        i = int(np.searchsorted(t, x))
        if i < 2 or i > len(t) - 3:
            return None
        return (h[i + 1] - h[i - 1]) / (t[i + 1] - t[i - 1])

    Z = np.arange(ZLO, ZHI + 1e-9, ZSTEP)
    fs = sorted(glob.glob(os.path.join(APRES, SITE, "pairs", "*.npz")))
    if not fs:
        sys.exit(f"no ApRES pairs under {APRES}/{SITE} - this side is local-only")
    DH, TR = [], []
    for f in fs:
        d = np.load(f)
        dt_days = float(d["dt_days"])
        if not (DT_LO < dt_days < DT_HI):
            continue
        t1 = datetime.fromisoformat(str(d["t1"]).replace(" ", "T"))
        if t1.tzinfo is None:
            t1 = t1.replace(tzinfo=timezone.utc)
        tr = tide_rate(dnum(t1) + dt_days / 2)
        if tr is None:
            continue
        # EVERY PAIR HAS ITS OWN DEPTH AXIS - all 101 points, but the grids
        # shift by up to 0.96 m against a 1.1-1.8 m step, so stacking by
        # index mixes depths that differ by most of a sample.
        z = np.asarray(d["fine_range_m"], float)
        y = np.asarray(d["fine_dh_m"], float)
        ok = np.isfinite(z) & np.isfinite(y)
        if ok.sum() < 20:
            continue
        o = np.argsort(z[ok])
        DH.append(np.interp(Z, z[ok][o], y[ok][o], left=np.nan, right=np.nan) / dt_days)
        TR.append(tr)
    DH, TR = np.array(DH), np.array(TR)
    A = np.full(Z.size, np.nan); S = np.full(Z.size, np.nan)
    for k in range(Z.size):
        y = DH[:, k]
        m = np.isfinite(y)
        if m.sum() < 10:
            continue
        x = TR[m] - TR[m].mean(); yy = y[m] - y[m].mean()
        sxx = (x ** 2).sum()
        if sxx == 0:
            continue
        a = (x * yy).sum() / sxx
        res = yy - a * x
        A[k] = a
        S[k] = math.sqrt((res ** 2).sum() / (m.sum() - 2) / sxx)
    return firn_to_true()(Z), 1e3 * A, 1e3 * S, DH.shape[0]


def radar_profile():
    rows = list(csv.DictReader(
        [l for l in open(RADAR_CSV) if not l.startswith("#")]))
    z = np.array([float(r["depth_m"]) for r in rows])
    v = np.array([float(r["pooled"]) for r in rows])
    e = np.array([float(r["formal"]) for r in rows])
    return z, v, e


# --------------------------------------------------------------------------
def legendre_fit(z, y, sig, K=ORDER):
    """Weighted Legendre fit of dh(z); strain is the analytic derivative.

    The covariance is scaled by the residual correlation length, measured
    as the lag at which the residual autocorrelation first falls below
    1/e. That is what turns a formally tiny error into an honest one.
    """
    ok = np.isfinite(z) & np.isfinite(y) & np.isfinite(sig) & (sig > 0)
    z, y, sig = z[ok], y[ok], sig[ok]
    z0, z1 = z.min(), z.max()
    x = 2 * (z - z0) / (z1 - z0) - 1.0
    P = np.zeros((len(x), K + 1)); dP = np.zeros((len(x), K + 1))
    P[:, 0] = 1.0
    if K >= 1:
        P[:, 1] = x; dP[:, 1] = 1.0
    for k in range(1, K):
        P[:, k + 1] = ((2 * k + 1) * x * P[:, k] - k * P[:, k - 1]) / (k + 1)
        dP[:, k + 1] = (2 * k + 1) * P[:, k] + dP[:, k - 1]
    W = np.diag(1 / sig ** 2)
    Ainv = np.linalg.inv(P.T @ W @ P)
    coef = Ainv @ (P.T @ W @ y)
    resid = y - P @ coef
    chi2 = float((resid ** 2 / sig ** 2).sum() / max(len(x) - (K + 1), 1))
    r = resid - resid.mean()
    ac = np.correlate(r, r, "full")[len(r) - 1:]
    ac = ac / ac[0] if ac[0] else ac
    kk = int(np.argmax(ac < 1 / np.e))
    L = float(max(kk, 1))
    cov = Ainv * max(chi2, 1.0) * L
    dxdz = 2.0 / (z1 - z0)
    strain = (dP @ coef) * dxdz * 100.0
    var = np.einsum("ij,jk,ik->i", dP, cov, dP) * dxdz ** 2 * 100.0 ** 2
    return dict(z=z, fit=P @ coef, strain=strain, sd=np.sqrt(var),
                chi2=chi2, L=L)


def main():
    grl_style.apply()
    zr, mr, er = radar_profile()
    za, ma, ea, npair = apres_profile()

    F = {}
    for nm, (z, m, e) in (("radar", (zr, mr, er)), ("ApRES", (za, ma, ea))):
        s = (z >= FIT_LO) & (z <= FIT_HI)
        F[nm] = legendre_fit(z[s], m[s], e[s])

    # the decomposition: curvature (gradient) vs depth-uniform level
    dec = {}
    for nm, f in F.items():
        g = np.polyfit(f["z"], f["strain"], 1)[0] * 100.0
        ge = math.hypot(f["sd"][0], f["sd"][-1]) / (f["z"][-1] - f["z"][0]) * 100.0
        i = int(np.argmin(abs(f["z"] - Z_REF)))
        dec[nm] = (g, ge, f["strain"][i], f["sd"][i])

    fig, ax = plt.subplots(
        1, 3, figsize=grl_style.figsize(grl_style.TWO_COL_MM, 0.42),
        width_ratios=[1.0, 1.0, 0.85], constrained_layout=True)

    # (a) measurement and model
    for nm, c, z, m, e in (("radar", RAD, zr, mr, er), ("ApRES", APR, za, ma, ea)):
        s = (z >= FIT_LO) & (z <= FIT_HI)
        ax[0].fill_betweenx(z[s], (m - e)[s], (m + e)[s], color=c, alpha=0.15, lw=0)
        ax[0].plot(m[s], z[s], "-", color=c, lw=0.7, alpha=0.5)
        ax[0].plot(F[nm]["fit"], F[nm]["z"], "-", color=c, lw=1.8,
                   label=f'{nm}  $\\chi^2$/dof {F[nm]["chi2"]:.1f}')
    ax[0].axvline(0, color="0.85", lw=0.7)
    ax[0].set_xlabel("Displacement per m of tide (mm)")
    ax[0].set_ylabel("Depth below surface (m)")
    ax[0].legend(loc="lower left", frameon=False)

    # (b) strain, the quantity both instruments can be asked for
    for nm, c in (("radar", RAD), ("ApRES", APR)):
        f = F[nm]
        ax[1].fill_betweenx(f["z"], f["strain"] - f["sd"], f["strain"] + f["sd"],
                            color=c, alpha=0.25, lw=0)
        ax[1].plot(f["strain"], f["z"], "-", color=c, lw=1.8, label=nm)
    ax[1].axvline(0, color="0.85", lw=0.7)
    ax[1].set_xlabel("Column strain (mm per 100 m per m tide)")
    ax[1].set_ylabel("Depth below surface (m)")
    ax[1].legend(loc="lower left", frameon=False)

    for a in ax[:2]:
        a.set_ylim(FIT_HI, FIT_LO)
        a.grid(alpha=0.15)

    # (c) what agrees and what does not
    items = [("bending curvature\n(strain gradient)", 0, "per 100 m"),
             (f"depth-uniform level\n(strain at {Z_REF:.0f} m)", 2, "")]
    y = np.array([1.0, 0.0])
    h = 0.26
    for k, (lab, idx, _u) in enumerate(items):
        for j, (nm, c) in enumerate((("radar", RAD), ("ApRES", APR))):
            v, e = dec[nm][idx], dec[nm][idx + 1]
            yy = y[k] + (h / 1.6 if j == 0 else -h / 1.6)
            ax[2].barh(yy, v, height=h, color=c, alpha=0.85)
            ax[2].errorbar(v, yy, xerr=e, color=INK, lw=1.0, capsize=2)
        vr, ere, va, ea_ = dec["radar"][idx], dec["radar"][idx + 1], \
            dec["ApRES"][idx], dec["ApRES"][idx + 1]
        ns = abs(vr - va) / math.hypot(ere, ea_)
        verdict = "consistent" if ns < 3 else "not consistent"
        ax[2].text(0.5, y[k] + 0.42, f"{ns:.1f}$\\sigma$ - {verdict}",
                   transform=ax[2].get_yaxis_transform(), ha="center",
                   fontsize=grl_style.SMALL_PT, color=SOFT)
    ax[2].axvline(0, color=INK, lw=0.8)
    ax[2].set_ylim(-0.55, 1.62)
    ax[2].set_yticks(y)
    ax[2].set_yticklabels([i[0] for i in items], fontsize=grl_style.SMALL_PT)
    ax[2].set_xlabel("mm per 100 m per m tide")
    ax[2].grid(axis="x", alpha=0.15)
    ax[2].set_axisbelow(True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, facecolor="white")
    print(f"wrote {OUT}")
    print(f"ApRES {SITE}: {npair} pairs, axis mapped {ZLO:.0f}-{ZHI:.0f} m "
          f"(reported) -> {za[0]:.0f}-{za[-1]:.0f} m (true)")
    for nm, f in F.items():
        print(f"  {nm:6s} chi2/dof {f['chi2']:.2f}, residual correlation "
              f"length {f['L']:.0f} samples")
    for lab, idx, unit in items:
        vr, ere = dec["radar"][idx], dec["radar"][idx + 1]
        va, ea_ = dec["ApRES"][idx], dec["ApRES"][idx + 1]
        ns = abs(vr - va) / math.hypot(ere, ea_)
        print(f"  {lab.replace(chr(10), ' '):44s} radar {vr:+6.2f}+/-{ere:4.2f}  "
              f"ApRES {va:+6.2f}+/-{ea_:4.2f}  {ns:4.1f} sigma")


if __name__ == "__main__":
    main()
