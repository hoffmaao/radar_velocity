#!/usr/bin/env python3
"""ApRES against the radar AT THE ApRES SITE: every estimate of the tidal response.

Two panels, no in-plot text - axis labels only, with identity and context
carried by the caption (by request):

  (a) ALL ApRES sites' half-hour strain rates against the CATS2008 tide
      rate: GA04 orange, GA01 aqua, GA05 yellow, GA10 magenta. The fitted
      line is GA04 only - the one site whose bed pick holds - and its slope
      is the rate-method number of record (-1.24 +/- 0.04 mm per metre of
      tide over the top 100 m, R2 = 0.75). GA01's visibly distinct cloud is
      the broken-bed-pick site behaving differently, which is why it is not
      pooled. One GA01 outlier (~430 ue/day) sits above the axis clip.

  (b) every estimate of the tidal response AT THE ApRES SITE, on one axis:
      the flexure model band (curvature goes as 1/L^2 with L unmeasured),
      the two ApRES GA04 estimates (rate method, and the retracted chained
      value kept so the history stays visible), the four calibrated radar
      lines' own blocks nearest that position, and their stack at two bin
      widths. The four radar lines share one blue with per-line markers -
      they are one instrument - so orange means ApRES in both panels.

WHY POSITION-MATCHED, AND WHAT CHANGED (19 Aug 2026). The ApRES site
coordinates were recovered, putting GA04 at 4.70 km along track and only
0.07 km off the line - i.e. at the GROUNDING-LINE END, not on the flat
floating section the earlier write-up assumed. This panel used to compare
ApRES against the radar LINE MEAN (-1.46 +/- 0.55), and they agreed
nicely. That agreement was an artefact of averaging: the profile runs
-3.7 to +2.5 mm/m along the line, so its mean is not the quantity a point
measurement should be compared with. Compared like for like at 4.70 km,
the radar reads +2.22 +/- 0.86 against ApRES -1.24 +/- 0.04 - opposite
signs, about 4 sigma apart. The line-mean result is not wrong and is not
discarded; it is simply a different quantity, and it still appears in the
strain-rate and error-budget figures.

HOW FAR TO TRUST THE RADAR END (diagnostics/far_end_check.m):
  - the SIGN is robust: all four lines independently give positive values
    at 4.65-4.75 km (+1.5 to +4.5), and a jackknife dropping each line in
    turn leaves +1.94 to +2.64, so no single line carries it. The
    radar-free GPS-curvature prediction also rises toward the grounding
    line, reaching about +0.8.
  - the MAGNITUDE is NOT resolved: the same stack gives +2.22 +/- 0.86 at
    500 m bins - interpolated between the 4.25 and 4.75 km bin centres
    that straddle the site - and +0.30 +/- 0.43 at 1 km, where no centre
    lies past 4.50 km so the figure is the 4.0-5.0 km bin CONTAINING the
    site rather than an interpolation. The two differ because the
    adjacent 4.25 km blocks are negative. The implied gradient (-0.34 to
    +2.51 over 500 m) is far sharper than the 1-3 km flexure scale, which
    is not physical for bending. Far-end blocks are also thinner than
    mid-line ones (valid samples -22%, fitted depth span -21%).
  - the DISAGREEMENT with ApRES survives either binning (4.0 and 3.6
    sigma), so it is not an artefact of the bin choice.
Both stack rows are drawn for exactly this reason.

Radar numbers are cited constants from diagnostics/far_end_check.m
(network products, 2026-08-19); the ApRES regressions are recomputed live
from the local EAGER_ApRES repo. Built locally in Python - the ApRES data
exists only on this machine.

Run from the repo root:  python3 scripts/figures/convergence.py
Writes figs/EAGER_2022_apres_radar_comparison.png
"""
import csv
import math
import os
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import grl_style
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# palette: validated project palette (see figure scripts / dataviz notes)
BLUE = "#2a78d6"; ORANGE = "#eb6834"; AQUA = "#1baf7a"
YELLOW = "#eda100"; MAGENTA = "#e87ba4"
# ink is pure black: axis labels, tick labels and tick marks all read black,
# matching the rest of the figure set. INK_SOFT stays for the genuinely
# recessive elements (the retracted row), which are not axis text.
INK = "#000000"; INK_SOFT = "#737373"; BAND = "#e6e6e2"

# ---- radar numbers AT THE ApRES SITE, cited from diagnostics/far_end_check.m
# GA04 projects to 4.70 km along track, 0.07 km off the line (tidal_response_map.m
# apres_along). These are the per-line blocks nearest that position - NOT line
# means, which is the whole point: the line mean of a profile running -3.7 to
# +2.5 is not the quantity a point measurement should be compared with.
APRES_ALONG_KM = 4.70
# All four radar rows share the project blue: they are four lines from ONE
# instrument and do not need four identities. That keeps orange meaning
# ApRES here as it does in panel (a). The markers stay per-line so the rows
# remain individually identifiable.
RADAR_AT_SITE = [               # (label, block km, A [mm/m], sigma, colour, marker)
    ("GL1", 4.75, +2.51, 1.33, BLUE, "s"),
    ("GL2", 4.70, +3.36, 2.53, BLUE, "^"),
    ("GL3", 4.66, +1.51, 1.62, BLUE, "D"),
    ("GL4", 4.65, +4.54, 2.89, BLUE, "v"),
]
RADAR_STACK_500 = (+2.22, 0.86)      # inverse-variance stack, 500 m bins,
                                     # interpolated between the 4.25 and
                                     # 4.75 km bin centres that straddle
                                     # the site
RADAR_STACK_1KM = (+0.30, 0.43)      # same stack at 1 km bins. There is no
                                     # centre past 4.50 km at that width, so
                                     # this is the 4.0-5.0 km bin CONTAINING
                                     # the site, not an interpolation - the
                                     # magnitude is NOT resolved; the sign is
MODEL_MID, MODEL_LO, MODEL_HI = -1.23, -2.2, -0.55   # L = 1.5..3 km
RETRACTED_CHAIN = 3.79               # apres_comparison.py chaining artefact,
                                     # also a GA04 quantity so it belongs here


_CATS = None


def cats_series():
    """The CATS2008 tide series over the ApRES window, parsed once and
    shared by every site's regression."""
    global _CATS
    if _CATS is None:
        tt, th = [], []
        fn = os.path.join(ROOT, "scripts/diagnostics/cats2008_apres_window.csv")
        with open(fn) as f:
            for r in csv.DictReader(f):
                tt.append(float(r["datenum"])); th.append(float(r["tide_m"]))
        _CATS = (tt, th)
    return _CATS


def apres_regression(site="GA04", cats=None):
    """One site's vsr vs tide rate; points and (for GA04) the fit, as in
    diagnostics/apres_rate_check.py.

    The CATS table is the same for every site, so the caller passes it in;
    the only thing this can fail to find is the site's own pair products."""
    tt, th = cats_series() if cats is None else cats

    def tide_rate(d):
        import bisect
        i = bisect.bisect_left(tt, d)
        if i < 2 or i > len(tt) - 3:
            return None
        return (th[i + 1] - th[i - 1]) / (tt[i + 1] - tt[i - 1])

    def dnum(dt):
        epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
        return 719529.0 + (dt - epoch).total_seconds() / 86400.0

    xs, ys = [], []
    src = os.path.expanduser(f"~/projects/EAGER_ApRES/results/{site}/pair_results.csv")
    with open(src) as f:
        for r in csv.DictReader(f):
            try:
                dt_days = float(r["dt_days"]); vsr = float(r["vsr_per_year_medfilt"])
            except (KeyError, ValueError):
                continue
            if not (0.01 < dt_days < 0.1) or not math.isfinite(vsr):
                continue
            t1 = datetime.fromisoformat(r["t1"].replace(" ", "T"))
            if t1.tzinfo is None:
                t1 = t1.replace(tzinfo=timezone.utc)
            tr = tide_rate(dnum(t1) + dt_days / 2)
            if tr is None:
                continue
            xs.append(tr); ys.append(vsr / 365.25)
    n = len(xs)
    if n < 3:
        # too few pairs for a fit (GA05/GA10 have handfuls); points only
        return xs, ys, float("nan"), float("nan"), float("nan"), 0.0, 0.0
    mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    A = sxy / sxx
    resid = [(y - my) - A * (x - mx) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(q * q for q in resid) / (n - 2) / sxx)
    r2 = 1 - sum(q * q for q in resid) / sum((y - my) ** 2 for y in ys)
    return xs, ys, A, se, r2, mx, my


def main():
    grl_style.apply()
    # every site with pair products; GA04 is the site of record (stable bed
    # pick), the others are context. Identity is carried by colour and told
    # in the caption - the figure itself carries no in-plot text beyond the
    # axis labels, by request.
    # GA04 is drawn LAST so the site of record (276 points) is not
    # overplotted by the context sites
    SITES = [("GA01", AQUA,   0.3, 12),
             ("GA05", YELLOW, 0.6, 20),
             ("GA10", MAGENTA, 0.6, 20),
             ("GA04", ORANGE, 0.5, 16)]
    # regress each site once and reuse; a site that contributes no points is
    # named on stdout, because with no in-plot text a silently missing
    # colour would leave the caption claiming a series that is not drawn.
    # Missing pair products and pair products that survive no filter are
    # both ways to draw nothing, so both are reported.
    cats = cats_series()
    fits = {}
    for site, *_ in SITES:
        try:
            fit = apres_regression(site, cats)
        except FileNotFoundError as e:
            print(f"warning: {site} has no pair products ({e.filename}) - "
                  f"absent from panel (a), but still named in the caption")
            continue
        if not fit[0]:
            print(f"warning: {site} has pair products but no pair survives "
                  f"the dt and tide-rate windows - absent from panel (a), "
                  f"but still named in the caption")
            continue
        fits[site] = fit
    if "GA04" not in fits:
        raise SystemExit("GA04 pair products are required for the fit of record")
    xs, ys, A, se, r2, mx, my = fits["GA04"]
    A_mm = 1e3 * A * 100
    se_mm = 1e3 * se * 100

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=grl_style.figsize(grl_style.TWO_COL_MM, 0.41),
        gridspec_kw={"width_ratios": [1.0, 1.05], "wspace": 0.60})
    fig.patch.set_facecolor("white")

    # ---- (a) all ApRES sites, rate vs tide rate ----------------------
    ax1.axhline(0, color="#bfbfbf", lw=1)
    ax1.axvline(0, color="#bfbfbf", lw=1)
    for site, col, alpha, size in SITES:
        if site not in fits:
            continue
        sx, sy, *_ = fits[site]
        ax1.scatter(sx, [1e6 * y for y in sy], s=size, color=col,
                    alpha=alpha, edgecolors="none", zorder=2)
    xf = [min(xs), max(xs)]
    ax1.plot(xf, [1e6 * (my + A * (x - mx)) for x in xf], color=INK, lw=2.2,
             zorder=3)
    # one GA01 point sits at ~430 ue/day (the broken-bed-pick site); the
    # axis holds the bulk of the data and the caption notes the clip
    ax1.set_ylim(-135, 145)
    ax1.set_xlabel("tide rate (m per day)", color=INK)
    ax1.set_ylabel("ApRES strain rate ($\\mu\\epsilon$ per day)", color=INK)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.tick_params(colors=INK)
    ax1.grid(alpha=0.18)

    # ---- (b) every estimate on one axis ------------------------------
    rows = []          # (y, label, value, err, colour, marker, bold, hollow)
    y = 0
    rows.append((y, "flexure model (L = 1.5$-$3 km)", MODEL_MID, None,
                 INK, None, False, False)); y += 1
    rows.append((y, "ApRES GA04, rate method", A_mm, se_mm,
                 ORANGE, "o", True, False)); y += 1
    rows.append((y, "ApRES GA04, chained (retracted)", RETRACTED_CHAIN, None,
                 INK_SOFT, "x", False, False)); y += 1.4
    for lbl, km, v, e, c, m in RADAR_AT_SITE:
        rows.append((y, f"radar {lbl} ({km:.2f} km)", v, e, c, m,
                     False, False)); y += 1
    rows.append((y, "radar stack, 500 m bins", RADAR_STACK_500[0],
                 RADAR_STACK_500[1], BLUE, "o", True, False)); y += 1
    rows.append((y, "radar stack, 1 km bins", RADAR_STACK_1KM[0],
                 RADAR_STACK_1KM[1], BLUE, "o", False, True)); y += 1

    ax2.axvspan(MODEL_LO, MODEL_HI, color=BAND, zorder=0)
    ax2.axvline(MODEL_MID, color=INK, lw=1.6, ls="--", zorder=1)
    ax2.axvline(0, color="#bfbfbf", lw=1, zorder=1)

    for yy, lbl, v, e, c, m, bold, hollow in rows:
        if m is None:
            continue
        if e is not None:
            ax2.plot([v - e, v + e], [yy, yy], color=c, lw=2.2 if bold else 1.4,
                     zorder=3)
        if m == "x":
            ax2.plot(v, yy, m, ms=9, color=c, mew=1.8, zorder=3)
            continue
        face = "white" if hollow else c
        ax2.plot(v, yy, m, ms=11 if bold else 8, mfc=face,
                 mec=c if hollow else "white", mew=1.4, color=c, zorder=4)

    ax2.set_yticks([r[0] for r in rows])
    ax2.set_yticklabels([r[1] for r in rows], fontsize=grl_style.BASE_PT, color=INK)
    ax2.invert_yaxis()
    ax2.set_ylim(rows[-1][0] + 0.9, -0.8)
    ax2.set_xlim(-4.2, 7.8)
    ax2.set_xlabel("column thickness change per metre of tide,\ntop 100 m (mm)",
                   color=INK)
    for sp in ("top", "right", "left"):
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(colors=INK, left=False)
    ax2.grid(axis="x", alpha=0.18)

    out = os.path.join(ROOT, "figs", "EAGER_2022_apres_radar_comparison.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")
    d500 = abs(RADAR_STACK_500[0] - A_mm) / math.hypot(RADAR_STACK_500[1], se_mm)
    d1km = abs(RADAR_STACK_1KM[0] - A_mm) / math.hypot(RADAR_STACK_1KM[1], se_mm)
    print(f"ApRES GA04: {A_mm:+.2f} +/- {se_mm:.2f} mm/m (R2={r2:.2f}) at "
          f"{APRES_ALONG_KM:.2f} km along track")
    print(f"radar at that position: {RADAR_STACK_500[0]:+.2f} +/- "
          f"{RADAR_STACK_500[1]:.2f} (500 m bins, {d500:.1f} sigma away), "
          f"{RADAR_STACK_1KM[0]:+.2f} +/- {RADAR_STACK_1KM[1]:.2f} "
          f"(1 km bins, {d1km:.1f} sigma away)")


if __name__ == "__main__":
    main()
