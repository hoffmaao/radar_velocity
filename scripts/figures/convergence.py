#!/usr/bin/env python3
"""The convergence figure: three independent estimates of the tidal response.

Two panels, no in-plot text - axis labels only, with identity and context
carried by the caption (by request):

  (a) ALL ApRES sites' half-hour strain rates against the CATS2008 tide
      rate: GA04 orange, GA01 aqua, GA05 yellow, GA10 magenta. The fitted
      line is GA04 only - the one site whose bed pick holds - and its slope
      is the rate-method number of record (-1.24 +/- 0.04 mm per metre of
      tide over the top 100 m, R2 = 0.75). GA01's visibly distinct cloud is
      the broken-bed-pick site behaving differently, which is why it is not
      pooled. One GA01 outlier (~430 ue/day) sits above the axis clip.

  (b) every estimate of the tidal response on one axis: the flexure model
      band (curvature goes as 1/L^2 with L unmeasured), ApRES GA04, the
      four calibrated radar lines from the network inversion, the pooled
      radar value, EAGER_2022 (hollow: uncalibrated and non-reproducible,
      excluded), and the retracted chained-analysis value as its own
      labelled row so the history stays visible.

Radar numbers are cited constants from diagnostics/line_means.m (network
products, 2026-08-17); the ApRES regressions are recomputed live from the
local EAGER_ApRES repo. Built locally in Python - the ApRES data exists
only on this machine.

Run from the repo root:  python3 scripts/figures/convergence.py
Writes figs/EAGER_2022_convergence.png
"""
import csv
import math
import os
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# palette: validated project palette (see figure scripts / dataviz notes)
BLUE = "#2a78d6"; ORANGE = "#eb6834"; AQUA = "#1baf7a"
YELLOW = "#eda100"; MAGENTA = "#e87ba4"
INK = "#333333"; INK_SOFT = "#737373"; BAND = "#e6e6e2"

# ---- radar numbers, cited from line_means.m (network, 500 m blocks) ----
RADAR_LINES = [                      # (label, A [mm/m], sigma, colour, marker)
    ("GL1", -1.21, 0.58, ORANGE,  "s"),
    ("GL2", -1.01, 0.85, AQUA,    "^"),
    ("GL3", -2.40, 0.78, YELLOW,  "D"),
    ("GL4", -1.09, 1.58, MAGENTA, "v"),
]
RADAR_POOLED = (-1.46, 0.55)         # sigma includes the 1.4x build excess
EAGER = (-4.72, 0.58)                # uncalibrated + non-reproducible
MODEL_MID, MODEL_LO, MODEL_HI = -1.23, -2.2, -0.55   # L = 1.5..3 km
RETRACTED_CHAIN = 3.79               # apres_comparison.py chaining artefact


def apres_regression(site="GA04"):
    """One site's vsr vs tide rate; points and (for GA04) the fit, as in
    diagnostics/apres_rate_check.py."""
    tt, th = [], []
    with open(os.path.join(ROOT, "scripts/diagnostics/cats2008_apres_window.csv")) as f:
        for r in csv.DictReader(f):
            tt.append(float(r["datenum"])); th.append(float(r["tide_m"]))

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
    # every site with pair products; GA04 is the site of record (stable bed
    # pick), the others are context. Identity is carried by colour and told
    # in the caption - the figure itself carries no in-plot text beyond the
    # axis labels, by request.
    SITES = [("GA04", ORANGE, 0.5, 16),
             ("GA01", AQUA,   0.3, 12),
             ("GA05", YELLOW, 0.6, 20),
             ("GA10", MAGENTA, 0.6, 20)]
    xs, ys, A, se, r2, mx, my = apres_regression("GA04")
    A_mm = 1e3 * A * 100
    se_mm = 1e3 * se * 100

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.2, 5.4), dpi=150,
        gridspec_kw={"width_ratios": [1.0, 1.05], "wspace": 0.60})
    fig.patch.set_facecolor("white")

    # ---- (a) all ApRES sites, rate vs tide rate ----------------------
    ax1.axhline(0, color="#bfbfbf", lw=1)
    ax1.axvline(0, color="#bfbfbf", lw=1)
    for site, col, alpha, size in SITES:
        try:
            sx, sy, *_ = apres_regression(site)
        except FileNotFoundError:
            continue
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
    ax1.tick_params(colors=INK_SOFT)
    ax1.grid(alpha=0.18)

    # ---- (b) every estimate on one axis ------------------------------
    rows = []          # (y, label, value, err, colour, marker, bold, hollow)
    y = 0
    rows.append((y, "flexure model (L = 1.5$-$3 km)", MODEL_MID, None,
                 INK, None, False, False)); y += 1
    rows.append((y, "ApRES GA04, rate method", A_mm, se_mm,
                 ORANGE, "o", True, False)); y += 1.4
    for lbl, v, e, c, m in RADAR_LINES:
        rows.append((y, f"radar {lbl}", v, e, c, m, False, False)); y += 1
    rows.append((y, "radar pooled (GL1$-$GL4)", RADAR_POOLED[0], RADAR_POOLED[1],
                 BLUE, "o", True, False)); y += 1.4
    rows.append((y, "EAGER_2022 (uncalibrated, excluded)", EAGER[0], EAGER[1],
                 INK_SOFT, "o", False, True)); y += 1.2
    rows.append((y, "chained ApRES analysis (retracted)", RETRACTED_CHAIN, None,
                 INK_SOFT, "x", False, False)); y += 1

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
    ax2.set_yticklabels([r[1] for r in rows], fontsize=9.5, color=INK)
    ax2.invert_yaxis()
    ax2.set_ylim(rows[-1][0] + 0.9, -0.8)
    ax2.set_xlim(-5.6, 5.0)
    ax2.set_xlabel("column thickness change per metre of tide,\ntop 100 m (mm)",
                   color=INK)
    for sp in ("top", "right", "left"):
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(colors=INK_SOFT, left=False)
    ax2.grid(axis="x", alpha=0.18)

    out = os.path.join(ROOT, "figs", "EAGER_2022_convergence.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")
    print(f"ApRES GA04: {A_mm:+.2f} +/- {se_mm:.2f} mm/m (R2={r2:.2f}) | "
          f"radar pooled {RADAR_POOLED[0]:+.2f} +/- {RADAR_POOLED[1]:.2f} | model {MODEL_MID:+.2f}")


if __name__ == "__main__":
    main()
