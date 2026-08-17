#!/usr/bin/env python3
"""The convergence figure: three independent estimates of the tidal response.

One page, two panels:

  (a) the ApRES evidence, shown as the measurement it actually is: 276
      half-hour strain rates from GA04's own pair products against the
      CATS2008 tide rate. If strain = A*tide the slope IS A - no chaining,
      no cumulative drift, nothing of ours between the instrument and the
      number.

  (b) every estimate of A on one axis, in mm of column-thickness change
      per metre of tide over the top 100 m: the thin-plate flexure model
      (drawn as a band, because its curvature goes as 1/L^2 and the
      flexure length is not measured), the ApRES value from (a), the four
      calibrated radar lines from the network inversion, and the pooled
      radar value. EAGER_2022 appears hollow - uncalibrated AND
      non-reproducible, excluded from the pool. The retracted chained
      ApRES analysis (+3.79) is marked faintly because pretending it never
      happened would misrepresent how the agreement was reached.

Radar numbers are from diagnostics/line_means.m on the network products
(2026-08-17 run) and are cited, not recomputed - this script has no server
access. The ApRES regression is recomputed live from the local repo so the
panel is evidence, not a quotation.

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


def apres_regression():
    """GA04 vsr vs tide rate; returns points and the fit, as in
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
    src = os.path.expanduser("~/projects/EAGER_ApRES/results/GA04/pair_results.csv")
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
    mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    A = sxy / sxx
    resid = [(y - my) - A * (x - mx) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(q * q for q in resid) / (n - 2) / sxx)
    r2 = 1 - sum(q * q for q in resid) / sum((y - my) ** 2 for y in ys)
    return xs, ys, A, se, r2, mx, my


def main():
    xs, ys, A, se, r2, mx, my = apres_regression()
    A_mm = 1e3 * A * 100          # strain per m -> mm per m over 100 m
    se_mm = 1e3 * se * 100

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.2, 5.4), dpi=150,
        gridspec_kw={"width_ratios": [1.0, 1.05], "wspace": 0.52})
    fig.patch.set_facecolor("white")

    # ---- (a) the ApRES evidence -------------------------------------
    ax1.axhline(0, color="#bfbfbf", lw=1)
    ax1.axvline(0, color="#bfbfbf", lw=1)
    ax1.scatter([x for x in xs], [1e6 * y for y in ys], s=14, color=ORANGE,
                alpha=0.45, edgecolors="none", zorder=2)
    xf = [min(xs), max(xs)]
    ax1.plot(xf, [1e6 * (my + A * (x - mx)) for x in xf], color=INK, lw=2.2,
             zorder=3)
    ax1.set_xlabel("tide rate (m per day)", color=INK)
    ax1.set_ylabel("ApRES strain rate ($\\mu\\epsilon$ per day)", color=INK)
    ax1.set_title("(a) ApRES GA04: their own half-hour strain rates\n"
                  f"follow the tide rate  (n={len(xs)}, $R^2$={r2:.2f})",
                  color=INK, fontsize=11)
    ax1.text(0.03, 0.04,
             f"slope: {A_mm:+.2f} $\\pm$ {se_mm:.2f} mm per metre of tide\n"
             "(over the top 100 m; rising tide $\\rightarrow$ column thins)",
             transform=ax1.transAxes, fontsize=9.5, color=INK,
             va="bottom")
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
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
                 INK_SOFT, "o", False, True)); y += 1

    ax2.axvspan(MODEL_LO, MODEL_HI, color=BAND, zorder=0)
    ax2.axvline(MODEL_MID, color=INK, lw=1.6, ls="--", zorder=1)
    ax2.axvline(0, color="#bfbfbf", lw=1, zorder=1)

    for yy, lbl, v, e, c, m, bold, hollow in rows:
        if m is None:
            continue
        if e is not None:
            ax2.plot([v - e, v + e], [yy, yy], color=c, lw=2.2 if bold else 1.4,
                     zorder=3)
        face = "white" if hollow else c
        ax2.plot(v, yy, m, ms=11 if bold else 8, mfc=face, mec=c if hollow else "white",
                 mew=1.4, color=c, zorder=4)

    # the retracted value, faint, so the history is visible
    yr = rows[-1][0] + 1.2
    ax2.plot(RETRACTED_CHAIN, yr, "x", ms=9, color=INK_SOFT, mew=1.8, zorder=3)
    ax2.text(RETRACTED_CHAIN - 0.35, yr, "chained ApRES analysis\n(retracted artefact)  ",
             fontsize=8, color=INK_SOFT, ha="right", va="center")

    ax2.set_yticks([r[0] for r in rows])
    ax2.set_yticklabels([r[1] for r in rows], fontsize=9.5, color=INK)
    ax2.invert_yaxis()
    ax2.set_ylim(yr + 1.6, -0.8)
    ax2.set_xlim(-5.6, 5.0)
    ax2.set_xlabel("column thickness change per metre of tide,\ntop 100 m (mm)",
                   color=INK)
    ax2.set_title("(b) Three independent routes, one answer", color=INK,
                  fontsize=11)
    ax2.text(MODEL_MID, rows[-1][0] + 0.55, "model", fontsize=8.5, color=INK,
             ha="center", va="top")
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)
    ax2.tick_params(colors=INK_SOFT, left=False)
    ax2.grid(axis="x", alpha=0.18)

    out = os.path.join(ROOT, "figs", "EAGER_2022_convergence.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")
    print(f"ApRES: {A_mm:+.2f} +/- {se_mm:.2f} mm/m | radar pooled "
          f"{RADAR_POOLED[0]:+.2f} +/- {RADAR_POOLED[1]:.2f} | model {MODEL_MID:+.2f}")


if __name__ == "__main__":
    main()
