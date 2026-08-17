#!/usr/bin/env python3
"""Method schematic - the concept figure for repeat-pass vertical strain.

Three panels:

  (a) THE SETTING. A ground-based radar is towed along the same ~5 km line
      thirteen times over four days, on a floating shelf that rides the
      tide and flexes against the grounding line. An ApRES sits on the
      same ice. The key physical fact the whole method rests on: the
      antenna rides the surface, so the tide itself is invisible to the
      radar - only deformation of the column below is observable.

  (b) THE OBSERVABLE. Two passes at different tide states, surface returns
      aligned (per-column coalignment): internal layers are displaced
      relative to the surface by an amount growing with depth. That
      displacement profile dh(z) IS the column strain.

  (c) THE ESTIMATION. Passes sample the tide at thirteen epochs; every
      pair is an interferogram, and the all-pairs network is inverted for
      per-epoch displacement with no reference pass. The joint fit
      x = a + b*t + c*tide separates the secular strain rate from the
      tidal response.

Run from the repo root:  python3 scripts/figures/method_schematic.py
Writes figs/EAGER_2022_method_schematic.png
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

BLUE = "#2a78d6"; ORANGE = "#eb6834"; AQUA = "#1baf7a"
YELLOW = "#eda100"; MAGENTA = "#e87ba4"
INK = "#333333"; INK_SOFT = "#737373"
ICE = "#eef3f8"; OCEAN = "#cfe3f5"; BED = "#d9cfc0"; LAYER = "#9db8d2"

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def bed(x):
    return 0.6 + max(0.0, x - 5.0) * 0.78


def flex(x, xgl=7.05, L=2.2):
    """Normalised tidal deflection of the floating surface: 1 far out,
    pinned to 0 at the grounding line."""
    if x >= xgl:
        return 0.0
    return 1.0 - math.exp(-(xgl - x) / L)


def panel_a(ax):
    xgl = 7.05
    xs = [i * 0.05 for i in range(201)]

    # bedrock
    ax.fill_between([0, 10], 0, [bed(0), bed(10)], color=BED, zorder=0)
    ax.fill_between(xs, 0, [bed(x) for x in xs], color=BED, zorder=1)
    # ocean under the floating part
    ax.fill_between([x for x in xs if x <= xgl],
                    [bed(x) for x in xs if x <= xgl], 2.2,
                    color=OCEAN, zorder=1)
    # ice body: floating base at 2.2, grounded base follows bed
    base = [max(2.2, bed(x)) for x in xs]
    ax.fill_between(xs, base, 4.0, color=ICE, zorder=2)
    ax.plot(xs, base, color=INK_SOFT, lw=1.0, zorder=3)
    ax.plot([0, 10], [4.0, 4.0], color=INK, lw=1.4, zorder=4)

    # internal layers
    for d in (0.55, 1.05, 1.55):
        ax.plot(xs, [4.0 - d for _ in xs], color=LAYER, lw=0.9, zorder=3)

    # tidal envelopes of the floating surface (exaggerated)
    xf = [x for x in xs if x <= 7.05]
    for sgn in (+1, -1):
        ax.plot(xf, [4.0 + sgn * 0.22 * flex(x) for x in xf],
                ls="--", color=BLUE, lw=1.1, zorder=4)
    ax.annotate("tide", xy=(0.55, 2.05), xytext=(0.55, 1.05), color=BLUE,
                ha="center", fontsize=9,
                arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.4))

    # grounding line
    ax.plot([xgl, xgl], [0.2, 2.25], ls=":", color=INK, lw=1.2, zorder=4)
    ax.text(xgl + 0.12, 0.35, "grounding line", fontsize=9, color=INK)

    # the radar sled + survey extent
    ax.add_patch(Rectangle((2.55, 4.02), 0.5, 0.22, fc=ORANGE, ec="white",
                           lw=0.8, zorder=6))
    ax.plot([2.8, 2.8], [4.24, 4.55], color=ORANGE, lw=1.4, zorder=6)
    ax.add_patch(FancyArrowPatch((1.0, 4.75), (5.6, 4.75),
                                 arrowstyle="<->", color=INK, lw=1.2,
                                 mutation_scale=12, zorder=6))
    ax.text(3.3, 4.9, "same 5 km line, 13 passes over 4 days",
            fontsize=9.5, color=INK, ha="center")

    # ApRES
    ax.plot([6.3, 6.3], [4.0, 4.35], color=MAGENTA, lw=2.0, zorder=6)
    ax.plot(6.3, 4.38, "v", ms=6, color=MAGENTA, zorder=6)
    ax.text(6.3, 4.52, "ApRES", fontsize=9, color=MAGENTA, ha="center")

    # column strain arrows within the shelf
    ax.add_patch(FancyArrowPatch((2.8, 3.62), (2.8, 3.30),
                                 arrowstyle="-|>", color=INK, lw=1.2,
                                 mutation_scale=13, zorder=6))
    ax.add_patch(FancyArrowPatch((2.8, 2.45), (2.8, 2.77),
                                 arrowstyle="-|>", color=INK, lw=1.2,
                                 mutation_scale=13, zorder=6))
    ax.text(3.0, 2.95, "column\nstrain", fontsize=8.5, color=INK, va="center")

    ax.text(1.0, 3.05, "floating shelf", fontsize=9, color=INK_SOFT)
    ax.text(8.35, 3.0, "grounded ice", fontsize=9, color=INK_SOFT)
    ax.text(1.0, 1.0, "ocean", fontsize=9, color=INK_SOFT)

    ax.text(0.1, 5.35, "(a) the setting: the radar rides the tide with the surface",
            fontsize=11, color=INK, weight="bold")
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.7)
    ax.axis("off")


def panel_b(ax):
    zs = [1.2, 3.1, 5.4]              # layer depths (drawing units)
    grow = [0.10, 0.28, 0.52]         # displacement growing with depth
    for k, (x0, col, lbl) in enumerate(
            [(1.6, BLUE, "pass at low tide"),
             (4.4, ORANGE, "pass at high tide")]):
        ax.plot([x0, x0], [0.4, 7.6], color=INK_SOFT, lw=0.8)
        # surface return, aligned between passes
        ax.plot([x0 - 0.35, x0 + 0.35], [7.2, 7.2], color=col, lw=3.0)
        for z, g in zip(zs, grow):
            # rising tide thins the column (the measured sign), so the
            # high-tide pass's layers sit SHALLOWER relative to the surface
            y = 7.2 - z + (g if k == 1 else 0.0)
            ax.plot([x0 - 0.28, x0 + 0.28], [y, y], color=col, lw=2.0)
        ax.text(x0, 7.95, lbl, fontsize=9, color=col, ha="center")
    # connectors showing the growing offset
    for z, g in zip(zs, grow):
        ax.annotate("", xy=(4.05, 7.2 - z + g), xytext=(1.95, 7.2 - z),
                    arrowprops=dict(arrowstyle="-", ls=":", color=INK_SOFT))
    ax.text(3.0, 6.55, "surface returns aligned\n(per-column coalignment)",
            fontsize=8.5, color=INK, ha="center")

    # dh(z) profile at right
    ax.plot([6.6, 6.6], [1.0, 7.2], color=INK_SOFT, lw=0.8)
    prof_x = [6.6 + 2.6 * (g / grow[-1]) for g in [0.0] + grow]
    prof_y = [7.2] + [7.2 - z for z in zs]
    ax.plot(prof_x, prof_y, "-o", color=INK, lw=1.8, ms=4)
    ax.text(8.1, 4.1, "dh(z)", fontsize=10, color=INK)
    ax.text(6.6, 0.35, "layer displacement relative to the surface;\n"
                       "strain of the column = dh(z) / z",
            fontsize=8.5, color=INK, ha="left")

    ax.text(0.1, 8.9, "(b) the observable: layers move relative to the surface",
            fontsize=11, color=INK, weight="bold")
    ax.set_xlim(0, 10); ax.set_ylim(0, 9.4)
    ax.axis("off")


def panel_c(ax):
    # tide over four days with the real-ish clustered pass times
    tdays = [0.10, 0.14, 0.55, 0.60, 1.10, 1.15, 1.62, 2.10, 2.15,
             2.62, 3.12, 3.20, 3.28]
    tide = lambda t: math.sin(2 * math.pi * t / 1.0) * 0.8
    ts = [i * 0.01 for i in range(0, 341)]
    ax.plot(ts, [tide(t) + 6.5 for t in ts], color=INK_SOFT, lw=1.2)
    pts = [(t, tide(t) + 6.5) for t in tdays]
    # all-pairs chords, very light: the network
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]],
                    color=BLUE, lw=0.6, alpha=0.22, zorder=2)
    for t, y in pts:
        ax.plot(t, y, "o", ms=6, mfc=ORANGE, mec="white", mew=0.8, zorder=3)
    ax.text(1.7, 8.15, "13 passes sample the tide;\nevery pair is an interferogram "
                       "(78) - inverted as one network, no reference pass",
            fontsize=8.5, color=INK, ha="center")
    ax.text(3.42, tide(3.42) + 6.75, "tide(t)", fontsize=9, color=INK_SOFT)

    # the joint fit
    ax.text(1.7, 4.55, r"$x_k = a + b\,t_k + c\,\mathrm{tide}(t_k)$",
            fontsize=12, color=INK, ha="center")
    ax.annotate("secular strain rate", xy=(1.35, 4.35), xytext=(0.35, 3.5),
                fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK_SOFT, lw=0.9))
    ax.annotate("tidal response", xy=(2.35, 4.35), xytext=(2.45, 3.5),
                fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK_SOFT, lw=0.9))

    ax.text(0.05, 8.9, "(c) the estimation: network + joint time-and-tide fit",
            fontsize=11, color=INK, weight="bold")
    ax.set_xlim(-0.1, 3.5); ax.set_ylim(2.8, 9.4)
    ax.axis("off")


def main():
    fig = plt.figure(figsize=(12.6, 8.6), dpi=150)
    fig.patch.set_facecolor("white")
    axa = fig.add_axes([0.03, 0.52, 0.94, 0.46])
    axb = fig.add_axes([0.03, 0.02, 0.46, 0.46])
    axc = fig.add_axes([0.53, 0.02, 0.44, 0.46])
    panel_a(axa)
    panel_b(axb)
    panel_c(axc)
    out = os.path.join(ROOT, "figs", "EAGER_2022_method_schematic.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
