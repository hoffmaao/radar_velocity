#!/usr/bin/env python3
"""Method schematic - the concept figure for repeat-pass vertical strain.

Three panels:

  (a) THE INTERFEROMETRY, after the "measuring englacial deformation
      with multipass radar" schematic in the Copenhagen vertical-velocity
      talk. Left: the two-pass geometry - the antenna rides the surface
      and the surface rides the tide, so the vertical baseline B_z IS
      the tide; nadir rays (no air path, no refraction bend) down to an
      internal layer displaced by dh between passes.
      Middle: the experiment - a layered floating column, the same line
      walked 13 times over four days, accumulation and vertical velocity
      acting on the column, the tide heaving it. Right: a REAL
      interferogram from the shelf, converted to vertical displacement
      relative to the surface - the pair with the STRONGEST CREDIBLE
      strain signal in the reprocessed all-pairs network products
      (scan of all 217 gate-passing CSARP_vvel_net pairs, ranked by rms
      dh at 150 m): EAGER_2022_GL3 pass 1 vs pass 11, 2.39 days and
      -0.34 m of tide apart, rms 9.9 mm at 150 m growing with depth into
      the grounding zone. The three GL1 short-interval pairs that rank
      higher were REJECTED as artefacts: their signal does not grow from
      zero at the surface (rms at 100 m >= rms at 150 m), and their
      per-column coalignment collapsed to 3/13 usable windows - the
      retracted-hinge mechanism, not strain.

  (b) THE OBSERVABLE. Two passes at different tide states, surface returns
      aligned (per-column coalignment): internal layers are displaced
      relative to the surface by an amount growing with depth. That
      displacement profile dh(z) IS the column strain.

  (c) THE ESTIMATION. The real GL3 pass times and tides (mean GPS platform
      elevation per pass): thirteen passes sample the tide, every pair is
      an interferogram, and the all-pairs network is inverted for
      per-epoch displacement with no reference pass. The joint fit
      x = a + b*t + c*tide separates the secular rate from the tidal
      response. The pair shown in (a) is highlighted.

Annotation labels are drawn on the figure, but the (a)/(b)/(c) panel
headers and caption lines are not - the lettering is added separately in
the slide or manuscript. The tide/network panel is labelled generically
(displacement vs time) so it reads as the estimation concept. The
embedded interferogram carries no title either - its pair metadata,
read from the npz so it always names the pair actually shown, prints to
stdout instead. Labels are placed clear of every line and arrow
(verified on zoomed crops).

Input: data/EAGER_2022_igram_schematic.npz, extracted on the CReSIS
server from the CSARP_multipass comp_mode 3 product through the production
chain (vdef.coalignPair with the vvel_defaults options, then
vdef.multilook [5 15]), cropped to the coherent column and decimated.
The npz records the product and ref/sec pass indices.

Run from the repo root:  python3 scripts/figures/method_schematic.py
Writes figs/EAGER_2022_method_schematic.png
"""
import datetime
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle

BLUE = "#2a78d6"; ORANGE = "#eb6834"; AQUA = "#1baf7a"
YELLOW = "#eda100"; MAGENTA = "#e87ba4"
INK = "#333333"; INK_SOFT = "#737373"
ICE = "#eef3f8"; OCEAN = "#cfe3f5"; BED = "#d9cfc0"; LAYER = "#9db8d2"

C_LIGHT = 299792458.0

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Diverging map for signed displacement: the house cool/warm pair around a
# neutral midpoint; masked (incoherent) samples in a flat gray that cannot
# be mistaken for zero.
DIV = LinearSegmentedColormap.from_list("vdef_div", [BLUE, "#f2f0ee", ORANGE])
DIV.set_bad("#c9c9c9")


def load_igram():
    fn = os.path.join(ROOT, "data", "EAGER_2022_igram_schematic.npz")
    return np.load(fn)


def boxmean(a, kt, kx):
    """Separable 'same'-support boxcar mean with edge correction."""
    def box1(m, k, axis):
        if k <= 1:
            return m
        m = np.moveaxis(m, axis, 0)
        n = m.shape[0]
        c = np.concatenate([np.zeros((1,) + m.shape[1:], m.dtype),
                            np.cumsum(m, axis=0)], axis=0)
        lo = np.clip(np.arange(n) - k // 2, 0, n)
        hi = np.clip(np.arange(n) + (k - k // 2), 0, n)
        cnt = (hi - lo).reshape((n,) + (1,) * (m.ndim - 1))
        s = (c[hi] - c[lo]) / cnt
        return np.moveaxis(s, 0, axis)
    return box1(box1(a, kt, 0), kx, 1)


def pass_label(t_posix):
    t = datetime.datetime.fromtimestamp(t_posix, datetime.timezone.utc)
    return t.strftime("%d %b %H:%M")


def panel_a_geometry(ax):
    """Two-pass geometry: the antenna rides the surface, which rides the
    tide, so the vertical baseline IS the tide; rays go straight down.

    Everything belonging to one epoch shares one colour: pass i (the
    reference, at the HIGHER tide in the embedded pair) is orange - solid
    surface, marker, label, solid layer; pass j is blue - dashed surface,
    dotted layer. Blue = lower tide matches panel (b)."""
    y_i, y_j = 7.4, 6.5      # surface height at each pass epoch
    # pass i sits right of the B_z label so its ray clears the text
    xi, xj = 3.7, 5.6        # antenna positions along the line

    # coordinate frame, top right: z up, x (along track) right, y into
    # the page
    ox, oy = 8.3, 9.0
    ax.add_patch(FancyArrowPatch((ox, oy), (ox, oy + 0.9), arrowstyle="-|>",
                                 color=INK, lw=1.1, mutation_scale=9))
    ax.add_patch(FancyArrowPatch((ox, oy), (ox + 0.9, oy), arrowstyle="-|>",
                                 color=INK, lw=1.1, mutation_scale=9))
    ax.plot(ox, oy, "o", ms=7, mfc="white", mec=INK, mew=1.0)
    ax.plot(ox, oy, "x", ms=4, color=INK)
    ax.text(ox + 0.05, oy + 0.95, "$z$", fontsize=9, color=INK)
    ax.text(ox + 1.05, oy + 0.1, "$x$", fontsize=9, color=INK)
    ax.text(ox - 0.3, oy - 0.15, "$y$", fontsize=9, color=INK, ha="right")

    # the surface at the two epochs: the whole shelf rides the tide
    ax.plot([0.3, 9.7], [y_i, y_i], color=ORANGE, lw=2.2)
    ax.plot([0.3, 9.7], [y_j, y_j], color=BLUE, lw=1.6, ls="--")
    ax.text(9.55, y_i + 0.25, "air", fontsize=8.5, color=INK_SOFT, ha="right")
    ax.text(9.55, y_j - 0.65, "ice", fontsize=8.5, color=INK_SOFT, ha="right")

    # vertical baseline = the tidal heave of the surface
    ax.add_patch(FancyArrowPatch((0.9, y_j), (0.9, y_i), arrowstyle="<|-|>",
                                 color=AQUA, lw=1.2, mutation_scale=8))
    ax.text(1.25, (y_i + y_j) / 2, "$B_z$ = tide", fontsize=9, color=AQUA,
            ha="left", va="center")

    # antennas ON the surface, one per epoch
    ax.plot(xi, y_i + 0.1, "s", ms=6, color=ORANGE)
    ax.text(xi + 0.4, y_i + 0.22, "pass $i$", fontsize=9, color=ORANGE,
            ha="left", va="bottom")
    ax.plot(xj, y_j + 0.1, "s", ms=6, color=BLUE)
    ax.text(xj + 0.4, y_j + 0.22, "pass $j$", fontsize=9, color=BLUE,
            ha="left", va="bottom")

    # nadir rays to the layer at each epoch
    lay = lambda x: 2.9 - (2.9 - 1.95) / (6.6 - 1.2) * (x - 1.2)
    ax.add_patch(FancyArrowPatch((xi, y_i), (xi, lay(xi) + 0.06),
                                 arrowstyle="-|>", color=AQUA, lw=1.3,
                                 mutation_scale=9))
    ax.add_patch(FancyArrowPatch((xj, y_j), (xj, lay(xj) - 0.55 + 0.06),
                                 arrowstyle="-|>", color=AQUA, lw=1.3,
                                 mutation_scale=9))

    # the layer at each epoch, coloured like its pass: solid orange at the
    # higher-tide epoch i, dotted blue at epoch j, lower with the surface
    ax.plot([1.2, 6.6], [2.9, 1.95], color=ORANGE, lw=1.6)
    xs = np.linspace(1.2, 6.6, 50)
    ax.plot(xs, lay(xs) - 0.55, ls=":", color=BLUE, lw=1.3)
    ax.text(1.15, 3.1, "layer", fontsize=9, color=ORANGE, ha="left",
            va="bottom")
    ax.add_patch(FancyArrowPatch((4.3, lay(4.3)), (4.3, lay(4.3) - 0.55),
                                 arrowstyle="<|-|>", color=INK, lw=1.1,
                                 mutation_scale=8))
    ax.text(4.38, lay(4.3) - 0.85, "$dh$", fontsize=9, color=INK, ha="left",
            va="top")

    ax.set_xlim(0, 10); ax.set_ylim(0.9, 10.4)
    ax.axis("off")


def panel_a_block(ax):
    """The experiment: a layered floating column walked 13 times."""
    X0, X1 = 1.2, 7.6           # front face extent
    SKX, SKY = 1.6, 0.9         # depth-direction skew

    def para(x0, x1, y0, y1, **kw):
        ax.add_patch(Polygon([(x0, y0), (x1, y0), (x1 + SKX, y1 + SKY),
                              (x0 + SKX, y1 + SKY)], closed=True, **kw))

    # ocean under the floating column
    ax.add_patch(Rectangle((X0, 1.9), X1 - X0, 0.7, fc=OCEAN, ec="none"))
    ax.add_patch(Polygon([(X1, 1.9), (X1, 2.6), (X1 + SKX, 3.5),
                          (X1 + SKX, 2.8)], closed=True, fc="#bcd6ee",
                         ec="none"))
    ax.text(X0 + 0.3, 2.22, "ocean", fontsize=8.5, color=INK_SOFT)

    # the ice column
    ax.add_patch(Rectangle((X0, 2.6), X1 - X0, 2.6, fc=ICE, ec=INK_SOFT,
                           lw=0.8))
    ax.add_patch(Polygon([(X1, 2.6), (X1, 5.2), (X1 + SKX, 6.1),
                          (X1 + SKX, 3.5)], closed=True, fc="#dde7f1",
                         ec=INK_SOFT, lw=0.8))
    para(X0, X1, 5.2, 5.2, fc="#f4f8fc", ec=INK_SOFT, lw=0.8)

    # internal layers as stacked sheets, topmost drawn last. Each sheet
    # gets front, top AND right-side faces so it reads as a closed slab
    # rather than an open box at the right-hand end.
    for k in range(4):
        y0 = 5.42 + 0.42 * k
        ax.add_patch(Rectangle((X0, y0), X1 - X0, 0.12, fc=ICE, ec=LAYER,
                               lw=0.8, zorder=3 + k))
        ax.add_patch(Polygon([(X1, y0), (X1, y0 + 0.12),
                              (X1 + SKX, y0 + 0.12 + SKY),
                              (X1 + SKX, y0 + SKY)], closed=True,
                             fc="#dfe9f3", ec=LAYER, lw=0.8, zorder=3 + k))
        para(X0, X1, y0 + 0.12, y0 + 0.12, fc="#f4f8fc", ec=LAYER, lw=0.8,
             zorder=3 + k)

    # the repeated passes on the surface sheet, both walking directions
    ytop = 5.42 + 0.42 * 3 + 0.12
    for f, col, flip in ((0.30, "#444444", False), (0.52, "#7a7a7a", True),
                         (0.74, "#ababab", False)):
        xa = X0 + f * SKX + 0.55
        xb = X1 + f * SKX - 0.55
        ya = ytop + f * SKY
        if flip:
            xa, xb = xb, xa
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, ya), arrowstyle="-|>",
                                     color=col, lw=1.3, mutation_scale=9,
                                     zorder=8))
        for xd in np.linspace(xa, xb, 6)[1:-1]:
            ax.plot(xd, ya, "o", ms=3, mfc="white", mec=col, mew=0.7,
                    zorder=9)
    # what moves the layers
    ax.plot(0.8, 8.0, "v", ms=8, color=YELLOW, zorder=5)
    ax.text(1.15, 8.0, "accumulation", fontsize=8.5, color=INK, va="center")
    ax.add_patch(FancyArrowPatch((0.55, 5.5), (0.55, 3.9), arrowstyle="-|>",
                                 color=MAGENTA, lw=2.4, mutation_scale=13))
    ax.text(0.42, 3.6, "vertical\nvelocity", fontsize=8.5, color=INK,
            ha="center", va="top")

    ax.set_xlim(0, 10); ax.set_ylim(1.2, 10.4)
    ax.axis("off")


def panel_a_igram(ax, cax, d):
    """The real pair, as a map of vertical displacement."""
    ph = np.asarray(d["phase"], dtype=float)
    coh = np.asarray(d["coherence"], dtype=float)
    x = d["x_km"]; z = d["depth_m"]; tw = d["twtt_us"]; nloc = d["n_local"]

    cpx = coh * np.exp(1j * ph)
    cpx[~np.isfinite(cpx)] = 0.0
    sm = boxmean(cpx, 3, 9)

    # surface reference, as in the production chain: the phase 50 ns below
    # the surface return defines zero for each column
    iref = int(np.argmin(np.abs(tw - 0.05)))
    refc = sm[iref - 1:iref + 2, :].mean(axis=0)
    refc = refc / np.maximum(np.abs(refc), 1e-12)
    phi = np.angle(sm * np.conj(refc)[None, :])

    # metres of layer displacement per phase cycle, at each depth
    cyc_m = C_LIGHT / (2.0 * float(d["fc"]) * nloc)
    dh_mm = -phi / (2 * np.pi) * cyc_m[:, None] * 1e3

    rows = np.isfinite(z) & (z <= 300.0)
    dh = np.where(coh > 0.3, dh_mm, np.nan)[rows]

    pcm = ax.pcolormesh(x, z[rows], np.ma.masked_invalid(dh), cmap=DIV,
                        vmin=-50, vmax=50, shading="nearest", rasterized=True)
    ax.set_ylim(300, 0)
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel("along track (km)", fontsize=9)
    ax.set_ylabel("depth below surface (m)", fontsize=9)
    ax.tick_params(labelsize=8, length=3, color=INK_SOFT)
    for s in ax.spines.values():
        s.set_color(INK_SOFT); s.set_linewidth(0.8)

    dtide = float(d["tide_sec"]) - float(d["tide_ref"])
    dt_days = abs(float(d["t_sec"]) - float(d["t_ref"])) / 86400.0
    prod = str(d["product"]) if "product" in getattr(d, "files", d) else "GL3"
    prod = prod.replace("EAGER_2022_", "")
    print("panel (a) pair (%s): %s vs %s, tide %+.2f m, %.2f days apart" %
          (prod, pass_label(float(d["t_sec"])), pass_label(float(d["t_ref"])),
           dtide, dt_days))

    cb = plt.colorbar(pcm, cax=cax)
    cb.set_label("layer displacement rel. surface, $dh$ (mm)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(INK_SOFT)
    cb.outline.set_linewidth(0.8)


ZS = [1.2, 3.1, 5.4]              # layer depths (drawing units)
GROW = [0.10, 0.28, 0.52]         # displacement growing with depth


def panel_b(ax):
    """The two-pass column sketch; the dh(z) response lives in its own
    formal axes (panel_b_profile)."""
    for k, (x0, col, lbl) in enumerate(
            [(1.6, BLUE, "pass at low tide"),
             (4.4, ORANGE, "pass at high tide")]):
        ax.plot([x0, x0], [0.4, 7.6], color=INK_SOFT, lw=0.8)
        # surface return, aligned between passes
        ax.plot([x0 - 0.35, x0 + 0.35], [7.2, 7.2], color=col, lw=3.0)
        for z, g in zip(ZS, GROW):
            # rising tide thins the column (the measured sign), so the
            # high-tide pass's layers sit SHALLOWER relative to the surface
            y = 7.2 - z + (g if k == 1 else 0.0)
            ax.plot([x0 - 0.28, x0 + 0.28], [y, y], color=col, lw=2.0)
        ax.text(x0, 7.95, lbl, fontsize=9, color=col, ha="center")
    # connectors showing the growing offset
    for z, g in zip(ZS, GROW):
        ax.annotate("", xy=(4.05, 7.2 - z + g), xytext=(1.95, 7.2 - z),
                    arrowprops=dict(arrowstyle="-", ls=":", color=INK_SOFT))
    ax.set_xlim(0.3, 5.8); ax.set_ylim(0, 8.6)
    ax.axis("off")


def panel_b_profile(ax):
    """The layer displacement response dh(z), on formal labelled axes.
    The values are the schematic ones from the column sketch, so the axes
    carry labels but no numeric ticks."""
    dh = [0.0] + GROW
    z = [0.0] + ZS
    ax.plot(dh, z, "-o", color=INK, lw=1.8, ms=4)
    ax.set_xlim(-0.04, GROW[-1] * 1.18)
    ax.set_ylim(ZS[-1] * 1.18, -0.15)               # depth increases down
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("layer displacement, dh(z)", fontsize=9, color=INK)
    ax.set_ylabel("depth below surface, z", fontsize=9, color=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_SOFT)
        ax.spines[side].set_linewidth(0.8)


def panel_c(ax, d):
    """The tide sampling and pair network, on formal labelled axes: real
    dates on x, metres of tide on y."""
    # the real sampling: per-pass times and tides (mean GPS platform
    # elevation) of the GL3 product; x is the day of December 2022 (UTC)
    dec1 = datetime.datetime(2022, 12, 1,
                             tzinfo=datetime.timezone.utc).timestamp()
    t = (d["t_pass"] - dec1) / 86400.0 + 1.0
    eta = d["tide_pass"] - d["tide_pass"].mean()

    # display curve: diurnal + semidiurnal harmonics through the 13 samples
    w1, w2 = 2 * np.pi, 4 * np.pi          # per day
    def design(tt):
        return np.column_stack([np.ones_like(tt), np.cos(w1 * tt),
                                np.sin(w1 * tt), np.cos(w2 * tt),
                                np.sin(w2 * tt)])
    coef, *_ = np.linalg.lstsq(design(t), eta, rcond=None)
    rms = float(np.sqrt(np.mean((design(t) @ coef - eta) ** 2)))
    print(f"panel (c) tide harmonics: rms misfit {rms*1000:.0f} mm")

    ts = np.linspace(t.min() - 0.08, t.max() + 0.3, 400)
    ax.plot(ts, design(ts) @ coef, color=INK_SOFT, lw=1.2, zorder=1)

    # all-pairs chords, very light: the network
    n = len(t)
    for i in range(n):
        for j in range(i + 1, n):
            ax.plot([t[i], t[j]], [eta[i], eta[j]],
                    color=BLUE, lw=0.6, alpha=0.22, zorder=2)
    # the pair embedded in panel (a); older extracts called ref_idx main_idx
    ref_key = "ref_idx" if "ref_idx" in getattr(d, "files", d) else "main_idx"
    ia, ja = int(d["sec_idx"]) - 1, int(d[ref_key]) - 1
    ax.plot([t[ia], t[ja]], [eta[ia], eta[ja]],
            color=ORANGE, lw=1.8, zorder=3)
    # point the callout at the later end of the chord, wherever the pair sits
    ie = ia if t[ia] > t[ja] else ja
    ax.annotate("the pair used in\ndisplacement figure",
                xy=(t[ie] + 0.03, eta[ie] + 0.01),
                xytext=(12.45, -0.55), fontsize=8, color=ORANGE, ha="center",
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9,
                                shrinkB=2))
    for ti, ei in zip(t, eta):
        ax.plot(ti, ei, "o", ms=6, mfc=ORANGE, mec="white", mew=0.8,
                zorder=4)
    ax.text(ts[-1] + 0.04, float((design(ts[-1:]) @ coef)[0]), "tide(t)",
            fontsize=9, color=INK_SOFT, va="center")

    # generic axes: the per-epoch displacement rides the tide in time, so
    # the panel reads as the estimation concept, not a tide gauge record.
    # The samples ARE the real GL3 pass times and tides; the labels stay
    # general and the ticks stay off.
    ax.set_xlim(9.3, 13.3)
    ax.set_ylim(-0.75, 0.75)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("time", fontsize=9, color=INK)
    ax.set_ylabel("displacement", fontsize=9, color=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_SOFT)
        ax.spines[side].set_linewidth(0.8)


def panel_c_formula(ax):
    """The joint fit, with its two terms named."""
    ax.text(5.0, 2.15, r"$x_k = a + b\,t_k + c\,\mathrm{tide}(t_k)$",
            fontsize=12, color=INK, ha="center")
    ax.annotate("secular strain rate", xy=(4.3, 1.75), xytext=(1.2, 0.35),
                fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK_SOFT, lw=0.9))
    ax.annotate("tidal response", xy=(6.3, 1.75), xytext=(7.2, 0.35),
                fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK_SOFT, lw=0.9))
    ax.set_xlim(0, 10); ax.set_ylim(0, 3)
    ax.axis("off")


def main():
    d = load_igram()

    fig = plt.figure(figsize=(12.6, 8.9), dpi=150)
    fig.patch.set_facecolor("white")

    # panel designators and header lines are deliberately absent: the
    # (a)/(b)/(c) lettering is added separately in the slide or manuscript
    ax_geo = fig.add_axes([0.02, 0.535, 0.21, 0.395])
    ax_blk = fig.add_axes([0.245, 0.545, 0.27, 0.375])
    ax_ig  = fig.add_axes([0.565, 0.60, 0.315, 0.30])
    cax    = fig.add_axes([0.888, 0.60, 0.011, 0.30])
    panel_a_geometry(ax_geo)
    panel_a_block(ax_blk)
    panel_a_igram(ax_ig, cax, d)

    axb = fig.add_axes([0.03, 0.125, 0.24, 0.42])
    axp = fig.add_axes([0.30, 0.165, 0.115, 0.30])
    axc = fig.add_axes([0.46, 0.28, 0.45, 0.245])
    axf = fig.add_axes([0.465, 0.125, 0.44, 0.115])
    panel_b(axb)
    panel_b_profile(axp)
    panel_c(axc, d)
    panel_c_formula(axf)

    out = os.path.join(ROOT, "figs", "EAGER_2022_method_schematic.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
