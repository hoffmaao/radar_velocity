#!/usr/bin/env python3
"""Elasticity schematic - the concept figure for the effective modulus.

The companion to method_schematic.py: that figure shows how repeat-pass
radar measures deformation, this one shows how the deformation becomes a
material property. Same palette, same conventions - no (a)/(b)/(c)
lettering, no titles inside the panels, the lettering and caption added
separately in the slide or manuscript.

  (a) THE SETTING AND THE BEAM. Ice grounded on the left, floating on the
      right, thickness h(x) from the tracked bed. The tide loads the
      underside; the ice bends. The model is an elastic beam of varying
      thickness on a hydrostatic foundation, clamped where flexure stops
      and flat where the shelf floats free.

  (b) THE DERIVATION, four steps from the beam equation to the modulus.
      Each step is a line of algebra with its terms named, so the reader
      can see exactly where E* enters and what it is divided by.

  (c) THE TWO OBSERVABLES, one beam. The surface admittance a(x) is the
      beam's DEFLECTION, normalised by a line mean, so it constrains the
      shape only. The englacial admittance dh(x) is its CURVATURE, in
      absolute units - and curvature amplitude goes as D^(-1/2), which is
      the amplitude equation the shape cannot supply. Drawn from the
      model at the fitted parameters.

  (d) THE ESTIMATION. Misfit over (E*, clamp) is a narrow diagonal
      valley: a stiffer beam clamped further landward fits nearly as
      well. The amplitude is eliminated in closed form, the two remaining
      parameters are searched, and the joint fit closes the valley.

Run from the repo root:  python3 scripts/figures/elasticity_schematic.py
Writes figs/EAGER_2022_elasticity_schematic.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon

BLUE = "#2a78d6"; ORANGE = "#eb6834"; AQUA = "#1baf7a"
YELLOW = "#eda100"; MAGENTA = "#e87ba4"
INK = "#000000"; INK_SOFT = "#737373"
ICE = "#eef3f8"; OCEAN = "#cfe3f5"; BED = "#d9cfc0"; LAYER = "#9db8d2"

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# The fitted values this project measured, used so the drawn curves are
# the real solution rather than a sketch.
E_STAR = 3.90e9          # Pa, joint GNSS + englacial strain, GL1/GL3/GL4
X0_STAR = -750.0         # m, fitted clamp position in the seaward frame
NU = 0.3
RHO_W = 1028.0
G = 9.81
H_ICE = 285.0            # m, tracked-bed thickness mid-line
ZR = 100.0               # m, the depth the column change is read at


def beam(x, E=E_STAR, h=H_ICE, A0=1.0):
    """Analytic clamped beam on a hydrostatic foundation.

    w = A0 [1 - exp(-x/l)(cos x/l + sin x/l)], l = (4D/rho_w g)^(1/4),
    the constant-thickness solution vdef.beamFlexure reproduces to 2e-4
    of A0 (test_beam_flexure check 1). Used here because the schematic
    wants the shape, not the varying-thickness correction.
    """
    D = E * h ** 3 / (12 * (1 - NU ** 2))
    lam = (4 * D / (RHO_W * G)) ** 0.25
    s = np.maximum(x, 0.0) / lam
    w = A0 * (1 - np.exp(-s) * (np.cos(s) + np.sin(s)))
    # curvature, analytic: w'' = A0 * (2/l^2) exp(-s)(cos s - sin s)
    d2w = A0 * (2 / lam ** 2) * np.exp(-s) * (np.cos(s) - np.sin(s))
    w = np.where(x < 0, 0.0, w)
    d2w = np.where(x < 0, 0.0, d2w)
    return w, d2w, lam


def panel_a(ax):
    """The setting: grounded ice, the flexure zone, the floating shelf,
    with the beam's boundary conditions marked where they apply."""
    x = np.linspace(-1.6, 5.0, 500)
    w, _, lam = beam(x * 1e3)
    surf_mean = 0.0
    amp = 0.34                                   # km of drawn deflection

    # ocean and bed
    ax.add_patch(Polygon([(-1.6, -3.05), (5.0, -3.05), (5.0, -1.55),
                          (-1.6, -1.55)], closed=True, fc=OCEAN, ec="none"))
    ax.add_patch(Polygon([(-1.6, -3.05), (-0.35, -3.05), (-0.35, -1.5),
                          (-1.6, -1.15)], closed=True, fc=BED, ec="none"))
    ax.text(-1.15, -2.45, "bed", fontsize=9, color=INK_SOFT)
    # ocean label sits low and left of the load annotation, which occupies
    # the middle of the water column
    ax.text(4.35, -2.9, "ocean", fontsize=9, color=INK_SOFT)

    # the ice: top and bottom surfaces both carry the deflection
    top = surf_mean + amp * w
    bot = top - 1.5
    ax.fill_between(x, bot, top, color=ICE, ec=INK_SOFT, lw=0.9)
    # an internal layer, the thing the radar actually tracks
    ax.plot(x, top - 0.62, "-", color=LAYER, lw=1.1)
    # label the layer over the ice, left of the h(x) arrow, so it neither
    # runs off the panel nor sits on the thickness annotation
    ax.text(0.95, top[np.argmin(abs(x - 0.95))] - 0.52, "internal layer",
            fontsize=8, color=LAYER, va="bottom")

    # thickness annotation
    ax.annotate("", xy=(3.5, top[np.argmin(abs(x - 3.5))]),
                xytext=(3.5, bot[np.argmin(abs(x - 3.5))]),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
    ax.text(3.58, -0.75, r"$h(x)$", fontsize=11, color=INK)

    # the clamp
    ax.plot([0, 0], [-1.75, 0.30], "-", color=INK, lw=2.2)
    ax.text(0.06, 0.40, "clamp: $w=0,\\ dw/dx=0$", fontsize=8.5, color=INK)
    ax.text(-1.5, 0.40, "grounded", fontsize=9, color=INK_SOFT)
    ax.text(3.05, 0.62, "free-floating:  $w = A_0$", fontsize=8.5, color=INK)

    # tidal loading arrows on the underside
    for xa in np.arange(0.45, 4.8, 0.62):
        ib = np.argmin(abs(x - xa))
        ax.add_patch(FancyArrowPatch((xa, bot[ib] - 0.62), (xa, bot[ib] - 0.10),
                                     arrowstyle="-|>", mutation_scale=9,
                                     color=BLUE, lw=1.2))
    ax.text(1.5, -2.62, r"hydrostatic load  $\rho_w g\,[A_0 - w(x)]$",
            fontsize=9.5, color=BLUE)

    # the deflection itself
    ax.annotate("", xy=(2.15, top[np.argmin(abs(x - 2.15))]),
                xytext=(2.15, surf_mean),
                arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.3))
    ax.text(2.22, 0.16, r"$w(x)$", fontsize=11, color=ORANGE)
    ax.plot(x, np.full_like(x, surf_mean), ":", color=INK_SOFT, lw=0.8)

    # flexural length
    ax.annotate("", xy=(0, -1.95), xytext=(lam / 1e3, -1.95),
                arrowprops=dict(arrowstyle="<->", color=INK_SOFT, lw=1.0))
    ax.text(lam / 2e3, -2.2, r"$\ell$", fontsize=10, color=INK_SOFT, ha="center")

    ax.set_xlim(-1.6, 5.0); ax.set_ylim(-3.05, 0.95)
    ax.axis("off")


def panel_b(ax):
    """The derivation: four lines from the loaded beam to the modulus."""
    y = 4.35
    dy = 1.02
    ax.text(0.0, y, r"$\dfrac{d^{2}}{dx^{2}}\!\left[D(x)\,"
                    r"\dfrac{d^{2}w}{dx^{2}}\right] = \rho_w g\,[A_0 - w(x)]$",
            fontsize=12.5, color=INK, va="center")
    ax.text(6.9, y, "beam on a\nhydrostatic foundation", fontsize=8.5,
            color=INK_SOFT, va="center")

    y -= dy
    ax.text(0.0, y, r"$D(x) = \dfrac{E^{*}\,h(x)^{3}}{12\,(1-\nu^{2})}$",
            fontsize=12.5, color=INK, va="center")
    ax.text(6.9, y, "the modulus enters\nONLY through $D$", fontsize=8.5,
            color=INK_SOFT, va="center")
    ax.annotate("", xy=(1.42, y - 0.30), xytext=(1.42, y - 0.62),
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=1.1))
    ax.text(1.52, y - 0.55, r"so $E^{*}$ is never separable from $h^{3}$:"
                            r"  $\;d\ln E^{*} = -3\,d\ln h$",
            fontsize=9, color=ORANGE, va="center")

    y -= 1.72
    ax.text(0.0, y, r"$\ell = \left(\dfrac{4D}{\rho_w g}\right)^{1/4}$",
            fontsize=12.5, color=INK, va="center")
    # mathtext will not parse \tfrac inside \left(...\right); x/\ell is
    # both legal and easier to read at this size
    ax.text(2.55, y, r"$\Rightarrow\ \ w(x) = A_0\left[1 - "
                     r"e^{-x/\ell}\left(\cos(x/\ell) + "
                     r"\sin(x/\ell)\right)\right]$",
            fontsize=11.5, color=INK, va="center")
    ax.text(6.9, y - 0.52, "flexural length: the SHAPE\ncarries $D$",
            fontsize=8.5, color=INK_SOFT, va="center")

    y -= 1.30
    ax.text(0.0, y, r"$\varepsilon_{zz}(x,z) = \dfrac{\nu}{1-\nu}\,"
                    r"(z_n - z)\,\dfrac{d^{2}w}{dx^{2}}$",
            fontsize=12.5, color=INK, va="center")
    ax.text(6.9, y, "bending strain:\nthe CURVATURE", fontsize=8.5,
            color=INK_SOFT, va="center")

    y -= dy
    ax.text(0.0, y, r"$\delta h(z_r) = \dfrac{\nu}{1-\nu}"
                    r"\left(z_n z_r - z_r^{2}/2\right)"
                    r"\dfrac{d^{2}w}{dx^{2}}$",
            fontsize=12.5, color=INK, va="center")
    ax.text(6.9, y, "what the radar\nmeasures, to $z_r$", fontsize=8.5,
            color=INK_SOFT, va="center")

    ax.set_xlim(-0.15, 9.6); ax.set_ylim(-0.55, 5.0)
    ax.axis("off")


def panel_c(ax_w, ax_s):
    """The two observables of one beam, at the fitted parameters."""
    x = np.linspace(-0.6, 4.6, 400)
    w, d2w, _ = beam(x * 1e3)

    # deflection, normalised by its own mean over the observed window -
    # exactly what the GPS line-mean normalisation does to a(x)
    obs = (x > 0.0) & (x < 4.3)
    a_x = w / np.mean(w[obs])
    ax_w.plot(x, a_x, "-", color=ORANGE, lw=2.0)
    ax_w.axhline(0, color=INK_SOFT, lw=0.7, ls=":")
    ax_w.axvline(0, color=INK, lw=1.6)
    ax_w.set_ylabel(r"$a(x)$  (line-mean units)", fontsize=9, color=ORANGE)
    ax_w.text(2.05, 0.30, "SHAPE only:\nthe normalisation\nremoves $A_0$",
              fontsize=8.5, color=ORANGE)
    ax_w.set_xlim(-0.6, 4.6); ax_w.set_ylim(-0.12, 1.55)
    ax_w.set_xticklabels([])
    ax_w.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax_w.spines[s].set_visible(False)

    # englacial: the curvature, absolute, in mm per metre of tide
    zn = H_ICE / 2.0
    k2 = (NU / (1 - NU)) * (zn * ZR - ZR ** 2 / 2)
    # no line-mean division here: unlike a(x) this observable is absolute,
    # which is the whole reason it supplies the amplitude equation
    dh = 1e3 * k2 * d2w
    ax_s.plot(x, dh, "-", color=BLUE, lw=2.0)
    ax_s.axhline(0, color=INK_SOFT, lw=0.7, ls=":")
    ax_s.axvline(0, color=INK, lw=1.6)
    ax_s.set_ylabel(r"$\delta h$  (mm per m tide)", fontsize=9, color=BLUE)
    ax_s.set_xlabel("seaward distance (km)", fontsize=9)
    ax_s.text(2.05, 0.55 * np.nanmax(dh),
              "ABSOLUTE:\namplitude $\\propto D^{-1/2}$", fontsize=8.5, color=BLUE)
    ax_s.set_xlim(-0.6, 4.6)
    ax_s.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax_s.spines[s].set_visible(False)


def panel_d(ax):
    """The misfit valley over (E*, clamp), and what closes it."""
    Eg = np.logspace(np.log10(0.6e9), np.log10(26e9), 150)
    x0 = np.linspace(-1500, 500, 150)
    EE, XX = np.meshgrid(Eg, x0)
    # The valley: E* and clamp trade off because a stiffer beam clamped
    # further landward reproduces the same curve over a short window. The
    # ridge follows constant l + x0, the analytic form of that trade-off,
    # anchored on the FITTED pair so the valley floor runs through the
    # solution the marker draws rather than past it.
    lam = (4 * (EE * H_ICE ** 3 / (12 * (1 - NU ** 2))) / (RHO_W * G)) ** 0.25
    lam_star = (4 * (E_STAR * H_ICE ** 3 / (12 * (1 - NU ** 2)))
                / (RHO_W * G)) ** 0.25
    J = (((lam + XX) - (lam_star + X0_STAR)) / 260.0) ** 2
    J += ((np.log10(EE / E_STAR)) / 1.15) ** 2 * 0.30
    ax.pcolormesh(np.log10(Eg / 1e9), x0 / 1e3, np.minimum(J, 9),
                  cmap="RdYlBu_r", shading="auto", rasterized=True)
    ax.contour(np.log10(Eg / 1e9), x0 / 1e3, J, levels=[1.0],
               colors=[INK], linewidths=1.1)
    ax.plot(np.log10(E_STAR / 1e9), X0_STAR / 1e3, "o", color="white",
            mec=INK, mew=1.4, ms=7)
    ax.text(np.log10(E_STAR / 1e9) + 0.06, X0_STAR / 1e3 + 0.05,
            r"$E^{*} = 3.9$ GPa", fontsize=9, color=INK)
    ax.text(-0.10, -1.34, "shape alone:\nthe valley stays open",
            fontsize=8.5, color=INK)
    ax.annotate("strain adds the\namplitude equation", xy=(0.60, -0.80),
                xytext=(0.72, -1.36), fontsize=8.5, color=BLUE,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.1))
    ax.set_xlabel(r"$\log_{10} E^{*}$  (GPa)", fontsize=9)
    ax.set_ylabel("clamp position (km)", fontsize=9)
    ax.tick_params(labelsize=8)


def main():
    fig = plt.figure(figsize=(12.6, 8.6), dpi=150)
    fig.patch.set_facecolor("white")

    # No panel lettering or headers: added in the slide or manuscript,
    # same convention as method_schematic.py
    ax_a = fig.add_axes([0.035, 0.545, 0.44, 0.40])
    panel_a(ax_a)

    ax_b = fig.add_axes([0.515, 0.535, 0.465, 0.415])
    panel_b(ax_b)

    ax_w = fig.add_axes([0.075, 0.285, 0.35, 0.185])
    ax_s = fig.add_axes([0.075, 0.075, 0.35, 0.185])
    panel_c(ax_w, ax_s)

    ax_d = fig.add_axes([0.575, 0.085, 0.36, 0.375])
    panel_d(ax_d)

    out = os.path.join(ROOT, "figs", "EAGER_2022_elasticity_schematic.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
