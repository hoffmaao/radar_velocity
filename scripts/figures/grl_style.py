"""AGU/GRL figure style, in one place, for the Python figures.

WHY THIS EXISTS. AGU's requirement is that lettering be close to 8 point
"at final print size" (6 point for sub/superscript), in Arial, Helvetica,
Times or Symbol, with figures 50-170 mm wide - one column 50-85 mm, two
column 105-170 mm. The trap is "at final print size": a figure authored
at 320 mm wide and then scaled into a 170 mm column has every font
multiplied by 0.53, so 12 point set in the script prints at 6.4 - under
the floor, and different in every figure because each was a different
width. Every figure in this repo was in that state.

THE FIX IS TO AUTHOR AT PRINT SIZE. Set the figure width in millimetres
to the column it will occupy, and set font sizes in points. Then what the
script says is what the page gets, no scaling step in between, and
consistency across figures is automatic rather than maintained by hand.

Raster output note: these figures are PNG, so type is rasterised and
AGU's "outlined, converted to curves, or embedded" clause does not bite.
Keep the DPI high enough that 8 point text survives - 400 for line art is
comfortable. If a vector submission is wanted later, the sizing here is
already correct and only the output call changes.

Usage:
    from grl_style import apply, figsize, TWO_COL_MM, BASE_PT
    apply()
    fig = plt.figure(figsize=figsize(TWO_COL_MM, 0.62))
"""
import matplotlib

ONE_COL_MM = 85.0     # AGU one-column figure: 50-85 mm
TWO_COL_MM = 170.0    # AGU two-column figure: 105-170 mm
MIN_MM = 50.0
MAX_MM = 170.0

BASE_PT = 8.0         # AGU: lettering close to 8 pt at final print size
SMALL_PT = 7.0        # dense in-plot annotation only; never the axis text
SUB_PT = 6.0          # AGU's floor, for sub/superscript

# Arial and Helvetica are both on AGU's accepted list; the rest are
# fallbacks so the figure still builds on a machine without them.
FONT_STACK = ['Helvetica', 'Arial', 'Nimbus Sans', 'Liberation Sans',
              'DejaVu Sans']

DPI = 400


def mm2in(x):
    return x / 25.4


def figsize(width_mm, aspect):
    """(width_in, height_in) for a figure `aspect` = height/width."""
    if not (MIN_MM <= width_mm <= MAX_MM):
        raise ValueError(
            f'width {width_mm} mm is outside AGU\'s 50-170 mm range')
    w = mm2in(width_mm)
    return (w, w * aspect)


def figsize_mm(width_mm, height_mm):
    """(width_in, height_in) from explicit millimetres."""
    if not (MIN_MM <= width_mm <= MAX_MM):
        raise ValueError(
            f'width {width_mm} mm is outside AGU\'s 50-170 mm range')
    return (mm2in(width_mm), mm2in(height_mm))


def apply(base_pt=BASE_PT):
    """Install the style. Call before creating any figure."""
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': FONT_STACK,
        'mathtext.fontset': 'dejavusans',
        'font.size': base_pt,
        'axes.titlesize': base_pt,
        'axes.labelsize': base_pt,
        'xtick.labelsize': base_pt,
        'ytick.labelsize': base_pt,
        'legend.fontsize': base_pt,
        'figure.titlesize': base_pt,
        # thinner furniture to match 8 pt text; the defaults are drawn for
        # a much larger canvas and look coarse at print size
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'grid.linewidth': 0.5,
        'lines.linewidth': 1.2,
        'patch.linewidth': 0.6,
        'legend.frameon': False,
        'savefig.dpi': DPI,
        'figure.dpi': DPI,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })
