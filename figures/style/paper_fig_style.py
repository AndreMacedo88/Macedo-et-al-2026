"""Shared figure style for the 2026 ASD paper panels (Python).

Single source of truth for the panel notebooks of Figures 1 and 2:
  * the cell-type palette, imported from `260601-paper_constants.py`
    (loaded via importlib because the filename starts with a digit);
  * Nature-compliant font sizes;
  * matplotlib rcParams tuned for editable-text SVG output;
  * `save_and_show()`, which writes each panel to SVG (editable text) + PNG
    and then renders the *on-disk* file inline, so the notebook cell shows
    exactly what was saved.

Each panel notebook does:
    import sys
    sys.path.insert(0, f"{REPO_ROOT}/figures/style")
    import paper_fig_style as S
    S.apply_rcparams()
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from IPython.display import SVG, Image, display

# ── locate this file's directory and load the digit-prefixed constants ──
_HERE = Path(__file__).resolve().parent


def _load_constants():
    path = _HERE / "260601-paper_constants.py"
    spec = importlib.util.spec_from_file_location("paper_constants", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_constants = _load_constants()

# corrected, project-matched cell-type palette (see the constants file)
paper_palette_celltypes = dict(_constants.paper_palette_celltypes)
# the full original palette, kept available for variations / other figures
mycolorscheme = dict(_constants.mycolorscheme)

# Figure-2 (Knoblich organoid) palettes
paper_palette_organoid_celltypes = dict(_constants.paper_palette_organoid_celltypes)
paper_palette_lineage = dict(_constants.paper_palette_lineage)
paper_palette_treatment = dict(_constants.paper_palette_treatment)

# ── cell-type ordering (progenitor -> neuron) ──────────────────────────
MOUSE_LINEAGE_ORDER = ["qNSC1", "qNSC2", "aNSC1", "aNSC2", "TAP", "NB"]
ORGANOID_EXC_ORDER = ["oRG", "RG", "ccRG", "IP", "L6_CThPN", "L56", "L4", "L23"]
ORGANOID_INH_ORDER = ["vRG", "ccvRG", "INP", "LGE_IN", "CGE_LGE_IN", "CGE_IN"]
ORGANOID_CELLTYPE_ORDER = ORGANOID_EXC_ORDER + ORGANOID_INH_ORDER + ["Astrocytes"]

# ── neutral / semantic colours ──────────────────────────────────────────
BG_GRAY = "#E0E0E0"  # non-lineage / background cells in UMAPs
PTIME_CMAP = "magma"  # pseudotime colormap
COL_SHORT = "#2166AC"  # 3'UTR shortening (volcano)
COL_LONG = "#B2182B"  # 3'UTR lengthening (volcano)
COL_NS = "lightgray"  # not significant

# ── Nature font sizes (pt) ──────────────────────────────────────────────
FONT_LABEL = 7  # axis labels (Nature max standard text)
FONT_TICK = 6  # tick labels (above 5 pt minimum)
FONT_LEGEND = 5.5  # legend text
FONT_LEGEND_TITLE = 6
FONT_ANNOT = 6  # in-panel annotations

# ── output root (the 3utr data tree, NOT the notebooks git repo) ────────
OUT_ROOT = Path(os.environ.get("M26_RESULTS_ROOT", Path(__file__).resolve().parents[2] / "results"))

# ── Nature column grid (single source of truth for panel widths) ─────────
# Page = 183 mm wide (Nature / Nature Medicine double column) x up to 247 mm.
# Half-width panels tile 2-up within the page (85 + gap + 85 ~< 183);
# full-width panels span it (174 + side margins). Heights are chosen per
# panel (content-driven); only the WIDTHS are fixed to keep a clean grid.
GRID_DPI = 300


def mm_to_px(mm):
    """Panel dimension in px at GRID_DPI; matplotlib figsize = px / GRID_DPI."""
    return round(mm / 25.4 * GRID_DPI)


def mm_to_in(mm):
    return mm / 25.4


HALF_W_MM, FULL_W_MM = 85.0, 174.0
HALF_W_PX = mm_to_px(HALF_W_MM)  # 1004 px -> 85.0 mm @ 300 dpi
FULL_W_PX = mm_to_px(FULL_W_MM)  # 2055 px -> 174.0 mm @ 300 dpi


def px_to_in(px):
    """Inches for a pixel dimension at GRID_DPI (matplotlib figsize helper)."""
    return px / GRID_DPI


# ── the tiling grid: more panel sizes than just half / full ────────────────
# Panels tile within the page separated by a fixed GAP_MM gutter, so any set of
# spans that adds up to a whole grid fills the extent EXACTLY:
#     HALF        + gap + HALF                  = FULL   (85 + 4 + 85 = 174)
#     QUARTER     + gap + THREEQUARTER          = FULL
#     THIRD       + gap + TWOTHIRD              = FULL
#     HALF        + gap + QUARTER + gap + QUARTER = FULL  (halves & quarters nest)
# quarter / half / full are ONE nested base-4 grid (2*QUARTER+gap = HALF,
# 2*HALF+gap = FULL), so their members mix-and-match freely on a row. Thirds are
# a SEPARATE base-3 grid: put thirds only next to thirds, never next to quarters.
# The SAME model drives heights, using the page's max height (FULL_H_MM) as the
# full extent, so quarter/third/... heights give a consistent vertical rhythm and
# also tile. To shuffle panels: just swap the *_W_PX / *_H_PX constant a panel
# uses (e.g. HALF_W_PX -> THREEQUARTER_W_PX). For any division not named below
# (fifths, sixths, ...) call grid_w_px(span, ncols) / grid_h_px(span, nrows).
GAP_MM = 4.0  # gutter between tiled panels (== FULL_W_MM - 2 * HALF_W_MM)
FULL_H_MM = 247.0  # Nature / Nature Medicine max figure height


def grid_len_mm(span, n, *, full, gap=GAP_MM):
    """Length (mm) of a panel spanning `span` of `n` equal tiles within `full` mm.

    Tiles are separated by `gap` mm gutters. `grid_len_mm(1, n, ...)` is one tile;
    a set of spans that sums to `n` (with the interleaving gaps) reconstitutes
    `full` exactly, which is what lets panels tile without leftover space.
    """
    col = (full - (n - 1) * gap) / n
    return span * col + (span - 1) * gap


def grid_w_mm(span, n, *, gap=GAP_MM):
    """Width in mm of a `span`-of-`n`-columns panel (page width = FULL_W_MM)."""
    return grid_len_mm(span, n, full=FULL_W_MM, gap=gap)


def grid_h_mm(span, n, *, gap=GAP_MM):
    """Height in mm of a `span`-of-`n`-rows panel (page height = FULL_H_MM)."""
    return grid_len_mm(span, n, full=FULL_H_MM, gap=gap)


def grid_w_px(span, n, **kw):
    return mm_to_px(grid_w_mm(span, n, **kw))


def grid_w_in(span, n, **kw):
    return mm_to_in(grid_w_mm(span, n, **kw))


def grid_h_px(span, n, **kw):
    return mm_to_px(grid_h_mm(span, n, **kw))


def grid_h_in(span, n, **kw):
    return mm_to_in(grid_h_mm(span, n, **kw))


# ── named column widths (mm) ── HALF_W_MM / FULL_W_MM defined above ────────
# EIGHTH is the next rung of the SAME nested base-2 ladder as quarter/half/full
# (2*EIGHTH + gap = QUARTER, 4*EIGHTH + 3*gap = HALF), so it mixes freely with them.
EIGHTH_W_MM = grid_w_mm(1, 8)  # 18.25 mm
QUARTER_W_MM = grid_w_mm(1, 4)  # 40.5 mm  (= grid_w_mm(1, 4))
THIRD_W_MM = grid_w_mm(1, 3)  # 55.33 mm
TWOTHIRD_W_MM = grid_w_mm(2, 3)  # 114.67 mm
THREEQUARTER_W_MM = grid_w_mm(3, 4)  # 129.5 mm

# ── named row heights (mm) ── fractions of the FULL page height ────────────
EIGHTH_H_MM = grid_h_mm(1, 8)  # 27.38 mm  -> half a quarter-height (2 stack -> QUARTER_H)
QUARTER_H_MM = grid_h_mm(1, 4)  # 58.75 mm  -> the go-to "short panel" height
THIRD_H_MM = grid_h_mm(1, 3)  # 79.67 mm  -> the go-to "medium panel" height
HALF_H_MM = grid_h_mm(1, 2)  # 121.5 mm
TWOTHIRD_H_MM = grid_h_mm(2, 3)  # 163.33 mm
THREEQUARTER_H_MM = grid_h_mm(3, 4)  # 184.25 mm

# ── the same widths / heights in px (matplotlib figsize = px / GRID_DPI) ───
EIGHTH_W_PX = mm_to_px(EIGHTH_W_MM)
QUARTER_W_PX = mm_to_px(QUARTER_W_MM)
THIRD_W_PX = mm_to_px(THIRD_W_MM)
TWOTHIRD_W_PX = mm_to_px(TWOTHIRD_W_MM)
THREEQUARTER_W_PX = mm_to_px(THREEQUARTER_W_MM)
# HALF_W_PX, FULL_W_PX defined above
EIGHTH_H_PX = mm_to_px(EIGHTH_H_MM)
QUARTER_H_PX = mm_to_px(QUARTER_H_MM)
THIRD_H_PX = mm_to_px(THIRD_H_MM)
HALF_H_PX = mm_to_px(HALF_H_MM)
TWOTHIRD_H_PX = mm_to_px(TWOTHIRD_H_MM)
THREEQUARTER_H_PX = mm_to_px(THREEQUARTER_H_MM)
FULL_H_PX = mm_to_px(FULL_H_MM)

# ── the same widths / heights in inches (cleanest figsize: no "/ 300") ─────
EIGHTH_W_IN = mm_to_in(EIGHTH_W_MM)
QUARTER_W_IN = mm_to_in(QUARTER_W_MM)
THIRD_W_IN = mm_to_in(THIRD_W_MM)
HALF_W_IN = mm_to_in(HALF_W_MM)
TWOTHIRD_W_IN = mm_to_in(TWOTHIRD_W_MM)
THREEQUARTER_W_IN = mm_to_in(THREEQUARTER_W_MM)
FULL_W_IN = mm_to_in(FULL_W_MM)
EIGHTH_H_IN = mm_to_in(EIGHTH_H_MM)
QUARTER_H_IN = mm_to_in(QUARTER_H_MM)
THIRD_H_IN = mm_to_in(THIRD_H_MM)
HALF_H_IN = mm_to_in(HALF_H_MM)
TWOTHIRD_H_IN = mm_to_in(TWOTHIRD_H_MM)
THREEQUARTER_H_IN = mm_to_in(THREEQUARTER_H_MM)
FULL_H_IN = mm_to_in(FULL_H_MM)

# consistency guardrails: the nested base-2 grid must reproduce the anchors
assert abs(grid_w_mm(2, 4) - HALF_W_MM) < 1e-9, "grid must reproduce HALF_W_MM"
assert abs(grid_w_mm(4, 4) - FULL_W_MM) < 1e-9, "grid must reproduce FULL_W_MM"
assert abs(grid_h_mm(2, 8) - QUARTER_H_MM) < 1e-9, "2 eighths + gap must equal a quarter"
assert abs(GAP_MM - (FULL_W_MM - 2 * HALF_W_MM)) < 1e-9, "gap must match anchors"


def apply_rcparams():
    """Set matplotlib rcParams for Nature-style, Inkscape-editable SVG."""
    mpl.rcParams.update(
        {
            # Nature/Nature Medicine want Arial/Helvetica. Arial is not installed;
            # Liberation Sans is its metric-identical open equivalent (so no reflow).
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Liberation Sans",
                "Arial",
                "Nimbus Sans",
                "DejaVu Sans",
            ],
            # render math in the SAME face as the text: Liberation Sans has every glyph we
            # use (Delta, psi, rho, prime, minus, x, >=), so Greek/variables/operators all
            # match the labels and the R panels -- no DejaVu/STIX mix.
            "mathtext.fontset": "custom",
            "mathtext.rm": "Liberation Sans",
            "mathtext.it": "Liberation Sans:italic",
            "mathtext.bf": "Liberation Sans:bold",
            "mathtext.cal": "Liberation Sans:italic",
            "font.size": FONT_LABEL,
            "axes.labelsize": FONT_LABEL,
            "axes.titlesize": FONT_LABEL,
            "xtick.labelsize": FONT_TICK,
            "ytick.labelsize": FONT_TICK,
            "legend.fontsize": FONT_LEGEND,
            "legend.title_fontsize": FONT_LEGEND_TITLE,
            "svg.fonttype": "none",  # keep text as <text> -> editable in Inkscape
            "pdf.fonttype": 42,  # TrueType (editable) if PDF is ever requested
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.bbox": "standard",  # NOT 'tight' -> preserve the exact figsize
            "figure.dpi": 150,
        }
    )


def palette_rgba(labels, palette, *, default=BG_GRAY, alpha=1.0):
    """(n, 4) RGBA array mapping each label via `palette` (missing -> `default`)."""
    return np.array(
        [to_rgba(palette.get(str(l), default), alpha=alpha) for l in labels],
        dtype=float,
    )


def scatter_random_order(
    ax,
    x,
    y,
    *,
    rgba=None,
    c=None,
    cmap=None,
    vmin=None,
    vmax=None,
    s=2.0,
    alpha=0.7,
    seed=0,
    rasterized=True,
    zorder=2,
    **kw,
):
    """Scatter points in a fixed-seed RANDOM order.

    Every UMAP in the paper must be plotted with cells in random order so that no
    cell-label group is systematically drawn on top of another (project rule). Pass
    EITHER per-point colours as `rgba` (an (n, 4) array, e.g. from `palette_rgba`) for
    categorical colouring, OR a per-point scalar `c` (+ `cmap`) for a continuous colour
    mapping. The permutation is seeded so the figure is reproducible.

    Returns the matplotlib PathCollection (so a colorbar can be attached when `c` is used).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    n = x.shape[0]
    order = np.random.default_rng(seed).permutation(n)
    ssize = np.asarray(s)[order] if np.ndim(s) > 0 else s
    if (rgba is None) == (c is None):
        raise ValueError("pass exactly one of `rgba` (categorical) or `c` (continuous)")
    if c is not None:
        return ax.scatter(
            x[order],
            y[order],
            c=np.asarray(c)[order],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            s=ssize,
            alpha=alpha,
            edgecolors="none",
            rasterized=rasterized,
            zorder=zorder,
            **kw,
        )
    return ax.scatter(
        x[order],
        y[order],
        c=np.asarray(rgba)[order],
        s=ssize,
        alpha=alpha,
        edgecolors="none",
        rasterized=rasterized,
        zorder=zorder,
        **kw,
    )


def benjamini_hochberg(pvals):
    """BH-adjusted p-values (NaNs preserved and excluded from the ranking)."""
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    ps = p[ok]
    n = ps.size
    order = np.argsort(ps)
    ranked = ps[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def permutation_calibrated_p(
    p_chisq,
    p_perm,
    *,
    at_floor,
    fit_xmin=0.3,
    fit_xmax=None,
    floor=1e-300,
):
    """Recalibrate chi-squared LRT p-values onto the within-sample permutation scale.

    The cell-level DM-LRT chi-squared p-values reach implausibly small values (e.g.
    -log10 p ~ 250) that statisticians distrust because each cell is not an independent
    sample (pseudoreplication). This fits a robust linear map on the -log10 scale
    (x = -log10 p_chisq, y = -log10 p_perm) over the genes where the within-sample
    permutation null is *informative* (not at its 1/(B+1) floor), then applies it to all
    genes to produce a "pseudoreplication-calibrated" p-value. The map only ever deflates
    (calibrated -log10 p is capped at the nominal value).

    IMPORTANT: the permutation p is floored at 1/(B+1), so the map is only *measured* up to
    ~-log10(min informative p_perm). Beyond `x_fit_max` the calibrated value is an
    EXTRAPOLATION and is flagged in the returned `extrapolated` mask. On this project's
    mouse data the measured slope is ~0.65 and lambda ~1.05-1.2 (mild), i.e. the chi-sq
    LRT is well calibrated against the permutation null over the measurable range.

    Parameters
    ----------
    p_chisq, p_perm : array-like of per-gene p-values.
    at_floor : boolean array; True where p_perm is at its permutation floor (e.g.
        n_perm_hits == 0) and so must be EXCLUDED from the calibration fit.
    fit_xmin, fit_xmax : restrict the -log10 p_chisq window used to fit the map.

    Returns dict(p_cal, neglog_cal, slope, intercept, x_fit_max, extrapolated, lam_gene).
    """
    x = -np.log10(np.clip(np.asarray(p_chisq, float), floor, 1))
    y = -np.log10(np.clip(np.asarray(p_perm, float), 1e-12, 1))
    at_floor = np.asarray(at_floor, dtype=bool)
    fit = np.isfinite(x) & np.isfinite(y) & (~at_floor) & (x >= fit_xmin)
    if fit_xmax is not None:
        fit &= x <= fit_xmax
    slope, intercept = np.polyfit(x[fit], y[fit], 1)
    x_fit_max = float(np.nanmax(x[fit]))
    y_cal = np.clip(np.minimum(intercept + slope * x, x), 0.0, None)
    return dict(
        p_cal=np.power(10.0, -y_cal),
        neglog_cal=y_cal,
        slope=float(slope),
        intercept=float(intercept),
        x_fit_max=x_fit_max,
        extrapolated=(x > x_fit_max),
        # per-gene inflation implied by the map (for reporting): 10^((1-slope/1)*x)-ish
        lam_gene=np.where(y > 0, x / np.maximum(y, 1e-9), np.nan),
    )


def panel_dir(panel, figure="Fig1"):
    """Return (and create) the output subdirectory for a panel."""
    d = OUT_ROOT / figure / panel
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_and_show(
    fig,
    panel,
    *,
    figure="Fig1",
    fig_dir=None,
    formats=("svg", "png"),
    dpi=300,
    show="png",
):
    """Save `fig` as <panel>.<fmt> and display the on-disk file inline.

    Parameters
    ----------
    fig : matplotlib Figure
    panel : str
        Panel name, e.g. "fig1a_umap". Also the output subdirectory name and
        the file stem. A variation is saved by passing e.g. "fig1a_umap__alt".
    figure : str
        Figure group subdir under the output root (default "Fig1").
    fig_dir : Path or None
        Override the output directory; default OUT_ROOT/<figure>/<panel>/.
    formats : tuple of str
        File formats to write. SVG keeps editable text; PNG is a raster mirror.
    show : {"png", "svg", None}
        Which saved file to render inline. PNG is pixel-faithful to disk; SVG
        is the editable artifact. None suppresses inline display.
    """
    fig_dir = Path(fig_dir) if fig_dir is not None else (OUT_ROOT / figure / panel)
    fig_dir.mkdir(parents=True, exist_ok=True)

    written = {}
    for ext in formats:
        fp = fig_dir / f"{panel}.{ext}"
        fig.savefig(fp, dpi=dpi, facecolor="white")
        written[ext] = fp

    w = round(fig.get_figwidth() * dpi)
    h = round(fig.get_figheight() * dpi)
    paths_str = ", ".join(str(p) for p in written.values())
    print(f"saved {panel}: {paths_str}  ({w}x{h} px @ {dpi} dpi)")

    # close so the figure is NOT auto-shown; we display the stored file instead
    plt.close(fig)
    if show == "png" and "png" in written:
        display(Image(filename=str(written["png"])))
    elif show == "svg" and "svg" in written:
        display(SVG(filename=str(written["svg"])))
    return None
