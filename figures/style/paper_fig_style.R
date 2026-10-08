# Shared figure style for the 2026 ASD paper panels (R).
#
# Sources the cell-type palette from 260601-paper_constants.R and adds
# Nature-compliant font sizes, a publication ggplot theme, and a
# save_and_show() helper that writes each panel to SVG (editable text, via
# svglite) + PNG and then renders the on-disk file inline so the cell shows
# exactly what was saved. If svglite is unavailable it falls back to an
# editable-text cairo_pdf.
#
# Each R panel notebook does:
#   source(file.path(REPO_ROOT, "figures", "style", "paper_fig_style.R"))

suppressPackageStartupMessages({
    library(ggplot2)
})

# locate this file's directory (when sourced); fall back to the known path
.paper_style_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) NA_character_
)
if (length(.paper_style_dir) != 1 || is.na(.paper_style_dir)) {
    .paper_style_dir <- file.path(Sys.getenv("M26_REPO_ROOT"), "figures", "style")
}
source(file.path(.paper_style_dir, "260601-paper_constants.R"))
# `paper_palette_celltypes` and `mycolorscheme` are now in scope.

# ── Nature font sizes (pt) ──────────────────────────────────────────────
FONT_LABEL <- 7
FONT_TICK <- 6
FONT_LEGEND <- 5.5
FONT_LEGEND_TITLE <- 6
FONT_ANNOT <- 6
LAYOUT_DPI <- 300

# ── output root (the 3utr data tree, not the notebooks git repo) ────────
OUT_ROOT <- Sys.getenv("M26_RESULTS_ROOT", file.path(.paper_style_dir, "..", "..", "results"))

# ── Nature column grid (single source of truth for panel widths) ─────────
# Page = 183 mm wide (Nature / Nature Medicine double column). Half-width
# panels tile 2-up; full-width spans the page. Only widths are fixed.
mm_to_in <- function(mm) mm / 25.4
HALF_W_MM <- 85.0
FULL_W_MM <- 174.0
HALF_W_IN <- HALF_W_MM / 25.4   # 3.3465 in -> 85 mm
FULL_W_IN <- FULL_W_MM / 25.4   # 6.8504 in -> 174 mm

# ── the tiling grid: more panel sizes than just half / full ────────────────
# (mirror of paper_fig_style.py; see that file for the full explanation).
# Panels tile within the page separated by a fixed GAP_MM gutter, so any set of
# spans that adds up to a whole grid fills the extent exactly:
#     half    + gap + half          = full   (85 + 4 + 85 = 174)
#     quarter + gap + threequarter  = full ;   third + gap + twothird = full
# quarter/half/full are one nested base-4 grid (mix freely on a row); thirds are
# a separate base-3 grid (put thirds only next to thirds). Heights use the page's
# max height (FULL_H_MM) as the full extent for a consistent vertical rhythm.
# To reshuffle a panel just swap the *_W_IN / *_H_IN constant it uses; for any
# division not named here call grid_w_in(span, ncols) / grid_h_in(span, nrows).
GAP_MM <- 4.0        # gutter between tiled panels (== FULL_W_MM - 2 * HALF_W_MM)
FULL_H_MM <- 247.0   # Nature / Nature Medicine max figure height

grid_len_mm <- function(span, n, full, gap = GAP_MM) {
    col <- (full - (n - 1) * gap) / n
    span * col + (span - 1) * gap
}
grid_w_mm <- function(span, n, gap = GAP_MM) grid_len_mm(span, n, FULL_W_MM, gap)
grid_h_mm <- function(span, n, gap = GAP_MM) grid_len_mm(span, n, FULL_H_MM, gap)
grid_w_in <- function(span, n, gap = GAP_MM) mm_to_in(grid_w_mm(span, n, gap))
grid_h_in <- function(span, n, gap = GAP_MM) mm_to_in(grid_h_mm(span, n, gap))

# named column widths (mm) -- HALF_W_MM / FULL_W_MM defined above
# Eighth is the next rung of the same nested base-2 ladder as quarter/half/full
# (2*eighth + gap = quarter), so it mixes freely with them.
EIGHTH_W_MM       <- grid_w_mm(1, 8)   # 18.25 mm
QUARTER_W_MM      <- grid_w_mm(1, 4)   # 40.5 mm
THIRD_W_MM        <- grid_w_mm(1, 3)   # 55.33 mm
TWOTHIRD_W_MM     <- grid_w_mm(2, 3)   # 114.67 mm
THREEQUARTER_W_MM <- grid_w_mm(3, 4)   # 129.5 mm
# named row heights (mm): fractions of the full page height
EIGHTH_H_MM       <- grid_h_mm(1, 8)   # 27.38 mm -> half a quarter (2 stack -> QUARTER_H)
QUARTER_H_MM      <- grid_h_mm(1, 4)   # 58.75 mm -> go-to "short panel" height
THIRD_H_MM        <- grid_h_mm(1, 3)   # 79.67 mm -> go-to "medium panel" height
HALF_H_MM         <- grid_h_mm(1, 2)   # 121.5 mm
TWOTHIRD_H_MM     <- grid_h_mm(2, 3)   # 163.33 mm
THREEQUARTER_H_MM <- grid_h_mm(3, 4)   # 184.25 mm

# the same in inches (ggsave / save_and_show() units)
EIGHTH_W_IN       <- mm_to_in(EIGHTH_W_MM)
QUARTER_W_IN      <- mm_to_in(QUARTER_W_MM)
THIRD_W_IN        <- mm_to_in(THIRD_W_MM)
TWOTHIRD_W_IN     <- mm_to_in(TWOTHIRD_W_MM)
THREEQUARTER_W_IN <- mm_to_in(THREEQUARTER_W_MM)
EIGHTH_H_IN       <- mm_to_in(EIGHTH_H_MM)
QUARTER_H_IN      <- mm_to_in(QUARTER_H_MM)
THIRD_H_IN        <- mm_to_in(THIRD_H_MM)
HALF_H_IN         <- mm_to_in(HALF_H_MM)
TWOTHIRD_H_IN     <- mm_to_in(TWOTHIRD_H_MM)
THREEQUARTER_H_IN <- mm_to_in(THREEQUARTER_H_MM)
FULL_H_IN         <- mm_to_in(FULL_H_MM)

# consistency checks: the nested base-2 grid must reproduce the anchors
stopifnot(
    abs(grid_w_mm(2, 4) - HALF_W_MM) < 1e-9,
    abs(grid_w_mm(4, 4) - FULL_W_MM) < 1e-9,
    abs(grid_h_mm(2, 8) - QUARTER_H_MM) < 1e-9,
    abs(GAP_MM - (FULL_W_MM - 2 * HALF_W_MM)) < 1e-9
)

# cell-type ordering (progenitor -> neuron) for Figure 2 (organoid)
ORGANOID_EXC_ORDER <- c("RG", "ccRG", "vRG", "oRG", "IP", "L23", "L4", "L56", "L6_CThPN")
ORGANOID_INH_ORDER <- c("ccvRG", "INP", "CGE_IN", "LGE_IN", "CGE_LGE_IN")
ORGANOID_CELLTYPE_ORDER <- c(ORGANOID_EXC_ORDER, ORGANOID_INH_ORDER, "Astrocytes")

# ── publication theme (panel titles removed; added in Inkscape assembly) ─
theme_paper <- theme_bw(base_size = FONT_LABEL, base_family = "Liberation Sans") +
    theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border       = element_rect(linewidth = 0.6),
        axis.text          = element_text(size = FONT_TICK, color = "black"),
        axis.title         = element_text(size = FONT_LABEL),
        plot.title         = element_blank(),
        legend.text        = element_text(size = FONT_LEGEND),
        legend.title       = element_text(size = FONT_LEGEND_TITLE),
        strip.text         = element_text(size = FONT_TICK + 1, face = "bold"),
        strip.background   = element_rect(fill = "#F0F0F0", linewidth = 0.4)
    )

.have_svglite <- requireNamespace("svglite", quietly = TRUE)

# Save `plot` as <panel>.svg (editable text) + .png, then display the stored
# file inline. width_in/height_in are the physical panel size in inches.
save_and_show <- function(plot, panel, width_in, height_in,
                          figure = "Fig1", fig_dir = NULL,
                          formats = c("svg", "png"), dpi = LAYOUT_DPI,
                          show = "png") {
    if (is.null(fig_dir)) fig_dir <- file.path(OUT_ROOT, figure, panel)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    svg_path <- file.path(fig_dir, paste0(panel, ".svg"))
    png_path <- file.path(fig_dir, paste0(panel, ".png"))

    for (ext in formats) {
        if (ext == "svg") {
            if (.have_svglite) {
                ggsave(svg_path, plot, width = width_in, height = height_in,
                       device = svglite::svglite)
            } else {
                pdf_path <- file.path(fig_dir, paste0(panel, ".pdf"))
                ggsave(pdf_path, plot, width = width_in, height = height_in,
                       device = cairo_pdf)
                message("svglite not installed: wrote editable-text PDF instead of SVG -> ",
                        pdf_path)
            }
        } else if (ext == "png") {
            ggsave(png_path, plot, width = width_in, height = height_in, dpi = dpi)
        }
    }
    cat(sprintf("saved %s to %s (%.2f x %.2f in)\n", panel, fig_dir, width_in, height_in))

    # render the stored file inline so the cell shows exactly what is on disk
    if (show == "png" && file.exists(png_path)) {
        IRdisplay::display_png(file = png_path)
    } else if (show == "svg" && .have_svglite && file.exists(svg_path)) {
        IRdisplay::display_svg(file = svg_path)
    }
    invisible(NULL)
}
