# Check out https://projects.susielu.com/viz-palette?mode=%22deuteranomaly%22&colors=%5B%22%2378b2ff%22%2C%22%239955ff%22%2C%22%230047ab%22%2C%22%23ff1e1d%22%2C%22%23e68339%22%2C%22%2384c64f%22%2C%22%23003333%22%5D
# for testing out colorblind-friendly palettes.

mycolorscheme = {
    "astrocyte (striatum)": "#78b2ff",
    "qNSC1": "#9955ff",
    "qNSC2": "#0047ab",
    "aNSC": "#ff1e1d",
    "TAP": "#e68339",
    "neuroblast": "#84c64f",
    "neuroblast/neuron (OB)": "#003333",
    "neuron": "#003333",
    "oligodendrocyte": "#ffae00",
    "off-target": "#8B4513",
    "OB": "#A6A6A6",
    "SVZ": "#D9D9D9",
}

# ──────────────────────────────────────────────────────────────────────
# paper_palette_celltypes — cell-type colours matched to THIS project's
# actual obs["celltype"] strings, for the 2026 ASD paper figures.
# Derived from `mycolorscheme` above with two corrections:
#   * the source palette had a single "aNSC"; this project has aNSC1 + aNSC2,
#     so aNSC1 keeps the original red and aNSC2 gets a darker shade.
#   * keys renamed to the data's exact strings ("neuroblast" -> "NB",
#     "neuron" -> "Neurons", "astrocyte (striatum)" -> "Astrocytes",
#     "oligodendrocyte" -> "Oligodendrocytes").
# Cell types present in the data but absent from the source palette
# (Microglia, Ependymal cells, Endothelial cells, OPC,
#  Olfactory Ensheathing cells, ROB_qNSC2) are intentionally omitted and
# fall back to the neutral non-lineage gray in the panels.
# ──────────────────────────────────────────────────────────────────────
paper_palette_celltypes = {
    # neurogenic lineage (coloured in Figure 1 panels A, C, E)
    "qNSC1": "#9955ff",
    "qNSC2": "#0047ab",
    "aNSC1": "#ff1e1d",
    "aNSC2": "#b3001b",  # darker shade of the aNSC red (this project has aNSC1 + aNSC2)
    "TAP": "#e68339",
    "NB": "#84c64f",  # source palette "neuroblast"
    # other cell types present in the data with a clear palette match
    "Neurons": "#003333",  # source palette "neuron"
    "Astrocytes": "#78b2ff",  # source palette "astrocyte (striatum)"
    "Oligodendrocytes": "#ffae00",  # source palette "oligodendrocyte"
}

# ──────────────────────────────────────────────────────────────────────
# Figure 2 (Knoblich human cortical organoids) palettes.
# Cell types are grouped by lineage and ordered progenitor -> neuron:
#   excitatory -> warm (YlOrBr ramp); inhibitory -> cool (BuPu ramp);
#   Astrocytes (excluded from both lineages) -> neutral gray.
# Keys match obs["celltype"] exactly.
# ──────────────────────────────────────────────────────────────────────
# Colours are a light->dark hue ramp that MUST track the corrected differentiation
# order (ORGANOID_EXC_ORDER / ORGANOID_INH_ORDER in paper_fig_style): lightest =
# progenitor, darkest = most-differentiated neuron. Keys are listed in that order so
# the ramp stays in sync; do not reorder keys without re-pairing the hex ramp.
paper_palette_organoid_celltypes = {
    # Excitatory lineage (progenitor -> neuron), light -> dark
    "oRG": "#fed778",
    "RG": "#fec34f",
    "ccRG": "#f98e23",
    "IP": "#ee7416",
    "L6_CThPN": "#db5d0b",
    "L56": "#c64902",
    "L4": "#a63a03",
    "L23": "#862e05",
    # Inhibitory lineage (progenitor -> neuron), light -> dark
    "vRG": "#b8cfe8",
    "ccvRG": "#9ab4d6",
    "INP": "#8c8cc1",
    "LGE_IN": "#8b60ac",
    "CGE_LGE_IN": "#863092",
    "CGE_IN": "#6e096a",
    # excluded from both lineages
    "Astrocytes": "#BDBDBD",
}

# Lineage-level colours (keys match obs["lineage"] + short aliases).
paper_palette_lineage = {
    "Excitatory neurons": "#cb4b02",
    "Inhibitory neurons": "#6e4b9e",
    "excitatory": "#cb4b02",
    "inhibitory": "#6e4b9e",
    "EXC": "#cb4b02",
    "INH": "#6e4b9e",
}

# CRISPR treatment arm (induced vs uninduced); orthogonal magenta/green.
paper_palette_treatment = {
    "induced": "#d01c8b",
    "uninduced": "#4dac26",
}

celltype_hex = [
    "#78b2ff",
    "#9955ff",
    "#0047ab",
    "#ff1e1d",
    "#e68339",
    "#84c64f",
    "#003333",
    "#ffae00",
]

lineage_hex = [
    "#78b2ff",
    "#9955ff",
    "#0047ab",
    "#ff1e1d",
    "#e68339",
    "#84c64f",
    "#003333",
]
