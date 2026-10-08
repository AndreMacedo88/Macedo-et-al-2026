mycolorscheme <- c(
    "astrocyte (striatum)" = "#78b2ff",
    "qNSC1" = "#9955ff",
    "qNSC2" = "#0047ab",
    "aNSC" = "#ff1e1d",
    "TAP" = "#e68339",
    "neuroblast" = "#84c64f",
    "neuroblast/neuron (OB)" = "#003333",
    "neuron" = "#003333",
    "oligodendrocyte" = "#ffae00",
    "off-target" = "#8B4513",
    "OB" = "#A6A6A6",
    "SVZ" = "#D9D9D9"
)

# ----------------------------------------------------------------------
# paper_palette_celltypes: cell-type colours matched to THIS project's
# actual obs[["celltype"]] strings, for the 2026 ASD paper figures.
# Derived from mycolorscheme above, with aNSC split into aNSC1 (original
# red) + aNSC2 (darker red), and keys renamed to the data's exact strings
# (neuroblast -> NB, neuron -> Neurons, astrocyte (striatum) -> Astrocytes,
# oligodendrocyte -> Oligodendrocytes). Data cell types absent from the
# source palette (Microglia, Ependymal cells, Endothelial cells, OPC,
# Olfactory Ensheathing cells, ROB_qNSC2) are omitted and fall back to the
# neutral non-lineage gray in the panels.
# ----------------------------------------------------------------------
paper_palette_celltypes <- c(
    # neurogenic lineage (coloured in Figure 1 panels A, C, E)
    "qNSC1" = "#9955ff",
    "qNSC2" = "#0047ab",
    "aNSC1" = "#ff1e1d",
    "aNSC2" = "#b3001b", # darker shade of the aNSC red (aNSC1 + aNSC2 in this project)
    "TAP"   = "#e68339",
    "NB"    = "#84c64f", # source palette "neuroblast"
    # other cell types present in the data with a clear palette match
    "Neurons"          = "#003333", # source palette "neuron"
    "Astrocytes"       = "#78b2ff", # source palette "astrocyte (striatum)"
    "Oligodendrocytes" = "#ffae00"  # source palette "oligodendrocyte"
)

# ----------------------------------------------------------------------
# Figure 2 (Knoblich human cortical organoids) palettes.
# excitatory cell types -> warm (YlOrBr); inhibitory -> cool (BuPu);
# Astrocytes -> gray. Keys match obs[["celltype"]] exactly.
# ----------------------------------------------------------------------
paper_palette_organoid_celltypes <- c(
    "RG" = "#fed778", "ccRG" = "#fec34f", "vRG" = "#fea937", "oRG" = "#f98e23",
    "IP" = "#ee7416", "L23" = "#db5d0b", "L4" = "#c64902", "L56" = "#a63a03",
    "L6_CThPN" = "#862e05",
    "ccvRG" = "#9ab4d6", "INP" = "#8c8cc1", "CGE_IN" = "#8b60ac",
    "LGE_IN" = "#863092", "CGE_LGE_IN" = "#6e096a",
    "Astrocytes" = "#BDBDBD"
)

paper_palette_lineage <- c(
    "Excitatory neurons" = "#cb4b02", "Inhibitory neurons" = "#6e4b9e",
    "excitatory" = "#cb4b02", "inhibitory" = "#6e4b9e",
    "EXC" = "#cb4b02", "INH" = "#6e4b9e"
)

paper_palette_treatment <- c(
    "induced" = "#d01c8b", "uninduced" = "#4dac26"
)
