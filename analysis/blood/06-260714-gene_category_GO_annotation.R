#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# GO-based functional-category annotation of the within-family DTU genes (S5h panel, Fig. 4D, 4E).
# -----------------------------------------------------------------------------
# The S5h lollipop originally coloured each gene by one hand-curated theme (a theme dict).
# That is incomplete: many genes belong to >1 category (e.g. CYTH1 is a GTPase/cytoskeleton
# ARF-GEF that is also immune - highly expressed in NK/T cells). This script replaces the
# hand curation with a reproducible, multi-membership assignment from Gene Ontology:
#   membership(gene, category) = TRUE iff the gene is annotated (propagated, GOALL = term +
#   all ancestors) to any anchor GO term of that category. A gene can be in several categories.
# Anchor GO terms per Fig. 4 theme (broad, standard roots so propagation is authoritative):
#   immune/myeloid       : GO:0002376 immune system process (BP)
#   RNA-binding/poly(A)   : GO:0003723 RNA binding (MF), GO:0006396 RNA processing (BP)
#   GTPase/cytoskeleton  : GO:0005856 cytoskeleton (CC), GO:0003924 GTPase activity (MF),
#                          GO:0007264 small GTPase mediated signal transduction (BP),
#                          GO:0031267 small GTPase binding (MF)
#   chromatin            : GO:0006325 chromatin organization (BP), GO:0000785 chromatin (CC)
# Also writes, per (gene,category), the gene's own (direct) GO terms that fall under the anchor
# (the evidence), so each attribution can be inspected. Genes in none -> "other".
# Env: conda apa-r. Output: DRIMSeq_allpairs_genewise/gene_category_membership.csv (+ _evidence.csv)
# =============================================================================
DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})
D <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/DRIMSeq_allpairs_genewise")
dtu <- read.csv(
  file.path(D, "genelevel_LOPO_tight.csv"),
  stringsAsFactors = FALSE
)
dtu <- dtu[!is.na(dtu$adj_pvalue) & dtu$adj_pvalue < 0.10, ] # annotate the FDR<0.1 set (plot uses <0.05)

clean_sym <- function(s) {
  s <- as.character(s)
  if (s == "GIMAP1-GIMAP5|GIMAP5") {
    return("GIMAP5")
  }
  sub("[|/].*$", "", sub("^[|,/ ]+", "", s))
}
genes <- unique(vapply(dtu$symbol, clean_sym, character(1)))
# H3-3B / RIPOR2 etc. are current HGNC symbols; map symbol, fall back to ALIAS for any miss.
entrez <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL"))
miss <- names(entrez)[is.na(entrez)]
if (length(miss)) {
  alias <- suppressWarnings(suppressMessages(mapIds(
    org.Hs.eg.db,
    miss,
    "ENTREZID",
    "ALIAS"
  )))
  entrez[miss] <- alias
}
message(
  "mapped ",
  sum(!is.na(entrez)),
  "/",
  length(genes),
  " genes to ENTREZ",
  if (any(is.na(entrez))) {
    paste0(
      " (unmapped: ",
      paste(names(entrez)[is.na(entrez)], collapse = ","),
      ")"
    )
  } else {
    ""
  }
)

CATS <- list(
  `immune/myeloid` = c("GO:0002376"),
  `RNA-binding/poly(A)` = c("GO:0003723", "GO:0006396"),
  `GTPase/cytoskeleton` = c(
    "GO:0005856",
    "GO:0003924",
    "GO:0007264",
    "GO:0031267"
  ),
  `chromatin` = c("GO:0006325", "GO:0000785")
)

ez <- entrez[!is.na(entrez)]
# GOALL = propagated (term + ancestors) -> authoritative membership
goall <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unname(ez),
  columns = "GOALL",
  keytype = "ENTREZID"
))
goall_by <- split(
  unique(goall[, c("ENTREZID", "GOALL")])$GOALL,
  unique(goall[, c("ENTREZID", "GOALL")])$ENTREZID
)
# direct GO terms (for evidence) + ancestor maps to test 'under anchor'
godir <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unname(ez),
  columns = c("GO", "ONTOLOGY"),
  keytype = "ENTREZID"
))
godir <- godir[!is.na(godir$GO), ]
ANC <- list(
  BP = as.list(GOBPANCESTOR),
  MF = as.list(GOMFANCESTOR),
  CC = as.list(GOCCANCESTOR)
)
term_name <- function(gos) {
  if (!length(gos)) {
    return(character(0))
  }
  unname(vapply(
    gos,
    function(g) {
      tt <- suppressWarnings(AnnotationDbi::Term(GOTERM[[g]]))
      if (is.null(tt)) g else tt
    },
    character(1)
  ))
}

memb <- data.frame(
  symbol = genes,
  entrez = entrez[genes],
  stringsAsFactors = FALSE
)
for (cn in names(CATS)) {
  memb[[cn]] <- FALSE
}
ev_rows <- list()
for (g in genes) {
  e <- entrez[[g]]
  if (is.na(e)) {
    next
  }
  gset <- goall_by[[as.character(e)]]
  gd <- godir[godir$ENTREZID == e, ]
  for (cn in names(CATS)) {
    anchors <- CATS[[cn]]
    is_member <- any(anchors %in% gset)
    memb[memb$symbol == g, cn] <- is_member
    if (is_member) {
      # evidence = this gene's direct terms that are (or descend from) an anchor
      supp <- character(0)
      for (i in seq_len(nrow(gd))) {
        t <- gd$GO[i]
        ont <- gd$ONTOLOGY[i]
        anc_i <- c(t, ANC[[ont]][[t]])
        if (any(anchors %in% anc_i)) supp <- c(supp, t)
      }
      supp <- unique(supp)
      ev_rows[[paste(g, cn)]] <- data.frame(
        symbol = g,
        category = cn,
        n_terms = length(supp),
        evidence_terms = paste(head(term_name(supp), 8), collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  }
}
catcols <- names(CATS)
memb$n_categories <- rowSums(memb[, catcols])
memb$categories <- apply(memb[, catcols], 1, function(x) {
  paste(catcols[x], collapse = ";")
})
memb$categories[memb$categories == ""] <- "other"

write.csv(memb, file.path(D, "gene_category_membership.csv"), row.names = FALSE)
ev <- do.call(rbind, ev_rows)
write.csv(ev, file.path(D, "gene_category_GO_evidence.csv"), row.names = FALSE)

message("\n== category membership (GO-propagated; multi-membership allowed) ==")
for (g in genes) {
  row <- memb[memb$symbol == g, ]
  message(sprintf("  %-10s [%d cat] %s", g, row$n_categories, row$categories))
}
message(
  "\nmulti-category genes: ",
  paste(memb$symbol[memb$n_categories >= 2], collapse = ", ")
)
message(
  "per-category totals: ",
  paste(sprintf("%s=%d", catcols, colSums(memb[, catcols])), collapse = "  ")
)
message("saved -> ", D, "/gene_category_membership.csv (+ _evidence.csv)")

# =============================================================================
# Curation: documented, paper-level refinement of the systematic GO calls.
# GO propagation (i) over-assigns via promiscuous single terms - the bare 'RNA binding' root from
# proteome-wide RNA-interactome screens, or lone spurious CC terms (centrosome / mitotic spindle /
# actin cytoskeleton attached to non-cytoskeletal proteins) - and (ii) under-assigns sparsely
# annotated genes. Each DROP/ADD below is explicit and has a stated reason (written to the curation log),
# so the final figure categories can be checked.
# =============================================================================
DROP <- list(
  c(
    "GRB2",
    "RNA-binding/poly(A)",
    "only the bare 'RNA binding' root (RNA-interactome capture); GRB2 is a Ras/immune-receptor signalling ADAPTOR, not an RBP"
  ),
  c(
    "VIM",
    "RNA-binding/poly(A)",
    "only 'double-stranded RNA binding'; vimentin is an intermediate-filament structural protein"
  ),
  c(
    "VIM",
    "immune/myeloid",
    "only 'cellular response to type II interferon'; peripheral for a structural protein"
  ),
  c(
    "APEX1",
    "GTPase/cytoskeleton",
    "only 'centrosome'; spurious for a DNA-repair/redox factor"
  ),
  c(
    "HNRNPC",
    "GTPase/cytoskeleton",
    "only 'actin cytoskeleton'; hnRNPC is a nuclear splicing RBP"
  ),
  c(
    "HNRNPC",
    "chromatin",
    "secondary; hnRNPC is canonically an RNA-binding/splicing protein"
  ),
  c(
    "CXCR2",
    "GTPase/cytoskeleton",
    "only 'mitotic spindle / microtubule cytoskeleton'; spurious for a chemokine GPCR"
  ),
  c(
    "KIF13B",
    "immune/myeloid",
    "only 'T cell activation'; KIF13B is a kinesin motor (cytoskeleton)"
  ),
  c(
    "YBX3",
    "GTPase/cytoskeleton",
    "only 'small GTPase binding'; YBX3 is a Y-box RNA-binding protein"
  ),
  c(
    "SNCA",
    "GTPase/cytoskeleton",
    "only 'actin cytoskeleton'; secondary for alpha-synuclein"
  ),
  c(
    "SNCA",
    "immune/myeloid",
    "GO terms are microglia/neuro-specific; tenuous in peripheral blood"
  )
)
ADD <- list(
  c(
    "GIMAP5",
    "immune/myeloid",
    "GTPase of Immunity-Associated Protein 5; essential for T-lymphocyte survival (canonical immune GTPase); GO annotation sparse"
  ),
  c(
    "GIMAP5",
    "GTPase/cytoskeleton",
    "AIG1-family GTP-binding protein (a GTPase by definition)"
  ),
  c(
    "MS4A6A",
    "immune/myeloid",
    "MS4A-family monocyte/macrophage/microglia surface marker; established myeloid immune gene; GO sparse"
  ),
  c(
    "CYTH1",
    "immune/myeloid",
    "cytohesin-1 activates LFA-1 integrin adhesion in lymphocytes and is enriched in NK/T cells (functional + expression immune role); GO lists only its ARF-GEF role"
  )
)
cur <- memb
log_rows <- list()
add_log <- function(g, cn, act, reason) {
  log_rows[[length(log_rows) + 1]] <<-
    data.frame(
      symbol = g,
      category = cn,
      action = act,
      reason = reason,
      stringsAsFactors = FALSE
    )
}
for (d in DROP) {
  if (d[1] %in% cur$symbol && isTRUE(cur[cur$symbol == d[1], d[2]])) {
    cur[cur$symbol == d[1], d[2]] <- FALSE
    add_log(d[1], d[2], "drop", d[3])
  }
}
for (d in ADD) {
  if (d[1] %in% cur$symbol && !isTRUE(cur[cur$symbol == d[1], d[2]])) {
    cur[cur$symbol == d[1], d[2]] <- TRUE
    add_log(d[1], d[2], "add", d[3])
  }
}
cur$n_categories <- rowSums(cur[, catcols])
cur$categories <- apply(cur[, catcols], 1, function(x) {
  paste(catcols[x], collapse = ";")
})
cur$categories[cur$categories == ""] <- "other"
write.csv(cur, file.path(D, "gene_category_curated.csv"), row.names = FALSE)
write.csv(
  do.call(rbind, log_rows),
  file.path(D, "gene_category_curation_log.csv"),
  row.names = FALSE
)
message("\n== CURATED (GO + documented edits) ==")
for (g in genes) {
  row <- cur[cur$symbol == g, ]
  message(sprintf("  %-10s [%d] %s", g, row$n_categories, row$categories))
}
message(
  "curated per-category totals: ",
  paste(sprintf("%s=%d", catcols, colSums(cur[, catcols])), collapse = "  ")
)
message(
  length(log_rows),
  " curation edits -> gene_category_curated.csv (+ _curation_log.csv)"
)
