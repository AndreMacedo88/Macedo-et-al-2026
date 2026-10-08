# Macedo et al. 2026: code for the 3'UTR alternative polyadenylation analyses

This repository contains the code that produces the analyses, numbers and figure panels of

> A. Lopes Martins Macedo et al., [[title]], [[journal or bioRxiv]] [[year, volume, DOI]].

It covers four data sets: single-cell 3'-end sequencing of the adult mouse ventricular-subventricular zone (Fig. 1, fig. S1), human cortical organoids from the CHOOSE screen (Fig. 2, figs. S2 and S3), post-mortem cortex (Fig. 3) and peripheral blood (Fig. 4, fig. S4).

Mapping of 3' ends and identification of poly(A) sites use the `utrcalling` package, which is distributed separately at [[GitHub URL of the tool]].

## Layout

```
preprocessing/   raw reads -> count matrices, cell annotation, pseudotime
  mouse_svz/       STAR + Cell Ranger, expression analysis, 3'-end workflows, site-level integration
  organoids/       download + Cell Ranger + 3'-end workflow, lineages and pseudotime
  blood/           STAR + 3'-end workflow, gene-level count table
analysis/        statistics behind the reported numbers
  mouse_svz/       Dirichlet-multinomial test, permutation and pseudobulk calibration, effect sizes, enrichment
  organoids/pseudotime/     developmental contrast (same tests as the mouse data)
  organoids/perturbation/   perturbation contrast, enrichment, overlap with the developmental contrast
  blood/           differential site usage, classifiers, within-family analyses
  pas_distance/    distance of poly(A) sites to the polyadenylation signal, mouse and organoids (Fig. 1B)
figures/         one notebook per figure panel (or group of panels)
  fig1/  fig2/development/  fig2/perturbation/  fig4/
  style/           shared plotting style and constants (Python and R)
src/             shared helper module (distance of sites to the polyadenylation signal)
environment/     conda environment files
```

Files are numbered in the order in which they run within each folder. Several notebook and panel names use the numbering of an earlier draft: `fig3*` notebooks in `figures/fig2/perturbation/` make panels of Fig. 2 and fig. S3, and `fig5*` notebooks in `figures/fig4/` make panels of Fig. 4 and fig. S4 (see the panel table below).

## Software

The analyses were run in three conda environments. We suggest recreating them from the files below, which list the package versions used for the paper; other installations of the same packages should work as well.

| Environment | File | Used for |
|---|---|---|
| `apa-python` | `environment/apa-python.yml` | Python notebooks and scripts (Python 3.14, scanpy, anndata, harmonypy, scikit-learn) |
| `apa-r` | `environment/apa-r.yml` | R scripts and R notebooks (R 4.4.3, DRIMSeq, clusterProfiler, ReactomePA, glmnet, pROC; IRkernel) |
| `apa-utrcalling` | `environment/apa-utrcalling.yml` | 3'-end workflows (install the `utrcalling` package from [[GitHub URL of the tool]] into this environment; the paper used version 2.5.3 for the mouse data and 3.1.10 for the organoid and blood data) |

```bash
conda env create -f environment/apa-python.yml
conda env create -f environment/apa-r.yml
conda env create -f environment/apa-utrcalling.yml
```

The preprocessing workflows also need Snakemake [[version]], STAR 2.7.9a and samtools 1.11 (created by Snakemake from the `envs/` files with `--use-conda`), and Cell Ranger 6.1.1 (mouse) and 8.0.1 (organoids), which are distributed by 10x Genomics.

## Data and paths

All code reads and writes data relative to one directory, given by an environment variable:

```bash
export M26_DATA_ROOT=/path/to/data            # required
export M26_RESULTS_ROOT=/path/to/results      # optional; figure panels (default: <repo>/results)
```

Below `M26_DATA_ROOT`, the code expects the directory structure of the server on which the analyses were run (for example `3utr/wildtypes_analysis/2026_analyses/utrcalling/` for the processed mouse data). Notebooks in `preprocessing/` and `analysis/` change into their original working directory below `M26_DATA_ROOT` in their first code cell, so relative paths resolve as they did originally. Figure notebooks write panels to `M26_RESULTS_ROOT/<Fig>/<panel>/`.

| Data | Source | Location below `M26_DATA_ROOT` |
|---|---|---|
| Mouse raw reads and processed matrices | GEO [[accession]] | `3utr/wildtypes_analysis/` [[list the deposited files and their paths]] |
| Mouse cell multiplexing (hashtag) assignment, batches 3 and 4 | [[deposited with the GEO record]] | `external/hto_assignment.csv` |
| Organoid raw reads and metadata | ArrayExpress E-MTAB-13148, Zenodo 8074791 (downloaded by `preprocessing/organoids/01_snakemake_workflow`) | `raw/organoids_knoblich/`, `3utr/organoids/knoblich/` |
| Blood raw reads, count matrices and participant metadata | controlled access, [[repository and accession]] | `raw/blood/`, `3utr/Blood/` |
| Reference genomes and annotations | Ensembl release 102 (GRCm38) and 104 (GRCh38); 10x Genomics `refdata-gex-mm10-2020-A` | `common_files/genomes/`, `common_files/annotations/`, `common_files/indices/`, `common_files/GRCh38_custom`, `references/` |
| Gene annotations | SFARI Gene (release 2026-03-22); Ensembl BioMart mouse-human orthologues | `common_files/annotations/` |

Participant-level tables (sample annotations, sibling-pair sheet) are not included. Templates with their column headers are provided (`*_TEMPLATE.csv`).

## Running the analyses

Activate the environment named in each step. Shell runners take their settings from environment variables, and their defaults are the values used for the paper.

### Mouse (Fig. 1, fig. S1)

1. `preprocessing/mouse_svz/01_wildtype_mice` (Snakemake; STAR alignment and Cell Ranger), then `02_AggWildtypes` (`cellranger aggr`, command in `_cmdline`).
2. Notebooks `03` to `05` (`apa-python`): quality control, integration, pseudotime and cell states.
3. `06_snakemake_pipeline` (Snakemake, `apa-utrcalling`): site calling on the gene-level 3'UTR axis, which gives the site count matrix used everywhere else.
4. Notebooks `07` to `10` (`apa-python`): lineage subset, integration on site counts and proportions, graph smoothing.
5. `analysis/mouse_svz/` in numbered order:
   - `01` (`apa-python`): effect sizes.
   - `02` (`apa-r`): Dirichlet-multinomial test. Run `Rscript 02_260401-dirichlet_multinomial_analysis_optimized_2.R <lineage h5ad> --method joint --spline_df 5 --bootstrap 199 --bootstrap_pval 0.1`.
   - `03` and `04` via their runners: permutation test and pseudobulk test.
   - `05` to `08`: smoothed W, composite ranking, comparison of the three tests, enrichment.
   - `09` and `10` (R): enrichment on the conservative gene set, curated gene sets.
   - `11` and `12`: neurons.
   - `13`: transcriptome vs 3'UTRome (Fig. 1H; it also runs the enrichment of fig. S1K).
   - `14`: confidence intervals of Fig. 1H.
6. Figure notebooks in `figures/fig1/`. Fig. 1B needs the distance analysis below.

### Organoids (Fig. 2, figs. S2 and S3)

1. `preprocessing/organoids/01_snakemake_workflow` (Snakemake): download, Cell Ranger, site calling (`utrcalling`, standard 3' chemistry).
2. Notebooks `02` to `08` (`apa-python`): metadata, lineages, pseudotime, embeddings, transfer of pseudotime to the site matrix.
3. Developmental contrast, `analysis/organoids/pseudotime/` in numbered order:
   - Runners `01`, `03` and `05`: Dirichlet-multinomial, permutation and pseudobulk tests.
   - `07` to `13`: comparison of tests, effect sizes, ranking, metrics, orthologues, ASD-risk genes.
   - `14` and `15`: enrichment.
   - `16`: W vs expression. Needs step 2 of the perturbation contrast first.
4. Perturbation contrast, `analysis/organoids/perturbation/`:
   - Runner `01`: DRIMSeq main effect, with the targeting-guide filter.
   - `03`: cell-type-adjusted contrast. Run `Rscript 03_drimseq_treatment_celltype_pseudobulk.R` with the options given in its header.
   - `04` to `10` in order.
5. Figure notebooks in `figures/fig2/`. Run `development/fig2b_proportions_umap` before `development/supp_3_developmental_landscape`, which reads its cached embedding.

### Distance of poly(A) sites to the polyadenylation signal (Fig. 1B)

`analysis/pas_distance/261007-pas_hexamer_distance_organoid_vs_mouse.ipynb` (`apa-python`; it also calls the `apa-utrcalling` environment) uses the final mouse and organoid site matrices and the GRCm38 and GRCh38 genome sequences. It writes the per-site distances that `figures/fig1/fig1b_pas_distance` plots and the organoid-mouse comparison given in the Methods.

### Blood (Fig. 4, fig. S4)

1. `preprocessing/blood/01-remapping.smk` (STAR), `02-utrmapping_genelevel.smk` (`apa-utrcalling`, gene level), notebook `03` (gene-level count table for quality control).
2. `analysis/blood/` in numbered order (R scripts in `apa-r`, Python in `apa-python`):
   - `04` runs through `04-run_repeated_cv.sh`.
   - `10` and `11` are the seed and initialisation checks of the within-family classifier.
3. Figure notebooks in `figures/fig4/`.

### Post-mortem cortex (Fig. 3)

[[Code for Fig. 3 is not included yet.]]

## Figure panels

| Panel | Notebook (`figures/` unless stated) |
|---|---|
| Fig. 1A, Fig. 4A | schematics, no code |
| Fig. 1B | `fig1/fig1b_pas_distance` |
| Fig. 1C | `fig1/fig1c_proportions_umap` |
| Fig. 1D | `fig1/fig1d_volcano_rhoW` |
| Fig. 1E, fig. S1E | `fig1/fig1e_exemplars` |
| Fig. 1F | `fig1/fig1f_bp_enrichment` |
| Fig. 1G | `fig1/fig1_neuron_global_w` |
| Fig. 1H | `analysis/mouse_svz/13_260701-correlations_over_pseudotime` |
| fig. S1A | `fig1/fig1a_umap` |
| fig. S1B, S1C | `fig1/supp_2_defensibility` |
| fig. S1D, S1F | `fig1/supp_3_metric_panel` |
| fig. S1G | `fig1/supp_8_exemplars_per_mouse` |
| fig. S1H, S1I | `fig1/supp_4_exemplars_genesets` |
| fig. S1J | `fig1/fig1g_disease_enrichment` |
| fig. S1K | `fig1/supp_12_expression_invisible_reactome` |
| Fig. 2A | `fig2/development/fig2b_proportions_umap` |
| Fig. 2B | `fig2/development/fig2c_pseudotime_volcano` |
| Fig. 2C, fig. S2H | `fig2/development/260701-exc_vs_inh_utrome` |
| Fig. 2D, fig. S2G, S2I | `fig2/development/fig2d_xspecies_exemplar` |
| Fig. 2E | `fig2/development/fig2f_pseudotime_enrichment` |
| Fig. 2F, fig. S3G | `fig2/perturbation/fig3a_treatment_volcano` |
| Fig. 2G, fig. S3I | `fig2/perturbation/fig3d_cross_experiment_overlap` |
| Fig. 2H | `fig2/perturbation/fig3b_convergence_ora` |
| fig. S2A | `fig2/development/supp_1_dataset_design` |
| fig. S2B, S2C | `fig2/development/supp_2_defensibility` |
| fig. S2D, S2E | `fig2/development/supp_3_developmental_landscape` |
| fig. S2F | `fig2/development/supp_4_xspecies_wexpr` |
| fig. S3A to S3D | `fig2/perturbation/supp_1_dataset_design` |
| fig. S3E, S3F | `fig2/perturbation/unassigned_panels` |
| fig. S3H | `fig2/perturbation/supp_3_full_enrichment` |
| Fig. 4B | `fig4/fig5b_dtu_volcano` |
| Fig. 4C | `fig4/fig5d_classifier` |
| Fig. 4D, 4E | `fig4/supp_8_within_family_drimseq` |
| Fig. 4F, 4G | `fig4/supp_9_lopo_within_family` |
| fig. S4A | `fig4/supp_1_qc` |
| fig. S4B | `fig4/fig5a_cohort_design` |
| fig. S4C | `fig4/fig5c_convergence` |
| fig. S4D | `fig4/supp_4_classifier_detail` |

Some notebooks also write variants of a panel that are not used in the paper.

## Licence

[[licence]]
