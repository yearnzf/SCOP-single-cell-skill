\---  
name: scop-singlecell  
description: >-  
  Use SCOP (a Seurat-based R toolkit) to run single-cell QC, marker/gene-set  
  heatmaps, SCOP-native expression figures and statistics, cell-cell  
  communication (CellChat backend), and SCENIC/regulon analysis with the pure  
  R/C++ backend (no Python required). Includes portable helper functions and  
  pitfalls learned from real projects.  
version: 1.1.0  
license: MIT  
tags:  
  \- single-cell  
  \- R  
  \- Seurat  
  \- SCOP  
  \- CellChat  
  \- SCENIC  
  \- bioinformatics  
\---

# SCOP single-cell skill

## What this is

A compact, GitHub-ready skill for analysing **Seurat** single-cell objects with  
**SCOP** (`scop`, the Seurat-based analysis toolkit). It documents the functions  
that actually work, provides reusable helpers, and records the pitfalls that cost  
hours in real analyses (device handling for network plots, SCENIC without Python,  
sample-level DE, heatmap sizing, etc.).

It is **data-agnostic**: point the workflows at your own `.rds` file and metadata  
columns.

## Full SCOP coverage

This skill catalogues **all 360 exported functions** of SCOP 0.9.2 in  
`references/function-catalog.md` (machine-readable: `references/scop_functions.csv`)  
and maps every functional domain to a workflow in `references/coverage-matrix.md`.

| Domain                                | Workflow                                                         | Status            |
| ------------------------------------- | ---------------------------------------------------------------- | ----------------- |
| Environment / database prep           | `scripts/run_database_prep.R`                                    | template          |
| Data I/O & conversion                 | `scripts/run_io_workflow.R`                                      | template          |
| QC & preprocessing                    | `scripts/run_qc_workflow.R`                                      | template          |
| Integration & batch correction        | `scripts/run_integration_workflow.R`                             | template          |
| Cell annotation & reference mapping   | `scripts/run_annotation_workflow.R`                              | template          |
| DE & markers                          | `scripts/run_de_enrichment_workflow.R`                           | template          |
| Enrichment / GSEA / GSVA / metabolism | `scripts/run_de_enrichment_workflow.R`                           | template          |
| Gene scoring                          | `scripts/run_scoring_composition_workflow.R`                     | template          |
| Cell composition                      | `scripts/run_scoring_composition_workflow.R`                     | template          |
| Trajectory / pseudotime / velocity    | `scripts/run_trajectory_workflow.R`                              | template          |
| Cell–cell communication               | `scripts/run_ccc_workflow.R`                                     | template          |
| SCENIC / TF                           | `scripts/run_scenic_workflow.R`                                  | template          |
| Spatial / deconvolution               | `scripts/run_spatial_workflow.R`, `run_deconvolution_workflow.R` | template          |
| ATAC / multiome                       | `scripts/run_atac_workflow.R`                                    | template          |
| Perturbation / regulatory networks    | `scripts/run_perturbation_workflow.R`                            | template          |
| Visualisation                         | `scripts/scop_helpers.R` + all workflows                         | recipe            |
| Per-gene expression + pseudobulk      | `scripts/run_expression_workflow.R`                              | **runnable demo** |
| Marker / gene-set heatmaps            | `scripts/run_heatmap_workflow.R`                                 | **runnable demo** |

"Template" means the script is complete but needs your data / external backend;  
check `args(scop::FunctionName)` because SCOP signatures evolve.

## When to use it

Use this skill when you need to:

- make publication-style **marker or gene-set heatmaps** from a Seurat object;
- compute **AUCell / gene-set scores** and visualise them;
- run **cell-cell communication** with SCOP's CellChat wrapper;
- infer **transcription-factor regulons** with SCENIC using only R/C++;
- make **per-gene expression figures** (violin, UMAP) and patient-paired  
  pseudobulk statistics for a small gene panel;
- avoid the common SCOP traps: blank network PDFs, `recordedplot`, arrow/psock  
  crashes, RSS that is always 1, fixed-size heatmaps that clip rows.

## Install

```r
# SCOP
pak::pkg_install("jinworks/SCOP")

# Optional but recommended
pak::pkg_install(c("Seurat", "ComplexHeatmap", "ggplot2", "patchwork", "DESeq2"))
# CellChat backend for CCC
pak::pkg_install("jinworks/CellChat")
```

If SCOP lives in a private library (e.g. `Rlibs/`), set:

```bash
export SCOP_LIBS="/absolute/path/to/Rlibs"
```

The workflows add `SCOP_LIBS` to `.libPaths()` automatically.

## Quick start

```bash
Rscript examples/make_demo_data.R      # synthetic 720-cell Seurat object
Rscript examples/quickstart.R          # SCOP per-gene figures + heatmap + pseudobulk
```

Outputs: `results/quickstart/` with one figure per gene (`individual/*.pdf/png`),  
a `GroupHeatmap`, `pseudobulk_DESeq2.csv`, and `gene_summary.csv`.

## Workflows

### 1. Portable setup

`scripts/scop_helpers.R` provides:

| Helper                 | Purpose                                                     |
| ---------------------- | ----------------------------------------------------------- |
| `scop_skill_init()`    | seed, data/results directories, optional `libs`             |
| `scop_two_group_pal()` | muted two-group palette                                     |
| `scop_get_layer()`     | Seurat v4/v5 `GetAssayData` compatibility                   |
| `save_scop_plot()`     | save ggplot/patchwork/Heatmap to PDF + PNG; unwraps `$plot` |
| `save_scop_baseplot()` | save base/circlize plots (SCOP `CCCNetworkPlot` etc.)       |
| `pseudobulk_deseq2()`  | patient/sample-paired pseudobulk DESeq2                     |
| `scop_gene_summary()`  | per-gene mean expression / % positive by group              |

Always source the helper relative to the script:

```r
source(file.path(script_dir, "scop_helpers.R"))
```

### 2. Expression figures and statistics

`scripts/run_expression_workflow.R`:

```r
INPUT_RDS   <- "data/annotated.rds"
GENES       <- c("PFKFB3", "SLC2A4", "SLC2A3", "HK2", "PKM", "LDHA", "VEGFA")
GROUP_COL   <- "condition"        # e.g. prox/dist
GROUP_LEVELS <- c("prox", "dist")
SAMPLE_COL  <- "sample_id"        # unique patient x condition
COVARIATES  <- c("patient")
```

Produces one SCOP `FeatureStatPlot` violin and one `FeatureDimPlot` UMAP per gene,  
a `GroupHeatmap` for the panel, and a pseudobulk DESeq2 table.

Key SCOP calls:

```r
FeatureStatPlot(srt, stat.by = gene, group.by = "condition", plot_type = "violin",
                comparisons = list(c("prox", "dist")),
                add_stat = "none", sig_label = "p.format")

FeatureDimPlot(srt, features = gene, reduction = "umap")

GroupHeatmap(srt, features = genes, group.by = "condition",
             exp_method = "zscore", limits = c(-2, 2))
```

### 3. Marker / gene-set heatmaps

`scripts/run_heatmap_workflow.R` wraps `GroupHeatmap`:

```r
GroupHeatmap(srt, features = marker_genes, group.by = "celltype",
             split.by = NULL, exp_method = "zscore", limits = c(-2, 2),
             width = 5, height = 6)
```

For gene-set scores, add `CellScoring(srt, method = "AUCell", features = gene_sets,
name = "scores")` first, then use `GroupHeatmap(..., assay = "scores")` or  
`FeatureStatPlot(..., stat.by = "scores")`.

### 4. Cell-cell communication (CellChat)

`scripts/run_ccc_workflow.R`:

```r
srt <- RunCellChat(srt, group.by = "celltype", species = "Homo_sapiens",
                   min.cells = 10, do.fast = FALSE)

save_scop_baseplot(CCCNetworkPlot(srt, plot_type = "circle"), "ccc/01_circle")
save_scop_plot(CCCHeatmap(srt, plot_type = "heatmap"), "ccc/02_heatmap")
save_scop_plot(CCCHeatmap(srt, plot_type = "dot",
                          display_by = "interaction", top_n = 25),
               "ccc/03_dot")
```

For two conditions:

```r
RunCellChat(srt, group.by = "celltype", group_column = "condition",
            group_cmp = c("prox", "dist"), species = "Homo_sapiens")
```

**Never** `print()` a `CCCNetworkPlot` object: it is a `recordedplot` and gives  
a blank page. Use `save_scop_baseplot()` (or open the device first and evaluate  
the call inside it).

### 5. SCENIC / TF regulons (no Python)

`scripts/run_scenic_workflow.R` uses SCOP's `backend = "cpp"` — pure R/C++,  
**no pySCENIC/conda required**:

```r
srt <- RunSCENIC(
  srt, assay = "RNA", layer = "counts",
  regulators = tfs, targets = targets,
  ranking_dbs = "data/cistarget/hg38_...feather",
  motif_annotations = "data/cistarget/motifs-...tbl",
  species = "Homo_sapiens", genome = "hg38",
  backend = "cpp", cores = 4, parallel_backend = "fork",
  min_regulon_size = 10, max_regulon_targets = 50
)
```

Then:

```r
SCENICPlot(srt, group.by = "celltype", plot_type = "rss_rank")$plot
SCENICPlot(srt, group.by = "celltype", plot_type = "rss_heatmap")$rss_matrix
GroupHeatmap(srt, features = rownames(srt@tools$SCENIC$scores),
             group.by = "celltype", assay = "scenic", layer = "data")
```

Optional: `RunDorothea()` for TF activity (needs `dorothea` + `decoupleR`),  
`RunCisTarget(backend = "r")` for R-only motif enrichment.

## Function map

`references/functions.md` lists the SCOP functions used here with key arguments.

## Pitfalls

Read `references/pitfalls.md` before debugging. The most common ones:

1. `CCCNetworkPlot` / `SCENICPlot` return **lists or recordedplot**, not always a  
   ggplot. Use the helpers or `$plot`.
2. Base/circlize plots + `ragg::agg_png` = blank PNG; use `grDevices::png`.
3. A bare expression passed to a helper is evaluated once (R promise); use  
   `substitute()` + `eval()` per device.
4. `cistarget2` with psock workers fails on arrow ALTREP feather data; use  
   `parallel_backend = "fork"` on macOS/Linux, `cores = 1` on Windows.
5. `RunSCENIC(backend = "cpp")` needs no Python; `backend = "python"` needs  
   `pyscenic==0.12.1` + `arboreto` + `ctxcore`.
6. SCENIC RSS must come from `SCENICPlot()`; do not hand-roll a per-regulon JS  
   formula (a scalar JS is always 0 → RSS always 1).
7. `GroupHeatmap` is fixed-size: device must be larger than body + annotations;  
   otherwise bottom rows/column labels are clipped.
8. `FeatureStatPlot` `add_stat` only accepts `"none"`, `"mean"`, `"median"`;  
   use `comparisons` + `pairwise_method` for tests.
9. Sample-level `RunDEtest` needs `sample_col`/`condition_col` and (for paired  
   edgeR) `design` + a named numeric `contrast`; low-count genes are filtered.
10. Cell-level Wilcoxon p-values are descriptive only — cells are not  
    independent replicates.

## Repository layout

```text
scop-skill/
├── SKILL.md
├── README.md
├── LICENSE
├── CITATION.cff
├── scripts/
│   ├── scop_helpers.R
│   ├── run_expression_workflow.R
│   ├── run_heatmap_workflow.R
│   ├── run_ccc_workflow.R
│   ├── run_scenic_workflow.R
│   ├── run_qc_workflow.R
│   ├── run_integration_workflow.R
│   ├── run_annotation_workflow.R
│   ├── run_de_enrichment_workflow.R
│   ├── run_scoring_composition_workflow.R
│   ├── run_trajectory_workflow.R
│   ├── run_spatial_workflow.R
│   ├── run_deconvolution_workflow.R
│   ├── run_atac_workflow.R
│   ├── run_perturbation_workflow.R
│   ├── run_io_workflow.R
│   ├── run_database_prep.R
│   └── build_function_catalog.R
├── examples/
│   ├── make_demo_data.R
│   ├── quickstart.R
│   └── demo_data.rds
├── references/
│   ├── function-catalog.md      # all 360 exports
│   ├── scop_functions.csv       # machine-readable catalog
│   ├── coverage-matrix.md       # domain -> workflow map
│   ├── functions.md
│   ├── pitfalls.md
│   ├── workflows.md
│   └── case-study-ec-expression.md
└── results/        (git-ignored)
```

## Data policy


The repository contains only a tiny synthetic demo object. Do not commit patient  
data, large `.rds` files, cisTarget databases (hundreds of MB), or generated  
figures. `data/` and `results/` are git-ignored except for `.gitkeep`.

## License

MIT — see `LICENSE`. SCOP and CellChat have their own licenses; cite them when  
using this skill.
