# SCOP single-cell skill

A GitHub-ready, data-agnostic **SCOP** workflow collection for Seurat single-cell
analysis: marker/gene-set heatmaps, per-gene expression figures, patient-paired
pseudobulk statistics, cell–cell communication (CellChat), and SCENIC/regulon
inference with a **pure R/C++ backend (no Python)**.

This repository distils a real project's SCOP analyses into reusable scripts,
helpers, and a list of pitfalls that are easy to hit and hard to diagnose.

## Why this exists

SCOP is powerful but has sharp edges:

- network plots return `recordedplot` and silently produce blank PDFs if printed;
- base/circlize plots break with `ragg::agg_png`;
- SCENIC's C++ backend avoids Python, but psock workers fail on arrow feather data;
- `GroupHeatmap` is fixed-size and clips rows/labels if the device is too small;
- RSS values can be wrong if calculated by hand.

This skill packages the working patterns so the next project does not repeat
those debugging sessions.

## Install

```r
pak::pkg_install("jinworks/SCOP")
pak::pkg_install(c("Seurat", "ComplexHeatmap", "ggplot2", "patchwork", "DESeq2"))
# optional CellChat backend
pak::pkg_install("jinworks/CellChat")
```

If SCOP is in a private library:

```bash
export SCOP_LIBS="/absolute/path/to/Rlibs"
```

## Quick start

```bash
git clone <this-repo>
cd scop-skill
Rscript examples/make_demo_data.R
Rscript examples/quickstart.R
```

The demo creates 720 synthetic cells (3 cell types, 3 patients, prox/dist) and
writes per-gene SCOP violins/UMAPs, a `GroupHeatmap`, a pseudobulk DESeq2 table,
and a per-gene summary.

## Push to GitHub

This folder is already a git repository with one initial commit:

```bash
git remote add origin git@github.com:<your-account>/scop-skill.git
git branch -M main
git push -u origin main
```

Replace `<your-account>` and the placeholder `repository-code` in `CITATION.cff`.

## Workflows

| Script | What it does |
|---|---|
| `scripts/run_expression_workflow.R` | per-gene `FeatureStatPlot` + `FeatureDimPlot`, `GroupHeatmap`, optional paired pseudobulk DESeq2 |
| `scripts/run_heatmap_workflow.R` | marker / gene-set heatmaps with `GroupHeatmap` |
| `scripts/run_ccc_workflow.R` | `RunCellChat` + circle/chord/embedding/heatmap/dot plots |
| `scripts/run_scenic_workflow.R` | `RunSCENIC(backend="cpp")` + `SCENICPlot` + regulon heatmap |
| `scripts/run_qc_workflow.R` | QC / normalise / PCA-UMAP / clustering |
| `scripts/run_integration_workflow.R` | integration, LISI, benchmark |
| `scripts/run_annotation_workflow.R` | SingleR / CellTypist / scmap / label transfer |
| `scripts/run_de_enrichment_workflow.R` | `RunDEtest` + ORA / GSEA / GSVA / metabolism |
| `scripts/run_scoring_composition_workflow.R` | AUCell scoring + composition tests |
| `scripts/run_trajectory_workflow.R` | Slingshot / Monocle / PAGA / Palantir / CellRank / scVelo |
| `scripts/run_spatial_workflow.R` | spatial QC, neighbourhood, gradients, networks |
| `scripts/run_deconvolution_workflow.R` | CIBERSORT / CARD / RCTD / SPOTlight / Cell2location |
| `scripts/run_atac_workflow.R` | ATAC QC, gene activity, chromVAR, links, motifs |
| `scripts/run_perturbation_workflow.R` | scTenifold, Scissor, Augur, scPagwas, VECTOR |
| `scripts/run_io_workflow.R` | Seurat ↔ h5ad / AnnData / loom / Giotto / SPE |
| `scripts/run_database_prep.R` | `Prepare*` annotation/enrichment/CCC/TF databases |
| `scripts/scop_helpers.R` | portable setup, saving, pseudobulk, summaries |

Edit the `CONFIG` block at the top of each workflow and run with
`Rscript scripts/<workflow>.R`.

### Example: per-gene figures

```r
p <- FeatureStatPlot(
  srt, stat.by = gene, group.by = "condition", plot_type = "violin",
  comparisons = list(c("prox", "dist")),
  add_stat = "none", sig_label = "p.format"
)
save_scop_plot(p, "results/individual/VEGFA_violin")
```

### Example: SCENIC without Python

```r
srt <- RunSCENIC(
  srt, assay = "RNA", layer = "counts",
  regulators = tfs, targets = targets,
  ranking_dbs = "data/cistarget/hg38_...feather",
  motif_annotations = "data/cistarget/motifs-...tbl",
  backend = "cpp", cores = 4, parallel_backend = "fork"
)
SCENICPlot(srt, group.by = "celltype", plot_type = "rss_rank")$plot
```

See `SKILL.md`, `references/workflows.md`, and the worked example in
`references/case-study-ec-expression.md`.

## Repository structure

```text
SKILL.md                 skill definition (for agent skill loaders)
README.md                this file
scripts/                 portable workflows and helpers
examples/                synthetic demo + quickstart
references/              function catalog, coverage matrix, function map,
                          pitfalls, workflow recipes, case study
results/                 generated outputs (git-ignored)
data/                    user data / cisTarget DB (git-ignored)
```

## Full SCOP coverage

All **360 exported functions** of SCOP 0.9.2 are catalogued in
`references/function-catalog.md` (CSV: `references/scop_functions.csv`), and every
domain is mapped to a workflow in `references/coverage-matrix.md`:

QC, integration, annotation, DE, enrichment, scoring, composition, trajectory,
cell–cell communication, SCENIC/TF, spatial/deconvolution, ATAC/multiome,
perturbation, I/O, database preparation, visualisation, and per-gene expression.

The templates are complete scripts; replace the CONFIG block and check
`args(scop::FunctionName)` for your SCOP version. The expression and heatmap
workflows run on the included demo.

## Data and databases

The repository ships only a **synthetic** demo object. Do not commit patient data.

SCENIC needs cisTarget databases (~300–400 MB), download from
<https://resources.aertslab.org/cistarget/> and point `RANKING_DBS` /
`MOTIF_ANNOTATIONS` at them. CellChat needs its species database via
`CellChatDB.human` / `CellChatDB.mouse`.

## Pitfalls

A non-exhaustive list (details in `references/pitfalls.md`):

- `CCCNetworkPlot` → `recordedplot`: use `save_scop_baseplot()`.
- base graphics + `ragg` → blank PNG: use `grDevices::png()`.
- psock + arrow ALTREP feather → serialization error: `parallel_backend = "fork"`
  (or `cores = 1` on Windows).
- `RunSCENIC(backend="cpp")` requires **no** pySCENIC; don't install Python just
  for SCOP's C++ path.
- RSS: use `SCENICPlot(...)$rss_matrix`, not a custom JS formula.
- `GroupHeatmap`: device must be bigger than the fixed body + annotations.
- cell-level Wilcoxon p-values are descriptive; use pseudobulk for inference.

## Citation

If you use this skill, cite SCOP and the underlying methods:

- SCOP: <https://github.com/jinworks/SCOP>
- CellChat: Jin et al., *Nature Communications* 2021
- SCENIC: Aibar et al., *Nature Methods* 2017
- pySCENIC/arboreto: Van de Sande et al., *Nature Protocols* 2020

## License

MIT. See `LICENSE`.
