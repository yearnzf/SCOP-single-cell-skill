# SCOP function map used by this skill

All functions below are from **scop 0.9.x** unless noted. Signatures change
between versions; check `args()` in your installed version.

## Setup / utilities

| Function | Purpose | Notes |
|---|---|---|
| `scop::PrepareEnv()` | create/configure Python env | only needed for Python backends (CellChat is R; SCENIC `cpp` needs no Python) |
| `scop::ListEnv()` | list Python envs | |
| `scop::check_python()` | check/install Python packages | used internally by Python backends |
| `scop::AnnotateFeatures()` | add TF / surfaceome / enzyme annotations | `db = "TF"`, `db = "CSPA"`, etc. |
| `scop::CellDimPlot()` | UMAP/TSNE coloured by metadata | |
| `scop::FeatureDimPlot()` | UMAP/TSNE coloured by gene/feature | one feature per call is easiest |
| `scop::CellStatPlot()` | cell-count/percent composition | `stat_type = "percent"`, `position = "stack"` |
| `scop::FeatureStatPlot()` | violin/box/bar of a feature by group | `comparisons`, `add_stat` (`none/mean/median`), `sig_label` |
| `scop::GroupHeatmap()` | grouped heatmap (celltype × condition etc.) | returns list with `$plot`; fixed-size |
| `scop::FeatureHeatmap()` | feature × group heatmap | aliases/variants |
| `scop::CellScoring()` | AUCell / gene-set scoring | `method = "AUCell"`, `features = list(...)`, `name = "scores"` |
| `scop::CellScoringPlot()` | visualise scores | |
| `scop::RunDEtest()` | DE testing (cell-level or sample-level) | `test.use = "wilcox"` / `"DESeq2"` / `"edgeR"` / `"limma"` |
| `scop::DEtestPlot()` | volcano / Manhattan / ring | |
| `scop::RunCellQC()` | per-cell QC metrics | |
| `scop::RunStandardWorkflow()` | normalise / HVG / scale / PCA / UMAP | check args in your version |
| `scop::RunDimsEstimate()` | estimate number of PCs | returns integer vector |
| `scop::RunUMAP2()` | UMAP wrapper | |

## Cell–cell communication (CellChat backend)

| Function | Purpose | Notes |
|---|---|---|
| `scop::RunCellChat()` | run CellChat via SCOP | `group.by`, `species`, `min.cells`, `do.fast`, `group_column`, `group_cmp` |
| `scop::CCCNetworkPlot()` | circle / chord / diff_network / embedding_network | returns `recordedplot`; **do not print** |
| `scop::CCCHeatmap()` | heatmap / dot / focused_heatmap / diff_heatmap | `dot` needs `display_by = "interaction"` |
| `scop::CCCStatPlot()` | bar / pathway_summary / role_scatter / comparison | |
| `scop::PrepareCCCDB()` | prepare custom CCC database | |

## SCENIC / TF

| Function | Purpose | Notes |
|---|---|---|
| `scop::RunSCENIC()` | GRNBoost2 → cisTarget → AUCell | `backend = "cpp"` is pure R/C++, no Python |
| `scop::RunCisTarget()` | motif enrichment only | `backend = "r"` uses RcisTarget |
| `scop::RunGENIE3()` | GENIE3 co-expression | R-only, slow |
| `scop::RunGRNBoost2()` | arboreto GRNBoost2 | Python backend |
| `scop::RunGRN()` | unified GRN interface | |
| `scop::SCENICPlot()` | 16 plot types (RSS, activity, network, overlap, …) | returns list; use `$plot`, `$rss_matrix`, `$rank_table` |
| `scop::RunDorothea()` | DoRothEA/CollecTRI TF activity | needs `dorothea` + `decoupleR` |
| `scop::DorotheaPlot()` | TF activity plots | |
| `scop::RunSCENICPlus()` | SCENIC+ (RNA+ATAC) | needs Python + multiome |

## Key argument notes

- `FeatureStatPlot(add_stat = ...)` only accepts `"none"`, `"mean"`, `"median"`.
  Statistical comparisons are controlled by `comparisons`, `pairwise_method`,
  `multiplegroup_comparisons`, `multiple_method`, `sig_label`.
- `GroupHeatmap(width, height)` and `SCENICPlot` / `CCCHeatmap` bodies are
  **fixed-size**; the graphics device must be larger than the body plus
  annotations/legends.
- `RunDEtest` sample-level methods need `sample_col` and `condition_col`; paired
  edgeR additionally needs `design` and a **named numeric** `contrast` vector.
- `RunCellChat(do.fast = TRUE)` requires `presto`; use `FALSE` otherwise.
