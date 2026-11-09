# SCOP coverage matrix

Every exported function of the installed SCOP (0.9.2, **360 functions**) is listed
in `references/function-catalog.md` / `references/scop_functions.csv`.
The table below maps each functional domain to a workflow, its status, and the
data/backends it needs.

Legend:

- **runnable demo** — runs end-to-end on the synthetic demo (`examples/`).
- **template** — runnable once you supply the input data/backends; signatures
  vary by SCOP version, so check `args()`.
- **recipe** — documented snippet for a set of related functions.

| Domain | SCOP functions | Workflow / file | Status | Requires |
|---|---|---|---|---|
| Environment & database prep | 40 | `scripts/run_database_prep.R` | template | network for `Prepare*` downloads |
| Data I/O & conversion | 16 | `scripts/run_io_workflow.R` | template | h5ad/loom/Giotto/SpatialExperiment files |
| QC & preprocessing | 40 | `scripts/run_qc_workflow.R` | template | raw counts |
| Integration & batch correction | 30 | `scripts/run_integration_workflow.R` | template | multiple samples |
| Cell annotation & reference mapping | 22 | `scripts/run_annotation_workflow.R` | template | reference object; CellTypist needs Python |
| Differential expression & markers | 10 | `scripts/run_de_enrichment_workflow.R` | template | annotated object |
| Enrichment / GSEA / GSVA / metabolism | 18 | `scripts/run_de_enrichment_workflow.R` | template | gene sets / SCOP DBs |
| Gene scoring & module scores | 5 | `scripts/run_scoring_composition_workflow.R` | template | gene sets |
| Cell composition | 7 | `scripts/run_scoring_composition_workflow.R` | template | cell-type labels |
| Trajectory, pseudotime & velocity | 26 | `scripts/run_trajectory_workflow.R` | template | UMAP/clusters; velocity needs spliced/unspliced + Python |
| Cell–cell communication | 23 | `scripts/run_ccc_workflow.R` | template | CellChat DB; CellphoneDB/LIANA need Python |
| SCENIC & transcription factors | 17 | `scripts/run_scenic_workflow.R` | template | cisTarget DB (~300–400 MB) |
| Spatial transcriptomics & deconvolution | 43 | `scripts/run_spatial_workflow.R`, `scripts/run_deconvolution_workflow.R` | template | spatial object / reference |
| ATAC & multiome | 13 | `scripts/run_atac_workflow.R` | template | ATAC fragments/peak matrix |
| Perturbation & regulatory networks | 11 | `scripts/run_perturbation_workflow.R` | template | processed object; some Python |
| Visualisation & themes | 26 | `scripts/scop_helpers.R` + all workflows | recipe | any Seurat object |
| Core utilities | 13 | `scripts/scop_helpers.R` | recipe | any Seurat object |
| Per-gene expression + pseudobulk | — | `scripts/run_expression_workflow.R` | **runnable demo** | Seurat object with UMAP |
| Marker / gene-set heatmaps | — | `scripts/run_heatmap_workflow.R` | **runnable demo** | Seurat object |
| End-to-end smoke test | — | `examples/quickstart.R` | **runnable demo** | none (synthetic data) |

## Notes by domain

### QC / preprocessing
`RunCellQC`, `RunStandardWorkflow`, `RunDimsEstimate`, `RunPCA`, `RunUMAP2`,
`FindNeighbors`, `FindClusters`, `RunDoublet*`, `RunDecontX`, `RunCellCycle*`,
`RunNMF`, `RunBANKSY`, `RunSmoothClust`. Signatures are version-sensitive; the
template marks the editable thresholds.

### Integration
`integration_scop()` coordinates multiple methods; individual wrappers include
`Harmony_integrate`, `CCA_integrate`, `RPCA_integrate`, `scVI_integrate`,
`fastMNN_integrate`, `BBKNN_integrate`, `Scanorama_integrate`, `Seurat_integrate`,
`LIGER_integrate`, `Conos_integrate`, `CSS_integrate`, `Coralysis_integrate`,
`GLUE_integrate`, `MultiMAP_integrate`, `scVI5_integrate`, `Harmony5_integrate`.
Some require Python/conda (`scVI`, `GLUE`, `Conos`, `LIGER`). Benchmark with
`RunIntegrationBenchmark()` and `RunLISI()`.

### Annotation
Reference-based: `RunSingleR`, `RunScmap`, `RunSciBet`, `RunKNNMap`,
`RunKNNPredict`, `RunSeuratMap`, `RunSymphonyMap`, `RunLabelTransfer`,
`RunReferenceMapping`. Python model: `RunCellTypist`/`TrainCellTypist`.
Malignant / CNV: `RunscMalignant*`, `RunCNV`/`CNVPlot`, `RunESTIMATE`.

### DE / enrichment
`RunDEtest()` supports cell-level (`wilcox`, `t`, `negbinom`, `poisson`, `LR`)
and sample-level (`DESeq2`, `edgeR`, `limma`, `dream`) tests; `DEtestPlot()`
draws volcano/Manhattan/ring. Enrichment: `RunEnrichment`, `RunGSEA`, `RunGSVA`,
`RunIREA`, `RunDynamicEnrichment`, `RunMetabolism`, `RunscFEA`.

### Trajectory / velocity
`RunSlingshot` + `LineagePlot` / `PseudotimeProjectionPlot`;
`RunMonocle2/3`, `RunPAGA`/`PAGAPlot`, `RunPalantir`,
`RunCellRank`/`RunCellRankTrends`/`RunCellRankEnrichment`,
`RunSCVELO`/`VelocityPlot`, `RunWOT`, `RunPHATE`, `RunVECTOR`.
`RunDynamicFeatures` + `DynamicPlot`/`DynamicHeatmap` for gene dynamics.

### CCC
R backends: `RunCellChat`, `RunCCC`, `RunSpatialCellChat`, `RunCOMMOT`,
`RunMDIC3`, `RunSpaTalk`. Python backends: `RunCellphoneDB`, `RunLIANA`,
`RunNichenetr`, `RunMultiNichenetr`, `RunSecAct*`.
Plots: `CCCNetworkPlot`, `CCCHeatmap`, `CCCStatPlot`.

### SCENIC / TF
`RunSCENIC(backend="cpp")` (no Python), `RunSCENICPlus` (multiome),
`RunCisTarget(backend="r")`, `RunGENIE3`, `RunGRNBoost2`, `RunGRN`,
`RunDorothea` (TF activity), `SCENICPlot`, `SCENICPlusPlot`, `DorotheaPlot`.

### Spatial
`ReadSpatialData`; `RunSpatialVariableFeatures`, `RunSpatialNeighborhood`,
`RunSpatialGradientFeatures`, `RunSpatialNetwork`, `RunSpatialIntegration`,
`RunSpatialBenchmark`, `RunSpatialEcoTyper`, `RunStatialKontextual`, `RunMistyR`.
Deconvolution: `RunCIBERSORT`, `RunCARD`, `RunRCTD`, `RunSPOTlight`,
`RunCell2location`, `RunBayesSpace`, `RunSTdeconvolve`, `RunCytoSPACE`,
`RunBisqueRNA`, `RunBayesPrism`, `RunscOMM`, `RunCell2fate`, `RunCSIDE`,
`RunMERINGUE`.

### ATAC / multiome
`ReadATACData`, `RunATACQC`, `FilterATACPeaks`, `CombineATACPeaks`,
`AnnotateATACPeaks`, `AddATACMotifs`, `RunGeneActivity`, `RunChromVAR`,
`RunLinkPeaks`, `RunMotifEnrichment`, `CoverageTrackPlot`, `RecountATAC`,
`ExportATACResults`, `ReadMultiomeData`, `RunSCENICPlus`.

### Perturbation / prioritisation
`RunscTenifoldKnk` (virtual KO), `RunscTenifoldNet`, `RunScissor`,
`RunAugur`, `RunscPagwas`, `RunVECTOR`, `RunFitDevo`, `RunMistyR`.

### Visualisation
All plotting functions return ggplot/patchwork/Heatmap or base graphics; use
`save_scop_plot()` for ggplot-like outputs and `save_scop_baseplot()` for
base/circlize outputs. `theme_scop` / `theme_spatial` are the themes.

## Maintenance check

After upgrading SCOP, regenerate the catalog and verify coverage:

```bash
SCOP_LIBS=/path/to/Rlibs Rscript scripts/build_function_catalog.R
Rscript -e 'd <- read.csv("references/scop_functions.csv"); stopifnot(all(table(d$category) > 0)); cat(nrow(d), "functions catalogued\n")'
```
