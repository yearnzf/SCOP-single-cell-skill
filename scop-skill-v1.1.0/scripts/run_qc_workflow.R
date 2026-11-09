# =============================================================================
# run_qc_workflow.R  (template)
# QC -> normalisation -> HVG -> PCA/UMAP -> clustering.
# Check `args(scop::RunCellQC)` / `args(scop::RunStandardWorkflow)` in your
# SCOP version; argument names have changed across releases.
# =============================================================================
INPUT_RDS <- "data/raw_counts.rds"
OUT_DIR   <- "results/qc"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
srt <- RunCellQC(srt, ...)                 # edit QC thresholds
srt <- RunStandardWorkflow(srt, ...)       # normalise / HVG / scale / PCA / UMAP
srt <- RunDimsEstimate(srt, ...)
srt <- RunUMAP2(srt, dims = 1:20)
srt <- FindNeighbors(srt, dims = 1:20)
srt <- FindClusters(srt, resolution = 0.5)

save_scop_plot(CellDimPlot(srt, group.by = "seurat_clusters"),
               file.path(OUT_DIR, "01_clusters"), 7, 6)
save_scop_plot(CellStatPlot(srt, stat.by = "seurat_clusters",
                            group.by = "celltype", stat_type = "percent",
                            position = "stack"),
               file.path(OUT_DIR, "02_composition"), 8, 5)
saveRDS(srt, file.path(OUT_DIR, "qc_processed.rds"))
message("Done: ", OUT_DIR)
