# =============================================================================
# run_trajectory_workflow.R  (template)
# Trajectory / pseudotime / velocity. Pick one method.
# =============================================================================
INPUT_RDS  <- "data/processed.rds"
OUT_DIR    <- "results/trajectory"
GROUP_BY   <- "celltype"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
# Choose one:
srt <- RunSlingshot(srt, group.by = GROUP_BY, ...)
srt <- RunMonocle3(srt, ...)
srt <- RunPAGA(srt, ...)
srt <- RunPalantir(srt, ...)
srt <- RunCellRank(srt, ...)
srt <- RunSCVELO(srt, ...)                 # needs spliced/unspliced layers + Python
save_scop_plot(LineagePlot(srt), file.path(OUT_DIR, "01_lineage"), 8, 6)
save_scop_plot(PseudotimeProjectionPlot(srt), file.path(OUT_DIR, "02_pseudotime"), 8, 6)
srt <- RunDynamicFeatures(srt, ...)
save_scop_plot(DynamicPlot(srt, ...), file.path(OUT_DIR, "03_dynamic_genes"), 9, 6)
save_scop_plot(DynamicHeatmap(srt, ...), file.path(OUT_DIR, "04_dynamic_heatmap"), 9, 8)
message("Done: ", OUT_DIR)
