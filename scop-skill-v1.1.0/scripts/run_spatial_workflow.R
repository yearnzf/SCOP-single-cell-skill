# =============================================================================
# run_spatial_workflow.R  (template)
# Spatial transcriptomics: QC, variable features, neighbourhood, gradient,
# network, integration and plots.
# =============================================================================
INPUT_RDS  <- "data/spatial.rds"      # object from ReadSpatialData()
OUT_DIR    <- "results/spatial"
GROUP_BY   <- "celltype"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
# srt <- ReadSpatialData(...)         # if starting from raw spatial data
srt <- RunSpotQC(srt, ...)
srt <- RunSpatialVariableFeatures(srt, ...)
srt <- RunSpatialNeighborhood(srt, ...)
srt <- RunSpatialGradientFeatures(srt, ...)
srt <- RunSpatialNetwork(srt, ...)
srt <- RunSpatialIntegration(srt, ...)  # if multiple slices
save_scop_plot(SpatialCellPlot(srt, group.by = GROUP_BY),
               file.path(OUT_DIR, "01_spatial_cells"), 8, 7)
save_scop_plot(SpatialNeighborhoodPlot(srt), file.path(OUT_DIR, "02_neighborhood"), 8, 7)
save_scop_plot(SpatialGradientPlot(srt), file.path(OUT_DIR, "03_gradient"), 8, 7)
save_scop_plot(SpatialNetworkPlot(srt), file.path(OUT_DIR, "04_network"), 8, 7)
message("Done: ", OUT_DIR)
