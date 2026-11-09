# =============================================================================
# run_deconvolution_workflow.R  (template)
# Bulk / spatial deconvolution with SCOP wrappers.
# =============================================================================
INPUT_RDS    <- "data/spatial.rds"    # spatial or bulk object
REFERENCE    <- "data/reference.rds"  # single-cell reference
OUT_DIR      <- "results/deconvolution"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
ref <- readRDS(REFERENCE)
# Choose one or more:
srt <- RunCIBERSORT(srt, ref, ...)
srt <- RunCARD(srt, ...)
srt <- RunRCTD(srt, ...)
srt <- RunSPOTlight(srt, ...)
srt <- RunCell2location(srt, ...)       # needs Python backend
srt <- RunSTdeconvolve(srt, ...)
save_scop_plot(DeconvolutionPlot(srt), file.path(OUT_DIR, "01_deconvolution"), 9, 6)
save_scop_plot(SpatialDeconvolutionPlot(srt), file.path(OUT_DIR, "02_spatial"), 8, 7)
message("Done: ", OUT_DIR)
