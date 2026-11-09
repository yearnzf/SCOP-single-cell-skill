# =============================================================================
# run_io_workflow.R  (template)
# Convert between Seurat, AnnData/h5ad, loom, Giotto, SpatialExperiment.
# =============================================================================
INPUT_RDS <- "data/annotated.rds"
OUT_DIR   <- "results/io"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

srt_to_h5ad(srt, file.path(OUT_DIR, "object.h5ad"))
srt_to_adata(srt, file.path(OUT_DIR, "object_adata"))
# srt_to_loom(srt, ...)   # if available in your SCOP version
# Readers:
# h5ad_to_srt("object.h5ad")
# loom_to_srt("object.loom")
# giotto_to_srt("giotto_object")
# spe_to_srt("spatial_experiment.rds")
message("Done: ", OUT_DIR)
