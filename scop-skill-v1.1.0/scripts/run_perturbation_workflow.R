# =============================================================================
# run_perturbation_workflow.R  (template)
# Perturbation / regulatory-network / priority methods.
# =============================================================================
INPUT_RDS  <- "data/processed.rds"
OUT_DIR    <- "results/perturbation"
GROUP_BY   <- "celltype"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
# Choose:
srt <- RunscTenifoldKnk(srt, ...)          # virtual knockout
save_scop_plot(scTenifoldKnkPlot(srt), file.path(OUT_DIR, "01_tenifold_knk"), 9, 6)
srt <- RunscTenifoldNet(srt, ...)
save_scop_plot(scTenifoldNetPlot(srt), file.path(OUT_DIR, "02_tenifold_net"), 9, 6)
srt <- RunScissor(srt, ...)                # phenotype priority
save_scop_plot(ScissorPlot(srt), file.path(OUT_DIR, "03_scissor"), 8, 6)
srt <- RunAugur(srt, ...)                  # cell-type prioritisation
srt <- RunscPagwas(srt, ...)
save_scop_plot(PlotscPagwas(srt), file.path(OUT_DIR, "04_pagwas"), 9, 6)
srt <- RunVECTOR(srt, ...)
save_scop_plot(VECTORPlot(srt), file.path(OUT_DIR, "05_vector"), 8, 6)
message("Done: ", OUT_DIR)
