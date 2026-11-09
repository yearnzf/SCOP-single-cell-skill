# =============================================================================
# run_atac_workflow.R  (template)
# ATAC / multiome: QC, gene activity, chromVAR, peak-gene links, motifs.
# =============================================================================
INPUT_RDS  <- "data/atac.rds"
OUT_DIR    <- "results/atac"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
# srt <- ReadATACData(...)            # if starting from fragments/peaks
srt <- RunATACQC(srt, ...)
srt <- RunGeneActivity(srt, ...)
srt <- RunChromVAR(srt, ...)
srt <- RunLinkPeaks(srt, ...)
srt <- RunMotifEnrichment(srt, ...)
save_scop_plot(CellDimPlot(srt, group.by = "celltype"),
               file.path(OUT_DIR, "01_clusters"), 7, 6)
save_scop_plot(CoverageTrackPlot(srt, ...), file.path(OUT_DIR, "02_coverage"), 9, 6)
message("Done: ", OUT_DIR)
