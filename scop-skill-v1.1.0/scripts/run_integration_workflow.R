# =============================================================================
# run_integration_workflow.R  (template)
# Batch integration and benchmark. Requires a list of Seurat objects.
# =============================================================================
INPUT_RDS_LIST <- c("data/sample1.rds", "data/sample2.rds")
OUT_DIR        <- "results/integration"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt_list <- lapply(INPUT_RDS_LIST, readRDS)
# integration_scop() orchestrates SCOP integration methods; see ?integration_scop
srt <- integration_scop(srt_list, ...)     # e.g. method = "Harmony"
srt <- RunIntegration(srt, ...)            # or run one method directly
srt <- RunLISI(srt, ...)                   # integration quality
srt <- RunIntegrationBenchmark(srt_list, ...)
save_scop_plot(IntegrationBenchmarkPlot(srt), file.path(OUT_DIR, "01_benchmark"), 9, 6)
save_scop_plot(CellDimPlot(srt, group.by = "batch"), file.path(OUT_DIR, "02_batch"), 7, 6)
message("Done: ", OUT_DIR)
