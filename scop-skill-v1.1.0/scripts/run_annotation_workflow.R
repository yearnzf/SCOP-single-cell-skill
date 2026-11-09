# =============================================================================
# run_annotation_workflow.R  (template)
# Cell-type annotation / reference mapping.
# Some methods need reference data or Python (CellTypist).
# =============================================================================
INPUT_RDS   <- "data/processed.rds"
OUT_DIR     <- "results/annotation"
REFERENCE   <- NULL     # Seurat object or SingleR reference
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
# Choose one or more:
srt <- RunSingleR(srt, ...)                # reference-based
srt <- RunCellTypist(srt, ...)             # needs Python + CellTypist models
srt <- RunSciBet(srt, ...)                 # needs reference
srt <- RunScmap(srt, ...)                  # needs reference
srt <- RunKNNMap(srt, ...)                 # label transfer
srt <- RunLabelTransfer(srt, ...)          # Seurat anchors
save_scop_plot(CellDimPlot(srt, group.by = "celltype"),
               file.path(OUT_DIR, "01_annotation"), 7, 6)
save_scop_plot(CellStatPlot(srt, stat.by = "celltype", group.by = "cluster",
                            stat_type = "percent", position = "stack"),
               file.path(OUT_DIR, "02_annotation_by_cluster"), 9, 5)
message("Done: ", OUT_DIR)
