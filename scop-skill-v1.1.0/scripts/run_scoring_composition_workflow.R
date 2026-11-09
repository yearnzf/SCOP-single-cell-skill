# =============================================================================
# run_scoring_composition_workflow.R  (template)
# Gene-set scoring (AUCell) + cell composition tests.
# =============================================================================
INPUT_RDS   <- "data/annotated.rds"
OUT_DIR     <- "results/scoring_composition"
GENE_SETS   <- list(Glycolysis = c("PFKFB3", "HK2", "PKM", "LDHA"),
                    Angiogenesis = c("VEGFA", "KDR", "ANGPT2"))
GROUP_BY    <- "condition"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
srt <- CellScoring(srt, method = "AUCell", features = GENE_SETS, name = "scores")
save_scop_plot(CellScoringPlot(srt, ...), file.path(OUT_DIR, "01_scores"), 9, 6)
save_scop_plot(CellStatPlot(srt, stat.by = "celltype", group.by = GROUP_BY,
                            stat_type = "percent", position = "stack"),
               file.path(OUT_DIR, "02_composition"), 8, 5)
srt <- RunProportionTest(srt, ...)
save_scop_plot(ProportionTestPlot(srt), file.path(OUT_DIR, "03_proportion_test"), 8, 5)
srt <- RunscCODA(srt, ...)
srt <- RunPropeller(srt, ...)
message("Done: ", OUT_DIR)
