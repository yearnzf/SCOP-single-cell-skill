# =============================================================================
# run_de_enrichment_workflow.R  (template)
# DE testing + enrichment / GSEA / GSVA / metabolism.
# =============================================================================
INPUT_RDS    <- "data/annotated.rds"
OUT_DIR      <- "results/de_enrichment"
GROUP_BY     <- "celltype"
GROUP1       <- "Disease"
GROUP2       <- "Control"
DE_TEST      <- "wilcox"       # or "DESeq2" / "edgeR" for sample-level
GENE_SETS    <- "data/gene_sets.gmt"
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({library(scop); library(Seurat)})

srt <- readRDS(INPUT_RDS)
srt <- RunDEtest(srt, group.by = GROUP_BY, group1 = GROUP1, group2 = GROUP2,
                 test.use = DE_TEST, only.pos = FALSE)
save_scop_plot(DEtestPlot(srt, group.by = GROUP_BY, test.use = DE_TEST,
                          group_use = GROUP1, plot_type = "volcano"),
               file.path(OUT_DIR, "01_volcano"), 7, 6)
srt <- RunEnrichment(srt, ...)             # ORA
save_scop_plot(EnrichmentPlot(srt, ...), file.path(OUT_DIR, "02_ora"), 9, 6)
srt <- RunGSEA(srt, ...)
save_scop_plot(GSEAPlot(srt, ...), file.path(OUT_DIR, "03_gsea"), 9, 6)
srt <- RunGSVA(srt, ...)
save_scop_plot(GSVAPlot(srt, ...), file.path(OUT_DIR, "04_gsva"), 9, 6)
srt <- RunMetabolism(srt, ...)
save_scop_plot(MetabolismPlot(srt, ...), file.path(OUT_DIR, "05_metabolism"), 9, 6)
message("Done: ", OUT_DIR)
