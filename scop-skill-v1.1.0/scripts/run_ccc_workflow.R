# =============================================================================
# run_ccc_workflow.R
# SCOP cell-cell communication (CellChat backend).
#
# Edit CONFIG, then:
#   Rscript scripts/run_ccc_workflow.R
#
# Requires: CellChat (>= 2.0), SCOP >= 0.9.
# For a two-condition comparison, set GROUP_COL="condition" and use
# `RunCellChat(..., group_column="condition", group_cmp=c("prox","dist"))`
# (see comments below).
# =============================================================================

INPUT_RDS  <- "data/annotated.rds"
GROUP_COL  <- "celltype"
SPECIES    <- "Homo_sapiens"       # "Mus_musculus" / "zebrafish"
OUT_DIR    <- "results/ccc"
MIN_CELLS  <- 10
DO_FAST    <- FALSE                # TRUE needs presto
# CONDITION_COL <- "condition"     # uncomment for dist vs prox comparison
# CONDITION_PAIR <- c("prox", "dist")

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
})
source(file.path(.script_dir, "scop_helpers.R"))

lib_env <- Sys.getenv("SCOP_LIBS", unset = "")
if (nzchar(lib_env)) {
  .libPaths(c(strsplit(lib_env, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
}
suppressPackageStartupMessages({
  library(scop)
  library(Seurat)
})

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("Install CellChat first, e.g. pak::pkg_install('jinworks/CellChat')")
}
if (!file.exists(INPUT_RDS)) stop("INPUT_RDS not found: ", INPUT_RDS)

cfg <- scop_skill_init(data_dir = dirname(INPUT_RDS), results_dir = OUT_DIR)
srt <- readRDS(INPUT_RDS)

srt <- RunCellChat(
  srt, group.by = GROUP_COL, species = SPECIES,
  min.cells = MIN_CELLS, do.fast = DO_FAST, verbose = TRUE
)
# For two-condition comparisons, replace the call above with:
# srt <- RunCellChat(srt, group.by = GROUP_COL, group_column = CONDITION_COL,
#                    group_cmp = CONDITION_PAIR, species = SPECIES,
#                    min.cells = MIN_CELLS, do.fast = DO_FAST)

# ---- network plots (base/circlize; must draw while a device is open) -------
save_scop_baseplot(
  CCCNetworkPlot(srt, method = "CellChat", plot_type = "circle"),
  file.path(cfg$results_dir, "01_network_circle"), 7, 7
)
save_scop_baseplot(
  CCCNetworkPlot(srt, method = "CellChat", plot_type = "chord"),
  file.path(cfg$results_dir, "02_network_chord"), 7, 7
)
save_scop_baseplot(
  CCCNetworkPlot(srt, method = "CellChat", plot_type = "embedding_network",
                 reduction = "umap"),
  file.path(cfg$results_dir, "03_embedding_network"), 8, 7
)

# ---- heatmaps / dot plots ------------------------------------------------
save_scop_plot(
  CCCHeatmap(srt, method = "CellChat", plot_type = "heatmap"),
  file.path(cfg$results_dir, "04_heatmap_sender_receiver"), 8, 7
)
save_scop_plot(
  CCCHeatmap(srt, method = "CellChat", plot_type = "dot",
             display_by = "interaction", top_n = 25),
  file.path(cfg$results_dir, "05_dot_interaction"), 12, 8
)
save_scop_plot(
  CCCStatPlot(srt, method = "CellChat", plot_type = "pathway_summary"),
  file.path(cfg$results_dir, "06_pathway_summary"), 10, 6
)

# ---- tables --------------------------------------------------------------
cc <- srt@tools$CellChat
if (!is.null(cc$long_table)) {
  write.csv(cc$long_table, file.path(cfg$results_dir, "CellChat_long_table.csv"),
            row.names = FALSE)
}
if (!is.null(cc$pair_table)) {
  write.csv(cc$pair_table, file.path(cfg$results_dir, "CellChat_pair_table.csv"),
            row.names = FALSE)
}
if (!is.null(cc$pathway_table)) {
  write.csv(cc$pathway_table, file.path(cfg$results_dir, "CellChat_pathway_table.csv"),
            row.names = FALSE)
}
message("Done: ", cfg$results_dir)
