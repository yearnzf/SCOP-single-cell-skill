# =============================================================================
# run_heatmap_workflow.R
# SCOP marker / gene-set heatmaps (GroupHeatmap + optional FeatureHeatmap).
#
# Edit CONFIG, then:
#   Rscript scripts/run_heatmap_workflow.R
# =============================================================================

INPUT_RDS   <- "data/annotated.rds"
OUT_DIR     <- "results/heatmaps"
GROUP_BY    <- "celltype"
SPLIT_BY    <- NULL                 # e.g. "condition" for celltype x condition
FEATURES    <- NULL                 # NULL -> use a character vector below
FEATURE_VEC <- c("PECAM1", "VWF", "KDR", "COL1A1", "DCN", "ACTA2", "MYH11")
EXP_METHOD  <- "zscore"             # "zscore" / "none" / "scale"
LIMITS      <- c(-2, 2)
CLUSTER_ROWS <- FALSE
CLUSTER_COLS <- FALSE
WIDTH_BODY  <- 5
HEIGHT_BODY <- 6

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
  library(ComplexHeatmap)
})

if (!file.exists(INPUT_RDS)) stop("INPUT_RDS not found: ", INPUT_RDS)
cfg <- scop_skill_init(data_dir = dirname(INPUT_RDS), results_dir = OUT_DIR)
srt <- readRDS(INPUT_RDS)
if (is.null(FEATURES)) FEATURES <- FEATURE_VEC
FEATURES <- intersect(FEATURES, rownames(srt))
if (!length(FEATURES)) stop("No requested features found in the object.")

h <- GroupHeatmap(
  srt, features = FEATURES, group.by = GROUP_BY, split.by = SPLIT_BY,
  exp_method = EXP_METHOD, limits = LIMITS,
  heatmap_palcolor = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  show_row_names = TRUE, row_names_side = "left",
  show_column_names = TRUE, column_names_rot = 45,
  cluster_rows = CLUSTER_ROWS, cluster_columns = CLUSTER_COLS,
  nlabel = 0, use_raster = FALSE, seed = cfg$seed,
  width = WIDTH_BODY, height = HEIGHT_BODY
)
save_scop_plot(h, file.path(cfg$results_dir, "01_GroupHeatmap"), 10, 9)
message("Done: ", cfg$results_dir)
