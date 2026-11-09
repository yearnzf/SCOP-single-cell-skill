# =============================================================================
# run_scenic_workflow.R
# SCOP SCENIC with the pure R/C++ backend (no Python / pySCENIC needed).
#
# Edit CONFIG, then:
#   Rscript scripts/run_scenic_workflow.R
#
# Required inputs:
#   - Seurat object with a count layer
#   - cisTarget ranking database (feather) and motif annotation table (.tbl)
#     Download from https://resources.aertslab.org/cistarget/
# =============================================================================

INPUT_RDS          <- "data/annotated.rds"
ASSAY              <- "RNA"
LAYER              <- "counts"
OUT_DIR            <- "results/scenic"
SPECIES            <- "Homo_sapiens"      # "Mus_musculus"
GENOME             <- "hg38"              # "mm10" for mouse
RANKING_DBS        <- c("data/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS  <- "data/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
REGULATORS         <- NULL                # NULL -> all TFs in motif annotation
TARGETS            <- NULL                # NULL -> VariableFeatures(object)
N_HVG              <- 2000
CORES              <- 4
MIN_EXPR_CELLS     <- 3
MIN_REGULON_SIZE   <- 10
MAX_REGULON_TARGETS <- 50

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
  library(RcisTarget)
  library(ComplexHeatmap)
})

if (!file.exists(INPUT_RDS)) stop("INPUT_RDS not found: ", INPUT_RDS)
if (!all(file.exists(RANKING_DBS))) stop("Ranking DB not found: ", RANKING_DBS)
if (!file.exists(MOTIF_ANNOTATIONS)) stop("Motif annotation not found: ", MOTIF_ANNOTATIONS)

cfg <- scop_skill_init(data_dir = dirname(INPUT_RDS), results_dir = OUT_DIR)
srt <- readRDS(INPUT_RDS)
srt@meta.data$celltype <- factor(as.character(srt@meta.data$celltype))

# Regulators: all TFs present in the motif annotation and in the object.
if (is.null(REGULATORS)) {
  motif_annot <- RcisTarget::importAnnotations(MOTIF_ANNOTATIONS)
  REGULATORS <- intersect(unique(motif_annot$TF), rownames(srt))
}
# Targets: HVG by default, plus all regulators (so low-expression TFs are kept).
if (is.null(TARGETS)) {
  srt <- FindVariableFeatures(srt, nfeatures = N_HVG, verbose = FALSE)
  TARGETS <- union(head(VariableFeatures(srt), N_HVG), REGULATORS)
}
TARGETS <- intersect(TARGETS, rownames(srt))

srt <- RunSCENIC(
  srt, assay = ASSAY, layer = LAYER,
  regulators = REGULATORS, targets = TARGETS,
  ranking_dbs = RANKING_DBS, motif_annotations = MOTIF_ANNOTATIONS,
  work_dir = file.path(cfg$results_dir, "scenic_work"),
  prefix = "scenic", species = SPECIES, genome = GENOME,
  backend = "cpp", cores = CORES, parallel_backend = "fork",
  min_expr_cells = MIN_EXPR_CELLS,
  min_regulon_size = MIN_REGULON_SIZE,
  max_regulon_targets = MAX_REGULON_TARGETS,
  n_rounds = 5000, seed = cfg$seed, force = FALSE, verbose = TRUE
)

# ---- tables --------------------------------------------------------------
sc <- srt@tools$SCENIC
write.csv(sc$regulons, file.path(cfg$results_dir, "regulons.csv"), row.names = FALSE)
write.csv(data.frame(regulon = rownames(sc$scores),
                     mean_activity = Matrix::rowMeans(sc$scores)),
          file.path(cfg$results_dir, "regulon_activity_summary.csv"), row.names = FALSE)

# ---- SCOP SCENIC figures -------------------------------------------------
save_scop_plot(
  SCENICPlot(srt, group.by = "celltype", plot_type = "rss_rank",
             top_n = 12, ncol = 2),
  file.path(cfg$results_dir, "01_RSS_rank"), 12, 8
)
save_scop_plot(
  SCENICPlot(srt, group.by = "celltype", plot_type = "rss_heatmap",
             top_n = 30, return_data = TRUE),
  file.path(cfg$results_dir, "02_RSS_heatmap"), 9, 9
)
save_scop_plot(
  SCENICPlot(srt, group.by = "celltype", plot_type = "activity_heatmap",
             top_n = 30),
  file.path(cfg$results_dir, "03_activity_heatmap"), 10, 9
)
save_scop_plot(
  SCENICPlot(srt, group.by = "celltype", plot_type = "network_graph",
             features = head(rownames(sc$scores), 1),
             max_targets = 20, label_nodes = "tfs", ncol = 1),
  file.path(cfg$results_dir, "04_regulon_network"), 9, 7
)

h <- GroupHeatmap(
  srt, features = head(rownames(sc$scores), 25),
  group.by = "celltype", assay = "scenic", layer = "data",
  exp_method = "zscore", limits = c(-2, 2),
  show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
  nlabel = 0, use_raster = FALSE, width = 5, height = 7
)
save_scop_plot(h, file.path(cfg$results_dir, "05_regulon_activity_heatmap"), 9, 9)

message("Done: ", cfg$results_dir)
