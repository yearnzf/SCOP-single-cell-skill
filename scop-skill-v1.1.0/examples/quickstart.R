# =============================================================================
# quickstart.R
# Minimal, runnable example: SCOP-native per-gene figures + heatmap + pseudobulk.
#   Rscript examples/make_demo_data.R
#   Rscript examples/quickstart.R
# =============================================================================
script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
})

# Optional local library path: set SCOP_LIBS=/path/to/Rlibs if SCOP is not in
# the default .libPaths().
lib_env <- Sys.getenv("SCOP_LIBS", unset = "")
lib_candidates <- c(
  if (nzchar(lib_env)) strsplit(lib_env, .Platform$path.sep, fixed = TRUE)[[1]],
  file.path(script_dir, "..", "Rlibs"),
  file.path(script_dir, "..", "..", "Rlibs")
)
lib_candidates <- lib_candidates[dir.exists(lib_candidates)]
if (length(lib_candidates)) .libPaths(c(lib_candidates, .libPaths()))

suppressPackageStartupMessages({
  library(scop)
  library(Seurat)
  library(ggplot2)
})

source(file.path(script_dir, "..", "scripts", "scop_helpers.R"))

cfg <- scop_skill_init(data_dir = "examples", results_dir = "results/quickstart")
srt <- readRDS(file.path(cfg$data_dir, "demo_data.rds"))
srt$condition <- factor(as.character(srt$condition), levels = c("prox", "dist"))
srt$sample_id <- paste(srt$patient, srt$condition, sep = "_")

genes <- c("PFKFB3", "SLC2A3", "VEGFA", "LDHA")
pal <- c(prox = "#5C7FA3", dist = "#C0604A")
out <- file.path(cfg$results_dir, "individual")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

for (g in genes) {
  p1 <- FeatureStatPlot(
    srt, stat.by = g, group.by = "condition", plot_type = "violin",
    comparisons = list(c("prox", "dist")),
    add_stat = "none", palcolor = pal, grid_major = FALSE,
    theme_use = "theme_classic", theme_args = list(base_size = 8),
    legend.position = "none", title = g,
    xlab = NULL, ylab = "log-normalized expression",
    sig_label = "p.format", sig_labelsize = 2.8
  ) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  save_scop_plot(p1, file.path(out, paste0(g, "_violin")), 3.4, 3.4)

  p2 <- FeatureDimPlot(
    srt, features = g, reduction = "umap",
    theme_use = "theme_classic", theme_args = list(base_size = 8),
    legend.position = "right", title = g
  ) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  save_scop_plot(p2, file.path(out, paste0(g, "_umap")), 4.0, 3.6)
}

# SCOP GroupHeatmap: genes x condition (all cells)
h <- GroupHeatmap(
  srt, features = genes, group.by = "condition",
  exp_method = "zscore", limits = c(-2, 2),
  heatmap_palcolor = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  show_row_names = TRUE, row_names_side = "left",
  show_column_names = TRUE, column_names_rot = 0,
  cluster_rows = FALSE, cluster_columns = FALSE, nlabel = 0,
  use_raster = FALSE, seed = 1, width = 4, height = 4
)
save_scop_plot(h, file.path(cfg$results_dir, "01_GroupHeatmap_genes_condition"), 7, 6)

# Patient-paired pseudobulk DESeq2 (3 patients x 2 conditions)
stats <- pseudobulk_deseq2(
  srt, genes = genes,
  group_col = "condition", group1 = "dist", group2 = "prox",
  sample_col = "sample_id", covariates = "patient",
  min_cells = 5, min_expr_cells = 4
)
write.csv(stats, file.path(cfg$results_dir, "pseudobulk_DESeq2.csv"), row.names = FALSE)

summ <- scop_gene_summary(srt, genes, group_col = "condition")
write.csv(summ, file.path(cfg$results_dir, "gene_summary.csv"), row.names = FALSE)

message("Quickstart complete: ", cfg$results_dir)
