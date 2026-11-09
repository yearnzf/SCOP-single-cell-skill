# =============================================================================
# run_expression_workflow.R
# SCOP-native expression figures + optional patient-paired pseudobulk DESeq2.
#
# Edit the CONFIG block, then:
#   Rscript scripts/run_expression_workflow.R
#
# Outputs (one figure per gene, no multi-panel):
#   results/expression/individual/<GENE>_violin.pdf/png
#   results/expression/individual/<GENE>_umap.pdf/png
#   results/expression/01_GroupHeatmap_genes_group.pdf/png
#   results/expression/pseudobulk_DESeq2.csv
#   results/expression/gene_summary.csv
# =============================================================================

# ------------------------------- CONFIG ------------------------------------
INPUT_RDS   <- "data/annotated.rds"          # Seurat object with a UMAP
OUT_DIR     <- "results/expression"
GENES       <- c("PFKFB3", "SLC2A4", "SLC2A3", "HK2", "PKM", "LDHA", "VEGFA")
GROUP_COL   <- "condition"                    # metadata column to compare
GROUP_LEVELS <- c("prox", "dist")             # reference first, test second
SAMPLE_COL  <- "sample_id"                    # unique patient x condition sample ID
COVARIATES  <- c("patient")                   # NULL for unpaired pseudobulk
REDUCTION   <- "umap"
ASSAY       <- "RNA"
GROUP_PAL   <- c(prox = "#5C7FA3", dist = "#C0604A")
# ---------------------------------------------------------------------------

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
  library(ggplot2)
})

if (!file.exists(INPUT_RDS)) stop("INPUT_RDS not found: ", INPUT_RDS)
cfg <- scop_skill_init(data_dir = dirname(INPUT_RDS), results_dir = OUT_DIR)
srt <- readRDS(INPUT_RDS)
if (!all(c(GROUP_COL, GENES) %in% c(colnames(srt@meta.data), rownames(srt)))) {
  stop("Check GROUP_COL / GENES: some are not present.")
}
srt@meta.data[[GROUP_COL]] <- factor(as.character(srt@meta.data[[GROUP_COL]]),
                                     levels = GROUP_LEVELS)
missing_genes <- setdiff(GENES, rownames(srt))
if (length(missing_genes)) {
  warning("Genes missing from object: ", paste(missing_genes, collapse = ", "))
  GENES <- intersect(GENES, rownames(srt))
}

# ---- per-gene SCOP figures -------------------------------------------------
out_ind <- file.path(cfg$results_dir, "individual")
dir.create(out_ind, recursive = TRUE, showWarnings = FALSE)
for (g in GENES) {
  p1 <- FeatureStatPlot(
    srt, stat.by = g, group.by = GROUP_COL, plot_type = "violin",
    comparisons = list(GROUP_LEVELS),
    add_stat = "none", palcolor = GROUP_PAL, grid_major = FALSE,
    theme_use = "theme_classic", theme_args = list(base_size = 8),
    legend.position = "none", title = g,
    xlab = NULL, ylab = "log-normalized expression",
    sig_label = "p.format", sig_labelsize = 2.8
  ) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  save_scop_plot(p1, file.path(out_ind, paste0(g, "_violin")), 3.4, 3.4)

  p2 <- FeatureDimPlot(
    srt, features = g, reduction = REDUCTION,
    theme_use = "theme_classic", theme_args = list(base_size = 8),
    legend.position = "right", title = g
  ) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  save_scop_plot(p2, file.path(out_ind, paste0(g, "_umap")), 4.0, 3.6)
}

# ---- one SCOP heatmap across all genes ------------------------------------
h <- GroupHeatmap(
  srt, features = GENES, group.by = GROUP_COL,
  exp_method = "zscore", limits = c(-2, 2),
  heatmap_palcolor = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  show_row_names = TRUE, row_names_side = "left",
  show_column_names = TRUE, column_names_rot = 0,
  cluster_rows = FALSE, cluster_columns = FALSE, nlabel = 0,
  use_raster = FALSE, seed = cfg$seed,
  width = 4, height = max(3, 0.5 * length(GENES))
)
save_scop_plot(h, file.path(cfg$results_dir, "01_GroupHeatmap_genes_group"),
               width = 7, height = max(5, 0.6 * length(GENES) + 3))

# ---- summary + optional pseudobulk ---------------------------------------
summ <- scop_gene_summary(srt, GENES, group_col = GROUP_COL, assay = ASSAY)
write.csv(summ, file.path(cfg$results_dir, "gene_summary.csv"), row.names = FALSE)

if (SAMPLE_COL %in% colnames(srt@meta.data) && requireNamespace("DESeq2", quietly = TRUE)) {
  stats <- pseudobulk_deseq2(
    srt, genes = GENES,
    group_col = GROUP_COL,
    group1 = GROUP_LEVELS[2], group2 = GROUP_LEVELS[1],
    sample_col = SAMPLE_COL, covariates = COVARIATES,
    assay = ASSAY, min_cells = 10, min_expr_cells = 4
  )
  write.csv(stats, file.path(cfg$results_dir, "pseudobulk_DESeq2.csv"),
            row.names = FALSE)
} else {
  message("Skipping pseudobulk: SAMPLE_COL missing or DESeq2 not installed.")
}

message("Done: ", cfg$results_dir)
