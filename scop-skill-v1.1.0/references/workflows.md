# Workflow recipes

Each recipe assumes a Seurat object `srt` and the helper file:

```r
source("scripts/scop_helpers.R")
suppressPackageStartupMessages({
  library(scop); library(Seurat); library(ggplot2)
})
```

## 0. Loading and checking

```r
srt <- readRDS("data/annotated.rds")
srt$condition <- factor(as.character(srt$condition), levels = c("control", "disease"))
stopifnot(all(c("celltype", "condition") %in% colnames(srt@meta.data)))
table(srt$celltype, srt$condition)
```

## 1. QC / standard workflow (optional)

SCOP provides QC and dimensional-reduction wrappers; exact arguments vary by
version, so check `args(scop::RunCellQC)` and `args(scop::RunStandardWorkflow)`
first.

```r
srt <- RunCellQC(srt, ...)
srt <- RunStandardWorkflow(srt, ...)
srt <- RunDimsEstimate(srt, ...)
srt <- RunUMAP2(srt, dims = 1:20)
CellDimPlot(srt, group.by = "celltype")
```

If your project already has UMAP/clusters, skip this step.

## 2. Marker heatmap

```r
markers <- c("PECAM1", "VWF", "KDR", "COL1A1", "DCN", "ACTA2", "MYH11")
h <- GroupHeatmap(
  srt, features = markers,
  group.by = "celltype", split.by = "condition",
  exp_method = "zscore", limits = c(-2, 2),
  show_row_names = TRUE, row_names_side = "left",
  show_column_names = TRUE, column_names_rot = 55,
  cluster_rows = FALSE, cluster_columns = FALSE, nlabel = 0,
  use_raster = FALSE, width = 5.5, height = 6
)
save_scop_plot(h, "results/marker_heatmap", 14, 11)
```

## 3. Gene-set scoring (AUCell)

```r
gene_sets <- list(
  Glycolysis = c("PFKFB3", "HK2", "PKM", "LDHA"),
  Angiogenesis = c("VEGFA", "KDR", "ANGPT2")
)
srt <- CellScoring(srt, method = "AUCell", features = gene_sets,
                   name = "scores")
GroupHeatmap(srt, features = names(gene_sets),
             group.by = "celltype", assay = "scores", layer = "data")
FeatureStatPlot(srt, stat.by = "scores", group.by = "condition",
                plot_type = "violin", comparisons = list(c("control", "disease")))
```

## 4. Cell–cell communication

```r
srt <- RunCellChat(srt, group.by = "celltype",
                   species = "Homo_sapiens",
                   min.cells = 10, do.fast = FALSE)

save_scop_baseplot(CCCNetworkPlot(srt, plot_type = "circle"),
                   "results/ccc/circle")
save_scop_baseplot(CCCNetworkPlot(srt, plot_type = "chord"),
                   "results/ccc/chord")
save_scop_plot(CCCHeatmap(srt, plot_type = "heatmap"),
               "results/ccc/heatmap")
save_scop_plot(CCCHeatmap(srt, plot_type = "dot",
                          display_by = "interaction", top_n = 25),
               "results/ccc/dot")
```

Two-condition comparison:

```r
srt <- RunCellChat(srt, group.by = "celltype",
                   group_column = "condition",
                   group_cmp = c("control", "disease"),
                   species = "Homo_sapiens")
save_scop_plot(CCCHeatmap(srt, plot_type = "diff_heatmap"),
               "results/ccc/diff_heatmap")
save_scop_baseplot(CCCNetworkPlot(srt, plot_type = "diff_network"),
                   "results/ccc/diff_network")
```

## 5. SCENIC / regulons (cpp backend)

Download cisTarget DBs first (hg38 v10 used here):

```text
hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

```r
motif_annot <- RcisTarget::importAnnotations("motifs-...tbl")
regulators <- intersect(unique(motif_annot$TF), rownames(srt))
srt <- FindVariableFeatures(srt, nfeatures = 2000)
targets <- union(head(VariableFeatures(srt), 2000), regulators)

srt <- RunSCENIC(
  srt, assay = "RNA", layer = "counts",
  regulators = regulators, targets = targets,
  ranking_dbs = "hg38_...feather",
  motif_annotations = "motifs-...tbl",
  species = "Homo_sapiens", genome = "hg38",
  backend = "cpp", cores = 4, parallel_backend = "fork",
  min_expr_cells = 3, min_regulon_size = 10,
  max_regulon_targets = 50, n_rounds = 5000
)

SCENICPlot(srt, group.by = "celltype", plot_type = "rss_rank", top_n = 12)$plot
SCENICPlot(srt, group.by = "celltype", plot_type = "rss_heatmap", top_n = 30)$plot
SCENICPlot(srt, group.by = "celltype", plot_type = "network_graph",
           features = rownames(srt@tools$SCENIC$scores)[1])$plot
GroupHeatmap(srt, features = rownames(srt@tools$SCENIC$scores)[1:25],
             group.by = "celltype", assay = "scenic", layer = "data")
```

## 6. Per-gene expression + patient-paired pseudobulk

```r
genes <- c("PFKFB3", "SLC2A3", "VEGFA", "LDHA")
srt$condition <- factor(as.character(srt$condition), levels = c("prox", "dist"))
srt$sample_id <- paste(srt$patient, srt$condition, sep = "_")

for (g in genes) {
  p1 <- FeatureStatPlot(
    srt, stat.by = g, group.by = "condition", plot_type = "violin",
    comparisons = list(c("prox", "dist")),
    add_stat = "none", sig_label = "p.format"
  )
  save_scop_plot(p1, file.path("results/individual", paste0(g, "_violin")),
                 3.4, 3.4)

  p2 <- FeatureDimPlot(srt, features = g, reduction = "umap")
  save_scop_plot(p2, file.path("results/individual", paste0(g, "_umap")),
                 4, 3.6)
}

stats <- pseudobulk_deseq2(
  srt, genes = genes,
  group_col = "condition", group1 = "dist", group2 = "prox",
  sample_col = "sample_id", covariates = "patient"
)
write.csv(stats, "results/individual/pseudobulk_DESeq2.csv", row.names = FALSE)
```

## 7. Single-cell RNA velocity / trajectory

SCOP also wraps trajectory and velocity methods (`RunSlingshot`,
`RunSCVELO`, `RunPAGA`, `RunPalantir`, `RunCellRank`, etc.). These need extra
Python/conda environments; they are outside the pure-R scope of this skill.
Follow the same pattern: run the SCOP wrapper, then save plots with
`save_scop_plot()` (ggplot) or `save_scop_baseplot()` (base graphics).
