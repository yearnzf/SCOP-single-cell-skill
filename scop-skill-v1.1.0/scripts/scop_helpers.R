# =============================================================================
# scop_helpers.R
# Small, dependency-light helpers shared by the SCOP workflows in this skill.
# =============================================================================

#' Initialise a SCOP analysis session
#'
#' @param data_dir Directory containing input RDS/CSV files.
#' @param results_dir Directory for figures/tables (created if missing).
#' @param libs Optional character vector of R library paths (e.g. a project Rlibs).
#' @param seed Random seed.
#' @return A list with normalised paths and the seed.
scop_skill_init <- function(data_dir = "data",
                            results_dir = "results",
                            libs = NULL,
                            seed = 1) {
  if (!is.null(libs)) .libPaths(c(libs, .libPaths()))
  options(stringsAsFactors = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  list(
    data_dir = normalizePath(data_dir, mustWork = FALSE),
    results_dir = normalizePath(results_dir, mustWork = FALSE),
    seed = seed
  )
}

#' Default two-group palette (muted blue / terracotta)
scop_two_group_pal <- function() {
  c(control = "#5C7FA3", disease = "#C0604A")
}

#' Get an assay layer across Seurat v4/v5
scop_get_layer <- function(srt, assay = "RNA", layer = "counts") {
  tryCatch(
    SeuratObject::GetAssayData(srt, assay = assay, layer = layer),
    error = function(e) SeuratObject::GetAssayData(srt, assay = assay, slot = layer)
  )
}

#' Save a SCOP ggplot / patchwork / Heatmap to PDF + PNG
#'
#' SCOP plotting functions return ggplot2/patchwork objects; some (e.g.
#' `GroupHeatmap`) return a list with `$plot`, and `SCENICPlot`/`CCCHeatmap`
#' also return lists. This helper normalises those cases.
save_scop_plot <- function(x, file, width = 6, height = 5, dpi = 300) {
  p <- if (is.list(x) && !is.null(x$plot)) x$plot else x
  if (is.null(p)) stop("No plot found in the object passed to save_scop_plot().")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(paste0(file, ".pdf"), width = width, height = height,
                 useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  grDevices::png(paste0(file, ".png"), width = width, height = height,
                 units = "in", res = dpi)
  print(p)
  grDevices::dev.off()
  invisible(NULL)
}

#' Save base/circlize graphics (e.g. SCOP `CCCNetworkPlot`)
#'
#' SCOP network plots return a `recordedplot`; `print()` gives a blank page.
#' The expression must be re-evaluated once per device, and the device must
#' already be open. Use this helper instead of `print(recordedplot)`.
#'
#' @param expr Bare expression that draws the plot, e.g.
#'   `save_scop_baseplot(CCCNetworkPlot(...), "ccc/circle")`.
save_scop_baseplot <- function(expr, file, width = 6, height = 6, dpi = 300) {
  expr <- substitute(expr)
  env <- parent.frame()
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(paste0(file, ".pdf"), width = width, height = height,
                 useDingbats = FALSE)
  eval(expr, envir = env)
  grDevices::dev.off()
  grDevices::png(paste0(file, ".png"), width = width, height = height,
                 units = "in", res = dpi)
  eval(expr, envir = env)
  grDevices::dev.off()
  invisible(NULL)
}

#' Patient/donor-paired pseudobulk DESeq2 for a small gene panel
#'
#' Aggregates raw counts by `sample_col` (one biological sample = one
#' patient × condition), fits `~ covariates + group_col`, and tests
#' `group1` vs `group2`. Intended for 3+ donors with paired conditions.
#'
#' @param srt Seurat object.
#' @param genes Character vector of genes to report.
#' @param group_col Metadata column with condition labels.
#' @param group1 Test group (e.g. "dist").
#' @param group2 Reference group (e.g. "prox").
#' @param sample_col Metadata column identifying independent samples
#'   (e.g. patient × condition sample ID).
#' @param covariates Optional metadata columns for the design, e.g. "patient".
#' @param min_cells Minimum cells per sample.
#' @param min_expr_cells Gene is testable if detected in >= this many samples.
#' @return data.frame with gene, baseMean, log2FoldChange (group1/group2),
#'   pvalue, padj, testable, note.
pseudobulk_deseq2 <- function(srt, genes, group_col, group1, group2,
                              sample_col, covariates = NULL,
                              assay = "RNA", layer = "counts",
                              min_cells = 10, min_expr_cells = 4) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("pseudobulk_deseq2() needs DESeq2: BiocManager::install('DESeq2')")
  }
  md <- srt@meta.data
  if (!all(c(group_col, sample_col, covariates) %in% colnames(md))) {
    stop("Missing metadata columns: ",
         paste(setdiff(c(group_col, sample_col, covariates), colnames(md)),
               collapse = ", "))
  }
  samples <- unique(as.character(md[[sample_col]]))
  ncell <- table(factor(md[[sample_col]], levels = samples))
  samples <- samples[ncell >= min_cells]
  if (length(samples) < 4) {
    stop("Need at least 4 samples after min_cells filtering.")
  }
  sm <- unique(md[, c(sample_col, group_col, covariates), drop = FALSE])
  sm <- sm[!duplicated(sm[[sample_col]]), , drop = FALSE]
  rownames(sm) <- sm[[sample_col]]
  sm <- sm[samples, , drop = FALSE]
  sm[[group_col]] <- factor(as.character(sm[[group_col]]),
                            levels = c(group2, group1))
  for (cv in covariates) sm[[cv]] <- factor(as.character(sm[[cv]]))
  if (any(table(sm[[group_col]]) < 2)) {
    stop("Need >=2 samples per group after filtering.")
  }

  counts <- scop_get_layer(srt, assay = assay, layer = layer)
  pb <- sapply(samples, function(s) {
    Matrix::rowSums(counts[, md[[sample_col]] == s, drop = FALSE])
  })
  pb <- round(as.matrix(pb))
  keep <- rowSums(pb > 0) >= min_expr_cells
  genes_use <- intersect(genes, rownames(pb)[keep])
  if (length(genes_use) == 0) {
    return(data.frame(gene = genes, baseMean = NA_real_,
                      log2FoldChange = NA_real_, pvalue = NA_real_,
                      padj = NA_real_, testable = FALSE,
                      note = "not detected in enough samples",
                      stringsAsFactors = FALSE))
  }
  rhs <- c(covariates, group_col)
  design <- stats::as.formula(paste("~", paste(rhs, collapse = " + ")))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pb[genes_use, , drop = FALSE],
    colData = sm,
    design = design
  )
  dds <- tryCatch(
    DESeq2::DESeq(dds, quiet = TRUE),
    error = function(e) {
      # Very small / low-dispersion datasets can fail parametric dispersion fit.
      tryCatch(DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE),
               error = function(e2) NULL)
    }
  )
  if (is.null(dds)) {
    # Last resort: use gene-wise dispersion estimates directly.
    dds <- tryCatch({
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
      SummarizedExperiment::mcols(dds)$dispersion <-
        SummarizedExperiment::mcols(dds)$dispGeneEst
      DESeq2::nbinomWaldTest(dds)
    }, error = function(e) NULL)
  }
  if (is.null(dds)) {
    return(data.frame(
      gene = genes, baseMean = NA_real_, log2FoldChange = NA_real_,
      pvalue = NA_real_, padj = NA_real_, testable = FALSE,
      note = "DESeq2 dispersion fit failed; use more samples/genes",
      stringsAsFactors = FALSE
    ))
  }
  res <- DESeq2::results(dds, contrast = c(group_col, group1, group2))
  out <- data.frame(
    gene = genes_use,
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    pvalue = res$pvalue,
    padj = res$padj,
    testable = TRUE,
    note = "",
    stringsAsFactors = FALSE
  )
  missing <- setdiff(genes, out$gene)
  if (length(missing)) {
    out <- rbind(out, data.frame(
      gene = missing, baseMean = NA_real_, log2FoldChange = NA_real_,
      pvalue = NA_real_, padj = NA_real_, testable = FALSE,
      note = "filtered by low detection", stringsAsFactors = FALSE
    ))
  }
  out[match(genes, out$gene), , drop = FALSE]
}

#' Per-gene summary table (mean expression and % positive by group)
scop_gene_summary <- function(srt, genes, group_col, assay = "RNA",
                              layer = "data") {
  expr <- scop_get_layer(srt, assay = assay, layer = layer)
  md <- srt@meta.data
  groups <- unique(as.character(md[[group_col]]))
  out <- do.call(rbind, lapply(genes, function(g) {
    do.call(rbind, lapply(groups, function(gr) {
      ii <- as.character(md[[group_col]]) == gr
      v <- as.numeric(expr[g, ii])
      data.frame(gene = g, group = gr, n_cells = sum(ii),
                 mean_expression = mean(v),
                 pct_positive = 100 * mean(v > 0),
                 stringsAsFactors = FALSE)
    }))
  }))
  out$group <- factor(out$group, levels = groups)
  out
}
