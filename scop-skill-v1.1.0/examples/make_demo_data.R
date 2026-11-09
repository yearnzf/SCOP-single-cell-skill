# =============================================================================
# make_demo_data.R
# Creates a tiny synthetic Seurat object used by examples/quickstart.R.
# Run from the repository root:
#   Rscript examples/make_demo_data.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(1)

genes <- c(
  "PFKFB3", "SLC2A4", "SLC2A3", "HK2", "PKM", "LDHA", "VEGFA",
  "PECAM1", "VWF", "KDR", "CDH5", "COL1A1", "DCN", "MYH11", "ACTA2",
  paste0("Gene", 1:55)
)
genes <- unique(genes)

celltypes <- c("Endothelial", "Fibroblast", "Smooth muscle")
conditions <- c("prox", "dist")
patients <- c("P1", "P2", "P3")

n_per <- 40
cells <- expand.grid(
  celltype = celltypes,
  condition = conditions,
  patient = patients,
  rep = seq_len(n_per),
  stringsAsFactors = FALSE
)
n <- nrow(cells)
cell_ids <- sprintf("Cell%03d", seq_len(n))

# Baseline expression: most genes low, a few markers high in specific lineages.
base_mu <- setNames(rep(2.5, length(genes)), genes)
base_mu[c("PECAM1", "VWF", "KDR", "CDH5")] <- 25
base_mu[c("COL1A1", "DCN")] <- 22
base_mu[c("MYH11", "ACTA2")] <- 22
base_mu[c("PFKFB3", "SLC2A3", "LDHA", "PKM", "VEGFA")] <- 12
base_mu[c("SLC2A4", "HK2")] <- 0.4

mu <- matrix(rep(base_mu, times = n), nrow = length(genes),
             dimnames = list(genes, cell_ids))
ct <- cells$celltype
mu["PECAM1", ct == "Endothelial"] <- mu["PECAM1", ct == "Endothelial"] * 4
mu["VWF", ct == "Endothelial"] <- mu["VWF", ct == "Endothelial"] * 4
mu["KDR", ct == "Endothelial"] <- mu["KDR", ct == "Endothelial"] * 3
mu["CDH5", ct == "Endothelial"] <- mu["CDH5", ct == "Endothelial"] * 3
mu["COL1A1", ct == "Fibroblast"] <- mu["COL1A1", ct == "Fibroblast"] * 5
mu["DCN", ct == "Fibroblast"] <- mu["DCN", ct == "Fibroblast"] * 5
mu["MYH11", ct == "Smooth muscle"] <- mu["MYH11", ct == "Smooth muscle"] * 5
mu["ACTA2", ct == "Smooth muscle"] <- mu["ACTA2", ct == "Smooth muscle"] * 5

# Disease (dist) effects in endothelial cells: VEGFA/SLC2A3 lower, PFKFB3 neutral.
is_ec <- ct == "Endothelial"
is_dist <- cells$condition == "dist"
mu["VEGFA", is_ec & is_dist] <- mu["VEGFA", is_ec & is_dist] * 0.35
mu["SLC2A3", is_ec & is_dist] <- mu["SLC2A3", is_ec & is_dist] * 0.55
mu["LDHA", is_ec & is_dist] <- mu["LDHA", is_ec & is_dist] * 0.7

# Add gene/cell-level overdispersion so pseudobulk DESeq2 has something to fit.
mu_noise <- mu * matrix(rgamma(length(mu), shape = 1.5, rate = 1.5),
                        nrow = nrow(mu), ncol = ncol(mu))
counts <- matrix(
  rnbinom(length(mu_noise), mu = as.vector(mu_noise), size = 3),
  nrow = nrow(mu_noise), dimnames = dimnames(mu_noise)
)
counts <- Matrix::Matrix(counts, sparse = TRUE)

meta <- data.frame(
  celltype = factor(ct, levels = celltypes),
  condition = factor(cells$condition, levels = conditions),
  patient = factor(cells$patient, levels = patients),
  row.names = cell_ids
)
meta$sample_id <- paste(meta$patient, meta$condition, sep = "_")

srt <- CreateSeuratObject(counts = counts, meta.data = meta, project = "scop-skill-demo")
srt <- NormalizeData(srt, verbose = FALSE)
srt <- FindVariableFeatures(srt, nfeatures = 60, verbose = FALSE)
srt <- ScaleData(srt, verbose = FALSE)
srt <- RunPCA(srt, npcs = 10, verbose = FALSE)
srt <- RunUMAP(srt, dims = 1:10, verbose = FALSE, seed.use = 1)

out <- file.path("examples", "demo_data.rds")
saveRDS(srt, out)
message("Wrote ", out, " (", ncol(srt), " cells x ", nrow(srt), " genes)")
