# =============================================================================
# build_function_catalog.R
# Regenerates references/function-catalog.md and references/scop_functions.csv
# from the installed SCOP package. Maintainer tool; not needed by end users.
#   SCOP_LIBS=/path/to/Rlibs Rscript scripts/build_function_catalog.R
# =============================================================================
lib <- Sys.getenv("SCOP_LIBS", unset = "")
if (nzchar(lib)) .libPaths(c(strsplit(lib, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
suppressPackageStartupMessages(library(scop))

ex <- sort(getNamespaceExports("scop"))
db <- tools::Rd_db("scop")

rd_tag <- function(rd, tag) {
  for (i in seq_along(rd)) {
    if (identical(attr(rd[[i]], "Rd_tag"), tag)) return(paste(unlist(rd[[i]]), collapse = ""))
  }
  ""
}
alias_map <- list()
for (nm in names(db)) {
  rd <- db[[nm]]
  ttl <- gsub("\\s+", " ", trimws(rd_tag(rd, "\\title")))
  alias_map[[sub("\\.Rd$", "", nm)]] <- ttl
  for (i in seq_along(rd)) {
    if (identical(attr(rd[[i]], "Rd_tag"), "\\alias")) {
      for (a in unlist(rd[[i]])) alias_map[[trimws(a)]] <- ttl
    }
  }
}
titles <- vapply(ex, function(f) if (!is.null(alias_map[[f]])) alias_map[[f]] else "", character(1))

patterns <- list(
  "SCENIC & transcription factors" = c(
    "^RunSCENIC", "^RunCisTarget", "^RunGENIE3", "^RunGRNBoost2$", "^RunGRN$",
    "^RunDorothea", "^DorotheaPlot", "^SCENICPlot", "^SCENICPlusPlot",
    "^PrepareTF$", "^PrepareJASPAR$", "^PrepareTRRUST$", "^PrepareTFLink$",
    "^AnnotateFeatures$", "^RunMotifEnrichment$", "^AddATACMotifs$"
  ),
  "Cell-cell communication" = c(
    "^RunCellChat", "^RunCCC", "^RunCellphoneDB", "^RunLIANA", "^RunNichenetr",
    "^RunMultiNichenetr", "^RunMDIC3", "^RunCOMMOT", "^RunSpatialCellChat",
    "^RunSpaTalk", "^RunSecAct", "^GetCCCObject$", "^ListCCCDB$",
    "^PrepareCCC", "^ccc_to_", "^CCC"
  ),
  "Spatial transcriptomics & deconvolution" = c(
    "^RunSpatial", "^Spatial", "^RunCARD", "^RunRCTD", "^RunCell2location",
    "^RunBayesSpace", "^RunSpaNorm", "^RunSTdeconvolve", "^RunSPOTlight",
    "^RunCytoSPACE", "^RunCSIDE", "^RunMERINGUE", "^RunStatial",
    "^ReadSpatialData", "^GetSpatial", "^SetSpatialImageAxes", "^SpatialCoordinates",
    "^RunDeconvolution", "^DeconvolutionPlot", "^RunCIBERSORT", "^RunBisqueRNA",
    "^RunBayesPrism", "^RunscOMM", "^RunCell2fate", "^RunSpotSweeper$"
  ),
  "ATAC & multiome" = c(
    "ATAC", "^RunChromVAR", "^RunGeneActivity", "^RunLinkPeaks",
    "^CoverageTrackPlot", "^ReadMultiomeData", "^RunSCENICPlus", "^SCENICPlusPlot"
  ),
  "Perturbation & regulatory networks" = c(
    "^RunscTenifold", "^scTenifold", "^RunScissor", "^ScissorPlot",
    "^RunAugur", "^RunscPagwas", "^PlotscPagwas", "^RunVECTOR$", "^VECTORPlot"
  ),
  "Trajectory, pseudotime & velocity" = c(
    "^RunSlingshot", "^LineagePlot", "^BranchStreamPlot", "^RunMonocle",
    "^RunPAGA", "^PAGAPlot", "^RunPalantir", "^PalantirTrajectoryPlot",
    "^RunCellRank", "^CellRankPlot", "^RunSCVELO", "^VelocityPlot",
    "^RunWOT", "^RunPHATE", "^PseudotimeProjectionPlot",
    "^RunCytoTRACE", "^CytoTRACEPlot", "^RuntAge", "^tAgePlot",
    "^RunFitDevo", "^FitDevoPlot", "^RunMistyR", "^MistyRPlot"
  ),
  "Enrichment, GSEA, GSVA & metabolism" = c(
    "^RunEnrichment", "^EnrichmentPlot", "^RunGSEA", "^GSEAPlot",
    "^RunGSVA", "^GSVAPlot", "^RunIREA", "^IREAPlot",
    "^RunDynamicEnrichment", "^RunDynamicFeatures", "^DynamicPlot", "^DynamicHeatmap",
    "^RunMetabolism", "^MetabolismPlot", "^RunscFEA", "^scFEA"
  ),
  "Differential expression & markers" = c(
    "^RunDEtest", "^DEtest", "^FindMarkers", "^FindAllMarkers",
    "^FindExpressedMarkers", "^VolcanoPlot", "^FoldChange", "^RunPermutation"
  ),
  "Gene scoring & module scores" = c(
    "^CellScoring", "^AddModuleScore", "^AverageExpression", "^AggregateExpression"
  ),
  "Cell composition" = c(
    "^CellStatPlot", "^RunProportionTest", "^ProportionTestPlot",
    "^RunscCODA", "^RunPropeller", "^RunMilo", "^RunRareQ"
  ),
  "Cell annotation & reference mapping" = c(
    "^RunSingleR", "^RunScmap", "^RunSciBet", "^RunCellTypist", "^TrainCellTypist",
    "^CellTypistModels", "^RunKNNMap", "^RunKNNPredict", "^RunLabelTransfer",
    "^RunSeuratMap", "^RunSymphonyMap", "^RunReferenceMapping",
    "^RunscMalignant", "^RunCNV", "^CNVPlot", "^RunESTIMATE",
    "^EstimateGenePlot", "^EstimateScorePlot", "^RenameClusters", "^RenameFeatures"
  ),
  "Integration & batch correction" = c(
    "_integrate$", "^integration_scop$", "^RunIntegration", "^RunLISI",
    "^SelectIntegrationFeatures", "^FindMultiModalNeighbors$", "^RunCoEmbedding$",
    "^RunCCA$", "^IntegrationBenchmarkPlot$"
  ),
  "Preprocessing & QC" = c(
    "^RunCellQC$", "^RunATACQC$", "^RunSpotQC$", "^RunDoublet", "^RunDecontX$",
    "^RunCellCycle", "^CellCycleScoring$", "^NormalizeData$",
    "^FindVariableFeatures$", "^ScaleData$", "^RunPCA$", "^RunDimsEstimate$",
    "^RunDimsReduction$", "^RunGLMPCA$", "^FindNeighbors$", "^FindClusters$",
    "^RunStandardWorkflow$", "^standard_scop$", "^RunSmoothClust$", "^RunBANKSY$",
    "^RunNMF$", "^RunHarmony2$", "^RunUMAP2?$", "^RunTriMap$", "^RunPaCMAP$",
    "^RunLargeVis$", "^RunMDS$", "^RunFR$", "^RunFWP$", "^RunPCAMap$",
    "^RunscDblFinder$", "^Runscds$", "^RunScrublet$", "^DimsEstimatePlot$",
    "^RunCSSMap$", "^RunDM$", "^RunmcRigor$", "^RunCHOIR$", "^RunMetaCell$"
  ),
  "Data I/O & conversion" = c(
    "^adata_to_srt$", "^h5ad_to_srt$", "^loom_to_", "^giotto_to_srt$",
    "^spe_to_srt$", "^srt_to_", "^ccc_to_", "^FetchH5$", "^FetchDataZero$",
    "^RecoverCounts$", "^ReadATACData$", "^ReadMultiomeData$", "^ReadSpatialData$",
    "^ExportATACResults$", "^CheckData", "^CheckATACData$", "^CheckDataList$",
    "^CheckDataMerge$", "^CheckDataType$", "^CombineATACPeaks$", "^FilterATACPeaks$"
  ),
  "Environment & database preparation" = c(
    "^Prepare", "^ListDB$", "^ListEnv$", "^RemoveEnv$", "^check_python$",
    "^env_requirements$", "^remove_python$", "^CycGenePrefetch$",
    "^CreateDataFile$", "^CreateMetaFile$", "^db_"
  ),
  "Visualisation & themes" = c(
    "^theme_", "^scop_logo$", "^CellDimPlot", "^FeatureDimPlot", "^CellCorHeatmap",
    "^CellDensityPlot", "^CellRankPlot", "^ClusterTreePlot", "^FerrisWheelPlot",
    "^GeneImmuneCorPlot", "^FeatureCorPlot", "^ProjectionPlot", "^GroupHeatmap$",
    "^FeatureHeatmap$", "^FeatureStatPlot$", "^DEtest.*Plot", "^EnrichmentPlot",
    "^GSEAPlot", "^GSVAPlot", "^LineagePlot", "^BranchStreamPlot", "^PAGAPlot",
    "^PalantirTrajectoryPlot", "^VelocityPlot", "^VECTORPlot", "^CytoTRACEPlot",
    "^RuntAge", "^tAgePlot", "^FitDevoPlot", "^MistyRPlot", "^PlotscPagwas",
    "^ScissorPlot", "^CNVPlot", "^Estimate.*Plot", "^CellScoringPlot",
    "^ProportionTestPlot", "^DeconvolutionPlot", "^Spatial.*Plot",
    "^Cell2locationPlot", "^COMMOTPlot", "^SpaTalkPlot", "^STdeconvolvePlot",
    "^StatialKontextualPlot", "^MetaCellPlot", "^NMFHeatmap", "^TACSPlot",
    "^ImmuneAbundancePlot", "^IntegrationBenchmarkPlot"
  ),
  "Core utilities" = c(
    "^srt_", "^DefaultReduction$", "^GetFeaturesData$", "^AddFeaturesData$",
    "^GetAssayData5$", "^ConvertHomologs$", "^GeneConvert$", "^GetSimilarFeatures$",
    "^RunSCExplorer$", "^RunGNIPLR$", "^%>%$", "^%\\|\\|%$"
  )
)

categorize <- function(f) {
  for (cat in names(patterns)) {
    if (any(vapply(patterns[[cat]], function(p) grepl(p, f, perl = TRUE), logical(1)))) return(cat)
  }
  "Other / utilities"
}
cat_vec <- vapply(ex, categorize, character(1))
out <- data.frame(function_name = ex, category = cat_vec, title = titles,
                  stringsAsFactors = FALSE)
dir.create("references", showWarnings = FALSE)
write.csv(out, "references/scop_functions.csv", row.names = FALSE)

lines <- c(
  "# SCOP function catalog",
  "",
  paste0("Generated from `scop` ", as.character(packageVersion("scop")),
         " (", length(ex), " exported functions)."),
  "",
  "Categories are heuristic; see `references/coverage-matrix.md` for workflows.",
  ""
)
for (cat in unique(out$category)) {
  d <- out[out$category == cat, , drop = FALSE]
  lines <- c(lines, paste0("## ", cat, " (", nrow(d), ")"), "",
             "| Function | Description |", "|---|---|")
  for (i in seq_len(nrow(d))) {
    desc <- ifelse(nzchar(d$title[i]), d$title[i], "—")
    lines <- c(lines, sprintf("| `%s` | %s |", d$function_name[i], desc))
  }
  lines <- c(lines, "")
}
writeLines(lines, "references/function-catalog.md")
message("Wrote references/function-catalog.md and references/scop_functions.csv")
print(table(out$category))
