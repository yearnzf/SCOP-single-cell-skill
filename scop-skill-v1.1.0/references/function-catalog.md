# SCOP function catalog

Generated from `scop` 0.9.2 (360 exported functions).

Categories are heuristic; see `references/coverage-matrix.md` for workflows.

## Core utilities (13)

| Function | Description |
|---|---|
| `%>%` | Objects exported from other packages |
| `%||%` | Objects exported from other packages |
| `AddFeaturesData` | Add features data |
| `ConvertHomologs` | Convert homologous gene symbols in expression objects |
| `DefaultReduction` | Find the default reduction name in a Seurat object |
| `GeneConvert` | Gene ID conversion function using biomart |
| `GetAssayData5` | Get expression data from Assay5 or Seurat object |
| `GetFeaturesData` | Get features data |
| `GetSimilarFeatures` | Find features with expression patterns similar to provided features |
| `RunGNIPLR` | Infer gene regulatory networks with GNIPLR |
| `RunSCExplorer` | Run SCExplorer |
| `srt_append` | Append a Seurat object to another |
| `srt_reorder` | Reorder idents by the gene expression |

## Data I/O & conversion (16)

| Function | Description |
|---|---|
| `adata_to_srt` | Convert an anndata object to a seurat object |
| `CheckDataList` | Check and preprocess a list of Seurat objects |
| `CheckDataMerge` | Check and preprocess a merged seurat object |
| `CheckDataType` | Check and report the type of data in Seurat object |
| `FetchDataZero` | FetchData but with zeroes for unavailable genes |
| `FetchH5` | Fetch data from the hdf5 file and returns a Seurat object |
| `giotto_to_srt` | Convert Giotto to Seurat |
| `h5ad_to_srt` | Read an .h5ad file and convert to a Seurat |
| `loom_to_adata` | Read a .loom file as an AnnData object |
| `loom_to_srt` | Read a .loom file and convert to a Seurat |
| `RecoverCounts` | Attempt to recover raw counts from the normalized matrix |
| `spe_to_srt` | Convert SpatialExperiment to Seurat |
| `srt_to_adata` | Convert a Seurat object to an AnnData object |
| `srt_to_giotto` | Convert Seurat to a native Giotto object |
| `srt_to_h5ad` | Convert a Seurat object to an .h5ad file |
| `srt_to_spe` | Convert Seurat to SpatialExperiment |

## SCENIC & transcription factors (17)

| Function | Description |
|---|---|
| `AddATACMotifs` | Add native motif matches and sequence statistics |
| `AnnotateFeatures` | Annotate Features |
| `DorotheaPlot` | Plot transcription factor activity |
| `PrepareJASPAR` | Prepare JASPAR databases |
| `PrepareTF` | Prepare TF databases |
| `PrepareTFLink` | Prepare TFLink databases |
| `PrepareTRRUST` | Prepare TRRUST databases |
| `RunCisTarget` | Run cisTarget motif enrichment on a GRN adjacency table |
| `RunDorothea` | Run transcription factor activity inference |
| `RunGENIE3` | Infer gene regulatory networks with GENIE3 |
| `RunGRN` | Infer gene regulatory networks with a selected backend method |
| `RunGRNBoost2` | Infer gene regulatory networks with GRNBoost2 |
| `RunMotifEnrichment` | Test motif enrichment within an explicit peak universe |
| `RunSCENIC` | Run SCENIC gene regulatory network analysis |
| `RunSCENICPlus` | Read SCENIC+ results or run an approximate SCENICPlus-compatible eGRN analysis |
| `SCENICPlot` | Plot SCENIC or SCENIC+ regulon activity |
| `SCENICPlusPlot` | Plot SCENIC+ eRegulon activity |

## Gene scoring & module scores (5)

| Function | Description |
|---|---|
| `AddModuleScore` | Add expression-program scores with a sparse native fast path |
| `AggregateExpression` | Aggregate expression by cell group |
| `AverageExpression` | Average expression by cell group |
| `CellScoring` | Cell scoring |
| `CellScoringPlot` | Plot CellScoring results from a Seurat object |

## ATAC & multiome (13)

| Function | Description |
|---|---|
| `AnnotateATACPeaks` | Annotate ATAC peaks with nearest genomic features |
| `CheckATACData` | Check an ATAC object without changing it |
| `CombineATACPeaks` | Combine explicit peak sets into one common interval universe |
| `CoverageTrackPlot` | Coverage track plot for ATAC data |
| `ExportATACResults` | Export an ATAC object and explicit result/resource tables |
| `FilterATACPeaks` | Filter an ATAC peak universe with explicit criteria |
| `ReadATACData` | Read raw ATAC counts into a ChromatinAssay |
| `ReadMultiomeData` | Read paired 10x RNA and ATAC counts |
| `RecountATAC` | Recount fragments over an explicit common ATAC peak universe |
| `RunATACQC` | Run scATAC quality control metrics |
| `RunChromVAR` | Compute native chromVAR deviations from an explicit sequence reference |
| `RunGeneActivity` | Quantify gene activity from indexed ATAC fragments |
| `RunLinkPeaks` | Link ATAC peaks to paired measured expression |

## Integration & batch correction (30)

| Function | Description |
|---|---|
| `BBKNN_integrate` | BBKNN integration function |
| `CCA_integrate` | Seurat v5 CCA integration |
| `ComBat_integrate` | ComBat integration function |
| `Conos_integrate` | Conos integration function |
| `Coralysis_integrate` | Coralysis integration function |
| `CSS_integrate` | CSS integration function |
| `fastMNN_integrate` | FastMNN integration function |
| `fastMNN5_integrate` | Seurat v5 fastMNN integration |
| `FindMultiModalNeighbors` | Seurat weighted-nearest-neighbor construction |
| `GLUE_integrate` | GLUE integration function |
| `Harmony_integrate` | Harmony integration function |
| `Harmony5_integrate` | Seurat v5 Harmony integration |
| `integration_scop` | Deprecated workflow entry points |
| `IntegrationBenchmarkPlot` | Plot integration benchmark results |
| `LIGER_integrate` | LIGER integration function |
| `MNN_integrate` | MNN integration function |
| `MultiMAP_integrate` | MultiMAP integration function |
| `RPCA_integrate` | Seurat v5 RPCA integration |
| `RunCCA` | Run canonical correlation analysis |
| `RunCoEmbedding` | Co-embed reference and query cells |
| `RunIntegration` | Integration workflow |
| `RunIntegrationBenchmark` | Benchmark single-cell integration methods |
| `RunLISI` | Compute LISI scores on a Seurat object |
| `Scanorama_integrate` | Scanorama integration function |
| `scVI_integrate` | ScVI integration function |
| `scVI5_integrate` | Seurat v5 scVI integration |
| `SelectIntegrationFeatures` | Select integration features with exact vectorized ranking |
| `Seurat_integrate` | Seurat integration function |
| `Uncorrected_integrate` | Uncorrected integration function |
| `WNN_integrate` | WNN integration function |

## Trajectory, pseudotime & velocity (26)

| Function | Description |
|---|---|
| `BranchStreamPlot` | Branch Stream Plot |
| `CellRankPlot` | Plot stored CellRank results |
| `CytoTRACEPlot` | Plot CytoTRACE 2 Results |
| `FitDevoPlot` | Plot FitDevo results |
| `LineagePlot` | Lineage Plot |
| `MistyRPlot` | Plot stored MISTy results |
| `PAGAPlot` | PAGA plot |
| `PalantirTrajectoryPlot` | Plot Palantir trajectories |
| `PseudotimeProjectionPlot` | Pseudotime Projection Plot |
| `RunCellRank` | Run CellRank analysis |
| `RunCellRankEnrichment` | Enrich CellRank trend modules |
| `RunCellRankTrends` | Fit and cluster CellRank lineage trends |
| `RunCytoTRACE` | Run CytoTRACE 2 |
| `RunFitDevo` | Run FitDevo developmental potential scoring |
| `RunMistyR` | Run mistyR multiview spatial modeling |
| `RunMonocle2` | Run Monocle2 analysis |
| `RunMonocle3` | Run Monocle3 analysis |
| `RunPAGA` | Run PAGA analysis |
| `RunPalantir` | Run Palantir analysis |
| `RunPHATE` | Run PHATE |
| `RunSCVELO` | Run scVelo workflow |
| `RunSlingshot` | RunSlingshot |
| `RuntAge` | Run tAge transcriptomic aging-clock prediction |
| `RunWOT` | Run WOT analysis |
| `tAgePlot` | Plot tAge transcriptomic aging-clock predictions |
| `VelocityPlot` | Velocity Plot |

## Cell-cell communication (23)

| Function | Description |
|---|---|
| `ccc_to_adata` | Convert CCC results to OmicVerse communication AnnData |
| `ccc_to_liana` | Convert CCC results to a LIANA-like table |
| `CCCHeatmap` | CCC heatmap and dot matrix plot |
| `CCCNetworkPlot` | CCC network and flow plots |
| `CCCStatPlot` | CCC statistical distribution and summary plots |
| `GetCCCObject` | Get an upstream CellChat-family object |
| `ListCCCDB` | List CCC ligand-receptor and prior-model databases |
| `PrepareCCCDB` | Prepare cell-cell communication databases |
| `RunCCC` | Run common cell-cell communication analyses |
| `RunCellChat` | Run CellChat analysis |
| `RunCellphoneDB` | Run CellphoneDB analysis |
| `RunCOMMOT` | Run official COMMOT spatial communication analysis |
| `RunLIANA` | Run LIANA cell-cell communication analysis |
| `RunMDIC3` | Run MDIC3 cell-cell communication inference |
| `RunMultiNichenetr` | Run MultiNicheNet analysis |
| `RunNichenetr` | Run NicheNet analysis |
| `RunSecAct` | Run SecAct secreted protein activity inference |
| `RunSecActCCC` | Run SecAct cell-cell communication analysis |
| `RunSecActPatternGenes` | Extract SecAct pattern-associated secreted proteins |
| `RunSecActSignalingPattern` | Run SecAct spatial signaling pattern analysis |
| `RunSecActVelocity` | Run SecAct spatial signaling velocity |
| `RunSpaTalk` | Run official SpaTalk spatial communication analysis |
| `RunSpatialCellChat` | Run Spatial CellChat analysis |

## Visualisation & themes (26)

| Function | Description |
|---|---|
| `Cell2locationPlot` | Plot cell2location spatial results |
| `CellCorHeatmap` | Cell Correlation Heatmap |
| `CellDensityPlot` | Cell density plot |
| `CellDimPlot` | Cell Dimensional Plot |
| `CellDimPlot3D` | 3D-Dimensional reduction plot for cell classification visualization. |
| `ClusterTreePlot` | Cluster tree plot |
| `COMMOTPlot` | Plot stored COMMOT results |
| `FeatureCorPlot` | Features correlation plot |
| `FeatureDimPlot` | Feature values on a 2D reduction |
| `FeatureDimPlot3D` | Feature values on a 3D reduction |
| `FeatureHeatmap` | Feature Heatmap |
| `FeatureStatPlot` | Feature statistical plots |
| `FerrisWheelPlot` | Ferris Wheel Plot |
| `GeneImmuneCorPlot` | Gene-immune correlation butterfly plot |
| `GroupHeatmap` | Group heatmap |
| `ImmuneAbundancePlot` | Immune abundance plots |
| `MetaCellPlot` | Visualize metacell partitions on a dimensionality reduction |
| `NMFHeatmap` | NMF similarity heatmap |
| `ProjectionPlot` | Projection Plot |
| `scop_logo` | scop logo |
| `SpaTalkPlot` | Plot stored SpaTalk results |
| `StatialKontextualPlot` | Plot stored Statial Kontextual results |
| `STdeconvolvePlot` | Plot STdeconvolve topic proportions |
| `TACSPlot` | Transcript-averaged cell scoring (TACS) |
| `theme_scop` | Default ggplot2 theme for scop plots |
| `theme_spatial` | ggplot2 theme for spatial plots |

## Preprocessing & QC (40)

| Function | Description |
|---|---|
| `CellCycleScoring` | Cell-cycle scoring with native sparse module scoring |
| `DimsEstimatePlot` | Dimension estimate diagnostic plot |
| `FindClusters` | Seurat clustering with transparent backend delegation |
| `FindNeighbors` | Find nearest neighbors |
| `FindVariableFeatures` | Find variable features |
| `NormalizeData` | Normalize a single-cell object |
| `RunBANKSY` | Run BANKSY spatial clustering |
| `RunCellCycle` | Run cell cycle scoring |
| `RunCellQC` | Cell-level quality control |
| `RunCHOIR` | Run CHOIR clustering |
| `RunCSSMap` | Single-cell reference mapping with CSS method |
| `RunDecontX` | Ambient RNA decontamination with decontX |
| `RunDimsEstimate` | Estimate useful dimensions from a reduction |
| `RunDimsReduction` | Run dimension reduction |
| `RunDM` | Run diffusion map (DM) |
| `RunDoubletCalling` | Doublet calling |
| `RunDoubletDetection` | Run doublet-calling with DoubletDetection |
| `RunFR` | Run Force-Directed Layout (Fruchterman-Reingold algorithm) |
| `RunFWP` | Run FWP feature-weight phenotype scoring |
| `RunGLMPCA` | Run generalized principal components analysis (GLMPCA) |
| `RunHarmony2` | Run Harmony algorithm |
| `RunLargeVis` | Run LargeVis (Dimensionality Reduction with a LargeVis-like method) |
| `RunmcRigor` | Run mcRigor metacell partition assessment |
| `RunMDS` | Run MDS (multi-dimensional scaling) |
| `RunMetaCell` | Run metacell partitioning for single-cell data |
| `RunNMF` | Run NMF |
| `RunPaCMAP` | Run PaCMAP |
| `RunPCA` | Run principal component analysis |
| `RunPCAMap` | Single-cell reference mapping with PCA method |
| `RunscDblFinder` | Run doublet-calling with scDblFinder |
| `Runscds` | Run doublet-calling with scds |
| `RunScrublet` | Run doublet-calling with Scrublet |
| `RunSmoothClust` | Run smoothclust spatial domain clustering |
| `RunSpotQC` | Run spot-level quality control |
| `RunStandardWorkflow` | Standard single-cell and spatial workflow |
| `RunTriMap` | Run TriMap |
| `RunUMAP` | Run UMAP |
| `RunUMAP2` | Run UMAP |
| `ScaleData` | Scale expression data |
| `standard_scop` | Deprecated workflow entry points |

## Cell composition (7)

| Function | Description |
|---|---|
| `CellStatPlot` | Statistical plot of cells |
| `ProportionTestPlot` | Proportion Test Plot |
| `RunMilo` | Milo differential abundance wrapper |
| `RunPropeller` | Propeller differential abundance wrapper |
| `RunProportionTest` | Proportion Test |
| `RunRareQ` | RareQ rare-cell population detection |
| `RunscCODA` | scCODA differential abundance |

## Cell annotation & reference mapping (22)

| Function | Description |
|---|---|
| `CellTypistModels` | Get available CellTypist models |
| `CNVPlot` | Plot SCOP CNV results |
| `EstimateGenePlot` | Gene and ESTIMATE score relationship plots |
| `EstimateScorePlot` | ESTIMATE score plots |
| `RenameClusters` | Rename clusters for the Seurat object |
| `RenameFeatures` | Rename features for the Seurat object |
| `RunCellTypist` | Run CellTypist cell type annotation |
| `RunCNV` | Run copy-number alteration inference |
| `RunESTIMATE` | Run ESTIMATE tumor microenvironment scoring |
| `RunKNNMap` | Single-cell reference mapping with KNN method |
| `RunKNNPredict` | Run KNN prediction |
| `RunLabelTransfer` | Transfer reference labels to query cells |
| `RunReferenceMapping` | Map query cells into a reference space |
| `RunSciBet` | Annotate single cells with a C++ SciBet-compatible classifier |
| `RunscMalignantFinder` | Run scMalignantFinder malignant cell identification |
| `RunscMalignantRegion` | Run scMalignantFinder malignant spatial region identification |
| `RunscMalignantStates` | Run scMalignantFinder cancer cell state scoring |
| `RunScmap` | Annotate single cells using scmap |
| `RunSeuratMap` | Single-cell reference mapping with Seurat method |
| `RunSingleR` | Annotate single cells using SingleR |
| `RunSymphonyMap` | Single-cell reference mapping with Symphony method |
| `TrainCellTypist` | Train a CellTypist model |

## Environment & database preparation (40)

| Function | Description |
|---|---|
| `check_python` | Check and install python packages |
| `CreateDataFile` | Create HDF5 data file from Seurat object |
| `CreateMetaFile` | Create Meta File in HDF5 format from Seurat object |
| `CycGenePrefetch` | Prefetch cell cycle genes |
| `db_DoubletDetection` | Deprecated doublet-calling entry points |
| `db_scDblFinder` | Deprecated doublet-calling entry points |
| `db_scds` | Deprecated doublet-calling entry points |
| `db_Scrublet` | Deprecated doublet-calling entry points |
| `env_requirements` | Python environment requirements |
| `ListDB` | List cached gene annotation databases |
| `ListEnv` | List conda-compatible Python environments |
| `PrepareCellChat` | Prepare CellChat databases |
| `PrepareCellTalk` | Prepare CellTalk databases |
| `PrepareChromosome` | Prepare Chromosome databases |
| `PrepareCORUM` | Prepare CORUM databases |
| `PrepareCSPA` | Prepare CSPA databases |
| `PrepareCustomDB` | Prepare custom gene annotation databases |
| `PrepareCytoTRACE2` | Prepare CytoTRACE2 model resources |
| `PrepareDB` | Prepare databases and reference resources |
| `PrepareDO` | Prepare DO databases |
| `PrepareENCODE` | Prepare ENCODE databases |
| `PrepareEnv` | Prepare the python environment |
| `PrepareEnzyme` | Prepare Enzyme databases |
| `PrepareGeneType` | Prepare GeneType databases |
| `PrepareGO` | Prepare GO databases |
| `PrepareHPO` | Prepare HPO databases |
| `PrepareHTFtarget` | Prepare hTFtarget databases |
| `PrepareIREA` | Prepare Immune Dictionary references |
| `PrepareKEGG` | Prepare KEGG databases |
| `PrepareMP` | Prepare MP databases |
| `PrepareMSigDB` | Prepare MSigDB databases |
| `PreparePFAM` | Prepare PFAM databases |
| `PrepareReactome` | Prepare Reactome databases |
| `PrepareSCExplorer` | Prepare Seurat objects for the SCExplorer |
| `PrepareSPRomeDB` | Prepare SPRomeDB databases |
| `PrepareSurfaceome` | Prepare Surfaceome databases |
| `PrepareVerSeDa` | Prepare VerSeDa databases |
| `PrepareWikiPathway` | Prepare WikiPathway databases |
| `remove_python` | Remove Python packages from a conda-compatible Python environment |
| `RemoveEnv` | Remove a conda-compatible Python environment |

## Spatial transcriptomics & deconvolution (43)

| Function | Description |
|---|---|
| `DeconvolutionPlot` | Plot deconvolution results |
| `GetSpatialDMResult` | Access stored SpatialDM results |
| `GetSpatialGraph` | Read or convert a stored spatial graph |
| `GetSpatialResult` | Access spatial analysis results |
| `ReadSpatialData` | Load and validate Visium, Visium HD, or Xenium data |
| `RunBayesSpace` | Run BayesSpace spatial clustering |
| `RunCARD` | Run CARD spatial deconvolution |
| `RunCell2fate` | Run Cell2fate RNA velocity analysis |
| `RunCell2location` | Run cell2location spatial deconvolution |
| `RunCIBERSORT` | Run CIBERSORT deconvolution |
| `RunCSIDE` | Run C-SIDE spatial differential expression |
| `RunCytoSPACE` | Run CytoSPACE spatial assignment |
| `RunDeconvolution` | Run bulk or pseudobulk deconvolution |
| `RunMERINGUE` | Run MERINGUE spatial autocorrelation analysis |
| `RunRCTD` | Run RCTD spatial deconvolution |
| `RunscOMM` | Run scOMM label prediction |
| `RunSpaNorm` | Run SpaNorm spatial normalization |
| `RunSpatialBenchmark` | Benchmark spatial domain clustering methods |
| `RunSpatialDM` | Run official SpatialDM spatial ligand-receptor association analysis |
| `RunSpatialEcoTyper` | Run SpatialEcoTyper spatial ecotype analysis |
| `RunSpatialGradientFeatures` | Run spatial gradient feature screening |
| `RunSpatialIntegration` | Run multi-sample spatial integration |
| `RunSpatialNeighborhood` | Run spatial neighborhood statistics |
| `RunSpatialNetwork` | Build a native spatial network |
| `RunSpatialVariableFeatures` | Run spatial variable feature detection |
| `RunSPOTlight` | Run SPOTlight spatial deconvolution |
| `RunSpotSweeper` | Run SpotSweeper spatial quality control |
| `RunStatialKontextual` | Run Statial Kontextual spatial relationships |
| `RunSTdeconvolve` | Run STdeconvolve reference-free spatial deconvolution |
| `SetSpatialImageAxes` | Persist the coordinate axis convention of a VisiumV2 image |
| `SpatialBenchmarkPlot` | Plot spatial domain clustering benchmarks |
| `SpatialCellChatPlot` | Plot stored SpatialCellChat results |
| `SpatialCellPlot` | Plot spatial cell boundaries |
| `SpatialCoordinates` | Read spatial coordinates with an explicit coordinate contract |
| `SpatialDeconvolutionPlot` | Plot stored spatial deconvolution proportions |
| `SpatialDMPlot` | Plot stored SpatialDM results |
| `SpatialGradientPlot` | Plot spatial gradient screening results |
| `SpatialIntegrationPlot` | Plot spatial integration results |
| `SpatialNeighborhoodPlot` | Spatial neighborhood plot |
| `SpatialNeighborhoodProfile` | Profile cell neighborhoods at several distances |
| `SpatialNetworkPlot` | Plot a native spatial network |
| `SpatialSpotPlot` | Spatial spot plot |
| `SpatialVariableFeaturePlot` | Plot spatial variable feature results |

## Differential expression & markers (10)

| Function | Description |
|---|---|
| `DEtestManhattanPlot` | DEtest Manhattan Plot |
| `DEtestPlot` | Differential Expression Test Plot |
| `DEtestRingPlot` | DEtest Ring Plot |
| `FindAllMarkers` | Find markers for all groups |
| `FindExpressedMarkers` | Find Expressed Markers |
| `FindMarkers` | Find markers between groups |
| `FoldChange` | Calculate expression fold change |
| `RunDEtest` | Differential gene test |
| `RunPermutation` | Permutation-based proportion test |
| `VolcanoPlot` | Volcano Plot |

## Enrichment, GSEA, GSVA & metabolism (18)

| Function | Description |
|---|---|
| `DynamicHeatmap` | Heatmap plot for dynamic features along lineages |
| `DynamicPlot` | Plot dynamic features across pseudotime |
| `EnrichmentPlot` | Enrichment Plot |
| `GSEAPlot` | GSEA Plot |
| `GSVAPlot` | Plots for GSVA (Gene Set Variation Analysis) |
| `IREAPlot` | Plot IREA response or polarization results |
| `MetabolismPlot` | Plots for metabolism pathway scoring |
| `RunDynamicEnrichment` | RunDynamicEnrichment |
| `RunDynamicFeatures` | Calculates dynamic features for lineages |
| `RunEnrichment` | Perform the enrichment analysis (over-representation) on the genes |
| `RunGSEA` | Perform the enrichment analysis (GSEA) on the genes |
| `RunGSVA` | Perform Gene Set Variation Analysis (GSVA) |
| `RunIREA` | Analyse cytokine responses or immune-cell polarization |
| `RunMetabolism` | Run metabolism pathway scoring |
| `RunscFEA` | Run scFEA flux estimation for a Seurat object |
| `scFEABalanceBarPlot` | Plot scFEA metabolite balance changes |
| `scFEAHeatmap` | Plot scFEA module flux heatmap |
| `scFEAVolcanoPlot` | Plot scFEA flux Cohen's d volcano plots |

## Perturbation & regulatory networks (11)

| Function | Description |
|---|---|
| `PlotscPagwas` | Plot scPagwas Scores |
| `RunAugur` | Prioritize perturbed cell types using Augur |
| `RunScissor` | Run Scissor phenotype-associated cell selection |
| `RunscPagwas` | Run scPagwas |
| `RunscTenifoldKnk` | Run scTenifoldKnk in-silico knockout analysis |
| `RunscTenifoldNet` | Run scTenifoldNet network comparison |
| `RunVECTOR` | Run VECTOR developmental direction inference |
| `ScissorPlot` | Plot Scissor results |
| `scTenifoldKnkPlot` | scTenifoldKnk Plot |
| `scTenifoldNetPlot` | scTenifoldNet Plot |
| `VECTORPlot` | Plot VECTOR results |

