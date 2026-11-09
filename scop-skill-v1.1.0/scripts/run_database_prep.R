# =============================================================================
# run_database_prep.R  (template)
# Prepare SCOP annotation / enrichment databases and Python envs.
# Run once; databases are cached under the SCOP data directory.
# =============================================================================
.libPaths(c(Sys.getenv("SCOP_LIBS", unset = ""), .libPaths()))
suppressPackageStartupMessages({library(scop)})

# Gene annotation / enrichment
PrepareDB(species = "Homo_sapiens", ...)
PrepareGO(species = "Homo_sapiens", ...)
PrepareKEGG(species = "Homo_sapiens", ...)
PrepareReactome(species = "Homo_sapiens", ...)
PrepareMSigDB(species = "Homo_sapiens", ...)
PrepareWikiPathway(species = "Homo_sapiens", ...)
PrepareCORUM(species = "Homo_sapiens", ...)

# TF / CCC / surfaceome
PrepareTF(species = "Homo_sapiens", ...)
PrepareJASPAR(species = "Homo_sapiens", ...)
PrepareTRRUST(species = "Homo_sapiens", ...)
PrepareTFLink(species = "Homo_sapiens", ...)
PrepareCCCDB(species = "Homo_sapiens", ...)
PrepareCellChat(species = "Homo_sapiens")
PrepareSurfaceome(species = "Homo_sapiens", ...)

# Python envs (optional; only for Python-backed methods)
ListEnv()
# PrepareEnv(modules = "scenic", ...)
# remove_python(...)
message("Database preparation finished.")
