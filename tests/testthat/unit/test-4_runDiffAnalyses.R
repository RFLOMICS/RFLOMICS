### ============================================================================
### [04_runDiffAnalysis]
### ----------------------------------------------------------------------------
# N. Bessoltane,
library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----
# load ecoseed data
data(ecoseed.mini.mae)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature"),
  "factorType"   = c("batch", "Bio")
)

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = ecoseed.mini.mae,
  omicsTypes  = c("RNAseq","metabolomics"),
  factorInfo  = factorInfo)

generateModelFormulae(MAE)

MAE <- setModelFormula(MAE, modelFormula = "~Repeat + temperature")
MAE <- setSelectedContrasts(MAE, contrastNames = selectedContrast)

MAE[["RNAtest"]]  <- 
  runDataProcessing(
    MAE[["RNAtest"]], 
    samples = colnames(MAE[["RNAtest"]][-1]), 
    lowCountFilter = list(method = "CPM"), 
    normalize = list(method = "TMM")
  )

MAE[["metatest"]] <- 
  runDataProcessing(
    MAE[["metatest"]], 
    transform = list(method = "log2")
  )

selectedContrast <- "(temperatureElevated - temperatureMedium)"

rna.se  <- MAE[["RNAtest"]]
meta.se <- MAE[["metatest"]]
  
test_that("runDiffAnalyis: RNAseq", {
  
  MAE1 <-
    runDiffAnalysis(
      MAE, SE.name = "RNAtest"
    )
  
  rna.se1 <-
    runDiffAnalysis(
      rna.se
    )
  
  getDiffSettings(rna.se1)
  
  meta.se1 <-
    runDiffAnalysis(
      meta.se
    )
  
  getDiffSettings(meta.se1)
  
})
