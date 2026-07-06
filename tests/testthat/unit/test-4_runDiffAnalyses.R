### ============================================================================
### [04_runDiffAnalysis]
### ----------------------------------------------------------------------------
# N. Bessoltane,
library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----
# load ecoseed data
data(ecoseed.mae)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorType"   = c("batch", "Bio", "Bio")
)

keep <- 
  colData(ecoseed.mae)$imbibition  != "EI" &
  colData(ecoseed.mae)$temperature != "Medium"

ecoseed.mae.mini <- ecoseed.mae[, keep]

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = ecoseed.mae.mini,
  omicsNames  = c("RNAtest", "metatest"),
  omicsTypes  = c("RNAseq","metabolomics"),
  factorInfo  = factorInfo)

MAE <- setModelFormula(MAE, modelFormula = "~Repeat + temperature + imbibition")

selectedContrast <- "(temperatureElevated - temperatureLow) in mean"
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


rna.se  <- MAE[["RNAtest"]]
meta.se <- MAE[["metatest"]]
  
test_that("runDiffAnalyis: RNAseq, all", {
  
  MAE1 <-
    runDiffAnalysis(
      MAE, SE.name = "RNAtest"
    )
  
  rna.se1 <-
    runDiffAnalysis(
      rna.se
    )
  
  expect_identical(MAE1[["RNAtest"]], rna.se1)
  
  diffSettings <- getDiffSettings(rna.se1)
  
  expect_identical(diffSettings$Model.formula, "~Repeat + temperature + imbibition")
  expect_identical(diffSettings$Contrasts.Sel$contrastName, selectedContrast)
  expect_identical(diffSettings$method, "edgeRglmfit")
  expect_null(diffSettings$validContrasts)
  
  
  # test contrasts input
  expect_error(
    runDiffAnalysis(
      rna.se, contrastNames = "toto"
    )
  )
  
  # "(imbibitionLI - imbibitionDS) in mean" no selected
  expect_no_error(
    runDiffAnalysis(
      rna.se, contrastNames = "(imbibitionLI - imbibitionDS) in mean"
    )
  )
  
  # validate
  rna.se1 <-
    validateContrasts(
      rna.se1, 
      analysisName = "all", 
      contrastNames = selectedContrast
    )
  diffSettings <- getDiffSettings(rna.se1)
  expect_equal(diffSettings$validContrasts, selectedContrast)
  
})


test_that("runDiffAnalyis: RNAseq, split", {
  
  rna.se1 <-
    runDiffAnalysis(
      rna.se, splitBy = "imbibition"
    )
  
  diffSettings <- getDiffSettings(rna.se1, analysisName = "imbibitionDS")
  
  expect_identical(diffSettings$Model.formula, "~Repeat + temperature")
  expect_identical(diffSettings$Contrasts.Sel$contrastName, 
                   "(temperatureElevated - temperatureLow)")
  expect_identical(diffSettings$method, "edgeRglmfit")
  expect_null(diffSettings$validContrasts)
  
  rna.se1 <-
    validateContrasts(
      rna.se1, 
      analysisName = "imbibitionDS", 
      contrastNames = "(temperatureElevated - temperatureLow)"
    )
  diffSettings <- getDiffSettings(rna.se1, analysisName = "imbibitionDS")
  expect_equal(diffSettings$validContrasts, "(temperatureElevated - temperatureLow)")
  
  expect_error(
    validateContrasts(
      rna.se1, 
      analysisName = "imbibitionDS", 
      contrastNames = selectedContrast
    )
  )

})



