### ============================================================================
### [03_runDataProcessing] unit tests
### ----------------------------------------------------------------------------
# N. Bessoltane

library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ---- 
# load ecoseed data
data(ecoseed.mae)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorType"   = c("batch", "Bio", "Bio")
)

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = ecoseed.mae,
  omicsNames  = c("RNAtest", "protetest", "metatest"),
  omicsTypes  = c("RNAseq","proteomics","metabolomics"),
  factorInfo  = factorInfo)

# set stat setting
## choice of model formulae
formulae <- generateModelFormulae(MAE)
MAE <- setModelFormula(MAE, modelFormula = formulae[[2]])

## list of all hypothesis grouped by contarst type (run on MAE or SE)
selectedContrasts <- generateExpressionContrast(MAE)$averaged
MAE <- setSelectedContrasts(MAE, contrastList = selectedContrasts)

## processed data
sampleToKeep <- colnames(MAE[["RNAtest"]])[-1]
MAE1 <- 
  runDataProcessing(MAE, SE.name = "RNAtest",
                    samples = sampleToKeep,
                    filterStrategy = "NbReplicates", 
                    cpmCutoff = 1,
                    normMethod = "TMM") |>
  runDataProcessing(SE.name = "protetest",
                    transformMethod = "none",
                    normMethod = "none") |> 
  runDataProcessing(SE.name = "metatest",
                    transformMethod = "log2",
                    normMethod = "median")

########### FUNCTIONS TESTS ###########

# ---- getProcessedData ----
test_that("getProcessedData returned value", {
  
  rna.S1  <- MAE1[["RNAtest"]]
  prot.S1 <- MAE1[["protetest"]]
  meta.S1 <- MAE1[["metatest"]]
  
  # getProcessedData-RflomicsMAE
  MAE2   <- MAE1 |>
    getProcessedData(SE.name = "RNAtest", norm = TRUE, log = TRUE) |>
    getProcessedData(SE.name = "protetest", norm = TRUE) |>
    getProcessedData(SE.name = "metatest", norm = TRUE)
  
  # getProcessedData-RflomicsSE
  rna.S2  <- getProcessedData(rna.S1, norm = TRUE, log = TRUE)
  prot.S2 <- getProcessedData(prot.S1, norm = TRUE)
  meta.S2 <- getProcessedData(meta.S1, norm = TRUE)
  
  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE2))
  expect_true("RflomicsSE"  %in% is(rna.S2))
  expect_true("RflomicsSE"  %in% is(prot.S2))
  expect_true("RflomicsSE"  %in% is(meta.S2))
  
  expect_identical(rna.S2,  MAE2[["RNAtest"]])
  expect_identical(prot.S2, MAE2[["protetest"]])
  expect_identical(meta.S2, MAE2[["metatest"]])
  
  # expected results
  expect_true( RFLOMICS:::.isFiltered(rna.S2))
  expect_false(RFLOMICS:::.isTransformed(rna.S2))
  expect_true( RFLOMICS:::.isNormalized(rna.S2))
  expect_identical(metadata(rna.S2)$DataProcessing$log, "log2")
  
  expect_false(RFLOMICS:::.isFiltered(prot.S2))
  expect_true(RFLOMICS:::.isTransformed(prot.S2))
  expect_true(RFLOMICS:::.isNormalized(prot.S2))  
  
  expect_false(RFLOMICS:::.isFiltered(meta.S2))
  expect_true(RFLOMICS:::.isTransformed(meta.S2))
  expect_true(RFLOMICS:::.isNormalized(meta.S2))
  
  # get Processing data from processed data
  expect_warning(getProcessedData(rna.S2, norm = TRUE, log = TRUE))
  expect_warning(getProcessedData(prot.S2, norm = TRUE))
  expect_warning(getProcessedData(meta.S2, norm = TRUE))
})

test_that("getProcessedData RNAseq case", {
  
  rna.S1  <- MAE1[["RNAtest"]]
  
  # get processed data
  rna.S2  <- getProcessedData(rna.S1, norm = TRUE, log = TRUE)
  
  expect(!identical(rna.S2, rna.S1), failure_message = "")
  expect(!identical(assay(rna.S2), assay(rna.S1)), failure_message = "")
  expect(!identical(metadata(rna.S2), metadata(rna.S1)), failure_message = "")
  expect(!identical(colnames(rna.S2), colnames(rna.S1)), failure_message = "")
  expect(!identical(names(rna.S2), names(rna.S1)), failure_message = "")
  
  ## filtering results
  expect_identical(colnames(rna.S2), sampleToKeep)
  expect_identical(sort(c(names(rna.S2), getFilteredFeatures(rna.S2))), names(rna.S1))
  expect_identical(getFilterSettings(rna.S2), getFilterSettings(rna.S1))
  expect_identical(getAnalysis(rna.S2, name = "PCAlist")$norm,
                   getAnalysis(rna.S1, name = "PCAlist")$norm)
  
  ## transformation results
  expect_identical(getAnalysis(rna.S2, "DataProcessing", "Transformation"), list())
  
  ## normalization results
  expect(!identical(getAnalysis(rna.S2, "DataProcessing", "Normalization"), 
                    getAnalysis(rna.S1, "DataProcessing", "Normalization")),
         failure_message = "")
  
  ## contrast 
  expect_equal(getSelectedContrasts(rna.S2), getSelectedContrasts(rna.S1))
})

test_that("getProcessedData arg", {
  
  rna.S1  <- MAE1[["RNAtest"]]
  prot.S1 <- MAE1[["protetest"]]
  meta.S1 <- MAE1[["metatest"]]
  
  # getProcessedData with default arg values
  rna.S2  <- getProcessedData(rna.S1)
  prot.S2 <- getProcessedData(prot.S1)
  meta.S2 <- getProcessedData(meta.S1)
  
  ## no change applied
  expect_false(identical(rna.S2, rna.S1))
  expect_identical(prot.S2, prot.S1)
  expect_identical(meta.S2, meta.S1)
  
  # get processed data with filter = TRUE
  rna.S2  <- getProcessedData(rna.S1,  filter = TRUE)
  prot.S2 <- getProcessedData(prot.S1, filter = TRUE)
  meta.S2 <- getProcessedData(meta.S1, filter = TRUE)
  # expected results
  expect_true( RFLOMICS:::.isFiltered(rna.S2))
  expect_false(RFLOMICS:::.isTransformed(rna.S2))
  expect_false( RFLOMICS:::.isNormalized(rna.S2))
  
  expect_false(RFLOMICS:::.isFiltered(prot.S2))
  expect_false(RFLOMICS:::.isTransformed(prot.S2))
  expect_false(RFLOMICS:::.isNormalized(prot.S2))  
  
  expect_false(RFLOMICS:::.isFiltered(meta.S2))
  expect_false(RFLOMICS:::.isTransformed(meta.S2))
  expect_false(RFLOMICS:::.isNormalized(meta.S2))
  
  # get processed data with trans = TRUE
  rna.S2  <- getProcessedData(rna.S1,  trans = TRUE)
  prot.S2 <- getProcessedData(prot.S1, trans = TRUE)
  meta.S2 <- getProcessedData(meta.S1, trans = TRUE)
  # expected results
  expect_true( RFLOMICS:::.isFiltered(rna.S2))
  expect_false(RFLOMICS:::.isTransformed(rna.S2))
  expect_false( RFLOMICS:::.isNormalized(rna.S2))
  
  expect_false(RFLOMICS:::.isFiltered(prot.S2))
  expect_true(RFLOMICS:::.isTransformed(prot.S2))
  expect_false(RFLOMICS:::.isNormalized(prot.S2))  
  
  expect_false(RFLOMICS:::.isFiltered(meta.S2))
  expect_true(RFLOMICS:::.isTransformed(meta.S2))
  expect_false(RFLOMICS:::.isNormalized(meta.S2))
  
  # get processed data with norm = TRUE
  rna.S2  <- getProcessedData(rna.S1,  norm = TRUE)
  prot.S2 <- getProcessedData(prot.S1, norm = TRUE)
  meta.S2 <- getProcessedData(meta.S1, norm = TRUE)
  # expected results
  expect_true( RFLOMICS:::.isFiltered(rna.S2))
  expect_false(RFLOMICS:::.isTransformed(rna.S2))
  expect_true( RFLOMICS:::.isNormalized(rna.S2))
  
  expect_false(RFLOMICS:::.isFiltered(prot.S2))
  expect_true(RFLOMICS:::.isTransformed(prot.S2))
  expect_true(RFLOMICS:::.isNormalized(prot.S2))  
  
  expect_false(RFLOMICS:::.isFiltered(meta.S2))
  expect_true(RFLOMICS:::.isTransformed(meta.S2))
  expect_true(RFLOMICS:::.isNormalized(meta.S2))
  
  # get processed data with log = TRUE
  rna.S2  <- getProcessedData(rna.S1,  log = TRUE)
  prot.S2 <- expect_warning(getProcessedData(prot.S1, log = TRUE))
  meta.S2 <- expect_warning(getProcessedData(meta.S1, log = TRUE))
  # expected results
  expect_false(RFLOMICS:::.isFiltered(rna.S2))
  expect_false(RFLOMICS:::.isTransformed(rna.S2))
  expect_false(RFLOMICS:::.isNormalized(rna.S2))
  
  expect_false(RFLOMICS:::.isFiltered(prot.S2))
  expect_false(RFLOMICS:::.isTransformed(prot.S2))
  expect_false(RFLOMICS:::.isNormalized(prot.S2))  
  
  expect_false(RFLOMICS:::.isFiltered(meta.S2))
  expect_false(RFLOMICS:::.isTransformed(meta.S2))
  expect_false(RFLOMICS:::.isNormalized(meta.S2))
})

test_that("getProcessedData / returned contrast", {
  
  # remove 1 samples
  rna.S1 <- runDataProcessing(MAE[["RNAtest"]],
                              samples = colnames(MAE[["RNAtest"]])[-1],
                              filterStrategy = "NbReplicates", 
                              cpmCutoff = 1,
                              normMethod = "TMM")
  
  # getProcessedData
  rna.S2 <- getProcessedData(rna.S1, filter = TRUE)
  
  contrast1 <- getSelectedContrasts(rna.S1)
  contrast2 <- getSelectedContrasts(rna.S2)
  
  expect_identical(metadata(rna.S1)$design$factorType,
                   metadata(rna.S2)$design$factorType)
  expect_identical(metadata(rna.S1)$design$Model.formula,
                   metadata(rna.S2)$design$Model.formula)
  
  expect_identical(contrast1, contrast2)
  
  # remove DS modality
  rna.S1 <- runDataProcessing(MAE[["RNAtest"]],
                              samples = colnames(MAE[["RNAtest"]])[-1:-9],
                              filterStrategy = "NbReplicates", 
                              cpmCutoff = 1,
                              normMethod = "TMM")
  
  # getProcessedData
  rna.S2 <- getProcessedData(rna.S1, filter = TRUE)
  
  contrast1 <- getSelectedContrasts(rna.S1)
  contrast2 <- getSelectedContrasts(rna.S2)
  
  expect_identical(metadata(rna.S1)$design$factorType,
                   metadata(rna.S2)$design$factorType)
  expect_identical(metadata(rna.S1)$design$Model.formula,
                   metadata(rna.S2)$design$Model.formula)
  
  expect_false(identical(contrast1, contrast2))
  expect_true(unique(contrast2$contrastName %in% contrast1$contrastName))
  
  # remove imbibition factor
  expect_error(runDataProcessing(MAE[["RNAtest"]],
                                 samples = colnames(MAE[["RNAtest"]])[-1:-18],
                                 filterStrategy = "NbReplicates", 
                                 cpmCutoff = 1,
                                 normMethod = "TMM"))
  
  # less than 2 rep
  expect_error(runDataProcessing(MAE[["RNAtest"]],
                                 samples = colnames(MAE[["RNAtest"]])[-1:-2],
                                 filterStrategy = "NbReplicates", 
                                 cpmCutoff = 1,
                                 normMethod = "TMM"))
  
  # no complete
  expect_error(runDataProcessing(MAE[["RNAtest"]],
                                 samples = colnames(MAE[["RNAtest"]])[-1:-3],
                                 filterStrategy = "NbReplicates", 
                                 cpmCutoff = 1,
                                 normMethod = "TMM"))
  
})
