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
MAE <- createRflomicsMAE(
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

########### FUNCTIONS TESTS ###########

# ---- runDataProcessing ----
test_that("runDataProcessing returned value", {

  MAE1    <- MAE
  sampleToKeep <- colnames(MAE[["RNAtest"]])[-1]

  ## method of RflomicsMAE class
  MAE1 <- runDataProcessing(MAE, SE.name = "RNAtest",
                            samples = sampleToKeep,
                            filterStrategy = "NbReplicates",
                            cpmCutoff = 1,
                            normMethod = "TMM")
  MAE1 <- runDataProcessing(MAE1, SE.name = "protetest",
                            transformMethod = "none",
                            normMethod = "none",
                            imputMethod = "MVI")
  MAE1 <- runDataProcessing(MAE1, SE.name = "metatest",
                            transformMethod = "log2",
                            normMethod = "median",
                            imputMethod = "MVI")

  ## returns a value of class RflomicsMAE/MAE
  expect_true("RflomicsMAE" %in% is(MAE1))

  ## Returns an object with the same length as the input
  expect_equal(length(MAE), length(MAE1))

  ## Doesn't modify the MAE@metadata
  expect_identical(metadata(MAE), metadata(MAE1))

  ## method of RflomicsSE class
  rna.S1  <- MAE[["RNAtest"]]
  rna.S1  <- runDataProcessing(rna.S1,
                               samples = sampleToKeep,
                               filterStrategy = "NbReplicates",
                               cpmCutoff = 1,
                               normMethod = "TMM")

  prot.S1 <- MAE[["protetest"]]
  prot.S1 <- runDataProcessing(prot.S1,
                               transformMethod = "none",
                               normMethod = "none",
                               imputMethod = "MVI")

  meta.S1 <- MAE[["metatest"]]
  meta.S1 <- runDataProcessing(meta.S1,
                               transformMethod = "log2",
                               normMethod = "median",
                               imputMethod = "MVI")

  ## returns a value of class RflomicsS
  expect_true("RflomicsSE" %in% is(rna.S1))
  expect_true("RflomicsSE" %in% is(prot.S1))
  expect_true("RflomicsSE" %in% is(meta.S1))

  ## we get the same results
  expect_identical(rna.S1,  MAE1[["RNAtest"]])
  expect_identical(prot.S1, MAE1[["protetest"]])
  expect_identical(meta.S1, MAE1[["metatest"]])

  ## Doesn't modify the colData
  expect_identical(colData(rna.S1),  colData(MAE[["RNAtest"]]))
  expect_identical(colData(prot.S1), colData(MAE[["protetest"]]))
  expect_identical(colData(meta.S1), colData(MAE[["metatest"]]))

  ## Doesn't modify the design
  expect_identical(metadata(rna.S1)$design$factorType,
                   metadata(MAE[["RNAtest"]])$design$factorType)
  expect_identical(metadata(prot.S1)$design$Model.formula,
                   metadata(MAE[["protetest"]])$design$Model.formula)
  expect_identical(metadata(meta.S1)$design$Contrasts.Sel,
                   metadata(MAE[["metatest"]])$design$Contrasts.Sel)

  ## Doesn't modify the matrix
  expect_identical(assay(rna.S1),  assay(MAE[["RNAtest"]]))
  expect_identical(assay(prot.S1), assay(MAE[["protetest"]]))
  expect_identical(assay(meta.S1), assay(MAE[["metatest"]]))

  ## modify the SE@metadata
  expect(!identical(metadata(MAE[["RNAtest"]]), metadata(rna.S1)),
         failure_message = "The two objects are not expected to be identical.")

  ### filtering settings: samples
  expect_identical(getSelectedSamples(rna.S1), sampleToKeep)
  expect_identical(getSelectedSamples(prot.S1), colnames(MAE[["protetest"]]))
  expect_identical(getSelectedSamples(meta.S1), colnames(MAE[["metatest"]]))
  ### filtering settings: features
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbReplicates", cpmCutoff = 1))
  expect_equal(
    getFilterSettings(prot.S1),
    list(method = "MVI", minValue = 6.226499, suppInfo = "missing value imputation"))
  expect_equal(
    getFilterSettings(meta.S1),
    list(method = "MVI", minValue = 1.98e-05, suppInfo = "missing value imputation"))
  ### Is it filtered ?
  expect_identical(RFLOMICS:::.isFiltered(rna.S1), FALSE)
  expect_identical(RFLOMICS:::.isFiltered(prot.S1), FALSE)
  expect_identical(RFLOMICS:::.isFiltered(meta.S1), FALSE)

  ### transformation settings
  expect_identical(getTransSettings(rna.S1), NULL)
  expect_identical(getTransSettings(prot.S1),
                   list(method = "none", suppInfo = "unknown"))
  expect_identical(getTransSettings(meta.S1),
                   list(method = "log2", suppInfo = NULL))
  ### Is it transformed ?
  expect_identical(RFLOMICS:::.isTransformed(rna.S1), FALSE)
  expect_identical(RFLOMICS:::.isTransformed(prot.S1), FALSE)
  expect_identical(RFLOMICS:::.isTransformed(meta.S1), FALSE)

  ### normalization settings
  expect_identical(getNormSettings(rna.S1),
                   list(method = "TMM", suppInfo = NULL))
  expect_identical(getNormSettings(prot.S1),
                   list(method = "none", suppInfo = "unknown"))
  expect_identical(getNormSettings(meta.S1),
                   list(method = "median", suppInfo = NULL))
  ### Is it normalized ?
  expect_identical(RFLOMICS:::.isNormalized(rna.S1), FALSE)
  expect_identical(RFLOMICS:::.isNormalized(prot.S1), FALSE)
  expect_identical(RFLOMICS:::.isNormalized(meta.S1), FALSE)

  ### add PCA for processed data (norm)
  expect(!is.null(getAnalysis(rna.S1,  "PCAlist", "norm")),
         failure_message = "We don't expect this object to be null.")
  expect(!is.null(getAnalysis(prot.S1,  "PCAlist", "norm")),
         failure_message = "We don't expect this object to be null.")
  expect(!is.null(getAnalysis(meta.S1,  "PCAlist", "norm")),
         failure_message = "We don't expect this object to be null.")

  ### reset analysis results
  expect_identical(getAnalysis(rna.S1,  "DiffExpAnal"),  list())
  expect_identical(getAnalysis(prot.S1, "DiffExpAnal"),  list())
  expect_identical(getAnalysis(meta.S1, "DiffExpAnal"),  list())
  expect_identical(getAnalysis(rna.S1,  "DiffExpEnrichAnal"),  list())
  expect_identical(getAnalysis(prot.S1, "DiffExpEnrichAnal"),  list())
  expect_identical(getAnalysis(meta.S1, "DiffExpEnrichAnal"),  list())
  expect_identical(getAnalysis(rna.S1,  "CoExpAnal"),  list())
  expect_identical(getAnalysis(prot.S1, "CoExpAnal"),  list())
  expect_identical(getAnalysis(meta.S1, "CoExpAnal"),  list())
  expect_identical(getAnalysis(rna.S1,  "CoExpEnrichAnal"),  list())
  expect_identical(getAnalysis(prot.S1, "CoExpEnrichAnal"),  list())
  expect_identical(getAnalysis(meta.S1, "CoExpEnrichAnal"),  list())

  ## reset integration analysis
  expect_identical(getAnalysis(MAE1, name = "IntegrationAnalysis"), list())

  expect_identical(getAnalysis(rna.S1,  "DiffExpAnal"),  list())



})

test_that("default values of runDataProcessing arguments", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  # run runDataProcessing with default value
  rna.S1  <- runDataProcessing(rna.S1)
  prot.S1 <- runDataProcessing(prot.S1)
  meta.S1 <- runDataProcessing(meta.S1)

  ### filtering settings: samples
  expect_identical(getSelectedSamples(rna.S1),  colnames(MAE[["RNAtest"]]))
  expect_identical(getSelectedSamples(prot.S1), colnames(MAE[["protetest"]]))
  expect_identical(getSelectedSamples(meta.S1), colnames(MAE[["metatest"]]))

  ### filtering settings: features
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbReplicates", cpmCutoff = 1))
  expect_equal(
    getFilterSettings(prot.S1),
    list(method = "MVI", minValue = 6.226499, suppInfo = "missing value imputation"))
  expect_equal(
    getFilterSettings(meta.S1),
    list(method = "MVI", minValue = 1.98e-05, suppInfo = "missing value imputation"))

  ### transformation settings
  expect_identical(getTransSettings(rna.S1), NULL)
  expect_identical(getTransSettings(prot.S1),
                   list(method = "log2", suppInfo = NULL))
  expect_identical(getTransSettings(meta.S1),
                   list(method = "log2", suppInfo = NULL))

  ### normalization settings
  expect_identical(getNormSettings(rna.S1),
                   list(method = "TMM", suppInfo = NULL))
  expect_identical(getNormSettings(prot.S1),
                   list(method = "median", suppInfo = NULL))
  expect_identical(getNormSettings(meta.S1),
                   list(method = "median", suppInfo = NULL))

})

test_that("Error/warning messages", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]
  sampleToKeep <- colnames(MAE[["RNAtest"]])[-1]

  # run runDataProcessing with incorrect argument
  # RNAseq data
  expect_error(runDataProcessing(rna.S1,
                                 samples = sampleToKeep,
                                 filterStrategy = "toto",
                                 cpmCutoff = 1,
                                 normMethod = "TMM"))
  expect_error(runDataProcessing(rna.S1,
                                 samples = sampleToKeep,
                                 filterStrategy = "NbReplicates",
                                 cpmCutoff = "toto",
                                 normMethod = "TMM"))
  expect_error(runDataProcessing(rna.S1,
                                 samples = sampleToKeep,
                                 filterStrategy = "NbReplicates",
                                 cpmCutoff = 1,
                                 normMethod = "toto"))
  expect_error(runDataProcessing(rna.S1,
                                 samples = sampleToKeep,
                                 filterStrategy = "NbReplicates",
                                 cpmCutoff = 1,
                                 normMethod = "median"))
  expect_error(runDataProcessing(rna.S1,
                                 samples = "toto",
                                 filterStrategy = "NbReplicates",
                                 cpmCutoff = 1,
                                 normMethod = "TMM"))
  expect_warning(runDataProcessing(rna.S1,
                                   samples = sampleToKeep,
                                   filterStrategy = "NbReplicates",
                                   cpmCutoff = 1,
                                   transformMethod = "none",
                                   normMethod = "TMM"))
  # proteomics data
  expect_error(
    runDataProcessing(prot.S1, transformMethod = "toto", normMethod = "none"))
  expect_error(
    runDataProcessing(prot.S1, transformMethod = "none", normMethod = "toto"))
  expect_warning(
    runDataProcessing(prot.S1, filterStrategy = "NbReplicates",
                      transformMethod = "none", normMethod = "none"))
  expect_warning(
    runDataProcessing(prot.S1, cpmCutoff = 1,
                      transformMethod = "none", normMethod = "none"))
})

# ---- runSampleFiltering ----
test_that("runSampleFiltering returned value", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  sampleToKeep <- colnames(MAE[["RNAtest"]])[-1]

  ## method of RflomicsMAE/RflomicsSE class
  MAE1   <- runSampleFiltering(MAE, SE.name = "RNAtest", samples = sampleToKeep)
  rna.S1 <- runSampleFiltering(rna.S1, samples = sampleToKeep)

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(rna.S1))
  expect_identical(rna.S1, MAE1[["RNAtest"]])
  expect_identical(assay(rna.S1), assay(MAE[["RNAtest"]]))
  expect_identical(colData(rna.S1), colData(MAE[["RNAtest"]]))

  ## filtering settings: features
  expect(!identical(metadata(rna.S1)$design$ExpDesign,
                    as.data.frame(colData(rna.S1))),
         failure_message = "")

  ## reset
  rna.S1 <- runDataProcessing(rna.S1)
  expect(!is.null(getFilterSettings(rna.S1)), failure_message = "")
  expect(!is.null(getNormSettings(rna.S1)), failure_message = "")
  rna.S1 <- runSampleFiltering(rna.S1, samples = sampleToKeep)
  expect_identical(getAnalysis(rna.S1, "DataProcessing", "featureFiltering"), list())
  expect_identical(getAnalysis(rna.S1, "DataProcessing", "Normalization"), list())
  expect_identical(getAnalysis(rna.S1, "DataProcessing", "Transformation"),list())
  expect_identical(getAnalysis(rna.S1, "PCAlist")$norm, NULL)

  ##
  rna.S1 <- runSampleFiltering(rna.S1)
  expect_identical(getSelectedSamples(rna.S1), colnames(rna.S1))
  expect_error(runSampleFiltering(rna.S1, samples = colnames(rna.S1)[-1:-3]))

  ##

})

# ---- runFeatureFiltering ----
test_that("runFeatureFiltering returned value", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]

  ## method of RflomicsMAE class
  MAE1   <- runFeatureFiltering(MAE, SE.name = "RNAtest",
                               filterMethod = "CPM",
                               filterStrategy = "NbReplicates",
                               cpmCutoff = 1)
  rna.S1 <- runFeatureFiltering(rna.S1,
                               filterMethod = "CPM",
                               filterStrategy = "NbReplicates",
                               cpmCutoff = 1)

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(rna.S1))
  expect_identical(rna.S1, MAE1[["RNAtest"]])

  ## filtering settings: features
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbReplicates", cpmCutoff = 1))
  expect_identical(RFLOMICS:::.isFiltered(rna.S1), FALSE)
  expect(length(getAnalysis(rna.S1, "DataProcessing", "featureFiltering")) != 0,
         failure_message = "")

  ## reset normalization
  expect_identical(getAnalysis(rna.S1, "DataProcessing", "Normalization"), list())
  expect_identical(getAnalysis(rna.S1, "DataProcessing", "Transformation"),list())
  expect_identical(rna.S1@metadata$DataProcessing$log, NULL)

})

test_that("default values of runFeatureFiltering arguments", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]

  ## method of RflomicsMAE class
  MAE1   <- runFeatureFiltering(MAE, SE.name = "RNAtest")
  rna.S1 <- runFeatureFiltering(rna.S1)

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(rna.S1))
  expect_identical(rna.S1, MAE1[["RNAtest"]])

  ## filtering settings: features
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbReplicates", cpmCutoff = 1))
  expect_identical(RFLOMICS:::.isFiltered(rna.S1), FALSE)
  expect(length(getAnalysis(rna.S1, "DataProcessing", "featureFiltering")) != 0,
         failure_message = "")

  # NULL arg
  rna.S1 <- runFeatureFiltering(rna.S1,
                               filterMethod = NULL,
                               filterStrategy = "NbReplicates",
                               cpmCutoff = 1)
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbReplicates", cpmCutoff = 1))

  rna.S1 <- runFeatureFiltering(rna.S1,
                               filterMethod = "CPM",
                               filterStrategy = NULL,
                               cpmCutoff = 5)
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbReplicates", cpmCutoff = 5))


  rna.S1 <- runFeatureFiltering(rna.S1,
                               filterMethod = "CPM",
                               filterStrategy = "NbConditions",
                               cpmCutoff = NULL)
  expect_identical(
    getFilterSettings(rna.S1),
    list(method = "CPM", filterStrategy = "NbConditions", cpmCutoff = 1))
})

test_that("Error/warning messages", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  # run runDataProcessing with incorrect argument
  # RNAseq data
  expect_error(runFeatureFiltering(rna.S1,
                                  filterMethod = "toto",
                                  filterStrategy = "NbReplicates",
                                  cpmCutoff = 1))
  expect_error(runFeatureFiltering(rna.S1,
                                  filterMethod = "CPM",
                                  filterStrategy = "toto",
                                  cpmCutoff = 1))
  expect_error(runFeatureFiltering(rna.S1,
                                  filterMethod = "CPM",
                                  filterStrategy = "NbConditions",
                                  cpmCutoff = "toto"))

  # proteomics/metabolics data
  ## Can't apply this method to omics types other than RNAseq.
  expect_no_error(runFeatureFiltering(prot.S1))
  expect_warning(runFeatureFiltering(prot.S1, filterMethod = "CMP"))

})

# ---- runTransformation ----
test_that("runTransformation returned value", {

  MAE1    <- MAE
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  ## method of RflomicsMAE class
  MAE1    <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none")

  MAE1    <-
    runFeatureFiltering(MAE1, SE.name = "metatest") |>
    runTransformData(SE.name = "metatest", transformMethod = "log2")

  prot.S1 <-
    runFeatureFiltering(prot.S1) |>
    runTransformData(transformMethod = "none")

  meta.S1 <-
    runFeatureFiltering(meta.S1) |>
    runTransformData(transformMethod = "log2")

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(prot.S1))
  expect_true("RflomicsSE" %in% is(meta.S1))
  expect_identical(prot.S1, MAE1[["protetest"]])
  expect_identical(meta.S1, MAE1[["metatest"]])

  ## transformation settings
  expect_identical(getTransSettings(prot.S1),
                   list(method = "none", suppInfo = "unknown"))
  expect_identical(RFLOMICS:::.isTransformed(prot.S1), FALSE)
  expect(length(getAnalysis(prot.S1, "DataProcessing", "Transformation")) != 0,
         failure_message = "")
  expect_identical(getTransSettings(meta.S1),
                   list(method = "log2", suppInfo = NULL))
  expect_identical(RFLOMICS:::.isTransformed(meta.S1), FALSE)
  expect(length(getAnalysis(meta.S1, "DataProcessing", "Transformation")) != 0,
         failure_message = "")

  ## reset normalization
  expect_identical(getAnalysis(prot.S1, "DataProcessing", "Normalization"), list())
  expect_identical(getAnalysis(meta.S1, "DataProcessing", "Normalization"), list())

})

test_that("default values of runTransformation arguments", {

  MAE1    <- MAE
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  ## method of RflomicsMAE class
  MAE1    <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest")

  MAE1    <-
    runFeatureFiltering(MAE1, SE.name = "metatest") |>
    runTransformData(SE.name = "metatest")

  prot.S1 <-
    runFeatureFiltering(prot.S1) |>
    runTransformData()

  meta.S1 <-
    runFeatureFiltering(meta.S1) |>
    runTransformData()

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(prot.S1))
  expect_true("RflomicsSE" %in% is(meta.S1))
  expect_identical(prot.S1, MAE1[["protetest"]])
  expect_identical(meta.S1, MAE1[["metatest"]])

  ## transformation settings
  expect_identical(getTransSettings(prot.S1),
                   list(method = "log2", suppInfo = NULL))
  expect_identical(RFLOMICS:::.isTransformed(prot.S1), FALSE)
  expect(length(getAnalysis(prot.S1, "DataProcessing", "Transformation")) != 0,
         failure_message = "")
  expect_identical(getTransSettings(meta.S1),
                   list(method = "log2", suppInfo = NULL))
  expect_identical(RFLOMICS:::.isTransformed(meta.S1), FALSE)
  expect(length(getAnalysis(meta.S1, "DataProcessing", "Transformation")) != 0,
         failure_message = "")

  # NULL arg
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  prot.S1 <-
    runFeatureFiltering(prot.S1) |>
    runTransformData(transformMethod = NULL)
  meta.S1 <-
    runFeatureFiltering(meta.S1) |>
    runTransformData(transformMethod = NULL)
  expect_identical(getTransSettings(prot.S1),
                   list(method = "log2", suppInfo = NULL))
  expect_identical(getTransSettings(meta.S1),
                   list(method = "log2", suppInfo = NULL))

})

test_that("Error/warning messages", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  # run with incorrect argument
  # RNAseq data
  expect_error(runTransformData(rna.S1))

  # proteomics/metabolics data
  expect_error(runTransformData(prot.S1, transformMethod = "toto"))
  expect_error(runTransformData(meta.S1, transformMethod = "toto"))

})

# ---- runNormalization ----
test_that("runNormalization returned value", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  ## method of RflomicsMAE class
  MAE1 <-
    runFeatureFiltering(MAE, SE.name = "RNAtest") |>
    runNormalization(SE.name = "RNAtest", normMethod = "TMM")

  MAE1 <-
    runFeatureFiltering(MAE1, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none") |>
    runNormalization(SE.name = "protetest", normMethod = "none")

  MAE1 <-
    runFeatureFiltering(MAE1, SE.name = "metatest") |>
    runTransformData(SE.name = "metatest", transformMethod = "log2") |>
    runNormalization(SE.name = "metatest", normMethod = "median")

  ## method of RflomicsSE class
  rna.S1  <-
    runFeatureFiltering(rna.S1) |>
    runNormalization(normMethod = "TMM")

  prot.S1 <-
    runFeatureFiltering(prot.S1) |>
    runTransformData(transformMethod = "none") |>
    runNormalization(normMethod = "none")

  meta.S1 <-
    runFeatureFiltering(meta.S1) |>
    runTransformData(transformMethod = "log2") |>
    runNormalization(normMethod = "median")

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(rna.S1))
  expect_true("RflomicsSE" %in% is(prot.S1))
  expect_true("RflomicsSE" %in% is(meta.S1))
  expect_identical(rna.S1,  MAE1[["RNAtest"]])
  expect_identical(prot.S1, MAE1[["protetest"]])
  expect_identical(meta.S1, MAE1[["metatest"]])

  ## normalization settings
  expect_identical(getNormSettings(rna.S1),
                   list(method = "TMM", suppInfo = NULL))
  expect_identical(RFLOMICS:::.isNormalized(rna.S1), FALSE)
  expect(length(getAnalysis(rna.S1, "DataProcessing", "Normalization")) != 0,
         failure_message = "")
  expect_identical(getNormSettings(prot.S1),
                   list(method = "none", suppInfo = "unknown"))
  expect_identical(RFLOMICS:::.isNormalized(prot.S1), FALSE)
  expect(length(getAnalysis(prot.S1, "DataProcessing", "Normalization")) != 0,
         failure_message = "")
  expect_identical(getNormSettings(meta.S1),
                   list(method = "median", suppInfo = NULL))
  expect_identical(RFLOMICS:::.isNormalized(meta.S1), FALSE)
  expect(length(getAnalysis(meta.S1, "DataProcessing", "Normalization")) != 0,
         failure_message = "")
})

test_that("default values of runTransformation arguments", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  ## method of RflomicsMAE class
  MAE1 <-
    runFeatureFiltering(MAE, SE.name = "RNAtest") |>
    runNormalization(SE.name = "RNAtest")

  MAE1 <-
    runFeatureFiltering(MAE1, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none") |>
    runNormalization(SE.name = "protetest")

  MAE1 <-
    runFeatureFiltering(MAE1, SE.name = "metatest") |>
    runTransformData(SE.name = "metatest", transformMethod = "log2") |>
    runNormalization(SE.name = "metatest")

  ## method of RflomicsSE class
  rna.S1  <-
    runFeatureFiltering(rna.S1) |>
    runNormalization()

  prot.S1 <-
    runFeatureFiltering(prot.S1) |>
    runTransformData(transformMethod = "none") |>
    runNormalization()

  meta.S1 <-
    runFeatureFiltering(meta.S1) |>
    runTransformData(transformMethod = "log2") |>
    runNormalization()

  ## returns a value of class RflomicsMAE/RflomicsSE
  expect_true("RflomicsMAE" %in% is(MAE1))
  expect_true("RflomicsSE" %in% is(rna.S1))
  expect_true("RflomicsSE" %in% is(prot.S1))
  expect_true("RflomicsSE" %in% is(meta.S1))
  expect_identical(rna.S1,  MAE1[["RNAtest"]])
  expect_identical(prot.S1, MAE1[["protetest"]])
  expect_identical(meta.S1, MAE1[["metatest"]])

  ## normalization settings
  expect_identical(getNormSettings(rna.S1),
                   list(method = "TMM", suppInfo = NULL))
  expect_identical(getNormSettings(prot.S1),
                   list(method = "median", suppInfo = NULL))
  expect_identical(getNormSettings(meta.S1),
                   list(method = "median", suppInfo = NULL))

  # NULL arg
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  rna.S1  <-
    runFeatureFiltering(MAE[["RNAtest"]]) |>
    runNormalization(normMethod = NULL)
  expect_identical(getNormSettings(rna.S1),
                   list(method = "TMM", suppInfo = NULL))

  prot.S1 <-
    runFeatureFiltering(MAE[["protetest"]]) |>
    runTransformData(transformMethod = "none") |>
    runNormalization(normMethod = NULL)
  expect_identical(getNormSettings(prot.S1),
                   list(method = "median", suppInfo = NULL))

  meta.S1 <-
    runFeatureFiltering(MAE[["metatest"]]) |>
    runTransformData(transformMethod = "log2") |>
    runNormalization(normMethod = NULL)
  expect_identical(getNormSettings(meta.S1),
                   list(method = "median", suppInfo = NULL))
})

test_that("Error/warning messages", {

  MAE1    <- MAE
  rna.S1  <- MAE[["RNAtest"]]
  prot.S1 <- MAE[["protetest"]]
  meta.S1 <- MAE[["metatest"]]

  # RNAseq data
  ## RNAseq data should be filtered (low count).
  expect_error(runNormalization(rna.S1))
  expect_error(
    runFeatureFiltering(rna.S1) |> runNormalization(normMethod = "toto"))

  # proteomics/metabolics data
  ## proteomics/metabolics data should be transformed before normalization.
  expect_error(runNormalization(prot.S1))
  expect_error(runNormalization(meta.S1))

  ## arg value not allowed for the parameter normMethod.
  ## Accepted values: median, totalSum, none
  expect_error(
    runTransformData(prot.S1, transformMethod = "none") |>
      runNormalization(normMethod = "toto"))
  expect_error(
    runTransformData(meta.S1, transformMethod = "log2") |>
      runNormalization(normMethod = "toto"))

})

test_that("Test explor plot", {
  
  sampleToKeep <- colnames(MAE[["RNAtest"]])[-1]
  MAE1 <- runDataProcessing(MAE, SE.name = "RNAtest",
                            samples = sampleToKeep,
                            filterStrategy = "NbReplicates",
                            cpmCutoff = 1,
                            normMethod = "TMM")
  MAE1 <- runDataProcessing(MAE1, SE.name = "protetest",
                            transformMethod = "none",
                            normMethod = "none",
                            imputMethod = "MVI")
  
  p <- plotLibrarySize(MAE1, SE.name = "RNAtest", raw = TRUE)
  expect_equal(is(p), "gg")
  expect_error(plotLibrarySize(MAE1, SE.name = "protetest"))
  
  p <- plotDataDistribution(MAE1, SE.name = "RNAtest", plot = "boxplot")
  expect_equal(is(p), "gg")
  
  p <- plotDataDistribution(MAE1, SE.name = "protetest", plot = "density")
  expect_equal(is(p), "gg")
  
  p <- plotOmicsPCA(MAE1, SE.name = "RNAtest")
  expect_equal(is(p), "gg")
  
  p <- plotExpDesignCompleteness(MAE1, omicName = "RNAtest")
  expect_equal(is(p), "gg")
})