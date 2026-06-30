### ============================================================================
### [03_runDataProcessing]
### ----------------------------------------------------------------------------
# N. Bessoltane

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

rna.se  <- MAE[["RNAtest"]]
meta.se <- MAE[["metatest"]]

test_that("runSampleFiltering: sample filtering", {
  
  # samples = NULL
  rna.se1  <- 
    runSampleFiltering(
      rna.se, samples = NULL
    )
  
  meta.se1  <- 
    runSampleFiltering(
      meta.se, samples = NULL
    )
  
  expect_identical(rna.se1, rna.se)
  expect_identical(meta.se1, meta.se)
  
  rna.se2  <- runSampleFiltering(rna.se)
  meta.se2 <- runSampleFiltering(meta.se)
  
  expect_identical(rna.se1, rna.se2)
  expect_identical(meta.se1, meta.se2)
  
  MAE1 <- MAE |>
    runSampleFiltering(
      SE.name = "RNAtest", 
      samples = NULL) |>
    runSampleFiltering(
      SE.name = "metatest", 
      samples = NULL)
  
  expect_identical(rna.se1, MAE1[["RNAtest"]])
  expect_identical(meta.se1, MAE1[["metatest"]])
  
  # samples = c("toto", "titi", "tata")
  expect_error(
    runSampleFiltering(
      rna.se, samples = c("toto", "titi", "tata")
    )
  )
  
  # valide values for "samples"
  rna.se1  <- 
    runSampleFiltering(
      rna.se, samples = colnames(rna.se)[-1]
    )
  
  meta.se1  <- 
    runSampleFiltering(
      meta.se, samples = colnames(meta.se)[-2]
    )
  
  expect_false(identical(assay(rna.se1), assay(rna.se)))
  expect_false(identical(assay(meta.se1), assay(meta.se)))
  
  expect_identical(colnames(rna.se1), colnames(rna.se)[-1])
  expect_identical(colnames(meta.se1), colnames(meta.se)[-2])
  
  expect_identical(getSelectedSamples(rna.se1), colnames(rna.se)[-1])
  expect_identical(getSelectedSamples(meta.se1), colnames(meta.se)[-2])
  
  # completness
  expect_error(
    runSampleFiltering(
      rna.se, samples = colnames(rna.se)[c(-1, -2)]
    )
  )
  
  expect_no_error(
    runSampleFiltering(
      rna.se, samples = colnames(rna.se)[c(-1, -4)]
    )
  )
  
  # compare with runSampleFiltering
  rna.se2  <- 
    runDataProcessing(
      rna.se, samples = colnames(rna.se)[-1]
    )
  
  meta.se2  <- 
    runDataProcessing(
      meta.se, samples = colnames(meta.se)[-2], 
    )
  
  expect_identical(assay(rna.se1), assay(rna.se2))
  expect_identical(assay(meta.se1), assay(meta.se2))
  
})

test_that("runFeatureFiltering: feature filtering", {
  
  # filtering -> NULL / default values
  expect_warning(
    rna.se1  <- 
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = 
        list(method           = NULL,
             strategy         = NULL,
             cpmCutoff        = NULL)
    )
    )
  
  expect_warning(
  meta.se1  <- 
    runFeatureFiltering(
      meta.se, 
      missingValueFilter = 
        list(method           = NULL,
             proportion       = NULL,
             nbCondition      = NULL)
    )
  )
  
  expect_identical(rna.se1, rna.se)
  expect_identical(meta.se1, meta.se)
  expect_false(RFLOMICS:::.isFiltered(rna.se1))
  expect_false(RFLOMICS:::.isFiltered(meta.se1))
  expect_null(getFilterSettings(rna.se1))
  expect_null(getFilterSettings(meta.se1))
  
  rna.se2  <- runFeatureFiltering(rna.se)
  meta.se2 <- runFeatureFiltering(meta.se)
  
  expect_false(identical(rna.se1, rna.se2))
  expect_false(identical(meta.se1, meta.se2))
  
  expect_warning(
    MAE1 <- MAE |>
      runFeatureFiltering(
        SE.name = "RNAtest", 
        lowCountFilter = 
          list(method           = NULL,
               strategy         = NULL,
               cpmCutoff        = NULL)
      ) |>
      runFeatureFiltering(
        SE.name = "metatest", 
        missingValueFilter = 
          list(method           = NULL,
               proportion       = NULL,
               nbCondition      = NULL)
      )
  )
  
  expect_identical(rna.se1, MAE1[["RNAtest"]])
  expect_identical(meta.se1, MAE1[["metatest"]])
  
  # filter -> "none"
  rna.se1  <- 
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "none")
    )
  
  meta.se1  <- 
    runFeatureFiltering(
      meta.se,
      missingValueFilter = list(method = "none")
    )
  
  expect_identical(assay(rna.se1), assay(rna.se))
  expect_identical(assay(meta.se1), assay(meta.se))
  expect_true(RFLOMICS:::.isFiltered(rna.se1))
  expect_true(RFLOMICS:::.isFiltered(meta.se1))
  expect_identical(getFilterSettings(rna.se1)$method, "none")
  expect_identical(getFilterSettings(meta.se1)$method, "none")

  # 
  expect_error(
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "toto")
    )
  )
  expect_error(
    runFeatureFiltering(
      meta.se, 
      missingValueFilter = list(method = "toto")
    )
  )
  
  # default values
  rna.se1  <- 
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "CPM")
    )
  expect_false(identical(assay(rna.se1), assay(rna.se)))
  expect_identical(getFilterSettings(rna.se1)$method, "CPM")
  expect_false(identical(getFilterSettings(rna.se1)$strategy, NULL))
  expect_false(identical(getFilterSettings(rna.se1)$cpmCutoff, NULL))
  expect_true(RFLOMICS:::.isFiltered(rna.se1))
  
  rna.se1  <- 
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "filterByExpr")
    )
  expect_false(identical(assay(rna.se1), assay(rna.se)))
  expect_identical(getFilterSettings(rna.se1)$method, "filterByExpr")
  expect_false(identical(getFilterSettings(rna.se1)$strategy, NULL))
  expect_true(identical(getFilterSettings(rna.se1)$cpmCutoff, NULL))
  
  # errors
  expect_error(
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "filterByExpr", strategy = "toto")
    )
  )
  expect_error(
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "CPM", cpmCutoff = -1)
    )
  )
  
  # ignore args
  expect_no_error(
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "filterByExpr", cpmCutoff = 1)
    )
  )
  expect_false(identical(getFilterSettings(rna.se1)$strategy, NULL))
  expect_true(identical(getFilterSettings(rna.se1)$cpmCutoff, NULL))
  
  # no missing values
  meta.se1  <- 
    runFeatureFiltering(
      meta.se,
      missingValueFilter = list(method = "GlobalFiltering")
    )
  expect_identical(assay(meta.se1), assay(meta.se))
  expect_identical(getFilterSettings(meta.se1)$method, "none")
  expect_true(RFLOMICS:::.isFiltered(meta.se1))
  
  # case with missing values
  # ?
  
  # compare with runFeatureFiltering
  rna.se1 <- 
    runFeatureFiltering(
      rna.se, 
      lowCountFilter = list(method = "filterByExpr")
    )
  
  rna.se2 <- 
    runDataProcessing(
      rna.se, 
      lowCountFilter = list(method = "filterByExpr")
    )
  
  expect_true(is(rna.se2, "RflomicsSE"))
  expect_identical(assay(rna.se1), assay(rna.se2))
  
})

test_that("runTransformData: transformation", {
  
  # trasnformation -> NULL / default values
  rna.se1  <- 
    runTransformData(
      rna.se, 
      method = NULL
    )
  
  meta.se1  <- 
    runTransformData(
      meta.se, 
      method = NULL
    )
  
  expect_true(is(rna.se1, "RflomicsSE"))
  expect_true(is(meta.se1, "RflomicsSE"))
  expect_false(identical(rna.se1, rna.se))
  expect_false(identical(meta.se1, meta.se))
  expect_true(RFLOMICS:::.isTransformed(rna.se1))
  expect_true(RFLOMICS:::.isTransformed(meta.se1))
  expect_false(identical(getTransSettings(rna.se1), NULL))
  expect_false(identical(getTransSettings(meta.se1), NULL))
  
  rna.se2  <- runTransformData(rna.se)
  meta.se2 <- runTransformData(meta.se)
  
  expect_identical(rna.se1, rna.se2)
  expect_identical(meta.se1, meta.se2)
  
  MAE1 <- MAE |>
    runTransformData(
      SE.name = "RNAtest", 
      method = NULL
    ) |>
    runTransformData(
      SE.name = "metatest", 
      method = NULL
    )
  
  expect_true(is(MAE1[["RNAtest"]], "RflomicsSE"))
  expect_true(is(MAE1[["metatest"]], "RflomicsSE"))
  expect_true(is(MAE1, "RflomicsMAE"))
  
  expect_identical(rna.se1, MAE1[["RNAtest"]])
  expect_identical(meta.se1, MAE1[["metatest"]])
  
  # filter -> "none"
  expect_no_error(
    runTransformData(
      rna.se, 
      method = "none"
    )
  )
  
  rna.se1  <- 
    runTransformData(
      rna.se,
      method = "none"
    )
  
  meta.se1  <- 
    runTransformData(
      meta.se,
      method = "none"
    )
  
  expect_identical(assay(rna.se1), assay(rna.se))
  expect_identical(assay(meta.se1), assay(meta.se))
  expect_true(RFLOMICS:::.isTransformed(rna.se1))
  expect_true(RFLOMICS:::.isTransformed(meta.se1))
  expect_identical(getTransSettings(rna.se1)$method, "none")
  expect_identical(getTransSettings(meta.se1)$method, "none")
  
  # 
  expect_error(
    runTransformData(
      rna.se, 
      method = "log10"
    )
  )
  expect_error(
    runTransformData(
      meta.se, 
      method = "toto"
    )
  )
  
  # 
  meta.se1  <- 
    runTransformData(
      meta.se, 
      method = "log10"
    )
  expect_false(identical(assay(meta.se1), assay(meta.se)))
  expect_identical(getTransSettings(meta.se1)$method, "log10")
  expect_true(RFLOMICS:::.isTransformed(meta.se1))
  
  meta.se1  <- 
    runTransformData(
      meta.se, 
      method = "log2"
    )
  expect_false(identical(assay(meta.se1), assay(meta.se)))
  expect_identical(getTransSettings(meta.se1)$method, "log2")
  expect_true(RFLOMICS:::.isTransformed(meta.se1))
  
  # case with missing values 
  # ?
  
  # comapre with runTransformData
  meta.se1  <- 
    runDataProcessing(
      meta.se, 
      transform = list(method = "log2")
    )
  
  meta.se2  <- 
    runTransformData(
      meta.se, 
      method = "log2"
    )
  
  expect_true(is(meta.se2, "RflomicsSE"))
  expect_identical(assay(meta.se1), assay(meta.se2))
  
  # data already transformed 
  expect_error(
    runTransformData(
      meta.se2, 
      method = "log2"
    )
  )
  
})

test_that("runNormalization: normalization", {
  
  # normalization -> NULL / default values
  rna.se1  <- 
    runNormalization(
      rna.se, 
      method = NULL
    )
  
  meta.se1  <- 
    runNormalization(
      meta.se, 
      method = NULL
    )
  
  expect_true(is(rna.se1, "RflomicsSE"))
  expect_true(is(meta.se1, "RflomicsSE"))
  expect_false(identical(rna.se1, rna.se))
  expect_false(identical(meta.se1, meta.se))
  expect_true(RFLOMICS:::.isNormalized(rna.se1))
  expect_true(RFLOMICS:::.isNormalized(meta.se1))
  expect_false(identical(getNormSettings(rna.se1), NULL))
  expect_false(identical(getNormSettings(meta.se1), NULL))
  
  rna.se2  <- runNormalization(rna.se)
  meta.se2 <- runNormalization(meta.se)
  
  expect_identical(rna.se1, rna.se2)
  expect_identical(meta.se1, meta.se2)
  
  MAE1 <- MAE |>
    runNormalization(
      SE.name = "RNAtest", 
      method = NULL
    ) |>
    runNormalization(
      SE.name = "metatest", 
      method = NULL
    )
  
  expect_true(is(MAE1[["RNAtest"]], "RflomicsSE"))
  expect_true(is(MAE1[["metatest"]], "RflomicsSE"))
  expect_true(is(MAE1, "RflomicsMAE"))
  
  expect_identical(rna.se1, MAE1[["RNAtest"]])
  expect_identical(meta.se1, MAE1[["metatest"]])
  
  # filter -> "none"
  rna.se1  <- 
    runNormalization(
      rna.se, 
      method = "none"
    )
  
  meta.se1  <- 
    runNormalization(
      meta.se,
      method = "none"
    )
  
  expect_identical(assay(rna.se1), assay(rna.se))
  expect_identical(assay(meta.se1), assay(meta.se))
  expect_true(RFLOMICS:::.isNormalized(rna.se1))
  expect_true(RFLOMICS:::.isNormalized(meta.se1))
  expect_identical(getNormSettings(rna.se1)$method, "none")
  expect_identical(getNormSettings(meta.se1)$method, "none")
  
  # 
  expect_error(
    runNormalization(
      meta.se, 
      method = "toto"
    )
  )
  
  # 
  rna.se1  <- 
    runNormalization(
      rna.se, 
      method = "TMM"
    )
  expect_false(identical(assay(rna.se1), assay(rna.se)))
  expect_identical(getNormSettings(rna.se1)$method, "TMM")
  expect_true(RFLOMICS:::.isNormalized(rna.se1))
  
  meta.se1  <- #  “median”, “totalSum”, “none”
    runNormalization(
      meta.se, 
      method = "median"
    )
  expect_false(identical(assay(meta.se1), assay(meta.se)))
  expect_identical(getNormSettings(meta.se1)$method, "median")
  expect_true(RFLOMICS:::.isNormalized(meta.se1))
  
  meta.se1  <- #  “median”, “totalSum”, “none”
    runNormalization(
      meta.se, 
      method = "totalSum"
    )
  expect_false(identical(assay(meta.se1), assay(meta.se)))
  expect_identical(getNormSettings(meta.se1)$method, "totalSum")
  expect_true(RFLOMICS:::.isNormalized(meta.se1))
  
  # comapre with runTransformData
  meta.se1  <- 
    runDataProcessing(
      meta.se, 
      normalize = list(method = "median")
    )
  
  meta.se2  <- 
    runNormalization(
      meta.se, 
      method = "median"
    )
  
  expect_true(is(meta.se2, "RflomicsSE"))
  expect_identical(assay(meta.se1), assay(meta.se2))
  
  # data already transformed 
  expect_error(
    runNormalization(
      meta.se2, 
      method = "median"
    )
  )
  
})

#test_that("runImputation imputation", { })

test_that("runOmicsPCA", { 
  
  rna.se1  <- 
    runOmicsPCA(
      rna.se, ncomp = 3
    )
  
  expect_true(is(rna.se1, "RflomicsSE"))
  expect_false(identical(rna.se1, rna.se))

  MAE  <- 
    runOmicsPCA(
      MAE, SE.name = "RNAtest", ncomp = 3
    )
  expect_identical(MAE[["RNAtest"]], rna.se1)
  
  rna.se2  <- 
    runDataProcessing(
      rna.se, ncomp = 3
    )

  expect_identical(metadata(rna.se2)[["PCA"]], metadata(rna.se1)[["PCA"]])
  
})

test_that("runDataProcessing", { 
  
  rna.se1  <- 
    runDataProcessing(
      rna.se, 
    )
  
  expect_identical(assay(rna.se1), assay(rna.se))
  expect_false(RFLOMICS:::.isFiltered(rna.se1))
  expect_false(RFLOMICS:::.isNormalized(rna.se1))
  expect_false(RFLOMICS:::.isTransformed(rna.se1))
  expect_false(RFLOMICS:::.isImputed(rna.se1))

  meta.se1  <- 
    runDataProcessing(
      meta.se, 
    )
  
  expect_identical(assay(meta.se1), assay(meta.se))
  expect_false(RFLOMICS:::.isFiltered(meta.se1))
  expect_false(RFLOMICS:::.isNormalized(meta.se1))
  expect_false(RFLOMICS:::.isTransformed(meta.se1))
  expect_false(RFLOMICS:::.isImputed(meta.se1))

})

# ---- visualisation ----
test_that("Test explor plot", {

  MAE1 <- MAE |>
    runDataProcessing(
      SE.name = "RNAtest", 
      lowCountFilter = list(method = "CPM"),
      normalize = list(method = "TMM")
    ) |>
    runDataProcessing(
      SE.name = "metatest", 
      #missingValueFilter = list(method = NULL),
      transform = list(method = "log2"),
      normalize = list(method = "median")
    )

  p <- plotLibrarySize(MAE1, SE.name = "RNAtest", raw = TRUE)
  expect(is(p, "gg"), "This plot is not ggplot")
  expect_error(plotLibrarySize(MAE1, SE.name = "metatest"))

  p <- plotDataDistribution(MAE1, SE.name = "RNAtest", plot = "boxplot")
  expect(is(p, "gg"), "This plot is not ggplot")

  p <- plotDataDistribution(MAE1, SE.name = "metatest", plot = "density")
  expect(is(p, "gg"), "This plot is not ggplot")

  p <- plotOmicsPCA(MAE1, SE.name = "RNAtest")
  expect(is(p, "gg"), "This plot is not ggplot")

  p <- plotExpDesignCompleteness(MAE1, omicName = "RNAtest")
  expect(is(p, "gg"), "This plot is not ggplot")
  
  p <- plotOmicsPCA(MAE1, SE.name = "metatest")
  expect(is(p, "gg"), "This plot is not ggplot") 
})


#checkExpDesignCompleteness

#miniRflomicsSE

