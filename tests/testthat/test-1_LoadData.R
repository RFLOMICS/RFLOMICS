library(RFLOMICS)
library(testthat)

# ---- Tests RFLOMICS constructor ----

# load data
# load ecoseed data
data(ecoseed)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorRef"    = c("rep1", "Low", "DS"),
  "factorType"   = c("batch", "Bio", "Bio"),
  "factorLevels" = c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI")
)

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = list(ecoseed.df$RNAtest, ecoseed.df$protetest, ecoseed.df$metatest),
  omicsNames  = c("RNAtest", "protetest", "metatest"),
  omicsTypes  = c("RNAseq","proteomics","metabolomics"),
  ExpDesign   = ecoseed.df$design,
  factorInfo   = factorInfo)
names(MAE) <- c("RNAtest", "protetest", "metatest")

# ---- test output class ----
test_that("FlomicsMultiAssay.constructor fonction return MultiAssayExperiment object", {
  
  # test type if MAE class
  expect_true("MultiAssayExperiment" %in% is(MAE))
  
  # test type if element of MAE class
  for (SE in names(MAE)) {
    expect_true("SummarizedExperiment" %in% is(MAE[[SE]]))
  }
})

test_that("test if RflomicsMAE / RflomicsSE", {
  expect_true(is(MAE, "RflomicsMAE"))
})

## ---- if input is mae object ----
test_that("test mae/se input", {
  
  MAE2 <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests",
    omicsData   = ecoseed.mae,
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    factorInfo   = factorInfo)
  names(MAE) <- c("RNAtest", "metatest", "protetest")
  
  expect_true(is(MAE2, "RflomicsMAE"))
  for (SE in names(MAE2)) {
    expect_true("SummarizedExperiment" %in% is(MAE2[[SE]]))
  }
  
})

## ---- if input is list of se ----
test_that("test mae/se input", {
  
  omicsData3 <- list(
    RNAtest = ecoseed.mae[["RNAtest"]], 
    metatest = ecoseed.mae[["metatest"]], 
    protetest = ecoseed.mae[["protetest"]]
  )
  
  factorInfo3 <- factorInfo[,c(-2, -4)]
  
  MAE3 <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests",
    omicsData   = omicsData3,
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = as.data.frame(colData(ecoseed.mae)),
    factorInfo  = factorInfo3)
  names(MAE) <- c("RNAtest", "metatest", "protetest")
  
  expect_true(is(MAE3, "RflomicsMAE"))
  for (SE in names(MAE3)) {
    expect_true("SummarizedExperiment" %in% is(MAE3[[SE]]))
  }
  
})

## ---- if input is list of se & data.frame ----
test_that("test mae/se input", {
  
  omicsData4 <- list(
    RNAtest = ecoseed.mae[["RNAtest"]], 
    metatest = ecoseed.df$metatest, 
    protetest = as.matrix(ecoseed.df$protetest)
  )
  
  MAE4 <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests",
    omicsData   = omicsData4,
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = ecoseed.df$design,
    factorInfo  = factorInfo)
  names(MAE) <- c("RNAtest", "metatest", "protetest")
  
  expect_true(is(MAE4, "RflomicsMAE"))
  for (SE in names(MAE4)) {
    expect_true("SummarizedExperiment" %in% is(MAE4[[SE]]))
  }
  
})


# ---- test getters ----
test_that("test getDatasetNames", {
  expect_true(all(getDatasetNames(MAE) %in% c("RNAtest", "metatest", "protetest")))
})

test_that("Design accessors", {
  expect_identical(getFactorTypes(MAE), metadata(MAE)$design$Factors.Type)
  
})

test_that("Factors types", {
  
  expect_identical(getBioFactors(MAE), c("temperature", "imbibition"))
  expect_identical(getBatchFactors(MAE), c("Repeat"))
  
})

# Test of internal function
test_that("Omics dictionnary", {
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "RNAtest"), 
                   list(variableName = "transcripts", valueType = "counts"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "metatest"), 
                   list(variableName = "metabolites", valueType = "XIC"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "protetest"), 
                   list(variableName = "proteins", valueType = "XIC"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "metatest"), 
                   RFLOMICS:::omicsDic(MAE[["metatest"]]))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "RNAtest"), 
                   RFLOMICS:::omicsDic(MAE[["RNAtest"]]))
  
  expect_error(RFLOMICS:::omicsDic(MAE))
  
})


# ---- Test order of samples ----
test_that("All omics data are ordred in same way", {
  
  expect_equal(colnames(MAE[[1]]), colnames(MAE[[2]]))
  expect_equal(colnames(MAE[[3]]), colnames(MAE[[2]]))
  expect_equal(colnames(MAE[[3]]), colnames(MAE[[2]]))
})

test_that("Test if samples in data matrix and rownames in design are ordered in same way", {
  
  expect_equal(colnames(MAE[[1]]), as.character(MAE[[1]]$samples))
  expect_equal(colnames(MAE[[2]]), as.character(MAE[[2]]$samples))
  expect_equal(colnames(MAE[[3]]), as.character(MAE[[3]]$samples))
})

# ---- Tests order of samples when samples are not all the same ----
test_that("Test if samples in data matrix and rownames in design are orderd in same way", {
  
  omicsData <- list(
    ecoseed.df$RNAtest,
    ecoseed.df$metatest,
    ecoseed.df$protetest)
  
  omicsData[[1]] <- omicsData[[1]][,-5]
  omicsData[[2]] <- omicsData[[2]][,-10]
  ExpDesign      <- ecoseed.df$design[-20, ]
  
  MAE <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests", 
    omicsData   = omicsData,
    omicsNames  = c("RNAtest", "metatest", "protetest"),
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = ExpDesign,
    factorInfo   = factorInfo)
  
  expect_equal(colnames(MAE[[1]]), as.vector(MAE[[1]]$samples))
  expect_equal(colnames(MAE[[2]]), as.vector(MAE[[2]]$samples))
  expect_equal(colnames(MAE[[3]]), as.vector(MAE[[3]]$samples))
})

