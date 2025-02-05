### ============================================================================
### [01_CreateRflomicsMAE] 
### ----------------------------------------------------------------------------
# N. Bessoltane, 

library(RFLOMICS)
library(testthat)

# ---- Tests RFLOMICS constructor ----

# load data
data(ecoseed.df)
data(ecoseed.mae)

# information about factors
factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorRef"    = c("rep1", "Low", "DS"),
  "factorType"   = c("batch", "Bio", "Bio"),
  "factorLevels" = c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI")
)

# create rflomicsMAE object
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = list(RNAtest = ecoseed.df$RNAtest, 
                     protetest = ecoseed.df$protetest, 
                     metatest = ecoseed.df$metatest),
  omicsTypes  = c("RNAseq","proteomics","metabolomics"),
  ExpDesign   = ecoseed.df$design,
  factorInfo  = factorInfo)

# ---- test output class ----
test_that("createRflomicsMAE returns RflomicsMAE/MultiAssayExperiment object", {
  
  # test type if MAE class
  expect_true("RflomicsMAE" %in% is(MAE))
  expect_true("MultiAssayExperiment" %in% is(MAE))
  
  # test type if element of MAE class
  for (SE in names(MAE)) {
    expect_true("RflomicsSE" %in% is(MAE[[SE]]))
    expect_true("SummarizedExperiment" %in% is(MAE[[SE]]))
  }
})

test_that("test MAE project name.", {
  
  # test type if MAE class
  expect_identical(getProjectName(MAE), "Tests")
})

test_that("test MAE names.", {
  
  # test type if MAE class
  expect_identical(names(MAE), c("RNAtest", "protetest", "metatest"))
})

test_that("test MAE omics types", {
  
  # test type if MAE class
  omicsType <- c("RNAseq","proteomics","metabolomics")
  names(omicsType) <- c("RNAtest", "protetest", "metatest")
  
  expect_identical(getOmicsTypes(MAE), omicsType)
})

## ---- colData ----
test_that("colData", {
  
  Repeat      <- factor(rep(c("rep1", "rep2", "rep3"), 9), 
                        levels =c("rep1", "rep2", "rep3"))
  imbibition  <- factor(c(rep("DS", 9), rep("EI", 9), rep("LI", 9)), 
                        levels =c("DS", "EI", "LI"))    
  temperature <- factor(rep(c(rep("Low",3), rep("Medium", 3), rep("Elevated", 3)), 3),  
                        levels =c("Low", "Medium", "Elevated"))
  Repeat      <- relevel(as.factor(Repeat),      ref="rep1")
  temperature <- relevel(as.factor(temperature), ref="Low")
  imbibition  <- relevel(as.factor(imbibition),  ref="DS")
  
  # groups
  groups.level <- character(0)
  # Boucles imbriquées pour créer les combinaisons
  for (v1 in c("Low", "Medium", "Elevated")) {
    for (v2 in c("DS", "EI", "LI")) {
      groups.level <- c(groups.level, paste(v1, v2, sep = "_"))
    }
  }
  groups <- factor(paste(temperature, imbibition, sep = "_"),
                   levels =groups.level)
  
  # samples
  samples.level <- character(0)
  # Boucles imbriquées pour créer les combinaisons
  for (v1 in groups.level) {
    for (v2 in 1:3) {
      samples.level <- c(samples.level, paste(v1, v2, sep = "_"))
    }
  }
  
  samples <- factor(paste(temperature, imbibition, sub("rep", "", Repeat), sep = "_"),
                    levels = samples.level)
  
  
  colData <- data.frame(Repeat      = Repeat, 
                        groups      = groups,
                        temperature = temperature,
                        imbibition  = imbibition,
                        samples     = samples)
  
  rownames(colData) <- samples
  
  expect_equal(as.data.frame(getDesignMat(MAE)), as.data.frame(colData))
})


test_that("test MAE metadtata", {
  
  # design
  expect_identical(getModelFormula(MAE), logical(0))
  expect_identical(getSelectedContrasts(MAE), data.frame())
  expect_identical(getFactorNames(MAE), factorInfo$factorName)
  
  # IntegrationAnalysis
  # expect_identical(getAnalysis(MAE, name = "IntegrationAnalysis"), list()) !!!!!!!
})


test_that("test SE metadtata", {
  
  for (SE in names(MAE)) {
    
    # design
    # expect_identical(getModelFormula(MAE[[SE]]), logical(0)) !!!!!!!!
    # expect_identical(getSelectedContrasts(MAE[[SE]]), data.frame()) !!!!!!!!
    expect_identical(getFactorNames(MAE[[SE]]), factorInfo$factorName)
    
    # DataProcessing
    if(getOmicsTypes(MAE[[SE]]) == "RNAseq"){
      expect(!is.null(MAE[[SE]]@metadata$DataProcessing$rowSumsZero), 
             failure_message = "This value should not be 0.")
    }
    else{
      expect_equal(MAE[[SE]]@metadata$DataProcessing$rowSumsZero, NULL)
    }
    expect_equal(as.vector(getSelectedSamples(MAE[[SE]])), colnames(MAE[[SE]]))
    expect_equal(getAnalysis(object  = MAE[[SE]], 
                             name    = "DataProcessing", 
                             subName = "featureFiltering"), 
                 list())
    expect_equal(getAnalysis(object  = MAE[[SE]], 
                             name    = "DataProcessing", 
                             subName = "Normalization"), 
                 list())
    expect_equal(getAnalysis(object  = MAE[[SE]], 
                             name    = "DataProcessing", 
                             subName = "Normalization"),
                 list())
    expect_equal(getAnalysis(object  = MAE[[SE]], 
                             name    = "DataProcessing", 
                             subName = "log"),
                 NULL)
  }
  
  # PCAlist
  expect_no_error(getAnalysis(MAE[[SE]], name = "PCAlist", subName = "raw"))
  expect_null(getAnalysis(MAE[[SE]], name = "PCAlist", subName = "norm"))
  
  # Analyses
  expect_equal(getAnalysis(MAE[[SE]], name = "DiffExpAnal"),       list())
  expect_equal(getAnalysis(MAE[[SE]], name = "CoExpAnal"),         list())
  expect_equal(getAnalysis(MAE[[SE]], name = "DiffExpEnrichAnal"), list())
  expect_equal(getAnalysis(MAE[[SE]], name = "CoExpEnrichAnal"),   list())
})


## ---- if input is mae object ----
test_that("test mae/se input", {
  
  MAE2 <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests",
    omicsData   = ecoseed.mae,
    omicsTypes  = c("RNAseq","proteomics","metabolomics"),
    factorInfo   = factorInfo)

  expect_true(is(MAE2, "RflomicsMAE"))
  for (SE in names(MAE2)) {
    expect_true("SummarizedExperiment" %in% is(MAE2[[SE]]))
  }

  expect_equal(colData(MAE), colData(MAE2))
  for(SE in names(MAE2)){
    expect_equal(MAE[[SE]], MAE2[[SE]])
  }
})

## ---- if input is list of se ----
test_that("test mae/se input", {
  
  omicsData3 <- list(
    RNAtest = ecoseed.mae[["RNAtest"]], 
    protetest = ecoseed.mae[["protetest"]],
    metatest = ecoseed.mae[["metatest"]]
  )
  
  factorInfo3 <- factorInfo[,c(-2, -4)]
  
  MAE3 <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests",
    omicsData   = omicsData3,
    omicsTypes  = c("RNAseq","proteomics","metabolomics"),
    ExpDesign   = as.data.frame(colData(ecoseed.mae)),
    factorInfo  = factorInfo3)

  expect_true(is(MAE3, "RflomicsMAE"))
  for (SE in names(MAE3)) {
    expect_true("SummarizedExperiment" %in% is(MAE3[[SE]]))
  }
  
  expect_equal(colData(MAE), colData(MAE3))
  for(SE in names(MAE3)){
    expect_equal(MAE[[SE]], MAE3[[SE]])
  }
})

## ---- if input is list of se & data.frame ----
test_that("test mae/se input", {
  
  omicsData4 <- list(
    RNAtest = ecoseed.mae[["RNAtest"]], 
    protetest = as.matrix(ecoseed.df$protetest),
    metatest = ecoseed.df$metatest
  )
  
  MAE4 <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests",
    omicsData   = omicsData4,
    omicsTypes  = c("RNAseq","proteomics","metabolomics"),
    ExpDesign   = ecoseed.df$design,
    factorInfo  = factorInfo)

  expect_true(is(MAE4, "RflomicsMAE"))
  for (SE in names(MAE4)) {
    expect_true("SummarizedExperiment" %in% is(MAE4[[SE]]))
  }
  
  expect_equal(colData(MAE), colData(MAE4))
  for(SE in names(MAE4)){
    expect_equal(MAE[[SE]], MAE4[[SE]])
  }
})


# ---- test getters ----
test_that("test getDatasetNames", {
  expect_true(all(getDatasetNames(MAE) %in% c("RNAtest", "metatest", "protetest")))
})

test_that("test getFactorTypes", {
  
  expect_identical(getFactorTypes(MAE), metadata(MAE)$design$Factors.Type)
  
  vec <- c("batch", "Bio", "Bio")
  names(vec) <- c("Repeat", "temperature", "imbibition" )
  expect_equal(getFactorTypes(MAE), vec)
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getFactorTypes(SE), vec)
})

test_that("Factors types", {
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  
  expect_identical(getBioFactors(MAE), c("temperature", "imbibition"))
  expect_equal(getBioFactors(SE), c("temperature", "imbibition"))
  
  expect_identical(getBatchFactors(MAE), c("Repeat"))
  expect_equal(getBatchFactors(SE), c("Repeat"))
  
  expect_null(getMetaFactors(MAE))
  expect_null(getMetaFactors(SE))
  
})

test_that("test getFactorModalities", {
  
  expect_equal(getFactorModalities(MAE, "imbibition"), c("DS", "EI", "LI"))
  expect_equal(getFactorModalities(MAE[["protetest"]], "imbibition"),  c("DS", "EI", "LI"))
})

# Test of internal function
test_that("Omics dictionnary", {
  expect_identical(RFLOMICS:::.omicsDic(MAE, SE.name = "RNAtest"), 
                   list(variableName = "transcript", valueType = "counts"))
  expect_identical(RFLOMICS:::.omicsDic(MAE, SE.name = "metatest"), 
                   list(variableName = "metabolite", valueType = "XIC"))
  expect_identical(RFLOMICS:::.omicsDic(MAE, SE.name = "protetest"), 
                   list(variableName = "protein", valueType = "XIC"))
  expect_identical(RFLOMICS:::.omicsDic(MAE, SE.name = "metatest"), 
                   RFLOMICS:::.omicsDic(MAE[["metatest"]]))
  expect_identical(RFLOMICS:::.omicsDic(MAE, SE.name = "RNAtest"), 
                   RFLOMICS:::.omicsDic(MAE[["RNAtest"]]))
  
  expect_error(RFLOMICS:::.omicsDic(MAE))
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
    RNAtest   = ecoseed.df$RNAtest,
    metatest  = ecoseed.df$metatest,
    protetest = ecoseed.df$protetest)
  
  omicsData[[1]] <- omicsData[[1]][,-5]
  omicsData[[2]] <- omicsData[[2]][,-10]
  ExpDesign      <- ecoseed.df$design[-20, ]
  
  MAE <- RFLOMICS::createRflomicsMAE(
    projectName = "Tests", 
    omicsData   = omicsData,
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = ExpDesign,
    factorInfo   = factorInfo)
  
  expect_equal(colnames(MAE[[1]]), as.vector(MAE[[1]]$samples))
  expect_equal(colnames(MAE[[2]]), as.vector(MAE[[2]]$samples))
  expect_equal(colnames(MAE[[3]]), as.vector(MAE[[3]]$samples))
})



test_that("Test check of NA in data", {
  
  omicsData <- list(
    RNAtest = ecoseed.df$RNAtest,
    metatest = ecoseed.df$metatest,
    protetest = ecoseed.df$protetest)
  
  omicsData[[1]][6,7] <- NA
  
  expect_no_error(RFLOMICS::createRflomicsMAE(
    projectName = "Tests", 
    omicsData   = omicsData,
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = ecoseed.df$design,
    factorInfo  = factorInfo))
  
  omicsData[[2]][6,7] <- -1
  
  expect_error(RFLOMICS::createRflomicsMAE(
    projectName = "Tests", 
    omicsData   = omicsData,
    omicsNames  = c("RNAtest", "metatest", "protetest"),
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = ecoseed.df$design,
    factorInfo  = factorInfo))
})


test_that("Test subRflomicsMAE", {
  
  miniMAE <- subRflomicsMAE(MAE, "protetest")
  expect_true("MultiAssayExperiment" %in% is(miniMAE))

  expect_equal(subRflomicsMAE(MAE), MAE)
  
  expect_null(subRflomicsMAE(MAE, "toto"))
})


test_that("readExpDesign", {
  
  design <- 
    data.frame(
      samples = paste0("S", 1:9),
      factor1 = rep(paste0("rep", 1:3),3),
      factor2 = c(rep("A",3), rep("B",3), rep("C",3))
    )
  
  fileName <- file.path(tempdir(), "design.txt")
  
  write.table(x = design, file = fileName, 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  design_bis <- RFLOMICS:::readExpDesign(fileName)
  expect_true("data.frame" %in% is(design_bis))
  expect_equal(names(design_bis), c("factor1", "factor2"))
})

test_that("readOmicsData", {
  
  data <- as.data.frame(matrix(sample(10:300, 100 * 10, replace = TRUE), 
                               nrow = 100, ncol = 10))
  
  data[[1]] <- paste0("gene", 1:100)
  
  names(data) <- c("genes", paste0("S", 1:9))
  
  fileName <- file.path(tempdir(), "design.txt")
  
  write.table(x = data, file = fileName, 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  data_bis <- RFLOMICS:::readOmicsData(fileName)
  expect_true("data.frame" %in% is(data_bis))
  expect_equal(names(data_bis), paste0("S", 1:9))
})


test_that("Test plot", {
  
  p <- plotConditionsOverview(MAE)
  expect_equal(is(p), "gg")
  
  p <- plotDataOverview(MAE)
  expect_equal(is(p), "gg")
  
  p <- plotDataOverview(MAE, realSize = TRUE)
  expect_equal(is(p), "gg")
})