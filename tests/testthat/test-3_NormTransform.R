### ============================================================================
### [03_normTransform]
### ----------------------------------------------------------------------------
# A. Hulot, N. Bessoltane

library(testthat)
library(RFLOMICS)

##### Checks on transformation and normalization of data through command line functions
##### Essentially testing the .checkTransNorm pipeline.

# ---- Construction of objects for the tests ----

# load ecoseed data
data(ecoseed.mae)
data(ecoseed.df)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorType"   = c("batch", "Bio", "Bio")
)

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = ecoseed.mae,
  omicsTypes  = c("RNAseq","proteomics","metabolomics"),
  factorInfo  = factorInfo)


# Comparison data:
protMat <- ecoseed.df$protetest
rnaMat <- ecoseed.df$RNAtest

rnaMat <- rnaMat[, match(colnames(assay(MAE[["RNAtest"]])), colnames(rnaMat))]
protMat <- protMat[, match(colnames(assay(MAE[["protetest"]])), colnames(protMat))]

# ---- Some functions ----

isTransformed <- function(object, SE.name) RFLOMICS:::.isTransformed(object[[SE.name]])
isNorm        <- function(object, SE.name) RFLOMICS:::.isNormalized(object[[SE.name]])

######################################-
########### FUNCTIONS TESTS ###########
######################################-

# ---- transformData, apply_transform ----

test_that("transformData and apply_transform yield expected results", {

  MAE2 <- MAE3 <- MAE4 <- MAE5 <- MAE5b <- MAE6 <- MAE
  rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))

  ####
  # --- no transformation, no modification asked, nothing is supposed to happen.


  MAE2 <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none")
  expect_equal(assay(MAE2[["protetest"]]), as.matrix(protMat))
  expect(!RFLOMICS:::.isTransformed(MAE2[["protetest"]]), failure_message = "It was transformed, it shouldn't be.")

  ####
  # --- right transformation, no assay modification.

  MAE5 <- MAE5b <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest" , transformMethod = "log2")
  expect_equal( assay(MAE5[["protetest"]]), as.matrix(protMat))
  expect(!RFLOMICS:::.isTransformed(MAE5[["protetest"]]), failure_message = "It was transformed, it shouldn't be.")

  ####
  # --- apply transformation:

  MAE5b[["protetest"]] <- RFLOMICS:::.applyTransformation(MAE5[["protetest"]])
  expect_equal( assay(MAE5b[["protetest"]]), as.matrix(log2(protMat + 10^-10)))
  expect(RFLOMICS:::.isTransformed(MAE5b[["protetest"]]), failure_message = "It wasn't transformed, it should be.")

})

# ---- RunNormalization, apply_norm ----

test_that("RunNormalization and apply_norm yield expected results", {

  MAE2 <- MAE3 <- MAE4 <- MAE5 <- MAE5b <- MAE6 <- MAE

  se.p <- getProcessedData(MAE[["RNAtest"]])
  rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(se.p))

  ####
  # --- Missing norm argument
  MAE2 <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none")
  MAE2 <-
    runFeatureFiltering(MAE2, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none") |>
    runNormalization(SE.name = "protetest" )
  expect_equal( assay(MAE2[["protetest"]]), as.matrix(protMat))
  expect( !isNorm(MAE2,"protetest"), failure_message = "It was normalized, it shouldn't be.")

  MAE2 <- runFeatureFiltering(MAE, SE.name = "RNAtest")
  MAE2 <- runNormalization(MAE2, SE.name = "RNAtest" )
  expect_equal( assay(MAE2[["RNAtest"]]), as.matrix(rnaSeqMat))
  expect(getNormSettings(MAE2[["RNAtest"]])$method == "TMM", failure_message = "TMM was not forced on RNAseq data")
  expect(!isNorm(MAE2,"RNAtest"), failure_message = "It was normalized, it shouldn't be.")

  ####
  # --- no norm, no modification asked, nothing is supposed to happen.
  MAE2 <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none")
  MAE2 <-
    runFeatureFiltering(MAE2, SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none") |>
    runNormalization(SE.name = "protetest",normMethod = "none")
  expect_equal( assay(MAE2[["protetest"]]), as.matrix(protMat))
  expect(!isNorm(MAE2,"protetest"), failure_message = "It was normalized, it shouldn't be.")

  ####
  # --- right transformation.

  MAE5 <- MAE5b <- MAE |>
    runFeatureFiltering(SE.name = "protetest") |>
    runTransformData(SE.name = "protetest", transformMethod = "none") |>
    runNormalization(SE.name = "protetest", normMethod = "median")
  # warning message: proteomics data should be transformed before normalization.
  expect_equal( assay(MAE5[["protetest"]]), as.matrix(protMat))
  expect(!isNorm(MAE5, "protetest"), failure_message = "It was normalized, it shouldn't be.")

  ####
  # --- apply Normalization:
  # warning message: proteomics data should be transformed before normalization.
  MAE5b[["protetest"]] <- RFLOMICS:::.applyTransformation(MAE5[["protetest"]])
  MAE5b[["protetest"]] <- RFLOMICS:::.applyNormalization(MAE5b[["protetest"]])
  # getNorm(MAE5b, "protetest") # median
  protMed <- apply(protMat, 2, FUN = function(vect) vect - median(vect))
  expect_equal( assay(MAE5b[["protetest"]]), as.matrix(protMed))
  expect(isNorm(MAE5b, "protetest"), failure_message = "It wasn't normalized.")

})



######################################-
########### PROTEO/METABO ############
######################################-

# ---- PROTEO - Check if yielding the expected result ----
# Test all possibles combination (so far) of data transformation and data normalization for proteomics data
# also true for metabolomics data.
# Also test results from PCA (raw and norm data)

test_that("Transformation and normalisation combination - proteomics", {

  # log10 is not an allowed value for the parameter transformMethod Accepted values: log2, none
  # casesMat <- expand.grid(c("none", "log2", "log10", "log1p", "squareroot"), c("none", "median", "totalSum"))
  casesMat <- expand.grid(c("none", "log2"), c("none", "median", "totalSum"))

  colnames(casesMat) <- c("Trans", "Norm")

  res_equal <- lapply(1:nrow(casesMat), FUN = function(i){

    case_vect <- casesMat[i,]
    # print(case_vect)
    # matrix version
    protMattransnorm <- protMat

    pca.raw <- FactoMineR::PCA(t(protMattransnorm), ncp = 5, graph = FALSE)

    expect_equal(pca.raw$eig, MAE[["protetest"]]@metadata$PCAlist$raw$eig)
    expect_equal(pca.raw$svd, MAE[["protetest"]]@metadata$PCAlist$raw$svd)
    expect_equal(pca.raw$ind, MAE[["protetest"]]@metadata$PCAlist$raw$ind)
    expect_equal(pca.raw$var, MAE[["protetest"]]@metadata$PCAlist$raw$var)
    # the call is obligatory different between the two

    protMattransnorm <- switch(as.character(case_vect[[1]]),
                               "none"       = protMattransnorm,
                               "log2"       = log2(protMattransnorm + 10^-10),
                               "log10"      = log10(protMattransnorm + 10^-10),
                               "log1p"      = log1p(protMattransnorm),
                               "squareroot" = sqrt(protMattransnorm)
    )

    protMattransnorm <- switch(as.character(case_vect[[2]]),
                               "none"     = protMattransnorm,
                               "median"   = apply(protMattransnorm, 2, FUN = function(vect) vect - median(vect)),
                               "totalSum" = apply(protMattransnorm, 2, FUN = function(vect) vect/sum(vect^2))
    )

    pca.norm <- FactoMineR::PCA(t(protMattransnorm), ncp = 5, graph = FALSE)

    # RFLOMICS version
    MAE2 <- MAE
    MAE2 <-
      runFeatureFiltering(MAE2, SE.name = "protetest") |>
      runTransformData(SE = "protetest", transformMethod = as.character(case_vect[[1]]))
    MAE2 <-
      runNormalization(MAE2, SE.name = "protetest", normMethod = as.character(case_vect[[2]]))

    MAE2 <- RFLOMICS::runOmicsPCA(MAE2, SE = "protetest")
    expect_equal(pca.norm$eig, MAE2[["protetest"]]@metadata$PCAlist$norm$eig)
    expect_equal(pca.norm$svd, MAE2[["protetest"]]@metadata$PCAlist$norm$svd)
    expect_equal(pca.norm$ind, MAE2[["protetest"]]@metadata$PCAlist$norm$ind)
    expect_equal(pca.norm$var, MAE2[["protetest"]]@metadata$PCAlist$norm$var)

  })

})

# ---- PROTEO - Check if transformation working without arguments ----
test_that("Transformation - no method - no modification", {

  MAE2 <-
    runFeatureFiltering(MAE, SE.name = "protetest") |>
    runTransformData(SE = "protetest")

  expect_error(runTransformData(MAE, SE = "protetest"))
  expect(getTransSettings(MAE2[["protetest"]])$method == "log2", failure_message = "The transformation is not log2 by default.")
  expect_equal( assay(MAE[["protetest"]]),  assay(MAE2[["protetest"]]))
})

######################################-
########### RNASEQ ####################
######################################-

# ---- RNAseq - Check if yielding the expected result ----

test_that("RNAseq - none + TMM + log2", {


  # RFLOMICS version
  MAE2 <- MAE
  # It is not recommended to transform RNAseq data.
  # MAE2 <- runTransformData(MAE2, SE = "RNAtest", transformMethod = "none")

  # matrix version (filtering is done in the FlomicsMultiAssay constructor)
  rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE2[["RNAtest"]]))

  pca.raw <- FactoMineR::PCA(t(log2(rnaSeqMat + 1)), ncp = 5, graph = FALSE)

  expect_equal(pca.raw$eig, MAE[["RNAtest"]]@metadata$PCAlist$raw$eig)
  expect_equal(pca.raw$svd, MAE[["RNAtest"]]@metadata$PCAlist$raw$svd)
  expect_equal(pca.raw$ind, MAE[["RNAtest"]]@metadata$PCAlist$raw$ind)
  expect_equal(pca.raw$var, MAE[["RNAtest"]]@metadata$PCAlist$raw$var)
  # the call is obligatory different between the two

  MAE2 <- MAE2 |>
    runFeatureFiltering(SE.name = "RNAtest") |>
    runNormalization(SE.name = "RNAtest", normMethod = "TMM") |>
    RFLOMICS::runOmicsPCA(SE.name = "RNAtest", raw = FALSE)
  se.p <- getProcessedData(MAE2[["RNAtest"]], filter = TRUE)

  # Manually transforming and normalizing rnaSeqMat
  rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(se.p))
  normFactors <- edgeR::calcNormFactors(rnaSeqMat, method = "TMM")
  libSize <-  colSums(rnaSeqMat)

  tnDat <- scale(rnaSeqMat, center = FALSE, scale = normFactors*libSize/mean(normFactors*libSize))

  tabMAE <- RFLOMICS:::getCoeffNorm(MAE2[["RNAtest"]])
  normFactorMAE <- tabMAE$norm.factors
  names(normFactorMAE) <- rownames(tabMAE)
  libsizeMAE <- tabMAE$lib.size
  names(libsizeMAE) <- rownames(tabMAE)

  expect_equal(normFactorMAE, normFactors)
  expect_equal(libsizeMAE, libSize)

  pca.norm <- FactoMineR::PCA(t(log2(tnDat+1)), ncp = 5, graph = FALSE)

  expect_equal(pca.norm$eig, MAE2[["RNAtest"]]@metadata$PCAlist$norm$eig, tolerance = 0)
  expect_equal(pca.norm$svd, MAE2[["RNAtest"]]@metadata$PCAlist$norm$svd, tolerance = 0)
  expect_equal(pca.norm$ind, MAE2[["RNAtest"]]@metadata$PCAlist$norm$ind, tolerance = 0)
  expect_equal(pca.norm$var, MAE2[["RNAtest"]]@metadata$PCAlist$norm$var, tolerance = 0)

})

# ---- RNAseq - behaviour of transforData and RunNormalization ---

test_that("RNAseq - correct behaviour of normalization and transformation",  {

  MAE2 <- MAE

  MAE2 <- MAE2 |>
    runFeatureFiltering(SE.name = "RNAtest") |>
    runNormalization(SE = "RNAtest")

  expect(getNormSettings(MAE2[["RNAtest"]])$method == "TMM",
         failure_message = "Normalization is not defaulted to TMM for RNAseq data.")
  expect(!RFLOMICS:::.isNormalized(MAE2[["RNAtest"]]),
         failure_message = "RunNormalization transformed the data when not asked to")


})

