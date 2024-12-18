### ============================================================================
### [05_coseq]
### ----------------------------------------------------------------------------
# A. Hulot,

library(testthat)
library(RFLOMICS)
library(coseq)

# ---- Construction of objects for the tests ----

# ---- Construction MAE RFLOMICS ready for differential analysis : ----
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

formulae <- generateModelFormulae( MAE)
MAE <- setModelFormula(MAE, formulae[[1]])

contrastList <- generateExpressionContrast(object = MAE) |>
  purrr::reduce(rbind) |>
  dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
                                "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))

MAE <- MAE |>
  setSelectedContrasts(contrastList)       |>
  # runDataProcessing(SE.name = "metatest", transformMethod = "log2",
  #                   normMethod = "totalSum") |>
  # runDiffAnalysis(SE.name = "metatest",  method = "limmalmFit")   |>
  # runDataProcessing(SE.name = "RNAtest", normMethod = "TMM")    |>
  # runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit")   |>
  runDataProcessing(SE.name = "protetest", transformMethod = "none",
                    normMethod = "median") |>
  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")

## ---- Construction of data tables differential analysis : ----

protMat <- ecoseed.df$protetest
# rnaMat <- ecoseed.df$RNAtest
# metMat <- ecoseed.df$metatest
condMat <- ecoseed.df$design

condMat$Repeat      <- factor(condMat$Repeat,
                              levels = c("rep1", "rep2", "rep3"))
condMat$imbibition  <- factor(condMat$imbibition,
                              levels = c("DS", "EI", "LI"))
condMat$temperature <- factor(condMat$temperature,
                              levels = c("Low", "Medium", "Elevated"))

condMat$Repeat      <- relevel(condMat$Repeat, ref = "rep1")
condMat$imbibition  <- relevel(condMat$imbibition, ref = "DS")
condMat$temperature <- relevel(condMat$temperature, ref = "Low")

orderNames <- rownames(condMat)

protMat <- protMat[match(orderNames, colnames(protMat))]
# rnaMat  <- rnaMat[match(orderNames, colnames(rnaMat))]
# metMat  <- metMat[match(orderNames, colnames(metMat))]

# Contrasts
design <- model.matrix(~Repeat + temperature + imbibition + temperature:imbibition, data = condMat)

# Not checking if the coefficients are ok in there.
# taking the ones computed by RFLOMICS functions.
# contrastsCoeff <- MAE@metadata$design@Contrasts.Coeff

## ---- Data Transformation ----
# rnaMat2 <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))
protMat2 <- apply(protMat, 2, FUN = function(vect) vect - median(vect))
# metMat2  <- apply(log2(metMat + 10^-10), 2, FUN = function(vect) vect/sum(vect^2))


# ---- Tests unitaires -----

test_that("Everything works as expected", {

  # MAE2 don't have any differential analysis results.
  expect_error(
    runCoExpression(
      object = MAE2, SE.name = "RNAtest", K = 2:10, replicates = 5, merge = "union",
      model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin",
      normFactors = "TMM"))

  # Wrong SE.name
  expect_error(
    runCoExpression(
      object = MAE, SE.name = "NA", K = 2:10, replicates = 5, merge = "union",
      model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin",
      normFactors = "TMM"))
})

## ---- Test: coseq.results.process ----
#
# test_that("When Median.min.rep doesn't correspond to a rep.ICL.min, an error message is returned",{
#
#   merge = "intersection"
#   K = 2:10
#   replicates = 2
#   iter <-  rep(K, each = replicates)
#   geneList <- getDEList(MAE[["RNAtest"]], operation = merge)
#
#   param.list = list(model = "Normal",
#                     GaussianModel = "Gaussian_pk_Lk_Ck",
#                     transformation = "arcsin",
#                     normFactors = "TMM",
#                     meanFilterCutoff = 50)
#
#   countMat <- SummarizedExperiment::assay(MAE[["RNAtest"]])[geneList,]
#
#   co <- capture.output(suppressMessages(
#     coseq.res.list <- lapply(1:replicates, function(x){
#
#       RFLOMICS:::.tryRflomics(
#         coseq::coseq(countMat, K = K, parallel = TRUE,
#                      model            = param.list[["model"]],
#                      transformation   = param.list[["transformation"]],
#                      meanFilterCutoff = param.list[["meanFilterCutoff"]],
#                      normFactors      = param.list[["normFactors"]],
#                      GaussianModel    = param.list[["GaussianModel"]],
#                      seed=4))
#     })))
#   names(coseq.res.list) <- c(1:replicates)
#
#   coseq.error.management <- RFLOMICS:::.coseq.error.manage(coseq.res.list = coseq.res.list,
#                                                            K = K,
#                                                            replicates = replicates)
#
#   CoExpAnal <-  RFLOMICS:::.coseq.results.process(coseqObjectList = coseq.error.management$coseq.res.list.values,
#                                                   K = K,
#                                                   conds = conds)
#
#   min.rep1 <- which.min(coseq::ICL(coseq.error.management$coseq.res.list.values[[1]]))
#   min.rep2 <- which.min(coseq::ICL(coseq.error.management$coseq.res.list.values[[2]]))
#
#   median <- unlist(lapply(1:length(K),
#                           function(i){
#                             median(coseq::ICL(coseq.error.management$coseq.res.list.values[[1]])[i],
#                                    coseq::ICL(coseq.error.management$coseq.res.list.values[[2]])[i])}))
#   names(median)=paste0("K=",K)
#   min.median <- which.min(median)
#
#   if((min.rep1 != min.median) & (min.rep2 != min.median)){
#     expect_equal(CoExpAnal[["error"]],"No min.median correspond to min.ICL.rep")
#   } else
#     expect_true(CoExpAnal[["results"]])
# })

## .....Test: coseq.error.managment

# test_that("For a given K, a likelihood equal to 0 is counted as failed job", {
#
#   merge = "intersection"
#   K = 2:10
#   replicates = 5
#
#   iter <-  rep(K, each = replicates)
#   geneList <- getDEList(MAE[["RNAtest"]], operation = merge)
#   countMat <- assay(MAE[["RNAtest"]])[geneList,][1:100,]
#
#   param.list = list(model = "Normal",
#                     GaussianModel = "Gaussian_pk_Lk_Ck",
#                     transformation = "arcsin",
#                     normFactors = "TMM",
#                     meanFilterCutoff = 50)
#
#   co <- capture.output(suppressMessages(
#     coseq.res.list <- lapply(1:replicates, function(x){
#       RFLOMICS:::.tryRflomics(
#         coseq::coseq(countMat, K = K, parallel = TRUE,
#                      model            = param.list[["model"]],
#                      transformation   = param.list[["transformation"]],
#                      meanFilterCutoff = param.list[["meanFilterCutoff"]],
#                      normFactors      = param.list[["normFactors"]],
#                      GaussianModel    = param.list[["GaussianModel"]],
#                      seed = x))
#     })))
#   names(coseq.res.list) <- c(1:replicates)
#
#   coseq.error.management <-
#     RFLOMICS:::.coseq.error.manage(coseq.res.list = coseq.res.list,
#                                    K = K,
#                                    replicates = replicates)
#
#   nbFailed <- dplyr::filter(coseq.error.management$jobs.tab.sum, K == "K=6")$n
#
#   K6 <-
#     unlist(
#       lapply(coseq.error.management$coseq.res.list.values,function(i){
#         coseq::likelihood(i)["K=6"]
#         })
#       )
#
#   nbK6eqNA <- sum(is.na(K6))
#   nbK6eq0 <- sum(K6 == 0, na.rm = TRUE)
#
#   expect_equal(nbFailed, nbK6eqNA + nbK6eq0)
#
# })

###########################################################
#######     FUNCTIONNAL TESTS                       #######
###########################################################

# ---- Tests seed -----


skip("Functionnal testing for Coseq skipped")
test_that("Two runs, same results - seed is working - RNAseq", {

  contrastNames <- c("(temperatureElevated - temperatureLow) in imbibitionDS",
                     "(imbibitionEI - imbibitionDS) in mean")
  MAE <- filterDiffAnalysis(MAE, SE.name = "RNAtest", logFC.cutoff = 2, p.adj.cutoff = 0.01)

  # getDEList(MAE[["RNAtest"]], contrasts = contrastNames, operation ="union")
  res1 <- runCoExpression(object = MAE, SE.name = "RNAtest",
                          K = 2:10, replicates = 5, merge = "union",
                          model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck",
                          transformation = "arcsin", normFactors = "TMM",
                          contrastNames = contrastNames
                          )

  res2 <- runCoExpression(object = MAE, SE.name = "RNAtest", K = 2:10,
                          replicates = 5, merge = "union",
                          model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck",
                          transformation = "arcsin", normFactors = "TMM",
                          contrastNames = contrastNames
                          )

  expect_identical(
    coseq::clusters(res1[["RNAtest"]]@metadata$CoExpAnal$results$coseqResults),
    coseq::clusters(res2[["RNAtest"]]@metadata$CoExpAnal$results$coseqResults))

})

test_that("Two runs, same results - seed is working - proteomics", {
    co <- capture.output(
        res1 <- runCoExpression(object = MAE, SE.name = "protetest", K = 2:20,
                                replicates = 5, merge = "union",
                                model = "Normal"))
    co <- capture.output(
        res2 <- runCoExpression(object = MAE, SE.name = "protetest", K = 2:20,
                                replicates = 5, merge = "union",
                                model = "Normal"))

    coexp_set <- getCoexpSettings(res2, SE.name = "protetest")
    expect_equal(coexp_set$K, 2:20)
    expect_equal(coexp_set$merge, "union")
    expect_equal(coexp_set$model, "Normal")

    expect_identical(
        coseq::clusters(
            getAnalysis(res1[["protetest"]],
                        name = "CoExpAnal",
                        subName = "results")$coseqResults),
        coseq::clusters(
            getAnalysis(res2[["protetest"]],
                        name = "CoExpAnal",
                        subName = "results")$coseqResults))

    # Check GaussianModel
    expect_equal(
        metadata(res1[["protetest"]])$CoExpAnal$results$coseqResults@metadata$GaussianModel,
        "Gaussian_pk_Lk_Bk")
    expect_equal(
        metadata(res2[["protetest"]])$CoExpAnal$results$coseqResults@metadata$GaussianModel,
        "Gaussian_pk_Lk_Bk")

    # Check transformation
    expect_equal(
        metadata(res1[["protetest"]])$CoExpAnal$results$coseqResults@transformation,
        "none")
    expect_equal(
        metadata(res2[["protetest"]])$CoExpAnal$results$coseqResults@transformation,
        "none")

    # Check normFactors
    expect_equal(
        metadata(res1[["protetest"]])$CoExpAnal$results$coseqResults@normFactors,
        rep(1, ncol(assay(res1[["protetest"]]))))

    expect_equal(
        metadata(res2[["protetest"]])$CoExpAnal$results$coseqResults@normFactors,
        rep(1, ncol(res2[["protetest"]])))
})

## ----- RNASeq -----
#
test_that("Coseq on RNAseq equivalence", {

  MAE <- filterDiffAnalysis(MAE, SE.name = "RNAtest",
                            logFC.cutoff = 0,
                            p.adj.cutoff = 0.05)

  # Parameters for the three analyses
  merge = "intersection"

  K = 2:10
  replicates = 3
  iter <-  rep(K, each = replicates)
  geneList <- getDEList(MAE[["RNAtest"]], operation = merge)

  param.list = list(model = "Normal",
                    GaussianModel = "Gaussian_pk_Lk_Ck",
                    transformation = "arcsin",
                    normFactors = "TMM",
                    meanFilterCutoff = 50)

  # RFLOMICS function:
  MAE <- runCoExpression(object = MAE, SE.name = "RNAtest",
                         K = K, replicates = replicates, merge = merge,
                         model = param.list$model,
                         GaussianModel = param.list$GaussianModel,
                         transformation = param.list$transformation,
                         normFactors = param.list$normFactors,
                         meanFilterCutoff = param.list$meanFilterCutoff
                         )

  clustersMAE <- coseq::clusters(MAE[["RNAtest"]]@metadata$CoExpAnal$results$coseqResults)
  resultsMAE  <- MAE[["RNAtest"]]@metadata$CoExpAnal$results$coseqResults
  tcountsMAE <- resultsMAE@tcounts
  yprofMAE <- resultsMAE@y_profiles

  # Test run coseqLocal only:
  countMat <- assay(MAE[["RNAtest"]])[geneList,]
  countMat <- countMat[, match(rownames(getDesignMat(MAE)), colnames(countMat))]
  identical(rownames(getDesignMat(MAE)), colnames(countMat),
            attrib.as.set = FALSE)

  test_rcl <- RFLOMICS:::.runCoseqLocal(countMat,
                                        conds = getDesignMat(MAE)$groups,
                                        K = K, replicates = replicates,
                                        param.list = param.list)


  clustersRCL <- coseq::clusters(test_rcl$coseqResults)
  tcountsRCL <- test_rcl$coseqResults@tcounts
  yprofRCL <- test_rcl$coseqResults@y_profiles

  # table(clustersRCL, clustersMAE)

  expect_equal(tcountsMAE, tcountsRCL)
  expect_equal(yprofRCL, yprofMAE)
  expect_identical(clustersRCL, clustersMAE)

  # Equivalent pipeline outside (supposed to be equivalent)

  countMat <- rnaMat[geneList,]
  countMat <- countMat[, match(rownames(getDesignMat(MAE)), colnames(countMat))]
  identical(rownames(getDesignMat(MAE)), colnames(countMat), attrib.as.set = FALSE)

  # set.seed(12345)
  co <- capture.output(suppressMessages(
    coseq_res <-  lapply(1:replicates, function(x){
      coseq::coseq(countMat, K = K,
                   model = param.list$model,
                   transformation = param.list$transformation,
                   GaussianModel = param.list$GaussianModel,
                   normFactors = param.list$normFactors,
                   meanFilterCutoff = param.list$meanFilterCutoff,
                   parallel = TRUE, seed = x)
    })
  ))
  names(coseq_res) <- c(1:replicates)
  s2 <- RFLOMICS:::.coseq.results.process(coseq_res, K,
                                          conds = getDesignMat(MAE)$group)

  max(unique(coseq::clusters(s2$coseqResults)))

  coseq.res <- s2$coseqResults

  clustersRES <- coseq::clusters(coseq.res)
  yprofRES <- coseq.res@y_profiles
  tcountsRES <- coseq.res@tcounts

  expect_equal(tcountsMAE, tcountsRES)
  expect_equal(yprofRES, yprofMAE)
  expect_identical(clustersRES, clustersMAE)

})

## ----- Proteomics -----

test_that("Coseq on Proteomics equivalence", {

    # Parameters for the three analyses
    merge = "union"
    geneList <- getDEList(MAE[["protetest"]], operation = merge)

    K = 2:20
    replicates = 3
    iter <-  rep(K, each = replicates)

    param.list = list(model = "Normal",
                      GaussianModel = "Gaussian_pk_Lk_Bk",
                      transformation = "none",
                      normFactors = "none",
                      meanFilterCutoff = NULL,
                      K = K,
                      replicates = replicates)

    # RFLOMICS function:
    co <- capture.output(
        MAE <- runCoExpression(object = MAE,
                               SE.name = "protetest",
                               contrastNames = NULL,
                               K = K, replicates = replicates,
                               merge = merge,
                               model = param.list$model,
                               GaussianModel = param.list$GaussianModel,
                               transformation = param.list$transformation,
                               normFactors = param.list$normFactors,
                               meanFilterCutoff = param.list$meanFilterCutoff
        ))

    clustersMAE <- coseq::clusters(metadata(MAE[["protetest"]])$CoExpAnal$results$coseqResults)
    resultsMAE  <- metadata(MAE[["protetest"]])$CoExpAnal$results$coseqResults
    tcountsMAE <- resultsMAE@tcounts
    yprofMAE <- resultsMAE@y_profiles

    # Test run coseqLocal only:

    countMat <- data.frame(assay(MAE[["protetest"]]))
    countMat <- apply(countMat, 2, FUN = function(x) x - median(x))
    countMat <- countMat[geneList,]
    countMat <- countMat[, match(rownames(getDesignMat(MAE)), colnames(countMat))]
    identical(rownames(getDesignMat(MAE)), colnames(countMat), attrib.as.set = FALSE)

    countMat2 <- t(apply(countMat, 1, function(x){scale(x, center = TRUE, scale = TRUE) }))
    colnames(countMat2) <- colnames(countMat)

    co <- capture.output(suppressWarnings(
        test_rcl <- RFLOMICS:::.runCoseqLocal(countMat2,
                                              param.list = param.list)))


    clustersRCL <- coseq::clusters(test_rcl$coseqResults)
    tcountsRCL <- test_rcl$coseqResults@tcounts
    yprofRCL <- test_rcl$coseqResults@y_profiles

    # table(clustersRCL, clustersMAE)

    expect_equal(tcountsMAE, tcountsRCL)
    expect_equal(data.frame(tcountsRCL), data.frame(countMat2)) # no transformation, no filtering
    expect_equal(yprofRCL, yprofMAE)
    expect_identical(clustersRCL, clustersMAE)

    # Equivalent pipeline outside (supposed to be equivalent)

    countMat <- protMat2[geneList,]
    countMat <- countMat[, match(rownames(getDesignMat(MAE)), colnames(countMat))]
    identical(rownames(getDesignMat(MAE)), colnames(countMat), attrib.as.set = FALSE)

    countMat2 <- t(apply(countMat , 1, function(x){scale(x, center = TRUE, scale = TRUE) }))
    colnames(countMat2) <- colnames(countMat)

    # set.seed(12345)
    co <- capture.output(suppressMessages(
        coseq_res <-  lapply(1:replicates, function(x){
            RFLOMICS:::.tryCatch_rflomics(coseq::coseq(countMat2, K = K,
                                                       model = param.list$model,
                                                       transformation = param.list$transformation,
                                                       GaussianModel = param.list$GaussianModel,
                                                       normFactors = param.list$normFactors,
                                                       meanFilterCutoff = param.list$meanFilterCutoff,
                                                       parallel = TRUE, seed = x))
        })
    ))
    names(coseq_res) <- paste0("K=", min(K), "-", max(K), "_", c(1:replicates))


    coseq.res <- RFLOMICS:::.coseq.results.process(coseq_res, K)$coseqResults

    clustersRES <- coseq::clusters(coseq.res)
    yprofRES <- coseq.res@y_profiles
    tcountsRES <- coseq.res@tcounts

    # table(clustersRES, clustersMAE)
    expect_equal(tcountsMAE, tcountsRES)
    expect_equal(data.frame(tcountsRES), data.frame(countMat2)) # no transformation, no filtering
    expect_equal(yprofRES, yprofMAE)
    expect_identical(clustersRES, clustersMAE)
})