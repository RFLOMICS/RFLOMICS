### ============================================================================
### [07_MOFA2] test data integration with mofa
### ----------------------------------------------------------------------------
# A. Hulot

library(testthat)
library(RFLOMICS)
library(MOFA2)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for integration analysis : ----
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
    runDataProcessing(SE.name = "metatest", transformMethod = "log2",
                      normMethod = "median") |>
    runDiffAnalysis(SE.name = "metatest",  method = "limmalmFit",
                    p.adj.cutoff = 0.2)   |>
    runDataProcessing(SE.name = "protetest", transformMethod = "none",
                      normMethod = "median") |>
    runDataProcessing(SE.name = "RNAtest", normMethod = "TMM")    |>
    runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit",
                    p.adj.cutoff = 0.2)

MAE0 <- MAE

###############################################################################
## This is to test equivalence, launch these after modifying the main methods #
###############################################################################




# ----- TESTS -----

test_that("prepareForIntegration is running", {

    # Check parameters
    expect_no_error(prepareForIntegration(MAE,
                                          omicsNames = NULL,
                                          variableLists = NULL))
    expect_no_error(prepareForIntegration(MAE,
                                          omicsNames = c("protetest", "metatest"),
                                          variableLists = NULL))
    expect_no_error(prepareForIntegration(MAE,
                                          omicsNames = c("protetest", "metatest"),
                                          variableLists = rownames(MAE)))

    # Check output classes
    mofaObj <- prepareForIntegration(MAE,
                                     omicsNames = c("protetest", "metatest"),
                                     variableLists = rownames(MAE))
    expect(is(object = mofaObj, class2 = "MOFA"),
           failure_message = "prepareForIntegration without arguments does not return a MOFA object")

    mixObj <- prepareForIntegration(MAE,
                                    omicsNames = c("protetest", "metatest"),
                                    variableLists = rownames(MAE),
                                    method = "mixOmics")

    expect(is(object = mixObj, class2 = "list"),
           failure_message = "prepareForIntegration without mixOmics does not return a list")


    expect_error(prepareForIntegration(MAE,
                                       omicsNames = c("protein", "metabolites"),
                                       variableLists = rownames(MAE)))

})


# Prepared Objects for rest of the tests
mofaObj <- prepareForIntegration(MAE,
                                 omicsNames = c("protetest", "metatest"),
                                 variableLists = rownames(MAE))
mixObj <- prepareForIntegration(MAE,
                                omicsNames = c("protetest", "metatest"),
                                variableLists = rownames(MAE),
                                method = "mixOmics")

# check configuration for MOFA in following tests
configMOFA <- FALSE
catchRes <- RFLOMICS:::.tryCatch_rflomics(
    runOmicsIntegration(MAE, preparedObject = mofaObj,
                        method = "MOFA", scale_views = TRUE,
                        maxiter = 1000, num_factors = 5))

if (!is.null(catchRes$result)) {
    configMOFA <- TRUE
} else if (!is.null(catchRes$error)) {

    grepRes <-
        grep(pattern = "mofapy|scipy", catchRes$error)

    if (length(grepRes) != 0 &&  grepRes == 1) {
        message("To use MOFA2, you need to correctly set up a python environment. ",
                "We recommend reading the README (Troubleshooting MOFA2 section).
                If scipy error, please check your mofapy2 version (>0.7.1). Encountered error:\n",
                catchRes$error)

    } else {
        configMOFA <- TRUE
        # other error messages are interesting to keep to continue debug.
    }

}

if (!configMOFA)  skip("Tests for MOFA integration skipped.")
# end check config

if (configMOFA) {

    test_that("RunOmics is running", {

        expect_error(runOmicsIntegration(MAE))

        # Checking parameters
        expect_error(runOmicsIntegration(MAE, preparedObject = mofaObj,
                                         method = "unknown"))
        expect_error(runOmicsIntegration(MAE, preparedObject = mixObj))
        expect_error(runOmicsIntegration(MAE, preparedObject = mixObj,
                                         method = "mixOmics",
                                         selectedResponse = "unknown"))

        expect_no_error(runOmicsIntegration(MAE, preparedObject = mixObj,
                                            method = "mixOmics"))
        expect_no_error(runOmicsIntegration(MAE, preparedObject = mofaObj,
                                            method = "MOFA"))
    })

    # Complete analyses for getters/summary
    MAE <- runOmicsIntegration(MAE, preparedObject = mofaObj,
                               method = "MOFA")
    MAE <- runOmicsIntegration(MAE, preparedObject = mixObj,
                               method = "mixOmics")


    test_that("Getters and setters for integration are running", {

        # getMixOmics
        expect_no_error(getMixOmics(MAE))
        expect_no_error(getMixOmics(MAE, response = "temperature"))
        expect_no_error(getMixOmics(MAE, response = "imbibition"))
        expect(!is.null(getMixOmics(MAE, response = "temperature")),
               failure_message = "Temperature results are NULL")
        expect_no_error(getMixOmics(MAE, response = "temperature",
                                    onlyResults = FALSE))
        expect_no_error(getMixOmicsSettings(MAE))

        # getMOFA
        expect_no_error(getMOFA(MAE))
        expect(is(getMOFA(MAE), "MOFA"),
               failure_message = "getMOFA does not return a MOFA object")
        expect_no_error(getMOFA(MAE, onlyResults = FALSE))

        expect_no_error(getMOFASettings(MAE))


    })

    test_that("Summary MixOmics", {

        # Uses internal function .getOneMORes

        expect_no_error(sumMixOmics(MAE))
        expect(is(sumMixOmics(MAE), "list"),
               failure_message = "sumMixOmics does not return a list")
        expect_no_error(sumMixOmics(MAE, selectedResponse = "temperature"))
        expect(is(sumMixOmics(MAE, selectedResponse = "temperature"), "matrix"),
               failure_message = "sumMixOmics for temperature is not a matrix")

    })


    test_that("Plot functions are working", {
        expect_no_error(RFLOMICS:::plotMOVarExp(MAE, selectedResponse = "temperature"))
        expect(is(RFLOMICS:::plotMOVarExp(MAE, selectedResponse = "temperature"),
                  "ggplot"),
               failure_message = "plotMOVarExp does not return a ggplot")

        expect_no_error(RFLOMICS:::.relationsMOFA(getMOFA(MAE)))
        expect(is(RFLOMICS:::.relationsMOFA(getMOFA(MAE)), "data.frame"),
               failure_message = ".relationsMOFA is not a data.frame")
    })

} else {
    test_that("RunOmics is running", {
        # Checking parameters
        expect_error(runOmicsIntegration(MAE, preparedObject = mixObj,
                                         method = "mixOmics",
                                         selectedResponse = "unknown"))

        expect_no_error(runOmicsIntegration(MAE, preparedObject = mixObj,
                                            method = "mixOmics"))
    })

    # Complete analyses for getters/summary
    MAE <- runOmicsIntegration(MAE, preparedObject = mixObj,
                               method = "mixOmics")


    test_that("Getters and setters for integration are running", {

        # getMixOmics
        expect_no_error(getMixOmics(MAE))
        expect_no_error(getMixOmics(MAE, response = "temperature"))
        expect_no_error(getMixOmics(MAE, response = "imbibition"))
        expect(!is.null(getMixOmics(MAE, response = "temperature")),
               failure_message = "Temperature results are NULL")
        expect_no_error(getMixOmics(MAE, response = "temperature",
                                    onlyResults = FALSE))
        expect_no_error(getMixOmicsSettings(MAE))
    })

    test_that("Summary MixOmics", {

        # Uses internal function .getOneMORes

        expect_no_error(sumMixOmics(MAE))
        expect(is(sumMixOmics(MAE), "list"),
               failure_message = "sumMixOmics does not return a list")
        expect_no_error(sumMixOmics(MAE, selectedResponse = "temperature"))
        expect(is(sumMixOmics(MAE, selectedResponse = "temperature"), "matrix"),
               failure_message = "sumMixOmics for temperature is not a matrix")

    })


    test_that("Plot functions are working", {
        expect_no_error(RFLOMICS:::plotMOVarExp(MAE, selectedResponse = "temperature"))
        expect(is(RFLOMICS:::plotMOVarExp(MAE, selectedResponse = "temperature"),
                  "ggplot"),
               failure_message = "plotMOVarExp does not return a ggplot")
    })

}

#######################################################################
## Equivalence testing, run these after modifying the main methods    #
#######################################################################

test_that("Equivalence", {
    skip("Functional test for Integration skipped")
    # Comparison data
    protMat <- ecoseed.df$protetest
    metMat <- ecoseed.df$metatest
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

    orderNames <- rownames(colData(MAE))
    condMat <- condMat[match(orderNames, rownames(condMat)),]
    protMat <- protMat[match(orderNames, colnames(protMat))]
    metMat  <- metMat[match(orderNames, colnames(metMat))]

    identical(orderNames, colnames(protMat), attrib.as.set = FALSE)
    identical(orderNames, colnames(metMat), attrib.as.set = FALSE)

    metMat2  <- apply(log2(metMat + 10^(-10)), 2,
                      FUN = function(vect) vect - median(vect))
    protMat2 <- apply(protMat, 2, FUN = function(vect) vect - median(vect))

    # Select data
    selectedData <- c("metatest", "protetest")
    selectMethod <- c("metatest" = "diff", "protetest" = "none")
    operation <- c("metatest" = "union", "protetest" = "union")

    variableList <-
        lapply(selectedData, function(set) {
            switch(
                selectMethod[[set]],
                "diff"  = getDEList(
                    object = MAE[[set]],
                    contrasts = getSelectedContrasts(MAE)$contrastName,
                    operation = operation[[set]]
                ),
                "none"  = names(MAE[[set]])
            )
        })
    names(variableList) <- selectedData


    protMat2 <- protMat2[variableList[["protetest"]],]
    metMat2 <- metMat2[variableList[["metatest"]],]

    # Transform
    metMat3 <- t(scale(t(metMat2), center = TRUE, scale = TRUE))
    protMat3 <- t(scale(t(protMat2), center = TRUE, scale = TRUE))
    # protMat3 <- protMat2
    # metMat3 <- metMat2

    # from .rbeFunction
    colBatch <- getBatchFactors(MAE)
    newFormula <-
        gsub(pattern = paste(paste(colBatch, "[+]"),  collapse = "|"),
             "", getModelFormula(MAE))
    colData <- getDesignMat(MAE)
    designToPreserve <-
        model.matrix(as.formula(newFormula), data = colData)

    metMat4 <- limma::removeBatchEffect(metMat3, batch = colData[, colBatch],
                                        design = designToPreserve)
    protMat4 <- limma::removeBatchEffect(protMat3, batch = colData[, colBatch],
                                         design = designToPreserve)


    MOFA.obj <- suppressWarnings(
        prepareForIntegration(object           = MAE,
                              omicsNames       = selectedData,
                              variableLists    = variableList,
                              method           = "MOFA",
                              transformData    = TRUE
        ))

    # ---- Equivalence after preparation : ----
    expect(is(MOFA.obj, "MOFA"), failure_message = "Prepared MAE is not a MOFA object")
    expect_equal(get_dimensions(MOFA.obj)$D, lengths(variableList))

    protRes <- MOFA.obj@data$protetest$group1

    expect_equal(dim(protMat4), dim(protRes))
    expect(identical(rownames(protMat4), rownames(protRes),
                     attrib.as.set = FALSE),
           failure_message = "proteins rownames are not identical")
    expect(identical(colnames(protMat4), colnames(protRes),
                     attrib.as.set = FALSE),
           failure_message = "proteins colnames are not identical")
    expect_identical(as.data.frame(protMat4), as.data.frame(protRes))

    metaRes <- MOFA.obj@data$metatest$group1

    expect_equal(dim(metMat4), dim(metaRes))
    expect(identical(rownames(metMat4), rownames(metaRes),
                     attrib.as.set = FALSE),
           failure_message = "metabolites rownames are not identical")
    expect(identical(colnames(metMat4), colnames(metaRes),
                     attrib.as.set = FALSE),
           failure_message = "metabolites colnames are not identical")
    expect_identical(as.data.frame(metMat4), as.data.frame(metaRes))

    # ---- Equivalence on results: ----

    catchRes <- RFLOMICS:::.tryCatch_rflomics(
        runOmicsIntegration(MAE, preparedObject = MOFA.obj,
                            method = "MOFA", scale_views = TRUE,
                            maxiter = 1000, num_factors = 5))

    if(!is.null(catchRes$result)){
        MAE3 <- catchRes$result

        # equivalence
        mofaobject <- create_mofa(data = list("protetest" = protMat4,
                                              "metatest" = metMat4),
                                  extract_metadata = TRUE)
        data_opts  <- get_default_data_options(mofaobject)
        model_opts <- get_default_model_options(mofaobject)
        train_opts <- get_default_training_options(mofaobject)

        data_opts$scale_views  <- TRUE
        train_opts$maxiter     <- 1000
        train_opts$verbose     <- FALSE
        model_opts$num_factors <- 5
        MOFAObject.untrained <- prepare_mofa(
            object           = mofaobject,
            data_options     = data_opts,
            model_options    = model_opts,
            training_options = train_opts
        )

        MOFAObject.trained <- suppressWarnings(
            run_mofa(MOFAObject.untrained,
                     use_basilisk = FALSE,
                     save_data = TRUE))

        resRFLOMICS <- get_factors(getMOFA(MAE3))$group1
        resEquivalence <- get_factors(MOFAObject.trained)$group1

        resRFLOMICSW <- get_weights(getMOFA(MAE3))$group1
        resEquivalenceW <- get_weights(MOFAObject.trained)$group1

        expect_equal(resRFLOMICSW, resEquivalenceW)
        expect_equal(resRFLOMICS, resEquivalence, tolerance = 10^-5)

    }
    else if(!is.null(catchRes$error)){

        grepRes <-
            grep(pattern = "mofapy", catchRes$error)

        if(length(grepRes) != 0 &&  grepRes == 1)
            warning("To use MOFA2, you need to correctly set up a python environment. ",
                    "We recommend reading the README (Troubleshooting MOFA2 section).")

        expect(length(grepRes) != 0 &&  grepRes == 1,
               failure_message = catchRes$error)

    }
})


