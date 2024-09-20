### ============================================================================
### [04_DE] 
### ----------------------------------------------------------------------------
# A. Hulot, 

library(testthat)
library(RFLOMICS)

#### Checks results of differential expression from edgeR and limma.
#### Using edger, limma and RFLOMICS functions to ensure everything is fine.

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


# ---- Construction of data tables differential analysis : ----

protMat <- ecoseed.df$protetest
rnaMat <- ecoseed.df$RNAtest
metMat <- ecoseed.df$metatest
condMat <- ecoseed.df$design

condMat$Repeat      <- factor(condMat$Repeat, levels = c("rep1", "rep2", "rep3"))
condMat$imbibition  <- factor(condMat$imbibition, levels = c("DS", "EI", "LI"))
condMat$temperature <- factor(condMat$temperature, levels = c("Low", "Medium", "Elevated"))

condMat$Repeat      <- relevel(condMat$Repeat, ref = "rep1")
condMat$imbibition  <- relevel(condMat$imbibition, ref = "DS")
condMat$temperature <- relevel(condMat$temperature, ref = "Low")

orderNames <- rownames(colData(MAE))
condMat <- condMat[match(orderNames, rownames(condMat)),]
protMat <- protMat[match(orderNames, colnames(protMat))]
rnaMat  <- rnaMat[match(orderNames, colnames(rnaMat))]
metMat  <- metMat[match(orderNames, colnames(metMat))]

identical(orderNames, colnames(protMat), attrib.as.set = FALSE)
identical(orderNames, colnames(rnaMat), attrib.as.set = FALSE)
identical(orderNames, colnames(metMat), attrib.as.set = FALSE)

# Contrasts
design <- model.matrix(~Repeat + temperature + imbibition + temperature:imbibition, data = condMat)

# Not checking if the coefficients are ok in there.
# taking the ones computed by RFLOMICS functions.

contrastsCoeff <- RFLOMICS:::.getContrastMatrixF(ExpDesign = condMat, factorBio = c("temperature", "imbibition"), modelFormula = "~Repeat + temperature + imbibition + temperature:imbibition",
                                                 contrastList = c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)", 
                                                                  "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3",
                                                                  "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))"))

MAE <- setSelectedContrasts(MAE, contrastList)

######################################-
########### FUNCTIONS TESTS ###########
######################################-

# ---- RNAseq : 

test_that("Differential analysis on RNAseq (counts) returns the same result within and outside of RFLOMICS pipeline", {
    
    ########################-
    ### RFLOMICS
    MAE[["RNAtest"]]@metadata$DiffExpAnal <- NULL
    
    MAE <- MAE |>
        runDataProcessing(SE.name = "RNAtest", normMethod = "TMM") |>
        runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit")
    
    ########################-
    ### equivalent pipeline
    se.pro <- getProcessedData(MAE[["RNAtest"]], filter = TRUE)
    rnaMat2 <- rnaMat %>%  dplyr::filter(rownames(.) %in% rownames(se.pro))
    expect_equal(dim(rnaMat2),  dim(se.pro))
    
    dge <- edgeR::DGEList(
      counts       = rnaMat2,
      group        = paste(condMat$temperature, condMat$imbibition, sep = "_"))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    
    
    rf.coef.norm <- as.data.frame(
      MAE[["RNAtest"]]@metadata$DataProcessing$Normalization$results$coefNorm)
    coef.norm <- as.data.frame(dge$samples)
  
    expect_equal(rf.coef.norm$norm.factors, coef.norm$norm.factors)
    expect_equal(rf.coef.norm$lib.size, coef.norm$lib.size)
    
    dge     <- edgeR::estimateGLMCommonDisp(dge, design = design)
    dge     <- edgeR::estimateGLMTrendedDisp(dge, design = design)
    dge     <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
    fit.RNA <- edgeR::glmFit(dge, design = design)
    
    ResGlm <-  lapply(1:nrow(contrastsCoeff), function(x){
        edgeR::glmLRT(fit.RNA, contrast = unlist(contrastsCoeff[x,]))
    })
    names(ResGlm) <- rownames(contrastsCoeff)
    
    ListResRNA <- lapply(ResGlm, function(x){
        res <- data.frame(edgeR::topTags(x, n = dim(x)[1])) %>%
            dplyr::rename("Abundance" = "logCPM", "pvalue" = "PValue", "Adj.pvalue" = "FDR")
        return(res)
    })
    names(ListResRNA) <- names(ResGlm)
    
    ## RNAseq
    DiffExpAnal <- getAnalysis(MAE[["RNAtest"]], 
                               name = "DiffExpAnal", 
                               subName = "results")
    
    MAESimple <- DiffExpAnal$DEF$`(temperatureElevated - temperatureLow) in imbibitionDS`
    MAEMean   <- DiffExpAnal$DEF$`(imbibitionEI - imbibitionDS) in mean`
    MAEInt    <- DiffExpAnal$DEF$`(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS`
    
    resSimple <- ListResRNA$`(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)`
    resMean   <- ListResRNA$`((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3`
    resInt    <- ListResRNA$`((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))`
    
    expect_equal(MAESimple, resSimple[match(rownames(MAESimple), rownames(resSimple)),])
    expect_equal(MAEMean, resMean)
    expect_equal(MAEInt, resInt)
    
})


# ---- Metabolomics : 

test_that("Diff Analysis on metabolomics returns the same result within and outside of RFLOMICS pipeline", {
    
    ########################-
    ### RFLOMICS
    
    MAE <- MAE |>
        runDataProcessing(SE.name = "metatest", transformMethod = "log2", 
                          normMethod = "totalSum")   |>
        runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")                       
    
    ########################-
    ### equivalent pipeline
    
    metMat2  <- apply(log2(metMat + 10^-10), 2, FUN = function(vect) vect/sum(vect^2))
    
    fitmet <- limma::lmFit(metMat2, design)
    ResGlmMet <-  lapply(1:nrow(contrastsCoeff), function(x){
        limma::contrasts.fit(fitmet, contrasts  = as.vector(unlist(contrastsCoeff[x,])))
    })
    names(ResGlmMet) <- rownames(contrastsCoeff)
    
    ListResMet <- lapply(ResGlmMet, function(x){
        fit2 <- limma::eBayes(x, robust = TRUE)
        res  <- limma::topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "AveExpr") %>% 
            dplyr::rename("Abundance" = "AveExpr","pvalue" = "P.Value", "Adj.pvalue" = "adj.P.Val")
        return(res)
    })
    
    # Tests results
    DiffExpAnal <- getAnalysis(MAE[["metatest"]], 
                               name = "DiffExpAnal", 
                               subName = "results")
    
    MAESimple <- DiffExpAnal$DEF$`(temperatureElevated - temperatureLow) in imbibitionDS`
    MAEMean   <- DiffExpAnal$DEF$`(imbibitionEI - imbibitionDS) in mean`
    MAEInt    <- DiffExpAnal$DEF$`(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS`
    
    resSimple <- ListResMet$`(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)`
    resMean   <- ListResMet$`((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3`
    resInt    <- ListResMet$`((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))`
    
    expect_equal(MAESimple, resSimple)
    expect_equal(MAEMean, resMean)
    expect_equal(MAEInt, resInt)
    
})

# ---- Proteomics : 

test_that("Diff Analysis on proteomics returns the same result within and outside of RFLOMICS pipeline", {
    
    ########################-
    ### RFLOMICS
    
    
    MAE <- MAE |>
        runDataProcessing(SE.name = "protetest", transformMethod = "none", 
                         normMethod = "median")    |>
        runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")
    
    
    ########################-
    ### equivalent pipeline
    
    protMat2 <- apply(protMat, 2, FUN = function(vect) vect - median(vect))
    
    fitprot <- limma::lmFit(protMat2, design)
    ResGlmProt <-  lapply(1:nrow(contrastsCoeff), function(x){
        limma::contrasts.fit(fitprot, contrasts  = as.vector(unlist(contrastsCoeff[x,])))
    })
    names(ResGlmProt) <- rownames(contrastsCoeff)
    
    ListResProt <- lapply(ResGlmProt, function(x){
        fit2 <- limma::eBayes(x, robust = TRUE)
        res  <- limma::topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "AveExpr") %>% 
            dplyr::rename("Abundance" = "AveExpr","pvalue" = "P.Value", "Adj.pvalue" = "adj.P.Val")
        return(res)
    })
    
    names(ListResProt) <- names(ResGlmProt)
    
    # Tests results
    DiffExpAnal <- getAnalysis(MAE[["protetest"]], 
                               name = "DiffExpAnal", 
                               subName = "results")
    
    MAESimple <- DiffExpAnal$DEF$`(temperatureElevated - temperatureLow) in imbibitionDS`
    MAEMean   <- DiffExpAnal$DEF$`(imbibitionEI - imbibitionDS) in mean`
    MAEInt    <- DiffExpAnal$DEF$`(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS`
    
    resSimple <- ListResProt$`(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)`
    resMean   <- ListResProt$`((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3`
    resInt    <- ListResProt$`((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))`
    
    expect_equal(as.matrix(MAESimple), as.matrix(resSimple[match(rownames(MAESimple), rownames(resSimple)),]))
    expect_equal(head(MAEMean), head(resMean))
    expect_equal(MAEInt, resInt)
})


