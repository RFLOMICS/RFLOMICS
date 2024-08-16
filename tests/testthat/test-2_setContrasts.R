library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in the code or package versions affect the expected results. 

# ---------------- RUN RFLOMICS ---------------

# load ecoseed data
data(ecoseed)

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


formulae <- generateModelFormulae(object = MAE) 
MAE <- setModelFormula(MAE, modelFormula = formulae[[1]])

# ----- 
test_that("generateExpressionContrast", {
    
    Contrasts.List.m <- generateExpressionContrast(MAE)
    Contrasts.List.f <- RFLOMICS:::.getExpressionContrastF(ExpDesign = colData(MAE), 
                                                           factorBio = c("temperature", "imbibition"), 
                                                           modelFormula = formulae[[1]])
    
    expect_equal(Contrasts.List.m, Contrasts.List.f)
})


