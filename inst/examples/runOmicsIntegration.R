# load ecoseed data
library(RFLOMICS)
data(ecoseed.mae)

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
names(MAE) <- c("RNAtest", "protetest", "metatest")

formulae <- generateModelFormulae( MAE) 
MAE <- setModelFormula(MAE, formulae[[1]])
contrastList <- Reduce(rbind, generateExpressionContrast(MAE)) 

MAE <- MAE |>
  setSelectedContrasts(contrastList[c(3,6,25)]) |>
  runDataProcessing(SE.name = "metatest", 
                    transformMethod = "log2",
                    normMethod = "median") |>
  runDataProcessing(SE.name = "protetest", 
                    transformMethod = "none",
                    normMethod = "median")   |>
  runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")     |>
  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")    

# Integration using MOFA
# Prepare mofa object:
mofaObj <- prepareForIntegration(MAE,
                                 omicsNames = c("protetest", "metatest"),
                                 variableLists = rownames(MAE), 
                                 method = "MOFA")

# Perform integration: 
# Not run: MAEtest <- runOmicsIntegration(MAE, 
#                                preparedObject = mofaObj, 
#                                method = "MOFA", num_factors = 5)

# Integration using MixOmics
mixObj <- prepareForIntegration(MAE,
                                omicsNames = c("protetest", "metatest"),
                                variableLists = rownames(MAE),
                                method = "mixOmics")
MAEtest <- runOmicsIntegration(MAE, preparedObject = mixObj, 
                               method = "mixOmics")


# Access mixOmics results:
getMixOmics(MAEtest, response = "temperature")
getMixOmicsSettings(MAEtest)
# mixOmics::plotIndiv(getMixOmics(MAEtest, response = "imbibition"))

# Access MOFA2 results:
# getMOFA(MAEtest)
# getMOFASettings(MAEtest)
# MOFA2::plot_variance_explained(getMOFA(MAEtest))

