library(RFLOMICS)

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
names(MAE) <- c("RNAtest", "protetest", "metatest")

formulae <- generateModelFormulae( MAE) 
MAE <- setModelFormula(MAE, formulae[[1]])
contrastList <- Reduce(rbind, generateExpressionContrast(MAE)) 

MAE <- MAE |>
    setSelectedContrasts(contrastList[c(3,6,25)]) |>
    runTransformData(SE.name = "metatest", transformMethod = "log2") |>
    runNormalization(SE.name = "metatest", normMethod = "median")    |>
    runNormalization(SE.name = "protetest", normMethod = "median")   |>
    runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")     |>
    runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")    

# Integration using MOFA
# Prepare mofa object:
mofaObj <- prepareForIntegration(MAE,
                                 omicsNames = c("protetest", "metatest"),
                                 variableLists = rownames(MAE),
                                 method = "MOFA")
class(mofaObj)

# Integration using MixOmics
mixObj <- prepareForIntegration(MAE,
                                omicsNames = c("protetest", "metatest"),
                                variableLists = rownames(MAE),
                                method = "mixOmics")
class(mixObj)