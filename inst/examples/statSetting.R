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

# generate all statistical model formulae
formulae <- generateModelFormulae(MAE)

# chose and set model formula to rflomicsMAE object
MAE <- setModelFormula(MAE, formulae[[1]])

# Generate expression of contrasts from chosen model
contrastList <- generateExpressionContrast(MAE, "averaged")

# Set the contrasts List and choose the first 3 contrasts of type averaged
MAE <- setSelectedContrasts(MAE, contrastList = contrastList[c(1,2,3),])

