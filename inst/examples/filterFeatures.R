# load ecoseed data
library(RFLOMICS)
data("ecoseed.mae")

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
                    normMethod = "median") |>
  runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")     |>
  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")    |>
  runCoExpression(SE.name = "protetest", K = 2:10, replicates = 5, 
                  merge = "union", min.data.size = 10)

# Define the selections options
selOpt = list("protetest" = c("cluster.1", "H3"), "metatest" = c("DE"))
MAE1 <- filterFeatures(MAE, selOpt)
MAE1

selOpt2 = list("protetest" = c("cluster.2", "H25"), metatest = c("DE"))
MAE2 <- filterFeatures(MAE, selOpt2,
                       type = c("metatest" = "intersection",
                                "protetest" = 'union'))
MAE2

