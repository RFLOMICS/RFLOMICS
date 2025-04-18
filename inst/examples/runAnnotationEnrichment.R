# load ecoseed data
library(RFLOMICS)
data(ecoseed.mae)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorType"   = c("batch", "Bio", "Bio")
)

# create rflomicsMAE object with ecoseed data
MAE <- createRflomicsMAE(
  projectName = "Tests",
  omicsData   = ecoseed.mae,
  omicsTypes  = c("RNAseq","proteomics","metabolomics"),
  factorInfo  = factorInfo)
#names(MAE) <- c("RNAtest", "protetest", "metatest")

# Set the statistical model and contrasts to test
formulae <- generateModelFormulae(MAE)
MAE <- setModelFormula(MAE, formulae[[1]])  

# Get the contrasts List and choose the first 3 contrasts of type averaged
contrastList <- generateExpressionContrast(MAE, "averaged")

MAE <- setSelectedContrasts(MAE, contrastList = contrastList[c(1, 2, 3),])

# Run the data preprocessing and perform the differential analysis 
MAE <- runDataProcessing(MAE, SE.name = "protetest",  
                         transformMethod = "log2",
                         normMethod = "median")

MAE <- runDiffAnalysis(MAE, 
                       SE.name = "protetest", 
                       method = "limmalmFit")

# Run GO annotation (enrichGO)
# Not run: need org.At.tair.db package
# MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
#                                list_args = list(OrgDb = "org.At.tair.db",
#                                                 keyType = "TAIR",
#                                                 pvalueCutoff = 0.05),
#                                from = "DiffExp", database = "GO",
#                                domain = "CC")

# Run KEGG annotation (enrichKEGG)
# need internet connection
# MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
#                                list_args = list(organism = "ath",
#                                                 keyType = "kegg",
#                                                 pvalueCutoff = 0.05),
#                                from = "DiffExp", database = "KEGG")


# Search for the pvalue cutoff:
# sumORA(MAE[["protetest"]], from = "DiffExp", database = "KEGG")

# need internet connection
# pathview package
# plotKEGG(MAE[["protetest"]], pathway_id = "ath00710", species = "ath",
#          contrastName = "cluster.4", from = "Coexp")

# From differential analysis proteins lists:
# plotClusterProfiler(MAE[["protetest"]],
#                     contrastName = "(temperatureElevated - temperatureMedium) in mean",
#                     database = "KEGG", from = "DiffExp",
#                     plotType = "heatplot", p.adj.cutoff = 0.05,
#                     domain = "no-domain")
# 
# plotEnrichComp(MAE[["protetest"]], from = "DiffExp",
#                database = "KEGG", matrixType = "FC")

# Get all results from KEGG on differential expression lists:
# results <- getEnrichRes(MAE[["protetest"]],
#                         from = "diffexp", database = "KEGG")

# Search for the pvalue cutoff:
# usedPvalue <- 
#     getEnrichPvalue(MAE[["protetest"]], from = "diffexp", database = "KEGG")
# settings <- 
#     getEnrichSettings(MAE[["protetest"]], from = "diffexp", database = "KEGG")

