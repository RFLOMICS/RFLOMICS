library(testthat)
library(RFLOMICS)

## Why is this code commented?
# Need either an internet connection or the file to test... 
# or org.at.tair

## ---- Construction MAE RFLOMICS ready for annotation analysis : ----
## load ecoseed data
data(ecoseed)

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
names(MAE) <- c("RNAtest", "protetest", "metatest")

formulae <- generateModelFormulae( MAE)
MAE <- setModelFormula(MAE, formulae[[1]])

contrastList <- generateExpressionContrast(object = MAE) |>
    purrr::reduce(rbind) |>
    dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                  "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
                                  "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))
MAE <- MAE |>
    setSelectedContrasts(contrastList = contrastList) |>
    runNormalization(SE.name = "protetest", normMethod = "median")  |>
    runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")

# |>
#     runCoExpression(SE.name = "RNAtest",
#                     K = 2:12,
#                     replicates = 2)

# ---- Annotation test function - DiffExpEnrichment ----
# 

#### Need org.At.tair.db to run
test_that("it's running from diffExpAnal - GO - proteomics", {

  # Selecting only one contrast
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE,
                                   nameList = "(temperatureElevated - temperatureLow) in imbibitionDS" , 
                                   SE.name = "protetest", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db",
                                                    keyType = "TAIR",
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))
  })

  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], 
                                   contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                   database = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO protetest from DiffExp) - There is no result in the enrichment metadata part.")

  # All contrasts
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db",
                                                    keyType = "TAIR",
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))

  })

  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS", database = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO protetest from DiffExp) - There is no result in the enrichment metadata part.")



})

# Need the custom file for this
# test_that("it's running from diffExpAnal - Custom - protetest", {
# 
#   df_custom <- vroom::vroom(file = paste0(system.file(package = "RFLOMICS"),
#                                           "/ExamplesFiles/ecoseed/AT_GOterm_EnsemblPlants.txt"),
#                             show_col_types = FALSE)
# 
#   MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "custom",
#                                  list_args = list(pvalueCutoff = 0.05),
#                                  col_term = "GO term accession",
#                                  col_gene = "Gene stable ID",
#                                  col_name = "GO term name",
#                                  col_domain = "GO domain",
#                                  annot = df_custom)
# 
#   expect(!is.null(MAE[["protetest"]]@metadata$DiffExpEnrichAnal$custom$summary), failure_message = "Custom proteomics didn't work (summary not present)")
# 
#   # Selecting only one contrast, custom file
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, nameList = "(temperatureMedium - temperatureLow) in imbibitionDS" , SE.name = "protetest", database = "custom",
#                                    list_args = list(pvalueCutoff = 0.05),
#                                    col_term = "GO term accession",
#                                    col_gene = "Gene stable ID",
#                                    col_name = "GO term name",
#                                    col_domain = "GO domain",
#                                    annot = df_custom)
#   })
# 
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], contrast = "(temperatureMedium - temperatureLow) in imbibitionDS" , database = "custom", domain = "biological_process")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from DiffExp) - There is no result in the enrichment metadata part.")
# 
#   #### Need org.At.tair.db to run
#   # All contrasts, GO database
#   # expect_no_error({
#   #   MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "GO",
#   #                                  list_args = list(OrgDb = "org.At.tair.db",
#   #                                                   keyType = "TAIR",
#   #                                                   pvalueCutoff = 0.05),
#   #                                  domain = c("BP", "MF", "CC"))
#   #
#   # })
# 
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS", 
#                                    database = "GO", domain = "BP")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from DiffExp) - There is no result in the enrichment metadata part.")
# 
#   # All contrasts, KEGG database
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "KEGG",
#                                    list_args = list(organism = "ath",
#                                                     keyType = "kegg",
#                                                     pvalueCutoff = 0.5))
# 
#   })
# 
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS", 
#                                    database = "KEGG")[["no-domain"]]
#     nrow(obj@result) > 0
#   }, failure_message = "(KEGG proteomics from DiffExp) - There is no result in the enrichment metadata part.")
# 
# })


# ---- Annotation test function - CoExpression enrichment ----

# #### Need org.At.tair.db to run
# test_that("it's running from CoExpAnal - GO - proteomics", {
# 
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, 
#                                    SE.name = "protetest", 
#                                    from = "CoExpAnal", 
#                                    database = "GO",
#                                    list_args = list(OrgDb = "org.At.tair.db",
#                                                     keyType = "TAIR",
#                                                     pvalueCutoff = 0.05),
#                                    domain = c("BP", "MF", "CC"))
# 
#   })
# 
# 
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], 
#                                    contrast = "cluster.1", 
#                                    database = "GO", 
#                                    domain = "BP", 
#                                    from = "CoExpAnal")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from CoExp) - There is no result in the enrichment metadata part.")
# 
# 
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
#                                    nameList = c("cluster.1", "cluster.2") ,
#                                    from = "CoExpAnal", 
#                                    database = "GO",
#                                    list_args = list(OrgDb = "org.At.tair.db",
#                                                     keyType = "TAIR",
#                                                     pvalueCutoff = 0.05),
#                                    domain = c("BP", "MF", "CC"))
# 
#   })
# 
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["protetest"]], 
#                                    contrast = "cluster.1", 
#                                    database = "GO", 
#                                    domain = "BP", from = "CoExpAnal")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from CoExp) - There is no result in the enrichment metadata part.")
# 
# })
