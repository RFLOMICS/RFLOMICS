### ============================================================================
### [06_clusterProfiler]
### ----------------------------------------------------------------------------
# A. Hulot,

library(testthat)
library(RFLOMICS)

## Why is this code commented?
# Need either an internet connection or the file to test...
# or org.at.tair

## ---- Construction MAE RFLOMICS ready for annotation analysis : ----
## load ecoseed data

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

formulae <- generateModelFormulae( MAE)
MAE <- setModelFormula(MAE, formulae[[1]])

contrastList <- generateExpressionContrast(object = MAE) |>
    purrr::reduce(rbind) |>
    dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)"))

MAE <- MAE |>
    setSelectedContrasts(contrastList = contrastList) |>
    runDataProcessing(SE.name = "protetest",
                      normMethod = "median",
                      transformMethod = "none")  |>
    runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")


########################################################
###                 TESTS WORKING                    ###
########################################################

# ---- Annotation test function - DiffExpEnrichment ----
#


#### Need org.At.tair.db to run

test_that("Automatically settings parameters or throwing errors
          - GO - Proteomics", {

              # automatically set settings if some are null

              skip_if_not_installed("org.At.tair.db")

              MAE[["protetest"]] <- runAnnotationEnrichment(
                  MAE[["protetest"]],
                  featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
                  database = "GO",
                  OrgDb = "org.At.tair.db",
                  pvalueCutoff = NULL,
                  qvalueCutoff = NULL,
                  minGSSize = NULL,
                  maxGSSize = NULL,
                  keyType = "TAIR",
                  domain = c("CC"),
                  from = "DiffExp")

              expect(
                  isTRUE(
                      getEnrichSettings(MAE[["protetest"]])$pvalueCutoff == 0.05 &&
                          getEnrichSettings(MAE[["protetest"]])$qvalueCutoff == 1 &&
                          getEnrichSettings(MAE[["protetest"]])$minGSSize == 10 &&
                          getEnrichSettings(MAE[["protetest"]])$maxGSSize == 500),
                  failure_message = "If null, one of qvalueCutoff (1),
          pvaluecutoff(0.05),
          minGsize (10) or maxGsize(500) are not set adequately")

              # error if NULL database
              expect_error(runAnnotationEnrichment(
                  MAE[["protetest"]],
                  featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
                  database = NULL,
                  OrgDb = "org.At.tair.db",
                  keyType = "TAIR",
                  domain = c("CC"),
                  from = "DiffExp"))

              # Error if domain is missing
              expect_error(runAnnotationEnrichment(
                  MAE[["protetest"]],
                  featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
                  database = "GO",
                  OrgDb = "org.At.tair.db",
                  keyType = "TAIR",
                  from = "DiffExp"))


              # Error if domain is not BP, CC, MF or ALL for GO
              expect_error(runAnnotationEnrichment(
                  MAE[["protetest"]],
                  featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
                  database = "GO",
                  domain = "domain", #
                  OrgDb = "org.At.tair.db",
                  keyType = "TAIR",
                  from = "DiffExp"))

              # Automatically setting domain for GO
              ## skip for long time running...
              # MAE[["protetest"]] <- runAnnotationEnrichment(
              #     MAE[["protetest"]],
              #     featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
              #     database = "GO",
              #     domain = "ALL",
              #     OrgDb = "org.At.tair.db",
              #     keyType = "TAIR",
              #     from = "DiffExp")

             #  expect(isTRUE(setequal(getEnrichSettings(MAE[["protetest"]])$domain,
             #                         c("BP", "CC", "MF"))),
             #         failure_message = "If domain is not set for GO, then domain
             # is not automatically set to BP, CC and MF")

              # Annotation not running if from is neither diffexp nor coexp
              expect_error(
                  runAnnotationEnrichment(
                      MAE[["protetest"]],
                      database = "GO",
                      domain = "ALL",
                      from = "Diffexpression",
                      OrgDb = "org.At.tair.db",
                      keyType = "TAIR"))

          })


test_that("Annotation not running if diffexp does not exist", {
    MAE_test <- MAE
    metadata(MAE_test[["protetest"]])$DiffExpAnal <- list()

    skip_if_not_installed("org.At.tair.db")
    expect_error(
        runAnnotationEnrichment(
            MAE_test[["protetest"]],
            database = "GO",
            domain = "ALL",
            OrgDb = "org.At.tair.db",
            keyType = "TAIR"))
})


test_that("Automatically filling parameters", {
    # automatic from to diffexp if null and no featureList
    skip_if_not_installed("org.At.tair.db")

    MAE[["protetest"]] <- runAnnotationEnrichment(
        MAE[["protetest"]],
        database = "GO",
        from = NULL,
        OrgDb = "org.At.tair.db",
        keyType = "TAIR",
        domain = c("CC"))

    expect(
        isTRUE(setequal(names(getEnrichRes(MAE[["protetest"]])),
                        getSelectedContrasts(MAE[["protetest"]])$contrastName)),
        failure_message = "featureNamesList is not automatically filled with contrasNames")

})


# ---- Getters tests ----


test_that("getEnrichRes - GO - proteomics is working", {
    skip_if_not_installed("org.At.tair.db")

    # Selecting only one contrast
    MAE <-
        runAnnotationEnrichment(
            MAE, SE.name = "protetest",
            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
            database = "GO",
            pvalueCutoff = 0.2,
            OrgDb = "org.At.tair.db",
            keyType = "TAIR",
            domain = c("CC"))

    # Identical to the one in "it's running from diffExpAnal - GO - proteomics"
    # # Getter from SE
    # expect({
    #     obj <- getEnrichRes(
    #         MAE[["protetest"]],
    #         featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
    #         database = "GO", domain = "BP")
    #     nrow(obj@result) > 0
    #
    # }, failure_message = "getEnrichRes from DiffExp SE - ",
    # "There is no result in the enrichment metadata part.")


    # Getter from MAE
    expect({
        obj <- getEnrichRes(
            MAE[["protetest"]],
            featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
            database = "GO", domain = "CC")
        nrow(obj@result) > 0

    }, failure_message = "getEnrichRes from DiffExp MAE -  ",
    "There is no result in the enrichment metadata part."
    )

    expect_no_error(getEnrichRes(MAE[["protetest"]]))

    expect_identical(getEnrichRes(MAE[["protetest"]]),
                     getEnrichRes(MAE[["protetest"]],
                                  from = "DiffExp", database = "GO"))

    expect_error(getEnrichRes(MAE))
    expect_error(getEnrichRes(MAE[["protetest"]], database = "data"))


    # SumORA testing
    expect_no_error(sumORA(MAE[["protetest"]], database = "GO",
                           featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS" ))

    expect(isTRUE(is.data.frame(sumORA(MAE[["protetest"]], database = "GO"))),
           failure_message = "sumORA without the featureListName is not a dataframe.")

    MAE <- runAnnotationEnrichment(
        MAE, SE.name = "protetest",
        database = "KEGG",
        organism = "ath",
        pvalueCutoff = 1)

    expect_no_error(getAnnotAnalysesSummary(MAE))
    expect_no_error(getAnnotAnalysesSummary(MAE, matrixType = "log2FC"))
    expect_no_error(getAnnotAnalysesSummary(MAE, matrixType = "GeneRatio"))
    expect_no_error(getAnnotAnalysesSummary(MAE, matrixType = "p.adjust"))
    expect_no_error(getAnnotAnalysesSummary(MAE, matrixType = "presence"))
    expect_error(getAnnotAnalysesSummary(MAE, matrixType = "type"))
    expect_error(getAnnotAnalysesSummary(MAE, from = "from"))

})



test_that("getEnrichSettings - GO - proteomics is working", {
    MAE <-
        runAnnotationEnrichment(
            MAE, SE.name = "protetest",
            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
            database = "GO",
            pvalueCutoff = 0.2,
            OrgDb = "org.At.tair.db",
            keyType = "TAIR",
            domain = c("CC"))

    expect_error(getEnrichSettings(MAE[["protetest"]], from = "Diffexp"))
})



# ---- KEGG tests (require internet connection) ----


test_that("KEGG annotation can be run", {

    # expects an organism name
    expect_error({
        runAnnotationEnrichment(
            MAE, SE.name = "protetest",
            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
            database = "KEGG",
            pvalueCutoff = 1)
    })

    # can be run if proper setting
    expect_no_error({
        runAnnotationEnrichment(
            MAE, SE.name = "protetest",
            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
            database = "KEGG",
            organism = "ath",
            pvalueCutoff = 1)
    })

    MAE <- runAnnotationEnrichment(
        MAE, SE.name = "protetest",
        featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
        database = "KEGG",
        organism = "ath",
        pvalueCutoff = 1)

    # automatically filling the keytype if none is provided:
    expect(setequal(getEnrichRes(MAE[["protetest"]], database = "KEGG")[[1]][["no-domain"]]@keytype,
                    "kegg"),
           failure_message = "KEGG enrichment: when providing no keytype,
           the keytype is not automatically kegg")

})

# ---- Custom tests ----

test_that("runEnrichment analysis with custom file is possible", {

    # don't run with empty argument
    expect_error(runAnnotationEnrichment(MAE[["protetest"]],
                 featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                 database = "custom"))

    # don't run with empty file
    expect_error(runAnnotationEnrichment(MAE[["protetest"]],
                            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                            database = "custom", annotation = data.frame()))

    # dummy annotation
    dummy_custom_annotation <- data.frame(
        "gene" = c("AT5G58070", "AT5G58070", "AT5G58070", "AT5G58070",
                   "AT4G04460", "AT3G08900", "AT3G08900",
                   "AT1G05510", "AT1G07080", "AT1G07140",
                   "AT1G13930", "AT1G14170", "AT1G14930",
                   "AT1G30360", "AT1G32060", "AT1G33120"),
        "term" = c("term1", "term2", "term3", "term4", "term2", "term2", "term3",
                   "term2", "term2", "term2", "term2", "term2", "term2", "term2",
                   "term2", "term2"),
        "name" = c("An imaginary term 1", "An imaginary term 2", "An imaginary term 3",
                   "An imaginary term 4", "An imaginary term 2", "An imaginary term 2",
                   "An imaginary term 3", "An imaginary term 2", "An imaginary term 2",
                   "An imaginary term 2", "An imaginary term 2", "An imaginary term 2",
                   "An imaginary term 2", "An imaginary term 2", "An imaginary term 2",
                   "An imaginary term 2"),
        "domain" = c("D1", "D2", "D1", "D1", "D2", "D2", "D1" , "D2", "D2", "D2",
                     "D2", "D2", "D2", "D2", "D2", "D2")
    )

    # can't run if gene or term is missing
    expect_error(runAnnotationEnrichment(MAE[["protetest"]],
                            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                            database = "custom", annotation =  dummy_custom_annotation[-1]))
    expect_error(runAnnotationEnrichment(MAE[["protetest"]],
                                         featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                         database = "custom",
                                         annotation =  dummy_custom_annotation[-2]))

    # run as expected if no name or domain or both
    expect_no_error(runAnnotationEnrichment(MAE[["protetest"]],
                                            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                            database = "custom",
                                            annotation =  dummy_custom_annotation[-3]))

    expect_no_error(runAnnotationEnrichment(MAE[["protetest"]],
                                            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                            database = "custom",
                                            annotation =  dummy_custom_annotation[-4]))

    expect_no_error(runAnnotationEnrichment(MAE[["protetest"]],
                                            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                            database = "custom",
                                            annotation =  dummy_custom_annotation[-c(3,4)]))

    expect_no_error(runAnnotationEnrichment(MAE[["protetest"]],
                            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS",
                            database = "custom",
                            annotation =  dummy_custom_annotation))

})

# ---- PlotClusterProfiler ----


test_that("plotClusterPRofiler - GO only", {
    skip_if_not_installed("org.At.tair.db")

    MAE <-
        runAnnotationEnrichment(
            MAE, SE.name = "protetest",
            featureList = "(temperatureElevated - temperatureLow) in imbibitionDS" ,
            database = "GO",
            pvalueCutoff = 1,
            OrgDb = "org.At.tair.db",
            keyType = "TAIR",
            domain = c("CC"))

    # trying to plot without indicating any feature list
    expect_error(plotClusterProfiler(MAE[["protetest"]], database = "GO"))

    # trying to plot without indicating any database
    expect_error(plotClusterProfiler(MAE[["protetest"]],
                                     featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS"))

    # trying to plot without indicating a domain name for GO
    expect_error(plotClusterProfiler(MAE[["protetest"]],
                        featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
                        database = "GO"))

    # trying to plot with a null a domain name for GO
    expect_error(plotClusterProfiler(MAE[["protetest"]],
                                     featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                     database = "GO", domain = NULL))

    # trying to plot
    expect_no_error(plotClusterProfiler(MAE[["protetest"]],
                                        featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                        database = "GO", domain = "CC"))

    plotToTest <- plotClusterProfiler(MAE[["protetest"]],
                                      featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                      database = "GO", domain = "CC")

    expect(is(plotToTest, "ggplot"),
           failure_message = "plotCusterProfiler did not return a ggplot object")

    expect(is(plotToTest, "enrichplotDot"),
           failure_message = "plotCusterProfiler did not return an enrichPlot")

    # trying to plot other types of plots
    expect_no_error(plotClusterProfiler(MAE[["protetest"]],
                                        featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                        database = "GO", domain = "CC",
                                        plotType = "heatplot"))

    expect_no_error(plotClusterProfiler(MAE[["protetest"]],
                                        featureListName = "(temperatureElevated - temperatureLow) in imbibitionDS",
                                        database = "GO", domain = "CC",
                                        plotType = "cnetplot"))
})


#---- PlotEnrichComp ----

test_that("plotEnrichComp - KEGG only", {

    contrastList <- generateExpressionContrast(object = MAE) |>
        purrr::reduce(rbind) |>
        dplyr::filter(tag %in% c("H1", "H2","H3"))

    MAE <- MAE |>
        setSelectedContrasts(contrastList = contrastList) |>
        runDataProcessing(SE.name = "protetest",
                          normMethod = "median",
                          transformMethod = "none")  |>
        runDiffAnalysis(SE.name = "protetest", method = "limmalmFit",
                        p.adj.cutoff = 0.05) |>
        runAnnotationEnrichment(
            SE.name = "protetest",
            featureList = NULL,
            database = "KEGG",
            pvalueCutoff = 1,
            organism = "ath",
            keyType = "kegg")

    expect_no_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG"))
    expect_no_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG", matrixType = "log2FC"))
    expect_no_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG", matrixType = "GeneRatio"))
    expect_no_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG", matrixType = "p.adjust"))
    expect_no_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG", matrixType = "presence"))
    expect_no_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG"))
    expect_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG", domain = "domain"))
    expect_error(plotEnrichComp(MAE[["protetest"]]))
    expect_error(plotEnrichComp(MAE[["protetest"]], database = "KEGG", matrixType = "type"))
    expect_error(plotEnrichComp(MAE[["protetest"]], database = "data", matrixType = "type"))

})


# --- CUSTOM FILE NEEDED ----
# test_that("it's running from diffExpAnal - Custom - protetest", {
#
#   df_custom <- vroom::vroom(file = paste0(system.file(package = "RFLOMICS"),
#                                           "/ExamplesFiles/ecoseed/AT_GOterm_EnsemblPlants.txt"),
#                             show_col_types = FALSE)
#
#   MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "custom",
#                                  pvalueCutoff = 0.05,
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
#                                    pvalueCutoff = 0.05,
#                                    col_term = "GO term accession",
#                                    col_gene = "Gene stable ID",
#                                    col_name = "GO term name",
#                                    col_domain = "GO domain",
#                                    annot = df_custom)
#   })
#
#   expect({
#     obj <- getEnrichRes(MAE[["protetest"]], contrast = "(temperatureMedium - temperatureLow) in imbibitionDS" , database = "custom", domain = "biological_process")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from DiffExp) - There is no result in the enrichment metadata part.")
#
#   #### Need org.At.tair.db to run
#   # All contrasts, GO database
#   # expect_no_error({
#   #   MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "GO",
#   #                                  OrgDb = "org.At.tair.db",
#   #                                                   keyType = "TAIR",
#   #                                                   pvalueCutoff = 0.05,
#   #                                  domain = c("BP", "MF", "CC"))
#   #
#   # })
#
#   expect({
#     obj <- getEnrichRes(MAE[["protetest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
#                                    database = "GO", domain = "BP")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from DiffExp) - There is no result in the enrichment metadata part.")
#
#   # All contrasts, KEGG database
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest", database = "KEGG",
#                                    organism = "ath",
#                                                     keyType = "kegg",
#                                                     pvalueCutoff = 0.5)
#
#   })
#
#   expect({
#     obj <- getEnrichRes(MAE[["protetest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
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
#                                    OrgDb = "org.At.tair.db",
#                                                     keyType = "TAIR",
#                                                     pvalueCutoff = 0.05,
#                                    domain = c("BP", "MF", "CC"))
#
#   })
#
#
#   expect({
#     obj <- getEnrichRes(MAE[["protetest"]],
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
#                                    OrgDb = "org.At.tair.db",
#                                                     keyType = "TAIR",
#                                                     pvalueCutoff = 0.05,
#                                    domain = c("BP", "MF", "CC"))
#
#   })
#
#   expect({
#     obj <- getEnrichRes(MAE[["protetest"]],
#                                    contrast = "cluster.1",
#                                    database = "GO",
#                                    domain = "BP", from = "CoExpAnal")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO proteomics from CoExp) - There is no result in the enrichment metadata part.")
#
# })
