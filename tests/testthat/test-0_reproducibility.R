### ============================================================================
### [00_reproducibility] Ensures the reproducibility of results in case of code 
### modifications or changes in package versions.
### ----------------------------------------------------------------------------
# N. Bessoltane, 

library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in 
# the code or package versions affect the expected results. 

# ---- Create RflomicsMAE object ----

# load ecoseed data
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

SE <- getRflomicsSE(MAE, "RNAtest")

sampleToKeep <- colnames(MAE[["protetest"]])[-7]

# ---- check design ----

test_that("test if RflomicsMAE / RflomicsSE", {
  
  check <- checkExpDesignCompleteness(MAE, omicName = "RNAtest")
  expect_false(check$error)
  expect_equal(check$message, "The experimental design is complete and balanced.")
})

## choice of model formulae
formulae <- generateModelFormulae(MAE)

MAE <- setModelFormula(MAE, modelFormula = formulae[[1]])

## list of all hypothesis grouped by contarst type (run on MAE or SE)
Contrasts.List <- generateExpressionContrast(MAE)

# utilisation en dehors de MAE ou SE


## choice of hypothesis
selectedContrasts <- rbind(Contrasts.List$simple[1:3,],
                           Contrasts.List$averaged[1:3,],
                           Contrasts.List$interaction[1:3,])
# selectedContrasts <- Contrasts.List$simple[1:4,]

MAE <- setSelectedContrasts(MAE, contrastList = selectedContrasts)

## data processing
## runDataProcessing uses the .raw names, not the SE name.
## creates the non .raw SE inside the MAE.
## Interface function, mostly
MAE <- MAE |> 
  runDataProcessing(SE.name = "RNAtest", samples = sampleToKeep,
                    filterStrategy = "NbReplicates", 
                    cpmCutoff = 1, normMethod = "TMM") |>
  runDataProcessing(SE.name = "protetest", samples = NULL,
                    normMethod = "none", transformMethod = "none") |>
  runDataProcessing(SE.name = "metatest", samples = NULL, 
                    normMethod = NULL, transformMethod = "log2")


test_that("Test generateReport", {
  expect_error(generateReport(MAE))

  co <- capture.output(
    suppressWarnings(
      expect_no_error(
        generateReport(
          MAE,
          archiveName = file.path(tempdir(), "archiv_report.tar.gz")
        )
      )
    )
  )
})

## diff analysis
MAE <- MAE |> 
  # runDiffAnalysis(SE.name = "RNAtest", p.adj.method = "BH", 
  #                 method = "edgeRglmfit", p.adj.cutoff = 0.05, 
  #                 logFC.cutoff = 0) |>
  runDiffAnalysis(SE.name = "protetest", p.adj.method = "BH", 
                  method = "limmalmFit",  p.adj.cutoff = 0.05, 
                  logFC.cutoff = 0) |>
  runDiffAnalysis(SE.name = "metatest", p.adj.method = "BH", 
                  method = "limmalmFit",  p.adj.cutoff = 0.05, 
                  logFC.cutoff = 0)

# MAE <- 
#   setValidContrasts(MAE, omicName = "RNAtest",
#                     contrastList = getSelectedContrasts(MAE[["RNAtest"]]))

MAE[["protetest"]] <- 
  setValidContrasts(MAE[["protetest"]], 
                    contrastList = getSelectedContrasts(MAE[["protetest"]])[1:6,])
MAE[["metatest"]] <- 
  setValidContrasts(MAE[["metatest"]], 
                    contrastList = getSelectedContrasts(MAE[["metatest"]])[1:6,])

## co expression
co <- capture.output(
  MAE <- MAE |> 
    runCoExpression(
      SE.name = "protetest",
      K = 2:10, 
      replicates = 5, 
      merge = "union", 
      model = "Normal", 
      min.data.size = 10)
) #|>


## ---- getRflomicsSE ----
test_that("test getRflomicsSE", {
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_true(is(SE, "RflomicsSE"))
  
  expect_null(getRflomicsSE(MAE))
  expect_null(getRflomicsSE(MAE, "toto"))
  
  expect_equal(getDatasetNames(SE), "RNAtest")
})

## ---- sub set of RflomicsMAE ----
test_that("test subRflomicsMAE", {
  
  miniMAE <- suppressWarnings(MAE[,, c("RNAtest", "metatest")])
  miniMAE.rf <- 
    suppressWarnings(subRflomicsMAE(MAE, c("RNAtest", "metatest")))
  expect_equal(miniMAE.rf, miniMAE)
})

# ---- contrasts ----

test_that("contrast", {
  
  Contrasts.names <- c("(temperatureMedium - temperatureLow) in imbibitionDS",                                                                                                                                                     
                       "(temperatureElevated - temperatureLow) in imbibitionDS",                                                                                                                                                   
                       "(temperatureElevated - temperatureMedium) in imbibitionDS",                                                                                                                                              
                       "(temperatureMedium - temperatureLow) in mean",          
                       "(temperatureElevated - temperatureLow) in mean",
                       "(temperatureElevated - temperatureMedium) in mean",
                       "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
                       "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
                       "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS")
  
  Contrasts.Coeff <- rbind(
    c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0, 0, 0.333333333333333, 0, 0.333333333333333, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0.333333333333333, 0, 0.333333333333333),
    c(0, 0, 0, -1, 1, 0, 0, -0.333333333333333, 0.333333333333333, -0.333333333333333, 0.333333333333333),
    c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
  ) %>% as.data.frame()
  
  names(Contrasts.Coeff) <- c("Intercept", "Repeatrep2", "Repeatrep3", 
                              "temperatureMedium", "temperatureElevated", 
                              "imbibitionEI", "imbibitionLI", 
                              "temperatureMedium:imbibitionEI", "temperatureElevated:imbibitionEI",
                              "temperatureMedium:imbibitionLI", "temperatureElevated:imbibitionLI")
  
  row.names(Contrasts.Coeff) <- Contrasts.names
  
  #expect_equal(getDiffSettings(MAE[["RNAtest"]])$contrastCoef, Contrasts.Coeff)
  expect_equal(getDiffSettings(MAE, "protetest")$contrastCoef, Contrasts.Coeff)
  expect_equal(getDiffSettings(MAE[["metatest"]])$contrastCoef, Contrasts.Coeff)
  
  expect_equal(
    checkExpDesignCompleteness(MAE[["RNAtest"]], 
                               sampleList = getSelectedSamples(MAE[["RNAtest"]]))$messages, 
    "The experimental design is complete but not balanced.")
  
  
})

# ---- processing ----

test_that("data processing", {
  
  ## low count Filtering
  # remplacer par le bon getter
  nb_lowGene <- length(getFilteredFeatures(MAE[["RNAtest"]]))
  expect_equal(nb_lowGene, 6725)
  expect_equal(getFilteredFeatures(MAE[["protetest"]]), NULL)
  expect_equal(getFilteredFeatures(MAE[["metatest"]]), NULL)
  
  ## sample filtering
  expect_equal(getSelectedSamples(MAE[["RNAtest"]]), sampleToKeep)
  
  ## transformation
  
  ## Normalisation
  
})


# ---- resetRflomicsMAE ----
test_that("resetRflomicsMAE", {
  
  ###use resetRflomicsMAE
  MAE0 <- RFLOMICS:::resetRflomicsMAE(MAE)
  
  expect_identical(getAnalysis(MAE0[[1]], "DiffExpAnal"),  list())
  expect_identical(getAnalysis(MAE0[[2]], "DiffExpAnal"),  list())
  expect_identical(getAnalysis(MAE0[[3]], "DiffExpAnal"),  list())
  
  expect_identical(getAnalysis(MAE0[[1]], "DiffExpEnrichAnal"),  list())
  expect_identical(getAnalysis(MAE0[[2]], "DiffExpEnrichAnal"),  list())
  expect_identical(getAnalysis(MAE0[[3]], "DiffExpEnrichAnal"),  list())
  
  expect_identical(getAnalysis(MAE0[[1]], "CoExpAnal"),  list())
  expect_identical(getAnalysis(MAE0[[2]], "CoExpAnal"),  list())
  expect_identical(getAnalysis(MAE0[[3]], "CoExpAnal"),  list())
  
  expect_identical(getAnalysis(MAE0[[1]], "CoExpEnrichAnal"),  list())
  expect_identical(getAnalysis(MAE0[[2]], "CoExpEnrichAnal"),  list())
  expect_identical(getAnalysis(MAE0[[3]], "CoExpEnrichAnal"),  list())
  
})

# ---- getAnalyzedDatasetNames ----
test_that("getAnalyzedDatasetNames", {
  
  names.list <- getAnalyzedDatasetNames(MAE)
  
  expect_identical(
    names(names.list),
    c("DataProcessing", "DiffExpAnal", "CoExpAnal"))
})

# ---- getLabs4plot ----
test_that("getLabs4plot", {
  
  Labs.list <- RFLOMICS:::getLabs4plot(MAE[[1]])
  expect_identical(Labs.list$title, "RNAtest: raw RNAseq data")
  expect_identical(Labs.list$x_lab, "RNAseq data")
  
  
  Labs.list <-  RFLOMICS:::getLabs4plot(MAE[[2]])
  expect_identical(Labs.list$title, "protetest: raw proteomics data")
  expect_identical(Labs.list$x_lab, "proteomics data")
})

# ---- rflomicsMAE2MAE ----
test_that("rflomicsMAE2MAE", {
  
  MAE_bis <- rflomicsMAE2MAE(MAE)
  expect_true(is(MAE_bis, "MultiAssayExperiment"))
  
})

test_that("Test diff plot", {
  
  co <- capture.output(
    p <- plotDiffAnalysis(
      MAE, SE.name = "protetest", 
      contrastName = "(temperatureMedium - temperatureLow) in mean"))
  expect_equal(names(p), c("MA.plot","Volcano.plot","Pvalue.hist" ))
  
  co <- capture.output(
    p <- plotHeatmapDesign(
      MAE, SE.name = "protetest", 
      contrastName = "(temperatureMedium - temperatureLow) in mean"))
  expect_true("HeatmapList" %in% is(p))
  
  co <- capture.output(
    p <- plotBoxplotDE(MAE, SE.name = "protetest", 
                       featureName = "AT1G01010", 
                       groupColor="groups",  raw = FALSE))
  expect_equal(is(p), "gg")
  
  co <- capture.output(
    p <- plotBoxplotDE(MAE, SE.name = "protetest", 
                       featureName = "", 
                       groupColor="groups",  raw = FALSE))
  expect_equal(is(p), "gg")
})


test_that("Test coseq plot", {
  
  p <- plotCoExpressionProfile(MAE, SE.name = "protetest") 
  expect_equal(is(p), "gg")
  
  p <- plotCoExpression(MAE, SE.name = "protetest") 
  expect_equal(names(p), c("profiles","boxplots","probapost_boxplots",
                           "probapost_barplots", "probapost_histogram",
                           "ICL", "logLike"))
  
  p <- plotCoseqContrasts(MAE, SE.name = "protetest")
  expect_equal(is(p), "gg")
  
  expect_equal(length(getCoexpClusters(MAE, SE.name = "protetest")), 6)
  
})

test_that("get summary analysis", {
  
  p <- getDiffAnalysesSummary(MAE, plot = TRUE)
  expect_equal(is(p), "gg")
  
  p <- getCoExpAnalysesSummary(MAE)
  expect_equal(is(p), "gg")
  
})

test_that("getters", {
  
  expect_equal(
    getDEList(
      MAE, SE.name = "protetest", 
      contrasts = "(temperatureElevated - temperatureLow) in imbibitionDS",
      operation = "intersection"),
    getDEList(
      MAE, SE.name = "protetest", 
      contrasts = "(temperatureElevated - temperatureLow) in imbibitionDS",
      operation = "union")
  )
  
  expect_equal(
    getValidContrasts(MAE, omicName="protetest"),
    metadata(MAE[["protetest"]])[["DiffExpAnal"]][["results"]][["Validcontrasts"]])
  
})

test_that(".generateEcoseedExampleData", {
  
  res <- RFLOMICS:::.generateEcoseedExampleData()
  expect_equal(names(res),
               c("projectName","ExpDesign","dF.List.ref","dF.Type.dFac",
                 "omicsNames","omicsTypes","omicsData"))
  
})

test_that(".getPackageInfo", {
  
  vers <- RFLOMICS:::.getPackageInfo(MAE, package = "coseq")
  expect_true(!is.null(vers))
})


test_that(".getKEGGRelease", {
  
  vers <- RFLOMICS:::.getKEGGRelease()
  
  if(!is.null(vers)){
    expect_true(!is.null(vers))
  }
  else{
    expect_null(vers)
  }
})

