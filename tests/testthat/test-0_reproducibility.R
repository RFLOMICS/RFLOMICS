### ============================================================================
### [00_reproducibility] Ensures the reproducibility of results in case of code 
### modifications or changes in package versions.
### ----------------------------------------------------------------------------
# N. Bessoltane, 


library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in the code or package versions affect the expected results. 

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

## diff analysis
MAE <- MAE |> 
  runDiffAnalysis(SE.name = "RNAtest", p.adj.method = "BH", 
                  method = "edgeRglmfit", p.adj.cutoff = 0.05, 
                  logFC.cutoff = 0) |>
  runDiffAnalysis(SE.name = "protetest", p.adj.method = "BH", 
                  method = "limmalmFit",  p.adj.cutoff = 0.05, 
                  logFC.cutoff = 0) |>
  runDiffAnalysis(SE.name = "metatest", p.adj.method = "BH", 
                  method = "limmalmFit",  p.adj.cutoff = 0.05, 
                  logFC.cutoff = 0)

MAE <- 
  setValidContrasts(MAE, omicName = "RNAtest",
                    contrastList = getSelectedContrasts(MAE[["RNAtest"]]))

MAE[["protetest"]] <- 
  setValidContrasts(MAE[["protetest"]], 
                    contrastList = getSelectedContrasts(MAE[["protetest"]])[1:6,])
MAE[["metatest"]] <- 
  setValidContrasts(MAE[["metatest"]], 
                    contrastList = getSelectedContrasts(MAE[["metatest"]])[1:6,])

## co expression
co <- capture.output(
  MAE <- MAE |> 
    #runCoExpression(SE.name = "RNAtest", 
    #contrastNames = "(temperatureMedium - temperatureLow) in imbibitionDS", 
    #K = 2:10, replicates = 5, merge = "union", model = "Normal", 
    #GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin",
    # normFactors = "TMM") |>
    runCoExpression(
      SE.name = "protetest",
      K = 2:10, 
      replicates = 5, 
      merge = "union", 
      model = "Normal", 
      min.data.size = 10)
) #|>
#runCoExpression(SE.name = "metatest", 
#contrastNames = "(temperatureMedium - temperatureLow) in imbibitionDS" ,
# K = 2:10, replicates = 5, merge = "union", model = "Normal")

# Enrichment
# Comment this when not locally used (need the org at tair package.)
MAE <- MAE |>
  # runAnnotationEnrichment(SE.name = "RNAtest", database = "GO",
  #                         domain = c("BP", "MF", "CC"),
  #                         list_args = list(OrgDb = "org.At.tair.db",
  #                                          keyType = "TAIR",
  #                                          pvalueCutoff = 0.05)) |>
  runAnnotationEnrichment(SE.name = "protetest", 
                          from = "CoExp",
                          database = "GO",
                          domain = c("BP"), 
                          pvalueCutoff = 0.05,
                          OrgDb = "org.At.tair.db", 
                          keyType = "TAIR")


# ---- test accessors ----
## ---- Rflomics class ----
test_that("test if RflomicsMAE / RflomicsSE", {
  expect_true(is(MAE, "RflomicsMAE"))
})

## ---- getDatasetNames ----
test_that("test getDatasetNames", {
  expect_true(all(getDatasetNames(MAE) %in% c("RNAtest", "metatest", "protetest")))
})

## ---- getOmicsTypes ----
test_that("test getOmicsTypes", {
  
  exp.res <- c("RNAseq", "proteomics", "metabolomics")
  names(exp.res) <- c("RNAtest", "protetest", "metatest")
  expect_equal(as.vector(getOmicsTypes(MAE)), as.vector(exp.res))
})

## ---- colData ----
test_that("colData", {
  
  Repeat      <- factor(rep(c("rep1", "rep2", "rep3"), 9), 
                        levels =c("rep1", "rep2", "rep3"))
  imbibition  <- factor(c(rep("DS", 9), rep("EI", 9), rep("LI", 9)), 
                        levels =c("DS", "EI", "LI"))    
  temperature <- factor(rep(c(rep("Low",3), rep("Medium", 3), rep("Elevated", 3)), 3),  
                        levels =c("Low", "Medium", "Elevated"))
  Repeat      <- relevel(as.factor(Repeat),      ref="rep1")
  temperature <- relevel(as.factor(temperature), ref="Low")
  imbibition  <- relevel(as.factor(imbibition),  ref="DS")
  
  # groups
  groups.level <- character(0)
  # Boucles imbriquées pour créer les combinaisons
  for (v1 in c("Low", "Medium", "Elevated")) {
    for (v2 in c("DS", "EI", "LI")) {
      groups.level <- c(groups.level, paste(v1, v2, sep = "_"))
    }
  }
  groups <- factor(paste(temperature, imbibition, sep = "_"),
                   levels =groups.level)
  
  # samples
  samples.level <- character(0)
  # Boucles imbriquées pour créer les combinaisons
  for (v1 in groups.level) {
    for (v2 in 1:3) {
      samples.level <- c(samples.level, paste(v1, v2, sep = "_"))
    }
  }
  
  samples <- factor(paste(temperature, imbibition, sub("rep", "", Repeat), sep = "_"),
                    levels = samples.level)
  
  
  colData <- data.frame(Repeat      = Repeat, 
                        groups      = groups,
                        temperature = temperature,
                        imbibition  = imbibition,
                        samples     = samples)
  
  rownames(colData) <- samples
  
  expect_equal(as.data.frame(getDesignMat(MAE)), as.data.frame(colData))
})


## ---- getRflomicsSE ----
test_that("test getRflomicsSE", {
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_true(is(SE, "RflomicsSE"))
  
  expect_null(getRflomicsSE(MAE))
  expect_null(getRflomicsSE(MAE, "toto"))
  
  expect_equal(getDatasetNames(SE), "RNAtest")
})


## ---- get design factors ----
test_that("test getFactorNames", {
  
  expect_equal(getFactorNames(MAE), c("Repeat", "temperature", "imbibition" ))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getFactorNames(SE), c("Repeat", "temperature", "imbibition" ))
})

test_that("test getFactorModalities", {
  
  expect_equal(getFactorModalities(MAE, "imbibition"), c("DS", "EI", "LI"))
  expect_equal(getFactorModalities(SE, "imbibition"),  c("DS", "EI", "LI"))
})



test_that("test getFactorTypes", {
  
  vec <- c("batch", "Bio", "Bio")
  names(vec) <- c("Repeat", "temperature", "imbibition" )
  expect_equal(getFactorTypes(MAE), vec)
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getFactorTypes(SE), vec)
})

test_that("test getBioFactors", {
  
  expect_equal(getBioFactors(MAE), c("temperature", "imbibition"))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getBioFactors(SE), c("temperature", "imbibition"))
})

test_that("test getBatchFactors", {
  
  expect_equal(getBatchFactors(MAE), c("Repeat"))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getBatchFactors(SE), c("Repeat"))
})

test_that("test getMetaFactors", {
  
  expect_null(getMetaFactors(MAE))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_null(getMetaFactors(SE))
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
  
  expect_equal(getDiffSettings(MAE[["RNAtest"]])$contrastCoef, Contrasts.Coeff)
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

# ---- diff analysis ----
# test_that("diff analysis", {
#     
#     ## RNAtest
#     stats.RNAseq <- data.frame(
#         All = c(2589, 5925, 1579, 5535, 9735, 4901, 1473, 5110,  817),
#         Up  = c(1323, 3017,  790, 2729, 4778, 2324,  568, 2453,  344),
#         Down= c(1266, 2908,  789, 2806, 4957, 2577,  905, 2657,  473)
#     ) %>% as.matrix()
#     rownames(stats.RNAseq) <- c(
#         "(temperatureMedium - temperatureLow) in imbibitionDS",                                                              
#         "(temperatureElevated - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureMedium) in imbibitionDS",
#         "(temperatureMedium - temperatureLow) in mean",
#         "(temperatureElevated - temperatureLow) in mean",
#         "(temperatureElevated - temperatureMedium) in mean",
#         "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS"
#     )
#     
#     expect_equal(MAE[["RNAtest"]]@metadata$DiffExpAnal$stats, stats.RNAseq)
#     
#     ## protetest
#     stats.prot <- data.frame(
#         All = c(36, 132,  53,  99, 238, 130,   0,   1,   0),
#         Up  = c(24,  67,  21,  59, 124,  47,   0,   1,   0),
#         Down= c(12,  65,  32,  40, 114,  83,   0,   0,   0)
#     ) %>% as.matrix()
#     rownames(stats.prot) <- c(
#         "(temperatureMedium - temperatureLow) in imbibitionDS",                                                              
#         "(temperatureElevated - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureMedium) in imbibitionDS",
#         "(temperatureMedium - temperatureLow) in mean",
#         "(temperatureElevated - temperatureLow) in mean",
#         "(temperatureElevated - temperatureMedium) in mean",
#         "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS"
#     )
#     
#     expect_equal(MAE[["protetest"]]@metadata$DiffExpAnal$stats, stats.prot)
#     
#     ## metatest
#     stats.met <- data.frame(
#         All = c(54, 67, 45, 73, 87, 69, 12, 18,  5),
#         Up  = c(41, 44, 23, 62, 57, 33,  6,  8,  1),
#         Down= c(13, 23, 22, 11, 30, 36,  6, 10,  4)
#     ) %>% as.matrix()
#     rownames(stats.met) <- c(
#         "(temperatureMedium - temperatureLow) in imbibitionDS",                                                              
#         "(temperatureElevated - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureMedium) in imbibitionDS",
#         "(temperatureMedium - temperatureLow) in mean",
#         "(temperatureElevated - temperatureLow) in mean",
#         "(temperatureElevated - temperatureMedium) in mean",
#         "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
#         "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS"
#     )
#     
#     expect_equal(MAE[["metatest"]]@metadata$DiffExpAnal$stats, stats.met)
# })

# ---- co-expression ----
# test_that("coseq analysis", {
#     
#     expect_equal(MAE[["RNAtest"]]@metadata$CoExpAnal$cluster.nb[[1]],   4)
#     expect_equal(MAE[["protetest"]]@metadata$CoExpAnal$cluster.nb[[1]], 3)
#     expect_equal(MAE[["metatest"]]@metadata$CoExpAnal$cluster.nb[[1]],  6) # false
# })

# ---- enrichment ----
# test_that("Enrichment analysis", {
#     
#     expect_true(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$GO$summary["BP"] == 1)
#     expect_true(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$GO$summary["MF"] == 6)
#     expect_true(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$GO$summary["CC"] == 11)
# })



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
    c("DataProcessing", "DiffExpAnal", "CoExpAnal", "CoExpEnrichAnal"))
  
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


test_that("Test explor plot", {
  
  p <- plotConditionsOverview(MAE)
  expect_equal(is(p), "gg")
  
  p <- plotDataOverview(MAE)
  expect_equal(is(p), "gg")
  
  p <- plotDataOverview(MAE, realSize = TRUE)
  expect_equal(is(p), "gg")
  
  p <- plotLibrarySize(MAE, SE.name = "RNAtest", raw = TRUE)
  expect_equal(is(p), "gg")
  expect_error(plotLibrarySize(MAE, SE.name = "protetest"))
  
  p <- plotDataDistribution(MAE, SE.name = "RNAtest", plot = "boxplot")
  expect_equal(is(p), "gg")
  
  p <- plotDataDistribution(MAE, SE.name = "protetest", plot = "density")
  expect_equal(is(p), "gg")
  
  p <- plotOmicsPCA(MAE, SE.name = "RNAtest")
  expect_equal(is(p), "gg")
  
  p <- plotExpDesignCompleteness(MAE, omicName = "RNAtest")
  expect_equal(is(p), "gg")
})


test_that("Test diff plot", {
  
  co <- capture.output(
    p <- plotDiffAnalysis(
      MAE, SE.name = "RNAtest", 
      contrastName = "(temperatureMedium - temperatureLow) in mean"))
  expect_equal(names(p), c("MA.plot","Volcano.plot","Pvalue.hist" ))
  
  co <- capture.output(
    p <- plotHeatmapDesign(
      MAE, SE.name = "RNAtest", 
      contrastName = "(temperatureMedium - temperatureLow) in mean"))
  expect_true("HeatmapList" %in% is(p))
  
  co <- capture.output(
    p <- plotBoxplotDE(MAE, SE.name = "RNAtest", 
                       featureName = "AT1G01010", 
                       groupColor="groups",  raw = FALSE))
  expect_equal(is(p), "gg")
  
  co <- capture.output(
    p <- plotBoxplotDE(MAE, SE.name = "RNAtest", 
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
    
  expect_equal(getDEList(
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

  

