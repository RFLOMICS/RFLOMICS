#' @import methods

#---- 00 common methods ----

setGeneric(
  name = "setElementToMetadata",
  def  = function(object,
                  name = NULL,
                  subName = NULL,
                  content = NULL)
    standardGeneric("setElementToMetadata")

)

setGeneric(
  name = "getAnalysis",
  def  = function(object,
                  name = NULL,
                  subName = NULL)
    standardGeneric("getAnalysis")

)

setGeneric(
  name = "resetRflomicsMAE",
  def  = function(object,
                  singleAnalyses = NULL,
                  multiAnalyses = NULL,
                  datasetNames = NULL)
    standardGeneric("resetRflomicsMAE")

)

setGeneric(
  name = "getAnalyzedDatasetNames",
  def  = function(object, analyses = NULL)
    standardGeneric("getAnalyzedDatasetNames")
)

setGeneric(
  name = "generateReport",
  def  = function(object,
                  reportName = NULL,
                  archiveName = NULL,
                  ...)
    standardGeneric("generateReport")

)

setGeneric(
  name = "getLabs4plot",
  def  = function(object, ...)
    standardGeneric("getLabs4plot")

)

setGeneric(
  name = "rflomicsMAE2MAE",
  def  = function(object,
                  raw = FALSE)
    standardGeneric("rflomicsMAE2MAE")

)


#---- 01 load data ----

setGeneric(
  name = "plotDataOverview",
  def  = function(object,
                  omicNames = NULL,
                  realSize = FALSE,
                  raw = FALSE,
                  completeCases = FALSE)
    standardGeneric("plotDataOverview")
)

setGeneric(
  name = "initRawRflomicsSE",
  def  = function(object, ...)
    standardGeneric("initRawRflomicsSE")
)

setGeneric(
  name = "getProjectName",
  def  = function(object)
    standardGeneric("getProjectName")
)

setGeneric(
  name = "subRflomicsMAE",
  def  = function(object, omicNames = NULL)
    standardGeneric("subRflomicsMAE")
)

setGeneric(
  name = "getFactorNames",
  def  = function(object)
    standardGeneric("getFactorNames")
)

setGeneric(
  name = "getFactorTypes",
  def  = function(object)
    standardGeneric("getFactorTypes")
)

setGeneric(
  name = "getBioFactors",
  def  = function(object)
    standardGeneric("getBioFactors")
)

setGeneric(
  name = "getBatchFactors",
  def  = function(object)
    standardGeneric("getBatchFactors")
)

setGeneric(
  name = "getMetaFactors",
  def  = function(object)
    standardGeneric("getMetaFactors")
)

setGeneric(
  name = "getDesignMat",
  def  = function(object)
    standardGeneric("getDesignMat")
)

setGeneric(
  name = "getDatasetNames",
  def  = function(object)
    standardGeneric("getDatasetNames")
)

setGeneric(
  name = "getOmicsTypes",
  def  = function(object)
    standardGeneric("getOmicsTypes")
)

setGeneric(
  name = "getRflomicsSE",
  def  = function(object, datasetName = NULL)
    standardGeneric("getRflomicsSE")
)

setGeneric(
  name = "getFactorModalities",
  def  = function(object, factorName)
    standardGeneric("getFactorModalities")
)

setGeneric(
  name = "plotConditionsOverview",
  def  = function(object, omicNames = NULL)
    standardGeneric("plotConditionsOverview")
)

#---- 02 stat setting ----

setGeneric(
  name = "generateModelFormulae",
  def  = function(object, updateModelFormula = NULL)
    standardGeneric("generateModelFormulae")
)

setGeneric(
  name = "setModelFormula",
  def  = function(object, modelFormula = NULL, ...)
    standardGeneric("setModelFormula")
)

setGeneric(
  name = "getModelFormula",
  def  = function(object, ...)
    standardGeneric("getModelFormula")
)

setGeneric(
  name = "updateModelFormula",
  def  = function(object, ...)
    standardGeneric("updateModelFormula")
)

setGeneric(
  name = "generateExpressionContrast",
  def  = function(object, contrastType = NULL, ...)
    standardGeneric("generateExpressionContrast")
)

setGeneric(
  name = "setSelectedContrasts",
  def  = function(object, contrastNames = NULL, ...)
    standardGeneric("setSelectedContrasts")
)

setGeneric(
  name = "getSelectedContrasts",
  def  = function(object, all = FALSE, ...)
    standardGeneric("getSelectedContrasts")
)

#---- 03 data processing ----

setGeneric(
  name = "subsetRflomicsSE",
  def  = function(object, bioFactor = NULL, level = NULL)
    standardGeneric("subsetRflomicsSE")
)

setGeneric(
  name = "runDataProcessing",
  def  = function(
    object,
    samples                 = NULL,
    MVencoding              = "NA",
    lowCountFilter = 
      list(method           = NULL,
           strategy         = NULL,
           cpmCutoff        = NULL),
    missingValueFilter = 
      list(
        method              = NULL,
        proportion          = NULL,
        nbCondition         = NULL),
    transform = list(method = NULL),
    normalize = list(method = NULL),
    impute    = list(method = NULL,
                     factor = NULL),
    ...)
  standardGeneric("runDataProcessing")
)

setGeneric(
  name = "filterLowAbundance",
  def  = function(
    object,
    method    = c("filterByExpr", "CPM", "none"),
    strategy  = NULL,
    cpmCutoff = 1)
    standardGeneric("filterLowAbundance")
)

setGeneric(
  name = "filterMissingValues",
  def  = function(
    object,
    method      = c("GlobalFiltering", "ConditionFiltering", "none"),
    proportion  = 0.5,
    nbCondition = 1)
    standardGeneric("filterMissingValues")
)

setGeneric(
  name = "runFeatureFiltering",
  def  = function(
    object,
    lowCountFilter = 
      list(method    = c("filterByExpr", "CPM", "none"),
           strategy  = NULL,
           cpmCutoff = 1),
    missingValueFilter    = 
      list(method      = c("GlobalFiltering", "ConditionFiltering", "none"),
           proportion  = 0.5,
           nbCondition = 1),
    ...)
    standardGeneric("runFeatureFiltering")
)

setGeneric(
  name = "runSampleFiltering",
  def  = function(object,
                  samples = NULL,
                  ...)
    standardGeneric("runSampleFiltering")
)

setGeneric(
  name = "runNormalization",
  def  = function(object,
                  method = NULL,
                  ...)
    standardGeneric("runNormalization")
)

setGeneric(
  name = "runTransformData",
  def  = function(object,
                  method = NULL,
                  ...)
    standardGeneric("runTransformData")
)


setGeneric(
  name = "runMVImputation",
  def  = function(object,
                  method = "minFeatureValue", 
                  factor = 0.1, ...)
    standardGeneric("runMVImputation")
)

setGeneric(
  name = "isProcessedData",
  def  = function(object,
                  filter = TRUE,
                  trans = TRUE,
                  norm = TRUE,
                  log = TRUE, ...)
    standardGeneric("isProcessedData")
)

setGeneric(
  name = "getTransSettings",
  def  = function(object, ...)
    standardGeneric("getTransSettings")
)

setGeneric(
  name = "getNormSettings",
  def  = function(object, ...)
    standardGeneric("getNormSettings")

)

setGeneric(
  name = "getImputSettings",
  def  = function(object, ...)
    standardGeneric("getImputSettings")
)

setGeneric(
  name = "getFilterSettings",
  def  = function(object, ...)
    standardGeneric("getFilterSettings")
)

setGeneric(
  name = "getFilteredFeatures",
  def  = function(object, ...)
    standardGeneric("getFilteredFeatures")
)

setGeneric(
  name = "getCoeffNorm",
  def  = function(object, ...)
    standardGeneric("getCoeffNorm")
)

setGeneric(
  name = "plotOmicsPCA",
  def  = function(object,
                  raw = c("raw", "norm"),
                  axes = c(1, 2),
                  groupColor = "groups",
                  ...)
    standardGeneric("plotOmicsPCA")
)

setGeneric(
  name = "runOmicsPCA",
  def  = function(object,
                  ncomp = 5,
                  ...)
    standardGeneric("runOmicsPCA")
)


setGeneric(
  name = "plotDataDistribution",
  def  = function(object, plot = "boxplot", raw = FALSE, ...)
    standardGeneric("plotDataDistribution")
)

setGeneric(
  name = "plotLibrarySize",
  def  = function(object, raw = FALSE, ...)
    standardGeneric("plotLibrarySize")
)

setGeneric(
  name = "plotMissingValues",
  def  = function(object, raw = FALSE, ...)
    standardGeneric("plotMissingValues")
)

setGeneric(
  name = "checkExpDesignCompleteness",
  def  = function(object, sampleList = NULL, ...)
    standardGeneric("checkExpDesignCompleteness")
)

setGeneric(
  name = "plotExpDesignCompleteness",
  def  = function(object, sampleList = NULL, ...)
    standardGeneric("plotExpDesignCompleteness")
)

setGeneric(
  name = "getSelectedSamples",
  def  = function(object, ...)
    standardGeneric("getSelectedSamples")
)


#---- 04 diff analysis ----

setGeneric(
  name = "runDiffAnalysis",
  def  = function(object,
                  contrastNames     = NULL,
                  method           = NULL,
                  p.adj.method     = "BH",
                  p.adj.cutoff     = 0.05,
                  logFC.cutoff     = 0,
                  splitByFactor    = NULL,
                  cmd              = FALSE,
                  ...)
    standardGeneric("runDiffAnalysis")
)

setGeneric(
  name = "generateContrastMatrix",
  def  = function(object, 
                  contrastNames = NULL, ...)
    standardGeneric("generateContrastMatrix")
)

setGeneric(
  name = "filterDiffAnalysis",
  def  = function(object,
                  analysisName = "all",
                  p.adj.cutoff = 0.05,
                  logFC.cutoff = 0,
                  ...)
    standardGeneric("filterDiffAnalysis")
)


setGeneric(
  name = "getDiffAnalysesSummary",
  def  = function(object,
                  analysisNames = "all",
                  plot = FALSE,
                  ylabelLength = 30,
                  nbMaxLabel = 20,
                  interface = FALSE)
    standardGeneric("getDiffAnalysesSummary")
)

## ----GET/SET---------

setGeneric(
  name = "getDiffStat",
  def  = function(object, analysisName = "all",...)
    standardGeneric("getDiffStat")
)

setGeneric(
  name = "getDEList",
  def  = function(object,
                  contrasts = NULL,
                  operation = "union",
                  analysisName = "all",
                  ...)
    standardGeneric("getDEList")
)

setGeneric(
  name = "getDEMatrix",
  def  = function(object, analysisName = "all", ...)
    standardGeneric("getDEMatrix")
)

setGeneric(
  name = "getDiffSettings",
  def  = function(object, analysisName = "all", ...)
    standardGeneric("getDiffSettings")
)


## ---- Plots ----

setGeneric(
  name = "plotDiffAnalysis",
  def  = function(object,
                  contrastName,
                  analysisName = "all",
                  typeofplots = c("MA.plot", "volcano", "histogram"),
                  ...)
    standardGeneric("plotDiffAnalysis")
)

setGeneric(
  name = "plotHeatmapDesign",
  def  = function(object,
                  contrastName,
                  analysisName = "all",
                  splitFactor="none",
                  title = "",
                  annotNames = NULL,
                  modalities = NULL,
                  drawArgs = list(),
                  heatmapArgs = list(),
                  ...)
    standardGeneric("plotHeatmapDesign")
)


setGeneric(
  name = "plotBoxplotDE",
  def  = function(object,
                  analysisName = "all",
                  featureName = NULL,
                  groupColor="groups",
                  raw = FALSE,
                  ...)
    standardGeneric("plotBoxplotDE")
)


setGeneric(
  name = "validateContrasts",
  def  = function(object,
                  analysisName  = "all",
                  contrastNames = NULL,
                  ...)
    standardGeneric("validateContrasts")
)

setGeneric(
  name = "getValidContrasts",
  def  = function(object, 
                  analysisName = "all",
                  contrastList = NULL, ...)
    standardGeneric("getValidContrasts")
)


#---- 05 co-expression ----

setGeneric(
  name = "runCoExpression",
  def  = function(object,
                  K                = 2:20,
                  replicates       = 5,
                  contrastNames    = NULL,
                  merge            = "union",
                  model            = "Normal",
                  GaussianModel    = NULL,
                  transformation   = NULL,
                  normFactors      = NULL,
                  meanFilterCutoff = NULL,
                  scale            = NULL,
                  min.data.size    = 100,
                  ...)
  standardGeneric("runCoExpression")
)

## ---- Plots ----

setGeneric(
  name = "plotCoExpressionProfile",
  def  = function(object, ...)
    standardGeneric("plotCoExpressionProfile")
)

setGeneric(
  name = "plotCoExpression",
  def  = function(object, ...)
    standardGeneric("plotCoExpression")
)

setGeneric(
  name = "plotCoseqContrasts",
  def  = function(object, ...)
    standardGeneric("plotCoseqContrasts")
)

## ---- GET/SET ----
setGeneric(
  name = "getCoexpSettings",
  def  = function(object, ...)
    standardGeneric("getCoexpSettings")
)

setGeneric(
  name = "getCoexpClusters",
  def  = function(object, clusterName = NULL, ...)
    standardGeneric("getCoexpClusters")
)

setGeneric(
  name = "getCoExpAnalysesSummary",
  def  = function(object, omicNames = NULL)
    standardGeneric("getCoExpAnalysesSummary")
)

#---- 06 annotation ----

setGeneric(
  name = "runAnnotationEnrichment",
  def  = function(object,
                  featureList = NULL,
                  from = "DiffExp",
                  universe = NULL,
                  database = "custom",
                  domain = "no-domain",
                  annotation = NULL,
                  OrgDb = NULL,
                  organism  = NULL,
                  keyType = NULL,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 1,
                  minGSSize = 10,
                  maxGSSize = 500,
                  ...)
  standardGeneric("runAnnotationEnrichment")
)

setGeneric(
  name = "plotClusterProfiler",
  def  = function(object,
                  featureListName = NULL,
                  database = NULL,
                  domain = "no-domain",
                  plotType = "dotplot",
                  showCategory = 15,
                  searchExpr = "",
                  nodeLabel = "all",
                  p.adj.cutoff = NULL,
                  ...)
    standardGeneric("plotClusterProfiler")

)


setGeneric(
  name = "plotEnrichComp",
  def = function(object,
                 from = "DiffExp",
                 database = NULL,
                 domain = "no-domain",
                 matrixType = "FC",
                 nClust = NULL,
                 ...)
    standardGeneric("plotEnrichComp")

)
setGeneric(
  name = "getEnrichRes",
  def  = function(object,
                  featureListName = NULL,
                  from = "DiffExp",
                  database = "GO",
                  domain = NULL, ...)

    standardGeneric("getEnrichRes")
)

setGeneric(
  name = "getEnrichSettings",
  def  = function(object,
                  from = "DiffExp",
                  database = "GO")

    standardGeneric("getEnrichSettings")
)

setGeneric(
  name = "sumORA",
  def  = function(object,
                  from = "DiffExp",
                  database = NULL,
                  featureListName = NULL)
    standardGeneric("sumORA")
)


setGeneric(
  name = "getAnnotAnalysesSummary",
  def  = function(object,
                  from         = "DiffExp",
                  listNames    = NULL,
                  omicNames    = NULL,
                  databaseList = NULL,
                  ...)
    standardGeneric("getAnnotAnalysesSummary")
)

#---- 07 integration ----

setGeneric(
  name = "prepareForIntegration",
  def  = function(object,
                  omicsNames = NULL,
                  rnaSeq_transfo = "limma (voom)",
                  variableLists = NULL,
                  group = NULL,
                  method = "MOFA",
                  cmd = FALSE,
                  ...)
    standardGeneric("prepareForIntegration")

)


setGeneric(
  name = "runOmicsIntegration",
  def  = function(object,
                  preparedObject = NULL,
                  method = "MOFA",
                  scale_views = FALSE,
                  maxiter = 1000,
                  num_factors = 10,
                  selectedResponse = NULL,
                  ncomp = 2,
                  link_datasets = 1,
                  link_response = 1,
                  sparsity = FALSE,
                  cases_to_try = 5,
                  cmd = FALSE,
                  ...)
  standardGeneric("runOmicsIntegration")

)

setGeneric(
  name = "getMixOmics",
  def  = function(object,
                  response = NULL,
                  onlyResults = TRUE)
    standardGeneric("getMixOmics")

)

setGeneric(
  name = "getMOFA",
  def  = function(object, onlyResults = TRUE)
    standardGeneric("getMOFA")

)

setGeneric(
  name = "getMOFASettings",
  def  = function(object)
    standardGeneric("getMOFASettings")

)

setGeneric(
  name = "setMOFA",
  def  = function(object, results = NULL)
    standardGeneric("setMOFA")

)

setGeneric(
  name = "setMixOmics",
  def  = function(object, results = NULL)
    standardGeneric("setMixOmics")

)

setGeneric(
  name = "sumMixOmics",
  def  = function(object, selectedResponse = NULL)
    standardGeneric("sumMixOmics")

)

setGeneric(
  name = ".getOneMORes",
  def  = function(object, selectedResponse)
    standardGeneric(".getOneMORes")

)


setGeneric(
  name = "getMixOmicsSettings",
  def  = function(object)
    standardGeneric("getMixOmicsSettings")

)

setGeneric(
  name = "plotMOVarExp",
  def  = function(object,
                  selectedResponse,
                  mode = NULL)
    standardGeneric("plotMOVarExp")

)





