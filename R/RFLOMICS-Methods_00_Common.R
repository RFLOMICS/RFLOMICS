### ============================================================================
### [00_common] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif

#' @import methods

# ---- resetRflomicsMAE ----
#' @title resetRflomicsMAE
#' @description
#' resetRflomicsMAE allows for initializing the object or initializing a
#' selection of results.
#' @param object An object of class \link{RflomicsMAE-class}
#' @param singleAnalyses vector of single omics analysis results names
#' (c("DataProcessing", "PCAlist", "DiffExpAnal",
#' "DiffExpEnrichAnal", "CoExpAnal", "CoExpEnrichAnal"))
#' @param multiAnalyses vector of multi omics analysis results names
#' (c("IntegrationAnalysis"))
#' @param datasetNames dataset name.
#' If dataset == NULL, all datasets will be reset
#' @return An object of class \link{RflomicsMAE-class}
#' @noRd
#' @keywords internal
setMethod(
  f = "resetRflomicsMAE",
  signature = "RflomicsMAE",
  definition = function(object,
                        singleAnalyses = NULL,
                        multiAnalyses  = NULL,
                        datasetNames   = NULL)
  {

    #default values
    all.datasets <- getDatasetNames(object)
    all.singleAnalyses <- c("DataProcessing",
                            "DiffExpAnal",
                            "DiffExpEnrichAnal",
                            "CoExpAnal",
                            "CoExpEnrichAnal")
    all.multiAnalyses <- c("IntegrationAnalysis")

    if(is.null(datasetNames) &&
       is.null(singleAnalyses) &&
       is.null(multiAnalyses)){
      datasetNames <- all.datasets
      singleAnalyses <- all.singleAnalyses
      multiAnalyses <- all.multiAnalyses
    }

    if(is.null(datasetNames) && !is.null(singleAnalyses))
      datasetNames <- all.datasets

    if(!is.null(datasetNames) && is.null(singleAnalyses))
      singleAnalyses <- all.singleAnalyses

    # for specific datasets
    if(!is.null(datasetNames)){

      datasetNames <- intersect(datasetNames, all.datasets)
      if(length(datasetNames) == 0)
        stop("The name of these data does not match the datasets available in the object.")

      singleAnalyses <- intersect(singleAnalyses, all.singleAnalyses)
      if(length(singleAnalyses) == 0)
        stop("The name of these analyses does not match the analyses available in the object.")

      # reset SO analysis
      for (data in datasetNames) {

        for (analysis in singleAnalyses) {

          metadata(object[[data]])[[analysis]] <-
            switch (
              analysis,
              "DataProcessing" = {
                genes_flt0 <-
                  metadata(object[[data]])[["DataProcessing"]][["rowSumsZero"]]

                list(
                  rowSumsZero      = genes_flt0,
                  selectedSamples  = colnames(object[[data]]),
                  featureFiltering = list(),
                  Normalization    = list(),
                  Transformation   = list(),
                  log = NULL)
              },
              {
                list()
              }
            )
        }
      }

      # reset MO analysis
      multiAnalyses <- all.multiAnalyses
    }

    # reset multi-omics analysis
    if(!is.null(multiAnalyses)){

      multiAnalyses <- intersect(multiAnalyses, all.multiAnalyses)
      if(length(multiAnalyses) == 0) stop("")

      for (analysis in multiAnalyses) {
        metadata(object)[[analysis]] <- list()
      }
    }

    return(object)
  })


# ---- generateReport ----
#' @title Generate RFLOMICS html report or archive
#' @description
#' This function is used to generate a html report from a
#' \link{RflomicsMAE-class} object or archive with results.
#' @param object a object of \link{RflomicsSE} class or
#' \link{RflomicsMAE-class} class.
#' @param reportName Name of the html report (default: date()_projectName.html).
#' @param archiveName name of archive with all analysis results
#' (default: date()_projectName.tar.gz).
#' @param tmpDir temporary directory (default: working directory)
#' @param ... other arguments to pass into the render function.
#' @return An html report or archive (tar.gz)
#' @importFrom rmarkdown render
#' @exportMethod generateReport
#' @rdname generateReport
#' @name generateReport
#' @aliases generateReport,RflomicsMAE-method
#' @example inst/examples/generateReport.R
setMethod(
  f          = "generateReport",
  signature  = "RflomicsMAE",
  definition = function(object,
                        reportName  = NULL,
                        archiveName = NULL,
                        tmpDir      = NULL,
                        ...) {

    # check analysis
    if(is.null(getAnalyzedDatasetNames(object)))
      stop("An exploratory analysis must be performed on at least one of the datasets.")

    if(is.null(reportName) && is.null(archiveName))
      stop("You must provide either a reportName or an archiveName.")

    projectName  <- getProjectName(object)

    # we need at least reportName or archiveName
    if(is.null(tmpDir)){
      if(!is.null(reportName))
        tmpDir <- dirname(reportName)
      else
        tmpDir <- dirname(archiveName)
    }

    # tmp dir
    if (file.access(tmpDir, 2) != 0)
      stop("No writing access in ", tmpDir)

    tmpDir <-
      file.path(tmpDir,
                paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName))
    dir.create(tmpDir, showWarnings = FALSE)

    # save MAE object in Rdata
    RDataName    <- paste0(projectName, ".MAE.RData")
    rflomics.MAE <- object
    save(rflomics.MAE, file = file.path(tmpDir, RDataName))

    # Set up parameters to pass to Rmd document
    param.list <-
      list(
        FEdata = file.path(tmpDir, RDataName),
        title  = paste0(projectName, " project"),
        rflomicsVersion = metadata(object)$rflomicsVersion,
        date = metadata(object)$date,
        outDir = tmpDir
      )

    # html name
    if(is.null(reportName))
      reportName <- file.path(
        tmpDir,
        paste0(format(Sys.time(), "%Y_%m_%d"), "_", projectName, ".html"))

    render(
      input             = system.file("RFLOMICSapp", "report.Rmd", package = "RFLOMICS"),
      output_file       = reportName,
      params            = param.list,
      knit_root_dir     = tmpDir,
      intermediates_dir = tmpDir,
      envir = new.env(parent = globalenv())
    )

    #Export results
    if (!is.null(archiveName)){

      # cp html in tmpDir
      file.copy(from = reportName, to = tmpDir)
      cmd <-
        paste0("tar -C ", dirname(tmpDir),
               " -czf ", archiveName, " ",
               basename(tmpDir))
      system(cmd)
      #message(cmd)

    } else{
      file.copy(from = reportName, to = dirname(tmpDir))
    }
    unlink(tmpDir, recursive = TRUE)
  }
)

## ---- get element from metadata slot from rflomicsSE/MAE ----

#' @title Get results from RFLOMICS object
#' @aliases getAnalysis,RflomicsMAE-method
#' @name getAnalysis
#' @rdname getAnalysis
#' @description
#' Get a specific analysis results from a Rflomcs MAE or a SE.
#' @return The analysis metadata slot (a list of results)
#' \itemize{
#'    \item getAnalysis: return list of results from a specific analysis.}
#' @param object The RflomicsMAE or RflomicsSE object from which to extract
#' the analysis.
#' @param name the name of element to add to metadata slot.
#' @param subName the name of sub element to add to metadata slot.
#' @exportMethod getAnalysis
setMethod(
  f = "getAnalysis",
  signature = "RflomicsMAE",
  definition = function(object, #SE.name = NULL,
                        name    = NULL,
                        subName = NULL){

    results <-
      .getAnalysis(object, name, subName)

    return(results)
  })

#' @rdname getAnalysis
#' @aliases getAnalysis,RflomicsSE-method
#' @name getAnalysis
#' @exportMethod getAnalysis
setMethod(
  f = "getAnalysis",
  signature = "RflomicsSE",
  definition = function(object,
                        name = NULL,
                        subName = NULL) {

    results <-
      .getAnalysis(object, name, subName)

    return(results)
  }
)

# ---- getAnalyzedDatasetNames ----
#' @rdname getAnalysis
#' @description
#' \itemize{
#'    \item getAnalyzedDatasetNames: return a list of performed analysis names.}
#' @param analyses vector of list of analysis name
#' @exportMethod getAnalyzedDatasetNames
#' @aliases getAnalyzedDatasetNames,RflomicsMAE-method
#' @name getAnalyzedDatasetNames
#' @examples
#' # See generateReport for an example that includes getAnalyzedDatasetNames
setMethod(
  f          = "getAnalyzedDatasetNames",
  signature  = "RflomicsMAE",
  definition = function(object, analyses = NULL) {

    all.analyses <- c("DataProcessing",
                      "DiffExpAnal", "DiffExpEnrichAnal",
                      "CoExpAnal", "CoExpEnrichAnal")

    if(is.null(analyses)) analyses <- all.analyses

    df.list <- list()
    for (dataset in getDatasetNames(object)) {

      if(is.null(object[[dataset]])) next

      for(analysis in analyses){

        if(length(metadata(object[[dataset]])[[analysis]]) == 0)
          next

        switch (
          analysis,
          "DataProcessing" = {
            if(isTRUE(metadata(object[[dataset]])[[analysis]]$done))
              df.list[[analysis]] <- c(df.list[[analysis]], dataset)
          },
          "DiffExpAnal" = {
            if(!is.null(getValidContrasts(object[[dataset]])))
              df.list[[analysis]] <- c(df.list[[analysis]], dataset)
          },
          "CoExpAnal" = {
            if(is.null(getAnalysis(object[[dataset]], name = analysis)$errors))
              df.list[[analysis]] <- c(df.list[[analysis]], dataset)
          },
          {
            for(db in names(metadata(object[[dataset]])[[analysis]])){
              df.list[[analysis]][[db]] <- c(df.list[[analysis]][[db]], dataset)
            }
          }
        )
      }
    }

    if(length(df.list) == 0) return(NULL)
    if(length(analyses) == 1) return(df.list[[1]])
    return(df.list)
  })

## ---- set element to metadata slot in rflomicsSE/MAE ----
#' @title setElementToMetadata
#' @description set element to metadata slot
#' @param object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results..
#' @param name the name of element to add to metadata slot.
#' @param subName the name of sub element to add to metadata slot.
#' @param content the content of element to add
#' @return An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}.
#' @keywords internal
#' @noRd
setMethod(
  f = "setElementToMetadata",
  signature = "RflomicsMAE",
  definition = function(object,
                        name = NULL,
                        subName = NULL,
                        content = NULL) {

    object <-
      .setElementToMetadata(object, name, subName, content)

    return(object)
  })

#' @keywords internal
#' @noRd
setMethod(
  f = "setElementToMetadata",
  signature = "RflomicsSE",
  definition = function(object,
                        name = NULL,
                        subName = NULL,
                        content = NULL) {

    object <-
      .setElementToMetadata(object, name, subName, content)

    return(object)
  })


## ---- getLabs4plot ----
#' @title getLabs4plot
#' @description get title, x label, y label for plot
#' @param object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results..
#' @param log .
#' @param processedData .
#' @param filtredData .
#' @return list of strings
#' @keywords internal
#' @noRd
setMethod(
  f = "getLabs4plot",
  signature = "RflomicsSE",
  definition = function(object) {

    labels <- list()
    omicsType <- getOmicsTypes(object)

    title <- NULL
    x_lab <- paste0(omicsType, " data")

    if(.isFiltered(object))
      title <- c(title,
                 paste0("filtred (", getFilterSettings(object)$method, ")"))

    if(.isTransformed(object) && getTransSettings(object)$method != "none"){
      method <- getTransSettings(object)$method
      title  <- c(title, paste0("transformed (", method, ")"))
    }

    if(.isNormalized(object) && getNormSettings(object)$method != "none"){
      method <- getNormSettings(object)$method
      title  <- c(title, paste0("normalized (", method, ")"))
    }

    if(!is.null(title)){
      title <- paste0(names(omicsType), ": ",
                      paste(title, collapse = " and "), " ",
                      omicsType, " data")
    } else {
      title <- paste0(names(omicsType), ": ",
                      "raw ", omicsType, " data")
    }

    if (!is.null(metadata(object)$DataProcessing$log)) {
      x_lab <- paste0(metadata(object)$DataProcessing$log, "(", omicsType, " data)")
    }

    labels$title <- title
    labels$x_lab <- x_lab

    return(labels)
  })


## ---- rflomicsMAE2MAE ----

#' @title convert rflomicsMAE to MAE
#' @aliases rflomicsMAE2MAE,RflomicsMAE-method
#' @name rflomicsMAE2MAE
#' @rdname rflomicsMAE2MAE
#' @description
#' Convert rflomicsMAE object to MultiAssayExperiment object.
#' @return object of MultiAssayExperiment class.
#' @param object The RflomicsMAE object to convert.
#' @param raw booleen. If TRUE raw omics data values.
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @exportMethod rflomicsMAE2MAE
setMethod(
  f = "rflomicsMAE2MAE",
  signature = "RflomicsMAE",
  definition = function(object, raw = FALSE){

    methods <- paste0(
      ""
    )

    # apply processing on matrix
    for(SE.name in names(object)){
      object[[SE.name]] <- getProcessedData(object[[SE.name]], norm = TRUE)
    }

    # contruct MAE from rflomicsMAE
    MAE <- MultiAssayExperiment(experiments = experiments(object),
                                colData     = colData(object),
                                sampleMap   = sampleMap(object),
                                metadata    =
                                  list(projectName = getProjectName(object),
                                       methods     = methods
                                  )
    )

    return(MAE)
  })

