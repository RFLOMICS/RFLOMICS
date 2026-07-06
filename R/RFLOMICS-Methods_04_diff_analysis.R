### ============================================================================
### [04_diff_analysis] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# D. Charif,
# N. Bessoltane,
# A. Hulot

##==== STAT METHOD ====

###==== METHOD to perform differential analysis ====

#' @title Run Differential Expression Analysis and process results
#' @name runDiffAnalysis
#' @aliases runDiffAnalysis,RflomicsSE-method
#' @description This is an interface method which run a differential
#' analysis on omics datasets stored in an object of class
#' \link{RflomicsSE} or \link{RflomicsMAE-class}. According to the type
#' of omics and to a list of contrasts, a differential analysis
#' is performed for each contrasts.
#' Two methods are available according to the type of object:
#' \itemize{
#' \item For RNAseq data:  the \code{glmFit} function/model of the
#' \code{edgeR} package is applied.
#' \item For proteomics and metabolomics data:  the \code{\link{lmFit}}
#' function/model of the \code{limma} package is applied.
#' }
#' @param object An object of class \link{RflomicsSE} or
#' class \link{RflomicsMAE-class}
#' @param SE.name SE.name the name of the dataset if the input object
#' is a \link{RflomicsMAE-class}
#' @param method A character vector giving the name of the differential
#' analysis method to run. Either "edgeRglmfit" or "limmalmFit".
#' @param contrastNames data.frame of contrast from generateExpressionContrast().
#' if NULL, it takes all selected contrasts.
#' @param p.adj.method The method chosen to adjust pvalue. Takes the same
#' values as the ones of adj.p.adjust method.
#' @param p.adj.cutoff The adjusted pvalue cut-off
#' @param logFC.cutoff The lof2FC cutoff
#' @param cmd Boolean. Used in the interface. If TRUE, print cmd for the user.
#' @param selectedModality apply differential analysis only on this modality.
#' @param ... Additional arguments.
#' @details
#' Functions and parameters used for RNAseq are those recommended in DiCoExpress
#' workflow (see the paper in reference).
#' Functions and parameters used for proteomics and metabolomics data are those
#' recommended in the (Efstathiou *et al.*, 2017)
#' @return A \link{RflomicsSE} or a \link{RflomicsMAE-class} object.
#' All the results are stored as a named list \code{DiffExpAnal}
#' in the metadata slot of a given \link{RflomicsSE} object.
#' Objects are:
#' \itemize{
#'   \item stats: data.frame giving a summary of the differential
#'  statistical analysis results by contrast:
#'  number of DE features, number of up and down regulated features
#'   \item setting: Parameters used for the differential analysis
#'  \item method: The method used for the differential analysis
#'  \item p.adj.method: The applied p-value correction method
#'   \item p.adj.cutoff: The cut-off applied for the adjusted p-value
#'   \item logFC.cutoff: The absolute log FC cut-off
#'   \item RawDEFres: a list giving for each contrast the raw results of
#'  the differential analysis method
#'   \item DEF: a list giving for each contrast a data.frame of non filtered
#'  differential expressed features with their statistics
#'   \item TopDEF: a list giving for each contrast a data.frame of
#'  differential expressed features ordered and filtered by p.adj.cutoff
#'  with their statistics
#'  \item mergeDEF: a data frame of 0 and 1 indicating for each features in row,
#'  if it is DE in a given contrasts in column
#'  \item contrasts: a data.table of the contrasts used for the differential
#'  analysis
#' }
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool
#' to process multifactorial RNAseq experiments from quality controls to
#' co-expression analysis through differential analysis based on contrasts
#' inside GLM models. Plant Methods 16, 68 (2020).
#'
#' Efstathiou G, Antonakis AN, Pavlopoulos GA, et al. ProteoSign:
#' an end-user online differential proteomics statistical analysis platform.
#' Nucleic Acids Res. 2017;45(W1):W300-W306.
#' @exportMethod runDiffAnalysis
#' @rdname runDiffAnalysis
#' @seealso \code{\link{getDiffSettings}}, \code{\link{getDEList}},
#'  \code{\link{getDEMatrix}}
#' @seealso \code{\link{plotDiffAnalysis}}, \code{\link{plotHeatmapDesign}},
#' \code{\link{plotBoxplotDE}}
#' @section Accessors:
#' @section Plots:
#' @example inst/examples/runDiffAnalysis.R
setMethod(
  f         = "runDiffAnalysis",
  signature = "RflomicsSE",
  definition = function(
    object,
    contrastNames    = NULL,
    method           = NULL,
    p.adj.method     = "BH",
    p.adj.cutoff     = 0.05,
    logFC.cutoff     = 0,
    splitByFactor    = NULL,
    cmd = FALSE,
    ...){
    
    # default methods
    default.methods <-
      switch (
        getOmicsTypes(object),
        "RNAseq" = "edgeRglmfit",
        "limmalmFit"
      )
    
    # check method
    if (is.null(method)) method <- default.methods
    if (isFALSE(method %in% default.methods))
      stop("The value '", method, "' is not supported for the argument 'method'.
           It is recommended to use the value '",default.methods[1],
           "' for '",getOmicsTypes(object), "' data")
    
    # check modelFormula
    if(is.null(getModelFormula(object)))
      stop("You need to set the model formula. (See setModelFormula())")
    
    # split data
    objects <- list()
    if(is.null(splitByFactor) || splitByFactor == "all"){
      
      splitLevels   <- "all"
      analysisNames <- "all" 
      objects[["all"]] <- object
    }
    else{
      bioFactors <- getBioFactors(object)
      
      if(!splitByFactor %in% bioFactors) 
        stop("No bio factor called ",splitByFactor," found.")
      
      splitLevels <- 
        getFactorModalities(object, factorName = splitByFactor)
      
      analysisNames <- vector()
      for(splitLevel in splitLevels){
        
        analysisName <- paste(splitByFactor, splitLevel, sep = "")
        
        objects[[analysisName]] <-
          subsetRflomicsSE(object, bioFactor = splitByFactor, level = splitLevel)
        
        analysisNames <- c(analysisNames, analysisName)
      }
    }
    
    # apply diff
    DiffExpAnals <- list()
    for(analysisName in analysisNames){
      
      object.min <- objects[[analysisName]]
      
      # check model
      modelFormula <- getModelFormula(object.min)
      if(is.null(modelFormula))
        stop("You need to set the model formula. (See setModelFormula())")
      
      # check contrast
      # If contrast names are provided, 
      # check whether they belong to the available contrasts.
      # Otherwise, use the selected contrasts if any are selected; 
      # if not, use all available contrasts.    
      allContrasts <- getSelectedContrasts(object.min, all = TRUE)
      if(nrow(allContrasts) == 0)
        stop("")
        
      if(is.null(contrastNames)){
        contrastNames <- 
          getSelectedContrasts(object.min, all = FALSE)$contrastName
        if(length(contrastNames) == 0)
          contrastNames <- allContrasts$contrastName
      }
      else
      {
        if(any(!contrastNames %in% allContrasts$contrastName))
          stop("The contrast names do not match the available contrasts.")
      }
      
      # define result output
      DiffExpAnal <- list(
        settings = list(),
        results  = list(),
        errors   = NULL
      )
      
      ## check completness
      Completeness <- checkExpDesignCompleteness(object.min)
      if (isTRUE(Completeness[["error"]])){
        DiffExpAnal[["errors"]] <- Completeness[["messages"]]
      
      }else{
        
        ## getcontrast
        contrastMatrix <-
          generateContrastMatrix(
            object.min, 
            contrastNames = contrastNames
            )
        
        ListRes <-
          switch(
            method,
            "edgeRglmfit" =
              .tryRflomics(
                .edgeRAnaDiff(
                  object          = object.min,
                  modelFormula    = modelFormula, 
                  Contrasts.Coeff = contrastMatrix,
                  FDR             = 1,
                  cmd             = cmd)
              ),
            "limmalmFit" =
              .tryRflomics(
                .limmaAnaDiff(
                  object          = object.min,
                  modelFormula    = modelFormula, 
                  Contrasts.Coeff = contrastMatrix,
                  p.adj.cutoff    = 1,
                  p.adj.method    = p.adj.method,
                  cmd             = cmd)
              )
          )
        
        if (!is.null(ListRes$error)) {
          DiffExpAnal[["errors"]] <- ListRes$error
          
        }else if (!is.null(ListRes$value)) {
          
          if (!is.null(ListRes$value[["RawDEFres"]]))
            DiffExpAnal[["results"]][["RawDEFres"]] <- ListRes$value[["RawDEFres"]]
          
          if (!is.null(ListRes$value[["ErrorList"]]))
            DiffExpAnal[["results"]][["runErrors"]] <- ListRes$value[["ErrorList"]]
          
          if (!is.null(ListRes$value[["DEF"]]))
            DiffExpAnal[["results"]][["DEF"]] <- ListRes$value[["DEF"]]
          
        }else{
          
          DiffExpAnal[["errors"]] <- "Something is not working correctly"
        }
      }
      # set settings
      DiffExpAnal[["settings"]][["method"]]           <- method
      DiffExpAnal[["settings"]][["Model.formula"]]    <- modelFormula
      selContrasts <- 
        allContrasts[allContrasts$contrastName %in% contrastNames,]
      DiffExpAnal[["settings"]][["Contrasts.Sel"]]    <- selContrasts
      DiffExpAnal[["settings"]][["contrastCoef"]]     <- contrastMatrix
      DiffExpAnal[["settings"]][["p.adj.method"]]     <- p.adj.method
      DiffExpAnal[["settings"]][["p.adj.cutoff"]]     <- p.adj.cutoff
      DiffExpAnal[["settings"]][["abs.logFC.cutoff"]] <- logFC.cutoff
      
      # 
      object <-
        setElementToMetadata(object  = object,
                             name    = "DiffExpAnal",
                             subName = analysisName,
                             content = DiffExpAnal)
      
      ## filtering
      object <-
        filterDiffAnalysis(object        = object,
                           analysisName  = analysisName,
                           p.adj.cutoff  = p.adj.cutoff,
                           logFC.cutoff  = logFC.cutoff)
    }
    
    # initiate analysis results
    object <-
      setElementToMetadata(object  = object,
                           name    = "DiffExpEnrichAnal",
                           content = list())
    object <-
      setElementToMetadata(object  = object,
                           name    = "CoExpAnal",
                           content = list())
    object <-
      setElementToMetadata(object  = object,
                           name    = "CoExpEnrichAnal",
                           content = list())

    return(object)
  })

#' @name runDiffAnalysis
#' @aliases runDiffAnalysis,RflomicsMAE-method
#' @rdname runDiffAnalysis
#' @exportMethod runDiffAnalysis
setMethod(
  f          = "runDiffAnalysis",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        contrastNames    = NULL,
                        method           = NULL,
                        p.adj.method     = "BH",
                        p.adj.cutoff     = 0.05,
                        logFC.cutoff     = 0,
                        splitByFactor    = NULL,
                        cmd              = FALSE,
                        ...){
    
    # all verifications are done in this method
    object[[SE.name]] <-
      runDiffAnalysis(object           = object[[SE.name]],
                      contrastNames    = contrastNames,
                      p.adj.method     = p.adj.method,
                      method           = method,
                      p.adj.cutoff     = p.adj.cutoff,
                      logFC.cutoff     = logFC.cutoff,
                      splitByFactor    = splitByFactor,
                      cmd = cmd)

    object <-
      setElementToMetadata(object,
                           name    = "IntegrationAnalysis",
                           content =  list())

    return(object)
  })

###==== METHOD to generateContrastMatrix ====

#' @rdname runDiffAnalysis
#' @name generateContrastMatrix
#' @aliases generateContrastMatrix,RflomicsSE-method
#' @description
#' \itemize{
#'    \item generateContrastMatrix:
#'  Defines contrast matrix or contrast list with contrast
#'  name and contrast coefficients}
#' @param contrastNames a data.frame of contrasts generated by
#' \link{generateExpressionContrast}
#' @return contrast matrix
#' @importFrom stats formula terms.formula
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
#' @noRd
#' @keywords internal
setMethod(
  f          = "generateContrastMatrix",
  signature  = "RflomicsSE",
  definition = function(object, 
                        contrastNames = NULL){

    modelFormula <- getModelFormula(object)
    if(is.null(modelFormula))
      stop("You need to set the model formula. (See generateModelFormulae() and 
           setModelFormula())")
      
    # get possible contrasts
    contrastList <- getSelectedContrasts(object, all = TRUE)
    if(nrow(contrastList) == 0)
      stop("No contrasts found in the object. (See generateModelFormulae() and 
           setModelFormula())")
    
    if(is.null(contrastNames)){
      contrastList  <- getSelectedContrasts(object, all = FALSE)
      
      if(nrow(contrastList) == 0)
        stop("No selected contrasts found in the object. (See setSelectedContrasts)")
    }
    else{
      contrastList <- contrastList[contrastList$contrastName %in% contrastNames,]
      if(nrow(contrastList) == 0)
        stop("None of the provided contrast names match the available contrasts.")
    }

    ExpDesign <- getDesignMat(object)
    factorBio <- getBioFactors(object)

    Contrasts.Coeff <-
      .getContrastMatrixF(ExpDesign    = ExpDesign,
                          factorBio    = factorBio,
                          modelFormula = modelFormula,
                          contrastList = contrastList$contrast)

    rownames(Contrasts.Coeff) <- contrastList$contrastName

    return(Contrasts.Coeff)
  })

###==== METHOD to filter differential analysis ====

#' @rdname runDiffAnalysis
#' @name filterDiffAnalysis
#' @aliases filterDiffAnalysis,RflomicsSE-method
#' @description
#' \itemize{
#'    \item filterDiffAnalysis: The filterDiffAnalysis method allows
#'    filtering the results of the differential analysis based on a
#'    new cutoff for p-value and fold change.}
#' @param p.adj.cutoff adjusted pvalue cutoff. Default is the parameter from
#' the differential analysis.
#' @param logFC.cutoff cutoff for absolute value of log2FC. Default is the
#' parameter from the differential analysis.
#' @param analysisName names of differential analyses.
#' @exportMethod filterDiffAnalysis
#' @importFrom data.table data.table
#' @importFrom purrr reduce
setMethod(
  f          = "filterDiffAnalysis",
  signature  = "RflomicsSE",
  definition = function(object,
                        analysisName = "all",
                        p.adj.cutoff = 0.05,
                        logFC.cutoff = 0){
    
    DiffExpAnal <- 
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)
    
    if (length(DiffExpAnal) == 0)
      stop("can't filter the DiffExpAnal object because it doesn't exist")
    
    if (is.null(p.adj.cutoff))
      p.adj.cutoff <- getDiffSettings(object[[SE.name]], analysisName = analysisName)$p.adj.cutoff
    
    if (is.null(logFC.cutoff))
      logFC.cutoff <- getDiffSettings(object[[SE.name]], analysisName = analysisName)$abs.logFC.cutoff
    
    if (is.null(DiffExpAnal[["results"]][["DEF"]]))
      stop("can't filter the DiffExpAnal object because it doesn't exist")
    
    # remplacera à terme les lignes ci-dessus
    DiffExpAnal[["settings"]][["p.adj.cutoff"]]      <- p.adj.cutoff
    DiffExpAnal[["settings"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
    DiffExpAnal[["results"]][["mergeDEF"]] <- NULL
    DiffExpAnal[["results"]][["TopDEF"]] <- NULL
    DiffExpAnal[["results"]][["stats"]] <- NULL
    DiffExpAnal[["results"]][["isNA"]] <- NULL
    DiffExpAnal[["errors"]] <- NULL
    
    ## TopDEF: Top differential expressed features
    stat.vec <- list()
    DEF_list <- data.frame(DEF = vector())
    for (x in names(DiffExpAnal[["results"]][["DEF"]])) {
      tab  <- DiffExpAnal[["results"]][["DEF"]][[x]]
      
      # Changing DEF if NA are encountered
      if (any(is.na(tab$Adj.pvalue))) {
        DiffExpAnal[["results"]][["DEF"]][[x]]  <- NULL
        DiffExpAnal[["results"]][["DEF"]][[x]]  <- tab <- tab[which(!is.na(tab$Adj.pvalue)),]
        DiffExpAnal[["results"]][["isNA"]][[x]] <- setdiff(rownames(object), rownames(DiffExpAnal[["results"]][["DEF"]][[x]]))
      }
      
      keep <- (tab$Adj.pvalue < p.adj.cutoff) & (abs(tab$logFC) > logFC.cutoff)
      tab  <- tab[keep,]

      if(nrow(tab) != 0){
        DiffExpAnal[["results"]][["TopDEF"]][[x]] <- tab
        
        # prepare to merge
        tmp      <- data.frame(DEF = rownames(tab))
        tmp[[x]] <- rep(1, length(rownames(tab)))
        DEF_list <- merge(DEF_list, tmp, by = "DEF", all = TRUE)
      }
      
      stat.vec[[x]] <- c("All"  = nrow(tab),
                         "Up"   = nrow(tab[tab$logFC >= 0,]),
                         "Down" = nrow(tab[tab$logFC <  0,])
      )
    }
    DEF_list[is.na(DEF_list)] <- 0
    
    if(nrow(DEF_list) != 0){
      DiffExpAnal[["results"]][["mergeDEF"]] <- DEF_list
    }else{
      DiffExpAnal[["errors"]] <-
        paste0("None of the contrasts yield any results.",
               " No differentially expressed features!")
    }
    
    DiffExpAnal[["results"]][["stats"]] <- do.call("rbind", stat.vec)
    
    object <-
      setElementToMetadata(object,
                           name    = "DiffExpAnal",
                           subName = analysisName,
                           content = DiffExpAnal)
    
    return(object)
  })

#' @rdname runDiffAnalysis
#' @name filterDiffAnalysis
#' @aliases filterDiffAnalysis,RflomicsMAE-method
#' @exportMethod filterDiffAnalysis
setMethod(
  f          = "filterDiffAnalysis",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        analysisNames = NULL,
                        p.adj.cutoff = 0.05,
                        logFC.cutoff = 0){

    if (!SE.name %in% names(object))
      stop(SE.name, " isn't the name of an experiment in ", object)

    object[[SE.name]] <-  filterDiffAnalysis(object = object[[SE.name]],
                                             analysisNames = analysisNames,
                                             p.adj.cutoff  = p.adj.cutoff,
                                             logFC.cutoff  = logFC.cutoff)

    return(object)
  })

###==== validateContrasts : (after differential analysis) ====

#' @rdname runDiffAnalysis
#' @name validateContrasts
#' @aliases validateContrasts,RflomicsSE-method
#' @description
#' \itemize{
#'    \item validateContrasts: Set the valid contrasts stored in \code{metadata} slot.}
#' @param contrastNames A data.frame of contrast
#' @exportMethod validateContrasts
setMethod(
  f          = "validateContrasts",
  signature  = "RflomicsSE",
  definition = function(object,
                        analysisName  = "all",
                        contrastNames = NULL){

    diffSettings <- getDiffSettings(object, analysisName = analysisName)
    if(is.null(diffSettings$method))
      stop("")
    
    unselectedContrasts <-
      contrastNames[!contrastNames %in% diffSettings$Contrasts.Sel$contrastName]

    if(length(unselectedContrasts) != 0)
      stop("These contrasts ", paste0(unselectedContrasts, collapse = ", "),
           " are not recognized.")
    
    DiffExpAnal <- 
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)
    DiffExpAnal[["settings"]][["validContrasts"]] <- contrastNames
    object <- 
      setElementToMetadata(
        object, 
        name = "DiffExpAnal", 
        subName = analysisName, 
        content = DiffExpAnal
      )

    return(object)
  })

#' @rdname runDiffAnalysis
#' @name validateContrasts
#' @aliases validateContrasts,RflomicsMAE-method
#' @param omicName a dataset name
#' @exportMethod validateContrasts
setMethod(
  f          = "validateContrasts",
  signature  = "RflomicsMAE",
  definition <- function(object,
                         omicName=NULL,
                         analysisName = "all",
                         contrastNames=NULL){

    if(!omicName %in% names(object))
      stop("This data name, ", omicName, ", does not exist in the your object")

    object[[omicName]] <-
      validateContrasts(object[[omicName]], 
                        analysisName = analysisName,
                        contrastList = contrastList)

    return(object)
  })

##==== GRAPHICAL METHOD ====

###==== Method to plot results of a differential analysis ====

#' @name plotDiffAnalysis
#' @aliases plotDiffAnalysis,RflomicsSE-method
#' @rdname runDiffAnalysis
#' @section Plots:
#' \itemize{
#'    \item plotDiffAnalysis method draws a MAplot, a volcano plot and the
#' p-values distribution from the results of a differential analysis.}
#' @param contrastName The contrastName for which the plots has to be drawn
#' @param typeofplots The plots you want to return. Default is all possible
#' plots: MA plot, Volcano plot and non adjusted pvalues histogram.
#' @param analysisName differiential analysis name
#' @exportMethod plotDiffAnalysis
#' @export
#' @examples
#' # See runDiffAnalysis for an example that includes plotDiffAnalysis
setMethod(
  f="plotDiffAnalysis",
  signature="RflomicsSE",
  definition <- function(object,
                         contrastName,
                         analysisName = "all",
                         typeofplots = c("MA.plot", "volcano", "histogram")){

    plots <- list()
    DiffExpAnal <-
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)

    resTable <- DiffExpAnal[["results"]][["DEF"]][[contrastName]]

    logFC.cutoff <- getDiffSettings(object, analysisName = analysisName)[["abs.logFC.cutoff"]]
    p.adj.cutoff <- getDiffSettings(object, analysisName = analysisName)[["p.adj.cutoff"]]

    if ("MA.plot" %in% typeofplots)
      plots[["MA.plot"]] <-
      .plotMA(data = resTable,
              p.adj.cutoff = p.adj.cutoff,
              logFC.cutoff = logFC.cutoff,
              contrastName=contrastName)
    if ("volcano" %in% typeofplots)
      plots[["Volcano.plot"]]   <-
      .plotVolcanoPlot(data = resTable,
                       p.adj.cutoff = p.adj.cutoff,
                       logFC.cutoff = logFC.cutoff,
                       contrastName=contrastName)
    if ("histogram" %in% typeofplots)
      plots[["Pvalue.hist"]]  <-
      .plotPValue(data =resTable, contrastName=contrastName)

    return(plots)
  })

#' @rdname runDiffAnalysis
#' @name plotDiffAnalysis
#' @aliases plotDiffAnalysis,RflomicsMAE-method
#' @exportMethod plotDiffAnalysis
setMethod(
  f          = "plotDiffAnalysis",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        contrastName,
                        analysisName = "all",
                        typeofplots = c("MA.plot", "volcano", "histogram")){

    return(plotDiffAnalysis(object = object[[SE.name]],
                            contrastName = contrastName,
                            analysisName = analysisName,
                            typeofplots = typeofplots))
  })

###==== Method to plot plotHeatmapDesign ====

#' @name plotHeatmapDesign
#' @aliases plotHeatmapDesign,RflomicsSE-method
#' @rdname runDiffAnalysis
#' @section Plots:
#' \itemize{
#'    \item plotHeatmapDesign method draws a heatmap from the results of a
#' differential analysis.}
#' @param contrastName The contrastName for which the MAplot has to be drawn
#' @param splitFactor characters. Default to none. Name of a feature in the
#' design matrix, splits the samples on the heatmap according to its modalities.
#' @param title characters. Title of the heatmap.
#' @param annotNames vector. Names of the annotations to keep in the Heatmap.
#' Default takes all available information.
#' @param modalities named list of vectors of modalities to subset and print
#' on the heatmap.
#' @param drawArgs,heatmapArgs  named lists. Any additional parameter passed
#' to ComplexHeatmap::Heatmap or ComplexHeatmap::draw
#' @exportMethod plotHeatmapDesign
#' @export
#' @importFrom tidyselect any_of
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importClassesFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importMethodsFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom grid gpar
#' @examples
#' # See runDiffAnalysis for an example that includes plotHeatmapDesign
setMethod(
  f          = "plotHeatmapDesign",
  signature  = "RflomicsSE",
  definition = function(object,
                        contrastName,
                        analysisName = "all",
                        splitFactor="none",
                        title = "",
                        annotNames = NULL,
                        modalities = NULL,
                        drawArgs = list(),
                        heatmapArgs = list()){
    

    Groups      <- getDesignMat(object)
    DiffExpAnal <- 
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)

    if (is.null(DiffExpAnal[["results"]][["TopDEF"]][[contrastName]])) {
      stop("There are no differentially expressed features")
    }

    resTable <- arrange(DiffExpAnal[["results"]][["TopDEF"]][[contrastName]],
                        Adj.pvalue)

    if (dim(resTable)[1] == 0) {
      stop("There are no differentially expressed features")
    }

    if (dim(resTable)[1] > 2000) {
      message("The number of differentially expressed variables is exceeding 2000, ",
              "only the first 2000 will be displayed")
      resTable <- resTable[seq_len(2000),]
      title <- ifelse(title == "",
                      paste0(title, "TOP 2000 DE ", .omicsDic(object)$variableName),
                      paste0(title, "\nTOP 2000 DE ", .omicsDic(object)$variableName))
    }

    # object2 <- .checkTransNorm(object, raw = FALSE)
    m.def  <- assay(object)

    m.def <- as.data.frame(m.def) %>%
      select(any_of(Groups$samples))

    # filter by DE
    m.def.filter <- subset(m.def,
                           rownames(m.def) %in% row.names(resTable))

    # normalize count

    # Center
    m.def.filter.center <- t(scale(t(m.def.filter), center = TRUE, scale = FALSE))

    # Annotations datatable
    df_annotation <- Groups %>% select(!samples & !groups)
    df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),]

    # Subset the dataset to print only interesting modalities
    if (!is.null(modalities)) {
      if (is.null(names(modalities))) {
        message("In heatmapPlot, modalities argument needs a named list. Not subsetting")
      }else{
        samplesToKeep <-
          Reduce(
            "intersect",
            lapply(seq_len(length(modalities)), function(i){
              col_nam <- names(modalities)[i]
              rownames(df_annotation[which(df_annotation[[col_nam]] %in% modalities[[i]]),])
            }
            ))

        df_annotation <-
          df_annotation[which(rownames(df_annotation) %in% samplesToKeep),]
        m.def.filter.center <-
          m.def.filter.center[, which(colnames(m.def.filter.center) %in% samplesToKeep)]

        df_annotation <-
          df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),]
      }
    }

    # Split management
    column_split.value <-
      if (splitFactor != "none") { df_annotation[, splitFactor] } else { NULL }

    # Select the right columns
    if (!is.null(annotNames)) {
      df_annotation <- df_annotation %>%
        select(any_of(annotNames))
    }

    # Color annotations
    nAnnot <- ncol(df_annotation)
    selectPal <- rownames(brewer.pal.info)[seq_len(nAnnot)]

    color_list <- lapply(seq_len(nAnnot),
                         FUN = function(i){
                           annot_vect <- unique(df_annotation[,i])

                           col_vect <-  colorRampPalette(
                             brewer.pal(n = min(length(annot_vect), 8),
                                        name = selectPal[i])
                           )(length(annot_vect))
                           names(col_vect) <- annot_vect
                           col_vect[!is.na(names(col_vect))]
                         })
    names(color_list) <- colnames(df_annotation)

    column_ha <- HeatmapAnnotation(df = df_annotation,
                                   col = color_list)

    namArg <- ifelse(getOmicsTypes(object) == "RNAseq",
                     "normalized counts", "XIC")

    # Arguments for Heatmap
    heatmapArgs <- c(
      list(matrix = m.def.filter.center,
           name = namArg,
           show_row_names = ifelse(dim(m.def.filter.center)[1] > 50, FALSE, TRUE),
           row_names_gp = gpar(fontsize = 8),
           column_names_gp = gpar(fontsize = 12),
           row_title_rot = 0 ,
           clustering_method_columns = "ward.D2",
           cluster_column_slices = FALSE,
           column_split = column_split.value,
           top_annotation = column_ha,
           column_title = title),
      heatmapArgs)

    # Arguments for drawing the heatmap
    drawArgs <- c(list(merge_legends = TRUE),
                  drawArgs)

    # Drawing heatmap in a null file to not plot it
    pdf(file = NULL)
    ha <- do.call(Heatmap, heatmapArgs)

    drawArgs$object <- ha
    ha <- do.call(draw, drawArgs)

    dev.off()

    return(ha)
  })

#' @name plotHeatmapDesign
#' @aliases plotHeatmapDesign,RflomicsMAE-method
#' @rdname runDiffAnalysis
#' @exportMethod plotHeatmapDesign
setMethod(
  f          = "plotHeatmapDesign",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        contrastName,
                        analysisName = "all",
                        splitFactor="none",
                        title = "", annotNames = NULL,
                        modalities = NULL,
                        drawArgs = list(),
                        heatmapArgs = list()){

    return(plotHeatmapDesign(object        = object[[SE.name]],
                             contrastName  = contrastName,
                             analysisName  = analysisName,
                             splitFactor   = splitFactor,
                             title         = title,
                             annotNames   = annotNames,
                             modalities   = modalities,
                             drawArgs     = drawArgs,
                             heatmapArgs  = heatmapArgs))

  })


###==== Method to plot BoxplotDE ====

#' @name plotBoxplotDE
#' @aliases plotBoxplotDE,RflomicsSE-method
#' @rdname runDiffAnalysis
#' @section Plots:
#' \itemize{
#'    \item plotBoxplotDE method draws a boxplot showing the expression of
#'    given differentially expressed feature.}
#' @param featureName variable name (gene/protein/metabolite name)
#' @param groupColor default to groups, indicate a variable in the design to
#' color the boxplots accordingly.
#' @param raw Boolean. Plot the raw data or the transformed ones (TRUE)
#' @exportMethod plotBoxplotDE
#' @examples
#' # See runDiffAnalysis for an example that includes plotBoxplotDE
setMethod(
  f          = "plotBoxplotDE",
  signature  = "RflomicsSE",
  definition = function(object,
                        featureName = NULL,
                        groupColor="groups",
                        raw = FALSE){

    # check variable name
    if (is.null(featureName) || featureName == "" || length(featureName) != 1) {
      message("set variable name")

      p <- ggplot() + theme_void() + ggtitle("set variable name")
      return(p)
    }

    if(isTRUE(raw)) object <- initRawRflomicsSE(object)

    Groups <- getDesignMat(object)

    # check presence of variable in SE
    object.DE <- tryCatch(object[featureName], error = function(e) e)
    if (!is.null(object.DE$message)) {
      message(object.DE$message)

      p <- ggplot() + theme_void() + ggtitle(object.DE$message)
      return(p)
    }

    Labels <- getLabs4plot(object)
    title <- paste0(Labels$title, "\n", "Feature: ", featureName)
    x_lab <- Labels$x_lab

    pseudo <- assay(object.DE)

    pseudo.gg <- melt(pseudo)
    colnames(pseudo.gg) <- c("featureName", "samples", "value")

    pseudo.gg <- pseudo.gg %>%
      full_join(Groups, by = "samples") %>%
      arrange(groups)

    pseudo.gg <- arrange(pseudo.gg, get(groupColor))

    pseudo.gg$groups <- factor(pseudo.gg$groups,
                               levels = unique(pseudo.gg$groups))

    p <-  ggplot(pseudo.gg,
                 aes(x = groups, y = value, label = featureName)) +
      geom_boxplot( aes(fill = get(groupColor))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = guide_legend(title = "condition")) +
      xlab("") +
      ylab(x_lab) +
      ggtitle(title)

    return(p)

  }
)

#' @rdname runDiffAnalysis
#' @name plotBoxplotDE
#' @aliases plotBoxplotDE,RflomicsMAE-method
#' @exportMethod plotBoxplotDE
setMethod(
  f          = "plotBoxplotDE",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        featureName = NULL,
                        groupColor = "groups",
                        raw = FALSE){

    plotBoxplotDE(object = object[[SE.name]],
                  featureName = featureName,
                  groupColor = groupColor, raw = raw)

  })




##==== ACCESSEUR: getteur and setteur ====

### ---- Get DE matrix from DiffExpAnalysis ----

#' @rdname runDiffAnalysis
#' @name getDEMatrix
#' @aliases getDEMatrix,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getDEMatrix: return a matrix of experimental design.}
#' @exportMethod getDEMatrix
#' @examples
#' # See runDiffAnalysis for an example that includes getDEMatrix
setMethod(
  f          = "getDEMatrix",
  signature  = "RflomicsSE",
  definition = function(object, analysisName = "all"){

    DiffRes <- getAnalysis(object,
                           name = "DiffExpAnal",
                           subName = analysisName)$results

    if (!is.null(DiffRes$mergeDEF))
      return(DiffRes$mergeDEF)

    else
      return(NULL)

  })

#' @rdname runDiffAnalysis
#' @name getDEMatrix
#' @aliases getDEMatrix,RflomicsMAE-method
#' @exportMethod getDEMatrix
setMethod(
  f          = "getDEMatrix",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name, analysisName = "all"){

    if(!omicName %in% names(object))
      stop("This data name, ", omicName,
           ", does not exist in the your object")

    getDEMatrix(object = object[[SE.name]], analysisName = analysisName)
  })

### ---- Get union or intersection from list of contrasts ----

# very similar to filter_DE_from_SE but returns a vector instead of a SE.

#' @rdname runDiffAnalysis
#' @name getDEList
#' @aliases getDEList,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getDEList: return a vector of union or intersection of differential
#'    expressed features from list of contrasts.}
#' @param contrasts Vector of characters, expect to be contrast names.
#' Default is null, the operation (union) is performed
#' on every contrasts found.
#' @param operation character. Either union or intersection.
#' Defines the operation to perform on the DE lists from the contrasts.
#' @exportMethod getDEList
#' @importFrom tidyselect starts_with any_of
#' @examples
#' # See runDiffAnalysis for an example that includes getDEList
setMethod(
  f          = "getDEList",
  signature  = "RflomicsSE",
  definition = function(object, 
                        contrasts = NULL, 
                        operation = "union", 
                        analysisName = "all"){

    validContrasts <- getValidContrasts(object)[["contrastName"]]

    if (is.null(validContrasts) || length(validContrasts) == 0){
      validContrasts <- names(getDEMatrix(object, analysisName = analysisName))[-1]

      if (is.null(validContrasts) || length(validContrasts) == 0)
        return(NULL)
    }

    if (is.null(contrasts) || length(contrasts) == 0){
      contrasts <- validContrasts
    }
    else{
      contrasts <- intersect(contrasts, validContrasts)
      if (is.null(contrasts) || length(contrasts) == 0)
        return(NULL)
    }

    if(is.null(getDEMatrix(object, analysisName = analysisName))) return(NULL)

    df_DE <- getDEMatrix(object, analysisName = analysisName) |>
      select(c("DEF", any_of(contrasts)))

    # if (is.null(df_DE) || nrow(df_DE) == 0 || ncol(df_DE) < 2)
    #   stop("")

    if (operation == "intersection") {
      DETab <- df_DE %>%
        mutate(SUMCOL = select(., -DEF) %>% rowSums(na.rm = TRUE)) %>%
        filter(SUMCOL == length(contrasts))

    } else {
      DETab <- df_DE %>%
        mutate(SUMCOL = select(., -DEF) %>% rowSums(na.rm = TRUE)) %>%
        filter(SUMCOL >= 1)
    }

    return(unique(DETab$DEF))
  })

#' @rdname runDiffAnalysis
#' @name getDEList
#' @aliases getDEList,RflomicsMAE-method
#' @exportMethod getDEList
setMethod(
  f          = "getDEList",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        contrasts = NULL,
                        operation = "union",
                        analysisName = "all"){

    if(!SE.name %in% names(object))
      stop("This data name, ", SE.name,
           ", does not exist in the object")

    getDEList(object = object[[SE.name]],
              contrasts = contrasts,
              operation = operation,
              analysisName = analysisName)
  })

### ---- Get diff setting ----
#' @rdname runDiffAnalysis
#' @name getDiffSettings
#' @aliases getDiffSettings,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getDiffSettings: return a list of differential expression analysis
#'    settings of a given omics dataset}
#' @exportMethod getDiffSettings
#' @examples
#' # See runDiffAnalysis for an example that includes getDiffSettings
setMethod(
  f          = "getDiffSettings",
  signature  = "RflomicsSE",
  definition = function(object, analysisName = "all"){
    DiffExpAnal <- 
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)
    return(DiffExpAnal[["settings"]])
  })

#' @rdname runDiffAnalysis
#' @name getDiffSettings
#' @aliases getDiffSettings,RflomicsMAE-method
#' @exportMethod getDiffSettings
setMethod(
  f          = "getDiffSettings",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name, analysisName = "all"){
    getDiffSettings(object = object[[SE.name]], analysisName = analysisName)
  })


### ---- get Valid Contrasts : (after differential analysis) ----

#' @rdname runDiffAnalysis
#' @name getValidContrasts
#' @aliases getValidContrasts,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getValidContrasts: return a data.frame of validated contrasts}
#' @exportMethod getValidContrasts
setMethod(
  f          = "getValidContrasts",
  signature  = "RflomicsSE",
  definition = function(object, analysisName = "all"){

    DiffExpAnal <- 
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)
    
    DiffExpAnal[["results"]][["Validcontrasts"]]
  })

#' @rdname runDiffAnalysis
#' @name getValidContrasts
#' @aliases getValidContrasts,RflomicsMAE-method
#' @param omicName a dataset name
#' @exportMethod getValidContrasts
setMethod(
  f          = "getValidContrasts",
  signature  = "RflomicsMAE",
  definition = function(object, omicName, analysisName = "all"){

    if(!omicName %in% names(object))
      stop("This data name, ", omicName,
           ", does not exist in the object")

    getValidContrasts(object[[omicName]], analysisName = analysisName)

  })


### ----- Get Summary for diffExpAnalysis : -----


#' @rdname runDiffAnalysis
#' @name getDiffStat
#' @aliases getDiffStat,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getDiffStat: Get summary table from diffExpAnalysis analysis}
#' @exportMethod getDiffStat
setMethod(
  f          = "getDiffStat",
  signature  = "RflomicsSE",
  definition = function(object, analysisName = "all"){

    DiffExpAnal <-
      getAnalysis(object, name = "DiffExpAnal", subName = analysisName)

    if(length(DiffExpAnal) == 0) return(NULL)

    if(is.null(DiffExpAnal[["results"]][["stats"]])) return(NULL)

    return(DiffExpAnal[["results"]][["stats"]])
  })

#' @rdname runDiffAnalysis
#' @name getDiffStat
#' @aliases getDiffStat,RflomicsMAE-method
#' @exportMethod getDiffStat
setMethod(
  f          = "getDiffStat",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name = NULL, 
                        analysisName = "all"){

    if(!SE.name %in% names(object))
      stop("This data name, ", SE.name,
           ", does not exist in the object")

      getDiffStat(object[[SE.name]], analysisName = analysisName)
  })

### ---- getDiffAnalysesSummary ----
#' @name getDiffAnalysesSummary
#' @aliases getDiffAnalysesSummary,RflomicsMAE-method
#' @section Accessors:
#' \itemize{
#'    \item getDiffAnalysesSummary: ...}
#' @param plot FALSE or TRUE
#' @param ylabelLength max length of the labels (characters)
#' @param nbMaxLabel number of labels to print
#' @param interface Boolean. Is this plot for the interface or commandline?
#' @param analysisNames differential analysis names to plot
#' @importFrom reshape2 melt
#' @importFrom stringr str_sub str_length
#' @exportMethod getDiffAnalysesSummary
#' @return a data.frame with differential analyses summary
#' @rdname runDiffAnalysis
setMethod(
  f          = "getDiffAnalysesSummary",
  signature  = "RflomicsMAE",
  definition = function(object, 
                        analysisNames = "all", 
                        plot = FALSE,
                        ylabelLength = 30,
                        nbMaxLabel = 20,
                        interface = FALSE){

    # DataProcessing
    omicNames <- 
      getAnalyzedDatasetNames(object, analyses = "DiffExpAnal")

    df.list <- list()
    for (dataset in omicNames) {
      for(analysisName in analysisNames){
        DiffExpAnal <- 
          getAnalysis(object[[dataset]], name = "DiffExpAnal", subName = analysisName)
        
        diffSettings <- 
          getDiffSettings(
            object[[dataset]], 
            analysisName = analysisName
          )
  
        #Validcontrasts <- getValidContrasts(object[[dataset]])
        if (is.null(diffSettings[["validContrasts"]]) || 
            length(diffSettings[["validContrasts"]]) == 0)
          Validcontrasts <- diffSettings[["Contrasts.Sel"]]
        
        ValidcontrastsNames <- Validcontrasts[,c("contrastName", "tag")]
        
        p.adj <- diffSettings[["p.adj.cutoff"]]
        logFC <- diffSettings[["abs.logFC.cutoff"]]
  
        df.list[[dataset]][[analysisName]] <-
          as.data.frame(DiffExpAnal[["results"]][["stats"]])[ValidcontrastsNames$contrastName,] %>%
          mutate(dataset = dataset, 
                 contrastNames = ValidcontrastsNames$contrastName, 
                 contrastTag   = ValidcontrastsNames$tag, 
                 analysisName  = analysisName, 
                 settings = paste0("(p.adj: ", p.adj, " & log2FC: ", logFC, ")"))
      }
    }

    if (length(df.list) == 0)
      return(NULL)


    df.list.flat <- unlist(df.list, recursive = FALSE)
    df <- Reduce(rbind, df.list.flat) |>
      melt(id.vars = c("dataset", "settings", "contrastNames", "contrastTag", 
                       "All", "analysisName"),
           value.name = "Up_Down") |>
      filter(All != 0, !is.na(All)) |>
      mutate(percent = Up_Down / All * 100)

    if (isFALSE(plot))
      return(df)

    # For the interface, the labels (contrastNames) were replaced
    # by tags if they exceed nbMaxLabel to simplify the figures.
    if (isTRUE(interface) && length(unique(df$contrasts)) > nbMaxLabel) {
      
      df$tabel <- df$contrastTag
    }else{

      df$tabel <-
          unlist(lapply(df$contrastNames, function(x){
          vec <- seq(1, nchar(x), ylabelLength)
          paste(substring(x, vec, vec + ylabelLength - 1) , collapse = "\n")
        }))
    }
    
    df$tabel <- factor(df$tabel, levels = unique(df$tabel))

    
    p <- ggplot(data = df, aes(y = tabel, x = percent, fill = variable)) +
      geom_col() +
      geom_text(aes(label = Up_Down), position = position_stack(vjust = 0.5)) +
      scale_x_continuous(breaks = seq(0, 100, 25),
                         labels = paste0(seq(0, 100, 25), "%")) +
      labs(fill = NULL, x = "", y="")
    
    if(length(unique(df$dataset)) == 1 && length(unique(df$analysisName)) > 1){
      
      p <- p + facet_wrap(~ analysisName, ncol = 4)
    }
    else{
      p <- p + facet_wrap(analysisName~ dataset+settings, ncol = 4)
    }

    return(p)
  }
)
