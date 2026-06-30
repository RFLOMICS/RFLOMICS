### ============================================================================
### [03_data_processing] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# D. Charif,
# N. Bessoltane,
# A. Hulot

##==== DATA PROCESSING METHOD ====

###==== runDataProcessing ====

#' @title Omic Data Exploratory and Preprocessing
#' @rdname runDataProcessing
#' @name runDataProcessing
#' @aliases runDataProcessing,RflomicsSE-method
#' @description
#' These functions applied a data processing (filtering, normalization
#' and/or transformation, PCA) on RNAseq, proteomics, or metabolomics data.
#'
#' runDataProcessing() calls the following functions:
#' @param object An object of class \link{RflomicsSE} or class \link{RflomicsSE}
#' @param samples samples to keep.
#' @param lowCountFilter low count filtering list of arguments.
#' @param missingValueFilter missed value filtering list of arguments
#' @param transform transformation list of arguments.
#' @param normalize normalization list of arguments.
#' @param impute imputation list of arguments.
#' @param ... additional arguments
#' @return An object of class \link{RflomicsSE} or class \link{RflomicsSE}
#' @exportMethod runDataProcessing
#' @seealso
#' \link{RflomicsMAE-class}
#' \link{RflomicsSE-class}
#' \link{getTransSettings}
#' \link{getFilterSettings}
#' \link{getFilteredFeatures}
#' \link{getCoeffNorm}
#' \link{getNormSettings}
#' \link{plotLibrarySize}
#' \link{plotDataDistribution}
#' \link{plotOmicsPCA}
#' @section Accessors:
#' @section Plots:
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al.
#' DiCoExpress: a tool to process multifactorial
#' RNAseq experiments from quality controls to co-expression analysis through
#' differential analysis based on contrasts inside GLM models.
#' Plant Methods 16, 68 (2020).
#'
#' @example inst/examples/runDataProcessing.R
setMethod(
  f          = "runDataProcessing",
  signature  = "RflomicsSE",
  definition = function(object,
                        samples                 = NULL,
                        MVencoding              = "NA",
                        lowCountFilter = 
                          list(method           = NULL,
                               strategy         = NULL,
                               cpmCutoff        = NULL),
                        missingValueFilter = 
                          list(method           = NULL,
                               proportion       = NULL,
                               nbCondition      = NULL),
                        transform = list(method = NULL),
                        normalize = list(method = NULL),
                        impute = 
                          list(method           = NULL,
                               factor           = NULL),
                        ...
                        ){
    
    object <- initRawRflomicsSE(object)
    done   <- NULL
    
    # keep selected samples
    if(!is.null(samples)){
      message("[RFLOMICS] #    => select samples... ", 
              getDatasetNames(object))
      object <- runSampleFiltering(object, samples)
      object <- updateModelFormula(object)
      done <- TRUE
    }
    
    # Remplacer les 0 par NA dans la matrice assay
    if(getOmicsTypes(object) %in% c("proteomics", "metabolomics") &&
       MVencoding == "0"){
      
      assayRaw <- assay(object)
      
      if(any(is.na(assayRaw)))
        stop("Contrary to what is indicated, the data contain missing values 
             encoded as NA.")
      
      
      assayRaw[assayRaw == 0] <- NA
      assay(object) <- assayRaw
    }
    

    # feature filtering
    if(!is.null(lowCountFilter$method) || !is.null(missingValueFilter$method)){
      message("[RFLOMICS] #    => feature filtering... ", 
              getDatasetNames(object))
      object <- runFeatureFiltering(object,
                                    lowCountFilter     = lowCountFilter,
                                    missingValueFilter = missingValueFilter)
    }
    
    # Run transformation...
    if(!is.null(transform$method)){
      if(getOmicsTypes(object) %in% c("proteomics", "metabolomics")){
        message("[RFLOMICS] #    => Data transformation (", 
                transform$method, ")... ", getDatasetNames(object))
        object <- runTransformData(object, method = transform$method)
        done <- TRUE
      }
    }


    # Run Normalisation
    if(!is.null(normalize$method)){
      message("[RFLOMICS] #    => Data normalization (", 
              normalize$method, ")... ", getDatasetNames(object))
      object <- runNormalization(object, method = normalize$method)
      done <- TRUE
    }
    
    # Run Imutation
    if(!is.null(normalize$method) && 
       getOmicsTypes(object) %in% c("proteomics", "metabolomics")){
      
      message("[RFLOMICS] #    => Data imputation (", 
              impute$method, ")... ", getDatasetNames(object))
      object <- runMVImputation(object, 
                                method = impute$method,
                                factor = impute$factor)
      done <- TRUE
    }
    
    # Run PCA for filtered & normalized data
    message("[RFLOMICS] #    => Computing PCA... ", 
            getDatasetNames(object))
    object <- runOmicsPCA(object, ...)

    # tag
    object <-
      setElementToMetadata(object,
                           name    = "DataProcessing",
                           subName = "done",
                           content =  done)

    # initiate analysis results
    for(anal in c("DiffExpAnal", "DiffExpEnrichAnal", "CoExpAnal", "CoExpEnrichAnal")){
      object <-
        setElementToMetadata(object,
                             name    = anal,
                             content =  list())
    }

    return(object)
  })

#' @rdname runDataProcessing
#' @name runDataProcessing
#' @aliases runDataProcessing,RflomicsMAE-method
#' @param SE.name SE.name the name of the dataset if the input object
#' is a \link{RflomicsMAE-class}
#' @exportMethod runDataProcessing
setMethod(
  f          = "runDataProcessing",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        samples                 = NULL,
                        MVencoding              = "NA",
                        lowCountFilter = 
                          list(method           = NULL,
                               strategy         = NULL,
                               cpmCutoff        = NULL),
                        missingValueFilter = 
                          list(method           = NULL,
                               proportion       = NULL,
                               nbCondition      = NULL),
                        transform = list(method = NULL),
                        normalize = list(method = NULL),
                        impute = 
                          list(method           = NULL,
                               factor           = NULL)
                        ){

    if (!SE.name %in% names(object))
      stop("SE name must be part of this list of names: ",
           getDatasetNames(object))

    SE.processed <-  
      runDataProcessing(object = object[[SE.name]],
                        samples            = samples,
                        lowCountFilter     = lowCountFilter,
                        missingValueFilter = missingValueFilter,
                        transform          = transform,
                        normalize          = normalize,
                        impute             = impute)

    object[[SE.name]] <- SE.processed

    object <-
      setElementToMetadata(object,
                           name    = "IntegrationAnalysis",
                           content =  list())
    
    return(object)
  })

###==== runSampleFiltering ====
# filtering per sample

#' @rdname runDataProcessing
#' @name runSampleFiltering
#' @aliases runSampleFiltering,RflomicsSE-method
#' @description
#' \itemize{
#' \item runSampleFiltering:
#'   This function applied sample filtering on an dataset.
#' }
#' @param samples list of sample names to keep (is NULL, do not change anything.)
#' @exportMethod runSampleFiltering
setMethod(
  f          = "runSampleFiltering",
  signature  = "RflomicsSE",
  definition = function(object, samples = NULL) {

    if(is.null(samples)) return(object)

    # check for samples overlap
    if(any(!samples %in% colnames(object)))
      stop("Some sample names are not part of the colnames of the object")

    # new design matrix
    object2 <- object[, samples]
    if(nrow(getDesignMat(object2)) == 0) stop("no samples in object!")

    # update colData after removing samples
    # and check if this removal affects the statistical model
    object2 <- .updateColData(object2)

    # check completness
    check.res <- checkExpDesignCompleteness(object2)
    if(check.res$error) stop(check.res$messages)
    message("[RFLOMICS] #       ", check.res$messages)
    
    return(object2)
  })


#' @rdname runDataProcessing
#' @aliases runSampleFiltering,RflomicsMAE-method
#' @exportMethod runSampleFiltering
setMethod(f          = "runSampleFiltering",
          signature  = "RflomicsMAE",
          definition = function(object,
                                SE.name,
                                samples=NULL) {

            object[[SE.name]] <-
              runSampleFiltering(object[[SE.name]], samples = samples)

            return(object)
          })


###==== runFeatureFiltering ====

# METHOD to filter data

#' @name runFeatureFiltering
#' @rdname runDataProcessing
#' @aliases runFeatureFiltering,RflomicsSE-method
#' 
#' @description
#' 
#' \itemize{
#'   \item \strong{runTransformData}: Performs feature filtering depending on 
#'   the omics data type. For \strong{RNA-seq data}, lowly expressed transcripts 
#'   are removed based on count statistics. For \strong{proteomics and metabolomics 
#'   data}, missing value filtering is applied.
#' }
#'
#' @param lowCountFilter A list specifying parameters for \strong{RNA-seq} 
#' low-count filtering.
#' \code{method} defines the filtering approach and supports two options:
#' \describe{
#'   \item{\code{"filterByExpr"} (default)}{
#'     Uses \code{edgeR::filterByExpr()}. In this case, \strong{strategy} must be
#'     \code{"groups"} (default and only supported value), and \strong{cpmCutoff} is not used
#'     (should be \code{NULL}).
#'   }
#'   \item{\code{"CPM"}}{
#'     Applies a CPM-based filtering approach. The \strong{strategy} parameter controls
#'     how features are retained: \code{"NbReplicates"} (default) keeps features expressed
#'     in at least a minimum number of replicates, while \code{"NbConditions"} keeps features
#'     expressed in a minimum number of conditions. The \strong{cpmCutoff} defines the CPM
#'     threshold used to determine expressed features (default: \code{1}).
#'   }
#'   \item{\code{"none"}}{
#'     ...
#'   }
#' }
#' @param missingValueFilter A list specifying parameters for missing value handling
#' in  \strong{proteomics and metabolomics data}.
#'
#' \strong{MVencoding} defines how missing values are encoded: \code{"NA"} or \code{"0"}.
#' If \code{"0"} is used, zero values are converted to \code{NA} in the processed data.
#'
#' \strong{method} defines the filtering approach and supports two options:
#' \describe{
#'   \item{\code{"GlobalFiltering"}}{
#'     Filters features based on the overall proportion of missing values.
#'     The \strong{proportion} parameter defines the maximum allowed proportion
#'     of missing values across all samples (default: \code{0.5}).
#'   }
#'   \item{\code{"ConditionFiltering"}}{
#'     Filters features based on missing values within conditions.
#'     The \strong{nbCondition} parameter defines the minimum number of conditions
#'     required (default: \code{1}), and \strong{proportion} defines the maximum
#'     allowed proportion of missing values per condition (default: \code{0.7}).
#'   }
#'   \item{\code{"none"}}{
#'     ...
#'   }
#' }
#' 
#' @exportMethod runFeatureFiltering
setMethod(
  f         = "runFeatureFiltering",
  signature = "RflomicsSE",
  definition = function(
    object,
    lowCountFilter = 
      list(method      = c("filterByExpr", "CPM", "none"),
           strategy    = NULL,
           cpmCutoff   = 1),
    missingValueFilter = 
      list(method      = c("GlobalFiltering", "ConditionFiltering", "none"),
           proportion  = 0.5,
           nbCondition = 1)
  ){
    
    # apply data processing
    object <- switch(
      getOmicsTypes(object),
      "RNAseq" = {

        # Filter low abundance
        do.call(filterLowAbundance, 
                c(list(object = object), lowCountFilter))
      },
      "proteomics" = {

        do.call(filterMissingValues, 
                c(list(object = object), missingValueFilter))
      },
      "metabolomics" = {
        
        do.call(filterMissingValues, 
                c(list(object = object), missingValueFilter))
      }
    )

    return(object)
  })

#' @rdname runDataProcessing
#' @aliases runFeatureFiltering,RflomicsMAE-method
#' @name runFeatureFiltering
#' @exportMethod runFeatureFiltering
setMethod(
  f         = "runFeatureFiltering",
  signature = "RflomicsMAE",
  definition = function(
    object, SE.name,
    lowCountFilter = 
      list(method      = c("filterByExpr", "CPM", "none"),
           strategy    = NULL,
           cpmCutoff   = 1),
    missingValueFilter = 
      list(method      = c("GlobalFiltering", "ConditionFiltering", "none"),
           proportion  = 0.5,
           nbCondition = 1)){

    if (!SE.name %in% names(object))
      stop("SE name must be part of this list of names: ",
           getDatasetNames(object))

    object[[SE.name]] <- 
      runFeatureFiltering(object             = object[[SE.name]],
                          lowCountFilter     = lowCountFilter,
                          missingValueFilter = missingValueFilter)
    
    return(object)
  })

###==== filterLowAbundance ====

# METHOD to filter data

#' @name filterLowAbundance
#' @description
#' \itemize{
#' \item filterLowAbundance:  This function aims at removing transcript from
#' the count data matrix of an omicsof type "RNAseq".
#' by applying filtering criterion described in reference.
#' }
#' @param method The filtering model ("CPM", "filterByExpr")
#' @param strategy The filtering strategy
#' ("NbConditions" or "NbReplicates") if method == "CPM"
#' @param cpmCutoff The CPM cutoff if method == "CPM".
#' @details
#' filterLowAbundance(): By default, gene/transcript with 0 count
#' are removed from the data. The function then two stategies of filtering are
#' proposed: 1) filterByExpr implemented in egdeR 2) computes the count per
#' million or read (CPM) for each gene in each sample and gives by
#' genes the number of sample(s) which are over the cpmCutoff
#' (NbOfsample_over_cpm).
#' Then Two filtering strategies are proposed:
#' \itemize{
#' \item NbConditions:  keep gene if the NbOfsample_over_cpm >= NbConditions
#' \item NbReplicates:  keep gene if the NbOfsample_over_cpm >= min(NbReplicates)
#' \item filterByExpr: the default filtering method implemented
#' in the edgeR filterByExpr() function.
#' }
#' @keywords internal
#' @noRd
#' @importFrom edgeR DGEList filterByExpr cpm
#' @seealso edgeR::filterByExpr
setMethod(
  f         = "filterLowAbundance",
  signature = "RflomicsSE",
  definition = function(
    object,
    method    = c("filterByExpr", "CPM", "none"),
    strategy  = NULL,
    cpmCutoff = 1){
    
    if (getOmicsTypes(object) != "RNAseq")
      stop("Can't apply filterLowAbundance to omics types other than RNAseq.")
    
    if(is.null(method)){
      warning("Filtering method name not specified.")
      return(object)
    } 
    
    if (.isFiltered(object))
      stop("Data is already filtered!")
    
    method   <- match.arg(method)

    assayRaw <- assay(object)
    Groups    <- colData(object)

    if(method == "filterByExpr"){
      # 
      suported.strategies <- c("groups")
      
      if(is.null(strategy)) strategy <- suported.strategies[1]
      if(!strategy %in% suported.strategies){
        stop(strategy, " should be one of ", suported.strategies)
      }
      
      dge  <- DGEList(counts = assayRaw, genes = rownames(assayRaw))
      
      keep <- switch (strategy,
        "groups" = filterByExpr(dge, group = Groups[["groups"]])
      )

      settings <-
        list(
          method    = method,
          strategy  = strategy,
          cpmCutoff = NULL)
      
    }
    else if(method == "CPM"){

      suported.strategies <- c("NbReplicates", "NbConditions")
      
      if(is.null(strategy)) strategy <- suported.strategies[1]
      if(!strategy %in% suported.strategies){
        stop(strategy, " should be one of ", 
             paste(suported.strategies, collapse = ", "))
      }

      if(is.null(cpmCutoff)) cpmCutoff <- 1
      if(!is.numeric(cpmCutoff) || cpmCutoff < 0)
        stop(cpmCutoff, " must be an integer value > 1")

      # filter cpm
      filter_cpm <- switch (strategy,
        "NbConditions" = length(unique(Groups$groups)),
        "NbReplicates" = min(table(Groups$groups))
      )

      keep <- rowSums(cpm(assayRaw) >= cpmCutoff) >= filter_cpm

      settings <-
        list(
          method    = method,
          strategy  = strategy,
          cpmCutoff = cpmCutoff)
      
    }
    else if(method == "none"){
      
      keep <- row.names(assayRaw)
      
      settings <-
        list(
          method    = method,
          strategy  = NULL,
          cpmCutoff = NULL)
    }

    # features to filtered
    object  <- object[keep]

    # output
    Filtering <- list(
      setting = settings,
      filtered = TRUE
    )

    message("[RFLOMICS] #       Step: Low counts Filtering...")
    message("[RFLOMICS] #       method: ", method,
            ", strategy: ", strategy, ", cpmCutoff: ", cpmCutoff)

    object <-
      setElementToMetadata(object,
                           name    = "DataProcessing",
                           subName = "featureFiltering",
                           content =  Filtering)

    return(object)
  })


###==== filterMissingValues ====

# METHOD to filter data

#' @name filterMissingValues
#' @description
#' \itemize{
#' \item filterMissingValues: This function aims to identify proteins or 
#' metabolites with a high proportion of missing values.
#' }
#' @param method Method used for missing value filtering 
#'        ("GlobalFiltering", "ConditionFiltering", or "none").
#' @param proportion Minimum proportion of samples without missing values 
#'        (used when method == "GlobalFiltering"), or per condition with 
#'        (when method == "ConditionFiltering"),
#' @param nbCondition Minimum number of conditions with at least 
#'        `proportion` proportion of non-missing values.
#' @keywords internal
#' @noRd
setMethod(
  f         = "filterMissingValues",
  signature = "RflomicsSE",
  definition = function(
    object,
    method      = c("GlobalFiltering", "ConditionFiltering", "none"),
    proportion  = 0.5,
    nbCondition = 1){
    
    if (!getOmicsTypes(object) %in% c("metabolomics", "proteomics"))
      stop("Can't apply filterMissingValues to omics types other than metabolomics and proteomics")
    
    if (.isFiltered(object))
      stop("Data is already filtered!")
    
    if(is.null(method)){
      warning("Filtering method name not specified.")
      return(object)
    } 
    
    method <- match.arg(method)
    
    assayRaw      <- assay(object)
    groups        <- getDesignMat(object)$groups
    names(groups) <- getDesignMat(object)$samples
    
    # no NA no filtering
    if(!any(is.na(assayRaw)) & method != "none"){
      message("No missing values were detected in the dataset.")
      method <- "none"
    }
    
    # high proportion of missing values
    if(method == "GlobalFiltering"){

        if(!is.numeric(proportion) | proportion > 1 | proportion < 0) 
          stop ("proportion must be a numeric value representing a proportion 
              between 0 and 1.")
        
        
        message("[RFLOMICS] #       Step: Missing Value Filtering...")
        message("[RFLOMICS] #       method: globalFilter",
                "; proportion: ", proportion)
        
        # proportion de valeurs présentes par ligne
        na_prop <- rowSums(!is.na(assayRaw)) / ncol(assayRaw)
        
        # Features à garder
        keep <- na_prop >= proportion
        
        # Features à retirer
        assayFilt <- assayRaw[keep,]
        
        settings <- list(
          method           = method,
          proportion       = proportion,
          nbCondition      = NULL
        )
    }
    else if(method == "ConditionFiltering"){
        
        if(is.null(nbCondition)) nbCondition <- 1
        if(is.null(proportion)) proportion <- 0.7
        
        if(!is.numeric(nbCondition) | nbCondition < 1 | 
           nbCondition > length(unique(getDesignMat(object2)$groups)))
          stop ("nbCondition must be an integer between 1 and number of condition : ", 
                length(unique(getDesignMat(object2)$groups)))
        
        if(!is.numeric(proportion) | proportion > 1 | proportion < 0) 
          stop ("proportion must be a numeric value representing 
              a proportion between 0 and 1.")
        
        message("[RFLOMICS] #       Step: Missing Value Filtering...")
        message("[RFLOMICS] #       method: ConditionFiltering",
                ", proportion: ", proportion,
                ", nbCondition: ", nbCondition)
        
        # Filtrage conditionnel
        keep <- apply(assayRaw, 1, function(row) {
          # proportion de non-NA par condition
          prop_non_na <- tapply(!is.na(row), groups, mean)
          # garder la ligne si au moins nbCondition conditions ont >= proportion
          sum(prop_non_na >= proportion) >= nbCondition
        })
        
        # Features à retirer
        assayFilt <- assayRaw[keep,]
        
        settings <- list(
          method           = method,
          proportion       = NULL,
          nbCondition      = nbCondition
        )
    }
    else if(method == "none"){
      
      message("[RFLOMICS] #       Step: Missing Value Filtering...")
      message("[RFLOMICS] #       method: none")
      
      assayFilt <- assayRaw
      
      settings <- 
        list(
          method           = method,
          proportion       = NULL,
          nbCondition      = NULL
        )
    }
    
    # output
    Filtering <- 
      list(
        setting  = settings,
        filtered = TRUE
      )
    
    object <-
      setElementToMetadata(object,
                           name    = "DataProcessing",
                           subName = "featureFiltering",
                           content =  Filtering)
    
    assay(object) <- assayFilt
    
    return(object)
  })

###==== runTransformData ====
# METHOD to transform data

#' @rdname runDataProcessing
#' @name runTransformData
#' @aliases runTransformData,RflomicsSE-method
#' @description
#' \itemize{
#' \item runTransformData:
#'    This function applied a transformation to the dataset. The transformation
#' method is chosen according to the dataset omicstype
#' (RNAseq: none, metabolomics/proteomics: log2 or log10)
#' }
#' @param method The transformation method to store in the metadata
#' @exportMethod runTransformData
setMethod(
  f          = "runTransformData",
  signature  = "RflomicsSE",
  definition = function(object,
                        method = NULL
                        ){
    
    if (.isTransformed(object)) {
      stop("The data were already transformed beforehand! method: ",
              getTransSettings(object)$method)
      return(object)
    }
    
    object <- 
      switch (
        getOmicsTypes(object),
        "RNAseq"       = .applyTrans_readcounts(object, method = method),
        "proteomics"   = .applyTrans_intensities(object, method = method),
        "metabolomics" = .applyTrans_intensities(object, method = method),
      )

    return(object)
  })

#' @rdname runDataProcessing
#' @name runTransformData
#' @aliases runTransformData,RflomicsMAE-method
#' @exportMethod runTransformData
setMethod(
  f          = "runTransformData",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        method = c("log10", "log1p", "log2", "squareroot", "none")
  ){

    object[[SE.name]] <-
      runTransformData(object[[SE.name]], method = method)

    return(object)

  })

###==== runNormalization ====
# METHOD to normalize data
# Function non generique pour les autres data

#' @rdname runDataProcessing
#' @name runNormalization
#' @aliases runNormalization,RflomicsSE-method
#' @description
#' \itemize{
#' \item runNormalization:
#'  This function applied a normalization on a dataset.
#' The normalization method is chosen according to the dataset omics type
#' (RNAseq: TMM, metabolomics/proteomics: median)
#' }
#' @param method Normalization method. Accepted values: TMM for RNAseq, and
#' median, totalSum, or none for proteomics and metabolomics data.
#' Default values: TMM for RNAseq data and median for proteomics and metabolomics
#' data
#' @return An object of class \link{RflomicsSE}
#' The applied normalization method and computed scaling factors
#' (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "method" and
#' "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a
#' \link{RflomicsSE} object.
#' @exportMethod runNormalization
setMethod(
  f          = "runNormalization",
  signature  = "RflomicsSE",
  definition = function(object,
                        method = NULL){
    
    if (.isNormalized(object)) {
      stop("The data were already normalized beforehand. Method: ",
              getNormSettings(object)$method)
      return(object)
    }
    
    # RNA-seq
    object <- 
      switch (
        getOmicsTypes(object),
        "RNAseq"       = .applyNorm_readcounts(object, method = method),
        "proteomics"   = .applyNorm_intensities(object, method = method),
        "metabolomics" = .applyNorm_intensities(object, method = method),
      )

    return(object)
  })

#' @rdname runDataProcessing
#' @name runNormalization
#' @aliases runNormalization,RflomicsMAE-method
#' @exportMethod runNormalization
setMethod(
  f          = "runNormalization",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        method = NULL){

    object[[SE.name]] <-
      runNormalization(object = object[[SE.name]],
                       method = method)
    return(object)
  })


###==== runMVImputation ====

# METHOD to filter data

#' @rdname runDataProcessing
#' @name runMVImputation
#' @aliases runMVImputation,RflomicsSE-method
#' @description
#' \itemize{
#' \item runMVImputation: Missing value imputation approach, applied to
#' proteomics and metabolomics data, replaces missing values (0 or NA) with
#' the minimum value among all non-zero values. Additionally, variables with
#' at least one condition group without any missing values are retained without
#' further filtering.
#' }
#' @param method The imputation method ("minFeatureValue") for proteomics and
#' metabolomics data.
#' @param factor factor
#' @exportMethod runMVImputation
setMethod(
  f         = "runMVImputation",
  signature = "RflomicsSE",
  definition = function(object, 
                        method = c("minFeatureValue", "none"), 
                        factor = 1.8){
    
    if(is.null(method)) return(object)
    
    method <- match.arg(method)
    
    if(!getOmicsTypes(object) %in% c("proteomics", "metabolomics"))
      stop("Can't apply data imputation on RNAseq data.")
    
    if(is.null(factor)) factor <- 1.8
    
    omics.df  <- assay(object)
    if(method == "minFeatureValue"){
        
        # -1.8 en log2 → /3.5 sans log2
        minVals <- min(omics.df, na.rm = TRUE) - 1.8
        
        omics.df[is.na(omics.df)] <- minVals
        
        # imputation
        Imputation <- 
          list(
            setting = list(method = method,
                           factor = factor),
            results = list(minVals = minVals, 
                           stat = colSums(is.na(omics.df))),
            imputed = TRUE)
    }
    else if(method == "none"){
        Imputation <- 
          list(
            setting = list(method = method,
                           factor = NULL),
            results = list(minVals = NULL, 
                           stat = colSums(is.na(omics.df))),
            imputed = TRUE)
    }
    
    message("[RFLOMICS] #       method: ", method)
    
    assay(object) <- omics.df
    object <- setElementToMetadata(object,
                                   name    = "DataProcessing",
                                   subName = "Imputation",
                                   content = Imputation)
    return(object)
  })

#' @rdname runDataProcessing
#' @name runMVImputation
#' @aliases runMVImputation,RflomicsMAE-method
#' @exportMethod runMVImputation
setMethod(
  f          = "runMVImputation",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        method = "minFeatureValue", 
                        factor = 1.8){
    
    object[[SE.name]] <-
      runMVImputation(object  = object[[SE.name]],
                      method  = method, 
                      factor  = factor)
    return(object)
  })

### ==== runOmicsPCA ====

#' @title runOmicsPCA
#' @name runOmicsPCA
#' @aliases runOmicsPCA,RflomicsSE-method
#' @description
#' \itemize{
#' \item runOmicsPCA:
#'  This function performs a principal component analysis on omic
#' data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a
#' "Normalization" slot is present in the metadata slot, then data are
#' normalized before running the PCA according to the indicated transform
#' method.
#' }
#' This function performs a principal component analysis on omic
#' data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a
#' "Normalization" slot is present in the metadata slot, then data are
#' normalized before running the PCA according to the indicated transform
#' method.
#' @param object An object of class \link{RflomicsSE-class}.
#' @param ncomp Number of components to compute. Default is 5.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runOmicsPCA
#' @importFrom FactoMineR PCA
#' @rdname runDataProcessing
#'
setMethod(
  f          = "runOmicsPCA",
  signature  = "RflomicsSE",
  definition = function(object, ncomp = 5) {
    
    pseudo  <- assay(object)
    
    if(getOmicsTypes(object) == "RNAseq"){
      
      pseudo  <- log2(pseudo + 1)
    }
    else if(getOmicsTypes(object) %in% c("proteomics", "metabolomics")){
      
      pseudo[is.na(pseudo)] <- 0
    }
    
    PCA.res <- PCA(t(pseudo), ncp = ncomp, graph = FALSE)

    object <- setElementToMetadata(object, name = "PCA", content = PCA.res)

    return(object)
  })

#' @rdname runDataProcessing
#' @aliases runOmicsPCA,RflomicsMAE-method
#' @name runOmicsPCA
#' @title runOmicsPCA
#' @param SE.name the name of the data the normalization have to be applied to.
#' @exportMethod runOmicsPCA
setMethod(f          = "runOmicsPCA",
          signature  = "RflomicsMAE",
          definition = function(object,
                                SE.name,
                                ncomp = 5) {

            object[[SE.name]] <- runOmicsPCA(object[[SE.name]], ncomp = ncomp)
            return(object)
          })

### ==== splitRflomicsSE ====
#' @name splitRflomicsSE
#' @description
#' \itemize{
#'    \item splitRflomicsSE...}
#' @param name description
#' @keywords internal
#' @noRd
setMethod(
  f         = "splitRflomicsSE",
  signature = "RflomicsSE",
  definition <- function(object, selectedModality = NULL){
    
    Target     <- getDesignMat(object)
    BioFactors <- getBioFactors(object)
    BioFactor  <- 
      BioFactors[unlist(lapply(BioFactors, function(x){grepl(x, selectedModality)}))]
    
    Modalities <- getFactorModalities(object, factorName = BioFactor)
    Modalitie  <- 
      Modalities[unlist(lapply(Modalities, function(x){grepl(x, selectedModality)}))]
    
    Ssamples   <- Target[Target[[BioFactor]] == Modalitie,]$samples
    
    object.f <- object[, Ssamples]
    object.f@colData[[BioFactor]] <- NULL
    for(i in colnames(object.f@colData)){
     object.f@colData[[i]] <- 
       factor(object.f@colData[[i]], 
              levels = unique(object.f@colData[[i]]))
    }
    
    object.f@metadata$DataProcessing$selectedSamples <- as.vector(Ssamples)
    object.f@metadata$design$factorType <- 
      object.f@metadata$design$factorType[names(object.f@metadata$design$factorType) != BioFactor]
    
    object.f@metadata$DataProcessing$Normalization$results$coefNorm <-
      object.f@metadata$DataProcessing$Normalization$results$coefNorm[as.vector(Ssamples),]
    
    return(object.f)
  })

## ---- checkExpDesignCompleteness ----

#' @name checkExpDesignCompleteness
#' @aliases checkExpDesignCompleteness,RflomicsSE-method
#' @rdname runDataProcessing
#' @description
#' \itemize{
#'    \item checkExpDesignCompleteness: return a string with message.
#'    This method checks some experimental design characteristics.
#'    A complete design (all combinations of factor modalities with at
#'    least 2 replicates for each have to be present) with
#'    at least one biological and one batch factors are required to use the
#'    RFLOMICS workflow.}
#' @exportMethod checkExpDesignCompleteness
#' @param sampleList list of samples to check.
#' @param raw booleen.
#' @exportMethod checkExpDesignCompleteness
setMethod(
  f         = "checkExpDesignCompleteness",
  signature = "RflomicsSE",
  definition <- function(object,
                         sampleList = NULL){

    if(!is.null(sampleList)) object <- object[,sampleList]

    output <- list()
    output[["error"]] <- FALSE

    # Only works with bio and batch factors for the rest of the function
    ExpDesign <- getDesignMat(object)
    bio.fact  <- getBioFactors(object)

    # check presence of bio factors
    if (!length(getBioFactors(object)) %in% seq_len(3)){
      output[["messages"]] <-
        "Error: You need at least 1 biological factor with at least 2 levels"
      output[["error"]] <- TRUE
      return(output)
    }
    # # check presence of bash factors
    # if (!length(getBatchFactors(object)) %in% c(1,2)){
    #   output[["messages"]] <-
    #     "Error: You need at least 1 batch factor with at least 2 replicates."
    #   output[["error"]] <- TRUE
    #   return(output)
    # }

    #remplacer le code ci-dessus par celui en bas
    group_count <- .countSamplesPerCondition(ExpDesign, bio.fact)

    # check presence of relicat / batch
    # check if design is complete
    # check if design is balanced
    # check nbr of replicates
    if(min(group_count$Count) == 0){

      output[["messages"]] <- "Error: The experimental design is not complete."
      output[["error"]]    <- TRUE
    }
    else if(min(group_count$Count) == 1){

      output[["messages"]] <-  "Error: You need at least 2 biological replicates."
      output[["error"]]    <- TRUE
    }
    else if(length(unique(group_count$Count)) != 1){

      output[["messages"]] <-
        "The experimental design is complete but not balanced."
    }
    else{
      output[["messages"]] <-
        "The experimental design is complete and balanced."
    }

    return(output)
  })

#' @rdname runDataProcessing
#' @name checkExpDesignCompleteness
#' @aliases checkExpDesignCompleteness,RflomicsMAE-method
#' @param omicName the name of the data the normalization have to be applied to.
#' @exportMethod checkExpDesignCompleteness
setMethod(f         = "checkExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){

            if (is.null(omicName)) stop("Argument omicName cannot be NULL.")

            SEObject <- getRflomicsSE(object, omicName)

            checkExpDesignCompleteness(SEObject, sampleList = sampleList)
          })

##==== ACCESSORS ====

###==== getTransSettings ====

# Get transformation parameters
#' @rdname runDataProcessing
#' @name getTransSettings
#' @aliases getTransSettings,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getTransSettings: return a list of transformation settings
#'    of a given omics dataset}
#' @exportMethod getTransSettings
#' @examples
#' # See runDataProcessing for an example that includes getTransSettings
setMethod(f          = "getTransSettings",
          signature  = "RflomicsSE",
          definition = function(object){

            getAnalysis(object,
                        name = "DataProcessing",
                        subName = "Transformation")$setting
          })

#' @rdname runDataProcessing
#' @name getTransSettings
#' @aliases getTransSettings,RflomicsMAE-method
#' @exportMethod getTransSettings
setMethod(f          = "getTransSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getTransSettings(object = object[[SE.name]])
          })

###==== getFilterSettings ====

# Get filtering parameters
#' @rdname runDataProcessing
#' @name getFilterSettings
#' @aliases getFilterSettings,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getFilterSettings: return a list the filtering settings of a given
#'    omics dataset}
#' @exportMethod getFilterSettings
#' @examples
#' # See runDataProcessing for an example that includes getFilterSettings
setMethod(f          = "getFilterSettings",
          signature  = "RflomicsSE",
          definition = function(object){

            getAnalysis(object,
                        name = "DataProcessing",
                        subName = "featureFiltering")$setting
          })

#' @rdname runDataProcessing
#' @name getFilterSettings
#' @aliases getFilterSettings,RflomicsMAE-method
#' @exportMethod getFilterSettings
setMethod(f          = "getFilterSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getFilterSettings(object = object[[SE.name]])
          })


###==== getImputSettings ====

# Get filtering parameters
#' @rdname runDataProcessing
#' @name getImputSettings
#' @aliases getImputSettings,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getImputSettings: return a list the imputation settings of a given
#'    omics dataset}
#' @exportMethod getImputSettings
#' @examples
#' # See runDataProcessing for an example that includes getImputSettings
setMethod(f          = "getImputSettings",
          signature  = "RflomicsSE",
          definition = function(object){
            
            getAnalysis(object,
                        name = "DataProcessing",
                        subName = "Imputation")$setting
          })

#' @rdname runDataProcessing
#' @name getImputSettings
#' @aliases getImputSettings,RflomicsMAE-method
#' @exportMethod getImputSettings
setMethod(f          = "getImputSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getImputSettings(object = object[[SE.name]])
          })

###==== getFilteredFeatures ====

# Get filtered features
#' @rdname runDataProcessing
#' @name getFilteredFeatures
#' @aliases getFilteredFeatures,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getFilteredFeatures: return a vector of filtered features of a given
#'    omics dataset}
#' @param count number of filtered features
#' @exportMethod getFilteredFeatures
#' @examples
#' # See runDataProcessing for an example that includes getFilteredFeatures
setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsSE",
          definition = function(object, count = FALSE){

            rawobject <- initRawRflomicsSE(object)
            filteredFeatures <- 
              intersect(rownames(rawobject), rownames(object))
            
            if(count) return(length(filteredFeatures))
            return(filteredFeatures)
          })

#' @rdname runDataProcessing
#' @exportMethod getFilteredFeatures
#' @name getFilteredFeatures
#' @aliases getFilteredFeatures,RflomicsMAE-method
setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsMAE",
          definition = function(object, count = FALSE,  SE.name){
            getFilteredFeatures(object = object[[SE.name]], count = count)
          })



###==== getSelectedSamples ====

# Get filtered samples
#' @rdname runDataProcessing
#' @name getSelectedSamples
#' @aliases getSelectedSamples,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getSelectedSamples: return a vector of selected samples of a given
#'    omics dataset}
#' @exportMethod getSelectedSamples
#' @examples
#' # See runDataProcessing for an example that includes getSelectedSamples
setMethod(f          = "getSelectedSamples",
          signature  = "RflomicsSE",
          definition = function(object){

            selectedSamples <- colnames(object)

            return(selectedSamples)
          })

#' @rdname runDataProcessing
#' @exportMethod getSelectedSamples
#' @name getSelectedSamples
#' @aliases getSelectedSamples,RflomicsMAE-method
setMethod(f          = "getSelectedSamples",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){

            getSelectedSamples(object[[SE.name]])
          })

###==== getCoeffNorm ====

#' @rdname runDataProcessing
#' @name getCoeffNorm
#' @aliases getCoeffNorm,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getCoeffNorm: return a named vector with normalization coefficients
#'    of a given omics dataset}
#' @exportMethod getCoeffNorm
#' @examples
#' # See runDataProcessing for an example that includes getCoeffNorm
setMethod(
  f          = "getCoeffNorm",
  signature  = "RflomicsSE",
  definition = function(object){

    metadata(object)[["DataProcessing"]][["Normalization"]][["results"]][["coefNorm"]]
  })

#' @rdname runDataProcessing
#' @name getCoeffNorm
#' @aliases getCoeffNorm,RflomicsMAE-method
#' @exportMethod getCoeffNorm

setMethod(f          = "getCoeffNorm",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){

            getCoeffNorm(object = object[[SE.name]])
          })

###==== getNormSettings ====
# Get normalizationparameters

#' @rdname runDataProcessing
#' @name getNormSettings
#' @aliases getNormSettings,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item getNormSettings: return a list of normalization settings
#'    of a given omics dataset}
#' @exportMethod getNormSettings
#' @examples
#' # See runDataProcessing for an example that includes getNormSettings
setMethod(f          = "getNormSettings",
          signature  = "RflomicsSE",
          definition = function(object){

            getAnalysis(object,
                        name = "DataProcessing",
                        subName = "Normalization")$setting
          })

#' @rdname runDataProcessing
#' @name getNormSettings
#' @aliases getNormSettings,RflomicsMAE-method
#' @exportMethod getNormSettings
setMethod(f          = "getNormSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getNormSettings(object = object[[SE.name]])
          })

##==== GRAPHICAL METHOD ====

###==== plotLibrarySize ====

#' @rdname runDataProcessing
#' @name plotLibrarySize
#' @aliases plotLibrarySize,RflomicsSE-method
#' @section Plots:
#' \itemize{
#'    \item plotLibrarySize: return barplot of library size by sample.}
#' @param raw a boolean
#' @exportMethod plotLibrarySize
#' @examples
#' # See runDataProcessing for an example that includes plotLibrarySize
setMethod(f          = "plotLibrarySize",
          signature  = "RflomicsSE",
          definition = function(object, raw = FALSE)
          {

            if (getOmicsTypes(object) != "RNAseq")
              stop("Data are not RNAseq!")


            if(raw) object <- initRawRflomicsSE(object)

            labels <- getLabs4plot(object)

            pseudo <- colSums(assay(object), na.rm = TRUE)

            Groups      <- getDesignMat(object)
            libSizeNorm <-
              full_join(Groups,
                        data.frame("value" = pseudo, "samples" = names(pseudo)),
                        by = "samples") |> arrange(groups)

            libSizeNorm$samples <-
              factor(libSizeNorm$samples, levels = levels(Groups$samples))

            p <-
              ggplot(libSizeNorm, aes(x = samples, y = value, fill = groups)) +
              geom_bar(stat = "identity" ) + 
              theme_bw() +
              theme(axis.text.x =  element_text(angle = 45, hjust = 1),
                    legend.position  = "none") +
              labs(x = "", y = "Total read count per sample") +
              ggtitle(labels$title)

            return(p)

          })

#' @rdname runDataProcessing
#' @name plotLibrarySize
#' @aliases plotLibrarySize,RflomicsMAE-method
#' @exportMethod plotLibrarySize
setMethod(f          = "plotLibrarySize",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name,
                                raw = FALSE){

            if (getOmicsTypes(object[[SE.name]]) == "RNAseq") {
              return(plotLibrarySize(object[[SE.name]],
                                     raw = raw))
            }else{
              stop("This function only applies to RNAseq data")
            }
          })

###==== plotDataDistribution ====

#' @rdname runDataProcessing
#' @name plotDataDistribution
#' @aliases plotDataDistribution,RflomicsSE-method
#' @description
#' \itemize{
#'    \item plotDataDistribution: return boxplot or density plot of expression
#'    or abundance distribution.}
#' @param plot plot type ("boxplot" or "density")
#' @param raw boolean. Plot the raw data or the transformed ones (raw = FALSE)
#' @exportMethod plotDataDistribution
#' @importFrom reshape2 melt
#' @examples
#' # See runDataProcessing for an example that includes plotDataDistribution
setMethod(
  f = "plotDataDistribution",
  signature = "RflomicsSE",
  definition = function(object, plot = "boxplot", raw = FALSE) {
    
    if(raw) object <- initRawRflomicsSE(object)
    
    pseudo <- assay(object)
    Groups <- getDesignMat(object)
    
    if(getOmicsTypes(object) == "RNAseq")
      pseudo <- log2(pseudo + 1)

    labels <- getLabs4plot(object)

    pseudo.gg <- pseudo %>% melt()
    colnames(pseudo.gg) <- c("features", "samples", "value")

    pseudo.gg <- pseudo.gg %>%
      full_join(Groups, by = "samples") %>%
      arrange(groups)

    pseudo.gg$samples <- factor(pseudo.gg$samples,
                                levels = unique(pseudo.gg$samples))
    switch(plot,
           "density" = {
             p <- ggplot(pseudo.gg) +  theme_bw() +
               geom_density( aes(x = value, group = samples, color = groups),
                             trim = FALSE) +
               xlab(labels$x_lab) +
               theme(legend.position = "none") +
               ggtitle(labels$title)
           },
           "boxplot" = {
             p <-  ggplot(pseudo.gg, aes(x = samples, y = value)) +
               geom_boxplot( aes(fill = groups), outlier.size = 0.3) +
               theme_bw() +
               theme(axis.text.x =  element_text(angle = 45, hjust = 1),
                     legend.position = "none",
                     plot.margin= margin(0.5,0.5,0.5,1,"cm")) +
               xlab("") +
               ylab(labels$x_lab) +
               ggtitle(labels$title) #+
             # geom_point(alpha = 1/100,size=0)
           }
    )

    return(p)
  }
)

#' @rdname runDataProcessing
#' @name plotDataDistribution
#' @aliases plotDataDistribution,RflomicsMAE-method
#' @exportMethod plotDataDistribution
setMethod(
  f = "plotDataDistribution",
  signature = "RflomicsMAE",
  definition = function(object, SE.name,
                        plot = "boxplot",
                        raw = FALSE) {
    plotDataDistribution(
      object = object[[SE.name]],
      plot = plot,
      raw = raw
    )
  }
)

### ---- plotOmicsPCA ----
#' @name plotOmicsPCA
#' @aliases plotOmicsPCA,RflomicsSE-method
#' @rdname runDataProcessing
#' @section Plots:
#' \itemize{
#'    \item plotOmicsPCA:
#' This function plot the factorial map from a PCA object stored
#' in a \link{RflomicsSE-class} object. By default, samples are
#' colored by groups (all combinations of level's factor)}
#' @param raw boolean. Does the pca have to be ran on raw data or transformed
#' @param axes A vector giving the two axis that have to be drawn for the
#' factorial map
#' @param groupColor All combination of level's factor
#' @importFrom FactoMineR coord.ellipse
#' @exportMethod plotOmicsPCA
#' @examples
#' # See runDataProcessing for an example that includes plotOmicsPCA
setMethod(
  f          = "plotOmicsPCA",
  signature  = "RflomicsSE",
  definition = function(object,
                        raw = TRUE,
                        axes = c(1, 2),
                        groupColor = "groups"){
    
    # define pca axis
    if (length(axes) != 2) axes <- c(1, 2)
    
    PC1 <- paste("Dim.", axes[1], sep = "")
    PC2 <- paste("Dim.", axes[2], sep = "")

    if (PC1 == PC2) PC2 <- PC1 + 1
    
    
    if(raw) object <- initRawRflomicsSE(object)
    
    # get pca score
    ExpDesign <- getDesignMat(object)
    score     <- as.data.frame(metadata(object)$PCA$ind$coord[, axes])
    score$samples <- row.names(score)
    score     <- right_join(score, ExpDesign, by = "samples")

    var1 <- round(metadata(object)$PCA$eig[axes, 2][1], digits = 3)
    var2 <- round(metadata(object)$PCA$eig[axes, 2][2], digits = 3)

    # get labels
    labels <- getLabs4plot(object)
    # plot
    p <- ggplot(score, aes(x = !!sym(PC1), y = !!sym(PC2), color = !!sym(groupColor)))  +
      geom_point(size = 2) +
      geom_text(aes(label = samples), size = 2, vjust = "inward", hjust = "inward") +
      xlab(paste(PC1, " (", var1, "%)", sep = "")) +
      ylab(paste(PC2, " (", var2, "%)", sep = "")) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      theme_bw() +
      theme(
        strip.text.x =  element_text(size = 8, face = "bold.italic"),
        strip.text.y =  element_text(size = 8, face = "bold.italic")
      ) + 
      ggtitle(labels$title)

    # ellipse corr
    aa <- select(score, all_of(groupColor), all_of(PC1), all_of(PC2))
    bb <- coord.ellipse(aa, bary = TRUE)
    p <- p + geom_polygon( data = bb$res,
      aes(x = !!sym(PC1), y = !!sym(PC2), fill = !!sym(groupColor)),
      show.legend = FALSE,
      alpha = 0.1
    )

    return(p)
  })

#' @rdname runDataProcessing
#' @name plotOmicsPCA
#' @aliases plotOmicsPCA,RflomicsMAE-method
#' @exportMethod plotOmicsPCA

setMethod(
  f          = "plotOmicsPCA",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        raw = FALSE,
                        axes = c(1, 2),
                        groupColor = "groups") {
    plotOmicsPCA(object[[SE.name]],
                 raw = raw,
                 axes = axes,
                 groupColor = groupColor)

  }
)

### ---- plotExpDesignCompleteness ----
#' @name plotExpDesignCompleteness
#' @aliases plotExpDesignCompleteness,RflomicsSE-method
#' @section Plots:
#' \itemize{
#'    \item plotExpDesignCompleteness:
#' This method checks that experimental design constraints are satisfied and
#' plot a summary of the design.
#' A complete design (all combinations of factor modalities with at least 2
#' replicates for each have to be present)
#' with at least one biological and one batch factors are required to use the
#' RFLOMICS workflow.}
#' @param sampleList list of samples to check.
#' @exportMethod plotExpDesignCompleteness
#' @rdname runDataProcessing
#' @examples
#' # See runDataProcessing for an example that includes plotExpDesignCompleteness
setMethod(f         = "plotExpDesignCompleteness",
          signature = "RflomicsSE",
          definition <- function(object, sampleList=NULL){

            # reduce object to sample list
            if(!is.null(sampleList))
              object <- object[,sampleList]

            check <- checkExpDesignCompleteness(object)

            # Only works with bio and batch factors for the rest of the
            # function
            ExpDesign <- getDesignMat(object)
            bio.fact <- getBioFactors(object)

            group_count <- .countSamplesPerCondition(ExpDesign, bio.fact)

            plot <- .plotExperimentalDesign(counts = group_count,
                                            message= check[["messages"]])
            return(plot)
          })

#' @param omicName a character string with the name of the dataset
#' @exportMethod plotExpDesignCompleteness
#' @aliases plotExpDesignCompleteness,RflomicsMAE-method
#' @rdname runDataProcessing
setMethod(f         = "plotExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){

            SEObject <- getRflomicsSE(object, omicName)

            plotExpDesignCompleteness(SEObject, sampleList = sampleList)
          })




###==== plotMissingValues ====

#' @rdname runDataProcessing
#' @name plotMissingValues
#' @aliases plotMissingValues,RflomicsSE-method
#' @description
#' \itemize{
#'    \item plotMissingValues: return barplot of % of messing values.
#' }
#' @param raw a boolean
#' @exportMethod plotMissingValues
#' @importFrom reshape2 melt
#' @examples
#' # See runDataProcessing for an example that includes plotMissingValues
setMethod(
  f          = "plotMissingValues",
  signature  = "RflomicsSE",
  definition = function(object, raw = FALSE)
  {

    if(raw) object <- initRawRflomicsSE(object)

    labels <- getLabs4plot(object)

    mat <- assay(object)

    df <-
      melt(mat) |>
      mutate(missed_value = ifelse(is.na(value), "yes (NA)", ifelse(value == 0, "yes (0)", "no"))) |>
      group_by(Var2, missed_value) |> count(name = "nb_missed") |>
      mutate(missed_p = nb_missed/nrow(mat)*100)

    p <- ggplot(df) +
      geom_col(aes(x = Var2, y = missed_p, fill = missed_value)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "", y = "% of missing values", fill = "Missed values") +
      ggtitle(labels$title)

    return(p)

  })

#' @rdname runDataProcessing
#' @name plotMissingValues
#' @aliases plotMissingValues,RflomicsMAE-method
#' @exportMethod plotMissingValues
setMethod(f          = "plotMissingValues",
          signature  = "RflomicsMAE",
          definition = function(object,
                                SE.name,
                                raw = FALSE){

            return(plotMissingValues(object[[SE.name]], raw = raw))
          })

## ---- CHECK ----

### ---- isProcessedData ----
#' @rdname runDataProcessing
#' @name isProcessedData
#' @param filter boolean. If TRUE, check if data is filtered (low counts/RNAseq)
#' @param trans boolean. If TRUE, check if data is transformed
#' @param norm boolean. If TRUE, check if data is normalized
#' @param log boolean. If TRUE, check if the data has been log-transformed
#' (RNAseq).
#' @aliases isProcessedData,RflomicsSE-method
#' @section Accessors:
#' \itemize{
#'    \item isProcessedData: return }
#' @exportMethod isProcessedData
setMethod(
  f          = "isProcessedData",
  signature  = "RflomicsSE",
  definition = function(object,
                        filter = FALSE,
                        trans = FALSE,
                        norm = FALSE,
                        log = FALSE){

    if(isFALSE(unique(filter, trans, norm))) {
      switch (
        getOmicsTypes(object),
        "RNAseq" = message(" is the ", getDatasetNames(object),
                           "data filtered and normalized"),
        message(" is the ", getDatasetNames(object),
                "data normalized and normalized")
      )
      norm = TRUE}

    if(!identical(colnames(object), getSelectedSamples(object))){
      message("The outlier samples have not been removed from the data.")
      return(FALSE)
    }

    switch (getOmicsTypes(object),
            "RNAseq" = {

              if(norm & isFALSE(.isNormalized(object))){
                message(getDatasetNames(object), " data have not been from the data.")
                return(FALSE)
              }
              if(filter & isFALSE(.isFiltered(object)))   return(FALSE)

              if(log) return(is.null(getAnalysis(object,
                                                 name = "DataProcessing",
                                                 subName = "log")))


            },
            {
              if(norm  & isFALSE(.isNormalized(object)))  return(FALSE)
              if(trans & isFALSE(.isTransformed(object))) return(FALSE)
            }
    )
    return(TRUE)
  })

#' @rdname runDataProcessing
#' @name isProcessedData
#' @aliases isProcessedData,RflomicsMAE-method
#' @exportMethod isProcessedData
setMethod(
  f          = "isProcessedData",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        filter = TRUE,
                        trans = TRUE,
                        norm = TRUE,
                        log = FALSE
  ){

    if (!SE.name %in% getDatasetNames(object)){
      stop("SE name must be part of this list of names: ",
           getDatasetNames(object))
    }

    results <-
      isProcessedData(object[[SE.name]],
                      filter = filter,
                      trans = trans,
                      norm = norm,
                      log = log)

    return(results)
  })
