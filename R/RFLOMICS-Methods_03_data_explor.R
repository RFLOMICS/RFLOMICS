### ============================================================================
### [03_data_processing] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# D. Charif, 
# N. Bessoltane, 
# A. Hulot

##==== DATA PROCESSING METHOD ====

###==== runDataProcessing ====

#' @title Data Exploratory and processing
#' @rdname runDataProcessing
#' @name runDataProcessing
#' @aliases runDataProcessing,RflomicsSE-method
#' @description
#' These functions applied a data processing (filtering, normalization 
#' and/or transformation, PCA) on RNAseq, proteomic, or metabolomic data.
#' 
#' runDataProcessing() calls the following functions:
#' @param object An object of class \link{RflomicsSE} or class \link{RflomicsSE} 
#' @param samples samples to keep.
#' @param filterStrategy strategy of RNAseq low count filtering. 
#' Mandatory for RNAseq data. Default value: "NbReplicates".
#' @param cpmCutoff CPM cutoff for RNAseq low count filtering.
#'  Mandatory for RNAseq data. Default value: 1.
#' @param normMethod of normalization. Mandatory for RNAseq data. 
#' Default value: RNAseq = TMM.
#' @param transformMethod method of transformation.
#' @param userNormMethod method used by user to normalize data.
#' @param userTransMethod method used by user to transforme data.
#' @return An object of class \link{RflomicsSE} or class \link{RflomicsSE} 
#' @exportMethod runDataProcessing
#' @seealso 
#' \link{RflomicsMAE-class} 
#' \link{RflomicsSE-class} 
#' \link{getProcessedData} 
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
                        samples=NULL, 
                        filterStrategy = NULL, 
                        cpmCutoff = NULL, 
                        transformMethod = NULL,
                        normMethod= NULL, 
                        imputMethod = NULL,
                        userTransMethod = "unknown",
                        userNormMethod = "unknown"){
    
    # apply data processing
    if(getOmicsTypes(object) == "RNAseq"){
      if(!is.null(transformMethod))
        warning(getOmicsTypes(object), " data don't need to be transformed.")
    }
    
    if(getOmicsTypes(object) != "RNAseq" &&
       (!is.null(filterStrategy) || !is.null(cpmCutoff)))
      warning(getOmicsTypes(object), " data don't need to be filtered")
    
    # keep selected samples
    if(is.null(samples)) samples <- colnames(object)
    message("[RFLOMICS] #    => select samples... ", 
            getDatasetNames(object))
    object <- runSampleFiltering(object, samples)
    
    # imputation
    message("[RFLOMICS] #    => feature filtering... ", 
            getDatasetNames(object))
    object <- runFeatureFiltering(object, 
                                  filterStrategy = filterStrategy, 
                                  cpmCutoff = cpmCutoff,
                                  imputMethod = imputMethod)
    
    # Run transformation...
    if(getOmicsTypes(object) != "RNAseq"){
      message("[RFLOMICS] #    => Data transformation... ", 
              getDatasetNames(object))
      object <- runTransformData(object, 
                                 transformMethod = transformMethod, 
                                 userTransMethod = userTransMethod) 
    }
    
    # Run Normalisation
    message("[RFLOMICS] #    => Data normalization... ", 
            getDatasetNames(object))
    object <- runNormalization(object, 
                               normMethod = normMethod,
                               userNormMethod = userNormMethod)
    
    # Run PCA for filtered & normalized data
    message("[RFLOMICS] #    => Compute PCA ", getDatasetNames(object))
    object <- runOmicsPCA(object, ncomp = 5, raw = FALSE)
    
    # tag
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing",
                           subName = "done",
                           content =  TRUE)
    
    # initiate analysis results
    object <- 
      setElementToMetadata(object, 
                           name    = "DiffExpAnal",
                           content =  list())
    object <- 
      setElementToMetadata(object, 
                           name    = "DiffExpEnrichAnal",
                           content =  list())    
    object <- 
      setElementToMetadata(object, 
                           name    = "CoExpAnal",
                           content =  list())
    object <- 
      setElementToMetadata(object, 
                           name    = "CoExpEnrichAnal",
                           content =  list())
    
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
                        samples=NULL, 
                        filterStrategy = NULL, 
                        cpmCutoff = NULL, 
                        transformMethod = NULL,
                        normMethod= NULL, 
                        userTransMethod = "unknown",
                        userNormMethod = "unknown"){
    
    if (!SE.name %in% names(object))
      stop("SE name must be part of this list of names: ",
           getDatasetNames(object))
    
    SE.processed <-  runDataProcessing(object = object[[SE.name]],
                                       samples = samples, 
                                       filterStrategy = filterStrategy, 
                                       cpmCutoff = cpmCutoff,
                                       normMethod= normMethod, 
                                       transformMethod = transformMethod,
                                       userTransMethod = userTransMethod,
                                       userNormMethod  = userNormMethod)
    
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
#' @exportMethod runSampleFiltering
setMethod(
  f          = "runSampleFiltering",
  signature  = "RflomicsSE",
  definition = function(object, samples = NULL) {
    
    if(is.null(samples)) samples <- colnames(object)
    
    # check if samples overlap with
    if(any(!samples %in% colnames(object)))
      stop("Some sample names are not part of the colnames of the object")
    
    # new design matrix
    object2 <- object[, samples]
    if(nrow(getDesignMat(object2)) == 0) stop("no samples in object!")
    
    # update colData after removing samples
    # and check if this removal affect the stat model
    object2 <- .updateColData(object2)
    
    # check completness
    check.res <- checkExpDesignCompleteness(object2)
    if(check.res$error) stop(check.res$messages)
    message("[RFLOMICS] #       ", check.res$messages)
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "selectedSamples",
                           content = samples)
    
    # initiate
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "featureFiltering", 
                           content =  list())
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "Transformation", 
                           content =  list())
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "Normalization", 
                           content =  list())
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "log", 
                           content =  NULL)
    
    object <- 
      setElementToMetadata(object, 
                           name    = "PCAlist", 
                           subName = "norm", 
                           content =  NULL)
    
    return(object)
    
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

#' @rdname runDataProcessing
#' @aliases runFeatureFiltering,RflomicsSE-method
#' @name runFeatureFiltering
#' @description 
#' \itemize{
#' \item runFeatureFiltering: This function allows filtering variables in omics 
#' data. In the case of RNA-seq data, it involves filtering out transcripts with 
#' low counts, while in the case of proteomics and metabolomics data, it applies 
#' the imputation procedure.
#' }
#' @param filterMethod The filtering model ("CPM") for RNAseq data.
#' @param filterStrategy The filtering strategy
#' ("NbConditions" or "NbReplicates") for RNAseq data.
#' @param cpmCutoff The CPM cutoff for RNAseq data.
#' @param imputMethod The imputation method ("MVI") for proteomics and 
#' metabolomics data.
#' @details
#' Low count filtering procedure: By default, transcript with 0 count 
#' are removed from the data. The function then computes the count per million 
#' or read (CPM) for each gene in each sample and gives by genes the number of 
#' sample(s) which are over the cpmCutoff (NbOfsample_over_cpm).
#' Then Two filtering strategies are proposed:
#' \itemize{
#' \item NbConditions:  keep gene if the NbOfsample_over_cpm >= NbConditions 
#' \item NbReplicates:  keep gene if the NbOfsample_over_cpm >= min(NbReplicat) 
#' \item filterByExpr: the default filtering method implemented 
#' in the edgeR filterByExpr() function. 
#' }
#' @details
#' Missing value imputation: This approach, applied to proteomics and 
#' metabolomics data, replaces missing values (0 or NA) with the minimum value 
#' among all non-zero values. Additionally, variables with at least one 
#' condition group without any missing values are retained without further filtering.
#' @exportMethod runFeatureFiltering
setMethod(
  f         = "runFeatureFiltering",
  signature = "RflomicsSE",
  definition = function(object, 
                        filterMethod   = NULL,
                        filterStrategy = NULL,
                        cpmCutoff      = NULL,
                        imputMethod    = NULL){
    
    # apply data processing
    object <-switch(
      getOmicsTypes(object),
      "RNAseq" = {
        
        if(!is.null(imputMethod))
          warning("we don't apply data imputation on ", 
                  getOmicsTypes(object), " data")
        
        # Filter low abundance
        filterLowAbundance(object         = object, 
                           filterMethod   = filterMethod,
                           filterStrategy = filterStrategy, 
                           cpmCutoff      = cpmCutoff)
      },
      {
        if(!is.null(filterMethod) || !is.null(filterStrategy) || 
           !is.null(cpmCutoff))
          warning("we don't apply low count filtering on ",
                  getOmicsTypes(object), " data.")
        
        # imputation
        dataImputation(object      = object, 
                       imputMethod = imputMethod) 
      }
    )
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "Transformation", 
                           content =  list())
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "Normalization", 
                           content =  list())
    
    object <- 
      setElementToMetadata(object, 
                           name    = "PCAlist", 
                           subName = "norm", 
                           content =  NULL)
    
    return(object)
  })

#' @rdname runDataProcessing
#' @aliases runFeatureFiltering,RflomicsMAE-method
#' @name runFeatureFiltering
#' @exportMethod runFeatureFiltering
setMethod(
  f         = "runFeatureFiltering",
  signature = "RflomicsMAE",
  definition = function(object, SE.name,
                        filterMethod = NULL,
                        filterStrategy = NULL,
                        cpmCutoff = NULL,
                        imputMethod = NULL){
    
    if (!SE.name %in% names(object))
      stop("SE name must be part of this list of names: ",
           getDatasetNames(object))
    
    SE.processed <- runFeatureFiltering(object = object[[SE.name]],
                                        filterMethod = filterMethod,
                                        filterStrategy = filterStrategy,
                                        cpmCutoff = cpmCutoff,
                                        imputMethod = imputMethod)
    
    object[[SE.name]] <- SE.processed
    
    return(object)
  })

###==== dataImputation ====

# METHOD to filter data

#' @name dataImputation
#' @description 
#' \itemize{
#' \item dataImputation: Missing value imputation approach, applied to 
#' proteomics and metabolomics data, replaces missing values (0 or NA) with 
#' the minimum value among all non-zero values. Additionally, variables with 
#' at least one condition group without any missing values are retained without 
#' further filtering.
#' }
#' @param imputMethod The imputation method ("MVI") for proteomics and 
#' metabolomics data.
#' @keywords internal
#' @noRd
setMethod(
  f         = "dataImputation",
  signature = "RflomicsSE",
  definition = function(object, imputMethod = "MVI"){
    
    if(getOmicsTypes(object) == "RNAseq")
      stop("Can't apply data imputation on RNAseq data.")
    
    if(imputMethod == "" || is.null(imputMethod))
      imputMethod <- "MVI"
    
    featureFiltering <- switch (
      imputMethod,
      "MVI" = {
        object2   <- getProcessedData(object)
        omics.df  <- assay(object2)
        design    <- getDesignMat(object2)
        min.value <- min(omics.df[omics.df != 0], na.rm = TRUE)
        
        # replace NA by 0 : done by createRflomicsMAE()
        # omics.df[is.na(omics.df)] <- 0
        # omics.df <- as.data.frame(lapply(omics.df, as.numeric))
        
        # filter-in variables with at leat 1 group with no 0
        keep.features <- Reduce(
          union, 
          lapply(unique(design$groups), function(x){ 
            
            samples.g <- design[design$groups == x,]$samples
            row.names(omics.df[Reduce(pmin, omics.df[samples.g]) > 0,])
          }))
        
        filteredFeatures <- setdiff(row.names(omics.df), keep.features)
        if(length(filteredFeatures) == 0) filteredFeatures <- NULL
        
        list(
          setting = list(method = "MVI", 
                         minValue = min.value,
                         suppInfo = "missing value imputation"),
          results = list(filteredFeatures = filteredFeatures),
          filtered = FALSE)
      },
      {
        stop("The ",imputMethod, 
             " method is not supported by rflomics for data imputation.")
      }
    )
    
    message("[RFLOMICS] #       approach: Data Imputation...")
    message("[RFLOMICS] #       method: ", imputMethod)
    
    object <- setElementToMetadata(object,
                                   name = "DataProcessing", 
                                   subName = "featureFiltering",
                                   content = featureFiltering)
    
    # replace 0 by min.value (applyFiltering)
    # omics.df[omics.df == 0] <- min.value
    
    return(object)
  })

###==== filterLowAbundance ====

# METHOD to filter data

#' @name filterLowAbundance
#' @description 
#' \itemize{
#' \item filterLowAbundance:  This function aims at removing transcript from 
#' the count data matrix of an omic of type "RNAseq".
#' by applying filtering criterion described in reference. 
#' }
#' @param filterMethod The filtering model ("CPM")
#' @param filterStrategy The filtering strategy 
#' ("NbConditions" or "NbReplicates")
#' @param cpmCutoff The CPM cutoff.
#' @details
#' filterLowAbundance(): By default, gene/transcript with 0 count 
#' are removed from the data. The function then
#' computes the count per million or read (CPM) for each gene 
#' in each sample and gives by
#' genes the number of sample(s) which are over the cpmCutoff 
#' (NbOfsample_over_cpm).
#' Then Two filtering strategies are proposed:
#' \itemize{
#' \item NbConditions:  keep gene if the NbOfsample_over_cpm >= NbConditions 
#' \item NbReplicates:  keep gene if the NbOfsample_over_cpm >= min(NbReplicat) 
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
  definition = function(object, 
                        filterMethod = "CPM", 
                        filterStrategy = "NbReplicates", 
                        cpmCutoff = 1){
    
    suported.strategies <- c("NbReplicates","NbConditions")
    
    if (getOmicsTypes(object) != "RNAseq")
      stop("Can't apply this method to omics types other than RNAseq.")
    
    if(.isFiltered(object)) 
      stop("Data is already filtered!")
    
    if (is.null(filterStrategy)) filterStrategy <- suported.strategies[1]
    if (isFALSE(filterStrategy %in% suported.strategies))
      stop("filterStrategy argument must be one of this tow options : ", 
           "NbReplicates or NbConditions")
    
    if(is.null(filterMethod)) filterMethod <- "CPM"
    if(filterMethod != "CPM")
      stop("filterMethod argument must be one of this tow options : CPM")
    
    if(is.null(cpmCutoff)) cpmCutoff <- 1
    if(!is.numeric(cpmCutoff) || cpmCutoff < 0)
      stop(cpmCutoff, " must be a integer value > 1")
    
    # filter outlier samples
    object2 <- getProcessedData(object)
    
    assayFilt  <- assay(object2)
    
    # nbr of genes with 0 count
    genes_flt0  <- object2[rowSums(assayFilt) <= 0, ]@NAMES
    
    # remove 0 count
    objectFilt  <- object2[rowSums(assayFilt)  > 0, ]
    assayFilt   <- assay(objectFilt)
    
    # filter cpm
    Groups       <- getDesignMat(object2)
    NbReplicate  <- table(Groups$groups)
    NbConditions <- length(unique(Groups$groups))
    
    # low count filtering
    keep <-
      switch(filterStrategy,
             "NbConditions" = { rowSums(cpm(assayFilt) >= cpmCutoff) >= NbConditions },
             "NbReplicates" = { rowSums(cpm(assayFilt) >= cpmCutoff) >= min(NbReplicate)},
             "filterByExpr" = 
               { 
                 dge <- DGEList(counts = assayFilt, genes = rownames(assayFilt))
                 filterByExpr(dge)
               }
      )
    
    # features to filtered
    genes_flt1  <- objectFilt[!keep]@NAMES
    
    # object <- objectFilt[keep]
    # output
    Filtering <- list(
      setting = list(method         = filterMethod,
                     filterStrategy = filterStrategy,
                     cpmCutoff      = cpmCutoff),
      results = list(filteredFeatures = c(genes_flt0, genes_flt1)),
      filtered = FALSE
    )
    
    message("[RFLOMICS] #       approach: Low counts Filtering...")
    message("[RFLOMICS] #       method: ", filterMethod, 
            ", strategy: ", filterStrategy, ", cpmCutoff: ", cpmCutoff)
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "featureFiltering", 
                           content =  Filtering)
    
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
#'    This function applied a transformation on dataset. The transformation 
#' method is chosen according to dataset omic type 
#' (RNaseq: none, metabolomics/proteomics: log2) 
#' }
#' @param transformMethod The transformation to store in the metadata
#' @exportMethod runTransformData
setMethod(
  f          = "runTransformData",
  signature  = "RflomicsSE",
  definition = function(object, 
                        transformMethod = "log2",
                        userTransMethod = "unknown"
  ){
    
    if (getOmicsTypes(object) == "RNAseq"){
      warning("It is not recommended to transform RNAseq data.")
      return(object)
    }
    
    if(.isTransformed(object)){
      warning("Data is already transformed!")
      return(object)
    }
    
    # accepted value for normMethod
    # default value : 1rst element
    default.methods <- 
      switch (getOmicsTypes(object),
              "proteomics"   = c("log2", "none"),
              "metabolomics" = c("log2", "none")
      )
    
    # check normMethod param
    if(is.null(transformMethod)) transformMethod <- default.methods[1]
    if(!transformMethod %in% default.methods)
      stop(transformMethod, 
           " is not an allowed value for the parameter transformMethod",
           " Accepted values: ", paste0(default.methods, collapse = ", "))
    
    
    if(userTransMethod == "" || is.null(userTransMethod))
      userTransMethod <- "unknown"
    
    userTransMethod <- 
      switch(transformMethod, "none" = userTransMethod, NULL)
    
    if(is.null(getFilterSettings(object)$method))
      stop("Before transforming the ", getOmicsTypes(object), " data, ", 
           "you must first run the feature filtering (missing value imputation).",
           " See ?runDataProcessing.")
    
    # output
    transformation <- 
      list(
        setting = list(method = transformMethod,
                       suppInfo = switch (transformMethod, 
                                          "none" = userTransMethod, 
                                          NULL)),
        results  = NULL,
        transformed = FALSE
      )
    
    message("[RFLOMICS] #       method: ", 
            switch(transformMethod, 
                   "none" = paste0("already transformed (", userTransMethod, ")"), 
                   transformMethod))
    
    object <- 
      setElementToMetadata(object,
                           name = "DataProcessing",
                           subName = "Transformation",
                           content = transformation)
    
    # initiate:
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "Normalization", 
                           content =  list())
    
    object <- 
      setElementToMetadata(object, 
                           name    = "PCAlist", 
                           subName = "norm", 
                           content =  NULL)
    
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
                        transformMethod = NULL,
                        userTransMethod = "unknown"
  ){
    
    object[[SE.name]] <- 
      runTransformData(object[[SE.name]], 
                       transformMethod = transformMethod,
                       userTransMethod = userTransMethod)
    
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
#' The normalisation method is chosen according to dataset omic type 
#' (RNAseq: TMM, metabolics/proteomics: median) 
#' }
#' @param normMethod Normalization method. Accepted values: TMM for RNAseq, and
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
                        normMethod = NULL,
                        userNormMethod = "unknown"
  ){
    
    # accepted value for normMethod
    # default value : 1rst element
    default.methods <- 
      switch (getOmicsTypes(object),
              "RNAseq"       = c("TMM"),
              "proteomics"   = c("median", "totalSum", "none"),
              "metabolomics" = c("median", "totalSum", "none")
      )
    
    if(.isNormalized(object)) stop("Data is already transformed!")
    
    # check normMethod param
    if(is.null(normMethod)) normMethod <- default.methods[1]
    if(!normMethod %in% default.methods)
      stop(normMethod, 
           " is not an allowed value for the parameter normMethod.",
           " Accepted values: ", paste0(default.methods, collapse = ", "))
    
    # check filtering status 
    if (is.null(getFilterSettings(object)$method))
      stop("Before transforming the ", getOmicsTypes(object), " data, ", 
           "you must first run the feature filtering. ",
           "See ?runDataProcessing.")
    
    # check if proteomics or metabolomics data are transformed 
    if(getOmicsTypes(object) != "RNAseq" && is.null(getTransSettings(object)$method))
      stop("Before normalizing the ", getOmicsTypes(object), " data, ", 
           "you must first run the transformation. See ?runDataProcessing.")
    
    # apply trans
    object2 <- getProcessedData(object, trans = TRUE)
    
    # calculation normalization coefficient
    coefNorm <- 
      switch(normMethod,
             "TMM"      = .tmmNormalization(object2),
             "median"   = .medianNormalization(object2),
             "totalSum" = .totalSumNormalization(object2),
             "none"     =  rep(1, ncol(assay(object2)))
      )
    
    # output
    Normalization <- list(
      setting = 
        list(method = normMethod,
             suppInfo = switch(normMethod, 
                               "none" = userNormMethod, 
                               NULL)),
      results = list(coefNorm = coefNorm),
      normalized = FALSE
    )
    
    message("[RFLOMICS] #       method: ", 
            switch(normMethod, 
                   "none" = paste0("already normalized (", userNormMethod, ")"), 
                   normMethod))
    
    object <- 
      setElementToMetadata(object,
                           name    = "DataProcessing",
                           subName = "Normalization",
                           content =  Normalization)
    
    # initiate:
    object <- 
      setElementToMetadata(object, 
                           name    = "PCAlist", 
                           subName = "norm", 
                           content =  NULL)
    
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
                        normMethod = NULL,
                        userNormMethod = "unknown"
  ){
    
    object[[SE.name]] <-  
      runNormalization(object         = object[[SE.name]],
                       normMethod     = normMethod,
                       userNormMethod = userNormMethod)
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
#' @param raw boolean. Does the pca have to be ran on raw data or transformed 
#' and normalized data? Default is FALSE, pca is ran on transformed and 
#' normalized data.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runOmicsPCA
#' @importFrom FactoMineR PCA
#' @rdname runDataProcessing
#'
setMethod(
  f          = "runOmicsPCA",
  signature  = "RflomicsSE",
  definition = function(object, ncomp = 5, raw = FALSE) {
    
    log <- ifelse(getOmicsTypes(object) == "RNAseq", TRUE, FALSE)
    
    pseudo  <- assay(getProcessedData(object, norm = !raw, log = log))
    
    PCAlist <- getAnalysis(object, name = "PCAlist")
    PCAlist[[ifelse(isTRUE(raw), "raw", "norm")]] <- 
      PCA(t(pseudo), ncp = ncomp, graph = FALSE)
    object <- setElementToMetadata(object, name = "PCAlist", content = PCAlist)
    
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
                                ncomp = 5,
                                raw = FALSE) {
            
            object[[SE.name]] <-  runOmicsPCA(object[[SE.name]],
                                              ncomp = ncomp,
                                              raw  = raw)
            return(object)
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
        "Error: You need at least 1 biological factor with at least 2 modalities."
      output[["error"]] <- TRUE
      return(output)
    }
    # check presence of bash factors
    if (!length(getBatchFactors(object)) %in% c(1,2)){ 
      output[["messages"]] <-  
        "Error: You need at least 1 batch factor with at least 2 replicats."
      output[["error"]] <- TRUE 
      return(output)
    }
    
    #remplacer le code ci-dessus par celui en bas
    group_count <- .countSamplesPerCondition(ExpDesign, bio.fact)
    
    # check presence of relicat / batch
    # check if design is complete
    # check if design is balanced
    # check nbr of replicats
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
            
            if(is.null(omicName)) stop("Arg omicName missed")
            
            SEObject <- getRflomicsSE(object, omicName)
            
            checkExpDesignCompleteness(SEObject, sampleList = sampleList)
          })

##==== ACCESSORS ====

###==== getProcessedData ====

#' @rdname runDataProcessing
#' @name getProcessedData
#' @param filter boolean. If TRUE, returned filtred (samples/features)
#' normalized data
#' @param trans boolean. If TRUE, returned transformed data
#' @param norm boolean. If TRUE, returned normalization data
#' @param log boolean. If TRUE, returned log10 matrix data. Only for RNAseq
#' @aliases getProcessedData,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item getProcessedData: return Rflomics object with a processed data  
#'    (filtering, normalization and/or transformation)}
#' @exportMethod getProcessedData
setMethod(
  f          = "getProcessedData",
  signature  = "RflomicsSE",
  definition = function(object,
                        filter = FALSE,
                        trans = FALSE,
                        norm = FALSE,
                        log = FALSE){
    
    # to apply normalization we must apply filtering and transformation
    if(norm)  filter = trans = TRUE
    # to apply transdormation we must apply filtering
    if(trans) filter = TRUE
    
    # filter samples
    object <- .applySampleFiltering(object)
    
    # filtering process
    if(filter){
      # filter features
      object <- .applyFeatureFiltering(object)
    }
    
    # transformation process
    if(trans & getOmicsTypes(object) %in% c("metabolomics", "proteomics")){
      object <- .applyTransformation(object)
    }
    
    
    # transformation process
    if(norm){
      object <- .applyNormalization(object)
    }
    
    
    # log
    if(log){
      if(getOmicsTypes(object) == "RNAseq")
        object <- .applyLog(object, log = "log2")
      else
        warning("Log is not recommended for ", getOmicsTypes(object), " data.")
    }
    
    return(object)
  })

#' @rdname runDataProcessing
#' @name getProcessedData
#' @aliases getProcessedData,RflomicsMAE-method
#' @exportMethod getProcessedData
setMethod(f          = "getProcessedData",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name,
                                filter = FALSE,
                                trans = FALSE,
                                norm = FALSE,
                                log = FALSE){
            
            if (!SE.name %in% getDatasetNames(object)){
              stop("SE name must be part of this list of names: ",
                   getDatasetNames(object))
            }
            
            object[[SE.name]] <- 
              getProcessedData(object[[SE.name]],
                               filter = filter,
                               trans = trans,
                               norm = norm,
                               log = log)
            
            return(object)
          })


###==== getTransSettings ====

# Get transformation parameters
#' @rdname runDataProcessing
#' @name getTransSettings
#' @aliases getTransSettings,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item getTransSettings: return a list of transformation settings 
#'    of a given omic dataset}
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
#'    omic dataset}
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

###==== getFilteredFeatures ====

# Get filtred features
#' @rdname runDataProcessing
#' @name getFilteredFeatures
#' @aliases getFilteredFeatures,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item getFilteredFeatures: return a vector of filtered features of a given 
#'    omic dataset}
#' @exportMethod getFilteredFeatures 
#' @examples
#' # See runDataProcessing for an example that includes getFilteredFeatures
setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsSE",
          definition = function(object){
            
            res <- getAnalysis(object, 
                               name = "DataProcessing",
                               subName = "featureFiltering")
            
            return(res[["results"]][["filteredFeatures"]])   
          })

#' @rdname runDataProcessing
#' @exportMethod getFilteredFeatures 
#' @name getFilteredFeatures
#' @aliases getFilteredFeatures,RflomicsMAE-method
setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getFilteredFeatures(object = object[[SE.name]])
          })



###==== getSelectedSamples ====

# Get filtred samples
#' @rdname runDataProcessing
#' @name getSelectedSamples
#' @aliases getSelectedSamples,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item getSelectedSamples: return a vector of selected samples of a given 
#'    omic dataset}
#' @exportMethod getSelectedSamples 
#' @examples
#' # See runDataProcessing for an example that includes getSelectedSamples
setMethod(f          = "getSelectedSamples",
          signature  = "RflomicsSE",
          definition = function(object){
            
            selectedSamples <- 
              getAnalysis(object, 
                          name = "DataProcessing", 
                          subName = "selectedSamples")
            
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
#'    of a given omic dataset}
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
#'    of a given omic dataset}
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
            
            
            if(isFALSE(raw))
              object <- getProcessedData(object, norm = TRUE)
            
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
    
    log <- ifelse(getOmicsTypes(object) == "RNAseq", TRUE, FALSE)
    
    if(isFALSE(raw))
      object <- getProcessedData(object, norm = TRUE, log = log)
    else
      object <- getProcessedData(object, log = log)
    
    labels <- getLabs4plot(object)
    
    # object2 <- .checkTransNorm(object, raw = raw)
    pseudo <- assay(object)
    Groups <- getDesignMat(object)
    
    # omicsType <- getOmicsTypes(object2)
    # 
    # x_lab <- paste0(omicsType, " data")
    # if (omicsType == "RNAseq") {
    #   x_lab <- paste0("log2(", omicsType, " data)")
    # }
    # 
    # # Raw data
    # if (raw) {
    #   title <- paste0(omicsType, " raw data")
    # } else {
    #   # Already normalized or transformed Data
    #   title <- paste0(omicsType, " data")
    #   
    #   if (.isTransformed(object2)) {
    #     title <- paste0("Transformed (", 
    #                     getTransSettings(object2)$method, ") ", 
    #                     title)
    #   }
    #   
    #   if (.isNormalizedalized(object2)) {
    #     title <- paste0(title, " - normalization: ", 
    #                     getNormSettings(object2)$method)
    #   }
    #   
    #   if (omicsType == "RNAseq") {
    #     x_lab <- paste0("log2(", omicsType, " data)")
    #   }
    # }
    
    pseudo.gg <- pseudo %>% melt()
    colnames(pseudo.gg) <- c("features", "samples", "value")
    
    pseudo.gg <- pseudo.gg %>%
      full_join(Groups, by = "samples") %>%
      arrange(groups)
    
    pseudo.gg$samples <- factor(pseudo.gg$samples, 
                                levels = unique(pseudo.gg$samples))
    
    switch(plot,
           "density" = {
             p <- ggplot(pseudo.gg) +
               geom_density( aes(x = value, group = samples, color = groups), 
                             trim = FALSE) +
               xlab(labels$x_lab) +
               theme(legend.position = "none") +
               ggtitle(labels$title)
           },
           "boxplot" = {
             p <-  ggplot(pseudo.gg, aes(x = samples, y = value)) +
               geom_boxplot( aes(fill = groups), outlier.size = 0.3) +
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
                        groupColor = "groups") 
  {
    
    # define pca axis
    if (length(axes) != 2)
      stop("PCA axes must be a vector of length 2")
    
    PC1 <- paste("Dim.", axes[1], sep = "")
    PC2 <- paste("Dim.", axes[2], sep = "")
    
    if (PC1 == PC2) PC2 <- PC1 + 1

    # get labels
    log <- ifelse(getOmicsTypes(object) == "RNAseq", TRUE, FALSE)
    object <- getProcessedData(object, norm = !raw, log = log)
    
    # get pca score
    ExpDesign <- getDesignMat(object)
    rawnorm <- ifelse(isTRUE(raw), "raw", "norm")
    score <- as.data.frame(metadata(object)$PCAlist[[rawnorm]]$ind$coord[, axes])
    score$samples <- row.names(score)
    score <- right_join(score, ExpDesign, by = "samples")
    
    var1 <- round(metadata(object)$PCAlist[[rawnorm]]$eig[axes, 2][1], digits = 3)
    var2 <- round(metadata(object)$PCAlist[[rawnorm]]$eig[axes, 2][2], digits = 3)
    
    # plot
    labels <- getLabs4plot(object)
    p <- ggplot(score, x = PC1, y = PC2, aes(color = groupColor))  +
      geom_point(size = 2) +
      geom_text(aes(label = samples), 
                size = 2, vjust = "inward", hjust = "inward") +
      xlab(paste(PC1, " (", var1, "%)", sep = "")) +
      ylab(paste(PC2, " (", var2, "%)", sep = "")) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      theme(
        strip.text.x =  element_text(size = 8, face = "bold.italic"),
        strip.text.y =  element_text(size = 8, face = "bold.italic")
      ) +
      ggtitle(labels$title)
    
    # ellipse corr
    aa <- select(score, all_of(groupColor), all_of(PC1), all_of(PC2))
    bb <- coord.ellipse(aa, bary = TRUE)
    p <- p + geom_polygon(
      data = bb$res,
      x = PC1, y = PC2, aes(fill = groupColor),
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
                                            message= check$message)
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

## ---- CHECK ----

### ---- isProcessedData ----
#' @rdname runDataProcessing
#' @name isProcessedData
#' @param filter boolean. If TRUE, check if data is filtred (low counts/RNAseq)
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
