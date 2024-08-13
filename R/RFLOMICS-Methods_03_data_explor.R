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
#' @param lowCountFiltering_strategy strategy of RNAseq low count filtering. 
#' Mandatory for RNAseq data. Default value: "NbReplicates".
#' @param lowCountFiltering_CPM_Cutoff CPM cutoff for RNAseq low count filtering.
#'  Mandatory for RNAseq data. Default value: 1.
#' @param normMethod of normalization. Mandatory for RNAseq data. 
#' Default value: RNAseq = TMM.
#' @param transformMethod method of transformation.
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
                        samples=colnames(object),
                        lowCountFiltering_strategy = "NbReplicates", 
                        lowCountFiltering_CPM_Cutoff = 1, 
                        normMethod = "none", 
                        transformMethod = "none")
  {
    
    # keep selected samples
    if(!is.null(samples))
      samples <- colnames(object)
    
    message("[RFLOMICS] #    => select samples...")
    object <- runSampleFiltering(object, samples)
    
    if(nrow(getDesignMat(object)) == 0) 
      stop("no samples in object!")
    
    # spported values:
    lowCountFiltering_strategy.sup     <- c("NbReplicates","NbConditions")
    transformation_method.sup          <- c("log2", "none")
    normalisation_method.abundance.sup <- c("median", "totalSum", "none")
    normalisation_method.count.sup     <- c("TMM")
    
    # apply data processing
    switch(
      object@metadata$omicType,
      "RNAseq" = {
        
        # Filter low abundance
        message("[RFLOMICS] #    => Low counts Filtering...")
        if(is.null(lowCountFiltering_strategy)   || 
           !lowCountFiltering_strategy %in% lowCountFiltering_strategy.sup) 
          stop("the low count filtering strategy : ", 
               lowCountFiltering_strategy, 
               " isn't supported by rflomics package. Supported values : ",  
               paste(lowCountFiltering_strategy.sup, collapse = ", "))
        
        if(is.null(lowCountFiltering_CPM_Cutoff) || 
           !is.numeric(lowCountFiltering_CPM_Cutoff)) 
          stop(lowCountFiltering_CPM_Cutoff, " must be a integer value > 1")
        
        SE.processed <- 
          filterLowAbundance(object = object, 
                             filterMethod= "CPM", 
                             filterStrategy = lowCountFiltering_strategy, 
                             cpmCutoff = lowCountFiltering_CPM_Cutoff)
        
        # Run Normalisation 
        message("[RFLOMICS] #    => Counts normalization...")
        if(is.null(normMethod) || normMethod != "TMM"){
          normMethod <- "TMM"
          warning("only ", normalisation_method.count.sup, 
                  " method is supported for ", 
                  getOmicsTypes(object), " normalisation.")
        }
        
        SE.processed <- 
          runNormalization(SE.processed, 
                           normMethod = normMethod)
      },
      {
        message("[RFLOMICS] #    => transformation data...")
        if(is.null(transformMethod)) transformMethod <- "none"
        if(! transformMethod %in% transformation_method.sup) 
          stop("the transformation method ", transformMethod,
               " is not support in rflomics package.Supported values : ", 
               paste(transformation_method.sup, collapse = ", "))
        SE.processed <- 
          runTransformData(object, 
                           transformMethod = transformMethod)
        
        message("[RFLOMICS] #    => Run normalization...")
        if(is.null(normMethod)) normMethod <- "none"
        if(! normMethod %in% normalisation_method.abundance.sup) 
          stop("the normalisation method ", normMethod,
               " is not support in rflomics package. Supported values : ", 
               paste(normalisation_method.abundance.sup, collapse = ", "))
        SE.processed <- 
          runNormalization(SE.processed, 
                           normMethod = normMethod)
      }
    )
    
    # Run PCA for filtered & normalized data 
    message("[RFLOMICS] #    => Compute PCA ")
    SE.processed <- runOmicsPCA(SE.processed, ncomp = 5, raw = FALSE)  
    SE.processed@metadata$DataProcessing[["done"]] <- TRUE
    
    return(SE.processed)
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
                        lowCountFiltering_strategy = "NbReplicates", 
                        lowCountFiltering_CPM_Cutoff = 1, 
                        normMethod= "none", 
                        transformMethod = "none"){
    
    if (!SE.name %in% names(object))
      stop("SE name must be part of this list of names: ",
           getDatasetNames(object))
    
    SE.processed <- 
      runDataProcessing(
        object = object[[SE.name]],
        samples = samples, 
        lowCountFiltering_strategy = lowCountFiltering_strategy, 
        lowCountFiltering_CPM_Cutoff = lowCountFiltering_CPM_Cutoff,
        normMethod= normMethod, 
        transformMethod = transformMethod)
    
    object[[SE.name]] <- SE.processed
    
    return(object)
  })

###==== filterLowAbundance ====

# METHOD to filter data

# Cette method est propre au RNASEQ => Est-ce que c'est vraiment ce que l'on souhaite ?
# Plutot qu'une fonction interface pour tous les omics ?
# Pourquoi ne pas avoir utilisée directement la fonction de edgeR ?
#' @rdname runDataProcessing
#' @aliases filterLowAbundance,RflomicsSE-method
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
#' @exportMethod filterLowAbundance
#' @importFrom edgeR DGEList filterByExpr cpm
#' @seealso edgeR::filterByExpr
setMethod(f         = "filterLowAbundance",
          signature = "RflomicsSE",
          definition = function(object, 
                                filterMethod= "CPM", 
                                filterStrategy = "NbConditions", 
                                cpmCutoff = 5){
            
            
            if (getOmicsTypes(object) != "RNAseq")
              stop("Can't apply this method to omics types other than RNAseq.")
            
            
            if (isFALSE(filterStrategy %in% c("NbReplicates","NbConditions")))
              stop("filterStrategy argument must be one of this tow options : ", 
                   "NbReplicates or NbConditions")
            
            if(filterMethod != "CPM")
              stop("filterMethod argument must be one of this tow options : ", 
                   "CPM")
            
            assayFilt  <- assay(object)
            
            # nbr of genes with 0 count
            genes_flt0  <- object[rowSums(assayFilt) <= 0, ]@NAMES
            
            # remove 0 count
            objectFilt  <- object[rowSums(assayFilt)  > 0, ]
            assayFilt   <- assay(objectFilt)
            
            # filter cpm
            Groups       <- getDesignMat(object)
            NbReplicate  <- table(Groups$groups)
            NbConditions <- length(unique(Groups$groups))
            
            # low count filtering
            keep <-
              switch(filterStrategy,
                     "NbConditions" = { 
                       rowSums(cpm(assayFilt) >= cpmCutoff) >= NbConditions },
                     "NbReplicates" = { 
                       rowSums(cpm(assayFilt) >= cpmCutoff) >= min(NbReplicate)},
                     "filterByExpr" = { 
                       dge <- DGEList(counts = assayFilt, 
                                      genes = rownames(assayFilt))
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
            
            object <- 
              setElementToMetadata(object, 
                                   name    = "DataProcessing", 
                                   subName = "featureFiltering", 
                                   content =  Filtering)
            
            return(object)
          })


#' @rdname runDataProcessing
#' @name filterLowAbundance
#' @aliases filterLowAbundance,RflomicsMAE-method
#' @exportMethod filterLowAbundance
setMethod(f          = "filterLowAbundance",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name, 
                                filterMethod= "CPM", 
                                filterStrategy = "NbConditions", 
                                cpmCutoff = 5){
            
            object[[SE.name]] <-  
              filterLowAbundance(object          = object[[SE.name]], 
                                 filterStrategy = filterStrategy,
                                 filterMethod   = filterMethod,
                                 cpmCutoff      = cpmCutoff)
            
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
#' @importFrom dplyr filter
setMethod(
  f          = "runSampleFiltering",
  signature  = "RflomicsSE",
  definition = function(object, samples = NULL) {
    
    if(is.null(samples))
      return(object)
    
    # check if samples overlap with
    if(any(!samples %in% colnames(object)))
      stop("Some sample names are not part of the colnames of the object")
    
    object <- 
      setElementToMetadata(object, 
                           name    = "DataProcessing", 
                           subName = "selectedSamples",
                           content = samples)
    
    
    
    # keep selected samples
    # keep only samples in data matrix, and colData
    SE.new <- object[, samples]
    
    SE.new <- .updateColData(SE.new)
    selectedContrasts  <- getSelectedContrasts(object)
    selectedContrasts2 <- updateSelectedContrasts(SE.new, selectedContrasts)

    object <- setSelectedContrasts(object, selectedContrasts2)
    
    return(object)
  })


#' @rdname runDataProcessing
#' @aliases runSampleFiltering,RflomicsMAE-method
#' @exportMethod runSampleFiltering
setMethod(f          = "runSampleFiltering",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name,
                                samples=NULL) {
            
            object[[SE.name]] <- 
              runSampleFiltering(object[[SE.name]], samples = samples)
            
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
#' @param transformMethod The transformation to store in the metadata or to 
#' store and apply if modifyAssay is TRUE.
#' @param modifyAssay Boolean. Do the transformation need to be applied on 
#' the data? The raw data will be replaced by the transformed ones.
#' @exportMethod runTransformData
setMethod(f          = "runTransformData",
          signature  = "RflomicsSE",
          definition = function(object, 
                                transformMethod = NULL, 
                                modifyAssay = FALSE){
            
            if (getOmicsTypes(object) == "RNAseq")
              stop("It is not recommended to transform RNAseq data.")
            
            # accepted value for normMethod
            # default value : 1rst element
            default.methods <- 
              switch (getOmicsTypes(object),
                      "proteomics"   = c("log2", "none"),
                      "metabolomics" = c("log2", "none")
              )
            
            # check normMethod param
            if(is.null(transformMethod)){
              transformMethod <- default.methods[1]
              message("The default method applied will be: ", transformMethod)
              
            }else if(!transformMethod %in% default.methods){
              stop(transformMethod, 
                   " is not an allowed value for the parameter transformMethod",
                   " Accepted values: ", paste0(default.methods, collapse = ", "))
            }
            
            # output
            transformation <- list(setting  = list(method = transformMethod), 
                                   results  = NULL,  
                                   transformed = FALSE)
            
            object <- 
              setElementToMetadata(object,
                                   name = "DataProcessing",
                                   subName = "Transformation",
                                   content = transformation)
            
            # @audrey je ne comprends pas l'interet de cette ligne de cmd
            if (modifyAssay) {
              object <-  .applyTransformation(object)
            }
            
            return(object)
          })

#' @rdname runDataProcessing
#' @name runTransformData
#' @aliases runTransformData,RflomicsMAE-method
#' @exportMethod runTransformData
setMethod(f          = "runTransformData",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name, 
                                transformMethod = NULL, 
                                modifyAssay = FALSE){
            
            object[[SE.name]] <-  runTransformData(object[[SE.name]], 
                                                   transformMethod = transformMethod, 
                                                   modifyAssay = modifyAssay)
            
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
#' @param modifyAssay Does the normalization have to be applied or just stored 
#' for later? Recommended it stays FALSE.
#' @return An object of class \link{RflomicsSE}
#' The applied normalization method and computed scaling factors 
#' (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "method" and 
#' "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a 
#' \link{RflomicsSE} object.
#' @exportMethod runNormalization
setMethod(f          = "runNormalization",
          signature  = "RflomicsSE",
          definition = function(object, 
                                normMethod = NULL, 
                                modifyAssay = FALSE){
            
            # accepted value for normMethod
            # default value : 1rst element
            default.methods <- 
              switch (getOmicsTypes(object),
                      "RNAseq"       = c("TMM"),
                      "proteomics"   = c("median", "totalSum", "none"),
                      "metabolomics" = c("median", "totalSum", "none")
              )
            
            # check normMethod param
            if(is.null(normMethod)){
              normMethod <- default.methods[1]
              message("The default method applied will be: ", normMethod)
              
            }else if(!normMethod %in% default.methods){
              stop(normMethod, 
                   " is not an allowed value for the parameter normMethod.",
                   " Accepted values: ", paste0(default.methods, collapse = ", "))
            }
            
            # check if proteomics or metabolomics data are transformed 
            if (getOmicsTypes(object) %in% c("proteomics", "metabolomics") &
                is.null(getTransSettings(object)$method))
              warning(getOmicsTypes(object), 
                      " data should be transformed before normalization.")
            
            if (getOmicsTypes(object) == "RNAseq" &
                is.null(getFilterSettings(object)$method))
              warning(getOmicsTypes(object), 
                      " data should be filtered (low count).")
            
            object2 <- getProcessedData(object, filter = TRUE, trans = TRUE)
            
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
              setting = list(method = normMethod),
              results = list(coefNorm = coefNorm),
              normalized = FALSE
            )
            
            object <- 
              setElementToMetadata(object,
                                   name    = "DataProcessing",
                                   subName = "Normalization",
                                   content =  Normalization)
            
            # @audrey, je ne vois pas l'interet de cette ligne. Nadia
            if (modifyAssay) object <- .applyNorm(object)
            
            return(object)
          })

#' @rdname runDataProcessing
#' @name runNormalization
#' @aliases runNormalization,RflomicsMAE-method
#' @exportMethod runNormalization
setMethod(f          = "runNormalization",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name, 
                                normMethod = NULL, 
                                modifyAssay = FALSE){
            
            object[[SE.name]] <-  
              runNormalization(object       = object[[SE.name]],
                               normMethod   = normMethod,
                               modifyAssay  = modifyAssay)
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
setMethod(f          = "runOmicsPCA",
          signature  = "RflomicsSE",
          definition = function(object, ncomp = 5, raw = FALSE) {
            
            log <- ifelse(getOmicsTypes(object) == "RNAseq", TRUE, FALSE)
            
            if(isFALSE(raw)){
              object2 <- getProcessedData(object, norm = TRUE, log = log)
              tag = "norm"
            }
            else{
              object2 <- getProcessedData(object, log = log)
              tag = "raw"
            }
            
            pseudo  <- assay(object2)
            
            object@metadata[["PCAlist"]][[tag]] <-
              PCA(t(pseudo), ncp = ncomp, graph = FALSE)
            
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
#' @exportMethod checkExpDesignCompleteness
setMethod(f         = "checkExpDesignCompleteness",
          signature = "RflomicsSE",
          definition <- function(object, sampleList=NULL){
            
            if(!is.null(sampleList))
              object <- object[,sampleList]
            
            output <- list()
            output[["error"]] <- FALSE
            
            # Only works with bio and batch factors for the rest of the function
            ExpDesign <- getDesignMat(object)
            bio.fact  <- getBioFactors(object)
            
            # check presence of bio factors
            if (!length(getBioFactors(object)) %in% seq_len(3)){ 
              output[["messages"]] <-  "Error: You need at least 1 biological factor with at least 2 modalities."
              output[["error"]]    <- TRUE
              return(output)
            }
            # check presence of bash factors
            if (!length(getBatchFactors(object)) %in% c(1,2)){ 
              output[["messages"]] <-  "Error: You need at least 1 batch factor with at least 2 replicats."
              output[["error"]]    <- TRUE 
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
              
              output[["messages"]] <- "The experimental design is complete but not balanced."
            }
            else{
              output[["messages"]] <- "The experimental design is complete and balanced."
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

###==== getOmicData ====

#' @rdname runDataProcessing
#' @name getOmicData
#' @param processedData boolean If TRUE, returned filtred data 
#' @param NormTrans boolean. If TRUE, apply normalisation and/or transformation
#' coefficients.
#' @param log boolean. If TRUE, apply log2 tranformation (for RNAseq data)
#' @aliases getOmicData,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item getOmicData: return a processed data matrix  
#'    (filtering, normalization and/or transformation)}
#' @exportMethod getOmicData
#' @importFrom SummarizedExperiment assay
setMethod(f          = "getOmicData",
          signature  = "RflomicsSE",
          definition = function(object, 
                                processedData = FALSE,
                                NormTrans = FALSE, 
                                log = FALSE,
                                ...){
            
            dataType <- getOmicsTypes(object)
            
            if(isTRUE(processedData)){
              
              omicData <- getProcessedData(object, 
                                           NormTrans = NormTrans)
            }else{
              
              omicData <- assay(object)
            }
            
            if(isTRUE(log)){
              
              omicData <- 
                switch(dataType,
                       "RNAseq" = log2(omicData + 1),
                       log2(omicData + 10^10)
                )
            }
            
            return(data.frame(omicData))
          })

#' @rdname runDataProcessing
#' @name getOmicData
#' @aliases getOmicData,RflomicsMAE-method
#' @exportMethod getOmicData
setMethod(f          = "getOmicData",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name,
                                processedData = FALSE,
                                NormTrans = FALSE, 
                                log = FALSE,
                                ...){
            
            if (!SE.name %in% getDatasetNames(object)){
              stop("SE name must be part of this list of names: ",
                   getDatasetNames(object))
            }
            
            pseudo <- getOmicData(object[[SE.name]],
                                  processedData = processedData,
                                  NormTrans = NormTrans, 
                                  log = log,
                                  ...)
            
            return(pseudo)
          })


###==== getProcessedData ====

#' @rdname runDataProcessing
#' @name getProcessedData
#' @param NormTrans boolean. If TRUE, returned Data are normalized 
#' and transformed.
#' @aliases getProcessedData,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item getProcessedData: return Rflomics object with a processed data  
#'    (filtering, normalization and/or transformation)}
#' @exportMethod getProcessedData
#' @importFrom SummarizedExperiment assay
setMethod(
  f          = "getProcessedData",
  signature  = "RflomicsSE",
  definition = function(object,
                        filter = FALSE,
                        trans = FALSE,
                        norm = FALSE,
                        log = FALSE){
    
    
    if(norm)  filter = trans = TRUE
    if(trans) filter = TRUE
    
    # filtering process
    if(filter){
      # filter samples
      selectedSamples <- getSelectedSamples(object)
      object <- object[, selectedSamples]
      object <- .updateColData(object)
      
      # filter features
      object <- .applyFiltering(object)
    }
    
    
    # transformation process
    if(trans && getOmicsTypes(object) %in% c("metabolomics", "proteomics"))
      object <- .applyTransformation(object)
    
    # transformation process
    if(norm)
      object <- .applyNorm(object)
    
    # log
    if(log && getOmicsTypes(object) == "RNAseq")
      object <- .applyLog(object, log = "log2")
    
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
            return(object@metadata$DataProcessing$Transformation$setting)   
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
            return(object@metadata$DataProcessing$featureFiltering$setting)   
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
            return(object@metadata$DataProcessing$featureFiltering$results$filteredFeatures)   
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
#' @aliases getFilteredSamples,RflomicsSE-method
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
setMethod(f          = "getCoeffNorm",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(metadata(object)[["DataProcessing"]][["Normalization"]][["results"]][["coefNorm"]])
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
            return(object@metadata$DataProcessing$Normalization$setting)
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



###==== set selected sample list ====

#' @rdname runDataProcessing
#' @name setSelectedSamples
#' @param samples sample names to keep
#' @aliases setSelectedSamples,RflomicsSE-method
#' @section Accessors: 
#' \itemize{
#'    \item setSelectedSamples: }
#' @exportMethod setSelectedSamples
#' @examples
#' # See runDataProcessing for an example that includes setSelectedSamples
setMethod(f          = "setSelectedSamples",
          signature  = "RflomicsSE",
          definition = function(object, samples = NULL){
            
            if(is.null(samples))
              return(object)
            
            # check if samples overlap with
            if(any(!samples %in% colnames(object)))
              stop("Some sample names are not part of the colnames of the object")
            
            object <- 
              setElementToMetadata(object, 
                                   name = "DataProcessing", 
                                   subName = "selectedSamples",
                                   content = samples)
            
            return(object)
          })

#' @rdname runDataProcessing
#' @name setSelectedSamples
#' @aliases setSelectedSamples,RflomicsMAE-method
#' @exportMethod setSelectedSamples 

setMethod(f          = "setSelectedSamples",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, samples = NULL){
            
            object[[SE.name]] <- 
              setSelectedSamples(object[[SE.name]], samples = samples)
            
            return(object)
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
#' @importFrom ggplot2 ggplot aes ggtitle element_text 
#' @importFrom ggplot2 theme labs ylab geom_bar
#' @importFrom dplyr full_join arrange
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
#' @importFrom ggplot2 ggplot geom_boxplot geom_density aes
#' @importFrom ggplot2 xlab theme element_text ylab margin ggtitle
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join arrange
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
    #   if (.isNorm(object2)) {
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
#' @param raw This argument indicates whether the scaled PCA has to be 
#' performed on raw [\sQuote{raw}] or normalized [\sQuote{norm}] data.
#' @param axes A vector giving the two axis that have to be drawn for the 
#' factorial map
#' @param groupColor All combination of level's factor
#' @importFrom dplyr mutate full_join select right_join
#' @importFrom FactoMineR coord.ellipse
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text aes xlab ylab 
#' geom_hline geom_vline geom_vline element_text ggtitle geom_polygon
#' @exportMethod plotOmicsPCA
#' @examples
#' # See runDataProcessing for an example that includes plotOmicsPCA
setMethod(
  f          = "plotOmicsPCA",
  signature  = "RflomicsSE",
  definition = function(object,
                        raw = c("raw", "norm"),
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
    if(raw != "raw")
      object <- getProcessedData(object, norm = TRUE, log = log)
    else
      object <- getProcessedData(object, log = log)
    
    labels <- getLabs4plot(object)
    
    # get pca score
    ExpDesign <- getDesignMat(object)
    score <- as.data.frame(metadata(object)$PCAlist[[raw]]$ind$coord[, axes])
    score$samples <- row.names(score)
    score <- right_join(score, ExpDesign, by = "samples")
    
    var1 <- round(metadata(object)$PCAlist[[raw]]$eig[axes, 2][1], digits = 3)
    var2 <- round(metadata(object)$PCAlist[[raw]]$eig[axes, 2][2], digits = 3)
    
    
    
    p <- ggplot(score, aes_string(x = PC1, y = PC2, color = groupColor))  +
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
    aa <- select(score, all_of(groupColor), PC1, PC2)
    bb <- coord.ellipse(aa, bary = TRUE)
    p <- p + geom_polygon(
      data = bb$res,
      aes_string(x = PC1, y = PC2, fill = groupColor),
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
                        raw = c("raw", "norm"),
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
