### ============================================================================
### [03_data_processing] internal functions
### ----------------------------------------------------------------------------
# N. Bessoltane
# D. Charif
# A. Hulot


# ---- sample filtering ----

# .applySampleFiltering
#' @title .applySampleFiltering
#'
#' @param object An object of class \link{RflomicsSE}
#' @description apply the filtering to the assay. Usually.
#' @keywords internal
#' @noRd
#'

.applySampleFiltering <- function(object) {
  
  selectedSamples <- getSelectedSamples(object)
  if(setequal(selectedSamples, colnames(object))) return(object)
  
  object <- object[, selectedSamples]
  object <- .updateColData(object)
  # update contrast list
  selectedContrasts  <- getSelectedContrasts(object)
  selectedContrasts2 <- updateSelectedContrasts(object, selectedContrasts)
  
  object <- setSelectedContrasts(object, selectedContrasts2)
  
  return(object)
}


# ---- features filtering - RNAseq data - low count - CMP ----

# .applyFeatureFiltering
#' @title .applyFeatureFiltering
#'
#' @param object An object of class \link{RflomicsSE}
#' @description apply the filtering to the assay. Usually.
#' @keywords internal
#' @noRd
#'
.applyFeatureFiltering <- function(object) {
  
  if (.isFiltered(object)) {
    warning("The data were already filterd beforehand! method: ",
            getFilterSettings(object)$method)    
    return(object)
  }
  
  if(is.null(getFilterSettings(object)$method)){
    warning("No filtering method specified, see ?runDataProcessing")
    return(object)
  }
  
  filtering.res <- 
    getAnalysis(object, 
                name = "DataProcessing", 
                subName = "featureFiltering")
  
  if(getFilterSettings(object)$method == "MVI"){
    omics.df  <- assay(object)
    # omics.df[is.na(omics.df)] <- 0
    # omics.df <- as.data.frame(lapply(omics.df, as.numeric))
    omics.df[omics.df == 0] <- getFilterSettings(object)$minValue
    assay(object) <- omics.df
  }
  
  filteredFeatures <- getFilteredFeatures(object)
  if(!is.null(filteredFeatures)){
    object <- 
      object[setdiff(names(object), filteredFeatures)]
    
    filtering.res[["filtered"]] <- TRUE
  }
  
  object <- setElementToMetadata(object,
                                 name    = "DataProcessing",
                                 subName = "featureFiltering",
                                 content = filtering.res)
  return(object)
}

# ---- transformation - prot/meta data ----

# .applyTransformation: apply the transformation method stored in 
# metadata(object)[["transform_method"]] and modify the assay.
#' @title apply_transformation
#'
#' @param object An object of class \link{RflomicsSE}
#' @keywords internal
#' @noRd
#'
.applyTransformation <- function(object) {
  
  if (.isTransformed(object)) {
    warning("The data were already transformed beforehand! method: ",
            getTransSettings(object)$method)
    return(object)
  }
  
  
  transform_method <- getTransSettings(object)$method
  if (is.null(transform_method)) {
    warning("No transformation method specified, see ?runDataProcessing")
    return(object)
  }
  
  assayTransform <- assay(object, withDimnames = TRUE)
  
  switch(transform_method,
         "log1p" = {
           assay(object) <- log1p(assayTransform)
         },
         "log2" = {
           assay(object) <- log2(assayTransform + 10^-10)
         },
         "log10" = {
           assay(object) <- log10(assayTransform + 10^-10)
         },
         "squareroot" = {
           assay(object) <- sqrt(assayTransform)
         },
         "none" = {
           assay(object) <- assayTransform
         },
         {
           stop("Could not recognize the transformation method. ",
                "No transformation applied. Please check your parameters.")
         } # default is none
  )
  
  metadata(object)[["DataProcessing"]][["Transformation"]][["transformed"]] <- 
    TRUE
  
  return(object)
}

# ---- normalization - RNAseq/prot/meta data ----

#' @description
#' .applyNormalization: apply the normalization method stored in 
#' metadata(object)[["Normalization"]] and modify the assay.
#' @title .applyNorm
#' @param object An object of class \link{RflomicsSE}
#' @description apply the normalization to the assay. Usually, after the transformation,
#' unless in the case of counts RNASeq data (TMM), where log2 is the second step.
#' @keywords internal
#' @noRd
#'
.applyNormalization <- function(object) {
  
  if (.isNormalized(object)) {
    warning("The data were already normalized beforehand! method: ",
            getNormSettings(object)$method)
    return(object)
  }
  
  norm_method <- getNormSettings(object)$method
  if (is.null(norm_method)) {
    warning("No normalization method specified, see ?runDataProcessing")
    return(object)
  }
  
  coefNorm <- getCoeffNorm(object)
  assayTransform <- assay(object)
  
  switch(norm_method,
         "median" = {
           assay(object) <- sweep(assayTransform, 2, coefNorm, "-")
         },
         "totalSum" = {
           assay(object) <- sweep(assayTransform, 2, coefNorm, "/")
         },
         "TMM" = {
           scales_factors <- coefNorm$norm.factors * coefNorm$lib.size
           assay(object) <- scale(assayTransform + 1, 
                                  center = FALSE, 
                                  scale = scales_factors)
         },
         "none" = {
           assay(object) <- assayTransform
         },
         {
           stop("Could not recognize the normalization method. ",
                "No normalization applied. Please check your parameters.")
         }
  )
  
  metadata(object)[["DataProcessing"]][["Normalization"]][["normalized"]] <- 
    TRUE
  
  return(object)
}

#' @title .medianNormalization
#' Interface to calculate the median normalization coefficient 
#' @param object rflomicsSE object
#' @return a data.frame with a row for each sample and columns group, lib.size 
#' and norm.factors containing the group labels, library sizes and normalization 
#' factors. Other columns can be optionally added to give more detailed sample 
#' information.
#' @keywords internal
#' @importFrom stats median
#' @noRd

.medianNormalization <- function(object){
  
  coef <- 
    apply(assay(object), 2, function(sample_vect) {median(sample_vect)})
  
  return(coef)
}

#' @title .tmmNormalization
#' Interface to the calcNormFactors function of the edgeR package  with the 
#' choosen TMM parameters as the normalization method
#' @param object rflomicsSE object
#' @return a data.frame with a row for each sample and columns group, lib.size 
#' and norm.factors containing the group labels, library sizes and normalization 
#' factors. Other columns can be optionally added to give more detailed sample 
#' information.
#' @keywords internal
#' @importFrom edgeR DGEList calcNormFactors
#' @noRd

.tmmNormalization <- function(object){
  
  groups <- getDesignMat(object)
  counts <- assay(object)
  
  dge <- edgeR::DGEList(counts=counts, group=groups$groups)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  nf  <- dge$samples
  return(nf)
}

#' @title .totalSumNormalization
#' Interface to calculate the totalSum normalization coefficient 
#' @param object rflomicsSE object
#' @return a data.frame with a row for each sample and columns group, lib.size 
#' and norm.factors containing the group labels, library sizes and normalization 
#' factors. Other columns can be optionally added to give more detailed sample 
#' information.
#' @keywords internal
#' @noRd

.totalSumNormalization <- function(object){
  
  coef <- 
    apply(assay(object), 2, function(sample_vect) {sum(sample_vect^2)})
  
  return(coef)
}

# ---- log trasnformation - RNAseq data ----

#' @title .applyLog
#'
#' @param object An object of class \link{RflomicsSE}
#' @param log log type
#' @description apply the log to the assay. Usually.
#' @keywords internal
#' @noRd
#'
.applyLog <- function(object, log = "log2") {
  
  if(getOmicsTypes(object) != "RNAseq") return(object)
  
  assay(object) <- 
    switch(log,
           "log2" = {
             if(.isNormalized(object))
               log2(assay(object))
             else
               log2(assay(object) + 1)
           }
    )
  
  metadata(object)[["DataProcessing"]][["log"]] <- "log2"
  
  return(object)
}

# ---- check data processing level ----

#' @title isFiltered, isNormalized, isTransformed,  
#'
#' @param object An object of class \link{RflomicsSE}
#' @description get if an assay has been transformed or normalized.
#' @keywords internal
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @noRd
#'
.isFiltered <- function(object) {
  featureFiltering <- 
    getAnalysis(object, name = "DataProcessing", subName = "featureFiltering")
  
  if(length(featureFiltering) == 0) return(FALSE)
  return(featureFiltering[["filtered"]])
}

.isTransformed <- function(object) {
  Transformation <- 
    getAnalysis(object, name = "DataProcessing", subName = "Transformation")
  
  if(length(Transformation) == 0) return(FALSE)
  Transformation[["transformed"]]
}

.isNormalized <- function(object) {
  Normalization <- 
    getAnalysis(object, name = "DataProcessing", subName = "Normalization")
  
  if(length(Normalization) == 0) return(FALSE)
  Normalization[["normalized"]]
}

# ---- update colData - levels ----

#' @title update colData after sample filtering
#' @param object rflomicsSE object
#' @return object rflomicsSE object
#' @keywords internal
#' @noRd
.updateColData <- function(object){
  
  colData.df <- as.data.frame(colData(object))
  
  for (factor in c(getBioFactors(object), getBatchFactors(object))){
    
    # if only one category remains after the filter, it's will be removed
    if (length(unique(colData.df[[factor]])) <= 1 ) {
      stop("The bio factor, ", factor, ", must have at least 2 levels.")
      # object[[factor]] <- NULL
      # colData.df[[factor]] <- NULL
      # factor.types <- getFactorTypes(object)
      # metadata(object)$design$factorType <- 
      #   factor.types[which(names(factor.types) != factor)]
      # # replace with setFactorTypes
    }
    else{
      F.levels <- levels(colData.df[[factor]])
      object[[factor]] <- 
        factor(colData.df[[factor]], 
               levels = intersect(F.levels, unique(colData.df[[factor]])))
    }
  }
  order_levels <- 
    with(colData.df, 
         do.call(order, 
                 colData.df[c(getBioFactors(object), getBatchFactors(object))]))
  object$samples <- 
    factor(object$samples, levels = unique(object$samples[order_levels]))
  object$groups  <- 
    factor(object$groups, levels = unique(object$groups[order_levels]))
  
  return(object)
}


## ---- countSamplesPerCondition: count nb of samples per condition to check completeness ----
#' @title countSamplesPerCondition
#' @param expDesign a data.frame with experimental design
#' @param bioFactors a vector of design bio factors
#' @return a data.frame with sample count per condition
#' @noRd
.countSamplesPerCondition <- function(expDesign, bioFactors) {
  
  #remplacer le code ci-dessus par celui en bas
  group_count <- group_by_at(expDesign, bioFactors) %>% 
    count(name = "Count")
  
  mod.fact <- lapply(names(group_count)[-ncol(group_count)], function(factor){
    unique(group_count[[factor]])
  }) 
  names(mod.fact) <- names(group_count)[-ncol(group_count)]
  
  full_join(expand.grid(mod.fact), group_count, by=bioFactors) %>% 
    mutate_at(.vars = "Count", .funs = function(x){ if_else(is.na(x), 0, x) }) %>%
    return()
}


## ---- plotExperimentalDesign ----
#' Plot the balance of data in an experimental design
#'
#' This function provides easy visualization of the balance of data in a data 
#' set given a specified experimental design. This function is useful for
#'  identifying missing data and other issues. 
#'  The core of this function is from the function ezDesign in the package ez.
#'
#' @param counts : the number of data in each cell of the design
#' @param cell_border_size : Numeric value specifying the size of 
#' the border seperating cells (0 specifies no border)
#'
#' @return A printable/modifiable ggplot2 object.
#' @keywords internal
#' @noRd
.plotExperimentalDesign <- function(counts, cell_border_size = 10, message=""){
  if (names(counts)[ncol(counts)] != "Count"){
    stop("the last column of the input data frame must be labelled Count")
  }
  if(ncol(counts) < 2){
    stop("data frame with less than 2 columns")
  }
  
  # #add color column
  # # #00BA38
  
  counts <- counts %>% 
    mutate(status = if_else(Count > 2 , "pass", 
                            if_else(Count == 2 , "warning", "error")))
  
  #list of factor names
  factors <- names(counts)[1:(dim(counts)[2]-2)]
  
  col.panel <- c("pass", "warning", "error")
  names(col.panel) <- c("#00BA38", "orange", "red")
  
  col.panel.u <- col.panel[col.panel %in% unique(counts$status)]
  
  switch (
    length(factors),
    "1" = { p <- ggplot(counts ,aes(x = !!sym(factors[1]), y = 1)) + 
      theme(axis.text.y = element_blank()) + ylab("") },
    "2" = { p <- ggplot(counts ,aes(x = !!sym(factors[1]), y = !!sym(factors[2]))) },
    "3" = {
      #get factor with min conditions -> to select for "facet_grid"
      factors.l <- lapply(factors, function(x){ length(unique(counts[[x]])) }) %>% unlist()
      names(factors.l) <- factors
      factor.min <- names(factors.l[factors.l == min(factors.l)][1])
      
      factors <- factors[factors != factor.min]
      
      #add column to rename facet_grid
      counts <- counts %>% mutate(grid = paste0(factor.min, "=",get(factor.min)))
      
      p <- ggplot(counts ,aes(x = !!sym(factors[1]), y = !!sym(factors[2]))) +
        facet_grid(grid~.) })
  
  p <- p + 
    geom_tile(aes(fill = status), color = "white",
              linewidth = 1, width = 1, height = 1) + 
    geom_text(aes(label = Count)) + 
    scale_fill_manual(values = names(col.panel.u), breaks = col.panel.u) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle(message)
  
  return(p)
}
