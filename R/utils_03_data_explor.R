### ============================================================================
### [03_data_processing] internal functions
### ----------------------------------------------------------------------------
# N. Bessoltane
# D. Charif
# A. Hulot

# ---- transformation ----

# .applyTrans_intensities: apply the transformation method.
#' @title applyTrans_intensities
#'
#' @param object An object of class \link{RflomicsSE}
#' @param method tranformation method
#' @keywords internal
#' @noRd
#'
.applyTrans_intensities <- function(object, 
                                   method = c("log10", "log1p", "log2", "squareroot", "none")) {
  
  method <- match.arg(method)
  
  assayRaw <- assay(object, withDimnames = TRUE)
  
  if (any(assayRaw < 0, na.rm = TRUE) &&
      method %in% c("log1p", "log2", "log10"))
    stop("Cannot use log transformation on negative values. Please check your data.")
  
  assay(object) <- 
    switch(
      method,
      "log2"       = log2(assayRaw + 10^-6*min(assayRaw[assayRaw != 0])),
      "log10"      = log10(assayRaw + 10^-6*min(assayRaw[assayRaw != 0])),
      #"squareroot" = assay(object) <- sqrt(assayRaw),
      "none"       = assay(object) <- assayRaw,
      {
        stop("Could not recognize the transformation method. ",
             "No transformation applied. Please check your parameters.")
      } # default is none
    )
  
  # output
  transformation <-
    list(
      setting = list(method = method),
      transformed = TRUE
    )
  
  message("[RFLOMICS] #       method: ", method)
  
  object <-
    setElementToMetadata(object,
                         name    = "DataProcessing",
                         subName = "Transformation",
                         content = transformation)
  return(object)
}


# .applyTrans_readcounts: apply the transformation method.
#' @title apply_transformation
#' @param object An object of class \link{RflomicsSE}
#' @param method tranformation method
#' @keywords internal
#' @noRd
#'
.applyTrans_readcounts <- function(object, 
                                   method = c("log2", "none")) {
  
  method <- match.arg(method)
  
  assayRaw <- assay(object, withDimnames = TRUE)
  
  if (any(assayRaw < 0, na.rm = TRUE) &&
      method %in% c("log1p", "log2", "log10"))
      stop("Cannot use log transformation on negative values. Please check your data.")

  assay(object) <- 
    switch(
      method,
      "log2" = log2(assayRaw + 1),
      "none" = assayRaw,
      {
        stop("Could not recognize the transformation method. ",
             "No transformation applied. Please check your parameters.")
      } # default is none
    )

  # output
  transformation <-
    list(
      setting = list(method = method),
      transformed = TRUE
    )
  
  message("[RFLOMICS] #       method: ", method)
  
  object <-
    setElementToMetadata(object,
                         name    = "DataProcessing",
                         subName = "Transformation",
                         content = transformation)

  return(object)
}

# ---- normalization - RNAseq/prot/meta data ----

#' @description
#' .applyNorm_readcounts: apply the normalization method
#' @title .applyNorm
#' @param object An object of class \link{RflomicsSE}
#' @description apply the normalization to the assay .
#' @keywords internal
#' @noRd
#'
.applyNorm_readcounts <- function(object, method = c("TMM", "none")) {
  
  method <- match.arg(method)
  
  assayRaw <- assay(object)
  
  if(method == "TMM"){
    coefNorm       <- .tmmNormalization(object)
    scales_factors <- coefNorm$norm.factors * coefNorm$lib.size
    scales_factors <- scales_factors / mean(scales_factors)
    assay(object)  <- scale(assayRaw, center = FALSE, scale = scales_factors)
  }
  else if(method == "none"){
    coefNorm <- rep(1, ncol(object))
  }

  # output
  Normalization <- list(
    setting = list(method = method),
    results = list(coefNorm = coefNorm),
    normalized = TRUE
  )
  
  message("[RFLOMICS] #       method: ", method)
  
  object <-
    setElementToMetadata(object,
                         name    = "DataProcessing",
                         subName = "Normalization",
                         content =  Normalization)

  return(object)
}

#' @description
#' .applyNorm_intensities: apply the normalization method
#' @title .applyNorm
#' @param object An object of class \link{RflomicsSE}
#' @description apply the normalization to the assay .
#' @keywords internal
#' @noRd
#'
.applyNorm_intensities <- function(object, method = c("median", "totalSum", "none")) {
  
  method <- match.arg(method)
  
  assayRaw <- assay(object)
  
  if(method == "median"){
    coefNorm <- .medianNormalization(object)
    assay(object) <- sweep(assayRaw, 2, coefNorm, "-")
    
  }
  else if(method == "totalSum"){
    coefNorm <- .totalSumNormalization(object)
    assay(object) <- sweep(assayRaw, 2, coefNorm, "/")
  }
  else if(method == "none"){
    coefNorm <- rep(1, ncol(assayRaw))
    assay(object) <- assayRaw
  }
  
  # output
  Normalization <- list(
    setting = list(method = method),
    results = list(coefNorm = coefNorm),
    normalized = TRUE
  )
  
  message("[RFLOMICS] #       method: ", method)
  
  object <-
    setElementToMetadata(object,
                         name    = "DataProcessing",
                         subName = "Normalization",
                         content =  Normalization)
  
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
    apply(assay(object), 2, function(sample_vect) {
      median(sample_vect, na.rm = TRUE)
      })

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
#' @importFrom edgeR DGEList normLibSizes
#' @noRd

.tmmNormalization <- function(object){

  groups <- getDesignMat(object)
  counts <- assay(object)

  dge <- DGEList(counts=counts, group=groups$groups)
  dge <- edgeR::normLibSizes(dge,method="TMM")
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
    apply(assay(object), 2, function(sample_vect) 
      {sum(sample_vect^2, na.rm = TRUE)
      })

  return(coef)
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

.isImputed <- function(object) {
  Imputation <-
    getAnalysis(object, name = "DataProcessing", subName = "Imputation")
  
  if(length(Imputation) == 0) return(FALSE)
  Imputation[["imputed"]]
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
#' @importFrom dplyr count group_by_at
#' @return a data.frame with sample count per condition
#' @noRd
.countSamplesPerCondition <- function(expDesign, bioFactors) {

  #remplacer le code ci-dessus par celui en bas
  group_count <- dplyr::group_by_at(expDesign, bioFactors) %>%
      dplyr::count(name = "Count")

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
