### ============================================================================
### [01_Load_Data] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif
# A. Hulot

# ---- ACCESSORS ----
## ---- getProjectName ----
#' @rdname RflomicsMAE-class
#' @name getProjectName
#' @aliases getProjectName,RflomicsMAE-method
#' @exportMethod getProjectName
#' @section Accessors:
#' \itemize{
#'   \item getProjectName: return a string with the name of the project}
setMethod(f          = "getProjectName",
          signature  = "RflomicsMAE",
          definition <- function(object){
            return(metadata(object)$projectName)
          })


## ---- getDesignMat:    get colData ----
#' @name getDesignMat
#' @aliases getDesignMat,RflomicsMAE-method
#' @exportMethod getDesignMat
#' @rdname RflomicsMAE-class
#' @section Accessors:
#' \itemize{
#'    \item getDesignMat: return a data.frame with experimental design.}
setMethod(f          = "getDesignMat",
          signature  = "RflomicsMAE",
          definition <- function(object){

            return(as.data.frame(colData(object)))
          })

#' @name getDesignMat,RflomicsSE-method
#' @exportMethod getDesignMat
#' @aliases getDesignMat,RflomicsSE-method
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getDesignMat: return a data.frame with experimental design.}
setMethod(f          = "getDesignMat",
          signature  = "RflomicsSE",
          definition <- function(object){

            return(as.data.frame(colData(object)))
          })

## ---- getDatasetNames: get experiment names from RflomicsMAE object ----

#' @name getDatasetNames
#' @aliases getDatasetNames,RflomicsMAE-method
#' @exportMethod getDatasetNames
#' @rdname RflomicsMAE-class
#' @section Accessors:
#' \itemize{
#'    \item getDatasetNames: return a vector with dataset names.}
setMethod(f          = "getDatasetNames",
          signature  = "RflomicsMAE",
          definition <- function(object){

            ExperimentNames <- unlist(metadata(object)$omicList)
            names(ExperimentNames) <- NULL

            return(ExperimentNames)
          })


#' @name getDatasetNames,RflomicsSE-method
#' @aliases getDatasetNames,RflomicsSE-method
#' @exportMethod getDatasetNames
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getDatasetNames: return a string with dataset name.}
setMethod(f          = "getDatasetNames",
          signature  = "RflomicsSE",
          definition <- function(object){

            return(names(metadata(object)$omicType))
          })

## ---- getOmicsTypes ----
# get omics types of experiments from RflomicsMAE object

#' @name getOmicsTypes
#' @aliases getOmicsTypes,RflomicsMAE-method
#' @exportMethod getOmicsTypes
#' @rdname RflomicsMAE-class
#' @importFrom purrr reduce
#' @section Accessors:
#' \itemize{
#'    \item getOmicsTypes: return a named vector with omics type of each
#'    dataset ("RNAseq", "proteomics", "metabolomics")}
setMethod(
  f          = "getOmicsTypes",
  signature  = "RflomicsMAE",
  definition <- function(object){

    OmicsTypes <- lapply(names(metadata(object)$omicList), function(x){

      datasetNames <- metadata(object)$omicList[[x]]
      datasetType  <- rep(x, length(datasetNames))
      names(datasetType) <- datasetNames

      return(datasetType)

    }) |> reduce(c)

    return(OmicsTypes)
  })


#' @name getOmicsTypes,RflomicsSE-method
#' @aliases getOmicsTypes,RflomicsSE-method
#' @exportMethod getOmicsTypes
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getOmicsTypes: return a named vector with omics type of dataset
#'    ("RNAseq", "proteomics", "metabolomics")}
setMethod(
  f          = "getOmicsTypes",
  signature  = "RflomicsSE",
  definition <- function(object){

    return(metadata(object)$omicType)
  })


## ---- getFactorNames:  get Factor names ----

#' @name getFactorNames
#' @exportMethod getFactorNames
#' @rdname RflomicsMAE-class
#' @aliases getFactorNames,RflomicsMAE-method
#' @section Accessors:
#' \itemize{
#'    \item getFactorNames: return a vector with the experimental factor names.}
setMethod(
  f          = "getFactorNames",
  signature  = "RflomicsMAE",
  definition <- function(object){

    return(names(metadata(object)$design$Factors.Type))
  })


#' @name getFactorNames,RflomicsSE-method
#' @aliases getFactorNames,RflomicsSE-method
#' @exportMethod getFactorNames
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getFactorNames: return a vector with the experimental factor names.}
setMethod(
  f          = "getFactorNames",
  signature  = "RflomicsSE",
  definition <- function(object){

    return(names(metadata(object)$design$factorType))
  })

## ---- getFactorTypes:  get Factor types ----

#' @name getFactorTypes
#' @exportMethod getFactorTypes
#' @aliases getFactorTypes,RflomicsMAE-method
#' @rdname RflomicsMAE-class
#' @section Accessors:
#' \itemize{
#'    \item getFactorTypes: return a named vector with experimental factor types
#'    ("bio", "batch" or "meta").}
setMethod(
  f          = "getFactorTypes",
  signature  = "RflomicsMAE",
  definition <- function(object){

    return(metadata(object)$design$Factors.Type)
  })


#' @name getFactorTypes,RflomicsSE-method
#' @exportMethod getFactorTypes
#' @aliases getFactorTypes,RflomicsSE-method
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getFactorTypes: return a named vector with experimental factor types
#'    ("bio", "batch" or "meta").}
setMethod(
  f          = "getFactorTypes",
  signature  = "RflomicsSE",
  definition <- function(object){
    return(metadata(object)$design$factorType)
  })


## ---- getBioFactors:   get bio factor ----

#' @name getBioFactors
#' @exportMethod getBioFactors
#' @rdname RflomicsMAE-class
#' @aliases getBioFactors,RflomicsMAE-method
#' @section Accessors:
#' \itemize{
#'    \item getBioFactors: return a vector with the biological factor names.}
setMethod(
  f          = "getBioFactors",
  signature  = "RflomicsMAE",
  definition <- function(object){

    factVect <- toupper(getFactorTypes(object))
    res <- names(factVect)[factVect == "BIO"]

    if(length(res) == 0) return(NULL)
    return(res)
  })

#' @name getBioFactors,RflomicsSE-method
#' @exportMethod getBioFactors
#' @aliases getBioFactors,RflomicsSE-method
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getBioFactors: return a vector with the biological factor names.}
setMethod(
  f          = "getBioFactors",
  signature  = "RflomicsSE",
  definition <- function(object){

    factVect <- toupper(getFactorTypes(object))
    res <- names(factVect)[factVect == "BIO"]

    if(length(res) == 0) return(NULL)
    return(res)
  })

## ---- getBatchFactors: get batch factor names ----

#' @name getBatchFactors
#' @exportMethod getBatchFactors
#' @aliases getBatchFactors,RflomicsMAE-method
#' @rdname RflomicsMAE-class
#' @section Accessors:
#' \itemize{
#'    \item getBatchFactors: return a vector with the batch factor names.}
setMethod(
  f          = "getBatchFactors",
  signature  = "RflomicsMAE",
  definition <- function(object){

    factVect <- toupper(getFactorTypes(object))
    res <- names(factVect)[factVect == "BATCH"]

    if(length(res) == 0) return(NULL)
    return(res)
  })


#' @name getBatchFactors,RflomicsSE-method
#' @aliases getBatchFactors,RflomicsSE-method
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getBatchFactors: return a vector with the batch factor names.}
setMethod(
  f          = "getBatchFactors",
  signature  = "RflomicsSE",
  definition <- function(object){

    factVect <- toupper(getFactorTypes(object))
    res <- names(factVect)[factVect == "BATCH"]

    if(length(res) == 0) return(NULL)
    return(res)
  })

## ---- getMetaFactors:  get meta Factor names ----

#' @name getMetaFactors
#' @aliases getMetaFactors,RflomicsMAE-method
#' @exportMethod getMetaFactors
#' @rdname RflomicsMAE-class
#' @section Accessors:
#' \itemize{
#'    \item getMetaFactors: return a vector with the metadata factor names.}
setMethod(
  f          = "getMetaFactors",
  signature  = "RflomicsMAE",
  definition <- function(object){

    factVect <- toupper(getFactorTypes(object))
    res <- names(factVect)[factVect == "META"]

    if(length(res) == 0) return(NULL)
    return(res)
  })


#' @name getMetaFactors,RflomicsSE-method
#' @aliases getMetaFactors,RflomicsSE-method
#' @exportMethod getMetaFactors
#' @rdname RflomicsSE-class
#' @section Accessors:
#' \itemize{
#'    \item getMetaFactors: return a vector with the metadata factor names.}
setMethod(
  f          = "getMetaFactors",
  signature  = "RflomicsSE",
  definition <- function(object){

    factVect <- toupper(getFactorTypes(object))
    res <- names(factVect)[factVect == "META"]

    if(length(res) == 0) return(NULL)
    return(res)
  })

## ---- getRflomicsSE:   get RflomicsSE object of one omic dataset ----

#' @name getRflomicsSE
#' @aliases getRflomicsSE,RflomicsMAE-method
#' @exportMethod getRflomicsSE
#' @rdname RflomicsMAE-class
#' @param datasetName the name of the RflomicsSE to retrieve
#' @section Accessors:
#' \itemize{
#'    \item getRflomicsSE: return a \link{RflomicsSE} object with selected
#'    dataset}
setMethod(
  f          = "getRflomicsSE",
  signature  = "RflomicsMAE",
  definition <- function(object,
                         datasetName = NULL){

    if(is.null(datasetName)) return(NULL)

    return(object[[datasetName]])
  })

## ---- getFactorModalities: ----

#' @name getFactorModalities
#' @aliases getFactorModalities,RflomicsMAE-method
#' @exportMethod getFactorModalities
#' @rdname RflomicsMAE-class
#' @param factorName factor name
#' @section Accessors:
#' \itemize{
#'    \item getFactorModalities: return a vector with the modality names of
#'    selected factor.}
setMethod(
  f          = "getFactorModalities",
  signature  = "RflomicsMAE",
  definition <- function(object, factorName){

    if(is.null(factorName)) return(NULL)
    if(!factorName %in% getFactorNames(object)) return(NULL)

    return(levels(getDesignMat(object)[[factorName]]))
  })

#' @name getFactorModalities,RflomicsSE-method
#' @aliases getFactorModalities,RflomicsSE-method
#' @exportMethod getFactorModalities
#' @rdname RflomicsSE-class
#' @param factorName factor name
#' @section Accessors:
#' \itemize{
#'    \item getFactorModalities: return a vector with the modality names of
#'    selected factor.}
setMethod(
  f          = "getFactorModalities",
  signature  = "RflomicsSE",
  definition <- function(object, factorName){

    if(is.null(factorName)) return(NULL)
    if(!factorName %in% getFactorNames(object)) return(NULL)

    return(levels(getDesignMat(object)[[factorName]]))
  })

## ---- subRflomicsMAE:  subset a RflomicsMAE from ----

#' @name subRflomicsMAE
#' @aliases subRflomicsMAE,RflomicsMAE-method
#' @rdname RflomicsMAE-class
#' @section Accessors:
#' \itemize{
#'    \item subRflomicsMAE: return a \link{RflomicsMAE-class} object with selected
#'    datasets.}
#' @param omicNames dataset name.
#' @exportMethod subRflomicsMAE
setMethod(
  f          = "subRflomicsMAE",
  signature  = "RflomicsMAE",
  definition <- function(object, omicNames = NULL){

    if(is.null(omicNames)) return(object)
    dataset.names <- names(object)

    if(!all(omicNames %in% dataset.names)) return(NULL)

    object.sub <- .tryCatch_rflomics(object[,, omicNames])
    if(!is.null(object.sub$error)) stop(object.sub$error)

    return(object.sub$result)
  })


# ---- PLOTS ----

## ---- plotDataOverview ----
#' @rdname RflomicsMAE-class
#' @aliases plotDataOverview,RflomicsMAE-method
#' @name plotDataOverview
#' @section Plots:
#' \itemize{
#'    \item plotDataOverview:
#' This function plot an overview of the loaded datasets displaying per sample
#' (n=number of entities (genes/metabolites/proteins); k=number of samples)}
#' @param omicNames a vector with dataset names
#' @param realSize booleen value, influence the display size
#' @param raw boolean. If TRUE, displays the raw data without any selection. If
#' FALSE, displays the data with removed samples.
#' @param complete.cases boolean. If true, only shows the complete cases of the
#' object.
#' @exportMethod plotDataOverview
#' @examples
#' # See createRflomicsMAE for an example that includes plotDataOverview
setMethod(
  f         = "plotDataOverview",
  signature = "RflomicsMAE",
  definition <- function(object,
                         omicNames=NULL,
                         realSize=FALSE,
                         raw = FALSE,
                         completeCases = FALSE){

    if(length(object) == 0) stop("object is NULL")

    object <- subRflomicsMAE(object, omicNames)


    if (is.null(object)) return(NULL)

    if (!raw) {
        for (SE.name in names(object)) {
            object[[SE.name]] <- object[[SE.name]][, getSelectedSamples(object, SE.name = SE.name)]
        }
    }

    if (completeCases) {
        object <- object[,complete.cases(object),]
    }

    Groups <- getDesignMat(object)

    nb_entities <-
      lapply(names(object), function(SE){
        dim(object[[SE]])[1] }) %>%
      unlist()
    names(nb_entities) <- names(object)

    data <-
      data.frame(nb_entities = nb_entities,
                 assay = names(nb_entities)) %>%
      full_join(data.frame(sampleMap(object)), by="assay") %>%
      mutate(y.axis = paste0(assay, "\n", "n=", nb_entities)) %>% arrange(primary)

    data$primary <- factor(data$primary, levels = levels(Groups$samples))

    nb_entities_ord <-
      select(data, y.axis, nb_entities) %>%
      unique() %>%
      arrange(desc(nb_entities))
    nb_entities_ord$nb_entities <- log(nb_entities_ord$nb_entities)
    tmp.vec <- c(0)
    breaks  <- vector()
    for(i in seq_len(length(nb_entities_ord$nb_entities))){
      tmp.vec[i+1] <- tmp.vec[i] + nb_entities_ord$nb_entities[i]
      breaks[i] <- tmp.vec[i] + nb_entities_ord$nb_entities[i]/2
    }

    switch (
      as.character(realSize),
      "TRUE"  = {
        p <- ggplot(data, aes(x=primary, y=log(nb_entities))) +
          geom_col(aes(fill = y.axis)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none",
                axis.text.y = element_text(hjust = 0)) +
          labs(x=paste0("Samples (k=", length(unique(sampleMap(object)$primary)), ")"), y="") +
          scale_y_continuous(breaks = (breaks), labels = nb_entities_ord$y.axis)

      },
      "FALSE" = {
        p <- ggplot(data, aes(x=primary, y=y.axis)) +
          geom_tile(aes(fill = y.axis), colour = "grey50") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.ticks = element_blank(), legend.position="none",
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(x=paste0("Samples (k=", length(unique(sampleMap(object)$primary)), ")"), y="")
      }
    )
    return(p)
  })

## ---- plotConditionsOverview ----

#' @rdname RflomicsMAE-class
#' @aliases plotConditionsOverview,RflomicsMAE-method
#' @name plotConditionsOverview
#' @section Plots:
#' \itemize{
#'    \item plotConditionsOverview:
#' A complete design and at least one biological and one batch factors are
#' required for using RFLOMICS workflow.}
#' @param omicNames a vector with dataset names
#' @importFrom purrr reduce
#' @exportMethod plotConditionsOverview
#' @examples
#' # See createRflomicsMAE for an example that includes plotConditionsOverview
setMethod(
  f         = "plotConditionsOverview",
  signature = "RflomicsMAE",
  definition <- function(object,
                         omicNames = NULL){

    # check presence of bio factors
    if (!length(getBioFactors(object)) %in% seq_len(3)){
      stop("No bio factor! or nbr of bio factors exceed 3!") }
    if (!length(getBatchFactors(object)) %in% c(1,2)){
      stop("No replicates found!") }

    BioFact <- getBioFactors(object)
    coldata <- getDesignMat(object) %>%
      mutate(samples=rownames(.))
    #coldata <- tibble::as_tibble(coldata)
    coldata <- sampleMap(object) %>% as.data.frame() %>%
      left_join(coldata, by = c("primary" = "samples"))

    all_combin_cond <- lapply(BioFact, function(x){
      df <- unique(coldata[x])
      rownames(df) <- seq_len(nrow(df))
      return(df)
    }) %>% reduce(merge)

    counts <- coldata %>% select(assay, all_of(BioFact)) %>%
      unique() %>%
      group_by_at(BioFact) %>% count(name = "Count") %>%
      right_join(all_combin_cond, by = BioFact) %>%
      mutate_at(.vars = "Count", .funs = function(x) {
        if_else(is.na(x), 0, x)})


    counts <- counts %>%
      mutate(status = if_else(Count == length(object) , "all_data",
                              if_else(Count == 0 , "no_data", "some_data")))

    #list of factor names
    factors <- names(counts)[seq_len(dim(counts)[2]-2)]

    col.panel <- c("all_data", "some_data", "no_data")
    names(col.panel) <- c("#00BA38", "orange", "red")

    col.panel.u <- col.panel[col.panel %in% unique(counts$status)]

    switch (
      length(factors),
      "1" = { p <- ggplot(counts, aes(x = !!sym(factors[1]), y = 1)) +
        theme(axis.text.y = element_blank()) + ylab("")
      },
      "2" = { p <- ggplot(counts, aes(x = !!sym(factors[1]), y = !!sym(factors[2])))
      },
      "3" = {
        #get factor with min conditions -> to select for "facet_grid"
        factors.l <-
          lapply(factors, function(x){
            length(unique(counts[[x]])) }) %>% unlist()
        names(factors.l) <- factors
        factor.min <- names(factors.l[factors.l == min(factors.l)][1])

        factors <- factors[factors != factor.min]

        #add column to rename facet_grid
        counts <- counts %>%
          mutate(grid = paste0(factor.min, "=",get(factor.min)))

        p <- ggplot(counts , aes(x = !!sym(factors[1]), y = !!sym(factors[2]))) +
          facet_grid(grid~.)
      }
    )

    p <- p + geom_tile(aes(fill = status),
                       color = "white", linewidth = 1,
                       width = 1, height = 1)  +
      geom_text(aes(label = Count)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x=element_text(angle=90, hjust=1),
            legend.position = "none")
    return(p)

  })
