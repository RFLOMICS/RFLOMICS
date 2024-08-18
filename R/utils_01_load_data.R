### ============================================================================
### [01_Load_Data] RflomicsMAE/SE constructors, functions, internal functions
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif

# ----  GLOBAL IMPORT & EXPORT ----
#' @importFrom dplyr mutate across if_else filter select
#' @importFrom stringr str_replace_all str_remove_all str_remove fixed str_split
#' @importFrom vroom vroom
#' @importFrom purrr reduce
#' @importFrom tidyr unite
#' @importFrom tidyselect where all_of
#' @importFrom magrittr "%>%" 

# ----  RflomicsMAE CLASS ----
## ---- createRflomicsMAE: create RflomicsMAE object from loaded data ----
#' @title RflomicsMAE class constructor
#' @description This function is a constructor for the class
#'  \link{RflomicsMAE-class}.
#' It initializes an object of class \link{RflomicsMAE-class}
#' from a list of omics datasets, a vector of dataset names, 
#' a vector of omics types, and an experimental design.
#' @param projectName Project name
#' @param omicsData list of omics dataset, list of data.frame, matrix or 
#' \link{SummarizedExperiment} objects, or \link{MultiAssayExperiment} object.
#' @param omicsNames vector of dataset names.
#' @param omicsTypes vector of dataset types: 
#' "RNAseq", "metabolomics", "proteomics"
#' @param ExpDesign a data.frame which describes the experimental design.
#' @param factorInfo a data.frame describing the experimental factors.
#' \itemize{
#' \item factorName: factor names 
#' \item factorRef: factor references 
#' \item factorType: factor type : "Bio", "batch", "Meta" 
#' \item factorLevels: levels of each factor with "," separation. 
#' }
#' @return An object of class \link{RflomicsMAE-class}
#' @seealso \link{RflomicsMAE-class}
#' @name createRflomicsMAE
#' @rdname createRflomicsMAE
#' @export
#' @example inst/examples/loadData.R
createRflomicsMAE <- function(projectName = NULL, 
                              omicsData   = NULL, 
                              omicsNames  = NULL, 
                              omicsTypes  = NULL, 
                              ExpDesign   = NULL, 
                              factorInfo  = NULL){
  
  # check arg
  ## => projectName
  if(is.null(projectName)) stop("projectName is mandatory.")
  projectName <- 
    str_replace_all(string = projectName, pattern = "[# /-]", replacement = "")
  
  ## => omicsData
  if(is.null(omicsData) || length(omicsData) == 0) 
    stop("the omicsData arg is mandatory.")
  
  ### => check class of omicsData
  omicsDataClass <- class(omicsData)
  if(!omicsDataClass %in% c("list", "MultiAssayExperiment"))
    stop("The omicsData argument must be of class list or MultiAssayExperiment.")

  ## => omicsNames
  if(is.null(omicsNames)){
    if(is.null(names(omicsData))) 
      stop("The omicsData must be named; otherwise, you need to use the omicsNames argument.")
    omicsNames <- names(omicsData)
  }
  else if(length(omicsNames) != length(omicsData))
    stop("the number of omicsData matrix must match the number of omicsNames.")
  
  omicsNames <- 
    str_replace_all(string = omicsNames, pattern = "[# /-]", replacement = "")
  if (isTRUE(any(duplicated(omicsNames))))
    stop("presence of duplicates in the omicsNames")
  
  names(omicsData) <- omicsNames
  
  ## => omicsData to list of data.frames
  omicsData.df <- list()
  for(dataName in names(omicsData)){
    
    omicsData.df[[dataName]] <-
      switch (
        class(omicsData[[dataName]])[1],
        "data.frame" = omicsData[[dataName]],
        "matrix" = as.data.frame(omicsData[[dataName]]),
        "SummarizedExperiment" = as.data.frame(assay(omicsData[[dataName]])),
        stop("The components of the omicsData object must be of class list, ", 
             "SummarizedExperiment, or matrix.")
      )
  }
  
  ## => omicsTypes
  if(is.null(omicsTypes))
    stop("the list of omicsTypes is mandatory.")
  
  if(length(omicsData.df) != length(omicsTypes))
    stop("the number of omicsData matrix must match the number of omicsTypes")
  
  if(isTRUE(any(!unique(omicsTypes) %in% c("RNAseq","metabolomics","proteomics"))))
    stop("omicsTypes must be part of RNAseq, metabolomics, or proteomics.")
  
  names(omicsTypes) <- omicsNames
  
  
  ## => ExpDesign
  if (is.null(ExpDesign)){
    ExpDesign <- 
      switch (omicsDataClass,
              "MultiAssayExperiment" = as.data.frame(colData(omicsData)),
              stop("the ExpDesign is mandatory.")
      )
  }
  if (nrow(ExpDesign) == 0 || ncol(ExpDesign) == 0)
    stop("the ExpDesign is mandatory.")
  
  designRownames <- 
    str_replace_all(string = rownames(ExpDesign), pattern = "[*# -/]", replacement = "")
  if (isTRUE(any(duplicated(designRownames))))
    stop("presence of duplicates in the ExpDesign colnames")
  
  rownames(ExpDesign) <- designRownames
  
  ## => factorInfo
  if (is.null(factorInfo)) stop("data.frame factorInfo is mandatory.")
  
  ### => factorInfo$factorName
  if (is.null(factorInfo$factorName)) stop("factorInfo$factorName is mandatory")
  
  if (any(!factorInfo$factorName %in% colnames(ExpDesign)))
    stop("factorInfo$factorName don't match ExpDesign colnames")
  
  ### => factorInfo$factorType
  if (is.null(factorInfo$factorType)) stop("factorInfo$factorType is mandatory.")
  
  if (any(!unique(factorInfo$factorType) %in% c("batch", "Bio", "Meta")))
    stop("factorInfo$factorType must be part of batch, Bio or Meta")
  
  factorBio   <- filter(factorInfo, factorType == "Bio")$factorName
  factorBatch <- filter(factorInfo, factorType == "batch")$factorName
  
  ## set ref and levels to ExpDesign
  refList <- vector()
  for (i in 1:nrow(factorInfo)){
    
    # set ref 
    if (!is.null(factorInfo$factorRef)){
      
      if(!factorInfo[i,]$factorRef %in% ExpDesign[[factorInfo[i,]$factorName]])
        stop("The factor ref: ", factorInfo[i,]$factorRef, " don't exist")
      
      ref <- factorInfo[i,]$factorRef
    }
    else{
      ref <- as.character(unique(ExpDesign[[factorInfo[i,]$factorName]])[1])
    }
    #ExpDesign <- ExpDesign[order(row.names(ExpDesign)), ]
    ExpDesign[[factorInfo[i,]$factorName]] <- 
      relevel(as.factor(ExpDesign[[factorInfo[i,]$factorName]]), ref=ref)
    
    refList <- c(refList, ref)
    
    # set level
    if (!is.null(factorInfo$factorLevels)){
      
      levels <- str_split(factorInfo[i,]$factorLevels, ",") |> 
        unlist() %>%
        str_remove(" ")
      if(any(!levels %in% ExpDesign[[factorInfo[i,]$factorName]]))
        stop("The factor levels: ", factorInfo[i,]$factorLevels, " don't exist")
      
      
      ExpDesign[[factorInfo[i,]$factorName]] <- 
        factor(ExpDesign[[factorInfo[i,]$factorName]], levels = levels)
    }
  }
  
  ## consctuct ExpDesign object
  names(refList)  <- factorInfo$factorName
  typeList <- factorInfo$factorType; names(typeList) <- factorInfo$factorName
  
  # Create the List.Factors list with the choosen level of 
  # reference for each factor
  names(typeList) <- names(ExpDesign)

  Design <- list(Factors.Type  = typeList, 
                 Model.formula = vector(), 
                 Contrasts.Sel = data.frame())

  ExpDesign   <- mutate(ExpDesign, samples=row.names(ExpDesign)) |>
    unite("groups", all_of(factorBio), sep = "_", remove = FALSE)
  
  order_levels <- 
    with(ExpDesign, do.call(order, ExpDesign[c(factorBio, factorBatch)]))
  ExpDesign$samples <- 
    factor(ExpDesign$samples, levels = unique(ExpDesign$samples[order_levels]))
  ExpDesign$groups <- 
    factor(ExpDesign$groups,  levels = unique(ExpDesign$groups[order_levels]))
  
  ## create SE object of each dataset
  SummarizedExperimentList <- list()
  listmap  <- list()
  omicList <- list()
  k <- 0
  
  for(data in omicsNames){
    
    k <- k+1
    omicType <- omicsTypes[data]
    
    RflomicsSE <- createRflomicsSE(omicsData.df[[data]], omicType, ExpDesign, typeList)
    
    #### run PCA for raw count
    SummarizedExperimentList[[data]] <- runOmicsPCA(RflomicsSE, raw = TRUE)

    # metadata for sampleMap for RflomicsMAE
    listmap[[data]] <- data.frame(
      primary = as.vector(SummarizedExperimentList[[data]]@colData$samples),
      colname = as.vector(SummarizedExperimentList[[data]]@colData$samples),
      stringsAsFactors = FALSE)
    
    colnames <- c(names(omicList[[omicType]]), k)
    omicList[[omicType]] <- c(omicList[[omicType]] ,data)
    names(omicList[[omicType]]) <- colnames
  }
  
  RfMAE <- RflomicsMAE(experiments = SummarizedExperimentList,
                       colData     = ExpDesign,
                       sampleMap   = listmap,
                       omicList    = omicList, 
                       projectName = projectName, 
                       design      = Design)
  
  # tag as raw data (le temps de trouver une solution pour 
  # ne pas faire co-exister les raw et les process)
  # names(RfMAE) <- paste(names(RfMAE), "raw", sep = ".")
  
  return(RfMAE)
}


## ---- RflomicsMAE: construct RflomicsMAE object ----
#' @title RflomicsMAE 
#' @description
#'  \link{RflomicsMAE-class} constructor.
#' @param experiments same as in MultiAssayExperiments.
#'  A list or experimentList of all combined experiments.
#' @param colData same as in MultiAssayExperiments. 
#' A DataFrame or data.frame of characteristics for all biological units
#' @param sampleMap same as in MultiAssayExperiments.
#' 	A DataFrame or data.frame of assay names, sample identifiers, 
#' 	and colname samples
#' @param omicList list of omics names in experiment lists.
#' @param projectName name of the project. This will be the name of the 
#' generated archive and report.
#' @param design experimental design. 
#' @param IntegrationAnalysis a list. Integration analysis results. 
#' Default is an empty list.
#' @seealso \link{MultiAssayExperiment}
#' @importFrom MultiAssayExperiment MultiAssayExperiment listToMap
#' @return An object of class \link{RflomicsMAE-class}
#' @name RflomicsMAE
#' @keywords internal
#' @noRd
RflomicsMAE <- function(experiments = ExperimentList(), 
                        colData     = S4Vectors::DataFrame(), 
                        sampleMap   = S4Vectors::DataFrame(
                          assay = factor(), 
                          primary = character(), 
                          colname = character()),
                        omicList       = list(),
                        projectName    = NULL,
                        design         = list(),
                        IntegrationAnalysis = list()){
  
  MAE <- NULL
  if (is(sampleMap, "DFrame") || is.data.frame(sampleMap)) {
    MAE <- MultiAssayExperiment(experiments, colData, sampleMap)
  } else if (is.list(sampleMap)) {
    MAE <- MultiAssayExperiment(experiments, colData, 
                                listToMap(sampleMap))
  } else{
    MAE <- MultiAssayExperiment(experiments, colData, sampleMap)
  }
  
  rflomicsMAE <- new("RflomicsMAE")
  for(slot in c("ExperimentList", "colData", "sampleMap", "drops")) {
    slot(rflomicsMAE, slot) <- slot(MAE, slot)
  }
  
  # Sys.setlocale('LC_TIME', 'C') # change for english
  # date <- format(Sys.time(), '%d %B %Y - %H:%M')
  # Sys.setlocale('LC_TIME') # return to default system time
  
  rflomicsMAE@metadata$omicList <- omicList
  rflomicsMAE@metadata$projectName <- projectName
  rflomicsMAE@metadata$design <- design
  rflomicsMAE@metadata$IntegrationAnalysis <- IntegrationAnalysis
  
  return(rflomicsMAE)
}

# ----  RflomicsSE CLASS ----
## ---- createRflomicsSE: create RflomicsSE object from loaded data ----
#' @title createRflomicsSE 
#' @description This function initializes an object of 
#' class \link{RflomicsSE} from a list of omics data.
#' @param omicsData omics dataset.
#' @return An object of class \link{RflomicsSE}
#' @name createRflomicsSE
#' @seealso \link{RflomicsSE-class}
#' @keywords internal
#' @noRd
createRflomicsSE <- function(omicData, omicType, ExpDesign, design){
  
  factorBio   <- names(design[design == "Bio"])
  factorBatch <- names(design[design == "batch"])
  
  # check overlap between design and data
  sample.intersect <- intersect(row.names(ExpDesign), colnames(omicData))
  if(length(sample.intersect) == 0) 
    stop("samples in omics data should match the names in experimental design")
  
  # select abundance from design table and reorder
  omicData <- select(omicData, all_of(sample.intersect))
  
  # remove row with sum == 0
  matrix <- as.matrix(omicData)
  # nbr of genes with 0 count
  genes_flt0  <- rownames(matrix[rowSums(matrix) <= 0, ])
  # remove 0 count
  matrix.filt  <- matrix[rowSums(matrix)  > 0, ]
  
  # create SE object
  colData   <- mutate(ExpDesign, samples=row.names(ExpDesign)) |>
    filter(samples %in% sample.intersect) |> 
    unite("groups", all_of(factorBio), sep = "_", remove = FALSE)
  
  for (factor in c(factorBio, factorBatch)){
    
    F.levels <- levels(colData[[factor]])
    colData[[factor]] <- factor(colData[[factor]], levels = intersect(F.levels, unique(colData[[factor]])))
  }
  
  order_levels <- with(colData, do.call(order, colData[c(factorBio, factorBatch)]))
  colData$samples <- factor(colData$samples, levels = unique(colData$samples[order_levels]))
  colData$groups  <- factor(colData$groups,  levels = unique(colData$groups[order_levels]))
  
  dataProcessing <- 
    list(rowSumsZero      = genes_flt0,
         selectedSamples  = colData$samples, 
         featureFiltering = list(), 
         Normalization    = list(), 
         Transformation   = list(),
         log = NULL)
  
  Design <- list(
    factorType = design[intersect(names(design), names(colData))],
    Model.formula = vector(),
    Contrasts.Sel = data.frame(),
    ExpDesign = DataFrame(colData))
  
  rflomicsSE <- 
    RflomicsSE(assays         = matrix.filt, 
               colData        = DataFrame(colData), 
               omicType       = omicType,
               design         = Design,
               DataProcessing = dataProcessing)
  
  return(rflomicsSE)
}

## ---- RflomicsSE: construct RflomicsSE object ----
#' @title RflomicsSE
#' @description
#' A constructor for the class \link{RflomicsSE}.
#' @param name description
#' @seealso SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @return An object of class \link{RflomicsSE}
#' @name RflomicsSE
#' @seealso \link{RflomicsSE-class}
#' @keywords internal
#' @noRd
RflomicsSE <- function(assays = NULL, colData = NULL, 
                       omicType = NULL,
                       design   = list() , 
                       DataProcessing = list() , 
                       PCAlist        = list() , 
                       DiffExpAnal    = list() ,      
                       CoExpAnal      = list() , 
                       DiffExpEnrichAnal = list() ,  
                       CoExpEnrichAnal   = list()){
  
  
  if(!is.null(assays))
    assays <- SimpleList(abundance = as.matrix(assays))
  
  SE <- SummarizedExperiment( assays  = assays, 
                              colData = colData)
  
  rflomicsSE <- new("RflomicsSE")
  for(slot in c("colData","assays","NAMES","elementMetadata")) {
    slot(rflomicsSE, slot) <- slot(SE, slot)
  }
  
  rflomicsSE@metadata <- 
    list("omicType" = omicType , 
         "design"   = design , 
         "DataProcessing" = DataProcessing , 
         "PCAlist"        = PCAlist , 
         "DiffExpAnal"    = DiffExpAnal ,      
         "CoExpAnal"      = CoExpAnal , 
         "DiffExpEnrichAnal" = DiffExpEnrichAnal ,  
         "CoExpEnrichAnal"   = CoExpEnrichAnal)
  
  return(rflomicsSE)
}

# ----  READ INPUR FILES ----
## ---- read_exp_design: read experimental design file ----
#' @title Read Experimental Design
#' @param file path to experimental design file
#' @return data.frame
#' @importFrom tidyselect where
#' @importFrom purrr reduce
#' @noRd
#' @keywords internal
readExpDesign <- function(file){
  
  if (missing(file)) {
    stop('Please provide a file path')
  }
  
  if(!file.exists(file))
  {
    stop(file, " don't exist!")
    return(NULL)
  }
  
  # read design and remove special characters
  # remove "_" from modality and factor names
  data <- vroom(file, delim = "\t", show_col_types = FALSE) %>%
    mutate(across(.cols = where(is.character), 
                  ~str_remove_all(.x, pattern = "[.,;:#@!?()$%&<>|=+-/]"))) %>%
    mutate(across(.cols = where(is.character), 
                  ~str_remove_all(.x, pattern = "[\\]\\[\'\"\ ]"))) %>%
    mutate(across(.cols = where(is.character), 
                  ~str_remove_all(.x, pattern = fixed("\\")))) %>% 
    mutate(across(.cols = c(-1), ~str_remove_all(.x, pattern = fixed("_")))) %>% 
    mutate(across(.cols = where(is.character), as.factor)) 
  
  names(data)  <- str_remove_all(string = names(data), pattern = "[.,;:#@!?()$%&<>|=+-/\\]\\[\'\"\ _]") %>%
    str_remove_all(., pattern = fixed("\\"))
  
  # check if there is duplication in sample names
  sample.dup <- as.vector(data[which(table(data[1]) > 1),1])[[1]]
  
  if (length(sample.dup) != 0) {
    
    stop("Duplicated sample names: ", paste0(sample.dup, collapse = ","))
  }
  
  # check if there is duplication in factor names
  # factor.dup <- as.vector(data[which(table(names(data[-1])) > 1),1])[[1]]
  factor.dup <- names(data[-1])[duplicated(names(data[-1]))]
  if (length(factor.dup) != 0) {
    stop("Duplicated factor name: ", paste0(factor.dup, collapse = ","))
  }
  
  # check if same name of moralities are used in diff factor
  mod.list <- sapply(names(data[-1]), function(x){ 
    unique(data[-1][[x]])
  }) %>% reduce(c)
  
  mod.dup <- mod.list[duplicated(mod.list)]
  if(length(mod.dup) != 0) {
    
    stop("Modality used in more than one factor: ", 
         paste0(mod.dup[1:10], collapse = ", "))
  }
  
  # warning if number of factors exceed n = 10
  n <- 10
  if (dim(data)[2]-1 >= n){
    
    data <- data[, 1:n]
    warning("Large number of columns! only the first ", n," will be displayed")
  }
  
  # check nbr of modality of the 5th fist columns
  index <- sapply(names(data[-1]), function(x){ if(length(unique(data[-1][[x]]))>n){ FALSE }else{ TRUE } })
  F.mod <- names(data[-1])[index]
  
  ratio <- length(F.mod)/length(names(data[-1]))
  
  if(ratio != 1)
  {
    warning("The select input contains a large number of options")
  }
  
  data            <- data.frame(data) 
  row.names(data) <- data[,1]
  data            <- data[,-1]
  return(data)
}


## ---- readOmicsData: read dataset matrix file ----
#' @title Read omics data 
#' @param file omics data matrix
#' @return data.frame
#' @importFrom vroom vroom
#' @noRd
#' @keywords internal
readOmicsData <- function(file){
  
  if(!file.exists(file))
  {
    stop(file, " don't exist!")
  }
  
  # read omics data and remove special characters
  data <- vroom(file, delim = "\t", show_col_types = FALSE)
  names(data)  <- str_remove_all(string = names(data), 
                                 pattern = "[.,;:#@!?()$%&<>|=+-/\\]\\[\'\"\ ]") %>%
    str_remove_all(., pattern = fixed("\\"))
  
  # check if there is duplication in sample names
  sample.dup <- as.vector(data[which(table(names(data[-1])) > 1),1])[[1]]
  
  if (length(sample.dup) !=0){
    
    stop("Duplicated sample names: ", paste0(sample.dup, collapse = ","))
  }
  
  # check if there is duplication in factor names
  entity.dup <- as.vector(data[which(table(data[1]) > 1),1])[[1]]
  
  if (length(entity.dup) !=0){
    
    stop("Duplicated feature names: ", paste0(entity.dup, collapse = ","))
  }
  
  data            <- data.frame(data) 
  row.names(data) <- data[,1]
  data            <- data[,-1]
  return(data)
}

## ---- checkSpecialCharacters: read dataset matrix file ----

# ----  INTERNAL FUNCTIONS ----
## ---- omicsDic: get variable name and type from omicstype ----
#' @title Omics Dictionary
#'
#' @param object a MAE object or a SE object (produced by Flomics). 
#' Expect to find a omicsType somewhere.
#' @param SE.name if object is a MAE, expect to find the experiment 
#' name from which the omics info has to be retrieved.
#' @return list of two elements: variableName and valueType.
#' @noRd
#' @keywords internal
.omicsDic <- function(object, SE.name = NULL){
  
  if (!is(object, "RflomicsSE") && !is(object, "RflomicsMAE")) {
    stop("Object must be a RflomicsSE or a RflomicsMAE, not a ",
         class(object))
  }
  
  if (is(object, "RflomicsMAE")) {
    if (missing(SE.name)) {
      stop("Please provide an Experiment name (SE.name).")
    }
    
    object <- object[[SE.name]]
  }
  
  omicsType <- getOmicsTypes(object)
  
  valReturn <- switch(omicsType,
                      "RNAseq"       =  list("variableName" = "transcripts",
                                             "valueType" = "counts"),
                      "proteomics"   =  list("variableName" = "proteins",
                                             "valueType" = "XIC"),
                      "metabolomics" =  list("variableName" = "metabolites",
                                             "valueType" = "XIC")
  )
  
  return(valReturn)
  
}

#' @title Omics Dictionary
#'
#' @param object a MAE object or a SE object (produced by Flomics). 
#' Expect to find a omicsType somewhere.
#' @param SE.name if object is a MAE, expect to find the experiment 
#' name from which the omics info has to be retrieved.
#' @return list of two elements: variableName and valueType.
#' @noRd
#' @keywords internal
omicsDic <- function(object, SE.name = NULL){
  
  if (!is(object, "RflomicsSE") && !is(object, "RflomicsMAE")) {
    stop("Object must be a RflomicsSE or a RflomicsMAE, not a ",
         class(object))
  }
  
  if (is(object, "RflomicsMAE")) {
    if (missing(SE.name)) {
      stop("Please provide an Experiment name (SE.name).")
    }
    
    object <- object[[SE.name]]
  }
  
  omicsType <- getOmicsTypes(object)
  
  valReturn <- switch(omicsType,
                      "RNAseq"       =  list("variableName" = "transcripts",
                                             "valueType" = "counts"),
                      "proteomics"   =  list("variableName" = "proteins",
                                             "valueType" = "XIC"),
                      "metabolomics" =  list("variableName" = "metabolites",
                                             "valueType" = "XIC")
  )
  
  return(valReturn)
  
}

## ---- checkNA: checks if there are NA/nan in the RflomicsSE assay ----
#' @title checkNA
#'
#' @param object An object of class \link{RflomicsSE}
#' @importFrom MultiAssayExperiment assay
#' @return boolean. if TRUE, NA/nan are detected in dataset matrix.
#' @keywords internal
#' @noRd
#'
.checkNA <- function(object) {
  NA_detect <- ifelse(any(is.na(assay(object))), TRUE, FALSE)
  return(NA_detect)
}

#' @title check_NA
#'
#' @param object An object of class \link{RflomicsSE}
#' @return boolean. if TRUE, NA/nan are detected in the SE::assay.
#' @keywords internal
#' @noRd
#'
check_NA <- function(object) {
  NA_detect <- ifelse(any(is.na(assay(object))), TRUE, FALSE)
  return(NA_detect)
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

# ---- PLOTS
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
#' @importFrom ggplot2 ggplot aes aes_string theme facet_grid labs ylab xlab
#' facet_grid element_blank geom_text scale_fill_manual geom_tile ggtitle
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
  
  switch (length(factors),
          "1" = { p <- ggplot(counts ,aes_string(x = factors[1], y = 1)) + 
            theme(axis.text.y = element_blank()) + ylab("") },
          "2" = { p <- ggplot(counts ,aes_string(x = factors[1], y = factors[2])) },
          "3" = {
            #get factor with min conditions -> to select for "facet_grid"
            factors.l <- lapply(factors, function(x){ length(unique(counts[[x]])) }) %>% unlist()
            names(factors.l) <- factors
            factor.min <- names(factors.l[factors.l == min(factors.l)][1])
            
            factors <- factors[factors != factor.min]
            
            #add column to rename facet_grid
            counts <- counts %>% mutate(grid = paste0(factor.min, "=",get(factor.min)))
            
            p <- ggplot(counts ,aes_string(x = factors[1], y = factors[2])) +
              facet_grid(grid~.) })
  
  p <- p + 
    geom_tile(aes(fill = status), color = "white",
              size = 1, width = 1, height = 1) + 
    geom_text(aes(label = Count)) + 
    scale_fill_manual(values = names(col.panel.u), breaks = col.panel.u) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle(message)
  
  return(p)
}

## ---- generateExampleData ----
#' generateExampleData
#' 
#' @return A list
#' @keywords internal
#' @noRd
.generateEcoseedExampleData <- function(){
  
  load(paste0(system.file(package = "RFLOMICS"),"/data/ecoseed.df.rda"))
  
  factorInfo <- data.frame(
    "factorName"   = c("Repeat", "temperature", "imbibition"),
    "factorRef"    = c("rep1", "Low", "DS"),
    "factorType"   = c("batch", "Bio", "Bio"),
    "factorLevels" = c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI")
  )
  
  ExpDesign <- ecoseed.df$design
  
  ExpDesign[["imbibition"]]  <- 
    factor(ExpDesign[["imbibition"]], 
           unlist(stringr::str_split(
             factorInfo$factorLevels[factorInfo$factorName == "imbibition"], ",")))
  
  ExpDesign[["temperature"]] <- 
    factor(ExpDesign[["temperature"]], 
           unlist(stringr::str_split(
             factorInfo$factorLevels[factorInfo$factorName == "temperature"], ",")))
  
  ExpDesign[["Repeat"]]      <- 
    factor(ExpDesign[["Repeat"]],
           unlist(stringr::str_split(
             factorInfo$factorLevels[factorInfo$factorName == "Repeat"], ",")))
  
  factorRef <- factorInfo$factorRef
  names(factorRef) <- factorInfo$factorName
  
  exampleData <- list(
    projectName   = "Ecoseed",
    ExpDesign     = ExpDesign,
    dF.List.ref   = factorRef,
    dF.Type.dFac  = factorInfo$factorType,
    omicsNames    = c("RNAseq.set1", "proteomics.set3", "metabolomics.set2"),
    omicsTypes    = c("RNAseq", "proteomics", "metabolomics"),
    omicsData     = 
      list("RNAseq.set1"       = ecoseed.df$RNAtest,
           "proteomics.set3"   = ecoseed.df$protetest,
           "metabolomics.set2" = ecoseed.df$metatest)
  )
  
  return(exampleData)
}

#' 