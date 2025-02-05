### ============================================================================
### [01_Load_Data] RflomicsMAE/SE constructors, functions, internal functions
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif

# ----  GLOBAL IMPORT & EXPORT ----
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
#' @param omicsData list of omics dataset: named list of data.frame, matrix or
#' \link{SummarizedExperiment} objects, or \link{MultiAssayExperiment} object.
#' @param omicsNames Vector of dataset names that we want to analyze.
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
  if (is.null(omicsData))
    stop("the omicsData arg is mandatory.")

  if (length(omicsData) == 0)
      stop("the omicsData arg is mandatory.")
  
  if(is.null(names(omicsData)))
    stop("The omicsData must be named.")

  ### => check class of omicsData
  omicsDataClass <- class(omicsData)
  if (!omicsDataClass %in% c("list", "MultiAssayExperiment"))
    stop("The omicsData argument must be of class list or MultiAssayExperiment.")

  ## => omicsNames
  if (is.null(omicsNames)){
    omicsNames <- names(omicsData)
  }
  else 
    if(any(!omicsNames %in% names(omicsData)))
      stop("the omicsNames values must match the names of omicsData object.")

  omicsNames <-
    str_replace_all(string = omicsNames, pattern = "[# /-]", replacement = "")
  if (isTRUE(any(duplicated(omicsNames))))
    stop("presence of duplicates in the omicsNames")

  ## => omicsData to list of data.frames
  omicsData.df <- list()
  for(dataName in omicsNames){
    
    omicsData.df[[dataName]] <-
      switch(
        class(omicsData[[dataName]])[1],
        "data.frame"           = omicsData[[dataName]],
        "matrix"               = as.data.frame(omicsData[[dataName]]),
        "SummarizedExperiment" = {
          tmp <- as.data.frame(assay(omicsData[[dataName]]))
          
          if(is(omicsData , "MultiAssayExperiment"))
            colnames(tmp) <- 
              sampleMap(omicsData)[sampleMap(omicsData)$assay == dataName,]$primary
          tmp
          },
        stop("The components of the omicsData object must be of class list, ",
             "SummarizedExperiment, or matrix.")
      )
    
    colnames(omicsData.df[[dataName]]) <- 
      str_replace_all(string = colnames(omicsData.df[[dataName]]), 
                      pattern = "[# /-]", replacement = "")
  }
  
  ## => omicsTypes
  if(is.null(omicsTypes))
    stop("the list of omicsTypes is mandatory.")

  if(length(omicsData.df) != length(omicsTypes))
    stop("the number of omicsData matrix must match the number of omicsTypes")

  if(isTRUE(any(!unique(omicsTypes) %in% c("RNAseq","metabolomics","proteomics"))))
    stop("omicsTypes must be part of RNAseq, metabolomics, or proteomics.")

  names(omicsTypes) <- omicsNames

  ## any(!is.na(A[A < 0])) valuer negative

  ## => check NA
  for(dataName in names(omicsData.df)){

    # replace NA by 0
    omicsData.df[[dataName]][is.na(omicsData.df[[dataName]])] <- 0
    rowNames <- row.names(omicsData.df[[dataName]])
    omicsData.df[[dataName]] <-
      as.data.frame(lapply(omicsData.df[[dataName]], as.numeric))
    row.names(omicsData.df[[dataName]]) <- rowNames

    # check negative values
    if(any(!is.na(omicsData.df[[dataName]][omicsData.df[[dataName]] < 0])))
      stop("The ",dataName, " data contains negative values")
  }

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

  metaFactors <- setdiff(colnames(ExpDesign), factorInfo$factorName)
  if(length(metaFactors) != 0){
    factorInfo <- data.frame(
      factorName = c(factorInfo$factorName, metaFactors),
      factorType = c(factorInfo$factorType, rep("Meta", length(metaFactors))))
  }
  
  factorBio   <- filter(factorInfo, factorType == "Bio")$factorName
  factorBatch <- filter(factorInfo, factorType == "batch")$factorName

  ## set ref and levels to ExpDesign
  # refList <- vector()
  for (i in 1:nrow(factorInfo)){
    
    if(factorInfo$factorType[i] == "Meta") next
    
    ExpDesign[[factorInfo[i,]$factorName]] <- 
      str_replace_all(string = ExpDesign[[factorInfo[i,]$factorName]], 
                      pattern = "[*# -/]", replacement = "")
    
    # set level
    if (!is.null(factorInfo$factorLevels)){
      levels <- str_split(factorInfo[i,]$factorLevels, ",") |>
        unlist() %>% str_remove(" ")
      
      if(any(!levels %in% ExpDesign[[factorInfo[i,]$factorName]]))
        stop("The factor levels: ", factorInfo[i,]$factorLevels, " don't exist") 
    }
    else{
      levels <- unique(ExpDesign[[factorInfo[i,]$factorName]])
    }
    ExpDesign[[factorInfo[i,]$factorName]] <-
      factor(ExpDesign[[factorInfo[i,]$factorName]], levels = levels)
    
    # set ref
    ref <- NULL
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
  }

  ## consctuct ExpDesign object
  #names(refList)  <- factorInfo$factorName
  typeList <- factorInfo$factorType
  names(typeList) <- factorInfo$factorName

  # Create the List.Factors list with the choosen level of
  # reference for each factor
  #names(typeList) <- names(ExpDesign)

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

    RflomicsSE <- 
      createRflomicsSE(
        omicData  = omicsData.df[[data]], 
        omicType  = omicType, 
        ExpDesign = ExpDesign, 
        design    = typeList)

    #### run PCA for raw count
    SummarizedExperimentList[[data]] <- runOmicsPCA(RflomicsSE, raw = TRUE)

    # metadata for sampleMap for RflomicsMAE
    listmap[[data]] <- data.frame(
      primary = as.vector(colData(SummarizedExperimentList[[data]])$samples),
      colname = as.vector(colData(SummarizedExperimentList[[data]])$samples),
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
                        projectName    = character(),
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

  # Sys.setlocale('LC_TIME', 'C') # change for english
  # date <- format(Sys.time(), '%d %B %Y - %H:%M')
  # Sys.setlocale('LC_TIME') # return to default system time

  metadata <- list(
    "omicList"            = omicList,
    "projectName"         = projectName,
    "design"              = design,
    "IntegrationAnalysis" = IntegrationAnalysis,
    "date"                = Sys.Date(),
    "sessionInfo"         = .writeSessionInfo(),
    "rflomicsVersion"     = packageVersion('RFLOMICS')
  )

  rflomicsMAE <- new("RflomicsMAE", metadata = metadata)
  for(slot in c("ExperimentList", "colData", "sampleMap", "drops")) {
    slot(rflomicsMAE, slot) <- slot(MAE, slot)
  }

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
#' @importFrom dplyr select
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
  matrixOmics <- as.matrix(omicData)
  # nbr of genes with 0 count
  genes_flt0  <- rownames(matrixOmics[rowSums(matrixOmics) <= 0, ])
  # remove 0 count
  matrix.filt  <- matrixOmics[rowSums(matrixOmics)  > 0, ]

  # check if transcriptomics, count matrixOmics
  if (omicType == "RNAseq" &&
      !is.integer(matrix.filt) &&
      !identical(matrix.filt, floor(matrix.filt))) {
    stop("OmicsType is RNAseq, expects counts. The omicData is not counts data.")
  }
  
  # create SE object
  colData   <- mutate(ExpDesign, samples = row.names(ExpDesign)) |>
    filter(samples %in% sample.intersect) |>
    unite("groups", all_of(factorBio), sep = "_", remove = FALSE)

  for (factor in c(factorBio, factorBatch)){

    F.levels <- levels(colData[[factor]])
    colData[[factor]] <- factor(colData[[factor]],
                                levels = intersect(F.levels, unique(colData[[factor]])))
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
    Contrasts.Sel = data.frame())

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

  metadata(rflomicsSE) <-
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
                      "RNAseq"       =  list("variableName" = "transcript",
                                             "valueType" = "counts"),
                      "proteomics"   =  list("variableName" = "protein",
                                             "valueType" = "XIC"),
                      "metabolomics" =  list("variableName" = "metabolite",
                                             "valueType" = "XIC")
  )

  return(valReturn)

}

## ---- generateExampleData ----
#' generateExampleData
#'
#' @return A list
#' @keywords internal
#' @noRd
.generateEcoseedExampleData <- function(){

  #load(paste0(system.file(package = "RFLOMICS"),"/data/ecoseed.df.rda"))
  data("ecoseed.df", package = "RFLOMICS", envir = environment())

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