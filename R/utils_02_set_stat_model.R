### ============================================================================
### [02_set_stat_model] function and internal function
### ----------------------------------------------------------------------------
# C. Payasant-Le-roux,
# N. Bessoltane
# D. Charif


#' @importFrom magrittr "%>%"
NULL

# ---- .generateModelFormulae ----
#' @title .generateModelFormulae
#' @description
#'
#' From a vector of character giving the name of the factors of an omics
#' experiment, and their type of effect: biological or batch, it returns
#' all models formulae that can be formulated in association with this factors.
#'  Batch effect factors do not appear in interaction terms with biological
#'  factor. Model formulae stop in second order interaction.
#'
#' @param FacBio a vector of character giving the name of the bio factors.
#' @param FacBatch a vector of character giving the name of the batch factors.
#' @importFrom stats update
#' @return a named list of object of class formula
#' @keywords internal
#' @noRd
#' @examples
#'
#' .generateModelFormulae(FacBio=c("Genotype","Temperature","Environment"),
#'                        FacBatch=c("Replicat"))
#' .generateModelFormulae(FacBio=c("Genotype","Temperature"),
#'                        FacBatch=c("Replicat"))
#' .generateModelFormulae(FacBio=c("Genotype"), FacBatch=c("Replicat"))
#'
#' .generateModelFormulae(FacBio=c("Genotype","Temperature"),
#'                        FacBatch=c("Replicat", "laboratory"))
#'
.generateModelFormulae <- function(FacBio=NULL, FacBatch=NULL){

  # Initialize
  formulae <- list()

  # Verify that nbr of bio factors are between 1 and 3.
  if(!length(FacBio) %in% seq_len(3)) stop(".... !")

  # Verify that nbr of batch factors are between 1 and 2.
  if(!length(FacBatch) %in% c(1,2)) stop(".... !")

  nFac <- length(FacBio)

  # get formulae without interation
  formulae[[1]] <- update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+",paste(FacBio,collapse="+"))),new=~.)

  # get formulae with interation if nbr of FacBio > 1
  if(nFac !=1)
    formulae[[2]] <- update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+","(",paste(FacBio,collapse="+"),")^2")),new=~.)

  formulae <- unlist(formulae)
  names(formulae) <- unlist(as.character(formulae))

  # Sort formulae

  formulae <- formulae[order(unlist(lapply(names(formulae),nchar)),decreasing=TRUE)]

  return(formulae)
}

# ---- contrastName2contrastDir ----

## contrastName to name of contrast directory
# Exemple:
# contrastName
# "(temperatureMedium - temperatureElevated) in imbibitionEI - (temperatureMedium - temperatureElevated) in imbibitionDS"
# contrastDir
# "temperatureMedium-temperatureElevated_in_imbibitionEI_vs_temperatureMedium-temperatureElevated_in_imbibitionDS"

#' Title
#'
#' @param a string: contrastName
#'
#' @return a string: contrastDir
#' @importFrom stringr str_replace_all str_remove_all
#' @keywords internal
#' @noRd
contrastName2contrastDir <- function(contrastName){
  # remplacement des comparaisons centrales
  tmp <- str_replace_all(contrastName,"[:blank:]-[:blank:]\\(","_vs_")
  # remplacement des comparaisons dans les parenthèses
  tmp <- str_replace_all(tmp,"[:blank:]-[:blank:]","-")
  # suppression des parenthèses
  tmp <- stringr::str_remove_all(tmp,c("\\(|\\)"))
  # remplcement des espaces par des _
  tmp <- str_replace_all(tmp,"[:blank:]","_")
  return(tmp)
}

# ---- .getExpressionContrast : function generating contrast expression devlopped by CPL ----

#' get contrast expression
#'
#' @param ExpDesign data.frame with only bio factors
#' @param factorBio vector of bio factors
#' @param modelFormula formula
#' @return a list of dataframe with all contrasts per type
#' @keywords internal
#' @noRd
.getExpressionContrastF <- function(ExpDesign, factorBio=NULL, modelFormula=NULL){

  # ExpDesign
  if(is.null(ExpDesign) || nrow(ExpDesign) == 0) stop("ExpDesign argument is mandatory.")

  # model formula
  if (is.null(modelFormula)) stop("modelFormula argument is mandatory.")
  if (is(modelFormula, "formula")) modelFormula <- paste(as.character(modelFormula), collapse = " ")
  modelFormula <- formula(modelFormula)

  # factorBio
  if (is.null(factorBio)) stop("factorBio argument is mandatory.")
  if (length(intersect(factorBio, names(ExpDesign))) == 0) stop("factorBio and names(ExpDesign) don't overlap!")

  # bio factor list in formulat
  labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")

  FactorBioInDesign <- intersect(labelsIntoDesign, factorBio)
  if (length(FactorBioInDesign) == 0) stop("factorBio and attributes (attr) of modelFormula don't overlap!")

  treatmentFactorsList <- lapply(FactorBioInDesign, function(x){(paste(x, levels(ExpDesign[[x]]), sep=""))})
  names(treatmentFactorsList) <- FactorBioInDesign

  interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)

  listOfContrastsDF <- list()
  # define all simple contrasts pairwise comparisons

  allSimpleContrast_df <- .defineAllSimpleContrasts(treatmentFactorsList)
  # if 1 factor or more than 1 + interaction
  if(length(treatmentFactorsList) == 1 || !isFALSE(interactionPresent)){

    listOfContrastsDF[["simple"]] <- select(allSimpleContrast_df, contrast, contrastName, groupComparison, type)
  }

  # define all simples contrast means
  # exists("allSimpleContrast_df", inherits = FALSE)
  if(length(treatmentFactorsList) != 1){
    allAveragedContrasts_df <- .define_averaged_contrasts(allSimpleContrast_df)
    listOfContrastsDF[["averaged"]] <-  select(allAveragedContrasts_df, contrast, contrastName, groupComparison, type)
  }

  # define all interaction contrasts
  if(length(treatmentFactorsList) != 1){
    if(interactionPresent){
      labelsIntoDesign            <- attr(terms.formula(modelFormula),"term.labels")
      labelOrder                  <- attr(terms.formula(modelFormula), "order")
      twoWayInteractionInDesign   <- labelsIntoDesign[which(labelOrder == 2)]
      groupInteractionToKeep      <- gsub(":", " vs ", twoWayInteractionInDesign)
      allInteractionsContrasts_df <- .defineAllInteractionContrasts(treatmentFactorsList, groupInteractionToKeep)

      listOfContrastsDF[["interaction"]] <- select(allInteractionsContrasts_df, contrast, contrastName, groupComparison, type)
    }
    #allInteractionsContrasts_df <- .defineAllInteractionContrasts(treatmentFactorsList)
    #listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
  }
  # choose the contrasts and rbind data frames of contrasts
  #selectedContrasts <- returnSelectedContrasts(listOfContrastsDF)

  tag <- 0
  for(C_type in names(listOfContrastsDF)){
    listOfContrastsDF[[C_type]]$tag <-
      paste0("H", (tag+1):(nrow(listOfContrastsDF[[C_type]])+tag))
    tag <- nrow(listOfContrastsDF[[C_type]])+tag
  }

  return(listOfContrastsDF)

}


# it is possible to define contrast combinations that are specifically suited
# to a particular experimental design and set of research questions

#Contrasts are used to test whether the levels of an effect are significantly
# different from one another. You can specify a contrast for each factor
# in the model. Contrasts represent linear combinations of the parameters.

# Comparison (or contrast) procedures are used to test more specific
# hypotheses about differences between means comparisons of the two drug
# groups to the control group (  a complex comparison)
# The contrasts are formed by applying a set of weights, called contrast
# coefficients, to the means
# A contrast is a test of the difference between the means of two groups
#  from the ANOVA. There are two categories of contrasts among the groups
#  tested by ANOVA, simple and complex.
# A simple contrast is a test of the difference between any two pairs,
# such as Experimental Group 1 and Control Group 2. A complex contrast
# is a test of the difference between combinations of groups.
# An example of a complex contrast is a test of the difference between
# a subgroup created by combining Experimental Groups 1, 2 and 4 combined,
# and a subgroup created by combining Control Groups 1 and 3

# pairwise comparisons between the groups
# biological_factors_list, interaction_in_model, model_matrix

## functions to define part of simple contrasts
#' define part of simple contrast data frame
#'
#' @param treatmentFactorsList list
#' @param i i
#' @param j j
#'
#' @return a dataframe of contrasts
#' @keywords internal
#' @importFrom utils combn
#' @importFrom data.table setDT setnames
#' @importFrom tidyr unite
#' @importFrom tidyselect all_of
#' @noRd
#' @author Christine Paysant-Le Roux
.define_partOfSimpleContrast_df <- function (treatmentFactorsList, i, j) {

  contrastPart <- fixFactor <- NULL

  comparisonPart <- treatmentFactorsList
  # combn(x,2) generate all combinations of the elements of x taken 2 at a time
  vectorFromCombn <- combn(treatmentFactorsList[[i]],2)[j,]
  comparisonPart[[i]] <- vectorFromCombn
  df_comparisonPart <- expand.grid(comparisonPart)
  setDT(df_comparisonPart)

  # paste all the column of the data table
  df_comparisonPart <- df_comparisonPart %>%
    unite(contrastPart, sep="_", remove=FALSE)

  #colnameFactor_i <- names(df_comparisonPart)[i]
  colnameFactor_i <- names(treatmentFactorsList)[i]

  #df_comparisonPart[, comparisonPart := df_comparisonPart[[colnameFactor_i]]]
  df_comparisonPart <- df_comparisonPart %>%
    mutate(comparisonPart = df_comparisonPart[[colnameFactor_i]])

  if(length(names(treatmentFactorsList)) != 1){
    colnamesToKeep <- setdiff(names(df_comparisonPart),c("contrastPart", "comparisonPart", colnameFactor_i))
    df_comparisonPart <- df_comparisonPart %>%
      unite(fixFactor, all_of(colnamesToKeep), sep = "_", remove = FALSE)
  }else{
    df_comparisonPart <- df_comparisonPart %>%
      mutate(fixFactor= NA)
  }


  colnamesToDelete <- names(treatmentFactorsList)
  df_comparisonPart <- df_comparisonPart %>%
    select(-all_of(colnamesToDelete))

  nameColumnContrast <- paste0("contrastPart", j)
  nameColumnComparison <- paste0("comparisonPart", j)
  data.table::setnames(df_comparisonPart,
                       c("contrastPart", "comparisonPart"),
                       c(nameColumnContrast, nameColumnComparison))
  return(df_comparisonPart)
}
#' compute a data table with all pairwise comparisons of one factor
#'
#' @param treatmentFactorsList list
#' @param i i
#' @noRd
#' @return a data frame with all simple contrasts
#' @keywords internal
#' @importFrom tidyselect all_of
#' @author Christine Paysant-Le Roux
.simpleContrastForOneFactor <- function (treatmentFactorsList, i){

  fixFactor <- groupComparison <- NULL

  contrastPart1 <- contrastPart2 <- contrastPart3 <- contrastPart4 <-  NULL
  comparisonPart1 <- comparisonPart2 <- NULL
  fixPart1 <- fixPart3 <- fixFactor1 <- fixFactor3 <- NULL
  comparisonPart3 <- comparisonPart4 <- NULL

  df_FirstComparisonPart <- .define_partOfSimpleContrast_df(treatmentFactorsList,i,2)
  df_FirstComparisonPart <- df_FirstComparisonPart %>% select(-fixFactor)

  df_SecondComparisonPart <- .define_partOfSimpleContrast_df(treatmentFactorsList,i,1)
  df_simpleContrasts_factor <- cbind(df_FirstComparisonPart, df_SecondComparisonPart)

  df_simpleContrasts_factor <- df_simpleContrasts_factor %>%
    mutate(contrast = paste0("(", contrastPart2,   " - ", contrastPart1, ")"),
           groupComparison = paste0("(", comparisonPart2, " - ", comparisonPart1, ")"))


  # case where fixFactor column is empty (NA inside)
  emptycolFixFactor <- unique(is.na(df_simpleContrasts_factor$fixFactor))
  if(emptycolFixFactor){
    df_simpleContrasts_factor <- df_simpleContrasts_factor %>%
      mutate(contrastName = groupComparison)

  } else {
    df_simpleContrasts_factor <- df_simpleContrasts_factor %>%
      mutate(contrastName = paste0(groupComparison, " in ", fixFactor ))
  }

  df_simpleContrasts_factor <- df_simpleContrasts_factor %>%
    mutate(type = "simple")

  colnamesToDelete <- c("contrastPart2", "comparisonPart2",
                        "contrastPart1", "comparisonPart1")

  df_simpleContrasts_factor <- df_simpleContrasts_factor %>%
    select(-all_of(colnamesToDelete)) %>%
    select("contrast", "groupComparison", "contrastName",
           "type", "fixFactor")



  return(df_simpleContrasts_factor)
}

#' define all simple contrasts
#'
#' @param treatmentFactorsList list of treatment factors
#'
#' @return a data frame with all simple contrasts
#' @importFrom data.table data.table
#' @keywords internal
#' @noRd
#' @author Christine Paysant-Le Roux
.defineAllSimpleContrasts <- function(treatmentFactorsList){
  # create a data table with 5 columns
  # empty data table
  allSimpleContrast_df <- data.table(contrast = character(),
                                     groupComparison = factor(),
                                     contrastName = character(),
                                     type = character(),
                                     fixFactor = factor())
  # create each data frame and rbind it to the allSimpleContrast_df
  for (i in seq_along(treatmentFactorsList)){
    dataTableToCreate <- .simpleContrastForOneFactor(treatmentFactorsList, i)
    allSimpleContrast_df <- rbind(allSimpleContrast_df, dataTableToCreate)
  }
  return(allSimpleContrast_df[])
}

#' Define averaged contrasts
#'
#' @param allSimpleContrast_df a data frame with all the simple contrasts
#' comparisons (test of the difference between any two pairs of groups)
#'
#' @return allAveragedContrasts_df : a data frame with all the averaged
#'  contrasts
#' @keywords internal
#' @noRd
#' @importFrom data.table data.table setcolorder
#' @author Christine Paysant-Le Roux
.define_averaged_contrasts <- function(allSimpleContrast_df){

  groupComparison <- contrast <- n <- fixFactor <- contrastName <- NULL

  allAveragedContrasts_df <- allSimpleContrast_df %>%
    group_by(groupComparison) %>%
    add_tally() %>%
    mutate(contrast= paste0(paste0("(", paste(contrast, collapse=" + ")),")/", n),
           meanIn  = paste(fixFactor, collapse=" + "),
           type    = "mean") %>%
    select(-contrastName, -fixFactor, -n) %>%
    unique() %>% data.table()

  allAveragedContrasts_df <- allAveragedContrasts_df %>%
    mutate(contrastName = paste(groupComparison, "mean", sep = " in "))

  setcolorder(allAveragedContrasts_df,
              c("contrast", "groupComparison", "contrastName",
                "type", "meanIn"))

  return(allAveragedContrasts_df[])
}
# define contrasts for interactions

#' define a data frame with part of interaction contrast
#'
#' @param treatmentFactorsList list
#' @param i i
#' @param j j
#' @param k k
#' @param row_i row_i
#' @param row_j row_j
#'
#' @return a dataframe with part of the interaction contrasts definition
#' @keywords internal
#' @importFrom data.table setDT setnames
#' @importFrom tidyr unite
#' @importFrom tidyselect all_of
#' @noRd
#' @author Christine Paysant-Le Roux
.define_partOfInteractionContrast_df <- function(
    treatmentFactorsList, i, j, k, row_i, row_j) {

  contrastPart_bis <- outsideGroup_bis <- fixFactor_bis <- NULL

  comparisonPart <- treatmentFactorsList
  # combn(x,2) generate all combinations of the elements of x
  # taken 2 at a time
  comparisonPart [[i]] <- combn(treatmentFactorsList[[i]],2)[row_i,]
  comparisonPart [[j]] <- combn(treatmentFactorsList[[j]],2)[row_j,]
  df_comparisonPart <- expand.grid(comparisonPart)
  setDT(df_comparisonPart)
  # paste all the column of the data table
  df_comparisonPart <- df_comparisonPart %>%
    unite(contrastPart_bis, names(df_comparisonPart),
          sep="_", remove=FALSE) %>%
    mutate(contrastPart = contrastPart_bis) %>%
    select(-contrastPart_bis)

  colnameFactor_i <- names(df_comparisonPart)[i]
  colnameFactor_j <- names(df_comparisonPart)[j]

  df_comparisonPart <- df_comparisonPart %>%
    mutate(comparisonPart = df_comparisonPart[[colnameFactor_i]]) %>%
    mutate(fixPart        = df_comparisonPart[[colnameFactor_j]])

  colnamesToKeep <- setdiff(names(treatmentFactorsList),c(colnameFactor_i, colnameFactor_j))
  if(length(colnamesToKeep)){

    df_comparisonPart <- df_comparisonPart %>%
      unite(outsideGroup_bis, all_of(colnamesToKeep), sep="_", remove=FALSE) %>%
      mutate(outsideGroup = outsideGroup_bis) %>%
      select(-outsideGroup_bis)
  }else{

    df_comparisonPart <- df_comparisonPart %>% mutate(outsideGroup = NA)
  }

  colnamesToKeep <- setdiff(names(df_comparisonPart),
                            c("contrastPart", "comparisonPart", "fixPart", "outsideGroup", colnameFactor_i))
  df_comparisonPart <- df_comparisonPart %>%
    unite(fixFactor_bis, all_of(colnamesToKeep), sep="_", remove=FALSE) %>%
    mutate(fixFactor = fixFactor_bis) %>%
    select(-fixFactor_bis)

  colnamesToDelete <- names(treatmentFactorsList)
  df_comparisonPart <- df_comparisonPart %>%
    select(-all_of(colnamesToDelete))


  nameColumnContrast <- paste0("contrastPart", k)
  nameColumnComparison <- paste0("comparisonPart", k)
  nameFixFactor <- paste0("fixFactor", k)
  namePartFixFactor <- paste0("fixPart", k)
  nameOutsideGroup <- paste0("outsideGroup", k)
  setnames(df_comparisonPart,
           c("contrastPart", "comparisonPart", "fixFactor",
             "fixPart", "outsideGroup"),
           c(nameColumnContrast, nameColumnComparison, nameFixFactor,
             namePartFixFactor, nameOutsideGroup))
  return(df_comparisonPart)
}

#' define interaction constrast for pairs of biological factors
#'
#' @param treatmentFactorsList list
#' @param i i
#' @param j j
#'
#' @return a dataframe
#' @importFrom tidyselect all_of
#' @importFrom data.table setnames
#' @keywords internal
#' @noRd
#' @author Christine Paysant-Le Roux
.defineInteractionConstrastForPairsOfFactors <- function(treatmentFactorsList,
                                                         i, j){

  contrastPart1 <- contrastPart2 <- contrastPart3 <- contrastPart4 <- NULL
  comparisonPart1 <- comparisonPart2 <- comparisonPart3 <- comparisonPart4 <- NULL
  fixPart1 <- fixPart3 <- fixFactor1 <- fixFactor3 <- NULL

  df_part1 <- .define_partOfInteractionContrast_df(treatmentFactorsList, i, j, 1, 2, 2)
  df_part2 <- .define_partOfInteractionContrast_df(treatmentFactorsList, i, j, 2, 1, 2)
  df_part3 <- .define_partOfInteractionContrast_df(treatmentFactorsList, i, j, 3, 2, 1)
  df_part4 <- .define_partOfInteractionContrast_df(treatmentFactorsList, i, j, 4, 1, 1)
  df_interactionContrasts <- cbind(df_part1, df_part2, df_part3, df_part4)
  df_interactionContrasts <- df_interactionContrasts %>%
    mutate(contrast = paste0("(", "(", contrastPart1, " - ", contrastPart2, ")"," - ",
                             "(", contrastPart3, " - ", contrastPart4, ")", ")"),
           groupComparison = paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " vs ",
                                    "(", fixPart1, " - ", fixPart3, ")"),
           contrastName  = paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " in ", fixFactor1, " - ",
                                  "(", comparisonPart3, " - ", comparisonPart4, ")", " in ", fixFactor3 ),
           type = "interaction")

  colnamesToDelete <- c("contrastPart1",  "comparisonPart1", "fixFactor1", "fixPart1", "outsideGroup1",
                        "contrastPart2", "comparisonPart2", "fixFactor2", "fixPart2", "outsideGroup2",
                        "contrastPart3", "comparisonPart3", "fixFactor3", "fixPart3", "outsideGroup3",
                        "contrastPart4", "comparisonPart4", "fixFactor4", "fixPart4")
  df_interactionContrasts <- df_interactionContrasts %>%
    select(-all_of(colnamesToDelete))

  setnames(df_interactionContrasts, "outsideGroup4", "outsideGroup")

  df_interactionContrasts <- df_interactionContrasts %>%
    mutate(groupInteraction = paste0(names(treatmentFactorsList)[i],
                                     " vs ", names(treatmentFactorsList)[j]))

  # if 3 factors bio / outsideGroup exist
  if(!is.null(df_interactionContrasts$outsideGroup)){

    ## add factor name of outsideGroup modality
    if(length(names(treatmentFactorsList)[-c(i,j)]) != 0){

      df_interactionContrasts <- df_interactionContrasts %>%
        mutate(outsideGroup = names(treatmentFactorsList)[-c(i,j)]) %>%
        group_by(outsideGroup, groupComparison) %>% add_tally() %>%
        mutate(contrast= paste0(paste0("(", paste(contrast, collapse=" + ")),")/", n),
               contrastName=paste0(groupComparison, " in ", outsideGroup)) %>%
        select(-n) %>% unique()
    }
  }

  return(df_interactionContrasts)
}

#' define all interaction contrasts
#'
#' @param treatmentFactorsList list
#' @param groupInteractionToKeep tokeep
#'
#' @return a dataframe with all the interaction contrasts
#' @importFrom data.table data.table
#' @keywords internal
#' @noRd
#' @author Christine Paysant-Le Roux
.defineAllInteractionContrasts <- function(treatmentFactorsList,
                                           groupInteractionToKeep = NULL){

  groupInteraction <- NULL

  allInteractionsContrasts_df <- data.table(contrast = character(),
                                            groupComparison = factor(),
                                            groupInteraction = character(),
                                            outsideGroup = character(),
                                            contrastName = character(),
                                            type = character())

  vecFori <- combn(length(treatmentFactorsList),2)[1,]
  vecForj <- combn(length(treatmentFactorsList),2)[2,]
  for (k in seq_along(vecFori)) {
    i <- vecFori[k]
    j <- vecForj[k]
    dataTableToCreate <- .defineInteractionConstrastForPairsOfFactors(treatmentFactorsList, i, j)
    allInteractionsContrasts_df <- rbind(allInteractionsContrasts_df,
                                         dataTableToCreate)
  }
  if(!missing(groupInteractionToKeep)){
    allInteractionsContrasts_df <- subset(allInteractionsContrasts_df,
                                          (groupInteraction %in% groupInteractionToKeep))
  }
  return(allInteractionsContrasts_df)
}


# ---- .getContrastMatrix : function generating contrast matrix devlopped by CPL ----

#' get contrast matrix
#'
#' @param ExpDesign data.frame with only bio factors
#' @param factorBio vector of bio factors
#' @param contrastList list of contrast expression
#' @param modelFormula formula
#' @return a list of dataframe with all contrast vectors
#' @keywords internal
#' @noRd
.getContrastMatrixF <- function(ExpDesign, factorBio, modelFormula, contrastList){

    modelFormula <- formula(paste(modelFormula, collapse = " "))

  # bio factor list in formula
  labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
  FactorBioInDesign <- intersect(factorBio, labelsIntoDesign)

  treatmentFactorsList <- lapply(FactorBioInDesign, function(x){paste(x, unique(ExpDesign[[x]]), sep="")})
  names(treatmentFactorsList) <- FactorBioInDesign

  treatmentCondenv <- new.env()

  interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
  isThreeOrderInteraction <- any(attr(terms.formula(modelFormula),"order") == 3)

  # get model matrix
  modelMatrix <- stats::model.matrix(modelFormula, data = ExpDesign)
  colnames(modelMatrix)[colnames(modelMatrix) == "(Intercept)"] <- "Intercept"
  # assign treatment conditions(group) to boolean vectors according to the design model matrix
  .assignVectorToGroups(treatmentFactorsList    = treatmentFactorsList,
                        modelMatrix             = modelMatrix,
                        interactionPresent      = interactionPresent,
                        isThreeOrderInteraction = isThreeOrderInteraction,
                        treatmentCondenv        = treatmentCondenv)
  # get the coefficient vector associated with each selected contrast
  colnamesGLMdesign <- colnames(modelMatrix)

  coefficientsMatrix <- sapply(contrastList, function(x)
    .returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))

  colnames(coefficientsMatrix) <- contrastList

  rownames(coefficientsMatrix) <- colnamesGLMdesign
  contrastMatrix <- as.data.frame(t(coefficientsMatrix))

  return(contrastMatrix)
}

#' compute group binary vector
#'
#' Compute a binary vector (a single vector of 0 and 1) returning the matched
#' group string(s) from a grepl match on the design model matrix colnames.
#'
#' @param biologicalGroups: factor giving group membership
#' (treatment condition associated to one sample)
#' @param colnamesMatrixDesign : vector giving the column names
#' of the model design matrix
#' @param interactionPresent: logical. If TRUE interaction is
#' include in the design model matrix
#'
#' @return a binary vector (a single vector of 0 and 1) returning
#' the matched group string(s) from a grepl match on the design model
#' matrix colnames
#' @keywords internal
#' @noRd
#' @author Christine Paysant-Le Roux
.computeGroupVector <- function(treatmentGroups, colnamesMatrixDesign,
                                interactionPresent,
                                isThreeOrderInteraction = isThreeOrderInteraction) {
  vectorLength <- length(colnamesMatrixDesign)
  groupVector <- rep(0, vectorLength)
  if(interactionPresent){
    simples <- unique(unlist(strsplit(treatmentGroups, "_")))
    order2interaction <- combn(simples,2, FUN=paste, collapse=':')
    toMatchList <- c(simples, order2interaction)
    if(isThreeOrderInteraction){
      order3interaction <- combn(simples,3, FUN=paste, collapse=':')
      toMatchList <- c(toMatchList, order3interaction)
    }
    toMatchList <- paste("^", toMatchList, "$", sep = "")
    #grepl return TRUE for each matched pattern
    pos <- grepl(paste(toMatchList, collapse = "|"), x= colnamesMatrixDesign)
    groupVector <- as.numeric(pos)
    # for intercept
    groupVector[1] <- 1
  } else {
    simples <- unlist(strsplit(treatmentGroups, "_"))
    toMatchList <- simples
    toMatchList <- paste("^", toMatchList, "$", sep = "")
    pos <- grepl(paste(toMatchList, collapse = "|"), x= colnamesMatrixDesign)
    groupVector <- as.numeric(pos)
    # for intercept
    groupVector[1] <- 1
  }
  return(groupVector)
}
#' Assign binary vector to groups
#'
#' @param treatmentFactorsList: list of treatment factors levels
#' @param interaction_in_model: logical. If TRUE interaction is include
#' in the design model matrix
#' @param modelMatrix: numeric matrix giving the design matrix of the GLM.
#' @param treatmentCondenv: the environment to use
#'
#' @return binary vector
#' @keywords internal
#' @noRd
#' @author Christine Paysant-Le Roux
.assignVectorToGroups <- function(treatmentFactorsList = treatmentFactorsList,
                                  modelMatrix = modelMatrix,
                                  interactionPresent = interactionPresent,
                                  isThreeOrderInteraction = isThreeOrderInteraction,
                                  treatmentCondenv = treatmentCondenv){
  # treatment conditions (group) compatible with colnames of the design model matrix
  treatmentGroups <- do.call(paste, c(expand.grid(treatmentFactorsList), sep = "_"))
  # assign binary vector to each group
  modelMatrixColnames <- colnames(modelMatrix)
  # treatmentCondenv <- new.env()
  binaryVectorsList <- lapply(treatmentGroups, function(x)
    .computeGroupVector(x, modelMatrixColnames, interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction))
  names(binaryVectorsList) <- treatmentGroups
  groupDF <- as.data.frame(binaryVectorsList, row.names = modelMatrixColnames)
  mapply(function(x, value) assign(x, value, pos = treatmentCondenv),
         treatmentGroups, binaryVectorsList)
  # http://r.789695.n4.nabble.com/Using-assign-with-mapply-td4681789.html
  # for (i in 1:n) assign(levels[i],
  #  indicatorForModelWithInteraction(levels[i],
  #  colnamesGLMdesign = colnamesGLMdesign), pos = levelsenv)
}

#' return contrast coefficients
#'
#' @param contrast contrast considered
#' @param colnamesGLMdesign colnames
#' @param treatmentCondenv: the environment to use
#' @noRd
#' @return the contrast vector
#' @keywords internal
#' @author Christine Paysant-Le Roux
.returnContrastCoefficients <- function(contrast,
                                        colnamesGLMdesign,
                                        treatmentCondenv){
  expression <- NULL
  if (!is.null(contrast)) {
    expression <- as.character(contrast)
    contrastVector <-rep(0,length(colnamesGLMdesign))
    ej <- parse(text = expression[1])
    # eval evaluates the expression argument in the environment
    #  specified by envir and returns the computed value
    contrastVector <- eval(ej, envir = treatmentCondenv)
  }
  return(contrastVector)
}

# ---- updateSelectedContrasts ----
#' updateSelectedContrasts
#'
#' @param object rflomicsMAE or rflomicsSE object
#' @param contrastList data.frame
#' @noRd
#' @return data.frame contrast
#' @keywords internal
updateSelectedContrasts <- function(object, contrastList=NULL){

  if (!is(object, "RflomicsSE") && !is(object, "RflomicsMAE"))
    stop("")

  allcontrasts <-
    generateExpressionContrast(object) |>
    reduce(rbind) |>
    filter(contrast %in% contrastList$contrast)

  return(allcontrasts)
}
