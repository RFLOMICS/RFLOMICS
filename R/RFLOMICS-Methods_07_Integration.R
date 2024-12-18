### ============================================================================
### [07_data_integration] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# A. Hulot

######################## METHODS FOR OMICS INTEGRATION ########################



# ---- prepareForIntegration ----
#' @title Preparation step for integration
#' @description This function transforms a RflomicsMAE produced by rflomics
#' into an untrained MOFA object or a list to use for mixOmics.
#' It checks for batch effects to correct them before integration.
#' It also transforms RNASeq counts data into continuous data using
#' \code{\link[limma]{voom}}.
#' This is the second step into the integration.
#' @param object An object of class \link{RflomicsMAE-class}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param omicsNames vector of characters strings,
#' referring to the names of the filtered table in 'object@ExperimentList'.
#' @param rnaSeq_transfo character string, only supports 'limma (voom)'
#' for now.
#' Transformation of the rnaSeq data from counts to continuous data.
#' @param variableLists list of variables to keep per dataset. Default is
#' keeping all features.
#' @param group Not implemented yet in the interface. Useful for MOFA2 run.
#' @param method one of MOFA or mixOmics.
#' Method for which the object is prepared.
#' @param transformData boolean.
#' Transform the data with the transform and normalization method?
#' Default is TRUE.
#' @param cmd used in the interface. Print cmd lines.
#' @return An untrained MOFA object or a list of dataset
#' @exportMethod prepareForIntegration
#' @rdname prepareForIntegration
#' @aliases prepareForIntegration
#' @importFrom MOFA2 create_mofa
#' @example inst/examples/prepareForIntegration.R
setMethod(
  f = "prepareForIntegration",
  signature = "RflomicsMAE",
  definition = function(object,
                        omicsNames = NULL,
                        rnaSeq_transfo = "limma (voom)",
                        variableLists = NULL,
                        group = NULL,
                        method = "MOFA",
                        transformData = TRUE,
                        cmd = FALSE) {
    method <- switch(
      toupper(method),
      "MIXOMICS" = "MixOmics",
      "MOFA"  = "MOFA",
      "MOFA2" = "MOFA",
      "MOFA+" = "MOFA"
    )

    # if no omicsNames we keep all SE
    if (is.null(omicsNames))
      omicsNames <- names(object)

    if (any(!omicsNames %in% names(object)))

    object <- object[, , omicsNames]

    # Checking for batch effects
    correct_batch <- FALSE
    ftypes <- getFactorTypes(object)

    if (any(ftypes == "batch")) {
      correct_batch <- TRUE
    }

    # Transformation before anything else, except for RNAseq data.
    if (transformData) {
      for (SEname in omicsNames) {
        if (getOmicsTypes(object[[SEname]]) != "RNAseq") {
          #object[[SEname]] <- .checkTransNorm(object[[SEname]])
          object[[SEname]] <- getProcessedData(object[[SEname]], norm = TRUE)
        }
        else{
          object[[SEname]] <- getProcessedData(object[[SEname]], filter = TRUE)
        }
      }
    }

    # On each selected omics, according to its type,
    # apply transformation if demanded.
    # Filter DE entities
    for (SEname in omicsNames) {
      SEobject <- object[[SEname]]
      omicsType <- getOmicsTypes(SEobject)

      list_args <- list(
        object = object,
        SEname = SEname,
        correctBatch = correct_batch,
        variableNames = variableLists[[SEname]],
        cmd = cmd
      )

      object <- switch(
        omicsType,
        "RNAseq" = {
          list_args$transformation <- rnaSeq_transfo
          do.call(".rnaseqRBETransform", list_args)
        },
        "proteomics" = do.call(".rbeTransform", list_args),
        "metabolomics" = do.call(".rbeTransform", list_args)
      )
    }

    # Check for duplicated features names across tables
    # MOFA add the entire view name if it's the case, MixOmics do not care.
    # For visualization purpose and coherence,
    # add .index at the end of duplicated variables.
    commonVarNames <- sum(duplicated(unlist(rownames(object))))
    if (commonVarNames > 0) {
      if (cmd) {
        message("[RFLOMICS] #   => Duplicated features names across tables,
                changing names for integration")
      }

      dupTab <- data.frame(
        "dataTable" = rep(names(object),
                          time = vapply(experiments(object), nrow, c(1))),
        "rownames" = unlist(rownames(object)),
        "dup" = duplicated(unlist(rownames(object))) +
          duplicated(unlist(rownames(object)), fromLast = TRUE)
      )
      dupTab <- dupTab[which(dupTab$dup == 1), ]
      omicstochange <- unique(dupTab$dataTable)

      res <- lapply(
        seq_len(length(object)),
        FUN = function(i) {
          SE.object <- object[[names(object)[i]]]
          if (names(object)[i] %in% omicstochange) {
            rownames(SE.object) <- paste(rownames(SE.object), i, sep = ".")
          }
          return(SE.object)
        }
      )
      names(res) <- names(object)

      object <- RflomicsMAE(
        experiments = res,
        colData   = colData(object),
        sampleMap = sampleMap(object),
        omicList    = metadata(object)$omicList,
        projectName = getProjectName(object),
        design      = metadata(object)$design
      )

      if (cmd) {
        message("[RFLOMICS] #   => Done replacing features names")
      }
    }

    # keep only columns corresponding to design factors
    # (remove samples and groups)
    colData(object) <- colData(object)[c(getBioFactors(object),
                                         getBatchFactors(object),
                                         getMetaFactors(object))]

    if (method == "MOFA") {

      MOFAObject <- create_mofa(object,
                                groups = group,
                                extract_metadata = TRUE)

      return(MOFAObject)

    } else if (method == "MixOmics") {
      # Common samples names:
      nsamp <- nrow(colData(object))
      object <- intersectColumns(object)

      if (nsamp != nrow(colData(object))) {
        warning("Removing ",
                nsamp - nrow(colData(object)),
                " samples not present in every experiment.")
      }

      MixOmicsObject <- list(blocks   = lapply(
        experiments(object),
        FUN = function(SE) {
          t(assay(SE))
        }
      ), metadata = colData(object))

      MixOmicsObject$blocks <- lapply(
        MixOmicsObject$blocks,
        FUN = function(mat) {
          mat[match(rownames(mat), rownames(MixOmicsObject$metadata)),]
        }
      )


      return(MixOmicsObject)
    }
  }
)


# ---- Run Omics integration ----

#' @title runOmicsIntegration
#' @description Runs the integration according to the selected method (MOFA or
#' mixOmics) and the settings given by the user. Requires to have the correct
#' entry format in preparedObject before running.
#' @param object An object of class \link{RflomicsMAE-class}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param preparedObject An untrained MOFA object or a list of dataset.
#' Usually a result of prepareForIntegration.
#' @param method one of MOFA or mixOmics.
#' Method for which the object is prepared.
#' @param scale_views boolean. If TRUE, scale each dataset to unit variance.
#' @param maxiter MOFA2 parameter.
#' Number of max iteration (otherwise stop when converged.)
#' @param num_factors MOFA2 parameter. The number of factor to compute.
#' @param selectedResponse character vector, used for mixOmics.
#' Response variables names for block.(s)plsda.
#' @param ncomp mixOmics parameter. Number of components to compute.
#' @param link_datasets mixOmics parameter.
#' Link between datasets in the computation.
#' @param link_response mixOmics parameter. Link between dataset and response.
#' @param sparsity boolean. Used to determine which mixOmics function to apply (either
#' block.plsda if FALSE or block.splsda if TRUE).
#' @param cases_to_try integer. If sparsity is set to TRUE, then cases_to_try
#' is used to determine the number of sets of variables to test for tuning.
#' @param cmd boolean. Used in the interface. If TRUE, print cmd in the console.
#' @param ... not in use at the moment
#' @return a RflomicsMAE object with the correct metadata slot filled with the
#' results and the settings.
#' @rdname runOmicsIntegration
#' @aliases runOmicsIntegration
#' @exportMethod runOmicsIntegration
#' @example inst/examples/runOmicsIntegration.R
#'
setMethod(
  f = "runOmicsIntegration",
  signature = "RflomicsMAE",
  definition = function(object,
                        preparedObject = NULL,
                        method = "MOFA",
                        scale_views = FALSE,
                        maxiter = 1000,
                        num_factors = 10,
                        selectedResponse = NULL,
                        ncomp = 2,
                        link_datasets = 1,
                        link_response = 1,
                        sparsity = FALSE,
                        cases_to_try = 5,
                        cmd = FALSE,
                        ...) {
    method <- switch(
      toupper(method),
      "MIXOMICS" = "MixOmics",
      "MOFA"  = "MOFA",
      "MOFA2" = "MOFA",
      "MOFA+" = "MOFA",
      {stop("This method is unrecognized")}
    )

    switch(toupper(method),
           "MOFA" = {

               if (!is(preparedObject, "MOFA")) {
                   stop("Using MOFA requires that preparedObject is a MOFA object.")
               }

               object <- setMOFA(object, NULL)

               if (cmd)
                   message("[RFLOMICS] #     => Running MOFA analysis")

               MOFA_run <- .runMOFAAnalysis(
                   object = preparedObject,
                   scale_views = scale_views,
                   maxiter = maxiter,
                   num_factors = num_factors
               )

               object <- setMOFA(
                   object,
                   list(
                       "MOFA_results" = MOFA_run$MOFAObject.trained,
                       "MOFA_untrained" = MOFA_run$MOFAObject.untrained,
                       "MOFA_messages" = MOFA_run$MOFA.messages,
                       "settings" = list(
                           scale_views = scale_views,
                           maxiter     = maxiter,
                           num_factors = num_factors,
                           selectData  = names(preparedObject@data)
                       )
                   )
               )

           },
           "MIXOMICS" = {

               if (!is(preparedObject, "list")) {
                   stop("Using mixOmics requires a list as preparedObject.")
               }

               object <- setMixOmics(object, NULL)

               if (cmd)
                   message("[RFLOMICS] #     => Running mixOmics analysis")

               if (is.null(selectedResponse))
                   selectedResponse <- getBioFactors(object)

               MixOmics_res <- lapply(
                   selectedResponse,
                   FUN = function(response_var) {
                       res_mixOmics <- .runMixOmicsAnalysis(
                           object = preparedObject,
                           selectedResponse = response_var,
                           scale_views = scale_views,
                           ncomp = ncomp,
                           link_datasets = link_datasets,
                           link_response = link_response,
                           sparsity = sparsity,
                           cases_to_try = cases_to_try
                       )

                       return(
                           list(
                               "MixOmics_tuning_results" = res_mixOmics$tuning_res,
                               "MixOmics_results"        = res_mixOmics$analysis_res
                           )
                       )
                   }
               )

               names(MixOmics_res) <- selectedResponse
               MixOmics_res$settings <- list(
                   scale_views      = scale_views,
                   ncomp            = ncomp,
                   sparsity         = sparsity,
                   cases_to_try     = cases_to_try,
                   selectedResponse = selectedResponse,
                   selectData  = names(preparedObject$blocks)
               )

               object <- setMixOmics(object, MixOmics_res)

           })

    return(object)
  }
)

# ----  Get a particular multi-omics result ----
#
#' @title Get a particular multi-omics result or settings.
#' @description
#' These methods are used to directly access the results of multi-omics
#' analyses or their settings, usually stored in the metadata of the
#' \link{RflomicsMAE-class} object. Setters are also available.
#'
#' @param object An object of class \link{RflomicsMAE-class}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param response a character giving the response variable to access
#' specifically.
#' @param onlyResults default return only the MixOmics or MOFA2 results.
#' If you want to access all information of the integration,
#' set onlyResults to FALSE.
#' In MixOmics case, works only when response is specified.
#' @return
#' For getters:
#' in getMixOmics, if response is NULL,
#' then all the mixOmics results are returned.
#' Otherwise, it gives the particular mixOmics result.
#' For MOFA, returns the untrained object and the trained object as a list.
#'
#' For setters: always returns a \link{RflomicsMAE-class} object.
#'
#' @exportMethod getMixOmics
#' @aliases getMixOmics
#' @rdname runOmicsIntegration

setMethod(
  f = "getMixOmics",
  signature = "RflomicsMAE",
  definition = function(object,
                        response = NULL,
                        onlyResults = TRUE) {
    toreturn <- metadata(object)[["IntegrationAnalysis"]][["mixOmics"]]

    if (is.null(toreturn)) {
      return(toreturn)
    }

    if (!is.null(response)) {
      toreturn <- toreturn[[response]]
      if (onlyResults)
        toreturn <- toreturn$MixOmics_results
      return(toreturn)
    } else{
      return(toreturn)
    }

  }
)

#' @rdname runOmicsIntegration
#' @exportMethod getMOFA
#' @aliases getMOFA
setMethod(
  f = "getMOFA",
  signature = "RflomicsMAE",
  definition = function(object, onlyResults = TRUE) {
    toreturn <- metadata(object)[["IntegrationAnalysis"]][["MOFA"]]
    if (onlyResults && !is.null(toreturn)) {
      toreturn <- toreturn[["MOFA_results"]]
    }

    return(toreturn)
  }
)
# ---- Get integration setting ----

#' @exportMethod getMOFASettings
#' @rdname runOmicsIntegration
#' @aliases getMOFASettings
setMethod(
  f = "getMOFASettings",
  signature = "RflomicsMAE",
  definition = function(object) {
    return(getMOFA(object, onlyResults = FALSE)$settings)
  }
)

#' @exportMethod getMixOmicsSettings
#' @rdname runOmicsIntegration
#' @aliases getMixOmicsSettings
setMethod(
  f = "getMixOmicsSettings",
  signature = "RflomicsMAE",
  definition = function(object) {
    return(metadata(object)[["IntegrationAnalysis"]][["mixOmics"]][["settings"]])
  }
)

# ---- Set Integration Results ----

#' @rdname runOmicsIntegration
#' @aliases setMOFA
#' @param results The MOFA or mixOmics results to set in the object.
#' If null, set to NULL.
#' @exportMethod setMOFA
setMethod(
  f = "setMOFA",
  signature = "RflomicsMAE",
  definition = function(object, results = NULL) {
    metadata(object)[["IntegrationAnalysis"]][["MOFA"]] <- results
    return(object)
  }
)

#' @rdname runOmicsIntegration
#' @exportMethod setMixOmics
#' @aliases setMixOmics
setMethod(
  f = "setMixOmics",
  signature = "RflomicsMAE",
  definition =  function(object, results = NULL) {
    metadata(object)[["IntegrationAnalysis"]][["mixOmics"]] <- results
    return(object)
  }
)



# ---- MixOmics summary ----

#' @title Get an overview of MixOmics integration results
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character.
#' Useful if MixOmics was run on several response variable.
#' If NULL, all variables are taken into account.
#' @return sumMixOmics: A data frame or a list of dataframe
#' (if selectedResponse is NULL) presenting the summary of mixOmics analyses.
#'
#' @rdname runOmicsIntegration
#' @aliases sumMixOmics
#' @exportMethod sumMixOmics

setMethod(
  f = "sumMixOmics",
  signature = "RflomicsMAE",
  definition =  function(object, selectedResponse = NULL) {
    if (is.null(metadata(object)$IntegrationAnalysis$mixOmics)) {
      stop("It seems this object has no mixOmics results.")
    }

    if (is.null(selectedResponse)) {
      posResponse <- names(metadata(object)$IntegrationAnalysis$mixOmics)
      posResponse <- posResponse[-which(posResponse == "settings")]
      res <- lapply(
        posResponse,
        FUN = function(selResponse) {
          .getOneMORes(object, selectedResponse = selResponse)
        }
      )
      names(res) <- posResponse
      return(res)
    } else {
      .getOneMORes(object, selectedResponse = selectedResponse)
    }
  }
)

#' @title get one MixOmics result
#' @description Get an overview of MixOmics integration results for
#' a specific response variable.
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string.
#' @return A data frame.
#' @keywords internal
#' @noRd

setMethod(
  f = ".getOneMORes",
  signature = "RflomicsMAE",
  definition =   function(object, selectedResponse) {
    res <- getMixOmics(object, response = selectedResponse)
    # Data_res <- res$MixOmics_results
    Data_res <- res

    df <- t(sapply(Data_res$X, dim))
    colnames(df) <- c("Ind", "Features")

    if (!is.null(res$MixOmics_tuning_results)) {
      df <- cbind(df, do.call("rbind", Data_res$keepX))
      colnames(df)[!colnames(df) %in% c("Ind", "Features")] <-
        paste("Comp", seq_len(length(Data_res$keepX[[1]])))
    }

    return(df)
  }
)

# ----- MixOmics: plot variance explained ----

#' @title plotMOVarExp
#'
#' @param object An object of class \link{RflomicsMAE-class}
#' @param selectedResponse a character string of the response variable to
#' consider
#' @param mode Can be NULL (default), "cumulative" or "comp".
#' Defines the type of graph to return
#' @return An object of class \link{RflomicsMAE-class}
#' @importFrom ggpubr ggarrange
#' @keywords internal
#' @noRd
#'
setMethod(
  f = "plotMOVarExp",
  signature = "RflomicsMAE",
  definition =   function(object, selectedResponse, mode = NULL) {
    if (is.null(getMixOmics(object,
                            response = NULL,
                            onlyResults = TRUE))) {
      stop("It seems this object has no mixOmics results.")
    }
    if (is.null(getMixOmics(object,
                            response = selectedResponse,
                            onlyResults = TRUE))) {
      stop("It seems you didn't run MixOmics on this particular variable.")
    }

    Data_res <- getMixOmics(object,
                            response = selectedResponse,
                            onlyResults = TRUE)
    gg_return <- NULL

    if (is.null(mode)) {
      gg_return <- ggarrange(.plot_MO_1(Data_res),
                             .plot_MO_2(Data_res),
                             ncol = 2)
    }
    else if (tolower(mode) == "cumulative") {
      gg_return <- .plot_MO_1(Data_res)
    }
    else if (tolower(mode) == "comp") {
      gg_return <- .plot_MO_2(Data_res)
    }

    return(gg_return)
  }
)
