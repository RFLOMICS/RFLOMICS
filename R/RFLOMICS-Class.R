### ============================================================================
### [RFLOMICS CLASS] accessors and plots
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif

#' @import methods

##==== RflomicsMAE Class ====

#' @name RflomicsMAE-class
#' @rdname RflomicsMAE-class
#' @title RflomicsMAE Class
#' @description
#' RflomicsMAE is a class that extends the \link{MultiAssayExperiment}
#' class by imposing a structure to the metadata slot. This class is used by
#' the Rflomics analysis workflow to store the experimental design, the settings
#' and results of a multi-omics integration analysis.
#' @param object An object of class \link{RflomicsMAE-class}
#' @section Slots:
#'  \itemize{
#'    \item ExperimentList:
#'      \itemize{
#'        \item A ExperimentList class object of \link{RflomicsSE} object
#'        for each assay dataset
#'        }
#'    \item colData: see \code{\link{MultiAssayExperiment}}
#'    \item sampleMap: see \code{\link{MultiAssayExperiment}}
#'    \item metadata:
#'      \itemize{
#'        \item projectName: string. Project name.
#'        \item omicList: list. Contains the list of omics datasets, with the
#'        type and name.
#'        \item design: The experimental design.
#'        \item IntegrationAnalysis: A list containing the multi-omics
#'        integration analysis settings and results.
#'        \item design: The experimental design
#'        \item sessionInfo:
#'        \item IntegrationAnalysis: A list containing the multi-omics
#'        integration analysis settings and results.
#'        }
#' }
#' @section Consructor:
#' \code{\link{createRflomicsMAE}}
#' @section Accessors:
#' @section Plots:
#' @section Methods:
#' \code{\link{generateModelFormulae}}
#' \code{\link{generateExpressionContrast}}
#' \code{\link{runDataProcessing}}
#' \code{\link{runDataProcessing}}
#' \code{\link{runDiffAnalysis}}
#' \code{\link{runCoExpression}}
#' \code{\link{runAnnotationEnrichment}}
#' @seealso \code{\link{MultiAssayExperiment}}
#' @aliases RflomicsMAE-class
#' @return A \code{RflomicsMAE} object.
#' @exportClass RflomicsMAE
#' @example inst/examples/loadData.R
setClass(
  Class    = "RflomicsMAE",
  contains = "MultiAssayExperiment",
  validity = function(object) {
    metadata <- metadata(object)

    # check if metadat is list()
    if (!is.list(metadata)) {
      return("The 'metadata' slot must be a list.")
    }

    # Vérification des éléments requis dans `metadata`
    required_elements <-
      c("omicList", "projectName", "design", "IntegrationAnalysis",
        "date","sessionInfo","rflomicsVersion")

    missing_elements <- setdiff(required_elements, names(metadata))
    if (length(missing_elements) > 0) {
      return(paste("The 'metadata' slot must contain the following elements:",
                   paste(missing_elements, collapse = ", ")))
    }

    # Vérification des types des éléments
    if (!is.character(metadata$projectName)) {
      return("The 'projectName' in 'metadata' must be a single string.")
    }
    if (!inherits(metadata$date, "Date")) {
      return("The 'date' in 'metadata' must be of class 'Date'.")
    }
    if (!inherits(metadata$rflomicsVersion, "package_version")) {
      return("The 'rflomicsVersion' in 'metadata' must be a package_version.")
    }
    if (!is.list(metadata$omicList)) {
      return("The 'omicList' in 'metadata' must be a list.")
    }
    if (!is.list(metadata$IntegrationAnalysis)) {
      return("The 'IntegrationAnalysis' in 'metadata' must be a list.")
    }
    if (!is.list(metadata$design)) {
      return("The 'design' in 'metadata' must be a list.")
    }

    TRUE
  }
)

##==== RflomicsSE Class ====

#' @title RflomicsSE Class
#' @name RflomicsSE-class
#' @rdname RflomicsSE-class
#' @description
#' RflomicsSE is a class that extends the \link{SummarizedExperiment} by imposing a structure
#' on the metadata slot. This class is used by the Rflomics analysis workflow to store the
#' experimental design, the settings and results of a single omic analysis.
#' The slot metadata is structured as follows:
#' @param object An object of class \link{RflomicsSE}
#' @section Slots:
#' See \link{SummarizedExperiment}
#'
#' The slot metadata is structured as follows:
#'  \itemize{
#'    \item omicType: the type of omics dataset
#'    \item design: experimental design
#'    \item DataProcessing: a list containing the data processing settings and results
#'    \item PCAlist: a list containing the PCA settings and results
#'    \item DiffExpAnal: a list containing the Differential Analysis settings and results
#'    \item CoExpAnal: a list containing the Coexpression Analysis settings and results
#'    \item DiffExpEnrichAnal: a list containing the enrichment analysis of the list of DE features settings and results
#'    \item CoExpEnrichAnal: a list containing the enrichment analysis of the list of co-expressed features settings and results
#'  }
#' @section Accessors:
#' @seealso \link{SummarizedExperiment}
#' @aliases RflomicsSE
#' @exportClass RflomicsSE
#' @return A \code{RflomicsSE} object.
setClass(
  Class    = "RflomicsSE",
  contains = "SummarizedExperiment",
  validity = function(object) {
    metadata <- metadata(object)

    # check if metadat is list()
    if (!is.list(metadata)) {
      return("The 'metadata' slot must be a list.")
    }

    # Vérification des éléments requis dans `metadata`
    required_elements <-
      c("omicType","design","DataProcessing","PCAlist","DiffExpAnal",
        "CoExpAnal","DiffExpEnrichAnal", "CoExpEnrichAnal")

    missing_elements <- setdiff(required_elements, names(metadata))
    if (length(missing_elements) > 0) {
      return(paste("The 'metadata' slot must contain the following elements:",
                   paste(missing_elements, collapse = ", ")))
    }

    TRUE
  }
)

