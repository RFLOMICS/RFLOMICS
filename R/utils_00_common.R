### ============================================================================
### [00_common] function and internal function
### ----------------------------------------------------------------------------
# N. Bessoltane,
# A. Hulot,
# D. Charif

# ---- INTERNAL FUNCTIONS ----
# ---- isContrastName : ----
#' @title Check if character vectors are contrasts Names
#' @param object a MAE object or a SE object (produced by Flomics).
#' If it's a RflomicsSE, expect to find
#'  a slot of differential analysis.
#' @param contrastName vector of characters.
#' @return boolean. TRUE if all of contrastName are indeed contrasts Names.
#' @noRd
#' @keywords internal
.isContrastName <- function(object, contrastName) {

  SelectedContrasts <- getSelectedContrasts(object)

  if(nrow(SelectedContrasts) == 0) return(NULL)

  SelectedContrasts <- SelectedContrasts$contrastName

  if(is.null(SelectedContrasts)) return(NULL)

  return(!any(!contrastName %in% SelectedContrasts))
}

# ---- isClusterName - Check if character vectors are tags Names : -----

#' @title Check if character vectors is a cluster name
#'
#' @param object a SE object (produced by Flomics). Expects to find
#'  a slot of coExpression analysis
#' @param clusterName vector of characters. For clusters, please
#' specify cluster.1, cluster.2, ... although 1,2,3 can work as well.
#' @return boolean. TRUE if all of tagName are indeed tags Names.
#' @noRd
#' @keywords internal
.isClusterName <- function(object, clusterName) {

  Cluster.list <- getCoexpClusters(object)
  if(is.null(Cluster.list)) return(NULL)

  return(!any(!clusterName %in% names(Cluster.list)))
}

# ---- .getOrigin - get origin of a particular name ----
#
#' @title get origin of a name given a rflomics MAE
#'
#' @param object a RflomicsSE, produced by rflomics
#' @param name name of the parameter to identify. For clusters, please
#' specify cluster.1, cluster.2, etc.
#' @return The origin of the name, one of Contrast, Tag or CoexCluster.
#' @noRd
#' @keywords internal

.getOrigin <- function(object, name) {

  if (.isContrastName(object, name)) return("DiffExp")
  if (.isClusterName(object, name)) return("CoExp")

  return("NoOriginFound")
}

# ---- tryCatch_rflomics - catch error, warning et message : -----

#' @title catch error, warning et message and results
#'
#' @param f function and its parameters
#' @return A list with 4 elements: message, warning, error and results
#' @noRd
#' @keywords internal
.tryCatch_rflomics <- function(f) {

  message_capture <- character()
  warning_capture <- character()

  result <- tryCatch({
    withCallingHandlers({
      output <- f
      list(result = output, messages = message_capture,
           warnings = warning_capture, error = NULL)
    }, message = function(m) {
      message_capture <<- c(message_capture, m$message)
      invokeRestart("muffleMessage")
    }, warning = function(w) {
      warning_capture <<- c(warning_capture, w$message)
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    list(result = NULL, messages = NULL, warnings = NULL, error = e$message)
  })

  # Inclure les messages et avertissements capturés dans le résultat
  if (length(message_capture) > 0) {
    result$messages <- message_capture
  }
  if (length(warning_capture) > 0) {
    result$warnings <- warning_capture
  }

  return(result)
}


#' @title .tryRflomics
#' @details
#' This function comes from
#' https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#  The author indicated that he merged Martins solution
#  (https://stackoverflow.com/a/4952908/2161065) and
#  the one from the R-help mailing list you get with demo(error.catching).
#' @param expr The expression that has been to be evaluated.
#' @return a named list
#' \itemize{
#' \item{\code{value:} }{The results of the expr evaluation or
#' NULL if an error occured }
#' \item{\code{warning:} }{warning message or NULL}
#' \item{\code{error:} }{error message or NULL}
#' }
#' @keywords internal
#' @noRd

.tryRflomics <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr,
             error    = function(e){ err <- e
             NULL
             }),
    warning = function(w){ warn <- w
    invokeRestart("muffleWarning")}
  )

  return(list(value = value,
              warning = warn,
              error = err)
  )

}

# ---- get package version  -----

#' @title get package version from info session stored in rflomicsMAE object
#'
#' @param rflomicsMAE rflomicsMAE object
#' @param package package name
#' @param info type of information about package (version, )
#' @return value of information
#' @noRd
#' @keywords internal
.getPackageInfo <- function(rflomicsMAE, package = NULL, info = "version") {

  sessionInfo <- metadata(rflomicsMAE)[["sessionInfo"]][["pkg.tab"]]

  if(is.null(package)) return(NULL)
  if(!info %in% c("version")) return(NULL)

  return(sessionInfo[[info]][sessionInfo[["package"]] == package])
}


## ---- writeSessionInfo ----
#' writeSessionInfo
#'
#' @return A list
#' @keywords internal
#' @importFrom reticulate py_config
#' @noRd
.writeSessionInfo <- function(){

  x <- sessionInfo()

  mkLabel <- function(L, n, type) {
    vers <- sapply(L[[n]], function(x) x[["Version"]])
    pkg <- sapply(L[[n]], function(x) x[["Package"]])
    data.frame( package=pkg, version=vers, type= type)
  }

  pkg <- list()

  if (!is.null(x$otherPkgs)) {
    pkg[[1]] <- mkLabel(x, "otherPkgs","attached")
  }
  if (!is.null(x$loadedOnly)) {
    pkg[[2]] <- mkLabel(x, "loadedOnly","loadedOnly")
  }

  return(list(
    "Rsession.conf" = c("R version"=x$R.version$version.string,
                        "Platform"= x$R.version$platform,
                        "OS"= x$running),
    "pkg.tab" = do.call(rbind, pkg),
    "python.conf" = py_config()))
}

## ---- getKEGGRelease ----
#' getKEGGRelease
#'
#' @return string
#' @keywords internal
#' @importFrom httr content GET
#' @importFrom stringr str_extract
#' @noRd
.getKEGGRelease <- function(){

  version <- NULL
  url <- "https://rest.kegg.jp/info/kegg"

  try.res  <- .tryCatch_rflomics(GET(url))
  response <-  try.res$result

  if(!is.null(response)){

    content <- content(response, "text")

    version <- str_extract(content, "Release \\d+\\.\\d+") |>
      str_remove(pattern = "Release ")
  }

  if(!is.null(try.res$warnings))
    message("warning: ", try.res$warnings)

  if(!is.null(try.res$error))
    message("error: ", try.res$error)


  return(version)
}

## ---- setElementToMetadata ----
#' .setElementToMetadata
#'
#' @param object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results..
#' @param name the name of element to add to metadata slot.
#' @param subName the name of sub element to add to metadata slot.
#' @param content the content of element to add
#' @return object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE}.
#' @keywords internal
#' @noRd
.setElementToMetadata <- function(object,
                                  name = NULL,
                                  subName = NULL,
                                  content = NULL){

  if(is.null(name))
    warning("Argument name is required")
  # if(is.null(content))
  #   warning("Argument content is required")

  if(!is.null(subName)){
    metadata(object)[[name]][[subName]] <- content
  } else{
    metadata(object)[[name]] <- content
  }

  # to do
  # Check the consistency of the object.
  # checkRflomicsMAE(object)

  return(object)
}


## ---- getAnalysis ----
#' getAnalysis
#'
#' @param object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results..
#' @param name the name of element to add to metadata slot.
#' @param subName the name of sub element to add to metadata slot.
#' @return list
#' @keywords internal
#' @noRd
.getAnalysis <- function(object,
                         name = NULL,
                         subName = NULL){

  results <- list()

  if(!name %in% names(metadata(object)))
    stop("There are no analysis results with this name:", name)

  results <- metadata(object)[[name]]

  if(!is.null(subName)){
    if(!subName %in% names(results))
      return(NULL)

    results <- results[[subName]]
  }

  return(results)
}

