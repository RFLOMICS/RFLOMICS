### ============================================================================
### [06_annot_analysis] function and internal function
### ----------------------------------------------------------------------------
# A. Hulot

######## INTERNAL - ANNOTATION CLUSTERPROFILER #########

#' @title see_pathview
#' @param ... Possible arguments for pathview function.
#' @return nothing. Plot on the currently opened device.
#' @keywords internal
#' @importFrom grid grid.raster
#' @importFrom stringr str_split
#' @importFrom pathview pathview
#' @importFrom png readPNG
#' @importFrom utils data
#' @noRd
#'
# Code from: https://stackoverflow.com/questions/60141841/
# how-to-get-pathview-plot-displayed-directly-rather-than-saving-as-a-file-in-r
# It deletes every file created by pathview
.see_pathview <- function(...) {

    if (!exists("bods")) {
        data("bods", package = "pathview")
    }

    msg <- capture.output(pathview(...), type = "message")
    msg <- grep("image file", msg, value = TRUE)
    filename <- sapply(strsplit(msg, " "), function(x)
        x[length(x)])
    if (length(filename) > 0) {
        img <- readPNG(filename)
        grid.raster(img)
        nam <- str_split(filename, "[.]")
        invisible(file.remove(filename))
        invisible(file.remove(paste0(nam[[1]][1], ".xml")))
        invisible(file.remove(paste0(nam[[1]][1], ".png")))
    }

    # rm(bods, envir = .GlobalEnv)

    return()
}

# ---- Valid URL ----
# From https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r

#' @title Test an url
#' @description
#' Called for pathview before trying to get the map.
#'
#' @param url the url to test
#' @param timeout timeout for open.connection
#' @return boolean. TRUE if accessible, FALSE if not
#' @noRd
#' @keywords internal

.validUrl <- function(url, timeout=2){
  con <- url(url)
  check <- suppressWarnings(try(
    open.connection(con,
                    open = "rt",
                    timeout = timeout),
    silent = TRUE)[1])
  suppressWarnings(try(close.connection(con), silent = TRUE))

  return(ifelse(is.null(check), TRUE, FALSE))
}

# ---- .CPR_message_processing ----

#' @title CPR message processing
#' @description
#' CPR message processing
#'
#' @param messagelist list of message clusterProfileR output
#' @return messagelist
#' @importFrom stringr str_detect
#' @noRd
#' @keywords internal
.CPR_message_processing <- function(messagelist){

  Expected_input_gene_ID <- NULL
  for(list in names(messagelist)){
    for(dom in names(messagelist[[list]])){

      if(str_detect(messagelist[[list]][[dom]][1], "Reading KEGG annotation online"))
        messagelist[[list]][[dom]] <- messagelist[[list]][[dom]][-1]

      if(str_detect(messagelist[[list]][[dom]][1], "No gene can be mapped")){

        if(is.null(Expected_input_gene_ID) &&
           str_detect(messagelist[[list]][[dom]][2], "Expected input gene ID")){
          Expected_input_gene_ID <- messagelist[[list]][[dom]][2]
        }
        messagelist[[list]][[dom]] <-
          c(messagelist[[list]][[dom]][1], Expected_input_gene_ID)
      }

      messagelist[[list]][[dom]] <- HTML(paste(messagelist[[list]][[dom]], collapse = "\n"))
    }
  }

  if(length(messagelist) == 0) return(NULL)
  return(messagelist)
}

