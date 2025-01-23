### ============================================================================
### [06_annot_analysis] function and internal function
### ----------------------------------------------------------------------------
# A. Hulot

######## INTERNAL - ANNOTATION CLUSTERPROFILER #########

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

