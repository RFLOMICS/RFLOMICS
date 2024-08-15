### ============================================================================
### [data] function and internal function
### ----------------------------------------------------------------------------
# D. Charif
# N. Bessoltane

#' @title Ecoseed project data
#' @description 
#' This dataset is provided by the EcoSeed project (FP7-KBBE; Impacts of 
#' Environmental Conditions on Seed Quality). that investigates the 
#' effect of seed production temperature on the germination potential of 
#' Arabidopsis thaliana.
#' 
#' This dataset is a multi-omics dataset composed of three data matrices:
#' transcriptomic (raw RNAseq read count data matrix), metabolomic and proteomic 
#' (relative abundance matrix as XIC).
#' 
#' These data are provided in 2 object: ecoseed.df and ecoseed.mae
#' 
#' @name ecoseed
#' @rdname ecoseed
#' @aliases ecoseed
#' @docType data
#' @usage data("ecoseed")
#' @details
#' @keywords datasets
#' @references FP7-KBBE; Impacts of Environmental Conditions on Seed Quality
#' @examples
#' data("ecoseed")
#' 
#' # list of data.frames
#' names(ecoseed.df)
#' head(ecoseed.df$design)
#' 
#' ecoseed.df$RNAtest[1:10, 1:10]
#' 
#' # MultAssayExperiment
#' ecoseed.mae
#' 
#' # An overview of the datasets contained in this package can be found in the vignette "Input Data".
"ecoseed"

#' @name ecoseed.df
#' @description 
#'  \itemize{
#'    \item ecoseed.df: a list of data.frame containing, 
#'      \itemize{
#'        \item design: a data.frame with experiment design,
#'        \item RNAtest: a data.frame with RNAseq data,
#'        \item protetest: a data.frame with proteomics data,
#'        \item metatest: a data.frame with metabolomics data
#'      }
#'  }
#'      
#' @rdname ecoseed
#' @aliases ecoseed.df
"ecoseed.df"

#' @name ecoseed.mae
#' @description 
#' \itemize{
#'    \item ecoseed.mae: a \link{MultiAssayExperiment} object, of RNAtest, 
#'    protetest and metatest data in \link{SummarizedExperiment}
#'    \itemize{
#'      \item ExperimentList class object of length 3: 
#'      \itemize{
#'        \item RNAtest: a \link{SummarizedExperiment} object with RNAseq data,
#'        \item protetest: a \link{SummarizedExperiment} object with proteomics data,
#'        \item metatest: a \link{SummarizedExperiment} object with metabolomics data
#'      }
#'      \item DataFrame with experiment design
#'      \item ...
#'    }
#'  }
#'      
#' @rdname ecoseed
#' @aliases ecoseed.mae
"ecoseed.mae"