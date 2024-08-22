### ============================================================================
### [data] function and internal function
### ----------------------------------------------------------------------------
# D. Charif
# N. Bessoltane

# ---- Data frame list ----

#' @title Ecoseed project data
#' @name ecoseed.df
#' @rdname ecoseed.df
#' @docType data
#' @keywords datasets
#' @usage data("ecoseed.df")  
#' @description 
#' This dataset is provided by the EcoSeed project (FP7-KBBE; Impacts of 
#' Environmental Conditions on Seed Quality). that investigates the 
#' effect of seed production temperature on the germination potential of 
#' Arabidopsis thaliana.
#' 
#' This dataset is a multi-omics dataset composed of three data matrices:
#' transcriptomics (raw RNAseq read count data matrix), metabolomics and proteomics
#' (relative abundance matrix as XIC).
#' 
#' These data are provided in 2 object: ecoseed.df and ecoseed.mae
#' @format A list of data.frame containing
#'      \itemize{
#'        \item design: a data.frame with experiment design,
#'        \item RNAtest: a data.frame with RNAseq data,
#'        \item protetest: a data.frame with proteomics data,
#'        \item metatest: a data.frame with metabolomics data
#'      }
#' @references FP7-KBBE; Impacts of Environmental Conditions on Seed Quality
#' @examples
#' data("ecoseed.df")
#' 
#' # list of data.frames
#' names(ecoseed.df)
#' head(ecoseed.df$design)
NULL

# ---- MAE ----

#' @title Ecoseed project data
#' @name ecoseed.mae
#' @rdname ecoseed.mae
#' @usage data("ecoseed.mae")    
#' @docType data
#' @keywords datasets
#' @description 
#' This dataset is provided by the EcoSeed project (FP7-KBBE; Impacts of 
#' Environmental Conditions on Seed Quality). that investigates the 
#' effect of seed production temperature on the germination potential of 
#' Arabidopsis thaliana.
#' 
#' This dataset is a multi-omics dataset composed of three data matrices:
#' transcriptomics (raw RNAseq read count data matrix), metabolomics and proteomics
#' (relative abundance matrix as XIC).
#' 
#' These data are provided in 2 object: ecoseed.df and ecoseed.mae
#' 
#' @format ecoseed.mae: a \link{MultiAssayExperiment} object, of RNAtest, 
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
#' @references FP7-KBBE; Impacts of Environmental Conditions on Seed Quality
#' @examples
#' data("ecoseed.mae")
NULL
