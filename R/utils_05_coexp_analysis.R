### ============================================================================
### [05_coExp_analysis] function and internal function
### ----------------------------------------------------------------------------
# D. Charif
# N. Bessoltane


################################### CO-EXPRESSION ##############################

#' @title coseq.results.process
#' @param coseqObjectList list of coseq object
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @importFrom stringr str_replace
#' @importFrom stats na.omit
#' @keywords internal
#' @noRd
.coseq.results.process <- function(coseq.res.list, K){

  CoExpAnal <- list()

  replicates <- length(coseq.res.list)

  detail.df <- data.frame()
  results.list <- list()
  for(x in names(coseq.res.list)){

    coseq.res <- coseq.res.list[[x]]

    # total failed
    if(!is.null(coseq.res$error)){
      detail.df <-
        rbind(detail.df,
              c(unlist(strsplit(x, "_")), "failed", coseq.res$error, "NA",
                "NA", "NA"))
    }

    if(!is.null(coseq.res$result)){

      logLike <- likelihood(coseq.res$result)
      ICL.res <- ICL(coseq.res$result)

      for(i in names(logLike)){
        detail.df <-
          rbind(detail.df,
                c(i, unlist(strsplit(x, "_"))[2],
                  ifelse(is.na(logLike[i]) | logLike[i] == 0,
                         "failed", "success"),
                  ifelse(is.na(logLike[i]) | logLike[i] == 0,
                         "The likelihood is equal to 0 or NA", NA),
                  ifelse(is.null(coseq.res$warnings), NA,
                         paste(unique(coseq.res$warnings), collapse = "; ")),
                  logLike[i], ICL.res[i]))
      }

      if(any(detail.df[,3] == "success")){
        results.list[[x]] <- coseq.res$result
      }
    }
  }
  names(detail.df) <-
    c("K", "n", "status", "errors", "warnings", "logLike", "ICL")

  # stats
  jobs.tab.sum <- detail.df |>
    group_by(K, status, errors, warnings) |>
    summarise(prop = round(n()/replicates*100),.groups = 'drop')

  CoExpAnal[["stats"]] <- jobs.tab.sum

  # warning
  CoExpAnal[["warnings"]] <-
    paste(unique(jobs.tab.sum$warnings[!is.na(jobs.tab.sum$warnings)]),
          collapse = "; ")

  if(any(detail.df$status == "success")){

    ICL.list <- list()

    #
    ICL.tab <-
      filter(detail.df, status == "success") %>%
      select(K, n, logLike, ICL) %>%
      mutate(ICL=as.numeric(ICL),
             logLike=as.numeric(logLike),
             n=as.numeric(n))
    ICL.list[["ICL.tab"]] <- ICL.tab

    # Summarize the table: by K, compute the median of the replicate's ICL.
    ICL.n <-
      group_by(ICL.tab, K) %>%
      summarize(median = median(as.numeric(ICL)),
                min = min(as.numeric(ICL)), n=n())
    ICL.list[["ICL.n"]] <- ICL.n

    nb_cluster <- ICL.n[ICL.n$median == min(ICL.n$median),]$K
    # Search for a replicate with a ICL min corresponding to the
    # K with the min median
    min_ICL <- ICL.n[ICL.n$K == nb_cluster,]$min
    min_ICL_rep <- filter(ICL.tab, K == nb_cluster, ICL == min_ICL)$n
    index <- paste0("K=", min(K), "-", max(K), "_", min_ICL_rep[1])
    coseq.res <- results.list[[index]]

    nb_cluster <- as.numeric(str_remove(string = nb_cluster, pattern = "K="))

    nK_success <- sum(ICL.n$n)

    # list of genes per cluster
    clusters_tmp <- clusters(coseq.res)
    if(length(unique(clusters_tmp)) != nb_cluster){
      CoExpAnal[["results"]] <- FALSE
      CoExpAnal[["error"]] <-
        "The optimal number of clusters does not correspond to the minimum ICL."
      return(CoExpAnal)
    }

    clusters <- lapply(seq_len(nb_cluster), function(i){
      names(clusters_tmp[clusters_tmp == i])
    })
    names(clusters) <- paste("cluster", seq_len(nb_cluster), sep = ".")

    #output
    CoExpAnal[["results"]]      <- TRUE
    CoExpAnal[["coseqResults"]] <- coseq.res
    CoExpAnal[["clusters"]]     <- clusters
    CoExpAnal[["cluster.nb"]]   <- nb_cluster
    CoExpAnal[["plots"]]        <- ICL.list

  } else{
    # Pb of convergence: if there is no K.min.rep which correspond to the
    #  median.min, return an error
    CoExpAnal[["results"]] <- FALSE

    CoExpAnal[["errors"]] <-
      paste(unique(jobs.tab.sum$errors[!is.na(jobs.tab.sum$errors)]),
            collapse = "; ")
  }

  return(CoExpAnal)
}

#' @title run Coseq for co-expression analysis on cluster
#' @param counts matrix
#' @param param.list list of coseq parameters
#' @return coseqResults
#' @keywords internal
#' @noRd
#'
.runCoseqLocal <- function(counts, param.list){

  replicates <- param.list[["replicates"]]
  K <- param.list[["K"]]
  iter <- rep(K, replicates)

  args <- list(
    object          = counts,
    model           = param.list[["model"]],
    transformation  = param.list[["transformation"]],
    meanFilterCutoff= param.list[["meanFilterCutoff"]],
    normFactors     = param.list[["normFactors"]],
    GaussianModel   = param.list[["GaussianModel"]],
    verbose         = FALSE
  )

  coseq.res.list <- lapply(seq_len(replicates), function(x){

    res <- .tryCatch_rflomics(
      do.call("coseq", c(args, list(K = K, seed = x, parallel = TRUE))))

    return(res)
  })
  names(coseq.res.list) <- paste0("K=",min(K),"-", max(K), "_", seq_len(replicates))

  coExpAnal <-
    .coseq.results.process(coseq.res.list = coseq.res.list, K = K)

  return(coExpAnal)
}
