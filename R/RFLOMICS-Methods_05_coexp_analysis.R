### ============================================================================
### [05_co-exp_analysis] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# D. Charif,
# N. Bessoltane,
# A. Hulot

##==== STAT METHOD ====

###==== METHOD runCoExpressions ====

#' @title Run CoExpression Analysis and process results
#' @name runCoExpression
#' @aliases runCoExpression,RflomicsSE-method
#' @description This method performs a co-expression/co-abundance analysis of
#' omic-data.
#' @details For now, only the coseq function of the coseq package is used.
#' For RNAseq data, parameters used are those recommended in DiCoExpress
#' workflow (see the reference). This parameters are: \code{model="normal"},
#' \code{transformation="arcsin"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="TMM"}, \code{meanFilterCutoff = 50}
#' For proteomic or metabolomic, data are scaled by protein or metabolite
#' to group them by expression profiles rather than by expression intensity.
#' After data scaling, recommended parameters (from \code{coseq} developers)
#' for co-expression analysis are:
#' \code{model="normal"}, \code{transformation="none"},
#' \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="none"},  \code{meanFilterCutoff = NULL}.
#' @return
#' An S4 object of class \link{RflomicsSE}
#' All the results are stored as a named list \code{CoExpAnal}
#'  in the metadata slot of a given \code{RflomicsSE} object. Objects are:
#' The runCoExpression method return several results, for \link{coseq}
#' method, objects are:
#' \itemize{
#' \item \code{settings:}  co-expression analysis settings. See \code{getCoexpSetting}
#' \item \code{results:}  boolean indicating if the co-expression analysis succeed
#' \itemize{
#' \item \code{coseqResults:} the raw results of \code{coseq}
#' \item \code{clusters:}  a List of clusters
#' \item \code{cluster.nb:}  The number of cluster
#' \item \code{plots:}  The plots of \code{coseq} results
#' \item \code{stats:}  A tibble summarising failed jobs: reason, propoif any
#' }
#' \item \code{errors:} error list.
#' }
#' @param object An object of class \link{RflomicsSE} or
#' class \link{RflomicsMAE-class}
#' @param SE.name SE.name the name of the dataset if the input object
#' is a \link{RflomicsMAE-class}
#' @param contrastNames names of the contrasts from which the DE entities
#' have to be taken. Can be NULL, in that case every contrasts from the differential
#' analysis are taken into consideration.
#' @param K Number of clusters (a single value or a vector of values)
#' @param replicates The number of iteration for each K.
#' @param model Type of mixture model to use \code{"Poisson"} or
#' \code{"normal"}. By default, it is the normal.
#' @param GaussianModel Type of \code{GaussianModel} to be used for the
#' Normal mixture model only. This parameters
#' is set to \code{"Gaussian_pk_Lk_Ck"} by default and doesn't have to
#' be changed except if an error message proposed
#' to try another model like \code{"Gaussian_pk_Lk_Bk"}.
#' @param transformation The transformation type to be used. By default,
#' it is the "arcsin" one.
#' @param normFactors The type of estimator to be used to normalize for
#' differences in library size.
#' By default, it is the "TMM" one.
#' @param merge \code{"union"} or \code{"intersection"}
#' @param meanFilterCutoff a cutoff to filter a gene with a mean expression lower than this value.
#' (only for RNAseq data, set to NULL for others).
#' @param scale Boolean. If TRUE scale all variables tounit variance.
#' @param min.data.size The minimum allowed number of variables (default: 100)
#' @param ... Additional arguments.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress:
#' a tool to process multifactorial RNAseq experiments from quality controls
#' to co-expression analysis through differential analysis based on contrasts
#' inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runCoExpression
#' @rdname runCoExpression
#' @section Accessors:
#' A set of getters and setters generic functions to access and
#' modify objects of the slot metadata of a \link{RflomicsMAE-class} object or
#' a \link{RflomicsMAE-class} object.
#' @section Plots:
#' A collection of functions for plotting results from omic analysis steps.
#' @seealso \code{\link[coseq]{coseq}}
#' @seealso \code{\link{createRflomicsMAE}}
#' @seealso \code{\link{generateModelFormulae}}
#' @seealso \code{\link{generateExpressionContrast}}
#' @seealso \code{\link{runDataProcessing}}
#' @seealso \code{\link{runDiffAnalysis}}
#' @example inst/examples/runCoExpression.R
setMethod(
  f = "runCoExpression",
  signature = "RflomicsSE",
  definition = function(object,
                        K                = 2:20,
                        replicates       = 5,
                        contrastNames    = NULL,
                        merge            = "union",
                        model            = "Normal",
                        GaussianModel    = NULL,
                        transformation   = NULL,
                        normFactors      = NULL,
                        meanFilterCutoff = NULL,
                        scale            = NULL,
                        min.data.size    = 100,
                        ...){

    # define result output
    CoExpAnal <- list(
      settings = list(),
      results  = list(),
      errors   = NULL
    )

    # default methods
    default.methods <-
      c(
        list(
          merge         = c("union", "intersection"),
          model         = "Normal",
          GaussianModel.all = c("Gaussian_pk_Lk_Ck", "Gaussian_pk_Lk_Bk")
        ),
        switch (
          getOmicsTypes(object),
          "RNAseq" = {
            list(
              transformation   = "arcsin",
              normFactors      = "TMM",
              meanFilterCutoff = 50,
              GaussianModel = "Gaussian_pk_Lk_Ck",
              scale = FALSE
            )
          },
          list(
            transformation   = "none",
            normFactors      = "none",
            meanFilterCutoff = NULL,
            GaussianModel = "Gaussian_pk_Lk_Bk",
            scale = TRUE
          )
        )
      )

    if (is.null(getDEMatrix(object))) {
      stop("Please run a differential analysis.
                       runCoExpression uses these results.")
    }

    validContrasts <- getValidContrasts(object)
    if(is.null(validContrasts) || nrow(validContrasts) == 0){
      validContrasts <- getSelectedContrasts(object)

      if(is.null(validContrasts) || nrow(validContrasts) == 0)
        stop("No defined contrasts")
    }

    if (is.null(contrastNames)){
      contrastNames <- validContrasts$contrastName
    }
    else{
      contrastNames <- intersect(contrastNames, validContrasts$contrastName)
      if(length(contrastNames) == 0)
        stop("No defined contrasts")
    }

    # check param
    if(!merge %in% default.methods[["merge"]])
      stop("Invalid value for argument 'merge'. The allowed values are ",
           paste(default.methods[["merge"]], collapse = " or "))
    if(!model %in% default.methods[["model"]])
      stop("Invalid value for argument 'model'. The allowed values are ",
           paste(default.methods[["model"]], collapse = " or "))

    if(is.null(GaussianModel)) GaussianModel <- default.methods[["GaussianModel"]]
    if(!GaussianModel %in% default.methods[["GaussianModel.all"]])
      stop("Invalid value for argument 'GaussianModel'. The allowed values are ",
           paste(default.methods[["GaussianModel.all"]], collapse = " or "))

    if(is.null(scale)) scale <- default.methods[["scale"]]
    if(!is.null(default.methods[["scale"]]) && scale != default.methods[["scale"]])
      stop("Invalid value for argument 'scale'. The allowed values are ",
           paste(default.methods[["scale"]], collapse = " or "))

    if(is.null(normFactors)) normFactors <- default.methods[["normFactors"]]
    if(!is.null(default.methods[["normFactors"]]) &&
       normFactors != default.methods[["normFactors"]])
      stop("Invalid value for argument 'normFactors'. The allowed values are ",
           paste(default.methods[["normFactors"]], collapse = " or "))

    if(is.null(transformation)) transformation <- default.methods[["transformation"]]
    if(!is.null(default.methods[["transformation"]]) &&
       transformation != default.methods[["transformation"]])
      stop("Invalid value for argument 'transformation'. The allowed values are ",
           paste(default.methods[["transformation"]], collapse = " or "))

    if(is.null(meanFilterCutoff)) meanFilterCutoff <- default.methods[["meanFilterCutoff"]]
    if(!is.null(default.methods[["meanFilterCutoff"]]) &&
       meanFilterCutoff != default.methods[["meanFilterCutoff"]])
      stop("Invalid value for argument 'meanFilterCutoff'. The allowed values are ",
           paste(default.methods[["meanFilterCutoff"]], collapse = " or "))

    message("[RFLOMICS] #     => Tool: coseq... ")

    CoExpAnal[["settings"]] <- list(
      "method" = "coseq",
      "model"            = model,
      "GaussianModel"    = GaussianModel,
      "transformation"   = transformation,
      "normFactors"      = normFactors,
      "meanFilterCutoff" = meanFilterCutoff,
      "contrastNames"  = contrastNames,
      "merge"       = merge,
      "replicates"       = replicates,
      "K"          = K,
      "scale"            = scale
    )
    names(CoExpAnal[["settings"]][["contrastNames"]]) <- contrastNames

    geneList <- getDEList(object = object,
                          contrasts = contrastNames,
                          operation = merge)

    if(length(geneList) < min.data.size)
      stop("The number of ", .omicsDic(object)$variableName,
           "s must be greater than ", min.data.size)

    # set default parameters based on data type
    counts <-
      switch(
        getOmicsTypes(object),
        "RNAseq" = {
          object2 <- getProcessedData(object, filter = TRUE)
          assay(object2)[geneList,]
        },
        {
          object2 <- getProcessedData(object, norm = TRUE)
          assay(object2)[geneList,]
        }
      )

    if(isTRUE(scale)){

      message("[RFLOMICS] # Scale each ", .omicsDic(object)$variableName,
              " (center = TRUE, scale = TRUE)")

      counts[] <-
        t(apply(counts,1,function(x){scale(x, center = TRUE, scale = TRUE) }))
    }

    # run coseq : on local machine or remote cluster
    used_function <- ".runCoseqLocal"
    # switch(
    #   as.character(clustermq),
    #   `FALSE` = ".runCoseqLocal",
    #   `TRUE`  = ".runCoseqClustermq"
    # )

    CoExpAnal[["results"]] <-
      do.call(used_function, list(counts, param.list = CoExpAnal[["settings"]]))

    if(!is.null(CoExpAnal[["results"]]$cluster.nb))
      message("[RFLOMICS] #   => Number of clusters: ",
              CoExpAnal[["results"]]$cluster.nb)

    if(!is.null(CoExpAnal[["results"]]$error))
      CoExpAnal[["errors"]] <- CoExpAnal[["results"]]$error

    object <-
      setElementToMetadata(object,
                           name    = "CoExpAnal",
                           content = CoExpAnal)
    return(object)
  })

#' @rdname runCoExpression
#' @name runCoExpression
#' @aliases runCoExpression,RflomicsMAE-method
#' @exportMethod runCoExpression
setMethod(
  f          = "runCoExpression",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        K                = 2:20,
                        replicates       = 5,
                        contrastNames    = NULL,
                        merge            = "union",
                        model            = "Normal",
                        GaussianModel    = NULL,
                        transformation   = NULL,
                        normFactors      = NULL,
                        meanFilterCutoff = NULL,
                        scale            = NULL,
                        min.data.size    = 100,
                        ...){

    if (!SE.name %in% names(object))
      stop(SE.name, " is not part of ", object)

    object[[SE.name]] <-
      runCoExpression(
        object           = object[[SE.name]],
        K                = K,
        replicates       = replicates,
        contrastNames    = contrastNames,
        merge            = merge,
        model            = model,
        GaussianModel    = GaussianModel,
        transformation   = transformation,
        normFactors      = normFactors,
        meanFilterCutoff = meanFilterCutoff,
        scale            = scale,
        min.data.size    = min.data.size,
        ...)

    return(object)
  })



##==== GRAPHICAL METHODS ====

### ---- getCoExpAnalysesSummary ----

#' @param omicNames the name of the experiment to summarize.
#' @section Plots:
#' \itemize{
#'    \item getCoExpAnalysesSummary: ...
#' }
#' @exportMethod getCoExpAnalysesSummary
#' @rdname runCoExpression
#' @name getCoExpAnalysesSummary
#' @aliases getCoExpAnalysesSummary,RflomicsMAE-method
setMethod(
  f = "getCoExpAnalysesSummary",
  signature = "RflomicsMAE",
  definition = function(object,
                        omicNames = NULL) {
    mean.y_profiles.list <- list()

    if (is.null(omicNames))
      omicNames <- getAnalyzedDatasetNames(object, "CoExpAnal")

    for (data in omicNames) {

      SE <- getProcessedData(object[[data]], filter = TRUE)

      CoExpAnal <- getAnalysis(SE, name = "CoExpAnal")$results

      Groups     <- getDesignMat(SE)
      cluster.nb <-
        CoExpAnal$cluster.nb
      coseq.res  <-
        CoExpAnal[["coseqResults"]]

      mean.y_profiles.list[[data]] <-
        lapply(seq_len(cluster.nb), function(cluster) {
          assays.data <- filter(as.data.frame(coseq.res@assays@data[[1]]),
                                get(paste0("Cluster_", cluster)) > 0.8)

          y_profiles.gg <-
            coseq.res@y_profiles[rownames(assays.data), ] %>%
            data.frame() %>%
            mutate(observations = rownames(.)) %>%
            melt(id = "observations", value.name = "y_profiles") %>%
            rename(samples = variable) %>%
            full_join(Groups , by = "samples")

          y_profiles.gg %>% group_by(groups) %>%
            summarise(mean = mean(y_profiles)) %>%
            mutate(cluster = paste0("cluster.", cluster))

        }) %>% reduce(rbind) %>% mutate(dataset = data)

    } %>% reduce(rbind)

    mean.y_profiles.gg <- reduce(mean.y_profiles.list, rbind)

    if (nrow(mean.y_profiles.gg) == 0)
      return(NULL)

    p <- ggplot(data = mean.y_profiles.gg, aes(x = groups, y = mean, group = 1)) +
      geom_line(aes(color = as.factor(cluster))) +
      geom_point(aes(color = as.factor(cluster))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      facet_grid(cols = vars(cluster), rows = vars(dataset), scales = "free") +
      labs(x = "Conditions", y = "Expression profiles mean")

    return(p)
  }
)

### ---- plotCoExpression ----
#' @name plotCoExpression
#' @aliases plotCoExpression,RflomicsSE-method
#' @section Plots:
#' \itemize{
#'    \item plotCoExpression:
#'    list plot of ICL, logLike and coseq object with min ICL
#' }
#' @exportMethod plotCoExpression
#' @rdname runCoExpression
setMethod(
  f="plotCoExpression",
  signature="RflomicsSE",
  definition <- function(object){

    object <- getProcessedData(object, filter = TRUE)

    CoExpAnal <- getAnalysis(object, name = "CoExpAnal")$results

    if(length(CoExpAnal) == 0) stop("No co-expression results!")

    Groups <- getDesignMat(object)

    coseq.res     <- CoExpAnal[["coseqResults"]]
    ICL.list      <- CoExpAnal[["plots"]]
    replicates    <- getCoexpSettings(object)$replicates
    K             <- getCoexpSettings(object)$K

    #### Plots
    ### plot ICL
    ICL.list[["ICL.tab"]] <-
      merge(ICL.list[["ICL.tab"]],
            expand.grid(K = paste("K", K, sep="="),
                        n = seq_len(replicates)), all = TRUE)

    ICL.list[["ICL.n"]] <-
      merge(ICL.list[["ICL.n"]],
            data.frame(K = paste("K", K, sep="=")), all = TRUE) %>%
      mutate(n = ifelse(is.na(n), 0, n))

    ICL.p <- ggplot(data = ICL.list[["ICL.tab"]]) +
      geom_boxplot(aes(x = K, y = ICL, group = K), na.rm = TRUE) +
      geom_text(data = ICL.list[["ICL.n"]],
                aes(x = seq_len(length(K)),
                    y = max(ICL.list[["ICL.tab"]]$ICL, na.rm = TRUE),
                    label = paste0("n=", n)),
                col = 'red', size = 4) +
      ylim(min(ICL.list[["ICL.tab"]]$ICL, na.rm = TRUE),
           max(ICL.list[["ICL.tab"]]$ICL, na.rm = TRUE))

    ### plot logLike
    logLike.p <- ggplot(data = ICL.list[["ICL.tab"]]) +
      geom_boxplot(aes(x = K, y = logLike, group = K), na.rm = TRUE) +
      geom_text(data = ICL.list[["ICL.n"]],
                aes(x = seq_len(length(K)),
                    y = max(ICL.list[["ICL.tab"]]$logLike, na.rm = TRUE),
                    label = paste0("n=", n)),
                col = 'red', size = 4) +
      ylim(min(ICL.list[["ICL.tab"]]$logLike, na.rm = TRUE),
           max(ICL.list[["ICL.tab"]]$logLike, na.rm = TRUE))

    ### coseq plots
    plot.coseq.res <- coseq::plot(coseq.res,
                           conds = Groups$groups,
                           collapse_reps = "average",
                           graphs = c("profiles",
                                      "boxplots",
                                      "probapost_boxplots",
                                      "probapost_barplots",
                                      "probapost_histogram"))


    return(c(plot.coseq.res, list("ICL" = ICL.p, "logLike" = logLike.p)))
  })

#' @rdname runCoExpression
#' @name plotCoExpression
#' @aliases plotCoExpression,RflomicsMAE-method
#' @exportMethod plotCoExpression
setMethod(f          = "plotCoExpression",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){

            plotCoExpression(object[[SE.name]])
          })

### ---- plotCoExpressionProfile ----

#' @param cluster cluster number
#' @param condition Default is group.
#' @param features Default is NULL.
#' @section Plots:
#' \itemize{
#'    \item plotCoExpressionProfile:
#'    ...
#' }
#' @exportMethod plotCoExpressionProfile
#' @importFrom reshape2 melt
#' @rdname runCoExpression
#' @name plotCoExpressionProfile
#' @aliases plotCoExpressionProfile,RflomicsSE-method
setMethod(
  f = "plotCoExpressionProfile",
  signature = "RflomicsSE",
  definition = function(object,
                        cluster = 1,
                        condition="groups",
                        features=NULL){

    object <- getProcessedData(object, filter = TRUE)
    Groups <- getDesignMat(object)

    CoExpAnal <- getAnalysis(object, name = "CoExpAnal")$results

    coseq.res  <- CoExpAnal[["coseqResults"]]

    assays.data <- as.data.frame(coseq.res@assays@data[[1]])
    assays.data <- assays.data[assays.data[[paste0("Cluster_",cluster)]] > 0.8,]

    y_profiles.gg <- coseq.res@y_profiles[rownames(assays.data),] %>%
      data.frame() %>%
      mutate(observations=rownames(.)) %>%
      melt(id="observations", value.name = "y_profiles") %>%
      rename(samples = variable) %>%
      full_join(Groups , by = "samples")

    y_profiles.gg <- arrange(y_profiles.gg, get(condition))
    y_profiles.gg$groups <- factor(y_profiles.gg$groups,
                                   levels = unique(y_profiles.gg$groups))

    p <- ggplot(data = y_profiles.gg,
                aes(x = groups, y = y_profiles)) +
      geom_boxplot(aes(fill = condition),
                   outlier.size = 0.3) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      xlab("Conditions") + ylab("Expression profiles") +
      ggtitle(paste0("Cluster: ",cluster, "; nb_observations: ",
                     dim(assays.data)[1]))

    if(!is.null(features)){

      df <- filter(y_profiles.gg, observations == features) %>%
        group_by(groups) %>%
        summarise(mean.y_profiles=mean(y_profiles))
      p <- p +
        geom_point(data = df,
                   aes(x = groups, y = mean.y_profiles),
                   color = "red", size = 2) +
        geom_line( data = df,
                   aes(x = groups, y = mean.y_profiles),
                   color = "red", group = 1) +
        ggtitle(paste0("Cluster: ",cluster,
                       "; nb_observations  ",
                       dim(assays.data)[1], "; red: ", features))
    }
    return(p)
  })

#' @name plotCoExpressionProfile
#' @aliases plotCoExpressionProfile,RflomicsMAE-method
#' @rdname runCoExpression
#' @exportMethod plotCoExpressionProfile
setMethod(
  f          = "plotCoExpressionProfile",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name,
                        cluster = 1,
                        condition="groups",
                        features=NULL){

    plotCoExpressionProfile(object = object[[SE.name]],
                            cluster = cluster,
                            condition = condition,
                            features = features)

  })

### ---- plotCoseqContrasts ----
#' @section Plots:
#' \itemize{
#'    \item plotCoseqContrasts: This function describes the composition of
#'    clusters according to the contrast to which the gene belongs
#' }
#' @export
#' @exportMethod plotCoseqContrasts
#' @importFrom tidyr pivot_longer
#' @name plotCoseqContrasts
#' @aliases plotCoseqContrasts,RflomicsSE-method
#' @rdname runCoExpression
setMethod(
  f = "plotCoseqContrasts",
  signature = "RflomicsSE",
  definition = function(object){

    # Get Selected contrasts for coexpression
    H <- getCoexpSettings(object)$contrastNames

    if(is.null(H) || length(H) < 2)
      return(NULL)

    CoExpAnal <- getAnalysis(object, name = "CoExpAnal")

    # Gene's repartition by clusters
    coseq.res  <- CoExpAnal[["results"]][["coseqResults"]]
    genesByclusters <-
      as.data.frame(ifelse(coseq.res@assays@data[[1]] > 0.8, 1, 0))
    genesByclusters <- rownames_to_column(genesByclusters, var = "DEF")

    # Gene's repartition by Contrasts
    genesByContrasts <-
      getDEMatrix(object) |>
      select("DEF",as.vector(H))

    # Pivot tab.clusters
    genesByclusters.piv <- pivot_longer(
      data = genesByclusters,
      cols = seq(2, (dim(genesByclusters)[2])),
      names_to = "C",
      values_to = "does.belong") |>
      filter(does.belong == 1) |>
      select(-does.belong)

    # Merge the table of Cluster and Contrast
    tab <- left_join(genesByclusters.piv,
                     as.data.frame(genesByContrasts), by = "DEF") |>
      select(-DEF) |>
      group_by(C) |>
      mutate(sum = rowSums(across(names(genesByContrasts)[-1])))

    # Summarize repartition for specific genes
    tab.spe <- filter(tab, sum == 1) |>
      select(-sum) |>
      pivot_longer(names(genesByContrasts)[-1], names_to = "H") |>
      filter(value == 1) |>
      group_by(C, H) |>
      count()

    # Summarize repartition for genes common to at least
    # 2 contrasts
    tab.com <- filter(tab, sum > 1) |>
      select(-sum) |>
      pivot_longer(names(genesByContrasts)[-1], names_to = "H") |>
      filter(value == 1) |>
      group_by(C, H) |>
      count() |> ungroup() |>
      mutate(H = "common") |> distinct()

    # Bind table and add prop
    tab <- rbind(tab.com,tab.spe) |>
      group_by(C) |>
      mutate(prop=(n/sum(n))*100)

    p <-  ggplot(tab,  aes(x = C, y = prop, fill = H)) +
      geom_bar(stat ="identity") +
      coord_flip() +
      labs(x = "", y = paste0("Proportion of ",
                              .omicsDic(object)$variableNamegenes)) +
      scale_fill_discrete(name = "",
                          breaks = c("common", as.vector(H)),
                          labels = c("commons to at least 2 contrasts",
                                     paste0("specifics to ", names(H)))) +
      theme(legend.position="bottom",legend.direction = "vertical")

    return(p)
  })

#' @rdname runCoExpression
#' @name plotCoseqContrasts
#' @aliases plotCoseqContrasts,RflomicsMAE-method
#' @exportMethod plotCoseqContrasts
setMethod(f          = "plotCoseqContrasts",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){

            plotCoseqContrasts(object = object[[SE.name]])

          })


##==== ACCESSORS ====


### ---- Get Co-exp setting ----
#' @exportMethod getCoexpSettings
#' @section Accessors:
#' \itemize{
#'    \item getCoexpSettings: Access to the co-expression analysis settings
#'    of a given omic dataset
#' }
#' @name getCoexpSettings
#' @aliases getCoexpSettings,RflomicsSE-method
#' @rdname runCoExpression
setMethod(f          = "getCoexpSettings",
          signature  = "RflomicsSE",
          definition = function(object){
            return(getAnalysis(object, name = "CoExpAnal")$settings)
          })

#' @rdname runCoExpression
#' @exportMethod getCoexpSettings
#' @name getCoexpSettings
#' @aliases getCoexpSettings,RflomicsMAE-method
setMethod(f          = "getCoexpSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getCoexpSettings(object = object[[SE.name]])
          })

### ---- getCoexpClusters ----

#' @section Accessors:
#' \itemize{
#'    \item getCoexpClusters: get members of a cluster.
#'  return The list of entities inside this cluster.
#' }
#' @exportMethod getCoexpClusters
#' @param clusterName name of the cluster
#' @rdname runCoExpression
#' @name getCoexpClusters
#' @aliases getCoexpClusters,RflomicsSE-method
setMethod(
  f          = "getCoexpClusters",
  signature  = "RflomicsSE",
  definition = function(object, clusterName = NULL) {

    CoExpAnal <- getAnalysis(object, name = "CoExpAnal")

    if(length(CoExpAnal) == 0) return(NULL)

    Clusters <- CoExpAnal[["results"]]$clusters

    if(length(Clusters) == 0) return(NULL)

    if(is.null(clusterName)) return(Clusters)

    if(any(!clusterName %in% names(Clusters))) return(NULL)

    Clusters <- Clusters[clusterName]

    if(length(Clusters) == 1)
      return(Clusters[[1]])

    return(Clusters)
  })


#' @exportMethod getCoexpClusters
#' @rdname runCoExpression
#' @name getCoexpClusters
#' @aliases getCoexpClusters,RflomicsMAE-method
setMethod(
  f          = "getCoexpClusters",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name, clusterName = NULL) {
    getCoexpClusters(object = object[[SE.name]],
                       clusterName = clusterName)
  })
