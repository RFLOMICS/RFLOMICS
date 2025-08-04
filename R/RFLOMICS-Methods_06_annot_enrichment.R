### ============================================================================
### [06_annot_analysis] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# A. Hulot,
# N. Bessoltane

##==== STAT METHOD ====

###==== METHOD runAnnotationEnrichment using CLUSTERPROFILER ====

#' @title Run Gene Enrichment Analysis and process results
#' @description This function performs over representation analysis (ORA) using
#' clusterprofiler functions. It can be used with custom annotation file
#' (via enricher), GO (enrichGO) or KEGG (enrichKEGG) annotations.
#' @param object An object of class \link{RflomicsSE} or
#' class \link{RflomicsMAE-class}
#' @param SE.name SE.name the name of the dataset if the input object
#' is a \link{RflomicsMAE-class}
#' @param featureList name of contrasts (tags or names) from which to extract DE
#' genes if from is DiffExpAnal.
#' @param from indicates if ListNames are from differential analysis results
#' (DiffExp) or from the co-expression analysis results (CoExp)
#' @param universe description
#' @param database is it a custom annotation, GO or KEGG annotations
#' @param domain subcatgory for the database (eg BP for GO)
#' @param annotation for custom annotation, a data frame of the annotation.
#' The data frame must contains at least two columns: gene and term, with the
#' omics name and the associated term id respectively. A column name
#' can be added with the full name of the term (if  term is not the full name
#' already). The column domain can be used to indicate either different databases
#' (grouped analyses of kegg and go for example) or different domains for
#' a single database (CC, MF and BP) for GO.
#' @param OrgDb OrgDb (with enrichGO)
#' @param organism supported organism listed in
#' 'https://www.genome.jp/kegg/catalog/org_list.html' (with enrichKEGG)
#' @param keyType keytype of input gene with enrichGO
#' (one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot' with enrichKEGG)
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant.
#' Tests must pass
#' i) pvalueCutoff on unadjusted pvalues,
#' ii) pvalueCutoff on adjusted pvalues and
#' iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param ... additionnal parameters for enrichKEGG, enrichGO, enricher.
#' @return A RflomicsMAE or a RflomicsSE, depending on the class of object
#' parameter. The enrichment results are added to the metadata slot, either
#' in DiffExpEnrichAnal or CoExpEnrichAnal.
#' @importFrom tidyselect all_of
#' @importFrom clusterProfiler enrichKEGG enrichGO enricher
#' @importFrom utils getFromNamespace
#' @exportMethod runAnnotationEnrichment
#' @rdname runAnnotationEnrichment
#' @name runAnnotationEnrichment
#' @aliases runAnnotationEnrichment,RflomicsSE-method
#' @section Accessors:
#' A set of getters and setters generic functions to access and
#' modify objects of the slot metadata of a \link{RflomicsMAE-class} object or
#' a \link{RflomicsMAE-class} object.
#' @section Plots:
#' A collection of functions for plotting results from omics analysis steps.
#' @example inst/examples/runAnnotationEnrichment.R
setMethod(
    f = "runAnnotationEnrichment",
    signature = "RflomicsSE",
    definition = function(object,
                          featureList = NULL,
                          from = "DiffExp",
                          universe = NULL,
                          database = "custom",
                          domain = "no-domain",
                          annotation = NULL,
                          OrgDb     = NULL,
                          organism  = NULL,
                          keyType   = NULL,
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 1,
                          minGSSize = 10,
                          maxGSSize = 500,
                          ...) {

        # define result output
        EnrichAnal <- list(
            settings = list(),
            results  = list(),
            errors   = NULL
        )

        param.list <- list()
        if (is.null(pvalueCutoff)) pvalueCutoff <- 0.05
        param.list[["pvalueCutoff"]] <- pvalueCutoff

        if (is.null(qvalueCutoff)) qvalueCutoff <- 1
        param.list[["qvalueCutoff"]] <- qvalueCutoff

        if (is.null(minGSSize)) minGSSize <- 10
        param.list[["minGSSize"]] <- minGSSize

        if (is.null(maxGSSize)) maxGSSize <- 500
        param.list[["maxGSSize"]] <- maxGSSize

        if (is.null(universe))
            param.list[["universe"]] <- names(getProcessedData(object, filter = TRUE))

        # check args
        ## database
        if (is.null(database)) stop("The 'database' argument is required.")

        ## domain / clusterProfileR fucntion
        switch(
            database,
            "GO" = {
                param.list[["OrgDb"]]   <- OrgDb
                param.list[["keyType"]] <- keyType

                func_to_use <- "enrichGO"
                if (is.null(domain) || "ALL" %in% domain) domain <- c("MF", "BP", "CC")
                if (any(!domain %in% c("MF", "BP", "CC")))
                    stop("The domain parameter is required, and its value must be in the following list:
               MF, CC, BP, ALL")
            },
            "KEGG"   = {

                if (is.null(organism))
                    stop("Please provide an organism name to use KEGG database (eg: ath, hsa")
                if (is.null(keyType)) keytype <- "kegg"

                param.list[["organism"]] <- organism
                param.list[["keyType"]]  <- keyType

                func_to_use <- "enrichKEGG"
                domain      <- "no-domain"
            },
            "custom" = {
                func_to_use <- "enricher"

                if (is.null(annotation)) {
                    stop("You need an annotation file for a custom enrichment")
                }

                if (nrow(annotation) == 0) {
                    stop("The custom file is empty (0 rows)")
                }

                if (any(!c("term", "gene") %in% colnames(annotation)))
                    stop("The data.frame 'annotation' must contain these two columns:",
                         " 'gene' and 'term'. At least one of them is missing.")

                if ("domain" %in% colnames(annotation)) {
                    domain <- unique(annotation[["domain"]])
                    domain <- domain[!domain %in% c("", " ", ".")]
                    domain <- domain[!is.na(domain)]

                    if (length(domain) > 10) {
                        domain <- "no-domain"
                        message("more than 10 domains were detected, switching to no domain analysis")
                    }
                }else{ domain <- "no-domain" }

                TERM2GENE <- list()
                TERM2NAME <- list()
                if (setequal(unique(domain), c("no-domain"))) {
                    TERM2GENE[["no-domain"]] <- list(
                        "term" = annotation[["term"]],
                        "gene" = annotation[["gene"]])
                    TERM2NAME[["no-domain"]] <- NA
                    if ("name" %in% colnames(annotation)) {
                        TERM2NAME[["no-domain"]] <- list(
                            "term" = annotation[["term"]],
                            "name" = annotation[["name"]])
                    }
                }else{
                    TERM2GENE <- lapply(domain, function(x){
                        annot1 <- unique(filter(annotation, domain == x)[,c("term", "gene")])
                        list(
                            "term" = annot1$term,
                            "gene" = annot1$gene)
                    })
                    names(TERM2GENE) <- domain

                    if ("name" %in% colnames(annotation)) {
                        TERM2NAME <- lapply(domain, function(x){

                            annot2 <- unique(filter(annotation, domain == x)[,c("term", "name")])
                            list(
                                "term" = annot2$term,
                                "name" = annot2$name)
                        })
                        names(TERM2NAME) <- domain
                    }
                }
            })

        ## lists of features
        if (!is.null(featureList)) from <- .getOrigin(object, featureList)
        if (is.null(from)) from <- "DiffExp"

        if (!from %in% c("CoExp", "DiffExp"))
            stop("None of the lists correspond to lists of differential features or co-expression clusters.")

        ## reset EnrichAnal
        object <- setElementToMetadata(object,
                                       name = paste0(from, "EnrichAnal"),
                                       subName = database,
                                       content = NULL)

        switch(
            from,
            "DiffExp" = {

                DiffExpAnal <- getAnalysis(object, name = "DiffExpAnal")

                if (length(DiffExpAnal) == 0)
                    stop("There is no differential expression analysis.")

                if (is.null(featureList)) {

                    featureList <- getSelectedContrasts(object)$contrastName
                    if (!is.null(getValidContrasts(object)))
                        featureList <- getValidContrasts(object)$contrastName
                }

                geneLists <-
                    lapply(featureList, function(contrastName) {
                        row.names(DiffExpAnal[["results"]][["TopDEF"]][[contrastName]])
                    })
                names(geneLists) <- featureList

            },
            "CoExp" = {

                geneLists <- getCoexpClusters(object)

                if (is.null(geneLists)) stop("There is no co-expression analysis.")

                if (!is.null(featureList))  geneLists <- geneLists[featureList]
            }
        )
        geneLists <- geneLists[lengths(geneLists) > 0]

        # for each list
        results_list <- list()
        overview_list <- list()
        errorMessages <- list()
        for (listname in names(geneLists)) {

            param.list$gene <- geneLists[[listname]]

            for (dom in domain) {

                switch(
                    database,
                    "GO" = {
                        param.list$ont <- dom
                    },
                    "custom" = {
                        param.list$TERM2GENE <- TERM2GENE[[dom]]
                        if ("name" %in% colnames(annotation)) {
                            param.list$TERM2NAME <- TERM2NAME[[dom]]
                        } else {
                            param.list$TERM2NAME <- NULL
                        }
                    }
                )

                # run clusterProfileR
                catchRes <- .tryCatch_rflomics(
                    do.call(getFromNamespace(func_to_use, ns = "clusterProfiler"),
                            param.list))


                # delete heavy slots
                if (!is.null(catchRes$result)) {

                    res1 <- catchRes$result
                    slot(res1, name = "geneSets", check = FALSE) <- list("removed")
                    slot(res1, name = "universe", check = FALSE) <- c("removed")

                    res   <- res1@result
                    res.n <- nrow(res[res$p.adjust < pvalueCutoff,])

                    results_list[[listname]][[dom]]  <- res1
                    overview_list[[listname]][[dom]] <- res.n

                }else{
                    errorMessages[[listname]][[dom]] <-
                        c(catchRes$message, catchRes$warning, catchRes$error)
                    overview_list[[listname]][[dom]] <- NA
                    results_list[[listname]][[dom]]  <- NULL
                }
            }
        }

        # generate summay
        if (length(overview_list) == 0) {
            EnrichAnal[["results"]][["summary"]] <- NULL
        } else {

            dt_res <- as.data.frame(do.call("rbind", overview_list))
            for (y in names(dt_res)) {
                dt_res[[y]] <- unlist(dt_res[[y]])
            }

            if (any(!is.na(dt_res))) {
                if (from == "DiffExp") {
                    dt_res <-
                        dt_res %>% mutate(Contrast = rownames(.)) %>%
                        relocate(Contrast)
                }else{
                    dt_res <-
                        dt_res %>% mutate(Cluster = rownames(.)) %>%
                        relocate(Cluster)
                }
                EnrichAnal[["results"]][["summary"]] <- dt_res
            }
        }

        EnrichAnal[["errors"]] <- .CPR_message_processing(errorMessages)

        results_list <- Filter(Negate(is.null), results_list)

        storedParam <-
            c("universe", "keyType", "pvalueCutoff",
              "qvalueCutoff", "OrgDb", "organism",
              "minGSSize", "maxGSSize")

        EnrichAnal[["settings"]] <-
            param.list[names(param.list) %in% storedParam]

        if(database == "KEGG")
            EnrichAnal[["settings"]][["keggRelease"]] <- .getKEGGRelease()

        if(!is.null(annotation))
            EnrichAnal[["settings"]][["annotation"]] <- annotation

        EnrichAnal[["settings"]] <-
            c(EnrichAnal[["settings"]], list("domain" = domain))

        EnrichAnal[["results"]][["enrichResult"]] <- results_list

        # if (database == "KEGG") rm(kegg_category, envir = .GlobalEnv)
        object <- setElementToMetadata(object,
                                       name = paste0(from, "EnrichAnal"),
                                       subName = database,
                                       content = EnrichAnal)
        return(object)
    })

#' @rdname runAnnotationEnrichment
#' @name runAnnotationEnrichment
#' @aliases runAnnotationEnrichment,RflomicsMAE-method
#' @exportMethod runAnnotationEnrichment
setMethod(
    f = "runAnnotationEnrichment",
    signature = "RflomicsMAE",
    definition = function(object,
                          SE.name,
                          featureList = NULL,
                          from = "DiffExp",
                          universe = NULL,
                          database = "custom",
                          domain = "no-domain",
                          annotation = NULL,
                          OrgDb = NULL,
                          organism  = NULL,
                          keyType = NULL,
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 1,
                          minGSSize = 10,
                          maxGSSize = 500,
                          ...){

        object[[SE.name]] <-
            runAnnotationEnrichment(
                object       = object[[SE.name]],
                featureList  = featureList,
                from         = from,
                universe     = universe,
                database     = database,
                domain       = domain,
                annotation   = annotation,
                OrgDb        = OrgDb,
                keyType      = keyType,
                organism     = organism,
                pvalueCutoff = pvalueCutoff,
                qvalueCutoff = qvalueCutoff,
                minGSSize    = minGSSize,
                maxGSSize    = maxGSSize,
                ...)

        return(object)
    }
)

##==== GRAPHICAL METHOD ====

###==== plotClusterProfiler ====

#' @section Plots:
#' \itemize{
#'    \item plotClusterProfiler: Plot a dotplot, a cnetplot or an heatplot, using enrichplot
#' package. It is a wrapper method destined for the RflomicsSE class..
#' }
#' @param featureListName the name of the contrast to consider.
#' For Co expression analysis, it is expected to be one of
#' "cluster.1", "cluster.2", etc.
#' @param database the database (GO, KEGG or custom)
#' @param domain if the database is GO, expect one of BP, MF or CC.
#' Default is NULL.
#' @param plotType type of plot. Define the function used inside.
#' One of dotplot, heatplot or cnetplot.
#' @param showCategory max number of terms to show.
#' @param searchExpr expression to search in the showCategory terms.
#' @param nodeLabel same as in enrichplot::cnetplot function, defines
# the labels on the graph. One of "all", "category" or "gene". Default is 'all'.
#' @param p.adj.cutoff pvalueCutoff to define the enrichment threshold.
# Default is the one find in the clusterprofiler results in object.
#' @param ... additionnal parameters for cnetplot, heatplot or
#' enrichplot functions.
#' @importFrom enrichplot cnetplot heatplot dotplot set_enrichplot_color
#' @importFrom ggrepel geom_label_repel
#' @exportMethod plotClusterProfiler
#' @rdname runAnnotationEnrichment
#' @name plotClusterProfiler
#' @aliases plotClusterProfiler,RflomicsSE-method
setMethod(
    f = "plotClusterProfiler",
    signature = "RflomicsSE",
    definition = function(object,
                          featureListName = NULL,
                          database = NULL,
                          domain = "no-domain",
                          plotType = "dotplot",
                          showCategory = 15,
                          searchExpr = "",
                          nodeLabel = "all",
                          p.adj.cutoff = NULL,
                          ...) {

        if (is.null(featureListName))
            stop("Please provide a featureListName (contrast name or coexpression cluster name)")

        dataPlot <- getEnrichRes(
            object,
            featureListName = featureListName,
            database = database
        )

        log2FC_vect <- NULL
        switch(.getOrigin(object, featureListName),
               "DiffExp" = {

                   DiffExpAnal <- getAnalysis(object, name = "DiffExpAnal")
                   inter <-
                       DiffExpAnal[["results"]][["TopDEF"]][[featureListName]]
                   log2FC_vect <- inter[["logFC"]]
                   names(log2FC_vect) <- rownames(inter)
               },
               "CoExp" = {},
               stop(
                   "Argument from is detected to be neither DiffExp nor CoExp.")

        )

        if (is.null(p.adj.cutoff)) {
            p.adj.cutoff <-
                getEnrichSettings(object,
                                  from = .getOrigin(object, featureListName),
                                  database = database)$pvalueCutoff
        }

        if (database == "GO") {
            if (is.null(domain)) {
                stop("database is GO, non null domain (BP, CC or MF) is expected.")
            } else if (setequal(domain, "no-domain")) {
                stop("database is GO, a domain name (BP, CC or MF) is expected")
            } else {
                dataPlot <- dataPlot[[domain]]
            }
        } else if (database == "custom" || database == "KEGG") {
            if (is.null(domain)) {
                domain <- "no-domain"
            }

            if (!domain %in% names(dataPlot)) {
                stop("domain is expected to be one of ",
                     paste(names(dataPlot), collapse = ","))
            } else {
                dataPlot <- dataPlot[[domain]]
            }
        }

        # Select categories to show
        dataTab <-
            dataPlot@result[which(dataPlot@result$p.adjust < p.adj.cutoff),]
        Categories <- dataTab$Description
        if (searchExpr != "") {
            Categories <-
                Categories[grep(toupper(searchExpr), toupper(Categories))]
        }
        NbtoPlot <- min(length(Categories), showCategory)
        if (NbtoPlot == 0)
            stop("There is no terms to show")

        Categories <- Categories[seq_len(NbtoPlot)]

        # Create the plot
        plotType <- tolower(plotType)
        returnplot <- NULL

        returnplot <-  switch(
            plotType,
            "cnetplot" = {
                cnetplot(
                    dataPlot,
                    showCategory = Categories,
                    color.params = list(foldChange = log2FC_vect),
                    node_label = nodeLabel,
                    ...) +
                    guides(colour = guide_colourbar(title = "log2FC"))
            },
            "heatplot" = {
                outgg <- heatplot(dataPlot,
                                  showCategory = Categories,
                                  foldChange = log2FC_vect,
                                  ...)
                outgg$scales$scales <- list()
                outgg + labs(fill = "log2FC") +
                    scale_fill_gradient2(
                        low = "blue",
                        mid = "white",
                        high = "red",
                        midpoint = 0,
                    ) +
                    theme(axis.text.y = element_text(size = 10))
            },
            {
                tryCatch(
                    dotplot(dataPlot,
                            showCategory = Categories,
                            ...),
                    error = function(e)
                        e,
                    warnings = function(w)
                        w
                )
            })

        return(returnplot)
    }
)


###==== plotEnrichComp ====

#' @section Plots:
#' \itemize{
#'    \item plotEnrichComp: plot an heatmap of all the enriched term found for a given
#' database and a given source (differential analysis or coexpression clusters).
#' Allow for the comparison of several enrichment results.
#' }
#' @param from indicates if the enrichment results are taken from differential
#' analysis results (DiffExp) or from the co-expression analysis results
#'  (CoExp)
#' @param database is it a custom annotation, GO or KEGG annotations
#' @param domain domain from the database (eg GO has three domains,
#' BP, CC and MF)
#' @param matrixType Heatmap matrix to plot, one of GeneRatio, p.adjust or
#' presence.
#' @param nClust number of separate cluster to plot on the heatmap, based o
#' n the clustering.
#' @param ... more arguments for ComplexHeatmap::Heatmap.
#' @importFrom reshape2 recast
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap ht_opt
#' @importFrom stats hclust dist
#' @importFrom stringr str_wrap
#' @exportMethod plotEnrichComp
#' @rdname runAnnotationEnrichment
#' @name plotEnrichComp
#' @aliases plotEnrichComp,RflomicsSE-method
setMethod(
    f = "plotEnrichComp",
    signature = "RflomicsSE",
    definition = function(object,
                          from = "DiffExp",
                          database = NULL,
                          domain = "no-domain",
                          matrixType = "FC",
                          nClust = NULL,
                          ...) {

        if (is.null(from))
            stop("The 'from' parameter is required.")

        if (!from %in% c("CoExp", "DiffExp"))
            stop("The 'from' parameter must be one of CoExp or DiffExp.")

        if (is.null(database))
            stop("The 'database' parameter is required.")

        allData <- getEnrichRes(object, from = from, database = database)

        if (length(allData) == 0) {
            stop("The selected database ", database,
                 " does not seem to exist in the object.")
        }

        allData <- allData[lengths(allData) > 0]
        if (length(allData) == 0) {
            stop("It seems there is no results here.")
        }

        pvalThresh <- allData[[1]][[1]]@pvalueCutoff

        domainPoss <- unique(unlist(lapply(allData, names)))
        if (missing(domain)) domain <- domainPoss
        if (is.null(domain))   domain <- domainPoss
        if (any(!domain %in% domainPoss)) {
            stop("Trying to select a domain that does not exist in the object.")
        }

        extract <-
            do.call("rbind", lapply(
                seq_len(length(allData)),
                FUN = function(i) {
                    allData[[i]] <- Filter(Negate(is.null), allData[[i]])
                    if (domain %in% names(allData[[i]])) {
                        cprRes <- allData[[i]][[domain]]
                        selectterms <- cprRes@result$p.adjust < cprRes@pvalueCutoff
                        nterms <- sum(selectterms)
                        if (nterms > 0) {
                            cprRes <- cprRes@result[selectterms, ]
                            cprRes$contrastName <- names(allData)[i]
                            cprRes$domain <- domain
                            return(cprRes)
                        }}
                }
            ))

        toKeep <- names(which(table(extract$ID) > 1))
        if (length(toKeep) == 0)
            stop("There is no common terms to show.")
        extract <- extract[extract$ID %in% toKeep, ]

        # handling description and ID
        extract$Description <-
            switch(
                database,
                "KEGG" = {
                    gsub(" - .*", "", extract$Description)
                },
                {
                    # if there is duplication in description
                    # not necessarily duplicated in the ID
                    # should not be grouped in the heatmap
                    if (sum(duplicated(extract$Description)) > 0) {
                        if (!identical(extract$Description,
                                       extract$ID)) {
                            posDup <-
                                unique(which(
                                    duplicated(extract$Description),
                                    duplicated(extract$Description, fromLast = TRUE)
                                ))
                            extract$Description[posDup] <-
                                paste0("(",
                                       extract$ID[posDup],
                                       ")",
                                       "\n",
                                       extract$Description[posDup])
                            extract$Description
                        } else{
                            extract$Description
                        }

                    } else {
                        extract$Description
                    }
                })

        extract$Description <-
            str_wrap(extract$Description, width = 20)
        extract$contrastName <-
            str_wrap(extract$contrastName, width = 30)

        extract$GeneRatio <- as.numeric(vapply(
            extract$GeneRatio,
            FUN = function(x)
                eval(parse(text = x)),
            FUN.VALUE = 1
        ))
        extract$BgRatio <- as.numeric(vapply(
            extract$BgRatio,
            FUN = function(x)
                eval(parse(text = x)),
            FUN.VALUE = 1
        ))
        extract$FC <- extract$GeneRatio / extract$BgRatio

        dat <- switch(
            matrixType,
            "GeneRatio" = {
                inter <- recast(
                    extract[, c("Description", "contrastName", "GeneRatio")],
                    Description ~ contrastName,
                    measure.var = "GeneRatio")
                rownames(inter) <- inter$Description
                # inter <- inter[, colnames(inter) != "ID"]
                inter <- inter[, colnames(inter) != "Description"]
                inter[is.na(inter)] <- 0
                inter
            },
            "p.adjust" = {
                inter <- recast(
                    extract[, c("Description", "contrastName", "p.adjust")],
                    Description ~ contrastName,
                    measure.var = "p.adjust")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "Description"]
                inter[is.na(inter)] <- 1
                inter
            },
            "presence" = {
                inter <- recast(
                    extract[, c("Description", "contrastName", "p.adjust")],
                    Description ~ contrastName,
                    measure.var = "p.adjust")
                rownames(inter) <- inter$Description
                inter <-
                    inter <- inter[, colnames(inter) != "Description"]
                inter[!is.na(inter)] <- 1
                inter[is.na(inter)]  <- 0
                inter
            },
            "FC" = {
                inter <- recast(
                    extract[, c("Description", "contrastName", "FC")],
                    Description ~ contrastName,
                    measure.var = "FC")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "Description"]
                inter[is.na(inter)] <- 0
                inter
            },
            "log2FC" = {
                inter <- recast(
                    extract[, c("Description", "contrastName", "FC")],
                    Description ~ contrastName,
                    measure.var = "FC")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "Description"]
                inter <- log2(inter)
                inter[is.infinite(as.matrix(inter))] <- 0
                # means FC is 0, shouldn't happen much...
                inter[is.na(inter)] <- 0
                # means it's not significant and not in the matrix.
                inter
            },
            stop("This matrix type is not supported. Please chose one of log2FC, FC, presence, p.adjust, GeneRatio")
        )

        if (nrow(dat) > 1) {
            switch(matrixType,
                   "presence" = {
                       hcPlot <- hclust(dist(dat, method = "binary"),
                                        method = "complete")
                       hcCol <-
                           hclust(dist(t(dat), method = "binary"),
                                  method = "complete")
                   },
                   {
                       hcPlot <- hclust(dist(dat, method = "euclidean"),
                                        method = "complete")
                       hcCol <-
                           hclust(dist(t(dat), method = "euclidean"),
                                  method = "complete")
                   })
        } else {
            hcPlot <- FALSE
        }

        colors <- switch(
            matrixType,
            "presence"   = {
                structure(c("white", "firebrick"), names = c("0", "1"))
            },
            "GeneRatio"  = {
                colorRamp2(c(0, max(dat)), c("white", "firebrick"))
            },
            "p.adjust"   = {
                colorRamp2(c(0, pvalThresh, 1),
                           c("firebrick", "white", "white"))
            },
            "FC"         = {
                colorRamp2(c(0, max(dat)), c("white", "firebrick"))
            },
            "log2FC"     = {
                colorRamp2(c(-max(abs(dat)), 0, max(abs(dat))),
                           c("blue", "white", "firebrick"))
            }
        )

        ht_opt(DENDROGRAM_PADDING = unit(0.1, "cm"))

        suppressWarnings(
            Heatmap(
                dat,
                col = colors,
                name = matrixType,
                cluster_columns = hcCol,
                show_column_dend = FALSE,
                cluster_rows = hcPlot,
                row_names_side = "left",
                column_names_rot = 30,
                rect_gp = gpar(col = "gray50", lwd = 0.5),
                width =  ncol(dat) * 5,
                height = nrow(dat) * 5,
                heatmap_legend_param = list(direction = "horizontal"),
                border = TRUE,
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = 10)
            )
        )

    }
)

##==== ACCESSORS ====

### ---- Get a particular enrichment result ----
#' @section Accessors:
#' \itemize{
#'    \item getEnrichRes: get a particular enrichment result.
#'    return enrichment results given in the form of lists of clusterprofiler
#'    results.
#' }
#' @param featureListName the contrast or cluster name on which the enrichment
#' was perform.
#' @param experiment if the object is a RflomicsMAE, then experiment is the
#' name of the RflomicsSE to look for.
#' @param from either diffexp or coexp, indicates where to search for the
#' results
#' @param database the database used for the enrichment (GO, KEGG or custom)
#' @param domain the subonology or subdomain for the database (eg CC, MF or
#' BP for GO.)
#' @param ... Not in use at the moment
#' @exportMethod getEnrichRes
#' @rdname runAnnotationEnrichment
#' @name getEnrichRes
#' @aliases getEnrichRes,RflomicsSE-method
setMethod(
    f          = "getEnrichRes",
    signature  = "RflomicsSE",
    definition = function(object,
                          featureListName = NULL,
                          from = "DiffExp",
                          database = "GO",
                          domain = NULL) {

        if (is.null(database) )
            stop("The 'database' parameter is required. It must be one of GO, KEGG or custom.")

        if (!database %in% c("GO", "KEGG", "custom"))
            stop("The 'database' parameter must be one of GO, KEGG or custom.")

        if (!is.null(featureListName)) from <- .getOrigin(object, featureListName)

        if (is.null(from))
            stop("The 'from' parameter is required. It must be one of DiffExp or CoExp")

        if (!from %in% c("CoExp", "DiffExp"))
            stop("The 'from' parameter must be one of DiffExp or CoExp")

        res <- getAnalysis(object, name = paste0(from, "EnrichAnal"),
                           subName = database)

        res <- res[["results"]][["enrichResult"]]

        if (!is.null(featureListName)){
            res <- res[[featureListName]]

            if (!is.null(domain))
                res <- res[[domain]]
        }

        return(res)
    })

#' @rdname runAnnotationEnrichment
#' @name getEnrichRes
#' @aliases getEnrichRes,RflomicsMAE-method
#' @exportMethod getEnrichRes
setMethod(
    f = "getEnrichRes",
    signature = "RflomicsMAE",
    definition = function(object,
                          experiment,
                          featureListName = NULL,
                          from = "DiffExp",
                          database = "GO",
                          domain = NULL) {

        if (missing(experiment)) {
            stop("Please indicate from which data you want to extract the enrichment results.")
        }

        res_return <- getEnrichRes(
            object[[experiment]],
            featureListName = featureListName,
            from = from,
            database = database,
            domain = domain
        )
    })

### ---- Get summary from ORA : ----

#' @section Accessors:
#' \itemize{
#'    \item sumORA: Get summary tables from ORA analyses -
#'    once an enrichment has been conducted.
#' }
#' @param database either NULL, GO, KEGG or custom.
#' if NULL, all tables are returned in a list.
#' @param from either DiffExp or CoExp.
#' @param featureListName the contrastName or clusterName to retrieve
#' the results from. If NULL, all results are returned.
#' @return a list of tables or a table
#' @exportMethod sumORA
#' @rdname runAnnotationEnrichment
#' @name sumORA
#' @aliases sumORA,RflomicsSE-method
setMethod(
    f = "sumORA",
    signature = "RflomicsSE",
    definition = function(object,
                          from = "DiffExp",
                          database = NULL,
                          featureListName = NULL) {

        if (!from %in% c("CoExp", "DiffExp"))
            stop("from argument must be one of DiffExp or CoExp.")

        listnames <- switch(from,
                            "DiffExp" = "Contrast",
                            "CoExp" = "Cluster")

        EnrichAnal <- getAnalysis(object,  name = paste0(from, "EnrichAnal"))
        if (is.null(database)) database <- names(EnrichAnal)

        list_res <-
            lapply(database, function(ontres) {
                interRes <- EnrichAnal[[ontres]]$results$summary
                if (!is.null(featureListName)) {
                    interRes <- interRes[which(interRes[[listnames]] == featureListName),]
                }
                interRes
            })
        names(list_res) <- database

        if (length(list_res) == 1) return(list_res[[1]])
        return(list_res)
    })

### ---- Get an enrichment arguments ----

#' @section Accessors:
#' \itemize{
#'    \item getEnrichSettings:
#'    get the settings of an enrichment analysis.
#' }
#' @param from where to search for the results (either CoExp or DiffExp)
#' @param database which database (GO, KEGG, custom...)
#' @return a list with all settings
#' @exportMethod getEnrichSettings
#' @rdname runAnnotationEnrichment
#' @name getEnrichSettings
#' @aliases getEnrichSettings,RflomicsSE-method
setMethod(
    f = "getEnrichSettings",
    signature = "RflomicsSE",
    definition = function(object,
                          from = "DiffExp",
                          database = "GO") {

        if (!database  %in% c("GO", "KEGG", "custom")) {
            stop(database, " is not a valid value. Choose one of GO, KEGG or custom.")
        }

        if (!from %in% c("CoExp", "DiffExp"))
            stop("The from parameter must be one of DiffExp or CoExp, not ", from)

        return(metadata(object)[[paste0(from, "EnrichAnal")]][[database]]$settings)
    }
)

### ---- getAnnotAnalysesSummary ----

#' @section Accessors:
#' \itemize{
#'    \item getAnnotAnalysesSummary:
#'    return A list of heatmaps, one for each ontology/domain.
#' }
#' @param from indicates if the enrichment results are taken from differential
#' analysis results (DiffExp) or from the co-expression analysis
#' results (CoExp)
#' @param matrixType Heatmap matrix to plot, one of GeneRatio, p.adjust
#' or presence.
#' @param ... more arguments for ComplexHeatmap::Heatmap.
#' @importFrom reshape2 recast
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap ht_opt
#' @importFrom stringr str_wrap
#' @exportMethod getAnnotAnalysesSummary
#' @rdname runAnnotationEnrichment
#' @name getAnnotAnalysesSummary
#' @aliases getAnnotAnalysesSummary,RflomicsMAE-method
setMethod(
    f = "getAnnotAnalysesSummary",
    signature = "RflomicsMAE",
    definition = function(object,
                          from = "DiffExp",
                          matrixType = "presence",
                          ...) {
        extract.list <- list()

        if (!matrixType %in% c("presence", "FC", "log2FC", "p.adjust", "GeneRatio"))
            stop("matrixType ", matrixType, " is not recognize. Please chose one of:",
                 " presence, FC, log2FC, p.adjust or GeneRatio")

        if (!from %in% c("CoExp", "DiffExp"))
            stop("The from parameter must be one of DiffExp or CoExp, not ", from)

        from <- paste0(from, "EnrichAnal")

        omicNames <- unique(unlist(getAnalyzedDatasetNames(object, from)))

        for (data in omicNames) {

            # for each database
            EnrichAnal <- getAnalysis(object[[data]], name = from)
            databases <- names(EnrichAnal)
            for (database in databases) {
                if (is.null(EnrichAnal[[database]]$results$enrichResult))
                    next

                clusterNames <-
                    names(EnrichAnal[[database]]$results$enrichResult)
                pvalThresh   <-
                    EnrichAnal[[database]]$settings$pvalueCutoff

                for (name in clusterNames) {
                    domains <-
                        names(EnrichAnal[[database]]$results$enrichResult[[name]])

                    for (dom in domains) {
                        cprRes <-
                            EnrichAnal[[database]]$results$enrichResult[[name]][[dom]]

                        if (is.null(cprRes)) next

                        cprRes <-
                            cprRes@result[cprRes@result$p.adjust < cprRes@pvalueCutoff, ]

                        if (nrow(cprRes) == 0) next

                        cprRes$contrastName <- name
                        cprRes$dataset      <- data

                        extract.list[[database]][[dom]] <-
                            rbind(extract.list[[database]][[dom]], cprRes)
                    }
                }
            }
        }


        p.list <- list()

        for (database in names(extract.list)) {
            for (dom in names(extract.list[[database]])) {
                extract <- extract.list[[database]][[dom]]

                extract$Description <-
                    str_wrap(extract$Description, width = 30)
                extract$contrastName <-
                    str_wrap(extract$contrastName, width = 30)

                extract$GeneRatio <-
                    as.numeric(vapply(
                        extract$GeneRatio,
                        FUN = function(x)
                            eval(parse(text = x)),
                        FUN.VALUE = 1
                    ))
                extract$BgRatio <-
                    as.numeric(vapply(
                        extract$BgRatio,
                        FUN = function(x)
                            eval(parse(text = x)),
                        FUN.VALUE = 1
                    ))
                extract$FC <- extract$GeneRatio / extract$BgRatio

                extract$contrastNameLabel <- extract$contrastName
                extract$contrastName <-
                    paste(extract$contrastName, extract$dataset, sep = "\n")

                split.df <- unique(extract[c("dataset", "contrastName", "contrastNameLabel")])
                split <- split.df$dataset
                names(split) <- split.df$contrastName

                if (nrow(extract) == 0)
                    next

                dat <- switch(
                    matrixType,
                    "GeneRatio" = {
                        inter <- recast(extract[, c("Description", "contrastName", "GeneRatio")],
                                        Description ~ contrastName,
                                        measure.var = "GeneRatio")
                        rownames(inter) <- inter$Description
                        #inter <- inter[, colnames(inter) != "Description"]
                        inter <- select(inter, -"Description")
                        inter[is.na(inter)] <- 0
                        inter
                    },
                    "p.adjust" = {
                        inter <-
                            recast(extract[, c("Description", "contrastName", "p.adjust")],
                                   Description ~ contrastName,
                                   measure.var = "p.adjust")
                        rownames(inter) <- inter$Description
                        #inter <- inter[, colnames(inter) != "Description"]
                        inter <- select(inter, -"Description")
                        inter[is.na(inter)] <- 1
                        inter
                    },
                    "presence" = {
                        inter <- recast(extract[, c("Description", "contrastName", "p.adjust")],
                                        Description ~ contrastName,
                                        measure.var = "p.adjust")
                        rownames(inter) <- inter$Description
                        #inter <- inter[, colnames(inter) != "Description"]
                        inter <- select(inter, -"Description")
                        inter[!is.na(inter)] <- 1
                        inter[is.na(inter)]  <- 0
                        inter
                    },
                    "FC" = {
                        inter <- recast(extract[, c("Description", "contrastName", "FC")],
                                        Description ~ contrastName,
                                        measure.var = "FC")
                        rownames(inter) <- inter$Description
                        #inter <- inter[, colnames(inter) != "Description"]
                        inter <- select(inter, -"Description")
                        inter[is.na(inter)] <- 0
                        inter
                    },
                    "log2FC" = {
                        inter <- recast(extract[, c("Description", "contrastName", "FC")],
                                        Description ~ contrastName,
                                        measure.var = "FC")
                        rownames(inter) <- inter$Description
                        #inter <- inter[, colnames(inter) != "Description"]
                        inter <- select(inter, -"Description")
                        inter <- log2(inter)
                        inter[is.infinite(as.matrix(inter))] <- 0
                        # means FC is 0, shouldn't happen much...
                        inter[is.na(inter)] <- 0
                        # means it's not significant and not in the matrix.
                        inter
                    }
                )

                colors <- switch(
                    matrixType,
                    "presence"   = {
                        structure(c("white", "firebrick"), names = c("0", "1"))
                    },
                    "GeneRatio"  = {
                        colorRamp2(c(0, max(dat)), c("white", "firebrick"))
                    },
                    "p.adjust"   = {
                        colorRamp2(c(0, pvalThresh, 1),
                                   c("firebrick", "white", "white"))
                    },
                    "FC"         = {
                        colorRamp2(c(0, max(dat)), c("white", "firebrick"))
                    },
                    "log2FC"     = {
                        colorRamp2(
                            c(-max(abs(dat)), 0, max(abs(dat))),
                            c("blue", "white", "firebrick"))
                    }
                )

                ht_opt(DENDROGRAM_PADDING = unit(0.1, "cm"))

                p.list[[database]][[dom]] <- suppressWarnings(
                    Heatmap(
                        t(dat),
                        col = colors,
                        name = matrixType,
                        row_split = split[names(dat)],
                        #cluster_columns = hcCol,
                        #cluster_rows = hcPlot,
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        row_names_side = "left",
                        #column_names_rot = 20,
                        row_labels = split.df[which(split.df$contrastName %in% names(dat)), ]$contrastNameLabel,
                        # column_names_rot = 0,
                        # column_names_centered = TRUE,
                        rect_gp = gpar(col = "gray80", lwd = 0.1),
                        width =  ncol(dat) * 5,
                        height = nrow(dat) * 5,
                        heatmap_legend_param = list(direction = "horizontal"),
                        border = TRUE,
                        column_names_gp = gpar(fontsize = 10),
                        row_names_gp = gpar(fontsize = 10)
                    )
                )


            }
        }
        if (length(p.list) == 0)
            return(NULL)
        return(p.list)
    }
)