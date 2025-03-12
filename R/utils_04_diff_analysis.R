### ============================================================================
### [04_diff_analysis] function and internal function
### ----------------------------------------------------------------------------
# D. Charif

# ---- Stat functions for differential analysis ----
## ---- edgeR ----


#' @title .edgeRAnaDiff
#'
#' @param object an object of class \link{RflomicsSE}
#' @param Contrasts.Coeff vector of coefficient, one for each contrast
#' @param FDR the FDR threshold
#' @param cmd if TRUE, verbose
#' @return A list of object of class \link{DGELRT}
#' @keywords internal
#' @importFrom stats model.matrix as.formula
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @noRd
#'
.edgeRAnaDiff <- function(object,
                          Contrasts.Coeff,
                          FDR = 1,
                          cmd = FALSE){

    modelFormula <- getModelFormula(object)
    if(length(modelFormula) == 0)
        stop("No model defined in the ", getDatasetNames(object), " object.")

    count_matrix <- assay(object)
    model_matrix <- model.matrix(as.formula(paste(modelFormula, collapse = " ")),
                                 data = getDesignMat(object))
    model_matrix <- model_matrix[colnames(object),]
    group        <- getCoeffNorm(object)$group
    lib.size     <- getCoeffNorm(object)$lib.size
    norm.factors <- getCoeffNorm(object)$norm.factors

    z <- y <- NULL

    ListRes <- list()

    # Construct the DGE obect
    dge <- DGEList(counts       = count_matrix,
                   group        = group,
                   lib.size     = lib.size,
                   norm.factors = norm.factors)

    # Run the model
    if (cmd) message("[RFLOMICS] [cmd] dge <- edgeR::estimateDisp(dge, design=model_matrix)")
    dge <- estimateDisp(dge, design=model_matrix)
    # if (cmd) message("[RFLOMICS] [cmd] dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)")
    # dge <- estimateGLMCommonDisp(dge, design=model_matrix)
    # if (cmd) message("[RFLOMICS] [cmd] dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)")
    # dge <- estimateGLMTrendedDisp(dge, design=model_matrix)
    # if (cmd) message("[RFLOMICS] [cmd] dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)")
    # dge <- estimateGLMTagwiseDisp(dge, design=model_matrix)
    if (cmd) message("[RFLOMICS] [cmd] fit.f <- edgeR::glmFit(dge,design=model_matrix)")
    fit.f <- glmFit(dge,design=model_matrix)


    if(cmd) message("[RFLOMICS] [cmd] apply model to each contrast")
    ResGlm <- lapply(rownames(Contrasts.Coeff), function(x){
        .tryRflomics(
            glmLRT(fit.f, contrast = unlist(Contrasts.Coeff[x,]))
        )
    })
    names(ResGlm) = rownames(Contrasts.Coeff)

    ListRes    <- list()
    error.list <- list()
    for(x in names(ResGlm)){

        if(!is.null(ResGlm[[x]]$error))
            error.list[[x]] <- ResGlm[[x]]$error

        if(!is.null(ResGlm[[x]]$value)){
            ListRes[["RawDEFres"]][[x]] <- ResGlm[[x]]$value
            res <- topTags(ResGlm[[x]]$value, n = dim(ResGlm[[x]]$value)[1])
            topDEF <- res$table[res$table$FDR <= FDR,]

            if(nrow(topDEF) == 0){

                error.list[[x]] <- "There are no results in table form."

            }else{
                ListRes[["DEF"]][[x]] <- topDEF
                ListRes[["DEF"]][[x]] <- dplyr::rename(topDEF,
                                                "Abundance"  = "logCPM",
                                                "pvalue"     = "PValue",
                                                "Adj.pvalue" = "FDR")
            }
        }
    }

    if(length(error.list) == 0)
        return(list(RawDEFres = ListRes[["RawDEFres"]], DEF = ListRes[["DEF"]]))
    else
        return(list(RawDEFres = ListRes[["RawDEFres"]], ErrorList = error.list))
}



## ---- limma ----

#' @title .limmaAnaDiff
#'
#' @param object an object of class \link{RflomicsSE}
#' @return A list
#' @keywords internal
#' @importFrom stats model.matrix as.formula
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @noRd
#'
.limmaAnaDiff <- function(object,
                          Contrasts.Coeff,
                          p.adj.cutoff = 1,
                          p.adj.method = "BH",
                          cmd = FALSE){

    count_matrix <- assay(object)
    modelFormula <- getModelFormula(object)
    if(length(modelFormula) == 0)
        stop("No model defined in the ", getDatasetNames(object), " object.")
    model_matrix <- model.matrix(as.formula(paste(modelFormula, collapse = " ")),
                                 data = getDesignMat(object))
    model_matrix <- model_matrix[colnames(object),]

    # Run the model
    if(cmd) message("[RFLOMICS] [cmd] fit linear model for each gene")
    fit <- lmFit(count_matrix, model_matrix)


    if(cmd) message("[RFLOMICS] [cmd] contrasts.fit(fit, contrasts = contrasts")

    ResGlm <-  lapply(rownames(Contrasts.Coeff), function(x){
        .tryRflomics(
            contrasts.fit(fit, contrasts = as.vector(unlist(Contrasts.Coeff[x,])))
        )
    })
    names(ResGlm) = rownames(Contrasts.Coeff)

    ListRes    <- list()
    error.list <- list()
    for(x in names(ResGlm)){
        # Construct a table of jobs summary

        if(!is.null(ResGlm[[x]]$error))
            error.list[[x]] <- ResGlm[[x]]$error

        if(!is.null(ResGlm[[x]]$value)){
            ListRes[["RawDEFres"]][[x]] <- ResGlm[[x]]$value

            fit2 <- eBayes(ResGlm[[x]]$value, robust=TRUE)
            res <- topTable(fit2, adjust.method = p.adj.method,
                            number=Inf, sort.by="AveExpr")
            topDEF <- res[res$adj.P.Val <= p.adj.cutoff,]

            if(nrow(topDEF) == 0){

                error.list[[x]] <- "There are no results in table form."

            }else{
                ListRes[["DEF"]][[x]] <- topDEF

                ListRes[["DEF"]][[x]] <- dplyr::rename(topDEF,
                                                "Abundance"="AveExpr",
                                                "pvalue"="P.Value",
                                                "Adj.pvalue"="adj.P.Val")
            }
        }
    }

    if(length(error.list) == 0)
        return(list(RawDEFres = ListRes[["RawDEFres"]], DEF = ListRes[["DEF"]]))
    else
        return(list(RawDEFres = ListRes[["RawDEFres"]], ErrorList = error.list))
}



# ---- Plot functions for differential analysis ----


#'.plotPValue
#'
#' @param data dataframe (ggplot2)
#' @param contrastName the contrast, useful for plot title
#' @return plot
#' @keywords internal
#' @noRd
.plotPValue <- function(data, contrastName = contrastName) {
    PValue <- NULL

    p <- ggplot(data = data) +
        geom_histogram(aes(x = pvalue), bins = 100) +
        labs(x = expression(p - value),
             y = "count",
             title = contrastName) +
        theme_bw(base_size = 10)

    return(p)
}



#' MA.plot
#'
#' @param data dataframe (ggplot2)
#' @param p.adj.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff |log2FC| cutoff (absolute value)
#' @param contrastName the contrast, useful for plot title
#' @return MA plot
#' @keywords internal
#' @importFrom ggpubr ggmaplot
#' @noRd

.plotMA <- function(data,
                    p.adj.cutoff,
                    logFC.cutoff,
                    contrastName = contrastName) {
    Abundance <- logFC <- Adj.pvalue <- NULL

    tmp <- select(data, "Abundance", "logFC", "Adj.pvalue") %>%
        dplyr::rename(.,
               baseMeanLog2 = Abundance,
               log2FoldChange = logFC,
               padj = Adj.pvalue
        )

    p <- ggmaplot(tmp,
                  main = contrastName,
                  fdr = p.adj.cutoff,
                  fc = 2 ^ logFC.cutoff,
                  size = 0.4,
                  ylab = bquote( ~ Log[2] ~ "fold change"),
                  xlab = bquote( ~ Log[2] ~ "mean expression"),
                  palette = c("#B31B21", "#1465AC", "grey30"),
                  select.top.method = c("padj", "fc"),
                  legend = "bottom",
                  top = 20,
                  font.label = c("plain", 7),
                  label.rectangle = TRUE,
                  font.legend = c(11, "plain", "black"),
                  font.main = c(11, "bold", "black"),
                  caption = paste(
                      "logFC cutoff=",
                      logFC.cutoff,
                      " and " ,
                      "FDR cutoff=",
                      p.adj.cutoff,
                      sep = ""
                  ),
                  ggtheme = theme_linedraw()
    )


    return(p)

}


#' Title
#'
#' @param data dataframe (ggplot2)
#' @param p.adj.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff log2FC cutoff (absolute value)
#' @param contrastName the contrast, useful for plot title
#' @return a volcano plot, made with the \link{EnhancedVolcano} package.
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @keywords internal
#' @noRd
#'
.plotVolcanoPlot <- function(data,
                             p.adj.cutoff,
                             logFC.cutoff,
                             contrastName) {
    # Find pvalue corresponding to the FDR cutoff for the plot
    # (mean between the last that passes the cutoff
    # and the first that is rejected to plot the line in the middle of
    #  the two points)
    # if pvalcutoff is 1 (no cutoff) no need to adjust

    if (p.adj.cutoff > 1) {
        stop("p.adj.cutoff must be between 0 and 1")
    }

    pval1 <- data$pvalue[data$Adj.pvalue < p.adj.cutoff] %>% last()
    pval2 <- data$pvalue[data$Adj.pvalue > p.adj.cutoff] %>% first()
    pvalCutoff <- (pval1 + pval2) / 2

    # If too low pvalues, unable to plot (error in if(d>0)...)
    # If drawconnectors is FALSE, it "works", with ylim being infinity,
    # it doesn't look like anything.
    # Modifiying the 0 pvalues to make sure it's working
    # default replacement in EnhancedVolcanoPlot
    nz_pval <- data$pvalue[data$pvalue != 0][1] * 10 ^ -1
    if (nz_pval == 0) {
        data$pvalue[data$pvalue == 0] <- data$pvalue[data$pvalue != 0][1]
    }

    Abundance <- logFC <- Adj.pvalue <- NULL
    p <- EnhancedVolcano(toptable = data,
                         lab = rownames(data),
                         x = 'logFC',
                         y = 'pvalue',
                         # pCutoff = p.adj.cutoff,
                         xlim = c(min(data[['logFC']], na.rm = TRUE) - 0.5,
                                  max(data[['logFC']], na.rm = TRUE) + 0.5),
                         pCutoff = pvalCutoff,
                         FCcutoff = logFC.cutoff,
                         axisLabSize = 14,
                         pointSize = 1.5,
                         labSize = 2.5,
                         title = contrastName,
                         titleLabSize = 20,
                         subtitle = "",
                         subtitleLabSize = 10,
                         caption = paste("logFC cutoff=",
                                         logFC.cutoff,
                                         " & " , "FDR cutoff=",
                                         p.adj.cutoff,
                                         sep = ""),
                         legendPosition = "bottom",
                         legendLabSize = 14,
                         legendIconSize = 3.5,
                         captionLabSize = 14,
                         col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                         colAlpha = 0.5,
                         drawConnectors = FALSE,
                         widthConnectors = 0.5,
                         max.overlaps = 15
    )

    return(p)
}



