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
#' @importFrom dplyr rename
#' @noRd
#'
.edgeRAnaDiff <- function(object,
                          modelFormula = NULL, 
                          Contrasts.Coeff,
                          FDR = 1,
                          cmd = FALSE){

    if(is.null(modelFormula))
      modelFormula <- getModelFormula(object)
  
    if(length(modelFormula) == 0)
        stop("No model defined in the ", getDatasetNames(object), " object.")

    count_matrix <- assay(object)
    model_matrix <- model.matrix(as.formula(paste(modelFormula, collapse = " ")),
                                 data = getDesignMat(object))
    model_matrix <- model_matrix[colnames(object),]
    #group        <- getCoeffNorm(object)$group
    #lib.size     <- getCoeffNorm(object)$lib.size
    #norm.factors <- getCoeffNorm(object)$norm.factors

    target       <- getDesignMat(object)
    #coeffNorm    <- getCoeffNorm(object)
    group        <- target$groups
    #lib.size     <- coeffNorm[coeffNorm$group %in% group,]$lib.size
    #norm.factors <- coeffNorm[coeffNorm$group %in% group,]$norm.factors
    
    z <- y <- NULL

    ListRes <- list()

    # Construct the DGE obect
    dge <- DGEList(counts       = count_matrix,
                   group        = group,
                   #lib.size     = lib.size,
                   norm.factors = rep(1, ncol(count_matrix)))

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
                ListRes[["DEF"]][[x]] <- rename(topDEF,
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
                          modelFormula = NULL, 
                          Contrasts.Coeff,
                          p.adj.cutoff = 1,
                          p.adj.method = "BH",
                          cmd = FALSE){

    if(is.null(modelFormula))
      modelFormula <- getModelFormula(object)
  
    count_matrix <- assay(object)
    
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
                  caption = bquote(log[2]*FC~cutoff~.(logFC.cutoff)~and~FDR~cutoff~.(p.adj.cutoff)),
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
#' @return a volcano plot
#' @importFrom ggrepel geom_text_repel
#' @keywords internal
#' @noRd
#'
.plotVolcanoPlot <- function(data,
                             p.adj.cutoff,
                             logFC.cutoff,
                             contrastName) {
    if (p.adj.cutoff > 1) {
        stop("p.adj.cutoff must be between 0 and 1")
    }

    data <- data[order(data[["pvalue"]], decreasing = FALSE),]

    pval1 <- data[["pvalue"]][data[["Adj.pvalue"]] < p.adj.cutoff][sum(data[["Adj.pvalue"]] < p.adj.cutoff)]
    pval2 <- data[["pvalue"]][data[["Adj.pvalue"]] > p.adj.cutoff][1]
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

    data[["Entity"]] <- rownames(data)
    data[["log2FC"]] <- data[["logFC"]]
    data[["-log10pvalue"]] <- -log10(data[["pvalue"]])
    data[["criteria"]] <- "none"
    data[["criteria"]][abs(data[["log2FC"]]) > logFC.cutoff] <- "log2FC"
    data[["criteria"]][data[["Adj.pvalue"]] < p.adj.cutoff] <- "Adjusted Pvalue"
    data[["criteria"]][data[["Adj.pvalue"]] < p.adj.cutoff & abs(data[["log2FC"]]) > logFC.cutoff] <- "Both"
    data[["Entity"]][data[["criteria"]] != "Both"] <- NA

    data <- data[order(data[["criteria"]]),]

    if (sum(data[["criteria"]] == "Both") > 2000) {
        data[["Entity"]][seq(from = 2000, to = nrow(data))] <- NA
    }

    data[["criteria"]][data[["criteria"]] == "Both"] <- "log2FC and Adj. Pvalue"

    p <- ggplot(data, aes(x = log2FC, y = `-log10pvalue`,
                          color = criteria, label = Entity)) +
        geom_point(size = 1.5, alpha = 0.5) +
        xlim(c(min(data[['logFC']], na.rm = TRUE) - 0.5,
               max(data[['logFC']], na.rm = TRUE) + 0.5)) +
        theme_bw() +
        xlab(bquote(~Log[2]~ "Fold Change")) +
        ylab(bquote(~-Log[10]~ "Pvalue")) +
        labs(caption = bquote(log[2]*FC~cutoff~.(logFC.cutoff)~and~FDR~cutoff~.(p.adj.cutoff)),
             title = contrastName) +
        theme(legend.position = "bottom",
              plot.subtitle = element_text(size = 10),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.title = element_blank()) +
        scale_color_manual(values = c("none" = 'grey30',
                                      "log2FC" = 'forestgreen',
                                      "Adjusted Pvalue" = 'royalblue',
                                      "log2FC and Adj. Pvalue" = 'red2')) +
        geom_vline(xintercept = logFC.cutoff, lty = "dashed", color = "gray28") +
        geom_vline(xintercept = -logFC.cutoff, lty = "dashed", color = "gray28") +
        geom_vline(xintercept = 0, lty = "dashed", color = "gray28") +
        geom_hline(yintercept = -log10(pvalCutoff), lty = "dashed", color = "gray28") +
        geom_text_repel(show.legend = FALSE, na.rm = TRUE)
    
    return(p)
}



#' fast_MA.plot
#'
#' @param data dataframe 
#' @param p.adj.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff |log2FC| cutoff (absolute value)
#' @param contrastName the contrast, useful for plot title
#' @return MA plot
#' @keywords internal
#' @noRd

.plotMA_fast <- function(data, p.adj.cutoff, logFC.cutoff, contrastName) {
  
  x <- data$Abundance
  y <- data$logFC
  p <- data$Adj.pvalue
  
  fc_thr <- logFC.cutoff
  sig <- p < p.adj.cutoff
  
  col <- rep("grey30", length(y))
  col[sig & y > fc_thr] <- "forestgreen"
  col[sig & y < -fc_thr] <- "red2"
  
  par(mar = c(7, 4, 4, 2) + 0.1)
  
  graphics::plot(
    x, y,
    pch = 16,
    cex = 0.5,
    col = col,
    main = contrastName,
    xlab = "Log2 mean expression",
    ylab = "Log2 fold change"
  )
  
  mtext(bquote(log[2]*FC~cutoff~.(logFC.cutoff)~and~FDR~cutoff~.(p.adj.cutoff)),
        side = 1,  
        line = 5,  
        adj = 0.5, 
        cex = 0.9)
  
  graphics::abline(h = c(-fc_thr, fc_thr), col = "black", lty = 2)
  
  ord <- order(p)
  top_idx <- ord[1:20]
  
  graphics::text(
    x[top_idx],
    y[top_idx],
    labels = rownames(data)[top_idx],
    cex = 0.5,
    pos = 3
  )
}



#' plotVolcano_fast
#'
#' @param data dataframe 
#' @param p.adj.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff log2FC cutoff (absolute value)
#' @param contrastName the contrast, useful for plot title
#' @return a volcano plot
#' @keywords internal
#' @noRd
#'
.plotVolcano_fast <- function(data,
                              p.adj.cutoff = 0.05,
                              logFC.cutoff = 1,
                              contrastName = "",
                              n_label = 20) {
  
  
  x <- data$logFC
  p <- data$pvalue
  padj <- data$Adj.pvalue
  
  ok <- is.finite(x) & is.finite(p) & is.finite(padj)
  
  x <- x[ok]
  p <- p[ok]
  padj <- padj[ok]
  rn <- rownames(data)[ok]
  
  # If too low pvalues, unable to plot (error in if(d>0)...)
  p[p == 0] <- min(p[p > 0], na.rm = TRUE) * 0.1
  
  y <- -log10(p)
  
  
  # significance (padj-based)
  sig <- padj < p.adj.cutoff
  
  col <- rep("grey80", length(x))
  col[sig & x > logFC.cutoff] <- "forestgreen"
  col[sig & x < -logFC.cutoff] <- "red2"
  
  # derive p-value cutoff equivalent to significant set
  pvalCutoff <- max(p[sig], na.rm = TRUE)
  
  hline <- -log10(pvalCutoff)
  
  
  # labels
  score <- y * abs(x)
  top <- order(score, decreasing = TRUE)
  top <- head(top, n_label)
  
  par(mar = c(7, 4, 4, 2) + 0.1)
  
  # PLOT
  graphics::plot(
    x, y,
    col = col,
    pch = 16,
    cex = 0.5,
    main = contrastName,
    xlab = "Log2 Fold Change",
    ylab = "-Log10 p-value"
  )
  
  mtext(bquote(log[2]*FC~cutoff~.(logFC.cutoff)~and~FDR~cutoff~.(p.adj.cutoff)),
        side = 1,   # marge du bas
        line = 5,   # ajuste la distance sous l’axe X
        adj = 0.5,  # centre
        cex = 0.9)
  
  # lines
  graphics::abline(
    v = c(-logFC.cutoff, logFC.cutoff),
    lty = 2,
    col = "grey50"
  )
  
  graphics::abline(
    h = hline,
    lty = 2,
    col = "grey50"
  )
  
  # labels
  if (length(top) > 0) {
    graphics::text(
      x[top],
      y[top],
      labels = rn[top],
      pos = ifelse(x[top] > 0, 4, 2),
      cex = 0.5,
      offset = 0.3
    )
  }
  invisible(NULL)
}

#'.plotPValue
#'
#' @param data dataframe 
#' @param contrastName the contrast, useful for plot title
#' @return plot
#' @keywords internal
#' @noRd
.plotPValue_fast <- function(data, contrastName = contrastName) {
  PValue <- NULL
  
  graphics::hist(
    data$pvalue,
    breaks = 100,
    col = "grey70",
    border = "white",
    main = contrastName,
    xlab = "p-value"
  )
}