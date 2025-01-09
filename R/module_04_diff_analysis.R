### ============================================================================
### [04_diff_analysis] shiny modules
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif,


#' @importFrom UpSetR upset
#' @importFrom DT formatSignif datatable renderDataTable styleInterval
#' formatStyle styleEqual
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom purrr reduce
#' @importFrom magrittr "%>%"

# ---- main module ----

## ---- UI function ----
DiffExpAnalysisUI <- function(id){

    #name space for id
    ns <- NS(id)

    tagList(
        fluidRow(
            uiOutput(ns("instruction"))),

        ### parametres for Diff Analysis
        fluidRow(
            column(3,
                   uiOutput(ns("DiffParamUI"))),
            column(9,
                   uiOutput(ns("ResultsMerge")),
                   uiOutput(ns("validateUI"))
            ),
            column(12,uiOutput(ns("ContrastsResults")))
        )
    )
}

## ---- SERVER function ----
DiffExpAnalysis <- function(input, output, session, dataset, rea.values){

    local.rea.values <-
        reactiveValues(
            p.adj.cutoff      = 0.05,
            abs.logFC.cutoff  = 0,
            selectedContrasts = NULL,
            DiffExpAnal = NULL)

    # list of tools for diff analysis
    MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")

    method <- switch(rea.values[[dataset]]$omicsType,
                     "RNAseq"       = MethodList[1],
                     "proteomics"   = MethodList[2],
                     "metabolomics" = MethodList[2])

    ### ---- instruction ----
    output$instruction <- renderUI({
        box(
            title = span(
                tagList(
                    icon("cogs"),
                    paste0("Differential analysis, using ", names(method)),
                    # a(names(method), href="https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf"),
                    tags$small("(Scroll down for instructions)")
                )
            ),
            solidHeader = TRUE,
            status = "warning",
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            diff_analysis_docs
        )
    })

    ### ---- DiffParamUI ----
    output$DiffParamUI <- renderUI({

        #we must run process before
        validate(
            need(rea.values[[dataset]]$process != FALSE,
                 "Please run \'Data exploration and pre-processing\'")
        )
        validate(
            need(!is.null(rea.values$Contrasts.Sel),
                 "Please run \'Data exploration and pre-processing\'")
        )
        #design must be complete
        validate(
            need(rea.values[[dataset]]$compCheck != FALSE,
                 metadata(session$userData$FlomicsMultiAssay)$completeCheck[["error"]])
        )
        local.rea.values$selectedContrasts <-
            getSelectedContrasts(
                getProcessedData(session$userData$FlomicsMultiAssay[[dataset]],
                                 filter = TRUE)
            )

        validate(
            need(nrow(local.rea.values$selectedContrasts) != 0,
                 message = "No contrast matches the sample selection")
        )

        contrastList <- local.rea.values$selectedContrasts$contrastName
        names(contrastList) <-
            paste0("[",local.rea.values$selectedContrasts$tag, "] ",
                   local.rea.values$selectedContrasts$contrastName)

        ## getcontrast
        box(
            title  = span(tagList(icon("sliders-h"), "  ", "Setting")),
            width  = 14,
            status = "warning",

            fluidRow(
                column(
                    12,
                    ## list of contrasts to test
                    pickerInput(
                        inputId  = session$ns("contrastList"),
                        label    = .addBSpopify(label = 'Selected contrasts:',
                                                content = "Contrasts/hypotheses on which to run the differential analysis. If you want to test all contrasts, select 'All'"),
                        choices  = contrastList),

                    # method for Diff analysis
                    selectInput(
                        inputId  = session$ns("AnaDiffMethod"),
                        label = .addBSpopify(label = 'Method:',
                                             content = "Differential analysis method. Fixed parameter according to omics type."),
                        choices  = method,
                        selected = method),

                    numericInput(
                        inputId = session$ns("p.adj.cutoff"),
                        label=.addBSpopify(label = 'Adjusted pvalue cutoff:',
                                           content = "The adjusted p-value cut-off. Pvalues are adjusted using Benjamini-Hochberg method."),
                        value=local.rea.values$p.adj.cutoff, min=0, max=1, 0.01),
                    numericInput(
                        inputId = session$ns("abs.logFC.cutoff"),
                        label=.addBSpopify(label = '|log2FC| cutoff:',
                                           content = "The absolute log2 FC cut-off"),
                        value=local.rea.values$abs.logFC.cutoff, min=0, max=100, 0.01),

                    # use of cluster. need setting step
                    # materialSwitch(
                    #   inputId = session$ns("clustermq"),
                    #   label=.addBSpopify(label = 'use remote Cluster:',
                    #                      content = "send calculation to the cluster"),
                    #   value = FALSE, status = "success"),

                    actionButton(inputId = session$ns("runAnaDiff"),
                                 label = "Run", class = "butt")#,
                )
            ))

    })

    # filter param
    output$validateUI <- renderUI({

        if (rea.values[[dataset]]$diffAnal == FALSE) return()
        if (is.null(rea.values[[dataset]]$DiffValidContrast) ||
            dim(rea.values[[dataset]]$DiffValidContrast)[1] == 0) return()

        fluidRow(
            column(width = 9),
            column(
                width = 3,
                actionButton(session$ns("validContrast"),"Validate", class="butt")#,
            )
        )
    })

    ##================================ RUN =======================================##
    ### ---- runAnaDiff ----
    ### run diff
    # Run the differential analysis for each contrast set
    # Filter
    #   -> return a dynamic user interface with a collapsible box for each contrast
    #         - Pvalue graph
    #         - MAplot
    #         - Table of the DE genes
    #   -> combine data : union or intersection
    observeEvent(input$runAnaDiff, {

        # list of chosen parameters
        param.list <- list(method        = input$AnaDiffMethod,
                           # clustermq   = input$clustermq,
                           p.adj.method  = "BH",
                           p.adj.cutoff  = input$p.adj.cutoff,
                           abs.logFC.cutoff = input$abs.logFC.cutoff)

        # Prevent multiple executions
        if(check_run_diff_execution(session$userData$FlomicsMultiAssay[[dataset]],
                                    param.list) == FALSE) return()

        # Initialization of reactive variables
        rea.values[[dataset]]$diffAnal   <- FALSE
        rea.values[[dataset]]$diffValid  <- FALSE
        rea.values[[dataset]]$coExpAnal  <- FALSE
        rea.values[[dataset]]$DiffExp  <- FALSE
        rea.values[[dataset]]$CoExp <- FALSE
        rea.values[[dataset]]$DiffValidContrast <- NULL

        #---- progress bar ----#
        progress <- shiny::Progress$new()
        progress$set(message = "Run Diff", value = 0)
        on.exit(progress$close())
        progress$inc(1/10, detail = "in progress...")
        #----------------------#

        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]

        # Run the analysis only if the 'diff' object is empty
        if (length(metadata(dataset.SE)$DiffExpAnal) == 0) {

            message("[RFLOMICS] # 04- Differential Analysis... ", dataset)
            # run diff analysis with selected method
            dataset.SE <-
                runDiffAnalysis(
                    object        = dataset.SE,
                    p.adj.method  = "BH",
                    method        = input$AnaDiffMethod,
                    # clustermq   = input$clustermq,
                    p.adj.cutoff  = input$p.adj.cutoff,
                    logFC.cutoff  = input$abs.logFC.cutoff,
                    cmd           = TRUE)

        }else{
            # If the differential analysis has already run, do not run it again
            message("[RFLOMICS] # 04 => Filtering differential analysis... ", dataset)
            ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering

            dataset.SE <-
                filterDiffAnalysis(
                    object        = dataset.SE,
                    p.adj.cutoff  = input$p.adj.cutoff,
                    logFC.cutoff  = input$abs.logFC.cutoff)
        }

        session$userData$FlomicsMultiAssay[[dataset]] <- dataset.SE

        rea.values[[dataset]]$diffAnal <- TRUE

        rea.values$datasetDiff <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "DiffExpAnal")
        rea.values$datasetDiffAnnot <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "DiffExpEnrichAnal")
        rea.values$datasetCoEx <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "CoExpAnal")
        rea.values$datasetCoExAnnot <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "CoExpEnrichAnal")

        #---- progress bar ----#
        progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
        #----------------------#

    }, ignoreInit = TRUE)

    ### ---- validContrast ----
    ### validate contrasts
    observeEvent(input$validContrast, {

        rea.values[[dataset]]$diffValid  <- FALSE
        rea.values[[dataset]]$coExpAnal  <- FALSE
        rea.values[[dataset]]$DiffExp  <- FALSE
        rea.values[[dataset]]$CoExp <- FALSE

        session$userData$FlomicsMultiAssay <-
            resetRflomicsMAE(session$userData$FlomicsMultiAssay,
                             datasetNames = dataset,
                             singleAnalyses = c("DiffExpEnrichAnal",
                                                "CoExpAnal",
                                                "CoExpEnrichAnal"),
                             multiAnalyses = c("IntegrationAnalysis"))

        session$userData$FlomicsMultiAssay[[dataset]] <-
            setValidContrasts(session$userData$FlomicsMultiAssay[[dataset]],
                              contrastList = rea.values[[dataset]]$DiffValidContrast)

        rea.values[[dataset]]$diffValid <- TRUE

        # reset reactive values
        rea.values$datasetDiff <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "DiffExpAnal")
        rea.values$datasetDiffAnnot <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "DiffExpEnrichAnal")
        rea.values$datasetCoEx <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "CoExpAnal")
        rea.values$datasetCoExAnnot <-
            getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                    analyses = "CoExpEnrichAnal")

    }, ignoreInit = TRUE)

    ##============================== DISPLAY =====================================##

    ### ---- display results ----
    # display results per contrast
    output$ContrastsResults <- renderUI({

        dataset.SE  <- session$userData$FlomicsMultiAssay[[dataset]]
        diffExpAnal <- getAnalysis(dataset.SE, name = "DiffExpAnal")

        if (rea.values[[dataset]]$diffAnal == FALSE ||
            is.null(diffExpAnal[["results"]][["DEF"]])) return()

        list(
            lapply(seq_len(nrow(local.rea.values$selectedContrasts)), function(i) {

                vect     <- unlist(local.rea.values$selectedContrasts[i,])
                stats    <- getDiffStat(dataset.SE)[vect["contrastName"],]
                DEList   <- getDEList(dataset.SE, contrasts = vect["contrastName"])

                diff.plots <-
                    plotDiffAnalysis(dataset.SE, contrastName = vect["contrastName"])
                stat <- paste0("[#DE: ", stats[["All"]], "]")

                # panel list
                # default panels (no DE)
                tabPanel.list <- list(

                    #### ---- pvalue plot ----
                    tabPanel(
                        title = "Pvalue's distribution",
                        tags$br(),
                        tags$i("You must have a look at the distribution of non-adjusted
            p-values to validate your analysis.
            The most desirable shape is a peak of p-values at 0 followed
            by a uniform (flat) distribution. If there is a peak in 1, consider
                   increasing the filtering threshold in the pre-processing
                   step."),
                        tags$br(),tags$hr(),tags$br(),
                        renderPlot({ diff.plots$Pvalue.hist })
                    ),

                    #### ---- MA plot ----
                    tabPanel(
                        title = "MA plot",
                        tags$br(),
                        tags$i(
                            paste0("It is expected that a majority of ",
                                   .omicsDic(dataset.SE)$variableName,
                                   " gather around 0. The red dots are the ",
                                   .omicsDic(dataset.SE)$variableName,
                                   " significantly over-expressed in the left factor's
                          level(s) in the contrast expression whereas ",
                                   "blue dots are ",
                                   .omicsDic(dataset.SE)$variableName,
                                   " significantly under-expressed in the rigth factor's
                          level(s) in the contrast expression.
                              Only the top 20 ",
                                   .omicsDic(dataset.SE)$variableName,
                                   " DE are labeled.")),
                        tags$br(),tags$hr(),tags$br(),
                        renderPlot({ diff.plots$MA.plot })
                    ),

                    #### ---- Volcano plot ----
                    tabPanel(
                        title = "Volcano plot",
                        tags$br(),
                        tags$i(
                            paste0("Red dots are ",
                                   .omicsDic(dataset.SE)$variableName," of interest: ",
                                   "displaying both large magnitude log2-fold-changes
                          (x axis) and high statistical significance (y axis)",
                                   "Only the top 20 ",
                                   .omicsDic(dataset.SE)$variableName,
                                   " DE are labeled.")),
                        tags$br(), tags$hr(), tags$br(),
                        renderPlot({ diff.plots$Volcano.plot }, height = 600))
                )

                # if DE
                ## PCA only if nbr of variables > 3
                if (length(DEList) > 3)
                    tabPanel.list <-
                    c(tabPanel.list,
                      list(
                          #### ---- PCA plot ----
                          tabPanel(
                              title = paste0("PCA on DE ", .omicsDic(dataset.SE)$variableName),
                              tags$br(),
                              .modVariablePCAUI(
                                  session$ns(paste0(vect["contrastName"],"-DE")))
                          )
                      )

                    )
                ## headmap, boxplot and table of DE only if DE nb > 0
                if (length(DEList) > 0){

                    tabPanel.list <-
                        c(tabPanel.list,
                          list(
                              #### ---- DEF table ----
                              tabPanel(
                                  title = "Table",
                                  tags$br(),
                                  .modDEtableUI(
                                      session$ns(paste0(vect["contrastName"],"-DE")))
                              ),

                              #### ---- Heatmap plot ----
                              tabPanel(
                                  title = "Heatmap",
                                  tags$br(),
                                  .modDEheatmapUI(
                                      session$ns(paste0(vect["contrastName"],"-DE")))
                              ),

                              #### ---- DE boxplot plot ----
                              tabPanel(
                                  title = "boxplots",
                                  tags$br(),
                                  .modVariableBoxplotUI(
                                      session$ns(paste0(vect["contrastName"],"-DE")))
                              )
                          )
                        )

                    stat <- paste0(
                        "[#DE: ", stats[["All"]], " ; ",
                        "Up: ",stats[["Up"]],
                        " (",round(stats[["Up"]]/stats[["All"]],2)*100," %)",
                        " ; ", "Down: ", stats[["Down"]],
                        " (",round(stats[["Down"]]/stats[["All"]],2)*100,"%)]"
                    )
                }

                #### ---- display panels ----
                # if error message specific to current contrast
                if(!is.null(
                    diffExpAnal[["results"]][["runErrors"]][[vect["contrastName"]]])){
                    fluidRow(
                        column(
                            width = 10,
                            box(
                                width=14,
                                status = "danger",
                                title =
                                    diffExpAnal[["results"]][["runErrors"]][[vect["contrastName"]]]
                            )
                        )
                    )
                }else{
                    # if no error
                    fluidRow(
                        column(
                            width = 10,
                            box(
                                width=14,
                                solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                status = ifelse(length(DEList) != 0, "success", "danger"),
                                title = tags$h5(
                                    paste0("[",vect["tag"], "] ", vect["contrastName"], " ", stat)),

                                do.call(what = tabsetPanel, args = tabPanel.list)
                            )
                        ),
                        column(
                            width = 2,
                            if (length(DEList) != 0){
                                checkboxInput(
                                    session$ns(paste0("checkbox_", vect[["contrastName"]])), "OK", value = TRUE)
                            }
                        )
                    )
                }
            })
        )
    })
    observe({

        dataset.SE  <- session$userData$FlomicsMultiAssay[[dataset]]
        diffExpAnal <- getAnalysis(dataset.SE, name = "DiffExpAnal")

        if(rea.values[[dataset]]$diffAnal == FALSE) return()
        if(is.null(diffExpAnal[["results"]][["TopDEF"]])) return()

        lapply(seq_len(length(rea.values$Contrasts.Sel$contrast)), function(i) {

            vect     <- unlist(rea.values$Contrasts.Sel[i,])

            # PCA axis for plot
            # update/adapt PCA axis
            callModule(UpdateRadioButtons, paste0(vect["contrastName"],"-diff"))

            DEList <- getDEList(dataset.SE, contrasts = vect["contrastName"])
            callModule(module = .modVariablePCA,
                       id = paste0(vect["contrastName"],"-DE"),
                       dataset.SE = dataset.SE,
                       id2 = vect["contrastName"],
                       variableList = DEList)

            callModule(module = .modDEtable,
                       id = paste0(vect["contrastName"],"-DE"),
                       dataset.SE = dataset.SE,
                       contrastName = vect["contrastName"])

            callModule(module = .modDEheatmap,
                       id = paste0(vect["contrastName"],"-DE"),
                       dataset.SE = dataset.SE,
                       contrastName = vect["contrastName"])

            callModule(module = .modVariableBoxplot,
                       id = paste0(vect["contrastName"],"-DE"),
                       dataset.SE = dataset.SE,
                       contrastName = vect["contrastName"],
                       variableList = DEList)

            # # update SelectizeInput for boxplot DE
            # DEList  <- rownames(diffExpAnal[["TopDEF"]][[vect["contrastName"]]])
            #
            # updateSelectizeInput(
            #   session = session,
            #   inputId = paste0(vect["contrastName"],"-DE"),
            #   choices = DEList,
            #   server = TRUE
            # )
        })
    })

    ### ---- summary ----
    # merge results on upset plot
    output$ResultsMerge <- renderUI({

        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
        diffExpAnal <- getAnalysis(dataset.SE, name = "DiffExpAnal")

        if (rea.values[[dataset]]$diffAnal == FALSE) return()

        if (!is.null(diffExpAnal[["errors"]])) {

            tagList(
                box(
                    width=14,
                    status = "danger",
                    solidHeader = TRUE,
                    title = diffExpAnal[["errors"]]
                )
            )

        } else if (!is.null(diffExpAnal[["results"]][["mergeDEF"]])){

            ### ---- barplot ----
            # 1rst panel
            tabPanel.list <- list(
                tabPanel(
                    title = "Summary",
                    renderPlot({
                        # get rflomicsMAE object with only 1
                        subMAE <-
                            subRflomicsMAE(session$userData$FlomicsMultiAssay, dataset)

                        contrasts <- getSelectedContrasts(subMAE[[dataset]])

                        subMAE <-
                            setValidContrasts(
                                object       = subMAE,
                                omicName     = dataset,
                                contrastList = getSelectedContrasts(subMAE[[dataset]]))

                        getDiffAnalysesSummary(subMAE, plot = TRUE, interface = TRUE)  +
                            theme(legend.text = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  strip.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = max(10, 20-nrow(contrasts))))
                    })
                )
            )

            ### ---- upset ----
            # 2nd panel (upset) if 2 validated contrasts or more
            DEF_mat <- diffExpAnal[["results"]][["mergeDEF"]]

            index <- sapply(names(DEF_mat)[-1], function(x){(
                input[[paste0("checkbox_",x)]])
            }) |> unlist()

            H_selected   <- names(index)[index]
            DEF_selected <- select(DEF_mat,any_of(H_selected))

            rea.values[[dataset]]$DiffValidContrast <-
                filter(local.rea.values$selectedContrasts, contrastName %in% H_selected)

            colnames(DEF_selected) <-
                rea.values[[dataset]]$DiffValidContrast$tag

            if (length(H_selected) > 1){

                tabPanel.list <- c(list(
                    tabPanel(
                        title = "Intersection between DE lists",
                        renderPlot({
                            upset(data = DEF_selected,
                                  sets = colnames(DEF_selected),
                                  order.by = "freq")})
                    )
                ),
                tabPanel.list)
            }

            tagList(
                box(
                    width=14,
                    status = "warning",
                    solidHeader = FALSE,

                    do.call(what = tabsetPanel, args = tabPanel.list)
                )
            )
        }
    })

    return(input)
}


# ---- sub modules ----

## ---- DE boxplot ----
.modVariableBoxplotUI <- function(id){

    #name space for id
    ns <- NS(id)
    uiOutput(ns("boxplotDEUI"))

}

.modVariableBoxplot <- function(input, output, session, dataset.SE, contrastName, variableList){

    output$boxplotDEUI <- renderUI(

        tagList(
            tags$i(paste0("Boxplot showing the expression/abundance profile of a selected DE ",
                          .omicsDic(dataset.SE)$variableName),
                   " colored by experimental factor's levels (see Levels radio buttons)."),

            tags$br(), tags$hr(), tags$br(),

            fluidRow(
                column(
                    width = 3,
                    # selectizeInput(
                    #   inputId = session$ns(paste0(contrastName,"-DE")),
                    #   label = paste0("Select ",.omicsDic(dataset.SE)$variableName,":"),
                    #   multiple = FALSE,
                    #   choices = NULL,
                    #   options = list(maxOptions = 1000)
                    # ),
                    selectizeInput(
                        inputId = session$ns(paste0(contrastName,"-DE")),
                        label = paste0(.omicsDic(dataset.SE)$variableName,":"),
                        multiple = FALSE,
                        choices = variableList,
                        options = list(maxOptions = 1000)
                    ),
                    radioButtons(
                        inputId = session$ns(paste0(contrastName,"-DEcondition")),
                        label = 'Levels:',
                        choices = c("groups",getBioFactors(dataset.SE)),
                        selected = getBioFactors(dataset.SE)[1]
                    )
                ),
                column(
                    width = 9,
                    renderPlot({

                        plotBoxplotDE(
                            object=dataset.SE,
                            featureName=input[[paste0(contrastName, "-DE")]],
                            groupColor=input[[paste0(contrastName, "-DEcondition")]])
                    })
                )
            )
        )
    )
}

## ---- DE heatmap ----
.modDEheatmapUI <- function(id){

    #name space for id
    ns <- NS(id)
    uiOutput(ns("heatmapDEUI"))

}

.modDEheatmap <- function(input, output, session, dataset.SE, contrastName){

    output$heatmapDEUI <- renderUI(

        tagList(
            tags$i(tags$p(
                paste0("Heatmap is performed on DE ", .omicsDic(dataset.SE)$variableName,
                       " expression data table which has been transformed by: ",
                       getTransSettings(dataset.SE)$method,
                       " method and normalized by: ",
                       getNormSettings(dataset.SE)$method , "  method."))),
            tags$i(tags$p(
                paste0("Clustering is independently performed on samples (row) and centered ",
                       .omicsDic(dataset.SE)$variableName,
                       " (column) using euclidian distance and complete aggregation method."))),
            tags$i(tags$p(
                " You can split the heatmap by the levels of a factor of interest
        (Levels radio buttons). You may also add or remove annotations to the heatmap by
        selecting biological/batch factors to display.")),
            tags$hr(),
            renderPlot({
                annot_arg <-
                    c(input[[paste0(contrastName, "-annotBio")]],
                      input[[paste0(contrastName, "-annotBatch")]])
                if (length(getMetaFactors(dataset.SE)) > 0) {
                    annot_arg <-
                        c(annot_arg, input[[paste0(contrastName, "-annotMeta")]])
                }

                plotHeatmapDesign(
                    object       = dataset.SE,
                    contrastName = contrastName,
                    splitFactor  = input[[paste0(contrastName, "-heat.condColorSelect")]],
                    annotNames   =  annot_arg)
            }),
            ## select cluster to plot
            column(
                6, radioButtons(
                    inputId  = session$ns(paste0(contrastName, "-heat.condColorSelect")),
                    label    = 'Levels:',
                    choices  = c("none", getBioFactors(dataset.SE)),
                    selected = "none", inline = TRUE)),

            ## select annotations to show
            column(
                width = 6 ,
                checkboxGroupInput(
                    inputId = session$ns(paste0(contrastName, "-annotBio")),
                    label = "Biological factors", inline = TRUE,
                    choices = getBioFactors(dataset.SE),
                    selected = getBioFactors(dataset.SE))),
            column(
                width = 6,
                checkboxGroupInput(
                    inputId = session$ns(paste0(contrastName, "-annotBatch")),
                    label = "Batch factors",  inline = TRUE,
                    choices = getBatchFactors(dataset.SE),
                    selected = NULL)),
            if (length(getMetaFactors(dataset.SE)) > 0) {
                column(
                    width = 6,
                    checkboxGroupInput(
                        inputId = session$ns(paste0(contrastName, "-annotMeta")),
                        label = "Metadata factors", inline = TRUE,
                        choices = getMetaFactors(dataset.SE),
                        selected = NULL))
            }
        )

    )
}


## ---- diff table ----
.modDEtableUI <- function(id){

    #name space for id
    ns <- NS(id)
    uiOutput(ns("tableDEUI"))

}

.modDEtable <- function(input, output, session, dataset.SE, contrastName){

    output$tableDEUI <- renderUI(

        tagList(
            tags$i("Table of results of the differential expression/abundance statistical analysis:"),
            tags$ul(
                tags$li(tags$i(paste0("Row names: ",.omicsDic(dataset.SE)$variableName," ID"))),
                tags$li(tags$i("logFC: log2 fold change")),
                tags$li(tags$i("Abundance: mean expression/abundance for the factor's levels")),
                tags$li(tags$i("t: t-statistic (limma-lmFit, proteomics/metabolomics only)")),
                tags$li(tags$i("pvalue: p-values")),
                tags$li(tags$i("Adj.value: adjusted p-value (BH)")),
                tags$li(tags$i("LR: likelihood ratio test (edgeR-glmLRT, transcriptomics RNAseq only)")),
                tags$li(tags$i("B: log-odds that the proteomics/metabolomics is differentially expressed (limma-topTable)")),
                tags$li(tags$i("Regulation: Up (green) or Down (red) regulated"))
            ),
            tags$hr(), tags$br(),
            ### DEF result table ###
            DT::renderDataTable({
                resTable <-
                    metadata(dataset.SE)$DiffExpAnal[["results"]][["TopDEF"]][[contrastName]]
                resTable$Regulation <- ifelse(resTable$logFC > 0, "Up", "Down")
                resTable %>% DT::datatable(
                    extensions = 'Buttons',
                    options = list(dom = 'lfrtipB',
                                   rownames = FALSE,
                                   pageLength = 10,
                                   buttons = c('csv', 'excel'),
                                   lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All")))) %>%
                    formatStyle(
                        'Regulation',
                        backgroundColor = DT::styleEqual(
                            c("Up", "Down"), c("#C7DCA7", c("#FFC5C5"))),
                        fontWeight = 'bold') %>%
                    formatSignif(columns = seq_len(dim(resTable)[2]), digits = 3)
            }, server = FALSE)
        )
    )
}


## ---- variable PCA ----
.modVariablePCAUI <- function(id){

    #name space for id
    ns <- NS(id)
    uiOutput(ns("pcaUI"))

}

.modVariablePCA <- function(input, output, session, dataset.SE, id2, variableList){

    output$pcaUI <- renderUI(

        tagList(
            tags$i(paste0("PCA plot of the differentially expressed  ",
                          .omicsDic(dataset.SE)$variableName,
                          ". It is expected that the separation between groups of
                       interest is better after the differential analysis.")),
            tags$br(), tags$hr(), tags$br(),
            fluidRow(
                column(
                    width = 12,
                    renderPlot({
                        # in case of data exploratory variableList: all features
                        # in case on diff analysis variableList: DE list
                        newDataset.SE <- runOmicsPCA(dataset.SE[variableList], ncomp = 5, raw = FALSE)

                        PC1.value <- as.numeric(
                            input[[paste0(id2,"-diff-Firstaxis")]][1])
                        PC2.value <- as.numeric(
                            input[[paste0(id2,"-diff-Secondaxis")]][1])
                        condGroup <-
                            input[[paste0(id2,"-pca.DE.condColorSelect")]][1]

                        plotOmicsPCA(newDataset.SE,
                                     raw = FALSE,
                                     axes = c(PC1.value, PC2.value),
                                     groupColor = condGroup)
                    })
                )
            ),
            fluidRow(
                column(
                    width = 6,
                    radioButtons(
                        inputId = session$ns(paste0(id2, "-pca.DE.condColorSelect")),
                        label = 'Levels:',
                        choices = c("groups",getBioFactors(dataset.SE)),
                        selected = "groups")),
                column(
                    width = 6,
                    UpdateRadioButtonsUI(session$ns(paste0(id2, "-diff"))))
            )
        )
    )
}


# ---- functions ----

## ----- check run diff execution ------
check_run_diff_execution <- function(object.SE, param.list = NULL){

    # filtering setting
    if (length(metadata(object.SE)[["DiffExpAnal"]]) == 0 ||
        !is.null(metadata(object.SE)[["DiffExpAnal"]])) return(TRUE)

    if(param.list$method != getDiffSettings(object.SE)$method) return(TRUE)
    if(param.list$p.adj.method  != getDiffSettings(object.SE)$p.adj.method) return(TRUE)
    if(param.list$p.adj.cutoff  != getDiffSettings(object.SE)$p.adj.cutoff) return(TRUE)
    if(param.list$abs.logFC.cutoff   != getDiffSettings(object.SE)$abs.logFC.cutoff)  return(TRUE)

    return(FALSE)
}
