### ============================================================================
### [00_common] commun modules
### ----------------------------------------------------------------------------
# N. Bessoltane,


# ---- update radio button ----
#' @keywords internal
#' @noRd
UpdateRadioButtonsUI <- function(id) {
    # Namespace for id
    ns <- NS(id)

    tagList(
        # Radio buttons for first axis
        radioButtons(
            inputId = ns("Firstaxis"),
            label = "Choice of PCs:",
            choices = list(
                "PC1" = 1,
                "PC2" = 2
            ),
            selected = 1,
            inline = TRUE
        ),

        # Radio buttons for second axis
        radioButtons(
            inputId = ns("Secondaxis"),
            label = "",
            choices = list(
                "PC2" = 2,
                "PC3" = 3
            ),
            selected = 2,
            inline = TRUE
        )
    )
}

#' @keywords internal
#' @noRd
UpdateRadioButtons <- function(input, output, session) {
    # Update Secondaxis when Firstaxis changes
    observeEvent(input$Firstaxis, {
        selected_first <- input$Firstaxis
        choices <- c("PC1" = 1, "PC2" = 2, "PC3" = 3)
        updateRadioButtons(session,
                           "Secondaxis",
                           choices = choices[-as.numeric(selected_first)],
                           selected = if (input$Secondaxis == selected_first) {
                               choices[-as.numeric(selected_first)][1]
                           } else {
                               input$Secondaxis
                           },
                           inline = TRUE)
    })

    # Update Firstaxis when Secondaxis changes
    observeEvent(input$Secondaxis, {
        selected_second <- input$Secondaxis
        choices <- c("PC1" = 1, "PC2" = 2, "PC3" = 3)
        updateRadioButtons(session,
                           "Firstaxis",
                           choices = choices[-as.numeric(selected_second)],
                           selected = if (input$Firstaxis == selected_second) {
                               choices[-as.numeric(selected_second)][1]
                           } else {
                               input$Firstaxis
                           },
                           inline = TRUE)
    })
}


#' @keywords internal
#' @noRd
RadioButtonsConditionUI <- function(id) {
    #name space for id
    ns <- NS(id)

    tagList(uiOutput(ns('condColor')),)
}

#' @keywords internal
#' @noRd
RadioButtonsCondition <- function(input, output, session, typeFact) {
    # select factors for color PCA plot
    output$condColor <- renderUI({
        factors <- getFactorTypes(session$userData$FlomicsMultiAssay)
        factors <- factors[factors %in% typeFact]
        condition <- names(factors)

        if (!any(typeFact %in% "Meta"))
            condition <- c("groups", condition)

        radioButtons(inputId = session$ns("condColorSelect"),
                     label = 'Levels:',
                     choices = condition,
                     selected = condition[1]
        )
    })
}

# ---- summary of all analysed data ----
#' @keywords internal
#' @noRd
.modSingleOmicAnalysesSummaryUI <- function(id) {
    ns <- NS(id)

    tagList(fluidPage(column(
        width = 12,
        fluidRow(uiOutput(ns("overView"))),
        fluidRow(uiOutput(ns("DiffSummary"))),
        fluidRow(uiOutput(ns("CoExSummary")))
    )))
}

#' @keywords internal
#' @noRd
.modSingleOmicAnalysesSummary <-
    function(input, output, session, rea.values) {

        # over view of dataset dimensions after processing
        output$overView <- renderUI({
            if (is.null(rea.values$datasetProcess))
                return()
            toto <<- session$userData$FlomicsMultiAssay
            box(title = "Dataset overview after data processing",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = FALSE,

                renderPlot({
                    plotDataOverview(
                        session$userData$FlomicsMultiAssay,
                        omicNames = rea.values$datasetProcess,
                        raw = FALSE
                    )
                })
            )
        })

        # summary of diff analysis on all dataset
        output$DiffSummary <- renderUI({

            res <- getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay)

            if (is.null(rea.values$datasetDiff))
                return()

            tabPanel.list <-
                list(
                    tabPanel(
                        title = "DE results",
                        renderPlot({
                            getDiffAnalysesSummary(
                                session$userData$FlomicsMultiAssay, plot = TRUE)
                        })
                    )
                )

            if (!is.null(rea.values$datasetDiffAnnot)) {

                H_tag  <- getSelectedContrasts(session$userData$FlomicsMultiAssay)

                tabPanel_annot.list <-
                    lapply(names(rea.values$datasetDiffAnnot), function(database) {

                        ListNames  <- vector()
                        domainList <- vector()
                        termNbr    <- 0
                        for(dataset in rea.values$datasetDiffAnnot[[database]]){

                            tmp <-
                                getAnalysis(session$userData$FlomicsMultiAssay[[dataset]],
                                            name = "DiffExpEnrichAnal",
                                            subName = database)

                            tmp <- tmp$results$summary %>%
                                mutate(sum = rowSums(across(where(is.numeric)))) %>%
                                filter(sum != 0)

                            ListNames <- unique(c(ListNames, row.names(tmp)))

                            tmp <- tmp[,-1]
                            tmp$sum <- NULL

                            domainList <- unique(c(domainList, names(colSums(tmp)[colSums(tmp) != 0])))
                            termNbr    <- termNbr + sum(colSums(tmp))
                        }
                        names(ListNames) <-
                            paste0("[",H_tag[H_tag$contrastName %in% ListNames,]$tag,"] ",
                                   ListNames)

                        tabPanel(
                            title = paste0("ORA results from ", database),
                            fluidRow(
                                #if (!identical(domainList, "no-domain")) {
                                column(
                                    width = 4,
                                    radioButtons(
                                        inputId  = session$ns(paste0(database, "-domain.diff")),
                                        label    = "Domain",
                                        choices  = domainList,
                                        selected = domainList[1],
                                        inline   = TRUE
                                    )
                                )
                                #}
                                ,column(
                                    width = 4,
                                    pickerInput(
                                        inputId  = session$ns(paste0(database, "-datasets.diff")),
                                        label    = "Dataset list:",
                                        choices  = rea.values$datasetDiffAnnot[[database]],
                                        selected = rea.values$datasetDiffAnnot[[database]],
                                        multiple = TRUE,
                                        options  = list(`actions-box` = TRUE)
                                    )
                                ),
                                column(
                                    width = 4,
                                    pickerInput(
                                        inputId  = session$ns(paste0(database, "-contrasts.diff")),
                                        label    = "Contrast list:",
                                        choices  = ListNames,
                                        selected = ListNames,
                                        multiple = TRUE,
                                        options  = list(`actions-box` = TRUE)
                                    )
                                )
                            ),
                            fluidRow(
                                column(
                                    width = 12,
                                    renderPlot({
                                        p.list <- getAnnotAnalysesSummary(
                                            session$userData$FlomicsMultiAssay,
                                            from = "DiffExp",
                                            databases = database,
                                            listNames = input[[paste0(database, "-contrasts.diff")]],
                                            omicNames = input[[paste0(database, "-datasets.diff")]]
                                        )

                                        p.list[[database]][[input[[paste0(database, "-domain.diff")]]]]
                                    },
                                    height = function() {min(1200, max(200, termNbr * 20))})
                                )
                            )
                        )
                    })

                tabPanel.list <- c(tabPanel.list, tabPanel_annot.list)
            }

            box(
                title = "Summary of differential expression analyses",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,
                tagList({

                    do.call(what = tabsetPanel, args = tabPanel.list)
                })
            )
        })

        # coexpression summary
        output$CoExSummary <- renderUI({

            if (is.null(rea.values$datasetCoEx))
                return()

            tabPanel.list <-
                list(
                    tabPanel(
                        title = "CoExp results",
                        renderPlot({
                            getCoExpAnalysesSummary(
                                session$userData$FlomicsMultiAssay)
                        })
                    )
                )

            if (!is.null(rea.values$datasetCoExAnnot)) {

                tabPanel_annot.list <-
                    lapply(names(rea.values$datasetCoExAnnot), function(database) {

                        # list of cluster
                        ListNames  <- vector()
                        domainList <- vector()
                        termNbr    <- 0
                        for(dataset in rea.values$datasetDiffAnnot[[database]]){

                            tmp <-
                                getAnalysis(session$userData$FlomicsMultiAssay[[dataset]],
                                            name = "CoExpEnrichAnal",
                                            subName = database)

                            tmp <- tmp$results$summary %>%
                                mutate(sum = rowSums(across(where(is.numeric)))) %>%
                                filter(sum != 0)

                            ListNames <- unique(c(ListNames, row.names(tmp)))

                            tmp <- tmp[,-1]
                            tmp$sum <- NULL

                            domainList <- unique(c(domainList, names(colSums(tmp)[colSums(tmp) != 0])))
                            termNbr    <- termNbr + sum(colSums(tmp))
                        }

                        tabPanel(
                            title = paste0("ORA results from ", database),
                            fluidRow(
                                column(
                                    width = 4,
                                    pickerInput(
                                        inputId  = session$ns(paste0(database, "-datasets.coex")),
                                        label    = "Dataset list:",
                                        choices  = rea.values$datasetCoExAnnot[[database]],
                                        selected = rea.values$datasetCoExAnnot[[database]],
                                        multiple = TRUE,
                                        options  = list(`actions-box` = TRUE)
                                    )
                                ),
                                column(
                                    width = 4,
                                    pickerInput(
                                        inputId  = session$ns(paste0(database, "-clusters.coex")),
                                        label    = "Cluster list:",
                                        choices  = ListNames,
                                        selected = ListNames,
                                        multiple = TRUE,
                                        options  = list(`actions-box` = TRUE)
                                    )
                                ),
                                column(
                                    width = 4,
                                    radioButtons(
                                        inputId = session$ns(paste0(database, "-domain.coex")),
                                        label = "Domain",
                                        choices = domainList,
                                        selected = domainList[1],
                                        inline = TRUE
                                    )
                                )
                            ),
                            fluidRow(
                                column(
                                    width = 12,
                                    renderPlot({

                                        p.list <- getAnnotAnalysesSummary(
                                            session$userData$FlomicsMultiAssay,
                                            from = "CoExp",
                                            databases = database,
                                            listNames = input[[paste0(database, "-clusters.coex")]],
                                            omicNames = input[[paste0(database, "-datasets.coex")]]
                                        )

                                        p.list[[database]][[input[[paste0(database, "-domain.coex")]]]]
                                    },
                                    height = function() {min(1200, max(200, termNbr * 20))})
                                )
                            )

                        )
                    })
                tabPanel.list <- c(tabPanel.list,tabPanel_annot.list)
            }

            box(
                title = "Summary of Co-expression analyses",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,

                do.call(what = tabsetPanel, args = tabPanel.list)
            )
        })
    }


# ---- selectizeModuleServer ----
# Module UI
.selectizeModuleUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("select_ui"))
    )
}

# Module Server
.selectizeModuleServer <- function(id, featureType, choices) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Dynamically generate the selectizeInput in the server.
        # Initialize choices to NULL to avoid loading everything on the client side.
        # choices = NULL
        output$select_ui <- renderUI({
            selectizeInput(
                inputId = ns("selectFeature"),
                label = paste0("Select DE ",featureType,":"),
                multiple = FALSE,
                choices = NULL,
                options = list(maxOptions = 1000)
            )
        })

        # Use server-side mode to handle large lists.
        updateSelectizeInput(
            session = session,
            inputId = "selectFeature",
            choices = choices,
            server = TRUE
        )

        # Return the reactive selection directly.
        return(reactive(input$selectFeature))
    })
}



# ---- Module load State ----

#' @keywords internal
#' @noRd
.modLoadStateUI <- function(id) {
    ns <- NS(id)

    textExp <- "On this page, you can restore previously saved (or bookmarked, in shiny language) states.
    When bookmarking, the app will automatically save the inputs and tables inside a folder named
    shiny_bookmarks/ (or tmp/shiny_bookmards/ depending on your operating system). If you have any state
    in this bookmark folder, they will appear here."

    tagList(
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "classesStyles.css")
        ),
        br(),
        box(
            title = "Loading and restoring previous bookmarked states",
            width = 12,
            status = "warning",
            tags$style(
                ".explain-p {
                            color: Gray;
                            text-justify: inter-word;
                            font-style: italic;
                            }"
            ),
            div(class = "explain-p", HTML(textExp)),
            br(),
            br(),
            uiOutput(ns("tableState"))
        )
    )
}

#' @keywords internal
#' @noRd
.modLoadState <-  function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        folders <- reactive({
            folder <- getwd()
            if (dir.exists(paste0(folder, "/shiny_bookmarks/"))) {
                folder <- paste0(folder, "/shiny_bookmarks/")
            } else if (dir.exists(paste0(folder, "/tmp/shiny_bookmarks/"))) {
                folder <- paste0(folder, "/tmp/shiny_bookmarks/")
            } else {
                return(NULL)
            }

            if (!is.null(folder)) {
                dirs <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
                if (length(dirs) == 0) return(NULL)

                df <- data.frame(
                    "State_Folder" = basename(dirs),
                    "Last_Modification" = file.info(dirs)$mtime,
                    stringsAsFactors = FALSE
                )
                df[order(df[["Last_Modification"]], decreasing = TRUE),]
            }


        })

        # Table of state folders available, with a button for each
        output$tableState <- renderUI({
            df <- folders()
            if (is.null(df)) {
                return(
                    renderText(expr = {"There is no shiny_bookmarks/ or tmp/shiny_bookmarks/ folder in your working directory"})
                )
            }

            df[["Choose"]] <- vapply(
                df[["State_Folder"]],
                function(x) {
                    as.character(actionButton(
                        ns(paste0("choose_", x)), "Choose",
                        onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})",
                                          ns("chosen"), x)
                    ))
                },
                character(1)
            )

            datatable(
                df,
                escape = FALSE,
                options = list(dom = "t", paging = FALSE),
                colnames = c("Stage Folder", "Last Modification", "Choose")
            )
        })

        # Observer for updating the url
        observeEvent(input$chosen, {
            chosenState <- input$chosen
            updateQueryString(
                paste0("?_state_id_=", chosenState),
                mode = "replace", session = session
            )
            session$reload()
        })
    })
}

# ---- Restore observer (used in server.R) ----
#' @keywords internal
#' @noRd
.updateValue <- function(
        session,
        rea.values,
        tabName,
        valueToObserve,
        nextValueToObserve,
        nextPatternToObserve) {

    rea.values$restoreValues[[nextValueToObserve]] <- .nextValues(nextPatternToObserve, rea.values$stateInput)
    rea.values$restoreValues$counter <- 1
    namAnalysis <- .namAnalysis(rea.values$restoreValues[[valueToObserve]], rea.values$restoreValues$counter)
    updateTabItems(session, inputId = "tabs", selected = namAnalysis)
    updateTabsetPanel(session, inputId = rea.values$restoreValues[[nextValueToObserve]][rea.values$restoreValues$counter], selected = tabName)
    shinyjs::js$showmenuItem(namAnalysis)

}

#' @keywords internal
#' @noRd
.processValue <- function(
        session,
        rea.values,
        valueToObserve,
        tabName
) {
    namAnalysis <- .namAnalysis(rea.values$restoreValues[[valueToObserve]], rea.values$restoreValues$counter)
    updateTabItems(session, inputId = "tabs", selected = namAnalysis)
    updateTabsetPanel(session, inputId = rea.values$restoreValues[[valueToObserve]][rea.values$restoreValues$counter], selected = tabName)
    shinyjs::js$showmenuItem("omicsanalysis")
    rea.values$restoreValues$counter <- rea.values$restoreValues$counter + 1
}

#' @keywords internal
#' @noRd
.initValue <- function(
        session,
        rea.values,
        valueToObserve,
        tabName
) {

    rea.values$restoreValues$counter <- 1
    namAnalysis <- .namAnalysis(rea.values$restoreValues[[valueToObserve]], rea.values$restoreValues$counter)
    # shinyjs::js$showmenuItem("omicsanalysis")
    updateTabItems(session, inputId = "tabs", selected = namAnalysis)
    updateTabsetPanel(session, inputId = rea.values$restoreValues[[valueToObserve]][rea.values$restoreValues$counter], selected = tabName)
    shinyjs::js$showmenuItem("omicsanalysis")
    rea.values$restoreValues$counter <- 2
}

#' @keywords internal
#' @noRd
.nextValues <- function(pattern, searchVector) {
    varinter <- searchVector[grep(pattern, names(searchVector))]
    return(gsub(pattern, "", names(varinter[which(varinter > 0)])))
}

#' @keywords internal
#' @noRd
.namAnalysis <- function(values, counter) {
    sub("([a-zA-Z]+)([0-9]+)", "\\1Analysis\\2", values[counter])
}

#' @keywords internal
#' @noRd
.observePreProcess <- function(session, rea.values){
    destroyPreProcessObserver <- FALSE
    observeEvent(rea.values$analysis, {# little difference here, can"t used the common observeRestore ?
        .initValue(session, rea.values, valueToObserve = "Processed", tabName = "Pre-processing")
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values$datasetProcess, {
        if (!rea.values$restoring) {
            destroyPreProcessObserver <- TRUE
        } else if (length(rea.values$restoreValues$Processed) >= rea.values$restoreValues$counter) {
            .processValue(session, rea.values, "Processed", tabName = "Pre-processing")
        } else if (length(rea.values$restoreValues$Processed) < rea.values$restoreValues$counter) {
            .updateValue(session, rea.values, tabName = "Pre-processing", "Processed", "Diff", "-validContrast$")
            destroyPreProcessObserver <- TRUE
        }
    }, ignoreInit = TRUE, once = destroyPreProcessObserver)

}


#' @keywords internal
#' @noRd
.observeRestore <- function(session, rea.values,
                            reavalToObserve,
                            valueToObserve, tabName,
                            nextValueToObserve,
                            nextPatternToObserve) {
    destroyThisObserver <- FALSE

    observeEvent(rea.values$restoreValues[[valueToObserve]], {
        .initValue(session, rea.values, valueToObserve = valueToObserve, tabName = tabName)
        if (length(rea.values$restoreValues[[valueToObserve]]) == 0) {
            .updateValue(session, rea.values, tabName = tabName,
                         valueToObserve = valueToObserve, nextValueToObserve = nextValueToObserve ,
                         nextPatternToObserve = nextPatternToObserve)
        }
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values[[reavalToObserve]], {
        if (!rea.values$restoring) {
            destroyThisObserver <- TRUE
        } else if (length(rea.values$restoreValues[[valueToObserve]]) >= rea.values$restoreValues$counter) {
            .processValue(session, rea.values, valueToObserve, tabName = tabName)
        } else if (length(rea.values$restoreValues[[valueToObserve]]) < rea.values$restoreValues$counter) {
            .updateValue(session, rea.values, tabName = tabName,
                         valueToObserve = valueToObserve, nextValueToObserve = nextValueToObserve ,
                         nextPatternToObserve = nextPatternToObserve)
            destroyThisObserver <- TRUE
        }
    }, ignoreInit = TRUE, once = destroyThisObserver)
}


#' @keywords internal
#' @noRd
.observeRestoreInte <- function(session, rea.values) {

    observeEvent(rea.values$restoreValues$mixOmics, {
        if (is.null(rea.values$stateInput[["mixomicsSetting-run_prep"]]) ||
            rea.values$stateInput[["mixomicsSetting-run_prep"]] < 1) {
            print("here, mixOmics bis")
            rea.values$restoreValues$MOFA <- "now2"
        }
        print(rea.values$restoreValues$mixOmics)
        shinyjs::js$showmenuItem("omicsintegration")
        updateTabItems(session, inputId = "tabs", selected = "Data Integration")
        updateTabItems(session, inputId = "tabs", selected = "withMixOmics")
        updateTabsetPanel(session, inputId = "integrationMenu-mixOmics", selected = "datavarsel")
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values$preparedfor_mixOmics, {
        print("Im in the prepared for mixOmics")
        shinyjs::js$showmenuItem("OmicsIntegration")
        updateTabItems(session, inputId = "tabs", selected = "Data Integration")
        updateTabItems(session, inputId = "tabs", selected = "withMixOmics")
        updateTabsetPanel(session, inputId = "integrationMenu-mixOmics", selected = "dataint")
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values$with_mixOmics, {
        rea.values$restoreValues$MOFA <- "now"
    }, ignoreInit = TRUE, once = TRUE)

    observe({
        if ("MOFASetting-run_prep" %in% names(rea.values$stateInput) &&
            rea.values$stateInput[["mofaSetting-run_prep"]] > 0) {
            rea.values$restoreValues$MOFA <- "now"
        }
    })

    observeEvent(rea.values$restoreValues$MOFA , {
        updateTabItems(session, inputId = "tabs", selected = "Data Integration")
        updateTabItems(session, inputId = "tabs", selected = "withMOFA")
        updateTabsetPanel(session, inputId = "integrationMenu", selected = "datavarsel")
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values$preparedfor_MOFA, {
        updateTabItems(session, inputId = "tabs", selected = "OmicsIntegration")
        updateTabItems(session, inputId = "tabs", selected = "withMOFA")
        updateTabsetPanel(session, inputId = "integrationMenu-MOFA", selected = "dataint")
    }, ignoreInit = TRUE, once = TRUE)

}

#' @keywords internal
#' @noRd
.observeDiff <- function(session, rea.values){
    destroyDiffObserver <- FALSE
    observeEvent(rea.values$restoreValues$Diff, {
        .initValue(session, rea.values, valueToObserve = "Diff", tabName = "Differential analysis")
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values$datasetDiff, {
        if (!rea.values$restoring) {
            destroyDiffObserver <- TRUE
        } else if (length(rea.values$restoreValues$Diff) >= rea.values$restoreValues$counter) {
            .processValue(session, rea.values, "Diff", tabName = "Differential analysis")
        } else if (length(rea.values$restoreValues$Diff) < rea.values$restoreValues$counter) {
            .updateValue(session, rea.values, tabName = "Differential analysis",
                         valueToObserve = "Diff", nextValueToObserve = "CoEx" ,
                         nextPatternToObserve = "-runCoSeq$")
            destroyDiffObserver <- TRUE
        }
    }, ignoreInit = TRUE, once = destroyDiffObserver)

}

#' @keywords internal
#' @noRd
.observeCoex <- function(session, rea.values){
    destroyCoExObserver <- FALSE
    observeEvent(rea.values$restoreValues$CoEx, {
        .initValue(session, rea.values, valueToObserve = "CoEx", tabName = "Co-expression analysis")
    }, ignoreInit = TRUE, once = TRUE)

    observeEvent(rea.values$datasetCoEx, {
        if (!rea.values$restoring) {
            destroyCoExObserver <- TRUE
        } else if (length(rea.values$restoreValues$CoEx) >= rea.values$restoreValues$counter) {
            .processValue(session, rea.values, "CoEx", tabName = "Co-expression analysis")
        } else if (length(rea.values$restoreValues$CoEx) < rea.values$restoreValues$counter) {
            .updateValue(session, rea.values, tabName = "Co-expression analysis",
                         valueToObserve = "CoEx", nextValueToObserve = "CustomAnnotDiff" ,
                         nextPatternToObserve = "-custom-DiffExpEnrichAnal-run$")
            destroyCoExObserver <- TRUE
        }
    }, ignoreInit = TRUE, once = destroyCoExObserver)

}