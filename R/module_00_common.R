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

            todisplay <- session$userData$FlomicsMultiAssay

            for (SE.name in names(todisplay)) {
                todisplay[[SE.name]] <-
                    todisplay[[SE.name]][, getSelectedSamples(todisplay, SE.name = SE.name)]
            }

            box(title = "Dataset overview after data processing",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = FALSE,

                renderPlot({
                    plotDataOverview(
                        todisplay,
                        omicNames = rea.values$datasetProcess
                    )
                })
            )
        })

        # summary of diff analysis on all dataset
        output$DiffSummary <- renderUI({
            if (is.null(rea.values$datasetDiff))
                return()

            box(
                title = "Summary of differential expression analyses",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,
                tagList({
                    tabPanel.list <- list(
                        tabPanel(title = "DE results",
                                 renderPlot({
                                     getDiffAnalysesSummary(
                                         session$userData$FlomicsMultiAssay, plot = TRUE)
                                 })))

                    p.list <- getAnnotAnalysesSummary(
                        session$userData$FlomicsMultiAssay,
                        from = "DiffExp",
                        matrixType = "presence"
                    )

                    if (!is.null(rea.values$datasetDiffAnnot)) {
                        tabPanel.list <-
                            c(tabPanel.list,
                              lapply(names(p.list), function(database) {
                                  tabPanel(
                                      title = paste0("ORA results from ", database),
                                      fluidRow(column(
                                          width = 12,
                                          radioButtons(
                                              inputId = session$ns(paste0(
                                                  database, "-domain.diff"
                                              )),
                                              label = "Domain",
                                              choices = names(p.list[[database]]),
                                              selected = names(p.list[[database]])[1],
                                              inline = TRUE
                                          )
                                      )),
                                      fluidRow(column(
                                          width = 12,
                                          renderPlot({
                                              p.list[[database]][[input[[paste0(database, "-domain.diff")]]]]
                                          }, height = function() {400*length(rea.values$datasetProcess)})
                                      ))
                                  )
                              })
                            )
                    }
                    do.call(what = tabsetPanel, args = tabPanel.list)
                })
            )
        })

        output$CoExSummary <- renderUI({
            if (is.null(rea.values$datasetCoEx))
                return()

            box(
                title = "Summary of Co-expression analyses",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,

                tagList({
                    tabPanel.list <-
                        list(
                            tabPanel(title = "CoExp results",
                                     renderPlot({
                                         getCoExpAnalysesSummary(
                                             session$userData$FlomicsMultiAssay)
                                     })
                            )
                        )

                    p.list <- getAnnotAnalysesSummary(
                        session$userData$FlomicsMultiAssay,
                        from = "CoExp",
                        matrixType = "presence"
                    )

                    if (!is.null(rea.values$datasetCoExAnnot)) {
                        tabPanel.list <-
                            c(tabPanel.list,
                              lapply(names(p.list), function(database) {
                                  tabPanel(
                                      title = paste0("ORA results from ", database),
                                      fluidRow(column(width = 12,
                                                      radioButtons(
                                                          inputId = session$ns(paste0(
                                                              database, "-domain.coex"
                                                          )),
                                                          label = "Domain",
                                                          choices = names(p.list[[database]]),
                                                          selected = names(p.list[[database]])[1],
                                                          inline = TRUE
                                                      )
                                      )),
                                      fluidRow(column(width = 12,
                                                      renderPlot({
                                                          p.list[[database]][[input[[paste0(database, "-domain.coex")]]]]
                                                      },
                                                      height = function() {500*length(rea.values$datasetProcess)},
                                                      width = "auto")
                                      ))
                                  )
                              }))
                    }

                    do.call(what = tabsetPanel, args = tabPanel.list)
                })
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
