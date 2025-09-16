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
