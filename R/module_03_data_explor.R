### ============================================================================
### [03_data_processing] shiny modules
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif, 

### UI
QCNormalizationTabUI <- function(id) {
  #name space for id
  ns <- NS(id)
  tagList(fluidRow(
    box(
      title = span(
        tagList(
          icon("filter"),
          "Data exploration and pre-processing",
          tags$small("(Scroll down for instructions)")
        )
      ),
      width = 12,
      status = "warning",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = TRUE,
      data_explor_docs
    )
  ),
  fluidRow(
    column(
      width = 4,
      box(
        width = 14,
        title = span(tagList(icon("sliders-h"), "  ", "Setting")),
        status = "warning",
        uiOutput(ns("selectSamplesUI")),
        plotOutput(ns("completenessUI")),
        hr(),
        uiOutput(ns("paramUI"))
      ),
      fluidRow(uiOutput(ns(
        "filtSummary1UI"
      ))),
      fluidRow(uiOutput(ns(
        "filtSummary2UI"
      )))
    ),
    
    column(8, uiOutput(
      ns("tabPanelUI")
    ))))
}


### server
QCNormalizationTab <-
  function(input, output, session, dataset, rea.values) {
    
    #---- sample list----
    output$selectSamplesUI <- renderUI({
      sampleList <-
        colnames(session$userData$FlomicsMultiAssay[[dataset]])
      pickerInput(
        inputId  = session$ns("selectSamples"),
        label    = 
          .addBSpopify(label = 'Samples list:', 
                       content = "Select samples to include in further analyses"),
        choices  = sampleList,
        options  = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = TRUE,
        selected = sampleList
      )
    })
    
    #---- Completeness----
    output$completenessUI <- renderPlot({
      plotExpDesignCompleteness(
        object = session$userData$FlomicsMultiAssay[[dataset]],
        sampleList = input$selectSamples
      )
    })
    
    #---- adapted parameters for each omics type ----
    output$paramUI <- renderUI({
      
      paramFeatureFilter <-
        switch (
          getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]]),
          "RNAseq" = {
            list(
              radioButtons(
                inputId  = session$ns("selectFilterMethod"),
                label    = 
                  .addBSpopify(
                    label = 'Low count filtering', 
                    content = "Genes with low counts will be removed based on count per million (cpm), accounting for the library size"),
                choices  =  list("CPM" = "CPM"),
                selected = "CPM"
              ),
              conditionalPanel(
                condition = 
                  paste0("input[\'", 
                         session$ns("selectFilterMethod"),
                         "\'] == \'CPM\'"),
                selectInput(
                  inputId  = session$ns("Filter_Strategy"),
                  label    = .addBSpopify(label = 'Strategy', 
                                          content = "Choose the strategy to filter genes based on count per million (cpm). Keep genes if the NbOfsample_over_cpm >= Strategy."),
                  choices  = c("NbConditions" = "NbConditions",
                               "NbReplicates" = "NbReplicates"),
                  selected = "NbReplicates"
                ),
                style="font-size:80%; font-family:Arial;"
              ),
              conditionalPanel(
                condition = 
                  paste0("input[\'", 
                         session$ns("selectFilterMethod"),
                         "\'] == \'CPM\'"),
                numericInput(
                  inputId = session$ns("FilterSeuil"),
                  label = .addBSpopify(label = 'CPM cut-off', 
                                       content = "Choose the cpm cut-off"),
                  value = 1, min = 1, max = 10, step = 1
                ),
                style="font-size:80%; font-family:Arial;"
              ),
              hr()
            )
          },
          {
            list(
              radioButtons(
                inputId  = session$ns("selectImputMethod"),
                label    = 
                  .addBSpopify(
                    label = 'Missing Value Imputation', 
                    content = paste0("Imputation method, cannot be changed for ",
                                     getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]]),
                                     " data.")),
                choices  =  list("MVI" = "MVI"),
                selected = "MVI"
              ),
              hr()
            )
          }
        )
      
      paramTransform <- 
        switch (
          getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]]),
          "RNAseq" = {},
          {
            list(
              radioButtons(
                inputId  = session$ns("dataTransform"),
                label    = 
                  .addBSpopify(label = 'Data Transformation', 
                               content = paste0("Choose transformation method. ", 
                                                "If the data is already transformed, ",
                                                "please specify the method, tool, or ",
                                                "platform used for the processing.")),
                choices  = c("log2" = "log2", 
                             "Already transformed" = "none"),
                selected = "log2"
              ),
              conditionalPanel(
                condition = 
                  paste0("input[\'", 
                         session$ns("dataTransform"),
                         "\'] == \'none\'"),
                textInput(inputId = session$ns("userTransMethod"), 
                          label = "by:", value = "unknown"),
                style="font-size:80%; font-family:Arial;"
              ),
              hr()
            )
          }
        )
      
      paramNormalization <-
        switch (
          getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]]),
          "RNAseq" = {
            list(
              radioButtons(
                inputId  = session$ns("selectNormMethod"),
                label    = 
                  .addBSpopify(label = 'Data Normalization', 
                               content = "Normalization method, cannot be changed for counts data."),
                choices  =  list("TMM (edgeR)" = "TMM"),
                selected = "TMM"
              ),
              hr()
            )
          },
          {
            list(       
              radioButtons(
                inputId = session$ns("selectProtMetNormMethod"),
                label   = 
                  .addBSpopify(label = 'Data Normalization', 
                               content = paste0("Choose normalization method. ",
                                                "If the data is already normalized, ",
                                                "please specify the method, tool, or ",
                                                "platform used for the processing.")),
                choices  = c("median" = "median",
                             "totalSum" = "totalSum",
                             "Already normalized" = "none"),
                selected = "median"
              ),
              conditionalPanel(
                condition = 
                  paste0("input[\'", 
                         session$ns("selectProtMetNormMethod"),
                         "\'] == \'none\'"),
                textInput(inputId = session$ns("userNormMethod"), 
                          label = "by:", value = "unknown")
              ),
              hr()
            )
          }
        )
      
      fluidRow(
        column(
          width = 12,
          paramFeatureFilter,
          paramTransform,
          paramNormalization,
          actionButton(session$ns("run"), "Run", class = "butt")
        )
      )
    })
    
    #### PCA axis for plot
    # update/adapt PCA axis
    callModule(UpdateRadioButtons, "factors")
    
    #### PCA for metadata axis for plot
    # update/adapt PCA axis
    callModule(UpdateRadioButtons, "meta")
    
    #---- dataset filtering summary----
    output$filtSummary1UI <- renderUI({
      SE.data  <-
        session$userData$FlomicsMultiAssay[[dataset]]
      
      tagList(
        box(
          title = length(names(SE.data)),
          width = 6,
          background = "maroon",
          paste0("Initial number of ", .omicsDic(SE.data)$variableName)
        ),
        box(
          title = length(colnames(SE.data)),
          width = 6,
          background = "light-blue",
          "Initial number of samples"
        )
      )
      
    })
    
    #---- Summary UI----
    output$filtSummary2UI <- renderUI({
      if (rea.values[[dataset]]$process == FALSE)
        return()
      
      SE.data  <- session$userData$FlomicsMultiAssay[[dataset]]
      
      tagList(
        box(
          title = length(names(SE.data)),
          width = 6,
          background = "fuchsia",
          paste0("Number of filtered ", 
                 .omicsDic(SE.data)$variableName)
        ),
        box(
          title = length(colnames(SE.data)),
          width = 6,
          background = "purple",
          "Number of filtered samples"
        )
      )
      
    })
    
    
    #---- Exploratory of Biological and Technical variability----
    output$tabPanelUI <- renderUI({
      if (rea.values$model == FALSE)
        return()
      
      MAE.data <- session$userData$FlomicsMultiAssay
      SE.data <- MAE.data[[dataset]]
      
      tabPanel.default.list <- list(
        tabPanel(
          "Distribution (boxplots)",
          tags$br(),
          tags$i(
            "It is expected to observe aligned boxplots/medians after running 
                        the pre-processing steps. Samples with shifted median may 
                        be outliers. Consider removing them."
          ),
          tags$br(),
          tags$hr(),
          uiOutput(session$ns("boxplotUI"))
        ),
        tabPanel(
          "Distribution (density)",
          tags$br(),
          tags$i(
            "It is expected to have a gaussian-like density distribution 
                      after running  the pre-processing steps.
                      If a second peak is observed at the beginning of the 
                      curve, this could indicate that the low, uninformative
                      values, have not been correctly filtered. 
                      In this case, you may 
                      increase the filtering threshold."
          ),
          tags$br(),
          tags$hr(),
          uiOutput(session$ns("CountDistUI"))
        ),
        tabPanel(
          "Principal component analysis",
          tags$br(),
          tags$i(
            "To observe the variability associated to each biological factor, 
               you can change the colour of samples according to 
               experimental factor (Factor button).
               It will help you to interpret the PCA axes. You may identify 
               outliers samples that drives the variability. Consider removing
               them from the analysis.
               In the best case scenario, biological replicates have to group 
               together with superposed ellipses. 
            If not, it may indicate batch effect."
          ),
          tags$br(),
          tags$hr(),
          uiOutput(session$ns("PCAcoordUI"))
        ),
        ### boxplot DE ###
        tabPanel(
          title = paste0("Feature profiles (boxplots)"),
          tags$br(),
          tags$i(paste0("Boxplot showing the expression profile of a selected ",
                        .omicsDic(SE.data)$variableName),
                 " colored by experimental factor's levels."),
          
          tags$br(), tags$hr(), tags$br(),
          uiOutput(session$ns("FeatureBoxplot"))
          
        )
      )
      
      tabPanel.list <- list()
      
      switch(
        getOmicsTypes(SE.data),
        "RNAseq" = {
          tabPanel.list <- c(list(
            tabPanel(
              "Library size",
              tags$br(),
              tags$i(
                "It is expected that the library sizes are equals or at least 
                close to each other after running the pre-processing steps."
              ),
              tags$br(),
              tags$hr(),
              uiOutput(session$ns("LibSizeUI"))
            )
          ), tabPanel.default.list)
          
        },
        "proteomics" = {
          tabPanel.list <- tabPanel.default.list
        },
        "metabolomics" = {
          tabPanel.list <- tabPanel.default.list
        }
      )
      
      # Exploratory of Biological and Technical variability
      box(
        width = 14,
        title = "Exploratory of Biological and Technical Variability",
        solidHeader = TRUE,
        status = "warning",
        do.call(what = tabsetPanel, args = tabPanel.list)
      )
    })
    
    #---- QC plot for raw/processed data----
    # library size plot only for RNAseq data
    output$LibSizeUI <- renderUI({
      plot <- renderPlot(
        plotLibrarySize(session$userData$FlomicsMultiAssay[[dataset]], raw = TRUE))
      
      if (rea.values[[dataset]]$process != FALSE) {
        plot <- list(renderPlot(
          plotLibrarySize(session$userData$FlomicsMultiAssay[[dataset]])
        ), plot)
      }
      
      return(plot)
    })
    
    # value (count/intensity) distribution (boxplot/density)
    output$boxplotUI <- renderUI({
      plot <- renderPlot(
        plotDataDistribution(
          session$userData$FlomicsMultiAssay[[dataset]],
          plot = "boxplot",
          raw = TRUE
        )
      )
      
      if (rea.values[[dataset]]$process != FALSE) {
        plot <- list(renderPlot(
          plotDataDistribution(
            session$userData$FlomicsMultiAssay[[dataset]],
            plot = "boxplot",
            raw = FALSE
          )
        ), plot)
      }
      
      return(plot)
    })
    
    output$CountDistUI <- renderUI({
      plot <- renderPlot(
        plotDataDistribution(
          session$userData$FlomicsMultiAssay[[dataset]],
          plot = "density",
          raw = TRUE
        )
      )
      
      if (rea.values[[dataset]]$process != FALSE) {
        plot <- list(renderPlot(
          plotDataDistribution(
            session$userData$FlomicsMultiAssay[[dataset]],
            plot = "density",
            raw = FALSE
          )
        ), plot)
      }
      
      return(plot)
    })
    
    # PCA plot
    output$PCAcoordUI <- renderUI({
      factors_type   <- getFactorTypes(session$userData$FlomicsMultiAssay)
      choices        <- names(factors_type)
      names(choices) <-
        paste(names(factors_type),  paste0("(", factors_type, ")"))
      
      ui <- list(fluidRow(
        tags$br(),
        column(
          width = 6,
          radioButtons(
            inputId = session$ns("PCA.factor.condition"),
            label = 'Factors:',
            choices = c("groups", choices),
            selected = "groups"
          )
        ),
        column(width = 6, UpdateRadioButtonsUI(session$ns("factors"))),
        tags$br(),
      ),
      plotOutput(session$ns("raw.PCAcoord")))
      
      if (rea.values[[dataset]]$process != FALSE) {
        ui <- list(plotOutput(session$ns("norm.PCAcoord")), ui)
      }
      
      return(ui)
    })
    
    output$raw.PCAcoord <- renderPlot({
      PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
      PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
      condGroup <- input$PCA.factor.condition
      
      plotOmicsPCA(
        session$userData$FlomicsMultiAssay[[dataset]],
        raw = "raw",
        axes = c(PC1.value, PC2.value),
        groupColor = condGroup
      )
      
    })
    output$norm.PCAcoord <- renderPlot({
      if (rea.values[[dataset]]$process == FALSE)
        return()
      
      PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
      PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
      condGroup <- input$PCA.factor.condition
      
      plotOmicsPCA(
        session$userData$FlomicsMultiAssay[[dataset]],
        raw = "norm",
        axes = c(PC1.value, PC2.value),
        groupColor = condGroup
      )
      
    })
    
    # Boxplot
    output$FeatureBoxplot <- renderUI({
      
      SE.data <- session$userData$FlomicsMultiAssay[[dataset]]
      
      ui <- list(
        plotOutput(session$ns("FeatureBoxplot.raw"))
      )
      
      if (rea.values[[dataset]]$process != FALSE) {
        ui <- list(plotOutput(session$ns("FeatureBoxplot.norm")), ui)
      }
      
      ui <- list(
        fluidRow(
          column(
            width = 6,
            selectizeInput(
              inputId = session$ns("DE"),
              label = paste0("Select ",.omicsDic(SE.data)$variableName,":"),
              multiple = FALSE,
              choices = names(SE.data),
              options = list(maxOptions = 1000)
            )
          ),
          column(
            width = 6,
            radioButtons(
              inputId = session$ns("DEcondition"),
              label = 'Levels:',
              choices = c("groups",getBioFactors(SE.data)),
              selected = getBioFactors(SE.data)[1]
            )
          )
        ), ui)
      return(ui)
    })
    
    output$FeatureBoxplot.raw <- renderPlot({
      
      plotBoxplotDE(
        object=session$userData$FlomicsMultiAssay[[dataset]], 
        featureName=input$DE, 
        groupColor=input$DEcondition, 
        raw = TRUE) 
    })
    output$FeatureBoxplot.norm <- renderPlot({
      
      plotBoxplotDE(
        object=session$userData$FlomicsMultiAssay[[dataset]], 
        featureName=input$DE, 
        groupColor=input$DEcondition, 
        raw = FALSE) 
    })
    
    #---- run preprocessing - Normalization/transformation, filtering...----
    observeEvent(input$run, {
      # check if input$selectSamples is empty
      if (is.null(input$selectSamples)) {
        showModal(modalDialog(title = "Error message", 
                              "Please select some samples to run the analysis."))
      }
      validate({
        need(!is.null(input$selectSamples), 
             message = "Please select some samples to run the analysis.")
      })
      
      # get parameters
      param.list <- switch(
        getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]]),
        "RNAseq" = {
          list(
            Filter_Strategy = input$Filter_Strategy,
            CPM_Cutoff = input$FilterSeuil,
            NormMethod = input$selectNormMethod
          )
        },
        {
          list(
            transform_method = input$dataTransform,
            NormMethod       = input$selectProtMetNormMethod,
            userNormMethod   = input$userNormMethod,
            userTransMethod  = input$userTransMethod
          )
        })
      param.list <-
        c(param.list, list(samples = input$selectSamples))
      
      if (check_run_process_execution(
        session$userData$FlomicsMultiAssay,
        dataset = dataset,
        param.list = param.list
      ) == FALSE &&
      rea.values[[dataset]]$process == TRUE)
        return()
      
      # re-initialize reactive values
      rea.values[[dataset]]$process   <- FALSE
      rea.values[[dataset]]$diffAnal  <- FALSE
      rea.values[[dataset]]$coExpAnal <- FALSE
      rea.values[[dataset]]$diffAnnot <- FALSE
      rea.values[[dataset]]$diffValid <- FALSE
      rea.values[[dataset]]$DiffValidContrast <- NULL
      
      message("[RFLOMICS] # 03- Data processing: ", dataset)
      
      catch.res <-
        .tryCatch_rflomics(runDataProcessing(
          object = session$userData$FlomicsMultiAssay,
          SE.name = dataset,
          samples = input$selectSamples,
          filterStrategy = param.list[["Filter_Strategy"]],
          cpmCutoff = param.list[["CPM_Cutoff"]],
          normMethod = param.list[["NormMethod"]],
          transformMethod = param.list[["transform_method"]],
          userNormMethod = param.list[["userNormMethod"]],
          userTransMethod = param.list[["userTransMethod"]]
        ))
      
      if(!is.null(catch.res$error))
        showModal(
          modalDialog(title = "Error message", catch.res$error))
      
      validate({
        need(is.null(catch.res$error), message = catch.res$error)
      })
      
      session$userData$FlomicsMultiAssay <- catch.res$result
      message(c(catch.res$messages, catch.res$warnings))
      
      rea.values[[dataset]]$process <- TRUE
      
      # re-initialize list of diff analized dataset
      rea.values$datasetProcess <- 
        getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                analyses = "DataProcessing")
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
    
    return(input)
    
  }


############## functions ###############

# ----- check run norm execution ------
check_run_process_execution <-
  function(object.MAE, dataset, param.list = NULL) {
    
    SE <- object.MAE[[dataset]]
    
    # check filtred samples
    if (!dplyr::setequal(getSelectedSamples(SE), param.list$samples))
      return(TRUE)
    
    switch(
      getOmicsTypes(SE),
      "RNAseq" = {
        # filtering setting
        if (is.null(getFilterSettings(SE)$filterStrategy) || 
            param.list$Filter_Strategy != getFilterSettings(SE)$filterStrategy)
          return(TRUE)
        if (is.null(getFilterSettings(SE)$cpmCutoff) ||
            param.list$CPM_Cutoff != getFilterSettings(SE)$cpmCutoff)
          return(TRUE)
      },
      {
        if (is.null(getTransSettings(SE)$method) ||
            param.list$transform_method != getTransSettings(SE)$method)
          return(TRUE)
        if (getTransSettings(SE)$method == "none" &&
            (is.null(getTransSettings(SE)$suppInfo) ||
             param.list$userTransMethod != getTransSettings(SE)$suppInfo))
          return(TRUE)
      }
    )
    
    if (is.null(getNormSettings(SE)$method) ||
        param.list$NormMethod != getNormSettings(SE)$method)
      return(TRUE)
    if (getNormSettings(SE)$method == "none" &&
        (!is.null(getNormSettings(SE)$suppInfo) &&
         param.list$userNormMethod != getNormSettings(SE)$suppInfo))
      return(TRUE)
    
    return(FALSE)
  }
