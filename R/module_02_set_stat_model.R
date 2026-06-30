### ============================================================================
### [02_set_stat_model] shiny modules
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif,

#' @importFrom shinyBS popify bsButton addPopover bsTooltip

.modGLMmodelUI <- function(id){

  ns <- NS(id)

  tagList(

    fluidRow(
      column(width= 12, uiOutput(ns("SetModelFormula")))
    ),
    fluidRow(
      column(width= 12, uiOutput(ns("SetContrasts")))
    ),
    fluidRow(
      column(width= 12, verbatimTextOutput(ns("printContrast")))
    ),
    tags$br()
  )
}


.modGLMmodel <- function(input, output, session, rea.values){

  # Construct the form to select the model
  output$SetModelFormula <- renderUI({

    validate(
      need(rea.values$loadData != FALSE, "Please load data")
    )

    box(status = "warning", width = 12, solidHeader = TRUE,
        title = "Select a model formula",

        tags$i("The proposed models are written based on the selected
                   biological and batch factors, with and
                   without interaction(s). Only the two orders interaction
                   terms between the biological factors are considered"),

        selectInput( inputId = session$ns("model.formulae"), label = "",
                     choices = rev(names(generateModelFormulae(session$userData$FlomicsMultiAssay))),
                     selectize = FALSE, size = 5),
        actionButton(session$ns("validModelFormula"),"Validate", class = "butt")
    )
  })

  # as soon as the "valid model formulae" button has been clicked
  # => The model formulae is set and the interface to select the contrasts appear
  observeEvent(input$validModelFormula, {

    rea.values$model          <- FALSE
    #rea.values$analysis       <- FALSE
    rea.values$Contrasts.Sel  <- NULL
    rea.values$datasetDiff    <- NULL
    #rea.values$datasetProcess <- NULL

    message("[RFLOMICS] # 02- Statistical setting...")
    message("[RFLOMICS] #    => model formula: ", input$model.formulae)

    # => Set the model formulae
    session$userData$FlomicsMultiAssay <-
      setModelFormula(session$userData$FlomicsMultiAssay, input$model.formulae)

    # => get list of expression contrast (hypothesis)
    # session$userData$FlomicsMultiAssay <- 
    #   generateExpressionContrast(session$userData$FlomicsMultiAssay)

    rea.values$model <- TRUE

  }, ignoreInit = TRUE)


  # => Get and Display all the contrasts
  #  => The contrasts have to be chosen
  output$SetContrasts <- renderUI({

    if (rea.values$model == FALSE) return()
    
    contrastList <-
      getSelectedContrasts(session$userData$FlomicsMultiAssay, all = TRUE)

    box(width=12, status = "warning", solidHeader = TRUE,
        title = "Select contrasts",
        tags$i("The proposed contrasts are computed according to the
                 selected model formula. For models without an interaction term,
                 all available contrasts are of the 'averaged' type.
                 You must choose the contrasts that
                 reflect the biological questions."),
        br(),
        br(),
        column(
          width = 12,
          lapply(unique(contrastList$type), function(contrastType) {

            contrastList_bis <- contrastList[contrastList$type == contrastType,]
            vect <- as.vector(contrastList_bis$contrast)
            names(vect) <-
              paste0("[",contrastList_bis$tag, "] ", 
                     as.vector(contrastList_bis$contrastName))

            pickerInput(
              inputId  = session$ns(paste0("ContrastType",contrastType)),
              label    = tags$span(style = "color: black;",
                                   paste0("Contrast type: ", contrastType)),
              choices  = vect,
              options  = list(`actions-box` = TRUE, size = 10,
                              `selected-text-format` = "count > 3"),
              multiple = TRUE,
              selected = NULL)
          })
        ),
        br(),
        actionButton(session$ns("validContrasts"), "Validate", class = "butt")
    )
  })

  # as soon as the "valid Contrasts" button has been clicked
  # => The selected contrasts are saved
  # => The load data item appears
  observeEvent(input$validContrasts, {
    
    #rea.values$analysis    <- FALSE
    rea.values$datasetDiff <- NULL
    #rea.values$datasetProcess <- NULL

    # reset analysis
    lapply(unlist(rea.values$datasetList), function(dataset){
      #rea.values[[dataset]]$process   <- FALSE
      rea.values[[dataset]]$diffAnal  <- FALSE
      rea.values[[dataset]]$coExpAnal <- FALSE
      rea.values[[dataset]]$DiffExp   <- FALSE
      rea.values[[dataset]]$diffValid <- FALSE
      rea.values[[dataset]]$DiffValidContrast <- NULL
    })
    
    #get list of selected contrast data frames with expression, name and type
    contrastList <- 
      getSelectedContrasts(session$userData$FlomicsMultiAssay, all = TRUE)
    
    contrast.sel.vec <- lapply(unique(contrastList$type), function(contrastType) {
      
      contrastList_bis <- contrastList[contrastList$type == contrastType,]
      contrastList_bis[contrastList_bis$contrast %in% 
                         input[[paste0("ContrastType",contrastType)]],]

    }) |> do.call(what = rbind)

    message("[RFLOMICS] #    => selected contrasts: ", nrow(contrast.sel.vec))

    # check if user has selected the contrasts to test
    if(nrow(contrast.sel.vec) == 0){

      showModal(modalDialog(title = "Error: no contrast selected",
                            "Please select at least one hypothesis (contrast) to test."))
    }

    ## continue only if message is true
    validate({
      need(nrow(contrast.sel.vec) != 0, message = "ok")
    })

    # define all the coefficients of selected contrasts and return a
    # contrast matrix with contrast sample name and associated coefficients
    session$userData$FlomicsMultiAssay <- 
      setSelectedContrasts(
        session$userData$FlomicsMultiAssay, 
        contrastNames = contrast.sel.vec$contrastName)
    rea.values$Contrasts.Sel <- contrast.sel.vec

    #rea.values$analysis <- TRUE

  }, ignoreInit = TRUE)

  return(input)
}

