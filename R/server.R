### ============================================================================
### [shiny server function]
### ----------------------------------------------------------------------------
# N. Bessoltane,

#' @importFrom shinyBS popify
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom magrittr "%>%"
#' @importFrom rmarkdown render
rflomicsServer <- function(input, output, session) {

    #### Increasing maximum possible size of loaded files (default is only 5MB)
    # https://stackoverflow.com/questions/18037737/how-to-change-maximum-upload-size-exceeded-restriction-in-shiny-and-save-user
    options(shiny.maxRequestSize = 100*1024^2) # 100 MB limit.

    # This is to get the desired menuItem selected initially.
    # selected=T seems not to work with a dynamic sidebarMenu.
    observeEvent(session, {
        updateTabItems(session = session, inputId = "tabs", selected = "coverPage")
    })

    jsCode <- 'shinyjs.hidemenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "none"; x.classList.remove("menu-open");};
shinyjs.showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "block"; x.classList.add("menu-open");};'

    #############################################
    # reactive value for reinitialisation of UIoutput
    rea.values <- reactiveValues(
        validate.status = 0,
        loadData = FALSE,
        model    = FALSE,
        analysis = FALSE,
        report   = FALSE,

        exampleData = NULL,

        datasetList     = NULL,
        contrastList    = NULL,
        Contrasts.Sel   = NULL,
        datasetProcess  = NULL,
        datasetDiff     = NULL, # list of dataset names with diff results
        datasetCoEx     = NULL,
        datasetDiffAnnot= NULL,
        datasetCoExAnnot= NULL,

        preparedfor_mixOmics = FALSE,
        preparedfor_MOFA = FALSE,

        with_mixOmics = FALSE,
        with_MOFA = FALSE,

        restoring       = FALSE
    )

    #############################################
    # dynamic sidebar menu #
    output$mysidebar <- renderUI({

        tagList(
            useShinyjs(),
            shinyjs::extendShinyjs(text = 'shinyjs.hidemenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "none"; x.classList.remove("menu-open");};
                         shinyjs.showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "block"; x.classList.add("menu-open");};',
                                   functions = c("hidemenuItem", "showmenuItem")),
            sidebarMenu(id = "tabs",
                        menuItem(text = "Welcome", tabName = "coverPage",
                                 icon = icon('dna'), selected = TRUE),
                        # menuItem(text = "Glossary page", tabName = "GlossaryPage",
                        #          icon = icon("address-book")),
                        menuItem(text = "Restore state", tabName = "loadState",
                                 icon = icon("upload")),
                        menuItem(text = "Load Data", tabName = "importData",
                                 icon = icon('download')),
                        menuItemOutput(outputId = "SetUpModelMenu"),
                        menuItemOutput(outputId = "omics"),
                        menuItemOutput(outputId = "omicsSumUI"),
                        menuItemOutput(outputId = "Integration")
            ),

            tags$br(),
            tags$br(),
            uiOutput("runReport"),
            tags$br(),
            tags$br(),
            uiOutput("downloadResults"),
            tags$br(),
            tags$br(),
            uiOutput("saveState"),

        )
    })

    # ---- ObserveEvent for restoring session ----

    reloadObserver <- observe({

        if (is.null(rea.values$restoring) || !rea.values$restoring) {
            reloadObserver$destroy()
        } else {
            observeEvent(rea.values$restoring,{
                if (rea.values$restoring)
                    updateTabItems(session, inputId = "tabs", selected = "importData")
            },
            once = TRUE, label = "restoreLoadData"
            )

            observeEvent(rea.values$loadData, {
                if (rea.values$restoring && rea.values$loadData && !rea.values$model)
                    updateTabItems(session, inputId = "tabs", selected = "SetUpModel")
            }, label = "restoreSetup",
            ignoreInit = TRUE, once = TRUE
            )

            observeEvent(rea.values$model, {
                rea.values$restoreValues <- list()
                processed <- rea.values$stateInput[grep("([0-9]+)-run$", names(rea.values$stateInput))]
                processed <- processed[which(processed > 0)]
                shinyjs::js$showmenuItem("omicsanalysis")
                rea.values$restoreValues$Processed <- gsub("-run$", "", names(processed))
                rea.values$restoreValues$counter <- 1
            }, ignoreInit = TRUE, once = TRUE)

            .observePreProcess(session, rea.values)
            .observeRestore(session, rea.values, reavalToObserve = "datasetDiff",
                            valueToObserve = "Diff", tabName = "Differential analysis",
                            nextValueToObserve = "CoEx", nextPatternToObserve = "-runCoSeq$")
            # .observeRestore(session, rea.values, reavalToObserve = "datasetCoEx",
            #                 valueToObserve = "CoEx", tabName =  "Co-expression analysis",
            #                 nextValueToObserve = "CustomAnnotDiff",
            #                 nextPatternToObserve = "-custom-DiffExpEnrichAnal-run$")

            .observeRestore(session, rea.values, reavalToObserve = "datasetCoEx",
                            valueToObserve = "CoEx", tabName =  "Co-expression analysis",
                            nextValueToObserve = "mixOmics",
                            nextPatternToObserve = "mixomicsSetting-run_prep$")

            .observeRestoreInte(session, rea.values)

            # TODO :
            # - Always go back to importData so that this can be properly loaded
            # updateTabItems(session, inputId = "tabs", selected = "importData")
            #
            # - Message the user at the end of the restoration:
            # showModal(modalDialog(title = "Restoration complete",
            #                       "Restored session is still in test, you may have to run all annotation manually, as well as data integration.
            #                   To do so, just click on the tabs, the restoration process will activate.
            #                   It is not advised to change anything in the load Data panel."))
            # - these two needs a new reactive value...
            #
            # Fixthis:
            # - Not able to navigate in annotation tabsetpanel (nested ones)
            # - Not able to make the connection between coexpression and integration
            #   although going from mixOmics to MOFA and between data
            #   selection and integration panels is done correctly
        }

    })

    #### Item for each omics #####
    # display omics Item
    output$omics <- renderMenu({

        validate({
            need(rea.values$analysis == TRUE, message="")
        })

        menuItem(text = "Omics Analysis", tabName = "OmicsAnalysis",
                 id = "omicsanalysis",
                 icon = icon('chart-area'),
                 shinyjs::useShinyjs(),
                 shinyjs::extendShinyjs(text = 'shinyjs.hidemenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "none"; x.classList.remove("menu-open");};
shinyjs.showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "block"; x.classList.add("menu-open");};', functions = c("hidemenuItem", "showmenuItem")),
                 list(
                     lapply(names(rea.values$datasetList), function(omics){

                         lapply(names(rea.values$datasetList[[omics]]), function(i){
                             menuSubItem(text = rea.values$datasetList[[omics]][[i]],
                                         tabName = paste0(omics, "Analysis", i))
                         })
                     })
                 )
        )
    })

    output$omicsSumUI <- renderMenu({

        validate(
            need(rea.values$analysis == TRUE && length(rea.values$datasetProcess) >= 2,
                 message = "")
        )

        menuItem(text = "Compare Analyses", tabName = "omicsSum",
                 icon = icon('code-compare')) # circle-nodes
    })

    #### Item for each data integration tools #####
    # display tool Item
    output$Integration <- renderMenu({

        validate({
            need(rea.values$analysis == TRUE && length(rea.values$datasetProcess) >= 2,
                 message = "")
        })

        menuItem(text = "Data Integration", tabName = "OmicsIntegration",
                 id = "omicsintegration",
                 shinyjs::useShinyjs(),
                 shinyjs::extendShinyjs(text = 'shinyjs.hidemenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "none"; x.classList.remove("menu-open");};
shinyjs.showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "block"; x.classList.add("menu-open");};', functions = c("hidemenuItem", "showmenuItem")),
                 icon = icon('network-wired'), startExpanded = FALSE,selected = FALSE,
                 menuSubItem(text = "with MOFA", tabName = "withMOFA" ),
                 menuSubItem(text = "with MixOmics", tabName = "withMixOmics")
        )
    })

    #### Item for report #####
    output$runReport <- renderUI({
        if(is.null(rea.values$datasetProcess)) return()

        # if(is.null(rea.values$datasetProcess) ||
        #    length(rea.values$datasetProcess) !=
        #    length(unlist(rea.values$datasetList))) return()

        column(
            width = 12,
            downloadButton(outputId = "report",
                           label = "Generate report", class = "butt"))
    })

    #### Item to download Results #####
    output$downloadResults <- renderUI({
        if(is.null(rea.values$datasetProcess)) return()

        column(
            width = 12,
            downloadButton(outputId = "download",
                           label = "Download results", class = "butt")
        )
    })

    #### Item to save state ####
    output$saveState <- renderUI({
        if (!rea.values$loadData) return()
        column(
            width = 12,
            bookmarkButton(label =  'Save State')
        )
    })

    #############################################
    # dynamic content #
    output$mycontent <- renderUI({

        items.list <- list()

        for(omics in c("RNAseq", "proteomics", "metabolomics")){

            items.list[[omics]] <-  lapply(1:10, function(i){
                tabItem(tabName = paste0(omics, "Analysis", i),
                        uiOutput(paste0(omics, "AnalysisUI", i)))
            })
        }

        itemsOmics <- purrr::reduce(items.list, c)

        items <- c(

            list(

                #### Cover Page        ####
                ###########################
                tabItem(tabName = "coverPage",
                        coverPageUI()

                ),
                #### Cover Page        ####
                ###########################
                # tabItem(tabName = "GlossaryPage",
                #         GlossaryPageUI()
                #
                # ),
                #### Load State        ####
                ###########################
                tabItem(tabName = "loadState",
                        .modLoadStateUI("staterefresh")
                ),
                #### Import data       ####
                ###########################
                tabItem(tabName = "importData",

                        .modLoadDataUI("data")
                ),

                #### Set Up statistical model & hypothesis ####
                ###############################################
                tabItem(tabName = "SetUpModel",

                        .modGLMmodelUI("model")
                )
            ),
            #### omics analysis ####
            ########################
            itemsOmics,
            list(

                #### analysis summary ####
                ########################
                tabItem(tabName = "omicsSum",
                        uiOutput(outputId = "omicsSum_UI")
                ),

                #### MOFA ####
                ########################
                tabItem(tabName = "withMOFA",
                        uiOutput(outputId = "withMOFA_UI")
                ),
                #### MixOmics ####
                ########################
                tabItem(tabName = "withMixOmics",
                        uiOutput(outputId = "withMixOmics_UI")
                )
            )

        )
        do.call(tabItems, items)
    })

    observe({

        lapply(names(rea.values$datasetList), function(omics){

            lapply(names(rea.values$datasetList[[omics]]), function(i){

                switch (omics,
                        "RNAseq" = {
                            output[[paste0("RNAseqAnalysisUI", i)]] <- renderUI({

                                tabsetPanel(
                                    shinyjs::useShinyjs(),
                                    shinyjs::extendShinyjs(text = jsCode, functions = c("hidemenuItem", "showmenuItem")),
                                    #### Data Exploratory & QC ####
                                    ###############################
                                    tabPanel(title = "Pre-processing",
                                             tags$br(),
                                             tags$br(),

                                             QCNormalizationTabUI(paste0("RNAseq",i))

                                    ),

                                    #### Diff analysis  ####
                                    ######################################
                                    tabPanel("Differential analysis",
                                             tags$br(),
                                             tags$br(),
                                             DiffExpAnalysisUI(paste0("RNAseq",i))
                                    ),
                                    #### Co-expression analysis  ####
                                    ######################################
                                    tabPanel("Co-expression analysis",
                                             tags$br(),
                                             tags$br(),
                                             CoSeqAnalysisUI(paste0("RNAseq",i))
                                             #verbatimTextOutput("Asuivre")
                                    ),
                                    #### enrichment analysis  CPR ####
                                    ######################################
                                    tabPanel("Annotation Enrichment",
                                             tags$br(),
                                             tags$br(),
                                             .modEnrichmentUI(paste0("RNAseq",i))
                                    ),
                                    id = paste0(omics, i)
                                )
                            })},
                        "proteomics" = {
                            output[[paste0("proteomicsAnalysisUI", i)]] <- renderUI({
                                tabsetPanel(
                                    shinyjs::useShinyjs(),
                                    shinyjs::extendShinyjs(text = jsCode, functions = c("hidemenuItem", "showmenuItem")),

                                    #### Data Exploratory & QC ####
                                    ###############################
                                    tabPanel("Pre-processing",
                                             tags$br(),
                                             tags$br(),

                                             QCNormalizationTabUI(paste0("proteomics",i))
                                    ),
                                    #### Diff analysis  ####
                                    ######################################
                                    tabPanel("Differential analysis",
                                             tags$br(),
                                             tags$br(),
                                             DiffExpAnalysisUI(paste0("proteomics",i))
                                    ),
                                    #### Co-expression analysis  ####
                                    ######################################
                                    tabPanel("Co-expression analysis",
                                             tags$br(),
                                             tags$br(),
                                             CoSeqAnalysisUI(paste0("proteomics",i))
                                    ),
                                    ### enrichment analysis CPR ####
                                    #####################################
                                    tabPanel("Annotation Enrichment",
                                             tags$br(),
                                             tags$br(),
                                             .modEnrichmentUI(paste0("proteomics",i))
                                    ),
                                    id = paste0(omics, i)
                                )
                            })},
                        "metabolomics" = {
                            output[[paste0("metabolomicsAnalysisUI", i)]] <- renderUI({

                                tabsetPanel(
                                    shinyjs::useShinyjs(),
                                    shinyjs::extendShinyjs(text = jsCode, functions = c("hidemenuItem", "showmenuItem")),

                                    #### Data Exploratory & QC ####
                                    ###############################
                                    tabPanel("Pre-processing",
                                             tags$br(),
                                             tags$br(),

                                             QCNormalizationTabUI(paste0("metabolomics",i))
                                    ),#### Diff analysis  ####
                                    ######################################
                                    tabPanel("Differential analysis",
                                             tags$br(),
                                             tags$br(),
                                             DiffExpAnalysisUI(paste0("metabolomics",i))
                                    ),
                                    #### Co-expression analysis  ####
                                    ######################################
                                    tabPanel("Co-expression analysis",
                                             tags$br(),
                                             tags$br(),
                                             CoSeqAnalysisUI(paste0("metabolomics",i))
                                    ),
                                    ### enrichment analysis CPR ####
                                    #####################################
                                    tabPanel("Annotation Enrichment",
                                             tags$br(),
                                             tags$br(),
                                             .modEnrichmentUI(paste0("metabolomics",i))
                                    ),
                                    id = paste0(omics, i)
                                )
                            })},
                )
            })
        })
    })

    #### analysis Summary ####
    ###############################
    output$omicsSum_UI <- renderUI({

        .modSingleOmicAnalysesSummaryUI("omics")
    })

    #### MOFA data integration ####
    ###############################
    output$withMOFA_UI <- renderUI({

        .modIntegrationAnalysisUI("mofaSetting", method = "MOFA")

    })

    #### MixOmics data integration ####
    ###################################
    output$withMixOmics_UI <- renderUI({

        .modIntegrationAnalysisUI("mixomicsSetting", method = "mixOmics")

    })

    ########################################################################
    ######################### MAIN #########################################

    ##########################################
    # Part0 : presentation page
    ##########################################

    # # dir
    # shinyDirChoose(input = input, id = 'dir0', roots = c(home = '~'))
    # output$filepaths <- renderPrint({parseDirPath(roots = c(home = '~'),
    # selection = input$dir0)})

    ##########################################
    # Part1 : load data
    ##########################################

    # load omics data and experimental design
    # set reference
    # set type of factor (bio/batch)
    # check design (complete and balanced)
    inputData <- callModule(.modLoadData, "data", rea.values)

    .modLoadState("staterefresh")

    ##########################################
    # Part2 : Set GLM model
    ##########################################

    # display set up model Item
    # if no error message
    inputModel <- list()
    observeEvent(rea.values[["loadData"]], {

        #continue only if message is true or warning
        validate({
            need(rea.values$validate.status == 0, message = "set design step failed")
        })

        # display design menu
        output$SetUpModelMenu <- renderMenu({

            validate({
                need(rea.values$loadData == TRUE, message = FALSE)
            })
            menuItem(text = "Experimental Design", tabName = "SetUpModel",
                     icon = icon('vials'))
        })

    }, ignoreInit = TRUE)

    # set GLM model
    # and select list of contrast to test
    {
        inputModel <- callModule(.modGLMmodel, "model", rea.values)
    }



    ##########################################
    # Part3 : ANALYSE
    ##########################################


    # for each omics data type
    # and for each dataser
    # if no error message
    observe({

        ##########################################
        # Part3 : Data Exploratory
        ##########################################
        lapply(names(rea.values$datasetList), function(omics){

            lapply(names(rea.values$datasetList[[omics]]), function(i){

                rea.values[[rea.values$datasetList[[omics]][[i]]]] <- reactiveValues(
                    process    = FALSE,
                    diffAnal   = FALSE,
                    diffValid  = FALSE,
                    coExpAnal  = FALSE,
                    diffAnnot  = FALSE,
                    coExpAnnot = FALSE,

                    compCheck  = TRUE,
                    message    = "",

                    DiffValidContrast = NULL,
                    CoExpClusterNames = NULL,
                    omicsType = omics
                )

                ##########################################
                # Part3 : Data Exploratory
                ##########################################
                inputNorm <- callModule(
                    module  = QCNormalizationTab, id = paste0(omics,i),
                    dataset = metadata(session$userData$FlomicsMultiAssay)$omicList[[omics]][[i]],
                    rea.values = rea.values)

                ##########################################
                # Part5 :  Diff Analysis
                ##########################################
                inputDiff <- callModule(
                    module  = DiffExpAnalysis, id = paste0(omics, i),
                    dataset = metadata(session$userData$FlomicsMultiAssay)$omicList[[omics]][[i]],
                    rea.values = rea.values)

                ##########################################
                # Part6 : Co-Expression Analysis
                ##########################################
                callModule(
                    module  = CoSeqAnalysis, id = paste0(omics, i),
                    dataset = metadata(session$userData$FlomicsMultiAssay)$omicList[[omics]][[i]],
                    rea.values = rea.values)

                ##########################################
                # Part7 : Enrichment Analysis CPR
                ##########################################
                callModule(
                    module  = .modEnrichment, id = paste0(omics, i),
                    dataset = metadata(session$userData$FlomicsMultiAssay)$omicList[[omics]][[i]],
                    rea.values = rea.values)

            })
        })

    })

    callModule(module = .modSingleOmicAnalysesSummary, id = "omics",
               rea.values = rea.values)
    callModule(module = .modIntegrationAnalysis, id = "mixomicsSetting",
               rea.values = rea.values, method = "mixOmics")
    callModule(module = .modIntegrationAnalysis, id = "mofaSetting",
               rea.values = rea.values, method = "MOFA")


    ##########################################
    # Part8 : RMD REPORT
    ##########################################

    output$report <- downloadHandler(
        # For PDF output, change this to "report.pdf"

        filename = function(){
            projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
            paste0(format(Sys.time(), "%Y_%m_%d"), "_", projectName, ".html")
        },
        content = function(file) {

            withProgress(
                message = 'Download in progress',
                detail = 'This may take a while...', value = 0, {

                    incProgress(0.2)
                    projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
                    outDir <- file.path(tempdir(),
                                        paste0(format(Sys.time(),"%Y_%m_%d"),"_",
                                               projectName))

                    dir.create(path = outDir, showWarnings=FALSE)

                    incProgress(0.3)
                    generateReport(object = session$userData$FlomicsMultiAssay,
                                   reportName = file)

                    incProgress(0.9)
                    owd <- setwd(tempdir())
                    on.exit(setwd(owd))
                })
        }
    )

    ##########################################
    # Part9 : Download results as an archive
    ##########################################

    output$download <- downloadHandler(
        # For PDF output, change this to "report.pdf"

        filename = function(){
            projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
            paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName, ".tar.gz")
        },
        content = function(file) {
            withProgress(
                message = 'Download in progress',
                detail = 'This may take a while...', value = 0, {

                    incProgress(0.2)
                    projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
                    outDir <-
                        file.path(tempdir(),
                                  paste0(format(Sys.time(),"%Y_%m_%d"), "_", projectName))

                    dir.create(path = outDir, showWarnings=FALSE)

                    incProgress(0.3)
                    generateReport(object = session$userData$FlomicsMultiAssay,
                                   archiveName = file)

                    incProgress(0.9)
                    owd <- setwd(tempdir())
                    on.exit(setwd(owd))

                })
        })


    # ---- Restore ----

    onRestore(function(state) {
        rea.values$stateDir <- state[["dir"]]
        rea.values$restoring <- TRUE
        rea.values$stateInput <- state[["input"]]
    })

    # onRestored(function(state) {
    #     message("RESTORED STATE Server restoration is completed")
    # })


}

