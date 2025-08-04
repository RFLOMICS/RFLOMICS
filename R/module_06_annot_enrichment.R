### ============================================================================
### [06_annotation] shiny modules
### ----------------------------------------------------------------------------
# N. Bessoltane,
# A. Hulot <- shiny functions

#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom DT datatable JS
#' @importFrom plotly renderPlotly

# ---- Module enrichment : main -----

#' @keywords internal
.modEnrichmentUI <- function(id) {
    options(shiny.maxRequestSize = 3000 * 1024 ^ 2)

    #name space for id
    ns <- NS(id)

    tagList(fluidRow(
        box(
            title = span(tagList(
                icon("chart-pie"),
                "Over Representation analysis (ClusterProfiler)",
                "   " ,
                tags$small("(Scroll down for instructions)")
            )),
            solidHeader = TRUE,
            status = "warning",
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            annot_analysis_docs
        )
    ),
    ## parameters + input data
    fluidRow(uiOutput(ns(
        "tabsetPanel_DB_UI"
    ))))
}

#' @importFrom data.table fread
#' @importFrom AnnotationDbi keytypes
#' @importFrom DT renderDataTable datatable
#' @keywords internal
.modEnrichment <-
    function(input,
             output,
             session,
             dataset,
             rea.values) {
        ns <- session$ns
        local.rea.values <- reactiveValues()

        # UI
        output$tabsetPanel_DB_UI <- renderUI({
            validate(
                need(
                    rea.values[[dataset]]$diffValid != FALSE,
                    "Please run a differential analysis and validate your choices."
                )
            )

            omicsType <-
                getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]])

            tabPanel.list <- list(tabPanel(title = "Custom annotation",
                                           .modEnrichmentDBUI(id = ns("custom"))))

            if (omicsType %in% c("RNAseq", "proteomics")) {
                tabPanel.list <- c(list(
                    tabPanel(title = "Gene Onlology database",
                             .modEnrichmentDBUI(id = ns("GO"))),
                    tabPanel(title = "KEGG database",
                             .modEnrichmentDBUI(id = ns("KEGG")))
                ),
                tabPanel.list)
            }
            do.call(what = tabsetPanel, args = tabPanel.list)
        })

        # SERVEUR
        callModule(
            module = .modEnrichmentDB,
            id = "GO",
            dataset = dataset,
            database = "GO",
            rea.values
        )
        callModule(
            module = .modEnrichmentDB,
            id = "KEGG",
            dataset = dataset,
            database = "KEGG",
            rea.values
        )
        callModule(
            module = .modEnrichmentDB,
            id = "custom",
            dataset = dataset,
            database = "custom",
            rea.values
        )

    }

# ---- Module enrichment : per database ----

#' @title .modEnrichmentDBUI
#' @keywords internal
#' @noRd
.modEnrichmentDBUI <- function(id) {
    #name space for id
    ns <- NS(id)

    tagList(
        br(),
        box(
            title = span(tagList(icon("sliders"), "  ", "Settings")),
            width = 3,
            status = "warning",
            uiOutput(ns("ParamDB_UI"))
        ),
        box(
            title = NULL,
            headerBorder = FALSE,
            width = 9,
            status = "warning",
            uiOutput(ns("tabsetRes_UI"))
        )
    )

}


#' @title .modEnrichmentDB
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table fread
#' @keywords internal
#' @noRd
.modEnrichmentDB <- function(input,
                             output,
                             session,
                             dataset,
                             database,
                             rea.values) {
    ns <- session$ns
    local.rea.values <- reactiveValues(dataPathAnnot = NULL)

    # UI

    ## setting
    output$ParamDB_UI <- renderUI({
        switch(
            database,
            "GO"     = {
                uiOutput(ns("AnnotParamGO_UI"))
            },
            "KEGG"   = {
                uiOutput(ns("AnnotParamKEGG_UI"))
            },
            "custom" = {
                uiOutput(ns("AnnotParamCustom_UI"))
            }
        )

    })

    ## display results of enrichment analysis for selected database
    output$tabsetRes_UI <- renderUI({
        ### for diff analysis results
        tabPanel.list <- list(
            tabPanel(
                "Enrichment from differential expression analysis results",
                br(),
                .modRunEnrichmentUI(ns("DiffExpEnrichAnal"))
            )
        )

        ### for coexp analysis results if exist
        if (rea.values[[dataset]]$coExpAnal) {
            tabPanel.list <- c(tabPanel.list,
                               list(
                                   tabPanel(
                                       "Enrichment from co-expression analysis results",
                                       br(),
                                       .modRunEnrichmentUI(ns("CoExpEnrichAnal"))
                                   )
                               ))
        }

        do.call(what = tabsetPanel, args = tabPanel.list)
    })

    ## GO settings
    output$AnnotParamGO_UI <-
        .annotParamGO(session, rea.values, dataset)

    ## KEGG settings
    output$AnnotParamKEGG_UI <-
        .annotParamKEGG(session, rea.values, dataset)

    ## custom settings
    output$AnnotParamCustom_UI <-
        .annotParamCustom(session, rea.values, dataset)

    ### select columns
    observeEvent(input$annotationFileCPR, {
        local.rea.values$dataPathAnnot <- NULL
        local.rea.values$dataPathAnnot <-
            input$annotationFileCPR$datapath
        output$selectColumnsCustom <-
            .annotFileColumns(session, local.rea.values,
                              dataset)

    })

    # SERVEUR
    callModule(
        module     = .modRunEnrichment,
        id         = "DiffExpEnrichAnal",
        dataset    = dataset,
        database   = database,
        listSource = "DiffExp",
        rea.values = rea.values,
        input2     = input,
        local.rea.values = local.rea.values
    )

    callModule(
        module     = .modRunEnrichment,
        id         = "CoExpEnrichAnal",
        dataset    = dataset,
        database   = database,
        listSource = "CoExp",
        rea.values = rea.values,
        input2     = input,
        local.rea.values = local.rea.values
    )

}

# ---- module enrichment : run per analysis type ----

#' @title .modRunEnrichmentUI
#' @keywords internal
#' @noRd
.modRunEnrichmentUI <- function(id) {
    #name space for id
    ns <- NS(id)

    tagList(uiOutput(ns("settings")),
            uiOutput(ns("enrichSummary")),
            uiOutput(ns("AnnotResults")))
}

#' @title .modRunEnrichment
#' @keywords internal
#' @importFrom DT renderDataTable datatable
#' @importFrom grid unit
#' @importFrom AnnotationDbi keytypes
#' @importMethodsFrom ComplexHeatmap draw
#' @noRd
.modRunEnrichment <-
    function(input,
             output,
             session,
             dataset,
             database,
             listSource,
             rea.values,
             input2,
             local.rea.values) {
        ns <- session$ns

        switch(listSource,
               "DiffExp" = {
                   title <- "Summary"
                   datasetList <- "datasetDiffAnnot"
               },
               "CoExp" = {
                   title <- "Summary"
                   datasetList <- "datasetCoExAnnot"
               })

        # UI
        ## setting
        output$settings <-
            .annotSettings(session, rea.values, dataset, listSource)

        ## display results
        output$enrichSummary <- .outEnrichSummary(session,
                                                  rea.values,
                                                  input,
                                                  listSource,
                                                  title,
                                                  dataset,
                                                  database)
        output$AnnotResults <- .outAnnotResults(session,
                                                input,
                                                input2,
                                                rea.values,
                                                listSource,
                                                dataset,
                                                database)

        # SERVER

        ## run enrichment
        # ---- run Annotation  ----
        observeEvent(input$run, {

            ### GO checks ###
            if (database == "GO") {
                condition <- input2$db.select != ""
                messCond <-
                    "Select a database for the analysis. If none is suggested in
      the list, please check and install the appropriate org.db package
      from Bioconductor."
                # check parameters
                if (!condition) {
                    showModal(modalDialog(title = "No database selected", messCond))
                }
                validate({
                    need(condition, messCond)
                })

                textToEval <-
                    paste0("keytypes(",
                           input2$db.select,
                           "::",
                           input2$db.select,
                           ")")
                accepted_keytypes <- eval(parse(text = textToEval))

                # Check keytypes for GO
                condition <- input2$keytype.go %in% accepted_keytypes
                messCond <- paste("Keytype must be one of: ",
                                  paste(accepted_keytypes, collapse = ", "),
                                  sep = " ")
                if (!condition) {
                    showModal(modalDialog(title = "Error in argument", messCond))
                }
                validate({
                    need(condition, message = messCond)
                })

                domain <- input2$ont.select
                if (input2$ont.select == "ALL")
                    domain <- c("MF", "BP", "CC")

                # set list args
                paramList <- list(
                    OrgDb        = input2$db.select,
                    keyType      = input2$keytype.go,
                    pvalueCutoff = input2$pValue_GO,
                    database     = database,
                    domain       = domain
                )
            }

            ### KEGG checks ###
            if (database == "KEGG") {
                # check param
                ## ID key
                condition <-
                    input2$keyTypeKEGG %in% c("kegg",
                                              "ncbi-geneid",
                                              "ncib-proteinid",
                                              "uniprot")
                messCond  <- "Please chose an appropriate
      KEGG keytype. It must be one of kegg, ncbi-geneid,
          ncbi-proteinid or uniprot"

                if (!condition) {
                    showModal(modalDialog(title = "Error in keytype", messCond))
                }
                validate({
                    need(condition, message = messCond)
                })

                ## code organism KEGG
                condition <- input2$organism != ""
                messCond  <-
                    "It seems you forgot to enter the KEGG organism code
      (three or more letters, eg: hsa for human or
      ath for arabidopsis thaliana)"
                if (!condition) {
                    showModal(modalDialog(title = "Missing organism code", messCond))
                }
                validate({
                    need(condition, message = messCond)
                })

                # set list args
                paramList <- list(
                    organism     = input2$organism,
                    keyType      = input2$keyTypeKEGG,
                    pvalueCutoff = input2$pValueKEGG,
                    database     = database,
                    domain       = "no-domain"
                )
            }

            ### Custom checks ###
            if (database == "custom") {
                # check param
                ## if custom annotation file
                condition <- !is.null(local.rea.values$dataPathAnnot)
                messCond  <- "Please load a custom annotation file."
                if (!condition) {
                    showModal(modalDialog(title = "Missing annotation file", messCond))
                }
                validate({
                    need(condition,  message = messCond)
                })

                # Check some custom elements and annotation file
                condition <-
                    input2$col_geneName != "" && input2$col_termID != ""
                messCond <-
                    "Please choose the columns names for the omics features names/ID and the ontology terms ID."
                if (!condition) {
                    showModal(modalDialog(title = "Missing mandatory settings", messCond))
                }
                validate({
                    need(condition, message = messCond)
                })

                # set param
                domain <- NULL
                annotation2 <- NULL
                col_domain_arg <- NULL

                annotation <-
                    fread(
                        file = local.rea.values$dataPathAnnot,
                        sep = "\t",
                        header = TRUE
                    )

                annotation <- rename(annotation,
                                     gene = input2$col_geneName,
                                     term = input2$col_termID)

                if (input2$col_domain != "")
                    annotation <- rename(annotation,
                                         domain = input2$col_domain)

                if (input2$col_termName != "")
                    annotation <- rename(annotation,
                                         name = input2$col_termName)

                ##
                paramList <- list(
                    pvalueCutoff = input2$pValue_custom,
                    database     = database,
                    annotation   = annotation
                )
            }

            # prevent multiple execution
            if (.checkRunORAExecution(session$userData$FlomicsMultiAssay[[dataset]],
                                      database   = database,
                                      fromSource = listSource,
                                      paramList  = paramList) == FALSE) return()

            #---- progress bar ----#
            progress <- Progress$new()
            progress$set(message = "Run annotation enrichment", value = 0)
            on.exit(progress$close())
            progress$inc(1 / 10, detail = "Checking everything")
            #-----------------------#

            # ---- Annotation on diff results: ----
            # run analysis
            messList <- switch(listSource,
                               "DiffExp" = "differential expression",
                               "CoExp" = "co-expression clusters")

            message("[RFLOMICS] # 06- ", database,
                    " Enrichment Analysis of ", messList ,
                    " lists from ", dataset)

            #---- progress bar ----#
            progress$inc(5 / 10, detail = paste("Run analyses ", 50, "%", sep = ""))
            #----------------------#

            rea.values[[dataset]][[listSource]] <- FALSE

            # run annotation diff analysis
            paramList <- c(
                paramList,
                list(
                    object = session$userData$FlomicsMultiAssay[[dataset]],
                    from = listSource))

            session$userData$FlomicsMultiAssay[[dataset]] <-
                do.call(runAnnotationEnrichment, paramList)

            rea.values[[dataset]][[listSource]] <- TRUE

            #---- progress bar ----#
            progress$inc(1, detail = paste("All done", 100, "%", sep = ""))
            #----------------------#

            rea.values[[datasetList]] <-
                getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                        analyses = paste0(listSource, "EnrichAnal"))

        }, ignoreInit = TRUE)
    }

# ---- check run enrichement analysis execution (à adapter) ----
#' @title .checkRunORAExecution
#' @noRd
#' @keywords internal
.checkRunORAExecution <-
    function(object.SE,
             database   = NULL,
             fromSource = NULL,
             paramList  = NULL) {

        dbRes <- getAnalysis(object.SE,
                             name = paste0(fromSource, "EnrichAnal"),
                             subName = database)

        if (is.null(dbRes)) {
            return(TRUE)
        }

        list_args <- dbRes$settings

        # check param
        # Do not appear in paramlist for custom ?
        if (!is.null(paramList$domain)) {
            if (!(setequal(paramList$domain, list_args$domain)))
                return(TRUE)
        }

        list_args_new <- paramList

        if (list_args_new$pvalueCutoff != list_args$pvalueCutoff)
            return(TRUE)

        switch(database,
               "GO" = {
                   if (list_args_new$OrgDb    != list_args$OrgDb)
                       return(TRUE)
                   if (list_args_new$keyType  != list_args$keyType)
                       return(TRUE)
               },
               "KEGG" = {
                   if (list_args_new$keyType  != list_args$keyType)
                       return(TRUE)
                   if (list_args_new$organism != list_args$organism)
                       return(TRUE)
               },
               "custom" = {

                   # if one analysis was launched with domain and the other without
                   # need to be able to launch it once again
                   if (!setequal(colnames(paramList$annotation),
                                 colnames(list_args$annotation)))
                       return(TRUE)

                   # Takes time
                   if (!setequal(paramList$annotation, list_args$annotation))
                       return(TRUE)
               }
        )

        return(FALSE)
    }



# ----- explain content ----
# functions for popify/bsshiny content


# ---- UI param functions ----
#' @noRd
#' @keywords internal
.annotSettings <- function(session, rea.values, dataset, listSource) {
    ns <- session$ns

    renderUI({
        validate(
            need(
                rea.values[[dataset]]$diffValid != FALSE,
                "Please run a differential analysis and validate the results before using this panel."
            )
        )

        ### lists
        ListNames <- switch(
            listSource,
            "DiffExp" = {
                ListNames.diff <- rea.values[[dataset]]$DiffValidContrast$contrastName
                names(ListNames.diff) <-
                    paste0("[",rea.values[[dataset]]$DiffValidContrast$tag,"] ",
                           ListNames.diff)
                ListNames.diff
            },
            "CoExp"   = rea.values[[dataset]]$CoExpClusterNames
        )

        tagList(column(
            width = 4,
            pickerInput(
                inputId  = ns("listToAnnot"),
                label    = "Selected lists:",
                choices  = ListNames)
        ),
        column(
            width = 6,
            br(),
            actionButton(
                inputId = ns("run"),
                label = "Run over-representation analysis",
                class = "butt"
            )
        ))
    })
}


#' @noRd
#' @keywords internal
.annotParamCustom <- function(session, rea.values, dataset) {
    ns <- session$ns

    renderUI({
        if (rea.values[[dataset]]$diffValid == FALSE)
            return()

        list(
            fileInput(
                inputId = ns("annotationFileCPR"),
                label = popify(
                    actionLink(inputId = "customAnnotFilepop",
                               label = "Annotation file"),
                    title = "Annotation File format",
                    content = paste0(
                        "<p>Required format:</p>",
                        "a tsv or csv file with at least two columns, ",
                        "one for omic features names and the other for",
                        " their associated terms."
                    ),
                    trigger = "click",
                    placement = "top"
                ),
                accept  = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            ),

            uiOutput(ns("selectColumnsCustom")),

            br(),
            ## alpha threshold
            numericInput(
                inputId = ns("pValue_custom"),
                label = "Adjusted p-value cutoff",
                value = 0.01 ,
                min = 0,
                max = 1,
                step = 0.01
            )
        )
    })
}

#' @noRd
#' @keywords internal
.annotFileColumns <- function(session, local.rea.values,
                              dataset) {
    ns <- session$ns
    dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
    varLabel <- .omicsDic(dataSE)$variableName
    varLabel <- paste0(toupper(substr(varLabel, 1, 1)),
                       substr(varLabel, 2, nchar(varLabel)))

    renderUI({
        annotation <- fread(
            file = local.rea.values$dataPathAnnot,
            sep = "\t",
            header = TRUE
        )

        column(
            width = 12,
            #br(),
            #hr(style = "border-top: 1px solid #000000;"),
            popify(
                actionLink(
                    inputId = ns(paste0("colSelAnnot")),
                    label = "Help",
                    icon = icon("question-circle")
                ),
                title = "Custom annotation",
                content = paste0("<p>Chose the appropriate columns names. ",
                                 "They are automatically derived from your file. ",
                                 "Items with a star are mandatory. </p>",
                                 "<p>- Term Name is used in case of a term name that is different from its id (eg GO:0009409 is the id, response to cold is the term name).</p>",
                                 "- Domain is used to differentiate databases or ontologies in the same file (eg: BP, CC and MF would be indicated in the domain column)."
                ),
                trigger = "click",
                placement = "top"
            ),

            # Select the right columns for the analysis
            pickerInput(
                inputId = ns("col_geneName"),
                label = paste0(varLabel, " ID: *"),
                choices = c("", colnames(annotation)),
                selected = ""
            ),
            pickerInput(
                inputId = ns("col_termID"),
                label = "Terms IDs: *",
                choices = c("", colnames(annotation)),
                selected = ""
            ),
            pickerInput(
                inputId = ns("col_termName"),
                label = "Term Names:",
                choices = c("", colnames(annotation)),
                selected = ""
            ),
            pickerInput(
                inputId = ns("col_domain"),
                label = "Domain:",
                choices = c("", colnames(annotation)),
                selected = ""
            ),
            hr(style = "border-top: 1px solid #000000;"),
        )
    })
}

#' @noRd
#' @keywords internal
.annotExFileColumns <- function(session, local.rea.values,
                                dataset) {
    ns <- session$ns
    dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
    varLabel <- .omicsDic(dataSE)$variableName
    varLabel <- paste0(toupper(substr(varLabel, 1, 1)),
                       substr(varLabel, 2, nchar(varLabel)))

    renderUI({
        column(
            width = 12,
            #br(),
            #hr(style = "border-top: 1px solid #000000;"),

            # Select the right columns for the analysis
            pickerInput(
                inputId = ns("col_geneName"),
                label = paste0(varLabel, " ID: *"),
                choices = "Gene stable ID",
                selected = "Gene stable ID"
            ),
            pickerInput(
                inputId = ns("col_termID"),
                label = "Terms IDs: *",
                choices = "GO term accession",
                selected = "GO term accession"
            ),
            pickerInput(
                inputId = ns("col_termName"),
                label = "Term Names:",
                choices = "GO term name",
                selected = "GO term name"
            ),
            pickerInput(
                inputId = ns("col_domain"),
                label = "Domain:",
                choices = "GO domain",
                selected = "GO domain"
            ),
            hr(style = "border-top: 1px solid #000000;"),
        )
    })
}

#' @noRd
#' @keywords internal
.annotParamGO <- function(session, rea.values, dataset) {
    ns <- session$ns

    renderUI({
        if (rea.values[[dataset]]$diffValid == FALSE)
            return()

        list(
            # select ontology
            pickerInput(
                inputId = ns("ont.select"),
                label = "Select Domain:",
                choices = c("ALL", "BP", "MF", "CC"),
                selected = "ALL",
                options = list(
                    `actions-box` = TRUE,
                    size = 10,
                    `selected-text-format` = "count > 3"
                )
            ),

            # select database
            pickerInput(
                inputId = ns("db.select"),
                label = span(
                    "Select Data Base:",
                    popify(
                        actionLink(ns("dbSel"), icon("question-circle")),
                        title = "",
                        content = paste0(
                            "The databases are automatically retrieved from your R library, searching for all package in the for of org.*db.",
                            " If none is found, please install the corresponding database using Bioconductor install manager."
                        ),
                        trigger = "click",
                        placement = "top"
                    )
                ),
                choices = c("", dir(.libPaths())[grep("org.*db", dir(.libPaths()))]),
                selected = "",
                options = list(
                    `actions-box` = TRUE,
                    size = 10,
                    `selected-text-format` = "count > 3"
                )
            ),

            # select keytype
            pickerInput(
                inputId = ns("keytype.go"),
                label = "Select id type:",
                choices = c(
                    "",
                    "ENTREZID",
                    "SYMBOL",
                    "TAIR",
                    "ENSEMBL",
                    "PMID",
                    "REFSEQ",
                    "GENENAME"
                ),
                selected = "",
                options = list(
                    `actions-box` = TRUE,
                    size = 10,
                    `selected-text-format` = "count > 3"
                )
            ),
            #),

            ## alpha threshold
            numericInput(
                inputId = ns("pValue_GO"),
                label = "Adjusted p-value",
                value = 0.01 ,
                min = 0,
                max = 1,
                step = 0.01
            )
        )

        #}
    })
}

#' @noRd
#' @keywords internal
.annotParamKEGG <- function(session, rea.values, dataset) {
    ns <- session$ns
    renderUI({
        if (rea.values[[dataset]]$diffValid == FALSE)
            return()

        list(
            pickerInput(
                inputId = ns("keyTypeKEGG"),
                label = "Select id type:",
                choices = c("", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
                selected = "",
                options = list(
                    `actions-box` = TRUE,
                    size = 10,
                    `selected-text-format` = "count > 3"
                )
            ),

            # type organism
            textInput(
                inputId = ns("organism"),
                label = popify(
                    actionLink(inputId = ns("colSelAnnot"),
                               label = "Organism code:"),
                    title = "KEGG Organism code",
                    content = paste0(
                        "Three or more letters indicating the species. ",
                        "For example: ",
                        "<ul><li>hsa for <i>Homo sapiens</i>,</li> ",
                        "<li>ath for <i>Arabidopsis thaliana</i>.",
                        "</li></ul>"
                    )
                ),
                value = "",
                width = NULL,
                placeholder = NULL
            ),

            ## alpha threshold
            numericInput(
                inputId = ns("pValueKEGG"),
                label = "Adjusted p-value",
                value = 0.01 ,
                min = 0,
                max = 1,
                step = 0.01
            )


        )
    })
}

# ---- main output functions -----
#' @noRd
#' @keywords internal
.outAnnotResults <- function(session,
                             input,
                             input2,
                             rea.values,
                             listSource,
                             dataset,
                             database) {
    ns <- session$ns

    # UI for all results
    renderUI({

        # Results:
        if (rea.values[[dataset]][[listSource]] == FALSE ||
            is.null(sumORA(session$userData$FlomicsMultiAssay[[dataset]],
                           from = listSource, database = database))) {
            return()
        }

        dataSE <- session$userData$FlomicsMultiAssay[[dataset]]

        varLabel0 <- .omicsDic(dataSE)$variableName

        # Explanation message for result table
        expMessage <-
            paste0(
                "This table present the results of the enrichment for the ",
                varLabel0,
                " list.",
                " <b>Description</b> gives the full name of  the term corresponding to the <b>ID</b> (if available).",
                " <b>GeneRatio</b> gives the number of ",
                varLabel0,
                " in the list that were also found in the term list.",
                " <b>BgRatio</b> stands for background ratio and gives the number of ",
                varLabel0,
                " of the term that were in the universe (number of annotated ",
                varLabel0,
                ").",
                " Pvalues, qvalues and adjusted pvalues were rounded to 3 digits."
            )

        # for KEGG database, 2 more columns:
        if (database == "KEGG") {
            expMessage <- paste0(
                expMessage,
                " <b>Category</b> and <b>subcategory</b> columns",
                " are given only for KEGG results. They are linked",
                " to this ontology."
            )
        }

        #foreach genes list selected (contrast)
        eres1 <-
            getEnrichRes(dataSE, from = listSource, database = database)


        lapply(names(eres1), function(listname) {
            eres <- getEnrichRes(
                dataSE,
                from = listSource,
                database = database,
                featureListName = listname
            )

            boxtitle <-
                switch (
                    listSource,
                    "DiffExp" = {
                        tag <-
                            rea.values[[dataset]]$DiffValidContrast[which(contrastName == listname),]$tag
                        paste0("[",tag, "] ", listname)
                    },
                    listname
                )

            sORA <- sumORA(dataSE, from = listSource,
                           database = database, featureListName = listname)

            cond <- as.numeric(unlist(sORA[-1], recursive = TRUE)) # warnings
            cond[is.na(cond)] <- 0
            if (any(cond > 0)) {

                eres <- Filter(Negate(is.null), eres)

                toKeep <- unlist(lapply(eres, FUN = function(enrichRes){
                    selectterms <- enrichRes@result$p.adjust < enrichRes@pvalueCutoff
                    return(sum(selectterms))
                }))
                choices <- names(eres[toKeep > 0])

                # ---- TabPanel plots ----
                tabPanel.list <- list(
                    tabPanel(
                        "DotPlot",
                        br(),
                        .outDotPlot(
                            input,
                            dataSE,
                            listname,
                            listSource,
                            "dotplot",
                            database
                        )
                    ),
                    tabPanel(
                        "Heatplot",
                        br(),
                        .outHeatPlot(
                            input,
                            dataSE,
                            listname,
                            listSource,
                            "heatplot",
                            database
                        )
                    ),
                    tabPanel(
                        "Cnetplot",
                        br(),
                        .outCnetPlot(
                            session,
                            input,
                            dataSE,
                            listname,
                            listSource,
                            "cnetplot",
                            database
                        )
                    ),
                    # ---- Tab Panel : Results table : ----
                    tabPanel(
                        "Result Table",
                        br(),
                        verticalLayout(
                            tags$style(
                                ".explain-p {
                                  color: Gray;
                                  text-justify: inter-word;
                                  font-style: italic;
                                }"
                            ),
                            div(class = "explain-p",
                                HTML(expMessage)),
                            hr(),

                            renderDataTable({
                                dataPlot <- getEnrichRes(
                                    object = dataSE,
                                    from = listSource,
                                    featureListName = listname,
                                    database = database,
                                    domain = input[[paste0(listname, "-domain")]]
                                )

                                pvalue <-
                                    getEnrichSettings(dataSE,
                                                      from = listSource,
                                                      database = database)$pvalueCutoff

                                datPlot <- dataPlot@result[dataPlot@result$p.adjust <  pvalue, ]
                                datPlot$pvalue <-
                                    round(datPlot$pvalue, 3)
                                datPlot$p.adjust <-
                                    round(datPlot$p.adjust, 3)
                                datPlot$qvalue <-
                                    round(datPlot$qvalue, 3)
                                datPlot$geneID <- unlist(lapply(
                                    datPlot$geneID,
                                    FUN = function(longString) {
                                        return(gsub("/", ", ", longString))
                                    }
                                ))

                                datatable(
                                    datPlot,
                                    rownames = FALSE,
                                    options = list(
                                        pageLength = 5,
                                        lengthMenu = c(5, 10, 15, 20),
                                        scrollX = TRUE
                                    )
                                )


                            }) # renderdatatable
                        )
                    ) # vertical layout

                )# TabsetPanel

                # ---- Tab Panel : only for KEGG, pathview : ----
                if (database == "KEGG") {
                    data <-   getEnrichRes(
                        object = dataSE,
                        featureListName = listname,
                        from = listSource,
                        database = "KEGG"
                    )[["no-domain"]]@result

                    pvalue <- getEnrichSettings(dataSE,
                                                from = listSource,
                                                database = database)$pvalueCutoff
                    mapChoices <-
                        sort(data$ID[data$p.adjust < pvalue])

                    tabPanel.list <- c(tabPanel.list,
                                       list(
                                           .outPathview(
                                               session,
                                               input,
                                               input2,
                                               data,
                                               dataSE,
                                               listSource,
                                               listname,
                                               mapChoices
                                           )
                                       ))

                } # if database KEGG

                # display results
                fluidRow(
                    box(
                        width = 12,
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        collapsed = TRUE,
                        status = "success",
                        title = boxtitle,

                        fluidRow(
                            column(
                                3,
                                # choice of domain
                                radioButtons(
                                    ns(paste0(listname, "-domain")),
                                    label = "Domain",
                                    choices = choices,
                                    selected = choices[1]
                                )
                            ),
                            column(2,
                                   renderUI({
                                       # number of terms
                                       dataPlot <- getEnrichRes(
                                           object = dataSE,
                                           from = listSource,
                                           featureListName = listname,
                                           database = database,
                                           domain = input[[paste0(listname, "-domain")]]
                                       )
                                       pvalue <-
                                           getEnrichSettings(dataSE,
                                                             from = listSource,
                                                             database = database)$pvalueCutoff
                                       datPlot <-
                                           dataPlot@result[dataPlot@result$p.adjust < pvalue, ]
                                       max_terms <- nrow(datPlot)

                                       numericInput(
                                           ns(paste0(
                                               listname, "-top.over"
                                           )),
                                           label = paste0("Top terms: (max: ",
                                                          max_terms, ")"),
                                           value = min(15, max_terms) ,
                                           min = 1,
                                           max = max_terms,
                                           step = 1
                                       )
                                   })),

                            column(
                                2,
                                # search term or gene
                                textInput(ns(
                                    paste0(listname, "-grep")
                                ),
                                label = "Search Expression")
                            ),
                            column(
                                2,
                                # rendering option
                                radioButtons(
                                    ns(paste0(listname, "-render")),
                                    label = "Rendering",
                                    choices = c("static", "interactive"),
                                    selected = "static"
                                )
                            ),
                            column(1,),
                            .popoverHelp(label = "")
                        ),
                        fluidRow(# display all tabsets
                            column(
                                width = 12,
                                do.call(what = tabsetPanel, args = tabPanel.list)
                            ))

                    ) # box
                ) # fluidrow
            } # if any result
            else return()
        })# lapply contrast names
    })
}

#' @noRd
#' @keywords internal
.popoverHelp <- function(label = "") {
    renderUI({
        column(1,
               br(),
               tags$style(
                   HTML(
                       "#classPop .popover {
              text-align: left;
              border-color: black;
              background-color: #d9edf7;
              width: 600px;
              max-width: 700px;
              color: black;
              opacity: 1;}
              "
                   )
               ),
               div(
                   id = "classPop",
                   span(
                       label,
                       popify(
                           h4(icon("question-circle"), "Help"),
                           title = "Parameters for enrichment results plots",
                           content = paste0(
                               "For each domain indicated (typically BP, CC and MF for GO enrichment), three graphs and the results ",
                               "table are displayed (all pvalues are rounded to 3 digits).",
                               " For more information on the plots, please check the clusterprofiler package vignette. ",
                               "All entries apply to the three graphs displayed. ",
                               "<p> You can choose the number of terms to be displayed, sorted by ascending order of ",
                               "adjusted pvalue.</p>",
                               "<p> The search expression allows you to display only the first terms containing the regular expression of interest. ",
                               "You can use regular expression patterns such as:",
                               "<ul><li> <b>||</b> for <i>or</i> statement; </li>",
                               "<li> <b>&</b> for <i>and</i> statement;</li>",
                               "</ul>"
                           ),
                           trigger = "click",
                           placement = "top"
                       ),
                       style = "color: #337ab7; border-color: #337ab7;"
                   )
               ))
    })

}

# ---- plots ----
#' @noRd
#' @keywords internal
.outHeatPlot <- function(input,
                         dataSE,
                         listname,
                         listSource,
                         plotType,
                         database) {
    varLabel0 <- .omicsDic(dataSE)$variableName
    renderUI({
        outplot <- tryCatch(
            plotClusterProfiler(
                dataSE,
                featureListName = listname,
                plotType = plotType,
                database = database,
                domain = input[[paste0(listname, "-domain")]],
                showCategory = input[[paste0(listname, "-top.over")]],
                searchExpr = input[[paste0(listname, "-grep")]]
            ),
            error = function(e) e,
            warning = function(w) w
        )

        plotExplain <- switch(
            listSource,
            "DiffExp" =  {
                paste0(
                    "This graph shows the top <b>",
                    input[[paste0(listname, "-top.over")]],
                    "</b> terms, arranged in alphabetical order. ",
                    " The tiles of <b>",
                    varLabel0,
                    "</b> present in the term are colored according to the log2 FC of the differential analysis of <b>",
                    listname,
                    "</b> . White tiles indicate the ",
                    varLabel0,
                    " is not present in the term."
                )
            },
            "CoExp" =
                paste0(
                    "This graph shows the top <b>",
                    input[[paste0(listname, "-top.over")]],
                    "</b> terms, arranged in alphabetical order. ",
                    "The color of the tile indicates the absence (white) or presence (black) of a ",
                    varLabel0,
                    " in the term."
                )
        )

        if (input[[paste0(listname, "-domain")]] == "no-domain") {
            subTitlePlot <- paste0("Database: ",
                                   database,
                                   " - ",
                                   input[[paste0(listname, "-top.over")]],
                                   " top enriched terms")
        } else {
            subTitlePlot <- paste0("Database: ",
                                   database,
                                   " ",
                                   input[[paste0(listname, "-domain")]],
                                   " - ",
                                   input[[paste0(listname, "-top.over")]],
                                   " top enriched terms")
        }


        if (is(outplot, "gg")) {

            switch(
                input[[paste0(listname, "-render")]],
                "static" = {

                    column(
                        width = 12,
                        tags$style(
                            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                        ),
                        div(class = "explain-p", HTML(plotExplain)),
                        hr(),
                        renderPlot({
                            outplot + labs(title = listname,
                                           subtitle = subTitlePlot)
                        })
                    )},
                "interactive" = {
                    column(
                        width = 12,
                        tags$style(
                            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                        ),
                        div(class = "explain-p", HTML(plotExplain)),
                        hr(),
                        renderPlotly({
                            ggplotly(outplot + labs(title = listname,
                                                    subtitle = subTitlePlot))
                        })
                    )
                })
        }else {
            renderText({outplot$message })
        }
    })
}

#' @noRd
#' @keywords internal
.outDotPlot <-
    function(input,
             dataSE,
             listname,
             listSource,
             plotType,
             database) {
        varLabel0 <- .omicsDic(dataSE)$variableName

        renderUI({
            outdot <- tryCatch(
                plotClusterProfiler(
                    dataSE,
                    featureListName = listname,
                    plotType = plotType,
                    database = database,
                    domain = input[[paste0(listname, "-domain")]],
                    showCategory = input[[paste0(listname, "-top.over")]],
                    searchExpr = input[[paste0(listname, "-grep")]]
                ),
                error = function(e) e,
                warning = function(w) w
            )

            plotExplain <- switch(
                listSource,
                "DiffExp" =  {
                    paste0(
                        "This graph shows the top <b>",
                        input[[paste0(listname, "-top.over")]],
                        "</b> terms, ordered by GeneRatio (number of ",
                        varLabel0,
                        " found differentially expressed, also called Count, relative to the total number of ",
                        varLabel0,
                        " for that term in the database). ",
                        "The dot color is linked to the adjusted pvalue.",
                        " Dot size is determined by the number of differentially expressed ",
                        varLabel0,
                        " in <b>",
                        listname,
                        "</b> and present in the term (Count)."
                    )
                },
                "CoExp" =
                    paste0(
                        "This graph shows the top <b>",
                        input[[paste0(listname, "-top.over")]],
                        "</b> terms, , ordered by GeneRatio (number of ",
                        varLabel0,
                        " found in the cluster , also called Count, relative to the total number of ",
                        varLabel0,
                        " for that term in the database). ",
                        "The dot color is linked to the adjusted pvalue.",
                        " Dot size is determined by the number of differentially expressed ",
                        varLabel0,
                        " in <b>",
                        listname,
                        "</b> and present in the term (Count)."
                    )
            )

            if (input[[paste0(listname, "-domain")]] == "no-domain") {
                subTitlePlot <- paste0("Database: ",
                                       database,
                                       " - ",
                                       input[[paste0(listname, "-top.over")]],
                                       " top enriched terms")
            } else {
                subTitlePlot <- paste0("Database: ",
                                       database,
                                       " ",
                                       input[[paste0(listname, "-domain")]],
                                       " - ",
                                       input[[paste0(listname, "-top.over")]],
                                       " top enriched terms")
            }

            if (is(outdot, "gg")) {
                switch(input[[paste0(listname, "-render")]],
                       "static" = {
                           column(
                               width = 12,
                               tags$style(
                                   ".explain-p {
                                    color: Gray;
                                    text-justify: inter-word;
                                    font-style: italic;
                                    }"
                               ),
                               div(class = "explain-p", HTML(plotExplain)),
                               hr(),

                               renderPlot({
                                   outdot +
                                       labs(title = listname,
                                            subtitle = subTitlePlot)
                               }))
                       },
                       "interactive" = {
                           column(
                               width = 12,
                               tags$style(
                                   ".explain-p {
                                    color: Gray;
                                    text-justify: inter-word;
                                    font-style: italic;
                                    }"
                               ),
                               div(class = "explain-p", HTML(plotExplain)),
                               hr(),
                               renderPlotly({
                                   ggplotly(outdot +
                                                labs(title = listname,
                                                     subtitle = subTitlePlot))
                               })
                           )
                       })

            }else
                renderText({
                    outdot$message
                })
        })
    }

#' @noRd
#' @keywords internal
.outCnetPlot <- function(session,
                         input,
                         dataSE,
                         listname,
                         listSource,
                         plotType,
                         database) {
    ns <- session$ns

    varLabel0 <- .omicsDic(dataSE)$variableName
    varLabel <- paste0(toupper(substr(varLabel0, 1, 1)),
                       substr(varLabel0, 2, nchar(varLabel0)))

    verticalLayout(renderUI({
        nodeLabelArg <- "none"
        if (input[[paste0(listname, "-genesLabels_cnet")]] &&
            input[[paste0(listname, "-termsLabels_cnet")]]) {
            nodeLabelArg <- "all"
        } else if (input[[paste0(listname, "-genesLabels_cnet")]] &&
                   !input[[paste0(listname, "-termsLabels_cnet")]]) {
            nodeLabelArg <- "gene"
        } else if (input[[paste0(listname, "-termsLabels_cnet")]] &&
                   !input[[paste0(listname, "-genesLabels_cnet")]]) {
            nodeLabelArg <- "category"
        }

        outcnet <- tryCatch(
            plotClusterProfiler(
                dataSE,
                featureListName = listname,
                # from = listSource,
                plotType = plotType,
                database = database,
                domain = input[[paste0(listname, "-domain")]],
                showCategory = input[[paste0(listname, "-top.over")]],
                searchExpr = input[[paste0(listname, "-grep")]],
                nodeLabel = nodeLabelArg
            ), warning = function(w) w,
            error = function(e) e
        )

        plotExplain <- switch(
            listSource,
            "DiffExp" =  {
                paste0(
                    "This graph represents the top <b>",
                    input[[paste0(listname, "-top.over")]],
                    "</b> terms (beige nodes), linked to their associated ",
                    varLabel0,
                    " (blue/red nodes).",
                    varLabel0,
                    " color of is determined by the log2FC derived from the differential analysis of <b>",
                    listname,
                    "</b> , with red indicating an up-regulated ",
                    varLabel0,
                    " and blue a down-regulated one.",
                    " The size of the term node correlates with the number of ",
                    varLabel0,
                    " to which they are linked."
                )
            },
            "CoExp" =
                paste0(
                    "This graph represents the top <b>",
                    input[[paste0(listname, "-top.over")]],
                    "</b> terms (beige nodes), linked to their associated ",
                    varLabel0,
                    ". ",
                    " The size of the term node correlates with the number of ",
                    varLabel0,
                    " to which they are linked."
                )
        )

        if (input[[paste0(listname, "-domain")]] == "no-domain") {
            subTitlePlot <- paste0("Database: ",
                                   database,
                                   " - ",
                                   input[[paste0(listname, "-top.over")]],
                                   " top enriched terms")
        } else {
            subTitlePlot <- paste0("Database: ",
                                   database,
                                   " ",
                                   input[[paste0(listname, "-domain")]],
                                   " - ",
                                   input[[paste0(listname, "-top.over")]],
                                   " top enriched terms")
        }

        if (is(outcnet, "gg")) {
            switch(
                input[[paste0(listname, "-render")]],
                "static" = {
                    column(
                        width = 12,
                        tags$style(
                            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                        ),
                        div(class = "explain-p", HTML(plotExplain)),
                        hr(),

                        renderPlot({
                            outcnet +
                                labs(title = listname,
                                     subtitle = subTitlePlot)

                        })
                    )
                },
                "interactive" = {
                    column(
                        width = 12,
                        tags$style(
                            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                        ),
                        div(class = "explain-p", HTML(
                            paste0(plotExplain, "
                                   Interactive version of this plot is not available yet."))),
                        hr(),
                        renderPlot({
                            outcnet +
                                labs(title = listname,
                                     subtitle = subTitlePlot)

                        })
                        # renderPlotly({
                        #     ggplotly(  outcnet +
                        #                    labs(title = listname,
                        #                         subtitle = subTitlePlot))
                        # })
                    )
                }
            )
        } else
            renderText({
                outcnet$message
            })

    }),
    fluidRow(column(
        2, checkboxInput(
            inputId = ns(paste0(listname,
                                "-genesLabels_cnet")),
            label = paste0(varLabel, " Labels"),
            value = FALSE
        )
    ),
    column(
        2,
        checkboxInput(
            inputId = ns(paste0(listname,
                                "-termsLabels_cnet")),
            label = "Terms Labels",
            value = TRUE
        )
    )))
}

#' @noRd
#' @keywords internal
.outPathview <- function(session,
                         input,
                         input2,
                         data,
                         dataSE,
                         listSource,
                         listname,
                         mapChoices) {
    ns <- session$ns

    varLabel0 <- .omicsDic(dataSE)$variableName
    pathviewExplain <-
        paste0(
            "<p> The link generated in this panel can be copy-pasted to access an interactive KEGG MAP. ",
            " On this MAP, you will find the ",
            varLabel0,
            "s of the list (contrast or cluster) highlighted."
        )

    tabPanel("KEGG Pathways results",
             br(),
             fluidRow(
                 column(
                     3,
                     selectInput(
                         inputId = ns(paste0(listname, "-MAP.sel")),
                         label = "Select map:",
                         choices = mapChoices,
                         multiple = FALSE,
                         selectize = FALSE,
                         size = 5,
                         selected = mapChoices[1]
                     ),
                 ),
                 column(
                     9,
                     h5("Link to interactive map online:"),
                     renderPrint({
                         link_to_map <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
                                               input[[paste0(listname, "-MAP.sel")]],
                                               "/",
                                               data[input[[paste0(listname, "-MAP.sel")]], "geneID"])
                         cat(paste0(link_to_map))
                     }),

                     tags$style(
                         ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                     ),
                     div(class = "explain-p", HTML(pathviewExplain)),
                 )
             ),
    )
}

# ---- Summary ----
#' @noRd
#' @keywords internal
.outEnrichSummary <- function(session,
                              rea.values,
                              input,
                              listSource,
                              title,
                              dataset,
                              database) {
    ns <- session$ns

    renderUI({

        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]

        if (rea.values[[dataset]][[listSource]] == FALSE ||
            is.null(sumORA(dataset.SE, listSource, database))) {
            return()
        }

        sORA <- sumORA(dataset.SE, from = listSource, database = database)
        errorMessages <-
            getAnalysis(dataset.SE,
                        name = paste0(listSource, "EnrichAnal"),
                        subName = database)$errors

        if (is.null(sORA) && is.null(errorMessages))
            return()

        options(DT.TOJSON_ARGS = list(na = 'string'))

        text.help <- c(paste0("0: no term was found significantly enriched for this ",
                              .omicsDic(dataset.SE)$variableName), " list")

        if(length(errorMessages) != 0){
            if(length(unique(unlist(errorMessages))) == 1){
                text.help <-
                    c(text.help, paste0("NA: ", unique(unlist(errorMessages))))

            }else{
                text.help <-
                    c(text.help, "NA: see the 'error messages tab'")
            }
        }

        # summary table of results
        tabPanel.list <- list(
            tabPanel(
                title = "Result table",
                br(),
                renderDataTable({
                    datatable(sORA, rownames = FALSE,
                              options = list( pageLength = 6, dom = 'tip'))
                }),
                br(),
                renderUI({
                    tagList(
                        HTML(paste(
                            "<div style='border: 1px solid black; padding: 10px; background-color: #f0f0f0; font-style: italic;'>",
                            paste(text.help, collapse = "<br>"),
                            "</div>"
                        ))
                    )
                })
            ))

        # summary enrichment plot
        if(any(as.numeric(unlist(sORA[-1], recursive = TRUE)) > 0, na.rm = TRUE))
            tabPanel.list <-
            c(tabPanel.list,
              list(
                  tabPanel(
                      title = "Comparison of results (Heatmap)",
                      br(),
                      .outCompResults(
                          session,
                          rea.values,
                          input,
                          listSource,
                          dataset,
                          database
                      )
                  )))

        # error message table
        if(!is.null(errorMessages) &&
           length(unique(unlist(errorMessages))) > 1)
            tabPanel.list <-
            c(tabPanel.list,
              list(
                  tabPanel(
                      title = "Error message",
                      br(),

                      renderDataTable({

                          errorMessages.df <- as.data.frame(unlist(errorMessages))
                          names(errorMessages.df) <- "error message"

                          datatable(errorMessages.df, rownames = TRUE, escape = 1,
                                    options = list( pageLength = 6, dom = 'tip'))
                      })
                  )))

        if (is.null(sORA))  {
            fluidRow(
                box(
                    width = 12,
                    solidHeader = TRUE,
                    collapsible = FALSE,
                    collapsed = FALSE,
                    status = "warning",
                    title = "There is no result for this enrichment analysis!",

                    if(length(errorMessages) != 0){

                        renderUI({
                            tagList(
                                HTML(paste(
                                    "<div style='border: 1px solid black; padding: 10px; background-color: #f0f0f0; font-style: italic;'>",
                                    unique(unlist(errorMessages)),
                                    "</div>"
                                ))
                            )
                        })
                    }
                )
            )
        } else {

            fluidRow(
                box(
                    width = 12,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    status = "primary",
                    title = title,

                    do.call(what = tabsetPanel, args = tabPanel.list)
                )
            )
        }
    })
}

#' @description draw the heatmap for the enrichment summary
#' @noRd
#' @keywords internal
.outCompResults <- function(session,
                            rea.values,
                            input,
                            listSource,
                            dataset,
                            database) {

    ns <- session$ns
    dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
    varLabel0 <- .omicsDic(dataSE)$variableName

    renderUI({
        if (rea.values[[dataset]][[listSource]] == FALSE ||
            is.null(sumORA(dataSE, listSource, database))) {
            return()
        }

        # Possible domains of ontology:
        possDomain <-  unique(unlist(lapply(
            getAnalysis(dataSE,
                        name = paste0(listSource, "EnrichAnal"),
                        subName = database)$results$enrichResult,
            FUN = function(x)
                names(x)
        )))

        # display results
        verticalLayout(fluidRow(column(
            6,
            radioButtons(
                inputId = ns(paste0(database, "-compDomain")),
                label = "Domain",
                choices = possDomain,
                selected = possDomain[1],
                inline = TRUE
            )
        ),
        column(
            6,
            radioButtons(
                inputId = ns(paste0(database, "-compType")),
                label = "Matrix Type",
                choices = c("presence", "FC"),
                selected = "FC",
                inline = TRUE
            )
        ),),
        fluidRow(
            column(12,
                   renderUI({
                       clustmeth <- switch(input[[paste0(database, "-compType")]],
                                           "presence" = "complete",
                                           "complete")
                       distmeth <-
                           switch(input[[paste0(database, "-compType")]],
                                  "presence" = "binary",
                                  "euclidean")
                       additionalMessage <- switch(
                           input[[paste0(database, "-compType")]],
                           "presence" = "Red tiles indicates the term has
                                 been found enriched using the list, white tile
                                 is an absent term. ",
                           "FC" = "Tiles are colored according to the
                           fold change (FC) of the enrichment. White tiles
                           denote the non-significance for a term in a list."
                       )


                       plotExplain <- switch(
                           listSource,
                           "DiffExp" =  {
                               paste0(
                                   "This graph shows all the enriched terms found with <b> ",
                                   database,
                                   "</b> database with differentially expressed ",
                                   varLabel0,
                                   " lists.",
                                   " Clustering of terms was computed using <b> ",
                                   distmeth,
                                   "</b> and <b>",
                                   clustmeth,
                                   "</b>. ",
                                   additionalMessage
                               )
                           },
                           "CoExp" =
                               paste0(
                                   "This graph shows all the enriched terms found with <b>",
                                   database,
                                   "</b> database with the co-expression clusters ",
                                   varLabel0,
                                   " lists.",
                                   " Clustering of terms was computed using <b> ",
                                   distmeth,
                                   "</b> and <b>",
                                   clustmeth,
                                   "</b>. ",
                                   additionalMessage
                               )
                       )

                       outHeatmap <- tryCatch({
                           outHeatmap <-  plotEnrichComp(
                               dataSE,
                               from = listSource,
                               database = database,
                               domain = input[[paste0(database, "-compDomain")]],
                               matrixType = input[[paste0(database, "-compType")]]
                           )
                       },
                       warning = function(warn)
                           warn,
                       error = function(err)
                           err)

                       if (is(outHeatmap, "Heatmap")) {
                           column(
                               width = 12,
                               tags$style(
                                   ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                               ),
                               hr(),
                               div(class = "explain-p", HTML(plotExplain)),
                               hr(),
                               renderPlot({
                                   draw(
                                       outHeatmap,
                                       heatmap_legend_side = "top",
                                       padding = unit(5, "mm"),
                                       gap = unit(2, "mm")
                                   )
                               },
                               width = "auto",
                               height = min(400 + nrow(outHeatmap@matrix) * 50, 1000))
                           )
                       } else {
                           renderText({
                               outHeatmap$message
                           })
                       }
                   })) ))
    })
}
