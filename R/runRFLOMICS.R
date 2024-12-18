### ============================================================================
### run interface
### ----------------------------------------------------------------------------
# D. Charif

#' @title Run RFLOMICS interface
#' @description running this function will open the shiny application.
#' Run the shiny application
#' @param ... More arguments to pass to shinyApp.
#' @return shinyApp
#' @importFrom shinyBS popify
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
# @import shiny
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' library(RFLOMICS)
#' ## Not run: runRFLOMICS()
#'
#'
runRFLOMICS <- function(...) {
  addResourcePath(prefix = 'www', system.file('RFLOMICSapp/www/', package = 'RFLOMICS'))
  shinyApp(ui = rflomicsUI, server = shinyServer(rflomicsServer), ...)
}
