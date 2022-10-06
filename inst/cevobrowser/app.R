library(cevomod)
library(cevoDatasets)
library(shiny)
library(shinydashboard)
library(tibble)


datasets <- readr::read_rds("/mnt/dane/projects_all/cancer_evolution/cancer_evolution/cevomod_analyses/data.Rds")


ui <- dashboardPage(
  dashboardHeader(title = "cevobrowser"),
  dashboardSidebar(
    radioButtons(
      "dataset_selection",
      label = "Dataset",
      choices = names(datasets),
      selected = names(datasets)[[1]]
    )
  ),
  dashboardBody(
    plotOutput("SFS")
  )
)


server <- function(input, output) {
  output$SFS <- renderPlot({
    dataset <- input$dataset_selection
    cd <- datasets[[dataset]]
    plot_SFS(cd)
  })
}


shinyApp(ui = ui, server = server)
