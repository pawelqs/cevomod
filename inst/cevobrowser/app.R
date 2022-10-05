library(cevomod)
library(cevoDatasets)
library(shiny)
library(shinydashboard)


datasets <- lst(
  AMLRO, OPUS_BRCA, OPUS_Larynx
)


ui <- dashboardPage(
  dashboardHeader(title = "cevobrowser"),
  dashboardSidebar(
    radioButtons(
      "dataset_selection",
      label = "Dataset",
      choices = c("OPUS_BRCA", "OPUS_Larynx", "AMLRO"),
      selected = "AMLRO"
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
