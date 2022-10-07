library(cevomod)
library(cevoDatasets)
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(tibble)

datasets <- readr::read_rds("/mnt/dane/projects_all/cancer_evolution/cancer_evolution/cevomod_analyses/data.Rds")
default_dataset <- names(datasets)[[1]]

header <- dashboardHeader(title = stringr::str_c("cevobrowser ", packageVersion("cevomod")))

sidebar <- dashboardSidebar(
  radioButtons(
    "dataset_selection",
    label = "Dataset",
    choices = names(datasets),
    selected = default_dataset
  ),
  sidebarMenu(
    menuItem("SFS", tabName = "SFS_tab", icon = icon("chart-simple")),
    menuItem("CNV", tabName = "CNV_tab", icon = icon("square-poll-horizontal")),
    menuItem("Models", tabName = "models_tab", icon = icon("chart-line"))
  )
)


SFS_tab <- tabItem(
  tabName = "SFS_tab",
  fluidRow(
    column(
      width = 8L,
      box(
        plotOutput("SFS_plot", height = "80vh"),
        height = "95vh",
        width = NULL,
      )
    ),
    column(
      width = 4L,
      valueBoxOutput("dataset_name", width = NULL),
      box(
        selectInput(
          "plot_type",
          label = "Plot",
          choices = c("SFS", "model")
        ),
        selectizeInput(
          "patients_list",
          label = "Select patients to show",
          choices = unique(datasets[[default_dataset]]$patient_id),
          multiple = TRUE,
          options = list(create = TRUE)
        ),
        height = "80vh",
        width = NULL
      )
    )
  )
)


CNV_tab <- tabItem(
  tabName = "CNV_tab",
  fluidRow(
    column(
      width = 8L,
      box(
        plotOutput("CNV_plot", height = "80vh"),
        height = "95vh",
        width = NULL,
      )
    ),
    column(
      width = 4L,
      box(
        uiOutput("cnv_plot_selector"),
        # sliderInput()
        height = "95vh",
        width = NULL
      )
    )
  )
)


models_tab <- tabItem(
  tabName = "models_tab",
  fluidRow(
    column(
      width = 8L,
      box(
        plotOutput("models_plot", height = "80vh"),
        height = "95vh",
        width = NULL,
      )
    ),
    column(
      width = 4L,
      box(
        height = "95vh",
        width = NULL
      )
    )
  )
)

body <- dashboardBody(
  tabItems(
    SFS_tab,
    CNV_tab,
    models_tab
  )
)


ui <- dashboardPage(header, sidebar, body)


server <- function(input, output) {

  rv <- reactiveValues(cd = datasets[[default_dataset]])

  observeEvent(input$dataset_selection, {
    rv$cd <- datasets[[input$dataset_selection]]
    updateSelectizeInput(
      getDefaultReactiveDomain(),
      "patients_list",
      choices = unique(rv$cd$metadata$patient_id)
    )
  })

  output$dataset_name <- renderValueBox({
    valueBox(
      value = rv$cd$name,
      subtitle = "Dataset name",
      icon = icon("database")
    )
  })

  output$n_patients <- renderValueBox({
    n_patients <- dplyr::n_distinct(rv$cd$metadata$patient_id)
    n_samples <- nrow(rv$cd$metadata)
    valueBox(
      value = stringr::str_c(n_patients, "/", n_samples),
      subtitle = "Patients/Samples",
      icon = icon("users")
    )
  })

  # output$n_mutations <- renderValueBox({
  #   n_patients <- dplyr::n_distinct(rv$cd$metadata$patient_id)
  #   n_samples <- nrow(rv$cd$metadata)
  #   valueBox(
  #     value = stringr::str_c(n_patients, "/", n_samples),
  #     subtitle = "Patients/Samples",
  #     icon = icon("users")
  #   )
  # })

  output$SFS_plot <- renderPlot({
    plot_SFS(rv$cd, geom = "bar") +
      ggplot2::facet_wrap(~sample_id, scales = "free") +
      hide_legend()
  })

  output$cnv_plot_selector <- renderUI({
    selectInput(
      "cnv_meta_field",
      label = "Plot",
      choices = get_CNVs_var_names(rv$cd)
    )
  })

  output$CNV_plot <- renderPlot({
    plot_CNV_heatmap(rv$cd, meta_field = input$cnv_meta_field)
  })

  output$models_plot <- renderPlot({
    plot_models(rv$cd) +
      hide_legend()
  })
}


shinyApp(ui = ui, server = server)
