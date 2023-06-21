library(cevomod)
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(tibble)
library(ggplot2)

theme_set(theme_minimal())
datasets <- test_data_fitted

# Header ----
header <- dashboardHeader(title = stringr::str_c("cevoapp ", packageVersion("cevomod")))


# Sidebar ----
sidebar <- dashboardSidebar(
  radioButtons(
    "dataset_selection",
    label = "Dataset",
    choices = names(datasets),
    selected = default_dataset
  ),
  selectizeInput(
    "patients_selection",
    label = "Select patients to show",
    choices = unique(datasets[[default_dataset]]$patient_id),
    multiple = TRUE,
    options = list(create = TRUE)
  ),
  sidebarMenu(
    menuItem("Overview", tabName = "overview_tab", icon = icon("magnifying-glass-chart")),
    menuItem("SFS", tabName = "SFS_tab", icon = icon("chart-simple")),
    menuItem("CNV", tabName = "CNV_tab", icon = icon("square-poll-horizontal")),
    menuItem("Models", tabName = "models_tab", icon = icon("chart-line")),
    menuItem("Residuals", tabName = "residuals_tab", icon = icon("database"))
  )
)


# Body ----

## Overview tab ----
overview_tab <- tabItem(
  tabName = "overview_tab",
  fluidRow(
    valueBoxOutput("dataset_name"),
    valueBoxOutput("n_patients"),
    valueBoxOutput("n_mutations")
  ),
  fluidRow(
    tabBox(
      id = "overview_plots_tabset", height = "35vh",
      tabPanel(
        "Sequencing Depth",
        plotOutput("overview_plot_DP", height = "30vh"),
        height = "35vh",
        width = 6L,
      ),
      tabPanel(
        "Private and shared mutations",
        plotOutput("overview_plot_private_shared", height = "30vh"),
        height = "35vh",
        width = 6L
      ),
      width = 6L
    ),
    box(
      # plotOutput("overview_plot_private_shared", height = "30vh"),
      height = "35vh",
      width = 6L
    )
  )
)


## SFS tab ----
SFS_tab <- tabItem(
  tabName = "SFS_tab",
  fluidRow(
    column(
      width = 9L,
      box(
        plotOutput("SFS_plot", height = "80vh"),
        height = "95vh",
        width = NULL,
      )
    ),
    column(
      width = 3L,
      box(
        selectInput(
          "plot_type",
          label = "Plot",
          choices = c("SFS", "model")
        ),
        height = "80vh",
        width = NULL
      )
    )
  )
)


## CNV tab ----
CNV_tab <- tabItem(
  tabName = "CNV_tab",
  fluidRow(
    column(
      width = 9L,
      box(
        plotOutput("CNV_plot", height = "80vh"),
        height = "95vh",
        width = NULL,
      )
    ),
    column(
      width = 3L,
      box(
        uiOutput("cnv_plot_selector"),
        height = "95vh",
        width = NULL
      )
    )
  )
)


## Models tab ----
models_tab <- tabItem(
  tabName = "models_tab",
  fluidRow(
    column(
      width = 9L,
      tabBox(
        id = "modelplots_tabset", height = "80vh",
        tabPanel(
          "SFS",
          plotOutput("models_SFS_plot", height = "80vh")
        ),
        tabPanel(
          "M(f) ~ 1/f",
          plotOutput("models_Mf_1f_plot", height = "80vh")
        ),
        tabPanel(
          "Cumulative tails",
          plotOutput("models_cum_tails_plot", height = "80vh")
        ),
        width = NULL
      )
    ),
    column(
      width = 3L,
      box(
        checkboxGroupInput(
          "model_layers_checkbox",
          "Select model layers:",
          choices = c("Neutral model", "Binomial Fit", "Clones", "Full model", "Drivers"),
          selected = c("Neutral model", "Clones", "Full model")
        ),
        height = "95vh",
        width = NULL
      )
    )
  )
)


## Residuals tab ----
residuals_tab <- tabItem(
  tabName = "residuals_tab",
  fluidRow(
    column(
      width = 9L,
      tabBox(
        id = "residualsplot_tabset", height = "80vh",
        tabPanel(
          "Sampling rate",
          plotOutput("sampling_rate_plot", height = "80vh")
        ),
        tabPanel(
          "Neutral model residuals",
          plotOutput("neutral_model_resid_plot", height = "80vh")
        ),
        tabPanel(
          "Full model residuals",
          plotOutput("full_model_resid_plot", height = "80vh")
        ),
        width = NULL
      )
    ),
    column(
      width = 3L,
      box(
        height = "95vh",
        width = NULL
      )
    )
  )
)


# Compose UI ----
body <- dashboardBody(
  tabItems(
    SFS_tab,
    CNV_tab,
    models_tab,
    residuals_tab,
    overview_tab
  )
)


ui <- dashboardPage(header, sidebar, body)


# Server ----
server <- function(input, output) {
  rv <- reactiveValues(cd = datasets[[default_dataset]])

  observeEvent(input$dataset_selection, {
    rv$cd <- datasets[[input$dataset_selection]]
    updateSelectizeInput(
      getDefaultReactiveDomain(),
      "patients_selection",
      choices = unique(rv$cd$metadata$patient_id)
    )
  })


  observeEvent(input$patients_selection, {
    rv$cd <- if (length(input$patients_selection) > 0) {
      datasets[[input$dataset_selection]] |>
        filter(patient_id %in% input$patients_selection)
    } else {
      datasets[[input$dataset_selection]]
    }
  })

  output$overview_plot_DP <- renderPlot({
    plot_sequencing_depth(rv$cd) +
      theme(axis.text.x = element_text(angle = 45))
  })

  output$overview_plot_private_shared <- renderPlot({
    plot_private_shared_mutations(rv$cd) +
      scale_fill_paletteer_d("PNWColors::Bay") +
      theme(axis.text.x = element_text(angle = 45))
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
      value = str_c(n_patients, "/", n_samples),
      subtitle = "Number of patients/samples",
      icon = icon("user")
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

  output$models_SFS_plot <- renderPlot({
    neutral_tail <- "Neutral model" %in% input$model_layers_checkbox
    binomial_fit <- "Binomial Fit" %in% input$model_layers_checkbox
    subclones <- "Clones" %in% input$model_layers_checkbox
    final_fit <- "Full model" %in% input$model_layers_checkbox
    show_drivers <- "Drivers" %in% input$model_layers_checkbox

    layers <- list(
      if (show_drivers) layer_mutations(drivers = rv$cd$cancer),
      hide_legend()
    )
    plot_models(
      rv$cd,
      neutral_tail = neutral_tail, binomial_layer = binomial_fit,
      subclones = subclones, final_fit = final_fit
    ) +
      layers
  })

  output$models_Mf_1f_plot <- renderPlot({
    neutral_lm_fitted <- "Neutral model" %in% input$model_layers_checkbox
    layers <- list(
      if (neutral_lm_fitted) {
        layer_lm_fits(rv$cd)
      }
    )
    plot_Mf_1f(rv$cd, from = 0.05, to = 0.5, scale = FALSE, mapping = ggplot2::aes(color = sample)) +
      layers +
      ggplot2::facet_wrap(~patient_id, scales = "free_y") +
      ggplot2::scale_color_brewer(palette = "Dark2")
  })

  output$models_cum_tails_plot <- renderPlot({
    plot_cumulative_tails(rv$cd) +
      ggplot2::facet_wrap(~patient_id) +
      ggplot2::aes(color = sample) +
      ggplot2::coord_cartesian(xlim = c(0.01, 1)) +
      ggplot2::scale_color_brewer(palette = "Dark2")
  })

  output$sampling_rate_plot <- renderPlot({
    plot_sampling_rate(rv$cd) +
      ggplot2::facet_wrap(~sample_id, scales = "free_y") +
      hide_legend()
  })

  output$neutral_model_resid_plot <- renderPlot({
    plot_residuals_powerlaw_model(rv$cd) +
      ggplot2::facet_wrap(~sample_id, scales = "free_y") +
      hide_legend()
  })

  output$full_model_resid_plot <- renderPlot({
    plot_residuals_full_model(rv$cd) +
      ggplot2::facet_wrap(~sample_id, scales = "free_y") +
      hide_legend()
  })
}


shinyApp(ui = ui, server = server)
