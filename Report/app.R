library(shiny)
library(shinydashboard)
library(dplyr)
library(readr)
library(here)
library(stringr)
library(PatientProfiles)
library(DT)
library(shinycssloaders)
library(shinyWidgets)
library(gt)
library(scales)
library(kableExtra)
library(tidyr)
library(stringr)
library(ggplot2)
library(fresh)
library(plotly)
#library(IncidencePrevalence)
#library(DiagrammeR)
#library(DiagrammeRsvg)

#load("mergedData.Rdata")

source(here("Report","functions.R"))

# Load data
Lasso <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/lasso_selectedfeatures_public_100k.csv") |>
  select(! database)
data_to <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/tableOne_summary_public_100k.csv") |>
  as_tibble() |> 
  visOmopResults::splitAdditional() |>
  visOmopResults::splitGroup() |>
  select(cohort_name, strata_level, variable_name, estimate_name, estimate_value, table, window) |>
  rename(age_group = strata_level)
data_to_tidy <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/tableOne_tidy_public_100k.csv")
data_to_sum <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/tableOne_SMD_public_100k.csv")

data_lsc <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/largeScale_summary_public_100k.csv") |>
  as_tibble() |> 
  visOmopResults::splitAdditional() |>
  visOmopResults::splitGroup() |>
  select(cohort_name, strata_level, variable_name, variable_level, estimate_name, estimate_value, table_name, variable_level) |>
  rename(age_group = strata_level)
data_lsc_tidy <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/largeScale_tidy_public_100k.csv")
data_lsc_sum <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/largeScale_SMD_public_100k.csv")

vac_attrition <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/vac_attrition_public_100k.csv")
unvac_attrition <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/unvac_attrition_public_100k.csv")

matched_vac_attrition <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/matched_vac_attrition_public_100k.csv")
matched_unvac_attrition <- read.csv("/home/AD_NDORMS/asanchezparada/HPVProject/Study/Results/matched_unvac_attrition_public_100k.csv")


# theme -----------------------
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#0c0e0c" 
  ),
  adminlte_sidebar(
    # width = "400px",
    dark_bg = "#78B7C5",
    dark_hover_bg = "#3B9AB2",
    dark_color = "white"
  ), 
  adminlte_global(
    content_bg = "#eaebea" 
  ),
  adminlte_vars(
    border_color = "#112446",
    active_link_hover_bg = "#FFF",
    active_link_hover_color = "#112446",
    active_link_hover_border_color = "#112446",
    link_hover_border_color = "#112446"
  )
)

# ui ------------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = "HPV Project"),
  ## menu -------------------------------------
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        text = "Background",
        tabName = "background"
      ),
      menuItem(
        text = "Study cohorts",
        tabName = "cohorts",
        menuSubItem(
          text = "Cohort Definitions",
          tabName = "definitions"
        ),
        menuSubItem(
          text = "Lasso Selected Features",
          tabName = "lasso"
        )
      ),
      menuItem(
        text = "Patient-level Characterisation",
        tabName = "char",
        menuSubItem(
          text = "Table One",
          tabName = "table_one"
        ),
        menuSubItem(
          text = "Large scale characteriscs",
          tabName = "lsc"
        )
      )
    )
  ),
  
  ## body ----
  dashboardBody(
    use_theme(mytheme),
    tabItems(
      # background  -------------------------------
      tabItem(
        tabName = "background",
        a(
          img(
            src = "https://raw.githubusercontent.com/darwin-eu/DrugUtilisation/main/man/figures/DrugUtilisation.png", 
            align = "right", width = "200px"
          ),
          href = "https://github.com/darwin-eu/DrugUtilisation/",
          target = "_blank"
        ),
        h3("xxxxx"),
        h5("DrugUtilisation is an R package that contains functions to conduct Patient level and Population level Drug Utilisation Studies in the OMOP common data model."),
        h5("In this study the main functionalities of the package are shown in a 'demonstration' study."),
        h5("The study shows: ... simvarstatin."),
        h5("The study was conducted in IQVIA Pharmetrics for academics (R).")
      ),
      # cohort attrition ------
      tabItem(
        tabName = "definitions",
        h4("Attrition for the study cohorts"),
        selectInput(inputId = "cohort",
                    label = "Select Cohort:",
                    choices = c("Unvaccinated cohort" = 1, "Vaccinated Cohort" = 2, "Matched Unvaccinated Cohort" = 3, "Matched Vaccinated Cohort" = 4)
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Tidy table",
            DT::dataTableOutput("attrition_tidy") |> withSpinner()
          ),
          tabPanel(
            "Diagram",
            #grVizOutput("attrition_diagram", width = "400px", height = "100%") %>% withSpinner()
            DT::dataTableOutput("attrition_tidy2") |> withSpinner()
          )
        )
      ),
      # Lasso ----------------------------------
      tabItem(
        tabName = "lasso",
        h4("Lasso Selected Features"),
        pickerInput(
          inputId = "year",
          label = "Select supopulation year: ",
          choices = c(Lasso[["year"]] |> unique() |> sort(), "all"),
          selected = "all",
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
          inline = TRUE
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Table",
            DT::DTOutput("lasso_table") |> withSpinner()
          )
        )
      ),
      
      # characteristics ----------------------------
      tabItem(
        tabName = "table_one",
        h4("Table One Characterisation"),
        selectors(
          data = data_to,
          prefix = "to",
          columns = c("cohort_name", "age_group", "variable_name", "estimate_name", "window"),
          default = list(age_group = "overall", estimate_name = c("count", "percentage", "median", "q25", "q75"))
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            DT::dataTableOutput("table_one_raw") |> withSpinner()
          ),
          tabPanel(
            "Tidy data",
            DT::dataTableOutput("table_one_tidy") |> withSpinner()
          ),
          tabPanel(
            "Summarised data",
            DT::dataTableOutput("table_one_sum") |> withSpinner()
          ),
          tabPanel(
            "SMD plot",
            plotlyOutput("to_smd") |> withSpinner()
          )
        )
      ),
      # lsc -----------------------------------
      tabItem(
        tabName = "lsc",
        h4("Large scale characteristics"),
        selectors(
          data = data_lsc,
          prefix = "lsc",
          columns = c("cohort_name", "age_group", "variable_name", "estimate_name", "table_name", "variable_level"),
          default = list(age_group = "overall", estimate_name = c("count", "percentage", "median", "q25", "q75"))
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            DT::dataTableOutput("lsc_raw") |> withSpinner()
          ),
          tabPanel(
            "Tidy data",
            DT::dataTableOutput("lsc_tidy") |> withSpinner()
          ),
          tabPanel(
            "Summarised data",
            DT::dataTableOutput("lsc_sum") |> withSpinner()
          ),
          tabPanel(
            "SMD plot",
            plotlyOutput("lsc_smd") |> withSpinner()
          )
        )
      )
      # end ---
    )
  )
)

# Server ---------------------------------------------------------------------
server <- function(input, output, session) {
  # deinitions: cohort attritions -------------
  
  data_of_interest <- reactive({
    if (input$cohort == 1) {
      unvac_attrition
    } else if (input$cohort == 2) {
      vac_attrition
    } else if (input$cohort == 3) {
      matched_unvac_attrition
    } else {
      matched_vac_attrition
    }
  })
  
  output$attrition_tidy <- DT::renderDataTable({
    data_of_interest()
  })
  
  
  # Lasso -----------------------------------
  # Select desired data
  
  # Filter data
  filtered_data <- reactive({
    print("ha")
    if (input$year != "all") {
      filter(Lasso, year == input$year)
    }else{
      Lasso
    }
  })
  
  output$lasso_table <- DT::renderDataTable({
    DT::datatable(filtered_data(), rownames = FALSE)
  })
  
  
  # Table one ------------------------------------
  # Computations
  TableOne_raw <- reactive({
    data_to|>
      filter(cohort_name %in% input$to_cohort_name, estimate_name %in% input$to_estimate_name, variable_name %in% input$to_variable_name, age_group %in% input$to_age_group, window %in% input$to_window)
  })  
  
  TableOne_tidy <- reactive({
    data_to_tidy |> 
      filter(estimate_name %in% input$to_estimate_name, variable_name %in% input$to_variable_name, age_group %in% input$to_age_group, window %in% input$to_window) |>
      select(! setdiff(c("total_unvac_matched_cohort", "total_vac_matched_cohort", "unvac_cohort", "vac_cohort"), input$to_cohort_name))
  })
  
  TableOne_sum <- reactive({
    data_to_sum |>
      filter(variable_name %in% input$to_variable_name, age_group %in% input$to_age_group, window %in% input$to_window)
    #      select(! c("Ctotal_unvac_matched_cohort", "total_vac_matched_cohort", "unvac_cohort", "vac_cohort"))
  })
  
  # Outputs
  output$table_one_raw <- DT::renderDataTable({
    DT::datatable(TableOne_raw(), options = list(scrollX = TRUE))
  })
  
  output$table_one_tidy <- DT::renderDataTable({
    DT::datatable(TableOne_tidy(), options = list(scrollX = TRUE))
  })
  
  output$table_one_sum <- DT::renderDataTable({
    DT::datatable(TableOne_sum(), options = list(scrollX = TRUE))
  })
  
  output$to_smd <- renderPlotly({
    plot_ly(data_to_sum[, c("original_smd", "matched_smd", "variable_name")], x = ~original_smd, y = ~matched_smd, text = ~variable_name, type = "scatter", mode = "markers") %>%
      layout(title = "SMD's for Table One Characterisation",
             xaxis = list(title = "Before Matching"),
             yaxis = list(title = "After Matching"), 
             shapes = list(
               list(
                 type = "line",
                 x0 = 0,
                 x1 = 1,
                 xref = "container",
                 y0 = 0.1,
                 y1 = 0.1,
                 line = list(color = "black")
               ),
               list(
                 type = "line",
                 x0 = 0,
                 x1 = 0.5,
                 xref = "container",
                 y0 = 0,
                 y1 = 0.5,
                 line = list(color = "red")
               )
             )
      )
  })
  
  # slc ---------------------------------
  LargeScale_raw <- reactive({
    data_lsc |>
      filter(cohort_name %in% input$lsc_cohort_name, estimate_name %in% input$lsc_estimate_name, variable_name %in% input$lsc_variable_name, table_name %in% input$lsc_table_name, age_group %in% input$lsc_age_group, variable_level %in% input$lsc_variable_level)
  })
  
  LargeScale_tidy <- reactive({
    data_lsc_tidy |>
      filter(estimate_name %in% input$lsc_estimate_name, variable_name %in% input$lsc_variable_name, table_name %in% input$lsc_table_name, variable_level %in% input$lsc_variable_level) |> # NO AGE GROUP: age_group %in% input$lsc_age_group) |>
      select(! setdiff(c("total_unvac_matched_cohort", "total_vac_matched_cohort", "unvac_cohort", "vac_cohort"), input$lsc_cohort_name))  # This not working
  })
  
  LargeScale_sum <- reactive({
    data_lsc_sum |>
      filter(variable_name %in% input$lsc_variable_name, table_name %in% input$lsc_table_name, variable_level %in% input$lsc_variable_level)
    #     select(! c("total_unvac_matched_cohort", "total_vac_matched_cohort", "unvac_cohort", "vac_cohort"))
  })
  
  output$lsc_raw <- DT::renderDataTable({
    datatable(LargeScale_raw(), options = list(scrollX = TRUE))
  })
  
  output$lsc_tidy <- DT::renderDataTable({
    datatable(LargeScale_tidy(), options = list(scrollX = TRUE))
  })
  
  output$lsc_sum <- DT:: renderDataTable({
    datatable(LargeScale_sum(), options = list(scrollX = TRUE))
  })
  
  output$lsc_smd <- renderPlotly({
    plot_ly(data_lsc_sum[, c("original_smd", "matched_smd", "variable_name")], x = ~original_smd, y = ~matched_smd, text = ~variable_name, type = "scatter", mode = "markers") |>
      layout(title = "SMD's for Large Scale Characterisation",
             xaxis = list(title = "Before Matching"),
             yaxis = list(title = "After Matching"),
             shapes = list(
               list(
                 type = "line",
                 x0 = 0,
                 x1 = 1,
                 xref = "container",
                 y0 = 0.1,
                 y1 = 0.1,
                 line = list(color = "black")
               ),
               list(
                 type = "line",
                 x0 = 0,
                 x1 = 0.5,
                 xref = "container",
                 y0 = 0,
                 y1 = 0.5,
                 line = list(color = "red")
               )
             )
      )
    
    
  })
  
  
  # end ---
}

shinyApp(ui = ui, server = server)