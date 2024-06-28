# Server ---------------------------------------------------------------------
server <- function(input, output, session) {
  # definitions: cohort attritions -------------
  
  attritions_tidy <- reactive({
    attritions |>
      filter(cohort_name == input$att_cohort_name)
    })

  output$attrition_tidy <- DT::renderDataTable({
    DT::datatable(attritions_tidy(), options = list(scrollX = TRUE))
  })
  
  output$attrition_diagram <- renderGrViz({
    cdi <- attritions |>
      filter(cohort_name == input$att_cohort_name) |>
      pull(cohort_definition_id) |> unique()
    
    plotCohortAttrition_mah(attritions |>
                              filter(cohort_name == input$att_cohort_name), cdi)
  })
  
  
  # lasso: matching conditions selected -----------------------------------
  # Select desired data
  
  #Filter data
  filtered_data <- reactive({
    lasso |>
      filter(year %in% input$lasso_year & study_group %in% input$lasso_study_group)
  })

  output$lasso_table <- DT::renderDataTable({
    DT::datatable(filtered_data(), rownames = FALSE)
  })


  # Table one ------------------------------------
  # Computations
  TableOne_raw <- reactive({
    characterisation_to|>
      filter(cohort_name %in% input$to_cohort_name, estimate_name %in% input$to_estimate_name, variable_name %in% input$to_variable_name, study_group %in% input$to_study_group, window %in% input$to_window)
  })

  reactive({
    print(input$to_study_group)
  })
  TableOne_tidy <- reactive({
    characterisation_to_smd[[input$to_study_group]] |>
      filter(variable_name %in% input$to_variable_name, variable_level %in% input$to_variable_level, window %in% input$to_window)
  })

  # Outputs
  output$table_one_raw <- DT::renderDataTable({
    DT::datatable(TableOne_raw(), options = list(scrollX = TRUE))
  })

  output$table_one_tidy <- DT::renderDataTable({
    DT::datatable(TableOne_tidy(), options = list(scrollX = TRUE))
  })

  # lsc ---------------------------------
  LargeScale_raw <- reactive({
    characterisation_lsc |>
      filter(estimate_name %in% input$lsc_estimate_name, variable_name %in% input$lsc_variable_name, variable_level %in% input$lsc_variable_level, study_group %in% input$lsc_study_group, cohort_name %in% input$lsc_cohort_name)
  })


  LargeScale_tidy <- reactive({
    characterisation_lsc_smd[[input$lsc_study_group]] |>
      filter(variable_name %in% input$lsc_variable_name, variable_level %in% input$lsc_variable_level)
  })

  output$lsc_raw <- DT::renderDataTable({
    datatable(LargeScale_raw(), options = list(scrollX = TRUE))
  })

  output$lsc_tidy <- DT::renderDataTable({
    datatable(LargeScale_tidy(), options = list(scrollX = TRUE))
  })

  output$lsc_smd <- renderPlotly({
    plot_ly(data = left_join(characterisation_lsc_smd[[input$lsc_study_group]] |>
                               dplyr::select(variable_name, variable_level, SMDs_x = SMDs),
                             characterisation_lsc_smd[[input$ylsc_study_group]] |>
                               dplyr::select(variable_name, variable_level, SMDs_y = SMDs),
                             by = c("variable_name", "variable_level")
                             ),
            x = ~SMDs_x, y = ~SMDs_y,
            text = ~variable_name, type = "scatter", 
            mode = "markers", 
            marker = list(symbol = "circle-open"), 
            name = "1:5 matched"
            ) |>
      layout(title = "SMD's for Large Scale Characterisation",
             xaxis = list(title = gsub("_", " ", input$lsc_study_group)),
             yaxis = list(title = gsub("_", " ", input$ylsc_study_group)),
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

  # OUTCOME MODEL ---------------

  IR_raw <- reactive({
    IR_results |>
      filter(study_group %in% input$ir_study_group, sex %in% input$ir_sex, group %in% input$ir_group,  outcome %in% input$ir_outcome, year %in% input$ir_year)
  })

  output$ir_raw <- DT::renderDataTable({
    datatable(IR_raw(), options = list(scrollX = TRUE))
  })

  IRR_raw <- reactive({
    IRR_results |>
      filter(study_group %in% input$irr_study_group, sex %in% input$irr_sex, outcome %in% input$irr_outcome, year %in% input$irr_year)
  })
  
    output$irr_raw <- DT::renderDataTable({
    datatable(IRR_raw(), options = list(scrollX = TRUE))
  })
    
  output$irr_forestplot <- renderPlot({
    name_file <- case_when(
      grepl("overall",input$irr_study_group) ~ "Vaccinated versus unvaccinated",
      grepl("1vs23",input$irr_study_group) ~ "1 versus more-than-1",
      grepl("2vs3",input$irr_study_group) ~ "2 versus more than-2"
    )
    ylimit <- case_when(
      grepl("overall",input$irr_study_group) ~ 3,
      grepl("1vs23",input$irr_study_group) ~ 10,
      grepl("2vs3",input$irr_study_group) ~ 10
    )
    ylimit <- max(ylimit)
    
    forestplot_data <- IRR_results |>
      filter(study_group %in% input$irr_study_group, sex %in% input$irr_sex, year %in% input$irr_year, outcome %in% input$irr_outcome)

    forestplot(
      df = forestplot_data,
      name = outcome,
      estimate = estimate,
      se = std.error,
      logodds = TRUE,
      colour = study_group,
      title = paste0(input$irr_sex," analyses: IRR at ",input$irr_year," years"),
      xlab = "Effect size",
      xlim = c(0.1,ylimit),
      
    ) +
    scale_color_manual(values = c("#fdbf69","#002147","#226584","#7eb3eb","#669900","#660099", "#660033", "#333300", "#003300")) +
    scale_fill_manual(values = c("#fdbf69","#002147","#226584","#7eb3eb","#669900","#660099", "#660033", "#333300", "#003300")) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5),  # Adjust the size as needed
      axis.title.x = element_text(size = 12),  # Adjust the font size for the x-axis title
      axis.title.y = element_text(size = 12),  # Adjust the font size for the y-axis title
      axis.text.x = element_text(size = 11),   # Adjust the font size for the x-axis labels
      axis.text.y = element_text(size = 11),    # Adjust the font size for the y-axis labels
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
      )

  }, height = 700)



  HR_raw <- reactive({
    HR_results |>
      filter(study_group %in% input$hr_study_group, sex %in% input$hr_sex, outcome %in% input$hr_outcome)
  })

  output$hr_raw <- DT::renderDataTable({
    datatable(HR_raw(), options = list(scrollX = TRUE))
  })
  
  output$hr_forestplot <- renderPlot({
    name_file <- case_when(
      grepl("overall",input$hr_study_group) ~ "Vaccinated versus unvaccinated",
      grepl("1vs23",input$hr_study_group) ~ "1 versus more-than-1",
      grepl("2vs3",input$hr_study_group) ~ "2 versus more than-2"
    )

    ylimit <- case_when(
      grepl("overall",input$hr_study_group) ~ 3,
      grepl("1vs23",input$hr_study_group) ~ 10,
      grepl("2vs3",input$hr_study_group) ~ 10
    )
    ylimit <- max(ylimit)
    
    forestplot_data <- HR_results |>
      filter(study_group %in% input$hr_study_group, sex %in% input$hr_sex, outcome %in% input$hr_outcome)

    forestplot(
      df = forestplot_data,
      name = outcome,
      estimate = estimate,
      se = std.error,
      logodds = TRUE,
      colour = study_group,
      title = paste0(input$hr_sex," analyses: HR estimate"),
      xlab = "Effect size",
      xlim = c(0.1,ylimit),
      
    ) +
      scale_color_manual(values = c("#fdbf69","#002147","#226584","#7eb3eb","#669900","#660099", "#660033", "#333300", "#003300")) +
      scale_fill_manual(values = c("#fdbf69","#002147","#226584","#7eb3eb","#669900","#660099", "#660033", "#333300", "#003300")) +
      theme(
        plot.title = element_text(size = 13, hjust = 0.5),  # Adjust the size as needed
        axis.title.x = element_text(size = 12),  # Adjust the font size for the x-axis title
        axis.title.y = element_text(size = 12),  # Adjust the font size for the y-axis title
        axis.text.x = element_text(size = 11),   # Adjust the font size for the x-axis labels
        axis.text.y = element_text(size = 11),    # Adjust the font size for the y-axis labels
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      )
    
  }, height = 700)
  
  #-----------

  # NEGATIVE CONTROL OUTCOME ---------------

  IR_nco_raw <- reactive({
    IR_nco |>
      filter(study_group %in% input$ir_nco_study_group, sex %in% input$ir_nco_sex, group %in% input$ir_nco_group, year %in% input$ir_nco_year, outcome %in% input$ir_nco_outcome)
  })
  
  output$ir_nco_raw <- DT::renderDataTable({
    datatable(IR_nco_raw(), options = list(scrollX = TRUE))
  })
  
  output$irr_nco_plot <- renderPlot({
    name_file <- case_when(
      grepl("overall",input$irr_nco_study_group) ~ "Vaccinated versus unvaccinated",
      grepl("1vs23",input$irr_nco_study_group) ~ "1 versus more-than-1",
      grepl("2vs3",input$irr_nco_study_group) ~ "2 versus more than-2"
    )
    ylimit <- case_when(
      grepl("overall",input$irr_nco_study_group) ~ 3,
      grepl("1vs23",input$irr_nco_study_group) ~ 10,
      grepl("2vs3",input$irr_nco_study_group) ~ 10
    )
    ylimit <- max(ylimit)
    
    forestplot_data <- IRR_nco |>
      filter(study_group %in% input$irr_nco_study_group, sex %in% input$irr_nco_sex, year %in% input$irr_nco_year) |>
      mutate(bias = case_when(
        lower >= 1 | upper <= 1 ~ "bias",
        .default = "no bias"
      ))
    
    forestplot(
      df = forestplot_data,
      name = outcome,
      estimate = estimate,
      se = std.error,
      logodds = TRUE,
      colour = bias,
      title = paste0(input$irr_nco_sex," analyses: IRR at ",input$irr_nco_year," years"),
      xlab = "Effect size",
      xlim = c(0.1,ylimit),
      
    ) +
      scale_color_manual(values = c("no bias" = "#002147", "bias" = "#fdbf69")) +
      scale_fill_manual(values = c("no bias" = "#002147", "bias" = "#fdbf69")) +
      theme(
        plot.title = element_text(size = 13, hjust = 0.5),  # Adjust the size as needed
        axis.title.x = element_text(size = 12),  # Adjust the font size for the x-axis title
        axis.title.y = element_text(size = 12),  # Adjust the font size for the y-axis title
        axis.text.x = element_text(size = 11),   # Adjust the font size for the x-axis labels
        axis.text.y = element_text(size = 11),    # Adjust the font size for the y-axis labels
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      )
    
  }, height = 700)

  output$hr_nco_plot <- renderPlot({
    name_file <- case_when(
      grepl("overall",input$hr_nco_study_group) ~ "Vaccinated versus unvaccinated",
      grepl("1vs23",input$hr_nco_study_group) ~ "1 versus more-than-1",
      grepl("2vs3",input$hr_nco_study_group) ~ "2 versus more than-2"
    )
    ylimit <- case_when(
      grepl("overall",input$hr_nco_study_group) ~ 3,
      grepl("1vs23",input$hr_nco_study_group) ~ 10,
      grepl("2vs3",input$hr_nco_study_group) ~ 10
    )
    ylimit <- max(ylimit)
    
    forestplot_data <- HR_nco |>
      filter(study_group %in% input$hr_nco_study_group, sex %in% input$hr_nco_sex) |>
      mutate(bias = case_when(
        lower >= 1 | upper <= 1 ~ "bias",
        .default = "no bias"
      ))
    
    forestplot(
      df = forestplot_data,
      name = outcome,
      estimate = estimate,
      se = std.error,
      logodds = TRUE,
      colour = bias,
      title = paste0(input$hr_nco_sex," analyses: HR estimate"),
      xlab = "Effect size",
      xlim = c(0.1,ylimit),
      
    ) +
      scale_color_manual(values = c("no bias" = "#002147", "bias" = "#fdbf69")) +
      scale_fill_manual(values = c("no bias" = "#002147", "bias" = "#fdbf69")) +
      theme(
        plot.title = element_text(size = 13, hjust = 0.5),  # Adjust the size as needed
        axis.title.x = element_text(size = 12),  # Adjust the font size for the x-axis title
        axis.title.y = element_text(size = 12),  # Adjust the font size for the y-axis title
        axis.text.x = element_text(size = 11),   # Adjust the font size for the x-axis labels
        axis.text.y = element_text(size = 11),    # Adjust the font size for the y-axis labels
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      )
    
  }, height = 700)
  # -----------
  
  # end ---
}
