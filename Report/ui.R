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
      ),
      menuItem(
        text = "Outcome Model",
        tabName = "out",
        menuSubItem(
          text = "IR",
          tabName = "ir"
        ),
        menuSubItem(
          text = "IRR",
          tabName = "irr"
        ),
        menuSubItem(
          text = "HR",
          tabName = "hr"
        )
      ),
      menuItem(
        text = "Negative Control Outcomes",
        tabName = "nco",
        menuSubItem(
          text = "IR",
          tabName = "ir_nco"
        ),
        menuSubItem(
          text = "IRR",
          tabName = "irr_nco"
        ),
        menuSubItem(
          text = "HR",
          tabName = "hr_nco"
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
        # a(
        #   img(
        #     src = "https://raw.githubusercontent.com/darwin-eu/DrugUtilisation/main/man/figures/DrugUtilisation.png", 
        #     align = "right", width = "200px"
        #   ),
        #   href = "https://github.com/darwin-eu/DrugUtilisation/",
        #   target = "_blank"
        # ),
        h3("Effectiveness of HPV vaccination in reducing HPV-related cancers and lesions"),
        h5("Background: HPV is the main cause of cervical cancer and is associated with many
other cancers and lesions. Since its recent implementation in the universal vaccination
programme, studies have reported evidence of the reduction of cervical cancers and high
grade CIN. However, those are ecological and low quality studies, and there is still lack
of evidence on the association between the vaccination and HPV-related rare cancers and
lesions."),
        h5("Objectives: I therefore aim to evaluate the effectiveness of HPV vaccination in reducing cervical, oropharyngeal, anal, vaginal, vulval and penile cancers, cervical pre-cancerous
lesions and warts in men and women in the UK. Additionally, it intends to compare the
effectiveness of different dose schedules in reducing the outcomes mentioned."),
        h5("Methods: This project is a cohort study. The data used for this study is from the
CPRD GOLD database, which collects primary care data from 674 practitions in the UK
that is later mapped into the OMOP CDM. Specifically, the data was extracted in July
2023 and has a total of 17.267.137 patients recorded. The research includes (1) an overall analysis for each sex comparing vaccinated and unvaccinated individuals, (2) a first
dose-stratified analysis for each sex comparing individuals vaccinated with one dose and
individuals vaccinated with more than one dose, and (3) a second dose-stratified analysis
for each sex comparing individuals vaccinated with two doses and individuals vaccinated
with more than two doses. The female population under study is required to have been
born in or after 1990 and the male population under study is required to have been born
in or after 2004. Additionally, all individuals need to be in observation from their 9th
birthday to their 18th or 15th birthday, respectively for women and men. To address any
potential confounding effects, people was exact matched by year of birth, region and year
of vaccination, and a yearly PS-matching was also applied using the nearest-neighbour
method with a 0,2 caliper (with PS computed from the Lasso regression selected variables). To address residual confounding, we implement the Negative control outcomes
(NCO) model to later calibrate biased results. The effectiveness of HPV vaccination is
estimated computing incidence rates and incidence rate ratios (using a Poisson regression)
at 5, 10 and 15 years, and a hazard ratio overall for each of the analyses."),
        h5("Results: Characterisation of the populations showed no sign of observed confounding
after the matching procedure in any of the female analyses carried out. I found evidence
of residual confounding in the female overall population, as suggested by the observed
associations between vaccination status and multiple NCOs after matching. Therefore,
calibration was applied to the results to adjust for the differences among vaccinated and
unvaccinated cohorts. The IRR at 15 years reported a lower incidence of CIN2/CIN3 in
the vaccinated cohort compared to the unvaccinated one (IRR = 0.33; 95% CI: 0.20-0.54).
The HR also reported a lower risk of developing these cervical pre-cancerous lesions for
the vaccinated individuals (HR = 0.36; 95% CI: 0.22-0.59). Regarding the dose-stratified
comparisons, the results show no significant association between different dose schedules
and the outcomes of interest. With respect to the male analyses, we do not have enough
power in the study to draw any conclusions."),
        h5("Conclusions: Our findings report a benefit of HPV vaccination suggesting over a 60% reduction of high-grade cervical precancerous lesions in women. The assessment of
the different vaccination schedules shows no benefit of more than one dose administration
after 15 years from vaccination.")
      ),
      # cohort attrition ------
      tabItem(
        tabName = "definitions",
        h4("Attrition for the study cohorts"),
        selectors(
          data = attritions,
          prefix = "att",
          columns = c("cohort_name"),
          multiple = FALSE,
          default = NULL
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Tidy table",
            DT::dataTableOutput("attrition_tidy") |> withSpinner()
          ),
          tabPanel(
            "Diagram",
            grVizOutput("attrition_diagram", width = "400px", height = "100%") %>% withSpinner()
          )
        )
      ),
      # Lasso ----------------------------------
      tabItem(
        tabName = "lasso",
        h4("Lasso Selected Features"),
        selectors(
          data = lasso,
          prefix = "lasso",
          columns = c("study_group", "year")
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
          data = characterisation_to,
          prefix = "to",
          columns = c("study_group"),
          multiple = FALSE
        ),
        selectors(
          data = characterisation_to,
          prefix = "to",
          columns = c("variable_name", "variable_level", "estimate_name", "window"),
          default = list(estimate_name = c("count", "percentage", "median", "q25", "q75")),
          multiple = TRUE
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            selectors(
              data = characterisation_to,
              prefix = "to",
              columns = c("cohort_name"),
              multiple = TRUE
            ),
            DT::dataTableOutput("table_one_raw") |> withSpinner()
          ),
          tabPanel(
            "Tidy data",
            DT::dataTableOutput("table_one_tidy") |> withSpinner()
          )
          # tabPanel(
          #   "SMD plot",
          #   plotlyOutput("to_smd") |> withSpinner()
          # )
        )
      ),
      # lsc -----------------------------------
      tabItem(
        tabName = "lsc",
        h4("Large scale characteristics"),
        selectors(
          data = characterisation_lsc,
          prefix = "lsc",
          columns = c("study_group"),
          multiple = FALSE
        ),
        selectors(
          data = characterisation_lsc,
          prefix = "lsc",
          columns = c("variable_name", "estimate_name", "variable_level"),
          default = list(estimate_name = c("count", "percentage", "median", "q25", "q75"))
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            selectors(
              data = characterisation_lsc,
              prefix = "lsc",
              columns = c("cohort_name"),
              multiple = TRUE
            ),
            DT::dataTableOutput("lsc_raw") |> withSpinner()
          ),
          tabPanel(
            "Tidy data",
            DT::dataTableOutput("lsc_tidy") |> withSpinner()
          ),
          tabPanel(
            "SMD plot",
            h6("Y-axis comparison"),
            selectors(
              data = characterisation_lsc,
              prefix = "ylsc",
              columns = c("study_group"),
              multiple = FALSE
            ),
            plotlyOutput("lsc_smd") |> withSpinner()
          )
        )
      ),
      # Outcome model ---------
      tabItem(
        tabName = "ir",
        h4("IR per 100.000 person-years"),
        selectors(
          data = IR_results,
          prefix = "ir",
          columns = c("study_group","sex","group", "outcome", "year"),
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            DT::dataTableOutput("ir_raw") |> withSpinner()
          )
        )
      ),
      tabItem(
        tabName = "irr",
        h4("Poisson regression (with offset: accounting for time to event)"),
        selectors(
          data = IRR_results,
          prefix = "irr",
          columns = c("sex","year"),
          multiple = FALSE
        ),
        selectors(
          data = IRR_results,
          prefix = "irr",
          columns = c("study_group", "outcome"),
          multiple = TRUE
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            DT::dataTableOutput("irr_raw") |> withSpinner()
          ),
          tabPanel(
            "Forest Plot",
            plotOutput("irr_forestplot") |> withSpinner()
          )
        )
      ),
      tabItem(
        tabName = "hr",
        h4("HR from Cox model"),
        selectors(
          data = HR_results,
          prefix = "hr",
          columns = c("sex"),
          multiple = FALSE
        ),
        selectors(
          data = HR_results,
          prefix = "hr",
          columns = c("study_group", "outcome"),
          multiple = TRUE
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            DT::dataTableOutput("hr_raw") |> withSpinner()
          ),
          tabPanel(
            "Forest Plot",
            plotOutput("hr_forestplot") |> withSpinner()
          )
        ),
      ),
      tabItem(
        tabName = "ir_nco",
        h4("IR per 100.000 person-years"),
        selectors(
          data = IR_nco,
          prefix = "ir_nco",
          columns = c("study_group", "sex", "group", "outcome", "year"),
          default = list(year = 15)
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw data",
            DT::dataTableOutput("ir_nco_raw") |> withSpinner()
          ),
        )
      ),
      tabItem(
        tabName = "irr_nco",
        h4("Poisson regression (with offset: accounting for time to event)"),
        selectors(
          data = IRR_nco,
          prefix = "irr_nco",
          columns = c("study_group","sex","year"),
          multiple = FALSE,
          default = list(year = 15)
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Forest Plot",
            plotOutput("irr_nco_plot") |> withSpinner()
          ),
        )
      ),
      tabItem(
        tabName = "hr_nco",
        h4("HR from Cox model"),
        selectors(
          data = HR_nco,
          prefix = "hr_nco",
          columns = c("study_group","sex"),
          multiple = FALSE
          ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Forest Plot",
            plotOutput("hr_nco_plot") |> withSpinner()
          ),
        )
       )
      # end ---
    )
  )
)