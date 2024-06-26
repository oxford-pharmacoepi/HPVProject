library(survival)
library(epiR)
#library(sjPlot)

cdm <- omopgenerics::bind(
  cdm$unvac_18f_coverage_0_any_dose |> mutate(treatment = 0, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(18)) |> dplyr::select(!date_of_birth),
  cdm$vac_18f_coverage_0_any_dose |> mutate(treatment = 1, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(18)) |> dplyr::select(!date_of_birth),
  cdm$matched_unvaccinated_coverage_0_18f_any_dose |> mutate(treatment = 0, group = "matched_overall"),
  cdm$matched_vaccinated_coverage_0_18f_any_dose |> mutate(treatment = 1, group = "matched_overall"),
  cdm$matched_unvac_cov_0_ratio_18f_any_dose |> mutate(treatment = 0, group = "matched_overall"),
  cdm$matched_vac_cov_0_ratio_18f_any_dose |> mutate(treatment = 1, group = "matched_overall"),
  cdm$not_matched_vaccinated_18f_coverage_0_1dose |> mutate(treatment = 0, group = "initial_1vs23") |> mutate(cohort_start_date = date_dose_1),
  cdm$not_matched_vaccinated_18f_coverage_0_23dose |> mutate(treatment = 1, group = "initial_1vs23") |> mutate(cohort_start_date = date_dose_1),
  cdm$not_matched_vaccinated_18f_coverage_0_2dose |> mutate(treatment = 0, group = "initial_2vs3") |> mutate(cohort_start_date = date_dose_1),
  cdm$not_matched_vaccinated_18f_coverage_0_3dose |> mutate(treatment = 1, group = "initial_2vs3") |> mutate(cohort_start_date = date_dose_1),
  cdm$matched_vaccinated_18f_coverage_0_1dose |> mutate(treatment = 0, group = "matched_1vs23"),
  cdm$matched_vaccinated_18f_coverage_0_23dose |> mutate(treatment = 1, group = "matched_1vs23"),
  cdm$matched_vaccinated_18f_coverage_0_2dose |> mutate(treatment = 0, group = "matched_2vs3"),
  cdm$matched_vaccinated_18f_coverage_0_3dose |> mutate(treatment = 1, group = "matched_2vs3"),
  name = "hpv_population_cohorts"
)

cdm <- omopgenerics::bind(
  cdm$unvac_15m_coverage_0_any_dose |> mutate(treatment = 0, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(15)) |> dplyr::select(!date_of_birth),
  cdm$vac_15m_coverage_0_any_dose |> mutate(treatment = 1, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(15)) |> dplyr::select(!date_of_birth),
  cdm$matched_unvac_cov_0_ratio_15m_any_dose |> mutate(treatment = 0, group = "matched_overall"),
  cdm$matched_vac_cov_0_ratio_15m_any_dose |> mutate(treatment = 1, group = "matched_overall"),
  cdm$not_matched_vaccinated_15m_coverage_0_1dose |> mutate(treatment = 0, group = "initial_1vs23") |> mutate(cohort_start_date = date_dose_1),
  cdm$not_matched_vaccinated_15m_coverage_0_23dose |> mutate(treatment = 1, group = "initial_1vs23") |> mutate(cohort_start_date = date_dose_1),
  cdm$not_matched_vaccinated_15m_coverage_0_2dose |> mutate(treatment = 0, group = "initial_2vs3") |> mutate(cohort_start_date = date_dose_1),
  cdm$not_matched_vaccinated_15m_coverage_0_3dose |> mutate(treatment = 1, group = "initial_2vs3") |> mutate(cohort_start_date = date_dose_1),
  cdm$matched_vaccinated_15m_coverage_0_1dose |> mutate(treatment = 0, group = "matched_1vs23"),
  cdm$matched_vaccinated_15m_coverage_0_23dose |> mutate(treatment = 1, group = "matched_1vs23"),
  cdm$matched_vaccinated_15m_coverage_0_2dose |> mutate(treatment = 0, group = "matched_2vs3"),
  cdm$matched_vaccinated_15m_coverage_0_3dose |> mutate(treatment = 1, group = "matched_2vs3"),
  name = "hpv_population_cohorts"
)

info(logger, "PREPARE OUTCOME TABLE")

# Real outcome cohorts -------------
sex <- "female"
type_outcomes <- "darwin"
resultsFolder <- paste0("Results_", type_outcomes, "_outcomes_",sex)
outcome_list <- settings(cdm[[paste0(type_outcomes, "_outcome_cohorts")]]) |> pull(cohort_name)
outcome_names_list <- settings(cdm[[paste0(type_outcomes, "_outcome_cohorts")]]) |> dplyr::select(cohort_definition_id, cohort_name) |> rename(outcome_name = cohort_name)

cdm$outcome_cohorts <- cdm[[paste0(type_outcomes, "_outcome_cohorts")]] |>
  compute(name = "outcome_cohorts", temporary = FALSE) |>
  left_join(outcome_names_list, by = "cohort_definition_id", copy = TRUE)

############################################## ----------------
full_outcome_table <- tibble("cohort_definition_id" = as.numeric(),
                             "group" = as.character(),
                             "subject_id" = as.numeric(),
                             "cohort_start_date" =  as.Date(character()),
                             "cohort_end_date" =  as.Date(character()), 
                             "cohort_year" = as.numeric(),
                             "treatment" = as.numeric(),
                             "pair_id" = as.numeric(),
                             "censored_date" = as.Date(character()),
                             "censor_reason" = as.character(),
                             
)


cdm <- cdm |>
  CDMConnector::insertTable(
    name = "full_outcome_table",
    table = full_outcome_table,
    overwrite = TRUE
  ) 

info(logger, "Prepare overall outcome table")

# INITIAL OVERALL ############---------
study_groups <- c("initial_overall", "initial_1vs23", "initial_2vs3")

for (study_group in study_groups) {
  cdm$total_initial_cohort <- cdm$hpv_population_cohorts |>
    filter(group == study_group) |>
    compute(name = "total_initial_cohort")
  
  print(0)
  # Create next dose
  next_dose <- cdm$total_initial_cohort |>
    filter(treatment == 0) |>
    addTableIntersectDate(tableName = "doses_allvac_cohort",
                          nameStyle = "next_dose_date",
                          window = c(1,Inf),
                          order = "first") |>
    dplyr::select(subject_id, next_dose_date)
  
  # Add next dose to outcome table
  cdm$partial_outcome_table <- cdm$total_initial_cohort |>
    left_join(next_dose, by = "subject_id") |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Add death date
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    left_join(cdm$death |>
                dplyr::select(person_id, death_date) |>
                rename(subject_id = person_id)) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Add event column: death, end obs, next dose
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    mutate(event_date = pmin(death_date, cohort_end_date, next_dose_date, na.rm = TRUE)) |>
    mutate(event_type = case_when(
      event_date == death_date ~ "death",
      event_date == cohort_end_date ~ "end of observation",
      event_date == next_dose_date ~ "next dose"
    )) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    mutate(censored_date = event_date) |>
    mutate(censor_reason = event_type) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    addCohortIntersectDate(targetCohortTable = "outcome_cohorts") |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Save for this group
  cdm$full_outcome_table <- cdm$full_outcome_table |>
    union_all(cdm$partial_outcome_table) |>
    dplyr::select(!c(next_dose_date, death_date, event_date, event_type)) |>
    compute(name = "full_outcome_table", temporary = FALSE)
}
# ------------

# MATCHED OVERALL ############---------
study_groups <- c("matched_overall", "matched_1vs23", "matched_2vs3")

for (study_group in study_groups) {
  cdm$total_matched_cohort <- cdm$hpv_population_cohorts |>
    filter(group == study_group) |>
    compute(name = "total_matched_cohort")
  
  print(0)
  # Create next dose
  next_dose <- cdm$total_matched_cohort |>
    filter(treatment == 0) |>
    addTableIntersectDate(tableName = "doses_allvac_cohort",
                          nameStyle = "next_dose_date",
                          window = c(1,Inf),
                          order = "first") |>
    dplyr::select(subject_id, next_dose_date)
  
  # Add next dose to outcome table
  cdm$partial_outcome_table <- cdm$total_matched_cohort |>
    left_join(next_dose, by = "subject_id") |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Add death date
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    left_join(cdm$death |>
                dplyr::select(person_id, death_date) |>
                rename(subject_id = person_id)) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Add event column: death, end obs, outcome, next dose
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    mutate(event_date = pmin(death_date, cohort_end_date, next_dose_date, na.rm = TRUE)) |>
    mutate(event_type = case_when(
      event_date == death_date ~ "death",
      event_date == cohort_end_date ~ "end of observation",
      event_date == next_dose_date ~ "next dose"
    )) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Add Censor pair column
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    group_by(pair_id) |>
    mutate(censored_date = min(event_date, na.rm = FALSE)) |>
    ungroup() |>
    mutate(censor_reason = case_when(
      censored_date == event_date ~ event_type,
      censored_date != event_date ~ "pair censored"
    )) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  cdm$partial_outcome_table <- cdm$partial_outcome_table |>
    addCohortIntersectDate(targetCohortTable = "outcome_cohorts") |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  # Save for this group
  cdm$full_outcome_table <- cdm$full_outcome_table |>
    union_all(cdm$partial_outcome_table) |>
    dplyr::select(!c(next_dose_date, death_date, event_date, event_type)) |>
    compute(name = "full_outcome_table", temporary = FALSE)
}
# ------------

# COMPUTE METRICS ---------------
info(logger, "COMPUTE METRICS: IR, IRR (Poisson Regression), Cox, Survival")

groups <- cdm$full_outcome_table |>
  distinct(group) |>
  pull()

# pull returns a vector
for (study_group in groups) {
  cdm$preoutcome_table <- cdm$full_outcome_table |>
    filter(group == study_group) |>
    compute(name = "preoutcome_table", temporary = FALSE)
  
  # PREPARE METRICS
  IR <- tibble()
  IRR <- tibble()
  Cox <- tibble()
  
  for (outcome in outcome_list) {
    cdm$outcome_table <- cdm$preoutcome_table |>
      dplyr::select(c(cohort_definition_id, group, subject_id, cohort_start_date, cohort_end_date, cohort_year, treatment, pair_id, censored_date, censor_reason, paste0(outcome, "_0_to_inf"))) |>
      mutate(time_to_event = pmin(censored_date, !!sym(paste0(outcome, "_0_to_inf")), na.rm = TRUE) - as.Date(cohort_start_date)) |>
      mutate(event_type = case_when(
        censored_date < !!sym(paste0(outcome, "_0_to_inf")) | is.na(!!sym(paste0(outcome, "_0_to_inf"))) ~ censor_reason,
        .default = outcome
      )) |>
      compute(name = "outcome_table", temporary = FALSE)
    
    for (ys in c(5, 10, 15, 33)) {
      # Vaccinated IR
      pop_vac <- cdm$outcome_table |> filter(treatment == 1) |>
        filter(time_to_event > 0) |>
        mutate(time_to_event = case_when(
          time_to_event >= 365.25*ys ~ 365.25*ys,
          time_to_event < 365.25*ys ~ time_to_event
        ))

      den_vac <- pop_vac |> pull(time_to_event) |> sum()/365.25

      num_vac <- pop_vac |>
        filter(event_type == outcome & time_to_event < 365.25*ys) |>
        tally() |>
        pull()

      data <- as.matrix(tibble(as.double(num_vac), as.double(den_vac)))

      ir <- epi.conf(
        data,
        ctype = "inc.rate",
        method = "exact",
        N = tally(pop_vac),
        design = 1,
        conf.level = 0.95
      ) * 100000

      IR <- IR |>
        rbind(as_tibble(ir) |> mutate(outcome = outcome, model = "IR", year = ys, group = "treated", n = num_vac))

      # Unvaccinated IR
      pop_unvac <- cdm$outcome_table |> filter(treatment == 0) |>
        filter(time_to_event > 0) |>
        mutate(time_to_event = case_when(
          time_to_event > 365.25*ys ~ 365.25*ys,
          time_to_event <= 365.25*ys ~ time_to_event
        ))

      den_unvac <- pop_unvac |> pull(time_to_event) |> sum()/365.25

      num_unvac <- pop_unvac |>
        filter(event_type == outcome & time_to_event < 365.25*ys) |>
        tally() |>
        pull()

      data <- as.matrix(tibble(as.double(num_unvac), as.double(den_unvac)))

      ir <- epi.conf(
        data,
        ctype = "inc.rate",
        method = "exact",
        N = tally(pop_unvac),
        design = 1,
        conf.level = 0.95
      ) * 100000

      IR <- IR |>
        rbind(as_tibble(ir) |> mutate(outcome = outcome, model = "IR", year = ys, group = "untreated", n = num_unvac))

      # Poisson Regression
      poisson_model <- glm(formula = outcome_flag ~ treatment, family = poisson(link = "log"), offset = log(time_to_event), data = cdm$outcome_table |>
                             mutate(time_to_event = case_when(
                               time_to_event > 365.25*ys ~ 365.25*ys,
                               time_to_event <= 365.25*ys ~ time_to_event
                             )) |>
                             mutate(outcome_flag = case_when(
                               event_type == outcome & time_to_event < 365.25*ys ~ 1,
                               .default = 0
                             )) |>
                             filter(time_to_event > 0)
      )

      poisson_coef <- poisson_model$coefficients[poisson_model$coefficients != "(Intercept)"]

      IRR <- IRR |>
        rbind(broom::tidy(poisson_model) |> mutate(outcome = outcome, model = "poisson", year = ys, irr = exp(estimate), lower = as.numeric(exp(estimate - 1.96*std.error)), upper = exp(estimate + 1.96*std.error)))
    }

    # Survival Study
    survivalplotdata <- tibble()
    
    for (exposure in c(0,1)) {
      info(logger, "SURVIVAL STUDY")
      info(logger, "Kaplan Meier plots")
      
      survival_data <- survfit(Surv(time_to_event, outcome_flag) ~ treatment, data = cdm$outcome_table |>
                                 filter(treatment == exposure) |>
                                 mutate(outcome_flag = case_when(
                                   event_type == outcome ~ 1,
                                   .default = 0
                                 ))) 
      
      survivalplotdata <- tibble(
        "time_to_event" = survival_data$time,
        "survival" = survival_data$surv,
        "lower" = survival_data$lower,
        "upper" = survival_data$upper
      )
      
      write.csv(survivalplotdata, file = paste0(resultsFolder, "/survivalplotdata_", study_group, "_", outcome,"_treatment_",exposure, ".csv"), row.names = FALSE)
    }

    info(logger, "loglog plot")
    png(paste0(resultsFolder, "/logPlot_", study_group, "_", outcome, ".png"))
    logplot <- plot(survival_data,
                    log = "xy",
                    main = "Kaplan Meier Plot",
                    xlab = "log(Time to event (days))", ylab = "log(Survival probability)",
                    #xlim = c(17, 37), ylim = c(0.75, 1),
                    conf.int = F, mark.time = F,
                    lty = 1:2)
    dev.off()
    
    info(logger, "Cox Regression")
    cox_model <- coxph(Surv(time_to_event, outcome_flag) ~ treatment, data = cdm$outcome_table|>
                         mutate(outcome_flag = case_when(
                           event_type == outcome ~ 1,
                           .default = 0
                         )))
    
    Cox <- Cox |> 
      rbind(broom::tidy(cox_model) |> 
              mutate(outcome = outcome, hr = exp(estimate), lower = exp(estimate - 1.96*std.error), upper = exp(estimate + 1.96*std.error))
      )
    
    print(45)
    #info(logger, "Check proportional Hazards Assumptions")
    # extract dummy variables (0/1 variables and set it upc in the table again): the cox.zph messed up the table
    # cdm$output_table <- cdm$output_table |>
    #   mutate(treatment = model.matrix(cox_model)[, treatment])
    # 
    # cox_model_2 <- coxph(Surv(time_to_event, outcome_flag) ~ treatment, data = cdm$outcome_table)
    #ph_check <- cox.zph(cox_model) # p > 0.05

  }
  
  write.csv(IR, file = paste0(resultsFolder, "/IR_", study_group, ".csv"), row.names = FALSE)
  write.csv(IRR, file = paste0(resultsFolder, "/IRR_", study_group, ".csv"), row.names = FALSE)
  write.csv(Cox, file = paste0(resultsFolder, "/Cox_", study_group, ".csv"), row.names = FALSE)
}

#--------------------