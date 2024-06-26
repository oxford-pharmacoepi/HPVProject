library(survival)
library(epiR)
#library(sjPlot)

# cdm <- cdmFromCon(
#   con = db,
#   cdmSchema = cdmSchema, 
#   writeSchema = writeSchema, 
#   cdmName = dbName, 
#   achillesSchema = achillesSchema, 
#   cohortTables = c("vac_cohort", "unvac_cohort", "total_vac_matched_cohort", "total_unvac_matched_cohort", "doses_allvac_cohort", "doses1_matched_cohort", "doses2_matched_cohort", "dose2_matched_cohort", "doses3_matched_cohort")
# )

cdm <- omopgenerics::bind(
  cdm$total_unvac_matched_cohort |> mutate(cohort_start_date = index_date) |> rename(treatment = vac_status) |> mutate(group = "overall"),
  cdm$total_vac_matched_cohort |> mutate(cohort_start_date = index_date) |> rename(treatment = vac_status) |> mutate(group = "overall"),
  cdm$dose1_matched_cohort |> mutate(cohort_start_date = index_date) |> mutate(group = "dose1vs23"),
  cdm$doses23_matched_cohort |> mutate(cohort_start_date = index_date) |> mutate(group = "dose1vs23"),
  cdm$dose2_matched_cohort |> mutate(cohort_start_date = index_date) |> mutate(group = "dose2vs3"),
  cdm$doses3_matched_cohort |> mutate(cohort_start_date = index_date) |> mutate(group = "dose2vs3"),
  name = "hpv_population_cohorts"
)

info(logger, "PREPARE OUTCOME TABLE")
# Fake outcome cohorts ################### ---------------
# Most common condition
# code <- cdm$total_matched_cohort |>
#   dplyr::select(subject_id) |>
#   right_join(cdm$condition_occurrence |>
#                dplyr::select(person_id, condition_concept_id) |>
#                rename(subject_id = person_id), by = "subject_id") |>
#   group_by(condition_concept_id) |>
#   tally() |>
#   filter(n == max(n)) |>
#   pull(condition_concept_id)
# 
# code_2 <- cdm$total_matched_cohort |>
#   dplyr::select(subject_id) |>
#   right_join(cdm$condition_occurrence |>
#                dplyr::select(person_id, condition_concept_id) |>
#                rename(subject_id = person_id), by = "subject_id") |>
#   group_by(condition_concept_id) |>
#   tally() |>
#   filter(! n == max(n)) |>
#   filter(! n == max(n)) |>
#   filter(! n == max(n)) |>
#   filter( n == max(n)) |>
#   pull(condition_concept_id)

# Flu cohort
# cdm <- generateConceptCohortSet(
#   cdm,
#   conceptSet = list("flu" = 4214962),
#   name = "flu"
# )
# 
# cdm <- generateConceptCohortSet(
#   cdm,
#   conceptSet = list("cough" = 254761),
#   name = "cough"
# )
# 
# cdm <- omopgenerics::bind(
#   cdm$flu,
#   cdm$cough,
#   name = "hpv_outcome_cohorts"
# )
#--------------------------------------

# Real outcome cohorts -------------
# info(logger, "GENERATE OUTOME COHORTS")
# cohort_json_dir <- "PhenotypeR/Cohorts/HPV_Outc_/"
# cohort_set <- read_cohort_set(cohort_json_dir)
# 
# cdm <-   generateCohortSet(cdm,
#                            cohort_set,
#                            name = "hpv_outcome_cohorts",
#                            computeAttrition = TRUE,
#                            overwrite = TRUE)


# --------------------------------
cdm$hpv_outcome_cohorts <- cdm$outcome_cohorts |> compute(name = "hpv_outcome_cohorts", temporary = FALSE)

outcome_list <- settings(cdm$hpv_outcome_cohorts) |> pull(cohort_name)
outcome_names_list <- settings(cdm$hpv_outcome_cohorts) |> dplyr::select(cohort_definition_id, cohort_name) |> rename(outcome_name = cohort_name)

cdm$hpv_outcome_cohorts <- cdm$hpv_outcome_cohorts |>
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
                             "year_of_birth" = as.numeric(),
                             "next_dose_date" = as.Date(character()),
                             "death_date" = as.Date(character()),
                             "outcome_date" = as.Date(character()),
                             "outcome_flag" = as.numeric(),
                             "outcome_name" = as.character(),
                             "event_date" = as.Date(character()),
                             "event_type" = as.character(),
                             "censored_date" = as.Date(character()),
                             "censor_reason" = as.character(),
                             "time_to_event" = as.numeric()
                             
)


cdm <- cdm |>
  CDMConnector::insertTable(
    name = "full_outcome_table",
    table = full_outcome_table,
    overwrite = TRUE
  ) 

info(logger, "Prepare overall outcome table")
# OVERALL ############---------
study_group <- "overall"

cdm$total_matched_cohort <- cdm$hpv_population_cohorts |>
  filter(group == study_group) |>
  dplyr::select(! c(vac_status, index_date)) |>
  compute(name = "total_matched_cohort")

print(0)
# Create next dose
next_dose <- cdm$total_matched_cohort |>
  filter(treatment == 0) |>
  left_join(cdm$doses_allvac_cohort |>
              dplyr::select(subject_id, cohort_start_date) |>
              rename(next_dose_date = cohort_start_date), by = "subject_id") |>
  group_by(subject_id) |>
  filter(next_dose_date == min(next_dose_date) | is.na(next_dose_date)) |>
  ungroup() |>
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

# Add outcome column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  left_join(cdm$hpv_outcome_cohorts |>
              dplyr::select(subject_id, cohort_start_date, outcome_name) |>
              rename(outcome_date = cohort_start_date), by = "subject_id") |>
  #filter(outcome_date >= cohort_start_date) |> #with our outomces this should not be necessary but should be checked that we do not outcomes before index date
  mutate(outcome_flag = case_when(
    is.na(outcome_date) ~ 0,
    !is.na(outcome_date) ~ 1
  )) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Add event column: death, end obs, outcome, next dose
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
  mutate(event_type = case_when(
    event_date == death_date ~ "death",
    event_date == cohort_end_date ~ "end of observation",
    event_date == outcome_date ~ "outcome",
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

# Add time to event column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  mutate(time_to_event = censored_date - cohort_start_date) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Save for this group
cdm$full_outcome_table <- cdm$full_outcome_table |>
  union_all(cdm$partial_outcome_table) |>
  compute(name = "full_outcome_table", temporary = FALSE)
# ------------

info(logger, "Prepare dose 1vs23 outcome table")
# DOSE stratas: 1 vs 23 ############ ---------
study_group <- "dose1vs23"

cdm$total_matched_cohort <- cdm$hpv_population_cohorts |>
  filter(group == study_group) |>
  compute(name = "total_matched_cohort") |>
  dplyr::select(! c(vac_status, index_date))

print(0)
# Create next dose
next_dose <- cdm$total_matched_cohort |>
  filter(treatment == 0) |>
  left_join(cdm$doses_allvac_cohort |>
              dplyr::select(subject_id, cohort_start_date) |>
              rename(next_dose_date = cohort_start_date), by = "subject_id") |>
  group_by(subject_id) |>
  filter(! next_dose_date == min(next_dose_date)) |>
  filter(next_dose_date == min(next_dose_date)  | is.na(next_dose_date)) |>
  ungroup() |>
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

# Add outcome column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  left_join(cdm$hpv_outcome_cohorts |>
              dplyr::select(subject_id, cohort_start_date, outcome_name) |>
              rename(outcome_date = cohort_start_date), by = "subject_id") |>
  filter(outcome_date >= cohort_start_date) |>
  mutate(outcome_flag = case_when(
    is.na(outcome_date) ~ 0,
    !is.na(outcome_date) ~ 1
  )) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Add event column: death, end obs, outcome, next dose
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
  compute(name = "partial_outcome_table") |>
  mutate(event_type = case_when(
    event_date == death_date ~ "death",
    event_date == cohort_end_date ~ "end of observation",
    event_date == outcome_date ~ "outcome",
    event_date == next_dose_date ~ "next dose"
  )) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Add Censor pair column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  group_by(pair_id) |>
  mutate(censored_date = min(event_date, na.rm = FALSE)) |>
  ungroup() |>
  compute(name = "partial_outcome_table") |>
  mutate(censor_reason = case_when(
    censored_date == event_date ~ event_type,
    censored_date != event_date ~ "pair censored"
  )) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Add time to event column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  mutate(time_to_event = censored_date - cohort_start_date) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Save for this group
cdm$full_outcome_table <- cdm$full_outcome_table |>
  union_all(cdm$partial_outcome_table) |>
  compute(name = "full_outcome_table", temporary = FALSE)
# -------------------

info(logger, "Prepare dose 2vs3 outcome table")
# DOSE stratas: 2 vs 3 ############ ---------------
study_group <- "dose2vs3"

cdm$total_matched_cohort <- cdm$hpv_population_cohorts |>
  filter(group == study_group) |>
  compute(name = "total_matched_cohort") |>
  dplyr::select(! c(vac_status, index_date))

print(0)
# Create next dose
next_dose <- cdm$total_matched_cohort |>
  filter(treatment == 0) |>
  left_join(cdm$doses_allvac_cohort |>
              dplyr::select(subject_id, cohort_start_date) |>
              rename(next_dose_date = cohort_start_date), by = "subject_id") |>
  group_by(subject_id) |>
  filter(! next_dose_date == min(next_dose_date)) |>
  filter(! next_dose_date == min(next_dose_date)) |>
  ungroup() |>
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

# Add outcome column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  left_join(cdm$hpv_outcome_cohorts |>
              dplyr::select(subject_id, cohort_start_date, outcome_name) |>
              rename(outcome_date = cohort_start_date), by = "subject_id") |>
  filter(outcome_date >= cohort_start_date) |>
  mutate(outcome_flag = case_when(
    is.na(outcome_date) ~ 0,
    !is.na(outcome_date) ~ 1
  )) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Add event column: death, end obs, outcome, next dose
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
  mutate(event_type = case_when(
    event_date == death_date ~ "death",
    event_date == cohort_end_date ~ "end of observation",
    event_date == outcome_date ~ "outcome",
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

# Add time to event column
cdm$partial_outcome_table <- cdm$partial_outcome_table |>
  mutate(time_to_event = censored_date - cohort_start_date) |>
  compute(name = "partial_outcome_table", temporary = FALSE)

# Save for this group
cdm$full_outcome_table <- cdm$full_outcome_table |>
  union_all(cdm$partial_outcome_table) |>
  compute(name = "full_outcome_table", temporary = FALSE)

# -----------------------

info(logger, "Prepare year of vaccination strata outcome table")
# YEAR OF VAC stratas ############ --------
cdm$total_matched_cohort <- cdm$full_outcome_table |>
  filter(group == "overall") |>
  compute(name = "total_matched_cohort", temporary = FALSE)

for (ys in c(2008:2023)) {
  print(ys)
  cdm$partial_outcome_table <- cdm$total_matched_cohort |>
    filter(cohort_year == ys) |>
    mutate(group = paste0("yov",ys)) |>
    compute(name = "partial_outcome_table", temporary = FALSE)
  
  #print(cdm$partial_outcome_table |> tally() |> pull())
  
  cdm$full_outcome_table <- cdm$full_outcome_table |>
    union_all(cdm$partial_outcome_table) |>
    compute(name = "full_outcome_table", temporary = FALSE)
}
# -----------------------

# COMPUTE METRICS ---------------
info(logger, "COMPUTE METRICS: IR, IRR (Poisson Regression), Cox, Survival")

groups <- cdm$full_outcome_table |>
  distinct(group) |>
  pull()

# pull returns a vector

# PREPARE METRICS
IR <- tibble("study_group" = as.character(),
             "treatment" = as.numeric(),
             "outcome" = as.character(),
             "years" = as.numeric(),
             "est" = as.numeric(),
             "lower" = as.numeric(),
             "upper" = as.numeric()
)

IRR <- tibble("study_group" = as.character(),
             "outcome" = as.character(),
             "years" = as.numeric(),
             "est" = as.numeric(),
             "lower" = as.numeric(),
             "upper" = as.numeric()
)

Cox <- tibbble("study_group" = as.character(),
          "outcome" = as.character(),
          "coef" = as.numeric(),
          "HR" = as.numeric(),
          "lower" = as.numeric(),
          "upper" = as.numeric(),
          "p-value" = as.numeric()
)

for (study_group in groups) {
  cdm$outcome_table <- cdm$full_outcome_table |>
    filter(group == study_group) |>
    compute(name = "outcome_table", temporary = FALSE)
  

  print(1)
  
  for (outcome in outcome_list) {
    #outcome <- "flu"
    print(2)
    for (ys in c(5, 10, 15)) {
      print(3)
      #ys <- 5
      # Vaccinated IR
      treatment <- 1
      pop_vac <- cdm$outcome_table |> filter(treatment == 1) |>
        mutate(time_to_event = case_when(
          time_to_event >= 365*ys ~ 365*ys,
          time_to_event < 365*ys ~ time_to_event
        ))
      
      den_vac <- pop_vac |> pull(time_to_event) |> sum()/365
      
      num_vac <- pop_vac |>
        filter(outcome_name == outcome & censor_reason == "outcome" & time_to_event < 365*ys) |>
        tally() |> pull()
      
      metric <- paste0("IR_",ys,"y_treatment1_",outcome)
      
      #IR <- as.numeric(num_vac*10^5/(den_vac))
      
      metric <- paste0("IR_CI_",ys,"y_treatment1_",outcome)
      
      data <- pop_vac |>
        mutate(outcome_flag = case_when(
          time_to_event <= 365*ys ~ outcome_flag,
          time_to_event > 365*ys ~ 0
        )) |>
        dplyr::select(outcome_flag, time_to_event) |>
        collect() |>
        as.matrix()
      
      data <- as.matrix(tibble(as.double(num_vac), as.double(den_vac)))
      
      IR_compute <- epi.conf(
        data,
        ctype = "inc.rate",
        method = "exact",
        N = tally(pop_vac),
        design = 1,
        conf.level = 0.95
      ) * 100000
      
      IR <- IR |>
        add_row("study_group" = study_group,
                "treatment" = 1,
                "outcome" = outcome,
                "years" = ys,
                "est" = IR_compute$est,
                "lower" = IR_compute$lower,
                "upper" = IR_compute$upper)
      # Unvaccinated IR
      treatment <- 0
      pop_unvac <- cdm$outcome_table |> filter(treatment == 0) |>
        mutate(time_to_event = case_when(
          time_to_event > 365*ys ~ 365*ys,
          time_to_event <= 365*ys ~ time_to_event
        ))
      
      den_unvac <- pop_unvac |> pull(time_to_event) |> sum()/365
      
      num_unvac <- pop_unvac |>
        filter(outcome_name == outcome & censor_reason == "outcome" & time_to_event < 365*ys) |>
        tally() |> pull()
      
      metric <- paste0("IR_",ys,"y_treatment0_", outcome)
      
      #IR[[metric]] <- num_unvac*10^5/(den_unvac)
      
      metric <- paste0("IR_CI_",ys,"y_treatment0_",outcome)
      
      data <- pop_unvac |>
        mutate(outcome_flag = case_when(
          time_to_event <= 365*ys ~ outcome_flag,
          time_to_event > 365*ys ~ 0
        )) |>
        dplyr::select(outcome_flag, time_to_event) |>
        collect() |>
        as.matrix()
      
      data <- as.matrix(tibble(as.double(num_unvac), as.double(den_unvac)))
      
      
      IR_compute <- epi.conf(
        data,
        ctype = "inc.rate",
        method = "exact",
        N = tally(pop_unvac),
        design = 1,
        conf.level = 0.95
      ) * 100000
      
      IR <- IR |>
        add_row("study_group" = study_group,
                "treatment" = 0,
                "outcome" = outcome,
                "years" = ys,
                "est" = IR_compute$est,
                "lower" = IR_compute$lower,
                "upper" = IR_compute$upper)
      
      # Check whether we have outcomes:
      check <- cdm$outcome_table |>
        mutate(outcome_flag = case_when(
          outcome_name == outcome & censor_reason == "outcome" & outcome_date <= cohort_start_date + years(ys) ~ 1,
          .default = 0 )) |>
        filter(outcome_flag == 1) |>
        count() |>
        pull()
      
      if (check > 0) {
      # Poisson Regression
      poisson_model <- glm(formula = outcome_flag ~ treatment, family = "poisson", data = cdm$outcome_table |>
                             mutate(outcome_flag = case_when(
                               outcome_name == outcome & censor_reason == "outcome" & outcome_date <= cohort_start_date + years(ys) ~ 1,
                               .default = 0
                             ))
      )
      poisson_coef <- poisson_model$coefficients[poisson_model$coefficients != "(Intercept)"]
      
      metric <- paste0("IRR_",ys,"y_", outcome)

      IRR <- IRR |>
        add_row("study_cohort" = study_cohort,
                "outcome" = outcome,
                "years" = ys,
                "est" = exp(poisson_coef),
                "lower" = NA,
                "upper" = NA)
      } else {
        metric <- paste0("IRR_",ys,"y_", outcome)
        IRR <- IRR |>
          add_row("study_cohort" = study_cohort,
                  "outcome" = outcome,
                  "years" = ys,
                  "est" = NA,
                  "lower" = NA,
                  "upper" = NA)
      }
    }
    
    print(4)  
    # Survival Study
    info(logger, "SURVIVAL STUDY")
    info(logger, "Kaplan Meier plots")
    
    survival_data <- survfit(Surv(time_to_event, outcome_flag) ~ treatment, data = cdm$outcome_table |>
                               mutate(outcome_flag = case_when(
                                 outcome_name == outcome & censor_reason == "outcome" ~ 1,
                                 .default = 0
                               ))
                             ) 
    
    png(paste0(resultsFolder, "/SurvivalPlot_", study_group, "_", outcome, ".png"))
    plot(survival_data,
         main = "Kaplan Meier Plot",
         xlab = "Time to event (days)", ylab = "Survival probability",
         #xlim = c(17, 37), ylim = c(0.75, 1),
         conf.int = T, mark.time = F,
         lty = 1:2) # 2 groups
    # Add c(1,3,4,2) in two places to re-order the legend to match
    # the ordering of the lines at 37 weeks
    legend(20, 0.40,
           # Inside c(), the ordering is the same as the order shown in surv.ex7.4
           c("Unvac", "Vac")[c(1,2)],
           lty = c(1,2),
           title = "Strata: ")
    dev.off()
    
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
                           outcome_name == outcome & censor_reason == "outcome" ~ 1,
                           .default = 0
                         )))
    summary_coxReg <- summary(cox_model)
    coef_coxReg <- summary_coxReg$coef
    # HRs, 95% CIs, p-values
    metric <- paste0("Cox_", outcome, )
   Cox <- Cox |>
      add_row("study_group" = study_group,
              "outcome" = outcome,
              "coef" = summary_coxReg$coef,
              "HR" = summary_coxReg$exp(coef),
              "lower" = summary_coxReg$lower,
              "upper" = summary_coxReg$upper,
              "p-value" = coef_coxReg[, "Pr(>|z|)"])
    
    
    #info(logger, "Forest plot")
    #sjPlot::plot_model(cox_model)
    print(45)
    #info(logger, "Check proportional Hazards Assumptions")
    # extract dummy variables (0/1 variables and set it upc in the table again): the cox.zph messed up the table
    # cdm$output_table <- cdm$output_table |>
    #   mutate(vac_status = model.matrix(cox_model)[, vac_status])
    # 
    # cox_model_2 <- coxph(Surv(time_to_event, outcome_flag) ~ vac_status, data = cdm$outcome_table)
    #ph_check <- cox.zph(cox_model) # p > 0.05
    
    #print(48)
    # provar random codi: flu
  }
  
  write.csv(unlist(IR), file = paste0(resultsFolder, "/IR_", study_group, ".csv"))
  write.csv(unlist(IRR), file = paste0(resultsFolder, "/IRR_", study_group, ".csv"))
  write.csv(unlist(Cox), file = paste0(resultsFolder, "/Cox_", study_group, ".csv"))
}
