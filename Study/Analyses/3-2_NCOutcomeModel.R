library(survival)
library(epiR)
library(EmpiricalCalibration)
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
  cdm$unvac_18f_coverage_0_any_dose |> mutate(treatment = 0, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(18)) |> dplyr::select(!date_of_birth),
  cdm$vac_18f_coverage_0_any_dose |> mutate(treatment = 1, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(18)) |> dplyr::select(!date_of_birth),
  # cdm$matched_unvaccinated_coverage_0_18f_any_dose |> mutate(treatment = 0, group = "matched_overall"),
  # cdm$matched_vaccinated_coverage_0_18f_any_dose |> mutate(treatment = 1, group = "matched_overall"),
  cdm$matched_unvac_cov_0_ratio_18f_any_dose |> mutate(treatment = 0, group = "matched_overall"),
  cdm$matched_vac_cov_0_ratio_18f_any_dose |> mutate(treatment = 1, group = "matched_overall"),
  cdm$not_matched_vaccinated_18f_coverage_0_1dose |> mutate(treatment = 0, group = "initial_1vs23"),
  cdm$not_matched_vaccinated_18f_coverage_0_23dose |> mutate(treatment = 1, group = "initial_1vs23"),
  cdm$not_matched_vaccinated_18f_coverage_0_2dose |> mutate(treatment = 0, group = "initial_2vs3"),
  cdm$not_matched_vaccinated_18f_coverage_0_3dose |> mutate(treatment = 1, group = "initial_2vs3"),
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
  cdm$not_matched_vaccinated_15m_coverage_0_1dose |> mutate(treatment = 0, group = "initial_1vs23"),
  cdm$not_matched_vaccinated_15m_coverage_0_23dose |> mutate(treatment = 1, group = "initial_1vs23"),
  cdm$not_matched_vaccinated_15m_coverage_0_2dose |> mutate(treatment = 0, group = "initial_2vs3"),
  cdm$not_matched_vaccinated_15m_coverage_0_3dose |> mutate(treatment = 1, group = "initial_2vs3"),
  cdm$matched_vaccinated_15m_coverage_0_1dose |> mutate(treatment = 0, group = "matched_1vs23"),
  cdm$matched_vaccinated_15m_coverage_0_23dose |> mutate(treatment = 1, group = "matched_1vs23"),
  cdm$matched_vaccinated_15m_coverage_0_2dose |> mutate(treatment = 0, group = "matched_2vs3"),
  cdm$matched_vaccinated_15m_coverage_0_3dose |> mutate(treatment = 1, group = "matched_2vs3"),
  name = "hpv_population_cohorts"
)

info(logger, "PREPARE OUTCOME TABLE")
# Negative Outcome Cohorts ################### ---------------
# Read the CSV file
#nco <- read.csv("Analyses/NCO.csv", header = TRUE) # Adjust the file path as needed

sex <- "female"
resultsFolder <- paste0("Results_NCO_",sex)
resultsFolder <- "Darwin_comparison_NCO"
outcome_list <- settings(cdm$nondarwin_nco_cohorts) |> pull(cohort_name)
outcome_names_list <- settings(cdm$nondarwin_nco_cohorts) |> dplyr::select(cohort_definition_id, cohort_name) |> rename(outcome_name = cohort_name)

cdm$outcome_cohorts <- cdm$nondarwin_nco_cohorts |>
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
# --------

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

# info(logger, "Prepare dose 1vs23 outcome table")
# # DOSE stratas: 1 vs 23 ############ ---------
# study_group <- "dose1vs23"
# 
# cdm$total_matched_cohort <- cdm$hpv_population_cohorts |>
#   filter(group == study_group) |>
#   compute(name = "total_matched_cohort")
# 
# print(0)
# # Create next dose
# next_dose <- cdm$total_matched_cohort |>
#   filter(treatment == 0) |>
#   left_join(cdm$doses_allvac_cohort |>
#               dplyr::select(subject_id, cohort_start_date) |>
#               rename(next_dose_date = cohort_start_date), by = "subject_id") |>
#   group_by(subject_id) |>
#   filter(! next_dose_date == min(next_dose_date)) |>
#   filter(next_dose_date == min(next_dose_date)  | is.na(next_dose_date)) |>
#   ungroup() |>
#   dplyr::select(subject_id, next_dose_date)
# 
# # Add next dose to outcome table
# cdm$partial_outcome_table <- cdm$total_matched_cohort |>
#   left_join(next_dose, by = "subject_id") |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add death date
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   left_join(cdm$death |>
#               dplyr::select(person_id, death_date) |>
#               rename(subject_id = person_id)) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add outcome column
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   left_join(cdm$outcome_cohorts |>
#               dplyr::select(subject_id, cohort_start_date, outcome_name) |>
#               rename(outcome_date = cohort_start_date), by = "subject_id") |>
#   filter(outcome_date >= cohort_start_date) |>
#   mutate(outcome_flag = case_when(
#     is.na(outcome_date) ~ 0,
#     !is.na(outcome_date) ~ 1
#   )) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add event column: death, end obs, outcome, next dose
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
#   compute(name = "partial_outcome_table") |>
#   mutate(event_type = case_when(
#     event_date == death_date ~ "death",
#     event_date == cohort_end_date ~ "end of observation",
#     event_date == outcome_date ~ "outcome",
#     event_date == next_dose_date ~ "next dose"
#   )) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add Censor pair column
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   group_by(pair_id) |>
#   mutate(censored_date = min(event_date, na.rm = FALSE)) |>
#   ungroup() |>
#   compute(name = "partial_outcome_table") |>
#   mutate(censor_reason = case_when(
#     censored_date == event_date ~ event_type,
#     censored_date != event_date ~ "pair censored"
#   )) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add time to event column
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   mutate(time_to_event = censored_date - cohort_start_date) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Save for this group
# cdm$full_outcome_table <- cdm$full_outcome_table |>
#   union_all(cdm$partial_outcome_table) |>
#   compute(name = "full_outcome_table", temporary = FALSE)
# 
# # -------------------
# 
# info(logger, "Prepare dose 2vs3 outcome table")
# # DOSE stratas: 2 vs 3 ############ ---------------
# study_group <- "dose2vs3"
# 
# cdm$total_matched_cohort <- cdm$hpv_population_cohorts |>
#   filter(group == study_group) |>
#   compute(name = "total_matched_cohort") |>
#   dplyr::select(! c(treatment, index_date))
# 
# print(0)
# # Create next dose
# next_dose <- cdm$total_matched_cohort |>
#   filter(treatment == 0) |>
#   left_join(cdm$doses_allvac_cohort |>
#               dplyr::select(subject_id, cohort_start_date) |>
#               rename(next_dose_date = cohort_start_date), by = "subject_id") |>
#   group_by(subject_id) |>
#   filter(! next_dose_date == min(next_dose_date)) |>
#   filter(! next_dose_date == min(next_dose_date)) |>
#   ungroup() |>
#   dplyr::select(subject_id, next_dose_date)
# 
# # Add next dose to outcome table
# cdm$partial_outcome_table <- cdm$total_matched_cohort |>
#   left_join(next_dose, by = "subject_id") |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add death date
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   left_join(cdm$death |>
#               dplyr::select(person_id, death_date) |>
#               rename(subject_id = person_id)) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add outcome column
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   left_join(cdm$outcome_cohorts |>
#               dplyr::select(subject_id, cohort_start_date, outcome_name) |>
#               rename(outcome_date = cohort_start_date), by = "subject_id") |>
#   filter(outcome_date >= cohort_start_date) |>
#   mutate(outcome_flag = case_when(
#     is.na(outcome_date) ~ 0,
#     !is.na(outcome_date) ~ 1
#   )) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add event column: death, end obs, outcome, next dose
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
#   mutate(event_type = case_when(
#     event_date == death_date ~ "death",
#     event_date == cohort_end_date ~ "end of observation",
#     event_date == outcome_date ~ "outcome",
#     event_date == next_dose_date ~ "next dose"
#   )) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add Censor pair column
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   group_by(pair_id) |>
#   mutate(censored_date = min(event_date, na.rm = TRUE)) |>
#   ungroup() |>
#   mutate(censor_reason = case_when(
#     censored_date == event_date ~ event_type,
#     censored_date != event_date ~ "pair censored"
#   )) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Add time to event column
# cdm$partial_outcome_table <- cdm$partial_outcome_table |>
#   mutate(time_to_event = censored_date - cohort_start_date) |>
#   compute(name = "partial_outcome_table", temporary = FALSE)
# 
# # Save for this group
# cdm$full_outcome_table <- cdm$full_outcome_table |>
#   union_all(cdm$partial_outcome_table) |>
#   compute(name = "full_outcome_table", temporary = FALSE)
# 
# 
# # -----------------------
# 
# info(logger, "Prepare year of vaccination strata outcome table")
# # YEAR OF VAC stratas ############ --------
# cdm$total_matched_cohort <- cdm$full_outcome_table |>
#   filter(group == "matched_overall") |>
#   mutate(cohort_year = year(cohort_start_date)) |> 
#   compute(name = "total_matched_cohort", temporary = FALSE)
# 
# for (ys in c(2008:2023)) {
#   print(ys)
#   cdm$partial_outcome_table <- cdm$total_matched_cohort |>
#     filter(cohort_year == ys) |>
#     mutate(group = paste0("yov",ys)) |>
#     compute(name = "partial_outcome_table", temporary = FALSE)
#   
#   cdm$full_outcome_table <- cdm$full_outcome_table |>
#     union_all(cdm$partial_outcome_table) |>
#     compute(name = "full_outcome_table", temporary = FALSE)
# }
# -----------------------

# COMPARISON PLUS AND MINUS 15 YEARS FOR MATCHED INDIVIDUALS

cdm$hpv_population_cohorts_split <- cdm$hpv_population_cohorts |>
  addAge(ageName = "age_at_vac") |>
  mutate(group = case_when(
    age_at_vac <= 15 ~ paste0(group,"_atv_0_15"),
    age_at_vac > 15 ~ paste0(group,"_atv_15_18"))
  )|>
  dplyr::select(- age_at_vac) |>
  compute(name = "hpv_population_cohorts_split", temporary = FALSE)


study_groups <- c("matched_overall_atv_0_15", "matched_overall_atv_15_18", "matched_1vs23_atv_0_15", "matched_1vs23_atv_15_18", "matched_2vs3_atv_0_15", "matched_2vs3_atv_15_18")

for (study_group in study_groups) {
  cdm$total_matched_cohort <- cdm$hpv_population_cohorts_split |>
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


# COMPUTE METRICS ---------------
info(logger, "COMPUTE METRICS: IR, IRR (Poisson Regression), Cox, Survival")

groups <- cdm$full_outcome_table |>
  distinct(group) |>
  pull()

groups <- c("matched_overall_atv_0_15", "matched_overall_atv_15_18", "matched_1vs23_atv_0_15", "matched_1vs23_atv_15_18", "matched_2vs3_atv_0_15", "matched_2vs3_atv_15_18")

# pull returns a vector
for (study_group in groups) {
  cdm$preoutcome_table <- cdm$full_outcome_table |>
    filter(group == study_group) |>
    compute(name = "preoutcome_table", temporary = FALSE)
  
  # PREPARE METRICS
  IR <- tibble()
  IRR <- tibble()
  Cox <- tibble()
  print(1)
  
  for (outcome in outcome_list) {
    cdm$outcome_table <- cdm$preoutcome_table |>
      dplyr::select(c(cohort_definition_id, group, subject_id, cohort_start_date, cohort_end_date, cohort_year, treatment, pair_id, censored_date, censor_reason, paste0(outcome, "_0_to_inf"))) |>
      mutate(time_to_event = pmin(censored_date, !!sym(paste0(outcome, "_0_to_inf")), na.rm = TRUE) - as.Date(cohort_start_date)) |>
      mutate(event_type = case_when(
        censored_date < !!sym(paste0(outcome, "_0_to_inf")) | is.na(!!sym(paste0(outcome, "_0_to_inf"))) ~ censor_reason,
        .default = outcome
      )) |>
      compute(name = "outcome_table", temporary = FALSE)
    
    print(2)
    for (ys in c(5, 10, 15, 33)) {
      print(3)
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
      
      metric <- paste0("IR_",ys,"y_treatment1_",outcome)
      
      ir1 <- as.numeric(num_vac*10^5/(den_vac))
      # IR[[metric]] <- ir1
      
      
      # metric <- paste0("IR_CI_",ys,"y_treatment1_",outcome)
      
      # data <- pop_vac |>
      #   mutate(outcome_flag = case_when(
      #     time_to_event <= 365.25*ys ~ outcome_flag,
      #     time_to_event > 365.25*ys ~ 0
      #   )) |>
      #   dplyr::select(outcome_flag, time_to_event) |>
      #   collect() |>
      #   as.matrix()
      
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
      
      metric <- paste0("IR_",ys,"y_treatment0_", outcome)
      
      ir0 <- num_unvac*10^5/(den_unvac)
      # IR[[metric]] <- ir0
      
      metric <- paste0("IR_CI_",ys,"y_treatment0_",outcome)
      
      # data <- pop_unvac |>
      #   mutate(outcome_flag = case_when(
      #     time_to_event < 365.25*ys ~ outcome_flag,
      #     time_to_event >= 365.25*ys ~ 0
      #   )) |>
      #   dplyr::select(outcome_flag, time_to_event) |>
      #   collect() |>
      #   as.matrix()
      
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
      # if (ir1 != 0 & ir0 != 0) {
      #   OR <- confint(poisson_model)
      #   metric <- paste0("CI_lower_IRR_",ys,"y_", outcome)
      #   IRR[[metric]] <- exp(OR[2])
      #   metric <- paste0("CI_higher_IRR_",ys,"y_", outcome)
      #   IRR[[metric]] <- exp(OR[4])
      #   
      # }
      # 
    }
    
    print(4)  
    # Survival Study
    info(logger, "SURVIVAL STUDY")
    info(logger, "Kaplan Meier plots")
    
    survival_data <- survfit(Surv(time_to_event, outcome_flag) ~ treatment, data = cdm$outcome_table |>
                               mutate(outcome_flag = case_when(
                                 event_type == outcome ~ 1,
                                 .default = 0
                               ))) 
    
    png(paste0(resultsFolder, "/SurvivalPlot_", study_group, "_", outcome, ".png"))
    plot(survival_data,
         main = "Kaplan Meier Plot",
         xlab = "Time to event (days)", ylab = "Survival probability",
         #xlim = c(17, 37), ylim = c(0.75, 1),
         ylim = c(0.90, 1),
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
                           event_type == outcome ~ 1,
                           .default = 0
                         )))
    # summary_coxReg <- summary(cox_model)
    # coef_coxReg <- summary_coxReg$coef
    # HRs, 95% CIs, p-values
    # metric <- paste0("Cox_", outcome)
    # Cox[[metric]] <- tibble("HR" = exp(coef_coxReg[, "coef"]),
    #                        exp(confint(cox_model)),
    #                        "p-value" = coef_coxReg[, "Pr(>|z|)"])
    Cox <- Cox |> 
      rbind(broom::tidy(cox_model) |> 
              mutate(outcome = outcome, hr = exp(estimate), lower = exp(estimate - 1.96*std.error), upper = exp(estimate + 1.96*std.error))
      )
    
    
    #info(logger, "Forest plot")
    #sjPlot::plot_model(cox_model)
    print(45)
    #info(logger, "Check proportional Hazards Assumptions")
    # extract dummy variables (0/1 variables and set it upc in the table again): the cox.zph messed up the table
    # cdm$output_table <- cdm$output_table |>
    #   mutate(treatment = model.matrix(cox_model)[, treatment])
    # 
    # cox_model_2 <- coxph(Surv(time_to_event, outcome_flag) ~ treatment, data = cdm$outcome_table)
    #ph_check <- cox.zph(cox_model) # p > 0.05
    
    #print(48)
    # provar random codi: flu
  }
  
  write.csv(IR, file = paste0(resultsFolder, "/IR_", study_group, ".csv"), row.names = FALSE)
  write.csv(IRR, file = paste0(resultsFolder, "/IRR_", study_group, ".csv"), row.names = FALSE)
  write.csv(Cox, file = paste0(resultsFolder, "/Cox_", study_group, ".csv"), row.names = FALSE)
}
# Plot negative control effect sizes 
plotForest(Cox$estimate, Cox$std.error, names = Cox$outcome)
# Fit the null distribution:  estimate the systematic error distribution
null <- fitNull(Cox$estimate, Cox$std.error)
null # SD < 25% indicates small  variability in the systematic error from one estimate to the next // mean > 0 indicates analysis positively biased
# Is the estimation of the systematic error distribution a good one?
# Test whether the calibrated p-value is truly calibrated, meaning the fraction of negative controls with a p-value below alpha is approximately the same as alpha
plotCalibration(Cox$estimate, Cox$std.error) # Calibrated p-value closer to the diagonal
# Plot the null distribution
plotCalibrationEffect(Cox$estimate, Cox$std.error, title = "NCO")
# Outcomes of interest
interest <- tibble(outcome = c("warts", "cin1"), estimate = c(0.118797844966718, 0.30989572445775), se = c(0.0661892435234304, 0.174767140651393))
p <- calibrateP(null, interest$estimate, interest$se) # if p < 0.05 --> significancy
model <- fitSystematicErrorModel(Cox$estimate, Cox$std.error, rep(0,40), estimateCovarianceMatrix = TRUE)
calibrated <- calibrateConfidenceInterval(interest$estimate, interest$se, model, ciWidth = 0.95)

plotCalibrationEffect(Cox$estimate,
                      Cox$std.error,
                      interest$estimate,
                      interest$se,
                      null, title = "Calibration effect")

#--------------------

