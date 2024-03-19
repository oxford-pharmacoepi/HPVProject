cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = dbName, 
  achillesSchema = achillesSchema, 
  cohortTables = c("vac_cohort", "unvac_cohort", "allvac_cohort", "total_vac_matched_cohort", "total_unvac_matched_cohort", "doses_allvac_cohort", "phenotyping_hpv_")
)

info(logger, "Prepare Outcome table")
# Instantiate all doses of vaccine
cohorts <- readCohortSet(path = here("Cohorts", "HIV_allvac"))
cdm <- generateCohortSet(cdm = cdm, cohortSet = cohorts, name = c("doses_allvac_cohort"))
# Put together matched cohorts
cdm$total_matched_cohort <- union_all(cdm$total_vac_matched_cohort, cdm$total_unvac_matched_cohort) |>
  compute(name = "total_matched_cohort", temporary = FALSE) |>
  select(! c(region, subclass, cohort_year, index_date)) 

# Create next dose
next_dose <- cdm$total_matched_cohort |>
  left_join(cdm$doses_allvac_cohort |>
               select(subject_id, cohort_start_date) |>
               rename(next_dose_date = cohort_start_date), by = "subject_id") |>
  group_by(subject_id) |>
  filter(next_dose_date > cohort_start_date | is.na(next_dose_date)) |>
  filter(next_dose_date == min(next_dose_date) | is.na(next_dose_date)) |>
  ungroup() |>
  select(subject_id, next_dose_date)

# Add next dose to outcome table
cdm$outcome_table <- cdm$total_matched_cohort |>
  left_join(next_dose, by = "subject_id") |>
  compute(name = "outcome_table", temporary = FALSE)

# Add death date
cdm$outcome_table <- cdm$outcome_table |>
  left_join(cdm$death |>
              select(person_id, death_date) |>
              rename(subject_id = person_id)) |>
  compute(name = "outcome_table", temporary = FALSE)

# Add outcome column
cdm$outcome_table <- cdm$outcome_table |>
  left_join(cdm$phenotyping_hpv_ |>
              select(subject_id, cohort_start_date) |>
              rename(outcome_date = cohort_start_date), by = "subject_id") |>
  mutate(outcome_flag = case_when(
    is.na(outcome_date) ~ 0,
    !is.na(outcome_date) ~ 1
    )) |>
  compute(name = "outcome_table", temporary = FALSE)

# Add event column: death, end obs, outcome, next dose
cdm$outcome_table <- cdm$outcome_table |>
  mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
  compute(name = "outcome_table", temporary = FALSE) |>
  mutate(event_type = case_when(
    event_date == death_date ~ "death",
    event_date == cohort_end_date ~ "end of observation",
    event_date == outcome_date ~ "outcome",
    event_date == next_dose_date ~ "next dose"
    )
  )

# Add Censor pair column
cdm$outcome_table <- cdm$outcome_table |>
  group_by(pair_id) |>
  mutate(censored_date = min(event_date, na.rm = FALSE)) |>
  ungroup() |>
  compute(name = "outcome_table", temporary = FALSE) |>
  mutate(censor_reason = case_when(
    censored_date == event_date ~ event_type,
    censored_date != event_date ~ "pair censored"
  ))

# Add time to event column
cdm$outcome_table <- cdm$outcome_table |>
  mutate(time_to_event = censored_date - cohort_start_date) |>
  compute(name = "outcome_table", temporary = FALSE)
