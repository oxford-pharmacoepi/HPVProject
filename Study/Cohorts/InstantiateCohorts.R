library(CDMConnector)
library(dplyr)
library(dbplyr)
library(here)
library(DBI)
library(duckdb)
library(tibble)
library(SqlRender)
library(DatabaseConnector)
library(remotes)
library(PatientProfiles)
library(CodelistGenerator)
library(tidyverse)
library(usethis)
library(RPostgres)

# Load vaccination records
if (instantiateVaccinationCohorts) {
  cohorts <- readCohortSet(path = here("Cohorts", "HPV_allvac"))
  cdm <- generateCohortSet(cdm = cdm, cohortSet = cohorts, name = c("doses_allvac_cohort"))
  
  cdm$firstdose_cohort <- cdm$doses_allvac_cohort |>
    group_by(subject_id) |>
    filter(cohort_start_date == min(cohort_start_date)) |>
    ungroup() |>
    compute(name = "firstdose_cohort", temporary = FALSE)
} else {
  cdm <- cdmFromCon(
    con = db,
    cdmSchema = cdmSchema, 
    writeSchema = writeSchema, 
    cdmName = server_dbi, 
    achillesSchema = achillesSchema,
    cohortTables = c("firstdose_cohort",
                     "doses_allvac_cohort")
  )
}

if (instantiateCohorts) {
  info(logger, "INSTANTANIATE COHORTS")
  
  # GENERAL COHORT RESTRICTING: sex, year of birth and time in observation ------------
  cdm$condition_cohort <- cdm$observation_period %>%
    # All individuals
    mutate(cohort_definition_id = 0) %>% # Add cohort id column
    rename(cohort_start_date = observation_period_start_date, cohort_end_date = observation_period_end_date, subject_id = person_id) %>% # rename according to the cohort vocabulary
    compute(name = "condition_cohort", temporary = FALSE) %>% # Save the result permanently
    newCohortTable(cohortSetRef = tibble(cohort_definition_id = 0, cohort_name = "condition_cohort")) %>%
    dplyr::select(!c(observation_period_id, period_type_concept_id)) %>% # eliminate unnecessary columns
    addSex() %>%
    filter(sex == sex_cohort) %>% # Filter to select sex
    compute(name = "condition_cohort", temporary = FALSE) %>%
    recordCohortAttrition(paste0("Restrict to ", sex_cohort)) %>%  # Record the changes made
    addDateOfBirth(cdm$observation_period) %>%
    addAge() %>%
    filter(year(date_of_birth) >= start_birth) %>%
    compute(name = "condition_cohort", temporary = FALSE) %>%
    recordCohortAttrition("Restrict to subjects born in 1990 or after") %>%
    mutate(!!paste0("date_",yy,"years") := as.Date(date_of_birth + years(yy))) %>%
    filter(cohort_end_date >= !!sym(paste0("date_",yy,"years"))) %>%
    compute(name = "condition_cohort", temporary = FALSE) %>%
    recordCohortAttrition("Restrict to subjects turning 18 in observation") |>
    mutate(date_9years = as.Date(date_of_birth + years(9))) %>%
    filter(cohort_start_date <= date_9years) %>%
    compute(name = "condition_cohort", temporary = FALSE) %>%
    recordCohortAttrition("Restrict to subjects in observation since 9yo")
  
  # Intersect the conditioned cohort to the cohort containing all vaccinated people  
  cdm$vac_status_cohort <- cdm$condition_cohort %>%
    addCohortIntersectFlag(targetCohortTable = "firstdose_cohort", window = c(-Inf, 0), indexDate = paste0("date_",yy,"years"), nameStyle = "vac_status") %>% #Intersection conditioned to vaccination before yyth birthday
    mutate(cohort_definition_id = vac_status + 1) %>% # Separate unvac and vac people into 2 cohorts (cdi = 1, 2)
    compute(name = "vac_status_cohort", temporary = FALSE)
  
  # Define unvaccinated cohort
  cdm[[paste0("unvac_",yy,ss,"_coverage_0_any_dose")]] <- cdm$vac_status_cohort %>%
    filter(cohort_definition_id == 1) %>%
    mutate(cohort_start_date = !!sym(paste0("date_",yy,"years"))) %>%
    compute(name = paste0("unvac_",yy,ss,"_coverage_0_any_dose"), temporary = FALSE) %>%
    newCohortTable(cohortSetRef = tibble(cohort_definition_id = 1, cohort_name = paste0("unvac_",yy,ss,"_coverage_0_any_dose")),
                   cohortAttritionRef = attrition(cdm$vac_status_cohort) |> mutate(cohort_definition_id = 1)) |>
    recordCohortAttrition(paste0("Restrict to unvaccinated at ",yy," yo"))
  
  # Define vaccinated cohort
  cdm[[paste0("vac_",yy,ss,"_coverage_0_any_dose")]] <- cdm$vac_status_cohort %>%
    filter(cohort_definition_id == 2) %>%
    mutate(cohort_start_date = !!sym(paste0("date_",yy,"years"))) %>%
    compute(name = paste0("vac_",yy,ss,"_coverage_0_any_dose"), temporary = FALSE) %>%
    newCohortTable(cohortSetRef = tibble(cohort_definition_id = 2, cohort_name = paste0("vac_",yy,ss,"_coverage_0_any_dose")),
                   cohortAttritionRef = attrition(cdm$vac_status_cohort) |> mutate(cohort_definition_id = 2)) |>
    recordCohortAttrition(paste0("Restrict to vaccinated at ",yy," yo"))
  # --------------------------------
  
  
} else {
  
  cdm <- cdmFromCon(
    con = db,
    cdmSchema = cdmSchema, 
    writeSchema = writeSchema, 
    cdmName = server_dbi, 
    achillesSchema = achillesSchema,
    cohortTables = c("firstdose_cohort",
                     "doses_allvac_cohort",
                     paste0("unvac_",yy,ss,"_coverage_0_any_dose"),
                     paste0("vac_",yy,ss,"_coverage_0_any_dose")
                     )
    )
  
}



