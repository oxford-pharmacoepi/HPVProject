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
cohorts <- readCohortSet(path = here("HPV_allvac_cohort.json"))
cdm <- generateCohortSet(cdm = cdm, cohortSet = cohorts, name = c("allvac_cohort"))


cdm$condition_cohort <- cdm$observation_period %>%
  # All individuals
  mutate(cohort_definition_id = 0) %>% # Add cohort id column
  rename(cohort_start_date = observation_period_start_date, cohort_end_date = observation_period_end_date, subject_id = person_id) %>% # rename according to the cohort vocabulary
  select(!c(observation_period_id, period_type_concept_id)) %>% # eliminate unnecessary columns
  compute(name = "condition_cohort", temporary = FALSE) %>% # Save the result permanently
  newCohortTable(cohortSetRef = tibble(cohort_definition_id = 0, cohort_name = "condition_cohort")) %>%
  addSex() %>%
  filter(sex == "Female") %>% # Filter to select only women
  compute(name = "condition_cohort", temporary = FALSE) %>%
  recordCohortAttrition("Select females") %>%  # Record the changes made
  addDateOfBirth(cdm$observation_period) %>%
  mutate(date_15years = date_of_birth + years(15)) %>%
  filter(cohort_start_date <= date_15years, cohort_end_date >= date_15years) %>%
  compute(name = "condition_cohort", temporary = FALSE) %>%
  recordCohortAttrition("Select subjects turning 15 in observation") %>%
  filter(date_of_birth >= as.POSIXct("2007-01-01") - years(15), date_of_birth <= as.POSIXct("2023-12-31") - years(15)) %>%
  compute(name = "condition_cohort", temporary = FALSE) %>%
  recordCohortAttrition("Select subjects turning 15 between 01/01/2007 and 31/12/2023") %>%
  addPriorObservation(indexDate = "date_15years", priorObservationName = "prior_observation_15y") %>%
  filter(prior_observation_15y >= 365) %>%
  compute(name = "condition_cohort", temporary = FALSE) %>%
  recordCohortAttrition("Select subjects with at least 365 days of prior observation") 

# Intersect the conditioned cohort to the cohort containing all vaccinated people  
cdm$vac_status_cohort <- cdm$condition_cohort %>%
  addCohortIntersectFlag(targetCohortTable = "allvac_cohort", window = c(-Inf, 0), indexDate = "date_15years", nameStyle = "vac_status") %>% #Intersection conditioned to vaccination before 15th birthday
  mutate(cohort_definition_id = vac_status + 1) %>% # Separate unvac and vac people into 2 cohorts (cdi = 1, 2)
  compute(name = "vac_status_cohort", temporary = FALSE)

# Define unvaccinated cohort
cdm$unvac_cohort <- cdm$vac_status_cohort %>%
  filter(cohort_definition_id == 1) %>%
  compute(name = "unvac_cohort", temporary = FALSE) %>%
  newCohortTable(
    cohortSetRef = tibble(cohort_definition_id = c(1), cohort_name = c("unvac_cohort")),
    cohortAttritionRef = attrition(cdm$vac_status_cohort) %>% mutate(cohort_definition_id = 1)
  ) %>% 
  recordCohortAttrition("Select not vaccinated")

# Define vaccinated cohort
cdm$vac_cohort <- cdm$vac_status_cohort %>%
  filter(cohort_definition_id == 2) %>%
  compute(name = "vac_cohort", temporary = FALSE) %>%
  newCohortTable(
    cohortSetRef = tibble(cohort_definition_id = c(2), cohort_name = c("vac_cohort")),
    cohortAttritionRef = attrition(cdm$condition_cohort) %>% mutate(cohort_definition_id = 2)
  ) %>% 
  recordCohortAttrition("Select vaccinated")

#listTables(con = db, schema = "results")
#cohortSet(cdm$vac_cohort) to see the name of the cohort
