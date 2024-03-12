library(tidyr)
library(dplyr)
library(dbplyr)
library(CDMConnector)
library(RPostgres)
library(PatientProfiles)
library(log4r)
library(visOmopResults)

# Details
log_file        <- here(resultsFolder, "log_CharacterisationData.txt")
logger          <- create.logger(logfile = log_file, level = "INFO")

# Prepare cohorts
cdm$vac_cohort <- cdm$vac_cohort |> 
  mutate(cohort_start_date = as.Date(date_15years)) |>
  select(! date_15years) |>
  compute(name = "vac_cohort", temporary = FALSE)

cdm$unvac_cohort <- cdm$unvac_cohort |>
  mutate(cohort_start_date = as.Date(date_15years)) |>
  select(! date_15years) |>
  compute(name = "unvac_cohort", temporary = FALSE)

# TABLE 1 CHARACTERISATION
info(logger, "TABLE ONE SUMMARY")
cdm <- omopgenerics::bind(
  cdm$vac_cohort |> select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"), 
  cdm$unvac_cohort |> select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  cdm$total_unvac_matched_cohort |> select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  cdm$total_vac_matched_cohort |> select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  name = "hpv_cohorts"
)

table_one_summary <- cdm$hpv_cohorts |> 
  select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addAge(cdm,ageGroup = list(c(0,15),c(16,19),c(20,Inf))) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
    strata = list("age_group"),
    ageGroup = list(
      c(0, 15), c(16, 19), c(20, Inf)
    ),
    cohortIntersect = list(
      "Medications" = list(
        targetCohortTable = "medications", value = "flag", window = list(c(-365,-1),c(-30,-1))
      ),
      "Conditions" = list(
        targetCohortTable = "conditions", value = "flag", window = list(c(-365,-1),c(-30,-1),c(-Inf,-1))
      ),
      "HIV_status" = list(
        targetCohortTable = "hiv_status", value = "flag", window = c(-Inf,-1)
      ),
      "vaccinations" = list(
        targetCohortTable = "vaccinations", value = "count", window = c(-Inf,-1)
      ),
      "papanicolau_smear_testing" = list(
        targetCohortTable = "papanicolau_smear_testing", value = c("count"), window = c(-Inf,-1)
      )
    ),
    #minCellCount = 5
  )

write.csv(table_one_summary,
          paste0(resultsFolder,"/tableOne_summary_",cdmSchema,".csv"), row.names = FALSE)

info(logger, "Tidy table one")
tableOne_tidy <- table_one_summary |> 
  as_tibble() |>
  splitAll() |>
  pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  select(! c(cdm_name, result_type, package_name, package_version))

write.csv(tableOne_tidy, paste0(resultsFolder,"/tableOne_tidy_",cdmSchema,".csv"), row.names = FALSE)
# Large Scale Characterisations
# We want to make characterisation at index date (= vac date = cohort_start_date)
info(logger, "Compute SMD for TableOne")
tableOne_red <- tableOne_tidy |> 
  filter(estimate_name == "percentage") |>
  mutate_all(~replace(., is.na(.), 0)) |>
  mutate(across(c(vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort), ~ as.numeric(.x)/100))

tableOne_smd <- tableOne_red |>
  mutate(original_smd = abs((vac_cohort - unvac_cohort)/sqrt((vac_cohort*(1-vac_cohort) + unvac_cohort*(1-unvac_cohort))/2))) |>
  mutate(matched_smd = abs((total_vac_matched_cohort - total_unvac_matched_cohort)/sqrt((total_vac_matched_cohort*(1-total_vac_matched_cohort) + total_unvac_matched_cohort*(1-total_unvac_matched_cohort))/2))) |>
  mutate_all(~replace(., is.na(.), 0))

write.csv(tableOne_smd, paste0(resultsFolder,"/tableOne_SMD_",cdmSchema,".csv"), row.names = FALSE)
# Vaccinated cohort
info(logger, "LARGE CHARACTERISATION SUMMARY")
slc <- summariseLargeScaleCharacteristics(cdm$hpv_cohorts, window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), eventInWindow = c("condition_occurrence","drug_exposure"))

write.csv(slc,
          paste0(resultsFolder,"/largeScale_summary_",cdmSchema,".csv"), row.names = FALSE)

info(logger, "Tidy Large Scale Characteristics")
slc_tidy <- slc |>
  as_tibble() |>
  splitAll() |>
  pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  select(! c(cdm_name, result_type, package_name, package_version))

write.csv(slc_tidy, paste0(resultsFolder,"/largeScale_tidy_",cdmSchema,".csv"), row.names = FALSE)

info(logger, "Compute SMD for LSC")
slc_red <- slc_tidy |> 
  filter(estimate_name == "percentage") |>
  mutate_all(~replace(., is.na(.), 0)) |>
  mutate(across(c(vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort), ~ as.numeric(.x)/100))

slc_smd <- slc_red |>
  mutate(original_smd = abs((vac_cohort - unvac_cohort)/sqrt((vac_cohort*(1-vac_cohort) + unvac_cohort*(1-unvac_cohort))/2))) |>
  mutate(matched_smd = abs((total_vac_matched_cohort - total_unvac_matched_cohort)/sqrt((total_vac_matched_cohort*(1-total_vac_matched_cohort) + total_unvac_matched_cohort*(1-total_unvac_matched_cohort))/2))) |>
  mutate_all(~replace(., is.na(.), 0))

write.csv(slc_smd, paste0(resultsFolder,"/largeScale_SMD_",cdmSchema,".csv"), row.names = FALSE)

