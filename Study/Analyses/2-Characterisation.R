library(tidyr)
library(dplyr)
library(dbplyr)
library(CDMConnector)
library(RPostgres)
library(PatientProfiles)
library(log4r)
library(visOmopResults)
library(data.table)

info(logger, "CHARACTERISATION")

# TABLE 1 CHARACTERISATION
info(logger, "TABLE ONE SUMMARY")
cdm <- omopgenerics::bind(
  cdm$vac_cohort |> mutate(cohort_start_date = date_15years) |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"), 
  cdm$unvac_cohort |> mutate(cohort_start_date = date_15years) |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  cdm$doses1_matched_cohort |> mutate(cohort_start_date = index_date) |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  cdm$doses2_matched_cohort |> mutate(cohort_start_date = index_date) |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  #cdm$total_unvac_matched_cohort |> mutate(cohort_start_date = index_date) |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  #cdm$total_vac_matched_cohort |> mutate(cohort_start_date = index_date)  |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date"),
  name = "hpv_cohorts"
)


table_one_summary <- cdm$hpv_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
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

write.csv(table_one_summary |> suppress(minCellCount = 5),
          paste0(resultsFolder,"/tableOne_summary_",cdmSchema,".csv"), row.names = FALSE)

info(logger, "Tidy table one")
tableOne_tidy <- table_one_summary |> 
  #as_tibble() |>
  splitAll() |>
  dplyr::mutate(across(c(age_group, window), ~replace(., is.na(.), "overall"))) |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_type, package_name, package_version)) |>
  dplyr::mutate(across(! c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0)))

tableOne_tidy_def <- tableOne_tidy |>
  filter(estimate_name == "count") |>
  mutate(across(! c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table), ~ as.numeric(.x))) |>
  mutate(across(! c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table), ~ replace(., . < 5, NA))) |>
  mutate(across(! c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table), ~ as.character(.x))) |>
  union_all(tableOne_tidy |>
              filter(estimate_name != "count"))

write.csv(tableOne_tidy_def, paste0(resultsFolder,"/tableOne_tidy_",cdmSchema,".csv"), row.names = FALSE)

# Large Scale Characterisations
# We want to make characterisation at index date (= vac date = cohort_start_date)
info(logger, "Compute SMD for TableOne")

tableOne_red <- tableOne_tidy |> 
  filter(estimate_name == "percentage") |>
  mutate(across(! c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table, window), ~ as.numeric(.x)/100))
  
tableOne_counts <- tableOne_tidy |>
  filter(estimate_name == "count") |>
    dplyr::select(! c(result_id, variable_level, estimate_name, estimate_type)) |>
    mutate(across(! c(age_group, variable_name, table, window), ~ as.numeric(.x))) |>
    mutate(across(! c(age_group, variable_name, table, window), ~ replace(., . < 5, NA)))
    #rename("Count VAC" = vac_cohort, "Count UNVAC" = unvac_cohort, "Count Matched VAC" = total_vac_matched_cohort, "Count Matched UNVAC" = total_unvac_matched_cohort, "Count 1 Dose Matched" = doses1_matched_cohort, "Count >1 Dose Matched" = doses2_matched_cohort)
#mutate("counts (UN, VAC, UM, VM)" = paste0(unvac_cohort, ", ", vac_cohort, ", ", total_unvac_matched_cohort, ", ", total_vac_matched_cohort)) |>
#dplyr::select(! c(vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort))

tableOne_smd <- tableOne_red |>
  mutate(original_smd = abs((vac_cohort - unvac_cohort)/sqrt((vac_cohort*(1-vac_cohort) + unvac_cohort*(1-unvac_cohort))/2))) |>
  #mutate(overall_matched_smd = abs((total_vac_matched_cohort - total_unvac_matched_cohort)/sqrt((total_vac_matched_cohort*(1-total_vac_matched_cohort) + total_unvac_matched_cohort*(1-total_unvac_matched_cohort))/2))) |>
  mutate(dose_matched_smd = abs((doses2_matched_cohort - doses1_matched_cohort)/sqrt((doses2_matched_cohort*(1-doses2_matched_cohort) + doses1_matched_cohort*(1-doses1_matched_cohort))/2))) |>
  mutate(across(c(original_smd, dose_matched_smd), ~replace(., is.na(.), 0))) |>
  # mutate(across(c(original_smd, overall_matched_smd, dose_matched_smd), ~replace(., is.na(.), 0))) |>
  #dplyr::select(c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table, window, original_smd, overall_matched_smd, dose_matched_smd)) |>
  dplyr::select(c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table, window, original_smd, dose_matched_smd)) |>
  left_join(tableOne_counts, by = c("variable_name","age_group","window","table")) |>
  dplyr::select(! c(estimate_name, estimate_type))


write.csv(tableOne_smd, paste0(resultsFolder,"/tableOne_SMD_",cdmSchema,".csv"), row.names = FALSE)

# to_help <- tableOne_tidy |>
#   filter(estimate_name == "count") |>
#   dplyr::select(age_group, variable_name, window, table, vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort) |>
#   left_join(tableOne_smd |> dplyr::select(age_group, variable_name, window, table, original_smd, matched_smd"Count VAC", "Count UNVAC", "Count Matched VAC", "Count Matched UNVAC"), 
#             by = c("variable_name","age_group","window","table"))
# to_help <- tableOne_tidy |>
#   filter(estimate_name == "count") |>
#   dplyr::select(age_group, variable_name, window, table, vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort) |>
#   left_join(tableOne_smd |> select(age_group, variable_name, window, table, original_smd, matched_smd, "counts (UN, VAC, UM, VM)"), 
#             by = c("variable_name","age_group","window","table"))

# Vaccinated cohort
info(logger, "LARGE CHARACTERISATION SUMMARY")

# matched_population <- union_all(cdm$total_unvac_matched_cohort, cdm$total_vac_matched_cohort) |>
#   rename(matched_strata = vac_status) |>
#   dplyr::select(subject_id, matched_strata)
# doses_population <- union_all(cdm$doses1_matched_cohort, cdm$doses2_matched_cohort) |>
#   rename(doses_strata = treatment) |>
#   dplyr::select(subject_id, doses_strata)
# 
# cdm <- CDMConnector::insertTable(table = total_population, name = "hpv_cohorts")
# cdm$hpv_cohorts <- union_all(cdm$unvac_cohort, cdm$vac_cohort) |> 
#   compute(name = "hpv_cohorts") |>
#   dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date","vac_status") |>
#   left_join(matched_population, by = "subject_id") |>
#   left_join(doses_population, by = "subject_id")

slc <- summariseLargeScaleCharacteristics(cdm$hpv_cohorts, 
                                          window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                          eventInWindow = c("condition_occurrence","drug_exposure"),
                                          minimumFrequency = 0
)

write.csv(slc |> suppress(minCellCount = 5),
          paste0(resultsFolder,"/largeScale_summary_",cdmSchema,".csv"), row.names = FALSE)

# ADDD AGE GROUP
info(logger, "Tidy Large Scale Characteristics")
slc_tidy <- slc |>
  # as_tibble() |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_type, package_name, package_version, type, analysis)) |> 
  dplyr::mutate(across(! c(result_id, variable_name, variable_level, estimate_name, estimate_type, table_name, concept), ~replace(., is.na(.), 0)))

slc_tidy_def <- slc_tidy |>
  filter(estimate_name == "count") |>
  mutate(across(! c(result_id, variable_name, variable_level, estimate_name, estimate_type, table_name, concept), ~ as.numeric(.x))) |>
  mutate(across(! c(result_id, variable_name, variable_level, estimate_name, estimate_type, table_name, concept), ~ replace(., . < 5, NA))) |>
  mutate(across(! c(result_id, variable_name, variable_level, estimate_name, estimate_type, table_name, concept), ~ as.character(.x))) |>
  union_all(slc_tidy |>
              filter(estimate_name != "count")) 

write.csv(slc_tidy_def, 
          paste0(resultsFolder,"/largeScale_tidy_",cdmSchema,".csv"), row.names = FALSE)

info(logger, "Compute SMD for LSC")
slc_red <- slc_tidy |> 
  filter(estimate_name == "percentage") |>
  mutate(across(! c(result_id, variable_name, variable_level, estimate_name, estimate_type, table_name, concept), ~ as.numeric(.x)/100))

slc_counts <- slc_tidy |>
  filter(estimate_name == "count") |>
  dplyr::select(! c(result_id, estimate_name, estimate_type)) |>
  mutate(across(! c(variable_name, variable_level, table_name, concept), ~ as.numeric(.x))) |>
  mutate(across(! c(variable_name, variable_level, table_name, concept), ~ replace(., . < 5, NA)))
#rename("Count VAC" = vac_cohort, "Count UNVAC" = unvac_cohort, "Count Matched VAC" = total_vac_matched_cohort, "Count Matched UNVAC" = total_unvac_matched_cohort, "Count 1 Dose Matched" = doses1_matched_cohort, "Count >1 Dose Matched" = doses2_matched_cohort)
#mutate("counts (UN, VAC, UM, VM)" = paste0(unvac_cohort, ", ", vac_cohort, ", ", total_unvac_matched_cohort, ", ", total_vac_matched_cohort)) |>
#dplyr::select(! c(vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort))

slc_smd <- slc_red |>
  mutate(original_smd = abs((vac_cohort - unvac_cohort)/sqrt((vac_cohort*(1-vac_cohort) + unvac_cohort*(1-unvac_cohort))/2))) |>
  #mutate(overall_matched_smd = abs((total_vac_matched_cohort - total_unvac_matched_cohort)/sqrt((total_vac_matched_cohort*(1-total_vac_matched_cohort) + total_unvac_matched_cohort*(1-total_unvac_matched_cohort))/2))) |>
  mutate(dose_matched_smd = abs((doses2_matched_cohort - doses1_matched_cohort)/sqrt((doses2_matched_cohort*(1-doses2_matched_cohort) + doses1_matched_cohort*(1-doses1_matched_cohort))/2))) |>
  mutate(across(c(original_smd, dose_matched_smd), ~replace(., is.na(.), 0))) |>
  # mutate(across(c(original_smd, overall_matched_smd, dose_matched_smd), ~replace(., is.na(.), 0))) |>
  #dplyr::select(c(result_id, age_group, variable_name, variable_level, estimate_name, estimate_type, table, window, original_smd, overall_matched_smd, dose_matched_smd)) |>
  dplyr::select(c(result_id, variable_name, variable_level, estimate_name, estimate_type, table_name, concept, original_smd, dose_matched_smd)) |>
  left_join(slc_counts, by = c("variable_name","variable_level", "table_name","concept")) |>
  dplyr::select(! c(estimate_name, estimate_type))

write.csv(slc_smd, paste0(resultsFolder,"/largeScale_SMD_",cdmSchema,".csv"), row.names = FALSE)

# slc_help <- slc_tidy |>
#   filter(estimate_name == "count") |>
#   dplyr::select(variable_name, variable_level, table_name, concept, vac_cohort, unvac_cohort, total_unvac_matched_cohort, total_vac_matched_cohort) |>
#   left_join(slc_smd |> dplyr::select(variable_name, variable_level, table_name, concept, original_smd, matched_smd, "Count VAC", "Count UNVAC", "Count Matched VAC", "Count Matched UNVAC"), 
#             by = c("variable_name","variable_level", "table_name","concept"))

# Save attritions
info(logger, "SAVE COHORT ATTRITIONS")

attritions <- attrition(cdm$hpv_cohorts) |>
  left_join(settings(cdm$hpv_cohorts), by = "cohort_definition_id") |>
  dplyr::select(! cohort_definition_id) |>
  setcolorder(c("cohort_name", "number_records", "number_subjects", "reason_id", "reason", "excluded_records", "excluded_subjects"))

write.csv(attritions, paste0(resultsFolder,"/study_attritions_",cdmSchema,".csv"), row.names = FALSE)
