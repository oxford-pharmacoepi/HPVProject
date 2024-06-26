to_Charac <- function(table, sex, index_date = FALSE) {
  if (sex == "f") {
    yy <- 18
  } else {
    yy <- 15
  }
  
  if (index_date) {
    to <- table |> 
      dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date")
  } else {
    to <- table |> 
      dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
      addDateOfBirth() |>
      mutate(cohort_start_date = date_of_birth + years(yy))
  }
  
  if (sex == "f") {
    print("female")
    to |>
      mutate(cohort_start_date = as.Date(cohort_start_date)) |>
      mutate(cohort_end_date = as.Date(cohort_end_date)) |>
      summariseCharacteristics(
        cohortIntersect = list(
          # "Medications" = list(
          #   targetCohortTable = "medications", value = "flag", window = list(c(-365,-1),c(-30,-1))
          # ),
          "Conditions" = list(
            targetCohortTable = "conditions", value = "flag", window = list(c(-Inf,-1))
          ),
          "HIV_status" = list(
            targetCohortTable = "hiv_status", value = "flag", window = c(-Inf,-1)
          ),
          "vaccinations" = list(
            targetCohortTable = "vaccinations", value = "count", window = c(-Inf,-1)
          ),
          "cervical_screening" = list(
            targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
          ),
          "smear_test" = list(
            targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
          )
        ),
      )
  } else {
    print("male")
    to |>
      mutate(cohort_start_date = as.Date(cohort_start_date)) |>
      mutate(cohort_end_date = as.Date(cohort_end_date)) |>
      summariseCharacteristics(
        cohortIntersect = list(
          # "Medications" = list(
          #   targetCohortTable = "medications", value = "flag", window = list(c(-365,-1),c(-30,-1))
          # ),
          "Conditions" = list(
            targetCohortTable = "conditions", value = "flag", window = list(c(-Inf,-1))
          ),
          "HIV_status" = list(
            targetCohortTable = "hiv_status", value = "flag", window = c(-Inf,-1)
          ),
          "vaccinations" = list(
            targetCohortTable = "vaccinations", value = "count", window = c(-Inf,-1)
          )
        ),
      )
  }
}

SMD_to <- function(to_table, target_cohort, comparator_cohort) {
  to_table |>
  visOmopResults::splitAll() |>
  dplyr::select(result_id, cdm_name, result_type, package_name, package_version, cohort_name, variable_name, variable_level, estimate_name, estimate_type, estimate_value, table, window, value) |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(all_of(c(target_cohort, comparator_cohort)), ~ as.numeric(.x)/100)) |>
  mutate(across(all_of(c(target_cohort, comparator_cohort)), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((!!sym(target_cohort) - !!sym(comparator_cohort))/sqrt((!!sym(target_cohort)*(1-!!sym(target_cohort)) + !!sym(comparator_cohort)*(1-!!sym(comparator_cohort)))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)
}

SMD_lsc <- function(lsc_table, target_cohort, comparator_cohort) {
  lsc_table |>
    visOmopResults::splitAll() |>
    dplyr::select(result_id, cdm_name, result_type, package_name, package_version, cohort_name, variable_name, variable_level, estimate_name, estimate_type, estimate_value, concept_id) |>
    tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
    dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version, overall)) |> 
    dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type), ~replace(., is.na(.), 0))) |>
    filter(estimate_name == "percentage") |>
    filter(! variable_name == "Sex") |>
    mutate(across(c(target_cohort, comparator_cohort), ~ as.numeric(.x)/100)) |>
    mutate(across(c(target_cohort, comparator_cohort), ~ as.numeric(.x))) |>
    mutate(SMDs = abs((!!sym(target_cohort) - !!sym(comparator_cohort))/sqrt((!!sym(target_cohort)*(1-!!sym(target_cohort)) + !!sym(comparator_cohort)*(1-!!sym(comparator_cohort)))/2))) |>
    arrange(desc(SMDs)) |>
    dplyr::select(! estimate_type)
}

bar_plot <- function(to_table, lsc_table, target_cohort, comparator_cohort, target_label, comparator_label, file_name) {
  to_table <- to_table |> dplyr::select(c(variable_name, variable_level, estimate_name, table, window, value, comparator_cohort, target_cohort, SMDs))
  lsc_table <- lsc_table |> dplyr::select(c(variable_name, variable_level, estimate_name, concept_id, comparator_cohort, target_cohort, SMDs))
  initial_graphic <- rbind(
    # specific
    to_table |> filter(variable_name %in% c("Hiv status", "Papanicolau smear testing")) |> mutate(type = "specific", variable_name = "Test") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
    to_table |> filter(variable_name == "cervical_screening") |> mutate(type = "specific") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
    to_table |> filter(variable_level == "Autoimmune disease") |> mutate(type = "specific") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
    #to_table |> filter(variable_level == "immunosupressants", window == "-365 to -1") |> mutate(type = "specific") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
    # general > 0.1
    lsc_table |> filter(variable_name %in% c("Conjunctivitis", "Pain in throat", "Eruption"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id),
    lsc_table |> filter(variable_name %in% c("tetanus toxoid vaccine, inactivated", "measles virus vaccine", "ibuprofen 20 MG/ML Oral Suspension"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id),
    # general < 0.1
    lsc_table |> filter(variable_name %in% c("Anxiety", "Cough", "Asthma daytime symptoms"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id),
    lsc_table |> filter(variable_name %in% c("amoxicillin 500 MG Oral Capsule", "penicillin V potassium 250 MG Oral Tablet", "Benzydamine 1.5 MG/ML Nasal Spray"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id)
  ) |>
    dplyr::select(! estimate_name)
  
  initial_graphic_m0.1 <- initial_graphic |>
    filter(SMDs <= 0.1) |>
    dplyr::select(variable_name, target_cohort, comparator_cohort) |>
    reshape2::melt(id.vars = "variable_name") |>
    ggplot(aes(x = variable_name, y = value, fill = variable)) + 
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values = setNames(c("#7EB3EB", "#FDBF69"), c(target_cohort, comparator_cohort)), 
                      labels = c(target_label, comparator_label)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    labs(title = "Initial characterisation: ASMD < 0.1", x = "Characteristics", y = "Proportion", fill = "Cohorts")
  
  ggsave(paste0(resultsFolder, "/", file_name, "_m0.1.png"), plot = initial_graphic_m0.1)
  
  initial_graphic_p0.1 <- initial_graphic |>
    filter(SMDs > 0.1) |>
    dplyr::select(variable_name, target_cohort, comparator_cohort) |>
    reshape2::melt(id.vars = "variable_name") |>
    ggplot(aes(x = variable_name, y = value, fill = variable)) + 
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values = setNames(c("#7EB3EB", "#FDBF69"), c(target_cohort, comparator_cohort)), 
                      labels = c(target_label, comparator_label)) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                         plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    labs(title = "Initial characterisation: ASMD > 0.1", x = "Characteristics", y = "Proportion", fill = "Cohorts")
  
  ggsave(paste0(resultsFolder, "/", file_name, "_p0.1.png"), plot = initial_graphic_p0.1)
  
}
