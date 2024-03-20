# Get data ----
# All potential outputs
library(dplyr)
library(here)
library(tidyr)
library(dplyr)
library(here)
library(readr)
library(fresh)
library(visOmopResults)
# Define result names
result_names <- c("cohort_definitions", "cohort_count", "code_counts", "cohort_overlap", 
                  "age_distribution", "time_distribution", "prevalence", "incidence", 
                  "index_events", "lsc_sample", "lsc_matched", "lsc_difference", "log",
                  "snapshot")

dataFolder <- "PhenoResults"
# Result files
result_files <- list.files(path = here(dataFolder), pattern = "\\.csv$", full.names = TRUE)

# Initialize data list with result_names
data <- vector("list", length(result_names))
names(data) <- result_names

# Loop through each file
for (file_path in result_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Check which result_name is in file_name
  for (resName in result_names) {
    if (grepl(resName, file_name)) {
      # Load CSV file
      csv_data <- read_csv(file_path)
      if (resName=="snapshot") {
        csv_data <- csv_data |> mutate(cdm_version=as.character(cdm_version))
      }
      # If data list already has data for this resName, append new data; otherwise, assign it
      if (!is.null(data[[resName]])) {
        data[[resName]] <- bind_rows(data[[resName]], csv_data)
      } else {
        data[[resName]] <- csv_data
      }
      
      # Break the loop once a match is found and assigned
      break
    }
  }
}
# Tranform data for shiny ----
# Orphan code counts
data$orphan_counts <- data$code_counts %>% 
  filter(strata_name == "recomendation") %>% 
  ungroup()  %>% 
  distinct() %>% 
  mutate(standard_concept_name = substr(additional_level, 1,
                                        unlist(gregexpr(';', additional_level))-2)) %>% 
  pivot_wider(names_from = variable_name, values_from = estimate_value) %>% 
  select("cdm_name", "cohort", "relationship_id",
         "standard_concept_id", "standard_concept_name", "Record count", "Person count")
# Code counts
data$code_counts <- data$code_counts %>% 
  filter(strata_name == "original_codes") %>% 
  ungroup()  %>% 
  distinct() %>% 
  mutate(standard_concept_name = substr(additional_level, 1,
                                        unlist(gregexpr(';', additional_level))-2)) %>% 
  pivot_wider(names_from = variable_name, values_from = estimate_value) %>% 
  select("cdm_name", "cohort", 
         "standard_concept_id", "standard_concept_name", "Record count", "Person count")
# Index events
data$index_events <- data$index_events %>% 
  pivot_wider(names_from = variable_name, values_from = estimate) %>% 
  select(cdm_name, cohort_name, codelist_name, group_name, group_level, standard_concept_id, standard_concept_name, 
         source_concept_name, source_concept_id, domain_id,  
         cdm_name, `Record count`, `Person count`)
# Cohort overlap

if(data$cohort_overlap %>% dplyr::tally() != 0){
data$cohort_overlap <- data$cohort_overlap %>%
  ungroup() %>%
  inner_join(data$cohort_count %>%
               select(cdm_name,
                      cohort_definition_id_x = cohort_definition_id,
                      cohort_name_x = cohort_name,
                      subject_counts_x = number_subjects),
             by = c("cdm_name", "cohort_definition_id_x")) %>%
  inner_join(data$cohort_count %>%
               select(cdm_name,
                      cohort_definition_id_y = cohort_definition_id,
                      cohort_name_y = cohort_name,
                      subject_counts_y = number_subjects),
             by = c("cdm_name", "cohort_definition_id_y")) %>%
  mutate(
    intersect_counts = as.integer(intersect_count)) %>%
  select(-intersect_count)}
# # Age distribution
data$age_distribution <- data$age_distribution %>%
  ungroup() %>%
  inner_join(data$cohort_count %>% select(cdm_name, cohort_definition_id, cohort_name)) %>%
  select(-cohort_definition_id)
# Time distribution
to_pivot <- colnames(data$time_distribution)[!colnames(data$time_distribution) %in%
                                               c("sex", "cohort_name", "cdm_name", "cohort_definition_id")]
data$time_distribution <- tibble(covariate = c("age", "prior_observation", "future_observation")) %>%
  left_join(
    data$time_distribution %>%
      ungroup() %>%
      inner_join(data$cohort_count %>% select(cdm_name, cohort_definition_id, cohort_name)) %>%
      select(-cohort_definition_id) %>%
      mutate(across(tidyr::everything(), ~ as.character(.x))) %>%
      pivot_longer(cols = to_pivot, values_to = "estimate_value") %>%
      mutate(
        estimate_type = case_when(
          grepl("Min", name) ~ "min",
          grepl("Max", name) ~ "max",
          grepl("Median", name) ~ "median",
          grepl("Mean", name) ~ "mean",
          grepl("Sd", name) ~ "sd",
        ),
        covariate = gsub("_Min|_Max|_Median|_Sd|_Mean", "", name),
        estimate_value = round(as.numeric(estimate_value), 3)
      ) %>%
      select(-name)
  ) %>%
  select(cdm_name, cohort_name, sex, covariate, estimate_type, estimate_value)
# LSC
data$lsc_table <- data$lsc_matched %>% 
  splitAdditional() %>%
  distinct() %>%
  mutate(
    estimate_name = paste0("matched_", estimate_name),
    estimate = as.numeric(estimate_value)
  ) %>% 
  pivot_wider(id_cols = c("variable_name", "variable_level", "group_level",
                          "table_name", "concept"),
    names_from = estimate_name, values_from = estimate) %>% 
  left_join(
    data$lsc_sample %>% 
      splitAdditional() %>%
      distinct() %>%
      mutate(
        estimate_name = paste0("sample_", estimate_name),
        estimate = as.numeric(estimate_value)
      ) %>% 
      pivot_wider(names_from = estimate_name, values_from = estimate)) %>% 
  mutate(
    difference_count = (sample_count - matched_count)/matched_count,
    difference_percentage = (sample_percentage - matched_percentage)/matched_percentage
  ) %>% 
  select(
    cdm_name, cohort_name = group_level, concept_name = variable_name, table_name, concept,
    window = variable_level, matched_count, matched_percentage, sample_count, sample_percentage, 
    difference_count, difference_percentage
  )
  
# Shiny theme ----
DUtheme <- create_theme(
  adminlte_color(
    light_blue = "#0c0e0c" 
  ),
  adminlte_sidebar(
    # width = "400px",
    dark_bg = "#78B7C5",
    dark_hover_bg = "#3B9AB2",
    dark_color = "white"
  ), 
  adminlte_global(
    content_bg = "#eaebea" 
  ),
  adminlte_vars(
    border_color = "#112446",
    active_link_hover_bg = "#FFF",
    active_link_hover_color = "#112446",
    active_link_hover_border_color = "#112446",
    link_hover_border_color = "#112446"
  )
)
