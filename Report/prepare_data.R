# Result folders
subFolders <- c("attritions", "lasso", "characterisation", "outcome_model_female_darwin", 
                "outcome_model_female_nondarwin", "outcome_model_male", "nco_female", 
                "nco_male", "calibration_female")
# Result files
for (sub in subFolders) {
  result_files <- list.files(path = here(dataFolder, sub), pattern = "\\.csv$", full.names = TRUE)
  assign(sub, vector("list"))

  for (file in result_files) {
    list <- get(sub)
    file_name <- tools::file_path_sans_ext(basename(file))
    list[[file_name]] <- read.csv(file)
    assign(sub, list)
  }
}
result_files <- list.files(path = here(dataFolder), pattern = "\\.csv$", full.names = TRUE)

# Attritions -----
# we want different files for different cohorts
names <- names(attritions)
for (name in names) {
  attritions[[name]] <- attritions[[name]] |> 
    dplyr::select(cohort_definition_id, number_records, number_subjects, reason_id, reason, excluded_records, excluded_subjects)
}

attritions <- bind_rows(attritions, .id = "cohort_name") |>
  mutate(cohort_name = gsub("att", "", 
                            gsub("_", " ", cohort_name)) |> 
           trimws())
# -------
# Lasso ------
names <- names(lasso)
for (name in names) {
  lasso[[name]] <- lasso[[name]] |> 
    dplyr::select(concept_id, coefficients, period = length_id, concept_name, year)
}

lasso <- bind_rows(lasso, .id = "study_group") |>
  mutate(study_group = gsub("lasso", "", 
                            gsub("female", "",
                                 gsub("cov0", "",
                                      gsub("_", " ", study_group)))
                            ) |> 
           trimws())
# ---------
# Characterisation -----
# Raw data
characterisation_lsc <- vector("list")
characterisation_to <- vector("list")
names <- names(characterisation)
for (name in names) {
  if (grepl("lsc",name)) {
    if (grepl("male",name)) {
      characterisation_lsc[[name]] <- characterisation[[name]] |>
        as_tibble() |>
        dplyr::select(group_name, group_level, strata_name, strata_level, variable_name, variable_level, estimate_name, estimate_type, estimate_value, additional_name, additional_level) |>
        visOmopResults::splitAll() |>
        mutate(sex = "male") |>
        filter(estimate_name %in% c("percentage", "count"))
      
    } else if (grepl("female",name)) {
      characterisation_lsc[[name]] <- characterisation[[name]] |>
        as_tibble() |>
        dplyr::select(group_name, group_level, strata_name, strata_level, variable_name, variable_level, estimate_name, estimate_type, estimate_value, additional_name, additional_level) |>
        visOmopResults::splitAll() |>
        mutate(sex = "female") |>
        filter(estimate_name %in% c("percentage", "count"))
      
    }
  } else if (grepl("to", name)) {
    if (grepl("male",name)) {
      characterisation_to[[name]] <- characterisation[[name]] |>
        as_tibble() |>
        dplyr::select(group_name, group_level, strata_name, strata_level, variable_name, variable_level, estimate_name, estimate_type, estimate_value, additional_name, additional_level) |>
        visOmopResults::splitAll() |>
        mutate(sex = "male") |>
        filter(estimate_name %in% c("percentage", "count"))
    } else if (grepl("female",name)) {
      characterisation_to[[name]] <- characterisation[[name]] |>
        as_tibble() |>
        dplyr::select(group_name, group_level, strata_name, strata_level, variable_name, variable_level, estimate_name, estimate_type, estimate_value, additional_name, additional_level) |>
        visOmopResults::splitAll() |>
        mutate(sex = "female") |>
        filter(estimate_name %in% c("percentage", "count"))
      
    }

  }
}

# Tidy table
names <- names(characterisation_lsc)
characterisation_lsc_smd <- vector("list")
characterisation_to_smd <- vector("list")

for (name in names) {
  if (grepl("female",name)) {
    yy <- 18
    ss <- "f"
    } else {
      yy <- 15
      ss <- "m"
      }
  if (grepl("overall",name) & grepl("initial", name)) {
    target <- paste0("vac_",yy,ss,"_coverage_0_any_dose")
    comparator <- paste0("unvac_",yy,ss,"_coverage_0_any_dose")
  } else if (grepl("overall",name) & grepl("matched", name)) {
    target <- paste0("matched_vac_cov_0_ratio_",yy,ss,"_any_dose")
    comparator <- paste0("matched_unvac_cov_0_ratio_",yy,ss,"_any_dose")
  } else if (grepl("1vs23",name) & grepl("initial", name)) {
    target <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_23dose")
    comparator <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_1dose")
  } else if (grepl("1vs23",name) & grepl("matched", name)) {
    target <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_23dose")
    comparator <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_1dose")
  } else if (grepl("2vs3",name) & grepl("initial", name)) {
    target <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_3dose")
    comparator <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_2dose")
  } else if (grepl("2vs3",name) & grepl("matched", name)) {
    target <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_3dose")
    comparator <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_2dose")
  }
  characterisation_lsc_smd[[name]] <- SMD_lsc(characterisation_lsc[[name]], target, comparator)
}

names <- names(characterisation_to)
for (name in names) {
  if (grepl("female",name)) {
    yy <- 18
    ss <- "f"
  } else {
    yy <- 15
    ss <- "m"
  }
  if (grepl("overall",name) & grepl("initial", name)) {
    target <- paste0("vac_",yy,ss,"_coverage_0_any_dose")
    comparator <- paste0("unvac_",yy,ss,"_coverage_0_any_dose")
  } else if (grepl("overall",name) & grepl("matched", name)) {
    target <- paste0("matched_vac_cov_0_ratio_",yy,ss,"_any_dose")
    comparator <- paste0("matched_unvac_cov_0_ratio_",yy,ss,"_any_dose")
  } else if (grepl("1vs23",name) & grepl("initial", name)) {
    target <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_23dose")
    comparator <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_1dose")
  } else if (grepl("1vs23",name) & grepl("matched", name)) {
    target <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_23dose")
    comparator <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_1dose")
  } else if (grepl("2vs3",name) & grepl("initial", name)) {
    target <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_3dose")
    comparator <- paste0("not_matched_vaccinated_",yy,ss,"_coverage_0_2dose")
  } else if (grepl("2vs3",name) & grepl("matched", name)) {
    target <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_3dose")
    comparator <- paste0("matched_vaccinated_",yy,ss,"_coverage_0_2dose")
  }
  characterisation_to_smd[[name]] <- SMD_to(characterisation_to[[name]], target, comparator)
}
  
# Save all in tibbles
characterisation_lsc <- bind_rows(characterisation_lsc, .id = "study_group")
  # mutate(study_group = gsub("lsc", "", 
  #                           gsub("_", " ", study_group)) |>
  #          trimws())

characterisation_to <- bind_rows(characterisation_to, .id = "study_group")
  # mutate(study_group = gsub("to", "", 
  #                           gsub("_", " ", study_group)) |>
  #          trimws())

# -------------

# Outcome model female -------
names_darwin <- names(outcome_model_female_darwin)
names_nondarwin <- names(outcome_model_female_nondarwin)
names_cal <- names(calibration_female)
IR_female_darwin <- vector("list")
IR_female_nondarwin <- vector("list")
IR_female_cal <- vector("list")
IRR_female_darwin <- vector("list")
IRR_female_nondarwin <- vector("list")
IRR_female_cal <- vector("list")
HR_female_darwin <- vector("list")
HR_female_nondarwin <- vector("list")
HR_female_cal <- vector("list")

for (name in names_darwin) {
  if (grepl("IR_",name)) {
    IR_female_darwin[[name]] <- outcome_model_female_darwin[[name]]
  } else if (grepl("IRR_", name)) {
    IRR_female_darwin[[name]] <- outcome_model_female_darwin[[name]] |> filter(term == "treatment")
  } else if (grepl("Cox_", name)) {
    HR_female_darwin[[name]] <- outcome_model_female_darwin[[name]] |> filter(term == "treatment")
  }
}
for (name in names_nondarwin) {
  if (grepl("IR_",name)) {
    IR_female_nondarwin[[name]] <- outcome_model_female_nondarwin[[name]]
  } else if (grepl("IRR_", name)) {
    IRR_female_nondarwin[[name]] <- outcome_model_female_nondarwin[[name]] |> filter(term == "treatment")
  } else if (grepl("Cox_", name)) {
    HR_female_nondarwin[[name]] <- outcome_model_female_nondarwin[[name]] |> filter(term == "treatment")
  }
}
for (name in names_cal) {
  if (grepl("IR_",name) & grepl("matched",name)) {
    IR_female_cal[[name]] <- calibration_female[[name]]
  } else if (grepl("IRR_",name) & grepl("matched",name)) {
    IRR_female_cal[[name]] <- calibration_female[[name]]
  } else if (grepl("HR_",name) & grepl("matched",name)) {
    HR_female_cal[[name]] <- calibration_female[[name]]
  }
}

IR_female_results <- rbind(
  bind_rows(IR_female_nondarwin, .id = "study_group") |>
  mutate(study_group = gsub("IR_", "",
                            gsub("_", " ", study_group))
         ),
  bind_rows(IR_female_darwin, .id = "study_group") |>
    mutate(study_group = gsub("IR_", "",
                              gsub("_", " ", study_group))
           )
  ) |>
  mutate(sex = "female")
  print("IR done")
  
IRR_female_results <- rbind(
  bind_rows(IRR_female_nondarwin, .id = "study_group") |>
    mutate(study_group = gsub("IRR", "",
                              gsub("_", " ", study_group))) |>
    dplyr::select(study_group, estimate, std.error, statistic, p.value, outcome, year, irr, lower, upper),
  bind_rows(IRR_female_darwin, .id = "study_group") |>
    mutate(study_group = gsub("IRR", "",
                              gsub("_", " ", study_group))) |>
    dplyr::select(study_group, estimate, std.error, statistic, p.value, outcome, year, irr, lower, upper), 
  bind_rows(IRR_female_cal, .id = "study_group") |>
    mutate(study_group = gsub("IRR", "",
                              gsub("_", " ", 
                                   gsub("matched", "calibrated", study_group)))) |>
    rename(outcome = name) |>
    mutate(year = str_extract(study_group, "\\d+$")) |>
    mutate(study_group = str_replace(study_group, "\\d+$", "")) |>
    dplyr::select(study_group, estimate, std.error,outcome, year, irr, lower, upper) |>
    mutate(statistic = NA, p.value = NA)
  ) |>
  mutate(sex = "female", year = as.integer(year))
print("IRR done")


HR_female_results <- rbind(
  bind_rows(HR_female_nondarwin, .id = "study_group") |>
    mutate(study_group = gsub("HR", "",
                              gsub("_", " ", study_group))) |>
    dplyr::select(study_group, estimate, std.error, statistic, p.value, outcome, hr, lower, upper),
  bind_rows(HR_female_darwin, .id = "study_group") |>
    mutate(study_group = gsub("HR", "",
                              gsub("_", " ", study_group))) |>
    dplyr::select(study_group, estimate, std.error, statistic, p.value, outcome,hr, lower, upper),
  bind_rows(HR_female_cal, .id = "study_group") |>
    mutate(study_group = gsub("HR", "",
                              gsub("_", " ", 
                                   gsub("matched", "calibrated", study_group)))) |>
    rename(outcome = name) |>
    mutate(study_group = str_replace(study_group, "\\d+$", "")) |>
    dplyr::select(study_group, estimate, std.error,outcome, hr, lower, upper) |>
    mutate(statistic = NA, p.value = NA)
  ) |>
  mutate(sex = "female")
print("HR done")
# ------------
# Outcome model for males ------
names <- names(outcome_model_male)
IR_male <- vector("list")
IRR_male <- vector("list")
HR_male <- vector("list")

for (name in names) {
  if (grepl("IR_",name)) {
    IR_male[[name]] <- outcome_model_male[[name]]
  } else if (grepl("IRR_", name)) {
    IRR_male[[name]] <- outcome_model_male[[name]] |> filter(term == "treatment")
  } else if (grepl("Cox_", name)) {
    HR_male[[name]] <- outcome_model_male[[name]] |> filter(term == "treatment")
  }
}

IR_male_results <- bind_rows(IR_male, .id = "study_group") |>
    mutate(study_group = gsub("IR_", "",
                              gsub("_", " ", study_group)),
           sex = "male"
           )


IRR_male_results <- bind_rows(IRR_male, .id = "study_group") |>
    mutate(study_group = gsub("IRR_", "",
                              gsub("_", " ", study_group)),
           sex = "male"
           ) |>
  dplyr::select(study_group, sex, estimate, std.error, statistic, p.value, outcome, year, irr, lower, upper)


HR_male_results <- bind_rows(HR_male, .id = "study_group") |>
    mutate(study_group = gsub("HR_", "",
                              gsub("_", " ", study_group)),
           sex = "male"
           ) |>
  dplyr::select(study_group, sex, estimate, std.error, statistic, p.value, outcome, hr, lower, upper)

# ----------
# Add all results -------
IR_results <- union_all(IR_female_results, IR_male_results) |>
  mutate(outcome = gsub("_", " ", outcome)) |>
  mutate(outcome = case_when(
    outcome == "hpv cin23 cohort" ~ "cin23",
    outcome == "hpv conization cohort broad" ~ "conization broad",
    outcome == "hpv conization cohort narrow" ~ "conization narrow",
    outcome == "hpv hpv screening cohort" ~ "cervical screening",
    outcome == "hpv inv cancer broad cohort" ~ "cervical cancer broad",
    outcome == "hpv inv cancer narrow cohort" ~ "cervical cancer narrow",
    outcome == "hpv smear test cohort" ~ "smear test",
    .default = outcome
  ))
IRR_results <- union_all(IRR_female_results, IRR_male_results) |>
  mutate(outcome = gsub("_", " ", outcome)) |>
  mutate(outcome = case_when(
    outcome == "hpv cin23 cohort" ~ "cin23",
    outcome == "hpv conization cohort broad" ~ "conization broad",
    outcome == "hpv conization cohort narrow" ~ "conization narrow",
    outcome == "hpv hpv screening cohort" ~ "cervical screening",
    outcome == "hpv inv cancer broad cohort" ~ "cervical cancer broad",
    outcome == "hpv inv cancer narrow cohort" ~ "cervical cancer narrow",
    outcome == "hpv smear test cohort" ~ "smear test",
    .default = outcome
  ))
HR_results <- union_all(HR_female_results, HR_male_results) |>
  mutate(outcome = gsub("_", " ", outcome)) |>
  mutate(outcome = case_when(
    outcome == "hpv cin23 cohort" ~ "cin23",
    outcome == "hpv conization cohort broad" ~ "conization broad",
    outcome == "hpv conization cohort narrow" ~ "conization narrow",
    outcome == "hpv hpv screening cohort" ~ "cervical screening",
    outcome == "hpv inv cancer broad cohort" ~ "cervical cancer broad",
    outcome == "hpv inv cancer narrow cohort" ~ "cervical cancer narrow",
    outcome == "hpv smear test cohort" ~ "smear test",
    .default = outcome
  ))
# ----------

# NCO female -----
names_nco <- names(nco_female)
IR_female_nco <- vector("list")
IRR_female_nco <- vector("list")
HR_female_nco <- vector("list")

for (name in names_nco) {
  if (grepl("IR_",name)) {
    IR_female_nco[[name]] <- nco_female[[name]]
  } else if (grepl("IRR_", name)) {
    IRR_female_nco[[name]] <- nco_female[[name]] |> filter(term == "treatment")
  } else if (grepl("Cox_", name)) {
    HR_female_nco[[name]] <- nco_female[[name]] |> filter(term == "treatment")
  }
}

IR_female_nco <- bind_rows(IR_female_nco, .id = "study_group") |>
  mutate(study_group = gsub("IR_", "",
                              gsub("_", " ", study_group))
         ) |>
  mutate(sex = "female")

IRR_female_nco <- bind_rows(IRR_female_nco, .id = "study_group") |>
  mutate(study_group = gsub("IRR_", "",
                            gsub("_", " ", study_group))
  ) |>
  mutate(sex = "female")

HR_female_nco <- bind_rows(HR_female_nco, .id = "study_group") |>
  mutate(study_group = gsub("HR_", "",
                            gsub("_", " ", study_group))
  ) |>
  mutate(sex = "female")
# ---------

# NCO male -------
names_nco <- names(nco_male)
IR_male_nco <- vector("list")
IRR_male_nco <- vector("list")
HR_male_nco <- vector("list")

for (name in names_nco) {
  if (grepl("IR_",name)) {
    IR_male_nco[[name]] <- nco_male[[name]]
  } else if (grepl("IRR_", name)) {
    IRR_male_nco[[name]] <- nco_male[[name]] |> filter(term == "treatment")
  } else if (grepl("Cox_", name)) {
    HR_male_nco[[name]] <- nco_male[[name]] |> filter(term == "treatment")
  }
}

IR_male_nco <- bind_rows(IR_male_nco, .id = "study_group") |>
  mutate(study_group = gsub("IR_", "",
                            gsub("_", " ", study_group))
  ) |>
  mutate(sex = "male")

IRR_male_nco <- bind_rows(IRR_male_nco, .id = "study_group") |>
  mutate(study_group = gsub("IRR_", "",
                            gsub("_", " ", study_group))
  ) |>
  mutate(sex = "male")

HR_male_nco <- bind_rows(HR_male_nco, .id = "study_group") |>
  mutate(study_group = gsub("HR_", "",
                            gsub("_", " ", study_group))
  ) |>
  mutate(sex = "male")
# --------------
# Add all NCO -------
IR_nco <- union_all(IR_female_nco, IR_male_nco) |>
  mutate(outcome = gsub("_", " ", outcome))

IRR_nco <- union_all(IRR_female_nco, IRR_male_nco) |>
  mutate(outcome = gsub("_", " ", outcome))

HR_nco <- union_all(HR_female_nco, HR_male_nco) |>
  mutate(outcome = gsub("_", " ", outcome))
# ----------
