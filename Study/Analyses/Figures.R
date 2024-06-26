library(reshape2)
library(ggplot2)
library(stringr)

# Create logger ---------------------
resultsFolder <- here("Figures")
loggerName <- gsub(":| |-", "_", paste0("log ", Sys.time(),".txt"))
logger <- create.logger()
logfile(logger) <- here(resultsFolder, loggerName)
level(logger) <- "INFO"
info(logger, "LOG CREATED")

source(here("Analyses", "Figures_functions.R"))

cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = dbName, 
  achillesSchema = achillesSchema,
  cohortTables = c("firstdose_cohort", 
                   "doses_allvac_cohort", 
                   # Coverage 0%
                   "unvac_18f_coverage_0_any_dose",
                   "vac_18f_coverage_0_any_dose",
                   "matched_vac_cov_0_ratio_18f_any_dose",
                   "matched_unvac_cov_0_ratio_18f_any_dose",
                   "unvac_15m_coverage_0_any_dose",
                   "vac_15m_coverage_0_any_dose",
                   # "vac_18f_coverage_0_1dose",
                   # "vac_18f_coverage_0_23dose",
                   # "vac_18f_coverage_0_2dose",
                   # "vac_18f_coverage_0_3dose",
                   # "matched_vac_18f_coverage_0_1dose",
                   # "matched_vac_18f_coverage_0_23dose",
                   # "matched_vac_18f_coverage_0_2dose",
                   # "matched_vac_18f_coverage_0_3dose",
                   # Coverage 60%
                   # "unvac_18f_coverage_60_any_dose",
                   # "vac_18f_coverage_60_any_dose",
                   # "matched_unvac_18f_coverage_60_any_dose",
                   # "matched_vac_18f_coverage_60_any_dose",
                   # "vac_18f_coverage_60_1dose",
                   # "vac_18f_coverage_60_23dose",
                   # "vac_18f_coverage_60_2dose",
                   # "vac_18f_coverage_60_3dose",
                   # "matched_vac_18f_coverage_60_1dose",
                   # "matched_vac_18f_coverage_60_23dose",
                   # "matched_vac_18f_coverage_60_2dose",
                   # "matched_vac_18f_coverage_60_3dose",
                   "medications",
                   "conditions",
                   "hiv_status",
                   "vaccinations",
                   "cervical_screening",
                   "smear_test",
                   "nondarwin_outcome_cohorts"
                  # "outcome_cohorts"
                   )
  
)

# histogram of the distribution of the age at first vac -----------
first_vac_data_males <- cdm$firstdose_cohort |>
  addSex() |>
  filter(sex == "Male") |>
  addAge() |>
  select(age) |>
  collect() |>
  table()

first_vac_data_fem <- cdm$firstdose_cohort |>
  addSex() |>
  filter(sex == "Female") |>
  addAge() |>
  select(age) |>
  collect() |>
  table()

hist_ageatfirstvac <- barplot(first_vac_data_fem,
                              main = "Histogram: age at first vaccination",
                              xlab = "Age", ylab = "Frequency",
                              xlim = c(5,40),   # Limit x-axis
                              ylim = c(0, max(first_vac_data) + 5000),  # Limit y-axis
                              col = "#FF9933"
)
barplot(first_vac_data_males, add = TRUE, col = "#006666")

text(x = hist_ageatfirstvac, y = first_vac_data_fem + 5500, labels = first_vac_data_fem)
text(x = hist_ageatfirstvac, y = first_vac_data_males + 2000, labels = first_vac_data_males)
legend("topright", legend = c("Women", "Men"), fill = c("#FF9933", "#006666"))
# ---------------------

# Characterisation of all vaccinated people --------------
table_one_summary <- cdm$firstdose_cohort |>
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addSex() |>
  addAge(cdm, ageGroup = list(c(0, 14), c(15, 16), c(17, 18))) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
    strata = list(c("age_group", "sex")),
    ageGroup = list(
      c(0, 14), c(15, 16), c(17, 18)
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
      "cervical_screening" = list(
        targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
      ),
      "smear_test" = list(
        targetCohortTable = "smear_test", value = c("count"), window = c(-Inf,-1)
      )
    ),
  )

formattable <- formatTable(table_one_summary,
                           formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                  "N" = "<count>",
                                                  "Mean (SD)" = "<mean> (<sd>)"),
                           header = c("group", "strata"),
                           split = c("group","strata", "additional"),
                           excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                           type = "gt")
# ------------

# CHARACTERISATIONS ---------
for (cov in c(0,60)) {
  cdm <- omopgenerics::bind(
    cdm[[paste0("vac_18f_coverage_",cov,"_any_dose")]],
    cdm[[paste0("unvac_18f_coverage_",cov,"_any_dose")]],
    name = "hpv_initial_overall_female_cohorts"
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("vac_15m_coverage_",cov,"_any_dose")]],
    cdm[[paste0("unvac_15m_coverage_",cov,"_any_dose")]],
    name = "hpv_initial_overall_male_cohorts"
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("matched_vac_cov_",cov,"_ratio_18f_any_dose")]],
    cdm[[paste0("matched_unvac_cov_",cov,"_ratio_18f_any_dose")]],
    name = "hpv_ratio_matched_overall_female_cohorts"
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("matched_vac_cov_",cov,"_ratio_15m_any_dose")]],
    cdm[[paste0("matched_unvac_cov_",cov,"_ratio_15m_any_dose")]],
    name = "hpv_ratio_matched_overall_male_cohorts"
  )
  
  # FEMALES
  yy <- 18
  ss <- "f"
  # Before matching overall
  to_initial_overall_female <- to_Charac(cdm$hpv_initial_overall_female_cohorts, sex = "f", index_date = FALSE)
  write.csv(to_initial_overall_female, paste0(resultsFolder, "/to_initial_overall_female_cov",cov,".csv"), row.names = FALSE)
  
  lsc_initial_overall_female <- summariseLargeScaleCharacteristics(cdm$hpv_initial_overall_female_cohorts |>
                                                                   addDateOfBirth() |>
                                                                   mutate(cohort_start_date = date_of_birth + years(yy)), 
                                                                 window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                 eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                 minimumFrequency = 0
  )
  write.csv(lsc_initial_overall_female, paste0(resultsFolder, "/lsc_initial_overall_female_cov",cov,".csv"), row.names = FALSE)
  # After matching overall
  to_ratio_matched_overall_female <- to_Charac(cdm$hpv_ratio_matched_overall_female_cohorts, sex = "f", index_date = TRUE)
  write.csv(to_ratio_matched_overall_female, paste0(resultsFolder, "/to_ratio_matched_overall_female_cov",cov,".csv"), row.names = FALSE)
  
  lsc_ratio_matched_overall_female <- summariseLargeScaleCharacteristics(cdm$hpv_ratio_matched_overall_female_cohorts, 
                                                                   window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                   eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                   minimumFrequency = 0
  )
  write.csv(lsc_ratio_matched_overall_female, paste0(resultsFolder, "/lsc_ratio_matched_overall_female_cov",cov,".csv"), row.names = FALSE)
  
  # MALE
  yy <- 15
  ss <- "m"
  # Before matching overall
  to_initial_overall_male <- to_Charac(cdm$hpv_initial_overall_male_cohorts, sex = "m", index_date = FALSE)
  write.csv(to_initial_overall_male, file = paste0(resultsFolder, "/to_initial_overall_male_cov",cov,".csv"), row.names = FALSE)
  
  slc_initial_overall_male <- summariseLargeScaleCharacteristics(cdm$hpv_initial_overall_male_cohorts |>
                                                                   addDateOfBirth() |>
                                                                   mutate(cohort_start_date = date_of_birth + years(yy)), 
                                                                 window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                 eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                 minimumFrequency = 0
  )
  write.csv(lsc_initial_overall_male, paste0(resultsFolder, "/lsc_initial_overall_male_cov",cov,".csv"), row.names = FALSE)
  # After matching overall
  to_ratio_matched_overall_male <- to_Charac(cdm$hpv_ratio_matched_overall_male_cohorts, sex = "m", index_date = TRUE)
  write.csv(to_ratio_matched_overall_male, paste0(resultsFolder, "/to_ratio_matched_overall_male_cov",cov,".csv"), row.names = FALSE)
  
  lsc_ratio_matched_overall_male <- summariseLargeScaleCharacteristics(cdm$hpv_ratio_matched_overall_male_cohorts, 
                                                                         window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                         eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                         minimumFrequency = 0
  )
  write.csv(lsc_ratio_matched_overall_male, file = paste0(resultsFolder, "/lsc_ratio_matched_overall_male_cov",cov,".csv"), row.names = FALSE)
  
}


# -------------

cov <- 0

to_initial_overall_female <- read.csv(paste0(resultsFolder, "/to_initial_overall_female_cov",cov,".csv")) |> as_tibble()
slc_initial_overall_female <- read.csv(paste0(resultsFolder, "/lsc_initial_overall_female_cov",cov,".csv")) |> as_tibble()
to_ratio_matched_overall_female <- read.csv(paste0(resultsFolder, "/to_ratio_matched_overall_female_cov",cov,".csv")) |> as_tibble()
slc_ratio_matched_overall_female <- read.csv(paste0(resultsFolder, "/lsc_ratio_matched_overall_female_cov",cov,".csv")) |> as_tibble()

to_initial_overall_male <- read.csv(paste0(resultsFolder, "/to_initial_overall_male_cov",cov,".csv")) |> as_tibble()
lsc_initial_overall_male <- read.csv(paste0(resultsFolder, "/lsc_initial_overall_male_cov",cov,".csv")) |> as_tibble()
# to_ratio_matched_overall_male <- read.csv(paste0(resultsFolder, "/to_ratio_matched_overall_female.csv")) |> as_tibble()
# slc_ratio_matched_overall_male <- read.csv(paste0(resultsFolder, "/lsc_ratio_matched_overall_female.csv")) |> as_tibble()

# Figures and plots -------
# FEMALES
# Before matching overall
target <- paste0("vac_18f_coverage_",cov,"_any_dose")
comparator <- paste0("unvac_18f_coverage_",cov,"_any_dose")

smd_to_initial_overall_female <- SMD_to(to_initial_overall_female, target, comparator)
write.csv(smd_to_initial_overall_female |>
            dplyr::select(variable_name, variable_level, table, window, value, smd_initial = SMDs), 
          paste0(resultsFolder, "/smd_to_initial_overall_female_cov",cov,".csv"))

gtTable(smd_to_initial_overall_female,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial Overall Female Population Characterisation",
        colsToMergeRows = "all_columns"
)

smd_lsc_initial_overall_female <- SMD_lsc(lsc_initial_overall_female, target, comparator)

write.csv(smd_lsc_initial_overall_female |>
            dplyr::select(variable_name, variable_level, target, comparator, smd_initial = SMDs), 
          paste0(resultsFolder, "/smd_lsc_initial_overall_female_cov",cov,".csv"))
# Initial plots
bar_plot(smd_to_initial_overall_female, smd_lsc_initial_overall_female, target, comparator, "Vaccinated", "Unvaccinated", "initial_graphic_female")

# After matching
target <- paste0("matched_vac_cov_",cov,"_ratio_18f_any_dose")
comparator <- paste0("matched_unvac_cov_",cov,"_ratio_18f_any_dose")

smd_to_ratio_matched_overall_female <- SMD_to(to_ratio_matched_overall_female, target, comparator)
write.csv(smd_to_ratio_matched_overall_female |>
            dplyr::select(variable_name, variable_level, table, window, value, smd_matched = SMDs), 
          paste0(resultsFolder, "/smd_to_ratio_matched_overall_female_cov",cov,".csv"))

gtTable(smd_to_ratio_matched_overall_female,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Ratio Matched Overall Female Population Characterisation",
        colsToMergeRows = "all_columns"
)

smd_lsc_ratio_matched_overall_female <- SMD_lsc(lsc_ratio_matched_overall_female, target, comparator)

write.csv(smd_lsc_ratio_matched_overall_female |>
            dplyr::select(variable_name, variable_level, target, comparator, smd_matched = SMDs), 
          paste0(resultsFolder, "/smd_lsc_ratio_matched_overall_female_cov",cov,".csv"))
# Matched plots
bar_plot(smd_to_ratio_matched_overall_female, smd_lsc_ratio_matched_overall_female, target, comparator, "Vaccinated", "Unvaccinated", "ratio_matched_graphic_female")

# Comparison plots
lsc_overall_female_graphic <- inner_join(smd_lsc_ratio_matched_overall_female |> rename(ASMD_matched = SMDs), 
                                         smd_lsc_initial_overall_female |> rename(ASMD_initial = SMDs)) |>
  ggplot(aes(x = ASMD_initial, y = ASMD_matched)) +
  geom_point(color = "#226584") +  # Scatter plot
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "#002147") +  # Horizontal line
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#002147") +  # Vertical line
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "#002147") +  # 1:1 line
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "AMSD before matching", y = "AMSD after matching", title = "Female population characterisation") +
  geom_rect(aes(ymin = 0.1, ymax = 0.2, xmin = 0, xmax = 0.6), fill = "#fdbf69", alpha = 0.01) + 
  geom_rect(aes(xmin = 0.1, xmax = 0.6, ymin = 0, ymax = 0.2), fill = "#c7dbf9", alpha = 0.01)
  


# MALES
# Before matching overall
target <- paste0("vac_15m_coverage_",cov,"_any_dose")
comparator <- paste0("unvac_15m_coverage_",cov,"_any_dose")

smd_to_initial_overall_male <- SMD_to(to_initial_overall_male, target, comparator)
write.csv(smd_to_initial_overall_male |>
            dplyr::select(variable_name, variable_level, table, window, value, smd_initial = SMDs), 
          paste0(resultsFolder, "/smd_to_initial_overall_male_cov",cov,".csv"))

gtTable(smd_to_initial_overall_male,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial Overall Male Population Characterisation",
        colsToMergeRows = "all_columns"
)

smd_lsc_initial_overall_male <- SMD_lsc(lsc_initial_overall_male, target, comparator)


write.csv(smd_lsc_initial_overall_male |>
            dplyr::select(variable_name, variable_level, target, comparator, smd_initial = SMDs), 
          paste0(resultsFolder, "/smd_lsc_initial_overall_male_cov",cov,".csv"))

# Initial plots
bar_plot(smd_to_initial_overall_male, smd_lsc_initial_overall_male, target, comparator, "Vaccinated", "Unvaccinated", "initial_graphic_male")

# -------------------

smd_to_overall_female <- full_join(smd_to_initial_overall_female |> rename(ASMD_initial = SMDs), smd_to_ratio_matched_overall_female |> rename(ASMD_matched = SMDs)) |>
  dplyr::select(variable_name, variable_level, estimate_name, window, vac_18f_coverage_0_any_dose, unvac_18f_coverage_0_any_dose, matched_vac_cov_0_ratio_18f_any_dose, matched_unvac_cov_0_ratio_18f_any_dose, ASMD_initial, ASMD_matched)

gtTable(smd_to_overall_female,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial Overall Male Population Characterisation",
        colsToMergeRows = "all_columns"
)


# Let's distinguish according to coverage -----


if (ss == "f") {
  to_initial_overall <- to_initial_overall_female
  slc_initial_overall <- slc_initial_overall_female
  to_ratio_matched_overall <- to_ratio_matched_overall_female
  slc_ratio_matched_overall <- slc_ratio_matched_overall_female
  group <- "female"
} else {
  to_initial_overall <- to_initial_overall_male
  slc_initial_overall <- slc_initial_overall_male
  to_ratio_matched_overall <- to_ratio_matched_overall_male
  slc_ratio_matched_overall <- slc_ratio_matched_overall_male
  group <- "male"
  
}

# Past Characterisation of Initial cohorts (only eligible criteria applied): index date at 18 years ------
cdm <- omopgenerics::bind(
  cdm[[paste0("vac_15m_coverage_",cov,"_any_dose")]],
  cdm[[paste0("unvac_15m_coverage_",cov,"_any_dose")]],
  name = "hpv_initial_overall_male_cohorts"
)
cdm <- omopgenerics::bind(
  cdm[[paste0("vac_18f_coverage_",cov,"_any_dose")]],
  cdm[[paste0("unvac_18f_coverage_",cov,"_any_dose")]],
  name = "hpv_initial_overall_female_cohorts"
)


to_initial_overall_female <- cdm$hpv_initial_overall_female_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addDateOfBirth() |>
  mutate(cohort_start_date = date_of_birth + years(18)) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
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
      "cervical_screening" = list(
        targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
      ),
      "smear_test" = list(
        targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
      )
    ),
  )
write.csv(to_initial_overall_female, paste0(resultsFolder, "/to_initial_overall_female.csv"))

to_initial_overall_male <- cdm$hpv_initial_overall_male_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addDateOfBirth() |>
  mutate(cohort_start_date = date_of_birth + years(15)) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
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
      )
    ),
  )

write.csv(to_initial_overall_male, paste0(resultsFolder, "/to_initial_overall_male.csv"))

to_initial_overall_format <- formatTable(to_initial_overall,
                                         formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                                "N" = "<count>",
                                                                "Mean (SD)" = "<mean> (<sd>)"),
                                         header = c("group", "strata"),
                                         split = c("group","strata", "additional"),
                                         excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                         type = "gt")

# Compute table with SMDs and Plot
aux <- to_initial_overall |>
  as_tibble() |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")), ~ as.numeric(.x)/100)) |>
  mutate(across(c(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((!!sym(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")) - !!sym(paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")))/sqrt((!!sym(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"))*(1-!!sym(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"))) + !!sym(paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose"))*(1-!!sym(paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose"))))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_to_initial_overall <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_initial = SMDs)

write.csv(smd_to_initial_overall, paste0(resultsFolder, "/smd_to_initial_overall_", group,".csv"))

gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial Overall Population Characterisation",
        colsToMergeRows = "all_columns"
)

slc_initial_overall_male <- summariseLargeScaleCharacteristics(cdm$hpv_initial_overall_male_cohorts |>
                                                                 addDateOfBirth() |>
                                                                 mutate(cohort_start_date = date_of_birth + years(18)), 
                                          window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                          eventInWindow = c("condition_occurrence","drug_exposure"),
                                          minimumFrequency = 0
)

write.csv(slc_initial_overall_male, paste0(resultsFolder, "/lsc_initial_cov",cov,"_overall_male.csv"))

aux_slc <- slc_initial_overall |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version, overall)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")), ~ as.numeric(.x)/100)) |>
  mutate(across(c(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((!!sym(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")) - !!sym(paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")))/sqrt((!!sym(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"))*(1-!!sym(paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"))) + !!sym(paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose"))*(1-!!sym(paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose"))))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_lsc_initial_overall <- aux_slc |>
  dplyr::select(variable_name, variable_level, paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose"), smd_initial = SMDs)

write.csv(smd_lsc_initial_overall, paste0(resultsFolder, "/smd_lsc_initial_overall_", group,".csv"))

# Compute graphics
# Specific hpv characteristics: hiv, smear_test, screening, (previous vaccines), autoimmune disease
# General health characteristics (>0.01 ASMD): anxiety, hormonal contraceptives, antibacterials sys, penumonia (-Inf,-1)
# General health characteristics (<0.01 ASMD): Opioids, lipic modifying agents, chronic liver disease, diabetes I

# Combining SLC
initial_graphic <- rbind(
  # specific
  aux |> filter(variable_name %in% c("Hiv status", "Papanicolau smear testing")) |> mutate(type = "specific", variable_name = "Test") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
  aux |> filter(variable_name == "cervical_screening") |> mutate(type = "specific") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
  aux |> filter(variable_level == "Autoimmune disease") |> mutate(type = "specific") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
  aux |> filter(variable_level == "immunosupressants", window == "-365 to -1") |> mutate(type = "specific") |> dplyr::select(! c(variable_name, table, value)) |> rename(variable_name = variable_level),
  # general > 0.1
  aux_slc |> filter(variable_name %in% c("Conjunctivitis", "Pain in throat", "Eruption"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id),
  aux_slc |> filter(variable_name %in% c("tetanus toxoid vaccine, inactivated", "measles virus vaccine", "ibuprofen 20 MG/ML Oral Suspension"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id),
  # general < 0.1
  aux_slc |> filter(variable_name %in% c("Anxiety", "Cough", "Asthma daytime symptoms"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id),
  aux_slc |> filter(variable_name %in% c("amoxicillin 500 MG Oral Capsule", "penicillin V potassium 250 MG Oral Tablet", "Benzydamine 1.5 MG/ML Nasal Spray"), variable_level == "-inf to -366") |> mutate(type = "general") |> rename(window = variable_level) |> dplyr::select(! concept_id)
) |>
  dplyr::select(! estimate_name)

vac_variable <- paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")
unvac_variable <- paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")

initial_graphic_m0.1 <- initial_graphic |>
  filter(SMDs <= 0.1) |>
  dplyr::select(variable_name, paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")) |>
  reshape2::melt(id.vars = "variable_name") |>
  ggplot(aes(x = variable_name, y = value, fill = variable)) + 
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = setNames(c("#7EB3EB", "#FDBF69"), c(vac_variable, unvac_variable)), 
                    labels = c("Vaccinated", "Unvaccinated")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Initial characterisation: ASMD < 0.1", x = "Characteristics", y = "Proportion", fill = "Cohorts")

initial_graphic_p0.1 <- initial_graphic |>
  filter(SMDs > 0.1) |>
  dplyr::select(variable_name, paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose"), paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")) |>
  reshape2::melt(id.vars = "variable_name") |>
  ggplot(aes(x = variable_name, y = value, fill = variable)) + 
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = setNames(c("#7EB3EB", "#FDBF69"), c(vac_variable, unvac_variable)), 
                    labels = c("Vaccinated", "Unvaccinated")) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Initial characterisation: ASMD > 0.1", x = "Characteristics", y = "Proportion", fill = "Cohorts")

# -----------------------

# Past Characterisation of vaccinated population by num of doses (only eligible criteria applied and dose separation): index date at 18 years and per dose count --------

to_initial_ndose <- cdm[[paste0("vac_18f_coverage_",cov,"_any_dose")]] |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addDateOfBirth() |>
  mutate(cohort_start_date = date_of_birth + years(18)) |>
  addCohortIntersectCount(
    targetCohortTable = "doses_allvac_cohort",
    nameStyle = "dose_count",
    indexDate = "cohort_start_date",
    window = c(-Inf,0)
  ) |>
  filter(dose_count < 4) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
    strata = list("dose_count"),
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
      "cervical_screening" = list(
        targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
      ),
      "smear_test" = list(
        targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
      )
    ),
  )

to_initial_ndose_format <- formatTable(to_initial_ndose,
                                       formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                              "N" = "<count>",
                                                              "Mean (SD)" = "<mean> (<sd>)"),
                                       header = c("group", "strata"),
                                       split = c("group","strata", "additional"),
                                       excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                       type = "gt")

# -----------------------

# Past Characterisation of Initial Dose cohorts (only eligible criteria and dose separation) -----
cdm <- omopgenerics::bind(
  cdm$vac_18f_coverage_0_1dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "1vs23dose"),
  cdm$vac_18f_coverage_0_23dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "1vs23dose"),
  cdm$vac_18f_coverage_0_2dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "2vs3dose"),
  cdm$vac_18f_coverage_0_3dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "2vs3dose"),
  name = "hpv_notmatched_dose_cohorts"
)

comparisons <- c("1vs23dose", "2vs3dose")
for (comp in comparisons) {
  to_doses <- cdm$hpv_notmatched_dose_cohorts |>
    filter(comparison == comp) |> 
    dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
    addDateOfBirth() |>
    mutate(cohort_start_date = date_of_birth + years(18)) |>
    mutate(cohort_start_date = as.Date(cohort_start_date)) |>
    mutate(cohort_end_date = as.Date(cohort_end_date)) |>
    summariseCharacteristics(
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
        "cervical_screening" = list(
          targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
        ),
        "smear_test" = list(
          targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
        )
      ),
    )
  assign(paste0("to_notmatched_dose_",comp),to_doses)
  
  to_doses_format <- formatTable(to_doses,
                                 formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                        "N" = "<count>",
                                                        "Mean (SD)" = "<mean> (<sd>)"),
                                 header = c("group", "strata"),
                                 split = c("group","strata", "additional"),
                                 excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                 type = "gt")
  assign(paste0("to_notmatched_dose_",comp,"_format"), to_doses_format)
}

aux <- to_notmatched_dose_1vs23dose |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(vac_18f_coverage_0_23dose, vac_18f_coverage_0_1dose), ~ as.numeric(.x)/100)) |>
  mutate(across(c(vac_18f_coverage_0_23dose, vac_18f_coverage_0_1dose), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((vac_18f_coverage_0_23dose - vac_18f_coverage_0_1dose)/sqrt((vac_18f_coverage_0_23dose*(1-vac_18f_coverage_0_23dose) + vac_18f_coverage_0_1dose*(1-vac_18f_coverage_0_1dose))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_initial_1vs23 <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_initial = SMDs)


gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial 1vs23 Dose Population Characterisation",
        colsToMergeRows = "all_columns"
)

aux <- to_notmatched_dose_2vs3dose |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(vac_18f_coverage_0_3dose, vac_18f_coverage_0_2dose), ~ as.numeric(.x)/100)) |>
  mutate(across(c(vac_18f_coverage_0_3dose, vac_18f_coverage_0_2dose), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((vac_18f_coverage_0_3dose - vac_18f_coverage_0_2dose)/sqrt((vac_18f_coverage_0_3dose*(1-vac_18f_coverage_0_3dose) + vac_18f_coverage_0_2dose*(1-vac_18f_coverage_0_2dose))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_initial_2vs3 <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_initial = SMDs)


gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial 2vs3 Dose Population Characterisation",
        colsToMergeRows = "all_columns"
)
# -----------------

# Past Characterisation of overall matched cohorts (index date is date of first vac or pair date) -------
cdm <- omopgenerics::bind(
  cdm[[paste0("matched_vaccinated_coverage_",cov,"_18f_any_dose")]],
  cdm[[paste0("matched_unvaccinated_coverage_",cov,"_18f_any_dose")]],
  name = "hpv_matched_overall_cohorts"
)

to_matched_overall <- cdm$hpv_matched_overall_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
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
      "cervical_screening" = list(
        targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
      ),
      "smear_test" = list(
        targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
      )
    ),
  )
write.csv(to_matched_overall, paste0(resultsFolder, "/to_matched_overall.csv"))
to_matched_overall_format <- formatTable(to_matched_overall,
                                         formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                                "N" = "<count>",
                                                                "Mean (SD)" = "<mean> (<sd>)"),
                                         header = c("group", "strata"),
                                         split = c("group","strata", "additional"),
                                         excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                         type = "gt")

aux <- to_matched_overall |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose), ~ as.numeric(.x)/100)) |>
  mutate(across(c(matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((matched_vaccinated_coverage_0_18f_any_dose - matched_unvaccinated_coverage_0_18f_any_dose)/sqrt((matched_vaccinated_coverage_0_18f_any_dose*(1-matched_vaccinated_coverage_0_18f_any_dose) + matched_unvaccinated_coverage_0_18f_any_dose*(1-matched_unvaccinated_coverage_0_18f_any_dose))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! c(estimate_type))

smd_matched_overall <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_matched = SMDs)

gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Matched Overall Population Characterisation",
        colsToMergeRows = "all_columns"
)

# Large Scale Characterisation
cdm <- omopgenerics::bind(
  cdm$unvac_18f_coverage_0_any_dose |> mutate(treatment = 0, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(18)) |> dplyr::select(!date_of_birth),
  cdm$vac_18f_coverage_0_any_dose |> mutate(treatment = 1, group = "initial_overall") |> addDateOfBirth() |> mutate(cohort_start_date = date_of_birth + years(18)) |> dplyr::select(!date_of_birth),
  cdm$matched_unvaccinated_coverage_0_18f_any_dose |> mutate(treatment = 0, group = "matched_overall"),
  cdm$matched_vaccinated_coverage_0_18f_any_dose |> mutate(treatment = 1, group = "matched_overall"),
  name = "hpv_matchedvsinitial_cohorts"
)
slc <- summariseLargeScaleCharacteristics(cdm$hpv_matchedvsinitial_cohorts, 
                                          window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                          eventInWindow = c("condition_occurrence","drug_exposure"),
                                          minimumFrequency = 0
)

write.csv(slc, paste0(resultsFolder, "/lsc_matched_overall.csv"))

slc_comp <- slc
slc_matchedvsinitial_overall <- slc_comp |>
  splitAll() |>
  filter(estimate_type == "percentage") |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_type, package_name, package_version, estimate_name, estimate_type)) |>
  mutate(across(! c(result_id, variable_name, variable_level, concept_id), ~ as.numeric(.x)/100)) |>
  mutate(across(! c(result_id, variable_name, variable_level, concept_id), ~ replace(., is.na(.), 0))) |>
  mutate(smd_initial = abs((vac_18f_coverage_0_any_dose - unvac_18f_coverage_0_any_dose)/sqrt((vac_18f_coverage_0_any_dose*(1-vac_18f_coverage_0_any_dose) + unvac_18f_coverage_0_any_dose*(1-unvac_18f_coverage_0_any_dose))/2))) |>
  mutate(smd_matched = abs((matched_vaccinated_coverage_0_18f_any_dose - matched_unvaccinated_coverage_0_18f_any_dose)/sqrt((matched_vaccinated_coverage_0_18f_any_dose*(1-matched_vaccinated_coverage_0_18f_any_dose) + matched_unvaccinated_coverage_0_18f_any_dose*(1-matched_unvaccinated_coverage_0_18f_any_dose))/2)))
  # mutate(across(! c(result_id, variable_name, variable_level, concept_id), ~ replace(., is.nan(.), 0)))
  
# Compute graphics
# Specific hpv characteristics: hiv, smear_test, screening, (previous vaccines), autoimmune disease
# General health characteristics (>0.01 ASMD): anxiety, hormonal contraceptives, antibacterials sys, penumonia (-Inf,-1)
# General health characteristics (<0.01 ASMD): Opioids, lipic modifying agents, chronic liver disease, diabetes I

matched_graphic <- rbind(
  # specific
  aux |> filter(variable_name %in% c("Hiv status", "Papanicolau smear testing")) |> mutate(type = "specific", variable_name = "Test"),
  #aux |> filter(variable_name == "cervical_screening") |> mutate(type = "specific")|>
  aux |> filter(variable_level == "Autoimmune disease") |> mutate(type = "specific"),
  aux |> filter(variable_level == "immunosupressants", window == "-365 to -1") |> mutate(type = "specific"),
  # general > 0.01
  aux |> filter(variable_level %in% c("Anxiety", "Asthma", "Pneumonia"), window == "-inf to -1") |> mutate(type = "general"),
  aux |> filter(variable_level %in% c("Hormonal contraceptices sys", "Antibacterials sys", "Beta blocking agents"), window == "-365 to -1") |> mutate(type = "general"),
  # general < 0.01
  aux |> filter(variable_level %in% c("Chronic liver disease", "Chronic kidney disease", "Inflammatory bowel disease"), window == "-inf to -1") |> mutate(type = "general"),
  aux |> filter(variable_level %in% c("Opioids", "Lipid modifying agents", "Diuretics"), window == "-365 to -1") |> mutate(type = "general")
) |>
  dplyr::select(variable_name, variable_level, matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose, SMDs)

matched_graphic_m0.01 <- matched_graphic |>
  filter(SMDs <= 0.01) |>
  dplyr::select(variable_level, matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose) |>
  reshape2::melt(id.vars = "variable_level") |>
  ggplot(aes(x = variable_level, y = value, fill = variable)) + 
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = c("matched_vaccinated_coverage_0_18f_any_dose" = "#0DB6C9", "matched_unvaccinated_coverage_0_18f_any_dose" = "#FF6123"), labels = c("Matched vaccinated", "Matched unvaccinated")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Characterisation after matching: ASMD < 0.01", x = "Characteristics", y = "Proportion", fill = "Cohorts")

matched_graphic_p0.01 <- matched_graphic |>
  filter(SMDs > 0.01) |>
  dplyr::select(variable_level, matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose) |>
  reshape2::melt(id.vars = "variable_level") |>
  ggplot(aes(x = variable_level, y = value, fill = variable)) + 
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = c("matched_vaccinated_coverage_0_18f_any_dose" = "#0DB6C9", "matched_unvaccinated_coverage_0_18f_any_dose" = "#FF6123"), labels = c("Matched vaccinated", "Matched unvaccinated")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Characterisation after matching: ASMD > 0.01", x = "Characteristics", y = "Proportion", fill = "Cohorts")

# SMD comparison graphic:
to_smd_graphic <- slc_matchedvsinitial_overall |> 
  inner_join(smd_matched_overall, by = c("variable_name", "variable_level", "table", "window", "value")) |>
  ggplot(aes(x = smd_initial, y = smd_matched)) +
  geom_point() +  # Scatter plot
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +  # Horizontal line
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # 1:1 line
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "AMSD Before matching", y = "AMSD after matching", title = "Population balance comparison") +
  xlim(0,0.12) +
  ylim(0,0.12)

lsc_smd_graphic <- slc_matchedvsinitial_overall |> 
  ggplot(aes(x = smd_initial, y = smd_matched)) +
  geom_point() +  # Scatter plot
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +  # Horizontal line
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # 1:1 line
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "AMSD Before matching", y = "AMSD after matching", title = "Population balance comparison") +
  xlim(0,0.6) +
  ylim(0,0.6)
# -----------------

# Past Characterisation of overall ratio 1:5 matched cohorts (index date is date of first vac or pair date) -------
cdm <- omopgenerics::bind(
  cdm[[paste0("matched_vac_cov_",cov,"_ratio_18f_any_dose")]],
  cdm[[paste0("matched_unvac_cov_",cov,"_ratio_18f_any_dose")]],
  name = "hpv_ratio_matched_overall_female_cohorts"
)

target <- paste0("matched_vac_cov_",cov,"_ratio_18f_any_dose")
comparator <- paste0("matched_unvac_cov_",cov,"_ratio_18f_any_dose")

to_ratio_matched_overall_female <- cdm$hpv_ratio_matched_overall_female_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
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
      "cervical_screening" = list(
        targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
      ),
      "smear_test" = list(
        targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
      )
    ),
  )

write.csv(to_ratio_matched_overall_female, paste0(resultsFolder, "/to_ratio_matched_overall_female.csv"))

to_ratio_matched_overall_format <- formatTable(to_ratio_matched_overall,
                                         formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                                "N" = "<count>",
                                                                "Mean (SD)" = "<mean> (<sd>)"),
                                         header = c("group", "strata"),
                                         split = c("group","strata", "additional"),
                                         excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                         type = "gt")

aux <- to_ratio_matched_overall |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(all_of(c(target, comparator)), ~ as.numeric(.x)/100)) |>
  mutate(across(all_of(c(target, comparator)), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((!!sym(target) - !!sym(comparator))/sqrt((!!sym(target)*(1-!!sym(target)) + !!sym(comparator)*(1-!!sym(comparator)))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_to_ratio_matched_overall <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_initial = SMDs)

write.csv(smd_to_ratio_matched_overall, paste0(resultsFolder, "/smd_to_ratio_matched_overall_", group,".csv"))

gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Initial Overall Population Characterisation",
        colsToMergeRows = "all_columns"
)


# Large Scale Characterisation

slc_ratio_matched_overall_female <- summariseLargeScaleCharacteristics(cdm$hpv_ratio_matched_overall_female_cohorts, 
                                          window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                          eventInWindow = c("condition_occurrence","drug_exposure"),
                                          minimumFrequency = 0
)
write.csv(slc_ratio_matched_overall_female, paste0(resultsFolder, "/lsc_ratio_matched_overall_female.csv"))

slc_ratio_matched_overall_male <- summariseLargeScaleCharacteristics(cdm$hpv_ratio_matched_overall_male_cohorts, 
                                                                       window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                       eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                       minimumFrequency = 0
)
write.csv(slc, paste0(resultsFolder, "/lsc_ratio_matched_overall_male.csv"))

aux_slc <- slc_ratio_matched_overall |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version, overall)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(target, comparator), ~ as.numeric(.x)/100)) |>
  mutate(across(c(target, comparator), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((!!sym(target) - !!sym(comparator))/sqrt((!!sym(target)*(1-!!sym(target)) + !!sym(comparator)*(1-!!sym(comparator)))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_lsc_ratio_matched_overall <- aux_slc |>
  dplyr::select(variable_name, variable_level, target, comparator, smd_initial = SMDs)

write.csv(smd_lsc_ratio_matched_overall, paste0(resultsFolder, "/smd_lsc_ratio_matched_overall_", group,".csv"))

slc_ratio_matched_overall <- slc |>
  splitAll() |>
  filter(estimate_type == "percentage") |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_type, package_name, package_version, estimate_name, estimate_type)) |>
  mutate(across(! c(result_id, variable_name, variable_level, concept_id), ~ as.numeric(.x)/100)) |>
  mutate(across(! c(result_id, variable_name, variable_level, concept_id), ~ replace(., is.na(.), 0))) |>
  mutate(smd_ratio_matched = abs((matched_vaccinated_coverage_0_ratio_18f_any_dose - matched_unvaccinated_coverage_0_ratio_18f_any_dose)/sqrt((matched_vaccinated_coverage_0_ratio_18f_any_dose*(1-matched_vaccinated_coverage_0_ratio_18f_any_dose) + matched_unvaccinated_coverage_0_ratio_18f_any_dose*(1-matched_unvaccinated_coverage_0_ratio_18f_any_dose))/2)))
# mutate(across(! c(result_id, variable_name, variable_level, concept_id), ~ replace(., is.nan(.), 0)))

# Compute graphics
# Specific hpv characteristics: hiv, smear_test, screening, (previous vaccines), autoimmune disease
# General health characteristics (>0.01 ASMD): anxiety, hormonal contraceptives, antibacterials sys, penumonia (-Inf,-1)
# General health characteristics (<0.01 ASMD): Opioids, lipic modifying agents, chronic liver disease, diabetes I

matched_graphic <- rbind(
  # specific
  aux |> filter(variable_name %in% c("Hiv status", "Papanicolau smear testing", "cervical_screening")) |> mutate(type = "specific", variable_name = "Test"),
  aux |> filter(variable_level == "Autoimmune disease") |> mutate(type = "specific"),
  aux |> filter(variable_level == "immunosupressants", window == "-365 to -1") |> mutate(type = "specific"),
  # general > 0.01
  aux |> filter(variable_level %in% c("Anxiety", "Asthma", "Pneumonia"), window == "-inf to -1") |> mutate(type = "general"),
  aux |> filter(variable_level %in% c("Hormonal contraceptices sys", "Antibacterials sys", "Beta blocking agents"), window == "-365 to -1") |> mutate(type = "general"),
  # general < 0.01
  aux |> filter(variable_level %in% c("Chronic liver disease", "Chronic kidney disease", "Inflammatory bowel disease"), window == "-inf to -1") |> mutate(type = "general"),
  aux |> filter(variable_level %in% c("Opioids", "Lipid modifying agents", "Diuretics"), window == "-365 to -1") |> mutate(type = "general")
) |>
  dplyr::select(variable_name, variable_level, matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose, SMDs)

ratio_matched_graphic_m0.01 <- ratio_matched_graphic |>
  filter(SMDs <= 0.01) |>
  dplyr::select(variable_level, matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose) |>
  reshape2::melt(id.vars = "variable_level") |>
  ggplot(aes(x = variable_level, y = value, fill = variable)) + 
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = c("matched_vaccinated_coverage_0_18f_any_dose" = "#0DB6C9", "matched_unvaccinated_coverage_0_18f_any_dose" = "#FF6123"), labels = c("Matched vaccinated", "Matched unvaccinated")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Characterisation after matching: ASMD < 0.01", x = "Characteristics", y = "Proportion", fill = "Cohorts")

ratio_matched_graphic_p0.01 <- ratio_matched_graphic |>
  filter(SMDs > 0.01) |>
  dplyr::select(variable_level, matched_vaccinated_coverage_0_18f_any_dose, matched_unvaccinated_coverage_0_18f_any_dose) |>
  reshape2::melt(id.vars = "variable_level") |>
  ggplot(aes(x = variable_level, y = value, fill = variable)) + 
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = c("matched_vaccinated_coverage_0_18f_any_dose" = "#0DB6C9", "matched_unvaccinated_coverage_0_18f_any_dose" = "#FF6123"), labels = c("Matched vaccinated", "Matched unvaccinated")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Characterisation after matching: ASMD > 0.01", x = "Characteristics", y = "Proportion", fill = "Cohorts")

# SMD comparison graphic:
to_smd_graphic <- slc_matchedvsinitial_overall |> 
  inner_join(smd_matched_overall, by = c("variable_name", "variable_level", "table", "window", "value")) |>
  ggplot(aes(x = smd_initial, y = smd_matched)) +
  geom_point() +  # Scatter plot
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +  # Horizontal line
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # 1:1 line
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "AMSD Before matching", y = "AMSD after matching", title = "Population balance comparison") +
  xlim(0,0.12) +
  ylim(0,0.12)

lsc_smd_graphic <- slc_matchedvsinitial_overall |> 
  ggplot(aes(x = smd_initial, y = smd_matched)) +
  geom_point() +  # Scatter plot
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +  # Horizontal line
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # 1:1 line
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "AMSD Before matching", y = "AMSD after matching", title = "Population balance comparison") +
  xlim(0,0.6) +
  ylim(0,0.6)
# -----------------

# Past Characterisation of dose matched cohorts (index date = first dose) ---------
cdm <- omopgenerics::bind(
  cdm[[paste0("matched_vac_18f_coverage_",cov,"_1dose")]] |> mutate(comparison = "1vs23dose"),
  cdm[[paste0("matched_vac_18f_coverage_",cov,"_23dose")]] |> mutate(comparison = "1vs23dose"),
  cdm[[paste0("matched_vac_18f_coverage_",cov,"_2dose")]] |> mutate(comparison = "2vs3dose"),
  cdm[[paste0("matched_vac_18f_coverage_",cov,"_3dose")]] |> mutate(comparison = "2vs3dose"),
  name = "hpv_matched_dose_cohorts"
)
comparisons <- c("1vs23dose", "2vs3dose")

for (comp in comparisons) {
  to_doses <- cdm$hpv_matched_dose_cohorts |>
    filter(comparison == comp) |> 
    dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
    mutate(cohort_start_date = as.Date(cohort_start_date)) |>
    mutate(cohort_end_date = as.Date(cohort_end_date)) |>
    summariseCharacteristics(
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
        "cervical_screening" = list(
          targetCohortTable = "cervical_screening", value = "count", window = c(-Inf,-1)
        ),
        "smear_test" = list(
          targetCohortTable = "smear_test", value = "flag", window = c(-Inf,-1)
        )
      ),
    )
  assign(paste0("to_matched_dose_",comp),to_doses)
  
  
 to_doses_format <- formatTable(to_doses,
                               formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                      "N" = "<count>",
                                                      "Mean (SD)" = "<mean> (<sd>)"),
                               header = c("group", "strata"),
                               split = c("group","strata", "additional"),
                               excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                               type = "gt")
  assign(paste0("to_matched_dose_",comp,"_format"), to_doses_format)
  
}

aux <- to_matched_dose_1vs23dose |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(matched_vac_18f_coverage_0_23dose, matched_vac_18f_coverage_0_1dose), ~ as.numeric(.x)/100)) |>
  mutate(across(c(matched_vac_18f_coverage_0_23dose, matched_vac_18f_coverage_0_1dose), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((matched_vac_18f_coverage_0_23dose - matched_vac_18f_coverage_0_1dose)/sqrt((matched_vac_18f_coverage_0_23dose*(1-matched_vac_18f_coverage_0_23dose) + matched_vac_18f_coverage_0_1dose*(1-matched_vac_18f_coverage_0_1dose))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_matched_1vs23 <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_matched = SMDs)


gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Matched 1vs23 Dose Population Characterisation",
        colsToMergeRows = "all_columns"
)

aux <- to_matched_dose_2vs3dose |>
  visOmopResults::splitAll() |>
  tidyr::pivot_wider(names_from = "cohort_name", values_from = "estimate_value") |>
  dplyr::select(! c(cdm_name, result_id, result_type, package_name, package_version)) |> 
  dplyr::mutate(across(! c(variable_name, variable_level, estimate_name, estimate_type, table), ~replace(., is.na(.), 0))) |>
  filter(estimate_name == "percentage") |>
  filter(! variable_name == "Sex") |>
  mutate(across(c(matched_vac_18f_coverage_0_3dose, matched_vac_18f_coverage_0_2dose), ~ as.numeric(.x)/100)) |>
  mutate(across(c(matched_vac_18f_coverage_0_3dose, matched_vac_18f_coverage_0_2dose), ~ as.numeric(.x))) |>
  mutate(SMDs = abs((matched_vac_18f_coverage_0_3dose - matched_vac_18f_coverage_0_2dose)/sqrt((matched_vac_18f_coverage_0_3dose*(1-matched_vac_18f_coverage_0_3dose) + matched_vac_18f_coverage_0_2dose*(1-matched_vac_18f_coverage_0_2dose))/2))) |>
  arrange(desc(SMDs)) |>
  dplyr::select(! estimate_type)

smd_matched_2vs3 <- aux |>
  dplyr::select(variable_name, variable_level, table, window, value, smd_matched = SMDs)


gtTable(aux,
        delim = "\n",
        style = "default",
        na = "-",
        title = "Matched 2vs3 Dose Population Characterisation",
        colsToMergeRows = "all_columns"
)
# -----------------------

# Future Characterisation of Initial cohorts (only eligible criteria applied): index date at 18 years ------
cdm <- omopgenerics::bind(
  cdm[[paste0("vac_18f_coverage_",cov,"_any_dose")]],
  cdm[[paste0("unvac_18f_coverage_",cov,"_any_dose")]],
  name = "hpv_initial_overall_cohorts"
)
cdm <- omopgenerics::bind(
  cdm[[paste0("vac_18m_cohort")]],
  cdm[[paste0("unvac_18m_cohort")]],
  name = "hpv_initial_overall_cohorts"
  )
# cdm <- omopgenerics::bind(
#   cdm[[paste0("vac_18f_cohort_mod")]],
#   cdm[[paste0("unvac_18f_cohort_mod")]],
#   name = "hpv_initial_overall_cohorts"
# )
future_to_initial_overall <- cdm$hpv_initial_overall_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addDateOfBirth() |>
  mutate(cohort_start_date = date_of_birth + years(18)) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  addAge(cdm,ageGroup = list(c(0, 14), c(15, 16), c(17, 18))) |>
  dplyr::select(! c(date_of_birth, age)) |>
  summariseCharacteristics(
    # strata = list("age_group"),
    cohortIntersect = list(
      "HIV_status" = list(
        targetCohortTable = "hiv_status", value = "flag", window = list(c(1,Inf))
      ),
      "cervical_screening" = list(
        targetCohortTable = "cervical_screening", value = "count", window = c(1, Inf)
      ),
      "outcomes" = list(
        targetCohortTable = "nondarwin_outcome_cohorts", value = "flag", window = list(c(1, Inf))
      )
    )
  )

future_to_initial_overall_format <- formatTable(future_to_initial_overall,
                                                formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                                       "N" = "<count>",
                                                                       "Mean (SD)" = "<mean> (<sd>)"),
                                                header = c("group", "strata"),
                                                split = c("group","strata", "additional"),
                                                excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                                type = "gt")

# -----------------------

# Future Characterisation of vaccinated population by num of doses (only eligible criteria applied and dose separation): index date at 18 years and per dose count --------
# paste0("vac_18f_coverage_",cov,"_any_dose")
future_to_initial_ndose <- cdm[["vac_18m_cohort"]] |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addDateOfBirth() |>
  mutate(cohort_start_date = date_of_birth + years(18)) |>
  addCohortIntersectCount(
    targetCohortTable = "doses_allvac_cohort",
    nameStyle = "dose_count",
    indexDate = "cohort_start_date",
    window = c(-Inf,0)
  ) |>
  #filter(dose_count < 4) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  summariseCharacteristics(
    strata = list("dose_count"),
    cohortIntersect = list(
      "HIV_status" = list(
        targetCohortTable = "hiv_status", value = "flag", window = c(1,Inf)
      ),
      # "cervical_screening" = list(
      #   targetCohortTable = "cervical_screening", value = "count", window = c(1, Inf)
      # ),
      "outcomes" = list(
        targetCohortTable = "outcome_cohorts", value = "flag", window = c(1, Inf)
      )
    ),
  )

future_to_initial_ndose_format <- formatTable(future_to_initial_ndose,
                                              formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                                     "N" = "<count>",
                                                                     "Mean (SD)" = "<mean> (<sd>)"),
                                              header = c("group", "strata"),
                                              split = c("group","strata", "additional"),
                                              excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                              type = "gt")

# -----------------------

# Future Characterisation of Initial Dose cohorts (only eligible criteria and dose separation) -----
cdm <- omopgenerics::bind(
  cdm$vac_18f_coverage_0_1dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "1vs23dose"),
  cdm$vac_18f_coverage_0_23dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "1vs23dose"),
  cdm$vac_18f_coverage_0_2dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "2vs3dose"),
  cdm$vac_18f_coverage_0_3dose |> dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |> mutate(comparison = "2vs3dose"),
  name = "hpv_notmatched_dose_cohorts"
)

comparisons <- c("1vs23dose", "2vs3dose")
for (comp in comparisons) {
  to_doses <- cdm$hpv_notmatched_dose_cohorts |>
    filter(comparison == comp) |> 
    dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
    addDateOfBirth() |>
    mutate(cohort_start_date = date_of_birth + years(18)) |>
    mutate(cohort_start_date = as.Date(cohort_start_date)) |>
    mutate(cohort_end_date = as.Date(cohort_end_date)) |>
    summariseCharacteristics(
      cohortIntersect = list(
        "HIV_status" = list(
          targetCohortTable = "hiv_status", value = "flag", window = c(1, Inf)
        ),
        "cervical_screening" = list(
          targetCohortTable = "cervical_screening", value = "count", window = c(1, Inf)
        ),
        "outcomes" = list(
          targetCohortTable = "outcome_cohorts", value = "flag", window = c(1, Inf)
        )
      ),
    )
  assign(paste0("future_to_notmatched_dose_",comp),to_doses)
  
  to_doses_format <- formatTable(to_doses,
                                 formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                        "N" = "<count>",
                                                        "Mean (SD)" = "<mean> (<sd>)"),
                                 header = c("group", "strata"),
                                 split = c("group","strata", "additional"),
                                 excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                 type = "gt")
  assign(paste0("future_to_notmatched_dose_",comp,"_format"), to_doses_format)
}
# -----------------

# Future Characterisation of Matched cohorts (only eligible criteria applied): index date at 18 years ------
cdm <- omopgenerics::bind(
  cdm[[paste0("matched_vaccinated_coverage_",cov,"_18f_any_dose")]],
  cdm[[paste0("matched_unvaccinated_coverage_",cov,"_18f_any_dose")]],
  name = "hpv_matched_overall_cohorts"
)

future_to_matched_overall <- cdm$hpv_matched_overall_cohorts |> 
  dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
  addAge(cdm,ageGroup = list(c(0, 14), c(15, 16), c(17, 18))) |>
  mutate(cohort_start_date = as.Date(cohort_start_date)) |>
  mutate(cohort_end_date = as.Date(cohort_end_date)) |>
  addAge(cdm,ageGroup = list(c(0, 14), c(15, 16), c(17, 18))) |>
  summariseCharacteristics(
    strata = list("age_group"),
    cohortIntersect = list(
      "HIV_status" = list(
        targetCohortTable = "hiv_status", value = "flag", window = list(c(1,Inf))
      ),
      # "cervical_screening" = list(
      #   targetCohortTable = "cervical_screening", value = "count", window = c(1, Inf)
      # ),
      "outcomes" = list(
        targetCohortTable = "nondarwin_outcome_cohorts", value = "flag", window = list(c(1, Inf))
      )
    )
  )

future_to_matched_overall_format <- formatTable(future_to_matched_overall,
                                                formatEstimateName = c("N%" = "<count> (<percentage>)",
                                                                       "N" = "<count>",
                                                                       "Mean (SD)" = "<mean> (<sd>)"),
                                                header = c("group", "strata"),
                                                split = c("group","strata", "additional"),
                                                excludeColumns = c("cdm_name", "result_id", "result_type", "package_name", "package_version", "estimate_type"),
                                                type = "gt")

# -----------------------

# Outcomes review
dates <- cdm$outcome_cohorts |> 
  filter(! cohort_definition_id == 6) |>
  filter(cohort_start_date > as.Date("2009-01-01")) |>
  addTableIntersectFlag(tableName = "doses_allvac_cohort",
                        indexDate = "cohort_start_date",
                        window = c(-Inf, 1),
                        nameStyle = "vaccination_flag"
  ) |>
  group_by(cohort_definition_id, vaccination_flag) |>
  tally() |>
  left_join(settings(cdm$outcome_cohorts) |> select(cohort_definition_id, cohort_name), by = "cohort_definition_id", copy = TRUE)

write.csv(dates, paste0(resultsFolder, "/outcomes_numbers_instudyperiod.csv"))

addAge() |> 
  mutate(year = year(cohort_start_date)) |>
  #filter(age < 33) |>
  pull(year) |> 
  table()
barplot(dates,
        widths = 2,
        las = 2)
# -----------