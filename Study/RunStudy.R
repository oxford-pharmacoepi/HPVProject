# Results folders
resultsFolder <- here("Results_June")
resultsFolder_lasso <- here(resultsFolder, "lasso")
resultsFolder_att <- here(resultsFolder, "attritions")
resultsFolder_char <- here(resultsFolder, "characterisation")
resultsFolder_nco_fem <- here(resultsFolder, "nco_female")
resultsFolder_nco_male <- here(resultsFolder, "nco_male")
resultsFolder_outcome_fem_darwin <- here(resultsFolder, "outcome_model_female_darwin")
resultsFolder_outcome_fem_nondarwin <- here(resultsFolder, "outcome_model_female_nondarwin")
resultsFolder_outcome_male <- here(resultsFolder, "outcome_model_male")

subDirs <- c(resultsFolder_lasso, resultsFolder_att, resultsFolder_char, resultsFolder_nco_fem, 
             resultsFolder_nco_male, resultsFolder_outcome_fem_darwin, resultsFolder_outcome_fem_nondarwin,
             resultsFolder_outcome_male)
for (dir in subDirs) {
  if (! file.exists(dir)){
    dir.create(dir)
  }
}

# create logger ----

loggerName <- gsub(":| |-", "_", paste0("log ", Sys.time(),".txt"))
logger <- create.logger()
logfile(logger) <- here(resultsFolder, loggerName)
level(logger) <- "INFO"
info(logger, "LOG CREATED")

# Start time
tic.clearlog()
tic.clear()
tic(msg = "HPVStudy Total run time: ")

# Parameters for analyses
instantiateVaccinationCohorts <- FALSE
instantiatePopulationCohorts <- FALSE
instantiateCharacteristics <- FALSE
instantiateOutcomeCohorts <- FALSE

# Study parameters
coverages <- 0
sex_cohorts <- c("Female", "Male")

# Code parameters
dbName <- "CPRD GOLD"

tic(msg = "Creation of CDM object: ")
# create cdm object ----
info(logger, "CREATING CDM OBJECT")
cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = server_dbi, 
  achillesSchema = achillesSchema,
  cohortTables = c(# Vaccination cohorts 
                   "firstdose_cohort",
                   "doses_allvac_cohort", 
                   # Female cohorts overall
                   "unvac_18f_coverage_0_any_dose",
                   "vac_18f_coverage_0_any_dose",
                   "matched_unvac_cov_0_ratio_18f_any_dose",
                   "matched_vac_cov_0_ratio_18f_any_dose",
                   # Female cohorts doses
                   "not_matched_vaccinated_18f_coverage_0_1dose",
                   "not_matched_vaccinated_18f_coverage_0_23dose",
                   "not_matched_vaccinated_18f_coverage_0_2dose",
                   "not_matched_vaccinated_18f_coverage_0_3dose",
                   "matched_vaccinated_18f_coverage_0_1dose",
                   "matched_vaccinated_18f_coverage_0_23dose",
                   "matched_vaccinated_18f_coverage_0_2dose",
                   "matched_vaccinated_18f_coverage_0_3dose",
                   # Male cohorts overall
                   "unvac_15m_coverage_0_any_dose",
                   "vac_15m_coverage_0_any_dose",
                   "matched_unvac_cov_0_ratio_15m_any_dose",
                   "matched_vac_cov_0_ratio_15m_any_dose",
                   # Male cohorts doses
                   "not_matched_vaccinated_15m_coverage_0_1dose",
                   "not_matched_vaccinated_15m_coverage_0_23dose",
                   "not_matched_vaccinated_15m_coverage_0_2dose",
                   "not_matched_vaccinated_15m_coverage_0_3dose",
                   "matched_vaccinated_15m_coverage_0_1dose",
                   "matched_vaccinated_15m_coverage_0_23dose",
                   "matched_vaccinated_15m_coverage_0_2dose",
                   "matched_vaccinated_15m_coverage_0_3dose",
                   # Covariates
                   "medications",
                   "conditions",
                   "hiv_status",
                   "vaccinations",
                   "cervical_screening",
                   "smear_test",
                   # Outcomes of interest
                   "darwin_outcome_cohorts",
                   "nondarwin_outcome_cohorts",
                   # NCO outcomes
                   "nondarwin_nco_cohorts"
                   )
  
)
info(logger, "CDM OBJECT CREATED")
toc(log = TRUE)

# create and export snapshot
# cdm_snapshot <- snapshot(cdm)
# write.csv(cdm_snapshot, here("Results", paste0(
#   "cdm_snapshot_", cdmName(cdm), "_" ,format(Sys.time(), "_%Y_%m_%d"), ".csv"
# )))

for (sex_cohort in sex_cohorts) {
  info(logger, paste0("Start study for ", sex_cohort, "cohort"))
  if (sex_cohort == "Male") {
    yy <- 15
    ss <- "m"
    start_birth <- 2004
  } else if (sex_cohort == "Female") {
    yy <- 18
    ss <- "f"
    start_birth <- 1990
  }

# instantiate necessary cohorts ----
tic(msg = "Instantiation of Cohorts: ")
info(logger, "INSTANTIATING STUDY COHORTS")
source(here("Cohorts", "InstantiateCohorts.R"))
info(logger, "STUDY COHORTS INSTANTIATED")
toc(log = TRUE)

# instantiate characteristics ----
tic(msg = "Instantiation of Characteristics: ")
info(logger, "INSTANTIATING STUDY CHARACTERISTICS")
source(here("Cohorts", "InstantiateCharacteristics.R"))
info(logger, "CHARACTERISTICS INSTANTIATED")
toc(log = TRUE)

# Overall matching ----
tic(msg = "Overall Matching: ")
info(logger, "RUN MATCHING")
source(here("Analyses", "1_0-OverallMatching_ratio.R"))
info(logger, "OVERALL MATCHING FINISHED")
toc(log = TRUE)

tic(msg = "Dose stratification Matching: ")
info(logger, "RUN DOSE MATCHING:")
source(here("Analyses", "1_3-AllDoseMatching_ratio.R"))
info(logger, "DOSE MATCHING FINISHED")
toc(log = TRUE)

tic(msg = "Characterisation: ")
info(logger, "RUN CHARACTERISATION")
source(here("Analyses", "2-Characterisation.R"))
info(logger, "CHARACTERISATION FINISHED")
toc(log = TRUE)

tic(msg = "Outcome model: ")
info(logger, "RUN OUTCOME MODEL")
source(here("Analyses", "3-OutcomeModel.R"))
info(logger, "OUTCOME MODEL FINISHED")
toc(log = TRUE)

tic(msg = "Outcome model: ")
info(logger, "RUN NCO")
source(here("Analyses", "3_2-NCOutcomeModel.R"))
info(logger, "NCO MODEL FINISHED")
toc(log = TRUE)
}

toc(log = TRUE)
tic.log(format = TRUE)
tic_log <- tic.log(format = TRUE)

# output$log <- tibble(cdm_name = input$cdmName, log = paste0(tic_log %>%  unlist(), collapse = "\n"))
# write_csv(output$log, here("Results", paste0(
#   "log_", cdmName(cdm), "_" ,format(Sys.time(), "_%Y_%m_%d"), ".csv"
# )))

# export results ----
info(logger, "EXPORTING RESULTS")
zip(
  zipfile = file.path(paste0(resultsFolder, "/Results_", cdmName(cdm), ".zip")),
  files = list.files(resultsFolder, full.names = TRUE)
)