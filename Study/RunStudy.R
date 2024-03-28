# create logger ----
resultsFolder <- here("Results")
loggerName <- gsub(":| |-", "_", paste0("log ", Sys.time(),".txt"))
logger <- create.logger()
logfile(logger) <- here(resultsFolder, loggerName)
level(logger) <- "INFO"
info(logger, "LOG CREATED")

# Start time
tic.clearlog()
tic.clear()
tic(msg = "HPVStudy Total run time: ")

# parameters ----
instantiateCohorts <- TRUE
instantiateCharacteristics <- FALSE

tic(msg = "Creation of CDM object: ")
# create cdm object ----
info(logger, "CREATING CDM OBJECT")
cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = dbName, 
  achillesSchema = achillesSchema,
  cohortTables = "allvac_cohort"

)
info(logger, "CDM OBJECT CREATED")
toc(log = TRUE)

# create and export snapshot
# cdm_snapshot <- snapshot(cdm)
# write.csv(cdm_snapshot, here("Results", paste0(
#   "cdm_snapshot_", cdmName(cdm), "_" ,format(Sys.time(), "_%Y_%m_%d"), ".csv"
# )))

tic(msg = "Instantiation of Cohorts: ")
# instantiate necessary cohorts ----
#info(logger, "INSTANTIATING STUDY COHORTS")
source(here("Cohorts", "InstantiateCohorts.R"))
info(logger, "STUDY COHORTS INSTANTIATED")
toc(log = TRUE)

# tic(msg = "Instantiation of Characteristics: ")
# info(logger, "INSTANTIATING STUDY CHARACTERISTICS")
# source(here("Cohorts", "InstantiateCharacteristics.R"))
# info(logger, "CHARACTERISTICS INSTANTIATED")
# toc(log = TRUE)

# tic(msg = "PhenotypeR: ")
# # run diagnostics ----
# info(logger, "RUN PHENOTYPER")
# source(here("PhenotypeR", "PhenotypeR.R"))
# info(logger, "PHENOTYPER FINISHED")
# toc(log = TRUE)

# tic(msg = "Matching: ")
# # run analyses ----
# info(logger, "RUN MATCHING")
# source(here("Analyses", "1-PSMatchingImplementation.R"))
# info(logger, "MATCHING FINISHED")
# toc(log = TRUE)
# 
# tic(msg = "Characterisation: ")
# info(logger, "RUN CHARACTERISATION")
# source(here("Analyses", "2-Characterisation.R"))
# info(logger, "CHARACTERISATION FINISHED")
# toc(log = TRUE)

#tic(msg = "Outcome model: ")
#info(logger, "RUN OUTCOME MODEL")
#source(here("Analyses", "3-OutcomeModel.R"))
#info(logger, "OUTCOME MODEL FINISHED")
#toc(log = TRUE)

toc(log = TRUE)
tic.log(format = TRUE)
tic_log <- tic.log(format = TRUE)

output$log <- tibble(cdm_name = input$cdmName, log = paste0(tic_log %>%  unlist(), collapse = "\n"))
write_csv(output$log, here("Results", paste0(
  "log_", cdmName(cdm), "_" ,format(Sys.time(), "_%Y_%m_%d"), ".csv"
)))

# export results ----
# info(logger, "EXPORTING RESULTS")
# zip(
#   zipfile = file.path(paste0(resultsFolder, "/Results_", cdmName(cdm), ".zip")),
#   files = list.files(resultsFolder, full.names = TRUE)
# )