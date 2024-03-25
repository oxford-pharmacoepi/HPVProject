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
instantiateCharacteristics <- TRUE

tic(msg = "Creation of CDM object: ")
# create cdm object ----
info(logger, "CREATING CDM OBJECT")
cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = dbName, 
  achillesSchema = achillesSchema,
<<<<<<< HEAD
  #cohortTables = c("vac_cohort", "unvac_cohort", "allvac_cohort", "vaccinations", "hiv_status", "conditions", "medications", "papanicolau_smear_testing", "cytology")
=======
<<<<<<< HEAD
  cohortTables = "allvac_cohort"
=======
  cohortTables = c("vac_cohort", "unvac_cohort", "allvac_cohort", "vaccinations", "hiv_status", "conditions", "medications", "papanicolau_smear_testing", "cytology")
>>>>>>> 5966c558c99a19585408c6c50ced1e975cbcf89b
>>>>>>> 3613a6e677e0a5296e76af66f8b57a30f68ffedb
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

tic(msg = "Instantiation of Characteristics: ")
info(logger, "INSTANTIATING STUDY CHARACTERISTICS")
source(here("Cohorts", "InstantiateCharacteristics.R"))
info(logger, "CHARACTERISTICS INSTANTIATED")
toc(log = TRUE)

# tic(msg = "PhenotypeR: ")
# # run diagnostics ----
# info(logger, "RUN PHENOTYPER")
# source(here("PhenotypeR", "PhenotypeR.R"))
# info(logger, "PHENOTYPER FINISHED")
# toc(log = TRUE)

tic(msg = "Matching: ")
# run analyses ----
info(logger, "RUN MATCHING")
source(here("Analyses", "1-PSMatchingImplementation.R"))
info(logger, "MATCHING FINISHED")
toc(log = TRUE)

tic(msg = "Characterisation: ")
info(logger, "RUN CHARACTERISATION")
source(here("Analyses", "2-Characterisation.R"))
info(logger, "CHARACTERISATION FINISHED")
toc(log = TRUE)

#tic(msg = "Outcome model: ")
#info(logger, "RUN OUTCOME MODEL")
#source(here("Analyses", "3-OutcomeModel.R"))
#info(logger, "OUTCOME MODEL FINISHED")
#toc(log = TRUE)


# export results ----
# info(logger, "EXPORTING RESULTS")
# zip(
#   zipfile = file.path(paste0(resultsFolder, "/Results_", cdmName(cdm), ".zip")),
#   files = list.files(resultsFolder, full.names = TRUE)
# )