# create logger ----
resultsFolder <- here("Results")
loggerName <- gsub(":| |-", "_", paste0("log ", Sys.time(),".txt"))
logger <- create.logger()
logfile(logger) <- here(resultsFolder, loggerName)
level(logger) <- "INFO"
info(logger, "LOG CREATED")

# parameters ----
instantiateCohorts <- TRUE

# create cdm object ----
info(logger, "CREATING CDM OBJECT")
cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = dbName, 
  achillesSchema = achillesSchema
)
info(logger, "CDM OBJECT CREATED")

# create and export snapshot

# instantiate necessary cohorts ----
#info(logger, "INSTANTIATING STUDY COHORTS")
source(here("Cohorts", "InstantiateCohorts.R"))
info(logger, "STUDY COHORTS INSTANTIATED")

info(logger, "INSTANTIATING STUDY CHARACTERISTICS")
source(here("Cohorts", "InstantiateCharacteristics.R"))
info(logger, "CHARACTERISTICS INSTANTIATED")

# run diagnostics ----
#info(logger, "RUN PHENOTYPER")
#source(here("PhenotypeR", "PhenotypeR.R"))
#info(logger, "PHENOTYPER FINISHED")

# run analyses ----
info(logger, "RUN MATCHING")
source(here("Analyses", "1-PSMatchingImplementation.R"))
info(logger, "MATCHING FINISHED")

info(logger, "RUN CHARACTERISATION")
source(here("Analyses", "2-Characterisation.R"))
info(logger, "CHARACTERISATION FINISHED")

# export results ----
info(logger, "EXPORTING RESULTS")
zip(
  zipfile = file.path(paste0(resultsFolder, "/Results_", cdmName(cdm), ".zip")),
  files = list.files(resultsFolder, full.names = TRUE)
)