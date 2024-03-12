library(tidyr)
library(CDMConnector)
library(RPostgres)
library(PatientProfiles)

# Details
log_file        <- here(resultsFolder, "log_instantaniationData.txt")
logger          <- create.logger(logfile = log_file, level = "INFO")

# instantiate medications
info(logger, "Instantiate medications")
codelistMedications <- codesFromConceptSet(here("MedicationsCharacteristics"), cdm)
cdm <- generateConceptCohortSet(
  cdm = cdm,
  name = "medications",
  conceptSet = codelistMedications
)

# instantiate conditions
info(logger, "Instantiate conditions")
codelistConditions <- codesFromCohort(path = here("ConditionsCharacteristics"), cdm = cdm)
autoimmune <- c("Type 1 Diabetes Mellitus", "Rheumatoid_arthritis", "Psoriasis", "Psoriatic Arthritis", "Multiple sclerosis", "Addison's disease", "Graves' disease", "Sjogren's syndrome", "Hashimoto thyroiditis", "Myasthenia gravis", "Vasculitis", "Pernicious anemia", "Celiac disease", "Scleroderma", "Sarcoidosis", "Ulcerative colitis", "Crohn's disease")
codelistConditions$Autoimmune_Disease <- codelistConditions[autoimmune] |> unlist() |> unname() |> unique()
codelistConditions[autoimmune] <- NULL
names(codelistConditions) <- gsub(" ", "_", tolower(names(codelistConditions)))
cdm <- generateConceptCohortSet(
  cdm = cdm, name = "conditions", conceptSet = codelistConditions,
  end = "observation_period_end_date", overwrite = TRUE
)

# instantiate HPV status
info(logger, "Instantiate HIV status")
HIVstatus <- readCohortSet(here("HIVStatus"))
cdm <- generateCohortSet(
  cdm = cdm,
  name = "hiv_status",
  cohortSet = HIVstatus,
)

# instantiate panicolau smear testing
info(logger, "Instantiate panicolau smear testing")
panicolauSmearTesting <- readCohortSet(here("PapanicolauSmearTest"))
cdm <- generateCohortSet(
  cdm = cdm,
  name = "papanicolau_smear_testing",
  cohortSet = panicolauSmearTesting,
)

# instantiate previous vaccinations
info(logger, "Instantiate previous vaccinations")
previousVaccinations <- readCohortSet(here("Vaccinations"))
cdm <- generateCohortSet(
  cdm = cdm,
  name = "vaccinations",
  cohortSet = previousVaccinations,
)

# instantiate cyinfo(logger, "Instantiate previous vaccinations")
cytology_results <- readCohortSet(here("Cytology_results"))
cdm <- generateCohortSet(
  cdm = cdm,
  name = "cytology",
  cohortSet = cytology_results,
)
