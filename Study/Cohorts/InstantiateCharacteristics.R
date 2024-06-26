library(tidyr)
library(CDMConnector)
library(RPostgres)
library(PatientProfiles)

if (instantiateCharacteristics) {
  info(logger, "INSTANTANIATE CHARACTERISTICS")
  # instantiate medications
  info(logger, "Instantiate medications")
  codelistMedications <- codesFromConceptSet(here("Cohorts", "MedicationsCharacteristics"), cdm)
  names(codelistMedications) <- gsub(" ", "_", tolower(names(codelistMedications)))
  # https://darwin-eu-dev.github.io/DrugUtilisation/reference/generateDrugUtilisationCohortSet.html
  cdm <- generateConceptCohortSet(
    cdm = cdm,
    name = "medications",
    conceptSet = codelistMedications,
    end = 0
  )
  
  # instantiate conditions
  info(logger, "Instantiate conditions")
  codelistConditions <- codesFromCohort(path = here("Cohorts", "ConditionsCharacteristics"), cdm = cdm)
  autoimmune <- c("Type 1 Diabetes Mellitus", "Rheumatoid_arthritis", "Psoriasis", "Psoriatic Arthritis", "Multiple sclerosis", "Addison's disease", "Graves' disease", "Sjogren's syndrome", "Hashimoto thyroiditis", "Myasthenia gravis", "Vasculitis", "Pernicious anemia", "Celiac disease", "Scleroderma", "Sarcoidosis", "Ulcerative colitis", "Crohn's disease")
  codelistConditions$Autoimmune_Disease <- codelistConditions[autoimmune] |> unlist() |> unname() |> unique()
  codelistConditions[autoimmune] <- NULL
  names(codelistConditions) <- gsub(" ", "_", tolower(names(codelistConditions)))
  cdm <- generateConceptCohortSet(
    cdm = cdm, name = "conditions", conceptSet = codelistConditions,
    end = 0, overwrite = TRUE
  )
  
  # instantiate HIV status
  info(logger, "Instantiate HIV status")
  HIVstatus <- readCohortSet(here("Cohorts","HIVStatus"))
  cdm <- generateCohortSet(
    cdm = cdm,
    name = "hiv_status",
    cohortSet = HIVstatus,
  )
  
  # instantiate papanicolau smear testing
  info(logger, "Instantiate panicolau smear testing")
  panicolauSmearTesting <- readCohortSet(here("Cohorts","PapanicolauSmearTest"))
  cdm <- generateCohortSet(
    cdm = cdm,
    name = "smear_test",
    cohortSet = panicolauSmearTesting,
  )
  
  # instantiate previous vaccinations
  info(logger, "Instantiate previous vaccinations")
  previousVaccinations <- CodelistGenerator::getATCCodes(cdm,
                                                         level = "ATC 2nd",
                                                         name  = "VACCINES")
  cdm <- CDMConnector::generateConceptCohortSet(
    cdm  = cdm,
    conceptSet = previousVaccinations,
    name = "vaccinations",
    limit = "all",
    end = 0,
    overwrite = TRUE
  )
  
  
  # instantiate cytology
  info(logger, "Instantiate cytology")
  cytology_results <- readCohortSet(here("Cohorts","Cytology_results"))
  cdm <- generateCohortSet(
    cdm = cdm,
    name = "cervical_screening",
    cohortSet = cytology_results,
  )
  
} else {
  cdm <- cdmFromCon(
    con = db,
    cdmSchema = cdmSchema, 
    writeSchema = writeSchema, 
    cdmName = server_dbi, 
    achillesSchema = achillesSchema,
    cohortTables = c("firstdose_cohort",
                     "doses_allvac_cohort",
                     paste0("unvac_",yy,ss,"_coverage_0_any_dose"),
                     paste0("vac_",yy,ss,"_coverage_0_any_dose"),
                     "medications",
                     "conditions",
                     "hiv_status",
                     "vaccinations",
                     "cervical_screening",
                     "smear_test"
                     )
    )
}