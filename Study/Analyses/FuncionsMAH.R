recordCohortAttrition2 <- function(cohort, reason, cohortId = NULL) {
  # get cohortId
  cohortId <- getCohortId(cohort, cohortId)
  # updateAttrition
  newAttrition <- updateAttrition(cohort, cohortId, reason) |>
    dplyr::as_tibble()
  # create cohort table
  cohort <- CDMConnector::newGeneratedCohortSet(
    cohortRef = cohort,
    cohortAttritionRef = newAttrition,
    cohortSetRef = CDMConnector::cohortSet(cohort),
    overwrite = TRUE
  )
  return(cohort)
}

getCohortId <- function(cohort, cohortId) {

  possibleCohortId <- CDMConnector::cohortSet(cohort) |> dplyr::pull("cohort_definition_id")

  if (is.null(cohortId)) {

    cohortId <- possibleCohortId

  } else if (!all(cohortId %in% possibleCohortId)) {

    cli::cli_abort("cohort_definition_id must be defined in the cohort_set.")

  }

  return(cohortId)

}

updateAttrition <- function(cohort, cohortId, reason) {

  oldAttrition <- CDMConnector::cohortAttrition(cohort)

  newRow <- CDMConnector::cohortCount(cohort) |>

    dplyr::select(

      "cohort_definition_id", "number_records", "number_subjects"

    ) |>

    dplyr::filter(.data$cohort_definition_id %in% .env$cohortId) |>

    dplyr::rename(

      "previous_records" = "number_records",

      "previous_subjects" = "number_subjects"

    ) |>

    dplyr::left_join(

      cohort |>

        dplyr::filter(.data$cohort_definition_id %in% .env$cohortId) |>

        dplyr::group_by(.data$cohort_definition_id) |>

        dplyr::summarise(

          "number_records" = dplyr::n(),

          "number_subjects" = dplyr::n_distinct(.data$subject_id),

          .groups = "drop"

        ) |>

        dplyr::collect(),

      by = "cohort_definition_id"

    ) |>

    dplyr::mutate(dplyr::across(

      dplyr::all_of(c("number_records", "number_subjects")),

      ~ dplyr::if_else(is.na(.x), as.integer(0), as.integer(.x))

    )) |>

    dplyr::mutate(

      "excluded_records" = .data$previous_records - .data$number_records,

      "excluded_subjects" = .data$previous_subjects - .data$number_subjects

    ) |>

    dplyr::inner_join(

      oldAttrition |>

        dplyr::select("cohort_definition_id", "reason_id") |>

        dplyr::filter(.data$cohort_definition_id %in% .env$cohortId) |>

        dplyr::group_by(.data$cohort_definition_id) |>

        dplyr::summarise(

          "reason_id" = max(.data$reason_id), .groups = "drop"

        ) |>

        dplyr::mutate(

          "reason_id" = .data$reason_id + 1, "reason" = .env$reason

        ),

      by = "cohort_definition_id"

    ) |>

    dplyr::select(

      "cohort_definition_id", "number_records", "number_subjects", "reason_id",

      "reason", "excluded_records", "excluded_subjects"

    )

  newAttrition <- oldAttrition |>

    dplyr::bind_rows(newRow)

  return(newAttrition)

}

addRegion <- function(x) {
  
  cdm <- attr(x, "cdm_reference")
  
  if (region_type == "char") {
    x |>
      dplyr::left_join(
        cdm$person |>
          dplyr::select("person_id", "care_site_id") |>
          dplyr::inner_join(
            cdm$care_site |> dplyr::select("care_site_id", "location_id"),
            by = "care_site_id"
          ) |>
          dplyr::inner_join(
            cdm$location |>
              dplyr::select("location_id", "region" = "location_source_value"),
            by = "location_id"
          ) |>
          dplyr::select("subject_id" = "person_id", "region"),
        by = "subject_id"
      ) |>
      CDMConnector::computeQuery()
  }else{
    x |>
      mutate(region = 1)
  }
}

addCytologyResults <- function(x, name){
  cdm <- attr(x, "cdm_reference")
  x |>
    PatientProfiles::addCohortIntersectFlag(
      targetCohortTable = "hpv_covariates",
      targetCohortId    = CDMConnector::settings(cdm[["hpv_covariates"]]) |> 
        dplyr::filter(cohort_name == name) |> 
        dplyr::pull(cohort_definition_id),
      indexDate  = "index_date",
      censorDate = NULL,
      targetStartDate = "cohort_start_date",
      targetEndDate   = "cohort_end_date",
      window = list(c(-Inf,-1)),
      nameStyle = "cytology_results"
    ) |>
    select(-starts_with("id_"))
}

addVisitsPriorYear <- function(x){
  x |>
    PatientProfiles::addTableIntersectCount(
      tableName  = "visit_occurrence",
      indexDate  = "index_date",
      censorDate = NULL,
      window     = list(c(-365,-1)),
      overlap    = TRUE,
      nameStyle  = "prior_visits"
    )
}

addHIV <- function(x, name){
  cdm <- attr(x, "cdm_reference")
  x |>
    PatientProfiles::addCohortIntersectFlag(
      targetCohortTable = "hpv_covariates",
      targetCohortId    = CDMConnector::settings(cdm[["hpv_covariates"]]) |> 
        dplyr::filter(cohort_name == name) |> 
        dplyr::pull(cohort_definition_id),
      indexDate  = "index_date",
      censorDate = NULL,
      targetStartDate = "cohort_start_date",
      targetEndDate   = "cohort_end_date",
      window = list(c(-Inf,-1)),
      nameStyle = "hiv_history"
    )
}

addPreviousPapanicolau <- function(x, name){
  cdm <- attr(x, "cdm_reference")
  x |>
    PatientProfiles::addCohortIntersectFlag(
      targetCohortTable = "hpv_covariates",
      targetCohortId    = CDMConnector::settings(cdm[["hpv_covariates"]]) |> 
        dplyr::filter(cohort_name == name) |> 
        dplyr::pull(cohort_definition_id),
      indexDate  = "index_date",
      censorDate = NULL,
      targetStartDate = "cohort_start_date",
      targetEndDate   = "cohort_end_date",
      window = list(c(-Inf,-1)),
      nameStyle = "papanicolau_smear_testing"
    )
}

addPreviousVaccinations <- function(x){
  cdm <- attr(x, "cdm_reference")
  codes <- CodelistGenerator::getATCCodes(cdm,
                                          level = "ATC 2nd",
                                          name  = "VACCINES")
  cdm <- CDMConnector::generateConceptCohortSet(
    cdm  = cdm,
    conceptSet = codes,
    name = "previous_vaccines",
    limit = "all",
    end = 0,
    overwrite = TRUE
  )
  attr(x,"cdm_reference") <- cdm
  x |> PatientProfiles::addCohortIntersectCount(
    targetCohortTable = "previous_vaccines",
    targetCohortId    = NULL,
    indexDate         = "index_date",
    censorDate = NULL,
    targetStartDate = "cohort_start_date",
    targetEndDate   = "cohort_end_date",
    window    = list(c(-Inf,-1)),
    nameStyle = "previous_vaccines")
}

#covariates <- readCohortSet(here("Cohorts","MAH_covariates"))
#cdm <- generateCohortSet(
#  cdm = cdm,
#  cohortSet = covariates,
#  overwrite = TRUE,
#  name = "hpv_covariates"
#)
