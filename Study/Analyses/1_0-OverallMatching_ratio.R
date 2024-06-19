library(tidyr)
library(CDMConnector)
library(MatchIt)
library(Matching)
library(RPostgres)
library(glmnet)
library(purrr)
library(log4r)
library(PatientProfiles)
library(here)
library(dplyr)
library(rlang)
source(here("Analyses","FuncionsMAH.R"))

info(logger = logger, "CREATE INITIAL POPULATION")

for (cov in coverages) {
  # Total Population
  vac_cohortname <- paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")
  unvac_cohortname <- paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")
  
  population <- union(cdm[[vac_cohortname]] |> mutate(vac_status = 1), cdm[[unvac_cohortname]] |> mutate(vac_status = 0)) |>
    addDateOfBirth() |>
    mutate(date_9years = date_of_birth + years(9)) |>
    mutate(!!paste0("date_",yy,"years") := date_of_birth + years(yy))
  
  size_pop <- population |> tally() |> pull()
  size_vacpop <- population |> filter(vac_status == 1) |> tally() |> pull()
  size_unvacpop <- population |> filter(vac_status == 0) |> tally() |> pull()
  info(logger = logger, paste0("Size Population = ", size_pop))
  info(logger = logger, paste0("Size Vacc Population = ", size_vacpop))
  info(logger = logger, paste0("Size Unvacc Population = ", size_unvacpop))
  
  name_control <- paste0("matched_unvac_cov_", cov,"_ratio_",yy,ss,"_any_dose")
  name_treatment <- paste0("matched_vac_cov_", cov,"_ratio_",yy,ss, "_any_dose")
  name_all <- paste0("overall_matched_coverage_", cov,"_ratio_",yy,ss,"_any_dose")
  
  # Total conditions and drugs
  Conditions <- cdm$condition_occurrence |>
    inner_join(population |> dplyr::select("person_id" = "subject_id", paste0("date_",yy,"years")), by = "person_id") |>
    compute() |>
    rename(subject_id = person_id, concept_id = condition_concept_id, occurrence_start_date = condition_start_date) |>
    left_join(cdm$concept |> dplyr::select(concept_id, concept_name), by = "concept_id") |>
    dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, paste0("date_",yy,"years")) |>
    mutate(occurrence_type = "condition") |>
    compute()
  
  Drugs <-  cdm$drug_exposure |>
    inner_join(population |> dplyr::select("person_id" = "subject_id", paste0("date_",yy,"years")), by = "person_id") |>
    compute() |>
    rename(subject_id = person_id, concept_id = drug_concept_id, occurrence_start_date = drug_exposure_start_date)|>
    left_join(cdm$concept |> dplyr::select(concept_id, concept_name), by = "concept_id") |>
    dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, paste0("date_",yy,"years")) |>
    mutate(occurrence_type = "drug") |>
    compute()
  
  ConditionsAndDrugs <- union_all(Conditions, Drugs) |>
    mutate(occurrence_flag = 1) |>
    filter(occurrence_start_date < !!sym(paste0("date_",yy,"years"))) |>
    dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, occurrence_type, occurrence_flag) |>
    compute()
  
  # Prepare covariates
  hiv_status <- population |>
    addCohortIntersectFlag(targetCohortTable = "hiv_status",
                           indexDate = paste0("date_",yy,"years"),
                           window = c(-Inf, -1),
                           nameStyle = "hiv_status_flag"
    ) |>
    dplyr::select(subject_id, hiv_status_flag) |>
    compute()
  
  papanicolau_smear_test <- population |>
    addCohortIntersectFlag(targetCohortTable = "smear_test",
                           indexDate = paste0("date_",yy,"years"),
                           window = c(-Inf, -1),
                           nameStyle = "papanicolau_smear_testing_flag"
    ) |>
    dplyr::select(subject_id, papanicolau_smear_testing_flag) |>
    compute()
  
  vaccinations <- population |>
    addCohortIntersectCount(targetCohortTable = "vaccinations",
                            indexDate = paste0("date_",yy,"years"),
                            window = c(-Inf, -1),
                            nameStyle = "vaccinations_flag"
    ) |>
    dplyr::select(subject_id, vaccinations_flag) |>
    compute()
  
  cytology <- population |>
    addCohortIntersectFlag(targetCohortTable = "cervical_screening",
                           indexDate = paste0("date_",yy,"years"),
                           window = c(-Inf, -1),
                           nameStyle = "hiv_status_flag"
    ) |>
    dplyr::select(subject_id, hiv_status_flag) |>
    compute()
  
  # Results file and cohorts
  lasso_selectedfeatures <- tibble("concept_id" = as.numeric(),
                                   "coefficients" = as.numeric(),
                                   "length_id" = as.character(),
                                   "concept_name" = as.character(),
                                   "year"     = as.numeric(),
                                   "database" = as.character()
  )
  
  total_matched_cohort <- tibble("cohort_definition_id" = as.numeric(),
                                 "subject_id" = as.numeric(),
                                 "cohort_year" = as.numeric(),
                                 "vac_status" = as.numeric(),
                                 "cohort_start_date" =  as.Date(character()),
                                 "cohort_end_date" =  as.Date(character()), 
                                 "pair_id" = as.numeric()
  )
  
  
  cdm <- cdm |>
    CDMConnector::insertTable(
      name = name_all,
      table = total_matched_cohort,
      overwrite = TRUE
    ) 
  
  # Subpopulation
  # initialisation of parameters
  last_pair <- 0
  
  for(sy in c(2008:2023)){#(year(studyEndDate)-year(studyStartDate))){
    # Execution time
    start_time <- Sys.time()
    
    print(sy)
    
    info(logger = logger, paste0("CREATE SUBPOPULATION ", sy))
    
    subpopulation <-  union(population |> 
                              filter(vac_status == 0,
                                     cohort_start_date <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                     cohort_end_date >= as.POSIXct(paste0(as.integer(sy),"-12-31")),
                                     date_9years <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                     !!sym(paste0("date_",yy,"years")) >= as.POSIXct(paste0(as.integer(sy),"-12-31"))), 
                            population |> 
                              filter(vac_status == 1, 
                                     cohort_start_date <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                     cohort_end_date >= as.POSIXct(paste0(as.integer(sy),"-12-31")),
                                     date_9years <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                     !!sym(paste0("date_",yy,"years")) >= as.POSIXct(paste0(as.integer(sy),"-12-31"))) |> 
                              inner_join(cdm$firstdose_cohort |> 
                                           filter(cohort_start_date >= as.POSIXct(paste0(as.integer(sy),"-01-01")), 
                                                  cohort_start_date <= as.POSIXct(paste0(as.integer(sy),"-12-31"))) |> 
                                           dplyr::select(subject_id)
                              )
    ) |> 
      compute()
    
    prev <- subpopulation |> 
      filter(vac_status == 0) |> 
      inner_join(cdm$firstdose_cohort |>
                   filter(cohort_start_date < as.POSIXct(paste0(as.integer(sy),"-01-01"))), by = "subject_id") |>
      collect()
    # THE PROBLEM IS THE NOT VACCINATED PEOPLE. I have people which are vaccinated after 15 years but before the index_date selected.
    # This could be ERRROR but this people will never be matched. SHould they be removed before lasso?
    
    if (nrow(prev) == 0) {
      check <- "OK!"
    }else{
      check <- "ERROR"
    }
    info(logger, "Check vac subpopulation is not already vaccinated: ", check)
    subpopulation |> tally()
    
    # Set new index date: 01/01/2008
    subpopulation <- subpopulation |> 
      mutate(index_date = as.POSIXct(paste0(as.integer(sy),"-01-01"))) |>
      compute()
    subpopulation |> tally()
    
    # Checks:
    # 365 prior information from index_date
    subpopulation <- subpopulation |> 
      addPriorObservation(indexDate = "index_date", priorObservationName = "prior_observation_indexdate") |> 
      filter(prior_observation_indexdate >= 365) |>
      compute()
    subpopulation |> tally()
    
    # Add year of birth
    subpopulation <- subpopulation |> 
      mutate(year_of_birth = year(date_of_birth)) |>
      compute()
    
    size_subpopulation <- as.numeric(count(distinct(subpopulation, subject_id)) |> compute() |> pull())
    size_vacsubpopulation <- as.numeric(count(distinct(subpopulation |> filter(vac_status == 1), subject_id)) |> compute() |> pull())
    size_unvacsubpopulation <- as.numeric(count(distinct(subpopulation |> filter(vac_status == 0), subject_id)) |> compute() |> pull())
    info(logger = logger, paste0("Size subpopulation = ", size_subpopulation))
    info(logger = logger, paste0("Size vaccinated subpopulation = ", size_vacsubpopulation))
    info(logger = logger, paste0("Size unvaccinated subpopulation = ", size_unvacsubpopulation))
    
    if (size_vacsubpopulation < 5 | size_unvacsubpopulation < 5) {
      info(logger = logger, "Not enough participants to proceed with the lasso / matching")
      population <- population |> 
        anti_join(subpopulation |> 
                    filter(vac_status == 1), 
                  by = "subject_id"
        )
      remaining_vac <- population |> filter(vac_status == 1) |> tally() |> pull()
      info(logger = logger, paste0("Remaining vaccinated subjects in population = ", remaining_vac))
    }else{
      
      # Subpopulation as Lasso Input
      sub_ConditionsAndDrugs <- subpopulation |> 
        dplyr::select(subject_id, index_date, year_of_birth, vac_status) |> 
        inner_join(
          # Occurrence Record before index date
          ConditionsAndDrugs, by = "subject_id"
        )
      size <- sub_ConditionsAndDrugs |> distinct(concept_id) |> tally() |> pull()
      info(logger, paste0("Conditions only for ",sy," subpopulation: ", size))
      
      sub_ConditionsAndDrugs <- sub_ConditionsAndDrugs |>
        filter(occurrence_start_date < index_date) %>%
        mutate(date_diff = !!datediff("index_date","occurrence_start_date", interval = "day")) |> 
        mutate(window = case_when(
          date_diff >= -31  & date_diff < 0        ~ "short",
          date_diff >= -365 & date_diff < -31      ~ "mid",
          date_diff < -365 & occurrence_type == "condition" ~ "long",
          TRUE ~ NA
        )) |>  
        filter(!is.na(window)) |>
        mutate(
          concept = paste0(window, "_", round(as.numeric(concept_id),0))) |>
        dplyr::select(subject_id, concept, concept_name, occurrence_flag) |>
        distinct() |>
        compute()
      
      size <- sub_ConditionsAndDrugs |> distinct(concept) |> tally() |> pull()
      info(logger, "Conditions after windows: ", size)
      
      lowfreq_oc <- sub_ConditionsAndDrugs |> 
        # Only those present with a frequency of > 0.5%
        group_by(concept) |> 
        tally() |>
        filter(n < 0.005*size_subpopulation) |> 
        dplyr::select(concept) |>
        compute()
      
      # Eliminate occurrences with low frequencies and generate a column for each occurrence
      sub_ConditionsAndDrugs <- sub_ConditionsAndDrugs |>  
        anti_join(lowfreq_oc, by = "concept") |>
        pivot_wider(id_cols = subject_id, names_from = concept, values_from = occurrence_flag, values_fill = 0) |> 
        compute()
      sub_ConditionsAndDrugs |> tally()
      
      # Add Drugs and Conditions, and covariates
      subpop_data <- subpopulation |> 
        dplyr::select(subject_id, index_date, year_of_birth, vac_status) |> 
        addVisitsPriorYear() |>
        addPreviousVaccinations() |>
        left_join(papanicolau_smear_test, by = "subject_id") |>
        left_join(cytology, by = "subject_id") |>
        left_join(hiv_status, by = "subject_id") |>
        left_join(sub_ConditionsAndDrugs, by = "subject_id") |> 
        compute()
      subpop_data |> tally()
      
      # Compute Lasso and extract Features of interest
      in_data <- subpop_data |> compute()
      
      # All conditions and drugs present in the subjects previous to the index_date
      x <- in_data |> 
        dplyr::select(! c(subject_id, index_date, vac_status, year_of_birth)) |> 
        as_tibble() |>
        mutate_all(~replace(., is.na(.), 0)) |> 
        data.matrix() 
      info(logger = logger, paste0("Number of selected features = ", ncol(x)))
      
      # Vac_status = 0/1
      y <- in_data |> 
        dplyr::select(vac_status) |> 
        compute() |> 
        as_tibble() |> 
        data.matrix()
      
      # Model definition
      lambdas <- 10^seq(2, -3, by = -.1) # magnitude order lambda
      best_model <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standarize = TRUE, nforlds = 5, family = "binomial")
      # how do we know this is the best lambda???????????????????????
      
      coef.lasso_reg <- coef(best_model, s = best_model$lambda.1se)
      selectedLassoFeatures <- names(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])
      selectedLassoFeatures <- selectedLassoFeatures[selectedLassoFeatures != "(Intercept)"]
      
      param <- length(selectedLassoFeatures)
      info(logger, "Lasso Selected features: ", param)
      
      
      if(!is.null(selectedLassoFeatures)){
        lasso_selectedfeatures <- lasso_selectedfeatures |> 
          union_all(
            tibble("concept_id" = names(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1]),
                   "coefficients" = as.numeric(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])
            ) |>
              filter(concept_id != "(Intercept)") |>
              mutate("length_id"  = gsub("_.*","",concept_id)) |>
              mutate("concept_id" = as.numeric(gsub(".*_","",concept_id))) |>
              inner_join(
                ConditionsAndDrugs |>
                  dplyr::select("concept_id","concept_name") |>
                  distinct() |>
                  as_tibble() |> as.data.frame(),
                by = "concept_id"
              ) |>
              mutate("year" = sy, "database" = dbName)
          )
      }
      
      dataMatching <- subpop_data |> 
          addRegion() |> 
          dplyr::select("vac_status", "subject_id", "year_of_birth", "region", all_of(selectedLassoFeatures)) |> 
          mutate(vac_status = case_when(
            vac_status == 0 ~ 1,
            vac_status == 1 ~ 0
          )) |>
          as_tibble()  |>  
          mutate_all(~replace(., is.na(.), 0)) |>
          compute()
      
      tryCatch(
        expr = {
          dataMatched <- matchit(
            vac_status ~ . -subject_id -year_of_birth -region, 
            data  = dataMatching,
            exact = c("year_of_birth","region"),
            caliper  = 0.2,
            method   = "nearest",
            distance = "glm",
            ratio = 5)
        },
        error = function(e){
          info(
            logger = logger, 
            paste("error on year:", sy, ", propensity scores not computed")
          )
        }
      )
              
    if(!"dataMatched" %in% ls()){
      info(logger, "Proceed with only exact matching for region and year of birth")
      dataMatched <- matchit(
        vac_status ~ year_of_birth + region, 
        data  = dataMatching,
        method = "exact",
        ratio = 5)
    } else {
      info(logger, "Proceed with PS and exact matching")
    }
  
    # Save matched cohort
    sub_matched <- as_tibble(match.data(dataMatched)) |> 
      dplyr::select("vac_status", "subject_id", "year_of_birth", "region", "subclass") |>
      mutate(vac_status = case_when(
        vac_status == 0 ~ 1,
        vac_status == 1 ~ 0
      )) |>
      compute()
      n_matched <- as.numeric(sub_matched |> tally() |> pull())
      info(logger = logger, paste0("Matched = ", n_matched))
      
      cdm <- cdm |>
        CDMConnector::insertTable(
          name = paste0("sub_",as.integer(sy),"_matched_cohort"),
          table = sub_matched
        )  
      
      cohort_name <- paste0("sub_",as.integer(sy),"_matched_cohort")
      cdm[[cohort_name]] <- cdm[[cohort_name]] |> 
        left_join(subpopulation |> 
                    dplyr::select(cohort_definition_id, subject_id, cohort_start_date, cohort_end_date, index_date), 
                  by = "subject_id"
        ) |>
        addCohortIntersectDate(
          targetCohortTable = "doses_allvac_cohort",
          indexDate = "cohort_start_date",
          window = c(0, Inf),
          nameStyle = "{value}_dose_1",
          order = "first"
        ) |>
        mutate(index_date = case_when(
          vac_status == 0 ~ index_date,
          vac_status == 1 ~ date_dose_1)) |>
        mutate(cohort_definition_id = as.integer(sy)) |>
        compute(name = cohort_name, temporary = FALSE) |>
        omopgenerics::newCohortTable(cohortSetRef = tibble(cohort_definition_id = as.integer(sy), 
                                                           cohort_name = paste0("sub_",as.integer(sy),"_matched_cohort")
        )
        )
      
      # Eliminate matched individuals from the population cohort
      population <- anti_join(population,cdm[[cohort_name]] |> 
                                dplyr::select(subject_id)
      ) |> 
        anti_join(subpopulation |> 
                    filter(vac_status == 1), 
                  by = "subject_id"
        )
      
      remaining_vac <- population |> filter(vac_status == 1) |> tally() |> pull()
      info(logger = logger, paste0("Remaining vaccinated subjects in population = ", remaining_vac))
      
      # Collect different years matched cohort into one table
      cdm[[name_all]] <- cdm[[name_all]] |> 
        union_all(cdm[[cohort_name]] |> 
                    mutate(cohort_start_date = index_date,
                           cohort_year = cohort_definition_id, 
                           pair_id = as.numeric(subclass) + last_pair, 
                           cohort_definition_id = vac_status + 10) |>
                    dplyr::select(! c(subclass, region, index_date))
        ) |> 
        compute()
      
      last_pair <- cdm[[name_all]] |> 
        dplyr::select(pair_id) |> 
        as.array() |> 
        collect() |> 
        max()
    }
    # Time execution
    end_time <- Sys.time()
    execution_time <- end_time - start_time
    print(execution_time)
  }
  
  info(logger, "RESULTS")
  n_total_matched <- cdm[[name_all]] |> tally() |> pull()
  info(logger, paste0("Number of total matched individuals = ", n_total_matched))
  info(logger, paste0("Number of vac matched individuals = ", n_total_matched/2))
  info(logger, paste0("Number of unvac matched individuals = ", n_total_matched/2))
  
  
  # Separate vac/unvac matched cohorts
  cdm[[name_control]] <- cdm[[name_all]] |> 
    filter(cohort_definition_id == 10) |> 
    left_join(cdm[[name_all]] |>
                 filter(cohort_definition_id == 11) |>
                 dplyr::select(pair_id, cohort_start_date) |> 
                 as_tibble() |>
                 group_by(pair_id) |>
                 mutate(index_date = as.Date.numeric(mean(as.numeric(as.Date(cohort_start_date))), origin = "1970-01-01")) |>
                 ungroup() |>
                 dplyr::select(pair_id, index_date) |>
                 distinct(), 
               by = "pair_id", copy = TRUE) |>
    mutate(cohort_start_date = index_date) |> 
    dplyr::select(cohort_definition_id, subject_id, cohort_start_date, cohort_end_date, pair_id) |>
    compute(name = name_control, temporary = FALSE) |>
    newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(10), cohort_name = c(name_control)),
                   cohortAttritionRef = attrition(cdm[[unvac_cohortname]]) %>% mutate(cohort_definition_id = 10)) |>
    recordCohortAttrition("Restrict to matched individuals")
  
  cdm[[name_treatment]] <- cdm[[name_all]] |> 
    filter(cohort_definition_id == 11) |>
    dplyr::select(cohort_definition_id, subject_id, cohort_start_date, cohort_end_date, pair_id) |>
    compute(name = name_treatment, temporary = FALSE) |>
    newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(11), cohort_name = c(name_treatment)),
                   cohortAttritionRef = attrition(cdm[[vac_cohortname]]) %>% mutate(cohort_definition_id = 11)) |>
    recordCohortAttrition("Restrict to matched individuals")
  
  # Write the results
  write.csv(lasso_selectedfeatures, here(resultsFolder_lasso, paste0("lasso_selectedfeatures_", yy,ss,"_overall_cov",cov,".csv")), row.names = FALSE)
}
