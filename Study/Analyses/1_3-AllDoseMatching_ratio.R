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
source(here("Analyses","FuncionsMAH.R"))

# source(here("Analyses","Functions 2.R"))
info(logger = logger, "CREATE INITIAL POPULATION")

for(cov in coverages) {
  info(logger, paste0("Coverage: ", cov))
  
  vac_cohortname <- paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")
  cdm$vac_cohort <- cdm[[vac_cohortname]] |>
    addDateOfBirth() |>
    mutate(date_9years = date_of_birth + years(9)) |>
    mutate(!!paste0("date_",yy,"years") := date_of_birth + years(yy)) |>
    compute(name = "vac_cohort", temporary = FALSE) |>
    newCohortTable(cohortSetRef = settings(cdm[[vac_cohortname]]) |> mutate(cohort_name == "vac_cohort"),
                   cohortAttritionRef = attrition(cdm[[vac_cohortname]]))
  
  # Run for the 2 comparisons
  for(comparison in c(1,2)) {
    info(logger, paste0("Comparison: ", comparison))
    # Total Population
    population <- cdm$vac_cohort %>%
      addCohortIntersectDate(
        targetCohortTable = "doses_allvac_cohort",
        indexDate = "cohort_start_date",
        window = c(0, Inf),
        nameStyle = "{value}_dose_1",
        order = "first"
      ) %>%
      addCohortIntersectCount(
        targetCohortTable = "doses_allvac_cohort",
        nameStyle = "dose_count",
        indexDate = paste0("date_",yy,"years"),
        window = c(-Inf,0)
      ) |>
      filter(dose_count == comparison || dose_count > comparison) |>
      mutate(treatment = case_when(
        dose_count == comparison ~ 0,
        dose_count > comparison ~ 1)
      )
    info(logger, "Cohort Names")
    name_control <- case_when(
      comparison == 1 ~ paste0("matched_vaccinated_",yy,ss,"_coverage_", cov, "_1dose"),
      comparison == 2 ~ paste0("matched_vaccinated_",yy,ss,"_coverage_", cov, "_2dose")
    )
    
    name_treatment <- case_when(
      comparison == 1 ~ paste0("matched_vaccinated_",yy,ss,"_coverage_", cov, "_23dose"),
      comparison == 2 ~ paste0("matched_vaccinated_",yy,ss,"_coverage_", cov, "_3dose")
    )
    
    name_all <- case_when(
      comparison == 1 ~ paste0("matched_vaccinated_",yy,ss,"_coverage_", cov, "_1vs23dose"),
      comparison == 2 ~ paste0("matched_vaccinated_",yy,ss,"_coverage_", cov, "_2vs3dose")
    )
    
    info(logger, "Instantiation of cohorts")
    # Unmatched dose cohorts
    cdm[[paste0("not_",name_control)]] <- population |> 
      filter(treatment == 0) |> 
      compute(name = paste0("not_",name_control), temporary = FALSE) |>
      newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(2), cohort_name = c(paste0("not_",name_control))),
                     cohortAttritionRef = attrition(cdm$vac_cohort) %>% mutate(cohort_definition_id = 2)) |>
      recordCohortAttrition(paste0("Restrict to ", comparison, " dose"))
    
    cdm[[paste0("not_",name_treatment)]] <- population |> 
      filter(treatment == 1) |>
      compute(name = paste0("not_",name_treatment), temporary = FALSE) |>
      newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(2), cohort_name = c(paste0("not_",name_treatment))),
                     cohortAttritionRef = attrition(cdm$vac_cohort) %>% mutate(cohort_definition_id = 2)) |>
      recordCohortAttrition(paste0("Restrict to more than ", comparison, " doses"))
    
    # Matched dose cohorts
    cdm[[name_control]] <- population |> 
      filter(treatment == 0) |> 
      compute(name = name_control, temporary = FALSE) |>
      newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(2), cohort_name = c(name_control)),
                     cohortAttritionRef = attrition(cdm$vac_cohort) %>% mutate(cohort_definition_id = 2)) |>
      recordCohortAttrition(paste0("Restrict to ", comparison, " dose"))
    
    cdm[[name_treatment]] <- population |> 
      filter(treatment == 1) |>
      compute(name = name_treatment, temporary = FALSE) |>
      newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(2), cohort_name = c(name_treatment)),
                     cohortAttritionRef = attrition(cdm$vac_cohort) %>% mutate(cohort_definition_id = 2)) |>
      recordCohortAttrition(paste0("Restrict to more than ", comparison, " doses"))
    
    size_pop <- population |> tally() |> pull()
    size_1vacpop <- population |> filter(treatment == 0) |> tally() |> pull()
    size_2vacpop <- population |> filter(treatment == 1) |> tally() |> pull()
    info(logger = logger, paste0("Size Population = ", size_pop))
    info(logger = logger, paste0("Size comparator dose Population = ", size_1vacpop))
    info(logger = logger, paste0("Size target doses Population = ", size_2vacpop))
    
    info(logger, "Instantiation of conditions and drugs")
    # Total conditions and drugs
    Conditions <- cdm$condition_occurrence |>
      inner_join(population |> dplyr::select("person_id" = "subject_id", date_dose_1), by = "person_id") |>
      compute() |>
      rename(subject_id = person_id, concept_id = condition_concept_id, occurrence_start_date = condition_start_date) |>
      left_join(cdm$concept |> dplyr::select(concept_id, concept_name), by = "concept_id") |>
      dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, date_dose_1) |>
      mutate(occurrence_type = "condition") |>
      compute()
    
    Drugs <-  cdm$drug_exposure |>
      inner_join(population |> dplyr::select("person_id" = "subject_id", date_dose_1), by = "person_id") |>
      compute() |>
      rename(subject_id = person_id, concept_id = drug_concept_id, occurrence_start_date = drug_exposure_start_date)|>
      left_join(cdm$concept |> dplyr::select(concept_id, concept_name), by = "concept_id") |>
      dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, date_dose_1) |>
      mutate(occurrence_type = "drug") |>
      compute()
    
    ConditionsAndDrugs <- union_all(Conditions, Drugs) |>
      mutate(occurrence_flag = 1) |>
      filter(occurrence_start_date < date_dose_1) |>
      dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, occurrence_type, occurrence_flag) |>
      compute()
    
    info(logger, "Prepare covariates")
    # Prepare covariates
    hiv_status <- population |>
      addCohortIntersectFlag(targetCohortTable = "hiv_status",
                             indexDate = "date_dose_1",
                             window = c(-Inf, -1),
                             nameStyle = "hiv_status_flag"
                             ) |>
      dplyr::select(subject_id, hiv_status_flag) |>
      compute()
    
    papanicolau_smear_test <- population |>
      addCohortIntersectFlag(targetCohortTable = "smear_test",
                             indexDate = "date_dose_1",
                             window = c(-Inf, -1),
                             nameStyle = "papanicolau_smear_testing_flag"
      ) |>
      dplyr::select(subject_id, papanicolau_smear_testing_flag) |>
      compute()
    
    vaccinations <- population |>
      addCohortIntersectCount(targetCohortTable = "vaccinations",
                             indexDate = "date_dose_1",
                             window = c(-Inf, -1),
                             nameStyle = "vaccinations_flag"
      ) |>
      dplyr::select(subject_id, vaccinations_flag) |>
      compute()
    
    cytology <- population |>
      addCohortIntersectFlag(targetCohortTable = "cervical_screening",
                             indexDate = "date_dose_1",
                             window = c(-Inf, -1),
                             nameStyle = "cytology_flag"
      ) |>
      dplyr::select(subject_id, cytology_flag) |>
      compute()
    
    # Results file and cohorts
    lasso_selectedfeatures <- tibble("concept_id" = as.numeric(),
                                     "coefficients" = as.numeric(),
                                     "length_id" = as.character(),
                                     "concept_name" = as.character(),
                                     "year"     = as.numeric(),
                                     "database" = as.character()
    )
    
    dose_matched_cohort <- tibble("cohort_definition_id" = as.numeric(),
                                  "subject_id" = as.numeric(),
                                  "cohort_year" = as.numeric(),
                                  "treatment" = as.numeric(),
                                  "cohort_start_date" =  as.Date(character()),
                                  "cohort_end_date" =  as.Date(character()), 
                                  "pair_id" = as.numeric()
    )
    
    
    cdm <- cdm |>
      CDMConnector::insertTable(
        name = name_all,
        table = dose_matched_cohort,
        overwrite = TRUE
      ) 
    
    # initialisation of parameters
    last_pair <- 0
    
    for(sy in 2008:2023){
      # Execution time
      start_time <- Sys.time()
      print(sy)
      info(logger = logger, paste0("CREATE SUBPOPULATION ", sy))
      
      subpopulation <-  population |> filter(cohort_start_date <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                             cohort_end_date >= as.POSIXct(paste0(as.integer(sy),"-12-31")),
                                             date_9years <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                             !!sym(paste0("date_",yy,"years")) >= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                             date_dose_1 >= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                             date_dose_1 <= as.POSIXct(paste0(as.integer(sy),"-12-31"))
      ) |> 
        rename(index_date = date_dose_1) |>
        compute()
      
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
      size_1vacsubpopulation <- as.numeric(count(distinct(subpopulation |> filter(treatment == 0), subject_id)) |> compute() |> pull())
      size_2vacsubpopulation <- as.numeric(count(distinct(subpopulation |> filter(treatment == 1), subject_id)) |> compute() |> pull())
      info(logger = logger, paste0("Size subpopulation = ", size_subpopulation))
      info(logger = logger, paste0("Size 1 dose subpopulation = ", size_1vacsubpopulation))
      info(logger = logger, paste0("Size 2 doses subpopulation = ", size_2vacsubpopulation))
      
      if (size_1vacsubpopulation < 5 | size_2vacsubpopulation < 5) {
        info(logger = logger, "Not enough participants to proceed with the lasso / matching")
        population <- population |> 
          anti_join(subpopulation, 
                    by = "subject_id"
          )
      }else{
        
        # Subpopulation as Lasso Input
        sub_ConditionsAndDrugs <- subpopulation |> 
          dplyr::select(subject_id, index_date, year_of_birth, treatment) |> 
          inner_join(
            # Occurrence Record before index date
            ConditionsAndDrugs, by = "subject_id"
          )
        
        size <- sub_ConditionsAndDrugs |> distinct(concept_id) |> tally() |> pull()
        info(logger, "Conditions only for subpopulation: ", size)
        
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
          dplyr::select(subject_id, index_date, year_of_birth, treatment) |> 
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
          dplyr::select(! c(subject_id, index_date, treatment, year_of_birth)) |> 
          as_tibble() |>
          mutate_all(~replace(., is.na(.), 0)) |> 
          data.matrix() 
        info(logger = logger, paste0("Number of selected features = ", ncol(x)))
        
        # Vac_status = 0/1
        y <- in_data |> 
          dplyr::select(treatment) |> 
          compute() |> 
          as_tibble() |> 
          data.matrix()
        
        # Model definition
        lambdas <- 10^seq(2, -3, by = -.1) # magnitude order lambda
        best_model <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standarize = TRUE, nfolds = 5, family = "binomial")
        
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
          dplyr::select("treatment", "subject_id", "year_of_birth", "region", all_of(selectedLassoFeatures)) |> 
          as_tibble()  |>  
          mutate_all(~replace(., is.na(.), 0)) |>
          compute()
        
        if("dataMatched" %in% ls()){
          rm(dataMatched)
        }
        
        tryCatch(
          expr = {
            dataMatched <- matchit(
              treatment ~ . -subject_id -year_of_birth -region, 
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
        
        if("dataMatched" %in% ls()){
          info(logger, "Proceed with PS and exact matching")
          # Save matched cohort
          sub_matched <- as_tibble(match.data(dataMatched)) |> 
            dplyr::select("treatment", "subject_id", "year_of_birth", "region", "subclass") |>
            compute()
          
        } else {
          
          tryCatch(
            expr = {
              dataMatched <- matchit(
                treatment ~ year_of_birth + region, 
                data  = dataMatching,
                method = "exact",
                ratio = 5)
            },
            error = function(e){
              info(
                logger = logger, 
                paste("error on year:", sy, ", no exact matched were found")
              )
            }
          )
          
          if("dataMatched" %in% ls()){
            info(logger, "Proceed with only exact matching for region and year of birth")
            # Save matched cohort
            sub_matched <- as_tibble(match.data(dataMatched)) |> 
              dplyr::select("treatment", "subject_id", "year_of_birth", "region", "subclass") |>
              compute()
            
          } else {
            sub_matched <- tibble("subject_id" = as.numeric(),
                                  "treatment" = as.numeric(),
                                  "subclass" = as.numeric(),
                                  "region" = as.character()
            )

          }

        }
        
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
          mutate(cohort_definition_id = as.integer(sy)) |>
          compute(name = cohort_name, temporary = FALSE) |>
          newCohortTable(cohortSetRef = tibble(cohort_definition_id = as.integer(sy), 
                                               cohort_name = paste0("sub_",as.integer(sy),"_matched_cohort")
                                               ),
                         cohortAttritionRef = attrition(cdm[[vac_cohortname]]) |> mutate(cohort_definition_id = as.integer(sy))
                         ) |>
          recordCohortAttrition(paste0("Year", sy, "matched"))
        
        # Eliminate matched individuals from the population cohort
        population <- anti_join(population, subpopulation) 
        remaining_pop <- population |> tally() |> pull()
        info(logger = logger, paste0("Remaining subjects in population = ", remaining_pop))
        
        # Collect different years matched cohort into one table
        cdm[[name_all]] <- cdm[[name_all]] |> 
          union_all(cdm[[cohort_name]] |> 
                      mutate(cohort_start_date = index_date,
                             cohort_year = cohort_definition_id, 
                             pair_id = as.numeric(subclass) + last_pair, 
                             cohort_definition_id = treatment + 10) |>
                      dplyr::select(! c(subclass, region, index_date))
          ) |> 
          compute(name = name_all, temporary = FALSE)
        
        tt <- cdm[[name_all]] |> 
          dplyr::select(pair_id) |> 
          as.array()
        
        if (is.na(nrow(tt))) {
          last_pair <- 0
        } else {
          last_pair <- max(tt)
        }
        
        print(last_pair)
      }
      # Time execution
      end_time <- Sys.time()
      execution_time <- end_time - start_time
      print(execution_time)
    }
    
    info(logger, "RESULTS")
    n_doses_matched <- cdm[[name_all]] |> tally() |> pull()
    info(logger, paste0("Number of doses matched individuals = ", n_doses_matched))
    info(logger, paste0("Number of 1 dose matched individuals = ", cdm[[name_all]] |> filter(treatment == 0) |> tally() |> pull()))
    info(logger, paste0("Number of 2 dose matched individuals = ", cdm[[name_all]] |> filter(treatment == 1) |> tally() |> pull()))
    
    
    # Separate 1 dose/2 dose matched cohorts
    cdm[[name_control]] <- cdm[[name_all]] |> 
      filter(cohort_definition_id == 10) |> 
      dplyr::select(cohort_definition_id, subject_id, cohort_start_date, cohort_end_date, pair_id) |>
      compute(name = name_control, temporary = FALSE) |>
      newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(10), cohort_name = c(name_control)),
                     cohortAttritionRef = attrition(cdm[[name_control]]) %>% mutate(cohort_definition_id = 10)) |>
      recordCohortAttrition("Restrict to matched individuals")
    
    cdm[[name_treatment]] <- cdm[[name_all]] |> 
      filter(cohort_definition_id == 11) |>
      dplyr::select(cohort_definition_id, subject_id, cohort_start_date, cohort_end_date, pair_id) |>
      compute(name = name_treatment, temporary = FALSE) |>
      newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(11), cohort_name = c(name_treatment)),
                     cohortAttritionRef = attrition(cdm[[name_treatment]]) %>% mutate(cohort_definition_id = 11)) |>
      recordCohortAttrition("Restrict to matched individuals")
    
    # Write the results
    write.csv(lasso_selectedfeatures,
              here(resultsFolder_lasso, paste0("lasso_selectedfeatures_", yy,ss,"_comparison_", comparison, "_cov",cov,".csv")), row.names = FALSE)
    
    # Save attritions
    info(logger, "SAVE COHORT ATTRITIONS")
    
    matched_dose1_attrition <- attrition(cdm[[name_control]])
    matched_doses23_attrition <- attrition(cdm[[name_treatment]])
    
    write.csv(matched_dose1_attrition, here(resultsFolder_att, paste0("/", name_control,"_cov",cov,"_attrition_",cdmSchema,".csv")), row.names = FALSE)
    write.csv(matched_doses23_attrition, here(resultsFolder_att, paste0("/", name_treatment,"_cov",cov,"_attrition_",cdmSchema,".csv")), row.names = FALSE)
    
  }
}
