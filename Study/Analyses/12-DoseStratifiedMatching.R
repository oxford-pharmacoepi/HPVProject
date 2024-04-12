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

#cohorts <- readCohortSet(path = here("Cohorts", "HIV_allvac"))
#cdm <- generateCohortSet(cdm = cdm, cohortSet = cohorts, name = c("doses_allvac_cohort"))

# Details
# log_file        <- here(resultsFolder, "log_DoseMatchingData2.txt")
# logger          <- create.logger(logfile = log_file, level = "INFO")
info(logger = logger, "CREATE INITIAL POPULATION")

# Total Population
population <- cdm$vac_cohort
doses <- cdm$doses_allvac_cohort |>
  left_join(population |>
              dplyr::select(subject_id, date_15years), by = "subject_id") |>
  filter(cohort_start_date <= date_15years) |>
  group_by(subject_id) |>
  count()

population <- population |>
  left_join(doses, by = "subject_id") |>
  mutate(num_doses = n) |>
  dplyr::select(! n) |>
  filter(num_doses == 2 || num_doses >= 3) |>
  mutate(treatment = case_when(
    num_doses == 2 ~ 0,
    num_doses > 2 ~ 1)
  )

size_pop <- population |> tally() |> pull()
size_2vacpop <- population |> filter(treatment == 0) |> tally() |> pull()
size_3vacpop <- population |> filter(treatment == 1) |> tally() |> pull()
info(logger = logger, paste0("Size Population = ", size_pop))
info(logger = logger, paste0("Size 2 dose Population = ", size_2vacpop))
info(logger = logger, paste0("Size 3 doses Population = ", size_3vacpop))

# Total conditions and drugs
Conditions <- cdm$condition_occurrence |> 
  rename(subject_id = person_id, concept_id = condition_concept_id, occurrence_start_date = condition_start_date) |>
  left_join(cdm$concept |> dplyr::select(concept_id, concept_name), by = "concept_id") |>
  dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date) |> 
  mutate(occurrence_type = "condition") |> 
  compute()

Drugs <-  cdm$drug_exposure |> 
  rename(subject_id = person_id, concept_id = drug_concept_id, occurrence_start_date = drug_exposure_start_date)|> 
  left_join(cdm$concept |> dplyr::select(concept_id, concept_name), by = "concept_id") |>
  dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date) |> 
  mutate(occurrence_type = "drug") |> 
  compute()

ConditionsAndDrugs <- union_all(Conditions, Drugs) |> 
  mutate(occurrence_flag = 1) |>
  inner_join(population, by = "subject_id") |>
  filter(occurrence_start_date < date_15years) |>
  dplyr::select(subject_id, concept_id, concept_name, occurrence_start_date, occurrence_type, occurrence_flag) |>
  compute()

# Prepare covariates
hiv_status <- population |>
  dplyr::select(subject_id, date_15years) |>
  inner_join(cdm$hiv_status |> 
               dplyr::select(subject_id, cohort_start_date) |> 
               rename(occurrence_start_date = cohort_start_date), by = "subject_id") |>
  filter(occurrence_start_date < date_15years) |>
  group_by(subject_id) |>
  filter(occurrence_start_date == min(occurrence_start_date)) |>
  ungroup() |>
  mutate(hiv_status = 1) |>
  dplyr::select(subject_id, hiv_status) |>
  compute()

papanicolau_smear_test <- population |>
  dplyr::select(subject_id, date_15years) |>
  inner_join(cdm$papanicolau_smear_testing |> 
               dplyr::select(subject_id, cohort_start_date) |> 
               rename(occurrence_start_date = cohort_start_date), by = "subject_id") |>
  filter(occurrence_start_date < date_15years) |>
  group_by(subject_id) |>
  filter(occurrence_start_date == min(occurrence_start_date)) |>
  ungroup() |>
  mutate(papanicolau = 1) |>
  dplyr::select(subject_id, papanicolau) |>
  compute()

vaccinations <- population |>
  dplyr::select(subject_id, date_15years) |>
  inner_join(cdm$vaccinations |> 
               dplyr::select(subject_id, cohort_start_date) |> 
               rename(occurrence_start_date = cohort_start_date), by = "subject_id") |>
  filter(occurrence_start_date < date_15years) |>
  group_by(subject_id) |>
  count() |>
  rename(prev_vaccinations = n) |>
  compute()

cytology <- population |>
  dplyr::select(subject_id, date_15years) |>
  inner_join(cdm$cytology |> 
               dplyr::select(subject_id, cohort_start_date) |> 
               rename(occurrence_start_date = cohort_start_date), by = "subject_id") |>
  filter(occurrence_start_date < date_15years) |>
  group_by(subject_id) |>
  filter(occurrence_start_date == min(occurrence_start_date)) |>
  ungroup() |>
  mutate(cytology = 1) |>
  dplyr::select(subject_id, cytology) |>
  compute()

# Results file and cohorts
lasso_selectedfeatures <- tibble("concept_id" = as.numeric(),
                                 "coefficients" = as.numeric(),
                                 "length_id" = as.character(),
                                 "concept_name" = as.character(),
                                 "year"     = as.numeric(),
                                 "database" = as.character()
)

dose_1vs2_matched_cohort <- tibble("cohort_definition_id" = as.numeric(),
                              "subject_id" = as.numeric(),
                              "cohort_year" = as.numeric(),
                              "vac_status" = as.numeric(1),
                              "cohort_start_date" =  as.Date(character()),
                              "cohort_end_date" =  as.Date(character()), 
                              "pair_id" = as.numeric()
)


cdm <- cdm |>
  CDMConnector::insertTable(
    name = "dose_1vs2_matched_cohort",
    table = dose_1vs2_matched_cohort,
    overwrite = TRUE
  ) 

# Subpopulation
# initialisation of parameters
last_pair <- 0

for(sy in 2008:2023){#(year(studyEndDate)-year(studyStartDate))){
  #sy <- 2009
  # Execution time
  start_time <- Sys.time()
  print(sy)
  info(logger = logger, paste0("CREATE SUBPOPULATION ", sy))
  
  subpopulation <-  population |> filter(cohort_start_date <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                         cohort_end_date >= as.POSIXct(paste0(as.integer(sy),"-12-31")),
                                         date_9years <= as.POSIXct(paste0(as.integer(sy),"-01-01")),
                                         date_15years >= as.POSIXct(paste0(as.integer(sy),"-01-01"))) |> 
    inner_join(cdm$firstdose_cohort |> 
                 filter(cohort_start_date >= as.POSIXct(paste0(as.integer(sy),"-01-01")), 
                        cohort_start_date <= as.POSIXct(paste0(as.integer(sy),"-12-31"))) |> 
                 rename(index_date = cohort_start_date) |>
                 dplyr::select(subject_id, index_date)
    ) |> 
    compute()
  
  
  # Set new index date: 01/01/2008
  # subpopulation <- subpopulation |> 
  #   mutate(index_date = as.POSIXct(paste0(as.integer(sy),"-01-01"))) |>
  #   compute()
  subpopulation |> tally()
  
  # Checks:
  # in observation at index date
  subpopulation <- subpopulation |> 
    filter(cohort_start_date <= index_date & cohort_end_date >= index_date) |>
    compute()
  subpopulation |> tally()
  
  # 365 prior information from index_date
  subpopulation <- subpopulation |> 
    addPriorObservation(indexDate = "index_date", priorObservationName = "prior_observation_indexdate") |> 
    filter(prior_observation_indexdate >= 365) |>
    compute()
  subpopulation |> tally()
  
  # Add year of birth and prior visits
  subpopulation <- subpopulation |> 
    inner_join(cdm$person |> 
                 dplyr::select(person_id, year_of_birth) |> 
                 rename(subject_id = person_id)
    ) |>
    addVisitsPriorYear() |>
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
      anti_join(subpopulation |> 
                  filter(treatment == 0), 
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
    best_model <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standarize = TRUE, nforlds = 5, family = "binomial")
    # how do we know this is the best lambda???????????????????????
    # Binomial??
    
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
    
    # Apply matching
    # Add Previous visits
    
    dataMatching <- subpop_data |> 
      dplyr::select("treatment", "subject_id", "year_of_birth", all_of(selectedLassoFeatures)) |> 
      as_tibble()  |>  
      mutate_all(~replace(., is.na(.), 0)) |>
      compute()
    
    Match
    if(is.null(selectedLassoFeatures)){
      print("is null")
      region_type <- "num"       # Matching no accepta "character" variables
      
      dataMatching <- subpop_data |>
        addRegion() |>
        dplyr::select("treatment", "subject_id", "year_of_birth", "region", all_of(selectedLassoFeatures)) |>
        as_tibble()  |>
        mutate_all(~replace(., is.na(.), 0)) |>
        compute()
      
      dataMatched <- Match(Tr = dataMatching$treatment,
                           X = subset(dataMatching, select = c(year_of_birth, region)),
                           M = 1,
                           exact = TRUE,
                           ties = FALSE,
                           replace = FALSE)
      
      # Add subclass
      treated <- dataMatching[dataMatched$index.treated, ] |>
        mutate(subclass = 1:dataMatched$wnobs)
      control <- dataMatching[dataMatched$index.control, ] |>
        mutate(subclass = 1:dataMatched$wnobs)
      
      # Save matched cohorts
      sub_matched <- union_all(treated, control)
      
      # dataMatched <- matchit(vac_status ~ . - subject_id,
      #                        data = dataMatching,
      #                        method = "exact",
      #                        ratio = 1,
      #                        ties = FALSE
      #                        )
    } else {
      print("not null")
      region_type <- "char"
      
      dataMatching <- subpop_data |> 
        addRegion() |> 
        dplyr::select("treatment", "subject_id", "year_of_birth", "region", all_of(selectedLassoFeatures)) |> 
        as_tibble()  |>  
        mutate_all(~replace(., is.na(.), 0)) |>
        compute()
      
      dataMatched <- matchit(treatment ~ . - subject_id,
                             data = dataMatching, 
                             exact = c("year_of_birth","region"), 
                             caliper = 0.2, 
                             method = "nearest", 
                             distance = "glm",
                             ratio = 1
      )
      # Save matched cohort
      sub_matched <- as_tibble(match.data(dataMatched)) |> 
        dplyr::select("treatment", "subject_id", "year_of_birth", "region", "subclass") |>
        compute()
      
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
      mutate(cohort_definition_id = as.integer(sy))
    
    newCohortTable(cdm[[cohort_name]], 
                   cohortSetRef = tibble(cohort_definition_id = as.integer(sy), 
                                         cohort_name = paste0("sub_",as.integer(sy),"_matched_cohort")
                   )
    )
    
    # Eliminate matched individuals from the population cohort
    population <- anti_join(population, subpopulation) 
    remaining_pop <- population |> tally() |> pull()
    info(logger = logger, paste0("Remaining subjects in population = ", remaining_pop))
    
    # Collect different years matched cohort into one table
    cdm$dose_1vs2_matched_cohort <- cdm$dose_1vs2_matched_cohort |> 
      union_all(cdm[[cohort_name]] |> 
                  mutate(cohort_year = cohort_definition_id, 
                         pair_id = as.numeric(subclass) + last_pair, 
                         cohort_definition_id = treatment + 10) |>
                  dplyr::select(! c(subclass, region))
      ) |> 
      compute()
    
    last_pair <- cdm$dose_1vs2_matched_cohort |> 
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
n_doses_matched <- cdm$dose_1vs2_matched_cohort |> tally() |> pull()
info(logger, paste0("Number of doses matched individuals = ", n_doses_matched))
info(logger, paste0("Number of 2 dose matched individuals = ", n_doses_matched/2))
info(logger, paste0("Number of 3 dose matched individuals = ", n_doses_matched/2))


# Separate 1 dose/2 dose matched cohorts
cdm$dose2_matched_cohort <- cdm$dose_1vs2_matched_cohort |> 
  filter(cohort_definition_id == 10) |> 
  compute(name = "dose2_matched_cohort", temporary = FALSE) |>
  newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(10), cohort_name = c("dose2_matched_cohort")),
                 cohortAttritionRef = attrition(cdm$vac_cohort) %>% mutate(cohort_definition_id = 10)) |>
  recordCohortAttrition("Exclude not matched individuals") |>
  compute(name = "dose2_matched_cohort", temporary = FALSE)

cdm$doses3_matched_cohort <- cdm$dose_1vs2_matched_cohort |> 
  filter(cohort_definition_id == 11) |>
  compute(name = "doses3_matched_cohort", temporary = FALSE) |>
  newCohortTable(cohortSetRef = tibble(cohort_definition_id = c(11), cohort_name = c("doses3_matched_cohort")),
                 cohortAttritionRef = attrition(cdm$vac_cohort) %>% mutate(cohort_definition_id = 11)) |>
  recordCohortAttrition("Exclude not matched individuals")

# Change cohort start date for the 15 years
# cdm$vac_cohort <- cdm$vac_cohort |>
#   mutate(cohort_start_date = as.Date(date_15years)) |>
#   dplyr::select(! date_15years) |>
#   compute(name = "vac_cohort", temporary = FALSE)
# 
# cdm$unvac_cohort <- cdm$unvac_cohort |>
#   mutate(cohort_start_date = as.Date(date_15years)) |>
#   dplyr::select(! date_15years) |>
#   compute(name = "unvac_cohort", temporary = FALSE)

# Write the results
write.csv(cdm$dose_1vs2_matched_cohort |> summary(),
          paste0(resultsFolder,"/attrition_dose2_matched_cohort_",cdmSchema,".csv"), row.names = FALSE)

write.csv(lasso_selectedfeatures,
          paste0(resultsFolder,"/lasso_selectedfeatures_dose2_",cdmSchema,".csv"), row.names = FALSE)

# Save attritions
info(logger, "SAVE COHORT ATTRITIONS")

matched_dose2_attrition <- attrition(cdm$dose2_matched_cohort)
matched_doses3_attrition <- attrition(cdm$doses3_matched_cohort)

write.csv(matched_doses23_attrition, paste0(resultsFolder,"/matched_dose2_attrition_",cdmSchema,".csv"), row.names = FALSE)
write.csv(matched_doses3_attrition, paste0(resultsFolder,"/matched_doses3_attrition_",cdmSchema,".csv"), row.names = FALSE)
