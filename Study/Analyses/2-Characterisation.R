library(reshape2)
library(ggplot2)
library(stringr)

source(here("Analyses", "Figures_functions.R"))

for (cov in coverages) {
  # Overall
  cdm <- omopgenerics::bind(
    cdm[[paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")]],
    cdm[[paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")]],
    name = paste0("hpv_crude_overall_",tolower(sex_cohort),"_cohorts")
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("vac_",yy,ss,"_coverage_",cov,"_any_dose")]],
    cdm[[paste0("unvac_",yy,ss,"_coverage_",cov,"_any_dose")]],
    name = paste0("hpv_matched_overall_",tolower(sex_cohort),"_cohorts")
  )

# Doses stratified
  cdm <- omopgenerics::bind(
    cdm[[paste0("not_matched_vaccinated_",yy,ss,"_cov",cov,"_1dose")]],
    cdm[[paste0("not_matched_vaccinated_",yy,ss,"_cov",cov,"_23dose")]],
    name = paste0("hpv_crude_1vs23dose_",tolower(sex_cohort),"_cohorts")
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("matched_vaccinated_",yy,ss,"_coverage_",cov,"_1dose")]],
    cdm[[paste0("matched_vaccinated_",yy,ss,"_coverage_",cov,"_23dose")]],
    name = paste0("hpv_matched_1vs23dose_",tolower(sex_cohort),"_cohorts")
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("not_matched_vaccinated_",yy,ss,"_coverage_",cov,"_2dose")]],
    cdm[[paste0("not_matched_vaccinated_",yy,ss,"_coverage_",cov,"_3dose")]],
    name = paste0("hpv_crude_2vs3dose_",tolower(sex_cohort),"_cohorts")
  )
  cdm <- omopgenerics::bind(
    cdm[[paste0("matched_vaccinated_",yy,ss,"_coverage_",cov,"_2dose")]],
    cdm[[paste0("matched_vaccinated_",yy,ss,"_coverage_",cov,"_3dose")]],
    name = paste0("hpv_matched_2vs3dose_",tolower(sex_cohort),"_cohorts")
  )

    # Before matching overall
    to_crude_overall_female <- to_Charac(cdm[[paste0("hpv_crude_overall_",tolower(sex_cohort),"_cohorts")]], sex = ss, index_date = TRUE)
    write.csv(to_crude_overall_female, here(resultsFolder_char, paste0("to_crude_overall_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    
    lsc_crude_overall_female <- summariseLargeScaleCharacteristics(cdm[[paste0("hpv_crude_overall_",tolower(sex_cohort),"_cohorts")]],
                                                                   window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                   eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                   minimumFrequency = 0
    )
    write.csv(lsc_crude_overall_female, here(resultsFolder_char,paste0("lsc_crude_overall_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    # After matching overall
    to_ratio_matched_overall_female <- to_Charac(cdm[[paste0("hpv_matched_overall_",tolower(sex_cohort),"_cohorts")]], sex = ss, index_date = TRUE)
    write.csv(to_ratio_matched_overall_female, here(resultsFolder_char, paste0("to_matched_overall_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    
    lsc_ratio_matched_overall_female <- summariseLargeScaleCharacteristics(cdm[[paste0("hpv_matched_overall_",tolower(sex_cohort),"_cohorts")]], 
                                                                           window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                           eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                           minimumFrequency = 0
    )
    write.csv(lsc_ratio_matched_overall_female, here(resultsFolder_char, paste0("lsc_matched_overall_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    # Before matching 1vs23 doses
    to_crude_1vs23dose_female <- to_Charac(cdm[[paste0("hpv_crude_1vs23dose_",tolower(sex_cohort),"_cohorts")]] |>
                                             mutate(cohort_start_date = date_dose_1), sex = ss, index_date = TRUE)
    write.csv(to_crude_1vs23dose_female, here(resultsFolder_char, paste0("to_crude_1vs23dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    
    lsc_crude_1vs23dose_female <- summariseLargeScaleCharacteristics(cdm[[paste0("hpv_crude_1vs23dose_",tolower(sex_cohort),"_cohorts")]] |>
                                                                       mutate(cohort_start_date = date_dose_1), 
                                                                     window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                     eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                     minimumFrequency = 0
    )
    write.csv(lsc_crude_1vs23dose_female, here(resultsFolder_char, paste0("lsc_crude_1vs23dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    # After matching 1vs23 doses
    to_matched_1vs23dose_female <- to_Charac(cdm[[paste0("hpv_matched_1vs23dose_",tolower(sex_cohort),"_cohorts")]], sex = ss, index_date = TRUE)
    write.csv(to_matched_1vs23dose_female, here(resultsFolder_char, paste0("to_matched_1vs23dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    
    lsc_matched_1vs23dose_female <- summariseLargeScaleCharacteristics(cdm[[paste0("hpv_matched_1vs23dose_",tolower(sex_cohort),"_cohorts")]],
                                                                       window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                       eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                       minimumFrequency = 0
    )
    write.csv(lsc_matched_1vs23dose_female, here(resultsFolder_char, paste0("lsc_matched_1vs23dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    # Before matching 2vs3 doses
    to_crude_2vs3dose_female <- to_Charac(cdm[[paste0("hpv_crude_2vs3dose_",tolower(sex_cohort),"_cohorts")]] |>
                                            mutate(cohort_start_date = date_dose_1), sex = ss, index_date = TRUE)
    write.csv(to_crude_2vs3dose_female, here(resultsFolder_char, paste0("to_crude_2vs3dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    
    lsc_crude_2vs3dose_female <- summariseLargeScaleCharacteristics(cdm[[paste0("hpv_crude_2vs3dose_",tolower(sex_cohort),"_cohorts")]] |>
                                                                      mutate(cohort_start_date = date_dose_1), 
                                                                    window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                    eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                    minimumFrequency = 0
    )
    write.csv(lsc_crude_2vs3dose_female, here(resultsFolder_char, paste0("lsc_crude_2vs3dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    # After matching 2vs3 doses
    to_matched_2vs3dose_female <- to_Charac(cdm[[paste0("hpv_matched_2vs3dose_",tolower(sex_cohort),"_cohorts")]], sex = ss, index_date = TRUE)
    write.csv(to_matched_2vs3dose_female, here(resultsFolder_char, paste0("to_matched_2vs3dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)
    
    lsc_matched_2vs3dose_female <- summariseLargeScaleCharacteristics(cdm[[paste0("hpv_matched_2vs3dose_",tolower(sex_cohort),"_cohorts")]],
                                                                      window = list(c(-Inf, -366), c(-365, -31), c(-31, -1)), 
                                                                      eventInWindow = c("condition_occurrence","drug_exposure"),
                                                                      minimumFrequency = 0
    )
    write.csv(lsc_matched_2vs3dose_female, here(resultsFolder_char, paste0("lsc_matched_2vs3dose_",tolower(sex_cohort),"_cov",cov,".csv")), row.names = FALSE)

  }