library(survival)

cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema, 
  writeSchema = writeSchema, 
  cdmName = dbName, 
  achillesSchema = achillesSchema, 
  cohortTables = c("vac_cohort", "unvac_cohort", "firstdose_cohort", "doses_allvac_cohort", "total_vac_matched_cohort", "total_unvac_matched_cohort", "phenotyping_hpv_")
)

info(logger, "PREPARE OUTCOME TABLE")

# Put together matched cohorts
cdm$total_matched_cohort <- union_all(cdm$total_vac_matched_cohort, cdm$total_unvac_matched_cohort) |>
  compute(name = "total_matched_cohort", temporary = FALSE) |>
  select(! c(region, subclass, cohort_year, index_date)) 

# Create next dose
next_dose <- cdm$total_matched_cohort |>
  left_join(cdm$doses_allvac_cohort |>
               select(subject_id, cohort_start_date) |>
               rename(next_dose_date = cohort_start_date), by = "subject_id") |>
  group_by(subject_id) |>
  filter(next_dose_date > cohort_start_date | is.na(next_dose_date)) |>
  filter(next_dose_date == min(next_dose_date) | is.na(next_dose_date)) |>
  ungroup() |>
  select(subject_id, next_dose_date)

# Add next dose to outcome table
cdm$outcome_table <- cdm$total_matched_cohort |>
  left_join(next_dose, by = "subject_id") |>
  compute(name = "outcome_table", temporary = FALSE)

# Add death date
cdm$outcome_table <- cdm$outcome_table |>
  left_join(cdm$death |>
              select(person_id, death_date) |>
              rename(subject_id = person_id)) |>
  compute(name = "outcome_table", temporary = FALSE)

# Add outcome column
cdm$outcome_table <- cdm$outcome_table |>
  left_join(cdm$phenotyping_hpv_ |>
              select(subject_id, cohort_start_date) |>
              rename(outcome_date = cohort_start_date), by = "subject_id") |>
  mutate(outcome_flag = case_when(
    is.na(outcome_date) ~ 0,
    !is.na(outcome_date) ~ 1
    )) |>
  compute(name = "outcome_table", temporary = FALSE)

# Add event column: death, end obs, outcome, next dose
cdm$outcome_table <- cdm$outcome_table |>
  mutate(event_date = pmin(death_date, cohort_end_date, outcome_date, next_dose_date, na.rm = TRUE)) |>
  compute(name = "outcome_table", temporary = FALSE) |>
  mutate(event_type = case_when(
    event_date == death_date ~ "death",
    event_date == cohort_end_date ~ "end of observation",
    event_date == outcome_date ~ "outcome",
    event_date == next_dose_date ~ "next dose"
    )
  )

# Add Censor pair column
cdm$outcome_table <- cdm$outcome_table |>
  group_by(pair_id) |>
  mutate(censored_date = min(event_date, na.rm = FALSE)) |>
  ungroup() |>
  compute(name = "outcome_table", temporary = FALSE) |>
  mutate(censor_reason = case_when(
    censored_date == event_date ~ event_type,
    censored_date != event_date ~ "pair censored"
  ))

# Add time to event column
cdm$outcome_table <- cdm$outcome_table |>
  mutate(time_to_event = censored_date - cohort_start_date) |>
  compute(name = "outcome_table", temporary = FALSE)

info(logger, "COMPUTE IR")

# Vaccinated rates
pop_5y_vac <- cdm$outcome_table |> filter(vac_status == 1) |>
  mutate(time_to_event = case_when(
    time_to_event > 365*5 ~ 365*5,
    time_to_event <= 365*5 ~ time_to_event
  ))

den_5y_vac <- pop_5y_vac |> pull(time_to_event) |> sum()/365

num_5y_vac <- pop_5y_vac |>
  filter(censor_reason == "outcome" & time_to_event < 5*365) |>
  tally() |> pull()


# Unvaccinated rates
pop_5y_unvac <- cdm$outcome_table |> filter(vac_status == 0) |>
  mutate(time_to_event = case_when(
    time_to_event > 365*5 ~ 365*5,
    time_to_event <= 365*5 ~ time_to_event
  ))

den_5y_unvac <- pop_5y_unvac |> pull(time_to_event) |> sum()/365

num_5y_unvac <- pop_5y_unvac |>
  filter(censor_reason == "outcome" & time_to_event < 5*365) |>
  tally() |> pull()
  
# IR 5 years
IR_5y_vac <- num_5y_vac*10^5/(den_5y_vac)
IR_5y_unvac <- num_5y_vac*10^5/(den_5y_vac)
  
info(logger, "COMPUTE IRR (Poisson Regression)")
# IRR 5 years
poisson_model_5y <- glm(formula = outcome_flag ~ vac_status, family = "poisson", data = cdm$outcome_table |>
                          #filter(censor_reason != "outcome" & time_to_event > 5*365) |>
                          mutate(outcome_flag = case_when(
                            outcome_date < cohort_start_date + years(5) ~ 1,
                            outcome_date > cohort_start_date + years(5) ~ 0
                          ))
                        )
poisson_coef_5y <- poisson_model_5y$coefficients[poisson_model_5y$coefficients != "(Intercept)"]
# not working because all values of outcome == 0

IRR <- exp(poisson_coef_5y)

# Survival Study
info(logger, "SURVIVAL STUDY")
info(logger, "Kaplan Meier plots")

survival_data <- survfit(Surv(time_to_event, outcome_flag) ~ vac_status, data = cdm$outcome_table) 

plot(survival_data,
     main = "Kaplan Meier Plot",
     xlab = "Time to event (days)", ylab = "Survival probability",
     #xlim = c(17, 37), ylim = c(0.75, 1),
     conf.int = F, mark.time = F,
     lty = 1:2) # 2 groups
# Add c(1,3,4,2) in two places to re-order the legend to match
# the ordering of the lines at 37 weeks
legend(20, 0.50,
       # Inside c(), the ordering is the same as the order shown in surv.ex7.4
       c("Vac", "Unvac")[c(1,2)],
       lty = c(1,2),
       title = "Strata: ")

info(logger, "Cox Regression")
cox_model <- coxph(Surv(time_to_event, outcome_flag) ~ vac_status, data = cdm$outcome_table)
summary_coxReg <- summary(cox_model)
coef_coxReg <- summary_coxReg$coef
# HRs, 95% CIs, p-values
cbind("HR" = exp(coef_coxReg[, "coef"]),
      exp(confint(cox_model)),
      "p-value" = coef_coxReg[, "Pr(>|z|)"])
