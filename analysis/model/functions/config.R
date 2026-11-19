# Study and nation settings for the 3-year analysis of the final SAP (v0.3, 30 October 2026, section 6.8).
# Change analysis settings here.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

STUDY <- list(
  cutoff          = as.Date("1958-09-01"),  # eligible at 65 (z = 1) if born on or after this date
  nip_start       = as.Date("2023-09-01"),  # index date for ineligible adults
  age_date        = as.Date("2023-08-31"),  # date for the cov_int_age check and the baseline table's age
  study_end       = as.Date("2026-08-31"),  # administrative study end
  induction_days  = 30L,                    # dementia on days 0 to 30 excludes the adult, HZ only from the HZ analysis
  bandwidths      = c("2" = 2, "1.5" = 1.5, "1" = 1),  # analysis bandwidths: years on either side of the cutoff
  primary_bandwidth = "2",
  outcomes        = c(dementia = "out_date_dementia", hz = "out_date_hz"),
  primary_outcome = "dementia",
  strata          = c("All", "Female", "Male"),
  incidence_years = c(1, 2, 3),             # model-based incidence this many years after the index date
  n_boot          = 1000L,                  # bootstrap samples for the effect of receipt
  max_boot_failed = 0.05,                   # CIs and SEs left missing above this share of failed samples
  seed            = 20261005L,
  min_dose_gap    = 28L,                    # days: a second RZV dose recorded sooner after the first is not counted
  min_flag_share  = 0.001,                  # comorbidity flags rarer than this define no subgroup (subgroup_flags())
  # Minimum events on each side of the cutoff for a national Weibull fit (bootstrap refits need 1), and minimum
  # adults with and without receipt on each side within delta's bandwidth for releasing the fitted uptakes. It
  # is one more than the disclosure threshold (SDC_DEFAULT), so that no released estimate rests on 7 or fewer.
  min_events_side = 8L
)

# The estimates of delta, the uptake discontinuity, each a regression of receipt by the end of follow-up on
# the running variable, separate each side of the cutoff, with a triangular kernel (uptake_fit() in fuzzy_rd.R).
#   main       SAP 6.5.1: at least one dose, local linear within rdrobust's MSE-optimal bandwidth, capped at
#              the analysis bandwidth
#   quadratic  SAP 6.5.1 and 6.6 (i): at least one dose, local quadratic over the whole analysis bandwidth
#   two_doses  SAP 6.6 (iii): two doses (derive.R), local linear as for main
DELTAS <- list(
  main      = list(receipt = "vacc_by_end",      p = 1L, bandwidth = "mserd"),
  quadratic = list(receipt = "vacc_by_end",      p = 2L, bandwidth = "analysis"),
  two_doses = list(receipt = "two_doses_by_end", p = 1L, bandwidth = "mserd")
)

# Primary disclosure control: suppress counts of 7 or fewer and statistics describing 7 or fewer people, and
# round larger counts to the nearest 5. Secondary disclosure, from combining released figures, is checked in
# each SDE's output review.
SDC_DEFAULT <- list(redact_at_or_below = 7L, round_to = 5L)

# data_end is the last date with complete linked data in each SDE. Follow-up ends by the earlier of data_end
# and study_end. Confirm data_end before the real analysis. run_nation.R --data-end sets it for one run.
NATIONS <- list(
  England = list(tre = "OpenSAFELY-TPP", bin = "month", data_end = as.Date("2026-08-31"), sdc = SDC_DEFAULT),
  Wales   = list(tre = "SAIL Databank",  bin = "week",  data_end = as.Date("2026-08-31"), sdc = SDC_DEFAULT)
)

# Baseline variables for the 3-year analysis (SAP Table 6).
# type: "cat" categorical, "bin" yes/no, "num" continuous (mean, SD, median and quartiles).
TABLE1_VARS <- list(
  # Age in fractional years, on age_date and at the index date (age_years() in derive.R)
  list(var = "age_aug2023",              type = "num",   label = "Age on 31 August 2023 (years)"),
  list(var = "age_index",                type = "num",   label = "Age at the index date (years)"),
  list(var = "cov_cat_sex",              type = "cat",   label = "Sex"),
  list(var = "cov_cat_ethnicity",        type = "cat",   label = "Ethnicity"),
  list(var = "cov_cat_imd",              type = "cat",   label = "Deprivation quintile (1 most deprived)"),
  list(var = "cov_cat_urban",            type = "cat",   label = "Urban or rural"),
  list(var = "cov_cat_region",           type = "cat",   label = "Region", nations = "England"),
  # Healthcare-use rates, annualised from the extract's counts over up to 2 years before the index date
  # (prior_rates() in derive.R)
  list(var = "pre_rate_hospitalisation", type = "num",   label = "Hospital admissions per year before the index date"),
  list(var = "pre_rate_gpvisit",         type = "num",   label = "Primary care visits per year before the index date")
)

# Optional comorbidity flags (SAP Table 6), each true if the adult has a record of the condition before the
# index date. The baseline table shows each flag the extract holds, and how many of them each adult has, and
# that number defines the subgroups of SAP 6.6 (ii) (comorbidity_subgroup() in describe.R). The log lists the
# flags the extract lacks.
COMORBIDITIES <- c(
  pre_bin_db           = "Diabetes",
  pre_bin_chol         = "High cholesterol",
  pre_bin_hypertension = "Hypertension",
  pre_bin_cvd          = "Cardiovascular disease (ischaemic heart disease, myocardial infarction, heart failure)",
  pre_bin_stroke       = "Stroke",
  pre_bin_depress      = "Depression",
  pre_bin_parkinson    = "Parkinson's disease",
  pre_bin_huntington   = "Huntington's disease",
  pre_bin_ms           = "Multiple sclerosis",
  pre_bin_mci          = "Mild cognitive impairment"
)

# Subgroups by the number of these comorbidities, for SAP 6.6 (ii) (comorbidity_subgroup() in describe.R).
COMORBIDITY_SUBGROUPS <- c("0 comorbidities", "1 comorbidity", "2 or more comorbidities")
