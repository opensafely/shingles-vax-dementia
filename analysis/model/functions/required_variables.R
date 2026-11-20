# Required extract variables for the SAP's 3-year analysis, their types and permitted categories, and the
# functions that read and check the extract. Names and types follow the dataset definition spreadsheet
# (OS_dataset_definition.xlsx). The trial run checks names, types and categories against the README.
#
# Extract conventions:
# - var_date_birth: the 15th of the recorded month (England) or the Thursday of the recorded week (Wales).
#   Month-start or Monday dates are converted. Any other date stops the run.
# - out_date_*: the first record from day 31 after the index date, when follow-up starts (SAP Table 3).
# - var_date_1vax_rzv: the first RZV dose on or after the earliest possible 65th birthday (the first day of
#   the recorded month or week, 65 years on), or 1 September 2023 for ineligible adults.
#   pre_bin_vax_shingles flags RZV or ZVL receipt before that date.
# - var_date_2vax_rzv: the first RZV dose at least 28 days (STUDY$min_dose_gap) after var_date_1vax_rzv, as
#   an earlier record is likely to duplicate the first dose.
# - pre_rate_hospitalisation and pre_rate_gpvisit: counts of admissions and visits, not rates, over the
#   2 years before the index date, truncated at var_date_lookback (the start of registration or of data
#   coverage). prior_rates() in derive.R annualises them.
# - pre_bin_dementia: dementia diagnosis or medication up to day 30 after the index date. Other pre_bin_* flags
#   refer to the period before the index date, and var_bin_* flags to the index date.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

REQUIRED_VARIABLES <- utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
name                      type         category   nations
var_date_birth            Date         Core       both
cov_int_age               Integer      Covariate  both
cov_cat_sex               Categorical  Inclusion  both
cov_cat_ethnicity         Categorical  Covariate  both
cov_cat_imd               Categorical  Covariate  both
cov_cat_urban             Categorical  Inclusion  both
cov_cat_region            Categorical  Inclusion  England
var_date_index            Date         Core       both
var_date_death            Date         Core       both
var_date_deregistered     Date         Core       both
var_date_residencyend     Date         Core       both
pre_bin_dementia          Binary       Exclusion  both
out_date_dementia         Date         Outcome    both
var_bin_residency         Binary       Inclusion  both
var_bin_registration      Binary       Inclusion  both
pre_bin_ic                Binary       Exclusion  both
pre_bin_vax_shingles      Binary       Exclusion  both
var_date_lookback         Date         Core       both
pre_rate_hospitalisation  Integer      Covariate  both
pre_rate_gpvisit          Integer      Covariate  both
out_date_hz               Date         Outcome    both
var_date_1vax_rzv         Date         Core       both
var_date_2vax_rzv         Date         Core       both
")

# Viral infections are reserved for the full study (SAP Table 6). Unused extract columns are discarded.
NOT_USED <- c("pre_bin_infection", "pre_int_com")   # pre_int_com: comorbidity_count() in describe.R counts the flags
ID_COLUMNS <- c("patient_id", "alf_pe", "ALF_PE")

# Permitted categorical values (SAP Table 6). Matching ignores case and surrounding spaces, and missing values
# are allowed. Other values cause an error.
CATEGORIES <- list(
  cov_cat_sex       = c("Female", "Male", "Intersex", "Unknown"),
  # Map the 16 or 18 census categories to these groups as the 2021 census does: Chinese under Asian, Gypsy or
  # Irish Traveller under White, and Arab under Other. Code missing or unstated ethnicity as "Unknown".
  cov_cat_ethnicity = c("White", "Asian or Asian British", "Black or Black British",
                        "Mixed", "Chinese or Other Ethnic Groups", "Unknown"),
  # Quintile 1 is most deprived (IMD in England, WIMD 2025 in Wales).
  cov_cat_imd       = c("1", "2", "3", "4", "5", "Unknown"),
  cov_cat_urban     = c("Urban", "Rural", "Unknown"),               # the first word of the ONS rural-urban classification
  # SAP labels. ehrQL names two of them differently: recode East to East of England, and
  # Yorkshire and The Humber to Yorkshire/Humberside.
  cov_cat_region    = c("North East", "North West", "Yorkshire and The Humber", "East Midlands", "West Midlands",
                       "East", "London", "South East", "South West", "Unknown")
)

# Optional: the comorbidity flags (COMORBIDITIES in functions/config.R) and the first dementia medication
# date on or after the index date (out_date_dementia_med). The latter contributes to the induction-period
# exclusion in cohort_criteria() (SAP 5.2.1.2).

read_extract <- function(path) {
  p <- tolower(path)
  d <- if (grepl("\\.csv(\\.gz)?$", p)) {
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
    } else {
      utils::read.csv(path, colClasses = "character", na.strings = c("", "NA"))
    }
  } else if (grepl("\\.(arrow|feather)$", p)) {
    arrow::read_feather(path)
  } else if (grepl("\\.parquet$", p)) {
    arrow::read_parquet(path)
  } else if (grepl("\\.rds$", p)) {
    readRDS(path)
  } else {
    stop("Unsupported extract format: ", path)
  }
  as.data.frame(d, stringsAsFactors = FALSE)
}

# Require this nation's variables, check and convert their types, and retain the permitted optional and
# identifier columns. Add missing other-nation variables as NA (region for Wales).
prepare_extract <- function(d, nation) {
  mine <- REQUIRED_VARIABLES$nations %in% c("both", nation)
  absent <- setdiff(REQUIRED_VARIABLES$name[mine], names(d))
  if (length(absent)) stop("The extract lacks: ", paste(absent, collapse = ", "))
  for (i in seq_len(nrow(REQUIRED_VARIABLES))) {
    v <- REQUIRED_VARIABLES$name[i]
    x <- if (v %in% names(d)) d[[v]] else rep(NA, nrow(d))
    d[[v]] <- switch(REQUIRED_VARIABLES$type[i],
                     Date = as_date(x, v), Integer = as_int(x, v),
                     Binary = as_flag(x, v), Categorical = as_cat(x, v))
  }
  extra <- intersect(names(COMORBIDITIES), names(d))               # optional comorbidity flags
  for (v in extra) d[[v]] <- as_flag(d[[v]], v)
  if ("out_date_dementia_med" %in% names(d)) {                     # optional dementia medication date
    d$out_date_dementia_med <- as_date(d$out_date_dementia_med, "out_date_dementia_med")
    extra <- c(extra, "out_date_dementia_med")
  }
  d[names(d) %in% c(REQUIRED_VARIABLES$name, extra, ID_COLUMNS)]
}

blank_to_na <- function(x) {
  v <- trimws(as.character(x))
  v[v %in% c("", "NA")] <- NA
  v
}

# Stop on values that cannot be read as YYYY-MM-DD dates. Read a date-time as its calendar date in its own time
# zone, so that midnight in British Summer Time stays on its day. No type check reports counts, because error
# messages can reach job logs that are not disclosure-checked.
as_date <- function(x, name) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(format(x, "%Y-%m-%d"), format = "%Y-%m-%d"))
  v <- blank_to_na(x)
  out <- as.Date(v, format = "%Y-%m-%d")
  bad <- sum(!is.na(v) & is.na(out))
  if (bad) stop(name, ": some values are not dates in YYYY-MM-DD form")
  out
}

as_int <- function(x, name) {
  if (is.numeric(x)) {
    if (any(x != round(x), na.rm = TRUE)) stop(name, ": non-integer values")
    return(as.integer(x))
  }
  v <- blank_to_na(x)
  out <- suppressWarnings(as.integer(v))
  bad <- sum(!is.na(v) & is.na(out))
  if (bad) stop(name, ": some values are not integers")
  out
}

# Accept the ehrQL T/F and SAIL 1/0 Boolean encodings, and also true/false, y/n and yes/no, ignoring case.
as_flag <- function(x, name) {
  if (is.logical(x)) return(x)
  v <- tolower(blank_to_na(x))
  out <- rep(NA, length(v))
  out[v %in% c("1", "t", "true", "y", "yes")] <- TRUE
  out[v %in% c("0", "f", "false", "n", "no")] <- FALSE
  bad <- sum(!is.na(v) & is.na(out))
  if (bad) stop(name, ": some values are not yes/no")
  out
}

# Return each value spelt as in CATEGORIES. On failure, report at most five invalid values, without counts.
as_cat <- function(x, name) {
  v <- blank_to_na(x)
  allowed <- CATEGORIES[[name]]
  if (is.null(allowed)) return(v)
  out <- allowed[match(tolower(v), tolower(allowed))]
  bad <- unique(v[!is.na(v) & is.na(out)])
  if (length(bad)) {
    stop(name, ": values other than ", paste0("\"", allowed, "\"", collapse = ", "), " or empty, such as ",
         paste0("\"", utils::head(bad, 5), "\"", collapse = ", "))
  }
  out
}
