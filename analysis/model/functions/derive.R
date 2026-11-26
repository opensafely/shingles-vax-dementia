# Derived variables: assumed birth dates, the RD running variable, eligibility, index dates, age, prior
# healthcare-use rates, outcome-specific follow-up and RZV receipt (SAP 5.2.1, 6.4.1, 6.5.1, 6.6, Tables 3
# and 6).
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

# Derive birth-date bins (bin_start, bin_width), assumed birth dates (dob), eligibility (z), index dates, the
# running variable (r_days, the SAP's x_i^c to within a day), the birth month for tables and figures
# (month_start) and age in completed years on 31 August 2023 (age_derived). England uses months. Wales uses
# Monday-starting weeks.
derive_birth <- function(d, bin, study) {
  dob <- d$var_date_birth
  start <- if (bin == "month") first_of_month(dob) else monday_on_or_before(dob)
  width <- if (bin == "month") days_in_month(start) else rep(7L, length(dob))
  middle <- start + if (bin == "month") 14L else 3L      # the 15th (SAP 5.2.1.1's example) or the Thursday
  # Each adult sits at the middle of the recorded month or week. A date at the start of its bin (TPP's
  # first of the month, SAIL's Monday) is moved there, and any other date stops the run.
  seen <- !is.na(dob)
  if (any(seen & dob != start & dob != middle)) {
    stop("var_date_birth holds dates that are neither the first day nor the middle (", if (bin == "month") "the 15th" else "the Thursday",
         ") of their ", bin, ": see the README, 'What the extract must contain'")
  }
  d$bin_start <- start
  d$bin_width <- width
  d$dob <- middle
  # Assign Welsh birth weeks to the month containing their Thursday for tables and figures by birth month, as
  # weekly counts would be too small to release. As the cutoff is a Monday and the first of a month, no month
  # mixes eligible and ineligible adults. Models and the density test retain the recorded weekly resolution.
  d$month_start <- first_of_month(d$dob)
  d$r_days <- as.integer(d$dob - study$cutoff)            # positive for eligible adults
  d$z <- as.integer(start >= study$cutoff)                # eligibility from the recorded bin
  if (any(seen & d$z != (d$dob >= study$cutoff))) stop("A recorded birth month or week straddles the cutoff")
  d$age_derived <- age_on(d$dob, study$age_date)
  index <- add_years(d$dob, 65L)                          # eligible: the assumed 65th birthday
  index[d$z %in% 0L] <- study$nip_start                   # ineligible: 1 September 2023
  d$index_date <- index
  attr(d, "dob_moved_to_middle") <- sum(seen & dob == start)    # dates that were bin starts
  d
}

# Include assumed birth dates in [cutoff - years, cutoff + years), with `years` rounded to whole calendar months.
in_bandwidth <- function(d, study, years) {
  m <- round(12 * years)
  !is.na(d$dob) & d$dob >= add_months(study$cutoff, -m) & d$dob < add_months(study$cutoff, m)
}

# Return age in completed years plus the elapsed fraction of the current birthday year.
# Unlike (date - dob) / 365.25, this gives exactly 65 on the 65th birthday.
age_years <- function(dob, on) {
  y <- age_on(dob, on)
  last <- add_years(dob, y)
  y + as.numeric(on - last) / as.numeric(add_years(dob, y + 1L) - last)
}

# Outcome-specific follow-up (SAP Table 3) starts after the 30-day induction period. Count time from the end
# of day 30 after the index date, so an outcome on day 31 has time 1 day, and an adult whose follow-up ends by
# day 30 is not at risk. Censor at death, deregistration, end of residency (var_date_residencyend) or the
# administrative end, the earlier of the study end and data_end (the end of complete data). An outcome on the
# censoring date counts as an event. Dementia does not censor HZ follow-up. Outcome dates start on day 31
# (README, convention 2). A date on or before day 30, which data_checks.csv counts, gives time 0 or less, so
# the adult is not at risk for that outcome.
follow_up <- function(d, outcome, study, data_end) {
  entry <- d$index_date + study$induction_days
  admin <- min(study$study_end, data_end)
  # A missing date means no death, deregistration or end of residency. The administrative end is always set.
  censor <- pmin(d$var_date_death, d$var_date_deregistered, d$var_date_residencyend, admin, na.rm = TRUE)
  ev <- d[[outcome]]
  event <- !is.na(ev) & ev <= censor
  end <- censor
  end[event] <- ev[event]
  time_days <- as.integer(end - entry)
  # Receipt for the fuzzy RD (SAP 6.5.1) is at least one recorded RZV dose on or before the end of follow-up
  # for this outcome (vacc_by_end). This reads 6.5.1's "relevant follow-up endpoint" as each adult's own end of
  # follow-up, not the administrative end for everyone. The baseline table counts doses by the administrative
  # end instead (vacc_by_date, two_doses_by_date). A dose before the assumed index date counts, as the true 65th
  # birthday can fall before the assumed one. A recorded second dose establishes receipt from its own date if the
  # first-dose date is missing (SAP 6.1 notes missing first doses). data_checks.csv counts both.
  dose <- pmin(d$var_date_1vax_rzv, d$var_date_2vax_rzv, na.rm = TRUE)
  # Two doses (SAP 6.6 (iii), and the baseline table) need two recorded doses at least min_dose_gap (28) days
  # apart, or a second dose with no first recorded. UKHSA recommends at least 8 weeks between the doses but
  # counts a second dose given at least 4 weeks after the first as valid, and advises repeating one given
  # sooner. A record less than 28 days after the first is therefore not counted as a second dose (SAP 6.1
  # notes that some records show two doses a short time apart). Doses recorded out of order are taken in
  # date order, as receipt above takes the earlier date. data_checks.csv counts both.
  dose2 <- pmax(d$var_date_1vax_rzv, d$var_date_2vax_rzv)        # the later date, when both are recorded
  dose2[(dose2 - dose < study$min_dose_gap) %in% TRUE] <- NA
  only2 <- is.na(d$var_date_1vax_rzv) & !is.na(d$var_date_2vax_rzv)
  dose2[only2] <- d$var_date_2vax_rzv[only2]
  data.frame(time_days = time_days,
             event = as.integer(event),
             at_risk = !is.na(time_days) & time_days > 0,
             vacc_by_end = as.integer(!is.na(dose) & dose <= end),
             two_doses_by_end = as.integer(!is.na(dose2) & dose2 <= end),
             vacc_by_date = as.integer(!is.na(dose) & dose <= admin),
             two_doses_by_date = as.integer(!is.na(dose2) & dose2 <= admin))
}

# Annualise prior hospital admissions and primary care visits (SAP Table 6) over the available look-back
# period, capped at two years, the same period the extract counts them over. Return missing rates when the
# look-back period is unavailable or not positive.
prior_rates <- function(d) {
  years <- pmin(2, as.numeric(d$index_date - d$var_date_lookback) / 365.25)
  years[is.na(years) | years <= 0] <- NA
  d$pre_rate_hospitalisation <- d$pre_rate_hospitalisation / years
  d$pre_rate_gpvisit <- d$pre_rate_gpvisit / years
  d
}
