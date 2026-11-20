#!/usr/bin/env Rscript
# Dementia incidence for the sample size calculation, run in each SDE before the main analysis (run_nation.R),
# on the same extract. Apply the main analysis's cohort selection (SAP 5.2.1) and dementia follow-up (SAP
# Table 3). Report the adults in the cohort on each side of the cutoff and, for ineligible adults only, the
# dementia events, person-years and crude rates per 1000 person-years by age. Ineligible adults are not offered
# RZV, so their rates estimate dementia incidence without the programme's vaccination at ages 65 to 69. No
# eligible adult's outcome is used, so the comparison of the two sides stays unseen until the main analysis.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)
#
# Usage:   Rscript sample_size_incidence.R --nation England|Wales --data FILE --out DIR [--data-end DATE]
# Options: --data-end DATE   Last date of complete linked data, as YYYY-MM-DD (default: the nation's data_end in
#                            functions/config.R).
#
# Input:   --data, the extract for run_nation.R (README, "What the extract must contain").
# Needs:   functions/*.R, and the R package arrow for a .arrow or .parquet extract.
#
# Outputs, in the --out folder (see README, "Dementia incidence for the sample size"):
#   sample_size_cohort.csv     adults in the analysis cohort, eligible and ineligible, at each bandwidth
#   sample_size_incidence.csv  ineligible adults: dementia events, person-years and rates by age
#   log.txt, session_info.txt

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
this_file <- gsub("~+~", " ", this_file, fixed = TRUE)          # Rscript encodes spaces in its own path
here <- if (length(this_file)) dirname(normalizePath(this_file[1])) else getwd()
for (f in sort(list.files(file.path(here, "functions"), pattern = "\\.R$", full.names = TRUE))) source(f)

opt <- parse_args(commandArgs(TRUE), list(nation = NA, data = NA, out = NA, "data-end" = NA))
nation <- opt$nation
if (!isTRUE(nation %in% names(NATIONS)) || is.na(opt$data) || is.na(opt$out)) {
  stop("Usage: Rscript sample_size_incidence.R --nation England|Wales --data FILE --out DIR [--data-end DATE]")
}
cfg <- NATIONS[[nation]]
if (!file.exists(opt$data[1])) stop("--data: no such file: ", opt$data[1])
if (!is.na(opt[["data-end"]][1])) {
  cfg$data_end <- as.Date(opt[["data-end"]][1], format = "%Y-%m-%d")
  if (is.na(cfg$data_end)) stop("--data-end must be a date as YYYY-MM-DD")
}
rule <- cfg$sdc
out <- opt$out
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(out, "log.txt"))
log_line("EPI-ZOSTER-110 dementia incidence for the sample size: ", nation, " (", cfg$tre, ")")
log_line("Follow-up to ", format(min(STUDY$study_end, cfg$data_end)))
log_line("Counts in this log are rounded to ", rule$round_to, ", and withheld at ", rule$redact_at_or_below, " or fewer")

# ---- Cohort selection (SAP 5.2.1), as in run_nation.R ----
d <- prepare_extract(read_extract(opt$data), nation)
d <- derive_birth(d, cfg$bin, STUDY)
fu <- follow_up(d, STUDY$outcomes[[STUDY$primary_outcome]], STUDY, cfg$data_end)
ff <- first_failed(cohort_criteria(d, nation, STUDY, fu))
co <- d[ff == 0L, , drop = FALSE]
fu <- fu[ff == 0L, , drop = FALSE]
log_line("Analysis cohort: ", sdc_count(nrow(co), rule), " adults")

# ---- Adults on each side of the cutoff ----
size <- do.call(rbind, lapply(names(STUDY$bandwidths), function(w) {
  inside <- in_bandwidth(co, STUDY, STUDY$bandwidths[[w]])
  data.frame(nation = nation, bandwidth = w, group = c("Eligible", "Ineligible"),
             n = sdc_count(c(sum(inside & co$z == 1L), sum(inside & co$z == 0L)), rule))
}))
write_table(size, out, "sample_size_cohort.csv")

# ---- Ineligible adults' dementia incidence by age (2-year bandwidth, the whole cohort) ----
ineligible <- co$z == 0L
write_table(incidence_by_age(co[ineligible, , drop = FALSE], fu[ineligible, , drop = FALSE], nation, STUDY, rule),
            out, "sample_size_incidence.csv")
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "session_info.txt"))
log_line("Done")
