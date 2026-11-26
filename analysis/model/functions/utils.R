# Shared utilities for disclosure control, dates, model fitting, logging, argument parsing and output.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

.log <- new.env()

log_open <- function(path) {
  .log$path <- path
  cat("", file = path)
}

log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S  "), ...)
  cat(msg, "\n", sep = "")
  if (!is.null(.log$path)) cat(msg, "\n", sep = "", file = .log$path, append = TRUE)
}

# ---- Disclosure control ----

# Suppress counts at or below the threshold. Round larger counts to the nearest specified base.
sdc_count <- function(x, rule) {
  out <- rule$round_to * round(x / rule$round_to)
  out[is.na(x) | x <= rule$redact_at_or_below] <- NA
  out
}

# Suppress summary statistics when the number of people they describe is at or below the threshold.
sdc_stat <- function(x, n, rule) ifelse(!is.na(n) & n > rule$redact_at_or_below, x, NA)

# ---- Dates ----

REF_MONDAY <- as.Date("1958-09-01")   # a Monday, and the eligibility cutoff

first_of_month <- function(d) as.Date(format(d, "%Y-%m-01"))
days_in_month <- function(month_start) as.integer(first_of_month(month_start + 31L) - month_start)
monday_on_or_before <- function(d) d - as.integer(d - REF_MONDAY) %% 7L

add_years <- function(d, n) {          # 29 February moves to 1 March in a common year
  lt <- as.POSIXlt(d)
  lt$year <- lt$year + n
  as.Date(lt)
}

add_months <- function(d, n) if (n == 0) d else seq(d, by = paste(n, "months"), length.out = 2)[2]

age_on <- function(dob, on) {
  as.integer(format(on, "%Y")) - as.integer(format(dob, "%Y")) - (format(on, "%m%d") < format(dob, "%m%d"))
}

# ---- Fitting and output ----

# Log fitting errors and return NULL. Log the distinct warnings and keep them in the "warnings" attribute.
fit_safely <- function(label, expr) {
  warns <- character(0)
  out <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      log_line("  ", label, ": not fitted (", conditionMessage(e), ")")
      NULL
    }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(warns)) log_line("  ", label, ": warning: ", paste(unique(warns), collapse = "; "))
  if (!is.null(out)) attr(out, "warnings") <- unique(warns)
  out
}

write_table <- function(df, dir, name) {
  utils::write.csv(df, file.path(dir, name), row.names = FALSE)
  log_line("  wrote ", name, " (", nrow(df), " rows)")
  invisible(df)
}

# Compute a combined MD5 checksum of the analysis code files from each file's name and MD5, ignoring folders.
# Drop carriage returns before hashing and write the list byte for byte, so that Windows line endings, in the
# files or from writing on Windows, give the same checksum.
code_checksum <- function(files) {
  f <- tempfile()
  on.exit(unlink(f))
  md5 <- vapply(files, function(x) {
    b <- readBin(x, "raw", file.size(x))
    writeBin(b[b != as.raw(13)], f)
    unname(tools::md5sum(f))
  }, "")
  writeBin(charToRaw(paste0(basename(files), " ", md5, "\n", collapse = "")), f)
  unname(tools::md5sum(f))
}

# ---- Argument parsing ----

# Parse --key value arguments (one or more values per key). Reject unknown options and options without a value.
parse_args <- function(args, defaults) {
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% names(defaults)) stop("Unknown argument: ", args[i])
    vals <- character(0)
    i <- i + 1
    while (i <= length(args) && !startsWith(args[i], "--")) {
      vals <- c(vals, args[i])
      i <- i + 1
    }
    if (!length(vals)) stop("--", key, " needs a value")
    out[[key]] <- vals
  }
  out
}
