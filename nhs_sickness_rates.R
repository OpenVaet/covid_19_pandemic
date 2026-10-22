# ============================================================================
# NHS Sickness Rate Analysis
# Breakpoints vs “observable norm” (2010–2018), seasonality-adjusted
# Run with RStudio's Source button or source("nhs_sickness_rates.R", encoding="UTF-8").
# Exports: england_sickness_results/run_<timestamp>_<id>/report.html plus
# source workbook, CSV tables, fitted models, plot, configuration and audits.
# Shared downloads are cached under england_sickness_results/cache/.
# This stops the run if input preparation fails; do not continue line by line.
# ============================================================================

# ---- Packages ----
# install.packages(c("dplyr", "tibble", "purrr", "ggplot2", "lubridate", "slider",
#                    "scales", "ggrepel", "strucchange", "readxl", "rvest", "curl", "jsonlite"))
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(lubridate)
library(slider)
library(scales)
library(ggrepel)
library(strucchange)

# ---- User-tunable parameters ----
# Each successful run is self-contained. The cache is shared across runs.
OUT_DIR <- getOption("nhs.output_dir", "england_sickness_results")
RUN_ANALYSIS <- getOption("nhs.run_analysis", TRUE)
MANUAL_CSV <- NULL         # Used only when DOWNLOAD_DATA is FALSE: Month, England
LOCAL_WORKBOOK <- getOption("nhs.local_workbook", NULL) # Optional offline .xlsx input

# Automatic input preparation (paths are relative to the working directory).
# TRUE: extract the national series from the downloaded/cached NHS workbook.
# FALSE: use MANUAL_CSV. No existing input file is overwritten.
DOWNLOAD_DATA <- TRUE
REFRESH_DOWNLOAD <- FALSE  # TRUE: fetch the workbook again, including revisions.
NHS_PUBLICATION_URL <- paste0(
  "https://digital.nhs.uk/data-and-information/publications/statistical/",
  "nhs-sickness-absence-rates/october-2025"
)
# NULL uses a verified direct link for October 2025, avoiding the publication
# webpage (which can return HTTP 403 to R). Other releases use page discovery.
# Alternatively, provide the desired workbook's direct .xlsx URL here.
NHS_WORKBOOK_URL <- NULL
# This publication's workbook contains April 2009 to October 2025 in Table 1.
# Its separate "Monthly data text file" CSV contains ONLY the release month.

BASELINE_START <- 2010
BASELINE_END   <- 2018
PROJ_START     <- 2019
PROJ_END       <- 2025

# Rolling window: (t-3, t-2, t-1, t, t+1, t+2)  => 6 points
ROLL_BEFORE <- 3
ROLL_AFTER  <- 2
ROLL_COMPLETE <- TRUE

# Breakpoints configuration
MAX_BREAKS <- 10
H_FRACTION <- 0.10    # minimum 10% of data per segment

# Breakpoint model choice:
# - "level": step changes in deviation (dev_100k_raw ~ 1)
# - "trend": changes in deviation trend (dev_100k_raw ~ month_num)
BP_MODEL_TYPE <- "level"

# ---- Helpers ----

safe_month <- function(x) {
  # Tries multiple month formats robustly; returns Date (first of month) or NA
  x <- as.character(x)
  
  # NHS month names are English, including on French-language computers.
  d <- suppressWarnings(my(x, locale = "C", quiet = TRUE))
  
  if (anyNA(d)) {
    d2 <- suppressWarnings(parse_date_time(
      x,
      orders = c(
        "my", "Ym", "Y-m", "Y/m",
        "bY", "Yb", "mY", "m/Y",
        "b-Y", "Y-b", "Ymd", "Y-m-d", "Y/m/d"
      ),
      tz = "UTC", locale = "C", quiet = TRUE
    ))
    d[is.na(d)] <- d2[is.na(d)]
  }
  
  # normalize to first day of month
  d <- floor_date(as.Date(d), unit = "month")
  as.Date(d)
}

get_optimal_breaks_bic <- function(bp_full_obj) {
  s <- summary(bp_full_obj)
  
  # BIC is typically provided as s$BIC
  if (!is.null(s$BIC)) {
    bic <- s$BIC
  } else if (!is.null(s$RSS) && is.matrix(s$RSS) && "BIC" %in% rownames(s$RSS)) {
    bic <- s$RSS["BIC", ]
  } else {
    stop("Could not extract BIC from breakpoint summary. Check strucchange version/output.")
  }
  
  idx <- which.min(bic)
  
  # Prefer named break counts if present
  if (!is.null(names(bic)) && nzchar(names(bic)[idx])) {
    return(as.integer(names(bic)[idx]))
  }
  
  # Fallback: assume entries correspond to 0,1,2,... breaks
  as.integer(idx - 1)
}

cap_max_breaks <- function(n, h, requested) {
  # maximum feasible breaks given minimum segment fraction h:
  # each segment must have at least floor(h*n) points; so max segments ≈ floor(1/h)
  # => max breaks ≈ floor(1/h) - 1
  max_by_h <- max(0, floor(1 / h) - 1)
  max_by_n <- max(0, n - 2)  # need at least 2 points to define any segmentation sensibly
  min(requested, max_by_h, max_by_n)
}

# ============================================================================
# Download and prepare the historical England series
# ============================================================================

download_nhs_file <- function(url, destination) {
  previous_timeout <- getOption("timeout")
  options(timeout = max(180, previous_timeout))
  on.exit(options(timeout = previous_timeout), add = TRUE)
  status <- tryCatch(
     {
      if (!requireNamespace("curl", quietly = TRUE)) stop("Install package curl")
      curl::curl_download(url, destination, quiet = FALSE,
        handle = curl::new_handle(followlocation = TRUE, timeout = 180))
      0L
    },
    error = function(e) stop("NHS download failed: ", conditionMessage(e),
                             "\nURL: ", url, call. = FALSE)
  )
  if (!identical(status, 0L) || !file.exists(destination) ||
      is.na(file.info(destination)$size) || file.info(destination)$size == 0) {
    stop("NHS download was unsuccessful or empty: ", url, call. = FALSE)
  }
}

find_nhs_workbook_url <- function(publication_url) {
  # Use the official file download independently of the HTML publication page.
  # Match the complete publication URL so changing releases cannot silently
  # select the October 2025 workbook.
  october_2025_page <- paste0(
    "https://digital.nhs.uk/data-and-information/publications/statistical/",
    "nhs-sickness-absence-rates/october-2025"
  )
  if (identical(sub("/+$", "", publication_url), october_2025_page)) {
    return(paste0(
      "https://files.digital.nhs.uk/83/6727AA/",
      "NHS%20Sickness%20Absence%20rates%2C%20October%202025.xlsx"
    ))
  }
  if (!requireNamespace("rvest", quietly = TRUE)) {
    stop('Install the download dependency: install.packages("rvest")', call. = FALSE)
  }
  page_path <- tempfile(fileext = ".html")
  on.exit(unlink(page_path), add = TRUE)
  download_nhs_file(publication_url, page_path)
  links <- rvest::html_elements(rvest::read_html(page_path), "a[href]")
  labels <- rvest::html_text2(links)
  hrefs <- rvest::html_attr(links, "href")
  # Exclude the distinct reason/staff-group workbooks on the same page.
  selected <- grepl("\\.xlsx($|[?#])", hrefs, ignore.case = TRUE) &
    grepl("NHS Sickness Absence Rates", labels, ignore.case = TRUE) &
    grepl("Monthly Tables", labels, ignore.case = TRUE) &
    !grepl("reason|benchmark", labels, ignore.case = TRUE)
  urls <- unique(xml2::url_absolute(hrefs[which(selected)], publication_url))
  if (length(urls) != 1L) {
    stop("Could not identify exactly one NHS monthly-rates workbook (found ",
         length(urls), "). Open NHS_PUBLICATION_URL and set NHS_WORKBOOK_URL ",
         "to the direct Monthly Tables .xlsx link.", call. = FALSE)
  }
  urls[[1]]
}

read_nhs_england_series <- function(workbook_path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop('Install the workbook dependency: install.packages("readxl")', call. = FALSE)
  }
  sheets <- readxl::excel_sheets(workbook_path)
  sheet <- sheets[tolower(trimws(sheets)) == "table 1"]
  if (length(sheet) != 1L) stop("The NHS workbook has no unique Table 1.", call. = FALSE)

  # Locate the header by its labels, rather than assume a fixed number of notes.
  preview <- readxl::read_excel(workbook_path, sheet = sheet, col_names = FALSE,
                               col_types = "text", n_max = 20,
                               .name_repair = "minimal")
  header_rows <- which(vapply(seq_len(nrow(preview)), function(i) {
    cells <- tolower(trimws(as.character(unlist(preview[i, ], use.names = FALSE))))
    all(c("month", "england") %in% cells)
  }, logical(1)))
  if (length(header_rows) != 1L) {
    stop("Could not locate the Month / England header in Table 1.", call. = FALSE)
  }
  header <- tolower(trimws(as.character(unlist(preview[header_rows, ], use.names = FALSE))))
  values <- readxl::read_excel(workbook_path, sheet = sheet, skip = header_rows,
                              col_names = FALSE, col_types = "text",
                              .name_repair = "minimal")
  months_text <- trimws(values[[match("month", header)]])
  rates_text <- trimws(values[[match("england", header)]])
  keep <- !(is.na(months_text) | months_text == "") |
          !(is.na(rates_text) | rates_text == "")
  months <- safe_month(months_text[keep])
  rates <- suppressWarnings(as.numeric(sub("%$", "", rates_text[keep])))
  if (!length(months) || anyNA(months) || any(!is.finite(rates)) ||
      any(rates < 0 | rates > 100) || anyDuplicated(months)) {
    stop("Invalid, missing or duplicate monthly observations in NHS Table 1. ",
         "The existing analysis CSV has not been replaced.", call. = FALSE)
  }
  order_rows <- order(months)
  months <- months[order_rows]
  rates <- rates[order_rows]
  expected <- seq(min(months), max(months), by = "month")
  baseline <- seq(as.Date(sprintf("%d-01-01", BASELINE_START)),
                  as.Date(sprintf("%d-12-01", BASELINE_END)), by = "month")
  if (length(months) != length(expected) || !all(baseline %in% months)) {
    stop("The NHS national series has gaps or does not cover the full baseline.",
         call. = FALSE)
  }
  # Excel stores 5.67 for 5.67%; do not multiply these values by 100.
  data.frame(Month = format(months, "%Y-%m-%d"), England = rates,
             stringsAsFactors = FALSE)
}

# Main function isolates each run and leaves incomplete runs clearly marked.
run_england_analysis <- function() {
  required <- c("readxl", "curl", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  run_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = OUT_DIR)
  dir.create(file.path(run_dir, "source"), recursive = TRUE)
  dir.create(file.path(run_dir, "plots"))
  NHS_DOWNLOAD_DIR <- file.path(OUT_DIR, "cache")
  dir.create(NHS_DOWNLOAD_DIR, showWarnings = FALSE)
  DATA_PATH <- file.path(run_dir, "source", "overall_sickness_rates.csv")
  completed <- FALSE
  on.exit(if (!completed) writeLines("Run incomplete. Do not use for comparison.", file.path(run_dir, "FAILED.txt")), add = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
  dput(list(publication = NHS_PUBLICATION_URL, workbook_override = NHS_WORKBOOK_URL,
    local_workbook = LOCAL_WORKBOOK, manual_csv = MANUAL_CSV, download = DOWNLOAD_DATA,
    refresh = REFRESH_DOWNLOAD, baseline_start = BASELINE_START, baseline_end = BASELINE_END,
    projection_end = PROJ_END, h = H_FRACTION, maximum_breaks = MAX_BREAKS,
    model_type = BP_MODEL_TYPE, baseline_time_unit = "calendar day"), file.path(run_dir, "configuration.R"))

if (!is.null(LOCAL_WORKBOOK)) {
  workbook_path <- normalizePath(LOCAL_WORKBOOK, mustWork = TRUE)
  national_series <- read_nhs_england_series(workbook_path)
  write.csv(national_series, DATA_PATH, row.names = FALSE)
} else if (DOWNLOAD_DATA) {
  dir.create(NHS_DOWNLOAD_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(DATA_PATH), recursive = TRUE, showWarnings = FALSE)
  release <- sub(".*/", "", sub("/+$", "", NHS_PUBLICATION_URL))
  release <- gsub("[^A-Za-z0-9_-]", "_", release)
  workbook_path <- file.path(NHS_DOWNLOAD_DIR, paste0("nhs_sickness_absence_", release, ".xlsx"))

  # An explicit URL override always downloads, avoiding reuse of another source.
  if (REFRESH_DOWNLOAD || !file.exists(workbook_path) || !is.null(NHS_WORKBOOK_URL)) {
    workbook_url <- if (is.null(NHS_WORKBOOK_URL)) {
      find_nhs_workbook_url(NHS_PUBLICATION_URL)
    } else NHS_WORKBOOK_URL
    temporary_workbook <- tempfile(tmpdir = NHS_DOWNLOAD_DIR, fileext = ".xlsx")
    # Validate the complete download before replacing either cached input.
    tryCatch({
      download_nhs_file(workbook_url, temporary_workbook)
      national_series <- read_nhs_england_series(temporary_workbook)
      if (!file.copy(temporary_workbook, workbook_path, overwrite = TRUE)) {
        stop("Could not save the NHS workbook: ", workbook_path, call. = FALSE)
      }
      writeLines(c(paste0("Publication: ", NHS_PUBLICATION_URL),
                   paste0("Workbook: ", workbook_url),
                   paste0("Downloaded UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE))),
                 paste0(workbook_path, ".source.txt"))
    }, finally = unlink(temporary_workbook))
  } else {
    message("Using cached NHS workbook: ", workbook_path,
            " (set REFRESH_DOWNLOAD <- TRUE to download again).")
    national_series <- read_nhs_england_series(workbook_path)
  }
  temporary_csv <- tempfile(tmpdir = dirname(DATA_PATH), fileext = ".csv")
  tryCatch({
    write.csv(national_series, temporary_csv, row.names = FALSE, na = "")
    if (!file.copy(temporary_csv, DATA_PATH, overwrite = TRUE)) {
      stop("Could not save the national series: ", DATA_PATH, call. = FALSE)
    }
  }, finally = unlink(temporary_csv))
  message(sprintf("Prepared %d monthly observations (%s to %s): %s",
                  nrow(national_series), national_series$Month[1],
                  tail(national_series$Month, 1), DATA_PATH))
}
if (is.null(LOCAL_WORKBOOK) && !DOWNLOAD_DATA) {
  if (is.null(MANUAL_CSV) || !file.exists(MANUAL_CSV))
    stop("Set MANUAL_CSV to an existing CSV with Month and England columns.", call. = FALSE)
  if (!file.copy(MANUAL_CSV, DATA_PATH, overwrite = TRUE)) stop("Could not copy MANUAL_CSV")
}
if (exists("workbook_path", inherits = FALSE)) {
  if (!file.copy(workbook_path, file.path(run_dir, "source", basename(workbook_path))))
    stop("Could not archive source workbook")
  if (file.exists(paste0(workbook_path, ".source.txt")))
    file.copy(paste0(workbook_path, ".source.txt"), file.path(run_dir, "source", "download_provenance.txt"))
}
if (!file.exists(DATA_PATH)) {
  stop("Input CSV not found: ", DATA_PATH,
       ". Set DOWNLOAD_DATA <- TRUE to download and prepare it.", call. = FALSE)
}

# ============================================================================
# Read data
# ============================================================================
raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE) |>
  as_tibble()

df <- raw |>
  rename(rate_pct = England) |>
  mutate(
    month    = safe_month(Month),
    rate_pct = as.numeric(rate_pct)
  ) |>
  arrange(month) |>
  mutate(
    rate_100k = rate_pct * 1000,                 # 1% = 1000 per 100k
    month_num = as.numeric(month),
    month_of_year = factor(month(month), levels = 1:12, labels = month.abb)
  )

# Data audit precedes every fit; missing values are never silently dropped.
expected_baseline <- seq(as.Date(sprintf("%d-01-01", BASELINE_START)),
  as.Date(sprintf("%d-12-01", BASELINE_END)), by = "month")
valid_dates <- df$month[!is.na(df$month)]
continuous <- length(valid_dates) > 0L &&
  identical(as.character(valid_dates), as.character(seq(min(valid_dates), max(valid_dates), by = "month")))
quality <- data.frame(check = c("Dates parse", "Unique months", "Rates finite and in 0-100",
  "No internal missing months", "Complete baseline"),
  passed = c(!anyNA(df$month), !anyDuplicated(df$month), all(is.finite(df$rate_pct) & df$rate_pct >= 0 & df$rate_pct <= 100),
    continuous, all(expected_baseline %in% df$month)))
write.csv(quality, file.path(run_dir, "quality_checks.csv"), row.names = FALSE)
if (!all(quality$passed)) stop("Input validation failed; see quality_checks.csv in ", run_dir, call. = FALSE)
# Keep the original all-available-months England analysis, including pre-2010 data.
# The comparison script separately applies a common 2010-2025 analysis window.
notes <- c(
  "England national monthly sickness absence rates, NHS Table 1. Values are percentages of available days lost, not percentages of people sick.",
  "Data coverage is the published workbook span; unavailable tail months are not imputed.",
  "Legacy baseline retained: rate_100k ~ elapsed calendar days + month-of-year, fitted in 2010-2018.",
  "Legacy breakpoints retain the full observed span and h=0.10; last-before and first-after dates are both exported.",
  "The comparison script re-fits monthly-index baselines and independent level breaks on the shared 2010-2025 window.",
  "Definition/recording changes cannot be inferred from rate values. Consult the source workbook notes and official quality documentation before interpreting breaks.",
  "Revision audit compares rate values with the preceding successful local snapshot only; the first run has no historical vintage for comparison.",
  "Original OLS mean/prediction intervals assume independent errors. Smoothed pointwise bounds are not a formal interval for a six-month average.",
  "Event labels in the original figure are annotations, not estimates of causal effects.")
writeLines(notes, file.path(run_dir, "READ_BEFORE_ANALYSIS.txt"))
old_path <- file.path(NHS_DOWNLOAD_DIR, "last_successful_series.csv")
old <- if (file.exists(old_path)) read.csv(old_path, stringsAsFactors = FALSE) else NULL
current <- data.frame(Month = format(df$month, "%Y-%m-%d"), England = df$rate_pct)
revisions <- data.frame(Month = character(), previous_rate_pct = numeric(), current_rate_pct = numeric(), change = character())
if (!is.null(old)) {
  revisions <- merge(old, current, by = "Month", all = TRUE, suffixes = c("_previous", "_current"))
  names(revisions)[2:3] <- c("previous_rate_pct", "current_rate_pct")
  revisions$change <- ifelse(is.na(revisions$previous_rate_pct), "added",
    ifelse(is.na(revisions$current_rate_pct), "removed", "revised"))
  revisions <- revisions[is.na(revisions$previous_rate_pct) | is.na(revisions$current_rate_pct) |
    abs(revisions$current_rate_pct - revisions$previous_rate_pct) > 1e-10, ]
}
write.csv(revisions, file.path(run_dir, "revision_changes.csv"), row.names = FALSE)
source_status <- data.frame(publication_url = NHS_PUBLICATION_URL,
  input_mode = if (!is.null(LOCAL_WORKBOOK)) "local workbook" else if (DOWNLOAD_DATA) "download/cache" else "manual CSV",
  source_md5 = unname(tools::md5sum(if (exists("workbook_path", inherits = FALSE)) workbook_path else DATA_PATH)),
  first_month = format(min(df$month), "%Y-%m"), last_month = format(max(df$month), "%Y-%m"),
  n_months = nrow(df), revision_comparison = if (is.null(old)) "No previous successful snapshot" else "Compared with previous successful snapshot",
  changed_rows = nrow(revisions), definition_review = "See workbook notes; no automated assertion of stable definitions")
write.csv(source_status, file.path(run_dir, "source_status.csv"), row.names = FALSE)

# ---- 6-month rolling average ----
df <- df |>
  mutate(
    roll6_100k = slide_dbl(
      rate_100k,
      ~ mean(.x, na.rm = TRUE),
      .before = ROLL_BEFORE,
      .after  = ROLL_AFTER,
      .complete = ROLL_COMPLETE
    )
  )

# ============================================================================
# Baseline model ("observable norm") 2010–2018: seasonality + time
# ============================================================================
df_fit <- df |>
  filter(year(month) >= BASELINE_START, year(month) <= BASELINE_END) |>
  filter(!is.na(rate_100k), !is.na(month_num), !is.na(month_of_year))

fit_2010_2018 <- lm(rate_100k ~ month_num + month_of_year, data = df_fit)

slope_2010_2018 <- unname(coef(fit_2010_2018)["month_num"]) * 365.25

# Predictions: mean + CI + prediction interval
pred_ci <- as_tibble(predict(fit_2010_2018, newdata = df, interval = "confidence", level = 0.95)) |>
  rename(norm_fit_100k = fit, norm_ci_low_100k = lwr, norm_ci_high_100k = upr)

pred_pi <- as_tibble(predict(fit_2010_2018, newdata = df, interval = "prediction", level = 0.95)) |>
  rename(norm_pi_low_100k = lwr, norm_pi_high_100k = upr) |>
  select(norm_pi_low_100k, norm_pi_high_100k)

df <- bind_cols(df, pred_ci, pred_pi) |>
  mutate(
    period = case_when(
      year(month) >= BASELINE_START & year(month) <= BASELINE_END ~ "fit_2010_2018",
      year(month) >= PROJ_START     & year(month) <= PROJ_END     ~ "proj_2019_2025",
      TRUE ~ "other"
    )
  )

# Roll the baseline expectation and bands to match rolling plotted series
df <- df |>
  mutate(
    norm_fit_roll6     = slide_dbl(norm_fit_100k, mean, .before = ROLL_BEFORE, .after = ROLL_AFTER, .complete = ROLL_COMPLETE),
    norm_pi_low_roll6  = slide_dbl(norm_pi_low_100k, mean, .before = ROLL_BEFORE, .after = ROLL_AFTER, .complete = ROLL_COMPLETE),
    norm_pi_high_roll6 = slide_dbl(norm_pi_high_100k, mean, .before = ROLL_BEFORE, .after = ROLL_AFTER, .complete = ROLL_COMPLETE)
  )

# Deviations vs baseline
df <- df |>
  mutate(
    dev_100k_raw  = if_else(!is.na(rate_100k)  & !is.na(norm_fit_100k),  rate_100k  - norm_fit_100k,  NA_real_),
    dev_100k_roll = if_else(!is.na(roll6_100k) & !is.na(norm_fit_roll6), roll6_100k - norm_fit_roll6, NA_real_),
    dev_pct_roll  = if_else(!is.na(roll6_100k) & !is.na(norm_fit_roll6) & norm_fit_roll6 != 0,
                            100 * (roll6_100k - norm_fit_roll6) / norm_fit_roll6,
                            NA_real_),
    dev_label_roll = if_else(!is.na(dev_pct_roll), sprintf("%+.1f%%", dev_pct_roll), NA_character_)
  )

# ============================================================================
# EXPORT: Monthly sickness rolling + smoothed prediction-band endpoints
#   Output aligns with monthly positivity format: "YYYY-MM,value"
# ============================================================================

OUT_SICKNESS_BANDS <- file.path(run_dir, "sickness_roll6_prediction_band.csv")

df_sickness_monthly <- df |>
  mutate(month_key = format(month, "%Y-%m")) |>
  filter(year(month) >= PROJ_START, year(month) <= PROJ_END) |>
  select(
    month_key,
    roll6_100k,
    rate_100k,
    norm_pi_low_roll6,
    norm_pi_high_roll6
  ) |>
  rename(
    month = month_key,
    sickness_rate_100k = rate_100k,
    sickness_roll6_100k = roll6_100k,
    sickness_pred_low_100k = norm_pi_low_roll6,
    sickness_pred_high_100k = norm_pi_high_roll6
  ) |>
  arrange(month)

write.csv(df_sickness_monthly, OUT_SICKNESS_BANDS, row.names = FALSE)
cat(sprintf("\nExported monthly sickness rolling + prediction band: %s\n", OUT_SICKNESS_BANDS))

print(head(df_sickness_monthly, 24))

# ============================================================================
# STRUCTURAL BREAKS (Bai–Perron) ON DEVIATION-FROM-NORM (RAW)
#   IMPORTANT FIX: refit a FULL breakpoints model with breaks=optimal_n
#   so confint() works (no vcov error).
# ============================================================================
cat("\n========================================\n")
cat("BAI-PERRON STRUCTURAL BREAK ANALYSIS\n")
cat("========================================\n")
cat(sprintf("Baseline norm: %d-%d (seasonality + time)\n", BASELINE_START, BASELINE_END))
cat("Breakpoints run on deviation-from-norm (RAW): dev_100k_raw\n")
cat(sprintf("Model type: %s\n\n", BP_MODEL_TYPE))

df_for_bp <- df |>
  filter(!is.na(dev_100k_raw), !is.na(month), !is.na(month_num)) |>
  arrange(month) |>
  mutate(row_idx = row_number())

bp_formula <- if (BP_MODEL_TYPE == "level") {
  dev_100k_raw ~ 1
} else {
  dev_100k_raw ~ month_num
}

# Cap max breaks to what is feasible (prevents "changed to 9" surprises)
max_breaks_eff <- cap_max_breaks(n = nrow(df_for_bp), h = H_FRACTION, requested = MAX_BREAKS)

bp_full <- breakpoints(
  bp_formula,
  data   = df_for_bp,
  h      = H_FRACTION,
  breaks = max_breaks_eff
)

optimal_n <- get_optimal_breaks_bic(bp_full)
cat(sprintf("BIC-optimal number of breaks: %d\n", optimal_n))

if (optimal_n > 0) {
  
  # REFIT as FULL model with breaks=optimal_n (this is the key fix)
  bp_opt_full <- breakpoints(
    bp_formula,
    data   = df_for_bp,
    h      = H_FRACTION,
    breaks = optimal_n
  )
  
  break_indices <- bp_opt_full$breakpoints
  
  all_breaks <- tibble(
    break_num = seq_along(break_indices),
    row_idx   = break_indices,
    month     = df_for_bp$month[break_indices]
  ) |>
    arrange(month) |>
    mutate(
      mean_before = map_dbl(row_idx, ~ mean(df_for_bp$dev_100k_raw[max(1, .x - 6):max(1, .x - 1)], na.rm = TRUE)),
      mean_after  = map_dbl(row_idx, ~ mean(df_for_bp$dev_100k_raw[min(nrow(df_for_bp), .x + 1):min(nrow(df_for_bp), .x + 6)], na.rm = TRUE)),
      direction = if_else(mean_after > mean_before, "up", "down")
    )
  
  # CI now works (because bp_opt_full is "breakpointsfull")
  bp_conf <- tryCatch(confint(bp_opt_full, level = 0.95), error = function(e) {
    warning("Breakpoint intervals unavailable: ", conditionMessage(e))
    list(confint = matrix(NA_real_, nrow = optimal_n, ncol = 3))
  })
  
  all_breaks <- all_breaks |>
    mutate(
      ci_lower_month = df_for_bp$month[pmax(1, bp_conf$confint[, 1])],
      ci_upper_month = df_for_bp$month[pmin(nrow(df_for_bp), bp_conf$confint[, 3])]
    )
  
  cat("\n----------------------------------------\n")
  cat("DETECTED STRUCTURAL BREAKS (on dev_100k_raw):\n")
  cat("----------------------------------------\n")
  for (i in 1:nrow(all_breaks)) {
    direction_symbol <- if_else(all_breaks$direction[i] == "up", "↑ UP", "↓ DOWN")
    cat(sprintf("  Break %d: %s  %s\n",
                i,
                format(all_breaks$month[i], "%B %Y"),
                direction_symbol))
    cat(sprintf("           95%% CI: [%s - %s]\n",
                format(all_breaks$ci_lower_month[i], "%b %Y"),
                format(all_breaks$ci_upper_month[i], "%b %Y")))
  }
  
  cat("\nSegment coefficients (optimal model):\n")
  print(coef(bp_opt_full))
  
} else {
  cat("\nNo structural breaks detected by BIC criterion.\n")
  all_breaks <- tibble(
    break_num = integer(),
    row_idx   = integer(),
    month     = as.Date(character()),
    direction = character(),
    ci_lower_month = as.Date(character()),
    ci_upper_month = as.Date(character())
  )
}

cat("========================================\n\n")

# ============================================================================
# LOCAL EXTREMA (for deviation labels)
# ============================================================================
cat("Finding all local highs and lows...\n")

df <- df |>
  mutate(
    is_local_max = roll6_100k > lag(roll6_100k) & roll6_100k > lead(roll6_100k),
    is_local_min = roll6_100k < lag(roll6_100k) & roll6_100k < lead(roll6_100k)
  )

local_extrema <- df |>
  filter(year(month) >= BASELINE_START, !is.na(roll6_100k), !is.na(norm_fit_roll6)) |>
  filter(is_local_max | is_local_min) |>
  mutate(
    extrema_type = if_else(is_local_max, "high", "low"),
    dev_label_extrema = sprintf("%+.1f%%", dev_pct_roll)
  ) |>
  select(month, roll6_100k, norm_fit_roll6, extrema_type, dev_pct_roll, dev_label_extrema)

cat(sprintf("Found %d local extrema from %d onwards\n", nrow(local_extrema), BASELINE_START))

# ============================================================================
# SEGMENT TRENDS (on roll6_100k) with correct boundaries
# ============================================================================
cat("\n========================================\n")
cat("SEGMENT TREND ANALYSIS (on roll6_100k)\n")
cat("========================================\n")

segment_trends <- list()

if (nrow(all_breaks) > 0) {
  break_dates <- sort(all_breaks$month)
  boundaries <- c(min(df$month, na.rm = TRUE), break_dates, max(df$month, na.rm = TRUE))
  
  for (i in seq_len(length(boundaries) - 1)) {
    start_date <- boundaries[i]
    end_date   <- boundaries[i + 1]
    
    seg_name <- if (i == 1) "pre_B1" else sprintf("B%d_to_B%d", i - 1, i)
    
    df_segment <- df |>
      filter(month >= start_date, month < end_date, !is.na(roll6_100k), !is.na(month_num))
    
    if (nrow(df_segment) >= 3) {
      fit_segment <- lm(roll6_100k ~ month_num, data = df_segment)
      slope <- unname(coef(fit_segment)[2]) * 365.25
      segment_trends[[seg_name]] <- list(start = start_date, end = end_date, slope = slope, fit = fit_segment)
      
      cat(sprintf("  %-10s (%s to %s): %+.1f per 100k/year\n",
                  seg_name, format(start_date, "%b %Y"), format(end_date, "%b %Y"), slope))
    }
  }
  
  # Add "Bk_to_end"
  k <- length(break_dates)
  df_k_end <- df |> filter(month >= break_dates[k], !is.na(roll6_100k), !is.na(month_num))
  if (nrow(df_k_end) >= 3) {
    fit_k_end <- lm(roll6_100k ~ month_num, data = df_k_end)
    slope_k_end <- unname(coef(fit_k_end)[2]) * 365.25
    segment_trends[[sprintf("B%d_to_end", k)]] <- list(
      start = break_dates[k],
      end   = max(df$month, na.rm = TRUE),
      slope = slope_k_end,
      fit   = fit_k_end
    )
    cat(sprintf("  %-10s (%s to %s): %+.1f per 100k/year\n",
                sprintf("B%d_to_end", k), format(break_dates[k], "%b %Y"),
                format(max(df$month, na.rm = TRUE), "%b %Y"), slope_k_end))
  }
}

cat("\nBaseline slope (2010-2018):\n")
cat(sprintf("  %+0.1f per 100k/year\n", slope_2010_2018))
cat("========================================\n\n")

# ============================================================================
# MATERIALIZE SEGMENT TRENDS INTO df
#   (creates trend_B1_B2, trend_B2_B3, trend_B3_end, trend_B1_end, trend_B2_end)
# NOTE: Uses base R pipe (|>) compatible code (no { } blocks on RHS).
# ============================================================================

# Always create the columns so ggplot filters won't error
df <- df |>
  mutate(
    trend_B1_B2 = NA_real_,
    trend_B2_B3 = NA_real_,
    trend_B3_end = NA_real_,
    trend_B1_end = NA_real_,
    trend_B2_end = NA_real_
  )

add_segment_trend <- function(df, fit, colname, start_inclusive, end_exclusive = NULL) {
  pred <- predict(fit, newdata = df)
  in_range <- if (is.null(end_exclusive)) {
    df$month >= start_inclusive
  } else {
    df$month >= start_inclusive & df$month < end_exclusive
  }
  df[[colname]] <- ifelse(in_range, pred, NA_real_)
  df
}

# Helper: fit roll6 trend for a date range (base-pipe compatible)
fit_roll_trend <- function(df, start_date, end_date = NULL) {
  d <- df |>
    filter(month >= start_date)

  if (!is.null(end_date)) {
    d <- d |>
      filter(month < end_date)
  }

  d <- d |>
    filter(!is.na(roll6_100k), !is.na(month_num))

  if (nrow(d) < 3) return(NULL)
  lm(roll6_100k ~ month_num, data = d)
}

# Helper: safely fetch a segment fit from segment_trends
get_segment_fit <- function(segment_trends, key) {
  if (is.null(segment_trends[[key]])) return(NULL)
  if (is.null(segment_trends[[key]]$fit)) return(NULL)
  segment_trends[[key]]$fit
}

# Build the trend columns if breaks exist
if (exists("all_breaks") && nrow(all_breaks) > 0) {

  # ---- B1->B2 ----
  if (nrow(all_breaks) >= 2) {
    fit_b1_b2 <- get_segment_fit(segment_trends, "B1_to_B2")
    if (is.null(fit_b1_b2)) {
      fit_b1_b2 <- fit_roll_trend(df, all_breaks$month[1], all_breaks$month[2])
    }
    if (!is.null(fit_b1_b2)) {
      df <- add_segment_trend(df, fit_b1_b2, "trend_B1_B2", all_breaks$month[1], all_breaks$month[2])
    }
  }

  # ---- B2->B3 and B3->end ----
  if (nrow(all_breaks) >= 3) {
    fit_b2_b3 <- get_segment_fit(segment_trends, "B2_to_B3")
    if (is.null(fit_b2_b3)) {
      fit_b2_b3 <- fit_roll_trend(df, all_breaks$month[2], all_breaks$month[3])
    }
    if (!is.null(fit_b2_b3)) {
      df <- add_segment_trend(df, fit_b2_b3, "trend_B2_B3", all_breaks$month[2], all_breaks$month[3])
    }

    # B3->end
    fit_b3_end <- get_segment_fit(segment_trends, "B3_to_end")
    if (is.null(fit_b3_end)) {
      fit_b3_end <- get_segment_fit(segment_trends, sprintf("B%d_to_end", 3))
    }
    if (is.null(fit_b3_end)) {
      fit_b3_end <- fit_roll_trend(df, all_breaks$month[3], NULL)
    }
    if (!is.null(fit_b3_end)) {
      df <- add_segment_trend(df, fit_b3_end, "trend_B3_end", all_breaks$month[3], NULL)
    }
  }

  # ---- B1->end (always) ----
  fit_b1_end <- fit_roll_trend(df, all_breaks$month[1], NULL)
  if (!is.null(fit_b1_end)) {
    df <- add_segment_trend(df, fit_b1_end, "trend_B1_end", all_breaks$month[1], NULL)
    segment_trends[["B1_to_end"]] <- list(
      start = all_breaks$month[1],
      end   = max(df$month, na.rm = TRUE),
      slope = unname(coef(fit_b1_end)[2]) * 365.25,
      fit   = fit_b1_end
    )
  }

  # ---- B2->end (if exists) ----
  if (nrow(all_breaks) >= 2) {
    fit_b2_end <- fit_roll_trend(df, all_breaks$month[2], NULL)
    if (!is.null(fit_b2_end)) {
      df <- add_segment_trend(df, fit_b2_end, "trend_B2_end", all_breaks$month[2], NULL)
      segment_trends[["B2_to_end"]] <- list(
        start = all_breaks$month[2],
        end   = max(df$month, na.rm = TRUE),
        slope = unname(coef(fit_b2_end)[2]) * 365.25,
        fit   = fit_b2_end
      )
    }
  }
}

# =====================================================================
# VISUALIZATION — Substack/print friendly
# (self-contained palette so it runs even if you execute this block alone)
# =====================================================================

# ---- Palette ----
color_raw       <- "grey50"
color_rolling   <- "#E69F00"
color_baseline  <- "#0072B2"
color_markers   <- "#D55E00"
color_breaks    <- "#CC79A7"
color_segment_1 <- "#7B3294"
color_segment_2 <- "#009E73"
color_segment_3 <- "#F0E442"
color_B1_end    <- "#7B3294"
color_B2_end    <- "#1B9E77"
color_deviation <- "#666666"

markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax & mandates", "Mandates dropped")
)

y_range <- range(c(df$rate_100k, df$roll6_100k, df$norm_fit_roll6), na.rm = TRUE)
y_top <- y_range[2]
y_bottom <- y_range[1]

markers <- markers |>
  mutate(y_label = if_else(row_number() %% 2 == 1, y_top * 0.985, y_top * 0.935))

# --- X-axis guides: keep only January separators (cleaner for print) ---
x_min <- floor_date(min(df$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df$month, na.rm = TRUE), unit = "month")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")
grid_jan <- tibble(x = x_breaks_jan)

p <- ggplot(df, aes(x = month)) +

  # Year separators (subtle, behind everything)
  geom_vline(
    data = grid_jan,
    aes(xintercept = x),
    color = "grey85",
    linewidth = 0.7
  ) +

  # Prediction band (projection period only)
  geom_ribbon(
    data = df |> filter(period == "proj_2019_2025"),
    aes(ymin = norm_pi_low_roll6, ymax = norm_pi_high_roll6),
    fill = color_baseline, alpha = 0.12
  ) +

  # Raw data (subdued)
  geom_line(aes(y = rate_100k), color = color_raw, linewidth = 0.9, alpha = 0.45) +

  # Rolling average (prominent)
  geom_line(aes(y = roll6_100k), color = color_rolling, linewidth = 2.0, na.rm = TRUE) +

  # Baseline expected (rolling): solid in fit period
  geom_line(
    data = df |> filter(period == "fit_2010_2018"),
    aes(y = norm_fit_roll6),
    color = color_baseline, linewidth = 1.8
  ) +

  # Baseline expected (rolling): dashed in projection
  geom_line(
    data = df |> filter(period == "proj_2019_2025"),
    aes(y = norm_fit_roll6),
    color = color_baseline, linewidth = 1.8, linetype = "dashed"
  ) +

  # Segment trend lines
  geom_line(
    data = df |> filter(!is.na(trend_B1_B2)),
    aes(y = trend_B1_B2),
    color = color_segment_1, linewidth = 1.8
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B2_B3)),
    aes(y = trend_B2_B3),
    color = color_segment_2, linewidth = 1.8
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B3_end)),
    aes(y = trend_B3_end),
    color = color_segment_3, linewidth = 1.8
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B1_end)),
    aes(y = trend_B1_end),
    color = color_B1_end, linewidth = 1.5, linetype = "dashed"
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B2_end)),
    aes(y = trend_B2_end),
    color = color_B2_end, linewidth = 1.5, linetype = "dotdash"
  ) +

  # Deviation labels at local extrema
  geom_label_repel(
    data = local_extrema,
    aes(x = month, y = roll6_100k, label = dev_label_extrema),
    color = color_deviation,
    fill = alpha("white", 0.88),
    size = 3.6,
    fontface = "bold",
    label.padding = unit(0.18, "lines"),
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.25, "lines"),
    segment.color = alpha(color_deviation, 0.6),
    segment.size = 0.35,
    direction = "y",
    nudge_y = if_else(local_extrema$extrema_type == "high", 120, -120),
    min.segment.length = 0,
    max.overlaps = 30,
    seed = 42
  ) +

  # Event markers
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dotted", linewidth = 1.0,
    color = color_markers, alpha = 0.75
  ) +
  geom_label(
    data = markers,
    aes(x = date, y = y_label, label = short_label),
    color = color_markers,
    fill = alpha("white", 0.92),
    size = 3.6,
    fontface = "bold",
    label.padding = unit(0.18, "lines"),
    label.size = 0.25
  ) +

  # Structural breaks
  {
    if (nrow(all_breaks) > 0) {
      list(
        geom_vline(
          data = all_breaks,
          aes(xintercept = month),
          linetype = "longdash", linewidth = 1.2,
          color = color_breaks, alpha = 0.85
        ),
        geom_label(
          data = all_breaks |> mutate(
            y_pos = if_else(row_number() %% 2 == 1,
                            y_bottom + (y_top - y_bottom) * 0.16,
                            y_bottom + (y_top - y_bottom) * 0.07),
            label_text = sprintf("B%d %s  %s",
                                 break_num,
                                 if_else(direction == "up", "↑", "↓"),
                                 format(month, "%b %Y"))
          ),
          aes(x = month, y = y_pos, label = label_text),
          color = color_breaks,
          fill = alpha("white", 0.93),
          size = 3.0,
          fontface = "bold",
          label.padding = unit(0.18, "lines"),
          lineheight = 0.95
        )
      )
    }
  } +

  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0.02, 0.10))
  ) +

  labs(
    title = "NHS Sickness Absence Rate (England): Breaks vs 2010–2018 Norm",
    subtitle = "Rolling (6m) vs seasonality-adjusted baseline + prediction band; breaks detected on deviation-from-norm",
    x = NULL,
    y = "Absence days per 100,000 available days",
    caption = {
      cap <- sprintf("%d breaks (BIC-optimal) | Baseline slope: %+.0f per 100k/year",
                     nrow(all_breaks), slope_2010_2018)
      if (!is.null(segment_trends[["B1_to_B2"]])) cap <- paste0(cap, sprintf(" | B1→B2: %+.0f/yr", segment_trends[["B1_to_B2"]]$slope))
      if (!is.null(segment_trends[["B2_to_B3"]])) cap <- paste0(cap, sprintf(" | B2→B3: %+.0f/yr", segment_trends[["B2_to_B3"]]$slope))
      if (!is.null(segment_trends[["B1_to_end"]])) cap <- paste0(cap, sprintf(" | B1→end: %+.0f/yr", segment_trends[["B1_to_end"]]$slope))
      cap
    }
  ) +

  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 13, color = "grey35", margin = margin(b = 12)),
    plot.caption = element_text(size = 11, color = "grey40", hjust = 0, margin = margin(t = 12)),
    axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.5),
    plot.margin = margin(t = 14, r = 18, b = 18, l = 14),
    legend.position = "none"
  )

if (interactive()) print(p)

ggsave(
  file.path(run_dir, "plots", "nhs_sickness_analysis_substack.png"),
  plot = p,
  width = 18,
  height = 10,
  dpi = 320,
  bg = "white"
)

cat("\nPlot saved in ", file.path(run_dir, "plots"), "\n")

# ============================================================================
# SUMMARY OUTPUT
# ============================================================================
cat("\n")
cat("========================================================================\n")
cat("                          ANALYSIS SUMMARY                               \n")
cat("========================================================================\n")
cat(sprintf("Structural breaks detected (BIC-optimal): %d\n", nrow(all_breaks)))
if (nrow(all_breaks) > 0) {
  for (i in 1:nrow(all_breaks)) {
    cat(sprintf("  Break %d: %s (%s)\n",
                i,
                format(all_breaks$month[i], "%B %Y"),
                if_else(all_breaks$direction[i] == "up", "upward", "downward")))
  }
}
cat("------------------------------------------------------------------------\n")
cat("TREND SLOPES (per 100k per year):\n")
cat(sprintf("  Baseline norm (2010-2018): %+7.1f\n", slope_2010_2018))

if (length(segment_trends) > 0) {
  for (nm in names(segment_trends)) {
    tr <- segment_trends[[nm]]
    # guard against divide-by-zero if slope baseline is ~0
    multiplier <- if (isTRUE(all.equal(slope_2010_2018, 0))) NA_real_ else tr$slope / slope_2010_2018
    cat(sprintf("  %-24s %+7.1f  (%s)\n",
                paste0(nm, ":"),
                tr$slope,
                ifelse(is.na(multiplier), "n/a vs baseline", sprintf("%.1fx baseline", multiplier))))
  }
}
cat("========================================================================\n")

# ---- Wales-compatible exports (rates in percent; deviations in percentage points) ----
canonical <- data.frame(organisation_id = "E", organisation = "England", staff_id = "XXX",
  staff_group = "All staff", date = df$month, month_label = format(df$month, "%Y-%m"),
  period_code = format(df$month, "%Ym%m"), year = year(df$month), month = month(df$month),
  month_index = (year(df$month) - BASELINE_START) * 12L + month(df$month) - 1L, rate_pct = df$rate_pct)
selected <- canonical$year >= BASELINE_START & canonical$year <= PROJ_END
observations <- canonical[selected, ]
write.csv(observations, file.path(run_dir, "monthly_observations_2010_2025.csv"), row.names = FALSE, na = "")
canonical$baseline_pct <- df$norm_fit_100k / 1000
canonical$baseline_mean_low95 <- df$norm_ci_low_100k / 1000
canonical$baseline_mean_high95 <- df$norm_ci_high_100k / 1000
canonical$deviation_pp <- df$dev_100k_raw / 1000
canonical$deviation_relative_pct <- 100 * canonical$deviation_pp / canonical$baseline_pct
write.csv(canonical[selected, ], file.path(run_dir, "baseline_and_deviations.csv"), row.names = FALSE, na = "")
write.csv(df, file.path(run_dir, "legacy_analysis_series.csv"), row.names = FALSE, na = "")
co <- as.data.frame(summary(fit_2010_2018)$coefficients)
co <- data.frame(term = rownames(co), co, row.names = NULL)
co$units <- "rate_100k; time coefficient per calendar day"
write.csv(co, file.path(run_dir, "baseline_coefficients.csv"), row.names = FALSE)
annual <- do.call(rbind, lapply(split(canonical[selected, ], canonical$year[selected]), function(z) {
  b <- mean(z$baseline_pct); a <- mean(z$rate_pct)
  data.frame(organisation = "England", staff_id = "XXX", staff_group = "All staff", year = z$year[1],
    n_months = nrow(z), complete_year = nrow(z) == 12L, mean_monthly_rate_pct = a,
    mean_monthly_baseline_pct = b, mean_deviation_pp = a-b, relative_to_baseline_pct = 100*(a-b)/b)
}))
write.csv(annual, file.path(run_dir, "annual_mean_monthly_rates.csv"), row.names = FALSE)
bp_export <- as.data.frame(all_breaks)
bp_export$last_month_before <- format(bp_export$month, "%Y-%m")
bp_export$first_month_after <- format(df_for_bp$month[bp_export$row_idx + 1L], "%Y-%m")
write.csv(bp_export, file.path(run_dir, "breakpoints_primary.csv"), row.names = FALSE, na = "")
bic <- as.numeric(AIC(bp_full, breaks = 0:max_breaks_eff, k = log(nrow(df_for_bp))))
write.csv(data.frame(n_breaks = 0:max_breaks_eff, BIC = bic, selected = (0:max_breaks_eff) == optimal_n),
  file.path(run_dir, "breakpoint_BIC.csv"), row.names = FALSE)
endpoints <- c(0L, all_breaks$row_idx, nrow(df_for_bp))
segments <- do.call(rbind, lapply(seq_len(length(endpoints)-1L), function(i) {
  z <- df_for_bp[seq.int(endpoints[i]+1L, endpoints[i+1L]), ]
  data.frame(segment = i, first_month = format(min(z$month), "%Y-%m"), last_month = format(max(z$month), "%Y-%m"),
    n_months = nrow(z), mean_deviation_pp = mean(z$dev_100k_raw)/1000)
}))
write.csv(segments, file.path(run_dir, "segments_primary.csv"), row.names = FALSE)
saveRDS(list(baseline = fit_2010_2018, breakpoint_full = bp_full, selected_breaks = optimal_n,
  segment_trends = segment_trends), file.path(run_dir, "fitted_models.rds"))
write.csv(data.frame(error = character()), file.path(run_dir, "model_errors.csv"), row.names = FALSE)
# Portable HTML: embed the plot so this report works when copied on its own.
html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE)
}
html_table <- function(x) {
  if (!nrow(x)) return("<p>No rows.</p>")
  for (j in seq_along(x)) if (is.numeric(x[[j]])) x[[j]] <- format(round(x[[j]], 3), trim = TRUE)
  paste0("<div class='scroll'><table><tr>", paste0("<th>", html_escape(names(x)), "</th>", collapse = ""),
    "</tr>", paste(apply(x, 1, function(r) paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>")), collapse = ""), "</table></div>")
}
png_path <- file.path(run_dir, "plots", "nhs_sickness_analysis_substack.png")
encoded <- jsonlite::base64_enc(readBin(png_path, "raw", n = file.info(png_path)$size))
html <- paste0("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>England NHS sickness absence</title><style>body{font:16px/1.6 system-ui,sans-serif;max-width:1200px;margin:40px auto;padding:0 24px;color:#192c3c}h1,h2{line-height:1.2}table{border-collapse:collapse;font-size:14px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}.scroll{overflow:auto}img{width:100%}li{margin-bottom:10px}</style><body>",
  "<h1>England NHS sickness absence</h1><p>National all-staff series. Published observations: ",
  source_status$first_month, " to ", source_status$last_month, ". Run: ", html_escape(basename(run_dir)), ".</p>",
  "<h2>Read before interpreting</h2><ul>", paste0("<li>", html_escape(notes), "</li>", collapse = ""), "</ul>",
  "<h2>Original England analysis</h2><img alt='England observed rate, baseline and detected breaks' src='data:image/png;base64,", encoded, "'>",
  "<h2>Annual means of monthly rates</h2><p>Unweighted monthly means; these are not denominator-weighted annual absence rates. Incomplete years are explicitly flagged.</p>", html_table(annual),
  "<h2>Independently detected breaks</h2><p>The legacy plot labels the last month before each break. The first month after it is also exported; the comparison report consistently uses that latter convention.</p>",
  html_table(bp_export), "<h2>Input checks</h2>", html_table(quality), "<h2>Source and revision audit</h2>", html_table(source_status),
  html_table(revisions), "<p>Source publication: <a href='", html_escape(NHS_PUBLICATION_URL), "'>NHS England monthly sickness absence</a>. Source workbook and retrieval provenance are included in source/.</p>",
  "<h2>Exports</h2><p>CSV observations, baseline/deviations, annual means, coefficients, BIC, breakpoints and segments; fitted_models.rds; configuration.R; sessionInfo.txt; original PNG. Use compare_england_wales.R for harmonised monthly baselines and common-window comparisons.</p></body></html>")
writeLines(html, file.path(run_dir, "report.html"), useBytes = TRUE)
write.csv(current, old_path, row.names = FALSE)
writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(run_dir, "RUN_COMPLETE"))
completed <- TRUE
message("England exports: ", normalizePath(run_dir), "\nOpen report.html in your browser.")
invisible(list(run_dir = normalizePath(run_dir), data = canonical, breakpoints = bp_export, baseline = fit_2010_2018))
}

if (isTRUE(RUN_ANALYSIS)) england_result <- run_england_analysis()
