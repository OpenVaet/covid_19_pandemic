# ============================================================================
# wales_positivity_rates.R
# Wales PCR positivity and testing volume, from Public Health Wales (PHW).
# Uses curl downloads, like the other Wales scripts. No httr2 / rlang dependency.
#
# UKHSA's England endpoints return no Wales data for these metrics. This script
# therefore downloads PHW's COVID-19 workbook and reads its "Tests by date" tab.
# PHW describes this tab as test counts, not deduplicated case/episode counts:
# https://phw.nhs.wales/foi-disclosures/foi-2023-203-covid-deaths-tests-and-incidences/
#
# IMPORTANT (checked 2026-10-23): the PHW workbook host timed out/returned HTTP
# 503 during development. Live workbook download and its current column layout
# could not be verified. The script stops clearly if the source is unavailable
# or the schema is ambiguous. It supports an explicit local workbook/CSV input.
# It never substitutes England data, archived placeholder zeros, or LFD counts.
#
# Calculations matching the England script's time aggregation:
#   positivity7(t) = 100 * sum(positive PCR tests, t-6:t) / sum(PCR tests, t-6:t)
#   tests7(t)      = sum(PCR tests, t-6:t)
#   Monthly outputs are arithmetic means of available DAILY rolling values.
#   Missing calendar days remain NA; a rolling value needs all seven days.
#   These are NOT monthly positivity ratios or monthly numbers of tests.
#   PHW and UKHSA source populations/testing practices may still differ.
#
# Outputs in data/wales/ (separate from the England files in data/):
#   tests_positivity_rates.csv       (month,value)
#   tests_total_7day_mean.csv        (month,value)
#   tests_positivity_and_volume.csv  (same three columns as uk_positivity_rates.R)
#   daily_tests_and_rolling.csv, monthly_coverage.csv, source_status.csv
#   source/PHW_COVID_Testing_Dashboard_Data.xlsx (cached source when applicable)
#
# Install missing packages once, in a fresh R session; on Windows:
#   install.packages(c("curl", "xml2", "readxl"), type = "binary")
# Linux/macOS: omit type="binary" if binaries are unavailable.
# Run from the same project working directory as the other scripts.
# ============================================================================

OUT_DIR <- getOption("wales.positivity.output_dir", file.path("data", "wales"))
RUN_ANALYSIS <- getOption("wales.positivity.run_analysis", TRUE)
REFRESH_DOWNLOAD <- getOption("wales.positivity.refresh", FALSE)
# NULL downloads/reuses the source. A path uses an explicitly supplied file.
LOCAL_WORKBOOK <- getOption("wales.positivity.local_workbook", NULL)
LOCAL_CSV <- getOption("wales.positivity.local_csv", NULL)
# Manual CSV: date (YYYY-MM-DD), positive_tests, total_tests. Optional
# negative_tests: when present, positive + negative is the PCR denominator.
# Supply national Wales daily PCR test counts, not people, cases, or episodes.
PHW_SOURCE_PAGE <- paste0("https://www2.nphs.wales.nhs.uk/CommunitySurveillanceDocs.nsf/",
  "61c1e930f9121fd080256f2a004937ed/0a14d73f334101f18025885b002a31fb?OpenDocument")
# Optional direct workbook URL, if PHW moves the attachment.
PHW_WORKBOOK_URL <- getOption("wales.positivity.workbook_url", NULL)
SHEET <- "Tests by date"
HEADER_ROW <- NULL                  # NULL finds a unique date header in rows 1:40
# Optional exact Excel column names if automatic matching is ambiguous.
DATE_COLUMN <- NULL
POSITIVE_COLUMN <- NULL
TOTAL_COLUMN <- NULL
NEGATIVE_COLUMN <- NULL
DOWNLOAD_TIMEOUT <- 90
CONNECT_TIMEOUT <- 20
DOWNLOAD_RETRIES <- 3L

need_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first in a fresh R session: ",
    'install.packages(c(', paste(sprintf('"%s"', missing), collapse = ", "), ')',
    if (.Platform$OS.type == "windows") ', type = "binary"' else "", ")", call. = FALSE)
}
normalise <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(trimws(x))))
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "",
                                              fileEncoding = "UTF-8")

fetch_file <- function(url, path) {
  need_packages("curl")
  problem <- ""
  for (attempt in seq_len(DOWNLOAD_RETRIES)) {
    if (file.exists(path)) unlink(path)
    h <- curl::new_handle(timeout = DOWNLOAD_TIMEOUT, connecttimeout = CONNECT_TIMEOUT,
      followlocation = TRUE, useragent = "R Wales PCR positivity research (public statistical download)")
    ans <- tryCatch(curl::curl_fetch_disk(url, path, handle = h), error = identity)
    if (!inherits(ans, "error") && ans$status_code == 200L &&
        file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0)
      return(invisible(ans))
    problem <- if (inherits(ans, "error")) conditionMessage(ans) else paste("HTTP", ans$status_code)
    if (file.exists(path)) unlink(path)
    if (!inherits(ans, "error") && ans$status_code >= 400L && ans$status_code < 500L &&
        !ans$status_code %in% c(408L, 429L)) break
    if (attempt < DOWNLOAD_RETRIES) {
      message("PHW download attempt ", attempt, " failed: ", problem, "; retrying.")
      Sys.sleep(attempt)
    }
  }
  stop("PHW download failed: ", problem, "\nURL: ", url,
    "\nNo substitute data were used. Download the PHW COVID-19 workbook manually and set ",
    "LOCAL_WORKBOOK, or supply a national daily PCR CSV through LOCAL_CSV.\n",
    "Source guidance: https://phw.nhs.wales/foi-disclosures/foi-2023-203-covid-deaths-tests-and-incidences/",
    call. = FALSE)
}

discover_workbook_url <- function() {
  if (!is.null(PHW_WORKBOOK_URL)) return(PHW_WORKBOOK_URL)
  need_packages("xml2")
  page <- tempfile(fileext = ".html")
  on.exit(unlink(page), add = TRUE)
  fetch_file(PHW_SOURCE_PAGE, page)
  doc <- xml2::read_html(page)
  links <- xml2::xml_attr(xml2::xml_find_all(doc, "//a[@href]"), "href")
  links <- links[!is.na(links)]
  links <- unique(links[grepl("PHW_COVID_Testing_Dashboard_Data\\.xlsx",
                              utils::URLdecode(links), ignore.case = TRUE)])
  if (length(links) != 1L) stop("Could not identify one PHW testing workbook link on ",
    PHW_SOURCE_PAGE, ".\nSet PHW_WORKBOOK_URL to the direct workbook URL or use LOCAL_WORKBOOK.",
    call. = FALSE)
  xml2::url_absolute(links, PHW_SOURCE_PAGE)
}

pick_column <- function(x, aliases, explicit = NULL, required = TRUE, label) {
  hits <- if (!is.null(explicit)) which(names(x) == explicit) else
    which(normalise(names(x)) %in% aliases)
  if (length(hits) == 1L) return(names(x)[hits])
  if (!required && !length(hits) && is.null(explicit)) return(NULL)
  stop("Cannot identify one ", label, " column. Set the matching *_COLUMN option.\nColumns: ",
       paste(names(x), collapse = "; "), call. = FALSE)
}
date_aliases <- c("date", "specimendate", "sampledate", "testdate", "dateofspecimen",
                  "dateofsample", "dateoftest")

parse_dates <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
  text <- trimws(as.character(x))
  out <- as.Date(rep(NA_character_, length(text)))
  serial <- !is.na(text) & grepl("^[0-9]{5}(\\.[0-9]+)?$", text)
  out[serial] <- as.Date(as.numeric(text[serial]), origin = "1899-12-30")
  for (fmt in c("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%d-%b-%Y")) {
    take <- is.na(out) & !is.na(text) & nzchar(text)
    out[take] <- suppressWarnings(as.Date(text[take], format = fmt))
  }
  out
}
parse_counts <- function(x, label) {
  text <- trimws(as.character(x))
  missing <- is.na(text) | text %in% c("", "NA", "N/A", "..", "-", "*", "[c]", "[x]") |
    grepl("^<\\s*[0-9]+$", text)
  cleaned <- gsub(",", "", text, fixed = TRUE)
  value <- suppressWarnings(as.numeric(cleaned))
  if (any(!missing & (!is.finite(value) | value < 0 | value != floor(value))))
    stop("Invalid daily test count in ", label, "; expected non-negative integers.", call. = FALSE)
  value[missing] <- NA_real_
  value
}

standardise_counts <- function(x) {
  # A national row must be selected explicitly if the sheet also holds regions.
  for (col in names(x)[normalise(names(x)) %in% c("area", "areaname", "geography", "location",
                                                 "healthboard", "localhealthboard")]) {
    values <- trimws(as.character(x[[col]]))
    national <- tolower(values) %in% c("wales", "all wales", "w92000004")
    if (any(national, na.rm = TRUE)) x <- x[!is.na(national) & national, , drop = FALSE] else
      stop("No explicit Wales total in geography column: ", col,
           ". Supply a national daily PCR CSV through LOCAL_CSV.", call. = FALSE)
  }
  for (col in names(x)[normalise(names(x)) %in% c("testtype", "testmethod")]) {
    is_pcr <- tolower(trimws(as.character(x[[col]]))) %in% c("pcr", "polymerase chain reaction")
    if (!any(is_pcr, na.rm = TRUE)) stop("No explicit PCR rows in ", col, call. = FALSE)
    x <- x[!is.na(is_pcr) & is_pcr, , drop = FALSE]
  }
  date_col <- pick_column(x, date_aliases, DATE_COLUMN, label = "date")
  positive_col <- pick_column(x, c("positivetests", "positive", "positives", "numberpositive",
    "numberofpositivetests", "numberoftestspositive", "positivepcrtests", "pcrpositives",
    "teststestedpositive", "positivecount"), POSITIVE_COLUMN, label = "positive PCR tests")
  total_col <- pick_column(x, c("totaltests", "tests", "testcount", "numberoftests",
    "teststotal", "numbertested", "totalpcrtests", "pcrtests", "total"), TOTAL_COLUMN,
    required = FALSE, label = "total PCR tests")
  negative_col <- pick_column(x, c("negativetests", "negative", "negatives", "numbernegative",
    "numberofnegativetests", "negativepcrtests", "pcrnegatives", "negativecount"),
    NEGATIVE_COLUMN, required = FALSE, label = "negative PCR tests")
  if (is.null(total_col) && is.null(negative_col)) stop(
    "Need total PCR tests or negative PCR tests as the denominator. Case/episode counts are not accepted.",
    call. = FALSE)
  raw_date <- x[[date_col]]
  keep <- !is.na(raw_date) & nzchar(trimws(as.character(raw_date)))
  x <- x[keep, , drop = FALSE]
  dates <- parse_dates(x[[date_col]])
  if (!length(dates) || anyNA(dates) || anyDuplicated(dates)) stop(
    "Invalid or duplicate specimen dates. Select one Wales PCR series; do not sum overlapping groups.",
    call. = FALSE)
  positive <- parse_counts(x[[positive_col]], positive_col)
  reported_total <- if (is.null(total_col)) rep(NA_real_, nrow(x)) else
    parse_counts(x[[total_col]], total_col)
  negative <- if (is.null(negative_col)) rep(NA_real_, nrow(x)) else
    parse_counts(x[[negative_col]], negative_col)
  total <- if (!is.null(negative_col)) positive + negative else reported_total
  if (any(positive > total, na.rm = TRUE)) stop("Positive PCR tests exceed the denominator.", call. = FALSE)
  if (any(!is.na(reported_total) & !is.na(total) & total > reported_total)) stop(
    "Positive plus negative tests exceed reported total tests; inspect the source columns.", call. = FALSE)
  daily <- data.frame(date = dates, positive_tests = positive, total_tests = total,
    negative_tests = negative, reported_total_tests = reported_total, source_row_present = TRUE)
  daily <- daily[order(daily$date), ]
  if (as.numeric(max(daily$date) - min(daily$date)) > 10000) stop("Unexpected date span.")
  daily <- merge(data.frame(date = seq(min(daily$date), max(daily$date), by = "day")), daily,
                 by = "date", all.x = TRUE, sort = TRUE)
  daily$source_row_present[is.na(daily$source_row_present)] <- FALSE
  if (!any(daily$total_tests > 0, na.rm = TRUE)) stop("No positive testing volume in the source.")
  attr(daily, "denominator") <- if (!is.null(negative_col)) "positive plus negative PCR tests" else
    paste0("PHW reported total PCR tests (", total_col, ")")
  daily
}

read_phw_workbook <- function(path) {
  need_packages("readxl")
  sheets <- readxl::excel_sheets(path)
  hits <- sheets[normalise(sheets) == normalise(SHEET)]
  if (length(hits) != 1L) stop("Expected sheet '", SHEET, "'. Available sheets: ",
    paste(sheets, collapse = "; "), ". Set SHEET explicitly; do not substitute episode/case tables.",
    call. = FALSE)
  header <- HEADER_ROW
  if (is.null(header)) {
    preview <- readxl::read_excel(path, sheet = hits, col_names = FALSE,
                                 col_types = "text", n_max = 40, .name_repair = "minimal")
    hit <- which(apply(preview, 1, function(r) any(normalise(r) %in%
      if (is.null(DATE_COLUMN)) date_aliases else normalise(DATE_COLUMN))))
    if (length(hit) != 1L) stop("Cannot identify one header row in '", hits,
      "'. Set HEADER_ROW explicitly.", call. = FALSE)
    header <- hit
  }
  x <- as.data.frame(readxl::read_excel(path, sheet = hits, skip = header - 1L,
                                      .name_repair = "minimal"))
  standardise_counts(x)
}

add_rolling <- function(daily) {
  rolling_sum <- function(x) vapply(seq_along(x), function(i) {
    if (i < 7L) return(NA_real_)
    z <- x[seq.int(i - 6L, i)]
    if (any(!is.finite(z))) NA_real_ else sum(z)
  }, numeric(1))
  daily$positive_tests_7day <- rolling_sum(daily$positive_tests)
  daily$tests_7day_total <- rolling_sum(daily$total_tests)
  daily$positivity_7day <- ifelse(daily$tests_7day_total > 0,
    100 * daily$positive_tests_7day / daily$tests_7day_total, NA_real_)
  daily
}
monthly_tables <- function(daily) {
  mean_available <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  groups <- split(daily, format(daily$date, "%Y-%m"))
  combined <- do.call(rbind, lapply(names(groups), function(month) {
    z <- groups[[month]]
    data.frame(month = month, positivity_mean_7day = mean_available(z$positivity_7day),
               tests_mean_7day_total = mean_available(z$tests_7day_total))
  }))
  coverage <- do.call(rbind, lapply(names(groups), function(month) {
    z <- groups[[month]]
    bounds <- seq(as.Date(paste0(month, "-01")), by = "month", length.out = 2L)
    expected <- as.integer(diff(bounds))
    data.frame(month = month, calendar_days = expected,
      source_days = sum(z$source_row_present), valid_daily_tests = sum(is.finite(z$total_tests)),
      valid_positivity_7day_days = sum(is.finite(z$positivity_7day)),
      valid_tests_7day_days = sum(is.finite(z$tests_7day_total)),
      complete_month = sum(is.finite(z$positivity_7day)) == expected &&
                       sum(is.finite(z$tests_7day_total)) == expected)
  }))
  list(combined = combined, coverage = coverage)
}

run_wales_positivity_rates <- function() {
  if (!is.null(LOCAL_WORKBOOK) && !is.null(LOCAL_CSV)) stop("Set only one local source.")
  dir.create(file.path(OUT_DIR, "source"), recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(OUT_DIR, "source", "PHW_COVID_Testing_Dashboard_Data.xlsx")
  source_url <- PHW_SOURCE_PAGE
  mode <- "cache"
  if (!is.null(LOCAL_CSV)) {
    source <- normalizePath(LOCAL_CSV, mustWork = TRUE)
    daily <- standardise_counts(read.csv(source, stringsAsFactors = FALSE, check.names = FALSE))
    mode <- "manual national PCR CSV"
    source_url <- NA_character_
  } else if (!is.null(LOCAL_WORKBOOK)) {
    source <- normalizePath(LOCAL_WORKBOOK, mustWork = TRUE)
    daily <- read_phw_workbook(source)
    mode <- "manual PHW workbook"
  } else if (file.exists(cache) && !REFRESH_DOWNLOAD) {
    source <- cache
    daily <- read_phw_workbook(source)
  } else {
    source_url <- discover_workbook_url()
    temporary <- tempfile(tmpdir = file.path(OUT_DIR, "source"), fileext = ".xlsx")
    on.exit(unlink(temporary), add = TRUE)
    fetch_file(source_url, temporary)
    daily <- read_phw_workbook(temporary)  # Validate before replacing the cached workbook.
    if (!file.copy(temporary, cache, overwrite = TRUE)) stop("Could not save the PHW source workbook.")
    source <- cache
    mode <- "download"
  }
  denominator <- attr(daily, "denominator")
  daily <- add_rolling(daily)
  if (!any(is.finite(daily$positivity_7day)) || !any(is.finite(daily$tests_7day_total)))
    stop("No complete seven-day windows with valid positive and total PCR test counts.", call. = FALSE)
  tables <- monthly_tables(daily)
  combined <- tables$combined
  # No monthly output is written until source validation and calculations succeed.
  write_csv(data.frame(month = combined$month, value = combined$positivity_mean_7day),
            file.path(OUT_DIR, "tests_positivity_rates.csv"))
  write_csv(data.frame(month = combined$month, value = combined$tests_mean_7day_total),
            file.path(OUT_DIR, "tests_total_7day_mean.csv"))
  write_csv(combined, file.path(OUT_DIR, "tests_positivity_and_volume.csv"))
  write_csv(daily, file.path(OUT_DIR, "daily_tests_and_rolling.csv"))
  write_csv(tables$coverage, file.path(OUT_DIR, "monthly_coverage.csv"))
  write_csv(data.frame(geography = "Wales", input_mode = mode, source_url = source_url,
    source_file = normalizePath(source, winslash = "/", mustWork = TRUE),
    source_md5 = unname(tools::md5sum(source)),
    processed_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    first_date = min(daily$date), last_date = max(daily$date), denominator = denominator,
    missing_calendar_days = sum(!daily$source_row_present),
    missing_test_values = sum(!is.finite(daily$total_tests)),
    monthly_aggregation = "Arithmetic mean of available daily trailing-7-day values",
    comparison_note = "PCR tests, not case/episode counts; source coverage and testing policies may differ from England"),
    file.path(OUT_DIR, "source_status.csv"))
  message("Wales PCR data: ", min(daily$date), " to ", max(daily$date),
    "; ", nrow(combined), " monthly rows.\nOutputs: ", normalizePath(OUT_DIR),
    "\nCheck monthly_coverage.csv before comparing incomplete months.")
  print(utils::head(combined, 12))
  invisible(list(daily = daily, monthly = combined, coverage = tables$coverage))
}

if (isTRUE(RUN_ANALYSIS)) run_wales_positivity_rates()
