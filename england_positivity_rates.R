# ============================================================================
# uk_positivity_rates.R
# UKHSA COVID-19 (Nation -> England)
# The requested uk_ filename retains the existing England-only data coverage.
# Download:
#   1) PCR positivity (7-day rolling)  -> monthly mean
#   2) PCR tests (7-day total count)   -> monthly mean
#
# Outputs (created only if missing):
#   - data/tests_positivity_rates.csv      (month,value)
#   - data/tests_total_7day_mean.csv       (month,value)
#   - data/tests_positivity_and_volume.csv (month,positivity_mean_7day,tests_mean_7day_total)
#
# Download method: curl::curl_fetch_disk(), as in wales_sickness_rates.R.
# No httr2 / rlang dependency and no package updates during execution.
# If required, install once in a fresh R session:
#   install.packages(c("curl", "jsonlite"), type = "binary")  # Windows
#   install.packages(c("curl", "jsonlite"))                  # Linux/macOS
#
# API / metric documentation:
# https://ukhsa-dashboard.data.gov.uk/access-our-data/examples
# https://ukhsa-dashboard.data.gov.uk/metrics-documentation/covid-19-testing-pcr-count-by-day
# Despite its name, PCRcountByDay is a trailing 7-day count. Its monthly mean
# below is NOT the total number of tests performed during the calendar month.
# ============================================================================

need_download_packages <- function() {
  pkgs <- c("curl", "jsonlite")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    command <- paste0("install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), ")",
      if (.Platform$OS.type == "windows") ', type = "binary"' else "", ")")
    stop("Missing or unloadable download packages: ", paste(missing, collapse = ", "),
      "\nRestart R, run: ", command, "\nThen run this script again.", call. = FALSE)
  }
}

base_api  <- "https://api.ukhsa-dashboard.data.gov.uk"
theme     <- "infectious_disease"
sub_theme <- "respiratory"
topic     <- "COVID-19"
geo_type  <- "Nation"
geography <- "England"

page_size <- 100
DOWNLOAD_TIMEOUT <- 120
CONNECT_TIMEOUT <- 20
DOWNLOAD_RETRIES <- 3L
RUN_ANALYSIS <- getOption("positivity.run_analysis", TRUE)

metric_url <- function(metric) {
  sprintf(
    "%s/themes/%s/sub_themes/%s/topics/%s/geography_types/%s/geographies/%s/metrics/%s",
    base_api, theme, sub_theme, topic, geo_type, geography, metric
  )
}

fetch_file <- function(url, path) {
  # Same transport as Wales; avoids RStudio's .rs.downloadFile wrapper.
  # TLS certificate verification remains enabled. No stale cache on failure.
  problem <- ""
  for (attempt in seq_len(DOWNLOAD_RETRIES)) {
    if (file.exists(path)) unlink(path)
    h <- curl::new_handle(timeout = DOWNLOAD_TIMEOUT,
      connecttimeout = CONNECT_TIMEOUT, followlocation = TRUE,
      useragent = "R UKHSA positivity research (public statistical download)")
    curl::handle_setheaders(h, Accept = "application/json")
    ans <- tryCatch(curl::curl_fetch_disk(url, path, handle = h), error = identity)
    if (!inherits(ans, "error") && ans$status_code == 200L &&
        file.exists(path) && !is.na(file.info(path)$size) &&
        file.info(path)$size > 0) return(invisible(ans))
    problem <- if (inherits(ans, "error")) conditionMessage(ans) else
      if (ans$status_code == 200L) "HTTP 200 with an empty response" else
        paste("HTTP", ans$status_code)
    if (file.exists(path)) unlink(path)
    # Retry connection failures, rate limits and temporary server failures.
    # Permanent access/request errors need attention rather than retries.
    if (!inherits(ans, "error") && ans$status_code >= 400L &&
        ans$status_code < 500L && !ans$status_code %in% c(408L, 429L)) break
    if (attempt < DOWNLOAD_RETRIES) {
      message("Download attempt ", attempt, " failed: ", problem, "; retrying.")
      Sys.sleep(attempt)
    }
  }
  stop("UKHSA download failed: ", problem, "\nURL: ", url,
    "\nNo incomplete response or stale cache was substituted.", call. = FALSE)
}

fetch_page <- function(url, page) {
  need_download_packages()
  query <- list(format = "json", page_size = page_size, page = page,
    age = "all", sex = "all", stratum = "default")
  encoded <- vapply(query, function(x) curl::curl_escape(as.character(x)), character(1))
  request_url <- paste0(url, if (grepl("?", url, fixed = TRUE)) "&" else "?",
    paste(paste0(names(encoded), "=", encoded), collapse = "&"))
  path <- tempfile("ukhsa_page_", fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  fetch_file(request_url, path)
  tryCatch(jsonlite::fromJSON(path, simplifyDataFrame = TRUE),
    error = function(e) stop("Invalid UKHSA JSON response: ", conditionMessage(e),
      "\nURL: ", request_url, call. = FALSE))
}

fetch_all_results <- function(url) {
  first <- fetch_page(url, 1)

  if (!is.list(first) || is.null(first[["count"]]) || is.null(first[["results"]])) {
    stop("Unexpected response structure from page 1 at: ", url)
  }

  count <- suppressWarnings(as.numeric(first[["count"]]))
  if (length(count) != 1L || is.na(count) || !is.finite(count) ||
      count < 1 || count != floor(count)) {
    stop("No records or invalid record count returned for: ", url, call. = FALSE)
  }
  total_pages <- ceiling(count / page_size)

  message("\n--- Fetching: ", url)
  message("Total records: ", count, " | Page size: ", page_size, " | Pages: ", total_pages)

  all_results <- vector("list", total_pages)
  if (!is.data.frame(first[["results"]]) || nrow(first[["results"]]) == 0L) {
    stop("Missing results on page 1 at: ", url, call. = FALSE)
  }
  all_results[[1]] <- first[["results"]]

  if (total_pages >= 2) {
    for (p in 2:total_pages) {
      dat <- fetch_page(url, p)
      res <- dat[["results"]]

      if (!is.data.frame(res) || nrow(res) == 0L ||
          !identical(names(res), names(first[["results"]])) ||
          !identical(dat[["count"]], first[["count"]])) {
        stop("Incomplete or changed API response on page ", p, " at: ", url,
          "\nNo partial series was exported. Retry the download.", call. = FALSE)
      }
      all_results[[p]] <- res

      if (p %% 10 == 0 || p == total_pages) {
        message("Fetched page ", p, " / ", total_pages)
      }
    }
  }

  all_results <- Filter(Negate(is.null), all_results)
  df <- do.call(rbind, all_results)

  if (!is.data.frame(df) || nrow(df) != count) {
    stop("Incomplete UKHSA download: expected ", count, " records, received ",
      nrow(df), ".\nURL: ", url, call. = FALSE)
  }
  required <- c("year", "month", "date", "metric_value")
  if (!all(required %in% names(df))) stop("Missing UKHSA columns: ",
    paste(setdiff(required, names(df)), collapse = ", "), call. = FALSE)
  if (anyNA(df$date) || anyDuplicated(df$date)) stop(
    "Missing or duplicate dates in the UKHSA series. Check the API filters before averaging.",
    call. = FALSE)

  df
}

to_monthly_mean <- function(df) {
  # expects df has year, month, metric_value, date
  df$date <- as.Date(df$date)
  df$month_key <- sprintf("%04d-%02d", df$year, df$month)

  monthly <- aggregate(
    metric_value ~ month_key,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  monthly <- monthly[order(monthly$month_key), ]
  names(monthly) <- c("month", "value")
  monthly
}

write_if_missing <- function(df, out_file) {
  if (file.exists(out_file)) {
    message("File already exists: ", out_file, " (skipping)")
    return(invisible(FALSE))
  }
  write.csv(df, out_file, row.names = FALSE)
  message("Wrote: ", out_file)
  invisible(TRUE)
}

# ---- Execution (cached CSVs retain the original paths and column names) ----
run_positivity_rates <- function() {
  dir.create("data", showWarnings = FALSE, recursive = TRUE)

  # ---- Metrics ----
  metric_positivity <- "COVID-19_testing_positivity7DayRolling"
  metric_tests      <- "COVID-19_testing_PCRcountByDay"

  out_pos <- file.path("data", "tests_positivity_rates.csv")
  out_tst <- file.path("data", "tests_total_7day_mean.csv")
  out_join <- file.path("data", "tests_positivity_and_volume.csv")

  # ---- Positivity ----
  pos_monthly <- NULL
  if (file.exists(out_pos)) {
    pos_monthly <- read.csv(out_pos, stringsAsFactors = FALSE)
  } else {
    pos_daily <- fetch_all_results(metric_url(metric_positivity))
    pos_monthly <- to_monthly_mean(pos_daily)
    write_if_missing(pos_monthly, out_pos)
  }

  # ---- Tests ----
  tst_monthly <- NULL
  if (file.exists(out_tst)) {
    tst_monthly <- read.csv(out_tst, stringsAsFactors = FALSE)
  } else {
    tst_daily <- fetch_all_results(metric_url(metric_tests))
    tst_monthly <- to_monthly_mean(tst_daily)
    write_if_missing(tst_monthly, out_tst)
  }

  # ---- Combined (month join) ----
  combined <- merge(
    transform(pos_monthly, positivity_mean_7day = value)[, c("month", "positivity_mean_7day")],
    transform(tst_monthly, tests_mean_7day_total = value)[, c("month", "tests_mean_7day_total")],
    by = "month",
    all = TRUE
  )

  combined <- combined[order(combined$month), ]

  write_if_missing(combined, out_join)

  message("\nPreview combined:")
  print(utils::head(combined, 24))
  invisible(combined)
}

if (isTRUE(RUN_ANALYSIS)) run_positivity_rates()
