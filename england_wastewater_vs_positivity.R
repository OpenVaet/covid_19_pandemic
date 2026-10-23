# England wastewater RNA versus COVID-19 PCR test positivity
# Run: source("england_wastewater_vs_positivity.R")
# Required packages (install once in a clean R session):
# install.packages(c("curl", "jsonlite", "readxl"), repos = "https://cloud.r-project.org")
# No httr2, automatic package installation, R Markdown or Pandoc is required.
#
# Uses the published population-weighted England series, not an unweighted
# average of sewage treatment works. The workbook contains earlier site data,
# but national averages were not published before June 2021 because coverage
# changed. Original workbook and API responses are retained for reproducibility.

RUN_ANALYSIS <- TRUE
OUTPUT_ROOT <- "england_wastewater_results"
REFRESH_DOWNLOADS <- FALSE  # TRUE downloads fresh sources; existing runs stay intact
LOCAL_WORKBOOK <- ""       # Optional path to the original EMHP XLSX workbook
LOCAL_POSITIVITY_CSV <- ""  # Optional DAILY CSV: date, positivity_pct (0--100)
# A local positivity CSV must contain the same England/all-age/PCR/7-day metric.
# Monthly CSVs from england_positivity_rates.R cannot be used as daily data.
LOCAL_POSITIVITY_SOURCE <- "" # Required provenance description when using a local CSV
BOOTSTRAP_REPLICATES <- 2000L
BOOTSTRAP_BLOCK_WEEKS <- 4L
RANDOM_SEED <- 20261023L
LAG_DAYS <- -21L:21L  # Positive k compares wastewater(t) with positivity(t+k)

WORKBOOK_URL <- paste0("https://assets.publishing.service.gov.uk/media/",
  "6278ebfc8fa8f5206f0d9c50/EMHP_SARS-CoV-2_RNA_Concentration_May_2022.xlsx")
PUBLICATION_URL <- paste0("https://www.gov.uk/government/publications/",
  "monitoring-of-sars-cov-2-rna-in-england-wastewater-monthly-statistics-15-july-2020-to-30-march-2022")
METHODS_URL <- paste0(PUBLICATION_URL,
  "/emhp-wastewater-monitoring-of-sars-cov-2-in-england-15-july-2020-to-30-march-2022")
POSITIVITY_DOC_URL <- paste0("https://ukhsa-dashboard.data.gov.uk/metrics-documentation/",
  "covid-19-testing-positivity-7-day-rolling")
POLICY_URL <- paste0("https://www.gov.uk/government/news/",
  "confirmatory-pcr-tests-to-be-temporarily-suspended-for-positive-lateral-flow-test-results")
METRIC <- "COVID-19_testing_positivity7DayRolling"
API_URL <- paste0("https://api.ukhsa-dashboard.data.gov.uk/themes/infectious_disease/",
  "sub_themes/respiratory/topics/COVID-19/geography_types/Nation/geographies/England/metrics/", METRIC)

fail <- function(...) stop(..., call. = FALSE)
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
need_packages <- function() {
  pkgs <- c("curl", "jsonlite", "readxl")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) fail("Missing packages: ", paste(missing, collapse = ", "),
    ". Restart R, install them with install.packages(), then run again.")
}

# Download to a temporary file, check HTTP status, then replace the destination.
# No HTML publication-page download is needed; the XLSX attachment is used directly.
fetch_file <- function(url, destination, attempts = 3L) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  part <- tempfile("download_", tmpdir = dirname(destination))
  on.exit(unlink(part), add = TRUE)
  for (attempt in seq_len(attempts)) {
    handle <- curl::new_handle(timeout = 120, connecttimeout = 25,
      followlocation = TRUE, useragent = "England-wastewater-research/1.0",
      httpheader = c(Accept = "*/*"))
    ans <- tryCatch(curl::curl_fetch_disk(url, part, handle = handle), error = identity)
    status <- if (inherits(ans, "error")) NA_integer_ else ans$status_code
    if (!is.na(status) && status >= 200 && status < 300 && file.info(part)$size > 0) {
      if (!file.copy(part, destination, overwrite = TRUE)) fail("Cannot save ", destination)
      return(invisible(destination))
    }
    why <- if (inherits(ans, "error")) conditionMessage(ans) else paste("HTTP", status)
    retry <- is.na(status) || status %in% c(408, 429) || status >= 500
    if (!retry || attempt == attempts) fail("Download failed: ", why, "\nURL: ", url,
      "\nNo substitute data were used. A verified local input can be supplied in the settings.")
    message("Download attempt ", attempt, " failed: ", why, "; retrying.")
    Sys.sleep(min(2 ^ attempt, 8))
  }
}

parse_numeric <- function(x, label) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | x == ""
  out <- suppressWarnings(as.numeric(x))
  if (any(!missing & !is.finite(out))) fail("Unrecognised numeric value in ", label,
    ": ", paste(unique(x[!missing & !is.finite(out)]), collapse = ", "))
  out[missing] <- NA_real_
  out
}

read_wastewater <- function(path) {
  sheet <- "Regional_concentrations"
  if (!sheet %in% readxl::excel_sheets(path)) fail("Workbook is missing sheet: ", sheet)
  x <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_names = FALSE,
    col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
  first <- tolower(trimws(x[[1]]))
  header <- which(first == "region name")
  england <- which(first == "england")
  if (length(header) != 1L || length(england) != 1L) fail("Unexpected regional workbook layout.")
  labels <- as.character(unlist(x[header, -1], use.names = FALSE))
  dates <- as.Date(labels, format = "%d/%m/%Y")
  if (anyNA(dates) || anyDuplicated(dates) || any(diff(dates) != 1))
    fail("Invalid, duplicate or nonconsecutive workbook date headings.")
  values <- parse_numeric(unlist(x[england, -1], use.names = FALSE), "England wastewater")
  if (any(values <= 0, na.rm = TRUE)) fail("Nonpositive published England RNA concentration.")
  if (min(dates) != as.Date("2021-06-01") || max(dates) != as.Date("2022-03-30"))
    fail("Workbook date range differs from this historical release; review the source first.")
  data.frame(date = dates, wastewater_gc_l = values)
}

validate_positivity <- function(x) {
  if (!all(c("date", "positivity_pct") %in% names(x)))
    fail("Positivity input requires DAILY columns date and positivity_pct; monthly data are not accepted.")
  dates <- suppressWarnings(as.Date(as.character(x$date), format = "%Y-%m-%d"))
  if (anyNA(dates) || anyDuplicated(dates)) fail("Invalid or duplicate positivity dates.")
  values <- parse_numeric(x$positivity_pct, "positivity")
  if (any(values < 0 | values > 100, na.rm = TRUE)) fail("Positivity must be a percentage between 0 and 100.")
  if (length(dates) < 30L || median(diff(sort(dates))) != 1)
    fail("Positivity input is not a daily time series.")
  out <- data.frame(date = dates, positivity_pct = values)
  out[order(out$date), ]
}

read_api_pages <- function(folder) {
  manifest <- file.path(folder, "manifest.json")
  if (!file.exists(manifest)) fail("Incomplete API cache: missing manifest.")
  m <- jsonlite::fromJSON(manifest)
  if (!identical(m$api_url, API_URL) || !identical(m$metric, METRIC))
    fail("API cache provenance does not match England PCR positivity.")
  n_pages <- as.integer(m$pages)
  if (length(n_pages) != 1L || is.na(n_pages) || n_pages < 1L) fail("Invalid API manifest.")
  pieces <- lapply(seq_len(n_pages), function(i) {
    file <- file.path(folder, sprintf("page_%03d.json", i))
    if (!file.exists(file)) fail("Incomplete API cache: ", file)
    obj <- jsonlite::fromJSON(file)
    if (!isTRUE(obj$count == m$count) || !is.data.frame(obj$results) || !nrow(obj$results))
      fail("API count changed or page is empty; refresh the full cache.")
    required <- c("date", "metric_value", "geography", "metric", "age", "sex", "stratum")
    if (!all(required %in% names(obj$results))) fail("Unexpected API columns.")
    y <- obj$results
    expected <- c(geography = "England", metric = METRIC, age = "all", sex = "all", stratum = "default")
    for (field in names(expected))
      if (anyNA(y[[field]]) || any(y[[field]] != expected[[field]]))
        fail("Unexpected API dimension: ", field)
    data.frame(date = y$date, positivity_pct = y$metric_value)
  })
  x <- do.call(rbind, pieces)
  if (nrow(x) != m$count) fail("API pagination incomplete: expected ", m$count, " rows, received ", nrow(x))
  list(data = validate_positivity(x), manifest = m)
}

get_positivity <- function(cache, refresh) {
  folder <- file.path(cache, "positivity_api")
  if (!refresh && file.exists(file.path(folder, "manifest.json"))) return(read_api_pages(folder))
  staging <- tempfile("api_pending_", tmpdir = cache)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  url <- function(i) paste0(API_URL, "?format=json&page_size=100&page=", i,
    "&age=all&sex=all&stratum=default")
  message("Downloading UKHSA daily England PCR positivity...")
  fetch_file(url(1), file.path(staging, "page_001.json"))
  first <- jsonlite::fromJSON(file.path(staging, "page_001.json"))
  count <- first$count
  if (length(count) != 1L || !is.finite(count) || count < 1L || count != as.integer(count))
    fail("Invalid UKHSA result count.")
  pages <- ceiling(count / 100)
  if (pages > 1000) fail("Unexpectedly large API response.")
  if (pages > 1L) for (i in 2:pages) {
    message("  Positivity page ", i, " / ", pages)
    fetch_file(url(i), file.path(staging, sprintf("page_%03d.json", i)))
  }
  jsonlite::write_json(list(api_url = API_URL, metric = METRIC, count = count,
    pages = pages, retrieved_utc = utc_now(), filters = "England; age=all; sex=all; stratum=default"),
    file.path(staging, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  checked <- read_api_pages(staging)
  # Only replace the reusable cache after every page has passed validation.
  if (dir.exists(folder)) unlink(folder, recursive = TRUE)
  if (!file.rename(staging, folder)) fail("Cannot commit API cache.")
  checked
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}

# Paired moving-block bootstrap on the full, regular weekly calendar, with gaps
# retained. Differences are computed BEFORE resampling. This is an approximate
# dependence-aware interval, not an IID test or a forecast prediction interval.
block_ci <- function(x, y, method, B, block) {
  n <- length(x)
  if (sum(is.finite(x) & is.finite(y)) < max(12L, 3L * block) || n < block)
    return(c(lower = NA_real_, upper = NA_real_))
  draws <- replicate(B, {
    starts <- sample.int(n - block + 1L, ceiling(n / block), replace = TRUE)
    idx <- head(unlist(lapply(starts, function(s) seq.int(s, length.out = block))), n)
    safe_cor(x[idx], y[idx], method)
  })
  if (sum(is.finite(draws)) < 0.9 * B) return(c(lower = NA_real_, upper = NA_real_))
  setNames(as.numeric(quantile(draws, c(0.025, 0.975), na.rm = TRUE)), c("lower", "upper"))
}

correlation_rows <- function(d, frequency, B = 0L, block = 4L) {
  w <- d$wastewater_gc_l; p <- d$positivity_pct
  changes_w <- c(NA_real_, diff(log10(w)))
  changes_p <- c(NA_real_, diff(p))
  series <- list(list(w, p, "pearson", "Levels: raw RNA vs positivity"),
    list(log10(w), p, "pearson", "Levels: log10 RNA vs positivity"),
    list(w, p, "spearman", "Levels: rank correlation"),
    list(changes_w, changes_p, "pearson", "Changes: log10 RNA vs positivity percentage points"))
  do.call(rbind, lapply(series, function(s) {
    ok <- is.finite(s[[1]]) & is.finite(s[[2]])
    ci <- if (B > 0) block_ci(s[[1]], s[[2]], s[[3]], B, block) else c(NA_real_, NA_real_)
    data.frame(frequency = frequency, comparison = s[[4]], method = s[[3]], n = sum(ok),
      correlation = safe_cor(s[[1]], s[[2]], s[[3]]), lower_95 = ci[1], upper_95 = ci[2],
      first_date = if (any(ok)) as.character(min(d$date[ok])) else NA_character_,
      last_date = if (any(ok)) as.character(max(d$date[ok])) else NA_character_, row.names = NULL)
  }))
}

lag_scan <- function(d, lags = LAG_DAYS) {
  w <- log10(d$wastewater_gc_l); p <- d$positivity_pct
  dw <- w - w[match(d$date - 7L, d$date)]
  dp <- p - p[match(d$date - 7L, d$date)]
  y <- vapply(lags, function(k) p[match(d$date + k, d$date)], numeric(nrow(d)))
  dy <- vapply(lags, function(k) dp[match(d$date + k, d$date)], numeric(nrow(d)))
  # Same wastewater dates for every lag within each analysis: fair comparison.
  common <- is.finite(w) & apply(is.finite(y), 1, all)
  common_change <- is.finite(dw) & apply(is.finite(dy), 1, all)
  do.call(rbind, lapply(seq_along(lags), function(i) data.frame(
    lag_days = lags[i], n_levels = sum(common), n_changes = sum(common_change),
    pearson_log_levels = safe_cor(w[common], y[common, i]),
    spearman_levels = safe_cor(w[common], y[common, i], "spearman"),
    pearson_7day_changes = safe_cor(dw[common_change], dy[common_change, i]))))
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}
html_table <- function(x) {
  for (i in seq_along(x)) if (is.numeric(x[[i]]))
    x[[i]] <- ifelse(is.na(x[[i]]), "", format(round(x[[i]], 3), trim = TRUE, scientific = FALSE))
  paste0("<div class='table'><table><thead><tr>",
    paste0("<th>", html_escape(names(x)), "</th>", collapse = ""),
    "</tr></thead><tbody>", paste(apply(x, 1, function(row)
      paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")), collapse = ""),
    "</tbody></table></div>")
}
embedded_png <- function(path, alt) {
  raw <- readBin(path, "raw", n = file.info(path)$size)
  paste0("<img alt='", html_escape(alt), "' src='data:image/png;base64,",
    jsonlite::base64_enc(raw), "'>")
}
png_plot <- function(file, draw, height = 1000) {
  grDevices::png(file, width = 1600, height = height, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

run_england_wastewater_vs_positivity <- function(output_root = OUTPUT_ROOT,
    refresh = REFRESH_DOWNLOADS, local_workbook = LOCAL_WORKBOOK,
    local_positivity_csv = LOCAL_POSITIVITY_CSV, local_positivity_source = LOCAL_POSITIVITY_SOURCE,
    bootstrap_replicates = BOOTSTRAP_REPLICATES) {
  need_packages()
  if (length(bootstrap_replicates) != 1L || !is.finite(bootstrap_replicates) || bootstrap_replicates < 200)
    fail("Use at least 200 bootstrap replicates (default 2000).")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(output_root, "cache")
  dir.create(cache, showWarnings = FALSE)
  run <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = output_root)
  dir.create(run)
  for (sub in c("data", "plots", "sources")) dir.create(file.path(run, sub))
  completed <- FALSE
  on.exit(if (!completed) writeLines("Run did not complete; do not use partial results.",
    file.path(run, "RUN_FAILED.txt")), add = TRUE)
  message("Writing report to: ", normalizePath(run, winslash = "/"))

  workbook <- file.path(cache, "EMHP_SARS-CoV-2_RNA_Concentration_May_2022.xlsx")
  if (nzchar(local_workbook)) {
    if (!file.exists(local_workbook)) fail("Local workbook not found: ", local_workbook)
    workbook <- normalizePath(local_workbook, winslash = "/")
    workbook_origin <- workbook
    workbook_retrieved <- "User-supplied file; original retrieval time unknown"
  } else {
    meta_file <- paste0(workbook, ".json")
    if (refresh || !file.exists(workbook) || !file.exists(meta_file)) {
      message("Downloading historical England wastewater workbook...")
      candidate <- tempfile(fileext = ".xlsx", tmpdir = cache)
      on.exit(unlink(candidate), add = TRUE)
      fetch_file(WORKBOOK_URL, candidate)
      read_wastewater(candidate)
      if (!file.copy(candidate, workbook, overwrite = TRUE)) fail("Cannot save workbook cache.")
      jsonlite::write_json(list(url = WORKBOOK_URL, retrieved_utc = utc_now()),
        meta_file, auto_unbox = TRUE, pretty = TRUE)
    }
    meta <- jsonlite::fromJSON(meta_file)
    if (!identical(meta$url, WORKBOOK_URL)) fail("Workbook cache provenance mismatch.")
    workbook_origin <- WORKBOOK_URL
    workbook_retrieved <- meta$retrieved_utc
  }
  ww <- read_wastewater(workbook)
  if (!file.copy(workbook, file.path(run, "sources", basename(workbook)))) fail("Cannot snapshot workbook.")

  if (nzchar(local_positivity_csv)) {
    if (!nzchar(local_positivity_source)) fail("Set LOCAL_POSITIVITY_SOURCE for the local CSV.")
    pos <- validate_positivity(read.csv(local_positivity_csv, stringsAsFactors = FALSE))
    pos_origin <- local_positivity_source
    pos_retrieved <- "User-supplied file; original retrieval time unknown"
    if (!file.copy(local_positivity_csv, file.path(run, "sources", "local_positivity.csv")))
      fail("Cannot snapshot local positivity CSV.")
  } else {
    result <- get_positivity(cache, refresh)
    pos <- result$data
    pos_origin <- API_URL
    pos_retrieved <- result$manifest$retrieved_utc
    dir.create(file.path(run, "sources", "positivity_api"))
    api_files <- list.files(file.path(cache, "positivity_api"), full.names = TRUE)
    if (!all(file.copy(api_files, file.path(run, "sources", "positivity_api")))) fail("Cannot snapshot API files.")
  }
  write_csv(ww, file.path(run, "data", "england_wastewater_daily.csv"))
  write_csv(pos, file.path(run, "data", "england_positivity_daily.csv"))
  daily <- merge(ww, pos, by = "date", all.x = TRUE, sort = TRUE)
  daily$paired <- is.finite(daily$wastewater_gc_l) & is.finite(daily$positivity_pct)
  if (sum(daily$paired) < 60L) fail("Fewer than 60 overlapping daily observations.")
  daily$log10_wastewater <- log10(daily$wastewater_gc_l)
  daily$policy_period <- ifelse(daily$date < as.Date("2022-01-11"), "Before 11 January 2022", "From 11 January 2022")
  write_csv(daily, file.path(run, "data", "paired_daily.csv"))
  write_csv(daily[!daily$paired, ], file.path(run, "data", "missing_dates.csv"))

  # Sunday endpoints retain each source's published trailing-seven-day value.
  # These have nonoverlapping time windows, although serial dependence remains.
  weekly <- daily[as.POSIXlt(daily$date)$wday == 0L, ]
  set.seed(RANDOM_SEED)
  weekly_stats <- correlation_rows(weekly, "Weekly Sunday endpoints", bootstrap_replicates, BOOTSTRAP_BLOCK_WEEKS)
  daily_stats <- correlation_rows(daily, "Daily (overlapping 7-day windows)")
  phase_stats <- do.call(rbind, lapply(unique(daily$policy_period), function(period)
    correlation_rows(daily[daily$policy_period == period, ], paste("Daily:", period))))
  # 11--16 January windows mix the two policy regimes; omit them from sensitivity.
  phase_stats <- rbind(phase_stats, correlation_rows(daily[daily$date >= as.Date("2022-01-17"), ],
    "Daily: post-policy, allowing 6-day window washout"))
  # Check that choosing Sunday did not determine the association.
  weekday_stats <- do.call(rbind, lapply(0:6, function(day)
    correlation_rows(daily[as.POSIXlt(daily$date)$wday == day, ],
      paste("Weekly endpoint weekday", day, "(Sunday=0)"))))
  block_sensitivity <- do.call(rbind, lapply(c(2L, 6L), function(b)
    transform(correlation_rows(weekly, paste("Weekly bootstrap block", b), bootstrap_replicates, b), block_weeks = b)))
  lag <- lag_scan(daily)
  for (item in c("weekly", "weekly_stats", "daily_stats", "phase_stats", "weekday_stats", "block_sensitivity", "lag"))
    write_csv(get(item), file.path(run, "data", paste0(item, ".csv")))

  quality <- data.frame(check = c("Workbook calendar days", "Missing wastewater days", "Missing positivity days",
    "Paired daily observations", "Scheduled Sunday endpoints", "Paired Sunday endpoints"),
    value = c(nrow(daily), sum(!is.finite(daily$wastewater_gc_l)), sum(!is.finite(daily$positivity_pct)),
      sum(daily$paired), nrow(weekly), sum(weekly$paired)))
  write_csv(quality, file.path(run, "data", "quality_checks.csv"))
  policy_date <- as.Date("2022-01-11")
  blue <- "#176B87"; orange <- "#C15728"
  timeline <- file.path(run, "plots", "time_series.png")
  png_plot(timeline, function() {
    par(mfrow = c(2, 1), mar = c(3, 7, 3, 2), mgp = c(4.7, 1, 0), las = 1)
    plot(daily$date, daily$wastewater_gc_l, type = "l", log = "y", col = blue, lwd = 2,
      xlab = "", ylab = "RNA gene copies/litre (log scale)", main = "England wastewater: published population-weighted 7-day mean")
    abline(v = policy_date, col = "grey50", lty = 3)
    plot(daily$date, daily$positivity_pct, type = "l", col = orange, lwd = 2,
      xlab = "Specimen / wastewater date", ylab = "PCR positivity (%)", main = "England PCR positivity: published trailing 7-day percentage")
    abline(v = policy_date, col = "grey50", lty = 3)
    mtext("Dotted line: 11 January 2022 testing policy change. Gaps are not interpolated.", side = 1, line = 2, cex = 0.75)
  })
  scatter <- file.path(run, "plots", "weekly_association.png")
  png_plot(scatter, function() {
    par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), las = 1)
    col <- ifelse(weekly$date < policy_date, blue, orange)
    plot(log10(weekly$wastewater_gc_l), weekly$positivity_pct, pch = 19, col = col,
      xlab = "log10 RNA gene copies/litre", ylab = "PCR positivity (%)", main = "Weekly levels (Sunday endpoints)")
    legend("topleft", c("Before 11 Jan 2022", "From 11 Jan 2022"), pch = 19, col = c(blue, orange), bty = "n", cex = 0.75)
    plot(c(NA, diff(log10(weekly$wastewater_gc_l))), c(NA, diff(weekly$positivity_pct)),
      pch = 19, col = col, xlab = "Weekly change in log10 RNA", ylab = "Weekly positivity change (percentage points)",
      main = "Weekly changes")
    abline(h = 0, v = 0, col = "grey70", lty = 3)
  }, height = 800)
  lagplot <- file.path(run, "plots", "lag_correlations.png")
  png_plot(lagplot, function() {
    par(mar = c(5, 5, 3, 2), las = 1)
    matplot(lag$lag_days, lag[c("pearson_log_levels", "spearman_levels", "pearson_7day_changes")],
      type = "l", lwd = 2, lty = c(1, 2, 1), col = c(blue, "#738C40", orange), ylim = c(-1, 1),
      xlab = "Lag in days: positive = wastewater precedes positivity", ylab = "Correlation",
      main = "Exploratory lag comparison on common dates")
    abline(v = 0, h = 0, lty = 3, col = "grey60")
    legend("bottomleft", c("Pearson: log10 levels", "Spearman: levels", "Pearson: 7-day changes"),
      col = c(blue, "#738C40", orange), lty = c(1, 2, 1), lwd = 2, bty = "n", cex = 0.8)
  }, height = 800)

  r <- weekly_stats$correlation
  f <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "not estimable")
  peak_idx <- if (any(is.finite(lag$pearson_log_levels))) which.max(lag$pearson_log_levels) else NA_integer_
  peak_text <- if (is.na(peak_idx)) "The lag analysis was not estimable." else sprintf(
    "The largest in-sample log-level Pearson correlation in the prespecified -21 to +21 day scan is %s at lag %+d days (%d matched dates). This is exploratory and does not establish a forecasting lead.",
    f(lag$pearson_log_levels[peak_idx]), lag$lag_days[peak_idx], lag$n_levels[peak_idx])
  conclusion <- paste0("On weekly Sunday endpoints, contemporaneous Pearson correlation is ", f(r[1]),
    " on the raw RNA scale and ", f(r[2]), " after log10 transformation; Spearman correlation is ", f(r[3]),
    ". Correlation between weekly changes in log10 RNA and positivity is ", f(r[4]), ". ",
    if (is.finite(r[2]) && is.finite(r[4]) && abs(r[4]) < abs(r[2]))
      "The weaker association in changes is consistent with shared longer-term patterns contributing to the level correlation. " else "",
    "These associations describe this historical period; they do not calibrate wastewater into infection prevalence or validate a predictive model.")
  sources <- data.frame(source = c("Wastewater workbook", "PCR positivity"),
    location = c(workbook_origin, pos_origin), retrieved_utc = c(workbook_retrieved, pos_retrieved))
  write_csv(sources, file.path(run, "sources", "provenance.csv"))
  files <- list.files(file.path(run, "sources"), full.names = TRUE, recursive = TRUE)
  write_csv(data.frame(file = substring(files, nchar(run) + 2L), md5 = unname(tools::md5sum(files))),
    file.path(run, "sources", "checksums.csv"))
  capture.output(sessionInfo(), file = file.path(run, "session_info.txt"))
  jsonlite::write_json(list(created_utc = utc_now(), seed = RANDOM_SEED,
    bootstrap_replicates = bootstrap_replicates, primary_block_weeks = BOOTSTRAP_BLOCK_WEEKS,
    lags_days = LAG_DAYS, missing_values = "Retained; no interpolation; complete pairs within each analysis",
    positivity = METRIC), file.path(run, "settings.json"), auto_unbox = TRUE, pretty = TRUE)

  report <- paste0("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>England wastewater and PCR positivity</title><style>",
    "body{font:16px/1.6 system-ui,sans-serif;color:#19313d;max-width:1100px;margin:40px auto;padding:0 24px}",
    "h1,h2{line-height:1.2}h2{margin-top:2em}img{width:100%;height:auto}.lead{background:#edf5f6;padding:20px;border-left:5px solid #176b87}",
    ".table{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:9px;text-align:left;border-bottom:1px solid #d8e2e6}",
    "th{background:#edf5f6}a{color:#176b87}.muted{color:#5b6d75}</style><body>",
    "<h1>England wastewater and COVID-19 PCR positivity</h1><p class='muted'>Historical observational comparison | Created ", utc_now(), "</p>",
    "<div class='lead'><strong>Conclusion</strong><p>", html_escape(conclusion), "</p></div>",
    "<p>Paired observations: ", min(daily$date[daily$paired]), " to ", max(daily$date[daily$paired]),
    ". The national workbook spans 1 June 2021 to 30 March 2022. Earlier site observations remain in the downloaded workbook but are not pooled into a new national estimate.</p>",
    "<h2>Data quality and alignment</h2>", html_table(quality),
    "<p>Wastewater is the published England population-weighted average of site-level trailing seven-day means, in N1 gene copies per litre. Positivity is the percentage of positive PCR tests among positive and negative PCR tests in the trailing seven days, by specimen date. Lateral-flow tests are excluded. Both are aligned by the dates in their source, without an extra rolling smoother. Missing dates are retained and never interpolated; see data/missing_dates.csv.</p>",
    embedded_png(timeline, "Wastewater concentrations and PCR positivity over time"),
    "<h2>Primary comparison: weekly endpoints</h2><p>Every Sunday is selected to avoid overlapping seven-day windows between adjacent endpoints. Dates and gaps stay on a regular weekly grid. First differences use adjacent calendar weeks; a missing endpoint cannot be bridged. Serial correlation can still remain.</p>",
    html_table(weekly_stats),
    "<p>Intervals are approximate paired moving-block bootstrap 95% percentile intervals, using ", bootstrap_replicates,
    " replicates and ", BOOTSTRAP_BLOCK_WEEKS, "-week blocks. Missing values are retained inside blocks. These intervals assume approximate stability of the dependence structure; epidemic waves and policy shifts can violate that assumption. They are not prediction intervals. No IID correlation p-values are reported.</p>",
    embedded_png(scatter, "Weekly level and change associations"),
    "<h2>Daily comparison</h2><p>Daily correlations are descriptive: adjacent values reuse six days of observations. The changes row here uses one-day changes in the published seven-day series.</p>", html_table(daily_stats),
    "<h2>Exploratory timing</h2><p>", html_escape(peak_text), "</p>",
    "<p>At lag k, RNA on date t is compared with positivity on date t+k. A positive lag means wastewater comes first. All lags use the same eligible wastewater dates within the levels analysis and within the changes analysis. Changes are differences over exactly seven calendar days. A lag peak after smoothing and an in-sample scan is not evidence of causality or validated early warning.</p>",
    embedded_png(lagplot, "Correlations at lags from minus 21 to plus 21 days"),
    "<h2>Testing-policy sensitivity</h2><p>From 11 January 2022, confirmatory PCR testing following a positive lateral-flow result was suspended for most people in England. This can change who enters the PCR denominator. The table splits at that date and also excludes the six mixed post-change windows. These are sensitivity comparisons, not estimates of the policy's causal effect.</p>",
    html_table(phase_stats),
    "<details><summary>Weekly endpoint and bootstrap sensitivity</summary><h3>Alternative endpoint weekdays</h3>", html_table(weekday_stats),
    "<h3>Alternative block lengths</h3>", html_table(block_sensitivity), "</details>",
    "<h2>Interpretation and limitations</h2><ul>",
    "<li>PCR positivity measures positivity among people tested, not population infection prevalence. Testing access, behaviour and policy can change the association.</li>",
    "<li>Wastewater catchments and their population weights change through time. Rainfall, shedding, variants, laboratory methods and sampling coverage also affect the signal. An ecological national association does not establish an individual-level relationship.</li>",
    "<li>The workbook uses 160 gc/L for site-level nondetects from June 2021. This script does not replace or recalculate the official national series. Earlier methods and detection limits differed.</li>",
    "<li>Limited Christmas sampling affects the wastewater series. The report preserves the workbook's missing values and the official published averages, without inventing data.</li>",
    "<li>The May 2022 workbook includes revisions to prior releases. The positivity API is retrieved at run time and may also be revised. Source snapshots, retrieval times and MD5 checksums are saved for reproducibility; this is not a real-time-vintage analysis.</li>",
    "<li>The short overlap cannot establish a stable relationship across later variants or later testing regimes. A forecasting claim would require held-out evaluation.</li></ul>",
    "<h2>Sources and reproducibility</h2>", html_table(sources),
    "<ul><li><a href='", PUBLICATION_URL, "'>UKHSA historical wastewater publication and workbook</a></li>",
    "<li><a href='", METHODS_URL, "'>Wastewater methods, missingness and revisions</a></li>",
    "<li><a href='", POSITIVITY_DOC_URL, "'>PCR positivity metric definition</a></li>",
    "<li><a href='", POLICY_URL, "'>January 2022 testing policy change</a></li></ul>",
    "<p>The run folder contains data/*.csv, plots/*.png, source snapshots, checksums, settings.json and session_info.txt. This HTML embeds its plots and can be shared on its own. RUN_COMPLETE.txt is written only after all outputs succeed.</p></body></html>")
  writeLines(report, file.path(run, "report.html"), useBytes = TRUE)
  writeLines(c(paste("Completed", utc_now()), conclusion), file.path(run, "RUN_COMPLETE.txt"))
  completed <- TRUE
  message(conclusion)
  message("Report: ", normalizePath(file.path(run, "report.html"), winslash = "/"))
  invisible(list(run_dir = normalizePath(run, winslash = "/"), weekly = weekly_stats,
    daily = daily_stats, lag = lag, quality = quality))
}

if (isTRUE(RUN_ANALYSIS)) run_england_wastewater_vs_positivity()
