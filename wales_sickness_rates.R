# Wales NHS sickness absence: monthly analysis, January 2010--December 2025
# Version 1.1, 22 October 2026; adds England-style figures. Compatible with R 4.3; no Pandoc required.
# Install once: install.packages(c("curl", "xml2", "strucchange", "sandwich", "jsonlite", "ggplot2", "ggrepel"))
# Run the WHOLE file: source("wales_sickness_rates.R", encoding = "UTF-8")
#
# READ BEFORE ANALYSIS
# Source checked 22 October 2026: all-staff and six named Wales-wide staff groups
# each have 192 numeric monthly values. "No occupation code" has only 89 months.
# The script rechecks this at every run. Annual/quarterly rows are NOT modelled.
# The official download form offers unformatted values with greater precision
# than its displayed 0.1-point rates. This script requests unformatted values.
# Outcome: FTE calendar days lost / available FTE calendar days, in percent.
# It is not the proportion of workers who were ill. Self-isolation is excluded.
# Staff classification and employer changes can affect comparisons: ambulance,
# administration and support groups deserve particular review. National totals
# avoid local boundary changes but not workforce-composition changes.
# Revisions flagged by the publisher are retained. A first snapshot cannot
# reconstruct revision history; later runs compare against the previous snapshot.
#
# METHODS
# Separate OLS baseline for each group: rate_pct ~ month_index + month_factor,
# fitted ONLY on January 2010--December 2018 (108 months), projected through 2025.
# Newey-West lag 12 covariance gives 95% baseline MEAN intervals, not prediction
# intervals. Monthly residuals (observed minus baseline) are NOT smoothed.
# Bai-Perron multiple level shifts: deviation_pp ~ 1 over 2010--2025. BIC selects
# 0..6 breaks independently for each group; no English dates are supplied.
# Primary minimum segment length: 24 months, matching the Dutch 8-quarter run.
# Sensitivities: 12 and 18 months. Compare them BEFORE interpreting onset dates:
# a minimum duration can force adjacent breaks apart and move estimated dates.
# These are level shifts in deviations, not changes in trend slope.
# Break-date intervals are approximate and conditional on the baseline and
# selected break count; they do not propagate baseline/selection uncertainty.
# Serial dependence, changing seasonality and baseline misspecification can
# affect BIC and break dates. HAC intervals do not make BIC serially robust.
# All staff and its components overlap; they are not independent replications.
# No pathogen-specific or causal interpretation follows from all-cause absence.
#
# SOURCES
# https://stats.gov.wales/en-GB/d4578da5-380f-439f-a056-8c5dd74d6509
# https://www.gov.wales/nhs-workforce-staff-directly-employed-nhs-wales-quality-report-html
# https://www.gov.wales/sickness-absence-nhs-quality-report-html
# https://cran.r-project.org/web/packages/strucchange/strucchange.pdf

# ---- Configuration ---------------------------------------------------------
OUT_DIR <- getOption("wales.output_dir", "wales_sickness_results")
SOURCE_MODE <- getOption("wales.source_mode", "refresh")          # "refresh", "cache", or "manual"
MANUAL_CSV <- "data/wales_sickness_raw.csv"
# Manual fallback: on the StatsWales page choose Full dataset, CSV, Unformatted
# numbers, Include reference codes and hierarchies = Yes, English. Save above.
# Then set SOURCE_MODE <- "manual". Never substitute the rolling annual series.
ORGANISATION_IDS <- "W"           # W = published Wales total; never sum rates
STAFF_IDS <- NULL                 # NULL = all named groups plus All staff
INCLUDE_UNCODED <- FALSE
START_YEAR <- 2010L
END_YEAR <- 2025L
BASELINE_END <- 2018L
MIN_SEGMENT_MONTHS <- 24L
SENSITIVITY_MONTHS <- c(12L, 18L)
MAX_BREAKS <- 6L
HAC_LAG <- 12L
MAKE_PLOTS <- TRUE
RUN_ANALYSIS <- getOption("wales.run_analysis", TRUE)              # FALSE loads functions without running

# England-style figures supplement the existing diagnostic figures.
MAKE_ENGLAND_STYLE_PLOTS <- TRUE
# Shared UK reference events; add verified country-specific rows if desired.
# Set STYLE_EVENTS <- STYLE_EVENTS[0, ] to hide all event markers.
STYLE_EVENTS <- data.frame(date=as.Date(c("2020-01-31","2020-12-08")),
  label=c("First UK cases","UK vaccine rollout"),
  source=c("https://www.gov.uk/government/news/cmo-confirms-cases-of-coronavirus-in-england",
    "https://www.gov.uk/government/news/uk-marks-one-year-since-deploying-worlds-first-covid-19-vaccine"))

DATASET_URL <- "https://stats.gov.wales/en-GB/d4578da5-380f-439f-a056-8c5dd74d6509"
QUALITY_URL <- "https://www.gov.wales/nhs-workforce-staff-directly-employed-nhs-wales-quality-report-html"

need_packages <- function() {
  pkgs <- c("curl", "xml2", "strucchange", "sandwich", "jsonlite")
  if (MAKE_ENGLAND_STYLE_PLOTS) pkgs <- c(pkgs, "ggplot2", "ggrepel")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages first: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}
write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
utc_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
text_file <- function(path) paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
normalise_text <- function(x) trimws(gsub("[[:space:]]+", " ", x))
require_cols <- function(x, columns) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop("Unexpected CSV schema: missing ", paste(absent, collapse = ", "),
    ". Download English CSV with reference codes and hierarchies included.", call. = FALSE)
}
fetch_file <- function(url, path, form = NULL) {
  # curl avoids RStudio's .rs.downloadFile wrapper; TLS verification stays on.
  problem <- ""
  for (attempt in seq_len(3L)) {
    h <- curl::new_handle(timeout = 120, connecttimeout = 20, followlocation = TRUE,
      useragent = "R Wales NHS sickness research (public statistical download)")
    if (!is.null(form)) {
      body <- paste(paste0(names(form), "=", vapply(form, curl::curl_escape, character(1))), collapse = "&")
      curl::handle_setopt(h, postfields = body)
      curl::handle_setheaders(h, "Content-Type" = "application/x-www-form-urlencoded")
    }
    ans <- tryCatch(curl::curl_fetch_disk(url, path, handle = h), error = identity)
    if (!inherits(ans, "error") && ans$status_code == 200L &&
        file.exists(path) && file.info(path)$size > 0) return(invisible(ans))
    problem <- if (inherits(ans, "error")) conditionMessage(ans) else paste("HTTP", ans$status_code)
    if (file.exists(path)) unlink(path)
    # Do not repeatedly retry permanent access or request errors.
    if (!inherits(ans, "error") && ans$status_code %in% c(400, 401, 403, 404, 405)) break
    if (attempt < 3L) Sys.sleep(attempt)
  }
  stop("Download failed: ", problem, "\nURL: ", url,
       "\nNo stale cache was substituted. See SOURCE_MODE and MANUAL_CSV at the top of the script.",
       call. = FALSE)
}
read_source_csv <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character(), fileEncoding = "UTF-8-BOM")
  names(x) <- sub("^\ufeff", "", names(x))
  require_cols(x, c("Data values", "Data description", "Staff group", "Staff group_reference",
    "Date", "Date_reference", "Organisation", "Organisation_reference", "Notes"))
  if (!nrow(x)) stop("Downloaded CSV is empty")
  if (any(trimws(x[["Data description"]]) != "Percentage of Staff Absent"))
    stop("Unexpected measure: inspect the downloaded file before analysis")
  x
}
load_source <- function(raw_dir, cache_dir) {
  cache <- file.path(cache_dir, "source.rds")
  csv <- file.path(raw_dir, "statswales_full.csv")
  if (SOURCE_MODE == "cache") {
    if (!file.exists(cache)) stop("No cache yet. Use refresh or manual mode first.")
    obj <- readRDS(cache)
    writeLines(obj$csv_text, csv, useBytes = TRUE)
    writeLines(obj$metadata, file.path(raw_dir, "metadata.csv"), useBytes = TRUE)
    writeLines(obj$quality_text, file.path(raw_dir, "quality_report.txt"), useBytes = TRUE)
    return(obj)
  }
  if (SOURCE_MODE == "manual") {
    if (!file.exists(MANUAL_CSV)) stop("Manual CSV not found: ", MANUAL_CSV)
    if (!file.copy(MANUAL_CSV, csv)) stop("Could not copy manual CSV")
  } else {
    # Verified public download form, not an invented API. No scraping of table pages.
    fetch_file(paste0(DATASET_URL, "/download"), csv, list(view_type = "unfiltered",
      format = "csv", view_choice = "raw", extended = "yes", download_language = "en-GB",
      selected_filter_options = "[]", action = "download"))
  }
  x <- read_source_csv(csv)
  meta <- ""; quality <- ""; meta_issues <- character()
  if (SOURCE_MODE == "refresh") {
    tryCatch({
      p <- file.path(raw_dir, "metadata.csv")
      fetch_file(paste0(DATASET_URL, "/download/metadata"), p)
      meta <- text_file(p)
    }, error = function(e) meta_issues <<- c(meta_issues, paste("Metadata:", conditionMessage(e))))
    tryCatch({
      p <- file.path(raw_dir, "quality_report.html")
      fetch_file(QUALITY_URL, p)
      doc <- xml2::read_html(p)
      main <- xml2::xml_find_first(doc, "//main")
      if (inherits(main, "xml_missing")) stop("Quality report main content not found")
      quality <- normalise_text(xml2::xml_text(main))
      writeLines(quality, file.path(raw_dir, "quality_report.txt"))
    }, error = function(e) meta_issues <<- c(meta_issues, paste("Quality report:", conditionMessage(e))))
  } else meta_issues <- "Manual import: source release date/metadata not automatically verified."
  obj <- list(data = x, csv_text = text_file(csv), metadata = meta, quality_text = quality,
    metadata_issues = meta_issues, retrieved_utc = utc_now(), origin_mode = SOURCE_MODE,
    source_url = DATASET_URL, csv_md5 = unname(tools::md5sum(csv)))
  obj
}

month_grid <- function() {
  dates <- seq(as.Date(sprintf("%d-01-01", START_YEAR)),
               as.Date(sprintf("%d-12-01", END_YEAR)), by = "month")
  data.frame(period_code = format(dates, "%Ym%m"), date = dates,
    month_label = format(dates, "%Y-%m"), year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")), month_index = seq_along(dates) - 1L)
}
parse_rate <- function(x) {
  clean <- trimws(gsub("\\[[^]]*\\]", "", x))
  # Keep suppression, blanks and non-numeric values missing, never zero.
  clean <- sub("%$", "", clean)
  valid <- grepl("^[+-]?[0-9]+([.][0-9]+)?$", trimws(clean))
  out <- rep(NA_real_, length(x)); out[valid] <- as.numeric(clean[valid]); out
}
audit_source <- function(obj) {
  x <- obj$data
  monthly <- grepl("^[0-9]{4}m(0[1-9]|1[0-2])$", x[["Date_reference"]])
  x <- x[monthly & x[["Date_reference"]] %in% month_grid()$period_code &
           x[["Organisation_reference"]] %in% ORGANISATION_IDS, , drop = FALSE]
  if (!nrow(x)) stop("No monthly observations for the requested organisations/window")
  if (!all(ORGANISATION_IDS %in% x[["Organisation_reference"]])) stop("An organisation ID was not found")
  d <- data.frame(organisation_id = x[["Organisation_reference"]], organisation = x$Organisation,
    staff_id = x[["Staff group_reference"]], staff_group = x[["Staff group"]],
    period_code = x[["Date_reference"]], rate_raw = x[["Data values"]], notes = x$Notes,
    stringsAsFactors = FALSE)
  d$rate_pct <- parse_rate(d$rate_raw)
  d$revised_flag <- grepl("revis|\\[r\\]", paste(d$notes, d$rate_raw), ignore.case = TRUE)
  d$provisional_flag <- grepl("provision|\\[p\\]", paste(d$notes, d$rate_raw), ignore.case = TRUE)
  d$row_present <- TRUE
  if (anyDuplicated(d[, c("organisation_id", "staff_id", "period_code")]))
    stop("Duplicate organisation/staff/month keys: inspect the raw CSV; no averaging performed")
  ids <- unique(d[, c("organisation_id", "organisation", "staff_id", "staff_group")])
  if (anyDuplicated(ids[, c("organisation_id", "staff_id")]))
    stop("Labels change within a staff/organisation code; review source definitions")
  if (!is.null(STAFF_IDS) && !all(STAFF_IDS %in% ids$staff_id)) stop("A requested STAFF_IDS value was not found")
  pieces <- lapply(seq_len(nrow(ids)), function(i) {
    z <- merge(month_grid(), d[d$organisation_id == ids$organisation_id[i] &
      d$staff_id == ids$staff_id[i], c("period_code", "rate_raw", "notes", "rate_pct",
      "revised_flag", "provisional_flag", "row_present")], by = "period_code", all.x = TRUE)
    z <- z[order(z$month_index), ]
    z$row_present[is.na(z$row_present)] <- FALSE
    z$revised_flag[is.na(z$revised_flag)] <- FALSE
    z$provisional_flag[is.na(z$provisional_flag)] <- FALSE
    z <- cbind(ids[rep(i, nrow(z)), ], z)
    z$selected <- if (is.null(STAFF_IDS)) INCLUDE_UNCODED | z$staff_group != "No occupation code" else z$staff_id %in% STAFF_IDS
    z
  })
  full <- do.call(rbind, pieces); rownames(full) <- NULL
  quality <- do.call(rbind, lapply(pieces, function(z) {
    invalid <- !is.finite(z$rate_pct) | z$rate_pct < 0 | z$rate_pct > 100
    missing <- is.na(z$rate_pct)
    eligible <- all(z$selected) && !any(invalid)
    data.frame(z[1, c("organisation_id", "organisation", "staff_id", "staff_group")],
      expected_months = nrow(z), rows_present = sum(z$row_present), numeric_months = sum(!missing),
      missing_months = sum(missing), missing_baseline = sum(missing & z$year <= BASELINE_END),
      invalid_numeric = sum(!missing & invalid), revised_flags = sum(z$revised_flag),
      provisional_flags = sum(z$provisional_flag), selected = z$selected[1], eligible = eligible,
      decision = if (!z$selected[1]) "Not selected (uncoded group by default)" else
        if (any(invalid)) "Excluded: incomplete/invalid; no imputation" else "Complete: eligible",
      missing_month_labels = paste(z$month_label[missing], collapse = "; "), row.names = NULL)
  }))
  list(data = full, quality = quality)
}
known_changes <- function() {
  # These workforce coding notes flag comparability risks; their quantitative
  # effects on this particular absence table are not assumed to be known.
  data.frame(period = c("2015 onward", "2015-2017", "April 2019", "2021", "COVID-19 period"),
    relevant_groups = c("Medical and dental", "Ambulance / support", "Ambulance / administration",
      "Local organisations", "All staff groups"),
    issue = c("GP trainee employment coverage changed under the lead-employer scheme.",
      "Support staff were reclassified as ambulance personnel.",
      "Call-taking and dispatch roles moved from administration into ambulance coding.",
      "Shared-services staff moved from Velindre to a separately reported employer.",
      "Self-isolation is excluded; separate pandemic management data use another method."),
    interpretation = "Review against original documentation; not an imposed statistical breakpoint.",
    source = QUALITY_URL, stringsAsFactors = FALSE)
}
revision_audit <- function(full, obj, cache_dir, run_dir) {
  previous <- file.path(cache_dir, "previous_snapshot.rds")
  vars <- c("rate_pct", "rate_raw", "notes", "row_present", "staff_group", "organisation")
  keys <- c("organisation_id", "staff_id", "period_code")
  scope <- list(ORGANISATION_IDS = sort(ORGANISATION_IDS), START_YEAR = START_YEAR, END_YEAR = END_YEAR)
  changes <- data.frame(organisation_id = character(), staff_id = character(), period_code = character(),
    field = character(), previous = character(), current = character())
  note <- "First snapshot: historical revisions cannot be assessed. Publisher revision flags are retained."
  if (file.exists(previous)) {
    old <- readRDS(previous)
    if (identical(scope, old$scope)) {
      z <- merge(old$data, full[, c(keys, vars)], by = keys, all = TRUE, suffixes = c("_old", "_new"))
      for (v in vars) {
        a <- as.character(z[[paste0(v, "_old")]]); b <- as.character(z[[paste0(v, "_new")]])
        same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
        ix <- which(!same)
        if (length(ix)) changes <- rbind(changes,
          data.frame(z[ix, keys], field = v, previous = a[ix], current = b[ix]))
      }
      note <- paste("Compared with snapshot imported", old$retrieved_utc,
        "- changed observation fields:", nrow(changes),
        "- metadata text changed:", !identical(obj$metadata, old$metadata),
        "- quality report text changed:", !identical(obj$quality_text, old$quality_text),
        ". Text changes require review; publication-date updates also count.")
      writeLines(old$metadata, file.path(run_dir, "previous_metadata.csv"), useBytes = TRUE)
      writeLines(old$quality_text, file.path(run_dir, "previous_quality_report.txt"), useBytes = TRUE)
    } else note <- "Scope differs from previous snapshot: no revision comparison; new reference snapshot saved."
  }
  write_csv(changes, file.path(run_dir, "revision_changes.csv"))
  snapshot <- list(data = full[, c(keys, vars)], scope = scope, metadata = obj$metadata,
    quality_text = obj$quality_text, retrieved_utc = obj$retrieved_utc)
  saveRDS(snapshot, file.path(run_dir, "source_snapshot.rds"))
  saveRDS(snapshot, previous)
  note
}

hac_vcov <- function(model) {
  effective_n <- NROW(sandwich::estfun(model))
  sandwich::NeweyWest(model, lag = min(HAC_LAG, effective_n - 1L),
    prewhite = FALSE, adjust = TRUE)
}
fit_baseline <- function(d) {
  d <- d[order(d$month_index), ]
  stopifnot(nrow(d) == nrow(month_grid()), all(is.finite(d$rate_pct)))
  d$month_factor <- factor(d$month, levels = 1:12)
  train <- d[d$year <= BASELINE_END, ]
  if (nrow(train) != 12L * (BASELINE_END - START_YEAR + 1L)) stop("Incomplete baseline")
  fit <- stats::lm(rate_pct ~ month_index + month_factor, data = train)
  V <- hac_vcov(fit)
  X <- stats::model.matrix(~ month_index + month_factor, data = d)
  d$baseline_pct <- as.numeric(stats::predict(fit, newdata = d))
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  d$baseline_mean_low95 <- d$baseline_pct - stats::qnorm(0.975) * se
  d$baseline_mean_high95 <- d$baseline_pct + stats::qnorm(0.975) * se
  d$deviation_pp <- d$rate_pct - d$baseline_pct
  d$deviation_relative_pct <- ifelse(d$baseline_pct > 0, 100 * d$deviation_pp / d$baseline_pct, NA_real_)
  co <- stats::coef(fit); cs <- sqrt(diag(V))
  coefficients <- data.frame(term = names(co), estimate = unname(co), se_HAC = unname(cs),
    low95_HAC = unname(co - stats::qnorm(0.975) * cs), high95_HAC = unname(co + stats::qnorm(0.975) * cs))
  annual <- do.call(rbind, lapply(split(d, d$year), function(y) {
    b <- mean(y$baseline_pct); actual <- mean(y$rate_pct)
    data.frame(year = y$year[1], n_months = nrow(y), mean_monthly_rate_pct = actual,
      mean_monthly_baseline_pct = b, mean_deviation_pp = actual - b,
      relative_to_baseline_pct = if (b > 0) 100 * (actual - b) / b else NA_real_)
  }))
  summary <- data.frame(n_months = nrow(d), n_baseline = nrow(train),
    baseline_slope_pp_per_year = unname(co["month_index"]) * 12,
    slope_low95_HAC = coefficients$low95_HAC[coefficients$term == "month_index"] * 12,
    slope_high95_HAC = coefficients$high95_HAC[coefficients$term == "month_index"] * 12,
    baseline_LjungBox_lag12_p = stats::Box.test(stats::residuals(fit), lag = 12, type = "Ljung-Box")$p.value,
    baseline_outside_0_100 = any(d$baseline_pct < 0 | d$baseline_pct > 100))
  list(data = d, coefficients = coefficients, annual = annual, summary = summary,
       fit = fit, vcov_HAC = V)
}
fit_breaks <- function(d, h) {
  max_k <- min(MAX_BREAKS, ceiling(nrow(d) / h) - 2L)
  if (max_k < 1L) stop("Minimum segment length leaves no admissible break model")
  bp <- strucchange::breakpoints(deviation_pp ~ 1, data = d, h = h, breaks = max_k)
  candidates <- 0:max_k
  bic <- as.numeric(stats::AIC(bp, breaks = candidates, k = log(nrow(d))))
  if (any(!is.finite(bic))) stop("Non-finite BIC: degenerate series requires review")
  k <- candidates[which.min(bic)]
  chosen <- if (k) strucchange::breakpoints(bp, breaks = k)$breakpoints else integer()
  # strucchange indexes the LAST observation in the preceding segment.
  segment <- findInterval(seq_len(nrow(d)) - 1L, chosen) + 1L
  means <- tapply(d$deviation_pp, segment, mean)
  fitted <- unname(means[as.character(segment)])
  breaks <- data.frame(min_segment_months = integer(), break_number = integer(),
    last_month_before = character(), first_month_after = character(), shift_pp = numeric(),
    ci95_first_after_low = character(), ci95_first_after_high = character(), ci_note = character())
  warn <- character()
  if (k) {
    ci_note <- "Approximate HAC CI; conditional on baseline and selected break count"
    ci <- tryCatch(withCallingHandlers(stats::confint(bp, breaks = k, level = 0.95,
      het.reg = FALSE, het.err = TRUE, vcov. = hac_vcov)$confint,
      warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) { ci_note <<- paste("CI unavailable:", conditionMessage(e)); NULL })
    bounds <- matrix(NA_real_, k, 2)
    if (!is.null(ci)) bounds[,] <- ci[, c(1, 3), drop = FALSE]
    labels <- function(j) {
      out <- rep(NA_character_, length(j)); ok <- is.finite(j) & j >= 1 & j <= nrow(d)
      out[ok] <- d$month_label[as.integer(j[ok])]; out
    }
    breaks <- data.frame(min_segment_months = h, break_number = seq_len(k),
      last_month_before = d$month_label[chosen], first_month_after = d$month_label[chosen + 1L],
      shift_pp = diff(as.numeric(means)), ci95_first_after_low = labels(bounds[, 1] + 1L),
      ci95_first_after_high = labels(bounds[, 2] + 1L), ci_note = ci_note)
    bad <- is.na(breaks$ci95_first_after_low) | is.na(breaks$ci95_first_after_high)
    breaks$ci_note[bad] <- paste(breaks$ci_note[bad], "Bounds unavailable/outside observed window.")
    if (length(warn)) breaks$ci_note <- paste(breaks$ci_note, paste(unique(warn), collapse = "; "))
  }
  endpoints <- c(0L, chosen, nrow(d))
  segments <- do.call(rbind, lapply(seq_len(length(endpoints) - 1L), function(i) {
    ix <- seq.int(endpoints[i] + 1L, endpoints[i + 1L])
    data.frame(min_segment_months = h, segment = i, start_month = d$month_label[min(ix)],
      end_month = d$month_label[max(ix)], n_months = length(ix), mean_deviation_pp = mean(d$deviation_pp[ix]))
  }))
  summary <- data.frame(min_segment_months = h, selected_breaks = k,
    at_candidate_limit = k == max_k, BIC_gain_over_zero = bic[1] - min(bic),
    segmented_LjungBox_lag12_p = stats::Box.test(d$deviation_pp - fitted, lag = 12, type = "Ljung-Box")$p.value,
    first_months_after = paste(d$month_label[chosen + 1L], collapse = "; "))
  list(breaks = breaks, segments = segments, summary = summary,
    bic = data.frame(min_segment_months = h, n_breaks = candidates, BIC = bic,
      delta_BIC = bic - min(bic), selected = candidates == k),
    segment = segment, fitted = fitted, full_model = bp, selected_k = k)
}
fit_group <- function(d) {
  base <- fit_baseline(d)
  hs <- unique(c(MIN_SEGMENT_MONTHS, SENSITIVITY_MONTHS))
  fits <- lapply(hs, function(h) fit_breaks(base$data, h)); names(fits) <- as.character(hs)
  primary <- fits[[as.character(MIN_SEGMENT_MONTHS)]]
  base$data$segment <- primary$segment
  base$data$segment_mean_deviation_pp <- primary$fitted
  list(data = base$data, coefficients = base$coefficients, annual = base$annual,
    summary = cbind(base$summary, primary$summary),
    breaks = do.call(rbind, lapply(fits, `[[`, "breaks")),
    segments = do.call(rbind, lapply(fits, `[[`, "segments")),
    bic = do.call(rbind, lapply(fits, `[[`, "bic")),
    sensitivity = do.call(rbind, lapply(fits, `[[`, "summary")),
    model = list(baseline = base$fit, baseline_vcov_HAC = base$vcov_HAC, break_fits = fits))
}

plot_group <- function(r) {
  d <- r$data; x <- d$year + (d$month - 1) / 12
  graphics::par(mfrow = c(3, 1), mar = c(4.3, 4.2, 2.5, 1), mgp = c(2.4, 0.7, 0))
  graphics::plot(x, d$rate_pct, type = "n", xlab = "Year", ylab = "Absence (%)",
    main = paste(d$organisation[1], "-", d$staff_group[1]), cex.main = 0.9,
    ylim = range(c(d$rate_pct, d$baseline_mean_low95, d$baseline_mean_high95)))
  graphics::polygon(c(x, rev(x)), c(d$baseline_mean_low95, rev(d$baseline_mean_high95)),
    col = "#0072B222", border = NA)
  graphics::lines(x, d$baseline_pct, col = "#0072B2", lty = 2, lwd = 1.4)
  graphics::lines(x, d$rate_pct, col = "#222222", lwd = 1.3)
  graphics::abline(v = BASELINE_END + 1, col = "grey60", lty = 3)
  graphics::legend("topleft", c("Observed monthly rate", "2010-2018 baseline; 95% mean CI"),
    col = c("#222222", "#0072B2"), lty = c(1, 2), bty = "n", cex = 0.75)
  graphics::plot(x, d$deviation_pp, type = "l", xlab = "Year", ylab = "Deviation (pp)",
    main = paste("Level shifts: minimum", MIN_SEGMENT_MONTHS, "months; BIC selects", r$summary$selected_breaks, "breaks"),
    cex.main = 0.85)
  graphics::abline(h = 0, col = "grey60", lty = 3)
  for (s in unique(d$segment)) {
    ix <- which(d$segment == s)
    graphics::segments(x[min(ix)], d$segment_mean_deviation_pp[ix[1]],
      x[max(ix)], d$segment_mean_deviation_pp[ix[1]], col = "#D55E00", lwd = 2)
  }
  edge <- which(diff(d$segment) != 0)
  if (length(edge)) graphics::abline(v = x[edge + 1L], col = "#D55E00", lty = 2)
  hs <- sort(unique(r$sensitivity$min_segment_months))
  graphics::plot(range(x), c(0.5, length(hs) + 0.5), type = "n", xlab = "First month after estimated break",
    ylab = "", yaxt = "n", main = "Sensitivity to minimum segment length", cex.main = 0.85)
  graphics::axis(2, seq_along(hs), labels = paste(hs, "months"), las = 1, cex.axis = 0.75)
  for (i in seq_along(hs)) {
    b <- r$breaks[r$breaks$min_segment_months == hs[i], ]
    graphics::abline(h = i, col = "grey90")
    if (nrow(b)) graphics::points(x[match(b$first_month_after, d$month_label)], rep(i, nrow(b)), pch = 19, col = "#D55E00")
  }
}
html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE)
}
html_table <- function(x, digits = 3) {
  if (!nrow(x)) return("<p>None.</p>")
  for (n in names(x)) if (is.numeric(x[[n]])) x[[n]] <- format(round(x[[n]], digits), trim = TRUE)
  cells <- apply(x, 1, function(z) paste0("<tr><td>", paste(html_escape(z), collapse = "</td><td>"), "</td></tr>"))
  paste0("<div class='scroll'><table><thead><tr><th>", paste(html_escape(names(x)), collapse = "</th><th>"),
    "</th></tr></thead><tbody>", paste(cells, collapse = "\n"), "</tbody></table></div>")
}
# ---- England-style presentation (embedded here; no companion file needed) ---
# Uses the existing fitted baseline and independently selected breaks. It does
# not change the baseline fit, breakpoint search, or the other diagnostic plots.
style_roll6 <- function(x) {
  # Match the original England window exactly: t-3, t-2, t-1, t, t+1, t+2.
  n <- length(x); out <- rep(NA_real_,n)
  if(n>=6L) for(i in seq.int(4L,n-2L)) {
    z <- x[seq.int(i-3L,i+2L)]
    if(all(is.finite(z)))out[i] <- mean(z)
  }
  out
}
style_prepare <- function(d,baseline_fit,break_table,country,staff_group="All staff",min_segment_months=24L) {
  required <- c("date","rate_pct","baseline_pct","month_index","month_factor")
  if(!all(required %in% names(d)))stop("England-style plot: missing model data columns")
  d <- d[order(d$date),];d$date <- as.Date(d$date)
  expected <- seq(min(d$date),max(d$date),by="month")
  if(!identical(as.character(d$date),as.character(expected)))stop("England-style plot requires contiguous monthly data")
  # Original England convention: ordinary pointwise OLS prediction limits,
  # then smooth the endpoints. This is a visual reference envelope, NOT a
  # calibrated prediction interval for a six-month average or a HAC mean CI.
  pi <- stats::predict(baseline_fit,newdata=d,interval="prediction",level=.95)
  s <- data.frame(country=country,staff_group=staff_group,date=d$date,
    month_index=d$month_index,rate_100k=d$rate_pct*1000,baseline_100k=d$baseline_pct*1000,
    rolling_100k=style_roll6(d$rate_pct*1000),baseline_roll6_100k=style_roll6(d$baseline_pct*1000),
    pi_low_pointwise_100k=pi[,"lwr"]*1000,pi_high_pointwise_100k=pi[,"upr"]*1000,
    pi_low_roll6_100k=style_roll6(pi[,"lwr"]*1000),pi_high_roll6_100k=style_roll6(pi[,"upr"]*1000))
  s$deviation_roll6_relative_pct <- ifelse(s$baseline_roll6_100k>0,
    100*(s$rolling_100k-s$baseline_roll6_100k)/s$baseline_roll6_100k,NA_real_)
  n <- nrow(s)
  before <- c(NA_real_,head(s$rolling_100k,-1L));after <- c(tail(s$rolling_100k,-1L),NA_real_)
  is_high <- s$rolling_100k>before & s$rolling_100k>after
  is_low <- s$rolling_100k<before & s$rolling_100k<after
  ix <- which((is_high|is_low) & is.finite(s$deviation_roll6_relative_pct))
  extrema <- s[ix,c("date","rolling_100k","deviation_roll6_relative_pct"),drop=FALSE]
  extrema$kind <- ifelse(is_high[ix],"high","low")
  extrema$label <- sprintf("%+.1f%%",extrema$deviation_roll6_relative_pct)
  extrema$label[abs(extrema$deviation_roll6_relative_pct)<.05] <- "0.0%"
  b <- break_table
  if(nrow(b)) {
    b$date <- as.Date(paste0(b$first_month_after,"-01"))
    if(anyNA(b$date) || anyDuplicated(b$date) || any(!b$date %in% s$date))stop("Invalid plot break dates")
    b <- b[order(b$date),];b$break_number <- seq_len(nrow(b))
    b$label <- paste0("B",b$break_number," ",ifelse(b$shift_pp>0,"up","down")," | ",format(b$date,"%Y-%m"))
  } else {
    b <- data.frame(break_number=integer(),first_month_after=character(),last_month_before=character(),
      shift_pp=numeric(),date=as.Date(character()),label=character())
  }
  # Descriptive straight lines fitted to the rolling series between detected
  # level breaks; their slopes are not a separate formal change-in-slope test.
  lines <- data.frame(date=as.Date(character()),value_100k=numeric(),id=character(),kind=character(),colour=character())
  trends <- data.frame(id=character(),kind=character(),start=character(),end=character(),n_months=integer(),slope_100k_year=numeric())
  trend_fit <- function(id,kind,start,end,colour) {
    z <- s[s$date>=start & s$date<=end & is.finite(s$rolling_100k),]
    if(nrow(z)<3L)return(invisible(NULL))
    model <- stats::lm(rolling_100k~month_index,data=z)
    lines <<- rbind(lines,data.frame(date=z$date,value_100k=as.numeric(predict(model,newdata=z)),id=id,kind=kind,colour=colour))
    trends <<- rbind(trends,data.frame(id=id,kind=kind,start=format(min(z$date),"%Y-%m"),end=format(max(z$date),"%Y-%m"),
      n_months=nrow(z),slope_100k_year=unname(coef(model)["month_index"])*12))
  }
  segment_colours <- c("#7B3294","#009E73","#C9AE00","#0072B2","#D55E00","#8A6D3B")
  if(nrow(b))for(i in seq_len(nrow(b))) {
    start <- b$date[i]
    # First-after convention: the next segment begins at its break month.
    end <- if(i<nrow(b))s$date[match(b$date[i+1L],s$date)-1L] else max(s$date)
    id <- if(i<nrow(b))paste0("B",i,"-B",i+1L) else paste0("B",i,"-end")
    trend_fit(id,"segment",start,end,segment_colours[(i-1L)%%length(segment_colours)+1L])
  }
  # Generalise the original B1-to-end / B2-to-end overlays to the first two
  # post-baseline breaks, so early Welsh baseline-period breaks remain explicit.
  tails <- head(which(b$date>=as.Date(sprintf("%d-01-01",BASELINE_END+1L))),2L)
  if(length(tails))for(j in seq_along(tails)) {
    i <- tails[j]
    trend_fit(paste0("B",i,"-end overlay"),paste0("tail",j),b$date[i],max(s$date),c("#7B3294","#1B9E77")[j])
  }
  limits <- range(c(s$rate_100k,s$rolling_100k,s$baseline_roll6_100k,
    s$pi_low_roll6_100k,s$pi_high_roll6_100k,lines$value_100k),finite=TRUE)
  list(data=s,breaks=b,extrema=extrema,lines=lines,trends=trends,country=country,
    staff_group=staff_group,min_segment_months=min_segment_months,limits=limits,
    baseline_slope_100k_year=unname(coef(baseline_fit)["month_index"])*12000)
}
style_axis_limits <- function(objects) {
  y <- range(unlist(lapply(objects,`[[`,"limits")),finite=TRUE)
  span <- max(diff(y),500)
  c(max(0,y[1]-.12*span),y[2]+.24*span)
}
style_plot <- function(obj,x_limits=range(obj$data$date),y_limits=style_axis_limits(list(obj)),events=STYLE_EVENTS) {
  s <- obj$data;b <- obj$breaks;extrema <- obj$extrema
  start_projection <- as.Date(sprintf("%d-01-01",BASELINE_END+1L))
  span <- diff(y_limits)
  extrema$nudge <- ifelse(extrema$kind=="high",.025*span,-.025*span)
  if(nrow(b))b$y_label <- y_limits[1]+span*rep(c(.075,.025),length.out=nrow(b))
  events <- events[events$date>=x_limits[1] & events$date<=x_limits[2],,drop=FALSE]
  if(nrow(events))events$y_label <- y_limits[2]-span*rep(c(.065,.14),length.out=nrow(events))
  post <- s[s$date>=start_projection,]
  train <- s[s$date<start_projection,]
  tick_start <- as.Date(paste0(format(x_limits[1],"%Y"),"-01-01"))
  ticks <- seq(tick_start,x_limits[2],by="6 months")
  ticks <- ticks[ticks>=x_limits[1]]
  month_labels <- function(x) paste0(month.abb[as.integer(format(x,"%m"))],"\n",format(x,"%Y"))
  caption <- sprintf("%d breaks (BIC) | Minimum segment: %d months | Baseline slope: %+.1f per 100k/year | Breaks label first month after",
    nrow(b),obj$min_segment_months,obj$baseline_slope_100k_year)
  recent <- obj$trends[obj$trends$kind=="segment" & obj$trends$end>=sprintf("%d-01",BASELINE_END+1L),]
  if(nrow(recent))caption <- paste0(caption,"\nSegment slopes: ",paste(sprintf("%s %+.0f/year",recent$id,recent$slope_100k_year),collapse=" | "))
  p <- ggplot2::ggplot(s,ggplot2::aes(x=date))+
    ggplot2::geom_vline(xintercept=as.numeric(start_projection),colour="grey65",linewidth=.45,linetype="dotted")+
    ggplot2::geom_ribbon(data=post,ggplot2::aes(ymin=pi_low_roll6_100k,ymax=pi_high_roll6_100k),
      fill="#0072B2",alpha=.12,na.rm=TRUE)+
    ggplot2::geom_line(ggplot2::aes(y=rate_100k),colour="grey55",linewidth=.65,alpha=.55)+
    ggplot2::geom_line(ggplot2::aes(y=rolling_100k),colour="#E69F00",linewidth=1.35,na.rm=TRUE)+
    ggplot2::geom_line(data=train,ggplot2::aes(y=baseline_roll6_100k),colour="#0072B2",linewidth=1.1,na.rm=TRUE)+
    ggplot2::geom_line(data=post,ggplot2::aes(y=baseline_roll6_100k),colour="#0072B2",linewidth=1.1,linetype="dashed",na.rm=TRUE)
  if(nrow(obj$lines))for(id in unique(obj$lines$id)) {
    z <- obj$lines[obj$lines$id==id,]
    lt <- switch(z$kind[1],segment="solid",tail1="dashed",tail2="dotdash")
    p <- p+ggplot2::geom_line(data=z,ggplot2::aes(x=date,y=value_100k),inherit.aes=FALSE,
      colour=z$colour[1],linetype=lt,linewidth=if(z$kind[1]=="segment")1.1 else .9)
  }
  if(nrow(extrema))p <- p+ggrepel::geom_label_repel(data=extrema,
    ggplot2::aes(x=date,y=rolling_100k,label=label),inherit.aes=FALSE,
    colour="#666666",fill="#FFFFFFE6",size=2.65,fontface="bold",label.size=.18,
    label.padding=grid::unit(.10,"lines"),box.padding=grid::unit(.18,"lines"),point.padding=grid::unit(.15,"lines"),
    segment.color="#66666670",segment.size=.25,direction="y",nudge_y=extrema$nudge,
    min.segment.length=0,max.overlaps=Inf,seed=42,max.iter=20000,max.time=3,
    xlim=as.numeric(x_limits),ylim=y_limits)
  if(nrow(events))p <- p+
    ggplot2::geom_vline(data=events,ggplot2::aes(xintercept=as.numeric(date)),
      colour="#D55E00",linetype="dotted",linewidth=.7,alpha=.7)+
    ggplot2::geom_label(data=events,ggplot2::aes(x=date,y=y_label,label=label),inherit.aes=FALSE,
      colour="#D55E00",fill="white",size=2.8,fontface="bold",label.padding=grid::unit(.12,"lines"),label.size=.2)
  if(nrow(b))p <- p+
    ggplot2::geom_vline(data=b,ggplot2::aes(xintercept=as.numeric(date)),
      colour="#CC79A7",linetype="longdash",linewidth=.85,alpha=.85)+
    ggplot2::geom_label(data=b,ggplot2::aes(x=date,y=y_label,label=label),inherit.aes=FALSE,
      colour="#B75C99",fill="white",size=2.45,fontface="bold",label.padding=grid::unit(.10,"lines"),label.size=.18)
  p+ggplot2::scale_x_date(breaks=ticks,labels=month_labels,expand=ggplot2::expansion(mult=0))+
    ggplot2::scale_y_continuous(labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE),expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_cartesian(xlim=x_limits,ylim=y_limits,clip="off")+
    ggplot2::labs(title=paste0("NHS Sickness Absence Rate (",obj$country,"): Breaks vs 2010-2018 Norm"),
      subtitle=paste0(obj$staff_group," | Rolling (6m), seasonal baseline and smoothed pointwise OLS prediction band; breaks on raw deviations"),
      x=NULL,y="Absence days per 100,000 available days",caption=caption)+
    ggplot2::theme_minimal(base_size=12)+
    ggplot2::theme(plot.title=ggplot2::element_text(size=17,face="bold",margin=ggplot2::margin(b=6)),
      plot.subtitle=ggplot2::element_text(size=10,colour="grey35",margin=ggplot2::margin(b=10)),
      plot.caption=ggplot2::element_text(size=9,colour="grey40",hjust=0,margin=ggplot2::margin(t=10)),
      axis.text.x=ggplot2::element_text(size=8),axis.text.y=ggplot2::element_text(size=10),
      axis.title.y=ggplot2::element_text(size=11,margin=ggplot2::margin(r=10)),
      panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),
      panel.grid.major.y=ggplot2::element_line(colour="grey88",linewidth=.35),
      plot.margin=ggplot2::margin(t=14,r=32,b=14,l=14),legend.position="none")
}
style_save_plot <- function(plot,png_path,pdf_path=NULL) {
  ggplot2::ggsave(png_path,plot=plot,width=18,height=10,dpi=220,bg="white",limitsize=FALSE)
  if(!is.null(pdf_path))ggplot2::ggsave(pdf_path,plot=plot,width=18,height=10,bg="white",limitsize=FALSE)
}
style_embed <- function(path,alt)paste0("<img alt='",html_escape(alt),"' src='data:image/png;base64,",
  jsonlite::base64_enc(readBin(path,"raw",n=file.info(path)$size)),"'>")
style_notes <- function()c(
  "England-style figures reproduce the original visual layers: raw monthly rates, a six-month average (t-3 through t+2), projected seasonal baseline, smoothed prediction envelope, local relative-deviation labels, independent break markers and descriptive segment/tail trends.",
  "The blue envelope is the original-style six-month smoothing of ordinary pointwise OLS 95% prediction limits. It is not a calibrated prediction interval for a rolling mean and does not replace the HAC baseline mean intervals in the other tables/diagnostic figures.",
  "Grey percentage labels compare the six-month observed average with the six-month baseline. Straight segment/tail lines are descriptive fits to overlapping rolling means, not separately tested trend breaks. A rolling value can span a detected break.",
  "Break dates label the first month in the new segment; the original England figure labelled the previous month. The existing models and break-date tables remain authoritative.",
  "Orange event markers are editable UK reference dates, not fitted effects or country-specific policy dates. English-only mandate annotations are not transferred to Wales.")
# ---- End England-style presentation ----------------------------------------

write_report <- function(results, quality, notes, run_dir) {
  parts <- c("<!doctype html><html lang='en'><meta charset='UTF-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
    "<title>Wales NHS sickness absence</title><style>body{font:16px/1.5 system-ui,sans-serif;max-width:1200px;margin:35px auto;padding:0 20px;color:#18232d}h1,h2{color:#174d70}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #d7dfe5;padding:7px;text-align:left}th{background:#edf4f8}.scroll{overflow:auto}img{width:100%;height:auto}li{margin:.4em 0}.note{background:#f1f6f9;padding:15px}a{color:#125d8c}</style>",
    "<h1>Wales NHS sickness absence, 2010-2025</h1>",
    paste0("<p>Run completed: ", utc_now(), ". Rates are percentages; deviations are percentage points (pp).</p>"),
    "<h2>Read before interpreting results</h2><div class='note'><ul>",
    paste0("<li>", html_escape(notes), "</li>"), "</ul></div>",
    paste0("<p>Sources: <a href='", DATASET_URL, "'>StatsWales dataset</a>; <a href='", QUALITY_URL, "'>quality report</a>.</p>"),
    "<p>Baseline: separate linear trend plus month-of-year effects fitted to 2010-2018 only. Breaks: BIC-selected level shifts in unsmoothed deviations, searched over 2010-2025. English dates are never imposed.</p>",
    "<h2>Data completeness</h2>", html_table(quality[, setdiff(names(quality), "missing_month_labels")]),
    "<p>See quality_checks.csv for missing months, known_definition_changes.csv for comparability flags, and revision_changes.csv for changes since the previous snapshot.</p>")
  for (i in seq_along(results)) {
    r <- results[[i]]
    parts <- c(parts, paste0("<h2>", html_escape(paste(r$data$organisation[1], "-", r$data$staff_group[1])), "</h2>"),
      "<h3>Baseline and primary model</h3>", html_table(r$summary),
      "<h3>Breakpoint sensitivity</h3>", html_table(r$sensitivity),
      "<p>Dates below identify the first month in the new segment; last_month_before is also supplied for comparison with other software conventions.</p>",
      html_table(r$breaks), "<h3>Recent yearly summaries</h3>", html_table(r$annual[r$annual$year >= 2018, ]),
      "<p>Yearly summaries are unweighted means of the 12 monthly rates, not the publisher's annual ratios. Staff denominators are not supplied here.</p>")
    if (MAKE_ENGLAND_STYLE_PLOTS) parts <- c(parts,
      "<h3>England-style analysis plot</h3>",
      "<p>Grey: monthly observations; orange: six-month average; blue: fitted/projected baseline and smoothed OLS prediction envelope; pink: independently estimated breaks. Relative-deviation labels and coloured descriptive trend lines follow the original England presentation.</p>",
      style_embed(file.path(run_dir, "plots", sprintf("england_style_%02d.png", i)),
        paste("England-style sickness absence analysis:", r$data$staff_group[1])))
    if (MAKE_PLOTS) parts <- c(parts, "<h3>Existing diagnostic figures</h3>",
      style_embed(file.path(run_dir, "plots", sprintf("series_%02d.png", i)),
        "Monthly rates, HAC baseline mean intervals, deviations and breakpoint sensitivity"))
  }
  writeLines(c(parts, "</html>"), file.path(run_dir, "report.html"), useBytes = TRUE)
}

main <- function() {
  need_packages()
  stopifnot(SOURCE_MODE %in% c("refresh", "cache", "manual"), START_YEAR == 2010L,
    END_YEAR == 2025L, BASELINE_END == 2018L, MIN_SEGMENT_MONTHS >= 12L,
    all(SENSITIVITY_MONTHS >= 12L), MAX_BREAKS >= 1L)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(OUT_DIR, "cache"); dir.create(cache_dir, showWarnings = FALSE)
  run_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC"), "_"), tmpdir = OUT_DIR)
  dir.create(run_dir); raw_dir <- file.path(run_dir, "source"); dir.create(raw_dir)
  message("Downloading/auditing Wales data (", SOURCE_MODE, ") ...")
  obj <- tryCatch(load_source(raw_dir, cache_dir), error = function(e) {
    writeLines(conditionMessage(e), file.path(run_dir, "DOWNLOAD_FAILED.txt"))
    stop(conditionMessage(e), call. = FALSE)
  })
  a <- audit_source(obj)
  write_csv(a$data, file.path(run_dir, "monthly_observations_2010_2025.csv"))
  write_csv(a$quality, file.path(run_dir, "quality_checks.csv"))
  write_csv(known_changes(), file.path(run_dir, "known_definition_changes.csv"))
  rev_note <- revision_audit(a$data, obj, cache_dir, run_dir)
  if (SOURCE_MODE != "cache") saveRDS(obj, file.path(cache_dir, "source.rds"))
  notes <- c(paste("Source mode:", SOURCE_MODE, "| original import:", obj$retrieved_utc,
    "| original mode:", obj$origin_mode, "| source CSV MD5:", obj$csv_md5),
    rev_note, obj$metadata_issues,
    "Complete monthly series only. Missing observations are not interpolated, and annual/quarterly rows are excluded.",
    "Unformatted source values are used. The published display rounds to 0.1 percentage point; this is not measurement accuracy.",
    "No provisional flag does not prove a figure is final or will never be revised.",
    "Staff reclassifications and changes in employment coverage require review; see known_definition_changes.csv and archived source notes.",
    "Existing diagnostic baseline bands are 95% HAC mean intervals. The additional England-style figures use a separate smoothed OLS prediction envelope; see the plot-specific notes. Long extrapolation depends on the linear-trend assumption.",
    "Compare 12-, 18- and 24-month segment constraints before interpreting break timing. A constraint can shift estimated breaks.",
    "Break intervals are approximate and conditional. BIC may overselect under autocorrelation; Ljung-Box values are descriptive diagnostics.",
    "The all-staff total overlaps with its components. Absence changes do not establish pathogen-specific or causal explanations.")
  writeLines(notes, file.path(run_dir, "READ_BEFORE_ANALYSIS.txt"), useBytes = TRUE)
  write_csv(data.frame(source_mode = SOURCE_MODE, original_import_utc = obj$retrieved_utc,
    origin_mode = obj$origin_mode, csv_md5 = obj$csv_md5, source_url = obj$source_url,
    metadata_available = nzchar(obj$metadata), quality_report_available = nzchar(obj$quality_text)),
    file.path(run_dir, "source_status.csv"))
  print(a$quality[, c("staff_group", "numeric_months", "missing_months", "eligible")], row.names = FALSE)
  eligible <- a$quality[a$quality$eligible, ]
  if (!nrow(eligible)) stop("No complete eligible series. Audit saved in ", run_dir, call. = FALSE)
  results <- list(); errors <- data.frame(organisation_id = character(), staff_id = character(), error = character())
  # Put all-staff total first in the report; do not infer it by adding group rates.
  eligible <- eligible[order(eligible$staff_group != "All staff", eligible$staff_group), ]
  for (i in seq_len(nrow(eligible))) {
    id <- eligible[i, ]; message("Fitting ", id$organisation, ": ", id$staff_group, " ...")
    d <- a$data[a$data$organisation_id == id$organisation_id & a$data$staff_id == id$staff_id, ]
    result <- tryCatch(fit_group(d), error = function(e) {
      errors <<- rbind(errors, data.frame(organisation_id = id$organisation_id, staff_id = id$staff_id,
        error = conditionMessage(e))); NULL
    })
    if (!is.null(result)) results[[paste(id$organisation_id, id$staff_id, sep = "_")]] <- result
  }
  write_csv(errors, file.path(run_dir, "model_errors.csv"))
  if (!length(results)) stop("No model succeeded. Inspect model_errors.csv in ", run_dir, call. = FALSE)
  combine <- function(name) do.call(rbind, lapply(results, function(r) {
    x <- r[[name]]
    if (name == "data") return(x)
    identity <- r$data[rep(1L, nrow(x)), c("organisation_id", "organisation", "staff_id", "staff_group"), drop = FALSE]
    cbind(identity, x)
  }))
  filenames <- c(data = "baseline_and_deviations.csv", coefficients = "baseline_coefficients.csv",
    annual = "annual_mean_monthly_rates.csv", summary = "staff_group_summary.csv",
    breaks = "breakpoints_all_sensitivities.csv", segments = "segments_all_sensitivities.csv",
    bic = "breakpoint_BIC.csv", sensitivity = "breakpoint_sensitivity.csv")
  for (name in names(filenames)) write_csv(combine(name), file.path(run_dir, filenames[[name]]))
  primary <- combine("breaks"); primary <- primary[primary$min_segment_months == MIN_SEGMENT_MONTHS, ]
  write_csv(primary, file.path(run_dir, "breakpoints_primary.csv"))
  saveRDS(lapply(results, `[[`, "model"), file.path(run_dir, "fitted_models.rds"))
  if (MAKE_PLOTS) {
    plot_dir <- file.path(run_dir, "plots"); dir.create(plot_dir)
    grDevices::pdf(file.path(run_dir, "staff_group_plots.pdf"), width = 11, height = 10, onefile = TRUE)
    tryCatch(invisible(lapply(results, plot_group)), finally = grDevices::dev.off())
    for (i in seq_along(results)) {
      grDevices::png(file.path(plot_dir, sprintf("series_%02d.png", i)), width = 1650, height = 1500, res = 150)
      tryCatch(plot_group(results[[i]]), finally = grDevices::dev.off())
    }
  }
  if (MAKE_ENGLAND_STYLE_PLOTS) {
    plot_dir <- file.path(run_dir, "plots"); dir.create(plot_dir, showWarnings = FALSE)
    style_objects <- style_plots <- list()
    for (i in seq_along(results)) {
      r <- results[[i]]
      b <- r$breaks[r$breaks$min_segment_months == MIN_SEGMENT_MONTHS, ]
      obj <- style_prepare(r$data, r$model$baseline, b, r$data$organisation[1], r$data$staff_group[1], MIN_SEGMENT_MONTHS)
      style_objects[[i]] <- obj; style_plots[[i]] <- style_plot(obj)
      style_save_plot(style_plots[[i]], file.path(plot_dir, sprintf("england_style_%02d.png", i)))
      if (r$data$organisation_id[1] == "W" && r$data$staff_id[1] == "XXX") {
        if (!file.copy(file.path(plot_dir, sprintf("england_style_%02d.png", i)),
          file.path(plot_dir, "wales_sickness_analysis_england_style.png"))) stop("Could not save all-staff plot")
      }
    }
    write_csv(do.call(rbind, lapply(style_objects, `[[`, "data")), file.path(run_dir, "england_style_plot_series.csv"))
    trend_tables <- lapply(seq_along(style_objects), function(i) {
      x <- style_objects[[i]]$trends
      cbind(country = rep(style_objects[[i]]$country, nrow(x)),
        staff_group = rep(style_objects[[i]]$staff_group, nrow(x)), x)
    })
    write_csv(do.call(rbind, trend_tables), file.path(run_dir, "england_style_descriptive_trends.csv"))
    write_csv(STYLE_EVENTS, file.path(run_dir, "england_style_event_annotations.csv"))
    grDevices::pdf(file.path(run_dir, "wales_england_style_plots.pdf"), width = 18, height = 10, onefile = TRUE)
    tryCatch(invisible(lapply(style_plots, print)), finally = grDevices::dev.off())
    notes <- c(notes, style_notes())
  }
  if (nrow(errors)) notes <- c(notes, "INCOMPLETE: some models failed; see model_errors.csv.")
  if (any(a$quality$selected & !a$quality$eligible)) notes <- c(notes, "INCOMPLETE: some selected series have gaps/invalid values and were excluded.")
  if (any(combine("sensitivity")$at_candidate_limit)) notes <- c(notes,
    paste0("Some fits reached the search limit of ", MAX_BREAKS,
      " breaks (at_candidate_limit = TRUE). These counts are capped; consider increasing MAX_BREAKS before interpreting those sensitivity fits."))
  if (anyNA(combine("breaks")$ci95_first_after_low) || anyNA(combine("breaks")$ci95_first_after_high))
    notes <- c(notes, "Some break-date confidence bounds are unavailable. Missing bounds are left blank and explained in ci_note.")
  writeLines(notes, file.path(run_dir, "READ_BEFORE_ANALYSIS.txt"), useBytes = TRUE)
  write_report(results, a$quality, notes, run_dir)
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  dput(list(SOURCE_MODE = SOURCE_MODE, ORGANISATION_IDS = ORGANISATION_IDS, STAFF_IDS = STAFF_IDS,
    INCLUDE_UNCODED = INCLUDE_UNCODED, START_YEAR = START_YEAR, END_YEAR = END_YEAR,
    BASELINE_END = BASELINE_END, MIN_SEGMENT_MONTHS = MIN_SEGMENT_MONTHS,
    SENSITIVITY_MONTHS = SENSITIVITY_MONTHS, MAX_BREAKS = MAX_BREAKS, HAC_LAG = HAC_LAG,
    MAKE_ENGLAND_STYLE_PLOTS = MAKE_ENGLAND_STYLE_PLOTS, STYLE_EVENTS = STYLE_EVENTS),
    file = file.path(run_dir, "configuration.R"))
  message("Saved ", length(results), " staff-group analyses in: ", normalizePath(run_dir, winslash = "/"),
    "\nOpen report.html. Primary break dates: breakpoints_primary.csv; sensitivities: breakpoint_sensitivity.csv.")
  if (nrow(errors) || any(a$quality$selected & !a$quality$eligible))
    warning("Incomplete requested analysis: inspect quality_checks.csv and model_errors.csv.", call. = FALSE)
  if (length(obj$metadata_issues)) warning("Source metadata was not fully verified; read source_status.csv and READ_BEFORE_ANALYSIS.txt.", call. = FALSE)
  invisible(list(run_dir = run_dir, quality = a$quality, results = results))
}
if (isTRUE(RUN_ANALYSIS)) analysis <- main()
