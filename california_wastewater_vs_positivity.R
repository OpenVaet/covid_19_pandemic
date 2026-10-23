# California wastewater RNA versus county COVID-19 NAAT test positivity
# San Francisco (06075), Santa Clara (06085), Sacramento (06067)
# Run in R/RStudio: source("california_wastewater_vs_positivity.R")
# Install once, in a clean R session if needed:
# install.packages(c("curl", "jsonlite"), repos = "https://cloud.r-project.org")
# No httr2, automatic package installation, R Markdown or Pandoc is required.
# The first download is approximately 185 MiB. Cached sources are reused.
#
# DESIGN, fixed before calculating correlations:
# * Historical CDC county NAAT positivity through 2023-05-10; no case-rate proxy.
# * Thursday-Wednesday weeks, >=2 distinct wastewater sampling days per week.
# * Separate plants, laboratories, N-gene targets and reported method metadata.
# * Primary concentration uses LOD/sqrt(2) for explicitly flagged nondetects.
# * Sensitivities: exclude nondetects; PMMoV normalization on identical weeks;
#   calendar-year estimates; 1/2/3 sampling-day thresholds; 2/4/8-week blocks.
# * Primary window matches the dates of the England national wastewater series.
#   Full-overlap results are additional; neither is a direct test of superiority
#   over England, whose geography, sampling and clinical testing differ.
# * Lag scans are exploratory retrospective associations, NOT forecast validation.

RUN_ANALYSIS <- TRUE
OUTPUT_ROOT <- "california_wastewater_results"
REFRESH_DOWNLOADS <- FALSE
LOCAL_WASTEWATER_CSV <- "" # Optional original California sample-level CSV
LOCAL_POSITIVITY_CSV <- "" # Optional CDC CSV, or county weekly export described below
LOCAL_POSITIVITY_SOURCE <- "" # Required description/URL for a local positivity CSV
# Local exports: county_fips, week_end, positivity_pct; OR your map script's
# geography_level, geography_id, week_end, positivity_pct. Values must be observed
# CDC county NAAT percentages for Thursday-Wednesday windows, not interpolated.
MIN_SAMPLING_DAYS <- 2L
MIN_PAIRED_WEEKS <- 20L
BOOTSTRAP_REPLICATES <- 1000L
BOOTSTRAP_BLOCK_WEEKS <- 4L
RANDOM_SEED <- 20261023L
LAG_WEEKS <- -3L:3L # Positive k: wastewater(t) versus positivity(t + k weeks)
ANALYSIS_START <- as.Date("2020-01-22")
ANALYSIS_END <- as.Date("2023-05-10")
COMMON_START <- as.Date("2021-06-01")
COMMON_END <- as.Date("2022-03-30")

COUNTIES <- c(`06075` = "San Francisco", `06085` = "Santa Clara", `06067` = "Sacramento")
CA_PACKAGE <- "b8c6ee3b-539d-4d62-8fa2-c7cd17c16656"
CA_LANDING <- "https://data.ca.gov/dataset/covid-19-wastewater-surveillance-data-california"
CA_API <- paste0("https://publish.data.ca.gov/api/3/action/package_show?id=", CA_PACKAGE)
CA_RESOURCE <- "16bb2698-c243-4b66-a6e8-4861ee66f8bf"
WW_URL <- paste0("https://publish.data.ca.gov/dataset/", CA_PACKAGE,
  "/resource/", CA_RESOURCE, "/download/master-covid-public.csv")
DICTIONARY_URL <- paste0("https://publish.data.ca.gov/dataset/", CA_PACKAGE,
  "/resource/e5b302ff-2a51-4033-9341-a6a1205075c5/",
  "download/nwss-data-dictionary_v3_1_0_20211208_metadata.csv")
CDC_LANDING <- paste0("https://data.cdc.gov/Public-Health-Surveillance/",
  "Weekly-COVID-19-County-Level-of-Community-Transmis/jgk8-6dpn")
CDC_API <- "https://data.cdc.gov/resource/jgk8-6dpn.csv"
CDC_METADATA <- "https://data.cdc.gov/api/views/jgk8-6dpn.json"

fail <- function(...) stop(..., call. = FALSE)
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
clean <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
lc <- function(x) tolower(clean(x))
yes <- function(x) lc(x) %in% c("yes", "true", "1")
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
num <- function(x, label) {
  s <- clean(x); out <- suppressWarnings(as.numeric(s))
  if (any(s != "" & !is.finite(out))) fail("Unrecognised numeric value in ", label)
  out[s == ""] <- NA_real_; out
}
mean_ok <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
day_end <- function(d) d + (3L - as.POSIXlt(d, tz = "UTC")$wday) %% 7L
weekly_diff <- function(x, d) x - x[match(d - 7L, d)]
md5 <- function(path) unname(tools::md5sum(path))
need_packages <- function() {
  p <- c("curl", "jsonlite")
  missing <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) fail("Missing packages: ", paste(missing, collapse = ", "),
    ". Restart R and install them with install.packages(), then run again.")
}

# HTTP status checked before saving. A failed refresh never replaces a good cache.
fetch <- function(url, destination, validate = function(p) TRUE) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  part <- tempfile("pending_", tmpdir = dirname(destination))
  on.exit(unlink(part), add = TRUE)
  for (attempt in 1:3) {
    h <- curl::new_handle(timeout = 600, connecttimeout = 25, followlocation = TRUE,
      useragent = "California-wastewater-research/1.0", httpheader = c(Accept = "*/*"))
    res <- tryCatch(curl::curl_fetch_disk(url, part, handle = h), error = identity)
    code <- if (inherits(res, "error")) NA_integer_ else res$status_code
    if (!is.na(code) && code >= 200 && code < 300 && file.info(part)$size > 0) {
      validate(part)
      if (!file.copy(part, destination, overwrite = TRUE)) fail("Cannot save ", destination)
      return(invisible(destination))
    }
    why <- if (inherits(res, "error")) conditionMessage(res) else paste("HTTP", code)
    retry <- is.na(code) || code %in% c(408, 429) || code >= 500
    if (!retry || attempt == 3) fail("Download failed: ", why, "\nURL: ", url,
      "\nNo substitute data used. Supply the original local CSV in the settings if necessary.")
    message("Attempt ", attempt, " failed: ", why, "; retrying.")
    Sys.sleep(2 ^ attempt)
  }
}
read_text_csv <- function(p) read.csv(p, colClasses = "character", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "", fileEncoding = "UTF-8-BOM")
check_cols <- function(x, fields) {
  if (!all(fields %in% names(x))) fail("Source schema changed; missing columns: ",
    paste(setdiff(fields, names(x)), collapse = ", "))
}
csv_validator <- function(fields) function(p) {
  h <- read.csv(p, nrows = 0L, check.names = FALSE); check_cols(h, fields); TRUE
}
json_validator <- function(p) { jsonlite::fromJSON(p); TRUE }
cached_source <- function(url, file, refresh, local = "", validate = function(p) TRUE) {
  if (nzchar(local)) {
    if (!file.exists(local)) fail("Local input not found: ", local)
    validate(local)
    return(list(path = normalizePath(local, winslash = "/"), url = "local input",
      acquired_utc = NA_character_, inspected_utc = utc_now(), md5 = md5(local),
      bytes = unname(file.info(local)$size), mode = "local"))
  }
  manifest <- paste0(file, ".source.json")
  if (refresh || !file.exists(file) || !file.exists(manifest)) {
    message("Downloading ", basename(file), "...")
    fetch(url, file, validate)
    m <- list(url = url, acquired_utc = utc_now(), md5 = md5(file), bytes = unname(file.info(file)$size))
    jsonlite::write_json(m, manifest, pretty = TRUE, auto_unbox = TRUE)
  }
  m <- jsonlite::fromJSON(manifest)
  if (!identical(m$url, url) || !identical(m$md5, md5(file)))
    fail("Cache provenance/hash mismatch: ", file, ". Refresh downloads or supply a local input.")
  validate(file)
  c(list(path = normalizePath(file, winslash = "/"), mode = "download/cache"), m)
}
query_url <- function(parameters) paste0(CDC_API, "?", paste(
  paste0(names(parameters), "=", vapply(parameters, utils::URLencode, character(1), reserved = TRUE)),
  collapse = "&"))

validate_positivity <- function(x) {
  if (all(c("fips_code", "date", "percent_test_results_reported") %in% names(x))) {
    if ("state_name" %in% names(x) && any(x$state_name != "California")) fail("Non-California CDC input.")
    fips <- x$fips_code; date <- substr(x$date, 1L, 10L); value <- x$percent_test_results_reported
  } else {
    if (all(c("geography_level", "geography_id") %in% names(x))) {
      x <- x[lc(x$geography_level) == "county", , drop = FALSE]
      x$county_fips <- x$geography_id
    }
    check_cols(x, c("county_fips", "week_end", "positivity_pct"))
    if ("positivity_status" %in% names(x)) {
      # The supplied map script uses reported / CDC missing or suppressed.
      # Reject non-observed substitute values.
      bad <- !is.na(x$positivity_pct) & clean(x$positivity_pct) != "" &
        !lc(x$positivity_status) %in% c("reported", "observed", "available")
      if (any(bad)) fail("Local positivity contains non-observed values; use the original CDC CSV.")
    }
    fips <- x$county_fips; date <- x$week_end; value <- x$positivity_pct
  }
  fips <- clean(fips)
  fips[grepl("^[0-9]{4}$", fips)] <- paste0("0", fips[grepl("^[0-9]{4}$", fips)])
  keep <- fips %in% names(COUNTIES)
  d <- as.Date(date[keep], format = "%Y-%m-%d")
  p <- num(value[keep], "positivity")
  out <- data.frame(county_fips = fips[keep], county = unname(COUNTIES[fips[keep]]),
    week_end = d, positivity_pct = p)
  if (anyNA(d) || any(day_end(d) != d) || anyDuplicated(out[c("county_fips", "week_end")]))
    fail("Positivity requires unique county/Wednesday week-end dates.")
  if (any(p < 0 | p > 100, na.rm = TRUE)) fail("Positivity must be percent, 0 to 100.")
  if (!all(names(COUNTIES) %in% out$county_fips)) fail("Positivity is missing a requested county.")
  out <- out[out$week_end >= ANALYSIS_START & out$week_end <= ANALYSIS_END, ]
  # Complete calendar, preserving missingness; no interpolation or zero imputation.
  grid <- expand.grid(county_fips = names(COUNTIES), week_end = seq(day_end(ANALYSIS_START),
    ANALYSIS_END, by = "7 days"), stringsAsFactors = FALSE)
  grid$county <- unname(COUNTIES[grid$county_fips])
  out <- merge(grid, out[c("county_fips", "week_end", "positivity_pct")], all.x = TRUE)
  out[order(out$county_fips, out$week_end), ]
}

WW_COLUMNS <- c("wwtp_name", "county_names", "population_served", "sample_id", "lab_id",
  "sample_collect_date", "sample_type", "sample_matrix", "sample_location",
  "concentration_method", "extraction_method", "pcr_type", "lod_sewage",
  "hum_frac_target_mic", "hum_frac_mic_conc", "hum_frac_mic_unit", "quality_flag",
  "qc_ignore", "dashboard_ignore", "analysis_ignore", "major_lab_method",
  "major_lab_method_desc", "pcr_target", "pcr_gene_target", "pcr_target_units",
  "pcr_target_avg_conc", "pcr_target_below_lod", "data_source", "site_id",
  "created_at", "county_treatmentplant", "label_name")
read_wastewater <- function(path) {
  h <- names(read.csv(path, nrows = 0L, check.names = FALSE))
  check_cols(setNames(as.data.frame(matrix(nrow = 0, ncol = length(h))), h), WW_COLUMNS)
  classes <- ifelse(h %in% WW_COLUMNS, "character", "NULL")
  x <- read.csv(path, colClasses = classes, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = "", fileEncoding = "UTF-8-BOM")
  x$source_row <- seq_len(nrow(x))
  x[] <- lapply(x, clean)
  plant_county <- sub(" county$", "", lc(x$county_treatmentplant))
  x <- x[plant_county %in% tolower(unname(COUNTIES)), , drop = FALSE]
  x$county <- unname(COUNTIES[match(sub(" county$", "", lc(x$county_treatmentplant)),
    tolower(unname(COUNTIES)))])
  x$county_fips <- names(COUNTIES)[match(x$county, COUNTIES)]
  x$date <- as.Date(x$sample_collect_date, format = "%Y-%m-%d")
  if (anyNA(x$date)) fail("Invalid wastewater sampling date.")
  x[x$date >= ANALYSIS_START & x$date <= ANALYSIS_END, , drop = FALSE]
}
linear_value <- function(value, units) {
  out <- rep(NA_real_, length(value)); u <- lc(units)
  direct <- u %in% c("copies/g dry sludge", "copies/l wastewater")
  out[direct] <- value[direct]
  logunits <- u == "log10 copies/l wastewater"
  out[logunits] <- 10 ^ value[logunits]
  out[!is.finite(out) | out <= 0] <- NA_real_; out
}
base_units <- function(x) sub("^log10 ", "", lc(x))

prepare_samples <- function(x) {
  # Disposition is hierarchical; every input row remains in the audit export.
  x$disposition <- "eligible"
  mark <- function(test, label) {
    test[is.na(test)] <- TRUE
    x$disposition[x$disposition == "eligible" & test] <<- label
  }
  mark(!lc(x$pcr_target) %in% c("sars-cov-2", "sars_cov_2", "sars-cov2") |
    !lc(x$pcr_gene_target) %in% c("n", "n1", "n2"), "other pathogen/target")
  # Upstream sites and metadata-free samples cannot safely be assigned to a plant.
  mark(lc(x$sample_location) != "wwtp", "not confirmed treatment-plant sample")
  metadata_fields <- c("wwtp_name", "lab_id", "sample_id", "sample_matrix",
    "sample_type", "pcr_target_units", "concentration_method", "extraction_method", "pcr_type")
  mark(Reduce(`|`, lapply(x[metadata_fields], function(z) clean(z) == "")), "incomplete method metadata")
  mark(!lc(x$pcr_target_units) %in% c("copies/g dry sludge", "copies/l wastewater",
    "log10 copies/l wastewater"), "unsupported concentration units")
  is_solids <- lc(x$sample_matrix) %in% c("primary sludge", "sludge")
  is_liquid <- lc(x$sample_matrix) == "raw wastewater"
  mark(!(is_solids & base_units(x$pcr_target_units) == "copies/g dry sludge" |
    is_liquid & base_units(x$pcr_target_units) == "copies/l wastewater"), "matrix/units mismatch")
  sample_key <- paste(x$wwtp_name, x$lab_id, x$sample_id, x$date, lc(x$pcr_gene_target), sep = "|")
  # Propagate an explicit quality/analysis exclusion across duplicate exports of a sample.
  # qc_ignore only removes entries from the QC dashboard; it is NOT a quality flag.
  bad_keys <- unique(sample_key[yes(x$quality_flag) | yes(x$analysis_ignore)])
  mark(sample_key %in% bad_keys, "explicit quality/analysis flag")
  # Do not merge inconsistent revisions. With no reliable revision ordering,
  # conflicting duplicate sample groups are all excluded and exported for review.
  duplicate_fields <- c(metadata_fields, "pcr_target_avg_conc", "lod_sewage",
    "hum_frac_mic_conc", "hum_frac_mic_unit", "pcr_target_below_lod", "major_lab_method")
  groups <- split(which(x$disposition == "eligible"), sample_key[x$disposition == "eligible"])
  for (ii in groups[lengths(groups) > 1L]) {
    z <- x[ii, duplicate_fields, drop = FALSE]
    if (nrow(unique(z)) == 1L) x$disposition[ii[-1L]] <- "duplicate identical measurement"
    else x$disposition[ii] <- "conflicting duplicate/revision"
  }
  # Parse only otherwise eligible observations: metadata-free duplicate exports
  # must not acquire units or values by undocumented forward/backward filling.
  x$linear_rna <- x$linear_lod <- x$linear_pmmov <- NA_real_
  eligible <- which(x$disposition == "eligible")
  x$linear_rna[eligible] <- linear_value(num(x$pcr_target_avg_conc[eligible], "RNA"), x$pcr_target_units[eligible])
  x$linear_lod[eligible] <- linear_value(num(x$lod_sewage[eligible], "LOD"), x$pcr_target_units[eligible])
  pmmov <- eligible[lc(x$hum_frac_target_mic[eligible]) == "pepper mild mottle virus" &
    base_units(x$hum_frac_mic_unit[eligible]) == base_units(x$pcr_target_units[eligible])]
  x$linear_pmmov[pmmov] <- linear_value(num(x$hum_frac_mic_conc[pmmov], "PMMoV"), x$hum_frac_mic_unit[pmmov])
  x$below_lod <- yes(x$pcr_target_below_lod)
  x$rna_primary <- x$linear_rna
  x$rna_primary[x$below_lod] <- x$linear_lod[x$below_lod] / sqrt(2)
  x$rna_detected_only <- x$linear_rna
  x$rna_detected_only[x$below_lod] <- NA_real_
  mark(!is.finite(x$rna_primary), "missing/invalid RNA or nondetect without valid LOD")
  x$pmmov_ratio <- x$rna_primary / x$linear_pmmov
  x$pmmov_ratio[!is.finite(x$pmmov_ratio)] <- NA_real_
  method_fields <- c("county_fips", "wwtp_name", "lab_id", "pcr_gene_target",
    "sample_matrix", "pcr_target_units", "sample_type", "concentration_method",
    "extraction_method", "pcr_type", "major_lab_method")
  signature <- do.call(paste, c(lapply(x[method_fields], lc), sep = " | "))
  signatures <- sort(unique(signature[x$disposition == "eligible"]))
  x$series_id <- ifelse(x$disposition == "eligible",
    sprintf("CA%02d", match(signature, signatures)), NA_character_)
  included <- x[x$disposition == "eligible", , drop = FALSE]
  if (!nrow(included)) fail("No eligible wastewater samples remain; inspect the source schema.")
  registry <- do.call(rbind, lapply(split(included, included$series_id), function(z) {
    first <- z[1L, c("series_id", method_fields), drop = FALSE]
    first$county <- z$county[1L]
    first$label <- paste(unique(z$label_name[nzchar(z$label_name)]), collapse = "; ")
    if (!nzchar(first$label)) first$label <- z$wwtp_name[1L]
    first$site_ids <- paste(sort(unique(z$site_id[nzchar(z$site_id)])), collapse = "; ")
    fips <- unique(unlist(regmatches(z$county_names, gregexpr("[0-9]{5}", z$county_names))))
    first$served_county_fips <- paste(sort(fips), collapse = "; ")
    first$geography_note <- if (!length(fips)) "catchment counties unknown" else
      if (any(fips != z$county_fips[1L])) "multi-county catchment; county positivity is a proxy" else
        "source lists one county; sewershed is not the whole county"
    pop <- suppressWarnings(as.numeric(z$population_served))
    first$population_min <- if (any(is.finite(pop))) min(pop, na.rm = TRUE) else NA_real_
    first$population_max <- if (any(is.finite(pop))) max(pop, na.rm = TRUE) else NA_real_
    first$first_sample <- min(z$date); first$last_sample <- max(z$date)
    first$samples <- nrow(z); first$nondetects <- sum(z$below_lod)
    first$missing_quality_flag <- sum(clean(z$quality_flag) == "")
    first$missing_nondetect_flag <- sum(clean(z$pcr_target_below_lod) == "")
    first$major_lab_method[first$major_lab_method == ""] <- "not reported"
    first
  }))
  list(audit = x, samples = included, registry = registry)
}

weekly_samples <- function(s, positivity) {
  # A sample is counted once; same-day samples get equal weight within that day.
  # Weekly means weight sampling days equally, not the number of laboratory rows.
  daily <- do.call(rbind, lapply(split(s, paste(s$series_id, s$date)), function(z) {
    matched <- is.finite(z$pmmov_ratio) & is.finite(z$rna_primary)
    data.frame(series_id = z$series_id[1L], county_fips = z$county_fips[1L], date = z$date[1L],
      rna = mean_ok(z$rna_primary), detected = mean_ok(z$rna_detected_only),
      ratio = mean_ok(z$pmmov_ratio[matched]), matched_rna = mean_ok(z$rna_primary[matched]))
  }))
  daily$week_end <- day_end(daily$date)
  weekly <- do.call(rbind, lapply(split(daily, paste(daily$series_id, daily$week_end)), function(z) {
    data.frame(series_id = z$series_id[1L], county_fips = z$county_fips[1L], week_end = z$week_end[1L],
      sampling_days = sum(is.finite(z$rna)), detected_days = sum(is.finite(z$detected)),
      matched_days = sum(is.finite(z$ratio)), rna_mean = mean_ok(z$rna),
      detected_mean = mean_ok(z$detected), ratio_mean = mean_ok(z$ratio),
      matched_rna_mean = mean_ok(z$matched_rna))
  }))
  weeks <- seq(day_end(ANALYSIS_START), ANALYSIS_END, by = "7 days")
  g <- merge(unique(s[c("series_id", "county_fips")]), data.frame(week_end = weeks), by = NULL)
  w <- merge(g, weekly, by = c("series_id", "county_fips", "week_end"), all.x = TRUE)
  for (k in c("sampling_days", "detected_days", "matched_days")) w[[k]][is.na(w[[k]])] <- 0L
  w <- merge(w, positivity[c("county_fips", "week_end", "positivity_pct")], all.x = TRUE)
  w <- w[order(w$series_id, w$week_end), ]; w$week_start <- w$week_end - 6L
  list(daily = daily, weekly = w)
}
safe_cor <- function(x, y, method = "pearson", minimum = MIN_PAIRED_WEEKS) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < minimum || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}
block_ci <- function(x, y, B, block, method = "pearson") {
  if (B < 1L || sum(is.finite(x) & is.finite(y)) < max(MIN_PAIRED_WEEKS, 3L * block) || length(x) < block)
    return(c(NA_real_, NA_real_))
  n <- length(x)
  # Preserve within-block gaps. Differences have already been computed on the
  # calendar; artificial joins between bootstrap blocks never create differences.
  draws <- replicate(B, {
    starts <- sample.int(n - block + 1L, ceiling(n / block), replace = TRUE)
    ii <- head(unlist(lapply(starts, function(a) seq.int(a, length.out = block))), n)
    safe_cor(x[ii], y[ii], method, minimum = 8L)
  })
  if (sum(is.finite(draws)) < 0.9 * B) return(c(NA_real_, NA_real_))
  as.numeric(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
}
trim_overlap <- function(d, x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(integer(0))
  seq.int(min(which(ok)), max(which(ok)))
}
estimate <- function(d, x, y, metric, period, B = 0L, block = BOOTSTRAP_BLOCK_WEEKS,
                     method = "pearson", variant = "primary") {
  ii <- trim_overlap(d, x, y)
  d <- d[ii, , drop = FALSE]; x <- x[ii]; y <- y[ii]
  ok <- is.finite(x) & is.finite(y)
  ci <- block_ci(x, y, B, block, method)
  data.frame(series_id = if (nrow(d)) d$series_id[1L] else NA_character_, period = period,
    variant = variant, metric = metric, method = method, n_weeks = sum(ok),
    calendar_weeks = length(x), coverage_pct = if (length(x)) 100 * mean(ok) else NA_real_,
    correlation = safe_cor(x, y, method), lower_95 = ci[1L], upper_95 = ci[2L],
    block_weeks = if (B > 0L) block else NA_integer_, bootstrap_replicates = B,
    first_week = if (any(ok)) as.character(min(d$week_end[ok])) else NA_character_,
    last_week = if (any(ok)) as.character(max(d$week_end[ok])) else NA_character_,
    status = if (sum(ok) < MIN_PAIRED_WEEKS) "insufficient paired weeks" else
      if (!is.finite(safe_cor(x, y, method))) "constant series" else "estimated", row.names = NULL)
}
period_data <- function(d, period) {
  if (period == "England dates") return(d[d$week_start >= COMMON_START & d$week_end <= COMMON_END, ])
  if (period == "Full overlap") return(d)
  yr <- as.integer(period)
  d[format(d$week_start, "%Y") == yr & format(d$week_end, "%Y") == yr, ]
}
primary_vectors <- function(d, minimum_days = MIN_SAMPLING_DAYS) {
  w <- d$rna_mean; w[d$sampling_days < minimum_days] <- NA_real_
  p <- d$positivity_pct
  list(raw = w, log = log10(w), p = p,
    dw = weekly_diff(log10(w), d$week_end), dp = weekly_diff(p, d$week_end))
}
lag_associations <- function(d, x, p, period, metric) {
  yy <- vapply(LAG_WEEKS, function(k) p[match(d$week_end + 7L * k, d$week_end)], numeric(nrow(d)))
  common <- is.finite(x) & apply(is.finite(yy), 1L, all)
  data.frame(series_id = rep(d$series_id[1L], length(LAG_WEEKS)), period = period,
    metric = metric, lag_weeks = LAG_WEEKS, n_weeks = sum(common),
    correlation = vapply(seq_along(LAG_WEEKS), function(j)
      safe_cor(x[common], yy[common, j]), numeric(1)),
    interpretation = "exploratory; common dates across lags; no forecast claim")
}
analyse <- function(weekly) {
  rows <- sensitivity <- lags <- list()
  add <- function(row, id) { row$series_id <- id; row }
  for (id in unique(weekly$series_id)) {
    message("Evaluating ", id, "...")
    all <- weekly[weekly$series_id == id, ]
    for (period in c("England dates", "Full overlap", "2020", "2021", "2022", "2023")) {
      d <- period_data(all, period); v <- primary_vectors(d)
      primary_period <- period %in% c("England dates", "Full overlap")
      B <- if (primary_period) BOOTSTRAP_REPLICATES else 0L
      specifications <- list(
        list(v$raw, v$p, "RNA level", "pearson"),
        list(v$log, v$p, "log10 RNA level", "pearson"),
        list(v$raw, v$p, "RNA rank", "spearman"),
        list(v$dw, v$dp, "weekly changes", "pearson"))
      for (a in specifications) rows[[length(rows) + 1L]] <- add(
        estimate(d, a[[1]], a[[2]], a[[3]], period, B = B, method = a[[4]]), id)
      if (!primary_period) next
      # Normalization comparison uses identical samples/days and identical weeks.
      keep <- d$matched_days >= MIN_SAMPLING_DAYS
      raw_matched <- ratio <- rep(NA_real_, nrow(d))
      raw_matched[keep] <- log10(d$matched_rna_mean[keep])
      ratio[keep] <- log10(d$ratio_mean[keep])
      detected <- log10(d$detected_mean); detected[d$detected_days < MIN_SAMPLING_DAYS] <- NA_real_
      for (a in list(list(raw_matched, "matched raw RNA"), list(ratio, "PMMoV ratio"),
                     list(detected, "exclude nondetects"))) {
        sensitivity[[length(sensitivity) + 1L]] <- add(estimate(d, a[[1]], v$p,
          "log10 RNA level", period, variant = a[[2]]), id)
      }
      for (threshold in c(1L, 3L)) {
        vv <- primary_vectors(d, threshold)
        sensitivity[[length(sensitivity) + 1L]] <- add(estimate(d, vv$log, vv$p,
          "log10 RNA level", period, variant = paste0("minimum ", threshold, " sampling days")), id)
      }
      for (block in c(2L, 8L)) sensitivity[[length(sensitivity) + 1L]] <- add(
        estimate(d, v$log, v$p, "log10 RNA level", period, B = B, block = block,
          variant = paste0(block, "-week bootstrap blocks")), id)
      # Each lag is tested on the same wastewater dates within this period.
      # Shift by date, never by row number after dropping missing weeks.
      for (metric in c("log10 RNA level", "weekly changes")) {
        xx <- if (metric == "log10 RNA level") v$log else v$dw
        pp <- if (metric == "log10 RNA level") v$p else v$dp
        lags[[length(lags) + 1L]] <- lag_associations(d, xx, pp, period, metric)
      }
    }
  }
  list(correlations = do.call(rbind, rows), sensitivities = do.call(rbind, sensitivity),
    lags = do.call(rbind, lags))
}

esc <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE)
}
html_table <- function(x, digits = 3L) {
  for (j in seq_along(x)) if (is.numeric(x[[j]])) x[[j]] <- ifelse(is.na(x[[j]]), "", formatC(x[[j]], digits = digits, format = "f"))
  paste0("<div class='scroll'><table><thead><tr>", paste0("<th>", esc(names(x)), "</th>", collapse = ""),
    "</tr></thead><tbody>", paste(apply(x, 1L, function(z) paste0("<tr>",
      paste0("<td>", esc(z), "</td>", collapse = ""), "</tr>")), collapse = ""), "</tbody></table></div>")
}
png_embed <- function(p, alt) {
  raw <- readBin(p, "raw", n = file.info(p)$size)
  paste0("<img alt='", esc(alt), "' src='data:image/png;base64,", jsonlite::base64_enc(raw), "'>")
}
plots <- function(weekly, registry, results, folder) {
  images <- character()
  palette <- c("San Francisco" = "#2563a6", "Santa Clara" = "#16816a", "Sacramento" = "#b76227")
  z <- merge(subset(results$correlations, period == "England dates" & metric == "log10 RNA level" &
    status == "estimated"), registry[c("series_id", "county", "label", "lab_id")], by = "series_id")
  if (nrow(z)) {
    p <- file.path(folder, "common_period_correlations.png")
    png(p, width = 1400, height = max(650, 95 * nrow(z)), res = 140)
    par(mar = c(5, 15, 4, 2)); y <- rev(seq_len(nrow(z)))
    plot(z$correlation, y, xlim = c(-1, 1), ylim = c(0.5, nrow(z) + 0.5), yaxt = "n",
      ylab = "", xlab = "Pearson r: log10 wastewater RNA versus county positivity",
      main = "England date window: 9 June 2021 to 30 March 2022", pch = 19,
      col = palette[z$county], cex = 1.2)
    axis(2, at = y, labels = paste(z$series_id, z$label, z$lab_id, paste0("n=", z$n_weeks)), las = 1, cex.axis = 0.8)
    segments(z$lower_95, y, z$upper_95, y, col = palette[z$county], lwd = 2)
    abline(v = 0, lty = 2, col = "grey55"); dev.off()
    images["common"] <- p
  }
  for (county in unname(COUNTIES)) {
    ids <- registry$series_id[registry$county == county]
    # Plot every eligible series, including segments too short for estimation.
    p <- file.path(folder, paste0(gsub(" ", "_", tolower(county)), "_series.png"))
    png(p, width = 1500, height = 420 * (length(ids) + 1L), res = 130)
    par(mfrow = c(length(ids) + 1L, 1L), mar = c(3, 5, 2.2, 1), oma = c(1, 0, 2, 0))
    d <- weekly[weekly$series_id == ids[1L], ]
    plot(d$week_end, d$positivity_pct, type = "l", lwd = 2, col = "#bc4143", xlab = "",
      ylab = "NAAT positive (%)", main = paste(county, "county clinical testing"), ylim = c(0, max(d$positivity_pct, na.rm = TRUE) * 1.05))
    for (id in ids) {
      d <- weekly[weekly$series_id == id, ]; reg <- registry[registry$series_id == id, ]
      w <- primary_vectors(d)$log
      plot(d$week_end, w, type = "o", pch = 16, cex = 0.3, col = palette[county], xlab = "",
        ylim = if (any(is.finite(w))) range(w, na.rm = TRUE) else c(0, 1),
        ylab = "log10 RNA", main = paste(id, reg$label, reg$lab_id,
          paste0("method: ", reg$major_lab_method), "|", reg$pcr_target_units))
      abline(v = as.numeric(c(COMMON_START, COMMON_END)), lty = 3, col = "grey55")
    }
    mtext("Separate scales and methods; gaps retained; Wednesday week ends", outer = TRUE, cex = 1)
    dev.off(); images[county] <- p
  }
  images
}

report <- function(prepared, weekly, positivity, results, images, sources, run) {
  reg <- prepared$registry
  headline <- merge(subset(results$correlations, period == "England dates" &
    metric == "log10 RNA level"), reg[c("series_id", "county", "label", "lab_id")], by = "series_id")
  conclusion <- vapply(unname(COUNTIES), function(county) {
    z <- headline[headline$county == county & headline$status == "estimated", ]
    if (!nrow(z)) return(paste0(county, ": no method segment meets the paired-week threshold in the England date window."))
    paste0(county, ": ", nrow(z), " eligible plant/method series; log-RNA Pearson correlations range from ",
      sprintf("%.2f", min(z$correlation)), " to ", sprintf("%.2f", max(z$correlation)),
      ". ", sum(z$lower_95 > 0, na.rm = TRUE), " of ", nrow(z),
      " have a pointwise 95% block-bootstrap interval wholly above zero. These are not independent county estimates.")
  }, character(1))
  dispositions <- as.data.frame(table(prepared$audit$county, prepared$audit$disposition), stringsAsFactors = FALSE)
  names(dispositions) <- c("county", "disposition", "rows"); dispositions <- dispositions[dispositions$rows > 0, ]
  pos_summary <- do.call(rbind, lapply(split(positivity, positivity$county), function(z) {
    ok <- is.finite(z$positivity_pct)
    data.frame(county = z$county[1L], observed_weeks = sum(ok), missing_weeks = sum(!ok),
      first_observed = if (any(ok)) as.character(min(z$week_end[ok])) else "none",
      last_observed = if (any(ok)) as.character(max(z$week_end[ok])) else "none")
  }))
  main_table <- merge(subset(results$correlations, period %in% c("England dates", "Full overlap") &
    metric %in% c("log10 RNA level", "weekly changes")), reg[c("series_id", "county", "label", "lab_id")])
  main_table <- main_table[order(main_table$county, main_table$series_id, main_table$period, main_table$metric),
    c("county", "series_id", "label", "lab_id", "period", "metric", "n_weeks", "coverage_pct",
      "first_week", "last_week", "correlation", "lower_95", "upper_95", "status")]
  normalized <- subset(results$sensitivities, variant %in% c("matched raw RNA", "PMMoV ratio") & period == "England dates")
  h <- c("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>California wastewater and COVID-19 positivity</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#172d3c;max-width:1200px;margin:auto;padding:32px}h1,h2{line-height:1.2}h2{margin-top:2em}p{max-width:1000px}.note{padding:18px;background:#eef4f8;border-left:5px solid #2563a6}.scroll{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:8px;border-bottom:1px solid #d9e2e8;text-align:left;vertical-align:top}th{background:#eef4f8;position:sticky;top:0}img{max-width:100%;height:auto}code{background:#eef4f8;padding:2px 4px}a{color:#175ba0}summary{cursor:pointer;font-weight:600}</style>",
    "<body><h1>California wastewater RNA versus COVID-19 test positivity</h1>",
    paste0("<p>San Francisco, Santa Clara and Sacramento. Generated ", esc(utc_now()), ".</p>"),
    "<div class='note'><strong>Question:</strong> Do measured wastewater concentrations track county NAAT positivity, and is the association stable across time, methods and normalization? This retrospective analysis does not establish a conversion from RNA to infections or validate an early-warning forecast.</div>",
    "<h2>Findings in the England date window</h2>",
    paste0("<p>", esc(conclusion), "</p>"),
    "<p>The range summarizes all estimable plant/method series, not a pooled county effect or a selection of the best result. Shared county positivity makes estimates dependent. Intervals are pointwise and approximate, without multiplicity correction. A positive level correlation can reflect shared epidemic waves; the weekly-change results below test a different, shorter-term association.</p>",
    if ("common" %in% names(images)) png_embed(images[["common"]], "Correlations and block-bootstrap intervals"),
    "<h2>Data checks before interpretation</h2>",
    "<p>Wastewater: California's historical Cal-SuWers sample archive. Clinical data: CDC revised weekly county archive jgk8-6dpn, field <code>percent_test_results_reported</code>, diagnostic and screening NAAT percentage positive. It is not a case-rate measure. This positivity archive ends on 10 May 2023; later wastewater is excluded.</p>",
    html_table(pos_summary, 0),
    "<p>There are no test-count denominators in this extraction: county percentages are never pooled and precision cannot be weighted by test volume. Testing eligibility, home-test uptake and reporting changed over the pandemic; positivity is not population prevalence. The revised archive is not a sequence of information available in real time.</p>",
    "<p><strong>Geography:</strong> plants are selected by the county containing the plant, then catchment counties are audited separately. San Francisco records can include San Mateo County in the catchment. Even a one-county catchment covers only part of that county. No county-wide wastewater composite is constructed; population metadata are shown rather than used to claim complete county coverage.</p>",
    "<p><strong>Method and revision checks:</strong> only broad SARS-CoV-2 N/N1/N2 targets at confirmed treatment plants are included. Upstream sites, incomplete method metadata, inconsistent matrix/units, explicit quality or analysis exclusions and conflicting duplicate samples are excluded. Identical measurement duplicates are counted once. Blank QC flags mean unknown, not certified quality. The <code>qc_ignore</code> field controls the QC report display and is not treated as a failure flag.</p>",
    "<p>The downloaded archive contains metadata-free duplicate exports, a dry-sludge-unit/raw-wastewater-matrix conflict for some Gilroy records, and a transition from an unreported major method ID to ID 1 for some VLT records in February 2023. The script does not silently repair these. Method-ID segments remain separate; a metadata transition need not mean the physical assay changed. No unique latest revision can be inferred from conflicting rows. Refreshing a source records a new hash; a single snapshot cannot reconstruct the full revision history.</p>",
    "<p>Many CAL4 sample identifiers have multiple differing concentrations and laboratory measurements. These could represent repeat analyses or revisions; the available identifiers do not resolve that distinction. The conservative duplicate rule removes such groups, substantially reducing CAL4 coverage. Those methods therefore may not support a common-window estimate. This exclusion is a data-identification limitation, not evidence that their correlation is low.</p>",
    html_table(dispositions, 0),
    "<details><summary>Included plant/method segments and catchment audit</summary>",
    html_table(reg[c("series_id", "county", "label", "lab_id", "pcr_gene_target", "sample_matrix", "pcr_target_units",
      "sample_type", "major_lab_method", "first_sample", "last_sample", "samples", "nondetects",
      "missing_quality_flag", "missing_nondetect_flag", "served_county_fips", "geography_note")], 0), "</details>",
    "<h2>Time alignment and statistical method</h2>",
    paste0("<p>Wastewater collection dates are assigned to Thursday-Wednesday weeks matching CDC's Wednesday endpoints. Replicate samples are averaged within day, then sampled days are averaged within week. At least ", MIN_SAMPLING_DAYS,
      " distinct days are required. Missing weeks remain missing. The primary common window contains complete weeks within 1 June 2021-30 March 2022 (endpoints 9 June 2021-30 March 2022). Each series also uses its full overlap with positivity. At least ", MIN_PAIRED_WEEKS, " paired weeks are required for estimation.</p>"),
    "<p>Log10 copies/L source values, including LODs and matching PMMoV values, are converted to the linear scale before averaging or dividing. Primary nondetect substitution is LOD/sqrt(2) in linear units, only for an explicit below-LOD flag and a valid LOD. The exclusion sensitivity removes those observations instead. Unflagged missing or nonpositive measurements are not invented.</p>",
    paste0("<p>Outputs include Pearson correlation on raw and log10 RNA, Spearman correlation, and Pearson correlation of consecutive-week log10-RNA changes with positivity changes in percentage points. Differences never bridge missing weeks. The main periods use paired moving-block bootstrap intervals (", BOOTSTRAP_REPLICATES,
      " draws, ", BOOTSTRAP_BLOCK_WEEKS, "-week blocks); 2- and 8-week blocks are sensitivity checks. Blocks preserve within-block missingness. These approximate intervals can remain unreliable with strong nonstationarity or long gaps; there are no IID p-values or causal claims.</p>"),
    "<h2>Correlations and weekly changes</h2>", html_table(main_table),
    "<h2>PMMoV normalization</h2><p>Both numerator and fecal marker must have compatible units in the same sample. Daily sample ratios are averaged into weekly ratios; the comparison RNA uses exactly those samples/days and weeks. This comparison is descriptive, not a search for whichever metric produces the largest correlation. The complete sensitivity CSV also includes sampling thresholds, nondetect exclusions and block lengths.</p>",
    html_table(normalized[c("series_id", "variant", "n_weeks", "correlation", "status")]),
    "<h2>Weekly trajectories</h2><p>Every included method segment is shown separately. Panels use their own concentration scales; no solid/liquid concentrations are spliced. Dotted vertical lines mark the England date window. The HTML embeds the images and can be viewed offline.</p>",
    unlist(lapply(unname(COUNTIES), function(county) c(paste0("<h3>", county, "</h3>"), png_embed(images[[county]], paste(county, "weekly trajectories"))))),
    "<h2>Exploratory lag analysis</h2><p>Positive lag k compares wastewater in week t with positivity in week t+k; negative lag means positivity precedes wastewater. All lags within a series/period use identical wastewater dates with complete clinical observations across the whole lag range. Levels and changes are exported separately. Do not describe the highest correlation as a proven lead time: autocorrelation, shared waves, testing delays and trying multiple lags can produce apparent leads. Sampling dates and retrospective positivity do not measure data-release delays.</p>",
    "<details><summary>Lag results</summary>", html_table(results$lags), "</details>",
    "<h2>Conclusion and limits of the England comparison</h2>",
    "<p>Use the common-window log-RNA estimates together with changes, calendar-year stability, missingness and PMMoV sensitivities. Stronger correlations in a California plant series would support a local surveillance association, not prove that US wastewater monitoring is intrinsically better than England's national aggregate. The England analysis uses a different geography and processing pipeline. A forecast claim would require a prespecified model, data available at each forecast origin and a held-out period.</p>",
    "<p>When the association is weak, the result concerns these measurements, populations and dates. It does not establish that wastewater is unrelated to infection burden or lacks surveillance value.</p>",
    "<h2>Exports and reproducibility</h2><p>CSV exports: sample audit with original row numbers and dispositions; eligible daily and complete weekly series; method/catchment registry; all primary and calendar-year correlations; sensitivities; exploratory lags; positivity; source manifest. The original statewide CSV stays in the reusable cache (or the supplied local path); the selected three-county snapshot is copied into each run as the sample audit. Keep the cache to preserve the complete source. No previous successful run is overwritten.</p>",
    "<p><a href='data/correlations.csv'>All correlations</a> · <a href='data/sensitivities.csv'>Sensitivities</a> · <a href='data/lag_scan.csv'>Lags</a> · <a href='data/series_registry.csv'>Method registry</a> · <a href='sources.csv'>Source manifest</a>. CSV links require the report folder; the HTML itself is self-contained.</p>",
    "<h2>Official sources</h2><ul>",
    paste0("<li><a href='", CA_LANDING, "'>California historical wastewater dataset</a>; <a href='", WW_URL, "'>sample CSV</a>; <a href='", DICTIONARY_URL, "'>NWSS metadata dictionary</a>.</li>"),
    paste0("<li><a href='", CDC_LANDING, "'>CDC revised weekly county archive</a>; <a href='", CDC_METADATA, "'>metric metadata and revision notes</a>.</li></ul>"),
    "<p>Source snapshots, URLs, acquisition timestamps and MD5 checksums are recorded below. MD5 is an integrity fingerprint, not a security guarantee.</p>",
    html_table(sources), "</body></html>")
  writeLines(h, file.path(run, "california_wastewater_vs_positivity.html"), useBytes = TRUE)
  writeLines(conclusion, file.path(run, "conclusions.txt"))
}

run_california_wastewater_vs_positivity <- function() {
  need_packages()
  if (BOOTSTRAP_REPLICATES < 100L || MIN_PAIRED_WEEKS < 12L || MIN_SAMPLING_DAYS < 1L ||
      MIN_SAMPLING_DAYS > 7L) fail("Invalid analysis settings.")
  if (nzchar(LOCAL_POSITIVITY_CSV) && !nzchar(LOCAL_POSITIVITY_SOURCE))
    fail("Set LOCAL_POSITIVITY_SOURCE to document the observed county NAAT metric and weekly window.")
  dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(OUTPUT_ROOT, "cache"); dir.create(cache, showWarnings = FALSE)
  run <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC"), "_"), tmpdir = OUTPUT_ROOT)
  dir.create(run); dir.create(file.path(run, "data")); dir.create(file.path(run, "plots"))
  dir.create(file.path(run, "source_metadata"))
  complete <- FALSE
  on.exit({
    capture.output(sessionInfo(), file = file.path(run, "session_info.txt"))
    if (!complete) writeLines(paste("Incomplete run", utc_now()), file.path(run, "RUN_FAILED.txt"))
  }, add = TRUE)
  set.seed(RANDOM_SEED)
  writeLines(capture.output(dput(list(counties = COUNTIES, analysis_start = ANALYSIS_START,
    analysis_end = ANALYSIS_END, common_start = COMMON_START, common_end = COMMON_END,
    minimum_days = MIN_SAMPLING_DAYS, minimum_paired_weeks = MIN_PAIRED_WEEKS,
    bootstrap_replicates = BOOTSTRAP_REPLICATES, bootstrap_block_weeks = BOOTSTRAP_BLOCK_WEEKS,
    seed = RANDOM_SEED, lag_weeks = LAG_WEEKS, refresh = REFRESH_DOWNLOADS,
    local_positivity_source = LOCAL_POSITIVITY_SOURCE))), file.path(run, "settings.R"))
  sources <- list()
  acquire <- function(label, url, name, local = "", validate = function(p) TRUE) {
    s <- cached_source(url, file.path(cache, name), REFRESH_DOWNLOADS, local, validate)
    sources[[label]] <<- data.frame(source = label, url = if (nzchar(local)) local else url,
      path = s$path, acquired_utc = if (is.null(s$acquired_utc)) NA_character_ else s$acquired_utc,
      inspected_utc = utc_now(), md5 = s$md5, bytes = s$bytes, mode = s$mode)
    s$path
  }
  ww_file <- acquire("California wastewater", WW_URL, "california_wastewater.csv",
    LOCAL_WASTEWATER_CSV, csv_validator(WW_COLUMNS))
  where <- paste0("fips_code in ('", paste(names(COUNTIES), collapse = "','"), "')")
  pos_url <- query_url(c(`$select` = "state_name,county_name,fips_code,date,percent_test_results_reported",
    `$where` = where, `$order` = "date,fips_code", `$limit` = "5000"))
  pos_file <- acquire("CDC county positivity", pos_url, "county_positivity.csv", LOCAL_POSITIVITY_CSV)
  raw_pos <- read_text_csv(pos_file)
  if (!nzchar(LOCAL_POSITIVITY_CSV)) {
    check_cols(raw_pos, c("fips_code", "date", "percent_test_results_reported"))
    count_url <- query_url(c(`$select` = "count(*) as n", `$where` = where))
    count_file <- acquire("CDC row count", count_url, "positivity_count.csv", validate = csv_validator("n"))
    expected <- num(read_text_csv(count_file)$n, "API count")
    if (length(expected) != 1L || is.na(expected) || expected >= 5000L || nrow(raw_pos) != expected)
      fail("CDC extraction is incomplete or exceeds the safety limit; review/paginate the API.")
  }
  for (a in list(c("California metadata", CA_API, "california_metadata.json"),
                 c("CDC metadata", CDC_METADATA, "cdc_metadata.json"),
                 c("NWSS dictionary", DICTIONARY_URL, "nwss_metadata.csv"))) {
    p <- acquire(a[1], a[2], a[3], validate = if (grepl("json$", a[3])) json_validator else function(p) TRUE)
    if (!file.copy(p, file.path(run, "source_metadata", a[3]))) fail("Cannot copy source metadata.")
  }
  positivity <- validate_positivity(raw_pos)
  if (any(vapply(split(positivity$positivity_pct, positivity$county_fips), function(x) sum(is.finite(x)), integer(1)) < MIN_PAIRED_WEEKS))
    fail("A county has too few observed clinical weeks.")
  message("Reading selected wastewater columns and auditing samples...")
  raw <- read_wastewater(ww_file)
  prepared <- prepare_samples(raw)
  weekly <- weekly_samples(prepared$samples, positivity)
  write_csv(prepared$audit, file.path(run, "data", "sample_audit.csv"))
  write_csv(prepared$registry, file.path(run, "data", "series_registry.csv"))
  write_csv(weekly$daily, file.path(run, "data", "daily_samples.csv"))
  write_csv(weekly$weekly, file.path(run, "data", "weekly_pairs.csv"))
  write_csv(positivity, file.path(run, "data", "county_positivity.csv"))
  results <- analyse(weekly$weekly)
  write_csv(results$correlations, file.path(run, "data", "correlations.csv"))
  write_csv(results$sensitivities, file.path(run, "data", "sensitivities.csv"))
  write_csv(results$lags, file.path(run, "data", "lag_scan.csv"))
  sources <- do.call(rbind, sources); rownames(sources) <- NULL
  write_csv(sources, file.path(run, "sources.csv"))
  images <- plots(weekly$weekly, prepared$registry, results, file.path(run, "plots"))
  report(prepared, weekly$weekly, positivity, results, images, sources, run)
  writeLines(paste("Complete", utc_now()), file.path(run, "RUN_COMPLETE.txt"))
  complete <- TRUE
  result <- normalizePath(file.path(run, "california_wastewater_vs_positivity.html"), winslash = "/")
  message("Completed report: ", result)
  invisible(list(report = result, directory = normalizePath(run, winslash = "/"), results = results,
    registry = prepared$registry))
}

if (isTRUE(RUN_ANALYSIS)) run_california_wastewater_vs_positivity()
