# Netherlands: quarterly sickness absence, 2010 Q1--2025 Q4
# Run with source("netherlands_sickness_rates.R", encoding = "UTF-8").
# First install: install.packages(c("curl", "jsonlite", "strucchange", "sandwich"))
#
# SOURCE REVIEW (22 October 2026; read before interpreting the results)
# * 80072ned's official CSV URL delivers a ZIP of semicolon-separated CSVs;
#   observations use decimal commas. All 36 industry/aggregate series have
#   64 numeric quarters in 2010--2025, no duplicate keys or missing values.
#   Three additional company-size series are extracted but excluded from models.
#   The 18 August 2026 release marks 2010--2025 final and 2026 provisional.
# * The supplied 24015NED query selects YEARS (JJ). This script selects KW
#   quarters and all AZW branches, retaining the requested B000125 when present.
#   The AZW API returned HTTP 406 during preparation: its current values,
#   completeness, release status and definitions could NOT be verified here.
#   They are downloaded and audited at runtime; unavailable sources are flagged.
# * CBS documents enterprise-group/register changes from 2010 that can affect
#   SBI/size classification. Treat 2010 as a boundary; do not splice pre-2010 data.
#   Rates measure days lost / available employee days, NOT people absent.
#   Long-term absence is included; ordinary maternity/childbirth leave excluded.
#   Employers may calculate using workdays or calendar days. Data are rounded
#   to 0.1 percentage point and are survey estimates (sector precision varies).
# * AZW broad health/welfare and CBS SBI Q cover the same aggregate. Sectors
#   and their aggregates overlap: do not sum rates or regard them as independent
#   samples. Definitions of narrower AZW and SBI groups need not coincide.
# * Final does not mean never revised. A first download cannot establish the
#   history of revisions. Subsequent runs compare values/status/definitions
#   against the previous successful snapshot; original source files are kept.
#
# METHODS
# Per sector: rate_pct ~ quarter_index + factor(quarter), fitted ONLY to the
# 36 quarters in 2010--2018; extrapolated unchanged through 2025.
# Bai--Perron level-shift model: deviation_pp ~ 1, on UNSMOOTHED observations
# over the full 2010--2025 window. BIC selects 0..6 breaks independently for
# each sector, with >=8 quarters per segment. No English or pandemic dates.
# These are shifts in the deviation level, not a model of changing slopes.
# HAC baseline mean intervals use Newey--West lag 4; they are NOT prediction
# intervals. Break-date intervals are conditional on the estimated baseline
# and selected break count. BIC is a Gaussian least-squares criterion; serial
# dependence, baseline misspecification and changing seasonal patterns can
# affect the selected breaks. No causal attribution or multiplicity claim.
#
# REFERENCES (official sources and package documentation)
# https://datasets.cbs.nl/CSV/CBS/nl/80072ned
# https://dataderden.cbs.nl/ODataApi/OData/24015NED
# https://www.cbs.nl/nl-nl/cijfers/detail/80072ned
# https://www.cbs.nl/nl-nl/onze-diensten/methoden/onderzoeksomschrijvingen/korte-onderzoeksbeschrijvingen/ziekteverzuim-kwartaalenquete--vanaf-1996--
# https://www.cbs.nl/nl-nl/nieuws/2022/10/ziekteverzuim-in-vierde-kwartaal-2021-verder-toegenomen/bedrijfstakken-en-branches
# https://cran.r-project.org/web/packages/strucchange/strucchange.pdf

# ---- Configuration ----------------------------------------------------------
OUT_DIR <- "netherlands_sickness_results"
SOURCE_MODE <- "refresh"          # "refresh" or explicit offline "cache"
SOURCES <- c("80072ned", "24015NED")
STRICT_SOURCES <- FALSE           # TRUE: abort modelling if any source failed
INCLUDE_COMPANY_SIZE <- FALSE
SECTOR_IDS <- NULL                # NULL = all; e.g. c("422400", "B000125")
START_YEAR <- 2010L
END_YEAR <- 2025L
BASELINE_END <- 2018L
MIN_SEGMENT_QUARTERS <- 8L
MAX_BREAKS <- 6L
HAC_LAG <- 4L
MAKE_PLOTS <- TRUE
RUN_ANALYSIS <- TRUE              # FALSE: load functions only

need_packages <- function() {
  pkgs <- c("curl", "jsonlite", "strucchange", "sandwich")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages first: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}
write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
require_cols <- function(x, cols, where) {
  if (!is.data.frame(x) || !all(cols %in% names(x)))
    stop("Unexpected schema in ", where, "; required: ", paste(cols, collapse = ", "))
}
get_col <- function(x, name, default = "") {
  if (name %in% names(x)) as.character(x[[name]]) else rep(default, nrow(x))
}
download_file <- function(url, dest) {
  h <- curl::new_handle(timeout = 60, connecttimeout = 20, followlocation = TRUE,
    useragent = "R CBS quarterly sickness analysis",
    httpheader = c(Accept = "application/json, application/zip, */*"))
  r <- curl::curl_fetch_disk(url, dest, handle = h)
  if (r$status_code != 200L || !file.exists(dest) || file.info(dest)$size == 0)
    stop("Download failed (HTTP ", r$status_code, "): ", url,
         ". No stale cache was substituted.", call. = FALSE)
  invisible(dest)
}
odata_get <- function(resource, raw_dir, query = "") {
  root <- "https://dataderden.cbs.nl/ODataApi/OData/24015NED/"
  url <- paste0(root, resource, query)
  pages <- list(); seen <- character()
  repeat {
    if (url %in% seen) stop("Repeated OData pagination URL: ", resource)
    seen <- c(seen, url)
    p <- file.path(raw_dir, sprintf("%s_%03d.json", resource, length(seen)))
    download_file(url, p)
    obj <- jsonlite::fromJSON(p, simplifyDataFrame = TRUE)
    value <- obj$value
    next_url <- obj[["odata.nextLink"]]
    if (is.null(next_url)) next_url <- obj[["@odata.nextLink"]]
    if (is.null(value) && !is.null(obj$d)) {
      value <- obj$d$results; next_url <- obj$d$`__next`
    }
    if (is.null(value)) stop("OData response has no data array: ", resource)
    if (length(value)) pages[[length(pages) + 1L]] <- as.data.frame(value)
    if (is.null(next_url) || !length(next_url) || !nzchar(next_url)) break
    if (startsWith(next_url, "/")) next_url <- paste0("https://dataderden.cbs.nl", next_url)
    if (!grepl("^https?://", next_url)) next_url <- paste0(root, next_url)
    url <- next_url
  }
  if (!length(pages)) return(data.frame())
  do.call(rbind, pages)
}
parse_rate <- function(x) {
  # Blank, dots, suppression symbols etc. remain NA, never zero.
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE)))
}
read_cbs_zip <- function(path) {
  z <- utils::unzip(path, list = TRUE)
  scratch <- tempfile("cbs_zip_")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  read_part <- function(name) {
    if (!name %in% z$Name) stop("Missing ZIP member: ", name)
    utils::unzip(path, files = name, exdir = scratch)
    txt <- readLines(file.path(scratch, name), encoding = "UTF-8", warn = FALSE)
    txt[1] <- sub("^\ufeff", "", txt[1])
    x <- utils::read.csv2(text = paste(c(txt, ""), collapse = "\n"),
                        colClasses = "character", na.strings = character(),
                        check.names = FALSE, stringsAsFactors = FALSE)
    names(x) <- sub("^\ufeff", "", names(x)); x
  }
  obs <- read_part("Observations.csv")
  codes <- read_part("BedrijfskenmerkenSBI2008Codes.csv")
  groups <- read_part("BedrijfskenmerkenSBI2008Groups.csv")
  periods <- read_part("PeriodenCodes.csv")
  measures <- read_part("MeasureCodes.csv")
  props <- read_part("Properties.csv")
  require_cols(obs, c("Measure", "BedrijfskenmerkenSBI2008", "Perioden", "Value"), "CBS observations")
  require_cols(codes, c("Identifier", "Title", "DimensionGroupId"), "CBS sector codes")
  require_cols(measures, c("Identifier", "Title", "Unit"), "CBS measures")
  m <- which(trimws(measures$Title) == "Ziekteverzuimpercentage")
  if (length(m) != 1L || measures$Unit[m] != "%") stop("Unexpected CBS measure/unit")
  obs <- obs[obs$Measure == measures$Identifier[m], , drop = FALSE]
  g <- trimws(groups$Title[match(codes$DimensionGroupId, groups$Id)])
  sectors <- data.frame(sector_id = trimws(codes$Identifier), sector = trimws(codes$Title),
    sector_group = g, company_size = g == "Bedrijfsgrootte", stringsAsFactors = FALSE)
  if (anyNA(sectors$company_size)) stop("Unmapped CBS sector group")
  p <- match(trimws(obs$Perioden), trimws(periods$Identifier))
  d <- data.frame(sector_id = trimws(obs$BedrijfskenmerkenSBI2008),
    period_code = trimws(obs$Perioden), rate_raw = obs$Value,
    value_attribute = get_col(obs, "ValueAttribute"),
    period_status = get_col(periods, "Status", "unknown")[p], stringsAsFactors = FALSE)
  list(data = d, sectors = sectors,
       metadata = list(properties = props, measures = measures, sectors = codes,
                       groups = groups, periods = periods),
       definitions = list(measures = measures, sectors = codes, groups = groups),
       source_url = "https://datasets.cbs.nl/CSV/CBS/nl/80072ned")
}
read_azw <- function(raw_dir) {
  info <- odata_get("TableInfos", raw_dir)
  codes <- odata_get("AZWBranches", raw_dir)
  periods <- odata_get("Perioden", raw_dir)
  properties <- odata_get("DataProperties", raw_dir)
  # No backslashes before $, and KW instead of the annual JJ filter.
  obs <- odata_get("TypedDataSet", raw_dir,
                  "?$filter=substringof('KW',Perioden)")
  require_cols(obs, c("AZWBranches", "Perioden", "Ziekteverzuimpercentage_1"), "AZW observations")
  require_cols(codes, c("Key", "Title"), "AZW sector codes")
  require_cols(periods, c("Key", "Title"), "AZW periods")
  require_cols(properties, c("Key", "Unit"), "AZW measures")
  measure <- properties[properties$Key == "Ziekteverzuimpercentage_1", , drop = FALSE]
  if (nrow(measure) != 1L || trimws(measure$Unit) != "%") stop("Unexpected AZW measure/unit")
  if (!"B000125" %in% trimws(codes$Key)) stop("Requested AZW branch B000125 is absent")
  sectors <- data.frame(sector_id = trimws(codes$Key), sector = trimws(codes$Title),
    sector_group = "AZW (see source definitions for aggregates)", company_size = FALSE)
  p <- match(trimws(obs$Perioden), trimws(periods$Key))
  status <- get_col(periods, "Status", "unknown")[p]
  status[is.na(status) | status == ""] <- "unknown"
  title <- periods$Title[p]
  status[!is.na(title) & grepl("\\*", title)] <- "Voorlopig (period title)"
  # Old OData tables may expose provisional status only in their description.
  desc <- paste(unlist(info), collapse = "\n")
  match_final <- regmatches(desc, regexec("tot en met ([0-9]{4}) zijn definitief", desc))[[1]]
  if (length(match_final) == 2L) {
    final_year <- as.integer(match_final[2])
    is_final <- as.integer(substr(obs$Perioden, 1, 4)) <= final_year & status == "unknown"
    status[is_final] <- "Definitief (table description)"
  }
  d <- data.frame(sector_id = trimws(obs$AZWBranches), period_code = trimws(obs$Perioden),
    rate_raw = as.character(obs$Ziekteverzuimpercentage_1), value_attribute = "not supplied by TypedDataSet",
    period_status = status, stringsAsFactors = FALSE)
  list(data = d, sectors = sectors,
       metadata = list(properties = info, measures = properties, sectors = codes, periods = periods),
       definitions = list(measures = properties, sectors = codes),
       source_url = "https://dataderden.cbs.nl/ODataApi/OData/24015NED/")
}
load_source <- function(id, raw_dir, cache_dir) {
  cache <- file.path(cache_dir, paste0(id, "_raw.rds"))
  if (SOURCE_MODE == "cache") {
    if (!file.exists(cache)) stop("No cached source: ", id, "; use refresh first")
    return(readRDS(cache))
  }
  if (id == "80072ned") {
    p <- file.path(raw_dir, "80072ned.zip")
    download_file("https://datasets.cbs.nl/CSV/CBS/nl/80072ned", p)
    obj <- read_cbs_zip(p)
  } else if (id == "24015NED") obj <- read_azw(raw_dir)
  else stop("Unsupported source: ", id)
  obj$retrieved_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  obj$source_files_md5 <- tools::md5sum(list.files(raw_dir, full.names = TRUE))
  # Cache only a successfully parsed complete source, preserving previous cache on failure.
  saveRDS(obj, cache)
  obj
}
quarter_grid <- function() {
  year <- rep(seq.int(START_YEAR, END_YEAR), each = 4L)
  q <- rep(1:4, END_YEAR - START_YEAR + 1L)
  data.frame(period_code = sprintf("%dKW%02d", year, q), year = year,
    quarter = q, quarter_label = paste0(year, "Q", q),
    date = as.Date(sprintf("%d-%02d-01", year, 3L * q - 2L)),
    quarter_index = seq_along(year) - 1L)
}
audit_source <- function(obj, id) {
  d <- obj$data
  grid <- quarter_grid()
  d <- d[grepl("^[0-9]{4}KW0[1-4]$", d$period_code) &
           d$period_code %in% grid$period_code, , drop = FALSE]
  if (anyDuplicated(paste(d$sector_id, d$period_code))) stop("Duplicate sector-quarter keys: ", id)
  if (any(!d$sector_id %in% obj$sectors$sector_id)) stop("Unknown sector code in ", id)
  if (anyDuplicated(obj$sectors$sector_id)) stop("Duplicate sector metadata in ", id)
  rows <- lapply(seq_len(nrow(obj$sectors)), function(i) {
    s <- obj$sectors[i, , drop = FALSE]; x <- d[d$sector_id == s$sector_id, , drop = FALSE]
    j <- match(grid$period_code, x$period_code)
    cbind(data.frame(dataset = id, sector_id = s$sector_id, sector = s$sector,
                    sector_group = s$sector_group, company_size = s$company_size), grid,
          data.frame(row_present = !is.na(j), rate_raw = x$rate_raw[j],
                     rate_pct = parse_rate(x$rate_raw[j]),
                     value_attribute = x$value_attribute[j], period_status = x$period_status[j]))
  })
  full <- do.call(rbind, rows)
  if (any(full$rate_pct < 0 | full$rate_pct > 100, na.rm = TRUE)) stop("Rate outside [0,100]: ", id)
  full$selected <- (INCLUDE_COMPANY_SIZE | !full$company_size) &
    (is.null(SECTOR_IDS) | full$sector_id %in% SECTOR_IDS)
  quality <- do.call(rbind, lapply(split(full, full$sector_id), function(x) {
    baseline <- x$year <= BASELINE_END
    data.frame(dataset = id, sector_id = x$sector_id[1], sector = x$sector[1],
      selected = x$selected[1], expected_quarters = nrow(grid), rows_present = sum(x$row_present),
      numeric_quarters = sum(!is.na(x$rate_pct)),
      missing_baseline = sum(is.na(x$rate_pct[baseline])),
      missing_total = sum(is.na(x$rate_pct)),
      provisional_quarters = sum(grepl("voorlop|provis", x$period_status, ignore.case = TRUE)),
      unknown_status = sum(is.na(x$period_status) | x$period_status == "unknown"),
      eligible = x$selected[1] && !anyNA(x$rate_pct),
      missing_quarters = paste(x$quarter_label[is.na(x$rate_pct)], collapse = ";"))
  }))
  list(data = full, quality = quality)
}
revision_audit <- function(full, obj, id, cache_dir, out_dir) {
  previous <- file.path(cache_dir, paste0(id, "_audit.rds"))
  vars <- c("rate_pct", "rate_raw", "row_present", "period_status", "value_attribute", "sector")
  cur <- full[, c("sector_id", "period_code", vars)]
  changes <- data.frame(sector_id = character(), period_code = character(), field = character(),
                        previous = character(), current = character())
  note <- "FIRST SNAPSHOT: historical revisions cannot be assessed."
  if (file.exists(previous)) {
    old <- readRDS(previous)
    z <- merge(old$data, cur, by = c("sector_id", "period_code"), all = TRUE,
               suffixes = c("_old", "_new"), sort = TRUE)
    for (v in vars) {
      a <- as.character(z[[paste0(v, "_old")]]); b <- as.character(z[[paste0(v, "_new")]])
      same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
      ix <- which(!same)
      if (length(ix)) changes <- rbind(changes,
        data.frame(z[ix, c("sector_id", "period_code")], field = v, previous = a[ix], current = b[ix]))
    }
    note <- paste0("Compared with snapshot retrieved ", old$retrieved_utc, ". Changed fields: ", nrow(changes),
      ". Measure/sector definitions changed: ", !identical(obj$definitions, old$definitions),
      ". Any metadata changed: ", !identical(obj$metadata, old$metadata),
      ". Metadata differences require review; routine publication updates also count.")
  }
  write_csv(changes, file.path(out_dir, paste0(id, "_revision_changes.csv")))
  snapshot <- list(data = cur, definitions = obj$definitions, metadata = obj$metadata,
                   retrieved_utc = obj$retrieved_utc)
  saveRDS(snapshot, file.path(out_dir, paste0(id, "_snapshot.rds")))
  saveRDS(snapshot, previous)
  note
}
hac_vcov <- function(model) {
  sandwich::NeweyWest(model, lag = min(HAC_LAG, stats::nobs(model) - 1L),
                     prewhite = FALSE, adjust = TRUE)
}
fit_sector <- function(d) {
  d <- d[order(d$quarter_index), ]
  stopifnot(nrow(d) == nrow(quarter_grid()), !anyNA(d$rate_pct))
  d$quarter_factor <- factor(d$quarter, levels = 1:4)
  train <- d[d$year <= BASELINE_END, ]
  fit <- stats::lm(rate_pct ~ quarter_index + quarter_factor, data = train)
  V <- hac_vcov(fit)
  X <- stats::model.matrix(~ quarter_index + quarter_factor, data = d)
  d$baseline_pct <- as.numeric(stats::predict(fit, newdata = d))
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  d$baseline_mean_low95 <- d$baseline_pct - stats::qnorm(0.975) * se
  d$baseline_mean_high95 <- d$baseline_pct + stats::qnorm(0.975) * se
  d$deviation_pp <- d$rate_pct - d$baseline_pct
  d$deviation_relative_pct <- ifelse(d$baseline_pct > 0, 100 * d$deviation_pp / d$baseline_pct, NA_real_)
  # strucchange uses ceiling(n/h)-2 as its maximum admissible break count.
  max_k <- min(MAX_BREAKS, ceiling(nrow(d) / MIN_SEGMENT_QUARTERS) - 2L)
  if (max_k < 1L) stop("Minimum segment size leaves no admissible break model")
  bp <- strucchange::breakpoints(deviation_pp ~ 1, data = d,
             h = MIN_SEGMENT_QUARTERS, breaks = max_k)
  candidates <- 0:max_k
  bic <- as.numeric(stats::AIC(bp, breaks = candidates, k = log(nrow(d))))
  if (any(!is.finite(bic))) stop("Non-finite BIC; degenerate series requires review")
  k <- candidates[which.min(bic)]
  chosen <- if (k) strucchange::breakpoints(bp, breaks = k)$breakpoints else integer()
  segment <- findInterval(seq_len(nrow(d)) - 1L, chosen) + 1L
  d$segment <- segment
  means <- tapply(d$deviation_pp, segment, mean)
  d$segment_mean_deviation_pp <- unname(means[as.character(segment)])
  breaks <- data.frame(break_number = integer(), last_quarter_before = character(),
    first_quarter_after = character(), shift_pp = numeric(),
    ci95_last_before_low = character(), ci95_last_before_high = character(),
    ci_note = character())
  if (k) {
    ci_note <- "HAC; conditional on baseline and selected number of breaks"
    ci <- tryCatch(stats::confint(bp, breaks = k, level = 0.95,
      het.reg = FALSE, het.err = TRUE, vcov. = hac_vcov)$confint,
      error = function(e) { ci_note <<- paste("CI unavailable:", conditionMessage(e)); NULL })
    bounds <- matrix(NA_real_, nrow = k, ncol = 2L)
    if (!is.null(ci)) bounds[,] <- ci[, c(1, 3), drop = FALSE]
    label <- function(j) {
      out <- rep(NA_character_, length(j))
      valid <- is.finite(j) & j >= 1 & j <= nrow(d)
      out[valid] <- d$quarter_label[as.integer(j[valid])]; out
    }
    breaks <- data.frame(break_number = seq_len(k), last_quarter_before = d$quarter_label[chosen],
      first_quarter_after = d$quarter_label[chosen + 1L], shift_pp = diff(as.numeric(means)),
      ci95_last_before_low = label(bounds[, 1]), ci95_last_before_high = label(bounds[, 2]),
      ci_note = ci_note)
    invalid <- is.na(breaks$ci95_last_before_low) | is.na(breaks$ci95_last_before_high)
    breaks$ci_note[invalid] <- paste(breaks$ci_note[invalid], "Bounds unavailable or outside observed window.")
  }
  co <- stats::coef(fit); co_se <- sqrt(diag(V))
  coefficients <- data.frame(term = names(co), estimate = unname(co), se_HAC = unname(co_se),
    low95_HAC = unname(co - stats::qnorm(0.975) * co_se), high95_HAC = unname(co + stats::qnorm(0.975) * co_se))
  summary <- data.frame(n_quarters = nrow(d), n_baseline = nrow(train), selected_breaks = k,
    baseline_slope_pp_per_year = unname(co["quarter_index"]) * 4,
    baseline_slope_low95 = coefficients$low95_HAC[coefficients$term == "quarter_index"] * 4,
    baseline_slope_high95 = coefficients$high95_HAC[coefficients$term == "quarter_index"] * 4,
    mean_2025_deviation_pp = mean(d$deviation_pp[d$year == 2025]),
    baseline_LjungBox_lag4_p = stats::Box.test(stats::residuals(fit), lag = 4, type = "Ljung-Box")$p.value,
    segmented_LjungBox_lag4_p = stats::Box.test(d$deviation_pp - d$segment_mean_deviation_pp,
                                              lag = 4, type = "Ljung-Box")$p.value)
  list(data = d, coefficients = coefficients, summary = summary, breaks = breaks,
       bic = data.frame(n_breaks = candidates, BIC = bic, selected = candidates == k),
       model = list(baseline = fit, baseline_vcov_HAC = V, breakpoints_full = bp, selected_k = k))
}
plot_sector <- function(result) {
  d <- result$data; x <- d$year + (d$quarter - 1) / 4
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  plot(x, d$rate_pct, type = "n", xlab = "Year", ylab = "Sickness absence (%)",
    main = paste(d$dataset[1], d$sector_id[1], d$sector[1]), cex.main = 0.85,
    ylim = range(c(d$rate_pct, d$baseline_mean_low95, d$baseline_mean_high95)))
  polygon(c(x, rev(x)), c(d$baseline_mean_low95, rev(d$baseline_mean_high95)),
          col = "#0072B222", border = NA)
  lines(x, d$baseline_pct, col = "#0072B2", lty = 2, lwd = 1.5)
  lines(x, d$rate_pct, col = "#222222", lwd = 1.5)
  abline(v = BASELINE_END + 1, col = "grey60", lty = 3)
  legend("topleft", c("Observed", "2010-2018 baseline; 95% mean CI"),
         col = c("#222222", "#0072B2"), lty = c(1, 2), bty = "n", cex = 0.8)
  plot(x, d$deviation_pp, type = "o", pch = 16, cex = 0.5,
       xlab = "Year", ylab = "Deviation (percentage points)",
       main = paste("Independent BIC selection:", result$summary$selected_breaks, "breaks"))
  abline(h = 0, col = "grey60", lty = 3)
  for (s in unique(d$segment)) {
    ix <- which(d$segment == s)
    segments(x[min(ix)], d$segment_mean_deviation_pp[ix[1]],
             x[max(ix)], d$segment_mean_deviation_pp[ix[1]], col = "#D55E00", lwd = 2)
  }
  edges <- which(diff(d$segment) != 0)
  if (length(edges)) abline(v = (x[edges] + x[edges + 1L]) / 2, col = "#D55E00", lty = 2)
}
main <- function() {
  need_packages()
  stopifnot(SOURCE_MODE %in% c("refresh", "cache"), START_YEAR == 2010L,
            END_YEAR == 2025L, BASELINE_END == 2018L,
            MIN_SEGMENT_QUARTERS >= 4, MAX_BREAKS >= 1)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(OUT_DIR, "cache")
  dir.create(cache_dir, showWarnings = FALSE)
  run_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC"), "_"), tmpdir = OUT_DIR)
  dir.create(run_dir)
  status <- data.frame(dataset = character(), state = character(), detail = character())
  all_data <- list(); all_quality <- list(); audit_notes <- character()
  for (id in SOURCES) {
    message("Downloading/auditing ", id, " (", SOURCE_MODE, ") ...")
    raw_dir <- file.path(run_dir, paste0("source_", id)); dir.create(raw_dir)
    tryCatch({
      obj <- load_source(id, raw_dir, cache_dir)
      jsonlite::write_json(obj$metadata, file.path(raw_dir, "metadata.json"), pretty = TRUE, na = "null")
      capture.output(str(obj$metadata$properties), file = file.path(raw_dir, "metadata_summary.txt"))
      a <- audit_source(obj, id)
      revision <- revision_audit(a$data, obj, id, cache_dir, run_dir)
      all_data[[id]] <- a$data; all_quality[[id]] <- a$quality
      status <- rbind(status, data.frame(dataset = id, state = "available", detail = paste("Retrieved", obj$retrieved_utc)))
      audit_notes <- c(audit_notes, paste(id, revision))
    }, error = function(e) {
      status <<- rbind(status, data.frame(dataset = id, state = "FAILED", detail = conditionMessage(e)))
      message("Source unavailable: ", id, ": ", conditionMessage(e))
    })
  }
  write_csv(status, file.path(run_dir, "source_status.csv"))
  writeLines(c("PRE-MODEL AUDIT", paste("Mode:", SOURCE_MODE), audit_notes,
    "Known classification boundary: CBS enterprise-group/register changes from 2010.",
    "Read source_*/metadata.json for current descriptions, definitions and release notes.",
    "Final figures can still be revised; first snapshots do not establish revision history.",
    "Missing series are excluded as whole series; no interpolation or annual-to-quarter allocation.",
    "AZW broad health/welfare and CBS SBI Q overlap. Rates must not be summed.",
    "HAC bands are baseline mean intervals, not future-observation prediction intervals.",
    "Breaks are exploratory level shifts. Baseline and break-count uncertainty is not fully propagated.",
    "Ljung-Box p-values are descriptive residual diagnostics, not post-selection formal tests."),
    file.path(run_dir, "READ_BEFORE_ANALYSIS.txt"))
  if (!length(all_data)) stop("No usable source; see ", file.path(run_dir, "source_status.csv"), call. = FALSE)
  full <- do.call(rbind, all_data); quality <- do.call(rbind, all_quality)
  write_csv(full, file.path(run_dir, "quarterly_observations_2010_2025.csv"))
  write_csv(quality, file.path(run_dir, "quality_checks.csv"))
  print(status, row.names = FALSE)
  message(sum(quality$eligible), " eligible series; ",
          sum(quality$selected & !quality$eligible), " selected series with gaps.")
  if (any(quality$selected & !quality$eligible))
    print(quality[quality$selected & !quality$eligible, ], row.names = FALSE)
  if (STRICT_SOURCES && any(status$state == "FAILED")) stop("Incomplete sources; audit saved in ", run_dir)
  eligible <- paste(quality$dataset[quality$eligible], quality$sector_id[quality$eligible], sep = "::")
  results <- list(); fit_errors <- data.frame(series = character(), error = character())
  split_data <- split(full, paste(full$dataset, full$sector_id, sep = "::"))
  for (key in eligible) {
    results[[key]] <- tryCatch(fit_sector(split_data[[key]]), error = function(e) {
      fit_errors <<- rbind(fit_errors, data.frame(series = key, error = conditionMessage(e))); NULL
    })
  }
  results <- Filter(Negate(is.null), results)
  write_csv(fit_errors, file.path(run_dir, "model_errors.csv"))
  if (!length(results)) stop("No eligible model fitted. See quality_checks.csv and model_errors.csv in ", run_dir)
  combine <- function(name) do.call(rbind, lapply(results, function(r) {
    x <- r[[name]]
    identity <- r$data[rep(1L, nrow(x)), c("dataset", "sector_id", "sector"), drop = FALSE]
    if (name == "data") x else cbind(identity, x)
  }))
  for (name in c("data", "coefficients", "summary", "breaks", "bic")) {
    filename <- switch(name, data = "baseline_and_deviations.csv", coefficients = "baseline_coefficients.csv",
      summary = "sector_summary.csv", breaks = "breakpoints.csv", bic = "breakpoint_BIC.csv")
    write_csv(combine(name), file.path(run_dir, filename))
  }
  saveRDS(lapply(results, `[[`, "model"), file.path(run_dir, "fitted_models.rds"))
  if (MAKE_PLOTS) {
    grDevices::pdf(file.path(run_dir, "sector_plots.pdf"), width = 11, height = 8.5, onefile = TRUE)
    tryCatch(invisible(lapply(results, plot_sector)), finally = grDevices::dev.off())
  }
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  dput(list(SOURCES = SOURCES, SOURCE_MODE = SOURCE_MODE, SECTOR_IDS = SECTOR_IDS,
    INCLUDE_COMPANY_SIZE = INCLUDE_COMPANY_SIZE, MIN_SEGMENT_QUARTERS = MIN_SEGMENT_QUARTERS,
    MAX_BREAKS = MAX_BREAKS, HAC_LAG = HAC_LAG), file = file.path(run_dir, "configuration.R"))
  message("Saved ", length(results), " sector models in: ", normalizePath(run_dir, winslash = "/"))
  incomplete <- any(status$state == "FAILED") || any(quality$selected & !quality$eligible) || nrow(fit_errors) > 0
  if (incomplete) warning("INCOMPLETE ANALYSIS: some requested sources/series were unavailable. Read source_status.csv, quality_checks.csv and model_errors.csv.", call. = FALSE)
  invisible(list(run_dir = run_dir, status = status, quality = quality, results = results))
}
if (isTRUE(RUN_ANALYSIS)) analysis <- main()
