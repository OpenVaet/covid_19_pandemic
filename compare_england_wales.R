# England-Wales NHS sickness absence comparison
# Version 1.1; adds matched England-style figures; R >= 4.3. No Pandoc, R Markdown or live website access required.
# Install once: install.packages(c("strucchange", "sandwich", "jsonlite", "ggplot2", "ggrepel"))
# First run nhs_sickness_rates.R and wales_sickness_rates.R, then:
# source("compare_england_wales.R", encoding = "UTF-8")
# Use RStudio's Source button; do not continue line by line after an error.
#
# Default: latest successful run in each results root below.
# To select exact inputs, set ENGLAND_RUN and WALES_RUN to a run folder OR ZIP.
# Windows paths: use forward slashes, e.g. "C:/Research/run_20261022.zip".
# Outputs: timestamped folder, self-contained report.html, CSVs, PNG/PDF plots,
# fitted models, input copies, checks, configuration, source metadata and hashes.
#
# The original England script uses elapsed calendar days and the full source
# window. This script deliberately REFITS BOTH baselines using month_index,
# restricts BOTH series to identical months, and re-estimates independent breaks.
# Baseline: OLS linear trend + month-of-year, January 2010-December 2018.
# Outcome: percent of available FTE calendar days lost; NOT percent of staff ill.
# Breaks: unsmoothed deviation level shifts; BIC over all feasible counts, including
# zero. Dates label FIRST month after each break. No country supplies the other's
# dates. h = 12/18/24 months shows the sensitivity to minimum segment duration.
# NW lag-12 intervals are conditional on the fitted baseline/chosen break count.
# HAC does not make BIC serially robust or account for model-selection uncertainty.
# Annual summaries are unweighted means of monthly rates, NOT weighted official
# annual absence rates. Partial years are labelled, and never compared as full years.
# All-cause absence, staff mix and recording practices cannot identify causal
# effects of infections, vaccination, mandates or any other intervention.

# ---- Configuration ---------------------------------------------------------
ENGLAND_ROOT <- "england_sickness_results"
WALES_ROOT <- "wales_sickness_results"
ENGLAND_RUN <- getOption("comparison.england_run", NULL)
WALES_RUN <- getOption("comparison.wales_run", NULL)
OUT_DIR <- getOption("comparison.output_dir", "england_wales_comparison_results")
START_YEAR <- 2010L
END_YEAR <- 2025L
BASELINE_END <- 2018L
HAC_LAG <- 12L
SEGMENT_MONTHS <- c(12L, 18L, 24L)
MAX_BREAKS <- NULL  # NULL searches all feasible counts; a numeric cap is recorded.
INCLUDE_WALES_STAFF_SUPPLEMENT <- TRUE
RUN_ANALYSIS <- getOption("comparison.run_analysis", TRUE)

# England-style figures supplement the existing diagnostic figures.
MAKE_ENGLAND_STYLE_PLOTS <- TRUE
# Shared UK reference events; add verified country-specific rows if desired.
# Set STYLE_EVENTS <- STYLE_EVENTS[0, ] to hide all event markers.
STYLE_EVENTS <- data.frame(date=as.Date(c("2020-01-31","2020-12-08")),
  label=c("First UK cases","UK vaccine rollout"),
  source=c("https://www.gov.uk/government/news/cmo-confirms-cases-of-coronavirus-in-england",
    "https://www.gov.uk/government/news/uk-marks-one-year-since-deploying-worlds-first-covid-19-vaccine"))
# Use the same existing sensitivity fit for both England-style panels.
STYLE_MIN_SEGMENT_MONTHS <- 24L  # Choose 12, 18 or 24; must be in SEGMENT_MONTHS.
# All rolling calculations use only common-window data, including the edges.

write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
text_file <- function(path) paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
observation_file <- "monthly_observations_2010_2025.csv"
is_complete_run <- function(path) {
  if (!dir.exists(path) || file.exists(file.path(path, "FAILED.txt")) ||
      !file.exists(file.path(path, observation_file))) return(FALSE)
  if (file.exists(file.path(path, "RUN_COMPLETE"))) return(TRUE)
  # Compatibility with the existing Wales script, which has no completion marker.
  if (!all(file.exists(file.path(path, c("report.html", "source_status.csv", "model_errors.csv"))))) return(FALSE)
  errors <- tryCatch(read.csv(file.path(path, "model_errors.csv")), error = function(e) NULL)
  !is.null(errors) && nrow(errors) == 0L
}
resolve_run <- function(spec, root, country) {
  if (is.null(spec)) {
    if (!dir.exists(root)) stop(country, ": results root not found: ", root,
      ". Run the country script first, or set ", toupper(country), "_RUN to a folder or ZIP.")
    dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[grepl("^run_", basename(dirs))]
    dirs <- dirs[vapply(dirs, is_complete_run, logical(1))]
    if (!length(dirs)) stop(country, ": no successful compatible run found in ", root)
    spec <- dirs[order(basename(dirs), decreasing = TRUE)][1]
  }
  original <- normalizePath(spec, winslash = "/", mustWork = TRUE)
  cleanup <- NULL
  if (!dir.exists(original)) {
    if (!grepl("\\.zip$", original, ignore.case = TRUE)) stop("Input must be a run folder or ZIP: ", original)
    entries <- utils::unzip(original, list = TRUE)$Name
    # Never extract absolute paths, parent traversal or ambiguous multiple runs.
    if (any(grepl("(^[/\\\\])|(^[A-Za-z]:)|((^|[/\\\\])\\.\\.([/\\\\]|$))", entries))) stop("Unsafe ZIP entry path")
    target <- entries[basename(entries) == observation_file]
    if (length(target) != 1L) stop(country, ": ZIP must contain exactly one run's observations")
    prefix <- dirname(target)
    cleanup <- tempfile("nhs_comparison_input_"); dir.create(cleanup)
    utils::unzip(original, exdir = cleanup)
    path <- file.path(cleanup, prefix)
  } else path <- original
  if (!is_complete_run(path)) {
    if (!is.null(cleanup)) unlink(cleanup, recursive = TRUE)
    stop(country, ": input is incomplete or contains model errors: ", original)
  }
  list(path = path, original = original, cleanup = cleanup)
}
prepare_series <- function(raw, country, staff = "XXX", end_date = as.Date(sprintf("%d-12-01", END_YEAR))) {
  required <- c("organisation_id", "staff_id", "date", "rate_pct")
  if (!all(required %in% names(raw))) stop(country, ": unexpected observations schema")
  org <- if (country == "England") "E" else "W"
  d <- raw[!is.na(raw$organisation_id) & !is.na(raw$staff_id) & raw$organisation_id == org & raw$staff_id == staff, , drop = FALSE]
  if (!nrow(d)) stop(country, ": no published national series for staff code ", staff)
  d$date <- as.Date(as.character(d$date), format = "%Y-%m-%d")
  if (anyNA(d$date) || any(format(d$date, "%d") != "01")) stop(country, ": invalid monthly date")
  d$rate_pct <- suppressWarnings(as.numeric(as.character(d$rate_pct)))
  d <- d[d$date >= as.Date(sprintf("%d-01-01", START_YEAR)) & d$date <= end_date, , drop = FALSE]
  if (!nrow(d) || anyDuplicated(d$date)) stop(country, ": empty or duplicate monthly series")
  d <- d[order(d$date), ]
  if (any(!is.finite(d$rate_pct)) || any(d$rate_pct < 0 | d$rate_pct > 100)) stop(country, ": missing/invalid rate")
  expected <- seq(as.Date(sprintf("%d-01-01", START_YEAR)), max(d$date), by = "month")
  if (!identical(as.character(d$date), as.character(expected))) stop(country, ": missing month(s) inside the observation window")
  baseline <- seq(as.Date(sprintf("%d-01-01", START_YEAR)), as.Date(sprintf("%d-12-01", BASELINE_END)), by = "month")
  if (!all(baseline %in% d$date)) stop(country, ": incomplete baseline")
  d$country <- country
  d$year <- as.integer(format(d$date, "%Y")); d$month <- as.integer(format(d$date, "%m"))
  d$month_index <- (d$year - START_YEAR)*12L + d$month - 1L
  d$month_label <- format(d$date, "%Y-%m")
  d
}
hac_vcov <- function(model) {
  effective_n <- NROW(sandwich::estfun(model))
  sandwich::NeweyWest(model, lag = min(HAC_LAG, effective_n - 1L),
    prewhite = FALSE, adjust = TRUE)
}
fit_baseline <- function(d) {
  d <- d[order(d$month_index), ]
  stopifnot(all(is.finite(d$rate_pct)))
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
  max_k <- ceiling(nrow(d) / h) - 2L
  if (!is.null(MAX_BREAKS)) max_k <- min(MAX_BREAKS, max_k)
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

plot_comparison <- function(series, breaks) {
  colors <- c(England = "#0067A5", Wales = "#C55223")
  graphics::par(mfrow = c(3, 1), mar = c(3.5, 4.5, 2.8, 1.2), mgp = c(2.6, .8, 0))
  x <- series$England$year + (series$England$month - 1)/12
  all_rates <- unlist(lapply(series, function(d) c(d$rate_pct, d$baseline_pct)))
  graphics::plot(range(x), range(all_rates), type = "n", xlab = "", ylab = "Absence (%)",
    main = "National all-staff monthly rates and 2010-2018 baselines")
  graphics::abline(v = BASELINE_END + 1, lty = 3, col = "grey65")
  for (country in names(series)) {
    d <- series[[country]]
    graphics::lines(x, d$rate_pct, col = colors[country], lwd = 1.4)
    graphics::lines(x, d$baseline_pct, col = colors[country], lty = 2, lwd = 1.1)
  }
  graphics::legend("topleft", c("England observed", "Wales observed", "Separate baseline (dashed)"),
    col = c(colors, "grey35"), lty = c(1, 1, 2), bty = "n", cex = .8)
  post <- series$England$year > BASELINE_END
  deviations <- unlist(lapply(series, function(d) d$deviation_pp[post]))
  graphics::plot(range(x[post]), range(c(0, deviations)), type = "n", xlab = "", ylab = "Deviation (pp)",
    main = "Observed minus each country's own baseline; unsmoothed")
  graphics::abline(h = 0, lty = 3, col = "grey65")
  for (country in names(series)) graphics::lines(x[post], series[[country]]$deviation_pp[post], col = colors[country], lwd = 1.6)
  graphics::legend("topleft", names(colors), col = colors, lty = 1, lwd = 1.6, bty = "n", cex = .8)
  hs <- sort(SEGMENT_MONTHS); positions <- seq_len(2*length(hs))
  labels <- unlist(lapply(hs, function(h) paste(c("England", "Wales"), "h =", h)))
  graphics::par(mar = c(3.5, 7, 2.8, 1.2))
  graphics::plot(range(x[post]), c(.5, length(positions)+.5), type = "n", xlab = "First month after break (post-2018 shown)",
    ylab = "", yaxt = "n", main = "Independently selected break dates: minimum-duration sensitivity")
  graphics::axis(2, positions, labels = labels, las = 1, cex.axis = .7)
  for (i in seq_along(hs)) for (j in 1:2) {
    country <- names(colors)[j]; y <- (i-1)*2+j
    b <- breaks[breaks$country == country & breaks$min_segment_months == hs[i] & breaks$first_month_after >= sprintf("%d-01", BASELINE_END+1L), ]
    graphics::abline(h = y, col = "grey90")
    if (nrow(b)) graphics::points(x[match(b$first_month_after, series$England$month_label)], rep(y, nrow(b)), pch = 19, col = colors[country])
  }
}
embed_image <- function(path, alt) {
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  paste0("<img alt='", html_escape(alt), "' src='data:image/png;base64,", jsonlite::base64_enc(bytes), "'>")
}
main <- function() {
  pkgs <- c("strucchange", "sandwich", "jsonlite")
  if (MAKE_ENGLAND_STYLE_PLOTS) pkgs <- c(pkgs, "ggplot2", "ggrepel")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  stopifnot(START_YEAR <= BASELINE_END, BASELINE_END < END_YEAR,
    all(SEGMENT_MONTHS >= 2L), all(SEGMENT_MONTHS == as.integer(SEGMENT_MONTHS)))
  if (MAKE_ENGLAND_STYLE_PLOTS && !STYLE_MIN_SEGMENT_MONTHS %in% SEGMENT_MONTHS)
    stop("STYLE_MIN_SEGMENT_MONTHS must be one of SEGMENT_MONTHS")
  if (!is.null(MAX_BREAKS) && (length(MAX_BREAKS) != 1 || !is.finite(MAX_BREAKS) || MAX_BREAKS < 1 || MAX_BREAKS %% 1 != 0))
    stop("MAX_BREAKS must be NULL or a positive integer")
  england <- resolve_run(ENGLAND_RUN, ENGLAND_ROOT, "England")
  on.exit(if (!is.null(england$cleanup)) unlink(england$cleanup, recursive = TRUE), add = TRUE)
  wales <- resolve_run(WALES_RUN, WALES_ROOT, "Wales")
  on.exit(if (!is.null(wales$cleanup)) unlink(wales$cleanup, recursive = TRUE), add = TRUE)
  inputs <- list(England = england, Wales = wales)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  run_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = OUT_DIR)
  dir.create(run_dir); dir.create(file.path(run_dir, "source")); dir.create(file.path(run_dir, "plots"))
  complete <- FALSE
  on.exit(if (!complete) writeLines("Comparison failed. Do not use partial results.", file.path(run_dir, "FAILED.txt")), add = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
  dput(list(england_input = england$original, wales_input = wales$original,
    start_year = START_YEAR, end_year = END_YEAR, baseline_end = BASELINE_END,
    HAC_lag = HAC_LAG, min_segments = SEGMENT_MONTHS, max_breaks = MAX_BREAKS,
    style_min_segment_months = STYLE_MIN_SEGMENT_MONTHS, style_events = STYLE_EVENTS,
    make_england_style_plots = MAKE_ENGLAND_STYLE_PLOTS,
    convention = "first month after", model = "independent level breaks on monthly baseline deviations"),
    file.path(run_dir, "configuration.R"))
  raw <- list(); national <- list(); source_rows <- list(); audit_notes <- list()
  for (country in names(inputs)) {
    input <- inputs[[country]]
    dest <- file.path(run_dir, "source", country); dir.create(dest)
    archive <- c(observation_file, "source_status.csv", "quality_checks.csv", "revision_changes.csv",
      "known_definition_changes.csv", "READ_BEFORE_ANALYSIS.txt", "configuration.R")
    archive <- archive[file.exists(file.path(input$path, archive))]
    if (!all(file.copy(file.path(input$path, archive), dest))) stop("Could not archive inputs")
    path <- file.path(input$path, observation_file)
    raw[[country]] <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    national[[country]] <- prepare_series(raw[[country]], country)
    d <- national[[country]]
    source_rows[[country]] <- data.frame(country = country, input_run = input$original,
      observation_md5 = unname(tools::md5sum(path)), first_month = d$month_label[1],
      last_month = tail(d$month_label, 1), available_months = nrow(d), missing_internal_months = 0L)
    notes_path <- file.path(input$path, "READ_BEFORE_ANALYSIS.txt")
    audit_notes[[country]] <- if (file.exists(notes_path)) text_file(notes_path) else "No source notes exported. Review source documentation."
  }
  source_info <- do.call(rbind, source_rows)
  # Validate each full input span before truncation. An inner join must not hide holes.
  common_end <- min(national$England$date[nrow(national$England)], national$Wales$date[nrow(national$Wales)])
  expected <- seq(as.Date(sprintf("%d-01-01", START_YEAR)), common_end, by = "month")
  if (common_end <= as.Date(sprintf("%d-12-01", BASELINE_END)))
    stop("The common window contains no post-baseline observations")
  for (country in names(national)) {
    national[[country]] <- national[[country]][national[[country]]$date <= common_end, ]
    stopifnot(identical(as.character(national[[country]]$date), as.character(expected)))
  }
  source_info$common_end <- format(common_end, "%Y-%m")
  source_info$common_months <- length(expected)
  source_info$months_excluded_after_common_end <- source_info$available_months - length(expected)
  write_csv(source_info, file.path(run_dir, "quality_checks.csv"))
  notes <- c(
    "Published national all-staff totals are used; staff-group percentages are never summed or averaged into national totals.",
    "Both countries are re-fitted on identical monthly dates: own 2010-2018 linear trend plus month-of-year effects. Exported country-script baselines/breakpoints are not reused.",
    "The comparison England baseline uses a monthly index; the original England analysis retains elapsed calendar days. Its original break window also starts before 2010. Small numerical/date differences are therefore expected.",
    "Breakpoints are estimated independently on unsmoothed deviations over the entire common window, including the baseline years. Early baseline-period breaks are reported in the full table.",
    "BIC selects zero through all feasible break counts unless MAX_BREAKS is explicitly capped. A 24-month duration restriction can move neighbouring estimated break dates; it does not identify the first abnormal month.",
    "Break dates label the first month of the new segment. Last-before dates and approximate conditional HAC intervals are also exported.",
    "Baseline mean intervals use Newey-West covariance with lag 12. Neither baseline nor break-date intervals include model-selection uncertainty. Serial correlation may still affect BIC.",
    "Annual summaries are unweighted means of available monthly rates. Matched partial-year rows cover exactly the same months. No extrapolation or missing-month imputation is performed.",
    "Rates measure days lost, not numbers of sick workers. Absolute excess is in percentage points; relative excess is 100 times mean excess divided by mean predicted baseline.",
    "Differences in baseline levels may reflect staff mix, employer coverage and recording. Administrative/occupation changes and revisions require source review; the comparison does not establish definition stability.",
    "Correlation is descriptive for autocorrelated series; no independent-month p-value is claimed. Concordant all-cause absence does not identify a cause or intervention effect.",
    "Input audits and revision files are copied with the inputs. This comparison does not refresh the original sources or reconstruct missing publication vintages.")
  writeLines(notes, file.path(run_dir, "READ_BEFORE_ANALYSIS.txt"))
  baselines <- lapply(national, fit_baseline)
  fitted <- lapply(baselines, `[[`, "data")
  tag <- function(x, country) cbind(country = rep(country, nrow(x)), x)
  collect <- function(key) do.call(rbind, lapply(names(baselines), function(cn) tag(baselines[[cn]][[key]], cn)))
  annual_long <- collect("annual")
  annual_long$complete_year <- annual_long$n_months == 12L
  write_csv(annual_long, file.path(run_dir, "annual_mean_monthly_rates.csv"))
  write_csv(collect("coefficients"), file.path(run_dir, "baseline_coefficients.csv"))
  write_csv(collect("summary"), file.path(run_dir, "baseline_diagnostics.csv"))
  common_columns <- c("country", "date", "month_label", "year", "month", "month_index", "rate_pct", "baseline_pct",
    "baseline_mean_low95", "baseline_mean_high95", "deviation_pp", "deviation_relative_pct")
  write_csv(do.call(rbind, lapply(fitted, function(d) d[, common_columns])), file.path(run_dir, "baseline_and_deviations.csv"))
  matched <- data.frame(month = fitted$England$month_label, England_rate_pct = fitted$England$rate_pct,
    Wales_rate_pct = fitted$Wales$rate_pct, England_baseline_pct = fitted$England$baseline_pct,
    Wales_baseline_pct = fitted$Wales$baseline_pct, England_deviation_pp = fitted$England$deviation_pp,
    Wales_deviation_pp = fitted$Wales$deviation_pp)
  matched$Wales_minus_England_rate_pp <- matched$Wales_rate_pct - matched$England_rate_pct
  matched$Wales_minus_England_deviation_pp <- matched$Wales_deviation_pp - matched$England_deviation_pp
  write_csv(matched, file.path(run_dir, "matched_monthly_comparison.csv"))
  annual <- merge(baselines$England$annual, baselines$Wales$annual, by = c("year", "n_months"), suffixes = c("_England", "_Wales"))
  annual$complete_year <- annual$n_months == 12L
  write_csv(annual, file.path(run_dir, "matched_annual_comparison.csv"))
  breaks_fits <- list(); breaks <- segments <- bics <- sensitivities <- list()
  for (country in names(fitted)) for (h in unique(SEGMENT_MONTHS)) {
    message("Fitting independent breaks: ", country, ", minimum segment ", h, " months")
    fit <- fit_breaks(fitted[[country]], h)
    key <- paste(country, h, sep = "_")
    breaks_fits[[key]] <- fit
    breaks[[key]] <- tag(fit$breaks, country); segments[[key]] <- tag(fit$segments, country)
    bics[[key]] <- tag(fit$bic, country); sensitivities[[key]] <- tag(fit$summary, country)
  }
  breaks <- do.call(rbind, breaks); segments <- do.call(rbind, segments)
  bics <- do.call(rbind, bics); sensitivities <- do.call(rbind, sensitivities)
  sensitivities$first_months_after_2018 <- vapply(seq_len(nrow(sensitivities)), function(i) {
    b <- breaks[breaks$country == sensitivities$country[i] & breaks$min_segment_months == sensitivities$min_segment_months[i], ]
    paste(b$first_month_after[b$first_month_after >= sprintf("%d-01", BASELINE_END+1L)], collapse = "; ")
  }, character(1))
  write_csv(breaks, file.path(run_dir, "breakpoints_all_sensitivities.csv"))
  write_csv(segments, file.path(run_dir, "segments_all_sensitivities.csv"))
  write_csv(bics, file.path(run_dir, "breakpoint_BIC.csv"))
  write_csv(sensitivities, file.path(run_dir, "breakpoint_sensitivity.csv"))
  post <- fitted$England$year > BASELINE_END
  correlations <- data.frame(window = c("All common months", "Post-baseline common months"),
    n_months = c(nrow(matched), sum(post)),
    Pearson_r_deviations = c(cor(matched$England_deviation_pp, matched$Wales_deviation_pp),
      cor(matched$England_deviation_pp[post], matched$Wales_deviation_pp[post])))
  write_csv(correlations, file.path(run_dir, "descriptive_correlations.csv"))
  # Optional Wales-only staff supplement uses each series' full available span.
  # It is explicitly separate from the matched all-staff country comparison.
  staff <- list(); staff_errors <- data.frame(staff_id = character(), note = character())
  if (INCLUDE_WALES_STAFF_SUPPLEMENT) {
    staff_ids <- unique(as.character(raw$Wales$staff_id[raw$Wales$organisation_id == "W"]))
    staff_ids <- staff_ids[!is.na(staff_ids) & nzchar(staff_ids)]
    for (id in staff_ids) {
      result <- tryCatch({
        d <- prepare_series(raw$Wales, "Wales", staff = id)
        fit <- fit_baseline(d); a <- fit$annual[fit$annual$n_months == 12L, ]
        if (!nrow(a)) stop("No complete calendar year")
        a <- a[which.max(a$year), ]
        data.frame(staff_id = id, staff_group = if ("staff_group" %in% names(d)) d$staff_group[1] else id, a)
      }, error = identity)
      if (inherits(result, "error")) staff_errors <- rbind(staff_errors, data.frame(staff_id = id, note = conditionMessage(result)))
      else staff[[id]] <- result
    }
  }
  staff <- if (length(staff)) do.call(rbind, staff) else data.frame()
  write_csv(staff, file.path(run_dir, "wales_staff_latest_full_year.csv"))
  write_csv(staff_errors, file.path(run_dir, "wales_staff_exclusions.csv"))
  saveRDS(list(baselines = baselines, break_fits = breaks_fits, inputs = source_info), file.path(run_dir, "fitted_models.rds"))
  png_path <- file.path(run_dir, "plots", "england_wales_comparison.png")
  grDevices::png(png_path, width = 1800, height = 1800, res = 160)
  tryCatch(plot_comparison(fitted, breaks), finally = grDevices::dev.off())
  grDevices::pdf(file.path(run_dir, "england_wales_comparison.pdf"), width = 11, height = 11)
  tryCatch(plot_comparison(fitted, breaks), finally = grDevices::dev.off())
  style_paths <- character()
  if (MAKE_ENGLAND_STYLE_PLOTS) {
    style_objects <- lapply(names(fitted), function(country) {
      b <- breaks_fits[[paste(country, STYLE_MIN_SEGMENT_MONTHS, sep = "_")]]$breaks
      style_prepare(fitted[[country]], baselines[[country]]$fit, b, country, "All staff", STYLE_MIN_SEGMENT_MONTHS)
    })
    names(style_objects) <- names(fitted)
    common_x <- range(expected)
    common_y <- style_axis_limits(style_objects)
    style_plots <- lapply(style_objects, style_plot, x_limits = common_x, y_limits = common_y)
    for (country in names(style_plots)) {
      style_paths[country] <- file.path(run_dir, "plots", paste0(tolower(country), "_england_style.png"))
      style_save_plot(style_plots[[country]], style_paths[country],
        file.path(run_dir, paste0(tolower(country), "_england_style.pdf")))
    }
    # Align panel columns exactly, including axis-label and margin widths.
    grobs <- lapply(style_plots, ggplot2::ggplotGrob)
    shared_widths <- do.call(grid::unit.pmax, lapply(grobs, `[[`, "widths"))
    for (i in seq_along(grobs)) grobs[[i]]$widths <- shared_widths
    draw_pair <- function() {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
      for (i in seq_along(grobs)) {
        grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
        grid::grid.draw(grobs[[i]]); grid::popViewport()
      }
      grid::popViewport()
    }
    pair_path <- file.path(run_dir, "plots", "england_wales_england_style_matched.png")
    grDevices::png(pair_path, width = 3960, height = 4400, res = 220)
    tryCatch(draw_pair(), finally = grDevices::dev.off())
    grDevices::pdf(file.path(run_dir, "england_wales_england_style_matched.pdf"), width = 18, height = 20)
    tryCatch(draw_pair(), finally = grDevices::dev.off())
    write_csv(do.call(rbind, lapply(style_objects, `[[`, "data")), file.path(run_dir, "england_style_plot_series.csv"))
    write_csv(do.call(rbind, lapply(style_objects, function(z) {
      cbind(country = rep(z$country, nrow(z$trends)), z$trends)
    })), file.path(run_dir, "england_style_descriptive_trends.csv"))
    write_csv(STYLE_EVENTS, file.path(run_dir, "england_style_event_annotations.csv"))
    write_csv(data.frame(first_month = format(common_x[1], "%Y-%m"), last_month = format(common_x[2], "%Y-%m"),
      y_min_100k = common_y[1], y_max_100k = common_y[2], min_segment_months = STYLE_MIN_SEGMENT_MONTHS,
      rolling_before = 3L, rolling_after = 2L, require_complete_window = TRUE), file.path(run_dir, "england_style_plot_settings.csv"))
    notes <- c(notes, style_notes(),
      "The two England-style comparison panels have identical date limits, y limits, smoothing windows, event annotations and minimum segment duration. Each uses its own matched-period fit and break dates, so they may differ from a country report fitted through a later endpoint.")
    writeLines(notes, file.path(run_dir, "READ_BEFORE_ANALYSIS.txt"))
  }
  last <- tail(annual, 1)
  pre <- vapply(fitted, function(d) mean(d$rate_pct[d$year <= BASELINE_END]), numeric(1))
  annual_table <- data.frame(Year = annual$year, Months = annual$n_months,
    Period = ifelse(annual$n_months == 12L, "Full year", paste0("Jan-", format(common_end, "%b"))),
    England_rate_pct = annual$mean_monthly_rate_pct_England, Wales_rate_pct = annual$mean_monthly_rate_pct_Wales,
    England_excess_pp = annual$mean_deviation_pp_England, Wales_excess_pp = annual$mean_deviation_pp_Wales,
    England_relative_pct = annual$relative_to_baseline_pct_England, Wales_relative_pct = annual$relative_to_baseline_pct_Wales)
  post_annual <- annual_table[annual_table$Year > BASELINE_END, ]
  complete_post <- post_annual[post_annual$Months == 12L, ]
  findings <- c(sprintf("In the 2010-2018 baseline period, Wales averaged %.2f%% and England %.2f%%: a %.2f percentage-point level difference.", pre["Wales"], pre["England"], pre["Wales"]-pre["England"]),
    sprintf("Over the %d post-baseline common months, the correlation of unsmoothed deviations from the separate baselines is %.3f (descriptive, not a causal test).", sum(post), correlations$Pearson_r_deviations[2]),
    sprintf("For %d (%d matched months), England averaged %.2f%%, %.2f pp above baseline (%+.1f%% relative); Wales averaged %.2f%%, %.2f pp above baseline (%+.1f%% relative).", last$year, last$n_months,
      last$mean_monthly_rate_pct_England, last$mean_deviation_pp_England, last$relative_to_baseline_pct_England,
      last$mean_monthly_rate_pct_Wales, last$mean_deviation_pp_Wales, last$relative_to_baseline_pct_Wales))
  if (nrow(complete_post)) findings <- c(findings, sprintf("Among complete post-baseline years, the largest average excess occurs in %d for England (%.2f pp) and %d for Wales (%.2f pp).",
    complete_post$Year[which.max(complete_post$England_excess_pp)], max(complete_post$England_excess_pp),
    complete_post$Year[which.max(complete_post$Wales_excess_pp)], max(complete_post$Wales_excess_pp)))
  writeLines(findings, file.path(run_dir, "findings.txt"))
  parts <- c("<!doctype html><html lang='en'><head><meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>England-Wales NHS sickness absence comparison</title>",
    "<style>body{font:16px/1.6 system-ui,sans-serif;max-width:1160px;margin:40px auto;padding:0 24px;color:#203040}h1{font-size:34px;line-height:1.15}h2{margin-top:40px;color:#144668}table{border-collapse:collapse;font-size:14px}td,th{padding:9px;border-bottom:1px solid #d9e1e8;text-align:left;vertical-align:top}th{background:#eaf1f7}.scroll{overflow:auto}.note{padding:18px 24px;background:#f0f5f8;border-left:4px solid #0067A5}img{width:100%;height:auto}li{margin:.5em 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}a{color:#0067A5}small{color:#526473}@media print{body{max-width:none;margin:0}h2{break-after:avoid}img,tr{break-inside:avoid}}</style></head><body>",
    "<h1>England and Wales: NHS sickness absence</h1>",
    paste0("<p>Matched monthly analysis: ", format(min(expected), "%Y-%m"), " to ", format(common_end, "%Y-%m"),
      " (", length(expected), " months). National all-staff totals.</p>"),
    "<div class='note'><strong>What the data show</strong><ul>", paste0("<li>", html_escape(findings), "</li>"), "</ul></div>",
    if (MAKE_ENGLAND_STYLE_PLOTS) c(
      "<h2>England-style plots over the same period</h2>",
      paste0("<p>Both panels use ", format(min(expected), "%Y-%m"), " to ", format(common_end, "%Y-%m"),
        ", identical vertical and horizontal scales, and a ", STYLE_MIN_SEGMENT_MONTHS,
        "-month minimum segment length. Breaks are independently estimated for each country.</p>"),
      "<p>Grey: raw monthly rate; orange: six-month average; blue: baseline and smoothed pointwise OLS prediction envelope; pink: break dates. Grey labels show relative deviation from the rolling baseline; coloured straight lines are descriptive segment/tail trends.</p>",
      style_embed(pair_path, "England and Wales in the original England plot style, aligned over identical months and axes"),
      "<p>The blue envelope reproduces the original presentation by smoothing pointwise OLS 95% prediction-limit endpoints. It is not a formal interval for the six-month average. The HAC baseline mean intervals remain in the other model outputs. Break labels use the first month after each break.</p>"
    ) else character(),
    "<h2>Observed rates and departures from baseline</h2>", embed_image(png_path, "Matched rates, baseline deviations and independent breakpoint sensitivity"),
    "<p>The higher Welsh rate is partly a pre-existing level difference. Percentage-point excess and relative excess answer different questions: the same absolute excess is a smaller percentage increase from a higher baseline.</p>",
    "<h2>Matched annual summaries</h2><p>Means of monthly rates; excess = observed minus own projected baseline. Relative excess uses the ratio of annual monthly means. Partial-year rows contain identical months for both countries.</p>",
    html_table(post_annual, digits = 3),
    "<details><summary>All years, including the baseline period</summary>", html_table(annual_table, digits = 3), "</details>",
    "<h2>Independently estimated breakpoints</h2><p>No English dates are imposed on Wales. The total count includes pre-2019 breaks; the last column below isolates post-2018 dates. All feasible counts are searched by default, including zero.</p>",
    html_table(sensitivities[, c("country", "min_segment_months", "selected_breaks", "at_candidate_limit", "first_months_after_2018")]),
    "<p>Do not read the earliest h=24 date as a precise onset: minimum segment lengths can force breaks apart. Compare all three specifications. Changing the end month can also move breaks; these dates are re-estimated on the common window, rather than copied from either country report.</p>",
    "<details><summary>All dates and approximate conditional 95% intervals</summary>", html_table(breaks), "</details>",
    "<h2>Data quality, coverage and provenance</h2>", html_table(source_info),
    "<p>Each series passed checks for duplicate months, internal gaps, missing/non-finite values, valid 0-100 rates and a complete 108-month baseline. Unavailable months beyond the common end are excluded explicitly, not imputed. Input source hashes identify the precise exported observations used.</p>",
    "<h2>Methods and interpretation</h2><ul>", paste0("<li>", html_escape(notes), "</li>"), "</ul>",
    "<details><summary>Baseline diagnostics</summary>", html_table(collect("summary")), "</details>",
    "<details><summary>Residual diagnostics and all break-count selections</summary>", html_table(sensitivities), "</details>",
    "<details><summary>Descriptive deviation correlations</summary>", html_table(correlations), "</details>")
  if (nrow(staff)) parts <- c(parts, "<h2>Wales-only staff supplement</h2><p>Latest complete year for each available staff series, using its own 2010-2018 baseline. Check the year column: this can extend beyond the matched country window above. These overlapping groups are not independent replications, and there is no matched English staff-group comparison here.</p>",
    html_table(staff, digits = 3), "<details><summary>Excluded staff series</summary>", html_table(staff_errors), "</details>")
  for (country in names(audit_notes)) parts <- c(parts, paste0("<details><summary>", country, " input interpretation notes</summary><pre>", html_escape(audit_notes[[country]]), "</pre></details>"))
  parts <- c(parts, "<h2>Sources and reproducibility</h2><ul>",
    "<li><a href='https://digital.nhs.uk/data-and-information/publications/statistical/nhs-sickness-absence-rates'>NHS England sickness absence statistics</a>: source release is recorded in the copied England source_status.csv.</li>",
    "<li><a href='https://stats.gov.wales/en-GB/d4578da5-380f-439f-a056-8c5dd74d6509'>StatsWales monthly sickness absence</a> and <a href='https://www.gov.wales/nhs-workforce-staff-directly-employed-nhs-wales-quality-report-html'>Welsh workforce quality report</a>.</li>",
    "<li><a href='https://cran.r-project.org/web/packages/strucchange/strucchange.pdf'>strucchange documentation</a>: multiple-break estimation and conditional intervals.</li></ul>",
    "<p>Exports include the matched England-style panels (individual PNG/PDF and aligned pair), their plotted values and descriptive trend slopes, plus matched monthly and annual CSVs, baseline coefficients and diagnostics, breakpoint dates/BIC/sensitivities, Wales staff supplement, PNG/PDF figure, fitted_models.rds, configuration.R, sessionInfo.txt and archived input tables/audits. This HTML embeds its figure and can be shared alone.</p>",
    paste0("<small>Generated ", format(Sys.time(), tz = "UTC", usetz = TRUE), "; run ", html_escape(basename(run_dir)), ".</small></body></html>"))
  writeLines(parts, file.path(run_dir, "report.html"), useBytes = TRUE)
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(run_dir, "RUN_COMPLETE"))
  complete <- TRUE
  message("Comparison saved: ", normalizePath(run_dir), "\nOpen report.html in your browser.")
  invisible(list(run_dir = normalizePath(run_dir), annual = annual, breaks = breaks, correlations = correlations))
}

if (isTRUE(RUN_ANALYSIS)) comparison_result <- main()
