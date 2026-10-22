# England, Wales and Netherlands: comparative sickness absence analysis
# R >= 4.3; standalone script; no Pandoc and no live source downloads required.
# Install once: install.packages(c("strucchange", "sandwich", "jsonlite"))
# Run the whole file: source("compare_england_wales_netherlands.R", encoding="UTF-8")
#
# INPUTS: completed exports from nhs_sickness_rates.R, wales_sickness_rates.R
# and netherlands_sickness_rates.R. Point *_RUN to a run folder or its ZIP.
# NULL selects the latest compatible completed run under the corresponding root.
# Windows example: WALES_RUN <- "C:/Research/run_20261022_105258_64e84c7f59d4.zip"
# No existing country script is sourced/executed, and no input file is modified.
#
# READ BEFORE ANALYSIS
# England is not the entire UK. England and Wales are separate NHS series; they
# are never pooled, since denominator weights and the other UK nations are absent.
# The Dutch broad care sector is not an exact workforce match to either NHS.
# Main Dutch hospital comparator B000114 excludes university medical centres;
# those are shown separately as 422800. GP and whole-economy series give context.
# Sector totals/components overlap and are not independent replications.
# Rates are days lost / available employee days, not percentages of people sick.
# NHS monthly averages are descriptive quarterly proxies, NOT official quarterly
# ratios: the required FTE-available-day denominators are not in these exports.
# CBS may use workdays/calendar days; rounded survey rates and scope differ.
# Revisions, definitions, source status and availability are retained in audits.
# A comparison of one snapshot cannot establish stable definitions or revision history.
#
# METHODS (same approach as the earlier England-Netherlands comparison)
# NHS: separate monthly OLS baselines in 2010-2018, linear month index plus
# month-of-year effects. Mean observations and predictions over COMPLETE quarters.
# Netherlands: separate quarterly OLS baselines in 2010-2018, linear quarter index
# plus quarter-of-year effects. All series are truncated to identical quarters.
# No interpolation, annual-to-quarter allocation or partial-quarter extrapolation.
# NW covariance: lag 12 for monthly NHS fits, lag 4 for quarterly Dutch fits.
# Quarterly NHS mean CIs use the averaged design matrix, not averaged CI endpoints.
# Independently estimated Bai-Perron LEVEL breaks on unsmoothed quarterly deviations;
# BIC compares zero through all feasible counts. Primary minimum = 8 quarters;
# sensitivities = 4 and 6 quarters. No English dates are imposed on other series.
# Dates label FIRST quarter after a break; no month-specific onset is inferred.
# HAC intervals are conditional on the baseline and selected count; they do not
# propagate selection uncertainty or make the Gaussian BIC serially robust.
# Annual summaries: means of matched quarters, explicitly flagging partial years.
# Correlations and contrasts are descriptive, not causal or independent-sample tests.

# ---- Configuration ---------------------------------------------------------
ENGLAND_ROOT <- "england_sickness_results"
WALES_ROOT <- "wales_sickness_results"
NETHERLANDS_ROOT <- "netherlands_sickness_results"
ENGLAND_RUN <- getOption("three_country.england_run", NULL)
WALES_RUN <- getOption("three_country.wales_run", NULL)
NETHERLANDS_RUN <- getOption("three_country.netherlands_run", NULL)
OUT_DIR <- getOption("three_country.output_dir", "england_wales_netherlands_results")
START_YEAR <- 2010L
END_YEAR <- 2025L
BASELINE_END <- 2018L
HAC_LAG <- 4L
PRIMARY_MIN_QUARTERS <- 8L
SEGMENT_QUARTERS <- c(4L, 6L, 8L)
MAX_BREAKS <- NULL             # NULL = all feasible counts, including zero
INCLUDE_OTHER_DUTCH_SECTORS <- TRUE
RUN_ANALYSIS <- getOption("three_country.run_analysis", TRUE)
# Explicit codes avoid ambiguous matching by translated sector names.
# These five core series are required; edit this table if deliberately changing scope.
DUTCH_SERIES <- data.frame(dataset = "24015NED",
  sector_id = c("B000125", "B000114", "422800", "800040", "B000127"),
  label = c("Netherlands: health and social care", "Netherlands: hospitals/specialist care",
    "Netherlands: university medical centres", "Netherlands: GPs/health centres", "Netherlands: whole economy"),
  short_label = c("NL broad care", "NL hospitals", "NL university centres", "NL GPs", "NL whole economy"),
  plot_main = c(TRUE, TRUE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE)

write_csv <- function(x, path) utils::write.csv(x, path, row.names=FALSE, na="", fileEncoding="UTF-8")
text_file <- function(path) paste(readLines(path, encoding="UTF-8", warn=FALSE), collapse="\n")
observation_file <- function(country) if (country == "Netherlands")
  "quarterly_observations_2010_2025.csv" else "monthly_observations_2010_2025.csv"
is_complete_run <- function(path, country) {
  if (!dir.exists(path) || file.exists(file.path(path,"FAILED.txt")) ||
      !file.exists(file.path(path,observation_file(country)))) return(FALSE)
  if (file.exists(file.path(path,"RUN_COMPLETE"))) return(TRUE)
  files <- if (country == "Netherlands") c("source_status.csv","model_errors.csv","sessionInfo.txt",
    "baseline_and_deviations.csv","configuration.R") else c("report.html","source_status.csv","model_errors.csv")
  if (!all(file.exists(file.path(path,files)))) return(FALSE)
  errors <- tryCatch(read.csv(file.path(path,"model_errors.csv")), error=function(e) NULL)
  !is.null(errors) && nrow(errors) == 0L
}
resolve_run <- function(spec, root, country) {
  if (is.null(spec)) {
    if (!dir.exists(root)) stop(country, ": missing results root ", root, "; set its *_RUN folder/ZIP path")
    dirs <- list.dirs(root, recursive=FALSE, full.names=TRUE)
    dirs <- dirs[grepl("^run_",basename(dirs))]
    dirs <- dirs[vapply(dirs,is_complete_run,logical(1),country=country)]
    if (!length(dirs)) stop(country, ": no completed compatible run under ",root)
    spec <- dirs[order(basename(dirs),decreasing=TRUE)][1]
  }
  original <- normalizePath(spec,winslash="/",mustWork=TRUE); cleanup <- NULL
  if (!dir.exists(original)) {
    if (!grepl("\\.zip$",original,ignore.case=TRUE)) stop("Input must be a run folder or ZIP: ",original)
    entries <- utils::unzip(original,list=TRUE)$Name
    if (any(grepl("(^[/\\\\])|(^[A-Za-z]:)|((^|[/\\\\])\\.\\.([/\\\\]|$))",entries))) stop("Unsafe ZIP entry path")
    target <- entries[basename(entries)==observation_file(country)]
    if (length(target)!=1L) stop(country,": ZIP must contain exactly one observation export")
    cleanup <- tempfile("three_country_input_"); dir.create(cleanup)
    # Only data and audit text are needed, never execute input R files or load RDS.
    members <- entries[grepl("\\.(csv|txt|json|html|R)$",entries,ignore.case=TRUE)]
    utils::unzip(original,files=members,exdir=cleanup)
    path <- file.path(cleanup,dirname(target))
  } else path <- original
  if (!is_complete_run(path,country)) {
    if (!is.null(cleanup)) unlink(cleanup,recursive=TRUE)
    stop(country,": incomplete/failed country export: ",original)
  }
  list(path=path,original=original,cleanup=cleanup)
}
quarter_fields <- function(dates) {
  yr <- as.integer(format(dates,"%Y")); q <- (as.integer(format(dates,"%m"))-1L)%/%3L+1L
  data.frame(year=yr,quarter=q,quarter_label=paste0(yr,"Q",q),
    date=as.Date(sprintf("%d-%02d-01",yr,3L*q-2L)),quarter_index=(yr-START_YEAR)*4L+q-1L)
}
validate_grid <- function(d, frequency, label) {
  if (!nrow(d) || anyNA(d$date) || anyDuplicated(d$date)) stop(label,": empty/invalid/duplicate dates")
  if (any(!is.finite(d$rate_pct)) || any(d$rate_pct<0 | d$rate_pct>100)) stop(label,": missing or invalid rate")
  expected <- seq(as.Date(sprintf("%d-01-01",START_YEAR)),max(d$date),by=if(frequency==12L) "month" else "3 months")
  if (!identical(as.character(d$date),as.character(expected))) stop(label,": missing period(s) inside observation window")
  if (sum(as.integer(format(d$date,"%Y"))<=BASELINE_END)!=frequency*(BASELINE_END-START_YEAR+1L))
    stop(label,": incomplete 2010-2018 baseline")
  invisible(TRUE)
}
prepare_nhs <- function(raw,country) {
  req <- c("organisation_id","staff_id","date","rate_pct")
  if (!all(req %in% names(raw))) stop(country,": unexpected monthly export columns")
  org <- if(country=="England") "E" else "W"
  d <- raw[which(raw$organisation_id==org & raw$staff_id=="XXX"),,drop=FALSE]
  d$date <- as.Date(as.character(d$date),format="%Y-%m-%d")
  if (anyNA(d$date) || any(format(d$date,"%d")!="01")) stop(country,": invalid monthly dates")
  d <- d[d$date>=as.Date(sprintf("%d-01-01",START_YEAR)) & d$date<=as.Date(sprintf("%d-12-01",END_YEAR)),]
  d$rate_pct <- suppressWarnings(as.numeric(as.character(d$rate_pct)))
  d <- d[order(d$date),]; validate_grid(d,12L,country)
  d$year <- as.integer(format(d$date,"%Y")); d$month <- as.integer(format(d$date,"%m"))
  d$month_index <- (d$year-START_YEAR)*12L+d$month-1L
  d$month_factor <- factor(d$month,levels=1:12)
  d
}
prepare_dutch <- function(raw,dataset,id) {
  req <- c("dataset","sector_id","period_code","rate_pct")
  if (!all(req %in% names(raw))) stop("Unexpected Dutch observation columns")
  d <- raw[which(tolower(raw$dataset)==tolower(dataset) & raw$sector_id==id),,drop=FALSE]
  if (!nrow(d)) stop("Dutch series absent: ",dataset,"::",id)
  if (anyNA(d$period_code) || any(!grepl("^[0-9]{4}KW0[1-4]$",d$period_code))) stop("Non-quarterly or invalid Dutch period code")
  d$year <- as.integer(substr(d$period_code,1,4)); d$quarter <- as.integer(substr(d$period_code,8,8))
  d <- d[d$year>=START_YEAR & d$year<=END_YEAR,]
  d$date <- as.Date(sprintf("%d-%02d-01",d$year,3L*d$quarter-2L))
  d$rate_pct <- suppressWarnings(as.numeric(as.character(d$rate_pct)))
  d <- d[order(d$date),]
  if ("row_present" %in% names(d)) {
    present <- tolower(as.character(d$row_present)) %in% c("true","1")
    if (!any(present)) stop("Dutch series has no published rows")
    # Permit only an explicitly unpublished tail, not a missing published value.
    last <- max(which(present))
    if (any(!present[seq_len(last)])) stop("Dutch series has an unpublished internal quarter")
    d <- d[seq_len(last),]
  }
  validate_grid(d,4L,paste(dataset,id))
  d$quarter_label <- paste0(d$year,"Q",d$quarter)
  d$quarter_index <- (d$year-START_YEAR)*4L+d$quarter-1L
  d$quarter_factor <- factor(d$quarter,levels=1:4)
  d
}
fit_native_baseline <- function(d,monthly=FALSE) {
  train <- d[d$year<=BASELINE_END,]
  form <- if(monthly) rate_pct~month_index+month_factor else rate_pct~quarter_index+quarter_factor
  fit <- stats::lm(form,data=train)
  V <- sandwich::NeweyWest(fit,lag=if(monthly)12L else HAC_LAG,prewhite=FALSE,adjust=TRUE)
  X <- stats::model.matrix(stats::delete.response(stats::terms(fit)),data=d)
  d$baseline_pct <- as.numeric(predict(fit,newdata=d))
  se <- sqrt(pmax(0,rowSums((X%*%V)*X)))
  d$baseline_mean_low95 <- d$baseline_pct-qnorm(.975)*se
  d$baseline_mean_high95 <- d$baseline_pct+qnorm(.975)*se
  d$deviation_pp <- d$rate_pct-d$baseline_pct
  d$deviation_relative_pct <- ifelse(d$baseline_pct>0,100*d$deviation_pp/d$baseline_pct,NA_real_)
  co <- coef(fit); cs <- sqrt(diag(V)); term <- if(monthly) "month_index" else "quarter_index"
  coefficients <- data.frame(term=names(co),estimate=unname(co),se_HAC=unname(cs),
    low95_HAC=unname(co-qnorm(.975)*cs),high95_HAC=unname(co+qnorm(.975)*cs))
  diagnostics <- data.frame(frequency=if(monthly)"monthly" else "quarterly", n_baseline=nrow(train),
    baseline_slope_pp_year=unname(co[term])*if(monthly)12 else 4,
    LjungBox_p=stats::Box.test(residuals(fit),lag=if(monthly)12 else 4,type="Ljung-Box")$p.value,
    baseline_outside_0_100=any(d$baseline_pct<0|d$baseline_pct>100))
  list(data=d,fit=fit,vcov=V,X=X,coefficients=coefficients,diagnostics=diagnostics)
}
aggregate_nhs <- function(fit,country) {
  d <- fit$data; fields <- quarter_fields(d$date)
  groups <- split(seq_len(nrow(d)),fields$quarter_label)
  audit <- do.call(rbind,lapply(groups,function(ix) data.frame(country=country,quarter=fields$quarter_label[ix[1]],
    n_months=length(ix),months=paste(format(d$date[ix],"%Y-%m"),collapse="; "),complete_quarter=length(ix)==3L)))
  groups <- groups[lengths(groups)==3L]
  result <- do.call(rbind,lapply(groups,function(ix) {
    Xmean <- colMeans(fit$X[ix,,drop=FALSE]); se <- sqrt(max(0,drop(t(Xmean)%*%fit$vcov%*%Xmean)))
    b <- mean(d$baseline_pct[ix]); actual <- mean(d$rate_pct[ix])
    data.frame(fields[ix[1],],rate_pct=actual,baseline_pct=b,
      baseline_mean_low95=b-qnorm(.975)*se,baseline_mean_high95=b+qnorm(.975)*se,
      deviation_pp=actual-b,deviation_relative_pct=if(b>0)100*(actual-b)/b else NA_real_)
  }))
  list(data=result,audit=audit)
}
hac_vcov <- function(model) {
  lag <- min(HAC_LAG,NROW(sandwich::estfun(model))-1L)
  # Same Bartlett weights as NeweyWest, without its trailing zero weight.
  # This avoids a harmless "more weights than observations" warning for
  # four-quarter segments while retaining the same covariance calculation.
  sandwich::vcovHAC(model,weights=1-(0:lag)/(lag+1),prewhite=FALSE,adjust=TRUE)
}
fit_breaks <- function(d, h) {
  max_k <- ceiling(nrow(d) / h) - 2L
  if (!is.null(MAX_BREAKS)) max_k <- min(MAX_BREAKS, max_k)
  if (max_k < 1L) stop("Minimum segment length leaves no admissible break model")
  model_warnings <- character()
  bp <- withCallingHandlers(
    strucchange::breakpoints(deviation_pp ~ 1, data = d, h = h, breaks = max_k),
    warning = function(w) {
      model_warnings <<- c(model_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  candidates <- 0:max_k
  bic <- as.numeric(stats::AIC(bp, breaks = candidates, k = log(nrow(d))))
  if (any(!is.finite(bic))) stop("Non-finite BIC: degenerate series requires review")
  k <- candidates[which.min(bic)]
  chosen <- if (k) strucchange::breakpoints(bp, breaks = k)$breakpoints else integer()
  # strucchange indexes the LAST observation in the preceding segment.
  segment <- findInterval(seq_len(nrow(d)) - 1L, chosen) + 1L
  means <- tapply(d$deviation_pp, segment, mean)
  fitted <- unname(means[as.character(segment)])
  breaks <- data.frame(min_segment_quarters = integer(), break_number = integer(),
    last_quarter_before = character(), first_quarter_after = character(), shift_pp = numeric(),
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
      out[ok] <- d$quarter_label[as.integer(j[ok])]; out
    }
    breaks <- data.frame(min_segment_quarters = h, break_number = seq_len(k),
      last_quarter_before = d$quarter_label[chosen], first_quarter_after = d$quarter_label[chosen + 1L],
      shift_pp = diff(as.numeric(means)), ci95_first_after_low = labels(bounds[, 1] + 1L),
      ci95_first_after_high = labels(bounds[, 2] + 1L), ci_note = ci_note)
    bad <- is.na(breaks$ci95_first_after_low) | is.na(breaks$ci95_first_after_high)
    breaks$ci_note[bad] <- paste(breaks$ci_note[bad], "Bounds unavailable/outside observed window.")
    if (length(warn)) breaks$ci_note <- paste(breaks$ci_note, paste(unique(warn), collapse = "; "))
  }
  endpoints <- c(0L, chosen, nrow(d))
  segments <- do.call(rbind, lapply(seq_len(length(endpoints) - 1L), function(i) {
    ix <- seq.int(endpoints[i] + 1L, endpoints[i + 1L])
    data.frame(min_segment_quarters = h, segment = i, start_quarter = d$quarter_label[min(ix)],
      end_quarter = d$quarter_label[max(ix)], n_quarters = length(ix), mean_deviation_pp = mean(d$deviation_pp[ix]))
  }))
  summary <- data.frame(min_segment_quarters = h, selected_breaks = k,
    diagnostic_warnings = paste(unique(c(model_warnings, warn)), collapse = "; "),
    at_candidate_limit = k == max_k, BIC_gain_over_zero = bic[1] - min(bic),
    segmented_LjungBox_lag4_p = stats::Box.test(d$deviation_pp - fitted, lag = 4, type = "Ljung-Box")$p.value,
    first_quarters_after = paste(d$quarter_label[chosen + 1L], collapse = "; "))
  list(breaks = breaks, segments = segments, summary = summary,
    bic = data.frame(min_segment_quarters = h, n_breaks = candidates, BIC = bic,
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


annual_summary <- function(d) do.call(rbind,lapply(split(d,d$year),function(z) {
  a <- mean(z$rate_pct); b <- mean(z$baseline_pct)
  data.frame(year=z$year[1],n_quarters=nrow(z),period=paste0("Q",min(z$quarter),"-Q",max(z$quarter)),
    complete_year=nrow(z)==4L,mean_rate_pct=a,mean_baseline_pct=b,mean_excess_pp=a-b,
    relative_excess_pct=if(b>0)100*(a-b)/b else NA_real_)
}))
add_identity <- function(x,key,registry) {
  id <- registry[match(key,registry$key),c("key","country","label","dataset","sector_id","role")]
  cbind(id[rep(1L,nrow(x)),,drop=FALSE],x)
}
embed_image <- function(path,alt) paste0("<img alt='",html_escape(alt),"' src='data:image/png;base64,",
  jsonlite::base64_enc(readBin(path,"raw",n=file.info(path)$size)),"'>")
plot_main <- function(results,registry) {
  keys <- registry$key[registry$plot_main]
  cols <- setNames(c("#0072B2","#8E44AD","#D55E00","#17836B",grDevices::hcl.colors(max(1,length(keys)-4),"Dark 3"))[seq_along(keys)],keys)
  graphics::par(mfrow=c(2,1),mar=c(4.3,4.8,3.2,1),mgp=c(2.8,.8,0))
  for (zoom in c(FALSE,TRUE)) {
    subset_data <- lapply(results[keys],function(d) if(zoom)d[d$year>=2018 & d$year<=2020,] else d)
    first <- subset_data[[1]]; x <- first$year+(first$quarter-1)/4
    ylim <- range(c(0,unlist(lapply(subset_data,`[[`,"deviation_pp"))))
    graphics::plot(range(x),ylim,type="n",xlab="",ylab="Above own baseline (percentage points)",
      xaxt=if(zoom)"n" else "s",main=if(zoom)"2018-2020 detail: complete quarterly observations" else "Sickness absence relative to each 2010-2018 baseline")
    graphics::abline(h=0,lty=2,col="grey50")
    graphics::abline(v=2019,lty=3,col="grey75")
    for (key in keys) {
      d <- subset_data[[key]]
      graphics::lines(x,d$deviation_pp,col=cols[key],lwd=2,lty=1,
        type=if(zoom)"o" else "l",pch=if(key=="Wales")17 else 16,cex=.6)
    }
    if (zoom) graphics::axis(1,at=x,labels=first$quarter_label,las=2,cex.axis=.75)
    else graphics::legend("topleft",registry$short_label[match(keys,registry$key)],col=cols,lty=1,lwd=2,bty="n",cex=.85)
  }
}
plot_breaks <- function(breaks,registry,results) {
  keys <- registry$key[registry$role %in% c("NHS","Dutch core")]
  graphics::par(mar=c(4,10.5,3,1),mgp=c(2.7,.8,0))
  first <- results[[keys[1]]]; x <- first$year+(first$quarter-1)/4
  graphics::plot(range(x),c(.5,length(keys)+.5),type="n",yaxt="n",xlab="First quarter after estimated break",ylab="",
    main=paste("Independent level breaks: minimum",PRIMARY_MIN_QUARTERS,"quarters"))
  graphics::axis(2,seq_along(keys),registry$short_label[match(keys,registry$key)],las=1,cex.axis=.8)
  for(i in seq_along(keys)) {
    graphics::abline(h=i,col="grey90")
    b <- breaks[breaks$key==keys[i] & breaks$min_segment_quarters==PRIMARY_MIN_QUARTERS,]
    if(nrow(b)) graphics::points(x[match(b$first_quarter_after,first$quarter_label)],rep(i,nrow(b)),pch=19,col="#0072B2")
  }
}
main <- function() {
  pkgs <- c("strucchange","sandwich","jsonlite")
  absent <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop("Install packages: ",paste(absent,collapse=", "))
  stopifnot(START_YEAR==2010L,BASELINE_END==2018L,END_YEAR==2025L,
    all(SEGMENT_QUARTERS>=4L),all(SEGMENT_QUARTERS%%1==0),PRIMARY_MIN_QUARTERS %in% SEGMENT_QUARTERS)
  if(!is.null(MAX_BREAKS) && (length(MAX_BREAKS)!=1L || !is.finite(MAX_BREAKS) || MAX_BREAKS<1 || MAX_BREAKS%%1!=0))
    stop("MAX_BREAKS must be NULL or a positive integer")
  if(!nrow(DUTCH_SERIES) || anyDuplicated(paste(DUTCH_SERIES$dataset,DUTCH_SERIES$sector_id))) stop("Invalid Dutch core series configuration")
  inputs <- list()
  on.exit(for(x in inputs) if(!is.null(x$cleanup)) unlink(x$cleanup,recursive=TRUE),add=TRUE)
  inputs$England <- resolve_run(ENGLAND_RUN,ENGLAND_ROOT,"England")
  inputs$Wales <- resolve_run(WALES_RUN,WALES_ROOT,"Wales")
  inputs$Netherlands <- resolve_run(NETHERLANDS_RUN,NETHERLANDS_ROOT,"Netherlands")
  dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  run_dir <- tempfile(paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S",tz="UTC"),"_"),tmpdir=OUT_DIR)
  dir.create(run_dir);dir.create(file.path(run_dir,"source"));dir.create(file.path(run_dir,"plots"))
  completed <- FALSE
  on.exit(if(!completed)writeLines("Run incomplete; do not use partial outputs.",file.path(run_dir,"FAILED.txt")),add=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(run_dir,"sessionInfo.txt"))
  dput(list(inputs=lapply(inputs,`[[`,"original"),START_YEAR=START_YEAR,END_YEAR=END_YEAR,BASELINE_END=BASELINE_END,
    Dutch_core=DUTCH_SERIES,include_other_Dutch=INCLUDE_OTHER_DUTCH_SECTORS,quarterly_h=SEGMENT_QUARTERS,
    primary_h=PRIMARY_MIN_QUARTERS,maximum_breaks=MAX_BREAKS,NHS_HAC_lag=12L,Dutch_HAC_lag=HAC_LAG,
    quarterly_method="equal mean of three monthly NHS rates and predictions; retain published Dutch quarterly rates"),
    file.path(run_dir,"configuration.R"))
  raw <- list(); source_info <- list(); source_notes <- list(); input_status <- list()
  for(country in names(inputs)) {
    src <- inputs[[country]]$path;dest <- file.path(run_dir,"source",country);dir.create(dest)
    # Preserve data, audits and metadata; model objects and binary originals remain in the country run.
    files <- list.files(src,recursive=TRUE,full.names=FALSE)
    archive <- files[grepl("(^|/)(monthly_observations_2010_2025.csv|quarterly_observations_2010_2025.csv|source_status.csv|quality_checks.csv|known_definition_changes.csv|READ_BEFORE_ANALYSIS.txt|configuration.R|metadata.json)$|revision_changes.csv$",files)]
    for(f in archive) {
      target <- file.path(dest,f);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      if(!file.copy(file.path(src,f),target))stop("Could not archive input: ",f)
    }
    path <- file.path(src,observation_file(country))
    raw[[country]] <- read.csv(path,stringsAsFactors=FALSE,check.names=FALSE,colClasses="character",na.strings="")
    source_info[[country]] <- data.frame(country=country,input_run=inputs[[country]]$original,
      observation_md5=unname(tools::md5sum(path)),input_rows=nrow(raw[[country]]))
    notes <- file.path(src,"READ_BEFORE_ANALYSIS.txt")
    source_notes[[country]] <- if(file.exists(notes))text_file(notes) else "No source notes exported"
    status <- file.path(src,"source_status.csv")
    input_status[[country]] <- if(file.exists(status))read.csv(status,stringsAsFactors=FALSE) else data.frame(note="No source status")
  }
  source_info <- do.call(rbind,source_info)
  write_csv(source_info,file.path(run_dir,"source_provenance.csv"))
  series <- list(England=prepare_nhs(raw$England,"England"),Wales=prepare_nhs(raw$Wales,"Wales"))
  registry <- data.frame(key=c("England","Wales"),country=c("England","Wales"),
    label=c("England: NHS all staff","Wales: NHS all staff"),short_label=c("England NHS","Wales NHS"),
    dataset="NHS",sector_id="XXX",role="NHS",frequency=12L,plot_main=TRUE,original_sector="All staff")
  exclusions <- data.frame(key=character(),reason=character())
  add_dutch <- function(dataset,id,label,short_label,role,plot_main) {
    key <- paste(dataset,id,sep="::");d <- prepare_dutch(raw$Netherlands,dataset,id)
    series[[key]] <<- d
    registry <<- rbind(registry,data.frame(key=key,country="Netherlands",label=label,short_label=short_label,
      dataset=dataset,sector_id=id,role=role,frequency=4L,plot_main=plot_main,
      original_sector=if("sector" %in% names(d))d$sector[1] else id))
  }
  for(i in seq_len(nrow(DUTCH_SERIES))) {
    s <- DUTCH_SERIES[i,]
    # Missing core series abort explicitly; a differently scoped sector is never substituted.
    add_dutch(s$dataset,s$sector_id,s$label,s$short_label,"Dutch core",s$plot_main)
  }
  if(INCLUDE_OTHER_DUTCH_SECTORS) {
    candidates <- unique(raw$Netherlands[,c("dataset","sector_id","sector")])
    for(i in seq_len(nrow(candidates))) {
      s <- candidates[i,];key <- paste(s$dataset,s$sector_id,sep="::")
      if(key %in% registry$key)next
      rows <- raw$Netherlands[which(raw$Netherlands$dataset==s$dataset & raw$Netherlands$sector_id==s$sector_id),]
      if("company_size" %in% names(rows) && any(tolower(rows$company_size)=="true",na.rm=TRUE))next
      if("selected" %in% names(rows) && !any(tolower(rows$selected)=="true",na.rm=TRUE))next
      # These two broad aggregates overlap across CBS/AZW. Drop duplicates only
      # when their complete quarter/value arrays are actually identical.
      alias <- switch(key,"80072ned::422400"="24015NED::B000125","80072ned::T001081"="24015NED::B000127",NULL)
      if(!is.null(alias) && alias %in% names(series)) {
        candidate <- tryCatch(prepare_dutch(raw$Netherlands,s$dataset,s$sector_id),error=identity)
        if(!inherits(candidate,"error") && identical(candidate$quarter_label,series[[alias]]$quarter_label) &&
           isTRUE(all.equal(candidate$rate_pct,series[[alias]]$rate_pct,tolerance=1e-12))) {
          exclusions <- rbind(exclusions,data.frame(key=key,reason=paste("Verified duplicate values/quarters of",alias)));next
        }
      }
      tryCatch(add_dutch(s$dataset,s$sector_id,paste0("Netherlands: ",s$sector),as.character(s$sector),"Dutch context",FALSE),
        error=function(e)exclusions <<- rbind(exclusions,data.frame(key=key,reason=conditionMessage(e))))
    }
  }
  # Determine complete-quarter coverage before fitting any model.
  availability <- do.call(rbind,lapply(registry$key,function(key) {
    d <- series[[key]];monthly <- registry$frequency[match(key,registry$key)]==12L
    if(monthly) {
      f <- quarter_fields(d$date);counts <- table(f$quarter_label)
      complete <- unique(f$date[f$quarter_label %in% names(counts)[counts==3L]])
    } else complete <- d$date
    data.frame(key=key,source_first=as.character(min(d$date)),source_last=as.character(max(d$date)),
      n_native_observations=nrow(d),last_complete_quarter_date=as.character(max(complete)),
      n_complete_quarters=length(complete))
  }))
  core_keys <- registry$key[registry$role %in% c("NHS","Dutch core")]
  common_end <- min(as.Date(availability$last_complete_quarter_date[availability$key %in% core_keys]))
  if(common_end<=as.Date(sprintf("%d-10-01",BASELINE_END)))stop("No post-baseline complete common quarter")
  expected <- seq(as.Date(sprintf("%d-01-01",START_YEAR)),common_end,by="3 months")
  # Optional context series do not shorten the primary common window.
  incomplete <- availability$key[as.Date(availability$last_complete_quarter_date)<common_end]
  if(length(incomplete)) {
    exclusions <- rbind(exclusions,data.frame(key=incomplete,reason="Context series does not reach common end; omitted as whole series"))
    registry <- registry[!registry$key %in% incomplete,];series[incomplete] <- NULL
  }
  availability$included <- availability$key %in% registry$key
  availability$common_end <- paste0(format(common_end,"%Y"),"Q",(as.integer(format(common_end,"%m"))-1L)%/%3L+1L)
  availability$matched_quarters <- ifelse(availability$included,length(expected),0L)
  availability$excluded_complete_tail_quarters <- ifelse(availability$included,availability$n_complete_quarters-length(expected),NA_integer_)
  write_csv(availability,file.path(run_dir,"quality_checks.csv"))
  write_csv(registry,file.path(run_dir,"series_registry.csv"))
  write_csv(exclusions,file.path(run_dir,"excluded_series.csv"))
  notes <- c(
    "England is not UK-wide. England and Wales are analysed separately and are never pooled without denominator weights.",
    "NHS scope differs from Dutch health/social care and from narrower branches. B000114 excludes university medical centres; 422800 is displayed separately. Whole-economy results are context, not a matched NHS comparator.",
    "The prior comparison's method is retained: monthly NHS baselines fitted in 2010-2018, then equal three-month averages; Dutch quarterly baselines fitted in 2010-2018. Baselines are re-estimated from raw rates rather than copied from country model exports.",
    "Both NHS baselines use a linear month index. This differs slightly from the original England script's calendar-day index. The common frequency and breakpoint search also differ from the original monthly England/Wales analyses.",
    "A complete quarter requires all three NHS months. Partial trailing quarters are dropped. The same complete quarters are used for every retained model; no interpolation or partial-quarter extrapolation occurs.",
    "Quarterly NHS values are unweighted averages of monthly rates, not official denominator-weighted quarterly rates. Available FTE-day denominators are not in these input exports. Do not interpret rate differences as counts of sick staff.",
    "CBS rates are rounded survey estimates; employers can use workdays or calendar days. Different recording rules, employer coverage and staff composition remain potential sources of between-country differences.",
    "Breaks are independent level shifts in quarterly observed-minus-baseline deviations over the full common window. BIC includes zero through all feasible counts unless MAX_BREAKS is explicitly capped; a minimum duration of 4/6/8 quarters is examined.",
    "Break dates identify the first quarter of the new segment, not the precise month of onset. Minimum segment duration, common endpoint and baseline assumptions can move the selected dates.",
    "Baseline mean CIs use Newey-West lag 12 for monthly NHS fits and lag 4 for Dutch quarterly fits. NHS quarterly CIs use the mean design vector. These are mean intervals, not prediction intervals.",
    "Approximate HAC break-date intervals are conditional on the baseline and selected break count. They do not include baseline/model-selection uncertainty; HAC does not make BIC robust to serial dependence.",
    "Annual summaries are unweighted means of matched quarters. A three-quarter 2025 result means January-September, not a full-year rate. Relative excess is 100 times mean excess divided by mean baseline.",
    "Overlapping Dutch branches and aggregates are not independent replications. Two known CBS/AZW aggregate duplicates are removed only if their dates and rates match exactly; exclusions are logged.",
    "Quarterly correlations and differences of baseline-adjusted rates are descriptive. No independent-observation p-values, causal effects, infection dates or pathogen attribution are inferred.",
    "Input source-status, definition notes, metadata and revision files are archived. This script performs no new source download and cannot reconstruct unseen revision vintages or establish unchanged definitions.")
  writeLines(notes,file.path(run_dir,"READ_BEFORE_ANALYSIS.txt"))
  results <- baselines <- list(); aggregation_audits <- list(); annual <- list()
  quarterly_cols <- c("year","quarter","quarter_label","date","quarter_index","rate_pct","baseline_pct",
    "baseline_mean_low95","baseline_mean_high95","deviation_pp","deviation_relative_pct")
  for(key in registry$key) {
    monthly <- registry$frequency[match(key,registry$key)]==12L
    base <- fit_native_baseline(series[[key]],monthly)
    baselines[[key]] <- base
    if(monthly) {
      agg <- aggregate_nhs(base,key);d <- agg$data;aggregation_audits[[key]] <- agg$audit
      write_csv(base$data,file.path(run_dir,paste0(tolower(key),"_monthly_baseline.csv")))
    } else d <- base$data
    d <- d[d$date<=common_end,quarterly_cols]
    stopifnot(identical(as.character(d$date),as.character(expected)))
    results[[key]] <- d;annual[[key]] <- add_identity(annual_summary(d),key,registry)
  }
  aggregation_audit <- do.call(rbind,aggregation_audits)
  aggregation_audit$used_in_common_window <- aggregation_audit$complete_quarter &
    aggregation_audit$quarter <= availability$common_end[1]
  write_csv(aggregation_audit,file.path(run_dir,"quarter_aggregation_audit.csv"))
  all_data <- do.call(rbind,lapply(names(results),function(key)add_identity(results[[key]],key,registry)))
  annual <- do.call(rbind,annual)
  write_csv(all_data,file.path(run_dir,"matched_quarterly_baseline_and_deviations.csv"))
  write_csv(annual,file.path(run_dir,"annual_mean_quarterly_rates.csv"))
  coefficient_table <- do.call(rbind,lapply(names(baselines),function(key)add_identity(baselines[[key]]$coefficients,key,registry)))
  diagnostic_table <- do.call(rbind,lapply(names(baselines),function(key)add_identity(baselines[[key]]$diagnostics,key,registry)))
  write_csv(coefficient_table,file.path(run_dir,"baseline_coefficients.csv"))
  write_csv(diagnostic_table,file.path(run_dir,"baseline_diagnostics.csv"))
  bfits <- bdates <- bsegments <- bscores <- bsummaries <- list()
  for(key in names(results)) {
    message("Independent quarterly breaks: ",key)
    for(h in unique(SEGMENT_QUARTERS)) {
      r <- fit_breaks(results[[key]],h);id <- paste(key,h,sep="__")
      bfits[[id]] <- r
      bdates[[id]] <- add_identity(r$breaks,key,registry);bsegments[[id]] <- add_identity(r$segments,key,registry)
      bscores[[id]] <- add_identity(r$bic,key,registry);bsummaries[[id]] <- add_identity(r$summary,key,registry)
    }
  }
  breaks <- do.call(rbind,bdates);segments <- do.call(rbind,bsegments)
  scores <- do.call(rbind,bscores);sensitivity <- do.call(rbind,bsummaries)
  sensitivity$post2018_first_quarters <- vapply(seq_len(nrow(sensitivity)),function(i) {
    b <- breaks[breaks$key==sensitivity$key[i] & breaks$min_segment_quarters==sensitivity$min_segment_quarters[i],]
    paste(b$first_quarter_after[b$first_quarter_after>paste0(BASELINE_END,"Q4")],collapse="; ")
  },character(1))
  write_csv(breaks,file.path(run_dir,"breakpoints_all_sensitivities.csv"))
  write_csv(breaks[breaks$min_segment_quarters==PRIMARY_MIN_QUARTERS,],file.path(run_dir,"breakpoints_primary.csv"))
  write_csv(sensitivity,file.path(run_dir,"breakpoint_sensitivity.csv"))
  write_csv(segments,file.path(run_dir,"segments_all_sensitivities.csv"))
  write_csv(scores,file.path(run_dir,"breakpoint_BIC.csv"))
  comparisons <- c(list(c("England","Wales")),unlist(lapply(registry$key[registry$country=="Netherlands"],
    function(key)list(c("England",key),c("Wales",key))),recursive=FALSE))
  correlations <- contrasts <- list()
  for(i in seq_along(comparisons)) {
    pair <- comparisons[[i]];a <- results[[pair[1]]];b <- results[[pair[2]]]
    for(window in c("All quarters","Post-2018")) {
      ix <- if(window=="All quarters")rep(TRUE,nrow(a)) else a$year>BASELINE_END
      correlations[[length(correlations)+1L]] <- data.frame(series_A=pair[1],series_B=pair[2],window=window,n_quarters=sum(ix),
        r_deviations=if(sd(a$deviation_pp[ix])>0 && sd(b$deviation_pp[ix])>0)cor(a$deviation_pp[ix],b$deviation_pp[ix]) else NA_real_)
    }
    for(year in sort(unique(a$year[a$year>BASELINE_END]))) {
      ix <- a$year==year
      contrasts[[length(contrasts)+1L]] <- data.frame(series_A=pair[1],series_B=pair[2],year=year,n_quarters=sum(ix),
        complete_year=sum(ix)==4L,rate_A_minus_B_pp=mean(a$rate_pct[ix]-b$rate_pct[ix]),
        excess_A_minus_B_pp=mean(a$deviation_pp[ix]-b$deviation_pp[ix]))
    }
  }
  correlations <- do.call(rbind,correlations);contrasts <- do.call(rbind,contrasts)
  write_csv(correlations,file.path(run_dir,"descriptive_correlations.csv"))
  write_csv(contrasts,file.path(run_dir,"baseline_adjusted_contrasts.csv"))
  detail <- all_data[all_data$key %in% core_keys & all_data$year>=2018 & all_data$year<=2020,]
  write_csv(detail,file.path(run_dir,"quarterly_detail_2018_2020.csv"))
  saveRDS(list(native_baselines=baselines,quarterly_break_models=bfits,source_provenance=source_info),file.path(run_dir,"fitted_models.rds"))
  png_main <- file.path(run_dir,"plots","comparison.png");png_breaks <- file.path(run_dir,"plots","breakpoints.png")
  grDevices::png(png_main,width=2000,height=1500,res=160)
  tryCatch(plot_main(results,registry),finally=grDevices::dev.off())
  grDevices::png(png_breaks,width=1900,height=850,res=160)
  tryCatch(plot_breaks(breaks,registry,results),finally=grDevices::dev.off())
  grDevices::pdf(file.path(run_dir,"comparison_plots.pdf"),width=12.5,height=9.5)
  tryCatch({plot_main(results,registry);plot_breaks(breaks,registry,results)},finally=grDevices::dev.off())
  last_year <- max(all_data$year)
  latest <- annual[annual$key %in% core_keys & annual$year==last_year,]
  history <- annual[annual$key %in% core_keys & annual$year>BASELINE_END,]
  findings <- c(paste0("The matched window is 2010Q1 to ",availability$common_end[1]," (",length(expected)," complete quarters). England and Wales remain separate national NHS series."),
    vapply(seq_len(nrow(latest)),function(i) {
      z <- latest[i,]
      sprintf("%s: %d %s mean %.2f%%; excess %+.2f percentage points (%+.1f%% relative to its own baseline).",z$label,z$year,z$period,z$mean_rate_pct,z$mean_excess_pp,z$relative_excess_pct)
    },character(1)))
  writeLines(findings,file.path(run_dir,"findings.txt"))
  core_registry <- registry[registry$key %in% core_keys,c("label","dataset","sector_id","original_sector")]
  display_annual <- function(x) x[,c("label","year","period","n_quarters","complete_year","mean_rate_pct","mean_baseline_pct","mean_excess_pp","relative_excess_pct")]
  parts <- c("<!doctype html><html lang='en'><head><meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>England, Wales and Netherlands: sickness absence</title>",
    "<style>body{font:16px/1.6 system-ui,sans-serif;color:#203040;max-width:1200px;margin:40px auto;padding:0 24px}h1{line-height:1.2;font-size:34px}h2{margin-top:38px;color:#174d70}.note{background:#eef4f8;border-left:4px solid #0072B2;padding:16px 24px}table{border-collapse:collapse;font-size:14px}td,th{border-bottom:1px solid #d9e1e8;padding:8px;text-align:left;vertical-align:top}th{background:#eaf1f7}.scroll{overflow:auto}img{width:100%;height:auto}li{margin:.5em 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}summary{cursor:pointer;font-weight:600;margin-top:20px}a{color:#0067A5}@media print{body{margin:0;max-width:none}h2{break-after:avoid}tr,img{break-inside:avoid}}</style></head><body>",
    "<h1>England, Wales and the Netherlands</h1><p>Comparative sickness absence, aligned to complete quarters. Separate 2010-2018 trend-and-seasonal baselines.</p>",
    "<div class='note'><strong>Scope</strong><p>England is the geographical coverage of the supplied NHS England series; it is not a UK-wide rate. Wales is shown separately. There is no pooled UK or England-Wales estimate.</p><ul>",
    paste0("<li>",html_escape(findings),"</li>"),"</ul></div>",
    "<h2>Comparison with the earlier England-Netherlands analysis</h2><p>The same monthly-baseline-to-quarterly-average approach is used for both NHS series. Dutch baselines use quarter-of-year effects. Adding Wales therefore gives a separate NHS comparison; it does not change the Dutch sector definitions.</p>",
    embed_image(png_main,"Quarterly baseline deviations and 2018-2020 detail for England, Wales and Dutch care sectors"),
    "<p>Quarterly NHS points are means of three monthly rates and predictions. The broad Dutch care aggregate includes social care; the hospital/specialist series excludes university medical centres. Scope and recording differences remain.</p>",
    "<h2>Latest matched year</h2><p>Check the period and complete_year columns. Q1-Q3 means January-September. These are averages of matched quarterly values, not official weighted annual rates.</p>",
    html_table(display_annual(latest),digits=3),
    "<details><summary>Core series: annual results from 2019</summary>",html_table(display_annual(history),digits=3),"</details>",
    "<details><summary>Core series: 2018-2020 quarterly detail</summary>",
    html_table(detail[,c("label","quarter_label","rate_pct","baseline_pct","deviation_pp")],digits=3),"</details>",
    "<h2>Independently estimated breakpoints</h2>",embed_image(png_breaks,"Independent primary break dates by NHS and Dutch core sector"),
    "<p>Dates are the first quarter after each estimated level break. The full common window, including baseline years, is searched. No English or pandemic dates are imposed. Counts include any pre-2019 breaks; the post-2018 column below isolates later dates.</p>",
    html_table(sensitivity[sensitivity$key %in% core_keys,c("label","min_segment_quarters","selected_breaks","at_candidate_limit","post2018_first_quarters")]),
    "<p>Changing the minimum segment duration or endpoint can move dates. An early estimated break is not evidence for a precise month of onset or a particular cause.</p>",
    "<details><summary>Core dates, level shifts and approximate conditional intervals</summary>",
    html_table(breaks[breaks$key %in% core_keys,c("label","min_segment_quarters","last_quarter_before","first_quarter_after","shift_pp","ci95_first_after_low","ci95_first_after_high","ci_note")]),"</details>",
    "<h2>Dutch sector context</h2><p>Each retained sector is fitted independently. Aggregate and component sectors overlap; this is descriptive sector coverage, not a count of independent replications.</p>",
    "<details><summary>Latest matched year: all retained Dutch sectors</summary>",
    html_table(display_annual(annual[annual$country=="Netherlands" & annual$year==last_year,]),digits=3),"</details>",
    "<details><summary>All retained sectors: primary break specification</summary>",
    html_table(sensitivity[sensitivity$min_segment_quarters==PRIMARY_MIN_QUARTERS,c("label","selected_breaks","first_quarters_after","BIC_gain_over_zero","segmented_LjungBox_lag4_p")]),"</details>",
    "<h2>Descriptive concordance</h2><p>Correlations use baseline-adjusted quarterly deviations over identical dates. Serial dependence precludes treating these quarters as independent replicates; no causal or significance claim is made.</p>",
    html_table(correlations[correlations$series_B %in% core_keys & correlations$window=="Post-2018",],digits=3),
    "<h2>Data checks, exclusions and provenance</h2>",html_table(core_registry),
    html_table(availability[availability$key %in% core_keys,]),
    "<p>Every retained series has a complete baseline, unique dates, no internal gaps and finite rates within 0-100. All retained models use the same quarters. Input checks precede fitting; optional context exclusions are logged.</p>",
    "<details><summary>Partial quarters and quarters excluded at the common endpoint</summary>",html_table(aggregation_audit[!aggregation_audit$used_in_common_window,]),"</details>",
    "<details><summary>Excluded or duplicate Dutch series</summary>",html_table(exclusions),"</details>",
    "<details><summary>Input locations and hashes</summary>",html_table(source_info),"</details>")
  for(country in names(input_status))parts <- c(parts,paste0("<details><summary>",country," source status and original audit notes</summary>"),
    html_table(input_status[[country]]),paste0("<pre>",html_escape(source_notes[[country]]),"</pre></details>"))
  parts <- c(parts,"<h2>Methods and limits</h2><ul>",paste0("<li>",html_escape(notes),"</li>"),"</ul>",
    "<details><summary>Native-frequency baseline diagnostics</summary>",html_table(diagnostic_table[,c("label","frequency","n_baseline","baseline_slope_pp_year","LjungBox_p","baseline_outside_0_100")]),"</details>",
    "<h2>Sources and exports</h2><ul>",
    "<li><a href='https://digital.nhs.uk/data-and-information/publications/statistical/nhs-sickness-absence-rates'>NHS England sickness absence statistics</a>.</li>",
    "<li><a href='https://stats.gov.wales/en-GB/d4578da5-380f-439f-a056-8c5dd74d6509'>StatsWales sickness absence</a> and <a href='https://www.gov.wales/nhs-workforce-staff-directly-employed-nhs-wales-quality-report-html'>Welsh workforce quality report</a>.</li>",
    "<li><a href='https://www.cbs.nl/nl-nl/cijfers/detail/80072ned'>CBS 80072ned</a> and <a href='https://dataderden.cbs.nl/ODataApi/OData/24015NED'>AZW 24015NED</a>; precise source vintages are recorded in the input exports.</li>",
    "<li><a href='https://cran.r-project.org/web/packages/strucchange/strucchange.pdf'>strucchange documentation</a>.</li></ul>",
    "<p>The run folder contains a standalone HTML report, PNG/PDF charts, matched quarterly data, native NHS monthly fits, annual summaries, 2018-2020 detail, coefficients, diagnostics, break dates/BIC/sensitivities, descriptive correlations and contrasts, models, configuration, session information, input copies and audits. No external assets are needed to view this HTML.</p>",
    paste0("<small>Generated ",format(Sys.time(),tz="UTC",usetz=TRUE),"; run ",html_escape(basename(run_dir)),".</small></body></html>"))
  writeLines(parts,file.path(run_dir,"report.html"),useBytes=TRUE)
  write_csv(data.frame(error=character()),file.path(run_dir,"model_errors.csv"))
  writeLines(format(Sys.time(),tz="UTC",usetz=TRUE),file.path(run_dir,"RUN_COMPLETE"))
  completed <- TRUE
  message("Comparison saved: ",normalizePath(run_dir),"\nOpen report.html in your browser.")
  invisible(list(run_dir=normalizePath(run_dir),data=all_data,annual=annual,breaks=breaks,registry=registry,availability=availability))
}
if(isTRUE(RUN_ANALYSIS)) three_country_result <- main()
