# California wastewater RNA -> county PCR positivity, high-testing validation
# Standalone: source("california_wastewater_positivity_prediction.R")
# Install once in a clean R session: install.packages(c("curl","jsonlite"))
# Historical public data; no automatic package installation, httr2 or Pandoc.
# First download is approximately 200 MB. Raw sources are cached, never replaced
# by proxy metrics. Outputs go to analysis/california_wastewater_prediction.
# All counties are screened. Plants/labs/methods remain separate, never pooled
# across incompatible concentration scales. The target is clinical positivity,
# not infection prevalence. All eligibility rules precede model evaluation.
# Version 1.0 | 23 October 2026

RUN_ANALYSIS <- getOption("capredict.autorun", TRUE)
OUTPUT_ROOT <- file.path("analysis", "california_wastewater_prediction")
REFRESH_DOWNLOADS <- FALSE
OPEN_HTML <- interactive()
LOCAL_WASTEWATER_CSV <- ""
LOCAL_CLINICAL_CSV <- ""
ANALYSIS_START <- as.Date("2020-02-01")
ANALYSIS_END <- as.Date("2023-05-10")
REFERENCE_START <- as.Date("2021-07-01")
REFERENCE_END <- as.Date("2021-12-31")
MIN_MONTH_REPORTS <- 3L
MIN_MONTH_SPAN_DAYS <- 14L
MIN_REFERENCE_MONTHS <- 4L
MIN_WEEKLY_TESTS <- 500
MIN_WEEKLY_TESTS_PER_100K <- 1000
PRIMARY_REFERENCE_FRACTION <- 0.75
SENSITIVITY_FRACTIONS <- c(0.5, 0.75, 1)
MIN_SAMPLING_DAYS <- 2L
MIN_CORRELATION_WEEKS <- 12L
MIN_TRAIN_WEEKS <- 16L
TUNING_WEEKS <- 8L
TEST_WEEKS <- 8L
EMBARGO_DAYS <- 14L
BOOTSTRAP_REPLICATES <- 999L
BLOCK_WEEKS <- 3L
RANDOM_SEED <- 20261023L
MAX_ACCEPTABLE_MAE_PP <- 3
MIN_MAE_SKILL <- 0.10
MAX_ACCEPTABLE_ABS_BIAS_PP <- 2
FOCUS_FIPS <- c("06075", "06085", "06067")
CA_LANDING <- "https://data.ca.gov/dataset/covid-19-wastewater-surveillance-data-california"
CA_API <- "https://publish.data.ca.gov/api/3/action/package_show?id=b8c6ee3b-539d-4d62-8fa2-c7cd17c16656"
WW_URL <- "https://publish.data.ca.gov/dataset/b8c6ee3b-539d-4d62-8fa2-c7cd17c16656/resource/16bb2698-c243-4b66-a6e8-4861ee66f8bf/download/master-covid-public.csv"
DICTIONARY_URL <- "https://publish.data.ca.gov/dataset/b8c6ee3b-539d-4d62-8fa2-c7cd17c16656/resource/e5b302ff-2a51-4033-9341-a6a1205075c5/download/nwss-data-dictionary_v3_1_0_20211208_metadata.csv"
CLINICAL_LANDING <- "https://data.chhs.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state"
CLINICAL_URL <- "https://data.chhs.ca.gov/dataset/f333528b-4d38-4814-bebb-12db1f10f535/resource/046cdd2b-31e5-4d34-9ed3-b48cdbc4be7a/download/covid19cases_test.csv"
CLINICAL_METADATA <- "https://data.chhs.ca.gov/api/3/action/package_show?id=covid-19-time-series-metrics-by-county-and-state"
CLINICAL_SCHEMA <- "https://data.chhs.ca.gov/api/3/action/datastore_search?resource_id=046cdd2b-31e5-4d34-9ed3-b48cdbc4be7a&limit=0"
# Official California county names and FIPS. Out-of-state/unknown rows are audited,
# never assigned to counties or mixed into denominators.
COUNTIES <- c(
  `06001` = "Alameda",
  `06003` = "Alpine",
  `06005` = "Amador",
  `06007` = "Butte",
  `06009` = "Calaveras",
  `06011` = "Colusa",
  `06013` = "Contra Costa",
  `06015` = "Del Norte",
  `06017` = "El Dorado",
  `06019` = "Fresno",
  `06021` = "Glenn",
  `06023` = "Humboldt",
  `06025` = "Imperial",
  `06027` = "Inyo",
  `06029` = "Kern",
  `06031` = "Kings",
  `06033` = "Lake",
  `06035` = "Lassen",
  `06037` = "Los Angeles",
  `06039` = "Madera",
  `06041` = "Marin",
  `06043` = "Mariposa",
  `06045` = "Mendocino",
  `06047` = "Merced",
  `06049` = "Modoc",
  `06051` = "Mono",
  `06053` = "Monterey",
  `06055` = "Napa",
  `06057` = "Nevada",
  `06059` = "Orange",
  `06061` = "Placer",
  `06063` = "Plumas",
  `06065` = "Riverside",
  `06067` = "Sacramento",
  `06069` = "San Benito",
  `06071` = "San Bernardino",
  `06073` = "San Diego",
  `06075` = "San Francisco",
  `06077` = "San Joaquin",
  `06079` = "San Luis Obispo",
  `06081` = "San Mateo",
  `06083` = "Santa Barbara",
  `06085` = "Santa Clara",
  `06087` = "Santa Cruz",
  `06089` = "Shasta",
  `06091` = "Sierra",
  `06093` = "Siskiyou",
  `06095` = "Solano",
  `06097` = "Sonoma",
  `06099` = "Stanislaus",
  `06101` = "Sutter",
  `06103` = "Tehama",
  `06105` = "Trinity",
  `06107` = "Tulare",
  `06109` = "Tuolumne",
  `06111` = "Ventura",
  `06113` = "Yolo",
  `06115` = "Yuba")

fail <- function(...) stop(..., call. = FALSE)
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
clean <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
lc <- function(x) tolower(clean(x))
yes <- function(x) lc(x) %in% c("yes", "true", "1")
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
num <- function(x, label) {
  s <- clean(x)
  # US thousands separators are accepted only in strictly grouped numeric forms.
  grouped <- grepl("^[+-]?[0-9]{1,3}(,[0-9]{3})+([.][0-9]+)?$", s)
  s[grouped] <- gsub(",", "", s[grouped], fixed=TRUE)
  s[s %in% c("None","NONE","NA","N/A")] <- ""
  out <- suppressWarnings(as.numeric(s))
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
WW_COLUMNS <- c("wwtp_name", "county_names", "population_served", "sample_id", "lab_id",
  "sample_collect_date", "sample_type", "sample_matrix", "sample_location",
  "concentration_method", "extraction_method", "pcr_type", "lod_sewage",
  "hum_frac_target_mic", "hum_frac_mic_conc", "hum_frac_mic_unit", "quality_flag",
  "qc_ignore", "dashboard_ignore", "analysis_ignore", "major_lab_method",
  "major_lab_method_desc", "pcr_target", "pcr_gene_target", "pcr_target_units",
  "pcr_target_avg_conc", "pcr_target_below_lod", "data_source", "site_id",
  "created_at", "county_treatmentplant", "label_name")
read_wastewater <- function(path) {
  h<-names(read.csv(path,nrows=0L,check.names=FALSE));check_cols(setNames(as.data.frame(matrix(nrow=0,ncol=length(h))),h),WW_COLUMNS)
  x<-read.csv(path,colClasses=ifelse(h%in%WW_COLUMNS,"character","NULL"),check.names=FALSE,stringsAsFactors=FALSE,na.strings="",fileEncoding="UTF-8-BOM")
  x$source_row<-seq_len(nrow(x));x[]<-lapply(x,clean)
  x$date<-as.Date(x$sample_collect_date,format="%Y-%m-%d")
  if(anyNA(x$date))fail("Invalid wastewater collection dates; inspect the source.")
  x<-x[x$date>=ANALYSIS_START & x$date<=ANALYSIS_END,,drop=FALSE]
  x$county_treatmentplant_original<-x$county_treatmentplant
  x$county_fips<-names(COUNTIES)[match(sub(" county$","",lc(x$county_treatmentplant)),tolower(COUNTIES))]
  x$county_assignment<-ifelse(is.na(x$county_fips),"unassigned","plant county field")
  tokens<-regmatches(x$county_names,gregexpr("[0-9]{5}",x$county_names))
  one<-vapply(tokens,function(z)length(unique(z))==1L,logical(1))
  inferred<-vapply(tokens,function(z)if(length(unique(z))==1L)unique(z)else "",character(1))
  use<-is.na(x$county_fips)&one&inferred%in%names(COUNTIES)
  x$county_fips[use]<-inferred[use];x$county_assignment[use]<-"single explicit catchment county"
  x$county<-unname(COUNTIES[x$county_fips])
  x$numeric_thousands_separator<-grepl(",",x$pcr_target_avg_conc,fixed=TRUE)
  x
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

wp_cor <- function(x,y,method="pearson") if(length(x)>=3L&&sd(x)>0&&sd(y)>0)unname(cor(x,y,method=method))else NA_real_
wp_month <- function(x)format(x,"%Y-%m")
wp_models <- list(RNA_log=target_logit~logw,RNA_spline=target_logit~splines::ns(logw,df=3),
  RNA_change=target_logit~logw+change_logw,AR_baseline=target_logit~prior_logit,RNA_plus_prior=target_logit~logw+prior_logit)
wp_fit <- function(name,d) {
  if(name=="training_mean")return(list(name=name,mean=mean(d$positivity)))
  if(name=="persistence")return(list(name=name))
  fit<-tryCatch(lm(wp_models[[name]],data=d),error=function(e)NULL)
  if(is.null(fit)||fit$rank<length(coef(fit))||any(!is.finite(coef(fit))))return(NULL)
  list(name=name,fit=fit)
}
wp_predict <- function(fit,d) {
  if(is.null(fit))return(rep(NA_real_,nrow(d)))
  if(fit$name=="training_mean")return(rep(fit$mean,nrow(d)))
  if(fit$name=="persistence")return(d$prior_positivity)
  tryCatch(100*plogis(predict(fit$fit,newdata=d)),error=function(e)rep(NA_real_,nrow(d)))
}
wp_metrics <- function(y,hat,baseline) {
  e<-hat-y;sse<-sum(e^2);sst<-sum((y-mean(y))^2);base_mae<-mean(abs(baseline-y))
  c(mae=mean(abs(e)),rmse=sqrt(mean(e^2)),bias=mean(e),r2=if(sst>0)1-sse/sst else NA_real_,
    pearson=wp_cor(hat,y),mae_skill=if(base_mae>0)1-mean(abs(e))/base_mae else NA_real_)
}
wp_bootstrap <- function(d,hat) {
  n<-nrow(d);L<-min(BLOCK_WEEKS,n)
  starts<-which(vapply(seq_len(n),function(i)i+L-1L<=n&&as.integer(d$target_end[min(n,i+L-1L)]-d$target_end[i])==7L*(L-1L),logical(1)))
  if(length(starts)<2L)return(c(mae_lo=NA,mae_hi=NA,skill_lo=NA,skill_hi=NA))
  covered<-unique(unlist(lapply(starts,function(i)i+seq_len(L)-1L)))
  # Never silently omit isolated test weeks from uncertainty calculations.
  if(length(covered)<n)return(c(mae_lo=NA,mae_hi=NA,skill_lo=NA,skill_hi=NA))
  b<-replicate(BOOTSTRAP_REPLICATES,{
    ii<-unlist(lapply(sample(starts,ceiling(n/L),replace=TRUE),function(i)i+seq_len(L)-1L))[seq_len(n)]
    m<-wp_metrics(d$positivity[ii],hat[ii],d$prior_positivity[ii]);m[c("mae","mae_skill")]
  })
  c(mae_lo=unname(quantile(b[1,],.025,na.rm=TRUE)),mae_hi=unname(quantile(b[1,],.975,na.rm=TRUE)),
    skill_lo=unname(quantile(b[2,],.025,na.rm=TRUE)),skill_hi=unname(quantile(b[2,],.975,na.rm=TRUE)))
}
wp_escape <- function(x) {x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
wp_table <- function(d,digits=3) {
  if(!nrow(d))return("<p>No rows meet these requirements.</p>")
  for(n in names(d))if(is.numeric(d[[n]]))d[[n]]<-ifelse(is.finite(d[[n]]),format(round(d[[n]],digits),trim=TRUE),"—")
  rows<-apply(d,1,function(x)paste0("<tr><td>",paste(wp_escape(x),collapse="</td><td>"),"</td></tr>"))
  paste0("<div class='scroll'><table><thead><tr><th>",paste(wp_escape(names(d)),collapse="</th><th>"),"</th></tr></thead><tbody>",paste(rows,collapse=""),"</tbody></table></div>")
}
wp_image <- function(path) paste0("<img src='data:image/png;base64,",gsub("[\r\n]","",jsonlite::base64_enc(readBin(path,"raw",n=file.info(path)$size))),"' alt='Surveillance analysis plot'>")

ca_clinical <- function(path,run) {
  x<-read_text_csv(path);check_cols(x,c('date','area','area_type','population','total_tests','positive_tests'))
  x$source_row<-seq_len(nrow(x));x$fips<-names(COUNTIES)[match(x$area,COUNTIES)]
  date_text<-clean(x$date);x$date<-as.Date(date_text,format='%Y-%m-%d')
  x$disposition<-ifelse(x$area_type!='County'|is.na(x$fips),'not one of the 58 California counties',ifelse(is.na(x$date),'missing date','eligible'))
  outside<-!is.na(x$date)&(x$date<ANALYSIS_START|x$date>ANALYSIS_END)
  x$disposition[x$disposition=='eligible'&outside]<-'outside analysis window'
  for(n in c('population','total_tests','positive_tests'))x[[n]]<-num(x[[n]],n)
  invalid<-!is.finite(x$total_tests)|!is.finite(x$positive_tests)|!is.finite(x$population)|x$population<=0|x$total_tests<0|x$positive_tests<0|x$positive_tests>x$total_tests
  x$disposition[x$disposition=='eligible'&invalid]<-'missing or invalid count/population'
  write_csv(x,file.path(run,'data','clinical_daily_audit.csv'))
  z<-x[x$disposition=='eligible',];if(anyDuplicated(z[c('fips','date')]))fail('Duplicate county specimen dates.')
  z$week_end<-day_end(z$date);groups<-split(seq_len(nrow(z)),paste(z$fips,z$week_end))
  vals<-t(vapply(groups,function(ii)c(days=length(unique(z$date[ii])),tests=sum(z$total_tests[ii]),positive=sum(z$positive_tests[ii]),population=median(z$population[ii]),population_versions=length(unique(z$population[ii]))),numeric(5)))
  first<-vapply(groups,`[`,integer(1),1L)
  w<-data.frame(fips=z$fips[first],week_end=z$week_end[first],vals,row.names=NULL)
  w$complete<-w$days==7L&w$population_versions==1L&w$tests>0
  w$positivity<-ifelse(w$complete,100*w$positive/w$tests,NA_real_)
  w$rate<-ifelse(w$complete,1e5*w$tests/w$population,NA_real_)
  w<-w[w$week_end<=ANALYSIS_END,];w$month<-wp_month(w$week_end)
  groups<-split(which(w$complete),paste(w$fips[w$complete],w$month[w$complete]))
  v<-t(vapply(groups,function(ii)c(reports=length(ii),span_days=as.integer(max(w$week_end[ii])-min(w$week_end[ii])),median_weekly_tests=median(w$tests[ii]),median_weekly_rate=median(w$rate[ii])),numeric(4)))
  j<-vapply(groups,`[`,integer(1),1L);m<-data.frame(fips=w$fips[j],month=w$month[j],v,row.names=NULL)
  m$adequate<-m$reports>=MIN_MONTH_REPORTS&m$span_days>=MIN_MONTH_SPAN_DAYS
  ref<-m[m$adequate&m$month>=wp_month(REFERENCE_START)&m$month<=wp_month(REFERENCE_END),]
  ref<-do.call(rbind,lapply(split(ref,ref$fips),function(a)data.frame(fips=a$fips[1],reference_months=nrow(a),reference_rate=median(a$median_weekly_rate))))
  m<-merge(m,ref,by='fips',all.x=TRUE,sort=FALSE);m$reference_ok<-!is.na(m$reference_months)&m$reference_months>=MIN_REFERENCE_MONTHS&is.finite(m$reference_rate)&m$reference_rate>0
  m$relative_volume<-m$median_weekly_rate/m$reference_rate
  write_csv(w,file.path(run,'data','clinical_weekly.csv'));write_csv(m,file.path(run,'data','monthly_testing.csv'))
  list(weekly=w,months=m,audit=x)
}
ca_qualify <- function(clinical,fraction) {
  m<-clinical$months;w<-clinical$weekly
  m$high_testing<-m$adequate&m$reference_ok&m$relative_volume>=fraction&m$median_weekly_tests>=MIN_WEEKLY_TESTS&m$median_weekly_rate>=MIN_WEEKLY_TESTS_PER_100K
  m$high_testing[is.na(m$high_testing)]<-FALSE
  j<-match(paste(w$fips,w$month),paste(m$fips,m$month))
  w$high_testing<-m$high_testing[j]&w$complete&w$tests>=MIN_WEEKLY_TESTS&w$rate>=MIN_WEEKLY_TESTS_PER_100K
  w$high_testing[is.na(w$high_testing)]<-FALSE
  list(weekly=w,months=m)
}
ca_weekly <- function(s) {
  # Average sample replicates within each date, then give sampled dates equal
  # weight within Thursday-Wednesday weeks. No interpolation of unsampled days.
  groups<-split(seq_len(nrow(s)),paste(s$series_id,s$date))
  v<-t(vapply(groups,function(ii){matched<-ii[is.finite(s$pmmov_ratio[ii])&is.finite(s$rna_primary[ii])];
    c(rna=mean_ok(s$rna_primary[ii]),detected=mean_ok(s$rna_detected_only[ii]),ratio=mean_ok(s$pmmov_ratio[matched]),matched_rna=mean_ok(s$rna_primary[matched]))},numeric(4)))
  j<-vapply(groups,`[`,integer(1),1L);d<-data.frame(series_id=s$series_id[j],county_fips=s$county_fips[j],date=s$date[j],v,row.names=NULL);d$week_end<-day_end(d$date)
  groups<-split(seq_len(nrow(d)),paste(d$series_id,d$week_end))
  v<-t(vapply(groups,function(ii)c(sampling_days=sum(is.finite(d$rna[ii])),detected_days=sum(is.finite(d$detected[ii])),matched_days=sum(is.finite(d$ratio[ii])),rna_mean=mean_ok(d$rna[ii]),detected_mean=mean_ok(d$detected[ii]),ratio_mean=mean_ok(d$ratio[ii]),matched_rna_mean=mean_ok(d$matched_rna[ii])),numeric(7)))
  j<-vapply(groups,`[`,integer(1),1L);w<-data.frame(series_id=d$series_id[j],county_fips=d$county_fips[j],week_end=d$week_end[j],v,row.names=NULL)
  w<-w[w$week_end<=ANALYSIS_END,];w<-w[order(w$series_id,w$week_end),]
  list(daily=d,weekly=w)
}
ca_pairs <- function(w,clinical,registry,lag,variant) {
  z<-w;z$fips<-z$county_fips;z$wastewater_end<-z$week_end;z$target_end<-z$week_end+7L*lag;z$target_start<-z$target_end-6L
  z$task<-if(lag==0)'same_week_estimation'else '7_day_sample_lead'
  value_col<-switch(variant,raw='rna_mean',raw_matched='matched_rna_mean',pmmov='ratio_mean',detected='detected_mean')
  day_col<-switch(variant,raw='sampling_days',raw_matched='matched_days',pmmov='matched_days',detected='detected_days')
  z$exposure<-z[[value_col]];z$adequate_sampling<-z[[day_col]]>=MIN_SAMPLING_DAYS&is.finite(z$exposure)&z$exposure>0
  j<-match(paste(z$series_id,z$week_end-7L),paste(z$series_id,z$week_end));z$previous_exposure<-z$exposure[j];z$previous_sampling<-z$adequate_sampling[j]
  k<-match(paste(z$fips,z$target_end),paste(clinical$fips,clinical$week_end));k1<-match(paste(z$fips,z$target_end-7L),paste(clinical$fips,clinical$week_end))
  z$positivity<-clinical$positivity[k];z$high_testing<-clinical$high_testing[k];z$weekly_tests<-clinical$tests[k];z$weekly_rate<-clinical$rate[k]
  z$prior_positivity<-clinical$positivity[k1];z$prior_high_testing<-clinical$high_testing[k1]
  z$geography_class<-registry$geography_class[match(z$series_id,registry$series_id)]
  z$pair_eligible<-z$adequate_sampling&is.finite(z$positivity)&z$high_testing&z$target_end<=ANALYSIS_END
  z$model_eligible<-z$pair_eligible&z$previous_sampling&is.finite(z$previous_exposure)&z$previous_exposure>0&is.finite(z$prior_positivity)&z$prior_high_testing
  for(n in c('pair_eligible','model_eligible','high_testing','prior_high_testing'))z[[n]][is.na(z[[n]])]<-FALSE
  z$logw<-log10(z$exposure);z$change_logw<-log10(z$exposure)-log10(z$previous_exposure)
  logit<-function(a)qlogis(pmin(.9995,pmax(.0005,a/100)))
  z$target_logit<-logit(z$positivity);z$prior_logit<-logit(z$prior_positivity);z$variant<-variant
  z
}
ca_split <- function(d) {
  n<-nrow(d);if(n<MIN_TRAIN_WEEKS+TUNING_WEEKS+TEST_WEEKS)return(NULL)
  test<-tail(seq_len(n),TEST_WEEKS);before<-which(d$target_end<min(d$target_end[test])-EMBARGO_DAYS)
  if(length(before)<TUNING_WEEKS)return(NULL)
  tune<-tail(before,TUNING_WEEKS);train<-which(d$target_end<min(d$target_end[tune])-EMBARGO_DAYS)
  if(length(train)<MIN_TRAIN_WEEKS)return(NULL)
  list(train=train,tune=tune,test=test)
}
ca_correlations <- function(z) {
  ans<-lapply(split(z,z$series_id),function(a){d<-a[a$pair_eligible,];d<-d[order(d$target_end),];n<-nrow(d);prior<-match(d$target_end-7L,d$target_end);ok<-!is.na(prior)
    data.frame(series_id=a$series_id[1],fips=a$fips[1],task=a$task[1],variant=a$variant[1],geography_class=a$geography_class[1],weeks=n,
      first=if(n)as.character(min(d$target_end))else NA_character_,last=if(n)as.character(max(d$target_end))else NA_character_,
      pearson_logrna=if(n>=MIN_CORRELATION_WEEKS)wp_cor(d$logw,d$positivity)else NA_real_,
      spearman=if(n>=MIN_CORRELATION_WEEKS)wp_cor(d$exposure,d$positivity,'spearman')else NA_real_,consecutive_changes=sum(ok),
      change_correlation=if(sum(ok)>=MIN_CORRELATION_WEEKS)wp_cor(d$logw[ok]-d$logw[prior[ok]],d$positivity[ok]-d$positivity[prior[ok]])else NA_real_)
  });do.call(rbind,ans)
}
ca_evaluate <- function(z,fraction) {
  summaries<-predictions<-scores<-audits<-partitions<-models<-list()
  for(id in unique(z$series_id)) {
    a<-z[z$series_id==id,];d<-a[a$model_eligible,];d<-d[order(d$target_end),];sp<-ca_split(d)
    au<-data.frame(series_id=id,fips=a$fips[1],fraction=fraction,variant=a$variant[1],task=a$task[1],geography_class=a$geography_class[1],high_testing_pairs=sum(a$pair_eligible),common_model_weeks=nrow(d),status='insufficient qualifying weeks',train_weeks=0L,tuning_weeks=0L,test_weeks=0L)
    if(is.null(sp)){audits[[id]]<-au;next}
    tr<-d[sp$train,];va<-d[sp$tune,];te<-d[sp$test,]
    names_fit<-c(names(wp_models),'training_mean','persistence');fits<-setNames(lapply(names_fit,function(n)wp_fit(n,tr)),names_fit)
    val<-sapply(fits,function(fit){hat<-wp_predict(fit,va);if(all(is.finite(hat)))mean(abs(hat-va$positivity))else Inf})
    ww<-c('RNA_log','RNA_spline','RNA_change');winner<-ww[which.min(val[ww])]
    if(!is.finite(val[winner])){au$status<-'model estimation failure';audits[[id]]<-au;next}
    dev<-rbind(tr,va);final<-setNames(lapply(names_fit,function(n)wp_fit(n,dev)),names_fit);hats<-sapply(final,wp_predict,d=te)
    if(any(!is.finite(hats[,winner]))){au$status<-'test prediction failure';audits[[id]]<-au;next}
    hat<-hats[,winner];metrics<-wp_metrics(te$positivity,hat,te$prior_positivity);ci<-wp_bootstrap(te,hat)
    passes<-is.finite(metrics['mae_skill'])&&metrics['mae']<=MAX_ACCEPTABLE_MAE_PP&&metrics['mae_skill']>=MIN_MAE_SKILL&&abs(metrics['bias'])<=MAX_ACCEPTABLE_ABS_BIAS_PP
    summaries[[id]]<-data.frame(series_id=id,fips=a$fips[1],fraction=fraction,variant=a$variant[1],task=a$task[1],geography_class=a$geography_class[1],winner=winner,
      train_weeks=nrow(tr),tuning_weeks=nrow(va),test_weeks=nrow(te),train_first=min(tr$target_end),train_last=max(tr$target_end),tune_first=min(va$target_end),tune_last=max(va$target_end),test_first=min(te$target_end),test_last=max(te$target_end),
      as.list(metrics),as.list(ci),passes_point_criteria=passes,positive_skill_interval=passes&&is.finite(ci['skill_lo'])&&ci['skill_lo']>0,
      extrapolation_weeks=sum(te$logw<min(dev$logw)|te$logw>max(dev$logw)),row.names=NULL)
    for(n in names_fit)if(all(is.finite(hats[,n])))scores[[paste(id,n)]]<-data.frame(series_id=id,fips=a$fips[1],fraction=fraction,variant=a$variant[1],task=a$task[1],model=n,selected_ww=n==winner,validation_mae=val[n],as.list(wp_metrics(te$positivity,hats[,n],te$prior_positivity)),row.names=NULL)
    pr<-te[c('series_id','fips','target_start','target_end','wastewater_end','exposure','positivity','prior_positivity','weekly_tests','weekly_rate','geography_class','task','variant')]
    pr$fraction<-fraction;pr$model<-winner;pr$prediction<-hat;pr$error_pp<-hat-pr$positivity;pr$within_development_exposure<-te$logw>=min(dev$logw)&te$logw<=max(dev$logw);predictions[[id]]<-pr
    dd<-d[c('series_id','fips','target_end','wastewater_end','exposure','positivity')];dd$partition<-'embargo';dd$partition[sp$train]<-'train';dd$partition[sp$tune]<-'tune';dd$partition[sp$test]<-'test';dd$fraction<-fraction;dd$variant<-a$variant[1];dd$task<-a$task[1];partitions[[id]]<-dd
    au$status<-'evaluated';au$train_weeks<-nrow(tr);au$tuning_weeks<-nrow(va);au$test_weeks<-nrow(te);audits[[id]]<-au
    models[[id]]<-list(selected_model=final[[winner]],exposure_range=range(dev$exposure),fraction=fraction,variant=a$variant[1],task=a$task[1],fips=a$fips[1],test_metrics=metrics)
  }
  bind<-function(x)if(length(x))do.call(rbind,x)else data.frame()
  list(summary=bind(summaries),predictions=bind(predictions),scores=bind(scores),audit=bind(audits),partitions=bind(partitions),models=models)
}

ca_label <- function(id,registry) {
  z<-registry[match(id,registry$series_id),]
  paste(z$label,z$lab_id,z$pcr_gene_target,paste0('[',id,']'))
}
ca_testing_plot <- function(f,clinical,path) {
  q<-ca_qualify(clinical,PRIMARY_REFERENCE_FRACTION);w<-q$weekly[q$weekly$fips==f,];m<-q$months[q$months$fips==f,]
  w<-w[order(w$week_end),];m<-m[order(m$month),];dates<-as.Date(paste0(m$month,'-15'));xr<-range(c(w$week_end,dates));threshold<-pmax(MIN_WEEKLY_TESTS_PER_100K,PRIMARY_REFERENCE_FRACTION*m$reference_rate)
  png(path,1400,850,res=130);on.exit(dev.off(),add=TRUE);par(mfrow=c(2,1),mar=c(4,5,3,1),las=1)
  plot(dates,m$median_weekly_rate,type='b',pch=16,col='#748894',xlim=xr,ylim=c(0,max(c(m$median_weekly_rate,threshold),na.rm=TRUE)*1.05),xlab='',ylab='Weekly PCR tests / 100,000',main=paste(COUNTIES[f],'| testing eligibility'))
  lines(dates,threshold,col='#ba513e',lty=2);points(dates[m$high_testing],m$median_weekly_rate[m$high_testing],col='#178477',pch=16)
  legend('topright',c('Monthly median','Primary threshold','Qualifying month'),col=c('#748894','#ba513e','#178477'),pch=c(16,NA,16),lty=c(1,2,NA),bty='n',cex=.8)
  plot(w$week_end,w$positivity,type='l',col='#b44157',xlim=xr,xlab='Specimen week ending Wednesday',ylab='PCR positivity (%)',main='Observed positivity; high-testing weeks highlighted')
  points(w$week_end[w$high_testing],w$positivity[w$high_testing],pch=16,col='#178477',cex=.6)
}
ca_prediction_plot <- function(id,task,r,registry,path,variant='raw') {
  d<-r$partitions[r$partitions$series_id==id&r$partitions$variant==variant&r$partitions$fraction==PRIMARY_REFERENCE_FRACTION&r$partitions$task==task,]
  p<-r$predictions[r$predictions$series_id==id&r$predictions$variant==variant&r$predictions$fraction==PRIMARY_REFERENCE_FRACTION&r$predictions$task==task,]
  d<-d[order(d$target_end),];p<-p[order(p$target_end),]
  png(path,1450,700,res=130);on.exit(dev.off(),add=TRUE);par(mfrow=c(1,2),mar=c(5,5,4,1),las=1)
  colors<-c(train='#87959e',tune='#c18b27',test='#1d8e83',embargo='#d4dbe0')
  # Points avoid drawing across gaps in the retained high-testing observations.
  plot(d$target_end,d$positivity,pch=16,col=colors[d$partition],xlab='Clinical week end',ylab='Positivity (%)',ylim=range(c(d$positivity,p$prediction)),main=paste(ca_label(id,registry),if(variant=='pmmov')'| PMMoV normalized'else '| raw RNA'),cex.main=.8,xaxt='n')
  axis.Date(1,at=pretty(d$target_end,n=4),format='%Y-%m')
  points(p$target_end,p$prediction,pch=1,col='#2769ad',cex=1.2)
  legend('topleft',c('Training','Tuning','Test observed','Test predicted'),col=c(colors[c('train','tune','test')],'#2769ad'),pch=c(16,16,16,1),bty='n',cex=.8)
  rr<-range(c(p$positivity,p$prediction,p$prior_positivity));plot(p$positivity,p$prediction,xlim=rr,ylim=rr,pch=16,col='#2769ad',xlab='Observed positivity (%)',ylab='Predicted positivity (%)',main=gsub('_',' ',task),cex.main=.9)
  abline(0,1,lty=2,col='#87959e');points(p$positivity,p$prior_positivity,pch=4,col='#b3703b')
  legend('topleft',c('Selected RNA model','Persistence'),col=c('#2769ad','#b3703b'),pch=c(16,4),bty='n',cex=.8)
}
ca_table <- function(d,digits=3) {
  if('task'%in%names(d))d$task<-gsub('_',' ',d$task)
  names(d)<-gsub('_',' ',names(d));wp_table(d,digits)
}
ca_report <- function(run,clinical,prepared,results,correlations,sources,coverage) {
  s<-results$summary;reg<-prepared$registry
  if(nrow(s)){s$county<-unname(COUNTIES[s$fips]);s$plant<-reg$label[match(s$series_id,reg$series_id)];s$lab<-reg$lab_id[match(s$series_id,reg$series_id)]}
  primary<-if(nrow(s))s[s$variant=='raw'&s$fraction==PRIMARY_REFERENCE_FRACTION&s$geography_class=='single county',]else data.frame()
  now<-if(nrow(primary))primary[primary$task=='same_week_estimation',]else data.frame()
  lead<-if(nrow(primary))primary[primary$task=='7_day_sample_lead',]else data.frame()
  conclusion<-if(!nrow(now))'No single-county plant/method series supports the full primary training, tuning and test design.'else paste0(nrow(now),' single-county plant/method series across ',length(unique(now$fips)),' counties support the primary same-week evaluation. The fixed point-performance criteria are met by ',sum(now$passes_point_criteria),' series, with ',sum(now$positive_skill_interval),' also having a positive lower bootstrap skill bound. Median held-out MAE across these series is ',round(median(now$mae),2),' percentage points. This is an exploratory assessment with eight held-out observations per series, not a validated statewide conversion to positivity.')
  au<-results$audit;groups<-split(seq_len(nrow(au)),paste(au$fraction,au$variant,au$task,au$geography_class))
  overview<-do.call(rbind,lapply(groups,function(ii){a<-au[ii,];j<-if(nrow(s))s$fraction==a$fraction[1]&s$variant==a$variant[1]&s$task==a$task[1]&s$geography_class==a$geography_class[1]else logical();b<-s[j,,drop=FALSE]
    data.frame(fraction=a$fraction[1],variant=a$variant[1],task=a$task[1],geography=a$geography_class[1],series_screened=nrow(a),series_evaluated=sum(a$status=='evaluated'),counties_evaluated=length(unique(b$fips)),median_MAE_pp=if(nrow(b))median(b$mae)else NA_real_,median_skill_pct=if(nrow(b))100*median(b$mae_skill,na.rm=TRUE)else NA_real_,pass_point_criteria=if(nrow(b))sum(b$passes_point_criteria)else 0L,positive_skill_bound=if(nrow(b))sum(b$positive_skill_interval)else 0L)}))
  overview<-overview[order(overview$fraction,overview$variant,overview$task,overview$geography),]
  show<-function(a)if(!nrow(a))a else a[c('series_id','county','plant','lab','task','winner','train_weeks','tuning_weeks','test_weeks','test_first','test_last','mae','bias','r2','mae_skill','skill_lo','skill_hi','passes_point_criteria','extrapolation_weeks')]
  extra<-if(nrow(s))s[s$variant=='raw'&s$fraction==PRIMARY_REFERENCE_FRACTION&s$geography_class!='single county',]else data.frame()
  ct<-correlations[correlations$variant=='raw'&correlations$fraction==PRIMARY_REFERENCE_FRACTION&correlations$task=='same_week_estimation',]
  ct$county<-unname(COUNTIES[ct$fips]);ct$plant<-reg$label[match(ct$series_id,reg$series_id)];ct<-ct[order(ct$county,ct$series_id),]
  disp<-as.data.frame(table(prepared$audit$disposition),stringsAsFactors=FALSE);names(disp)<-c('disposition','rows')
  focus<-list();for(f in intersect(FOCUS_FIPS,clinical$weekly$fips)){path<-file.path(run,'plots',paste0('testing_',f,'.png'));ca_testing_plot(f,clinical,path);focus[[f]]<-wp_image(path)}
  plots<-character();plot_rows<-if(nrow(s))s[s$variant=='raw'&s$fraction==PRIMARY_REFERENCE_FRACTION,]else data.frame()
  if(nrow(plot_rows))for(i in seq_len(nrow(plot_rows))){id<-plot_rows$series_id[i];task<-plot_rows$task[i];path<-file.path(run,'plots',paste0(id,'_',task,'.png'));ca_prediction_plot(id,task,results,reg,path);plots<-c(plots,paste0("<details><summary>",wp_escape(paste(ca_label(id,reg),task,plot_rows$geography_class[i],sep=' · ')),"</summary>",wp_image(path),"</details>"))}
  # Correlations are descriptive across all eligible weeks, never model selection.
  cc<-ct[ct$geography_class=='single county'&is.finite(ct$pearson_logrna),]
  cor_text<-if(nrow(cc))paste0(nrow(cc),' single-county series across ',length(unique(cc$fips)),' counties have at least ',MIN_CORRELATION_WEEKS,' eligible paired weeks. Their median descriptive Pearson correlation between log10 RNA and positivity is ',round(median(cc$pearson_logrna),3),'.')else 'No series meets the descriptive correlation minimum.'
  raw_sens<-if(nrow(s))s[s$variant=='raw'&s$fraction!=PRIMARY_REFERENCE_FRACTION,]else data.frame()
  variant_sens<-if(nrow(s))s[s$variant!='raw',]else data.frame()
  candidates<-if(nrow(s))s[s$variant=='pmmov'&s$fraction==PRIMARY_REFERENCE_FRACTION&s$geography_class=='single county'&s$task=='same_week_estimation'&s$passes_point_criteria,]else data.frame()
  candidate_text<-if(nrow(candidates))paste0('The prespecified PMMoV sensitivity identifies ',nrow(candidates),if(nrow(candidates)==1)' local candidate meeting the point-performance criteria. A positive lower skill bound is present for 'else ' local candidates meeting the point-performance criteria. A positive lower skill bound is present for ',sum(candidates$positive_skill_interval),' of these candidates. These candidates are highlighted after reviewing the test results and require new external validation; they do not replace the primary raw-RNA analysis.')else 'No PMMoV-normalized same-week model meets the point-performance criteria under the primary volume rule.'
  candidate_images<-character()
  if(nrow(candidates))for(i in seq_len(nrow(candidates))){id<-candidates$series_id[i];path<-file.path(run,'plots',paste0(id,'_pmmov_candidate.png'));ca_prediction_plot(id,'same_week_estimation',results,reg,path,'pmmov');candidate_images<-c(candidate_images,wp_image(path))}
  # On matched raw/PMMoV samples, confirm identical test weeks before comparing.
  matched<-results$predictions[results$predictions$variant%in%c('raw_matched','pmmov'),]
  comparison<-data.frame()
  if(nrow(matched)){
    mm<-merge(matched[matched$variant=='raw_matched',c('series_id','task','target_end','positivity','prediction')],matched[matched$variant=='pmmov',c('series_id','task','target_end','positivity','prediction')],by=c('series_id','task','target_end'),suffixes=c('_raw','_norm'))
    if(nrow(mm))comparison<-do.call(rbind,lapply(split(mm,paste(mm$series_id,mm$task)),function(x)data.frame(series_id=x$series_id[1],task=x$task[1],same_test_weeks=nrow(x),raw_MAE=mean(abs(x$prediction_raw-x$positivity_raw)),PMMoV_MAE=mean(abs(x$prediction_norm-x$positivity_norm)))))
  }
  write_csv(comparison,file.path(run,'data','normalization_comparison.csv'))
  h<-paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>California wastewater prediction under high testing</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#193447;background:#f4f7f9;margin:0}main{max-width:1300px;margin:auto;padding:30px}h1,h2{line-height:1.2}h1{font-size:33px}h2{margin-top:34px}.card{background:#fff;border:1px solid #dbe5ea;border-radius:10px;padding:22px;margin:20px 0}.verdict{border-left:6px solid #168476}.scroll{overflow:auto;max-height:650px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px;border-bottom:1px solid #dce5ea;text-align:left;white-space:nowrap}th{position:sticky;top:0;background:#ebf2f6}img{max-width:100%;height:auto}details{background:white;margin:12px 0;padding:15px;border:1px solid #dbe5ea;border-radius:8px}summary{cursor:pointer;font-weight:600}a{color:#176a91}.small{font-size:13px;color:#576d79}li{margin:8px 0}code{background:#e9f0f4;padding:2px 4px}</style></head><body><main>",
    "<p class='small'>California · Historical evaluation · ",utc_now(),"</p><h1>Can wastewater predict positivity during high testing?</h1><div class='card verdict'><strong>Primary conclusion</strong><p>",conclusion,"</p><p>PMMoV normalization is assessed separately below. Any local candidate identified in that sensitivity requires further validation.</p></div>",
    "<p>All 58 California counties are screened. The wastewater source contains earlier measurements than the national WVAL archive, allowing an analysis beginning in 2020. Plants, laboratories, targets and method segments are evaluated separately. The target is reported county PCR test positivity, not infection prevalence or a statewide infection count.</p>",
    "<h2>What changed from the national analysis</h2><ul><li>Measured RNA concentrations replace retrospectively standardized WVALs. No CDC category or color threshold enters the models.</li><li>CDPH's archived county PCR counts supply both positivity (100 × weekly positive tests / weekly total tests) and testing volume. This avoids mixing a CDC percentage with a separate HHS denominator. It also means these numerical results concern CDPH PCR positivity, not exactly the previous CDC NAAT target.</li><li>Both sources use the same Thursday–Wednesday weeks. Same-week estimates are retrospective; the separate seven-day sample-lead task uses wastewater from the preceding week. Collection dates do not establish operational availability or publication lead time.</li></ul>",
    "<h2>Coverage across California</h2>",ca_table(coverage,0),
    "<p>Counties with absent data or insufficient overlap remain in the audit. A plant covers a sewershed, often only part of a county. Multiple methods in one county share the same clinical outcome and are not independent replications. There is no statewide pooled wastewater concentration, pooled positivity or population-representative statewide performance claim.</p>",
    "<h2>Primary held-out results: single-county catchments</h2><p>MAE and bias are percentage points. Skill = 1 − model MAE / persistence MAE; positive values favor wastewater. The selected model is chosen on tuning data only. The untouched test is scored once after refitting on training plus tuning.</p>",ca_table(show(primary)),
    "<h2>Shared/uncertain catchments: exploratory results</h2><p>These are excluded from the primary conclusion, but retained to show sites such as San Francisco. County positivity is a geographic proxy here; the source does not provide population fractions that would justify a catchment-weighted clinical target.</p>",ca_table(show(extra)),
    "<h2>Prespecified eligibility and validation</h2><div class='card'><ul>",
    "<li>A clinical week needs all seven specimen dates with valid positive/total PCR counts and a consistent population denominator. Zero testing, missing days, missing dates and invalid counts are excluded. Tests are not deduplicated people. Population is the 2020 California Department of Finance estimate supplied in the clinical file.</li>",
    "<li>The reference is the median of monthly median complete-week testing rates in July–December 2021; at least ",MIN_REFERENCE_MONTHS," adequately observed reference months are required. A month needs at least ",MIN_MONTH_REPORTS," complete weeks whose endpoints span ",MIN_MONTH_SPAN_DAYS," days. Its median testing must reach ",MIN_WEEKLY_TESTS," weekly tests, ",MIN_WEEKLY_TESTS_PER_100K," per 100,000, and ",PRIMARY_REFERENCE_FRACTION*100,"% of its reference. Every retained target and preceding clinical week must also satisfy the absolute floors and its month's high-testing rule.</li>",
    "<li>The 50%, 75% and 100% reference thresholds are separate sensitivity analyses, not a search for the most favorable outcome. Calendar months are assigned by week end. Eligibility uses the full month's retrospective volume; it is a conditional historical analysis, not a deployable prospective eligibility classifier. High volume is an operational rule, not proof of representative testing.</li>",
    "<li>At least ",MIN_SAMPLING_DAYS," distinct wastewater sampling days are required per week. Same-day replicates are averaged first, then sampling days receive equal weight. Missing weeks are never interpolated. Log10 RNA enters the models; source log10 concentrations are converted to linear units before averaging.</li>",
    "<li>Within each plant/method, the last ",TEST_WEEKS," qualifying weeks form the untouched test, the preceding ",TUNING_WEEKS," choose the model, and at least ",MIN_TRAIN_WEEKS," earlier weeks train it. A ",EMBARGO_DAYS,"-day calendar embargo separates partitions. Test periods may contain gaps and differ across series. No random week split is used.</li>",
    "<li>Wastewater-only candidates: linear log10 RNA, a natural spline with three degrees of freedom, and log10 RNA plus its previous-week change. The response is logit positivity, clipped at 0.05%/99.95% only for transformation. Selection minimizes tuning MAE. Mean positivity, persistence, autoregression and RNA plus prior positivity are comparison models, all scored on identical weeks within an evaluation.</li>",
    "<li>For the seven-day sample-lead task, persistence uses positivity from the wastewater collection week. For same-week estimation, it uses the previous clinical week. Clinical reporting delays are unknown; the comparison concerns sample-date information, not a live forecasting system.</li>",
    "<li>Exploratory pass criteria: MAE ≤ ",MAX_ACCEPTABLE_MAE_PP," points, ≥ ",MIN_MAE_SKILL*100,"% improvement over persistence, and |bias| ≤ ",MAX_ACCEPTABLE_ABS_BIAS_PP," points. Paired moving blocks of ",BLOCK_WEEKS," consecutive weeks and ",BOOTSTRAP_REPLICATES," replicates give nominal 95% uncertainty for MAE and skill. Intervals are omitted if any test observation cannot belong to a block. Eight test observations provide limited precision; these are pointwise, unadjusted intervals, not individual prediction intervals.</li></ul></div>",
    "<h2>Testing-volume sensitivity overview</h2>",ca_table(overview[overview$variant=='raw',]),
    "<details><summary>Held-out results at other volume thresholds</summary>",ca_table(raw_sens),"</details>",
    "<h2>PMMoV and nondetect sensitivity</h2><p>The primary uses raw RNA. PMMoV is a fecal marker; dividing RNA by its concentration can adjust for some variation in fecal strength. Normalization is compared to raw RNA on exactly the same samples, dates, eligible weeks and partitions. Each gets its own tuning-selected model. Nondetect exclusion is a separate sensitivity and can change the observation set. The primary specification remains fixed regardless of which sensitivity looks better.</p>",ca_table(comparison),
    "<div class='card'><strong>Local candidates for further validation</strong><p>",candidate_text,"</p>",ca_table(show(candidates)),paste(candidate_images,collapse=''),"</div>",
    "<details><summary>All normalization and nondetect results</summary>",ca_table(variant_sens),"</details>",
    "<h2>Correlation is descriptive</h2><p>",cor_text," These correlations use all high-testing paired weeks, including eventual test observations, and are never used to choose models. Consecutive-week change correlations do not bridge gaps. A high correlation does not establish low prediction error.</p>",ca_table(ct),
    "<h2>Testing histories in the three original counties</h2>",paste(focus,collapse=''),
    "<h2>Untouched-test visual checks</h2><p>Every evaluable series at the primary volume threshold is available below, including shared-catchment sensitivity results. Points avoid drawing lines through excluded weeks.</p>",paste(plots,collapse=''),
    "<h2>Source quality and exclusions</h2>",ca_table(disp,0),
    "<p>Only broad SARS-CoV-2 N/N1/N2 targets are used, consistent with the earlier California analysis. Variant-specific and other pathogen targets are excluded. Treatment-plant samples need laboratory, matrix, units and method metadata. Explicit quality/analysis flags are propagated across duplicate samples. Identical duplicates are counted once; conflicting sample identifiers are excluded because revision order is not established. The qc_ignore field controls dashboard display and is not a quality failure.</p>",
    "<p>Explicit nondetects use LOD/√2 only when a valid detection limit is reported. Missing or nonpositive unflagged RNA is not imputed. Strictly grouped US thousands separators are parsed while original strings remain in the audit. A missing plant-county field is resolved only from a single explicit served-county FIPS in that same source row. Missing or conflicting geography remains outside the primary cohort. Laboratory-method ID changes remain separate segments even if they may represent metadata changes rather than assay changes.</p>",
    "<details><summary>Plant, laboratory, method and catchment registry</summary>",ca_table(reg),"</details>",
    "<h2>Interpretation and limits</h2><p>Earlier RNA data increase the overlap with high testing. They do not eliminate changes in test-seeking, screening, home testing, repeat testing, variants, shedding or sewershed coverage. The clinical and wastewater files are revised historical snapshots. Unlike WVAL, raw RNA does not use a future reference baseline, but these files still cannot reconstruct exactly what was available at each historical forecast origin.</p>",
    "<p>A model passing this short test is a candidate for additional validation across independent waves and sites. Applying a 2020–2022 calibration to later years would estimate a hypothetical positivity under that testing regime, not measured later positivity or prevalence. This script therefore produces no post-testing-collapse positivity reconstructions.</p>",
    "<h2>Sources and reproducibility</h2><ul><li><a href='",CA_LANDING,"'>California historical Cal-SuWers wastewater archive</a> and <a href='",DICTIONARY_URL,"'>NWSS metadata dictionary</a>.</li><li><a href='",CLINICAL_LANDING,"'>CDPH archived county PCR counts</a>; <a href='",CLINICAL_SCHEMA,"'>field definitions</a>. Clinical data ceased routine updates in December 2023. Counts are assigned by specimen collection date and county of residence; unknown dates and jurisdictions are audited separately.</li><li><a href='https://wwwnc.cdc.gov/eid/article/29/8/22-1634_article'>McManus et al., 2023</a>: wastewater forecasting in Denmark explicitly considered testing changes and variants; its target was incidence, not positivity.</li><li><a href='https://wwwnc.cdc.gov/eid/article/30/10/24-0637_article'>Lininger et al., 2024</a>: Oregon wastewater/clinical associations support surveillance value, while association alone does not establish predictive calibration.</li></ul>",
    "<p>The script exports sample/daily/weekly audits, clinical counts, all eligibility decisions, model partitions, tuning scores, held-out predictions, model objects, settings and source checksums. Raw downloads remain in the reusable cache. Each run has a new folder; the HTML embeds all figures and works independently.</p>",ca_table(sources),"</main></body></html>")
  writeLines(enc2utf8(h),file.path(run,'index.html'),useBytes=TRUE);writeLines(conclusion,file.path(run,'conclusion.txt'));write_csv(overview,file.path(run,'evaluation_summary.csv'))
  invisible(list(conclusion=conclusion,overview=overview))
}

run_california_wastewater_positivity_prediction <- function() {
  need_packages()
  if(!PRIMARY_REFERENCE_FRACTION%in%SENSITIVITY_FRACTIONS||MIN_TRAIN_WEEKS<16||TUNING_WEEKS<8||TEST_WEEKS<8||EMBARGO_DAYS<14||BOOTSTRAP_REPLICATES<99||MIN_SAMPLING_DAYS<1||MIN_SAMPLING_DAYS>7)fail('Invalid analysis settings.')
  set.seed(RANDOM_SEED);dir.create(OUTPUT_ROOT,recursive=TRUE,showWarnings=FALSE);cache<-file.path(OUTPUT_ROOT,'cache');dir.create(cache,showWarnings=FALSE)
  run<-tempfile(paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S',tz='UTC'),'_'),tmpdir=OUTPUT_ROOT);dir.create(run)
  for(n in c('data','plots','source_metadata'))dir.create(file.path(run,n))
  complete<-FALSE;on.exit({capture.output(sessionInfo(),file=file.path(run,'session_info.txt'));if(!complete)writeLines('INCOMPLETE: inspect the R error.',file.path(run,'RUN_FAILED.txt'))},add=TRUE)
  config<-ls(envir=environment(run_california_wastewater_positivity_prediction),pattern='^[A-Z][A-Z_0-9]+$');dput(mget(config,envir=environment(run_california_wastewater_positivity_prediction)),file.path(run,'settings.R'))
  sources<-list()
  acquire<-function(label,url,name,local='',validate=function(p)TRUE){a<-cached_source(url,file.path(cache,name),REFRESH_DOWNLOADS,local,validate);sources[[label]]<<-data.frame(source=label,url=url,path=a$path,acquired_utc=if(is.null(a$acquired_utc))NA_character_ else a$acquired_utc,inspected_utc=utc_now(),md5=a$md5,bytes=a$bytes,mode=a$mode);a$path}
  ww_file<-acquire('Cal-SuWers raw samples',WW_URL,'california_wastewater.csv',LOCAL_WASTEWATER_CSV,csv_validator(WW_COLUMNS))
  clinical_file<-acquire('CDPH county PCR counts',CLINICAL_URL,'clinical_daily.csv',LOCAL_CLINICAL_CSV,csv_validator(c('date','area','area_type','population','total_tests','positive_tests')))
  for(a in list(c('Wastewater metadata',CA_API,'wastewater_metadata.json'),c('Clinical metadata',CLINICAL_METADATA,'clinical_metadata.json'),c('Clinical field definitions',CLINICAL_SCHEMA,'clinical_schema.json'),c('NWSS metadata dictionary',DICTIONARY_URL,'nwss_dictionary.csv'))){p<-acquire(a[1],a[2],a[3],validate=if(grepl('json$',a[3]))json_validator else function(p)TRUE);file.copy(p,file.path(run,'source_metadata',a[3]))}
  message('Preparing matched PCR positivity and testing-volume histories...');clinical<-ca_clinical(clinical_file,run)
  message('Auditing California wastewater samples...');raw<-read_wastewater(ww_file);unassigned<-raw[is.na(raw$county_fips),];write_csv(unassigned,file.path(run,'data','unassigned_wastewater_rows.csv'))
  prepared<-prepare_samples(raw[!is.na(raw$county_fips),]);reg<-prepared$registry
  reg$geography_class<-ifelse(grepl('^source lists one county',reg$geography_note),'single county',ifelse(grepl('^multi-county',reg$geography_note),'shared catchment','unknown catchment'))
  prepared$registry<-reg
  write_csv(prepared$audit,file.path(run,'data','sample_audit.csv'));write_csv(reg,file.path(run,'data','series_registry.csv'))
  message('Aggregating dates into matched weekly windows...');weekly<-ca_weekly(prepared$samples);write_csv(weekly$daily,file.path(run,'data','wastewater_daily.csv'));write_csv(weekly$weekly,file.path(run,'data','wastewater_weekly.csv'))
  all_results<-all_cor<-list();coverage<-data.frame(fips=names(COUNTIES),county=unname(COUNTIES),source_rows=0L,eligible_samples=0L,method_series=0L,single_county_series=0L)
  for(i in seq_len(nrow(coverage))){f<-coverage$fips[i];coverage$source_rows[i]<-sum(raw$county_fips==f,na.rm=TRUE);coverage$eligible_samples[i]<-sum(prepared$samples$county_fips==f);coverage$method_series[i]<-sum(reg$county_fips==f);coverage$single_county_series[i]<-sum(reg$county_fips==f&reg$geography_class=='single county')}
  write_csv(coverage,file.path(run,'data','county_coverage.csv'))
  scenarios<-data.frame(fraction=SENSITIVITY_FRACTIONS,variant='raw');scenarios<-rbind(scenarios,data.frame(fraction=PRIMARY_REFERENCE_FRACTION,variant=c('raw_matched','pmmov','detected')))
  for(i in seq_len(nrow(scenarios))){frac<-scenarios$fraction[i];variant<-scenarios$variant[i];q<-ca_qualify(clinical,frac)
    if(variant=='raw')write_csv(q$months,file.path(run,'data',paste0('monthly_eligibility_',frac,'.csv')))
    for(lag in 0:1){p<-ca_pairs(weekly$weekly,q$weekly,reg,lag,variant);tag<-paste(frac,variant,lag,sep='_');message('Evaluating ',tag)
      write_csv(p,file.path(run,'data',paste0('pairs_',tag,'.csv')));cr<-ca_correlations(p);cr$fraction<-frac;all_cor[[tag]]<-cr;all_results[[tag]]<-ca_evaluate(p,frac)}
  }
  bind<-function(n){x<-lapply(all_results,`[[`,n);x<-x[vapply(x,nrow,integer(1))>0];if(length(x))do.call(rbind,x)else data.frame()}
  results<-setNames(lapply(c('summary','predictions','scores','audit','partitions'),bind),c('summary','predictions','scores','audit','partitions'))
  correlations<-do.call(rbind,all_cor);for(n in names(results))write_csv(results[[n]],file.path(run,'data',paste0(n,'.csv')));write_csv(correlations,file.path(run,'data','correlations.csv'))
  saveRDS(lapply(all_results,`[[`,'models'),file.path(run,'fitted_models.rds'));sources<-do.call(rbind,sources);write_csv(sources,file.path(run,'sources.csv'))
  report<-ca_report(run,clinical,prepared,results,correlations,sources,coverage)
  complete<-TRUE;writeLines(paste('COMPLETE',utc_now()),file.path(run,'RUN_COMPLETE.txt'));message(report$conclusion);message('Report: ',normalizePath(file.path(run,'index.html'),winslash='/'))
  if(OPEN_HTML)utils::browseURL(file.path(run,'index.html'))
  invisible(list(directory=normalizePath(run,winslash='/'),clinical=clinical,registry=reg,results=results,correlations=correlations,coverage=coverage))
}
if(isTRUE(RUN_ANALYSIS))run_california_wastewater_positivity_prediction()
