# County wastewater -> clinical positivity: high-testing evaluation
# Version 1.0 | 23 October 2026
# Run after us_wastewater_positivity_map.R:
# install.packages(c("curl", "jsonlite"))  # once, in a clean R session
# source("us_wastewater_positivity_prediction.R")
#
# The target is reported county NAAT positivity under a high-testing regime,
# NOT infection prevalence. High test volume does not establish random sampling.
# HHS counts are an eligibility indicator, not the exact denominator of the
# revised CDC positivity percentages. Never reconstruct positivity from these
# two sources or use HHS counts as binomial weights for the CDC percentages.
# WVAL is numeric, not a category. Current WVAL history was retrospectively
# recalibrated by CDC; this evaluates a revised-data relationship, not a live
# forecasting system with historical data-release availability.

RUN_ANALYSIS <- getOption("wwpredict.autorun", TRUE)
MAP_RUN_DIR <- ""  # empty: newest complete run under MAP_OUTPUT_ROOT
MAP_OUTPUT_ROOT <- file.path("analysis", "county_wastewater_map")
OUTPUT_ROOT <- file.path("analysis", "wastewater_positivity_prediction")
LOCAL_HHS_ARCHIVE_DIR <- "" # optional directory containing original di4u-7yu6_*.csv files
LOCAL_HHS_INDEX <- "" # optional original HHS archive-index JSON
REFRESH_DOWNLOADS <- FALSE
OPEN_HTML <- interactive()

REFERENCE_START <- as.Date("2021-07-01")
REFERENCE_END <- as.Date("2021-12-31")
ANALYSIS_START <- as.Date("2022-01-01")
ANALYSIS_END <- as.Date("2023-05-10")
MIN_MONTH_REPORTS <- 3L
MIN_MONTH_SPAN_DAYS <- 14L
MIN_REFERENCE_MONTHS <- 4L
MIN_WEEKLY_TESTS <- 500
MIN_WEEKLY_TESTS_PER_100K <- 1000
PRIMARY_REFERENCE_FRACTION <- 0.75
SENSITIVITY_FRACTIONS <- c(0.50, 0.75, 1.00)
MAX_TESTING_REPORT_AGE_DAYS <- 7L

EXCLUDE_SHARED_CATCHMENTS <- TRUE
PANEL_SELECTION_WEEKS <- 12L # choose modal site set in first 12 source calendar weeks
MIN_PANEL_SELECTION_OBSERVATIONS <- 6L
MIN_TRAIN_WEEKS <- 16L
TUNING_WEEKS <- 8L
TEST_WEEKS <- 8L
EMBARGO_DAYS <- 14L # two unused calendar weeks between partitions
MIN_CORRELATION_WEEKS <- 12L
BOOTSTRAP_REPLICATES <- 499L
BLOCK_WEEKS <- 3L
RANDOM_SEED <- 24023L
# Explicit exploratory performance criteria, not clinical validation thresholds.
MAX_ACCEPTABLE_MAE_PP <- 3
MIN_MAE_SKILL <- 0.10
MAX_ACCEPTABLE_ABS_BIAS_PP <- 2
FOCUS_FIPS <- c("06075", "06085", "06067")

HHS_ARCHIVE_API <- "https://healthdata.gov/resource/pcq6-asfs.json"
HHS_METADATA_URL <- "https://healthdata.gov/api/views/pcq6-asfs.json"
HHS_VOLUME_METADATA_URL <- "https://healthdata.gov/api/views/di4u-7yu6.json"

wp_stop <- function(...) stop(...,call.=FALSE)
wp_csv <- function(x,p) write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
wp_read <- function(p) read.csv(p,colClasses="character",check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("","NA"),fileEncoding="UTF-8-BOM")
wp_need <- function(x,n) if(!all(n%in%names(x)))wp_stop("Missing columns: ",paste(setdiff(n,names(x)),collapse=", "))
wp_num <- function(x) {
  s<-trimws(as.character(x));v<-suppressWarnings(as.numeric(s));missing<-is.na(s)|s==""
  if(any(!missing & !is.finite(v)))wp_stop("Invalid numeric input: ",paste(head(unique(s[!missing&!is.finite(v)])),collapse=", "))
  v[missing]<-NA_real_;v
}
wp_date <- function(x) {
  s<-substr(x,1,10);a<-as.Date(s,format="%Y-%m-%d");bad<-is.na(a)
  a[bad]<-as.Date(s[bad],format="%m/%d/%Y")
  if(anyNA(a))wp_stop("Invalid dates.");a
}
wp_fips <- function(x) {
  if(any(!grepl("^[0-9]{1,5}$",x)))wp_stop("Invalid county FIPS.")
  sprintf("%05d",as.integer(x))
}
wp_month <- function(x) format(x,"%Y-%m")
wp_md5 <- function(p) unname(tools::md5sum(p))
wp_utc <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
wp_query <- function(u,p) paste0(u,"?",paste(vapply(seq_along(p),function(i)paste0(utils::URLencode(names(p)[i],reserved=TRUE),"=",utils::URLencode(p[[i]],reserved=TRUE)),""),collapse="&"))
wp_fetch <- function(url,dest,refresh=REFRESH_DOWNLOADS) {
  if(file.exists(dest)&&!refresh)return(dest)
  tmp<-tempfile(tmpdir=dirname(dest));on.exit(unlink(tmp),add=TRUE)
  for(k in 1:3) {
    h<-curl::new_handle(timeout=180,connecttimeout=25,followlocation=TRUE,useragent="County-positivity-evaluation/1.0")
    z<-tryCatch(curl::curl_fetch_disk(url,tmp,h),error=identity)
    code<-if(inherits(z,"error"))NA_integer_ else z$status_code
    if(!is.na(code)&&code>=200&&code<300&&file.info(tmp)$size>0) {
      first<-readLines(tmp,n=1,warn=FALSE)
      if(grepl("^\\s*<(!doctype|html)",first,ignore.case=TRUE))wp_stop("HTML returned instead of source data: ",url)
      if(!file.copy(tmp,dest,overwrite=TRUE))wp_stop("Cannot save download.")
      return(dest)
    }
    if(k==3L||(!is.na(code)&&code<500&&!code%in%c(408,429)))wp_stop("Download failed: ",if(inherits(z,"error"))conditionMessage(z)else code,"\n",url)
    Sys.sleep(2^k)
  }
}
wp_inputs <- function(run) {
  root<-MAP_RUN_DIR
  required<-c("wastewater_county_weekly.csv","positivity_county_weekly.csv","boundary_counties.csv")
  if(!nzchar(root)) {
    dirs<-list.dirs(MAP_OUTPUT_ROOT,recursive=FALSE,full.names=TRUE)
    dirs<-dirs[vapply(dirs,function(d)file.exists(file.path(d,"RUN_COMPLETE.txt"))&&all(file.exists(file.path(d,"data",required))),logical(1))]
    if(!length(dirs))wp_stop("Run us_wastewater_positivity_map.R first, or set MAP_RUN_DIR to its completed output folder.")
    root<-sort(dirs,decreasing=TRUE)[1]
  }
  paths<-file.path(root,"data",required)
  if(any(!file.exists(paths)))wp_stop("Incomplete input map folder: ",root)
  w<-wp_read(paths[1]);p<-wp_read(paths[2]);a<-wp_read(paths[3])
  wp_need(w,c("fips","week_end","wval","sites","shared_sites","site_ids"));wp_need(p,c("fips","week_end","positivity"));wp_need(a,c("fips","name","state"))
  w$fips<-wp_fips(w$fips);p$fips<-wp_fips(p$fips);a$fips<-wp_fips(a$fips)
  w$week_end<-wp_date(w$week_end);p$week_end<-wp_date(p$week_end)
  for(n in c("wval","sites","shared_sites"))w[[n]]<-wp_num(w[[n]])
  p$positivity<-wp_num(p$positivity)
  if(anyDuplicated(w[c("fips","week_end")])||anyDuplicated(p[c("fips","week_end")]))wp_stop("Duplicate county/week inputs.")
  if(any(as.POSIXlt(w$week_end)$wday!=6)||any(as.POSIXlt(p$week_end)$wday!=3))wp_stop("Unexpected source reporting weeks.")
  if(any(w$wval<0,na.rm=TRUE)||any(p$positivity<0|p$positivity>100,na.rm=TRUE))wp_stop("Out-of-range source values.")
  if(any(!w$fips%in%a$fips)||any(!p$fips%in%a$fips))wp_stop("Unknown input county FIPS.")
  wp_csv(data.frame(file=normalizePath(paths,winslash="/"),md5=unname(tools::md5sum(paths))),file.path(run,"sources","map_inputs.csv"))
  w<-w[w$week_end>=ANALYSIS_START-28&w$week_end<=ANALYSIS_END+3,]
  p<-p[p$week_end>=ANALYSIS_START-28&p$week_end<=ANALYSIS_END,]
  list(w=w,p=p,areas=a,root=root)
}

wp_testing <- function(cache,run,areas) {
  index_url<-wp_query(HHS_ARCHIVE_API,list(`$select`="update_date,archive_link,rows",`$where`="update_date >= '2021-06-01T00:00:00' and update_date < '2023-06-01T00:00:00'",`$order`="update_date",`$limit`="5000"))
  index_file<-if(nzchar(LOCAL_HHS_INDEX))LOCAL_HHS_INDEX else wp_fetch(index_url,file.path(cache,"hhs_archive_index.json"))
  index<-jsonlite::fromJSON(index_file,simplifyVector=FALSE)
  if(!length(index)||length(index)>=5000)wp_stop("Unexpected HHS archive index; pagination may be required.")
  wp_csv(data.frame(source="HHS archive index",url=index_url,md5=wp_md5(index_file)),file.path(run,"sources","index_manifest.csv"))
  file.copy(index_file,file.path(run,"sources","hhs_archive_index.json"))
  for(item in list(c(HHS_METADATA_URL,"hhs_archive_metadata.json"),c(HHS_VOLUME_METADATA_URL,"hhs_volume_metadata.json"))) {
    q<-wp_fetch(item[1],file.path(cache,item[2]));file.copy(q,file.path(run,"sources",item[2]))
  }
  dest<-file.path(cache,"hhs_archive");dir.create(dest,showWarnings=FALSE)
  required<-c("fips","county","state","date","total_tests_last_7_days","total_tests_per_100k_last_7_days")
  pieces<-manifest<-missing_dates<-vector("list",length(index))
  for(i in seq_along(index)) {
    z<-index[[i]];u<-z$archive_link$url
    if(length(u)!=1||!grepl("^https://us-dhhs-aa[.]s3[.]us-east-2[.]amazonaws[.]com/di4u-7yu6_[^/]+[.]csv$",u))wp_stop("Unexpected archive URL; inspect index.")
    name<-basename(u);q<-if(nzchar(LOCAL_HHS_ARCHIVE_DIR))file.path(LOCAL_HHS_ARCHIVE_DIR,name)else wp_fetch(u,file.path(dest,name))
    if(!file.exists(q))wp_stop("Missing original HHS snapshot: ",q)
    x<-wp_read(q);wp_need(x,required)
    if(nrow(x)!=as.integer(z$rows))wp_stop("HHS archive row-count mismatch: ",name)
    x<-x[required]
    missing_date<-is.na(x$date)|!nzchar(trimws(x$date))
    missing_dates[[i]]<-data.frame(archive_file=name,rows_missing_date=sum(missing_date))
    x<-x[!missing_date,,drop=FALSE];x$fips<-wp_fips(x$fips);x$date<-wp_date(x$date)
    x<-x[x$fips%in%areas$fips,];x$archive_update<-rep(z$update_date,nrow(x));x$archive_file<-rep(name,nrow(x))
    pieces[[i]]<-x;manifest[[i]]<-data.frame(file=name,url=u,md5=wp_md5(q),rows=as.integer(z$rows),update=z$update_date)
    if(i==1L||i%%20L==0L||i==length(index))message("Testing archive ",i," / ",length(index))
  }
  q<-do.call(rbind,pieces);m<-do.call(rbind,manifest);wp_csv(m,file.path(run,"sources","hhs_snapshot_manifest.csv"))
  wp_csv(do.call(rbind,missing_dates),file.path(run,"data","testing_missing_dates_audit.csv"))
  q$tests<-wp_num(q$total_tests_last_7_days);q$rate<-wp_num(q$total_tests_per_100k_last_7_days)
  bad<-(!is.na(q$tests)&q$tests<0)|(!is.na(q$rate)&q$rate<0)
  wp_csv(q[bad,],file.path(run,"data","invalid_testing_rows.csv"));q$tests[bad]<-NA;q$rate[bad]<-NA
  # Latest archived revision for each county/report date; retain revision audit.
  q<-q[order(q$archive_update,decreasing=TRUE),];key<-paste(q$fips,q$date)
  duplicate<-duplicated(key)|duplicated(key,fromLast=TRUE)
  wp_csv(q[duplicate,],file.path(run,"data","testing_repeated_report_dates.csv"))
  q<-q[!duplicated(key),];q<-q[order(q$fips,q$date),];q$month<-wp_month(q$date)
  wp_csv(q,file.path(run,"data","testing_report_snapshots.csv"))
  groups<-split(seq_len(nrow(q)),paste(q$fips,q$month))
  message('Summarising testing volume by county/month')
  values<-t(vapply(groups,function(ii) {
    ok<-ii[is.finite(q$rate[ii])&is.finite(q$tests[ii])]
    c(reports=length(ok),span_days=if(length(ok))as.integer(max(q$date[ok])-min(q$date[ok]))else 0L,
      median_weekly_tests=if(length(ok))median(q$tests[ok])else NA_real_,median_weekly_tests_per_100k=if(length(ok))median(q$rate[ok])else NA_real_)
  },numeric(4)))
  first<-vapply(groups,`[`,integer(1),1L)
  monthly<-data.frame(fips=q$fips[first],month=q$month[first],values,row.names=NULL)
  monthly$adequate_reports<-monthly$reports>=MIN_MONTH_REPORTS&monthly$span_days>=MIN_MONTH_SPAN_DAYS
  ref<-monthly[monthly$month>=wp_month(REFERENCE_START)&monthly$month<=wp_month(REFERENCE_END)&monthly$adequate_reports,]
  ref<-do.call(rbind,lapply(split(ref,ref$fips),function(x)data.frame(fips=x$fips[1],reference_months=nrow(x),reference_weekly_rate=median(x$median_weekly_tests_per_100k))))
  monthly<-merge(monthly,ref,by="fips",all.x=TRUE,sort=FALSE)
  monthly$reference_ok<-!is.na(monthly$reference_months)&monthly$reference_months>=MIN_REFERENCE_MONTHS&is.finite(monthly$reference_weekly_rate)&monthly$reference_weekly_rate>0
  monthly$relative_volume<-monthly$median_weekly_tests_per_100k/monthly$reference_weekly_rate
  wp_csv(monthly,file.path(run,"data","monthly_testing.csv"))
  list(reports=q,months=monthly)
}
wp_qualify <- function(testing,p,fraction) {
  m<-testing$months
  m$high_testing<-m$adequate_reports&m$reference_ok&is.finite(m$relative_volume)&m$relative_volume>=fraction&
    m$median_weekly_tests>=MIN_WEEKLY_TESTS&m$median_weekly_tests_per_100k>=MIN_WEEKLY_TESTS_PER_100K
  m$high_testing[is.na(m$high_testing)]<-FALSE
  z<-p;z$month<-wp_month(z$week_end);hit<-match(paste(z$fips,z$month),paste(m$fips,m$month))
  z$month_high<-m$high_testing[hit];z$month_high[is.na(z$month_high)]<-FALSE
  z$volume_report_date<-as.Date(NA);z$weekly_tests<-z$weekly_rate<-NA_real_
  pg<-split(seq_len(nrow(z)),z$fips);rg<-split(seq_len(nrow(testing$reports)),testing$reports$fips);r<-testing$reports
  for(f in intersect(names(pg),names(rg))) {
    ii<-pg[[f]];rr<-rg[[f]];rr<-rr[order(r$date[rr])]
    j<-findInterval(as.numeric(z$week_end[ii]),as.numeric(r$date[rr]));ok<-j>0
    z$volume_report_date[ii[ok]]<-r$date[rr[j[ok]]];z$weekly_tests[ii[ok]]<-r$tests[rr[j[ok]]];z$weekly_rate[ii[ok]]<-r$rate[rr[j[ok]]]
  }
  z$report_age_days<-as.integer(z$week_end-z$volume_report_date)
  z$high_testing<-z$month_high&!is.na(z$report_age_days)&z$report_age_days<=MAX_TESTING_REPORT_AGE_DAYS&
    is.finite(z$weekly_tests)&z$weekly_tests>=MIN_WEEKLY_TESTS&is.finite(z$weekly_rate)&z$weekly_rate>=MIN_WEEKLY_TESTS_PER_100K
  z$high_testing[is.na(z$high_testing)]<-FALSE
  list(p=z,m=m)
}
wp_panels <- function(w) {
  rows<-lapply(split(w,w$fips),function(x) {
    x<-x[order(x$week_end),];first<-min(x$week_end)
    z<-x[x$week_end<=first+7L*(PANEL_SELECTION_WEEKS-1L)&is.finite(x$wval)&!is.na(x$site_ids)&nzchar(x$site_ids),]
    if(EXCLUDE_SHARED_CATCHMENTS)z<-z[z$shared_sites==0,]
    if(nrow(z)<MIN_PANEL_SELECTION_OBSERVATIONS)return(NULL)
    tab<-sort(table(z$site_ids),decreasing=TRUE);chosen<-names(tab)[1]
    data.frame(fips=x$fips[1],site_ids=chosen,panel_selection_end=first+7L*(PANEL_SELECTION_WEEKS-1L),selection_observations=nrow(z),modal_observations=as.integer(tab[1]))
  });z<-do.call(rbind,rows);if(is.null(z))wp_stop("No usable fixed site panels.");rownames(z)<-NULL;z
}
wp_pairs <- function(w,p,panels,offset) {
  z<-p[p$week_end>=ANALYSIS_START&p$week_end<=ANALYSIS_END,]
  z$target_end<-z$week_end;z$target_start<-z$target_end-6L;z$wastewater_end<-z$target_end-offset;z$wastewater_start<-z$wastewater_end-6L
  wk<-paste(w$fips,w$week_end);j<-match(paste(z$fips,z$wastewater_end),wk);j1<-match(paste(z$fips,z$wastewater_end-7L),wk)
  panel<-match(z$fips,panels$fips);z$panel_selection_end<-panels$panel_selection_end[panel]
  z$wval<-w$wval[j];z$previous_wval<-w$wval[j1];z$site_ids<-w$site_ids[j]
  z$stable_panel<-!is.na(j)&!is.na(j1)&!is.na(panel)&w$site_ids[j]==panels$site_ids[panel]&w$site_ids[j1]==panels$site_ids[panel]
  if(EXCLUDE_SHARED_CATCHMENTS)z$stable_panel<-z$stable_panel&w$shared_sites[j]==0&w$shared_sites[j1]==0
  z$stable_panel[is.na(z$stable_panel)]<-FALSE
  # In the lead task, last clinical value available by wastewater collection end
  # is target-14 days; target-7 would already be in the future.
  prior_gap<-if(offset>0)14L else 7L
  k<-match(paste(z$fips,z$target_end-prior_gap),paste(p$fips,p$week_end));z$prior_positivity<-p$positivity[k];z$prior_high_testing<-p$high_testing[k]
  z$prior_high_testing[is.na(z$prior_high_testing)]<-FALSE
  z$pair_eligible<-z$high_testing&z$stable_panel&is.finite(z$wval)&is.finite(z$positivity)
  z$model_eligible<-z$pair_eligible&is.finite(z$previous_wval)&is.finite(z$prior_positivity)&z$prior_high_testing
  z$logw<-log1p(z$wval);z$change_logw<-log1p(z$wval)-log1p(z$previous_wval)
  logit<-function(v)qlogis(pmin(.9995,pmax(.0005,v/100)))
  z$target_logit<-logit(z$positivity);z$prior_logit<-logit(z$prior_positivity)
  z$task<-if(offset<0)"contemporaneous_estimation"else"11_day_sample_lead"
  z
}
wp_cor <- function(x,y,method="pearson") if(length(x)>=3L&&sd(x)>0&&sd(y)>0)unname(cor(x,y,method=method))else NA_real_
wp_correlations <- function(z) {
  do.call(rbind,lapply(split(z[z$pair_eligible,],z$fips[z$pair_eligible]),function(d) {
    d<-d[order(d$target_end),];prior<-match(d$target_end-7L,d$target_end);ok<-!is.na(prior)
    data.frame(fips=d$fips[1],weeks=nrow(d),first=min(d$target_end),last=max(d$target_end),
      pearson_logw=if(nrow(d)>=MIN_CORRELATION_WEEKS)wp_cor(d$logw,d$positivity)else NA_real_,
      spearman=if(nrow(d)>=MIN_CORRELATION_WEEKS)wp_cor(d$wval,d$positivity,"spearman")else NA_real_,
      consecutive_changes=sum(ok),change_correlation=if(sum(ok)>=MIN_CORRELATION_WEEKS)wp_cor(d$logw[ok]-d$logw[prior[ok]],d$positivity[ok]-d$positivity[prior[ok]])else NA_real_)
  }))
}
wp_partition <- function(d) {
  d<-d[order(d$target_end),];n<-nrow(d)
  if(n<MIN_TRAIN_WEEKS+TUNING_WEEKS+TEST_WEEKS)return(NULL)
  test<-tail(seq_len(n),TEST_WEEKS)
  before_test<-which(d$target_end<min(d$target_end[test])-EMBARGO_DAYS)
  if(length(before_test)<TUNING_WEEKS)return(NULL)
  tune<-tail(before_test,TUNING_WEEKS)
  train<-which(d$target_end<min(d$target_end[tune])-EMBARGO_DAYS)
  if(length(train)<MIN_TRAIN_WEEKS)return(NULL)
  if(any(d$panel_selection_end[train]>max(d$wastewater_end[train])))return(NULL)
  list(train=train,tune=tune,test=test)
}
wp_models <- list(WVAL_log=target_logit~logw,WVAL_spline=target_logit~splines::ns(logw,df=3),
  WVAL_change=target_logit~logw+change_logw,AR_baseline=target_logit~prior_logit,WVAL_plus_prior=target_logit~logw+prior_logit)
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
wp_evaluate <- function(z,fraction) {
  summaries<-predictions<-scores<-audits<-models<-list();k<-0L
  groups<-split(z,z$fips)
  for(f in names(groups)) {
    county<-groups[[f]];d<-county[county$model_eligible,];d<-d[order(d$target_end),];partition<-wp_partition(d)
    status<-if(!any(county$stable_panel))'no eligible fixed site panel'else if(!any(county$pair_eligible))'no high-testing overlap'else 'insufficient high-testing weeks'
    audit<-data.frame(fips=f,fraction=fraction,task=z$task[1],source_high_testing_weeks=sum(county$high_testing),
      high_testing_with_wastewater=sum(county$high_testing&is.finite(county$wval)&is.finite(county$positivity)),
      high_testing_pairs=sum(county$pair_eligible),common_model_weeks=nrow(d),
      status=status,train_weeks=0L,tuning_weeks=0L,test_weeks=0L)
    if(is.null(partition)){audits[[f]]<-audit;next}
    tr<-d[partition$train,];va<-d[partition$tune,];te<-d[partition$test,]
    names_to_fit<-c(names(wp_models),"training_mean","persistence")
    fits<-setNames(lapply(names_to_fit,function(n)wp_fit(n,tr)),names_to_fit)
    validation<-sapply(fits,function(fit){hat<-wp_predict(fit,va);if(all(is.finite(hat)))mean(abs(hat-va$positivity))else Inf})
    ww_names<-c("WVAL_log","WVAL_spline","WVAL_change");winner<-ww_names[which.min(validation[ww_names])]
    if(!is.finite(validation[winner])){audit$status<-"model estimation failure";audits[[f]]<-audit;next}
    # Model selection ends here. The final test values have not been used.
    development<-rbind(tr,va);refits<-setNames(lapply(names_to_fit,function(n)wp_fit(n,development)),names_to_fit)
    hats<-sapply(refits,wp_predict,d=te)
    if(any(!is.finite(hats[,winner]))){audit$status<-"test prediction failure";audits[[f]]<-audit;next}
    k<-k+1L;hat<-hats[,winner];metrics<-wp_metrics(te$positivity,hat,te$prior_positivity);interval<-wp_bootstrap(te,hat)
    supported<-is.finite(metrics['mae_skill'])&&metrics['mae']<=MAX_ACCEPTABLE_MAE_PP&&metrics['mae_skill']>=MIN_MAE_SKILL&&abs(metrics['bias'])<=MAX_ACCEPTABLE_ABS_BIAS_PP
    precision_supported<-supported&&is.finite(interval['skill_lo'])&&interval['skill_lo']>0
    summaries[[k]]<-data.frame(fips=f,fraction=fraction,task=z$task[1],winner=winner,train_weeks=nrow(tr),tuning_weeks=nrow(va),test_weeks=nrow(te),
      train_first=min(tr$target_end),train_last=max(tr$target_end),tune_first=min(va$target_end),tune_last=max(va$target_end),test_first=min(te$target_end),test_last=max(te$target_end),
      as.list(metrics),as.list(interval),passes_point_criteria=supported,positive_skill_interval=precision_supported,
      extrapolation_weeks=sum(te$logw<min(development$logw)|te$logw>max(development$logw)),row.names=NULL)
    for(n in names_to_fit)if(all(is.finite(hats[,n])))scores[[paste(f,n)]]<-data.frame(fips=f,fraction=fraction,task=z$task[1],model=n,selected_ww=n==winner,
      validation_mae=validation[n],as.list(wp_metrics(te$positivity,hats[,n],te$prior_positivity)),row.names=NULL)
    out<-te[c("fips","target_start","target_end","wastewater_start","wastewater_end","wval","positivity","prior_positivity","site_ids","weekly_tests","weekly_rate")]
    out$fraction<-fraction;out$task<-z$task[1];out$model<-winner;out$prediction<-hat;out$error_pp<-hat-out$positivity
    out$within_development_wval_range<-te$logw>=min(development$logw)&te$logw<=max(development$logw);predictions[[k]]<-out
    audit$status<-"evaluated";audit$train_weeks<-nrow(tr);audit$tuning_weeks<-nrow(va);audit$test_weeks<-nrow(te);audits[[f]]<-audit
    models[[paste(f,fraction,z$task[1],sep="_")]]<-list(fit=refits[[winner]],range=range(development$wval),site_ids=unique(te$site_ids),test_metrics=metrics)
  }
  bind<-function(x)if(length(x))do.call(rbind,x)else data.frame()
  list(summary=bind(summaries),predictions=bind(predictions),scores=bind(scores),audit=bind(audits),models=models)
}

wp_escape <- function(x) {x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
wp_table <- function(d,digits=3) {
  if(!nrow(d))return("<p>No rows meet these requirements.</p>")
  for(n in names(d))if(is.numeric(d[[n]]))d[[n]]<-ifelse(is.finite(d[[n]]),format(round(d[[n]],digits),trim=TRUE),"—")
  rows<-apply(d,1,function(x)paste0("<tr><td>",paste(wp_escape(x),collapse="</td><td>"),"</td></tr>"))
  paste0("<div class='scroll'><table><thead><tr><th>",paste(wp_escape(names(d)),collapse="</th><th>"),"</th></tr></thead><tbody>",paste(rows,collapse=""),"</tbody></table></div>")
}
wp_image <- function(path) paste0("<img src='data:image/png;base64,",gsub("[\r\n]","",jsonlite::base64_enc(readBin(path,"raw",n=file.info(path)$size))),"' alt='Surveillance analysis plot'>")
wp_county_plot <- function(f,w,p,months,results,areas,path) {
  a<-areas[match(f,areas$fips),];clinical<-p[p$fips==f,];water<-w[w$fips==f,];m<-months[months$fips==f,]
  clinical<-clinical[order(clinical$week_end),];water<-water[order(water$week_end),];m<-m[order(m$month),]
  grDevices::png(path,1400,1200,res=130);on.exit(dev.off(),add=TRUE)
  par(mfrow=c(3,1),mar=c(4.5,5,3,1),las=1)
  dates<-as.Date(paste0(m$month,"-15"));threshold<-pmax(MIN_WEEKLY_TESTS_PER_100K,PRIMARY_REFERENCE_FRACTION*m$reference_weekly_rate)
  xrange<-range(c(dates,clinical$week_end,water$week_end),na.rm=TRUE)
  ylim<-c(0,max(c(m$median_weekly_tests_per_100k,threshold),na.rm=TRUE)*1.05)
  plot(dates,m$median_weekly_tests_per_100k,type='b',pch=16,col='#6a7f8a',ylim=ylim,xlim=xrange,xlab='',ylab='Weekly tests / 100,000',main=paste(a$name,a$state,'| monthly testing eligibility'))
  lines(dates,threshold,col='#b44835',lty=2);ok<-m$high_testing;points(dates[ok],m$median_weekly_tests_per_100k[ok],pch=16,col='#167b72',cex=1.3)
  legend('topright',c('Monthly median of report snapshots','Prespecified threshold','Qualifying month'),col=c('#6a7f8a','#b44835','#167b72'),lty=c(1,2,NA),pch=c(16,NA,16),bty='n',cex=.8)
  plot(clinical$week_end,clinical$positivity,type='l',col='#b3414c',xlim=xrange,xlab='',ylab='Positivity (%)',main='Observed clinical positivity; high-testing weeks highlighted')
  ok<-clinical$high_testing;points(clinical$week_end[ok],clinical$positivity[ok],pch=16,col='#167b72',cex=.7)
  pr<-results$predictions
  if(nrow(pr)){pr<-pr[pr$fips==f&pr$task=='contemporaneous_estimation'&pr$fraction==PRIMARY_REFERENCE_FRACTION,];if(nrow(pr))points(pr$target_end,pr$prediction,pch=1,col='#245eaa',cex=1.1)}
  if(nrow(pr))legend('topright',c('Observed','High-testing week','Held-out wastewater estimate'),col=c('#b3414c','#167b72','#245eaa'),lty=c(1,NA,NA),pch=c(NA,16,1),bty='n',cex=.8)else legend('topright',c('Observed','High-testing week'),col=c('#b3414c','#167b72'),lty=c(1,NA),pch=c(NA,16),bty='n',cex=.8)
  plot(water$week_end,log1p(water$wval),type='l',col='#236e90',xlim=xrange,xlab='Source week end',ylab='log(1 + numeric WVAL)',main='Wastewater county median | original observations, before eligibility filtering')
}
wp_report <- function(run,inputs,testing,primary_p,primary_months,correlations,results,settings) {
  s<-results$summary;au<-results$audit;c0<-correlations[correlations$fraction==PRIMARY_REFERENCE_FRACTION&correlations$task=='contemporaneous_estimation',]
  selected<-if(nrow(s))s[s$fraction==PRIMARY_REFERENCE_FRACTION&s$task=='contemporaneous_estimation',]else data.frame()
  enough_cor<-c0[is.finite(c0$pearson_logw),]
  conclusion<-if(!nrow(selected))paste0("No county has enough qualifying observations for the prespecified training, tuning and held-out test. A validated wastewater-to-positivity predictor cannot be established from this archive under the primary high-testing rule.")else
    paste0(nrow(selected)," counties support the full evaluation. Median held-out MAE is ",round(median(selected$mae),2)," percentage points; median skill against persistence is ",round(100*median(selected$mae_skill,na.rm=TRUE),1),"%. ",sum(selected$passes_point_criteria)," meet the prespecified point-performance criteria; ",sum(selected$positive_skill_interval)," also have a positive lower bootstrap skill bound. These are exploratory, local results, not validation for 2026 or population prevalence.")
  summary_rows<-do.call(rbind,lapply(split(au,paste(au$fraction,au$task)),function(x){
    y<-if(nrow(s))s[s$fraction==x$fraction[1]&s$task==x$task[1],]else data.frame()
    data.frame(reference_fraction=x$fraction[1],task=x$task[1],counties_reviewed=nrow(x),counties_evaluated=sum(x$status=='evaluated'),
      median_test_mae=if(nrow(y))median(y$mae)else NA,median_skill_pct=if(nrow(y))100*median(y$mae_skill,na.rm=TRUE)else NA,
      pass_point_criteria=if(nrow(y))sum(y$passes_point_criteria)else 0)
  }))
  selected_show<-selected
  if(nrow(selected_show)){selected_show$county<-inputs$areas$name[match(selected_show$fips,inputs$areas$fips)];selected_show$state<-inputs$areas$state[match(selected_show$fips,inputs$areas$fips)];selected_show<-selected_show[c('fips','county','state','winner','train_weeks','tuning_weeks','test_weeks','mae','rmse','bias','r2','mae_skill','skill_lo','skill_hi','passes_point_criteria')]}
  sensitivity_show<-if(nrow(s))s[s$fraction!=PRIMARY_REFERENCE_FRACTION,]else data.frame()
  if(nrow(sensitivity_show)){sensitivity_show$county<-inputs$areas$name[match(sensitivity_show$fips,inputs$areas$fips)];sensitivity_show$state<-inputs$areas$state[match(sensitivity_show$fips,inputs$areas$fips)];sensitivity_show<-sensitivity_show[c('fips','county','state','fraction','task','winner','mae','bias','mae_skill','skill_lo','skill_hi','passes_point_criteria')]}
  cshow<-merge(c0,inputs$areas[c('fips','name','state')],by='fips',all.x=TRUE);cshow<-cshow[order(-cshow$weeks,cshow$fips),]
  focus_audit<-au[au$fips%in%FOCUS_FIPS&au$fraction==PRIMARY_REFERENCE_FRACTION&au$task=='contemporaneous_estimation',c('fips','high_testing_with_wastewater','high_testing_pairs','common_model_weeks','status')]
  focus_audit<-merge(inputs$areas[inputs$areas$fips%in%FOCUS_FIPS,c('fips','name','state')],focus_audit,by='fips',all.x=TRUE)
  focus<-intersect(FOCUS_FIPS,unique(inputs$w$fips));images<-character()
  for(f in focus){q<-file.path(run,'plots',paste0('county_',f,'.png'));wp_county_plot(f,inputs$w,primary_p,primary_months,results,inputs$areas,q);images<-c(images,wp_image(q))}
  fmt<-function(x)if(length(x)&&is.finite(x))format(round(x,3),trim=TRUE)else 'not estimable'
  h<-paste0("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Wastewater as a predictor of positivity</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#18374a;background:#f4f7f9;margin:0}main{max-width:1300px;margin:auto;padding:32px}h1,h2{line-height:1.2}h1{font-size:32px}h2{margin-top:36px}.card{background:white;border:1px solid #dce4e8;border-radius:10px;padding:22px;margin:18px 0}.verdict{border-left:6px solid #27798c}.scroll{overflow:auto;max-height:650px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px;border-bottom:1px solid #dbe4e8;text-align:left;white-space:nowrap}th{position:sticky;top:0;background:#edf3f6}img{max-width:100%;height:auto}a{color:#126b91}code{background:#e9eff3;padding:2px 5px}.small{font-size:13px;color:#536b7a}li{margin-bottom:8px}</style><main>",
    "<p class='small'>US county surveillance · ",wp_utc(),"</p><h1>Can wastewater predict clinical test positivity?</h1>",
    "<div class='card verdict'><strong>Primary conclusion</strong><p>",conclusion,"</p></div>",
    "<p>The target is the percentage of reported county NAAT tests that are positive during a high-testing regime. It is not the percentage of residents infected. We use numeric WVALs throughout; CDC category labels and color thresholds are irrelevant to these models.</p>",
    "<h2>What the high-testing restriction leaves</h2><p>",nrow(enough_cor)," counties have at least ",MIN_CORRELATION_WEEKS," eligible paired weeks for a descriptive correlation. Across these counties, median Pearson correlation between log(1+WVAL) and positivity is ",fmt(median(enough_cor$pearson_logw,na.rm=TRUE)),"; median Spearman correlation is ",fmt(median(enough_cor$spearman,na.rm=TRUE)),". These correlations use eligible observations from the whole period and are not held-out prediction scores. The change correlations require consecutive source weeks and never bridge gaps.</p>",
    wp_table(summary_rows),
    "<h2>Eligibility fixed before examining the relationship</h2><div class='card'><ul>",
    "<li>For each county, the reference is the median of monthly median weekly testing rates during July–December 2021. At least ",MIN_REFERENCE_MONTHS," adequately observed reference months are required.</li>",
    "<li>A month needs at least ",MIN_MONTH_REPORTS," distinct HHS report dates spanning ",MIN_MONTH_SPAN_DAYS," days. Its median weekly testing must be at least ",MIN_WEEKLY_TESTS," tests and ",MIN_WEEKLY_TESTS_PER_100K," tests per 100,000 residents, and at least ",100*PRIMARY_REFERENCE_FRACTION,"% of its own reference rate.</li>",
    "<li>Each retained clinical week must also have a preceding HHS volume report no more than ",MAX_TESTING_REPORT_AGE_DAYS," days old meeting both absolute floors. Missing volume is ineligible. These are explicit analytic thresholds, not CDC rules or proof that testing is representative.</li>",
    "<li>Sensitivity analyses repeat the entire evaluation at 50%, 75%, and 100% of reference volume. They are reported separately; we do not select a threshold because it produces the highest correlation.</li>",
    "<li>HHS report snapshots describe overlapping seven-day measures. We take medians, never sum them into monthly test counts. The report-date schedule differs from the revised positivity archive. HHS volume is an eligibility indicator, not the matched denominator of the CDC percentage, so no binomial weighting is used.</li>",
    "<li>A fixed county site set is selected from the first ",PANEL_SELECTION_WEEKS," source calendar weeks using availability only. The modal set is retained; other site combinations are excluded. Weeks involving shared catchments are excluded by default. These rules reduce changing-site artifacts and duplicate catchment information.</li>",
    "<li>Duplicate county/report dates use the latest archived revision, with all repeats exported for audit. A source snapshot with missing dates is excluded and counted in testing_missing_dates_audit.csv; no report date is guessed from its download timestamp.</li></ul></div>",
    "<h2>Prediction design</h2><ul><li><strong>Contemporaneous estimation:</strong> wastewater Sunday–Saturday is paired with clinical Thursday–Wednesday ending three days earlier. The windows overlap for four days. This is a retrospective estimate, not a same-window measurement or early-warning forecast.</li>",
    "<li><strong>Sample-date lead:</strong> wastewater is paired with the clinical week ending 11 days later, starting five days after wastewater collection week end. Actual publication delays are unavailable; this is not validated operational lead time.</li>",
    "<li>Within each county, the last ",TEST_WEEKS," eligible weeks form an untouched test set; the preceding ",TUNING_WEEKS," eligible weeks select the wastewater model; earlier weeks form training (at least ",MIN_TRAIN_WEEKS,"). A ",EMBARGO_DAYS,"-day embargo separates partitions. No random week splitting is used. The previous clinical value used by all comparison models must itself be from a high-testing week.</li>",
    "<li>Wastewater-only candidates: linear log(1+WVAL), a 3-df natural spline, and log(1+WVAL) plus its previous-week change. Response: logit positivity, clipped only at 0.05%/99.95% for transformation. The tuning set selects by MAE; training+tuning are then refitted once and tested. Training mean, last-observed positivity, an autoregression, and wastewater plus prior positivity are comparisons.</li>",
    "<li>In the 11-day lead task the persistence input is clinical positivity 14 days before the target; using target-minus-seven would use a clinical value from after wastewater collection. All models are compared on exactly the same retained weeks.</li>",
    "<li>MAE/RMSE/bias are in percentage points. Skill = 1 − model MAE / persistence MAE; positive is better. R² uses held-out observations and may be negative. Nominal 95% intervals use paired moving blocks of ",BLOCK_WEEKS," consecutive weeks and ",BOOTSTRAP_REPLICATES," replicates; short tests provide imprecise, exploratory uncertainty. Intervals are omitted if gaps prevent every test observation from belonging to a consecutive block. They are not individual prediction intervals and are not adjusted for testing many counties.</li>",
    "<li>Exploratory pass criteria: MAE ≤ ",MAX_ACCEPTABLE_MAE_PP," points, at least ",100*MIN_MAE_SKILL,"% MAE improvement over persistence, and |bias| ≤ ",MAX_ACCEPTABLE_ABS_BIAS_PP," points. A positive lower bootstrap skill bound is reported separately. Thresholds are fixed, not tuned to make a model pass.</li></ul>",
    "<h2>Held-out primary results</h2>",wp_table(selected_show),
    "<h2>Held-out sensitivity results</h2><p>These use different eligibility rules and do not replace the primary analysis. A positive point estimate of skill whose interval includes zero does not establish an improvement over persistence. All test periods contain only eight observations; fragmented periods can preclude a block-bootstrap interval.</p>",wp_table(sensitivity_show),
    "<h2>California counties reviewed previously</h2><p>Paired weeks meet the high-testing and fixed-site rules. Common model weeks additionally need a previous high-testing clinical observation and a previous-week wastewater value. Full evaluation needs at least 16 training, eight tuning and eight test weeks, with calendar gaps between them.</p>",wp_table(focus_audit),
    "<h2>County testing and surveillance histories</h2>",paste(images,collapse=''),
    "<h2>Descriptive county correlations</h2>",wp_table(cshow),
    "<h2>Limits and how to extend this analysis</h2><div class='card'><p>The national WVAL archive begins in January 2022, after much of the mass-testing period. A strong in-sample association over one wave cannot supply the missing training and validation years. If few counties qualify, earlier raw wastewater concentrations (for example California's 2021 archive) are the appropriate extension; lowering eligibility requirements just to obtain a favorable model is not.</p>",
    "<p>High testing reduces one source of instability but does not remove differences in who seeks tests, repeated tests, home testing, laboratory reporting, sewershed coverage, variants, shedding or weather. Positivity can remain useful for trends after mass testing declines; it is not rendered universally useless on one date. A future estimate from this calibration would be a hypothetical positivity under the reference testing regime, not observed positivity or population prevalence.</p>",
    "<p>CDC retrospectively recalibrates the numeric WVAL history using later data. Chronological holdout protects the regression fitting, but cannot remove that upstream look-ahead. Consequently these results concern the current revised historical series; a real-time forecast needs archived as-of WVALs and publication timestamps. No 2026 positivity values are generated or inserted into the map by this script.</p></div>",
    "<h2>Sources and supporting research</h2><ul>",
    "<li><a href='https://data.cdc.gov/d/atcp-73re'>CDC site-level wastewater WVAL archive</a>; <a href='https://data.cdc.gov/d/jgk8-6dpn'>CDC revised county positivity archive</a>.</li>",
    "<li><a href='https://healthdata.gov/d/di4u-7yu6'>HHS county community profiles</a> and <a href='https://healthdata.gov/d/pcq6-asfs'>version archive</a>: original county testing count/rate snapshots, not state substitutions.</li>",
    "<li><a href='https://wwwnc.cdc.gov/eid/article/29/8/22-1634_article'>McManus et al. (2023), Predicting COVID-19 Incidence Using Wastewater Surveillance Data, Denmark, 2021–2022</a>: modeling explicitly accounts for changing testing and variants; incidence is a different target from positivity.</li>",
    "<li><a href='https://wwwnc.cdc.gov/eid/article/30/10/24-0637_article'>Lininger et al. (2024), Oregon wastewater and respiratory tests</a>: strong statewide SARS-CoV-2 rank correlation, which by itself does not establish county prediction accuracy.</li>",
    "<li><a href='https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2026.1827546/full'>Kalangara et al. (2026), comparison of US laboratory surveillance sources</a>: declining volumes and differing positivity levels coexist with useful trend agreement.</li>",
    "<li><a href='https://www.cdc.gov/wastewater/about/wval.html'>CDC WVAL methods and historical recalibration</a>.</li></ul>",
    "<h2>Exports</h2><p>All reviewed counties, exclusions, testing months, paired data, tuning scores, untouched-test predictions and model objects are exported beside this report. Raw HHS snapshots remain in the reusable cache; source URLs and checksums are retained per run. Preserve that cache for reproducibility. The HTML embeds its figures and works independently.</p><p class='small'>Input map run: ",wp_escape(inputs$root),"</p></main></html>")
  writeLines(enc2utf8(h),file.path(run,'index.html'),useBytes=TRUE)
  list(conclusion=conclusion,summary=summary_rows)
}

run_us_wastewater_positivity_prediction <- function() {
  miss<-c('curl','jsonlite')[!vapply(c('curl','jsonlite'),requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss))wp_stop('Install in a clean R session: ',paste(miss,collapse=', '))
  if(!PRIMARY_REFERENCE_FRACTION%in%SENSITIVITY_FRACTIONS||MIN_TRAIN_WEEKS<PANEL_SELECTION_WEEKS||EMBARGO_DAYS<14||BOOTSTRAP_REPLICATES<99)wp_stop('Invalid analysis settings.')
  set.seed(RANDOM_SEED);dir.create(OUTPUT_ROOT,recursive=TRUE,showWarnings=FALSE);cache<-file.path(OUTPUT_ROOT,'cache');dir.create(cache,showWarnings=FALSE)
  run<-tempfile(paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S',tz='UTC'),'_'),tmpdir=OUTPUT_ROOT);dir.create(run)
  for(n in c('sources','data','plots'))dir.create(file.path(run,n))
  complete<-FALSE;on.exit({capture.output(sessionInfo(),file=file.path(run,'session_info.txt'));if(!complete)writeLines('INCOMPLETE: inspect R error.',file.path(run,'RUN_FAILED.txt'))},add=TRUE)
  settings<-mget(c('REFERENCE_START','REFERENCE_END','ANALYSIS_START','ANALYSIS_END','MIN_MONTH_REPORTS','MIN_MONTH_SPAN_DAYS','MIN_REFERENCE_MONTHS','MIN_WEEKLY_TESTS','MIN_WEEKLY_TESTS_PER_100K','PRIMARY_REFERENCE_FRACTION','SENSITIVITY_FRACTIONS','MAX_TESTING_REPORT_AGE_DAYS','EXCLUDE_SHARED_CATCHMENTS','PANEL_SELECTION_WEEKS','MIN_PANEL_SELECTION_OBSERVATIONS','MIN_CORRELATION_WEEKS','FOCUS_FIPS','MIN_TRAIN_WEEKS','TUNING_WEEKS','TEST_WEEKS','EMBARGO_DAYS','BOOTSTRAP_REPLICATES','BLOCK_WEEKS','RANDOM_SEED','MAX_ACCEPTABLE_MAE_PP','MIN_MAE_SKILL','MAX_ACCEPTABLE_ABS_BIAS_PP'),envir=environment(run_us_wastewater_positivity_prediction))
  dput(settings,file.path(run,'settings.R'))
  inputs<-wp_inputs(run);testing<-wp_testing(cache,run,inputs$areas);panels<-wp_panels(inputs$w);wp_csv(panels,file.path(run,'data','fixed_site_panels.csv'))
  all_results<-all_cor<-list();primary_p<-primary_months<-NULL
  for(frac in SENSITIVITY_FRACTIONS) {
    qualified<-wp_qualify(testing,inputs$p,frac)
    wp_csv(qualified$m,file.path(run,'data',paste0('monthly_eligibility_',frac,'.csv')))
    if(frac==PRIMARY_REFERENCE_FRACTION){primary_p<-qualified$p;primary_months<-qualified$m}
    for(offset in c(-3L,11L)) {
      pair<-wp_pairs(inputs$w,qualified$p,panels,offset);tag<-paste(frac,pair$task[1],sep='_')
      message('Evaluating ',tag)
      wp_csv(pair,file.path(run,'data',paste0('pairs_',tag,'.csv')))
      cr<-wp_correlations(pair)
      if(!is.null(cr)&&nrow(cr)){cr$fraction<-frac;cr$task<-pair$task[1];all_cor[[tag]]<-cr}
      all_results[[tag]]<-wp_evaluate(pair,frac)
    }
  }
  bind<-function(name){x<-lapply(all_results,`[[`,name);x<-x[vapply(x,nrow,integer(1))>0];if(length(x))do.call(rbind,x)else data.frame()}
  results<-list(summary=bind('summary'),predictions=bind('predictions'),scores=bind('scores'),audit=bind('audit'))
  correlations<-if(length(all_cor))do.call(rbind,all_cor)else data.frame(fips=character(),fraction=numeric(),task=character(),pearson_logw=numeric(),spearman=numeric())
  for(n in names(results))wp_csv(results[[n]],file.path(run,'data',paste0(n,'.csv')))
  wp_csv(correlations,file.path(run,'data','descriptive_correlations.csv'))
  saveRDS(unlist(lapply(all_results,`[[`,'models'),recursive=FALSE),file.path(run,'fitted_models.rds'))
  report<-wp_report(run,inputs,testing,primary_p,primary_months,correlations,results,settings)
  wp_csv(report$summary,file.path(run,'evaluation_summary.csv'))
  writeLines(report$conclusion,file.path(run,'conclusion.txt'))
  complete<-TRUE;writeLines(paste('COMPLETE',wp_utc()),file.path(run,'RUN_COMPLETE.txt'))
  message(report$conclusion);message('Report: ',normalizePath(file.path(run,'index.html'),winslash='/'))
  if(OPEN_HTML)utils::browseURL(file.path(run,'index.html'))
  invisible(list(directory=normalizePath(run,winslash='/'),results=results,correlations=correlations))
}
if(isTRUE(RUN_ANALYSIS))run_us_wastewater_positivity_prediction()
