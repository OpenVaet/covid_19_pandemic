# England and Wales: out-of-sample validation of monthly sickness-rate baselines
# Version 1.0 | R >= 4.3 | Outcome: sickness absence rate, NOT mortality.
# Install once: install.packages(c("forecast", "ggplot2", "jsonlite"))
# Run the country scripts first, then Source this whole file in RStudio.
# Default inputs: latest successful country exports. Exact run folders or ZIPs
# may be specified below. Windows example: "C:/Research/run_20261022.zip".
# No internet, Pandoc, R Markdown or separate helper script is needed to run.
# Literature and exact model/interval definitions are included in report.html.
#
# Retrospective forecasting experiment: final/revised source values, not vintage
# real-time releases. All data after 2018 are discarded BEFORE forecasting.
# Monthly protocol assumes the previous month's value is already available.
# Actual publication delays are not simulated; this is a statistical hindcast.
# Selection is locked using 2013-2016; 2017-2018 is a later evaluation period.
# No guarantee of future coverage follows from covering every historical month.

ENGLAND_ROOT <- "england_sickness_results"
WALES_ROOT <- "wales_sickness_results"
ENGLAND_RUN <- getOption("baseline.england_run", NULL)
WALES_RUN <- getOption("baseline.wales_run", NULL)
OUT_DIR <- getOption("baseline.output_dir", "sickness_baseline_validation_results")
RUN_ANALYSIS <- getOption("baseline.run_analysis", TRUE)
MONTHLY_INFORMATION_LAG <- 1L
LEVELS <- c(50, 80, 90, 95, 98)
START_DATE <- as.Date("2010-01-01")
FIRST_TARGET <- as.Date("2011-01-01")
COMMON_START <- as.Date("2013-01-01")
SELECTION_END <- as.Date("2016-12-01")
END_DATE <- as.Date("2018-12-01")
# Percentages below are percentage points of available FTE calendar days lost.
# Optional practical point-error threshold; not a statistical significance cutoff.
POINT_TOLERANCE_PP <- 0.25
SEED <- 20261022L

write_csv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
escape_html <- function(x) {
  x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE); x <- gsub(">","&gt;",x,fixed=TRUE)
  gsub('"',"&quot;",x,fixed=TRUE)
}
bind_rows <- function(xs) { xs <- Filter(function(x)!is.null(x)&&nrow(x)>0,xs); if(length(xs)) do.call(rbind,xs) else data.frame() }
obs_file <- "monthly_observations_2010_2025.csv"
is_complete_run <- function(p) {
  if (!dir.exists(p)||file.exists(file.path(p,"FAILED.txt"))||!file.exists(file.path(p,obs_file))) return(FALSE)
  if (file.exists(file.path(p,"RUN_COMPLETE"))) return(TRUE)
  f <- file.path(p,"model_errors.csv")
  file.exists(file.path(p,"report.html")) && file.exists(f) &&
    isTRUE(tryCatch(nrow(read.csv(f))==0L,error=function(e)FALSE))
}
resolve_run <- function(spec,root,country) {
  if(is.null(spec)) {
    ds <- list.dirs(root,recursive=FALSE,full.names=TRUE)
    ds <- ds[grepl("^run_",basename(ds)) & vapply(ds,is_complete_run,logical(1))]
    if(!length(ds)) stop(country,": no successful exports in ",root,". Run the country script or set ",toupper(country),"_RUN to a folder or ZIP.")
    spec <- ds[order(basename(ds),decreasing=TRUE)][1]
  }
  original <- normalizePath(spec,winslash="/",mustWork=TRUE); cleanup <- NULL
  if(dir.exists(original)) p <- original else {
    if(!grepl("\\.zip$",original,ignore.case=TRUE)) stop("Input must be a run folder or ZIP")
    entries <- utils::unzip(original,list=TRUE)$Name
    if(any(grepl("(^[/\\\\])|(^[A-Za-z]:)|((^|[/\\\\])\\.\\.([/\\\\]|$))",entries))) stop("Unsafe ZIP entry path")
    target <- entries[basename(entries)==obs_file]
    if(length(target)!=1L) stop(country,": ZIP must contain exactly one observation file")
    cleanup <- tempfile("baseline_input_"); dir.create(cleanup)
    utils::unzip(original,exdir=cleanup); p <- file.path(cleanup,dirname(target))
  }
  if(!is_complete_run(p)) stop(country,": input export is incomplete: ",original)
  list(path=p,original=original,cleanup=cleanup)
}
prepare_series <- function(raw,country) {
  req <- c("organisation_id","staff_id","date","rate_pct")
  if(!all(req %in% names(raw))) stop(country,": unexpected observation schema")
  code <- if(country=="England") "E" else "W"
  d <- raw[!is.na(raw$organisation_id)&!is.na(raw$staff_id)&raw$organisation_id==code&raw$staff_id=="XXX",req]
  d$date <- as.Date(as.character(d$date),format="%Y-%m-%d")
  if(anyNA(d$date)) stop(country,": invalid dates")
  d <- d[d$date>=START_DATE & d$date<=END_DATE,]
  d <- d[order(d$date),]; d$rate_pct <- suppressWarnings(as.numeric(as.character(d$rate_pct)))
  if(anyDuplicated(d$date)||!identical(as.character(d$date),as.character(seq(START_DATE,END_DATE,by="month")))) stop(country,": require every month January 2010-December 2018, exactly once")
  if(any(!is.finite(d$rate_pct))||any(d$rate_pct<0|d$rate_pct>100)) stop(country,": missing or invalid rate; no imputation performed")
  data.frame(country=country,date=d$date,rate_pct=d$rate_pct,t=seq_len(nrow(d)),
    year=as.integer(format(d$date,"%Y")),month=as.integer(format(d$date,"%m")))
}

model_catalogue <- function() {
  rows <- list()
  add <- function(id,label,family,window,min_n,details) {
    rows[[length(rows)+1L]] <<- data.frame(model=id,label=label,family=family,window=window,min_n=min_n,details=details)
  }
  add("snaive","Last year's same month","snaive",0,12,
      "Point forecast repeats latest same calendar month. Normal PI uses RMS annual differences from available history; PI unavailable until 12 such differences (24 training months).")
  add("smooth3","Last-year cyclic 3-month smooth","smooth3",0,12,
      "Average adjacent calendar months in the last 12-month cycle, wrapping within that observed cycle. Normal PI uses RMS past errors at the SAME forecast horizon; at least 12 errors per horizon. Not centred smoothing of future target data.")
  for(w in c(12,24,36,0)) {
    tag <- if(w==0) "exp" else as.character(w); label <- if(w==0) "expanding" else paste0(w," months")
    add(paste0("harmonic_",tag),paste0("Smooth seasonal curve / ",label),"harmonic",w,if(w==0)36 else w,
        "OLS intercept plus annual and semiannual sine/cosine terms (K=2); no trend. Student-t prediction intervals include new-observation noise and coefficient uncertainty; iid Gaussian errors assumed.")
  }
  for(fam in c("seasonal_mean","seasonal_trend","harmonic_trend")) for(w in c(24,36,0)) {
    tag <- if(w==0) "exp" else as.character(w); label <- if(w==0) "expanding" else paste0(w," months")
    name <- switch(fam,seasonal_mean="Calendar-month averages",seasonal_trend="Trend + calendar month",harmonic_trend="Trend + smooth seasonality")
    detail <- switch(fam,seasonal_mean="OLS calendar-month factors; equivalent to averaging each month's past observations. Common residual variance pooled across months.",seasonal_trend="OLS linear time trend plus 11 calendar-month indicators.",harmonic_trend="OLS linear time trend plus annual and semiannual sine/cosine terms (K=2).")
    add(paste0(fam,"_",tag),paste0(name," / ",label),fam,w,if(w==0)36 else w,
        paste(detail,"Student-t prediction intervals; iid Gaussian errors assumed."))
  }
  for(fam in c("ets","sarima","harmonic_ar1")) for(w in c(36,0)) {
    tag <- if(w==0) "exp" else as.character(w); label <- if(w==0) "expanding" else "36 months"
    name <- switch(fam,ets="Additive ETS",sarima="Automatic seasonal ARIMA",harmonic_ar1="Smooth seasonality + AR(1)")
    detail <- switch(fam,ets="forecast::ets(model='AZZ', additive.only=TRUE, ic='aicc'): additive errors; trend/seasonality and damping selected on training data only.",sarima="forecast::auto.arima: seasonal period 12, p/q<=2, P/Q<=1, d<=1, D<=1, stepwise=TRUE, approximation=FALSE, AICc. Search is bounded, not exhaustive.",harmonic_ar1="Annual and semiannual Fourier regressors with fixed AR(1) errors, intercept, no trend; estimated by maximum likelihood.")
    add(paste0(fam,"_",tag),paste0(name," / ",label),fam,w,36,
        paste(detail,"Model-conditional normal prediction intervals; parameter/model-selection uncertainty is not fully propagated."))
  }
  bind_rows(rows)
}
features <- function(t) data.frame(t=t,month=factor((t-1L)%%12L+1L,levels=1:12),
  s1=sin(2*pi*t/12),c1=cos(2*pi*t/12),s2=sin(4*pi*t/12),c2=cos(4*pi*t/12))
cyclic_prediction <- function(y,h=1L) {
  n <- length(y); pos <- ((h-1L)%%12L)+1L; cycle <- tail(y,12)
  vapply(pos,function(j)mean(cycle[((c(j-1L,j,j+1L)-1L)%%12L)+1L]),numeric(1))
}
predict_model <- function(history,spec,h) {
  # The function receives HISTORY ONLY. No access to target outcomes or later data.
  if(nrow(history)<spec$min_n) return(NULL)
  d <- if(spec$window>0) tail(history,spec$window) else history
  y <- d$rate_pct; future_t <- max(d$t)+seq_len(h); fam <- spec$family
  lo <- hi <- matrix(NA_real_,h,length(LEVELS)); scale <- rep(NA_real_,h)
  train_x <- features(d$t); train_x$y <- y; new_x <- features(future_t)
  detail <- ""; residuals <- numeric(); df <- NA_real_
  if(fam %in% c("snaive","smooth3")) {
    if(fam=="snaive") {
      center <- rep(tail(y,12),length.out=h)
      residuals <- if(length(y)>12) y[13:length(y)]-y[1:(length(y)-12)] else numeric()
      detail <- "Seasonal naive; variance from past annual differences"
    } else {
      center <- cyclic_prediction(y,seq_len(h))
      # Match the uncertainty calibration to the actual forecast horizon.
      # Annual-frozen h=12 errors must not be calibrated using only h=1 errors.
      for(j in seq_len(h)) {
        past_origins <- if(length(y)-j>=12L)seq.int(12L,length(y)-j) else integer()
        err <- vapply(past_origins,function(o)y[o+j]-cyclic_prediction(y[seq_len(o)],j),numeric(1))
        if(length(err)>=12L)scale[j] <- sqrt(mean(err^2))
      }
      detail <- "Cyclic 3-month last-year curve; scale from historical same-horizon errors"
    }
    if(fam=="snaive" && length(residuals)>=12) {
      scale <- rep(sqrt(mean(residuals^2)),h)*sqrt(floor((seq_len(h)-1)/12)+1)
    }
    for(k in seq_along(LEVELS)) { q <- qnorm((1+LEVELS[k]/100)/2); lo[,k] <- center-q*scale; hi[,k] <- center+q*scale }
  } else if(fam %in% c("harmonic","seasonal_mean","seasonal_trend","harmonic_trend")) {
    form <- switch(fam,harmonic=y~s1+c1+s2+c2,seasonal_mean=y~month,
      seasonal_trend=y~t+month,harmonic_trend=y~t+s1+c1+s2+c2)
    fit <- stats::lm(form,data=train_x)
    if(fit$rank!=length(coef(fit))||df.residual(fit)<3) stop("Rank-deficient or insufficient regression degrees of freedom")
    center <- as.numeric(predict(fit,newdata=new_x)); df <- df.residual(fit)
    for(k in seq_along(LEVELS)) {
      z <- predict(fit,newdata=new_x,interval="prediction",level=LEVELS[k]/100)
      lo[,k] <- z[,"lwr"]; hi[,k] <- z[,"upr"]
    }
    scale <- (hi[,1]-lo[,1])/(2*qt((1+LEVELS[1]/100)/2,df))
    residuals <- residuals(fit); detail <- paste(deparse(form),collapse=" ")
  } else {
    ts_y <- ts(y,start=c(d$year[1],d$month[1]),frequency=12)
    if(fam=="ets") fit <- forecast::ets(ts_y,model="AZZ",additive.only=TRUE,ic="aicc")
    if(fam=="sarima") fit <- forecast::auto.arima(ts_y,seasonal=TRUE,ic="aicc",max.p=2,max.q=2,
      max.P=1,max.Q=1,max.d=1,max.D=1,start.p=1,start.q=1,start.P=1,start.Q=1,
      stepwise=TRUE,approximation=FALSE,allowdrift=TRUE,allowmean=TRUE)
    if(fam=="harmonic_ar1") {
      xr <- as.matrix(train_x[,c("s1","c1","s2","c2")])
      fit <- forecast::Arima(ts_y,order=c(1,0,0),seasonal=c(0,0,0),xreg=xr,include.mean=TRUE,method="ML")
      if(!is.null(fit$code)&&fit$code!=0) stop("AR(1) optimiser did not converge")
      fc <- forecast::forecast(fit,h=h,xreg=as.matrix(new_x[,colnames(xr)]),level=LEVELS)
    } else fc <- forecast::forecast(fit,h=h,level=LEVELS)
    if(!is.null(fit$code)&&fit$code!=0) stop("ARIMA optimiser did not converge")
    center <- as.numeric(fc$mean); lo <- as.matrix(fc$lower); hi <- as.matrix(fc$upper)
    scale <- (hi[,1]-lo[,1])/(2*qnorm((1+LEVELS[1]/100)/2))
    residuals <- as.numeric(residuals(fit)); detail <- fc$method
  }
  if(length(center)!=h||any(!is.finite(center))) stop("Invalid point forecast")
  valid_pi <- rowSums(is.finite(lo)&is.finite(hi))==length(LEVELS)
  if(any(valid_pi & (rowSums(lo>hi,na.rm=TRUE)>0))) stop("Reversed prediction bounds")
  z <- data.frame(h=seq_len(h),prediction=center,scale=scale,df=df,pi_available=valid_pi,
    n_train=nrow(d),train_start=min(d$date),train_end=max(d$date),fit_description=detail)
  for(k in seq_along(LEVELS)) { z[[paste0("lower",LEVELS[k])]] <- lo[,k]; z[[paste0("upper",LEVELS[k])]] <- hi[,k] }
  z
}

make_forecasts <- function(series,catalogue) {
  rows <- issues <- list(); count <- 0L
  for(country in names(series)) for(protocol in c("monthly_update","annual_frozen")) {
    d <- series[[country]]
    origins <- if(protocol=="monthly_update") which(d$date>=FIRST_TARGET)-MONTHLY_INFORMATION_LAG else which(d$month==12 & d$year<2018)
    origins <- sort(unique(origins[origins>=1L]))
    for(o in origins) {
      targets <- if(protocol=="monthly_update") o+MONTHLY_INFORMATION_LAG else o+seq_len(12)
      targets <- targets[targets<=nrow(d) & d$date[targets]>=FIRST_TARGET]
      if(!length(targets)) next
      hs <- targets-o; history <- d[seq_len(o),]
      for(m in seq_len(nrow(catalogue))) {
        spec <- catalogue[m,]; warns <- character()
        fit <- tryCatch(withCallingHandlers(predict_model(history,spec,max(hs)),warning=function(w) {
          warns <<- c(warns,conditionMessage(w)); invokeRestart("muffleWarning")
        }),error=function(e)e)
        if(inherits(fit,"error") || length(warns)) issues[[length(issues)+1L]] <- data.frame(country=country,protocol=protocol,
          origin=as.character(d$date[o]),model=spec$model,severity=if(inherits(fit,"error"))"error" else "warning",
          message=if(inherits(fit,"error"))conditionMessage(fit) else paste(unique(warns),collapse="; "))
        if(is.null(fit)||inherits(fit,"error")) next
        z <- fit[hs,]; z$country <- country; z$protocol <- protocol; z$model <- spec$model
        z$date <- d$date[targets]; z$actual <- d$rate_pct[targets]; z$origin <- d$date[o]
        stopifnot(all(z$train_end<z$date),all(z$date<=END_DATE))
        rows[[length(rows)+1L]] <- z
      }
      count <- count+1L
      if(count%%12L==0) message(country," / ",protocol," / origin ",d$date[o])
    }
  }
  forecasts <- bind_rows(rows)
  issues <- if(length(issues))bind_rows(issues) else data.frame(country=character(),protocol=character(),origin=character(),model=character(),severity=character(),message=character())
  list(forecasts=forecasts,issues=issues)
}
score_forecasts <- function(d) {
  d$error <- d$actual-d$prediction; d$absolute_error <- abs(d$error)
  d$squared_error <- d$error^2; d$within_tolerance <- d$absolute_error<=POINT_TOLERANCE_PP
  wsum <- 0.5*d$absolute_error
  for(lev in LEVELS) {
    lo <- d[[paste0("lower",lev)]]; hi <- d[[paste0("upper",lev)]]; a <- 1-lev/100
    d[[paste0("covered",lev)]] <- d$actual>=lo & d$actual<=hi
    d[[paste0("width",lev)]] <- hi-lo
    is <- hi-lo+2/a*(lo-d$actual)*(d$actual<lo)+2/a*(d$actual-hi)*(d$actual>hi)
    d[[paste0("interval_score",lev)]] <- is; wsum <- wsum+a/2*is
  }
  d$wis <- wsum/(length(LEVELS)+0.5)
  d$year <- as.integer(format(d$date,"%Y")); d
}
summarise_rows <- function(d,expected,period) {
  pi <- is.finite(d$wis); n <- nrow(d); np <- sum(pi)
  out <- data.frame(period=period,n_forecasts=n,n_expected=expected,n_intervals=np,
    complete_points=n==expected,complete_intervals=np==expected,
    MAE=mean(d$absolute_error),RMSE=sqrt(mean(d$squared_error)),bias=mean(d$error),
    point_within_tolerance_pct=100*mean(d$within_tolerance),WIS=if(np)mean(d$wis[pi]) else NA_real_)
  for(lev in LEVELS) {
    ok <- !is.na(d[[paste0("covered",lev)]]); covered <- sum(d[[paste0("covered",lev)]],na.rm=TRUE)
    out[[paste0("coverage",lev)]] <- if(sum(ok))100*covered/sum(ok) else NA_real_
    out[[paste0("covered_n",lev)]] <- covered
    out[[paste0("misses",lev)]] <- sum(ok)-covered
    out[[paste0("width",lev)]] <- if(sum(ok))mean(d[[paste0("width",lev)]][ok]) else NA_real_
    out[[paste0("all_covered",lev)]] <- sum(ok)==expected && covered==expected
  }
  out
}
summary_table <- function(d,catalogue) {
  periods <- list(early_2011_2012=c("2011-01-01","2012-12-01"),selection_2013_2016=c("2013-01-01","2016-12-01"),
    holdout_2017_2018=c("2017-01-01","2018-12-01"),common_2013_2018=c("2013-01-01","2018-12-01"),available_2011_2018=c("2011-01-01","2018-12-01"))
  rows <- list()
  for(p in names(periods)) for(pr in unique(d$protocol)) for(cn in c("England","Wales","Pooled")) for(id in catalogue$model) {
    z <- d[d$protocol==pr&d$model==id&d$date>=as.Date(periods[[p]][1])&d$date<=as.Date(periods[[p]][2]),]
    if(cn!="Pooled") z <- z[z$country==cn,]
    if(!nrow(z)) next
    expected <- length(seq(as.Date(periods[[p]][1]),as.Date(periods[[p]][2]),by="month"))*if(cn=="Pooled")2 else 1
    row <- summarise_rows(z,expected,p); row$country <- cn; row$protocol <- pr; row$model <- id
    row$first_target <- min(z$date); row$last_target <- max(z$date)
    rows[[length(rows)+1L]] <- row
  }
  x <- bind_rows(rows); x$label <- catalogue$label[match(x$model,catalogue$model)]; x
}
select_models <- function(s) {
  rows <- list()
  for(pr in unique(s$protocol)) for(cn in c("England","Wales","Pooled")) {
    z <- s[s$protocol==pr&s$country==cn&s$period=="selection_2013_2016"&s$complete_intervals&is.finite(s$WIS),]
    if(!nrow(z)) stop("No complete selection-period candidates for ",cn," / ",pr)
    z <- z[order(z$WIS,z$MAE,z$model),]
    best <- z[1,]; point <- z[order(z$MAE,z$WIS),][1,]
    rows[[length(rows)+1L]] <- data.frame(country=cn,protocol=pr,model=best$model,label=best$label,
      development_WIS=best$WIS,development_MAE=best$MAE,point_accuracy_winner=point$model,
      point_accuracy_winner_MAE=point$MAE)
  }
  bind_rows(rows)
}

by_year_table <- function(d) {
  rows <- lapply(split(d,interaction(d$country,d$protocol,d$model,d$year,drop=TRUE)),function(z) {
    out <- summarise_rows(z,12,paste0("year_",z$year[1])); out$country <- z$country[1]
    out$protocol <- z$protocol[1]; out$model <- z$model[1]; out$year <- z$year[1]; out
  }); bind_rows(rows)
}
diagnostics_table <- function(d,selected) {
  rows <- list()
  for(i in which(selected$country!="Pooled")) {
    sel <- selected[i,]; z <- d[d$country==sel$country & d$protocol==sel$protocol & d$model==sel$model & d$date>=COMMON_START,]
    z <- z[order(z$date),]; e <- z$error
    lb <- if(length(e)>24)Box.test(e,lag=12,type="Ljung-Box",fitdf=0)$p.value else NA_real_
    rows[[length(rows)+1L]] <- data.frame(country=sel$country,protocol=sel$protocol,model=sel$model,
      n=length(e),error_acf1=if(length(e)>2)acf(e,plot=FALSE,lag.max=1)$acf[2] else NA_real_,
      Ljung_Box_lag12_p_descriptive=lb,n_absent_intervals=sum(!z$pi_available),
      outside_physical_range=sum(z$lower98<0|z$upper98>100,na.rm=TRUE))
  }; bind_rows(rows)
}
coverage_by_month <- function(d) {
  d <- d[d$date>=COMMON_START,]; d$month <- as.integer(format(d$date,"%m"))
  rows <- lapply(split(d,interaction(d$country,d$protocol,d$model,d$month,drop=TRUE)),function(z) {
    out <- summarise_rows(z,6,"calendar_month_2013_2018"); out$country <- z$country[1]
    out$protocol <- z$protocol[1]; out$model <- z$model[1]; out$month <- z$month[1]; out
  }); bind_rows(rows)
}
html_table <- function(d,digits=3) {
  if(!nrow(d))return("<p>No eligible results.</p>")
  cols <- lapply(d,function(x) {
    if(is.logical(x))return(ifelse(is.na(x),"—",ifelse(x,"YES","NO")))
    if(is.numeric(x)) {
      nd <- if(all(x[is.finite(x)]==round(x[is.finite(x)])))0L else digits
      return(ifelse(is.finite(x),formatC(x,format="f",digits=nd),"—"))
    }
    ifelse(is.na(x),"—",as.character(x))
  })
  rows <- vapply(seq_len(nrow(d)),function(i)paste0("<tr>",paste0("<td>",escape_html(vapply(cols,`[`,character(1),i)),"</td>",collapse=""),"</tr>"),character(1))
  paste0('<div class="table-wrap"><table><thead><tr>',paste0("<th>",escape_html(names(d)),"</th>",collapse=""),
    "</tr></thead><tbody>",paste(rows,collapse="\n"),"</tbody></table></div>")
}
plot_theme <- function() ggplot2::theme_minimal(base_size=12)+ggplot2::theme(
  panel.grid.minor=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"),
  legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
save_figure <- function(p,name,run,width=12,height=7) {
  png <- file.path(run,"plots",paste0(name,".png"))
  ggplot2::ggsave(png,p,width=width,height=height,dpi=170,bg="white")
  ggplot2::ggsave(file.path(run,"plots",paste0(name,".pdf")),p,width=width,height=height)
  png
}
make_plots <- function(d,s,selected,run) {
  files <- list()
  for(pr in unique(d$protocol)) {
    title <- if(pr=="monthly_update")"Monthly update: one month ahead" else "Annual freeze: predict the next calendar year"
    ds <- bind_rows(lapply(which(selected$protocol==pr & selected$country!="Pooled"),function(i) {
      q <- selected[i,]; z <- d[d$country==q$country & d$protocol==pr & d$model==q$model,]
      z$panel <- paste0(q$country," | ",q$label);z
    }))
    p <- ggplot2::ggplot(ds,ggplot2::aes(date,actual))+
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower98,ymax=upper98),fill="#b7cbe5",alpha=.5,na.rm=TRUE)+
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower95,ymax=upper95),fill="#77a7d5",alpha=.5,na.rm=TRUE)+
      ggplot2::geom_line(ggplot2::aes(y=prediction,color="Forecast"),linewidth=.55)+
      ggplot2::geom_line(ggplot2::aes(color="Observed"),linewidth=.45)+
      ggplot2::geom_point(data=ds[!is.na(ds$covered95)&!ds$covered95,],color="#bf3a34",size=1.7)+
      ggplot2::geom_vline(xintercept=as.numeric(as.Date("2017-01-01")),linetype="dashed",color="#466078")+
      ggplot2::facet_wrap(~panel,ncol=1,scales="fixed")+
      ggplot2::scale_color_manual(values=c(Forecast="#14689d",Observed="#202932"))+
      ggplot2::scale_x_date(date_breaks="1 year",date_labels="%Y")+
      ggplot2::labs(title=title,subtitle="Model choices locked using 2013-2016; dashed line starts the 2017-2018 holdout",
        x=NULL,y="Sickness absence rate (%)",color=NULL,
        caption="Dark band: 95% prediction interval. Light band: 98%. Red dots: observations outside 95%. Raw monthly targets, not smoothed observations.")+plot_theme()
    files[[paste0(pr,"_forecasts")]] <- save_figure(p,paste0(pr,"_forecasts"),run,13,8)
    rank <- s[s$protocol==pr&s$country!="Pooled" & s$period%in%c("selection_2013_2016","holdout_2017_2018")&s$complete_intervals,]
    rank$period <- factor(rank$period,levels=c("selection_2013_2016","holdout_2017_2018"),labels=c("Selection: 2013-2016","Holdout: 2017-2018"))
    ordered <- aggregate(WIS~label,rank,mean);ordered <- ordered$label[order(ordered$WIS,decreasing=TRUE)]
    rank$label <- factor(rank$label,levels=ordered)
    p <- ggplot2::ggplot(rank,ggplot2::aes(WIS,label,color=period))+ggplot2::geom_point(size=2.3)+
      ggplot2::facet_wrap(~country,nrow=1)+ggplot2::scale_color_manual(values=c("#65758b","#007f7a"))+
      ggplot2::labs(title=paste0("Probabilistic accuracy | ",title),subtitle="Lower weighted interval score is better. Holdout results are shown for audit; they do not change the locked choices.",
        x="Mean WIS (percentage points)",y=NULL,color=NULL)+plot_theme()
    files[[paste0(pr,"_ranking")]] <- save_figure(p,paste0(pr,"_ranking"),run,14,9)
  }
  selected_rows <- bind_rows(lapply(which(selected$country!="Pooled"),function(i) {
    q <- selected[i,];d[d$country==q$country & d$protocol==q$protocol & d$model==q$model & d$date>SELECTION_END,]
  }))
  calibration <- bind_rows(lapply(split(selected_rows,interaction(selected_rows$country,selected_rows$protocol,drop=TRUE)),function(z) {
    data.frame(country=z$country[1],protocol=z$protocol[1],nominal=LEVELS,
      observed=vapply(LEVELS,function(l)100*mean(z[[paste0("covered",l)]],na.rm=TRUE),numeric(1)))
  }))
  p <- ggplot2::ggplot(calibration,ggplot2::aes(nominal,observed,color=country))+
    ggplot2::geom_abline(slope=1,intercept=0,color="#a3adb8",linetype="dashed")+
    ggplot2::geom_line()+ggplot2::geom_point(size=2)+ggplot2::facet_wrap(~protocol)+
    ggplot2::scale_color_manual(values=c(England="#176fa1",Wales="#007f74"))+
    ggplot2::scale_x_continuous(breaks=LEVELS)+ggplot2::coord_cartesian(ylim=c(35,100))+
    ggplot2::labs(title="Do the selected intervals cover as often as their labels suggest?",subtitle="Holdout only: 24 months per country. Coverage changes in steps of 4.17 percentage points.",
      x="Nominal prediction-interval coverage (%)",y="Observed coverage (%)",color=NULL)+plot_theme()
  files$calibration <- save_figure(p,"holdout_calibration",run,12,5)
  files
}
embed_image <- function(p,alt) paste0('<img loading="lazy" alt="',escape_html(alt),'" src="data:image/png;base64,',jsonlite::base64_enc(readBin(p,"raw",n=file.info(p)$size)),'">')

literature <- function() data.frame(
  source=c("Asghar et al. (2021), BMJ Open 11:e053885","Hyndman & Athanasopoulos, Forecasting: Principles and Practice, ch. 5",
    "Hyndman & Khandakar (2008), Journal of Statistical Software 27(3)","Hyndman & Athanasopoulos, ch. 7",
    "Hyndman & Athanasopoulos, ch. 10","Bracher et al. (2021), PLOS Computational Biology 17:e1008618",
    "Bracher et al. (2022), correction","forecast package: ETS and auto.arima documentation"),
  url=c("https://doi.org/10.1136/bmjopen-2021-053885","https://otexts.com/fpp3/toolbox.html",
    "https://doi.org/10.18637/jss.v027.i03","https://otexts.com/fpp3/regression.html",
    "https://otexts.com/fpp3/dynamic-regression.html","https://doi.org/10.1371/journal.pcbi.1008618",
    "https://doi.org/10.1371/journal.pcbi.1010592","https://pkg.robjhyndman.com/forecast/reference/auto.arima.html"),
  contribution=c("Directly relevant NHS ambulance study used monthly 2009–2018 data and ARIMA/SARIMA; seasonality was useful for forecasting. This motivates a candidate, not a presumption of superiority for national all-staff rates.",
    "Rolling-origin evaluation keeps future outcomes out of model fitting. Seasonal naive forecasts provide a demanding simple benchmark. Future observations require prediction intervals, not confidence intervals for their means.",
    "Introduces automatic ARIMA and exponential-smoothing state-space forecasting. We refit each candidate using each origin's training history.",
    "Calendar-month effects and smooth Fourier terms offer transparent seasonal regression baselines. Linear trends are optional candidates, rather than assumed improvements.",
    "Regression with autocorrelated errors can preserve smooth seasonality while accounting for persistence. We include a deliberately small Fourier-plus-AR(1) model.",
    "Weighted interval score combines sharpness and penalties for missed observations. We score all candidates at the same five interval levels and include the median absolute-error term; smaller WIS is better.",
    "The correction concerns example WIS values and figures; it is linked alongside the original scoring paper.",
    "Documents training-only AICc selection and the bounded stepwise ARIMA search. Software versions are saved; package updates may change selected internal specifications."))
conclusion_html <- function(s,selected) {
  paragraphs <- character()
  for(i in seq_len(nrow(selected))) {
    q <- selected[i,]; h <- s[s$country==q$country&s$protocol==q$protocol&s$model==q$model&s$period=="holdout_2017_2018",]
    b <- s[s$country==q$country&s$protocol==q$protocol&s$model=="snaive"&s$period=="holdout_2017_2018",]
    protocol <- if(q$protocol=="monthly_update")"monthly update" else "annual freeze"
    if(!nrow(h)||!h$complete_intervals) {
      paragraphs <- c(paragraphs,paste0("<p><b>",q$country," / ",protocol,": ",escape_html(q$label),"</b> was selected, but did not produce a complete holdout. No validated recommendation is issued.</p>"));next
    }
    improvement <- if(nrow(b)&&b$complete_intervals&&b$WIS>0)100*(1-h$WIS/b$WIS) else NA_real_
    assessment <- if(is.finite(improvement)&&improvement>0)paste0("Its holdout WIS was ",sprintf("%.1f",improvement),"% lower than the last-year same-month benchmark.") else
      if(is.finite(improvement))paste0("It did not improve on the benchmark in the holdout (WIS ",sprintf("%.1f",-improvement),"% higher); the selection advantage did not carry forward.") else "Benchmark comparison unavailable."
    calibration <- if(h$coverage95<90)paste0(" <b>Calibration concern:</b> observed coverage was only ",sprintf("%.1f%%",h$coverage95),
      " for nominal 95% intervals. This score winner should not be treated as providing validated 95% bounds; recalibration would need a fresh evaluation period.") else
      " The small holdout does not establish exact 95% or 98% calibration."
    paragraphs <- c(paragraphs,paste0("<p><b>",q$country," / ",protocol,": ",escape_html(q$label),".</b> Selected by 2013–2016 WIS. In 2017–2018: MAE <b>",sprintf("%.3f",h$MAE)," percentage points</b>, WIS ",sprintf("%.3f",h$WIS),
      "; 95% interval covered <b>",h$covered_n95,"/",h$n_expected,"</b> observations (mean width ",sprintf("%.3f",h$width95)," pp), and 98% covered <b>",h$covered_n98,"/",h$n_expected,"</b> (width ",sprintf("%.3f",h$width98)," pp). ",assessment,calibration,"</p>"))
  }
  paste(paragraphs,collapse="\n")
}
write_report <- function(run,d,s,selected,years,diagnostics,catalogue,issues,manifest,plots) {
  tab <- function(x)html_table(x)
  section <- function(title,body)paste0("<section><h2>",title,"</h2>",body,"</section>")
  refs <- literature()
  ref_html <- paste(vapply(seq_len(nrow(refs)),function(i)paste0('<li><a href="',refs$url[i],'">',escape_html(refs$source[i]),"</a>. ",escape_html(refs$contribution[i]),"</li>"),character(1)),collapse="")
  selected_display <- selected[,c("country","protocol","label","development_WIS","development_MAE","point_accuracy_winner")]
  names(selected_display) <- c("Series","Protocol","Selected predictor","Selection WIS","Selection MAE","Best point model ID")
  body <- section("What this analysis concludes",paste0(
    '<div class="callout">The strongest predictor depends on the country and when the forecast is made. These are the choices selected before the final two evaluation years. A lower overall score does not guarantee calibrated 95% or 98% intervals: the coverage results are part of the conclusion.</div>',
    conclusion_html(s,selected),"<p>The pooled choice applies the same model specification separately to both countries and averages their scores equally. It does not combine their rates or estimate a UK rate.</p>"))
  body <- paste0(body,section("Forecasting question and fair comparison",paste0(
    "<p>Outcome: the <b>raw monthly national all-staff sickness absence rate</b> (% of available FTE calendar days lost). The smoothed six-month curves and 2010–2018 fitted baselines in the earlier reports are not inputs: using them would introduce future information. Each country's model uses only that country's preceding observations.</p>",
    "<div class=grid><div><h3>Monthly update</h3><p>At each month, fit using data through the previous month and predict one month ahead. A 12-month window moves each month. This measures statistical predictability assuming the previous value is available.</p></div><div><h3>Annual freeze</h3><p>Fit through December, then forecast January–December without using any observations from that target year. A one-year baseline therefore uses only the preceding calendar year; forecast horizons are 1–12 months.</p></div></div>",
    "<p>Forecasts begin in <b>2011</b> where a method has enough history. Three-year models first qualify in <b>2013</b>. The fair comparison uses exactly <b>2013–2018</b>; shorter-window models are also scored over the same months. Selection uses <b>2013–2016</b>, and the selected model ID is frozen for <b>2017–2018</b>. Parameters are still refitted at each legitimate origin. Annual forecasts have only four selection years and two holdout years, so apparent differences can be fragile.</p>",
    "<p>Early-period and available-window tables are descriptive and never used to rank methods against different target months. A method must supply every selection-period prediction and interval to be eligible. Failures are logged rather than replaced by easier models or silently omitted.</p>",tab(selected_display))))
  body <- paste0(body,section("How accuracy and uncertainty are judged",paste0(
    "<p><b>Primary selection criterion:</b> mean weighted interval score (WIS), with the same 50%, 80%, 90%, 95% and 98% central prediction intervals for every method. MAE breaks exact WIS ties. A separate best-MAE model is shown because point accuracy and interval quality can favour different models. RMSE, bias, interval width and coverage are also exported.</p>",
    "<p>For interval [L,U] with tail probability α, IS = (U−L) + (2/α)(L−y)1(y&lt;L) + (2/α)(y−U)1(y&gt;U). WIS = [0.5|y−m| + Σ(α/2)IS]/(K+0.5), with K=5. It penalises both unnecessarily broad intervals and missed observations. Widths and scores are in percentage points.</p>",
    "<p>A 95% prediction interval concerns a future observation, whereas a confidence interval for the mean is narrower and answers a different question. Normal or Student-t assumptions are model-specific; none of these intervals is guaranteed to be calibrated. ARIMA/ETS intervals condition on fitted parameters and selected model structure. OLS intervals assume independent residuals, which may be unrealistic; out-of-sample coverage checks that assumption's consequences.</p>",
    "<p><b>Covering every month is not the optimisation target.</b> Even if independent, correctly calibrated 95% intervals were used, the chance of covering all 96 months would be 0.95<sup>96</sup> ≈ 0.7%; for 98%, 0.98<sup>96</sup> ≈ 14.4%. Serial dependence changes those calculations. The report nevertheless lists all-month coverage explicitly. Choosing a level after seeing misses does not create a simultaneous guarantee.</p>",
    "<p>The user-facing 95% and 98% comparisons are reported together; the script never chooses 98% merely to obtain a perfect retrospective count. Point errors within ±",POINT_TOLERANCE_PP," pp are an editable practical tolerance, not a significance test.</p>")))
  body <- paste0(body,section("Literature informing the candidate set",paste0("<ol>",ref_html,"</ol>",
    "<p>Application to these two national series is our modelling choice, not a result established by the cited ambulance study. Neural networks and other highly flexible learners are not included in this small 108-month dataset. Adding external predictors would require historical values and publication timing; none are assumed available.</p>")))
  examples <- s[s$country!="Pooled" & ((s$model=="harmonic_12" & s$period=="available_2011_2018") |
    (s$model%in%c("harmonic_36","seasonal_mean_36") & s$period=="common_2013_2018")),]
  body <- paste0(body,section("Your one-year and three-year examples, tested",paste0(
    "<p>The one-year smooth baseline fits a periodic curve with annual and semiannual waves to the previous 12 months. The three-year smooth version uses the same form over 36 months; the calendar-month average is a separate three-year alternative. This table answers the coverage examples directly; ranking still uses identical target months.</p>",
    tab(examples[,c("country","protocol","label","first_target","last_target","n_expected","MAE","covered_n95","all_covered95","width95","covered_n98","all_covered98","width98")]))))
  for(pr in unique(d$protocol)) {
    title <- if(pr=="monthly_update")"Monthly update" else "Annual freeze"
    body <- paste0(body,section(paste0(title,": selected forecasts"),paste0(embed_image(plots[[paste0(pr,"_forecasts")]],paste(title,"observations and forecasts")),
      "<p>The forecast model displayed was selected using 2013–2016. Earlier predictions from that model are retrospective descriptions, not an assertion that it would have been selected in 2011. Bands are individual-month prediction intervals; they are not simultaneous bands for the whole path.</p>")),
      section(paste0(title,": candidate comparison"),embed_image(plots[[paste0(pr,"_ranking")]],paste(title,"candidate WIS"))))
    for(period in c("selection_2013_2016","holdout_2017_2018","common_2013_2018")) {
      z <- s[s$protocol==pr & s$period==period,]
      z <- z[order(z$country,!z$complete_intervals,z$WIS),]
      view <- z[,c("country","label","n_forecasts","n_intervals","n_expected","complete_intervals","MAE","RMSE","WIS","coverage95","width95","coverage98","width98")]
      label <- switch(period,selection_2013_2016="Selection 2013–2016",holdout_2017_2018="Holdout 2017–2018 — audit, no reselection",common_2013_2018="Full common period 2013–2018 — retrospective ranking")
      body <- paste0(body,'<details><summary>',title,": ",label,"</summary>",tab(view),"</details>")
    }
  }
  body <- paste0(body,section("95% versus 98%: did every month fall inside?",paste0(
    "<p>Counts below are from the common 2013–2018 window (72 months per country). YES means every required observation had an interval and fell within it. It does not mean the point predictions were exact.</p>",
    tab(s[s$period=="common_2013_2018" & s$country!="Pooled",c("country","protocol","label","n_expected","n_intervals","misses95","all_covered95","width95","misses98","all_covered98","width98")]))))
  body <- paste0(body,section("Calibration and stability",paste0(embed_image(plots$calibration,"Holdout nominal versus observed interval coverage"),
    "<p>There are only 24 holdout observations per country: a single miss changes coverage by 4.17 percentage points. Nominal 95% versus 98% cannot be finely distinguished in such a small sample. The two countries and adjacent months are correlated; the pooled result is not 48 independent trials.</p>",
    "<p>The following forecast-error diagnostics are descriptive across 2013–2018, including the selection period. Ljung–Box p-values are not adjusted for model selection, and annual forecast errors are dependent within a forecast year. Use them as prompts for refinement, not confirmatory tests.</p>",tab(diagnostics))))
  sel_years <- bind_rows(lapply(which(selected$country!="Pooled"),function(i) {
    q<-selected[i,];years[years$country==q$country & years$protocol==q$protocol & years$model==q$model,]
  }))
  body <- paste0(body,'<details><summary>Selected predictors: results by target year</summary>',tab(sel_years[,c("country","protocol","model","year","n_forecasts","n_intervals","MAE","WIS","coverage95","coverage98")]),"</details>")
  early <- s[s$country!="Pooled" & s$period%in%c("early_2011_2012","available_2011_2018"),]
  body <- paste0(body,'<details><summary>Early years and model-specific available windows (not a fair leaderboard)</summary>',
    "<p>Rows with fewer than the expected number of intervals cannot establish all-period coverage. Same-month naive and cyclic smoothing have point forecasts from 2011, but no intervals until at least 12 past forecasting errors are available. Cyclic smoothing uses errors at the same forecast horizon, so long-horizon intervals become available later. A smooth harmonic regression can estimate a pooled residual variance from a single year, at the cost of stronger shape and error assumptions.</p>",
    tab(early[,c("period","country","protocol","label","first_target","last_target","n_forecasts","n_intervals","n_expected","MAE","WIS","coverage95","coverage98","all_covered95","all_covered98")]),"</details>")
  body <- paste0(body,section("Exact candidates and implementation",tab(catalogue)),
    section("Data provenance and limits",paste0(tab(manifest),
      "<p>Both countries supply a complete 108-month national all-staff series for January 2010–December 2018; no missing months, duplicate months or invalid rates are accepted. No values are imputed and no outliers are removed. The input is an existing official-data export, not a new web download. Original CSVs, metadata, definition-change notes and revision-check files, where present, are copied to the run's source folder.</p>",
      "<p>This is a <b>retrospective pseudo-out-of-sample analysis of the supplied data vintage</b>. Original publication vintages and historical release delays are unavailable: revisions may improve historical predictability, and the previous month might not have been published at the forecast date. Identical data and software versions are needed for exact replication.</p>",
      "<p>England and Wales have different workforce composition and recording systems. Each is modelled separately. Recorded definition changes are preserved from the input runs; no automatic continuity correction is inferred. Rates are rounded aggregate administrative measures, not independent person-level observations.</p>",
      "<p>The supplied Wales notes flag GP-trainee employment coverage from 2015 and staff-category reclassification during 2015–2017. Reclassification between staff groups may cancel in the all-staff aggregate, whereas coverage changes need not; this is not a proven correction or an imposed breakpoint. These historical features remain part of the observed forecasting task.</p>",
      "<p>Good one-month-ahead forecasting does not establish a stable counterfactual for detecting a future structural break: rolling models can adapt to and absorb a persistent increase. Selection over a small candidate set and a short historical period remains uncertain. Results support conclusions only for the assessed pre-2019 forecasting task.</p>")),
    '<details><summary>Model warnings and failures (',nrow(issues),")</summary>",tab(issues),"</details>",
    section("Reproducing this report",paste0("<p>Run <code>validate_sickness_baselines.R</code> after creating the England and Wales exports. Set <code>ENGLAND_RUN</code> and <code>WALES_RUN</code> to exact folders or ZIPs to reproduce a chosen snapshot. Dependencies: forecast, ggplot2, jsonlite. No HTML rendering tool is required.</p>",
      "<p>Exports include all month-level forecasts and intervals, model catalogue, full score tables, yearly and calendar-month diagnostics, locked selections, missed-observation tables, source copies/hashes, configuration and session information. Plots are exported as PNG and PDF. This HTML embeds its figures and can be moved on its own.</p>")))
  css <- "body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;color:#20303d;background:#f4f7fa;margin:0;line-height:1.6}main{max-width:1280px;margin:auto;background:white;padding:48px}h1{font-size:36px;line-height:1.2;max-width:950px}h2{font-size:25px;color:#123e55;margin-top:0}h3{margin:0}section{margin:40px 0;padding-bottom:20px;border-bottom:1px solid #d9e1e8}.eyebrow{color:#097d78;font-weight:700;letter-spacing:.1em;text-transform:uppercase;font-size:12px}.subtitle{font-size:18px;color:#5b6d7a}.callout{background:#e8f4f2;border-left:5px solid #087e76;padding:18px 22px;margin-bottom:20px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:24px}.grid>div{background:#f1f5fa;padding:22px}.table-wrap{overflow:auto;max-height:700px;border:1px solid #dbe3e8;margin:18px 0}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:9px 12px;border-bottom:1px solid #dde4ea;text-align:left;vertical-align:top}th{background:#183f56;color:white;position:sticky;top:0;white-space:nowrap}tr:nth-child(even){background:#f3f6f9}img{width:100%;height:auto;display:block}details{padding:16px;border:1px solid #d7e0e7;border-radius:6px;margin:16px 0}summary{font-weight:600;cursor:pointer}a{color:#006f9c}code{background:#eef2f6;padding:2px 5px}footer{font-size:12px;color:#667d8b}@media(max-width:750px){main{padding:20px}.grid{grid-template-columns:1fr}h1{font-size:28px}}@media print{body{background:white}main{padding:0;max-width:none}.table-wrap{max-height:none}details{break-inside:avoid}img{max-height:850px;object-fit:contain}}"
  html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Which baseline best predicts monthly NHS sickness absence?</title><style>',css,
    '</style></head><body><main><div class="eyebrow">England &amp; Wales · Forecast validation</div><h1>Which baseline best predicts monthly sickness absence?</h1><p class="subtitle">2011–2018 forecast targets · fair comparison from 2013 · raw monthly observations</p>',
    body,"<footer>Generated ",format(Sys.time(),tz="UTC",usetz=TRUE)," · Script version 1.0 · Outcome is sickness absence, not mortality.</footer></main></body></html>")
  writeLines(html,file.path(run,"report.html"),useBytes=TRUE)
}

run_analysis <- function() {
  needed <- c("forecast","ggplot2","jsonlite")
  missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop("Install required packages first: install.packages(c(",paste(sprintf('"%s"',missing),collapse=","),"))",call.=FALSE)
  stopifnot(identical(LEVELS,c(50,80,90,95,98)),MONTHLY_INFORMATION_LAG==1L)
  set.seed(SEED); dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  run <- tempfile(paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=OUT_DIR)
  dir.create(run);dir.create(file.path(run,"plots"));dir.create(file.path(run,"source"))
  complete <- FALSE; inputs <- list()
  on.exit({for(x in inputs)if(!is.null(x$cleanup))unlink(x$cleanup,recursive=TRUE)
    if(!complete)writeLines("Run failed. Do not use partial output for conclusions.",file.path(run,"FAILED.txt"))},add=TRUE)
  inputs$England <- resolve_run(ENGLAND_RUN,ENGLAND_ROOT,"England")
  inputs$Wales <- resolve_run(WALES_RUN,WALES_ROOT,"Wales")
  series <- list(); manifest <- list()
  for(cn in names(inputs)) {
    p <- inputs[[cn]]$path; f <- file.path(p,obs_file)
    raw <- read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
    series[[cn]] <- prepare_series(raw,cn)
    dest <- file.path(run,"source",cn);dir.create(dest)
    keep <- c(obs_file,"source_status.csv","known_definition_changes.csv","revision_changes.csv","quality_checks.csv","configuration.R","READ_BEFORE_ANALYSIS.txt","sessionInfo.txt")
    keep <- keep[file.exists(file.path(p,keep))];file.copy(file.path(p,keep),dest)
    manifest[[cn]] <- data.frame(country=cn,input=inputs[[cn]]$original,
      observations_md5=unname(tools::md5sum(f)),n_months=nrow(series[[cn]]),start=START_DATE,end=END_DATE,
      missing_months=0,duplicate_months=0,missing_rates=0,metadata_files=paste(keep,collapse="; "))
  }
  manifest <- bind_rows(manifest);catalogue <- model_catalogue()
  write_csv(bind_rows(series),file.path(run,"input_series_2010_2018.csv"))
  write_csv(manifest,file.path(run,"input_manifest.csv"));write_csv(catalogue,file.path(run,"model_catalogue.csv"))
  message("Forecasting ",nrow(catalogue)," candidate specifications independently for each country and protocol...")
  result <- make_forecasts(series,catalogue);d <- score_forecasts(result$forecasts)
  if(anyDuplicated(d[,c("country","protocol","model","date")]))stop("Duplicate forecast keys")
  s <- summary_table(d,catalogue); selected <- select_models(s)
  years <- by_year_table(d); diagnostic <- diagnostics_table(d,selected)
  write_csv(d,file.path(run,"forecasts_and_intervals.csv"));write_csv(s,file.path(run,"model_scores.csv"))
  write_csv(selected,file.path(run,"locked_model_choices.csv"));write_csv(years,file.path(run,"scores_by_year.csv"))
  write_csv(coverage_by_month(d),file.path(run,"scores_by_calendar_month.csv"))
  write_csv(diagnostic,file.path(run,"selected_error_diagnostics.csv"));write_csv(result$issues,file.path(run,"model_issues.csv"))
  misses <- d[(!is.na(d$covered95)&!d$covered95)|(!is.na(d$covered98)&!d$covered98),]
  write_csv(misses,file.path(run,"months_outside_intervals.csv"));write_csv(literature(),file.path(run,"literature.csv"))
  config <- list(script_version="1.0",start=START_DATE,first_target=FIRST_TARGET,common_start=COMMON_START,
    selection_end=SELECTION_END,end=END_DATE,levels=LEVELS,point_tolerance_pp=POINT_TOLERANCE_PP,
    seed=SEED,monthly_information_lag=MONTHLY_INFORMATION_LAG,models=catalogue,
    selection="Complete 2013-2016 WIS; then MAE; no holdout reselection",input_paths=lapply(inputs,`[[`,"original"))
  dput(config,file=file.path(run,"configuration.R"));writeLines(capture.output(sessionInfo()),file.path(run,"sessionInfo.txt"))
  saveRDS(list(config=config,series=series,forecasts=d,scores=s,selected=selected),file.path(run,"validation_results.rds"))
  plots <- make_plots(d,s,selected,run)
  write_report(run,d,s,selected,years,diagnostic,catalogue,result$issues,manifest,plots)
  writeLines("All requested analysis and exports completed.",file.path(run,"RUN_COMPLETE"))
  complete <- TRUE
  message("Complete: ",normalizePath(file.path(run,"report.html"),winslash="/"))
  invisible(list(run_dir=normalizePath(run,winslash="/"),scores=s,selected=selected))
}

if(isTRUE(RUN_ANALYSIS)) baseline_validation <- run_analysis()
