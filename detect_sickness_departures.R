# England and Wales: frozen-baseline monitoring, January 2019-December 2020
# Version 1.0 | R >= 4.3 | Outcome: monthly sickness absence rate, NOT mortality.
# Install once: install.packages(c("forecast", "ggplot2", "jsonlite"))
# Source this WHOLE script in RStudio after running the England and Wales scripts.
# Set ENGLAND_RUN and WALES_RUN below to exact export folders or ZIPs if needed.
# No download, Pandoc, R Markdown, or separate helper script is required.
# Outputs: timestamped run directory, self-contained HTML, CSVs, PNG/PDF, RDS.
# Allow several minutes for simulation/refitting. All random seeds are recorded.
#
# Prespecified primary model: K=2 Fourier seasonality + intercept + AR(1) errors.
# Fit separately to 2010-2018, then FREEZE through 2019-2020 (no outcome updates).
# Historical validation: 24-month paths from Dec 2012, 2013, 2014 and Dec 2016.
# Development paths overlap; Dec 2016 -> 2017-2018 is a separate audit path.
# Primary upward CUSUM threshold targets a 5% chance of ANY alarm across BOTH
# countries and ALL 24 months, under the fitted joint Gaussian AR(1) null.
# Bootstrap refitting approximates parameter-estimation uncertainty; no exact
# conditional false-alarm guarantee is claimed. Independent simulations audit it.
# A pointwise 98% crossing, a CUSUM alarm, and an estimated onset are different.
# Onset bootstrap range is exploratory, conditional on a single positive step
# and the selected alarm window, not a selection-adjusted confidence interval.
# Source data are current/revised vintages: this is NOT a real-time release replay.

ENGLAND_ROOT <- "england_sickness_results"
WALES_ROOT <- "wales_sickness_results"
ENGLAND_RUN <- getOption("monitor.england_run", NULL)
WALES_RUN <- getOption("monitor.wales_run", NULL)
OUT_DIR <- getOption("monitor.output_dir", "sickness_departure_results")
RUN_ANALYSIS <- getOption("monitor.run_analysis", TRUE)
START_DATE <- as.Date("2010-01-01")
TRAIN_END <- as.Date("2018-12-01")
END_DATE <- as.Date("2020-12-01")
HORIZON <- 24L
LEVELS <- c(50, 80, 90, 95, 98)
POINT_TOLERANCE_PP <- 0.25
SEED <- getOption("monitor.seed", 20261022L)
B_CALIBRATION <- getOption("monitor.bootstrap_calibration", 3000L)
B_VERIFICATION <- getOption("monitor.bootstrap_verification", 3000L)
B_HISTORICAL <- getOption("monitor.bootstrap_historical", 1000L)
B_HISTORICAL_VERIFY <- getOption("monitor.bootstrap_historical_verify", 1000L)
# Additional robustness check: paired circular blocks of fitted innovations.
# It is reported separately and does not change the prespecified primary rule.
RUN_BLOCK_SENSITIVITY <- getOption("monitor.block_sensitivity", TRUE)
B_BLOCK <- getOption("monitor.bootstrap_block", 2000L)
BLOCK_LENGTH <- 12L
FAMILY_ALPHA <- 0.05
CUSUM_K <- 0.5
ONSET_MIN_POST <- 2L
# Changing these settings defines a different analysis; do so before inspecting
# 2019-2020. 95%/98% pointwise bands are descriptive, not repeated-test control.

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
  if (file.exists(file.path(p,"RUN_COMPLETE"))||file.exists(file.path(p,"INPUT_VERIFIED"))) return(TRUE)
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
  if(anyDuplicated(d$date)||!identical(as.character(d$date),as.character(seq(START_DATE,END_DATE,by="month")))) stop(country,": require every month January 2010-December 2020, exactly once")
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

PRIMARY_ID <- "harmonic_ar1_exp"
COUNTRIES <- c("England","Wales")
FOURIER_NAMES <- c("s1","c1","s2","c2")
primary_fit <- function(history,h=HORIZON,verify_formula=FALSE) {
  if(nrow(history)<36L) stop("Primary model needs at least 36 training months")
  if(any(diff(history$t)!=1L)||any(!is.finite(history$rate_pct))) stop("Invalid primary training series")
  xr <- as.matrix(features(history$t)[,FOURIER_NAMES])
  xn <- as.matrix(features(max(history$t)+seq_len(h))[,FOURIER_NAMES])
  fit <- forecast::Arima(ts(history$rate_pct,frequency=12),order=c(1,0,0),
    xreg=xr,include.mean=TRUE,method="ML")
  if(!is.null(fit$code)&&fit$code!=0L) stop("Primary AR(1) optimiser did not converge")
  co <- coef(fit); phi <- unname(co["ar1"])
  if(!is.finite(phi)||abs(phi)>=1||!is.finite(fit$sigma2)||fit$sigma2<=0) stop("Invalid stationary primary model")
  mu <- as.numeric(co["intercept"]+xr%*%co[FOURIER_NAMES])
  mu_future <- as.numeric(co["intercept"]+xn%*%co[FOURIER_NAMES])
  center <- mu_future+phi^seq_len(h)*(tail(history$rate_pct,1)-tail(mu,1))
  se <- sqrt(fit$sigma2*cumsum(phi^(2*(seq_len(h)-1L))))
  if(verify_formula) {
    fc <- forecast::forecast(fit,h=h,xreg=xn,level=95)
    ref_se <- as.numeric(fc$upper-fc$mean)/qnorm(.975)
    if(max(abs(center-as.numeric(fc$mean)))>1e-7||max(abs(se-ref_se))>1e-7)
      stop("Closed-form AR(1) forecasts disagree with forecast::forecast")
  }
  list(fit=fit,phi=phi,sigma=sqrt(fit$sigma2),mu=mu,mu_future=mu_future,
    mean=center,se=se,history=history,innovations=as.numeric(residuals(fit)))
}

fit_pair <- function(histories) {
  if(!identical(histories$England$date,histories$Wales$date)) stop("Training calendars differ")
  lapply(histories,primary_fit,verify_formula=TRUE)
}

bootstrap_pair <- function(fits,B,seed,label="bootstrap",method="parametric",block_length=BLOCK_LENGTH) {
  # Draw complete paired training + future paths under the fitted null, refit
  # each synthetic training series, then studentize its 24-step forecast errors.
  # No actual target outcomes are inputs to this function.
  set.seed(seed)
  n <- nrow(fits[[1]]$history); H <- length(fits[[1]]$mean)
  phi <- vapply(fits,`[[`,numeric(1),"phi")
  sig <- vapply(fits,`[[`,numeric(1),"sigma")
  innovations <- do.call(cbind,lapply(fits,function(z)z$innovations[-1L]))
  rho <- cor(innovations)[1,2]
  if(!is.finite(rho)||abs(rho)>=.999) stop("Singular cross-country innovation correlation")
  Sigma <- outer(sig,sig)*matrix(c(1,rho,rho,1),2)
  Gamma <- Sigma/(1-outer(phi,phi))
  chol_e <- chol(Sigma); chol_u <- chol(Gamma)
  if(!method%in%c("parametric","paired_block"))stop("Unknown bootstrap method")
  residual_pool <- sweep(innovations,2,colMeans(innovations),`-`)
  residual_pool <- sweep(residual_pool,2,apply(residual_pool,2,sd),`/`)
  residual_pool <- sweep(residual_pool,2,sig,`*`)
  paired_blocks <- function(N) {
    starts <- sample.int(nrow(residual_pool),ceiling(N/block_length),replace=TRUE)
    ix <- unlist(lapply(starts,function(s)((s-1L+seq_len(block_length)-1L)%%nrow(residual_pool))+1L),use.names=FALSE)
    residual_pool[ix[seq_len(N)],,drop=FALSE]
  }
  mu <- do.call(cbind,lapply(fits,function(z)c(z$mu,z$mu_future)))
  z <- array(NA_real_,c(B,H,2),dimnames=list(NULL,NULL,COUNTRIES))
  failures <- character(); b <- 0L; attempts <- 0L
  max_attempts <- B+max(20L,ceiling(.03*B))
  while(b<B && attempts<max_attempts) {
    attempts <- attempts+1L
    burn <- if(method=="paired_block")200L else 0L;N <- n+H+burn
    e <- if(method=="paired_block")paired_blocks(N) else matrix(rnorm(N*2),N,2)%*%chol_e
    u <- matrix(0,N,2)
    previous <- if(method=="paired_block")c(0,0) else as.numeric(matrix(rnorm(2),1)%*%chol_u)
    for(j in seq_len(N)) { u[j,] <- phi*previous+e[j,]; previous <- u[j,] }
    u <- u[burn+seq_len(n+H),,drop=FALSE]
    y <- mu+u
    draw <- tryCatch(withCallingHandlers({
      ans <- matrix(NA_real_,H,2)
      for(cc in seq_along(COUNTRIES)) {
        history <- fits[[cc]]$history; history$rate_pct <- y[seq_len(n),cc]
        bf <- primary_fit(history,H)
        ans[,cc] <- (y[n+seq_len(H),cc]-bf$mean)/bf$se
      }
      if(any(!is.finite(ans))) stop("Nonfinite bootstrap error")
      ans
    },warning=function(w)stop(conditionMessage(w))),error=function(e)e)
    if(inherits(draw,"error")) failures <- c(failures,conditionMessage(draw)) else {
      b <- b+1L; z[b,,] <- draw
    }
    if(b>0 && b%%500L==0L && !inherits(draw,"error")) message(label,": ",b,"/",B," draws")
  }
  if(b<B||length(failures)>max(5L,.02*attempts)) stop(label,": excessive bootstrap refit failures (",length(failures),"/",attempts,")")
  list(z=z,rho=rho,failures=failures,attempts=attempts,seed=seed,B=B,method=method,
    block_length=if(method=="paired_block")block_length else NA_integer_)
}

cusum_path <- function(z,k=CUSUM_K) {
  out <- numeric(length(z)); previous <- 0
  for(i in seq_along(z)) { previous <- max(0,previous+z[i]-k);out[i] <- previous }
  out
}
bootstrap_statistics <- function(z) {
  C <- array(0,dim(z),dimnames=dimnames(z))
  for(cc in seq_len(dim(z)[3])) for(b in seq_len(dim(z)[1])) C[b,,cc] <- cusum_path(z[b,,cc])
  list(cusum=C,max_cusum=apply(C,1,max),max_z=apply(z,1,max))
}
calibrate_monitor <- function(calibration,verification) {
  cs <- bootstrap_statistics(calibration$z); vs <- bootstrap_statistics(verification$z)
  # Conservative empirical order-statistic quantile; finite Monte Carlo
  # uncertainty is audited on a separate, independently seeded set of draws.
  limit <- unname(quantile(cs$max_cusum,1-FAMILY_ALPHA,type=1))
  upper <- unname(quantile(cs$max_z,1-FAMILY_ALPHA,type=1))
  audit <- bind_rows(lapply(c("CUSUM primary","Single-month upper envelope (alternative)"),function(rule) {
    alarm <- if(rule=="CUSUM primary") vs$max_cusum>limit else vs$max_z>upper
    bt <- binom.test(sum(alarm),length(alarm))
    data.frame(rule=rule,target_probability=FAMILY_ALPHA,observed_probability=mean(alarm),
      alarms=sum(alarm),simulated_paths=length(alarm),mc_lower95=bt$conf.int[1],mc_upper95=bt$conf.int[2],
      threshold=if(rule=="CUSUM primary")limit else upper)
  }))
  qs <- lapply(LEVELS,function(lev) {
    a <- (1-lev/100)/2
    list(lower=apply(calibration$z,c(2,3),quantile,probs=a,type=8),
      upper=apply(calibration$z,c(2,3),quantile,probs=1-a,type=8))
  });names(qs) <- as.character(LEVELS)
  list(cusum_limit=limit,upper_limit=upper,quantiles=qs,audit=audit)
}

apply_monitor <- function(fits,targets,calibration,origin,role) {
  bind_rows(lapply(seq_along(COUNTRIES),function(cc) {
    cn <- COUNTRIES[cc]; f <- fits[[cc]]; target <- targets[[cc]]
    stopifnot(nrow(target)==HORIZON,all(target$date>max(f$history$date)))
    z <- data.frame(country=cn,origin=as.Date(origin),role=role,date=target$date,h=seq_len(HORIZON),
      actual=target$rate_pct,prediction=f$mean,scale=f$se,model=PRIMARY_ID,
      train_start=min(f$history$date),train_end=max(f$history$date),n_train=nrow(f$history),pi_available=TRUE)
    for(lev in LEVELS) {
      q <- calibration$quantiles[[as.character(lev)]]
      z[[paste0("lower",lev)]] <- f$mean+f$se*q$lower[,cc]
      z[[paste0("upper",lev)]] <- f$mean+f$se*q$upper[,cc]
      z[[paste0("conditional_lower",lev)]] <- f$mean-qnorm((1+lev/100)/2)*f$se
      z[[paste0("conditional_upper",lev)]] <- f$mean+qnorm((1+lev/100)/2)*f$se
    }
    z$standardized_error <- (z$actual-z$prediction)/z$scale
    z$cusum <- cusum_path(z$standardized_error)
    z$cusum_limit <- calibration$cusum_limit
    z$cusum_above_limit <- z$cusum>z$cusum_limit
    z$joint_upper_envelope <- z$prediction+calibration$upper_limit*z$scale
    z$above_joint_envelope <- z$actual>z$joint_upper_envelope
    z$above98 <- z$actual>z$upper98; z$below98 <- z$actual<z$lower98
    z$outside98 <- z$above98|z$below98
    score_forecasts(z)
  }))
}
first_date <- function(d,flag) {
  ix <- which(!is.na(flag)&flag)
  if(length(ix))as.character(d[min(ix)]) else NA_character_
}
departure_table <- function(d) bind_rows(lapply(split(d,interaction(d$country,d$origin,drop=TRUE)),function(z) {
  z <- z[order(z$date),]
  data.frame(country=z$country[1],origin=z$origin[1],role=z$role[1],
    first_exit98=first_date(z$date,z$outside98),first_above98=first_date(z$date,z$above98),
    first_below98=first_date(z$date,z$below98),first_primary_cusum=first_date(z$date,z$cusum_above_limit),
    first_joint_upper_envelope=first_date(z$date,z$above_joint_envelope),
    months_above98=sum(z$above98),months_below98=sum(z$below98),
    primary_cusum_limit=z$cusum_limit[1])
}))

threshold_mc_sensitivity <- function(d,bootstrap,label) {
  # Binomial order-statistic interval for the simulated null's 95th quantile.
  # This is Monte Carlo precision, NOT a confidence interval for true onset.
  maxima <- sort(bootstrap_statistics(bootstrap$z)$max_cusum);B <- length(maxima)
  ranks <- pmax(1L,pmin(B,qbinom(c(.025,.975),B,1-FAMILY_ALPHA)+c(0L,1L)))
  bounds <- maxima[ranks]
  bind_rows(lapply(COUNTRIES,function(cn) {
    z <- d[d$country==cn,];z <- z[order(z$date),]
    data.frame(country=cn,calibration=label,draws=B,threshold=z$cusum_limit[1],
      threshold_mc_lower95=bounds[1],threshold_mc_upper95=bounds[2],
      estimated_alarm=first_date(z$date,z$cusum_above_limit),
      earliest_mc95=first_date(z$date,z$cusum>bounds[1]),
      latest_mc95=first_date(z$date,z$cusum>bounds[2]))
  }))
}

gls_step <- function(errors,precision,min_post=ONSET_MIN_POST) {
  n <- length(errors)
  if(n<min_post) return(NULL)
  candidates <- seq_len(n-min_post+1L)
  G <- outer(seq_len(n),candidates,`>=`)*1
  numer <- as.numeric(crossprod(G,precision%*%errors))
  denom <- colSums(G*(precision%*%G))
  delta <- pmax(numer/denom,0); improvement <- pmax(numer,0)^2/denom
  best <- which.max(improvement)
  list(onset=candidates[best],delta=delta[best],improvement=improvement[best],
    candidates=data.frame(onset_h=candidates,step_pp=delta,gls_improvement=improvement))
}
estimate_onsets <- function(monitored,bootstrap,scope="through_first_alarm") {
  # Exploratory dating, separate from the calibrated prospective alarm.
  # Noise paths come from pre-2019-only simulations. The fitted positive step
  # is estimated retrospectively; no step is inserted into the forecast model.
  summaries <- profiles <- distributions <- list()
  for(cc in seq_along(COUNTRIES)) {
    cn <- COUNTRIES[cc]; z <- monitored[monitored$country==cn,];z <- z[order(z$date),]
    alarm <- which(z$cusum_above_limit)
    n <- if(scope=="through_first_alarm") { if(length(alarm))min(alarm) else 0L } else nrow(z)
    base <- data.frame(country=cn,scope=scope,window_end=if(n)as.character(z$date[n]) else NA_character_,
      n_months=n,onset=NA_character_,range_lower95=NA_character_,range_upper95=NA_character_,
      step_pp=NA_real_,left_boundary=NA,zero_step_bootstrap_pct=NA_real_,
      boundary_bootstrap_pct=NA_real_,status="No primary alarm; onset not estimated")
    if(n>=ONSET_MIN_POST) {
      E <- sweep(bootstrap$z[,seq_len(n),cc,drop=FALSE][,,1],2,z$scale[seq_len(n)],`*`)
      V <- cov(E); ridge <- 1e-8*max(diag(V)); P <- solve(V+diag(ridge,n))
      fit <- gls_step(z$error[seq_len(n)],P)
      if(fit$delta<=0) base$status <- "Positive step not supported in this window" else {
        step <- fit$delta*(seq_len(n)>=fit$onset)
        # Centre simulated errors; conditional uncertainty under the fitted step.
        noise <- sweep(E,2,colMeans(E),`-`)
        boot <- t(vapply(seq_len(nrow(noise)),function(b) {
          f <- gls_step(noise[b,]+step,P);c(onset_h=f$onset,step_pp=f$delta)
        },numeric(2)))
        positive <- boot[,"step_pp"]>0
        bounds <- if(any(positive))as.integer(quantile(boot[positive,"onset_h"],c(.025,.975),type=1)) else c(NA_integer_,NA_integer_)
        base$onset <- as.character(z$date[fit$onset]);base$range_lower95 <- as.character(z$date[bounds[1]])
        base$range_upper95 <- as.character(z$date[bounds[2]]);base$step_pp <- fit$delta
        base$left_boundary <- fit$onset==1L;base$zero_step_bootstrap_pct <- 100*mean(!positive)
        base$boundary_bootstrap_pct <- 100*mean(boot[positive,"onset_h"]==1L)
        base$status <- "Exploratory single-step estimate; range conditional on model and fixed window"
        pr <- fit$candidates;pr$country <- cn;pr$scope <- scope;pr$date <- z$date[pr$onset_h]
        profiles[[length(profiles)+1L]] <- pr
        bd <- as.data.frame(boot);bd$country <- cn;bd$scope <- scope;bd$draw <- seq_len(nrow(bd));bd$date <- z$date[bd$onset_h]
        distributions[[length(distributions)+1L]] <- bd
      }
    } else if(n>0) base$status <- "Alarm too early to estimate a step with two post-onset months"
    summaries[[length(summaries)+1L]] <- base
  }
  list(summary=bind_rows(summaries),profiles=bind_rows(profiles),draws=bind_rows(distributions))
}

benchmark_paths <- function(series,catalogue,folds) {
  out <- issues <- list()
  for(i in seq_len(nrow(folds))) for(cn in COUNTRIES) {
    origin <- folds$origin[i];d <- series[[cn]];history <- d[d$date<=origin,]
    target <- head(d[d$date>origin,],HORIZON)
    stopifnot(nrow(target)==HORIZON,all(target$date<=TRAIN_END))
    for(m in seq_len(nrow(catalogue))) {
      spec <- catalogue[m,];warns <- character()
      pred <- tryCatch(withCallingHandlers(predict_model(history,spec,HORIZON),warning=function(w) {
        warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")
      }),error=function(e)e)
      if(inherits(pred,"error")||length(warns)) issues[[length(issues)+1L]] <- data.frame(country=cn,origin=origin,model=spec$model,
        severity=if(inherits(pred,"error"))"error" else "warning",message=if(inherits(pred,"error"))conditionMessage(pred) else paste(unique(warns),collapse="; "))
      if(is.null(pred)||inherits(pred,"error"))next
      pred$country <- cn;pred$origin <- origin;pred$role <- folds$role[i];pred$model <- spec$model
      pred$date <- target$date;pred$actual <- target$rate_pct
      out[[length(out)+1L]] <- score_forecasts(pred)
    }
    message("24-month benchmark: ",cn," from ",origin)
  }
  list(forecasts=bind_rows(out),issues=if(length(issues))bind_rows(issues) else
    data.frame(country=character(),origin=as.Date(character()),model=character(),severity=character(),message=character()))
}
benchmark_summary <- function(d,catalogue) {
  out <- list()
  for(role in unique(d$role))for(cn in c(COUNTRIES,"Pooled"))for(id in catalogue$model) {
    z <- d[d$role==role & d$model==id,];if(cn!="Pooled")z <- z[z$country==cn,]
    if(!nrow(z))next
    expected <- HORIZON*if(role=="development")3L else 1L
    if(cn=="Pooled")expected <- expected*2L
    s <- summarise_rows(z,expected,role);s$country <- cn;s$model <- id
    s$label <- catalogue$label[match(id,catalogue$model)]
    s$unique_months_per_country <- length(unique(z$date));out[[length(out)+1L]] <- s
  }
  ans <- bind_rows(out);ans$rank <- NA_integer_
  for(role in unique(ans$period))for(cn in unique(ans$country)) {
    ix <- which(ans$period==role & ans$country==cn & ans$complete_intervals)
    ans$rank[ix] <- rank(ans$WIS[ix],ties.method="min")
  }
  ans
}

primary_scores <- function(d) bind_rows(lapply(split(d,interaction(d$country,d$origin,drop=TRUE)),function(z) {
  s <- summarise_rows(z,HORIZON,z$role[1]);s$country <- z$country[1];s$origin <- z$origin[1]
  s$first_target <- min(z$date);s$last_target <- max(z$date)
  s$cusum_alarm <- any(z$cusum_above_limit);s$upper_envelope_alarm <- any(z$above_joint_envelope);s
}))

fit_diagnostics <- function(fits,origin) bind_rows(lapply(names(fits),function(cn) {
  f <- fits[[cn]];r <- f$innovations[-1L]
  data.frame(country=cn,origin=as.Date(origin),n_train=nrow(f$history),ar1=f$phi,
    innovation_sd_pp=f$sigma,innovation_acf1=as.numeric(acf(r,plot=FALSE,lag.max=1)$acf[2]),
    ljung_box12_p=Box.test(r,lag=12,type="Ljung-Box",fitdf=1)$p.value,
    shapiro_p=shapiro.test(r)$p.value,
    note="Descriptive model checks; bootstrap validity depends on Gaussian stationary AR(1) assumptions")
}))

make_monitor_plots <- function(series,final,historical,scores,departures,run) {
  files <- list();cols <- c("Forecast"="#176c9e","Observed"="#1d2937")
  actual <- bind_rows(lapply(series,function(d)d[d$date>=as.Date("2016-01-01"),]))
  marks <- departures[!is.na(departures$first_primary_cusum),];marks$alarm_date <- as.Date(marks$first_primary_cusum)
  p <- ggplot2::ggplot()+
    ggplot2::geom_ribbon(data=final,ggplot2::aes(date,ymin=lower98,ymax=upper98),fill="#c9dcec",alpha=.8)+
    ggplot2::geom_ribbon(data=final,ggplot2::aes(date,ymin=lower95,ymax=upper95),fill="#91b9d6",alpha=.8)+
    ggplot2::geom_line(data=actual,ggplot2::aes(date,rate_pct,color="Observed"),linewidth=.5)+
    ggplot2::geom_line(data=final,ggplot2::aes(date,prediction,color="Forecast"),linewidth=.6)+
    ggplot2::geom_point(data=final[final$above98,],ggplot2::aes(date,actual),color="#b96d0a",size=1.8)+
    ggplot2::geom_vline(xintercept=as.numeric(as.Date("2019-01-01")),linetype="dashed",color="#627488")+
    ggplot2::geom_vline(data=marks,ggplot2::aes(xintercept=as.numeric(alarm_date)),color="#b62d36",linetype="dotted")+
    ggplot2::facet_wrap(~country,ncol=1,scales="fixed")+ggplot2::scale_color_manual(values=cols)+
    ggplot2::scale_x_date(date_breaks="1 year",date_labels="%Y")+
    ggplot2::labs(title="England and Wales | the frozen 2019-2020 projection",
      subtitle="Separate models fitted through December 2018; no refitting or updating with 2019-2020 outcomes",
      x=NULL,y="Sickness absence rate (%)",color=NULL,
      caption="Dark/light bands: bootstrap 95%/98% pointwise prediction intervals. Amber dots: above 98%. Red dotted line: first primary CUSUM alarm.")+
    plot_theme()
  files$projection <- save_figure(p,"frozen_projection_2019_2020",run,13,8)
  p <- ggplot2::ggplot(final,ggplot2::aes(date,cusum))+
    ggplot2::geom_line(color="#176c9e",linewidth=.75)+
    ggplot2::geom_line(ggplot2::aes(y=cusum_limit),color="#b62d36",linetype="dashed")+
    ggplot2::geom_vline(data=marks,ggplot2::aes(xintercept=as.numeric(alarm_date)),color="#b62d36",linetype="dotted")+
    ggplot2::facet_wrap(~country,ncol=1,scales="fixed")+
    ggplot2::scale_x_date(breaks=seq(as.Date("2019-01-01"),END_DATE,by="3 months"),date_labels="%b %Y")+
    ggplot2::labs(title="Accumulated upward deviations from the frozen baseline",
      subtitle="Primary CUSUM: one threshold calibrated jointly over 24 months and both countries",
      x=NULL,y="CUSUM statistic",
      caption="Dashed red: simulated 5% family false-alarm threshold. A large isolated observation can also trigger CUSUM; an alarm is not proof of persistence.")+
    plot_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30,hjust=1))
  files$cusum <- save_figure(p,"cusum_2019_2020",run,13,7)
  hold <- historical[historical$role=="audit_2017_2018",]
  p <- ggplot2::ggplot(hold,ggplot2::aes(date,actual))+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower98,ymax=upper98),fill="#c9dcec")+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower95,ymax=upper95),fill="#91b9d6")+
    ggplot2::geom_line(ggplot2::aes(color="Observed"),linewidth=.6)+
    ggplot2::geom_line(ggplot2::aes(y=prediction,color="Forecast"),linewidth=.6)+
    ggplot2::geom_point(data=hold[hold$outside98,],color="#b96d0a",size=2)+
    ggplot2::facet_wrap(~country,ncol=1,scales="fixed")+ggplot2::scale_color_manual(values=cols)+
    ggplot2::scale_x_date(breaks=seq(as.Date("2017-01-01"),TRAIN_END,by="3 months"),date_labels="%b %Y")+
    ggplot2::labs(title="Earlier audit | a 24-month frozen forecast",
      subtitle="Fit through December 2016, predict January 2017-December 2018 without a December 2017 refit",
      x=NULL,y="Sickness absence rate (%)",color=NULL,
      caption="Bootstrap 95%/98% pointwise prediction intervals. Amber dots: outside 98%. This is not the earlier annually refitted validation.")+
    plot_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30,hjust=1))
  files$audit <- save_figure(p,"historical_24_month_audit",run,13,8)
  rank <- scores[scores$country=="Pooled" & scores$period=="development" & scores$complete_intervals,]
  rank <- rank[order(rank$WIS),];rank <- head(rank,10)
  rank$label <- factor(rank$label,levels=rev(rank$label))
  p <- ggplot2::ggplot(rank,ggplot2::aes(WIS,label,fill=model==PRIMARY_ID))+
    ggplot2::geom_col(width=.65)+ggplot2::scale_fill_manual(values=c("TRUE"="#176c9e","FALSE"="#aab9c8"),guide="none")+
    ggplot2::labs(title="24-month development benchmark | ten lowest interval scores",
      subtitle="Equal country weights; overlapping paths from December 2012, 2013 and 2014. Blue: prespecified primary model.",
      x="Mean weighted interval score (lower is better)",y=NULL,
      caption="Benchmark candidates use their native intervals. Primary monitoring bands additionally use refit bootstrap calibration.")+
    plot_theme()
  files$ranking <- save_figure(p,"benchmark_24_month_ranking",run,13,6)
  files
}

report_html <- function(run,final,historical,departures,onsets,scores,primary,audits,diagnostics,quality,manifest,plots,issues,definition_notes,block_sensitivity=NULL,threshold_sensitivity=NULL) {
  image_tag <- function(path,alt) paste0('<figure><img alt="',escape_html(alt),'" src="data:image/png;base64,',
    jsonlite::base64_enc(readBin(path,"raw",n=file.info(path)$size)),'"></figure>')
  date_label <- function(x) ifelse(is.na(x),"No crossing / alarm",format(as.Date(x),"%B %Y"))
  dt <- departures[,c("country","first_exit98","first_above98","first_primary_cusum","first_joint_upper_envelope")]
  for(j in 2:ncol(dt))dt[[j]] <- date_label(dt[[j]])
  names(dt) <- c("Country","First outside 98% (either side)","First above 98%","First primary CUSUM alarm","First joint upper-envelope crossing")
  robustness_top <- ""
  if(!is.null(block_sensitivity)) {
    bd <- departure_table(block_sensitivity$monitoring)
    bi <- match(departures$country,bd$country)
    dt[[5]] <- date_label(bd$first_primary_cusum[bi]);names(dt)[5] <- "Block sensitivity: first CUSUM alarm"
    shared <- pmax(as.Date(departures$first_primary_cusum),as.Date(bd$first_primary_cusum[bi]))
    shared_text <- paste(paste0(departures$country,": ",date_label(as.character(shared))),collapse="; ")
    robustness_top <- paste0('<div class="callout caution"><strong>Check date sensitivity.</strong> A paired-block resampling check is reported alongside the Gaussian AR(1) calibration. First month by which both calibrations have signalled: ',
      escape_html(shared_text),', using the estimated thresholds. This is agreement between two analyses, not proof of the true onset month. Earlier dates depend on calibration assumptions; finite-simulation ranges are shown below.</div>')
  }
  mc_html <- ""
  if(!is.null(threshold_sensitivity)) {
    mt <- threshold_sensitivity[,c("country","calibration","estimated_alarm","earliest_mc95","latest_mc95")]
    for(j in 3:5)mt[[j]] <- date_label(mt[[j]])
    names(mt) <- c("Country","Calibration","Alarm at estimated threshold","MC range: earliest","MC range: latest")
    mc_html <- paste0('<h3>How precise is the simulated alarm month?</h3>',html_table(mt),
      '<p>Each threshold is estimated from a finite simulation sample. A binomial order-statistic 95% interval for that threshold is translated into the earliest and latest crossing on the observed path. These ranges describe <strong>Monte Carlo precision only</strong>, conditional on the fitted model and data. They are not confidence intervals for the true change date and do not cover model misspecification.</p>')
  }
  onset <- onsets[onsets$scope=="through_first_alarm",c("country","onset","range_lower95","range_upper95","step_pp","left_boundary","status")]
  names(onset) <- c("Country","Estimated step onset","Conditional range start","Conditional range end","Step (pp)","Left boundary?","Interpretation")
  historical_dates <- departure_table(historical)
  ah <- primary[primary$period=="audit_2017_2018",]
  dev <- scores[scores$country=="Pooled" & scores$period=="development" & scores$complete_intervals,]
  dev <- dev[order(dev$WIS),];best <- dev[1,];prim <- dev[dev$model==PRIMARY_ID,]
  winner_text <- if(best$model==PRIMARY_ID) "The prespecified model has the lowest pooled development WIS among the complete 24-month candidates. This does not establish that it is optimal for detecting changes." else
    paste0("The prespecified model is not the lowest-WIS candidate for 24-month forecasting: ",escape_html(best$label),
      " ranks first in the pooled development comparison. The monitoring model remains fixed as agreed; do not describe it as the empirically best 24-month predictor.")
  audit_warning <- if(any(ah$cusum_alarm))
    "The primary rule also alarms in the 2017-2018 audit. Earlier alarms may reflect real historical departures or model misspecification; the empirical audit does not support treating every alarm as a new, unusual 2019-2020 event." else
    "Neither country has a primary CUSUM alarm in the 2017-2018 audit. One earlier two-year audit is still too little to establish a reliable long-run false-alarm rate."
  earlier <- historical_dates[historical_dates$role=="development" & !is.na(historical_dates$first_primary_cusum),]
  development_warning <- if(nrow(earlier))paste0(" Earlier development paths do alarm: ",
    paste(paste0(earlier$country," in ",format(as.Date(earlier$first_primary_cusum),"%B %Y"),
      " (fit through ",format(earlier$origin,"%B %Y"),")"),collapse="; "),
    ". These overlapping historical periods are not known to be unchanged processes; those alarms cannot simply be counted as confirmed false positives.") else
    " No primary alarms occur in the three earlier development paths; these overlapping paths do not establish a real-world false-alarm rate."
  summary_rows <- ah[,c("country","MAE","RMSE","WIS","covered_n95","covered_n98","width95","width98","cusum_alarm")]
  names(summary_rows) <- c("Country","MAE (pp)","RMSE (pp)","WIS","Covered /24 at 95%","Covered /24 at 98%","Mean 95% width (pp)","Mean 98% width (pp)","CUSUM alarm?")
  bt <- scores[scores$country=="Pooled",c("period","rank","label","n_forecasts","n_intervals","MAE","WIS","coverage95","coverage98")]
  bt <- bt[order(bt$period,bt$rank,na.last=TRUE),]
  has_notes <- nrow(definition_notes)>0
  final_diag <- diagnostics[diagnostics$origin==TRAIN_END,]
  diagnostic_note <- if(any(final_diag$ljung_box12_p<.05|final_diag$shapiro_p<.05))
    '<div class="callout caution">At least one final-model residual diagnostic has p &lt; 0.05. Remaining dependence or departures from Gaussian errors can weaken the simulation calibration. These checks are descriptive, and the simulated false-alarm target should not be read as a verified real-data probability.</div>' else
    '<p>No final-model residual diagnostic listed below has p &lt; 0.05. This does not prove the Gaussian stationary model is correct.</p>'
  sensitivity_html <- if(is.null(block_sensitivity))'<p>The paired-block sensitivity check was disabled in this run.</p>' else {
    sd <- departure_table(block_sensitivity$monitoring)[,c("country","first_above98","first_primary_cusum")]
    for(j in 2:ncol(sd))sd[[j]] <- date_label(sd[[j]])
    names(sd) <- c("Country","Block sensitivity: first above 98%","Block sensitivity: first CUSUM alarm")
    paste0('<h3>Sensitivity to the error assumptions</h3><p>As a separate check, resample paired, centred fitted innovations in circular blocks of ',
      BLOCK_LENGTH,' months, preserve within-block and cross-country dependence, and refit both models on each simulated training sample. A 200-month burn-in initializes the synthetic process. Both the ',
      B_BLOCK,' calibration and ',B_BLOCK,' independent verification draws use only 2010–2018 residuals. The model form and CUSUM reference value remain unchanged.</p>',
      html_table(sd),'<p>This check allows some residual serial structure and non-Gaussian tails beyond the parametric AR(1) null. Its reliability is limited by nine years of data and the chosen block length; it is not an assumption-free correction. It is reported separately, not selected according to which 2019–2020 date it produces. Sensitivity dates are model-dependent and are not combined with the primary rule.</p>',
      '<details><summary>Independent paired-block simulation audit</summary>',html_table(block_sensitivity$calibration$audit),'</details>')
  }
  body <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>England &amp; Wales — frozen sickness-rate monitoring</title><style>',
    'body{font:16px/1.58 system-ui,-apple-system,Segoe UI,sans-serif;color:#203244;background:#f4f6f8;margin:0}main{max-width:1160px;margin:32px auto;background:white;padding:40px 48px;border-radius:12px}h1{font-size:36px;line-height:1.15;margin:8px 0 18px}h2{font-size:25px;margin-top:44px;color:#143c59}h3{margin-top:28px}p{max-width:100ch}.eyebrow{letter-spacing:.12em;text-transform:uppercase;color:#176c9e;font-size:12px;font-weight:700}.lead{font-size:19px}.callout{background:#edf4f9;border-left:4px solid #176c9e;padding:16px 20px;margin:24px 0}.caution{background:#fff6e8;border-color:#b96d0a}.table-wrap{overflow:auto;margin:18px 0}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:9px 11px;text-align:left;border-bottom:1px solid #dce4eb;vertical-align:top}th{background:#edf2f6;font-weight:650}tr:nth-child(even){background:#f8fafb}figure{margin:24px 0}img{width:100%;height:auto}a{color:#176c9e}code{background:#edf2f6;padding:2px 5px;border-radius:3px}details{margin:20px 0}summary{cursor:pointer;font-weight:650}small,.muted{color:#5a6c7e}@media(max-width:700px){main{margin:0;padding:24px 18px;border-radius:0}h1{font-size:29px}}@media print{main{margin:0;padding:0}details{display:block}table{font-size:10px}h2{break-after:avoid}figure{break-inside:avoid}}',
    '</style></head><body><main><div class="eyebrow">Monthly sickness absence · England and Wales</div>',
    '<h1>When do observations leave the frozen projection?</h1><p class="lead">A common seasonal model, fitted separately in each country through December 2018, projects January 2019–December 2020. The report separates individual unusual months, accumulated upward deviations, and retrospective estimates of onset.</p>',
    '<p class="muted">Generated ',escape_html(format(Sys.time(),tz="Europe/Paris",usetz=TRUE)),
    '. Outcome: percentage of available full-time-equivalent calendar days lost to sickness. All differences and interval widths are percentage points (pp).</p>',
    '<div class="callout"><strong>Primary decision rule:</strong> an upward CUSUM with its threshold fixed before examining 2019–2020. Under the fitted joint null, calibration targets a 5% probability of at least one alarm across both countries and all 24 months. This is a model-based target, not a guaranteed real-world error rate.</div>',
    '<h2>Departure dates</h2>',html_table(dt),robustness_top,mc_html,
    '<p><strong>A prediction-band crossing is not a breakpoint.</strong> The 98% band applies separately to each month, on both sides. The primary CUSUM accumulates positive standardized forecast errors. Its alarm month is when evidence crosses the threshold; it is not automatically the month the process changed. The joint upper envelope is a separately calibrated alternative for a single large month, not a second component of the primary rule. Combining either rule into one decision would require new calibration.</p>',
    image_tag(plots$projection,"Frozen 2019-2020 projection and observations for England and Wales"),
    image_tag(plots$cusum,"CUSUM paths and primary threshold"),
    '<h2>Estimated onset: a separate, exploratory calculation</h2>',html_table(onset),
    '<p>For each country with an alarm, fit one positive, persistent level step to the forecast deviations from January 2019 through that country’s first CUSUM alarm. Candidate onsets are searched independently, with at least two post-onset months. Generalized least squares uses the bootstrap forecast-error covariance. The range contains the middle 95% of onset estimates from simulations around that fitted step, conditional on a positive fitted step.</p>',
    '<div class="callout caution">These are <strong>conditional bootstrap ranges, not selection-adjusted 95% confidence intervals</strong>. They do not account for choosing the alarm window after observing the data, uncertainty over the shape or number of changes, or data-definition revisions. A January 2019 boundary estimate means the change may have started before monitoring began. A large isolated shock can trigger CUSUM and make a persistent-step model a poor description.</div>',
    '<p>The CSVs also show a full-2019–2020 single-step sensitivity analysis. That analysis uses later observations and may be dominated by a later large shock; it should not replace the first-alarm date or be described as prospectively detected onset. No dates are imposed from England onto Wales.</p>',
    '<h2>Does the baseline survive a true two-year forecast?</h2>',
    '<p>Three development origins are used: December 2012 → 2013–2014, December 2013 → 2014–2015, and December 2014 → 2015–2016. These yield 72 forecast instances per country but only 48 unique target months. They overlap and are not independent. The later audit fits once through December 2016 and predicts all of 2017–2018; it does not refit after 2017. No development path includes 2017 or 2018.</p>',
    '<p>The primary model was prespecified from the earlier analysis. The audit is separated computationally, but those years were already examined in previous work; it is not a newly untouched prospective test. The final 2010–2018 refit is allowed only after the monitoring design is fixed.</p>',
    html_table(summary_rows),'<div class="callout caution">',audit_warning,development_warning,'</div>',
    image_tag(plots$audit,"Two-year forecast audit for 2017-2018"),
    '<h3>Comparison with the previous 21 candidate baselines</h3><p>',winner_text,
    ' Ranking uses weighted interval score (WIS), which rewards narrow, well-calibrated predictive distributions; lower is better. A model does not win simply by using a wider interval. Countries receive equal weight. Candidates with incomplete interval predictions are not ranked. Benchmark candidates use their native interval methods; the primary monitoring and audit tables use additional bootstrap refitting.</p>',
    image_tag(plots$ranking,"Pooled development WIS ranking"),
    '<details><summary>All pooled benchmark scores and historical alarm dates</summary>',html_table(bt),html_table(historical_dates),'</details>',
    '<h2>Calibration and assumptions</h2>',
    '<p>The model is <code>rate(t) = intercept + annual Fourier terms + semiannual Fourier terms + AR(1) error</code>. It has no linear trend. Its structure is shared; coefficients and error variance are estimated separately. Seasonal effects, coefficients, and state are frozen at the end of 2018. Forecast uncertainty still changes with horizon. No observed 2019–2020 value updates the projections or thresholds.</p>',
    '<p>Paired Gaussian AR(1) training and future paths are simulated using the innovation correlation estimated from the two training series. Both models are refitted on each synthetic training sample. Studentized forecast errors generate pointwise bands and the maximum CUSUM distribution. This approximates parameter-estimation uncertainty and within- and cross-country dependence. It assumes stationary, homoscedastic Gaussian errors and stable seasonality; it does not cover model misspecification or every possible data-generating process.</p>',
    '<p><code>S(0)=0; S(h)=max[0, S(h−1)+(observed(h)−forecast(h))/conditional_SE(h)−0.5]</code>. The threshold is the empirical 95th percentile of the largest S across all 24 horizons and both countries. There is no reset after the first alarm and no rolling replacement of the baseline. The independent simulation audit below uses a different seed; its confidence limits describe Monte Carlo uncertainty only.</p>',
    html_table(audits[audits$role=="final_2019_2020",]),
    '<p>Final calibration draws: ',B_CALIBRATION,'; independent verification draws: ',B_VERIFICATION,
    '. Each historical origin uses ',B_HISTORICAL,' calibration and ',B_HISTORICAL_VERIFY,' independent verification draws. The seed schedule and any refit failures are exported. Nominal targets may differ from the independent Monte Carlo estimate and from empirical historical behavior.</p>',
    diagnostic_note,sensitivity_html,'<details><summary>Model diagnostics and historical simulation checks</summary>',html_table(diagnostics),html_table(audits[audits$role!="final_2019_2020",]),'</details>',
    '<h2>Data checks and interpretation</h2>',html_table(quality),
    '<p>The script requires every month from January 2010 through December 2020 exactly once for each all-staff national total, with finite values between 0 and 100. It stops on gaps or duplicate months and does not impute. Observations after 2020 are discarded before the analysis; actual 2019–2020 values are supplied only to the evaluation stage. Input files, their MD5 hashes, revision flags, quality notes, and available definition-change notes are retained.</p>',
    '<p>These are current/revised data, not archived releases as known at each historical date. Publication lags are not replayed. An alarm date is therefore the month of the observation, not its actual announcement date. Staffing composition, reporting practices, and definition changes can affect rates. Aggregation can also hide changes within staff groups.</p>',
    if(has_notes)paste0('<details><summary>Definition notes supplied with the country exports</summary>',html_table(definition_notes),'</details>') else '',
    '<p><strong>Interpretation:</strong> an excess relative to this counterfactual baseline does not establish its cause. A monitoring alarm is not evidence that COVID-19, a particular policy, or a specific workplace exposure caused the change. Sickness absence is not mortality.</p>',
    '<details><summary>Input provenance and model-fitting issues</summary>',html_table(manifest),html_table(issues),'</details>',
    '<h2>Reproduce and inspect</h2><p>Set <code>ENGLAND_RUN</code> and <code>WALES_RUN</code> to the corresponding country export folders or ZIP files, then source the entire script. Defaults locate the latest successful export under <code>england_sickness_results</code> and <code>wales_sickness_results</code>. Each successful run has <code>RUN_COMPLETE</code>, a self-contained <code>report.html</code>, a recorded configuration and R session, input copies, monthly forecasts, crossing dates, onset profiles and draws, validation scores, simulation checks, and PNG/PDF figures. No Pandoc or R Markdown is needed.</p>',
    '<h2>Method references</h2><ol>',
    '<li>Hyndman &amp; Athanasopoulos, <a href="https://otexts.com/fpp3/toolbox.html#time-series-cross-validation">Forecasting: Principles and Practice, time-series cross-validation</a>. Forecast horizons must match the intended use; later observations must not enter training.</li>',
    '<li>R Core, <a href="https://stat.ethz.ch/R-manual/R-devel/library/stats/html/arima.html">ARIMA regression</a> and <a href="https://stat.ethz.ch/R-manual/R-devel/library/stats/html/predict.arima.html">ARIMA prediction</a>. Native prediction standard errors omit parameter-estimation uncertainty, motivating the additional refit simulation.</li>',
    '<li>NIST/SEMATECH, <a href="https://www.itl.nist.gov/div898/handbook/pmc/section3/pmc323.htm">CUSUM control charts</a>. Accumulation helps detect smaller sustained departures. This script calibrates the finite-horizon statistic by dependent simulation instead of using independent-observation tabulated limits.</li>',
    '<li>Gandy &amp; Kvaløy (2013), <a href="https://doi.org/10.1002/sjos.12006">Guaranteed conditional performance of control charts via bootstrap methods</a>. Estimation uncertainty affects monitoring performance. This script uses a simpler refit-bootstrap approximation and does not implement or claim that paper’s conditional performance guarantee.</li>',
    '<li>Bracher et al. (2021), <a href="https://doi.org/10.1371/journal.pcbi.1008618">Evaluating epidemic forecasts in an interval format</a>. Weighted interval score evaluates calibration and sharpness jointly.</li>',
    '<li>R boot package, <a href="https://stat.ethz.ch/R-manual/R-devel/library/boot/html/tsboot.html">Bootstrapping of time series</a>. Documents block resampling and reconstructing a process from resampled residuals; the paired circular-block implementation in this script is written explicitly for the two-country sensitivity check.</li>',
    '</ol><p class="muted">This report and its figures are self-contained. Source URLs and files are preserved in the country exports. All model and alarm choices must be fixed before reusing the procedure on a new monitoring period.</p></main></body></html>')
  writeLines(body,file.path(run,"report.html"),useBytes=TRUE)
}

run_analysis <- function() {
  for(p in c("forecast","ggplot2","jsonlite"))if(!requireNamespace(p,quietly=TRUE))stop("Install missing package: ",p)
  Bs <- c(B_CALIBRATION,B_VERIFICATION,B_HISTORICAL,B_HISTORICAL_VERIFY,if(RUN_BLOCK_SENSITIVITY)B_BLOCK)
  if(any(!is.finite(Bs)|Bs<200L|Bs!=floor(Bs)))stop("Each bootstrap count must be an integer >=200; defaults are recommended")
  if(HORIZON!=24L||CUSUM_K<0||FAMILY_ALPHA<=0||FAMILY_ALPHA>=.5)stop("Invalid monitoring configuration")
  dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  run <- file.path(OUT_DIR,paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
  if(dir.exists(run))stop("Output run already exists; choose another output directory")
  dir.create(run);dir.create(file.path(run,"plots"));dir.create(file.path(run,"source"))
  success <- FALSE;on.exit(if(!success)writeLines("Run failed; do not treat partial outputs as complete.",file.path(run,"FAILED.txt")),add=TRUE)
  config_names <- c("ENGLAND_ROOT","WALES_ROOT","ENGLAND_RUN","WALES_RUN","OUT_DIR","START_DATE","TRAIN_END","END_DATE",
    "HORIZON","LEVELS","SEED","B_CALIBRATION","B_VERIFICATION","B_HISTORICAL","B_HISTORICAL_VERIFY","RUN_BLOCK_SENSITIVITY","B_BLOCK","BLOCK_LENGTH","FAMILY_ALPHA","CUSUM_K","ONSET_MIN_POST")
  dput(mget(config_names,envir=environment(run_analysis)),file.path(run,"configuration.R"))
  inputs <- list(England=resolve_run(ENGLAND_RUN,ENGLAND_ROOT,"England"),Wales=resolve_run(WALES_RUN,WALES_ROOT,"Wales"))
  on.exit(for(x in inputs)if(!is.null(x$cleanup))unlink(x$cleanup,recursive=TRUE),add=TRUE)
  series <- list();quality <- manifest <- definitions <- list()
  for(cn in COUNTRIES) {
    path <- file.path(inputs[[cn]]$path,obs_file);raw <- read.csv(path,stringsAsFactors=FALSE)
    series[[cn]] <- prepare_series(raw,cn)
    code <- if(cn=="England")"E" else "W"
    selected <- raw[raw$organisation_id==code & raw$staff_id=="XXX",]
    dates <- as.Date(as.character(selected$date));selected <- selected[dates>=START_DATE & dates<=END_DATE,]
    flag_count <- function(name)if(name%in%names(selected))sum(as.logical(selected[[name]]),na.rm=TRUE) else NA_integer_
    quality[[cn]] <- data.frame(country=cn,months=nrow(series[[cn]]),training_months=sum(series[[cn]]$date<=TRAIN_END),
      monitoring_months=sum(series[[cn]]$date>TRAIN_END),missing_rates=sum(!is.finite(series[[cn]]$rate_pct)),
      duplicate_months=anyDuplicated(series[[cn]]$date),revised_flags=flag_count("revised_flag"),provisional_flags=flag_count("provisional_flag"),
      vintage="Current export; original release vintages and publication delays not reconstructed")
    dest <- file.path(run,"source",cn);dir.create(dest)
    existing <- list.files(inputs[[cn]]$path,full.names=TRUE,recursive=FALSE)
    copy <- existing[basename(existing)%in%c(obs_file,"RUN_COMPLETE","configuration.R","sessionInfo.txt","READ_BEFORE_ANALYSIS.txt","source_status.csv",
      "known_definition_changes.csv","quality_checks.csv","revision_changes.csv","data_dictionary.csv","input_manifest.csv","source_manifest.csv")]
    if(any(!file.copy(copy,dest,overwrite=FALSE)))stop("Could not preserve source files for ",cn)
    writeLines("Validated observation snapshot for detect_sickness_departures.R; original country export was checked before copying.",file.path(dest,"INPUT_VERIFIED"))
    manifest[[cn]] <- data.frame(country=cn,input=inputs[[cn]]$original,observations_file=obs_file,
      md5=unname(tools::md5sum(path)),bytes=file.info(path)$size,copied_metadata_files=length(copy)-1L)
    notes <- file.path(inputs[[cn]]$path,"known_definition_changes.csv")
    if(file.exists(notes)) { nd <- read.csv(notes,stringsAsFactors=FALSE);nd$country <- cn;definitions[[cn]] <- nd }
  }
  quality <- bind_rows(quality);manifest <- bind_rows(manifest);definitions <- bind_rows(definitions)
  write_csv(bind_rows(series),file.path(run,"input_series_2010_2020.csv"));write_csv(quality,file.path(run,"data_quality.csv"))
  write_csv(manifest,file.path(run,"input_manifest.csv"));write_csv(definitions,file.path(run,"definition_notes.csv"))
  folds <- data.frame(origin=as.Date(c("2012-12-01","2013-12-01","2014-12-01","2016-12-01")),
    role=c(rep("development",3),"audit_2017_2018"))
  write_csv(folds,file.path(run,"historical_origins.csv"))
  catalogue <- model_catalogue();write_csv(catalogue,file.path(run,"model_catalogue.csv"))
  benchmark <- benchmark_paths(series,catalogue,folds)
  scores <- benchmark_summary(benchmark$forecasts,catalogue)
  write_csv(benchmark$forecasts,file.path(run,"benchmark_forecasts.csv"));write_csv(scores,file.path(run,"benchmark_scores.csv"))
  write_csv(benchmark$issues,file.path(run,"model_issues.csv"))
  historical <- audits <- diagnostics <- boot_logs <- seeds <- list()
  record_boot <- function(b,origin,role,purpose) {
    boot_logs[[length(boot_logs)+1L]] <<- data.frame(origin=origin,role=role,purpose=purpose,
      draws=b$B,attempts=b$attempts,failures=length(b$failures),messages=paste(unique(b$failures),collapse="; "),innovation_correlation=b$rho,
      method=b$method,block_length=b$block_length)
    seeds[[length(seeds)+1L]] <<- data.frame(origin=origin,role=role,purpose=purpose,seed=b$seed)
  }
  for(i in seq_len(nrow(folds)+1L)) {
    is_final <- i>nrow(folds)
    origin <- if(is_final)TRAIN_END else folds$origin[i]
    role <- if(is_final)"final_2019_2020" else folds$role[i]
    histories <- lapply(series,function(d)d[d$date<=origin,])
    fits <- fit_pair(histories)
    diagnostics[[i]] <- fit_diagnostics(fits,origin)
    message("Bootstrap monitoring from ",origin," (",role,")")
    b <- bootstrap_pair(fits,if(is_final)B_CALIBRATION else B_HISTORICAL,SEED+i*1000L,"calibration")
    v <- bootstrap_pair(fits,if(is_final)B_VERIFICATION else B_HISTORICAL_VERIFY,SEED+i*1000L+1L,"verification")
    record_boot(b,origin,role,"calibration");record_boot(v,origin,role,"verification")
    cal <- calibrate_monitor(b,v);a <- cal$audit;a$origin <- origin;a$role <- role;audits[[i]] <- a
    # Target observations first enter here, AFTER fit and calibration have ended.
    targets <- lapply(series,function(d)head(d[d$date>origin,],HORIZON))
    monitored <- apply_monitor(fits,targets,cal,origin,role)
    if(is_final) {
      final <- monitored;final_bootstrap <- b;final_verification <- v;final_fits <- fits;final_calibration <- cal
    } else historical[[i]] <- monitored
  }
  historical <- bind_rows(historical);audits <- bind_rows(audits);diagnostics <- bind_rows(diagnostics)
  block_sensitivity <- NULL
  if(isTRUE(RUN_BLOCK_SENSITIVITY)) {
    message("Final-model sensitivity: paired 12-month residual blocks")
    bb <- bootstrap_pair(final_fits,B_BLOCK,SEED+6000L,"block calibration",method="paired_block")
    bv <- bootstrap_pair(final_fits,B_BLOCK,SEED+6001L,"block verification",method="paired_block")
    record_boot(bb,TRAIN_END,"block_sensitivity","calibration");record_boot(bv,TRAIN_END,"block_sensitivity","verification")
    bc <- calibrate_monitor(bb,bv)
    bm <- apply_monitor(final_fits,lapply(series,function(d)d[d$date>TRAIN_END,]),bc,TRAIN_END,"block_sensitivity")
    block_sensitivity <- list(bootstrap=bb,verification=bv,calibration=bc,monitoring=bm)
    write_csv(bm,file.path(run,"block_sensitivity_monitoring.csv"));write_csv(departure_table(bm),file.path(run,"block_sensitivity_dates.csv"))
    write_csv(bc$audit,file.path(run,"block_sensitivity_verification.csv"))
  }
  threshold_sensitivity <- threshold_mc_sensitivity(final,final_bootstrap,"Parametric primary")
  if(!is.null(block_sensitivity))threshold_sensitivity <- rbind(threshold_sensitivity,
    threshold_mc_sensitivity(block_sensitivity$monitoring,block_sensitivity$bootstrap,"Paired-block sensitivity"))
  departures <- departure_table(final);primary <- primary_scores(historical)
  onset_alarm <- estimate_onsets(final,final_bootstrap,"through_first_alarm")
  onset_full <- estimate_onsets(final,final_bootstrap,"full_2019_2020")
  onsets <- rbind(onset_alarm$summary,onset_full$summary)
  exports <- list(monthly_monitoring=final,departure_dates=departures,onset_estimates=onsets,
    onset_profiles=bind_rows(list(onset_alarm$profiles,onset_full$profiles)),
    onset_bootstrap_draws=bind_rows(list(onset_alarm$draws,onset_full$draws)),
    historical_monitoring=historical,historical_departures=departure_table(historical),historical_primary_scores=primary,
    calibration_verification=audits,threshold_monte_carlo=threshold_sensitivity,model_diagnostics=diagnostics,bootstrap_log=bind_rows(boot_logs),random_seeds=bind_rows(seeds))
  for(nm in names(exports))write_csv(exports[[nm]],file.path(run,paste0(nm,".csv")))
  plots <- make_monitor_plots(series,final,historical,scores,departures,run)
  report_html(run,final,historical,departures,onsets,scores,primary,audits,diagnostics,quality,manifest,plots,benchmark$issues,definitions,block_sensitivity,threshold_sensitivity)
  saveRDS(list(fits=final_fits,calibration=final_calibration,bootstrap=final_bootstrap,verification=final_verification,
    final=final,historical=historical,onsets=onsets,benchmark_scores=scores,primary_scores=primary,block_sensitivity=block_sensitivity,
    threshold_sensitivity=threshold_sensitivity),file.path(run,"analysis_objects.rds"))
  writeLines(capture.output(sessionInfo()),file.path(run,"sessionInfo.txt"))
  writeLines(c("Common seasonal AR(1) form; separate country fits and onset searches.",
    "All final model training ends December 2018. No observed monitoring outcomes enter calibration.",
    "2017-2018 audit is a 24-month frozen path, not two refitted annual paths.",
    "Independent calibration/verification seeds are in random_seeds.csv.",
    "Closed-form primary forecasts checked against forecast::forecast at every historical and final fit.",
    "Bootstrap onset ranges are exploratory and not adjusted for alarm-window selection.",
    "RUN_COMPLETE denotes computation completed, not proof that the model assumptions hold."),file.path(run,"READ_ME.txt"))
  writeLines("Completed",file.path(run,"RUN_COMPLETE"));success <- TRUE
  writeLines(normalizePath(run,winslash="/"),file.path(OUT_DIR,"latest.txt"))
  message("Completed. Open: ",normalizePath(file.path(run,"report.html"),winslash="/"))
  invisible(run)
}

if(isTRUE(RUN_ANALYSIS)) run_analysis()
