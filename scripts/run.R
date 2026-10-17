root <- normalizePath(commandArgs(TRUE)[1],winslash='/')
dir.create(file.path(root,'results'),showWarnings=FALSE)
dir.create(file.path(root,'logs'),showWarnings=FALSE)
suppressPackageStartupMessages({library(forecast);library(tsfknn);library(imputeTS);library(tsoutliers);library(jsonlite);library(zoo)})
source(file.path(root,'inputs/original/NNAR-tuning-script.R'))
source(file.path(root,'inputs/original/ro-eval-NNAR-script.R'))
raw <- ts(read.csv(file.path(root,'inputs/original/Blood_data_Tema.csv'))$QTY_DEMANDED,start=c(2013,1),frequency=12)
warnings <- character(); fallbacks <- list()
cw <- function(expr) withCallingHandlers(expr,warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
prefix <- function(x,n) ts(as.numeric(x)[seq_len(n)],start=c(2013,1),frequency=12)
mape <- function(pred,truth) mean(abs((as.numeric(pred)-as.numeric(truth))/as.numeric(truth)),na.rm=TRUE)*100
metrics <- function(pred,truth) c(MAPE=mape(pred,truth),MAE=mean(abs(pred-truth),na.rm=TRUE),RMSE=sqrt(mean((pred-truth)^2,na.rm=TRUE)))
prepare <- function(x,n,cache) {
  key <- as.character(n)
  if(!exists(key,cache,inherits=FALSE)) {
    imp <- cw(na_kalman(prefix(x,n),model='StructTS',smooth=TRUE,nit=-1))
    fit <- tryCatch(cw(tso(imp,types=c('AO','LS','TC'),maxit.iloop=10,tsmethod='auto.arima')),error=function(e) {
      fallbacks[[length(fallbacks)+1L]] <<- list(n=n,error=conditionMessage(e),action='ML with allowdrift=FALSE, ic=bic')
      cw(tso(imp,types=c('AO','LS','TC'),maxit.iloop=10,tsmethod='auto.arima',args.tsmethod=list(allowdrift=FALSE,ic='bic',method='ML')))
    })
    assign(key,fit$yadj,cache)
  }
  get(key,cache,inherits=FALSE)
}
clean_cache <- new.env(parent=emptyenv()); grid_cache <- new.env(parent=emptyenv()); nn_cache <- new.env(parent=emptyenv())
global <- prepare(raw,length(raw),clean_cache)
stopifnot(all(!is.na(raw[76:93])),max(abs(raw[76:93]-global[76:93]))<1e-10)
nnpars <- expand.grid(repeats=c(25,30,35),p=1:3)[,c('p','repeats')]
kn <- function(train,h,lag,k) as.numeric(knn_forecasting(train,h=h,lags=lag,k=k,msas='MIMO',cf='mean',transform='none')$prediction)
get_grid <- function(x,n,h,mode,method,seed=1,cc=clean_cache,gc=grid_cache,nc=nn_cache,g=global) {
  key <- paste(mode,method,n,h,seed,sep='_')
  if (!exists(key,gc,inherits=FALSE)) {
    train <- if(mode=='G') prefix(g,n) else prepare(x,n,cc)
    if(method=='KNN') {
      pars <- do.call(rbind,lapply(1:6,function(lag) data.frame(lag=lag,k=seq_len(min(52,n-lag-h+1)))))
      pred <- vapply(seq_len(nrow(pars)),function(j)kn(train,h,pars$lag[j],pars$k[j]),numeric(h))
      ids <- paste(pars$lag,pars$k,sep=':')
    } else {
      pars <- nnpars; ids <- paste(pars$p,pars$repeats,sep=':')
      pred <- vapply(seq_len(nrow(pars)),function(j) {
        fitkey <- paste(mode,n,pars$p[j],pars$repeats[j],seed,sep='_')
        if(!exists(fitkey,nc,inherits=FALSE)) {
          set.seed(seed)
          fit <- cw(nnetar(train,p=pars$p[j],P=0,lambda='auto',repeats=pars$repeats[j]))
          pp <- as.numeric(cw(forecast(fit,h=18,PI=FALSE)$mean))
          assign(fitkey,pp,nc)
        }
        get(fitkey,nc,inherits=FALSE)[seq_len(h)]
      },numeric(h))
    }
    stopifnot(all(is.finite(pred)))
    assign(key,list(pars=pars,pred=pred,ids=ids),gc)
  }
  get(key,gc,inherits=FALSE)
}
choose <- function(x,n,h,H,mode,method,seed=1,cc=clean_cache,gc=grid_cache,nc=nn_cache,g=global) {
  ends <- n-H-c(12,6,0)
  stopifnot(max(ends+H)<=n,min(ends)>=24)
  folds <- lapply(ends,function(end) {
    z <- get_grid(x,end,H,mode,method,seed,cc,gc,nc,g)
    z$loss <- apply(z$pred,2,mape,truth=as.numeric(x)[end+seq_len(H)])
    z
  })
  outer <- get_grid(x,n,h,mode,method,seed,cc,gc,nc,g)
  ids <- Reduce(intersect,c(lapply(folds,function(z)z$ids),list(outer$ids)))
  losses <- vapply(folds,function(z)z$loss[match(ids,z$ids)],numeric(length(ids)))
  avg <- rowMeans(losses); stopifnot(all(is.finite(avg)))
  sel <- which.min(avg)
  pp <- outer$pred[,match(ids,outer$ids),drop=FALSE]
  list(pred=pp,ids=ids,selected=sel,validation=avg[sel],ends=ends)
}

pi_checks <- fromJSON(file.path(root,'results/nnar-point-equivalence.json'))
stopifnot(nrow(pi_checks)==9L,all(pi_checks$pi_delta==0),all(pi_checks$horizon_delta==0))
pi_changes <- 0L
replace_pi <- function(expr) {
 if(!is.call(expr)) return(expr)
 if(identical(expr[[1]],quote(forecast::forecast)) && 'PI' %in% names(expr)) {
  stopifnot(identical(expr[['PI']],TRUE));expr[['PI']] <- FALSE;pi_changes <<- pi_changes+1L
 }
 for(i in seq_along(expr)) if(i>1L && is.call(expr[[i]])) expr[[i]] <- replace_pi(expr[[i]])
 expr
}
for(nm in c('NNAR_tuning_forecast','ro_eval_NNAR_forecast')) {
 fn <- get(nm);body(fn) <- replace_pi(body(fn));assign(nm,fn)
}
stopifnot(pi_changes==2L)

rows <- list();predictions <- list();selections <- list();count <- 0L
add <- function(n,method,mode,selection,rule,seed,pred) {
  truth <- as.numeric(raw)[(n+1):length(raw)]
  count <<- count+1L
  rows[[count]] <<- data.frame(origin=n,horizon=length(truth),method=method,preprocessing=mode,selection=selection,rule=rule,seed=seed,t(metrics(as.numeric(pred),truth)))
  predictions[[count]] <<- data.frame(origin=n,target=n+seq_along(truth),method=method,preprocessing=mode,selection=selection,rule=rule,seed=seed,truth=truth,prediction=as.numeric(pred))
}
for(n in 75:91) {
  h <- length(raw)-n
  for(method in c('KNN','NNAR')) for(seed in if(method=='NNAR') 1:5 else 1) {
    for(rule in if(seed==1) c('matched','6','12','18') else 'matched') {
      H <- if(rule=='matched') h else as.integer(rule)
      for(mode in c('G','P')) {
        fit <- choose(raw,n,h,H,mode,method,seed)
        truth <- as.numeric(raw)[(n+1):length(raw)]
        errors <- apply(fit$pred,2,mape,truth=truth)
        oracle <- which.min(errors)
        stopifnot(errors[oracle]<=errors[fit$selected]+1e-10)
        add(n,method,mode,'O',rule,seed,fit$pred[,oracle])
        add(n,method,mode,'V',rule,seed,fit$pred[,fit$selected])
        selections[[length(selections)+1L]] <- data.frame(origin=n,method=method,mode=mode,rule=rule,seed=seed,inner_H=H,inner_ends=paste(fit$ends,collapse=';'),candidates=length(fit$ids),selected=fit$ids[fit$selected],oracle=fit$ids[oracle],validation_MAPE=fit$validation)
      }
    }
  }
  train <- prepare(raw,n,clean_cache)
  for(mode in c('G','P')) {
    tr <- if(mode=='G') prefix(global,n) else train
    add(n,'ARIMA',mode,'none','none',1,cw(forecast(auto.arima(tr,seasonal=FALSE),h=h)$mean))
  }
  add(n,'Mean12','P','none','none',1,rep(mean(tail(train,12)),h))
  add(n,'Naive','P','none','none',1,rep(tail(train,1),h))
  add(n,'SeasonalNaive','P','none','none',1,rep(tail(train,12),length.out=h))
  checkpoint <- list(errors=do.call(rbind,rows),predictions=do.call(rbind,predictions),selections=do.call(rbind,selections))
  saveRDS(checkpoint,file.path(root,'results/factorial-checkpoint.rds'))
  cat('Completed origin',n,'of 91;',count,'metric rows\n');flush.console()
}
out <- checkpoint
out$summary <- do.call(rbind,lapply(split(out$errors,interaction(out$errors[,c('method','preprocessing','selection','rule','seed')],drop=TRUE)),function(z) data.frame(z[1,c('method','preprocessing','selection','rule','seed')],median_MAPE=median(z$MAPE),mean_MAPE=mean(z$MAPE),median_MAE=median(z$MAE),median_RMSE=median(z$RMSE))))
rownames(out$summary)<-NULL
saveRDS(out,file.path(root,'results/factorial.rds'))
write_json(out$summary,file.path(root,'results/factorial-summary.json'),pretty=TRUE,digits=12)
write_json(out$errors,file.path(root,'results/factorial-errors.json'),pretty=TRUE,digits=12)
write_json(out$predictions,file.path(root,'results/factorial-predictions.json'),pretty=TRUE,digits=12)
write_json(out$selections,file.path(root,'results/factorial-selections.json'),pretty=TRUE,digits=12)
saveRDS(list(clean=as.list(clean_cache),nn=as.list(nn_cache)),file.path(root,'results/model-cache.rds'))
cat('Running original NNAR control flow with PI disabled; see README.md\n');flush.console()
original_nnar <- cw(ro_eval_NNAR_forecast(global,18,seed.number=1))
write_json(original_nnar$NNAR_errors,file.path(root,'results/original-nnar-point-only-runtime.json'),pretty=TRUE,digits=12)
saveRDS(original_nnar,file.path(root,'results/original-nnar-point-only-runtime.rds'))
recomputed <- subset(out$errors,method=='NNAR' & preprocessing=='G' & selection=='O' & rule=='matched' & seed==1)
nn_delta <- max(abs(recomputed$MAPE-original_nnar$NNAR_errors$MAPE))
stopifnot(nn_delta<1e-10)
cat('Checking future invariance using fresh caches\n');flush.console()
invariance <- list()
for(n in c(75,83,91)) for(method in c('KNN','NNAR')) {
  h <- length(raw)-n
  altered <- raw;altered[(n+1):length(raw)]<-raw[(n+1):length(raw)]*20+500
  a <- choose(raw,n,h,h,'P',method,1)
  b <- choose(altered,n,h,h,'P',method,1,new.env(parent=emptyenv()),new.env(parent=emptyenv()),new.env(parent=emptyenv()))
  delta <- max(abs(a$pred-b$pred)); same <- a$ids[a$selected]==b$ids[b$selected]
  stopifnot(delta==0,same)
  invariance[[length(invariance)+1L]]<-data.frame(origin=n,method=method,max_abs_delta=delta,same_selection=same)
}
write_json(list(nnar_original_point_only_max_difference=nn_delta, interval_simulation_omitted=TRUE, protocol_amendment="README.md",invariance=do.call(rbind,invariance),raw_outer_targets_unchanged=TRUE),file.path(root,'results/factorial-verification.json'),pretty=TRUE,digits=12)
write_json(list(warnings=unique(warnings),warning_count=length(warnings),fallbacks=fallbacks),file.path(root,'logs/factorial-diagnostics.json'),pretty=TRUE)
capture.output(sessionInfo(),file=file.path(root,'results/session-info.txt'))
cat('Complete\n');print(subset(out$summary,rule %in% c('matched','none') & seed==1))
