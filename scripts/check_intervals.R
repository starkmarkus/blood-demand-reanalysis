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

cat('Checking NNAR point forecasts against original PI=TRUE option\n');flush.console()
pi_checks <- lapply(seq_len(nrow(nnpars)),function(j) {
  set.seed(1)
  fit <- cw(nnetar(prefix(global,75),p=nnpars$p[j],P=0,lambda='auto',repeats=nnpars$repeats[j]))
  withpi <- cw(forecast(fit,PI=TRUE,h=18,level=c(95),interval='c'))
  without <- cw(forecast(fit,PI=FALSE,h=18))
  shorter <- cw(forecast(fit,PI=FALSE,h=6))
  data.frame(p=nnpars$p[j],repeats=nnpars$repeats[j],pi_delta=max(abs(withpi$mean-without$mean)),horizon_delta=max(abs(shorter$mean-head(without$mean,6))))
})
pi_checks <- do.call(rbind,pi_checks);stopifnot(max(pi_checks$pi_delta)<1e-10,max(pi_checks$horizon_delta)<1e-10)
write_json(pi_checks,file.path(root,'results/nnar-point-equivalence.json'),pretty=TRUE,digits=12)

