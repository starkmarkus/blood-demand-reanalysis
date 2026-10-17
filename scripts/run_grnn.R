args <- commandArgs(TRUE)
root <- normalizePath(args[1],winslash='/')
exprs <- parse(file.path(root,'scripts/run.R'))
for(ex in exprs) {
 if(is.call(ex) && identical(ex[[1]],as.name('<-')) && identical(ex[[2]],as.name('pi_checks'))) break
 eval(ex,envir=.GlobalEnv)
}
if(length(args)>=2L && file.exists(args[2])) {
 old <- readRDS(args[2]);invisible(list2env(old$clean,envir=clean_cache))
}
for(f in list.files(file.path(root,'inputs/grnn/tsfgrnn/R'),pattern='[.]R$',full.names=TRUE))source(f)
stopifnot(identical(formals(grnn_forecasting)$msas,quote(c('recursive','MIMO'))))
grnn_base <- grnn_forecasting
grnn_cache <- new.env(parent=emptyenv())
grnn_forecasting <- function(timeS,h,lags=NULL,sigma=NULL,msas=c('recursive','MIMO'),scale=TRUE) {
 msas <- match.arg(msas)
 if(!is.null(sigma))return(grnn_base(timeS,h,lags,sigma,msas,scale))
 key <- paste(c(format(as.numeric(timeS),digits=17),h,lags,msas,scale),collapse='|')
 if(!exists(key,grnn_cache,inherits=FALSE))assign(key,cw(grnn_base(timeS,h,lags,sigma,msas,scale)),grnn_cache)
 get(key,grnn_cache,inherits=FALSE)
}
source(file.path(root,'inputs/original/GRNN-lag-best-script.R'))
source(file.path(root,'inputs/original/ro-eval-GRNN-script.R'))
cat('Running unchanged original GRNN forward control flow\n');flush.console()
original <- cw(ro_eval_GRNN_forecast(global,18))
write_json(original$GRNN_errors,file.path(root,'results/grnn-original-runtime.json'),pretty=TRUE,digits=12)
published <- fromJSON(file.path(root,'reference/grnn-published-errors.json'))
historical_delta <- max(abs(as.matrix(original$GRNN_errors)-as.matrix(published[,c('RMSE','MAE','MAPE')])))
cat('Original GRNN median MAPE',median(original$GRNN_errors$MAPE),'max rounded-table delta',historical_delta,'\n');flush.console()
grid <- function(x,n,h,mode,strategy,cc=clean_cache) {
 tr <- if(mode=='G')prefix(global,n) else prepare(x,n,cc)
 vapply(1:4,function(lag)as.numeric(grnn_forecasting(tr,h,lags=lag,msas=strategy,scale=FALSE,sigma=NULL)$prediction),numeric(h))
}
choose_grnn <- function(x,n,mode,strategy,cc=clean_cache) {
 h <- 93-n; ends <- n-h-c(12,6,0)
 inner <- lapply(ends,function(end)apply(grid(x,end,h,mode,strategy,cc),2,mape,truth=as.numeric(x)[end+seq_len(h)]))
 loss <- rowMeans(do.call(cbind,inner));stopifnot(all(is.finite(loss)))
 pred <- grid(x,n,h,mode,strategy,cc)
 outer <- apply(pred,2,mape,truth=as.numeric(x)[n+seq_len(h)])
 list(pred=pred,selected=which.min(loss),oracle=which.min(outer),ends=ends,loss=loss)
}
errors <- list();predictions <- list();selections <- list();control <- list()
add <- function(n,mode,strategy,selection,pred) {
 h <- 93-n;truth <- as.numeric(raw)[n+seq_len(h)]
 errors[[length(errors)+1L]] <<- data.frame(origin=n,preprocessing=mode,strategy=strategy,selection=selection,t(metrics(pred,truth)))
 predictions[[length(predictions)+1L]] <<- data.frame(origin=n,target=n+seq_len(h),preprocessing=mode,strategy=strategy,selection=selection,truth=truth,prediction=pred)
}
for(n in 75:91) {
 h <- 93-n
 a <- grid(raw,n,h,'G','MIMO');chosen <- which.min(apply(a,2,mape,truth=as.numeric(raw)[n+seq_len(h)]))
 b <- grid(raw,n,h,'G','recursive')
 original_pred <- as.numeric(original$GRNN_forecasts[[n-74]])
 stopifnot(max(abs(original_pred-b[,chosen]))<1e-10)
 add(n,'G','mixed','test-MIMO',original_pred)
 control[[length(control)+1L]] <- data.frame(origin=n,selected_MIMO_lag=chosen,recursive_oracle_lag=which.min(apply(b,2,mape,truth=as.numeric(raw)[n+seq_len(h)])),point_delta=max(abs(original_pred-b[,chosen])))
 for(mode in c('G','P'))for(strategy in c('MIMO','recursive')) {
  fit <- choose_grnn(raw,n,mode,strategy)
  add(n,mode,strategy,'O',fit$pred[,fit$oracle]);add(n,mode,strategy,'V',fit$pred[,fit$selected])
  selections[[length(selections)+1L]] <- data.frame(origin=n,preprocessing=mode,strategy=strategy,inner_ends=paste(fit$ends,collapse=';'),selected=fit$selected,oracle=fit$oracle)
 }
 checkpoint <- list(errors=do.call(rbind,errors),predictions=do.call(rbind,predictions),selections=do.call(rbind,selections),control=do.call(rbind,control))
 saveRDS(checkpoint,file.path(root,'results/grnn-checkpoint.rds'))
 saveRDS(as.list(grnn_cache),file.path(root,'results/grnn-cache.rds'))
 cat('GRNN completed origin',n,'of 91\n');flush.console()
}
out <- checkpoint
out$summary <- do.call(rbind,lapply(split(out$errors,interaction(out$errors[,c('preprocessing','strategy','selection')],drop=TRUE)),function(z)data.frame(z[1,c('preprocessing','strategy','selection')],median_MAPE=median(z$MAPE),mean_MAPE=mean(z$MAPE),median_MAE=median(z$MAE),median_RMSE=median(z$RMSE))))
rownames(out$summary)<-NULL
invariance <- list()
for(n in c(75,83,91))for(strategy in c('MIMO','recursive')) {
 a <- choose_grnn(raw,n,'P',strategy)
 altered <- raw;altered[(n+1):93]<-raw[(n+1):93]*20+500
 saved_cache <- grnn_cache;grnn_cache <- new.env(parent=emptyenv())
 b <- choose_grnn(altered,n,'P',strategy,new.env(parent=emptyenv()))
 grnn_cache <- saved_cache
 delta <- max(abs(a$pred-b$pred));same <- identical(a$selected,b$selected)
 stopifnot(delta==0,same)
 invariance[[length(invariance)+1L]] <- data.frame(origin=n,strategy=strategy,max_forecast_delta=delta,same_selection=same)
}
for(name in c('errors','predictions','selections','summary','control'))write_json(out[[name]],file.path(root,paste0('results/grnn-',name,'.json')),pretty=TRUE,digits=12)
write_json(do.call(rbind,invariance),file.path(root,'results/grnn-invariance.json'),pretty=TRUE,digits=12)
write_json(list(source_version='0.1.0',source_loaded=TRUE,original_version_unrecorded=TRUE,original_median_MAPE=median(original$GRNN_errors$MAPE),max_historical_rounded_metric_difference=historical_delta,original_search_strategy='MIMO',original_refit_strategy='recursive',warnings=unique(warnings),fallbacks=fallbacks),file.path(root,'results/grnn-verification.json'),pretty=TRUE,digits=12,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(root,'results/grnn-session-info.txt'))
print(out$summary)
