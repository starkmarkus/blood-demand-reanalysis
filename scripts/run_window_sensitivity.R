args <- commandArgs(TRUE)
root <- normalizePath(args[1],winslash='/')
# Load the unchanged computation definitions, stopping before the original experiment.
exprs <- parse(file.path(root,'scripts/run.R'))
for(ex in exprs) {
 if(is.call(ex) && identical(ex[[1]],as.name('<-')) && identical(ex[[2]],as.name('pi_checks'))) break
 eval(ex,envir=.GlobalEnv)
}
if(length(args)>=2L && file.exists(args[2])) {
 old <- readRDS(args[2])
 list2env(old$clean,envir=clean_cache);list2env(old$nn,envir=nn_cache)
}
designs <- expand.grid(windows=2:4,spacing=c(3,6))
ends_for <- function(n,h,j,s) n-h-rev(0:(j-1))*s
choose_windows <- function(x,n,method,j,s,ids,cc=clean_cache,gc=grid_cache,nc=nn_cache) {
 h <- 93-n;ends <- ends_for(n,h,j,s)
 stopifnot(min(ends)>=24,max(ends+h)<=n)
 folds <- lapply(ends,function(end) {
  z <- get_grid(x,end,h,'P',method,1,cc,gc,nc)
  stopifnot(all(ids %in% z$ids))
  apply(z$pred[,match(ids,z$ids),drop=FALSE],2,mape,truth=as.numeric(x)[end+seq_len(h)])
 })
 losses <- do.call(cbind,folds); avg <- rowMeans(losses)
 stopifnot(all(is.finite(avg)))
 z <- get_grid(x,n,h,'P',method,1,cc,gc,nc)
 pred <- z$pred[,match(ids,z$ids),drop=FALSE]
 selected <- which.min(avg)
 oracle <- which.min(apply(pred,2,mape,truth=as.numeric(x)[n+seq_len(h)]))
 list(pred=pred,selected=selected,oracle=oracle,ends=ends,losses=avg)
}
errors <- list();predictions <- list();selections <- list();checks <- list()
reference <- fromJSON(file.path(root,'results/factorial-errors.json'))
for(n in 75:91) {
 h <- 93-n
 for(method in c('KNN','NNAR')) {
  allends <- unique(unlist(lapply(seq_len(nrow(designs)),function(i)ends_for(n,h,designs$windows[i],designs$spacing[i]))))
  grids <- lapply(c(allends,n),function(end)get_grid(raw,end,h,'P',method))
  ids <- Reduce(intersect,lapply(grids,function(z)z$ids))
  native <- choose(raw,n,h,h,'P',method,1)
  # Explicit column matching avoids data-mask ambiguity.
  prior <- reference[reference$origin==n & reference$method==method & reference$preprocessing=='P' & reference$selection=='V' & reference$rule=='matched' & reference$seed==1,]
  delta <- abs(mape(native$pred[,native$selected],as.numeric(raw)[n+seq_len(h)])-prior$MAPE)
  stopifnot(nrow(prior)==1L,delta<1e-9)
  checks[[length(checks)+1L]] <- data.frame(origin=n,method=method,native_primary_delta=delta)
  for(i in seq_len(nrow(designs))) {
   j <- designs$windows[i];s <- designs$spacing[i]
   fit <- choose_windows(raw,n,method,j,s,ids)
   for(selection in c('O','V')) {
    id <- if(selection=='O')fit$oracle else fit$selected
    pred <- fit$pred[,id];truth <- as.numeric(raw)[n+seq_len(h)]
    errors[[length(errors)+1L]] <- data.frame(origin=n,method=method,windows=j,spacing=s,selection=selection,t(metrics(pred,truth)))
    predictions[[length(predictions)+1L]] <- data.frame(origin=n,target=n+seq_len(h),method=method,windows=j,spacing=s,selection=selection,truth=truth,prediction=pred)
   }
   selections[[length(selections)+1L]] <- data.frame(origin=n,method=method,windows=j,spacing=s,inner_ends=paste(fit$ends,collapse=';'),candidates=length(ids),candidate_ids=paste(ids,collapse=';'),selected=ids[fit$selected],oracle=ids[fit$oracle])
  }
 }
 checkpoint <- list(errors=do.call(rbind,errors),predictions=do.call(rbind,predictions),selections=do.call(rbind,selections),checks=do.call(rbind,checks))
 saveRDS(checkpoint,file.path(root,'results/window-checkpoint.rds'))
 cat('Window sensitivity completed origin',n,'of 91\n');flush.console()
}
out <- checkpoint
out$summary <- do.call(rbind,lapply(split(out$errors,interaction(out$errors[,c('method','windows','spacing','selection')],drop=TRUE)),function(z)data.frame(z[1,c('method','windows','spacing','selection')],median_MAPE=median(z$MAPE),mean_MAPE=mean(z$MAPE),median_MAE=median(z$MAE),median_RMSE=median(z$RMSE))))
rownames(out$summary)<-NULL
invariance <- list()
for(n in c(75,83,91)) for(method in c('KNN','NNAR')) {
 h <- 93-n; allends <- unique(unlist(lapply(seq_len(nrow(designs)),function(i)ends_for(n,h,designs$windows[i],designs$spacing[i]))))
 ids <- Reduce(intersect,lapply(c(allends,n),function(end)get_grid(raw,end,h,'P',method)$ids))
 altered <- raw;altered[(n+1):93] <- raw[(n+1):93]*20+500
 # The largest window design exercises all four historical scoring windows.
 a <- choose_windows(raw,n,method,4,6,ids)
 b <- choose_windows(altered,n,method,4,6,ids,new.env(parent=emptyenv()),new.env(parent=emptyenv()),new.env(parent=emptyenv()))
 delta <- max(abs(a$pred-b$pred)); same <- identical(a$selected,b$selected)
 stopifnot(delta==0,same)
 invariance[[length(invariance)+1L]] <- data.frame(origin=n,method=method,max_forecast_delta=delta,same_selection=same)
}
for(name in c('errors','predictions','selections','summary','checks'))write_json(out[[name]],file.path(root,paste0('results/window-',name,'.json')),pretty=TRUE,digits=12)
write_json(do.call(rbind,invariance),file.path(root,'results/window-invariance.json'),pretty=TRUE,digits=12)
write_json(list(warnings=unique(warnings),fallbacks=fallbacks),file.path(root,'results/window-diagnostics.json'),pretty=TRUE)
saveRDS(list(clean=as.list(clean_cache),nn=as.list(nn_cache)),file.path(root,'results/window-cache.rds'))
capture.output(sessionInfo(),file=file.path(root,'results/window-session-info.txt'))
print(out$summary)
