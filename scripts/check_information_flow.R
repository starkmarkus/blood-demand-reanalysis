root <- normalizePath(commandArgs(TRUE)[1],winslash='/')
exprs <- parse(file.path(root,'scripts/run.R'))
for(ex in exprs) {
 if(is.call(ex) && identical(ex[[1]],as.name('<-')) && identical(ex[[2]],as.name('pi_checks'))) break
 eval(ex,envir=.GlobalEnv)
}
# Refit the prefix rather than inferring preprocessing from stored forecast errors.
clean <- prepare(raw,77,new.env(parent=emptyenv()))
preds <- fromJSON(file.path(root,'results/factorial-predictions.json'))
baseline <- list()
for(method in c('Mean12','Naive','SeasonalNaive')) {
 z <- preds[preds$origin==77 & preds$method==method,]
 p <- switch(method,Mean12=rep(mean(tail(clean,12)),16),Naive=rep(tail(clean,1),16),SeasonalNaive=rep(tail(clean,12),length.out=16))
 delta <- max(abs(z$prediction-as.numeric(p)))
 stopifnot(nrow(z)==16L,delta<1e-9)
 baseline[[method]] <- list(max_forecast_difference=delta)
}
source(file.path(root,'inputs/original/Lag-normalise-invert-scaling-script.R'))
# Mirror only the original LSTM's deterministic data preparation and inversion inputs.
# No network is fitted and no numerical LSTM performance claim is made.
lstm_inputs <- function(series,origin) {
 horizon <- length(series)-origin
 diffed <- diff(series,differences=1)
 supervised <- lag_transform(diffed,1)
 N <- length(series)-1;n <- N-horizon
 train <- supervised[1:n,];test <- supervised[(n+1):N,]
 # Inspect availability before scaling; no claim of reproducing the historical
 # normalization/Keras runtime is made by this deterministic check.
 list(train=as.numeric(as.matrix(train)),x=as.numeric(test[,1]),offset=as.numeric(series[n+seq_len(horizon)]))
}
n <- 75
a <- lstm_inputs(global,n)
altered <- global;altered[(n+1):93] <- global[(n+1):93]*20+500
b <- lstm_inputs(altered,n)
stopifnot(identical(a$train,b$train),a$x[1]==b$x[1],a$offset[1]==b$offset[1],any(a$x[-1]!=b$x[-1]),any(a$offset[-1]!=b$offset[-1]))
map <- data.frame(step=1:18,target_index=76:93,latest_observed_input_index=75:92,available_at_fixed_origin=75:92<=75,x_original=a$x,x_perturbed=b$x,offset_original=a$offset,offset_perturbed=b$offset)
write_json(list(prefix77=list(raw_last=as.numeric(raw[77]),prepared_last=as.numeric(clean[77]),checks=baseline),lstm=list(scope='Unscaled lagged inputs and inversion offsets only; no normalization or LSTM runtime reproduction',origin=75,unchanged_training=TRUE,changed_input_steps=which(a$x!=b$x),changed_inversion_steps=which(a$offset!=b$offset),availability=map)),file.path(root,'results/information-flow-checks.json'),pretty=TRUE,digits=12,auto_unbox=TRUE)
cat('Prefix-77 baselines and LSTM information-availability checks passed.\n')
