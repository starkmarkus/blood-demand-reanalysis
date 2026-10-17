root <- normalizePath(commandArgs(TRUE)[1],winslash='/')

library(jsonlite)
d <- fromJSON(file.path(root,'results/factorial-summary.json'))
out <- file.path(root,'figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
spec <- data.frame(prep=c('G','P',rep('P',12)),
 selection=c('O','O',rep('V',8),rep('none',4)),
 rule=c('matched','matched','matched','6','12','18',rep('matched',4),rep('none',4)),
 seed=c(rep(1,6),2:5,rep(1,4)),
 label=c('Published setup (G-O)','Prefix / test selection (P-O)',
 'Historical: matched H, seed 1','Historical: H = 6, seed 1',
 'Historical: H = 12, seed 1','Historical: H = 18, seed 1',
 'Historical: matched H, seed 2','Historical: matched H, seed 3',
 'Historical: matched H, seed 4','Historical: matched H, seed 5',
 'ARIMA (prefix)','12-month mean (prefix)','Naive (prefix)','Seasonal naive (prefix)'))
values <- list()
for (method in c('KNN','NNAR')) {
 vals <- rep(NA_real_,14)
 for (i in 1:14) {
   if(method=='KNN' && i %in% 7:10) next
   m <- if(i>10) c('ARIMA','Mean12','Naive','SeasonalNaive')[i-10] else method
   r <- d[d$method==m & d$preprocessing==spec$prep[i] & d$selection==spec$selection[i] & d$rule==spec$rule[i] & d$seed==spec$seed[i],]
   stopifnot(nrow(r)==1L)
   vals[i] <- r$median_MAPE
 }
 values[[method]] <- vals
}
draw <- function() {
 layout(matrix(c(1,2,3),1,3),widths=c(3.2,2.5,2.5))
 par(family='sans',cex=1.0,mgp=c(2.1,.5,0),tcl=-.25)
 par(mar=c(4.4,.1,2.4,.1))
 plot(NA,xlim=c(0,1),ylim=c(.5,15.2),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
 text(.98,14:1,labels=spec$label,adj=1,cex=.93)
 for (j in 1:2) {
   method <- c('KNN','NNAR')[j]
   par(mar=c(4.4,.5,2.4,.6))
   plot(NA,xlim=c(10,30.5),ylim=c(.5,15.2),xlab='Median MAPE (%)',ylab='',
        axes=FALSE,main=method,xaxs='i',yaxs='i')
   y <- 14:1
   abline(v=seq(10,30,5),col='#e0e0e0',lwd=.7)
   abline(h=c(12.5,4.5),col='#999999',lty=3)
   axis(1,at=seq(10,25,5))
   vals <- values[[method]]
   col <- c(rep('#a65a24',2),rep('#20577b',8),rep('#333333',4))
   points(vals,y,pch=c(rep(5,2),rep(16,8),rep(15,4)),col=col,cex=1.0)
   text(30.5,y,labels=ifelse(is.na(vals),'--',sprintf('%.2f',vals)),adj=1,cex=.85,col=col)
   text(30.5,14.9,'Value',adj=1,cex=.85)
 }
}
pdf(file.path(out,'Figure_1.pdf'),width=8.2,height=6.2,useDingbats=FALSE)
draw();dev.off()
png(file.path(out,'Figure_1.png'),width=2460,height=1860,res=300)
draw();dev.off()
write_json(list(rows=spec,values=values),file.path(out,'figure-data.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=15)
