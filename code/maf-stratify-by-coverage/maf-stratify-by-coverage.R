make4PanelPlot <- function(dat,main="") {
    layout(matrix(1:4,nrow=2),wid=c(1,3),hei=c(3,1))
    par(mar=c(2,3.8,1,0),mgp=c(2.4,0.5,0),tcl=-0.25,cex.axis=cex,cex.lab=cex)
    r <- hist(dat$iaaf,br=br,plot=F)
    cols <- rep(grey(0.6),length(r$breaks))
    cols[1] <- "red"
    out <- barplot(r$counts/sum(r$counts),horiz=TRUE,yaxs="i",tcl=0,ylab="Mutation AF: IPSCs",axes=F,col=cols)
    axis(1,at=c(0,0.1,0.2),tcl=-0.25)
    axis(2,at=quantile(out,seq(0,1,0.2)),labels=seq(0,1,0.2),las=2)
    plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F)
    text(-0.4,0.5,main,cex=2.1,xpd=NA)
    par(mar=c(2,0.5,1,1),mgp=c(1.7,0.5,0))
    smoothScatter(dat$faaf,dat$iaaf, colramp = Lab.palette,yaxs="i",xaxs="i",xlab="",ylab="",axes=F,xlim=c(0,1),ylim=c(0,1))
    lines(c(0,afThresh),c(afThresh,afThresh),lty=2,col="red")
    lines(c(afThresh,afThresh),c(0,afThresh),lty=2,col="red")
    r <- hist(dat$faaf,br=br,plot=F)
    par(mar=c(3,0.5,0,1),mgp=c(1.7,0.5,0))
    cols <- rep(grey(0.6),length(r$breaks))
    cols[1] <- "red"
    out <- barplot(r$counts/sum(r$counts),xaxs="i",tcl=0,xlab="Mutation AF: Fibroblasts",axes=F,col=cols)
    axis(1,at=quantile(out,seq(0,1,0.2)),labels=seq(0,1,0.2))
    axis(2,at=c(0,0.2,0.35),las=2)
    legend("topright",c("AF=0"),fill="red",bty="n",cex=1.8)
}

wes <- readRDS("Data/mut-files/rdas/mut-wes-fdr5.collapseDinucs.unfiltered.rda")
wgs <- readRDS("Data/mut-files/rdas/mut-wgs-fdr5.collapseDinucs.unfiltered.rda")
outDir <- "plots/maf-heatmaps/stratified-by-coverage/"

p1 <- paste(wes$chr,wes$pos,sep=":")
p2 <- paste(wgs$chr,wgs$pos,sep=":")
wgs <- wgs[!p2%in%p1,]
m <- rbind(wes,wgs)
p <- paste(m$chr,m$pos,m$ips,sep=":")
stopifnot(any(duplicated(p)))

fcov <- m$fref + m$falt
icov <- m$iref + m$ialt

null <- sapply(seq(30,60,5),function(i) {
    dat <- m[fcov > i & icov > i,]
    outFile <- paste0(outDir,"/maf-heatmap-",i,"X.png")
    message("Plotting ",outFile)
    png(file=outFile,res=250,units="in",hei=5,wid=5)
    make4PanelPlot(dat=dat,main=paste0(i,"X"))
    dev.off()
})

p1 <- paste(wes$chr,wes$pos,sep=":")
p2 <- paste(wgs$chr,wgs$pos,sep=":")
wgs <- wgs[!p2%in%p1,]

fcov <- wes$fref + wes$falt
icov <- wes$iref + wes$ialt
wesMax <- mean(c(fcov,icov))*3
sum(fcov > wesMax | icov > wesMax)

fcov <- wgs$fref + wgs$falt
icov <- wgs$iref + wgs$ialt
wgsMax <- mean(c(fcov,icov))*3
sum(fcov > wgsMax | icov > wgsMax)

## Check if there are sites present in multiple donors
p1 <- paste(wes$chr,wes$pos,wes$fib,sep=":")
any(duplicated(p1))
p1 <- unique(paste(wes$chr,wes$pos,wes$fib,sep=":")) ## This collapses mutations from the same donor (which are legit) but keeps any that are replicated across donors (because these will have a different donor ID)
m1 <- t(matrix(unlist(strsplit(p1,":")),nrow=3))
any(duplicated(p1))

p1 <- paste(wgs$chr,wgs$pos,wgs$fib,sep=":")
any(duplicated(p1))
p1 <- unique(paste(wgs$chr,wgs$pos,wgs$fib,sep=":")) ## This collapses mutations from the same donor (which are legit) but keeps any that are replicated across donors (because these will have a different donor ID)
m1 <- t(matrix(unlist(strsplit(p1,":")),nrow=3))
any(duplicated(p1))

