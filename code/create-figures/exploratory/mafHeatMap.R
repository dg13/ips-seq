source("code/utils/col-rb.R")

cargs <- commandArgs(trail=TRUE)
stopifnot(length(cargs)==2)
inFile <- cargs[1]
outFile <- cargs[2]

lwd <- 1.7
afThresh <- 0.7

Lab.palette <- colorRampPalette(col.rb, space = "Lab")
br <- c(0,seq(0,1,0.02))
cex <- 1.7
message("Reading ",inFile)
m <- readRDS(inFile)
message("Plotting to ",outFile)
png(file=outFile,res=250,units="in",hei=5,wid=5)
layout(matrix(1:4,nrow=2),wid=c(1,3),hei=c(3,1))
par(mar=c(2,3.8,1,0),mgp=c(2.4,0.5,0),tcl=-0.25,cex.axis=cex,cex.lab=cex)
r <- hist(m$iaaf,br=br,plot=F)
cols <- rep(grey(0.6),length(r$breaks))
cols[1] <- "red"
out <- barplot(r$counts/sum(r$counts),horiz=TRUE,yaxs="i",tcl=0,ylab="Mutation AF: IPSCs",axes=F,col=cols)
axis(1,at=c(0,0.1,0.2),tcl=-0.25)
axis(2,at=quantile(out,seq(0,1,0.2)),labels=seq(0,1,0.2),las=2)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",axes=F)
par(mar=c(2,0.5,1,1),mgp=c(1.7,0.5,0))
smoothScatter(m$faaf,m$iaaf, colramp = Lab.palette,yaxs="i",xaxs="i",xlab="",ylab="",axes=F,xlim=c(0,1),ylim=c(0,1))
lines(c(0,afThresh),c(afThresh,afThresh),lty=2,col="red")
lines(c(afThresh,afThresh),c(0,afThresh),lty=2,col="red")
r <- hist(m$faaf,br=br,plot=F)
par(mar=c(3,0.5,0,1),mgp=c(1.7,0.5,0))
cols <- rep(grey(0.6),length(r$breaks))
cols[1] <- "red"
out <- barplot(r$counts/sum(r$counts),xaxs="i",tcl=0,xlab="Mutation AF: Fibroblasts",axes=F,col=cols)
axis(1,at=quantile(out,seq(0,1,0.2)),labels=seq(0,1,0.2))
axis(2,at=c(0,0.2,0.35),las=2)
legend("topright",c("AF=0"),fill="red",bty="n",cex=1.8)
dev.off()
