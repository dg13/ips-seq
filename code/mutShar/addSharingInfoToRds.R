checkOverlapWithShared <- function(l,sharingDat) {
    p1 <- paste(l$chr,l$pos,sep=":")
    p2 <- sharingDat[sharingDat$donor==convertToFriendly(l[1,"fib"]),]
    p2 <- paste(p2[,1],p2[,2],sep=":")
    p1%in%p2
}

cargs <- commandArgs(trail=TRUE)
stopifnot(length(cargs)==3)

inFile <- cargs[1]
pairingFile <- cargs[2]
outFile <- cargs[3]

source("code/utils/funcs.R")

## inFile <- "Data/newCalls2/rdas/mut-wgs-fdr5.collapseDinucs.rda"
## pairingFile <- "Data/mutShar/paired.all.txt.gz"
message("Reading pairing info from ",pairingFile)
p <- read.table(pairingFile,col.names=c("chr","pos","donor"))
message("Reading ",inFile)
x <- readRDS(inFile)
message("Splitting into single lines and double lines per donor")
x1lines <- x[!convertToFriendly(x$fib)%in%p$donor,]
x2lines <- x[convertToFriendly(x$fib)%in%p$donor,]
x2lines <- x2lines[order(x2lines$fib),]
d <- split(x2lines,x2lines$fib)
d <- d[order(names(d))]
message("Searching for overlap with sites found in ",pairingFile)
sh <- lapply(d,checkOverlapWithShared,sharingDat=p)
x1lines$shar <- NA
x2lines$shar <- unlist(sh)
pa <- x2lines[x2lines$shar,"pos"]
stopifnot(all(pa%in%p$pos)) ## Sanity check 
message("Concatenating all data")
out <- rbind(x1lines,x2lines)
saveRDS(out,file=outFile)
