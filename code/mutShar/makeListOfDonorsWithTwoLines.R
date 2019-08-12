cargs <- commandArgs(trail=TRUE)
stopifnot(length(cargs)==2)
inFile <- cargs[1]
outFile <- cargs[2]
message("Reading ",inFile)
l <- readRDS(inFile)
message("Writing to ",outFile)
write.table(unique(l$m2l$donor),outFile,col=F,row=F,qu=F)
