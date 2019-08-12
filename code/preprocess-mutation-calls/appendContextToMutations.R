minMuts <- 10

source("code/utils/funcs.R")
cargs <- commandArgs(trail=TRUE)
if(length(cargs)!=3)
    stop("Usage: R --args <inFile> <fastaFile> <outFile>")
inFile <- cargs[1]
if(!file.exists(inFile))
    stop("Error - infile ",inFile," not found")
faFile <- cargs[2]
outFile <- cargs[3]

ignoreMissing <- FALSE

message("Reading ",inFile)
d <- readRDS(inFile)
d <- filterMutTable(d)

message("Reading DNA context from ",faFile)
fa <- read.table(faFile,as.is=TRUE)
fa <- t(matrix(fa[,1],nrow=2))

out <- appendMutationContext(dat=d,fasta=fa,ignoreMissing=ignoreMissing)
message("Writing to ",outFile)
saveRDS(out,outFile)

