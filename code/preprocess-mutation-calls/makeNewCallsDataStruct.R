addAlleleFrequencies <- function(dat) {
    if(!("altIsMinor"%in%colnames(x))) {
        message("Adding allele frequencies, altIsMinor flag")
        o1 <- data.frame(iraf=dat$iref/(dat$iref + dat$ialt),iaaf=dat$ialt/(dat$iref + dat$ialt),fraf=dat$fref/(dat$fref + dat$falt),faaf=dat$falt/(dat$fref + dat$falt))
        dat <- cbind(dat,o1)
        dat <- cbind(dat,dat$fraf > dat$faaf)
        colnames(dat)[ncol(dat)] <- "altIsMinor"
    }
    dat
}

maskSite <- function(i,mut,cnv,ips) {
    cni <- cnv[cnv$line==ips[i],]
    ret <- FALSE
    if(nrow(cni)>0) {
        if(any(mut[i,"chr"]==cni$chr & mut[i,"pos"]>=cni$start & mut[i,"pos"]<=cni$end))
            ret <- TRUE
    }
    ret
}

markSitesInCnvs <- function(dat) {
    message("Reading CNV location file ",cnvFile)
    ips <- gsub("HPSI.*-","",dat$ips)
    cn <- read.table(cnvFile,sep="\t",he=T,strings=F)
    dat$inCnv <- sapply(1:nrow(dat),maskSite,mut=dat,cnv=cn,ips=ips)
    dat
}

doCnvOverlap <- function(line,dat,cn,bedtoolsPath="/software/team170/bedtools2/bin/bedtools") {
    message("Checking CNV overlap for ",line)
    tempdir(check=TRUE)
    bedFile1 <- tempfile()
    bedFile2 <- tempfile()
    bedFile3 <- tempfile()
    out <- dat[dat$ips==line,c("chr","pos","pos")]
    write.table(out,file=bedFile1,row=F,col=F,qu=F,sep="\t")
    write.table(cn[which(cn$line==convertToFriendly(line)),c("chr","start","end")],file=bedFile2,sep="\t",row=F,col=F,qu=F)
    cmd <- paste("sort -k 1,1n -k 2,2n",bedFile2,"| ",bedtoolsPath," merge -i stdin >",bedFile3)
    ## message(cmd)
    system(cmd)
    cmd <- paste(bedtoolsPath,"intersect -a ",bedFile1,"-b",bedFile3)
    ## message(cmd)
    r <- system(cmd,intern=TRUE)
    ret <- ""
    if(length(r)>0) {
        m <- t(matrix(unlist(strsplit(r,"\t")),nrow=3))
        ret <- which(dat$pos%in%m[,3] & dat$ips==line)
    }
    system(paste("rm ",bedFile1,bedFile2,bedFile3))
    ret
}

markSitesInCnvs2 <- function(dat,bedtoolsPath="/software/team170/bedtools2/bin/bedtools") {
    message("Reading CNV location file ",cnvFile)
    cn <- read.table(cnvFile,sep="\t",he=T,strings=F)
    cn <- cn[cn$chr!="X" & cn$chr!="Y",]
    lines <- unique(dat$ips)
    lines <- lines[convertToFriendly(lines)%in%cn$line]
    rownumbers <- sapply(lines,doCnvOverlap,dat=dat,cn=cn,bedtoolsPath=bedtoolsPath)
    rownumbers <- na.omit(as.numeric(unique(unlist(rownumbers))))
    dat$inCnv <- FALSE
    dat[unique(unlist(rownumbers)),"inCnv"] <- TRUE
    dat
}

addExacUk1kgAf <- function(dat) {
    dat$exFreq <- as.numeric(gsub("^\\.$",0,gsub("EXAC=","",dat$exac)))
    dat$ukFreq <- as.numeric(gsub("^\\.$",0,gsub("1KG=","",dat$uk1kg)))
    dat
}

cleanOrphanDinucs <- function(dat) {
    dat <- dat[order(dat$ips,dat$chr,dat$pos),]
    dinucs <- dat[dat$dinuc,]
    findOrphan <- function(i) { ret <- NA; if(dinucs[i-1,"pos"]+1!=dinucs[i,"pos"] & dinucs[i+1,"pos"]-1!=dinucs[i,"pos"]) { ret <- i }; ret }
    idx <- na.omit(sapply(2:(nrow(dinucs)-1),findOrphan))
    message("Found ",length(idx)," orphan dinucs")
    if(dinucs[1,"pos"]+1!=dinucs[2,"pos"])
        idx <- c(1,idx)
    if(dinucs[nrow(dinucs),"pos"]-1!=dinucs[(nrow(dinucs)-1),"pos"])
        idx <- c(idx,nrow(dinucs))
    dinucs <- dinucs[-idx,]
    checkDinucOrdering(dinucs)
    dat <- dat[!dat$dinuc,]
    dat <- rbind(dinucs,dat)
    dat <- dat[order(dat$ips,dat$chr,dat$pos),]
    dat
}

source("code/utils/funcs.R")
source("code/preprocess-mutation-calls/funcs.R")
cnvFile <- "Data/cnv/atab_20161018.tsv"
if(!file.exists(cnvFile))
    stop("CNV location file ",cnvFile," not found")

cargs <- commandArgs(trail=TRUE)
if(length(cargs)<2)
    stop("Usage: R --args <inFile> <outFile>")
inFile <- cargs[1]
outFile <- cargs[2]
if(!file.exists(inFile))
    stop("Error - infile ",inFile," not found")
unfilteredOutFile <- gsub(".rd[as]$",".unfiltered.rda",outFile)

## Start processing
popAfThresh <- 0.001
afThresh <- 0.7
x <- readRawMutFileV4(inFile)
message("Removing ",sum(x$chr>22)," sex chromosome sites")
x <- x[x$chr<23,]
message("Removing ",sum(grepl(",",x$alt))," sites with multiple ALT alleles")
x <- x[grep(",",x$alt,invert=TRUE),]
x <- addAlleleFrequencies(dat=x)
x$ndon <- markSegSites2(dat=x,f=inFile)
x <- markSitesInCnvs2(dat=x,bedtoolsPath="/software/team170/bedtools2/bin/bedtools")
# x <- addExacUk1kgAf(dat=x)
message("Removing ",sum(x$exFreq > popAfThresh | x$ukFreq > popAfThresh)," sites with AF > ",popAfThresh," in ExAC or UK1KG")
x <- x[x$exFreq < popAfThresh & x$ukFreq < popAfThresh,]
message("Removing ",sum(x$ndon!=1)," sites found in > 1 donor")
x <- x[x$ndon==1,]
message("Removing ",sum(x$inCnv)," sites overlapping CNVs")
x <- x[!x$inCnv,]
message("Writing allele frequency unfiltered file to ",unfilteredOutFile)
saveRDS(x,file=unfilteredOutFile)
message("Removing ",sum(x$faaf>afThresh | x$iaaf>afThresh)," sites with AF > ",afThresh)
x <- x[x$faaf < afThresh & x$iaaf < afThresh,]
x$dinuc <- x$dinuc>1 ## change to logical
if(any(x[seq(1,nrow(x)-1,2),"pos"]+1==x[seq(2,nrow(x),2),"pos"])) {
    message("Adjacent mutations detected - assuming dinucleotides not collapsed in this file")
    message("Removing orphan dinucs")
    x <- cleanOrphanDinucs(dat=x)
}
message("Writing to ",outFile)
saveRDS(x,file=outFile)
