makeSharedLineDataStruct <- function(dat) {
    m <- t(matrix(unlist(strsplit(names(table(dat$ips)),"[-_]")),nrow=3))
    colnames(m) <- c("batch","donor","line")
    m <- as.data.frame(m,stringsAsFactors=FALSE)
    d <- names(table(m$donor)[table(m$donor)==2])
    m <- m[m$donor%in%d,]
    ## Remove lines with double subscripts
    d <- m[grep("[0-9][0-9]",m$line),"donor"]
    m2l <- m[!(m$donor%in%d),]
    m2l <- cbind(m2l,paste(m2l[,1],"-",m2l[,2],"_",m2l[,3],sep=""),stringsAsFactors=FALSE)
    m2l <- cbind(m2l,paste(m2l[,2],"_",m2l[,3],sep=""),stringsAsFactors=FALSE)
    colnames(m2l) <- c("batch","donor","line","name","friendly")
    m2l <- cbind(m2l,sapply(1:nrow(m2l),function(i) { sum(dat$ips==m2l[i,"name"]) }),stringsAsFactors=FALSE)
    colnames(m2l) <- c("batch","donor","line","name","friendly","count")
    list(m2l=m2l)
}

cargs <- commandArgs(trail=TRUE)
inFile <- cargs[1]
outFile <- cargs[2]

message("Reading ",inFile)
x <- readRDS(inFile)
out <- makeSharedLineDataStruct(dat=x)
message("Writing to ",outFile)
saveRDS(out,file=outFile)
