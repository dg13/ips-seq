readRawMutFileV3 <- function(f) {
    colNames <- c("chr","pos","fib","ips","pval","ref","alt","fref","falt","iref","ialt","ukFreq","exFreq","clinsig","sift","pphen","hgnc","ensGene","csq","dinuc")
    message("Reading raw mutation calls from ",f)
    ret <- read.table(f,sep="\t",stringsAsFactors=FALSE,col.names=colNames)
    ## ret$dinuc <- as.logical(ret$dinuc-1)
    if(sum(ret$chr=="X" | ret$chr=="Y")>0) {
        message("Setting X and Y to 23, 24")
        ret[ret$chr=="X","chr"] <- 23
        ret[ret$chr=="Y","chr"] <- 24
    }
    ret$chr <- as.numeric(as.character(ret$chr))
    ret <- ret[order(ret$chr),]
    ret
}

readRawMutFileV4 <- function(f) {
    colNames <- c("chr","pos","fib","ips","pval","ref","alt","fref","falt","iref","ialt","ukFreq","exFreq","clinsig","sift","pphen","hgnc","ensGene","csq","dinuc")
    message("Reading raw mutation calls from ",f)
    ret <- read.table(f,sep="\t",stringsAsFactors=FALSE,col.names=colNames)
    ## ret$dinuc <- as.logical(ret$dinuc-1)
    if(sum(ret$chr=="X" | ret$chr=="Y")>0) {
        message("Setting X and Y to 23, 24")
        ret[ret$chr=="X","chr"] <- 23
        ret[ret$chr=="Y","chr"] <- 24
    }
    ret$chr <- as.numeric(as.character(ret$chr))
    ret[grep(",",ret$ukFreq),"ukFreq"] <- as.character(unlist(lapply(d <- strsplit(ret[grep(",",ret$ukFreq),"ukFreq"],","),function(l) { sort(as.numeric(l[l!="."]))[1] })))
    ret[grep(",",ret$exFreq),"exFreq"] <- as.character(unlist(lapply(d <- strsplit(ret[grep(",",ret$exFreq),"exFreq"],","),function(l) { sort(as.numeric(l[l!="."]))[1] })))
    ret[ret$exFreq==".","exFreq"] <- 0
    ret[ret$ukFreq==".","ukFreq"] <- 0
    ret$exFreq <- as.numeric(ret$exFreq)
    ret$ukFreq <- as.numeric(ret$ukFreq)    
    ret <- ret[order(ret$chr),]
    ret
}

readRawMutFileV2 <- function(f,fill=FALSE) {
    rawColNames1 <- c("ft","ips","fib","chr","pos","iref","ialt","fref","falt","pval","ref","alt","csq","info","exac","uk1kg","dinuc")
    rawColNames2 <- c("ft","ips","fib","chr","pos","iref","ialt","fref","falt","pval","ref","alt","csq","exac","uk1kg","dinuc")
    message("Reading raw mutation calls from ",f)
    ret <- read.table(f,sep="\t",as.is=TRUE,fill=fill)
    if(fill)
        ret <- ret[-(nrow(ret)),]
    if(ncol(ret)==17) {
        colnames(ret) <- rawColNames1
    } else if (ncol(ret)==16) {
        colnames(ret) <- rawColNames2
    } else {
        stop("Error - unexpected column number ",ncol(ret))
    }
    ## ret$dinuc <- as.logical(ret$dinuc-1)
    if(sum(ret$chr=="X" | ret$chr=="Y")>0) {
        message("Setting X and Y to 23, 24")
        ret[ret$chr=="X","chr"] <- 23
        ret[ret$chr=="Y","chr"] <- 24
    }
    ret$chr <- as.numeric(as.character(ret$chr))
    ret <- ret[order(ret$chr),]
    ret
}

readRawMutFile <- function(f,fill=FALSE) {
    rawColNames <- c("ft","ips","fib","chr","pos","iref","ialt","fref","falt","pval","ref","alt","csq","info","exac","uk1kg")
    message("Reading raw mutation calls from ",f)
    ret <- read.table(f,sep="\t",as.is=TRUE,fill=fill)
    if(fill)
        ret <- ret[-(nrow(ret)),]
    if(ncol(ret)!=length(rawColNames))
        stop("Error in assigning columns")
    colnames(ret) <- rawColNames
    ret
}

getDonorCountVec <- function(dat) {
    d1 <- split(dat,dat$chr)
    lapply(d1,function(l) {
               message("Splitting by position for chr ",l[1,"chr"])
               l1 <- split(l,l$pos)
               donVec <- unlist(lapply(l1,function(l2) { length(table(l2$fib)) }))
               donVec[match(l$pos,names(donVec))]
           })
}

markSegSites <- function(dat) {
    i <- getDonorCountVec(dat=dat)
    i <- unlist(i)
    ### Below is some hacky lines to get the list of donor counts ("i") in the same order as the original data frame ("dat")
    i <- cbind(t(matrix(unlist(strsplit(names(i),"\\.")),nrow=2)),i)
    i <- data.frame(chrx=as.numeric(i[,1]),posx=as.numeric(i[,2]),ndon=as.numeric(i[,3]))
    i <- i[order(i$chr,i$pos),]
    ret <- dat[order(dat$chr,dat$pos),]
    ret <- cbind(ret,i)
    if(sum(ret$chr!=ret$chrx & ret$pos!=ret$posx)>0)
        stop("Error in position ordering for segregating sites count")
    ret <- ret[order(ret$ips,ret$chr,ret$pos),]
    ret <- ret[,!(colnames(ret)%in%c("chrx","posx"))]
    ret
}


markSegSites2 <- function(dat,f) {
    cmd <- paste("./code/preprocess-mutation-calls/countNumberDonors.pl ",f,sep="")
    segSites <- system(cmd,intern=TRUE)
    m <- t(matrix(unlist(strsplit(segSites,"\t")),nrow=3))
    p1 <- paste(m[,1],m[,2],sep=":")
    p2 <- paste(dat[,1],dat[,2],sep=":")
    ret <- rep(1,nrow(dat))
    i <- p2%in%p1
    ## message("Found ",sum(i)," sites segregating in > 1 donor")
    ret[i] <- 2
    ret
}
