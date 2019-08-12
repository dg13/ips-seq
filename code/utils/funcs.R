removeDinucSecondBase <- function(dat) {
    xd <- dat[dat$dinuc,]
    xd <- xd[seq(1,nrow(xd),2),]
    ret1 <- dat[!dat$dinuc,] ## Add the SNVs back in
    ret <- rbind(ret1,xd)
    ret <- ret[order(ret$ips,ret$chr,ret$pos),]
    ret
}

convertToFriendly <- function(n) {
    gsub("HPSI.*-","",n)
}

convertToDonor <- function(n) {
    gsub("_[0-9]*$","",gsub("HPSI.*-","",n))
}

getVarDetMatrix <- function(dat,version="mine",field="info") {
    if(field=="info") {
        ret <- getVarDetMatrixVep(dat=dat,version=version)
    }
    if(field=="csq") {
        ret <- getVarDetMatrixCsq(dat=dat,version=version)
    }
    ret
}

getVarDetMatrixVep <- function(dat,version="mine") {
    fields <- strsplit(dat[,"info"],"\\|")
    len <- unlist(lapply(fields,length))
    tab <- table(len)
    if(length(tab)>1)
        stop("Error in info string formatting - fields have variable numbers of fields")
    ncol <- as.numeric(names(tab))
    x1 <- t(matrix(unlist(fields),nrow=ncol))
    if(version=="mine") {
        colnames(x1) <- c("varDetails","varConsequences","impact","geneSymbol","ensGene","featureType","ensTrans","bioType","exon","intron","hgvsc","hgvsp","cdnaPos","cdsPos","protPos","aa","codons","existVar","dist","strand","flags","geneSymbolSrc","hgnc")
    } else if(version=="petr") {
          colnames(x1) <- c("varDetails","varConsequences","impact","geneSymbol","ensGene","featureType","ensTrans","bioType","exon","intron","hgvsc","hgvsp","cdnaPos","cdsPos","protPos","aa","codons","existVar","dist","strand","geneSymbolSrc","hgnc","sift")
    }
    as.data.frame(x1,stringsAsFactors=FALSE)
}

getVarDetMatrixCsq <- function(dat,version=1) {
    dat[dat$csq==".","csq"] <- "NA:NA:.:NA:NA"
    dat <- selectMostDeleteriousFromMultiFieldRows(dat=dat,version=version) ## Deals with rows that have multi-consequence outcomes (e.g. missense and splice in the same gene)
    fields <- strsplit(dat$csq,":")
    len <- unlist(lapply(fields,length))
    tab <- table(len)
    if(length(tab)>1)
        stop("Error in info string formatting - fields have variable numbers of fields")
    ncol <- as.numeric(names(tab))
    m <- t(matrix(unlist(fields),nrow=ncol))
    data.frame(gene=m[,1],tx=m[,2],csq=m[,3],aa=m[,4],pos=m[,5],stringsAsFactors=FALSE)
}

## Functions to mark MNPs

markDinucs <- function(dat) {
    message("Marking dinucs")
    dat <- cbind(dat,FALSE)
    colnames(dat)[ncol(dat)] <- "dinuc"
    dat <- dat[order(dat$ips,dat$chr,dat$pos),]
    rl <- makeRunLengthEncodingStruct(dat)
    pos <- rl[rl[,2]==1 & rl[,3]==1,1]
    message("Found ",length(pos)," dinucs")
    pos <- sort(c(pos,pos+1))
    dat[pos,"dinuc"] <- TRUE ## this works even if pos is length==0
    rl <- makeRunLengthEncodingStruct(dat)
    if(sum(rl[,2]>3 & rl[,3]==1)>0)
        stop("Error >quadnucleotide change found")
    pos <- rl[rl[,2]==3 & rl[,3]==1,1]  ## marks quadnucs
    if(length(pos)>0) {
        message("Removing ",length(pos)," quadnuc runs")
        pos <- sort(c(pos-2,pos-1,pos,pos+1))
        dat <- dat[-(pos),]
    }
    rl <- makeRunLengthEncodingStruct(dat)    
    pos <- rl[rl[,2]==2 & rl[,3]==1,1]  ## marks trinucs
    if(length(pos)>0) {
        message("Removing ",length(pos)," trinuc runs")
        pos <- sort(c(pos-1,pos,pos+1))
        dat <- dat[-(pos),]
    }
    dat
}

## Functions to select the most deleterious site in dinucleotides or multi consequence variants

## Wrappers
selectMostDeleteriousFromMultiFieldRows <- function(dat,version=1) {
    ## Wrapper function to go between different csq formats
    message(version)
    if(version==1) {
        ret <- selectMostDeleteriousFromMultiFieldRowsV1(dat)
    }
    if(version==2) {
        ret <- selectMostDeleteriousFromMultiFieldRowsV2(dat)
    }
    ret
}

selectMostDeleteriousFromMultiFieldRowsV1 <- function(dat) {
    fields <- strsplit(dat[,"csq"],":")
    len <- unlist(lapply(fields,length))
    i <- which(len!=5)
    if(length(i)>0) {
        outcomes <- sapply(i,chooseMostDeleteriousVersion,dat=dat)
        dat[i,"csq"] <- outcomes
    } else {
          message("No multiple outcome fields found")
      }
    dat
}

selectMostDeleteriousFromMultiFieldRowsV2 <- function(dat) {
    fields <- strsplit(dat[,"csq"],",")
    len <- unlist(lapply(fields,length))
    i <- which(len>1)
    if(length(i)>0) {
        outcomes <- sapply(i,chooseMostDeleteriousVersion,dat=dat,version=2)
        dat[i,"csq"] <- outcomes
    } else {
          message("No multiple outcome fields found")
      }
    dat
}

chooseMostDeleteriousVersion <- function(j,dat,version=1) {
    ## Wrapper function to go between different csq formats
    if(version==1) {
        ret <- chooseMostDeleteriousVersionV1(j,dat)
    }
    if(version==2) {
        ret <- chooseMostDeleteriousVersionV2(j,dat)
    }
}

chooseMostDeleteriousVersionV1 <- function(j,dat,field1="\\|",field2=":") {
    d1 <- unlist(strsplit(dat[j,"csq"],field1))
    m <- t(matrix(unlist(strsplit(d1,field2)),ncol=length(d1)))
    m <- data.frame(gene=m[,1],tx=m[,2],csq=m[,3],aa=m[,4],pos=m[,5],stringsAsFactors=FALSE)
    csqCodes <- getCsqCodes(m$csq)
    m <- m[order(csqCodes,decreasing=TRUE),]
    paste(m[1,],collapse=":")
}

handleNoncoding <- function(d1,field1,field2) {
    i <- grep("(non_coding|intron|splice|utr)",d1)
    ret <- data.frame(gene=NA,tx=NA,csq=NA,aa=NA,pos=NA,biotype=NA,stringsAsFactors=FALSE)
    if(length(i)>0) {
        d1 <- d1[i]
        m <- t(matrix(unlist(strsplit(d1,field2)),ncol=length(d1)))
        ret <- data.frame(gene=m[,2],tx=m[,3],csq=m[,1],aa=NA,pos=NA,biotype=m[,4],stringsAsFactors=FALSE)
    }
    ret
}

handleCoding <- function(d1,field1,field2) {
    i <- grep("(synonymous|stop_retained|missense|start_lost|stop_lost|stop_gained)",d1)
    ret <- data.frame(gene=NA,tx=NA,csq=NA,aa=NA,pos=NA,biotype=NA,stringsAsFactors=FALSE)
    if(length(i)>0) {
        d1 <- d1[i]
        m <- t(matrix(unlist(strsplit(d1,field2)),ncol=length(d1)))
        ret <- data.frame(gene=m[,2],tx=m[,3],csq=m[,1],aa=m[,6],pos=m[,7],biotype=m[,4],stringsAsFactors=FALSE)
    }
    ret
}

handleCompound <- function(d1,field1,field2) {
    i <- grep("@",d1)
    d1 <- d1[i]
    ret <- data.frame(gene=NA,tx=NA,csq=NA,aa=NA,pos=NA,biotype=NA,stringsAsFactors=FALSE)
    if(length(i)>0) {
        m <- t(matrix(unlist(strsplit(d1,field2)),ncol=length(d1)))
        ret <- data.frame(gene=NA,tx=NA,csq=NA,aa=NA,pos=m[,1],biotype=NA,stringsAsFactors=FALSE)
    }
    ret
}

chooseMostDeleteriousVersionV2 <- function(j,dat,field1=",",field2="\\|") {
    print(paste("Here",j))
    d1 <- unlist(strsplit(dat[j,"csq"],field1))
    print(dat[j,])
    m1 <- handleNoncoding(d1,field1=field1,field2=field2)
    m2 <- handleCoding(d1,field1=field1,field2=field2)
    m3 <- handleCompound(d1,field1=field1,field2=field2)
    m <- rbind(m1,m2,m3)
    print(paste("Here 1"))
    m <- m[!is.na(m$gene),]
    print(m$csq)
    m <- m[m$biotype=="protein_coding" | m$biotype=="lincRNA",]
    print(paste("Here 2"))
    csqCodes <- getCsqCodesV2(m$csq)
    m <- m[order(csqCodes,decreasing=TRUE),]
    paste(m[1,],collapse="|")
}

## m <- t(matrix(unlist(strsplit(d1,field2)),ncol=length(d1)))
## if(ncol(m)==4) { 
##     m <- data.frame(gene=m[,2],tx=NA,csq=m[,1],aa=NA,pos=NA,stringsAsFactors=FALSE)
## } else if(ncol(m)==7) {
##     m <- data.frame(gene=m[,2],tx=m[,3],csq=m[,1],aa=m[,6],pos=m[,7],stringsAsFactors=FALSE)
## } else if(ncol(m)==1) {
##     m <- data.frame(gene=NA,tx=NA,csq=NA,aa=NA,pos=m[,1],stringsAsFactors=FALSE)
## } else {
##     stop("Syntax error")
## }

getCsqCodesV2 <- function(vec) {
    vec[is.na(vec)] <- "."
    i <- grep("^(\\.|non_coding|intron|synonymous|stop_retained|3_prime_utr|5_prime_utr|splice_region|missense|splice_donor|splice_acceptor|start_lost|stop_lost|stop_gained)$",vec,invert=TRUE)
    if(length(i)>0)
        stop("Unknown field - \"",vec[i[1]],"\" found in ",length(i)," locations of CSQ field")
    vec[is.na(vec)] <- 0 ## NA set to 0 impact
    vec <- gsub("\\.",0,vec)
    vec <- gsub("non_coding",0,vec)
    vec <- gsub("intron",1,vec)
    vec <- gsub("synonymous",2,vec)
    vec <- gsub("stop_retained",2,vec)
    vec <- gsub("3_prime_utr",3,vec)
    vec <- gsub("5_prime_utr",4,vec)
    vec <- gsub("splice_region",5,vec)
    vec <- gsub("missense",6,vec)
    vec <- gsub("splice_donor",7.5,vec)
    vec <- gsub("splice_acceptor",7.5,vec)
    vec <- gsub("start_lost",9,vec)
    vec <- gsub("stop_lost",10,vec)
    vec <- gsub("stop_gained",11,vec)
    as.numeric(vec)
}

getCsqCodes <- function(vec) {
    i <- grep("^(\\.|synonymous_variant|missense_variant|stop_lost|start_stop_splice|stop_gained)$",vec,invert=TRUE)
    if(length(i)>0)
        stop("Unknown field - \"",vec[i[1]],"\" found in ",length(i)," locations of CSQ field")
    vec[is.na(vec)] <- 0
    vec <- gsub("\\.",0,vec)
    vec <- gsub("synonymous_variant",1,vec)
    vec <- gsub("missense_variant",2,vec)
    vec <- gsub("stop_lost",3,vec)
    vec <- gsub("start_stop_splice",4,vec)
    vec <- gsub("stop_gained",5,vec)
    as.numeric(vec)
}

selectBaseFromDinuc <- function(i,csqCodes) {
    i1 <- csqCodes[i]
    i2 <- csqCodes[i+1]
    if(i1<i2) i+1 else i;
}



####################################################################################################
## Old functions kept to not break legacy code - most up to date functions above

checkIsDinuc <- function(i,dat) {
    ## if(!i%%100)
    ##     message("chr ",dat[1,"chr"]," row ",i)
    ## 
    sum(dat$pos==dat[i,"pos"]-1 | dat$pos==dat[i,"pos"]+1)>0
}

checkDinucsChrByChr <- function(dat) {
    message("Marking dinucs in ",dat[1,"ips"])
    d <- split(dat,dat$chr) ## split into chrs
    lapply(d,markTriDinucs)
    ## unlist(lapply(d,function(l) { sapply(1:nrow(l),checkIsDinuc,dat=l) }))
}

checkDinucsSampleBySample <- function(dat) {
    d <- split(dat,dat$ips)
    unlist(lapply(d,checkDinucsChrByChr))
}

makeRunLengthEncodingStruct <- function(dat) {
    rl <- rle(diff(dat$pos)==1) ## RLE of boolean vector indicating locations of single position differences in position - i.e. di/trinucs
    rl1 <- cbind(rl$lengths,rl$values)
    rl1 <- cbind(cumsum(rl1[,1]),rl1) ## Cumsum gives row positions of runs
    rl1
}

markTriDinucs <- function(dat) {
    rl <- rle(diff(dat$pos)==1) ## RLE of boolean vector indicating locations of single position differences in position - i.e. di/trinucs
    rl1 <- cbind(rl$lengths,rl$values)
    rl1 <- cbind(cumsum(rl1[,1]),rl1) ## Cumsum gives index positions of runs of TRUE
    if(sum(rl1[,2]>2 & rl1[,3]==1)>0)
        stop("Error - > trinucleotide change found in sample ",dat$ips[1])
    pos <- rl1[rl1[,2]==2 & rl1[,3]==1,1]  ## marks trinucs
    if(length(pos)>0) {
        message("Removing ",length(pos)," trinuc runs in sample ",dat$ips[1],", chromosome ",dat$chr[1])
        pos <- sort(c(pos-1,pos,pos+1))
        dat <- dat[-(pos),]
    }
    pos <- rl1[rl1[,2]==1 & rl1[,3]==1,1]
    ##     message("Found ",length(pos)," dinucs in sample ",dat$ips[1],", chromosome ",dat$chr[1])
    pos <- sort(c(pos,pos+1))
    dat[pos,"dinuc"] <- TRUE ## this works even if pos is length==0
    dat
}

setOrder <- function(orderBy="alpha") {
    ret <- NA
    if(orderBy=="mut") {
        ret <- cbind(rep(c("A","C","G","T"),each=4),"[",c(paste(rep("C",48),c(rep(">A",16),rep(">G",16),rep(">T",16)),sep=""),paste(rep("T",48),c(rep(">A",16),rep(">C",16),rep(">G",16)),sep="")),"]",rep(c("A","C","G","T")))
    }
    if(orderBy=="alpha") {
        ret <- cbind(rep(c("A","C","G","T"),each=24),"[",rep(c(paste(rep("C",12),c(rep(">A",4),rep(">G",4),rep(">T",4)),sep=""),paste(rep("T",12),c(rep(">A",4),rep(">C",4),rep(">G",4)),sep="")),2),"]",rep(c("A","C","G","T")))
    }
    apply(ret,1,paste,collapse="")
}

orderCounts <- function(vec,orderBy="alpha") {
    cats <- setOrder(orderBy=orderBy)
    match(cats,vec)
}


## 

getShared <- function(l) {
    message("Processing ",l[1,"fib"])
    p <- paste(l$chr,l$pos,sep=":")
    ret <- rep(FALSE,length(p))
    ret[p%in%names(which(table(p)>1))] <- TRUE
    ret
}

addShared <- function(mut,minCounts=0,removeSingleLines=FALSE) {
    n1 <- mut$falt + mut$fref
    n2 <- mut$ialt + mut$iref
    mut <- mut[n1 & n2 >= minCounts,]
    dat <- split(mut,mut$fib)
    n <- unlist(lapply(dat,function(l) { length(unique(l$ips)) }))
    oneLine <- dat[which(n==1)]
    message("Found ",length(oneLine)," donors with a single line")
    dat <- dat[which(n==2)]
    message("Found ",length(dat)," donors with a single line")
    shar <- lapply(dat,getShared)
    dat <- do.call(rbind,dat)
    dat$shar <- unlist(shar)
    if(!removeSingleLines) {
        oneLine <- do.call(rbind,oneLine)
        oneLine$shar <- FALSE
        dat <- rbind(dat,oneLine)
    }
    dat
}

checkDinucOrdering <- function(dat) {
    i <- seq(1,nrow(dat)-1,2)
    iplus1 <- seq(2,nrow(dat),2)
    stopifnot(all(dat[i,"pos"]+1==dat[iplus1,"pos"]))
}


callAdjacentMuts <- function(dat) {
    ## Function to subset data frame to only adjacent mutations
    ind <- which(diff(dat[,2])==1)
    ind <- sort(c(ind,ind+1))
    message("Found ",length(ind)," adjacent mutations")
    dat[ind,]
}

callNonAdjacentMuts <- function(dat) {
    ## Function to subset data frame to only adjacent mutations
    ind <- which(diff(dat[,2])==1)
    ind <- sort(c(ind,ind+1))
    ind <- which(!(1:nrow(dat)%in%ind))
    message("Found ",length(ind)," non-adjacent mutations")
    dat[ind,]
}

filterMutTable <- function(dat) {
    message("Found ",length(table(dat$fib))," donors")
    message("Found ",length(table(dat$ips))," lines")
    message("Removing ",sum(dat$chr==23)," X chromosome mutations")
    dat <- dat[dat$chr<23,] ## Remove X chromosome
    message("Removing ",sum(dat$dinuc)/2," dinucleotide mutations")
    dat <- dat[!dat$dinuc,]
    i <- grepl(",",d$alt)
    message("Removing ",sum(i)," multi-allelics")
    dat <- dat[!i,]
    tab <- table(dat$ips)
    message("Removing ",sum(tab<minMuts)," donors with < ",minMuts," mutations")
    nam <- names(tab[tab>minMuts])
    dat <- dat[dat$ips%in%nam,]
    message("Retaining ",length(table(dat$fib))," donors")
    message("Retaining ",length(table(dat$ips))," lines")
    dat
}

appendMutationContext <- function(dat,fasta,format="cosmic",ignoreMissing=FALSE) {
    s <- paste(">",dat$chr,":",dat$pos,"-",dat$pos,sep="")
    if(sum(!(s%in%fasta[,1]))>0) {
        if(!ignoreMissing)
            stop("Error - sequence(s): \n",paste(s[which(!(s%in%fasta[,1]))],collapse="\n"),"\nnot found in ",faFile)
        else {
            message("**Warning**\nIgnoring ",sum(!(s%in%fasta[,1]))," records missing from ",faFile)
            toKeep <- s%in%fasta[,1]
            s <- s[toKeep]
            dat <- dat[toKeep,]
        }
    }
    fasta <- fasta[match(s,fasta[,1]),]
    dat <- data.frame(dat,trinuc=fasta[,2],stringsAsFactors=FALSE)
    dat <- data.frame(dat,pat="[N>N]",stringsAsFactors=FALSE)
    m <- t(matrix(unlist(strsplit(dat$trinuc,"")),nrow=3))
    if(sum(m[,2]!=dat$ref)>0)
        stop("Error in ref base assignment")
    m <- makeRevCompMat(m)
    dat[(dat$ref=="C" & dat$alt=="A")|(dat$ref=="G" & dat$alt=="T"),"pat"] <- "[C>A]"
    dat[(dat$ref=="C" & dat$alt=="G")|(dat$ref=="G" & dat$alt=="C"),"pat"] <- "[C>G]"
    dat[(dat$ref=="C" & dat$alt=="T")|(dat$ref=="G" & dat$alt=="A"),"pat"] <- "[C>T]"
    dat[(dat$ref=="T" & dat$alt=="A")|(dat$ref=="A" & dat$alt=="T"),"pat"] <- "[T>A]"
    dat[(dat$ref=="T" & dat$alt=="C")|(dat$ref=="A" & dat$alt=="G"),"pat"] <- "[T>C]"
    dat[(dat$ref=="T" & dat$alt=="G")|(dat$ref=="A" & dat$alt=="C"),"pat"] <- "[T>G]"
    if(sum(dat$pat=="[N>N]")>0) {
        print(tab[tab$pat=="[N>N]",])
        stop("Error - some mutations (multi-allelic, N, ?) not ATGC->ATGC")
    }
    data.frame(dat,trinucPat=paste(m[,1],dat$pat,m[,3],sep=""),stringsAsFactors=FALSE)
}

makeRevCompMat <- function(mat) {
    if(sum(grepl("[^ATGC]",mat)>0))
       stop("Error - non ATGC bases found")
    t(apply(mat,1,revCompRow))
}

revCompRow <- function(r) {
    ret <- NA
    if(r[2]!="C" & r[2]!="T") {
        ret <- rev(unlist(lapply(r,revComp)))
    } else {
          ret <- r
      }
    ret
}

revComp <- function(base) {
    if(base=="A") { ret <- "T" }
    if(base=="C") { ret <- "G" }
    if(base=="G") { ret <- "C" }
    if(base=="T") { ret <- "A" }
    ret
}
