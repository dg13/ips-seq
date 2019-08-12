inFile1 <- "Data/fitnessImpact/wgs.all.sift-pphen.txt.gz"
inFile2 <- "Data/fitnessImpact/wes.all.sift-pphen.txt.gz"
inFile3 <- "Data/fitnessImpact/hwes.all.sift-pphen.txt.gz"
message("Reading ",inFile1)
p <- read.table("Data/fitnessImpact/wgs.all.sift-pphen.txt.gz",col.names=c("chr","pos","uk10k","exac","vepGene","vepImpact","csqGene","csqImpact","sift","pphen"))
saveRDS(p,file="Data/fitnessImpact/wgs.all.sift-pphen.rds")
message("Reading ",inFile2)
p <- read.table("Data/fitnessImpact/wes.all.sift-pphen.txt.gz",col.names=c("chr","pos","uk10k","exac","vepGene","vepImpact","csqGene","csqImpact","sift","pphen"))
saveRDS(p,file="Data/fitnessImpact/wes.all.sift-pphen.rds")
message("Reading ",inFile3)
p <- read.table("Data/fitnessImpact/hwes.all.sift-pphen.txt.gz",col.names=c("chr","pos","uk10k","exac","vepGene","vepImpact","csqGene","csqImpact","sift","pphen"))
saveRDS(p,file="Data/fitnessImpact/hwes.all.sift-pphen.rds")
