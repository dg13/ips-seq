#!/bin/sh

. ./code/utils/variables.sh

./code/mutShar/concatenateSharedCalls.sh Data/mutShar/paired-wes/ Data/mutShar/paired-wes.all.txt.gz
./code/mutShar/concatenateSharedCalls.sh Data/mutShar/paired-hwes/ Data/mutShar/paired-hwes.all.txt.gz
./code/mutShar/concatenateSharedCalls.sh Data/mutShar/paired-wgs/ Data/mutShar/paired-wgs.all.txt.gz
$RSCRIPT code/mutShar/addSharingInfoToRds.R Data/mut-files/rdas/mut-wgs-fdr5.bcsqFilter.rda Data/mutShar/paired-wgs.all.txt.gz Data/mut-files/rdas/mut-wgs-fdr5.bcsqFilter.wsharinfo.rda
$RSCRIPT code/mutShar/addSharingInfoToRds.R Data/mut-files/rdas/mut-wgs-fdr1.bcsqFilter.rda Data/mutShar/paired-wgs.all.txt.gz Data/mut-files/rdas/mut-wgs-fdr1.bcsqFilter.wsharinfo.rda
$RSCRIPT code/mutShar/addSharingInfoToRds.R Data/mut-files/rdas/mut-wes-fdr5.bcsqFilter.rda Data/mutShar/paired-wes.all.txt.gz Data/mut-files/rdas/mut-wes-fdr5.bcsqFilter.wsharinfo.rda
$RSCRIPT code/mutShar/addSharingInfoToRds.R Data/mut-files/rdas/mut-hwes-fdr5.bcsqFilter.rda Data/mutShar/paired-hwes.all.txt.gz Data/mut-files/rdas/mut-hwes-fdr5.bcsqFilter.wsharinfo.rda
