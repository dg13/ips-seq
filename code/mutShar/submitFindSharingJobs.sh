#!/usr/bin/env bash

. ./code/utils/variables.sh

FDR=1
message "$0: Running sharing comparison at FDR threshold $FDR (i.e. mutation has to be < than this pvalue in *either* line)"
./code/mutShar/doAllFindSharedMuts.sh Data/mutShar/wes-2donors.txt Data/mut-files/raw/exomes.ft.txt.gz 4 Data/mutShar/paired-wes `awk '$1=="'"$FDR"'" { print $2 }' Data/mut-files/fdr-thresholds/fdr-wes.txt` farm clobber
./code/mutShar/doAllFindSharedMuts.sh Data/mutShar/hwes-2donors.txt Data/mut-files/raw/high-vs-low-exomes.ft.txt.gz 4 Data/mutShar/paired-hwes `awk '$1=="'"$FDR"'" { print $2 }' Data/mut-files/fdr-thresholds/fdr-hwes.txt` farm clobber
./code/mutShar/doAllFindSharedMuts.sh Data/mutShar/wgs-2donors.txt Data/mut-files/raw/wgs.396.ft.txt.gz 4 Data/mutShar/paired-wgs `awk '$1=="'"$FDR"'" { print $2 }' Data/mut-files/fdr-thresholds/fdr-wgs.txt` farm clobber
