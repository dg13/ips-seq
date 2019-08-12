#!/usr/bin/env bash

## Code to create data structs for both lines from the same donor, mutations filtered only on whether they are significant in one or other line

. ./code/utils/variables.sh

## Check CL options
N_STANDARD_ARGS=5
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <wgsRawCallFile> <donorName> <nCores> <outDir> <pvalThresh>" 1>&2
    echo "$0: Code to create data structs for both lines from the same donor, mutations filtered only on whether they are significant in one or other line" 1>&2
    exit
fi
INFILE=$1
checkFileExists $INFILE
DONOR=$2
NCORES=$3
OUTDIR=$4
PVAL=$5

## Execute
BUFFER=6G
message "$0: gunzip -c $INFILE | grep $DONOR | sort --parallel=$NCORES --buffer-size=$BUFFER  -k 1,1 -k 2,2n |  ./code/mutShar/findShared.pl --pval=$PVAL stdin | gzip > $OUTDIR/${DONOR}.fdr1.txt.gz"
gunzip -c $INFILE | grep $DONOR | sort --parallel=$NCORES --buffer-size=$BUFFER  -k 1,1 -k 2,2n |  ./code/mutShar/findShared.pl --pval=$PVAL stdin | gzip > $OUTDIR/${DONOR}.txt.gz
