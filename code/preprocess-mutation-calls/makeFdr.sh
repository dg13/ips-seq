#!/bin/sh -e

## Code to sort the raw calls and produce FDR thresholded sets
. ./code/utils/funcs2.sh

## Check CL options
N_STANDARD_ARGS=4
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <inFile> <outDir> <fdrFile> <stem>" 1>&2
    echo "$0: Code to sort the raw calls and produce FDR thresholded sets" 1>&2
    exit
fi
INFILE=$1
OUTDIR=$2
FDR=$3
STEM=$4
checkFileExists $FDR
checkFileExists $INFILE
mkdirIfDoesntExist $OUTDIR

## Execute
message "$0: Creating FDR1 set"
message "$0: ./code/preprocess-mutation-calls/makeFdrCallset.pl $INFILE $FDR --fdr 1 | gzip > $OUTDIR/${STEM}-fdr1.txt.gz"
./code/preprocess-mutation-calls/makeFdrCallset.pl $INFILE $FDR --fdr 1 | gzip > $OUTDIR/${STEM}-fdr1.txt.gz
message "$0: Creating FDR5 set"
message "$0: ./code/preprocess-mutation-calls/makeFdrCallset.pl $INFILE $FDR --fdr 5 | gzip > $OUTDIR/${STEM}-fdr5.txt.gz"
./code/preprocess-mutation-calls/makeFdrCallset.pl $INFILE $FDR --fdr 5 | gzip > $OUTDIR/${STEM}-fdr5.txt.gz
message "$0: Creating FDR10 set"
message "$0: ./code/preprocess-mutation-calls/makeFdrCallset.pl $INFILE $FDR --fdr 10 | gzip > $OUTDIR/${STEM}-fdr10.txt.gz"
./code/preprocess-mutation-calls/makeFdrCallset.pl $INFILE $FDR --fdr 10 | gzip > $OUTDIR/${STEM}-fdr10.txt.gz
