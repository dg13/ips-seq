#!/bin/sh -e

## Code to collapse dinucs into a single record based on deleteriousness

. ./code/utils/funcs2.sh

## Check CL options
N_STANDARD_ARGS=3
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <inDir> <outDir> <stem>" 1>&2
    echo "$0: Code to collapse dinucs into a single record based on deleteriousness" 1>&2
    exit
fi
INDIR=$1
OUTDIR=$2
STEM=$3
mkdirIfDoesntExist $OUTDIR
export PERL5LIB=$PERL5LIB:`pwd`/code

## Execute
message "$0: Filtering FDR1 set -> ${OUTDIR}/${STEM}-fdr1.txt.gz"
message "$0: ./code/preprocess-mutation-calls/collapseDinucs.pl ${INDIR}/${STEM}-fdr1.txt.gz | gzip > ${OUTDIR}/${STEM}-fdr1.txt.gz"
./code/preprocess-mutation-calls/collapseDinucs.pl ${INDIR}/${STEM}-fdr1.txt.gz | gzip > ${OUTDIR}/${STEM}-fdr1.txt.gz
message "$0: Filtering FDR5 set -> ${OUTDIR}/${STEM}-fdr5.txt.gz"
message "$0: ./code/preprocess-mutation-calls/collapseDinucs.pl ${INDIR}/${STEM}-fdr5.txt.gz | gzip > ${OUTDIR}/${STEM}-fdr5.txt.gz"
./code/preprocess-mutation-calls/collapseDinucs.pl ${INDIR}/${STEM}-fdr5.txt.gz | gzip > ${OUTDIR}/${STEM}-fdr5.txt.gz
message "$0: Filtering FDR10 set -> ${OUTDIR}/${STEM}-fdr10.txt.gz"
message "$0: ./code/preprocess-mutation-calls/collapseDinucs.pl ${INDIR}/${STEM}-fdr10.txt.gz | gzip > ${OUTDIR}/${STEM}-fdr10.txt.gz"
./code/preprocess-mutation-calls/collapseDinucs.pl ${INDIR}/${STEM}-fdr10.txt.gz | gzip > ${OUTDIR}/${STEM}-fdr10.txt.gz
