#!/bin/sh -e

## Wrapper for making rdas
. ./code/utils/funcs2.sh

## Check CL options
N_STANDARD_ARGS=2
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <inFile> <outFile>" 1>&2
    echo "$0: Wrapper for making rdas" 1>&2
    exit
fi
INFILE=$1
OUTFILE=$2
checkFileExists $INFILE

## Execute
message "$0: Rscript code/preprocess-mutation-calls/makeNewCallsDataStruct.R $INFILE $OUTFILE"
Rscript code/preprocess-mutation-calls/makeNewCallsDataStruct.R $INFILE $OUTFILE
