#!/bin/sh -e

## Code to put all the shared mutation locations into a single file

SCRDIR=$(dirname $0)
if [ ! -e $SCRDIR/../utils/variables.sh ]
then
     echo "$0: Error $SCRDIR/../utils/variables.sh not found" 1>&2
     exit 1
fi
source $SCRDIR/../utils/variables.sh

## Check CL options
N_STANDARD_ARGS=2
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <inDir> <outFile>" 1>&2
    echo "$0: Code to put all the shared mutation locations into a single file" 1>&2
    exit
fi
INDIR=$1
checkDirExists $INDIR
OUTFILE=$2

## Execute
for FILE in `ls $INDIR/*.txt.gz`;
do
    DONOR=`echo $( basename $FILE ) | tr '\\.' ' ' | awk '{ print $1 }'`;
    message "Processing $DONOR";
    gunzip -c $FILE | awk -vOFS="\t" '{ print $1,$2,"'"$DONOR"'" }';
done | gzip -f > $OUTFILE
