#!/usr/bin/env bash

## Submit wrapper to create files for pairs of lines from the same donor

. ./code/utils/variables.sh

MODE="debug"
CLOBBER="noclobber"

## Check CL options
N_STANDARD_ARGS=5
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "$0: Submit wrapper to create files for pairs of lines from the same donor" 1>&2
    echo "Usage: $0 <donorList> <wgsFile> <ncores> <outDir> <pvalThresh>" 1>&2
    exit
fi
LISTFILE=$1
WGSFILE=$2
NCORES=$3
OUTDIR=$4
PVAL=$5
mkdirIfDoesntExist $OUTDIR
MAXMEMPERTHREAD=2

## Check optional mode and clobber args
if [ $# -gt `expr $N_STANDARD_ARGS` ]
then
    checkModeArg $6
    MODE=$6
fi

if [ $# -gt `expr $N_STANDARD_ARGS + 1` ]
then
    checkClobberArg $7
    CLOBBER=$7
fi

## Execute
cat $LISTFILE | while read DONOR
do
    QUEUE=normal
    MEM=`echo "$MAXMEMPERTHREAD * $NCORES" | bc`
    JOBNAME="get.${DONOR}.pairedLines"
    OUTFILE="$OUTDIR/paired.${DONOR}.txt.gz"
    NCORES2=`echo $NCORES + 2 | bc -l` ## request more cores, because running shell pipeline which can bump over 100% when running concurrently
    OSTR=$(oStrM2 $JOBNAME $QUEUE $MEM $NCORES2)
    RSTR=$(rStrM $MEM)
    CMD="./code/mutShar/findSharedMuts.sh $WGSFILE $DONOR $NCORES $OUTDIR $PVAL"
    if [ $MODE = "debug" ]
    then
	if [ $CLOBBER = "noclobber" ]
	then
	    if [ -e $OUTFILE ]
	    then
		message "$0: $OUTFILE found. Skipping..."
	    else
		echo "bsub $OSTR -R\"$RSTR\" $CMD"
	    fi
	else
	    echo "bsub $OSTR -R\"$RSTR\" $CMD"
	fi
    else
	if [ $CLOBBER = "noclobber" ]
	then
	    if [ -e $OUTFILE ]
	    then
		message "$0: $OUTFILE found. Skipping..."
	    else
		bsub $OSTR -R"$RSTR" $CMD
	    fi
	else
	    bsub $OSTR -R"$RSTR" $CMD
	fi
    fi
done
