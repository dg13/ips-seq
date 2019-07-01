#!/bin/sh -e

## Wrapper for creation of R data structs from filtered call sets
. ./code/utils/funcs2.sh

MODE="debug"
CLOBBER="noclobber"

## Check CL options
N_STANDARD_ARGS=3
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <inDir> <outDir> <stem>" 1>&2
    echo "$0: Wrapper for creation of R data structs from filtered call sets" 1>&2
    exit
fi
INDIR=$1
OUTDIR=$2
checkDirExists $INDIR
mkdirIfDoesntExist $OUTDIR
STEM=$3

## Check optional mode and clobber args
if [ $# -gt `expr $N_STANDARD_ARGS` ]
then
    checkModeArg $4
    MODE=$4
fi

if [ $# -gt `expr $N_STANDARD_ARGS + 1` ]
then
    checkClobberArg $5
    CLOBBER=$5
fi

## Execute
for FDR in 1 5 10
do
    for SUBDIR in bcsqFilter collapseDinucs
    do
	QUEUE=normal
	MEM=10
	JOBNAME="makeRda.$STEM.$FDR.$SUBDIR.$STEM"
	OSTR=$(oStr2 $JOBNAME $QUEUE $MEM)
	RSTR=$(rStr $MEM)
	INFILE="$INDIR/$SUBDIR/mut-${STEM}-fdr${FDR}.txt.gz"
	if [ ! -e $INFILE ]
	then
	    message "$0: Warning input file $INFILE is missing, skipping this submission"
	    SAVEMODE=$MODE
	    MODE="debug"
	fi
	OUTFILE="$OUTDIR/mut-${STEM}-fdr${FDR}.${SUBDIR}.rda"
	CMD="./code/preprocess-mutation-calls/makeNewCallsDataStruct.sh $INFILE $OUTFILE"
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
	if [ ! -e $INFILE ]
	then
	    MODE=$SAVEMODE
	fi
    done
done
