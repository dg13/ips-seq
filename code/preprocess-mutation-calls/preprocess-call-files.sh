#!/bin/sh -e

## Wrapper to run all the basic filtering for a call set

SCRDIR=$(dirname $0)
if [ ! -e $SCRDIR/../utils/variables.sh ]
then
     echo "$0: Error $SCRDIR/../utils/variables.sh not found" 1>&2
     exit 1
fi
source $SCRDIR/../utils/variables.sh

## Functions
doSort() {
    INF=$1
    OUTF=$2
    PVALTHRESH=0.05
    BUFF=${NCORES}G
    message "$0: Creating call file from $INF weak filtering on pvalue (<$PVALTHRESH), then sorted by line, chr then position, output to $OUTF"
    message "$0:gunzip -c $INF | awk '\$1!=\"#\" && \$5<$PVALTHRESH' | sort --parallel=$NCORES --buffer-size=$BUFF -k 4,4 -k 1,1n -k 2,2n | gzip > $OUTF"
    gunzip -c $INF | awk '$1!="#" && $5<'$PVALTHRESH'' | sort --parallel=$NCORES --buffer-size=$BUFF -k 4,4 -k 1,1n -k 2,2n | gzip > $OUTF
}

doRle() {
    INF=$1
    OUTF=$2
    OUTF=`echo $OUTF | sed 's/\.gz$//'`
    message "$0: Doing RLE for $INF"
    message "$0: R --vanilla --slave --args $INF $OUTF < code/newCalls/makeRleOfPositions-v2.R"
    R --vanilla --slave --args $INF $OUTF < code/newCalls/makeRleOfPositions-v2.R
    message "$0: gzip'ing $OUTF"
    gzip -f $OUTF
}

doMarkDinucs() {
    SORTEDFILE=$1
    RLEFILE=$2
    OUTF=$3
    message "$0: Doing dinuc marking, > 2 mutated nucs removed for $INF"
    message "$0: ./code/newCalls/markDinucs.pl $SORTEDFILE $RLEFILE  | awk '\$NF<3' | gzip > $OUTF"
    ./code/newCalls/markDinucs.pl $SORTEDFILE $RLEFILE  | awk '$NF<3' | gzip > $OUTF ## $NF<3 filters any line where last field is >2, removing multinucleotide runs
}

## Check CL options
N_STANDARD_ARGS=4
echo "$0: Checking CL options" 1>&2
if [ $# -lt $N_STANDARD_ARGS ]
then
    echo "Usage: $0 <rawCallFile> <outDir> <stem> <step>" 1>&2
    echo "$0: Wrapper to run all the basic filtering for the WES data set" 1>&2
    exit
fi
INFILE=$1
checkFileExists $INFILE
MSTROUTDIR=$2
mkdirIfDoesntExist $MSTROUTDIR
STEM=$3
STEP=$4
R=/software/R-3.1.2/bin/R
NCORES=16
if [ $# -gt $N_STANDARD_ARGS ]
then
    NCORES=$5
fi

## Execute
## 1. Sort by line, chromosome and then position
OUTDIR=$MSTROUTDIR/sorted
SORTEDFILE=$OUTDIR/sorted-${STEM}.txt.gz
if [ $STEP -le 0 ]
then
    mkdirIfDoesntExist $OUTDIR
    doSort $INFILE $SORTEDFILE
fi

## 2. Create run length encoding of position diffs to identify runs of positions separated by 1 base
OUTDIR=$MSTROUTDIR/sorted
RLEFILE=$OUTDIR/rl-${STEM}.txt.gz
if [ $STEP -le 1 ]
then
    mkdirIfDoesntExist $OUTDIR
    doRle $SORTEDFILE $RLEFILE
fi

## 3. Mark n-nucleotide runs and remove trinucs
OUTDIR=$MSTROUTDIR/dinucsMarked
DINUCFILE=$OUTDIR/dinucsMarked-${STEM}.txt.gz
if [ $STEP -le 2 ]
then
    mkdirIfDoesntExist $OUTDIR
    doMarkDinucs $SORTEDFILE $RLEFILE $DINUCFILE
fi

## 4. Make FDR call sets
FDRDIR=$MSTROUTDIR/fdrCallSets
if [ $STEP -le 3 ]
then
    FDRTHRESH=$SCRDIR/fdr-${STEM}.txt
    checkFileExists $FDRTHRESH ## Needs to be created by hand
    mkdirIfDoesntExist $FDRDIR
    message "$0: Doing FDR call sets"
    message "$0: ./code/newCalls/makeFdr.sh $DINUCFILE $OUTDIR $THRESHFILE mut-${STEM}"
    ./code/newCalls/makeFdr-v2.sh $DINUCFILE $FDRDIR $FDRTHRESH mut-${STEM}
fi

## 5. Select highest impact mutation in BCSQ field
BCSFDIR=$MSTROUTDIR/bcsqFilter
if [ $STEP -le 4 ]
then
    checkDirExists $FDRDIR
    mkdirIfDoesntExist $BCSFDIR
    message "$0: Running BCSQ field filtering"
    message "$0: ./code/newCalls/filterBcsqField.sh $FDRDIR $OUTDIR mut-${STEM}"
    ./code/newCalls/filterBcsqField-v2.sh $FDRDIR $BCSFDIR mut-${STEM}
fi

## 6. Collapse dinucs
DINUCDIR=$MSTROUTDIR/collapseDinucs
if [ $STEP -le 4 ]
then
    checkDirExists $BCSFDIR
    mkdirIfDoesntExist $DINUCDIR
    message "$0: Running dinucleotide collapsing"
    message "$0: ./code/newCalls/collapseDinucs.sh $BCSFDIR $DINUCDIR mut-${STEM}"
    ./code/newCalls/collapseDinucs-v2.sh $BCSFDIR $DINUCDIR mut-${STEM}
fi    
