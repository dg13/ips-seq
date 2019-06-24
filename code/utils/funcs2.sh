checkFileExists() ## Standard file check function
{
    if [ -z $1 ]
    then
	message "$0: File arg not set in checkFileExists"
	exit
    fi
    if [ ! -e $1 ]
    then
	echo "$0: File $1 not found" 1>&2
	exit 1
    fi
    if [ -d $1 ]
    then
	echo "$0: $1 is a directory" 1>&2
	exit 1
    fi
}

checkDirExists() ## Standard dir check function
{
    if [ -z $1 ]
    then
	echo "$0: Argument to checkDirExists is zero length" 1>&2
	exit 1
    fi
    if [ ! -d $1 ]
    then
	echo "$0: Directory $1 not found or exists but is not a directory" 1>&2
	exit 1
    fi
}

getTmpFile()
{
    ARCH=`uname`
    if [ $ARCH == "Darwin" ]
    then
	EXCLUSIVETMPFILENAME123456789=`mktemp -t tmp`
    else
	EXCLUSIVETMPFILENAME123456789=`mktemp`
    fi
    echo "$EXCLUSIVETMPFILENAME123456789"
}

getTmpDir()
{
    ARCH=`uname`
    if [ $ARCH == "Darwin" ]
    then
	EXCLUSIVETMPDIRNAME123456789=`mktemp -d -t tmp`
    else
	EXCLUSIVETMPDIRNAME123456789=`mktemp -d`
    fi
    echo "$EXCLUSIVETMPDIRNAME123456789"
}

checkDirNotEmpty()
{
    N=`ls $1 | wc -l`
    if [ $N -eq 0 ]
    then
	echo "$0: Error - directory \"$1\" is empty" 1>&2
	exit 1
    fi
}

checkOnFarm()
{
    ONFARM=`echo $HOSTNAME | egrep "^(farm|bc|darwin|gen1)" | wc -l`
    if [ $ONFARM -ne 1 ]
    then
	echo "$0: $HOSTNAME does not appear to be on farm" 1>&2
	exit 1
    fi
}

onFarm()
{
    ONFARM=`echo $HOSTNAME | egrep "^(farm|bc|darwin|gen1)" | wc -l`
    if [ $ONFARM -eq 1 ]
    then
	echo 1
    else
	echo 0
    fi
}

checkFileNotEmpty()
{
    if [ ! -s $1 ]
    then
	echo "$0: file $1 is empty" 1>&2
	exit 1
    fi
}

mkdirIfDoesntExist()
{
    if [ -z $1 ]
    then
	message "$0: No argument given to mkdirIfDoesntExist"
	exit 1
    fi
    if [ -f $1 ]
    then
	echo "$0: File $1 exists and is a regular file" 1>&2
	exit 1
    fi
    if [ ! -e $1 ]
    then
	echo "$0: Creating dir $1" 1>&2
	echo "mkdir -p $1"
	## echo "PATH: $PATH"
	mkdir -p $1
    fi
}

oStr() 
{
    JOBN=$1
    QUEUE=$2
    MEM=$3
    INB=`echo $MEM \* 1000000 | bc -l | awk '{ printf "%d\n",$1 }'`
    echo "-M$INB -q $QUEUE -J $JOBN -o farmOut/$JOBN.%J.stdout -e farmOut/$JOBN.%J.stderr"
}

oStr2() 
{
    JOBN=$1
    QUEUE=$2
    MEM=$3
    if [ -z "$JOBN" ] ||  [ -z "$QUEUE" ] ||  [ -z "$MEM" ]
    then
	message "$0: argument not set in oStr2"
	exit
    fi
    INKB=`echo $MEM \* 1000 | bc -l | awk '{ printf "%d\n",$1 }'`
    echo "-M$INKB -q $QUEUE -J $JOBN -o farmOut/$JOBN.%J.stdout -e farmOut/$JOBN.%J.stderr"
}

oStrM()
{
    JOBN=$1
    QUEUE=$2
    NCPU=$3
    MEM=$4
    if [ -z "$JOBN" ] ||  [ -z "$QUEUE" ] ||  [ -z "$MEM" ] || [ -z $NCPU ]
    then
	message "$0: argument not set in oStrM"
	exit
    fi
    INB=`echo $MEM \* 1000000 | bc -l | awk '{ printf "%d\n",$1 }'`
    echo "-M$INB -q $QUEUE -J $JOBN -o farmOut/$JOBN.%J.stdout -e farmOut/$JOBN.%J.stderr -n$NCPU"
}

oStrM2()
{
    JOBN=$1
    QUEUE=$2
    MEM=$3
    NCPU=$4
    if [ -z "$JOBN" ] ||  [ -z "$QUEUE" ] ||  [ -z "$MEM" ] || [ -z $NCPU ]
    then
	message "$0: argument not set in oStrM2"
	exit
    fi
    INKB=`echo $MEM \* 1000 | bc -l | awk '{ printf "%d\n",$1 }'`
    echo "-M$INKB -q $QUEUE -J $JOBN -o farmOut/$JOBN.%J.stdout -e farmOut/$JOBN.%J.stderr -n$NCPU"
}

rStr() 
{
    MEM=$1 ## Memory requirements in Gb
    if [ -z "$MEM" ]
    then
	message "$0: argument not set in rStr"
	exit
    fi
    INKB=`echo $MEM \* 1000 | bc -l | awk '{ printf "%d\n",$1 }'`
    echo "select[mem>$INKB] rusage[mem=$INKB]"
}

rStrM()
{
    MEM=$1 ## Memory requirements in Gb
    INKB=`echo $MEM \* 1000 | bc -l | awk '{ printf "%d\n",$1 }'`
    echo "select[mem>$INKB] rusage[mem=$INKB] span[hosts=1]"
}

isGzip()
{
    ISGZIP=`echo $1 | grep gz | wc -l`
    if [ $ISGZIP -eq 1 ]
    then
	RET=1
    else
	RET=0
    fi
    echo $RET
}

checkFileDoesntExist()
{
    if [ -e $1 ]
    then
	echo "$0: File $1 found" 1>&2
	exit 1
    fi
}

getChrSize()
{
    CHRSTR="chr"$CHR
    CHRSIZE=`awk '$1=="'$CHRSTR'" { print $2 }' $CHROMSIZEFILE`
    echo $CHRSIZE
}

checkClobberArg()
{
    echo "$0: Reading clobber option \"$1\" from CL"
    if [ -z $1 ]
    then
	message "$0: Clobber arg not set (length zero)"
	exit 1
    fi
    if [ $1 != "clobber" ] && [ $1 != "noclobber" ]
    then
	message "$0: Clobber argument $1 not recognised (must be clobber|noclobber)"
	exit 1
    fi
}

checkModeArg()
{
    echo "$0: Reading mode option \"$1\" from CL" 1>&2
    if [ -z $1 ]
    then
	message "$0: Mode arg not set (length zero)"
	exit 1
    fi
    if [ $1 != "debug" ] && [ $1 != "farm" ] && [ $1 != "local" ]
    then
	message "$0: Mode argument $1 not recognised (must be farm|debug)"
	exit 1
    fi
}

checkInPath()
{
    EXE=$1
    command -v $EXE >/dev/null 2>&1 || { echo >&2 "I require $EXE but it's not installed.  Aborting.";exit 1;  }
}

checkDirNotEmpty()
{
    DIR=$1
    if [ ! "$(ls -A $DIR)" ]
    then
	echo "$0: Directory $DIR is empty" 1>&2
	exit 1
    fi
}

isDirEmpty()
{
    DIR=$1
    if [ ! "$(ls -A $DIR)" ]
    then
	RET=0 ## Empty
    else
	RET=1 ## Not empty
    fi
    echo $RET
}

assertIsSetNotEmpty() 
{
    [[ ! ${!1} && ${!1-_} ]] && {
        echo "$0: $1 is not set, aborting." >&2
        exit 1
    }
    : "${!1:? "$1 is empty, aborting."}"
}

message()
{
    echo $1 1>&2
}
