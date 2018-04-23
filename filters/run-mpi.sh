#!/bin/bash

usage()
{
    local name=$(basename $0)
    echo "USAGE: IM={1,2,3,4} $name <filter> <nr-cores>"
    echo "( filter={blurc,thresc} )"
    echo "( by default, 1 node -> 16 cores )"
    echo -e "\nExample: IM=1 $name thresc 4\tApply the threshold filter to image 1 with 4 cores."
}

test -z $IM && 
{
    echo "Error: \`IM\` was not set"
    usage && exit 1
}

PROGRAM=$1

test -z $PROGRAM &&
{
    echo "Error: \`filter\` was not set"
    usage && exit 1
}

if [[ $PROGRAM != "thresc" ]] && [[ $PROGRAM != "blurc" ]]; then 
    echo "Error: filter can only be \"blurc\" or \"thresc\""
    usage
    exit 1
else
    FILESTATS="results-$(eval echo $PROGRAM)-mpi.txt"
fi

NR=$2

test -z $NR &&
{
    echo "Error: \`nr-threads\` was not set"
    usage && exit 1
}

# Now check if programs are built
if [[ ! -f $PROGRAM ]]; then
    echo -e "Please make first:\tmake MPI=1"
    exit 1
fi

ARGS=""
if [[ $PROGRAM == "thresc" ]]; then
    ARGS="../imgs/im$IM.ppm out$(eval echo $IM)-$PROGRAM-$(eval echo $NR)mpinodes.ppm"
else
    ARGS="10 ../imgs/im$IM.ppm out$(eval echo $IM)-$PROGRAM-$(eval echo $NR)mpinodes.ppm"
fi
echo "Running \"$PROGRAM\" on image $IM on $NR MPI cores! (Defaulting ARGS: '$ARGS')"

# Check if our runtime statistics file exists
test -f $FILESTATS ||
{
    echo -e "filter\timg\tcores\telapsed_time" > $FILESTATS
}

# Run the program and filter out the elapsed time
ELAPSED_TIME=$(salloc -N$NR mpprun ./$PROGRAM $ARGS | grep -E 'secs$' | cut -d' ' -f3)
# Prepend the number of threads used
ELAPSED_TIME="$PROGRAM-img$IM-$NR-$ELAPSED_TIME"
echo "$ELAPSED_TIME" | sed -e 's/-/\t/g' >> $FILESTATS
