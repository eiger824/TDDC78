#!/bin/bash

die()
{
    echo "$@"
    echo "USAGE: $0 <nproc> <progname> <args...>"
    exit 1
}

test $# -ge 2 || die "Not enough input arguments"

# First arg: number of processes
echo $1 | grep -qE '^[0-9]*$'
test $? -eq 0 || die "Wrong process number"
NPROC=$1

shift

# Second arg: does program exist?
test -f $1 || die "Program does not exist"

# Finally: run the program
echo "Running '$@' under \`mpirun' wrapper"
mpirun -np $NPROC $@

