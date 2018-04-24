#!/bin/bash

rm -f thresc-pthreads-results.png
rm -f blurc-pthreads-results.png

if [[ $1 != "pthreads" ]] && [[ $1 != "mpi" ]]; then
    echo "USAGE: $(basename $0) pthreads|mpi"
    exit 1
fi

octave plots-$1.m
