#!/bin/bash

if [ $(basename $PWD) != "filters" ]; then
    echo "Please run this script from the filters directory."
    exit 1
fi

rm -f thresc-pthreads-results.png
rm -f blurc-pthreads-results.png

if [[ $1 != "pthreads" ]] && [[ $1 != "mpi" ]]; then
    echo "USAGE: $(basename $0) pthreads|mpi"
    exit 1
fi

octave scripts/plots-$1.m
