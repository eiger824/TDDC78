#!/bin/bash

#
# Filename:			get_latest.sh
#
# Author:			Santiago Pagola
# Brief:			Secure-Copies files from NSC to the local computer
# Last modified:	sön 15 apr 2018 13:52:50 CEST
#

test -z $1 &&
{
    echo "USAGE: USERNAME=<user> $(basename $0) FILE"
    exit 1
}

FILE=$1

test -z $USERNAME &&
{
    echo "USERNAME is empty"
    echo "USAGE: USERNAME=<user> $(basename $0) FILE"
    exit 2
}

scp $USERNAME@triolith.nsc.liu.se:~/TDDC78/filters/$FILE .
