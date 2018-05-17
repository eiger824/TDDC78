#!/bin/bash

rm -f output.txt

test ! -f output.txt && touch output.txt

for i in `seq 1 16`; do
	echo "Running with $i threads"
	# Run it and store the result in a file
    OUTPUT=$(OMP_NUM_THREADS=$i ./laplsolv-parallel 2>&1)
	ELAPSED_TIME=$(echo $OUTPUT | cut -d" " -f8)
    ITERATIONS=$(echo $OUTPUT | cut -d" " -f12)
	TEMPERATURE=$(echo $OUTPUT | cut -d" " -f18)
	echo -e "$i\t$ELAPSED_TIME\t$ITERATIONS\t$TEMPERATURE" >> output.txt
done
cat output.txt
