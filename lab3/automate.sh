#!/bin/bash

rm -f output.txt

test ! -f output.txt && touch output.txt

echo -e "Threads\tTime\t\t\tIter\tTemp" >> output.txt
for i in 1 2 4 8 9 10 11 12 14 16 32; do
	echo "Running with $i threads"
	# Run it and store the result in a file
    OUTPUT=$(OMP_NUM_THREADS=$i ./laplsolv-parallel 2>&1)
	ELAPSED_TIME=$(echo $OUTPUT | cut -d" " -f2)
	ITERATIONS=$(echo $OUTPUT | cut -d" " -f6)
	TEMPERATURE=$(echo $OUTPUT | cut -d" " -f12)
	echo -e "$i\t$ELAPSED_TIME\t$ITERATIONS\t$TEMPERATURE\n" >> output.txt
done
cat output.txt
