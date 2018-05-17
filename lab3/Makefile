FORTCC  = gfortran

all: laplsolv-sequential laplsolv-parallel
everything: laplsolv-sequential laplsolv-parallel run plots

laplsolv-sequential: laplsolv-sequential.f90
	${FORTCC} -fopenmp $^ -o $@

laplsolv-parallel: laplsolv-parallel.f90
	${FORTCC} -fopenmp $^ -o $@

run:
	./automate.sh

plots: output.txt
	octave plots/plots-openmp.m 2>&1 > /dev/null
	eog exec-times.png &

clean:
	rm -f *~ laplsolv-parallel laplsolv-sequential *.txt *.png
