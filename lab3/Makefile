FORTCC  = gfortran

all: laplsolv-sequential laplsolv-parallel run plots

laplsolv-sequential: laplsolv-sequential.f90
	${FORTCC} -fopenmp $^ -o $@

laplsolv-parallel: laplsolv-parallel.f90
	${FORTCC} -fopenmp $^ -o $@

run:
	./automate.sh

plots:
	octave plots/plots-openmp.m

clean:
	rm -f *~ laplsolv-parallel laplsolv-sequential laplsolv-parallel-new output.txt
