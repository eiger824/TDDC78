FORTCC  = gfortran
CC 		= gcc

all: laplsolv-sequential laplsolv-parallel laplsolv-parallel-new

laplsolv-sequential: laplsolv-sequential.f90
	${FORTCC} -fopenmp $^ -o $@

laplsolv-parallel: laplsolv-parallel.f90
	${FORTCC} -fopenmp $^ -o $@

laplsolv-parallel-new: laplsolv-parallel.c
	${CC} -g -fopenmp -Wall -Wextra -std=c99 $< -o $@

clean:
	rm -f *~ laplsolv-parallel laplsolv-sequential laplsolv-parallel-new output.txt
