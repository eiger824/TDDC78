all: laplsolv-sequential laplsolv-parallel

laplsolv-sequential: laplsolv-sequential.f90
	ifort -fopenmp $^ -o $@

laplsolv-parallel: laplsolv-parallel.f90
	ifort -fopenmp $^ -o $@

clean:
	rm -f *~ laplsolv-parallel laplsolv-sequential
