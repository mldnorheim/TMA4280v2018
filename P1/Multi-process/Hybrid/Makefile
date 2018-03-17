# Run hybrid
hybrid: hybrid.cpp
	mpiCC hybrid.cpp -lm -o hybrid -fopenmp
	mpirun -np 4 hybrid 4 100 # include argument for which of zeta/mach we want to run?
#                    nthreads n
