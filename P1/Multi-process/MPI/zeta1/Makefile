# Run mpizeta1 (Q5) - mpirun -np 8 zeta1   # 8 = nprocs = 2^N where N = 3
mpizeta1: zeta1.cpp
	mpiCC -o zeta1 zeta1.cpp
	mpirun -np 4 zeta1   # 2 = nprocs = 2^N where N = 1

# Run error-duration-plot for zeta1
errdurzeta1: errdurzeta1.cpp
	mpiCC -o errdurzeta1 errdurzeta1.cpp -lm
	mpirun -np 8 errdurzeta1

# Run mpimach1 (Q5)
mpimach1: mach1.cpp
	mpiCC -o mach1 mach1.cpp
	mpirun -np 2 mach1   # 2 = nprocs = 2^N where N = 1
	mpirun -np 8 mach1   # 8 = nprocs = 2^N where N = 3

# Run error-duration-plot for mach1 
errdurmach1: errdurmach1.cpp
	mpiCC -o errdurmach1 errdurmach1.cpp -lm
	mpirun -np 8 errdurmach1

# Run plotting of error and timings (Q5)
#plottings: 

# Run mpireduc (Q7)
mpireduc: reduc.cpp
	mpiCC -o reduc reduc.cpp -lm
	mpirun -np 4 reduc 0 0
	mpirun -np 4 reduc 0 1
	mpirun -np 4 reduc 1 0
	mpirun -np 4 reduc 1 1

# Run mpireduc00 (Q7)
mpireduc00: reduc.cpp
	mpiCC -o reduc reduc.cpp -lm
	mpirun -np 4 reduc 0 0

# Run mpireduc01 (Q7)
mpireduc01: reduc.cpp
	mpiCC -o reduc reduc.cpp -lm
	mpirun -np 4 reduc 0 1

# Run mpireduc10 (Q7)
mpireduc10: $(reduc.cpp)
	mpiCC -o reduc reduc.cpp -lm
	mpirun -np 4 reduc 1 0

# Run mpireduc11 (Q7)
mpireduc11: reduc.cpp
	mpiCC -o reduc reduc.cpp -lm
	mpirun -np 4 reduc 1 1

