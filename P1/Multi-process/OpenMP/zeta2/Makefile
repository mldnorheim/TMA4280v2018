# Run ompzeta2 (Q8)
ompzeta2: zeta2.cpp
	g++ -fopenmp zeta2.cpp -lm -o zeta2
	./zeta2 4 100     
      # ./zeta2 nthreads n

# Run error-duration-plot for zeta2 (Q8)
errdurzeta2: errdurzeta2.cpp
	g++ -fopenmp errdurzeta2.cpp -lm -o errdurzeta2
	./errdurzeta2 2

# Run ompmach2 (Q8)
ompmach2: mach2.cpp
	g++ -fopenmp mach2.cpp -lm -o mach2
	./mach2 4 100     
      # ./mach2 nthreads n

# Run error-duration-plot for mach2 (Q8)
errdurmach2: errdurmach2.cpp
	g++ -fopenmp errdurmach2.cpp -lm -o errdurmach2
	./errdurmach2 2
