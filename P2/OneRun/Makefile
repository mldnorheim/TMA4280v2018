# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
MPICC=mpicxx
MPIFC=mpif90
TARGET=poisson
FSRC=fst.f
CSRC=poisson.cpp
LIBS=-lm
OPENMP?=0
OPTFLAGS=-O3
#OPENMP=0
ifeq ($(strip $(OPENMP)),1)
OPTFLAGS+=-fopenmp
endif
CFLAGS=$(OPTFLAGS)
FFLAGS=$(OPTFLAGS)


all:	$(SRC)
	$(MPIFC) $(FFLAGS) -c fst.f
	$(MPICC) $(CFLAGS) -c poisson.cpp
	$(MPICC) $(CFLAGS) -o $(TARGET) fst.o poisson.o $(LIBS) 

clean:
	rm -f *.o $(TARGET)
