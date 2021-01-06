MPICC=mpicc
LDFLAGS= -lm -llammps_intelmpi -lfftw3 -openmp -lgslcblas
BLAH=-std=c++0x

INC=-I ~/path_to-lammps/src
LIB=-L ~/path_to-lammps/src
OFLAGS=-o Lammps_drivensystem_run.out
SOURCE=Lammps_drivensystem_run.cpp

all:
	$(MPICC) $(BLAH) $(INC) $(LIB) $(SOURCE) $(LDFLAGS) $(OFLAGS)


