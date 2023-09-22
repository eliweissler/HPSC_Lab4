
# -----------------------------------------------------
# On local Ubuntu
# -----------------------------------------------------
#
# Be sure to use MPI compiler wrappers.
# Whenever compiling an MPI program, you should use the MPI wrappers:
#
#    C - mpicc
#    C++ - mpiCC, mpicxx, mpic++
#    FORTRAN - mpifort, mpif77, mpif90
#
# These wrappers do all of the dirty work for you of making sure that
# all of the appropriate compiler flags, libraries, include directories,
# library directories, etc. are included when you compile your program.


# -----------------------------------------------------
# On Summit:
# -----------------------------------------------------
#
# module load intel
# module load impi

# -----------------------------------------------------
# Set compiler based on platform
# -----------------------------------------------------

# Default compiler is Intel MPI compiler

COMPILER := mpicc
UBUNTU   := ubuntu

# Change if on ubuntu

ifeq ($(strip $(LOGNAME)),scott)
COMPILER := mpic++ 
endif

# -----------------------------------------------------
# Make esPIC
# -----------------------------------------------------


esPIC: esPIC.cpp esPIC.h gauss_seidel.h plotter.h mpiInfo.h particles.h particle_plotter.h LaplacianOnGrid.h
	$(COMPILER) esPIC.cpp -g -o esPIC -std=c++11 -fopenmp

