# Makefile for MPI dummy library. 
# Must be edited for a specific machine. Does NOT read in 
# the make.def file of NPB 3.4
FC = f90
CC = cc
AR = ar
CP = cp
F08 = def

# Enable if either Cray or IBM: (no such flag for most machines: see wtime.h)
# MACHINE	=	-DCRAY
# MACHINE	=	-DIBM

libmpi.a: mpi_dummy.o mpi_dummy_c.o wtime.o
	$(AR) r libmpi.a mpi_dummy.o mpi_dummy_c.o wtime.o

mpi_dummy.o: mpi_dummy.f90 mpif.h
	$(FC) -c $(FFLAGS) mpi_dummy.f90
# For a Cray C90, try:
#	cf90 -dp -c $(FFLAGS) mpi_dummy.f90
# For an IBM 590, try:
#	xlf90 -c $(FFLAGS) mpi_dummy.f90

mpi_dummy.f90: mpi_dummy_${F08}.f90
	$(CP) -f mpi_dummy_${F08}.f90 $@

mpi_dummy_c.o: mpi_dummy_c.c mpi.h
	$(CC) -c ${MACHINE} -o mpi_dummy_c.o mpi_dummy_c.c

wtime.o: wtime.c
# For most machines or CRAY or IBM
	$(CC) -c ${MACHINE} wtime.c
# For a precise timer on an SGI Power Challenge, try:
#	$(CC) -o wtime.o -c wtime_sgi64.c

test: test.f90
	$(FC) -o test -I. test.f90 -L. -lmpi



clean: 
	- rm -f *~ *.o *.mod
	- rm -f test libmpi.a mpi_dummy.f90
