# Standalone tools that live under NPB-MPI/common (not NAS benchmark objects).
# From this directory:  make
# From NPB3.4-MPI/:     make probe_counter_freq

GCC ?= gcc
MPICC ?= mpicc
CFLAGS ?= -Wall -Wextra -std=c99 -O2
LDFLAGS ?= -lm

.PHONY: all probe_counter_freq clean

all: probe_counter_freq

probe_counter_freq: probe_counter_freq.x probe_counter_freq-mpi.x

probe_counter_freq.x: probe_counter_freq.c
	$(GCC) $(CFLAGS) -D_GNU_SOURCE -o $@ $< $(LDFLAGS)

probe_counter_freq-mpi.x: probe_counter_freq.c
	$(MPICC) $(CFLAGS) -D_GNU_SOURCE -DUSE_MPI -o $@ $< $(LDFLAGS)

clean:
	rm -f probe_counter_freq.x probe_counter_freq-mpi.x
