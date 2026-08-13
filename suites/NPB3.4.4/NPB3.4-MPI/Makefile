SHELL=/bin/sh
CLASS=U
SUBTYPE=
VERSION=
F08=def
MFLAGS=CLASS=$(CLASS) F08=$(F08)
SFILE=config/suite.def

default: header
	@ sys/print_instructions

BT: bt
bt: header
	cd BT; $(MAKE) $(MFLAGS) SUBTYPE=$(SUBTYPE) VERSION=$(VERSION)

SP: sp
sp: header
	cd SP; $(MAKE) $(MFLAGS)

LU: lu
lu: header
	cd LU; $(MAKE) $(MFLAGS) VERSION=$(VERSION)

MG: mg
mg: header
	cd MG; $(MAKE) $(MFLAGS)

FT: ft
ft: header
	cd FT; $(MAKE) $(MFLAGS)

IS: is
is: header
	cd IS; $(MAKE) $(MFLAGS)

CG: cg
cg: header
	cd CG; $(MAKE) $(MFLAGS)

EP: ep
ep: header
	cd EP; $(MAKE) $(MFLAGS)

DT: dt
dt: header
	cd DT; $(MAKE) $(MFLAGS)

# Awk script courtesy cmg@cray.com, modified by Haoqiang Jin
suite:
	@ awk -f sys/suite.awk SMAKE=$(MAKE) $(SFILE) | $(SHELL)


# It would be nice to make clean in each subdirectory (the targets
# are defined) but on a really clean system this will won't work
# because those makefiles need config/make.def
clean:
	- rm -f core *~ */core */*~
	- rm -f */*.o */*.mod */*.obj */*.exe */npbparams.h */mpinpb.*
	- rm -f MPI_dummy/test MPI_dummy/libmpi.a MPI_dummy/mpi_dummy.f90
	- rm -f sys/setparams sys/makesuite sys/setparams.h
	- rm -f btio.*.out*

# After editing tacvar.conf — drop generated measure lib + timer glue, then rebuild.
tacvar_clean:
	- rm -f common/tacvar_npb.o common/tacvar_npb_tf.o common/c_timers.o
	- rm -rf config/tacvar_build_*

filt:
	mkdir -p bin
	cc -O2 -std=c99 -Wall -o bin/filt.x common/filt.c -lm

# nspg uses TACVAR timer macros from generated config (same OUTDIR as benchmarks).
TACVAR_ROOT    ?= $(abspath ../../../src/measure)
TACVAR_CONF    ?= tacvar.conf
TACVAR_OUTDIR  ?= config/tacvar_build_$(shell uname -m)
PAPI_HOME      ?= /home/hpckey/01-App/papi-7.1.0
NSPG_MPICC     ?= mpicc

nspg:
	mkdir -p bin $(TACVAR_OUTDIR)
	python3 $(TACVAR_ROOT)/tools/gen_config.py --conf $(TACVAR_CONF) \
		--compiler "$(NSPG_MPICC)" --consumer npb-mpi --outdir $(TACVAR_OUTDIR) \
		--measure-root $(TACVAR_ROOT) \
		--papi-inc $(PAPI_HOME)/include --papi-lib $(PAPI_HOME)/lib
	$(NSPG_MPICC) -O2 -std=c99 -Wall -D_GNU_SOURCE \
		`sed -n 's/^TACVAR_MEASURE_CFLAGS := //p' $(TACVAR_OUTDIR)/tacvar_generated_config.mk` \
		-I$(TACVAR_OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(TACVAR_OUTDIR)/tacvar_generated_config.h \
		-o bin/test_nspg.x common/test_nspg.c \
		`sed -n 's/^TACVAR_MEASURE_LDFLAGS := //p' $(TACVAR_OUTDIR)/tacvar_generated_config.mk`

# NAS: drop site config. Next `make <bench>` copies make.def.template → make.def.
veryclean: clean tacvar_clean
	- rm -f config/make.def config/suite.def
	- rm -f bin/sp.* bin/lu.* bin/mg.* bin/ft.* bin/bt.* bin/is.* 
	- rm -f bin/ep.* bin/cg.* bin/dt.* bin/filt.x bin/test_nspg.x

header:
	@ sys/print_header

