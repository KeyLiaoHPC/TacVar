# Makefile for libtacvar_measure.a
# Usage:
#   make CONF=path/to/tacvar.conf CONSUMER=npb-mpi OUTDIR=build/npb-mpi
#   make CONF=... CONSUMER=lmbench CC=gcc OUTDIR=build/lmbench

MEASURE_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CONF         ?= $(MEASURE_ROOT)/../..//suites/lmbench/tacvar.conf
CONSUMER     ?= lmbench
OUTDIR       ?= $(MEASURE_ROOT)/build/$(CONSUMER)
CC           ?= gcc
CFLAGS       ?= -std=c99 -Wall -Wextra -O2 -D_GNU_SOURCE
PAPI_HOME    ?= /home/hpckey/01-App/papi-7.1.0
PAPI_INC     ?= $(PAPI_HOME)/include
PAPI_LIB     ?= $(PAPI_HOME)/lib
GEN          := $(MEASURE_ROOT)/tools/gen_config.py

.PHONY: all clean gen lib

all: lib

gen: $(OUTDIR)/tacvar_generated_config.h $(OUTDIR)/tacvar_generated_config.mk

$(OUTDIR)/tacvar_generated_config.h $(OUTDIR)/tacvar_generated_config.mk: $(CONF) $(GEN)
	@mkdir -p $(OUTDIR)
	python3 $(GEN) --conf $(CONF) --compiler $(CC) --consumer $(CONSUMER) \
		--outdir $(OUTDIR) --measure-root $(MEASURE_ROOT) \
		--papi-inc $(PAPI_INC) --papi-lib $(PAPI_LIB)

-include $(OUTDIR)/tacvar_generated_config.mk

lib: gen $(OUTDIR)/libtacvar_measure.a

# Object build rules after include so TACVAR_MEASURE_* are available.
$(OUTDIR)/tacvar_measure.o: $(MEASURE_ROOT)/tacvar_measure.c \
		$(OUTDIR)/tacvar_generated_config.h
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(OUTDIR)/tacvar_generated_config.h \
		-c $< -o $@

$(OUTDIR)/tacvar_csv.o: $(MEASURE_ROOT)/tacvar_csv.c \
		$(OUTDIR)/tacvar_generated_config.h
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(OUTDIR)/tacvar_generated_config.h \
		-c $< -o $@

$(OUTDIR)/timer_util.o: $(MEASURE_ROOT)/timers/timer_util.c
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) -c $< -o $@

$(OUTDIR)/counter_perf_event.o: $(MEASURE_ROOT)/counters/counter_perf_event.c \
		$(OUTDIR)/tacvar_generated_config.h
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(OUTDIR)/tacvar_generated_config.h \
		-c $< -o $@

$(OUTDIR)/counter_papi.o: $(MEASURE_ROOT)/counters/counter_papi.c \
		$(OUTDIR)/tacvar_generated_config.h
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(OUTDIR)/tacvar_generated_config.h \
		-c $< -o $@

$(OUTDIR)/counter_rdpmc_x86.o: $(MEASURE_ROOT)/counters/counter_rdpmc_x86.c \
		$(OUTDIR)/tacvar_generated_config.h
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(OUTDIR)/tacvar_generated_config.h \
		-c $< -o $@

$(OUTDIR)/counter_pmu_arm.o: $(MEASURE_ROOT)/counters/counter_pmu_arm.c \
		$(OUTDIR)/tacvar_generated_config.h
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -I$(OUTDIR) \
		-DTACVAR_HAS_GENERATED_CONFIG \
		-include $(OUTDIR)/tacvar_generated_config.h \
		-c $< -o $@

$(OUTDIR)/x86_events.o: $(MEASURE_ROOT)/events/x86_events.c
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -c $< -o $@

$(OUTDIR)/armv8_events.o: $(MEASURE_ROOT)/events/armv8_events.c
	$(CC) $(CFLAGS) $(TACVAR_MEASURE_CFLAGS) -c $< -o $@

$(OUTDIR)/libtacvar_measure.a: $(addprefix $(OUTDIR)/,$(TACVAR_MEASURE_OBJS))
	rm -f $@
	$(AR) cr $@ $^
	ranlib $@ || true

clean:
	rm -rf $(OUTDIR)
