
SHELL   = /bin/bash
CC     ?= gcc
FUTC   ?= futhark
CFLAGS ?= -std=c99 -fopenmp -O3 -march=native

C_LIBS    = -std=c99 -lm
MC_LIBS   = -std=c99 -lm -lpthread
OCL_LIBS  = -lm -lOpenCL
CUDA_LIBS = -lm -lcuda -lcudart -lnvrtc

FUTTYPE=c
LIBS=$(MC_LIBS)

all: futtop_c

futtop: futtop.c libmultigrid.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

futtopExplicit: futtopExplicit.c libexplicit.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

futtopNonlinear: futtopNonlinear.c libnonlinear.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

futtopNonlinearMMA: futtopNonlinearMMA.c libnonlinear.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

computeCompliance: computeCompliance.c libcompliance.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

libmultigrid.c: libmultigrid.fut src/*.fut
	$(FUTC) $(FUTTYPE) --library $<

libexplicit.c: libexplicit.fut src/*.fut
	$(FUTC) $(FUTTYPE) --library $<

libnonlinear.c: libnonlinear.fut src/*.fut
	$(FUTC) $(FUTTYPE) --library $<

libcompliance.c: libcompliance.fut src/*.fut
	$(FUTC) $(FUTTYPE) --library $<

clean:
	rm -f $$(ls src/*.fut | cut -f1 -d'.')
	rm -f futtop futtop_library.c futtop_library.h
	rm -f src/*.c

bench:
	$(FUTC) bench --backend=opencl src/*.fut
