FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

HAS_INEQ:=$(shell pkg-config --variable=ineq edipack2)
GLOB_INC:=$(shell pkg-config --cflags scifor edipack2)

OBJS=ASSERTING.o COMMON.o

ifeq ($(HAS_INEQ),true)
all: $(OBJS) NORMAL HYBRID REPLICA GENERAL INEQ
else
all: $(OBJS) NORMAL HYBRID REPLICA GENERAL
endif

NORMAL HYBRID REPLICA GENERAL: 
	make -C ./$@_NORMAL
	make -C ./$@_SUPERC
	make -C ./$@_NONSU2

INEQ:
	make -C ./$@_NORMAL_NORMAL

.f90.o:
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o 

