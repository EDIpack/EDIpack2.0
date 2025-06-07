FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

TSTDIR=src
TSTBIN=bin

HAS_INEQ:=$(shell pkg-config --variable=ineq edipack)
GLOB_INC:=$(shell pkg-config --cflags scifor edipack)

OBJS=ASSERTING.o COMMON.o
.SUFFIXES: .f90
.PHONY: clean

ifeq ($(HAS_INEQ),true)
all: bin ineq_normal_normal
else
all: bin
endif

bin: \
normal_normal normal_superc normal_nonsu2 \
hybrid_normal hybrid_superc hybrid_nonsu2 \
replica_normal replica_superc replica_nonsu2 \
general_normal general_superc general_nonsu2


normal_normal:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/NORMAL/NORMAL
	echo NORMAL/NORMAL >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/NORMAL_NORMAL
	@echo ""

normal_superc:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/NORMAL/SUPERC
	echo NORMAL/SUPERC >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/NORMAL_SUPERC
	@echo ""

normal_nonsu2:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/NORMAL/NONSU2
	echo NORMAL/NONSU2 >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/NORMAL_NONSU2
	@echo ""

hybrid_normal:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/HYBRID/NORMAL
	echo HYBRID/NORMAL >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/HYBRID_NORMAL
	@echo ""

hybrid_superc:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/HYBRID/SUPERC
	echo HYBRID/SUPERC >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/HYBRID_SUPERC
	@echo ""

hybrid_nonsu2:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/HYBRID/NONSU2
	echo HYBRID/NONSU2 >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/HYBRID_NONSU2
	@echo ""

replica_normal:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/REPLICA/NORMAL
	echo REPLICA/NORMAL >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/REPLICA_NORMAL
	@echo ""

replica_superc:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/REPLICA/SUPERC
	echo REPLICA/SUPERC >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/REPLICA_SUPERC
	@echo ""

replica_nonsu2:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/REPLICA/NONSU2
	echo REPLICA/NONSU2 >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/REPLICA_NONSU2
	@echo ""

general_normal:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/GENERAL/NORMAL
	echo GENERAL/NORMAL >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/GENERAL_NORMAL
	@echo ""

general_superc:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/GENERAL/SUPERC
	echo GENERAL/SUPERC >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/GENERAL_SUPERC
	@echo ""

general_nonsu2:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/GENERAL/NONSU2
	echo GENERAL/NONSU2 >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/GENERAL_NONSU2
	@echo ""

ineq_normal_normal:
	@echo "Compile:" $@
	mkdir -p $(TSTBIN)/INEQ/NORMAL/NORMAL
	echo INEQ/NORMAL/NORMAL >> $(TSTBIN)/list_dir
	make -C $(TSTDIR)/INEQ_NORMAL_NORMAL
	@echo ""


test: 
	./test.sh



.f90.o:
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}



clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
