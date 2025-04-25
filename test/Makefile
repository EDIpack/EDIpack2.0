TSTDIR=src
TSTBIN=bin

HAS_INEQ:=$(shell pkg-config --variable=ineq edipack2)

ifeq ($(HAS_INEQ),true)
all: bin ineq
else
all: bin
endif
	make -C $(TSTDIR)

bin:
	mkdir -p $(TSTBIN)
####### BATH_TYPE FOLDERS
	mkdir -p $(TSTBIN)/NORMAL
	mkdir -p $(TSTBIN)/HYBRID
	mkdir -p $(TSTBIN)/REPLICA
	mkdir -p $(TSTBIN)/GENERAL
####### ED_MODE FOLDERS
	mkdir -p $(TSTBIN)/NORMAL/NORMAL
	mkdir -p $(TSTBIN)/NORMAL/SUPERC
	mkdir -p $(TSTBIN)/NORMAL/NONSU2	
	mkdir -p $(TSTBIN)/HYBRID/NORMAL
	mkdir -p $(TSTBIN)/HYBRID/SUPERC
	mkdir -p $(TSTBIN)/HYBRID/NONSU2
	mkdir -p $(TSTBIN)/REPLICA/NORMAL
	mkdir -p $(TSTBIN)/REPLICA/SUPERC
	mkdir -p $(TSTBIN)/REPLICA/NONSU2	
	mkdir -p $(TSTBIN)/GENERAL/NORMAL
	mkdir -p $(TSTBIN)/GENERAL/SUPERC
	mkdir -p $(TSTBIN)/GENERAL/NONSU2
####### INSERT FOLDER TO BE TESTED IN LIST_DIR
	echo NORMAL/NORMAL >> $(TSTBIN)/list_dir
	echo NORMAL/SUPERC >> $(TSTBIN)/list_dir
	echo NORMAL/NONSU2 >> $(TSTBIN)/list_dir
	echo HYBRID/NORMAL  >> $(TSTBIN)/list_dir
	echo HYBRID/SUPERC  >> $(TSTBIN)/list_dir
	echo HYBRID/NONSU2  >> $(TSTBIN)/list_dir
	echo REPLICA/NORMAL >> $(TSTBIN)/list_dir
	echo REPLICA/SUPERC >> $(TSTBIN)/list_dir
	echo REPLICA/NONSU2 >> $(TSTBIN)/list_dir
	echo GENERAL/NORMAL >> $(TSTBIN)/list_dir
	echo GENERAL/SUPERC >> $(TSTBIN)/list_dir
	echo GENERAL/NONSU2 >> $(TSTBIN)/list_dir


ineq:
	mkdir -p $(TSTBIN)
	mkdir -p $(TSTBIN)/INEQ
####### INEQ BATH_TYPE FOLDERS
	mkdir -p $(TSTBIN)/INEQ/NORMAL
	mkdir -p $(TSTBIN)/INEQ/NORMAL/NORMAL
####### INSERT FOLDER TO BE TESTED IN LIST_DIR
	echo INEQ/NORMAL/NORMAL >> $(TSTBIN)/list_dir


test: 
	./test.sh



clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
