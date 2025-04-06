#TO BE CHANGED BY USER:
EXE=hm_bethe
CC=mpic++
PLAT=gnu
DIREXE=.

# LIBRARIES TO BE INCLUDED
LIB_ED=edipack2_cbinding

FLAG = -O2 -D_MPI

#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################

ifdef LIB_ED
GLOB_INC+=$(shell pkg-config --cflags ${LIB_ED})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_ED})
endif





all:
	@echo ""
	$(call colorecho,"compiling $(EXE) ", 6)
	@echo ""
	$(CC) $(FLAG) $(EXE).cpp -o $(DIREXE)/$(notdir $(EXE)) ${GLOB_INC} ${GLOB_LIB}

clean: 
	@echo "Cleaning:"
	@rm -f *.o *~
	@rm -fv  $(DIREXE)/$(EXE)



#########################################################################
