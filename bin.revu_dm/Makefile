#Makefile 

###############################################################################
# Include definitions 
###############################################################################
include ../include.mk_dm
include ./objects.mk.$(RAMS_VERSION)

###############################################################################
# Compiler commands
###############################################################################
INCLUDES  = $(UTILS_INCS) $(PAR_INCS) $(HDF5_INCS)
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)
C_COMMAND = $(C_COMP) -c $(C_OPTS) -D$(CMACH) $(INCLUDES) \
            $(PAR_DEFS) $(HDF5_DEFS)

###############################################################################
# File extension rules
###############################################################################
.SUFFIXES:
.SUFFIXES: .F90 .f90 .f .c .o .a

# f90 rule.
.f90.a:
	@echo ""
	$(F_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.f90=.o)
	rm -f $(<F:.f90=.o)

# F90 rule.
.F90.a:
	@echo ""
	$(F_COMMAND) -D$(CMACH) $<
	$(ARCHIVE) $@ $(<F:.F90=.o)
	rm -f $(<F:.F90=.o)

# c rule.
.c.a:
	@echo ""
	$(C_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.c=.o)
	rm -f $(<F:.c=.o)

###############################################################################
# Define archive and executable names
###############################################################################
EXE=$(BASE)-$(RAMS_VERSION)
ARC=$(BASE)-$(RAMS_VERSION).a

# Define targets.

all: $(EXE)

$(EXE): $(ARC) $(MAIN) FORCE
	@echo ""
	$(F_COMP) -o $(EXE) $(MAIN_OBJ) $(LOADER_OPTS) $(ARC) \
	$(HDF5_LIBS) $(LIBS) $(PAR_LIBS)
	rm -f $(MAIN:.f90=.f) *.o
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""
#	if [ ! -d ../run ]; then mkdir ../run; fi;
#	ln -fs `pwd`/$(EXE) ../run/$(BASE)
#	cp $(EXE) ../run

$(MAIN): FORCE
	@echo ""
	$(F_COMMAND) $@

$(ARC): $(OBJ)

FORCE:

check: FORCE
	@echo ""
	check

install:
	@echo ""
#	if [ ! -d ../run ]; then mkdir ../run; fi;
#	ln -fs `pwd`/$(EXE) ../run/$(BASE)
#	cp $(EXE) ../run
	@echo ""

clean:
	@echo ""
	rm -f $(ARC) $(EXE) $(BASE) *.o *.mod *.f
#	rm -f ../run/$(BASE) ../run/$(EXE)
	@echo ""

# Include dependencies (if this does not work, try using gnu make or copy the
# contents of the file dep_file here).

#include dep_file.mk
