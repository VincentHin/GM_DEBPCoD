# Generic Makefile for the compilation of any EBT program on either laptop or
# the lisa cluster. Both the dynamic and bifurcation version of the model can
# be compiled. These programs will be called "<model>" and "<model>_bif",
# respectively.
#
# Possible targets to make:
#
# 	make <model>
#	make <model_bif
#	make all			(makes all possible EBT target programs)
#	make both			(makes first <model> and <model>_bif)
#	make cleanebt			(cleans up all EBT programs)
#	make cleancvf			(cleans up all CVF files)
#	make allclean			(combines previous two)
#
#	make PROBLEM=model1 both	(makes <model1> and <model1>_bif)
#	make PROBLEM=model1 cleanone	(cleans up <model1> and <model1>_bif)
#
#==============================================================================
# Tayloring
# The following settings might need adaptation for specific uses

# EXTRAINCLUDES = MeasureBifstats.c

# EXTRAINCLUDEDIRS = $(HOME)/programs/include $(HOME)/programs/lib /opt/local/include
EXTRAINCLUDEDIRS = /opt/local/include

# EXTRALIBDIRS = $(HOME)/programs/lib /opt/local/lib
EXTRALIBDIRS = /opt/local/lib

# EXTRALIBS = fftw3

TMPDIR   = /tmp

# Determining the names of all EBT programs and CVF files in the current directory
# All .c files in the current directory are considered EBT programs

ALLCFILES = $(wildcard *.c)
ALLPROBLEMS = $(patsubst %.c,%,$(ALLCFILES))
ALLCVF = $(wildcard *.cvf)

# Alternatively, if not all .c files in the current directory are EBT programs
# or not all CVF file should be clean by 'make allclean' and 'make cleancvf'
# specify the targets here:

# ALLPROBLEMS = <model1> <model2>
# ALLCVF = <run1.cvf> <run2.cvf>

#==============================================================================
# Default configuration flags

# If present use icc, otherwise gcc

ICC=$(shell which icc)
ifeq ("$(ICC)", "")
INTEL_CC ?= 0
else
INTEL_CC ?= 1
endif

# Compile for debugging when DEBUG="debug" or DEBUG="1"

MAKEDEBUG=0
ifeq ("$(DEBUG)", "debug")
MAKEDEBUG=1
INTEL_CC=0
endif
ifeq ("$(DEBUG)", "1")
MAKEDEBUG=1
INTEL_CC=0
endif

#==============================================================================
# Find the ebtclean program

EBTCLEAN = $(shell which ebtclean)
ifeq ("$(EBTCLEAN)", "")
SYSTEM = $(shell uname)
ifeq ("$(SYSTEM)", "Darwin")			# This is Mac OS specific
EBTCLEAN = $(patsubst %/Resources,%/MacOs/ebtclean, $(EBTPATH))
endif
endif

#==============================================================================

# Compilation settings for Intel C++ compiler
ifeq ("$(INTEL_CC)", "1")
CC	 = icc
WARN     = -w1 -Wcheck

# Optimization settings for Intel C++ compiler
# OPTFLAGS = -O2
# OPTFLAGS = -O3 -ipo -static -xHost
# OPTFLAGS = -O3 -ipo -static -axavx,sse4.2,sse3
# Appropriate definition to run on all cpu3 cores of LISA
OPTFLAGS = -O3 -ipo -static -xCORE-AVX-I

# Linker settings for Intel C++ compiler
LD.c     = xild -r

#  Compilation settings for GCC compiler
else
CC	?= clang
WARN     = -Wall -Wpointer-arith -Wcast-qual -Wcast-align
OPTFLAGS = -O3 -march=native
LD.c     = ld -r
endif

# Debugging settings (mostly to be used on laptop via "debug on")
ifeq ("$(MAKEDEBUG)", "1")
CC       = /opt/local/bin/clang
CCFLAGS  = -g -DDEBUG=1 $(WARN) $(SPECDEFS)
TMPDIR   = .
CLEANCMD = echo ' '; echo Object files of module library preserved for debugging; echo ' '
else
CCFLAGS  = $(OPTFLAGS)  $(WARN) $(SPECDEFS)
CLEANCMD = echo ' '; echo Cleaning object files of module library; rm -f $(EBTOBJS); echo ' '
endif

# Generic compilation and linking flags
EBTDIR      = $(EBTPATH)/integrator
INCLUDEDIRS = . $(EBTDIR) $(EBTDIR)/Odesolvers $(EXTRAINCLUDEDIRS)
CFLAGS	    = $(CCFLAGS) -DPROBLEMFILE="<$(PROBLEM).h>" $(patsubst %, -I%, $(INCLUDEDIRS))
LDFLAGS     = $(patsubst %, -L%, $(EXTRALIBDIRS))
LIBS        = -lm $(patsubst %, -l%, $(EXTRALIBS))

# Names of the EBT modules linked into the program-independent library
EBTMODS  =  ebtinit ebtmain ebtcohrt ebttint ebtutils ebtstop
EBTOBJS  =  $(TMPDIR)/ebtinit.o $(TMPDIR)/ebtmain.o  $(TMPDIR)/ebtcohrt.o \
	    $(TMPDIR)/ebttint.o $(TMPDIR)/ebtutils.o $(TMPDIR)/ebtstop.o

#==============================================================================

# The following targets are valid if PROBLEM is not defined

ifeq ("$(PROBLEM)", "")
all:
	for I in $(ALLPROBLEMS) ; do $(MAKE) PROBLEM=$${I} both ; done

both:
	$(MAKE) PROBLEM=$(word 1,$(ALLPROBLEMS)) both

allclean: cleanebt cleancvf

cleanebt:
	@echo "Cleaning up everything...."
	@for I in $(ALLPROBLEMS) ; do $(MAKE) PROBLEM=$${I} cleanone ; done

cleancvf:
	@echo "Cleaning up $(ALLCVF)"
	@for I in $(patsubst %.cvf, %, $(ALLCVF)) ; do $(EBTCLEAN) -f $${I} >/dev/null ; done

# Re-invoke make but now with PROBLEM defined and the same target

$(ALLPROBLEMS) $(patsubst %,%_bif,$(ALLPROBLEMS))::
	$(MAKE) PROBLEM=$(subst _bif,,$@) $@

#==============================================================================

# The following targets are valid if PROBLEM is defined

else
PROBLEMUC=$(shell echo $(PROBLEM) | tr "[:lower:]" "[:upper:]")
PROBLEM_BIF = $(PROBLEM)_bif

# Generate the list of additional dependency files

EXTRADEPENDS = $(foreach FILE, $(EXTRAINCLUDES), $(shell ls -1 $(patsubst %, %/$(FILE), $(INCLUDEDIRS)) 2>/dev/null | head -1))

#==============================================================================

# The dependencies of the executables

both: $(PROBLEM) $(PROBLEM_BIF)

$(PROBLEM):     SPECDEFS = -DBIFURCATION=0 -DPROGRAMNAME=$(PROBLEMUC)
$(PROBLEM_BIF): SPECDEFS = -DBIFURCATION=1 -DPROGRAMNAME=$(PROBLEMUC)

$(PROBLEM): $(PROBLEM).o $(PROBLEM).lib.o
	$(LINK.c) $(LDFLAGS) -o $@ $(PROBLEM).o     $(PROBLEM).lib.o     $(LIBS)

$(PROBLEM_BIF): $(PROBLEM_BIF).o $(PROBLEM_BIF).lib.o
	$(LINK.c) $(LDFLAGS) -o $@ $(PROBLEM_BIF).o $(PROBLEM_BIF).lib.o $(LIBS)

#==============================================================================

# The dependencies of the problem-specific object files

$(PROBLEM).o: $(PROBLEM).c $(PROBLEM).h $(EXTRADEPENDS)
	$(COMPILE.c) -o $@ $(PROBLEM).c

$(PROBLEM_BIF).o: $(PROBLEM).c $(PROBLEM).h $(EXTRADEPENDS)
	$(COMPILE.c) -o $@ $(PROBLEM).c

#==============================================================================

# The dependencies of the problem object library file

$(PROBLEM).lib.o: $(PROBLEM).h
	@${MAKE} -s -k methtest 2>/dev/null
	for I in $(EBTMODS) ; do $(COMPILE.c) -o $(TMPDIR)/$${I}.o $(EBTDIR)/$${I}.c ; done
	$(LD.c) -o $@ $(EBTOBJS)
	@${CLEANCMD}

$(PROBLEM_BIF).lib.o: $(PROBLEM).h
	@${MAKE} -s -k methtest 2>/dev/null
	for I in $(EBTMODS) ; do $(COMPILE.c) -o $(TMPDIR)/$${I}.o $(EBTDIR)/$${I}.c ; done
	$(LD.c) -o $@ $(EBTOBJS)
	@${CLEANCMD}

#==============================================================================

# The dependencies of some additional targets

cleanone:
	@echo "Cleaning up $(PROBLEM) and $(PROBLEM_BIF)"
	@rm -f  $(PROBLEM)      $(PROBLEM).o        $(PROBLEM).lib.o     $(PROBLEM)module.*
	@rm -f  $(PROBLEM_BIF)  $(PROBLEM_BIF).o    $(PROBLEM_BIF).lib.o $(PROBLEM_BIF)module.*
	@rm -Rf $(PROBLEM).dSYM $(PROBLEM_BIF).dSYM
	@for I in $(EBTMODS) ; do rm -f ./$${I}.o $(TMPDIR)/$${I}.o ; done

methtest:
	@echo ' '
	@$(COMPILE.c) -o /dev/null $(EBTDIR)/methtest.c 2>&1 | grep '#error' |\
		sed 's/.* #error //'
	@echo ' '

endif

#==============================================================================
