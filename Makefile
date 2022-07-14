# Makefile for compiling Multipack.
#
# Usages:
#	make
#	make install
#	make doc
#	make MODULES="<list of modules>"
#
# (Rudimentary) Original by Travis Oliphant
# Rewritten by Pearu Peterson <pearu@ioc.ee>
#
# $Revision: 1.22 $
# $Date: 1999/06/14 13:24:25 $

PREFIX:=$(shell python -c "import sys;print sys.prefix")

# Users can modify below:
FC=f77
CC=cc
FFLAGS=-O
FSHARED=-shared 
OPT=-O3
CFLAGS=$(OPT) -I$(PREFIX)/include/python1.5/

# List of modules to be included for compilation (if not specified in
# command line):
ifeq (,$(MODULES))
MODULES:=minpack odepack quadpack fitpack
endif

# User provided libraries:
LINKLIBS:=
    # by default: use only Multipack provided libraries
    # '!' will force skipping Multipack library. For example,
    # '!linpack_lite' means that 'linpack_lite' will not be
    # compiled. Instead, it is assumed that user has provided the
    # corresponing replacement, say 'linpack'. 
    # Note also that user replacement should follow the skipped
    # library. Libraries (such as blas, lapack, etc) that are used by
    # others, should be always latest.
    # Examples:
#LINKLIBS:=minpack odepack quadpack ddierckx \
#	!lapack !linpack_lite !blas linpack lapack blas_intel
#LINKLIBS:=!blas blas_intel

# Where to look for user provided libraries (they are searched first):
LINKLIBSDIR:=
# Example:
LINKLIBSDIR:=/numeric/netlib/lib

# Installation stuff:
INSTALLDIR = $(PREFIX)/lib/python1.5/site-packages/multipack
# Note that 'make install' will try to create a file multipack.pth to
# the directory ../$(INSTALLDIR)
INSTALL = install -c
INSTALLDATA = install -c -m 644

#**************************************************
# The front end user should not modify what follows
#**************************************************
# Run setmodules.py (it solves dependences):
SETMODS:=$(shell ./setmodules.py $(MODULES) : $(LINKLIBS))

# Existing modules:
MODULES:=$(patsubst m.%,%,$(filter m.%,$(SETMODS)))

# The corresponding python files (to be imported to Multipack.py):
PYFILES:=$(patsubst p.%,%,$(filter p.%,$(SETMODS)))
PYFILENAMES:=$(patsubst %.py,%,$(PYFILES))

# Multipack provided libraries that need to be compiled:
MLINKLIBS:=$(patsubst l.%,%,$(filter l.%,$(SETMODS)))

# All libraries to be used in linking:
LLINKLIBS:=$(patsubst ll.%,%,$(filter ll.%,$(SETMODS)))

# Where to look for Multipack provided libraries:
LINKLIBSDIR+=./lib

# Shared libraries of the modules:
MODSO:=$(patsubst %,_%module.so,$(MODULES))

# ./setmodules.py will create `m_multipackmodule.c'
all: Multipack debug

Multipack: $(MODSO) Multipack.py
Multipack.py: _Multipack.py
	cat _Multipack.py > Multipack.py # header of Multipack.py
	echo -e "$(patsubst %,\nfrom % import *,$(PYFILENAMES))" >> Multipack.py
$(patsubst %,lib/m_%module.c,$(MODULES)): lib/m_%module.c : __%.c multipack.h
	./setmodules.py -m -q $* > /dev/null
	mv m_$*module.c lib
$(MODSO): _%module.so : lib/m_%module.o $(MLINKLIBS)
	$(FC) $(FFLAGS) -shared -o $@ lib/m_$*module.o \
	$(patsubst %,-L%,$(LINKLIBSDIR)) \
	$(patsubst %,-l%,$(LLINKLIBS))
	@rm -f Multipack.py
# Compile Multipack provided libraries:
.PHONY: $(MLINKLIBS) doc
$(MLINKLIBS):
	@echo "Compiling library: $(@)"
	cd $(@) && make FC=$(FC) FFLAGS="$(FFLAGS)" LIB=$(@)
	mv $(@)/lib$(@).a lib

# If special treatment is necessary, define them here:
blas:
	@echo "Compiling library: blas"
	cd lapack/BLAS/SRC && make
	mv lapack/$(@).a lib/lib$(@).a
linpack_lite:
	@echo "Compiling library: linpack"
	cd linpack && make FC=$(FC) FFLAGS="$(FFLAGS)"
	mv linpack/lib$(@).a lib
lapack:
	@echo "Compiling library: lapack"
	cd lapack && make pack FORTRAN=$(FC) OPTS="$(FFLAGS)"
	mv lapack/$(@).a lib/lib$(@).a
# eof special treatments

install:
	$(INSTALL) -d $(INSTALLDIR)
	$(INSTALL) $(MODSO) $(INSTALLDIR)
	$(INSTALLDATA) Multipack.py common_routines.py $(PYFILES) $(INSTALLDIR)
	cd $(INSTALLDIR) && echo "multipack" > ../multipack.pth

DOCMESS = "\
*************************************************\n\
Please find Multipack User's Guide in files:\n\
\tdoc/usersguide.dvi\n\
\tdoc/usersguide.ps\n\
\tdoc/usersguide.html (only if you have tth)\n\
*************************************************\n\
"
doc:
	cd doc && make
	@echo -e $(DOCMESS)
clean:
	rm *.o *.a *.so *.pyc

distclean:
	rm -f *.o *.so
	rm -f minpack/*.o minpack/*.a
	rm -f odepack/*.o odepack/*.a
	rm -f quadpack/*.o quadpack/*.a
	rm -f blas/*.o blas/*.a
	rm -f lapack/SRC/*.a lapack/BLAS/SRC/*.a
	rm -f lapack/SRC/*.o lapack/BLAS/SRC/*.o
	rm -f mach/*.o mach/*.a
	rm -f linpack/*.o linpack/*.a
	rm -f ddierckx/*.o ddierckx/*.a
	rm -f lib/*.a lib/*.o lib/*.c
	find . "(" -name '*~' -o -name '*.pyc' -o -name '*.pyo' \
	-o -name 'ref_*.pkl' ")" -exec rm {} \;
debug:
	@echo "MODULES=$(MODULES)"
	@echo "LINKLIBS=$(LINKLIBS)"
	@echo "LINKLIBSDIR=$(LINKLIBSDIR)"
	@echo "MLINKLIBS=$(MLINKLIBS)"
	@echo "LLINKLIBS=$(LLINKLIBS)"
	@echo "PYFILES=$(PYFILES)"
	@echo "PYFILENAMES=$(PYFILENAMES)"
	@echo "INSTALLDIR=$(INSTALLDIR)"
	@echo "INSTALL=$(INSTALL)"
