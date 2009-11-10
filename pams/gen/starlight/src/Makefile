#########################################################################
#
# Makefile for gammagamma monte carlo
#
#
########################################################################
#
#choose ONE of the following two lines, depending on operating system
include Makefile.Linux
#include Makefile.SunOS
#
BIN= ../bin
#

#SLIB = `cernlib mathlib packlib kernlib jetset74`
SLIB = `cernlib mathlib jetset74`

#SLIB = `cernlib -v pro mathlib packlib kernlib jetset74`
#
#CGLIB=`cernlib -v new -G mathlib` 
#
# 	Libraries
#
#
#CERNLIB=/cern/pro/lib
#LIBS= -L$(CERNLIB) -lmathlib -ljetset74 -lkernlib
#
# Additional suffixes
#
.SUFFIXES: .F .cdf  .car .cra
#
#      Executable Name
EXE = starlight
#
#
#          Fortran Sources
#
FSOURCES = starlight.f \
	   vladimir.f \
	   diffLum_vm.f \
	   flux.f \
	   formf.f \
	   nrbw.f \
	   sigma_A.f \
	   sigmagp.f \
	   vmomenta.f\
	   sigmavm.f\
	   rws.f\
	   t.f\
	   sigmavmw.f\
	   setConst.f\
	   momenta.f\
	   decayEvent.f\
	   taudecay.f\
	   thetalep.f\
	   transform.f\
	   thephi.f\
	   input.f\
           newparam.f\
	   diffLum_2gamma.f\
	   readDiffLum.f\
	   sigmacalc.f\
	   sigmadelta.f\
	   sigma2.f\
	   sigmui.f\
	   pickw.f\
	   picky.f\
	   tablecalc.f\
	   pp.f\
	   pperpdist.f\
	   writejetsetText.f\
	   writejetsetGSTARtext.f\
	   writejetsetNtuple.f\
	   writeText.f\
	   writeGSTARtext.f\
	   writeNtuple.f\
	   jtog.f\
	   twodecay.f\
	   pickwy_vm.f\
	   heptup.f\
	   photonbreakup.f\
	   hadronbreakup.f\
	   pttablegen.f\
	   vmpt.f\
	   vmsigmapt.f\
	   nofe.f\
           radmul.F\
	   PofB.f



#
#      C Sources
#
CSOURCES =
#
#
#
# Objects list
#
OBJECTS= $(FSOURCES:.f=.o) \
         $(CSOURCES:.c=.o)
#
#
#
#---------
main:
	make $(BIN)/$(EXE)
#
# Target definitions
#
#
$(BIN)/$(EXE):  $(OBJECTS)
# use the following line to compile on Solaris
#	$(FC) -o $(BIN)/$(EXE) $(OBJECTS) $(SLIB) -L/opt/WS5.0/lib -lM77 -lM77 -lF77 -lsunmath
#use the following line to compile on Linux
	$(FC) -o $(BIN)/$(EXE) $(OBJECTS) $(SLIB) -L/usr/pgi/linux86/lib -lpgftnrtl -lpgc -L/usr/local/pkg/egcs-1.1.2-lhc/lib/gcc-lib/i686-pc-linux-gnu/egcs-2.91.66 -lg2c 
#
#
#    Second-level dependencies
#
.F.o:    
	$(FC) -c -o $*.o $*.F                              

#  Cleanup

	rm -rf $(BIN)/starlight.dat










