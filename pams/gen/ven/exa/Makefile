###########################################################################
#
# Makefile for analysis of generator ZEBRA files
# Syntax: make FF=filename  , where filename is e.g. dndy_ana (w/o suffix!)
#
# Modified for STAR use:  Matt Bloomer 3/5/95
###########################################################################
#
# 	Package and STAR Library definitions
#
PKG_HOME:=$(dir $(shell pwd))
STAR_DIR=$(STAR_ROOT)/$(STAR_LEVEL)
#
#       Library definitions
#
LIBS =  `cernlib -v pro kernlib` \
	-legz \
        `cernlib -v pro packlib mathlib pawlib` \
        `cernlib -v pro geant321`
#
#	Library directories (mainly for libegz.a)
#
LIBDIRS = $(PKG_HOME)/newlib/$(STAR_ARCH) $(PKG_HOME)/lib/$(STAR_ARCH) \
	 $(STAR_LIB)/dev/sim/lib/
LLIBDIRS = $(addprefix -L,$(LIBDIRS)) 
#
#	Include files
#
INCLUDES = -I$(PKG_HOME)/inc -I$(STAR_DIR)/sim/inc
#
# 	Compile options:  use same as for lib/Makefile
#
include flags.mak
#
#	Executable directory
#
OUTDIR = $(PKG_HOME)/exa/$(STAR_ARCH)
#
#	Special treatment for AIX:  won't work for AIX 
#

.PHONY: usage clean test

usage:
	@ echo 
	@ echo Each target is the name of an executable.  Available targets are:
	@ echo
	@ echo " dndy_ana: rapidity distributions by species"
	@ echo " dump:	dump of egz data to screen
	@ echo " spectra: rapidity distributions by species (the right way)"
	@ echo
	@ echo "Make executable by typing make (target)"
	@ echo "N.B. not working for AIX"
	@ echo

% : %.F
	$(FC) $(INCLUDES) $(FFLAGS) $@.F read_zebra_evt.F \
	$(LLIBDIRS) $(LIBS) -o $(OUTDIR)/$@

clean:
	rm *.o
