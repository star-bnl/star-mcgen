c     File:- event.inc
c     =================
c
c     This file is intended to be used by the analysis programs.
c     The information from the ZEBRA banks is being unpacked into
c     the common block /event/ with the primary information as
c     read by the input and with secondary information as requested
c     by setup parameters.
c
c	maxptl  : maximum number of particles pr event.
c	nptls	: number of particles in event
c	nnproj  : # neutron spectators (in proj)
c	npproj  : # proton  spectators (in proj)
c	naproj  : # spectators (in proj)
c	nntarg  : # neutron spectators (in target)
c	nptarg  : # proton  spectators (in target)
c	natarg  : # spectators (in target)
c       ntry    : # of tries to get an interaction since last call
c                 This number is needed to get a cross section from
c		        the number of attempts + the impact parameter
c		        range.
c
c	px	: 4 momentum of each particle 
c	py	: 
c	pz	:
c	energy	:
c	idptl	: ID for particle (as given by event generator)
c	ioptl 	: Parent sequence number ( 0 if no parent)
c
c	pt	: perpendicular momentum (sqrt(px**2+py**2)
c	y	: rapidity
c	pt      : transverse momentum (GeV/c)
c	mass	: mass in GeV of particle
c	theta	: polar angle theta in radians
c	phi	: azimuthal angle in radians
c	gtype	: GEANT particle ID.
c
      INTEGER MAXPTL
      PARAMETER (MAXPTL=30000)
      INTEGER NPTLS,NTRY,NNPROJ,NPPROJ,NAPROJ,NNTARG,NPTARG,NATARG
     +,IDPTL(MAXPTL),IOPTL(MAXPTL),GTYPE(MAXPTL)
      REAL BIMEVT,PX(MAXPTL),PY(MAXPTL),PZ(MAXPTL),ENERGY(MAXPTL)
     +,MASS(MAXPTL),THETA(MAXPTL),PHI(MAXPTL),Y(MAXPTL),PT(MAXPTL)
c
      COMMON/EGEVNT/NPTLS,BIMEVT,NNTARG,NPTARG,NATARG,NNPROJ,NPPROJ
     +,NAPROJ,NTRY,PX,PY,PZ,ENERGY,IDPTL,IOPTL,MASS,THETA,PHI,Y,GTYPE
     +,PT
c
c	Reaction description
c
c	atarg  : Target Mass Number
c	ztarg  : Target Charge Number
c	aproj  : Projectile Mass number
c	atarg  : Projectile Charge Number
c	Bmin   : minimum impact parameter (normally 0)
c	Bmax   : maximum impact parameter
c	Nbeam  : Number of beam particles thrown into the
c	          area defined by bmin->bmax
c	
       INTEGER ATARG,ZTARG,APROJ,ZPROJ,NBEAM
       REAL SQRTS,BMIN,BMAX
       COMMON/EGRXN/ATARG,ZTARG,APROJ,ZPROJ,SQRTS,BMIN,BMAX,NBEAM
c
