	subroutine setConst

	implicit NONE

	include 'const.inc'
	include 'global.inc'
	include 'D2LParam.inc'
	include 'inputp.inc'
	include 'bw.inc'

	double precision Wmin_default

c       define physical constants and masses
	hbarc = .197327053
	hbarcmev=hbarc*1000.
      	pi = 3.141592654
      	alpha = 1/137.0359895
      	mp = 0.93827231
      	mpi = 0.13956995
      	mK = 0.493677
      	mel = 0.00051099907
      	mmu = 0.105658389
      	mtau = 1.777

C  Generic radius first, then specific

        RNuc = 1.2 * A**(1.0/3.0)
	if (Z .eq. 79) RNuc=6.38
	if (Z .eq. 82) RNuc=6.62


c	unless otherwise defined later, default is Wmin=0
	Wmin_default = 0.

C  define constants for specific final states

c       define masses, widths and spins
        if(ip.eq.11) then 
	   mass = mel
c	Wmin is settable in  input file-- default is 0.01; Wmin up
c	to 0.15 is safe for Summer 2000 triggering for e+e- pairs
	   Wmin_default = 0.01
	   spin = 0.5
	endif
        if(ip.eq.13) then
	   mass = mmu
	   spin = 0.5
	   Wmin_default = 2.*mmu
	endif
        if(ip.eq.15) then
	   mass = 1.777
	   spin = 0.5
	   Wmin_default = 2*mtau
	endif
        if(ip.eq.115) then
	   mass = 1.318
           width =  1.04 * 10.0**(-6.)
	   spin = 2.
	endif
        if(ip.eq.221) then
	   mass = 0.54730
           width = 1. * 10.0**(-6.)
	   spin = 0.
	endif
        if(ip.eq.225) then
	   mass = 1.2754
           width = 2.6 * 10.0**(-6.)
	   spin = 2.
	endif
        if(ip.eq.331) then
	   mass = 0.95777
           width = 5. * 10.0**(-6.)
	   spin = 0.
	endif
        if(ip.eq.335) then
	   mass = 1.525
           width = .1 * 10.0**(-6.)
	   spin = 2.
	endif
        if(ip.eq.441) then
	   mass = 2.980
           width = 6.3 * 10.0**(-6.)
	   spin = 0.
	endif
        if(ip.eq.10221) then
	   mass = 0.980
           width = .56 * 10.0**(-6.)
	   spin = 0.
	endif
        if(ip.eq.33) then
	   mass = 1.540
           width = .1 * 10.0**(-6.)
	   spin = 2.
	endif
       	if (ip.eq.113) then
           mass = 0.7685
           width = 0.1507
           spin = 1.
           bslope=11.0
           f2o4pi=2.02
           Wmin_default = 2.*mpi
           Wmax = mass + 5.*width

C Breit-Wigner parameters (no direct pipi)

           ANORM =-2.75 
	   BNORM=0.
	   C=1.0

      	endif

C  rho0+direct pi+pi-.  maximum W set by Wmax

       	if (ip.eq.913) then
           mass = 0.7685
           width = 0.1507
           spin = 1.
           bslope=11.0
           f2o4pi=2.02
           Wmin_default = 2.*mpi
C  use the same (ZEUS) 1.5 GeV maximum mass
           Wmax = mass+5.*width

C  Breit-Wigner parameters, with direct pipi

	   ANORM =-2.75 
	   BNORM = 1.84
	   C =1.0
      	endif

      	if (ip.eq.223) then
           mass = 0.78194
           width = 0.00843
           spin = 1.
           bslope=10.0
           f2o4pi=23.13
           Wmin _default= mass - 5.*width
           Wmax = mass + 5.*width
C artificial Breit-Wigner parameters
	   ANORM=-2.75
	   BNORM=0.
	   C=1.0
      	endif
      	if (ip.eq.333) then
           width = 0.00443
           mass = 1.019413
           spin = 1.
           bslope=7.0
           f2o4pi=13.71
           Wmin_default = 2.*mK
           Wmax = mass + 5.*width
C artificial Breit-Wigner parameters
	   ANORM=-2.75
	   BNORM=0.
	   C=1.0
      	endif
       if (ip.eq.443) then
           mass = 3.09688
           width = 0.000087
           spin = 1.
           bslope=4.0
           f2o4pi=10.45
           Wmin_default = mass - 5.*width
           Wmax = mass + 5.*width
C artificial Breit-Wigner parameters - no direct pions
	   ANORM=-2.75
	   BNORM=0.
	   C=1.0
      	endif

c	set Wmin equal to the default values if Wmin is set 
c	to -1 in the input file
	if(Wmin.eq.-1) Wmin = Wmin_default

c	define luminosities, etc.
       	IF(Z.eq.79) THEN
           lum=2.
           Q0=0.060
           rho0=0.159407
      	ELSEIF(Z.eq.53)THEN
           lum=27.
           Q0=0.069
           rho0=0.161626
      	ELSEIF(Z.eq.29)THEN
           lum=95.
           Q0=0.087
           rho0=0.166878
      	ELSEIF(Z.eq.14)THEN
           lum=440.
           Q0=0.115
           rho0=0.177128
      	ELSEIF(Z.eq.8)THEN
           lum=980.
           Q0=0.138
           rho0=0.188459
      	ELSEIF(Z.eq.82)THEN
c       (for LHC)
           lum=1.
           Q0=0.059
           rho0=0.159176
      	ELSEIF(Z.eq.20)THEN
c       (for LHC)
           lum=20000.
           Q0=0.102
           rho0=0.171907
      	ELSE
           write (*,*) 'Warning: Luminosity  and density not'
	   write(*,*) 'defined for this projectile!'
           lum=1.
C  take 'typical' density
	   rho0=0.16
C  Q0 isn't used
	   Q0=0.
      	ENDIF

c	calculate Ep

        Ep= gamma_em*mp

C  rapidity range now comes from input file
C  n.b. Ymax is used for both 2-photon and gamma-p Ymin not used for 2 photon

           Ymin  =  -Ymax

C  find maximum photon energy (for VM only, in GeV, lab frame 
C  use beam energy, nuclear size based cutoff

	   EgMax=4.*gamma_em*hbarc/RNuc

	RETURN
	END



