	subroutine setConst

	implicit NONE

	include 'const.inc'
	include 'global.inc'
	include 'D2LParam.inc'
	include 'range.inc'
	include 'inputp.inc'
	include 'bw.inc'

c       define constants
	hbarc = .197327053
      	pi = 3.141592654
      	alpha = 1/137.0359895
      	mp = 0.93827231
      	mpi = 0.13956995
      	mK = 0.493677
      	mel = 0.00051099907
      	mmu = 0.105658389
      	mtau = 1.777
        RNuc = 1.2 * A**(1.0/3.0)

c       define masses, widths and spins
        if(ip.eq.11) then 
	   mass = 0.00051099907
	endif
        if(ip.eq.13) then
	   mass = 0.105658389
	   spin = 0.5
	endif
        if(ip.eq.15) then
	   mass = 1.777
	   spin = 0.5
	endif
        if(ip.eq.115) then
	   mass = 1.3181
           width =  1.04 * 10.0**(-6.)
	   spin = 2.
	endif
        if(ip.eq.221) then
	   mass = 0.54745
           width = 1. * 10.0**(-6.)
	   spin = 0.
	endif
        if(ip.eq.225) then
	   mass = 1.275
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
           Wmin = 2.*mpi
           Wtop = mass + 5.*width
      	endif
      	if (ip.eq.223) then
           mass = 0.78194
           width = 0.00843
           spin = 1.
           bslope=10.0
           f2o4pi=23.13
           Wmin = mass - 5.*width
           Wtop = mass + 5.*width
      	endif
      	if (ip.eq.333) then
           width = 0.00443
           mass = 1.019413
           spin = 1.
           bslope=7.0
           f2o4pi=13.71
           Wmin = 2.*mK
           Wtop = mass + 5.*width
      	endif
       if (ip.eq.443) then
           mass = 3.09688
           width = 0.087
           spin = 1.
           bslope=4.0
           f2o4pi=10.45
           Wmin = mass - 5.*width
           Wtop = mass + 5.*width
      	endif

c	define luminosities, etc.
       	IF(Z.eq.79) THEN
           lum=2.
           gamma_ta=108.4
           Q0=0.060
           rho0=0.159407
      	ELSEIF(Z.eq.53)THEN
           lum=27.
           gamma_ta=112.9
           Q0=0.069
           rho0=0.161626
      	ELSEIF(Z.eq.29)THEN
           lum=95.
           gamma_ta=124.5
           Q0=0.087
           rho0=0.166878
      	ELSEIF(Z.eq.14)THEN
           lum=440.
           gamma_ta=135.2
           Q0=0.115
           rho0=0.177128
      	ELSEIF(Z.eq.8)THEN
           lum=980.
           gamma_ta=135.2
           Q0=0.138
           rho0=0.188459
      	ELSEIF(Z.eq.82)THEN
c       (for LHC)
           lum=1.
           gamma_ta=2940.
           Q0=0.059
           rho0=0.159176
      	ELSEIF(Z.eq.20)THEN
c       (for LHC)
           lum=20000.
           gamma_ta=3730.
           Q0=0.102
           rho0=0.171907
      	ELSE
           write(*,*) 'ERROR: Luminosity not defined for this 
     &		projectile!'
           STOP 98
      	ENDIF

c	calculate Ep
        Ep= gamma_ta*mp

c	define rapidities, max photon energy
      	IF( Z.eq.79 .or. Z.eq.53 .or. Z.eq.29 .or. Z.eq.14 .or. 
     &		Z.eq.8 )THEN
C       >> At RHIC
           Ymin  =  -5.6
           Ytop  =   5.6
           EgMax =  25.0
      	ENDIF
      	IF( Z.eq.82 .or. Z.eq.20 )THEN
C       >> At LHC
           Ymin  =  -8.9
           Ytop  =   8.9
           EgMax = 600.0
      	ENDIF

c	define constants for Breit-Wigner normalization
	ANORM =-2.75 
	BNORM = 1.84
	BNORM_0 =0.0
	C =1.0
	RETURN
	END
