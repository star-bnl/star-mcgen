c	2/3/2000	JS
c	commented out, for the short term, the lines which make it not
c	write out if the values of Egamma go outside the prescribed range

c	This code calculates a table of differential luminosity as a
c	function of W (CM energy) and Y.  dN/dWdY is calculated as 
c	dN/dWdY = (photon energy)*flux*( cross section)*(Breit-Wigner width).  
c	For reference, see STAR Note 386.

      SUBROUTINE diffLum_vm

      IMPLICIT NONE

      include 'inputp.inc'
      include 'D2LParam.inc'
      include 'const.inc'
      include 'global.inc'
      include 'range.inc'
      include 'bw.inc'
      DOUBLE PRECISION formf,flux,sigmagp,nrbw,sigma_A
      DOUBLE PRECISION Av,Wgp,cs,cvma
      DOUBLE PRECISION W,dW,dY
      DOUBLE PRECISION Egamma,y
      DOUBLE PRECISION t,tmin,tmax
      DOUBLE PRECISION csgA,int
      DOUBLE PRECISION bwnorm,testint,dndWdY
      DOUBLE PRECISION xg(1:5),ag(1:5)
      DOUBLE PRECISION ax,bx
      INTEGER          I,J,K,NGAUSS, Z, A


C     DATA FOR GAUSS INTEGRATION
      DATA xg/0.1488743390,0.4333953941,0.6794095683,0.8650633667,
     $        0.9739065285/
      DATA ag/0.2955242247,0.2692667193,0.2190863625,0.1494513492,
     $        0.0666713443/
      NGAUSS = 5

c     Write parameters for this calculation to starlight.dat so they
c     are there for future reference.
      open (unit=20,file='starlight.dat',status='unknown')
      write (20,*) Z
      write (20,*) A
      write (20,*) gamma_em
      write (20,*) Wtop
      write (20,*) numw
      write (20,*) Ytop
      write (20,*) numy
      write (20,*) gg_or_gP

      WRITE(*,*) 'Generating W/Y map for Monte Carlo ...'

      dW = (Wtop-Wmin)/DFLOAT(numw)
      dY  = (Ytop-Ymin)/DFLOAT(numy)      

C     Normalize the Breit-Wigner Distribution
      testint=0.0

c     Write the values of W used in the calculation to starlight.dat.
      DO 50 I=0,numw-1
        W = Wmin + DFLOAT(I)*dW + 0.5*dW
        testint = testint + nrbw(W,ANORM,BNORM_0,C)*dW
	write(20,*) W 
 50   CONTINUE

c     Write the values of Y used in the calculation to starlight.dat.
      DO 60 I=0,numy-1
        Y = Ymin + DFLOAT(I)*dY + 0.5*dY
	write(20,*) Y
 60   CONTINUE

      bwnorm = 1./testint
      WRITE(*,*) 'BW Norm:',bwnorm

      Eth=0.5*(((Wmin+mp)*(Wmin+mp)-mp*mp)/(Ep+DSQRT(Ep*Ep-mp*mp)))

      int=0.
      DO 102 I=0,numw-1

        W = Wmin + DFLOAT(I)*dW + 0.5*dW

        DO 101 J=0,numy-1
          Y = Ymin + DFLOAT(J)*dY + 0.5*dY
          Egamma = 0.5*W*DEXP(Y)

          IF(Egamma.lt.Eth) then
c	GOTO 101
c	this is just a kluge to keep it from not writing out when
c	outside range-- have to think about the best way to handle this
	Egamma = Eth
	endif
          IF(Egamma.gt.EgMax) then
c	GOTO 101
c	this is just a kluge to keep it from not writing out when
c	outside range-- have to think about the best way to handle this
	Egamma = EgMax
	endif

C       Find gamma-proton CM energy
          Wgp=DSQRT(2.*Egamma*(Ep+DSQRT(Ep*Ep-mp*mp))+mp*mp)

C       Calculate V.M.+proton cross section
	  cs=DSQRT(16.*pi*f2o4pi*bslope*hbarc*hbarc*
     &     sigmagp(Wgp)/alpha)

C       Calculate V.M.+Nucleus cross section
          cvma=sigma_A(cs)

C       Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
          Av=(alpha*cvma*cvma)/(16.*pi*f2o4pi*hbarc*hbarc)

          tmin   = ( (W**2)/(4.*Egamma*gamma_ta) )**2
          tmax   = tmin + 0.25
          ax     = 0.5*(tmax-tmin)
          bx     = 0.5*(tmax+tmin)
          csgA   = 0.
          DO 100 K=1,NGAUSS

            t    = ax*xg(K)+bx
            csgA = csgA + ag(K)*formf(t)*formf(t)
            t    = ax*(-xg(K))+bx
            csgA = csgA + ag(K)*formf(t)*formf(t)

 100      CONTINUE
          csgA = 0.5*(tmax-tmin)*csgA
          csgA = Av*csgA

          dndWdY = Egamma*flux(Egamma)*csgA*nrbw(W,ANORM,BNORM_0,
     &	bwnorm)
	  write(20,*) dndWdY

 101    CONTINUE
 102  CONTINUE

      close (unit=20)

      RETURN
      END
