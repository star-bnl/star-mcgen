c     This subroutine calculates the vector meson cross section assuming
c     a narrow resonance.  For reference, see STAR Note 386.

      SUBROUTINE sigmavm

      IMPLICIT NONE

      include 'range.inc'
      include 'global.inc'
      include 'const.inc'
      DOUBLE PRECISION formf,flux,sigmagp,sigma_A
      DOUBLE PRECISION Av,Wgp,cs,cvma
      DOUBLE PRECISION W,Ymin,Ytop,dY
      DOUBLE PRECISION y1,y2,y12,ega1,ega2,ega12
      DOUBLE PRECISION t,tmin,tmax
      DOUBLE PRECISION csgA1,csgA2,csgA12,int,dR,rate
      DOUBLE PRECISION tmp
      DOUBLE PRECISION xg(1:5),ag(1:5)
      DOUBLE PRECISION ax,bx
      INTEGER          J,K,NY,NGAUSS

C     >> DATA FOR GAUSS INTEGRATION
      DATA xg/0.1488743390,0.4333953941,0.6794095683,0.8650633667,
     $        0.9739065285/
      DATA ag/0.2955242247,0.2692667193,0.2190863625,0.1494513492,
     $        0.0666713443/
      NGAUSS = 5

      NY   =  1200
      dY   = (Ytop-Ymin)/DFLOAT(NY)


      WRITE(*,*) 'Using Narrow Resonance ...'

      W = mass
      Eth=0.5*(((W+mp)*(W+mp)-mp*mp)/(Ep+DSQRT(Ep*Ep-mp*mp)))
      WRITE(*,*) 'gamma+nucleon  Threshold:',Eth
      int=0.

      tmp = 0.0
      DO 200 J=0,NY-1

        y1  = Ymin + DFLOAT(J)*dY
        y2  = Ymin + DFLOAT(J+1)*dY
        y12 = 0.5*(y1+y2) 

        ega1  = 0.5*W*DEXP(y1)
        ega2  = 0.5*W*DEXP(y2)  
        ega12 = 0.5*W*DEXP(y12)

        IF(ega1.lt.Eth)GOTO 200
        IF(ega2.gt.EgMax)GOTO 200

        IF(J.eq.0)THEN
C       >> 1st Point (Calculated only the first time)     =====>>>

C       >> Find gamma-proton CM energy
        Wgp=DSQRT(2.*ega1*(Ep+DSQRT(Ep*Ep-mp*mp))+mp*mp)
C       >> Calculate V.M.+proton cross section
        cs=DSQRT(16.*pi*f2o4pi*bslope*hbarc*hbarc*sigmagp(Wgp)/alpha)
C       >> Calculate V.M.+Nucleus cross section
        cvma=sigma_A(cs)
C       >> Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
        Av=(alpha*cvma*cvma)/(16.*pi*f2o4pi*hbarc*hbarc)

        tmin  = ( (W**2)/(4.*ega1*gamma_ta) )**2
        tmax  = tmin + 0.25
        ax    = 0.5*(tmax-tmin)
        bx    = 0.5*(tmax+tmin)
        csgA1  = 0.
        DO 50 K=1,NGAUSS

          t     = ax*xg(K)+bx
          csgA1 = csgA1 + ag(K)*formf(t)*formf(t)
          t     = ax*(-xg(K))+bx
          csgA1 = csgA1 + ag(K)*formf(t)*formf(t)

 50    CONTINUE
        csgA1 = Av*0.5*(tmax-tmin)*csgA1
        ELSE
        csgA1 = csgA2
        ENDIF

C       >> Middle Point                      =====>>>

C       >> Find gamma-proton CM energy
        Wgp=DSQRT(2.*ega12*(Ep+DSQRT(Ep*Ep-mp*mp))+mp*mp)
C       >> Calculate V.M.+proton cross section
        cs=DSQRT(16.*pi*f2o4pi*bslope*hbarc*hbarc*sigmagp(Wgp)/alpha)
C       >> Calculate V.M.+Nucleus cross section
        cvma=sigma_A(cs)
C       >> Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
        Av=(alpha*cvma*cvma)/(16.*pi*f2o4pi*hbarc*hbarc)

        tmin   = ( (W**2)/(4.*ega12*gamma_ta) )**2
        tmax   = tmin + 0.25
        ax     = 0.5*(tmax-tmin)
        bx     = 0.5*(tmax+tmin)
        csgA12 = 0.
        DO 60 K=1,NGAUSS

          t      = ax*xg(K)+bx
          csgA12 = csgA12 + ag(K)*formf(t)*formf(t)
          t      = ax*(-xg(K))+bx
          csgA12 = csgA12 + ag(K)*formf(t)*formf(t)

 60    CONTINUE
        csgA12 = 0.5*(tmax-tmin)*csgA12
        csgA12 = Av*csgA12


C       >> Second Point                      =====>>>

C       >> Find gamma-proton CM energy
        Wgp=DSQRT(2.*ega2*(Ep+DSQRT(Ep*Ep-mp*mp))+mp*mp)
C       >> Calculate V.M.+proton cross section
        cs=DSQRT(16.*pi*f2o4pi*bslope*hbarc*hbarc*sigmagp(Wgp)/alpha)
C       >> Calculate V.M.+Nucleus cross section
        cvma=sigma_A(cs)
C       >> Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
        Av=(alpha*cvma*cvma)/(16.*pi*f2o4pi*hbarc*hbarc)

        tmin  = ( (W**2)/(4.*ega2*gamma_ta) )**2
        tmax  = tmin + 0.25
        ax    = 0.5*(tmax-tmin)
        bx    = 0.5*(tmax+tmin)
        csgA2 = 0.
        DO 70 K=1,NGAUSS

          t     = ax*xg(K)+bx
          csgA2 = csgA2 + ag(K)*formf(t)*formf(t)
          t     = ax*(-xg(K))+bx
          csgA2 = csgA2 + ag(K)*formf(t)*formf(t)

 70    CONTINUE
        csgA2 = 0.5*(tmax-tmin)*csgA2
        csgA2 = Av*csgA2


C       >> Sum the contribution for this W,Y. 
        dR  = ega1*flux(ega1)*csgA1
        dR  = dR + 4.*ega12*flux(ega12)*csgA12
        dR  = dR + ega2*flux(ega2)*csgA2
        tmp = tmp+2.*dR*(dY/6.)
        dR  = dR*(dY/6.)
C       >> The 2 accounts for the 2 beams
        dR  = 2.*dR
        int = int+dR

 200    CONTINUE

      rate=lum*int
      WRITE(*,*) 'Cross section (mb):',10.*int
      WRITE(*,*) 'Production rate   :',rate,' Hz'

      RETURN
      END
