      SUBROUTINE HADRONBREAKUP(PBREAKUP,b,zp,ap,gamma)

C  This subroutine calculates the number of hadronic breakup
C  reactions at impact parameter b for a nucleus with Z=zp, A=ap
C  Woods-Saxon density profile is assumed

C  To normalize to an absolute  probablity P=1-exp(-pbreakup)

      IMPLICIT REAL*8(A-H,O-Z)

      integer IFIRST

      DIMENSION DEN2(20001),DEN1(20001)

C  this subroutine calculates the probability for hadronic nuclear breakup
C  at an impact parameter b for a nuclear with charge Z, atomic number A

      SAVE IFIRST
      DATA IFIRST /0/

      IF (IFIRST .NE. 0) GOTO 100

C  Initialize

C Integration delta x, delta z

      IFIRST=1
      DELL=.05
      DELR=.01

C  use two sigma_NN's.  52 mb at RHIC, 88 mb at LHC
C  gamma is in cm system

      SIGNN=5.2
      IF (gamma .GT. 500.) SIGNN=8.8


      if (zp .eq. 79) THEN
         R1=6.38
         A1=0.535
      elseif (zp .eq. 82) THEN
         R1=6.62
         A1=0.535
      ELSE
         WRITE(6,2)zp
 2       FORMAT(' NUCLEAR PARAMETERS NOT DEFINED FOR Z.', 
     *'USING DEFAULTS.  ZP =',F7.2)
         R1=1.2*A**(1./3.)
         A1=0.535
      ENDIF

      SIGNN=5.0

      WRITE(6,12)R1,A1,SIGNN
 12   FORMAT(' Nuclear density R= ',F7.4,' fm.  thick= ',F7.4,
     *' Sigma_NN= ',F7.4)
 
      R2=R1
      RHO1=ap
      RHO2=RHO1
      NZ1=((R1+5.)/DELR)
      NR1=NZ1
      NZ2=((R2+5.)/DELR)
      NR2=NZ2
      RR1=-DELR
      NY=((R1+5.)/DELL)
      NX=2*NY
      DO 47 IR1=1,NR1
      DEN1(IR1)=0.
      RR1=RR1+DELR
      Z1=-DELR/2.
      DO 45 IZ1=1,NZ1
      Z1=Z1+DELR
      RSQ=RR1*RR1+Z1*Z1
   45 DEN1(IR1)=DEN1(IR1)+1./(1.+EXP((SQRT(RSQ)-R1)/A1))
      DEN1(IR1)=DEN1(IR1)*2.*DELR
   47 DEN2(IR1)=DEN1(IR1)
      AN1=0.
      RR1=0.
      DO 61 IR1=1,NR1
      RR1=RR1+DELR
   61 AN1=AN1+RR1*DEN1(IR1)*DELR*2.*3.141592654
      AN2=AN1

      delo=.05
c .1 to turn mb into fm^2


C  Calculate breakup probability here

 100  PBREAKUP=0.
      if (b .gt. 25.) return

      Y=-.5*DELL
      DO 11 IY=1,NY
      Y=Y+DELL
      X=-DELL*FLOAT(NY+1)
      DO 11 IX=1,NX
      X=X+DELL
      XB=B-X
      RPU=SQRT(X*X+Y*Y)
      IRUP=(RPU/DELR)+1
      RTU=SQRT(XB*XB+Y*Y)
      IRUT=(RTU/DELR)+1
      T1=DEN2(IRUT)*RHO2/AN2
      T2=DEN1(IRUP)*RHO1/AN1
c  Eq.(6) BCW, Baltz, Chasman, White, Nucl. Inst. & Methods A 417, 1 (1998):
 11   PBREAKUP=PBREAKUP+2.*T1*(1.-EXP(-SIGNN*T2))*DELL*DELL
      RETURN
      END
