c     This subroutine calculates momentum and energy of vector meson
c     given W and Y.

      SUBROUTINE VMOMENTA(W,Y,E,px,py,pz,tcheck)

      IMPLICIT NONE

      include 'const.inc'
      include 'global.inc'
      include 'D2LParam.inc'
      include 'inputp.inc'
      include 'range.inc'
      REAL W,Y,dW,dY
      REAL ran,xtest,xt
      REAL Egam,Epom,tmin,pt1,pt2,phi1,phi2
      REAL px1,py1,px2,py2
      REAL E,pz,py,px,pt
      double precision formf,t1,t2
      INTEGER tcheck

      dW = (Wtop-Wmin)/DFLOAT(numw)

      dY  = (Ytop-Ymin)/DFLOAT(numy)      
      
C       >> Find Egam,Epom in CM frame
        Egam = 0.5*W*EXP(Y)
        Epom = 0.5*W*EXP(-Y)
C       >> Draw t1 according to form factor
 202    xt = ran(ISEED)
        t1  = 0.5*xt
C       >> Check tmin
        tmin = (Egam/gamma_em)**2
        IF(tmin.gt.0.5)THEN
          WRITE(*,*) 'WARNING: tmin=',tmin
          WRITE(*,*) 'Will pick a new W,Y'
	  tcheck = 1
	  return
        ENDIF
        IF( t1 .lt. tmin )GOTO 202
        xtest = ran(ISEED)
        IF( xtest .gt. formf(t1)*formf(t1) )GOTO 202
        pt1 = SQRT( t1 - tmin )
        phi1 = 6.28*ran(ISEED)

C       >> Draw t2 according to form factor
 203    xt = ran(ISEED)
        t2  = 0.5*xt
C       >> Check tmin
        tmin = (Epom/gamma_em)**2
        IF(tmin.gt.0.5)THEN
          WRITE(*,*) 'WARNING: tmin=',tmin
          WRITE(*,*) 'Will pick a new W,Y'
          tcheck = 1
	  return
        ENDIF
        IF( t2 .lt. tmin )GOTO 203
        xtest = ran(ISEED)
        IF( xtest .gt. formf(t2)*formf(t2) )GOTO 203
        pt2 = SQRT( t2 - tmin )
        phi2 = 6.28*ran(ISEED)

C       >> Compute px1,py1,px2,py2
        px1 = pt1*COS(phi1)
        py1 = pt1*SIN(phi1)
        px2 = pt2*COS(phi2)
        py2 = pt2*SIN(phi2)

        px = px1 + px2
        py = py1 + py2

C       Compute vector sum Pt = Pt1 + Pt2 to find pt for the vector meson
        pt = SQRT( px**2 + py**2 )
        
c	I guess W is the mass of the vector meson (not necessarily
c	on-mass-shell), and E is the energy
        E  = SQRT(W**2+pt**2)*COSH(Y)
        pz = SQRT(W**2+pt**2)*SINH(Y)

      return
      END
