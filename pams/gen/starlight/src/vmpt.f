c     This subroutine calculates momentum and energy of vector meson
c     given W and Y, including interference. 
C     It gets the pt distribution from a lookup table.  

      SUBROUTINE VMPT(W,Y,E,px,py,pz,tcheck)

      IMPLICIT NONE

      include 'const.inc'
      include 'global.inc'
      include 'D2LParam.inc'
      include 'inputp.inc'

      include 'Ftable.inc'
      include 'pttable.inc'

      REAL W,Y,dW,dY
      REAL ran,theta,xpt,ypt

      REAL yfract, ptfract,pt1,pt2,yleft

      REAL E,pz,py,px,pt


      INTEGER tcheck, IY,IPT,j

      dW = (Wmax-Wmin)/DFLOAT(numw)

      dY  = (Ymax-Ymin)/DFLOAT(numy)      
      
C  Y is already fixed; choose a pt
C  Follow the approavh in pickwy.f

C in   fptarray(IY,pt) IY=1 corresponds to Y=0, IY=numy/2 corresponds to +y

      IY=INT(ABS(Y)/dY)+1
      if (IY .GT. numy/2) THEN
C         WRITE(6,22)IY,numy
C 22      FORMAT(' Error in VMPT IY>numy/2.  IY, numy=',I6,2X,I6)
         IY=numy/2
      ENDIF

C      yfract=ABS(Y)-(2.*IY-numy)*dY
       yleft=ABS(Y)-(IY-1)*dY
       yfract=yleft*dY


C      if (fract .lt. 0. .or. fract .gt. 1) write(6,23)fract,Y,IY
C 23   format(' error in vmpt.  fract=',E11.3,'  Y,IY=',F7.4,2X,I5)

C       write(6,111)Y,IY,ptmax,dpt
C  111  format(' VMPT. Y,IY,PTMAX,dpt=',F7.4,2X,I5,2X,F7.4,2X,F7.4)

 50   xpt=ran(iseed)

      do 55 j=1,NPT+1
         if (xpt .lt. fptarray(IY,j)) goto 60 
 55   continue
 60   continue
      
C  now do linear interpolation - start with extremes

      if (j .eq. 1) then
         pt1=xpt/fptarray(IY,J)*dpt/2.
         goto 80
      endif
      if (j .eq. NPT+1) then
         pt1=(ptmax-dpt/2.) + 
     *   dpt/2.*(xpt-fptarray(IY,J))/(1.-fptarray(IY,J))
         goto 80
      endif

C  we're in the middle

      ptfract=(xpt-fptarray(IY,J))/(fptarray(IY,J+1)-fptarray(IY,J))
      pt1=j*dpt+ptfract*dpt

C  at an extreme in y?
      if (IY.eq. numy/2) then
         pt=pt1
         goto 120
      endif

 80   continue

C  interpolate in y repeat for next fractional y bin

      do 85 j=1,NPT+1
         if (xpt .lt. fptarray(IY+1,j)) goto 90 
 85   continue
 90   continue
      
C  now do linear interpolation - start with extremes

      if (j .eq. 1) then
         pt2=xpt/fptarray(IY+1,J)*dpt/2.
         goto 100
      endif
      if (j .eq. NPT+1) then
         pt2=(ptmax-dpt/2.) + 
     *   dpt/2.*(xpt-fptarray(IY+1,J))/(1.-fptarray(IY+1,J))
         goto 100
      endif

C  we're in the middle

      ptfract=(xpt-fptarray(IY+1,J))
     */(fptarray(IY+1,J+1)-fptarray(IY+1,J))
      pt2=j*dpt+ptfract*dpt

 100  continue

C  now interpolate in y

      pt=yfract*pt2+(1-yfract)*pt1

 120  continue

C  we have a pt
 
      theta=ran(iseed)*2.*pi
      px=pt*cos(theta)
      py=pt*sin(theta)

C     write(6,*)xpt,theta,px,py

c	I guess W is the mass of the vector meson (not necessarily
c	on-mass-shell), and E is the energy

        E  = SQRT(W**2+pt**2)*COSH(Y)
        pz = SQRT(W**2+pt**2)*SINH(Y)
c	randomly choose to make pz negative 50% of the time
	if (ran(iseed).ge.0.5) pz = -pz
      return
      END



