C     This subroutine calculates px py and pz given w and y 

      subroutine momenta(w,y,E,px,py,pz)

      implicit NONE

      include 'D2LParam.inc'
      include 'const.inc'
      include 'tables.inc'
      include 'inputp.inc'
      real y,px,py,pz,E1,E2,E,w,signpx,pp1,pp2,pt
      real anglepp1,anglepp2,pp,ran

c     E1 and E2 are for the two photons in the CM frame
      E1 = w * exp(y) / 2.
      E2 = w * exp(-y) / 2.
c	I think this way to calculate pz is wrong... it assumes pt=0
c	I'll do it in the way I think is right after getting pt
c      pz =  E1 - E2
      signpx = ran(iseed)
c      call ranmar(signpx,1)
      if (signpx.gt.0.5) pz = - pz


c     calculate px and py
c	pp gives the pt of each photon; multiply by cos(phi) or sin(phi)
c	to get x and y components-- phi is random between 0 and 2*pi
	anglepp1 = ran(iseed)
	anglepp2 = ran(iseed)
c      call ranmar(anglepp1,1)
c      call ranmar(anglepp2,1)
      pp1 = pp(E1)
      pp2 = pp(E2)
      px = pp1*cos(2.*pi*anglepp1)+pp2*cos(2.*pi*anglepp2)
      py = pp1*sin(2.*pi*anglepp1)+pp2*sin(2.*pi*anglepp2)

C       Compute vector sum Pt = Pt1 + Pt2 to find pt for the produced
c	particle
        pt = SQRT( px**2 + py**2 )

c       W is the mass of the produced particle (not necessarily
c       on-mass-shell); now compute its energy and pz
        E  = SQRT(W**2+pt**2)*COSH(Y)
        pz = SQRT(W**2+pt**2)*SINH(Y)

      return
      end

