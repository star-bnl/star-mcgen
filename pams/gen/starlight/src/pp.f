c     this function calculates the perpendicular momentum distribution
c     function
c     gives a randomly chosen pp when given a photon energy
c     uses a binary seach assuming a monotonically increasing function

      real function pp(E)

      implicit NONE

      include 'const.inc'
      include 'D2LParam.inc'
      include 'inputp.inc'
      real E
      real ranx,ptest,dp,pdis,pperpdist
      real u0,u1,ran
      integer i

      u0 =2.*(E*RNuc/(hbarc*gamma_em))**2
      u1 =2.*(1.5)**2
c      call ranmar(ran,1)
      ranx = ran(iseed) * pperpdist(u0,u1)

c     ptest is a dimensionless pperp
      ptest = 0.75
      dp = ptest/2.

      do 50 i = 1,10
        u1 =2.*(ptest)**2
        pdis = pperpdist(u0,u1)
        if (ranx.lt.pdis) then
          ptest = ptest - dp
        else
          ptest = ptest + dp
        endif
        dp = dp /2.
 50   continue

      pp = ptest * hbarc / RNuc

      return
      end
