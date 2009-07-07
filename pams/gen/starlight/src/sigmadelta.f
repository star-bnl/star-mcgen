c     	This subroutine carries out a delta function cross section. 
c	calculation.  For reference, see STAR Note 243, Eq. 8.

      subroutine sigmadelta

      implicit NONE

      include 'global.inc'
      include 'inputp.inc'
      include 'const.inc'
      include 'Ftable.inc'
      include 'sig.inc'
      integer j,ivalw,i
      real remainw,sum

      wdelta = mass

CVM: multipy all farray() by f_max

c     calculate the differential cross section and place in the sigma table
      do 300 i  = 1,numw
      do 400 j = 1, numy


C  Changed '8' to '4' in the following equation, to fix the integration of
C  a delta function. Now, it matches the standard literature (cf. Eq. 67 of
C  G. Baur et al, Phys. Rep. 364, 359 (2002).  STAR Note 243 gives the 
C  incorrect '4'

        sigma (i,j) = (spin*2.+1.)*
     &   4*pi**2.*width/(mass**3.)*f_max*farray(i,j)*hbarc**2./100.

 400  continue
 300  continue

c     find the index, i, for the value of w just less than the mass
c     because we want to use the value from the sigma table that has w =
c     mass
      do 100 i = 1,numw
        if(mass.gt.warray(i)) ivalw = i
 100  continue

      remainw = (mass-warray(ivalw)) / (warray(ivalw+1)-warray(ivalw))
      ivalwd = ivalw
      remainwd = remainw

c     if we are interested rho pairs at threshold
c     then just set sigma to be 100 nb
      if (ip.eq.33) then
        sum  = 0
        do 500 j = 1,numy-1
          sum = sum + 2.0 * (yarray(j+1) - yarray(j)) *
     &            100.0 *  10.0**(-9.) * (.1/mass) *
     & ( (1.-remainw)*f_max*(farray(ivalw  ,j)+farray(ivalw  ,j+1))/2. +
     &   (   remainw)*f_max*(farray(ivalw+1,j)+farray(ivalw+1,j+1))/2. )
 500    continue
      else

c       sum to find the total cross section
c       the two is to account for the fact that we integrate
c       over just 1/2 of the rapidity range
        sum = 0.
        do 600 j = 1,numy-1
          sum =  sum + 2. *  (yarray(j+1) - yarray(j)) *
     &    ( (1.-remainw)*(sigma(ivalw  ,j)+sigma(ivalw  ,j+1))/2. +
     &      (   remainw)*(sigma(ivalw+1,j)+sigma(ivalw+1,j+1))/2. )
 600    continue
      endif

      print *,'the total crossection is:',sum, ' barns.'

      return
      end
