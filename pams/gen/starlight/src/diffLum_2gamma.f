c     This subroutine computes Coulomb/Hadronic breakup probabilities, and then does diff. lum. calculation
      subroutine diffLum_2gamma

      implicit NONE

      include 'inputp.inc'
      include 'const.inc'
      include 'D2LParam.inc'
      Double Precision OldNorm
      integer i,j
      real xlum,wmev,sum,D2LDMDY
      real w(1000),y(1000)

      Normalize = 1./SQRT(1.*numw*numy) !if your grid is very fine, you'll want high accuracy --> small Normalize
      OldNorm   = Normalize
c     xlum is the luminosity as a function of w and y

c     Write to starlight.dat the set of values for w used in the
c     calculation.
      do 175 i = 1,numw
        w(i)  = Wmin + (wmax-Wmin)/float(numw)*float(i)
        write(20,*) w(i)
 175  continue

c     Write to starlight.dat the set of values for y used in the
c     calculation.
      do 200 i = 1,numy
        y(i) = ymax*float(i-1)/float(numy)
        write(20,*) y(i)
 200  continue

c     For each (w,y) pair, calculate the differential luminosity and write
c     it to starlight.dat.
      do 400 i = 1,numw
        sum = 0.
        do 300 j = 1,numy
          wmev = w(i) *1000.

C     Convert from photon flux dN/dW to Lorentz invariant
C     photon number WdN/dW 

           xlum =  wmev * D2LDMDY(wmev,y(j))
           write(20,*) xlum 
           if (j.eq.1) OldNorm = Normalize ! save value of Integral for each new W(i) and Y(1)
 300    continue
        Normalize = OldNorm
 400  continue

      return
      end
