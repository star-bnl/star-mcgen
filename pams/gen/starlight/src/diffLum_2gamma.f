c     This subroutine 
      subroutine diffLum_2gamma

      implicit NONE

      include 'inputp.inc'
      include 'D2LParam.inc'
      integer i,j
      real xlum,wmev,sum,D2LDMDY
      real w(1000),y(1000)

c     xlum is the luminosity as a function of w and y

c     Write the parameters being used for the calculation to 
c     starlight.dat so that they are there for future reference.
      open (unit=20,file='starlight.dat',status='unknown')
      write (20,*) Z
      write (20,*) A
      write (20,*) gamma_em
      write (20,*) wmax
      write (20,*) numw
      write (20,*) ymax
      write (20,*) numy
      write (20,*) gg_or_gP

c     Write to starlight.dat the set of values for w used in the
c     calculation.
      do 175 i = 1,numw
        w(i)  = wmax/float(numw)*float(i)
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
c	I don't understand why this is multiplied by wmev.. but at least
c	multiplying by wmev instead w(i) correctly accounts for the units
           xlum =  wmev * D2LDMDY(wmev,y(j))
           write(20,*) xlum
 300    continue
 400  continue

      close (unit=20)
      return
      end
