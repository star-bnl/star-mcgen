c     This function calculates cross section for two-particle decay.
c     For reference, see STAR Note 243, Eq. 9.

      subroutine sigma2

      implicit NONE

      include 'inputp.inc'
      include 'Ftable.inc'
      include 'sig.inc'
      include 'const.inc'
      integer i,j
      real sum,sigmui

c     calculate the two-lepton differential cross section
c     the 100 is to convert to barns
c     the two is to account for the fact that we integrate only
c     over one half of the rapidity range
      do 200 i  = 1,numw
      do 100 j = 1, numy
        sigma (i,j) = 2. * sigmui(warray(i))* farray(i,j) / 100.
c	write(*,*) 'i,j,sigma',i,j,sigma(i,j)
c	write(*,*) 'w,f',warray(i),farray(i,j)
 100  continue
 200  continue

c     calculate the total two-lepton cross section
      sum = 0.
      do 400 i  = 1,numw-1
      do 300 j = 1, numy-1
        sum = sum + 2.* (   ( sigma(i,j)+sigma(i+1,j)+
     &                 sigma(i,j+1)+sigma(i+1,j+1) )/4. *
     &                (yarray(j+1)-yarray(j))   *
     &                (warray(i+1)-warray(i))    /
     &                ((warray(i+1)+warray(i))/2.)  )

 300  continue
 400  continue


      if (ip.eq.11)
     &     print *,'The total electron crossection is:',sum

      if (ip.eq.13)
     &     print *,'The total muon crossection is:',sum

      if (ip.eq.15)
     &     print *,'The total tau crossection is:',sum

      return
      end
