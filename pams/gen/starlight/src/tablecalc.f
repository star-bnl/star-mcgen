C     this subroutine calculates the tables that are used
c     elsewhere in the montecarlo
c     the tauon decay is taken from V-A theory, 1 - 1/3 cos(theta)
c     the energy of each of the two leptons in tau decay
c     is calculated using formula 10.35 given
c     in introduction to elementary particles by D. griffiths
c     which assmes that the mass of the electron is 0.
c     the highest energy attainable by the electron in such a system is
c     .5 * mass of the tau

      subroutine tablecalc

      implicit NONE 

      include 'const.inc'
      include 'tables.inc'
      integer i
      real E,theta

      tautolangle(0) = 0.
      dgammade(0) = 0.

      do 100 i= 1,100
c     calculate energy of tau decay
        E = float(i)/100. * .5 * mtau
        dgammade(i) = dgammade(i-1) + E**2. * (1. - 4.*E/(3.*mtau))

c     calculate angles for tau
        theta = pi * float(i) / 100.
        tautolangle(i) = tautolangle(i-1) + (1 + .3333333 * cos(theta))
 100  continue


      return
      end

