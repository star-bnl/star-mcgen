C     This function calculates the cross-section as a function of
c     angle for a given W and Y, for the production of two muons.
c     (or tauons)
C     expression is taken from Brodsky et al. PRD 1971, 1532
c     equation 5.7
C     factor that are not dependant on theta are scrapped, so the
c     the absolute crosssections given by this function are inaccurate
c     here we are working in the CM frame of the photons and the last
c     term is 0

      real function thetalep (W,theta)

      implicit NONE

      include 'global.inc'
      real W,theta,W1sq
      real moverw

      W1sq =  (W / 2.)**2
      moverw = mass**2 / W1sq

      thetalep = 2. + 4.*(1-moverw)*
     &          ((1.-moverw)*sin(theta)**2*cos(theta)**2. + moverw) /
     &           (1. - (1. - moverw)*cos(theta)**2)**2


       return
       end

