C     carries out a lorentz transform of the frame.  (Not a
C     boost!)
      subroutine transform (betax,betay,betaz,E,px,py,pz,iFbadevent)

      implicit NONE

      double precision betax,betay,betaz,beta,gamma,gob
      real px,py,pz,E,E0,px0,py0,pz0
      integer iFbadevent

      E0 = E
      px0 = px
      py0 = py
      pz0 = pz

      beta = sqrt(betax**2 + betay**2 + betaz**2)
c	write(*,*) 'beta',beta
	if(beta.ge.1.0) iFbadevent = 1
      gamma = 1./sqrt(1. - beta**2)
      gob = (gamma-1)/beta**2
c	write(*,*) 'beta,gob',beta,gob
      E = gamma*(E0 - betax*px0 - betay*py0 - betaz*pz0)

      px = -gamma*betax*E0     + (1. + gob*betax**2)*px0
     &    +  gob*betax*betay*py0 + gob*betax*betaz*pz0

      py = -gamma*betay*E0          + gob*betay*betax*px0
     &    +  (1. + gob*betay**2)*py0 + gob*betay*betaz*pz0

      pz = -gamma*betaz*E0     +  gob*betaz*betax*px0
     &    +  gob*betaz*betay*py0 + (1. + gob*betaz**2)*pz0

      return
      end

