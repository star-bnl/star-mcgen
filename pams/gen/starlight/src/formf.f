c     This function calculates the nuclear form factor assuming a hard
c     sphere + Yukawa.

      DOUBLE PRECISION FUNCTION formf(t)

      IMPLICIT NONE

      include 'D2LParam.inc'
      include 'const.inc'
      DOUBLE PRECISION t
      Double Precision q,sph,yuk,a0,R,r0

C     >> Use Parameterization from FRITIOF
      r0=1.16*(1.-1.16*A**(-2./3.))
      R=r0*A**(1./3.)

      q   = DSQRT(t)

      sph = DSIN(q*R/hbarc) - (q*R/hbarc)*DCOS(q*R/hbarc)
      sph = ((3.*hbarc*hbarc*hbarc)/(q*q*q*r0*r0*r0))*sph
      sph = sph/DFLOAT(A)

      a0  = 0.70
      yuk = 1./( 1. + ((a0*a0*q*q)/(hbarc*hbarc)) )
      formf = sph*yuk

      RETURN
      END
