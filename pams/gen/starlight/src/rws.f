C     Wood-Saxon Nuclear Density Distribution for Nucleus
C     with mass number A.

      DOUBLE PRECISION FUNCTION rws(r)

      IMPLICIT NONE

      include 'D2LParam.inc'
      DOUBLE PRECISION r
      DOUBLE PRECISION r0,RadNuc,C

C     >> Approximately r0=1.12 (R=r0*A**1/3), C=0.53
C     >> Use Parameterization from FRITIOF
      r0=1.16*(1.-1.16*A**(-2./3.))
      RadNuc=r0*A**(1./3.)
C     >> Use constant skin-thickness, C=0.53 fm
      C=0.53

      rws=1.0/( 1. + DEXP( (r-RadNuc)/C ) )

      RETURN
      END
