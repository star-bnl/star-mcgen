      DOUBLE PRECISION FUNCTION nofe(Egamma,bimp)
C
C     Function for the calculation of the "photon density".
C     nofe = numberofgammas/(energy*area)     
C     Assume beta=1.0 and gamma>>1, i.e. neglect the
C     (1/gamma**2)*K0(X) term.
C
      IMPLICIT NONE

      DOUBLE PRECISION Egamma,bimp

      INCLUDE 'D2LParam.inc'
      INCLUDE 'const.inc'

      DOUBLE PRECISION DbesK1
      DOUBLE PRECISION X,factor1,factor2,factor3
C
      X=(bimp*Egamma)/(gamma_em*hbarc)
      factor1=(DFLOAT(Z*Z)*alpha)/(pi*pi)
      factor2=1./(Egamma*bimp*bimp)
      factor3=X*X*DbesK1(X)*DbesK1(X)
      nofe=factor1*factor2*factor3
C
      RETURN
      END
