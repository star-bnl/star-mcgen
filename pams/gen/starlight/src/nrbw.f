      DOUBLE PRECISION FUNCTION nrbw(W,ANORM,B,C)

C     Relativistic Breit-Wigner according to J.D. Jackson,
C     Nuovo Cimento 34, 6692 (1964), with nonresonant term. A is the strength
C     of the resnonant term and b the strength of the non-resonant
C     term. C is an overall normalization. 

      IMPLICIT NONE

      include 'global.inc'
      include 'const.inc'
      include 'inputp.inc'
      DOUBLE PRECISION W,ANORM,B,C
      DOUBLE PRECISION ppi,ppi0,GammaPrim,rat
      DOUBLE PRECISION aa,bb,cc

C  width depends on energy - Jackson Eq. A.2

C  if below threshold, then return 0.  Added 5/3/2001 SRK
C  0.5% extra added for safety margin
      
      if (W .LT. 2.01*mpi) THEN
         nrbw=0.
         return
      ENDIF
      ppi=DSQRT( (W/2.)**2 - mpi*mpi )
      ppi0=0.358


C handle phi-->K+K- properly
      if (ip .eq. 333) then
         if (W .LT. 2.*mK) THEN
            nrbw=0.
            return
         ENDIF
         ppi=DSQRT( (W/2.)**2- mK*mK)
         ppi0=DSQRT( (mass/2)**2-mK*mK)
      endif

      rat=ppi/ppi0
      GammaPrim=width*(mass/W)*rat*rat*rat

      aa=ANORM*DSQRT(GammaPrim*mass*W)
      bb=W*W-mass*mass
      cc=mass*GammaPrim

C  real part^2

       nrbw = ( (aa*bb)/(bb*bb+cc*cc) + B)**2

C imaginary part^2

       nrbw = nrbw + ( (aa*cc)/(bb*bb+cc*cc) )**2


C  Alternative, a simple, no-background BW, following J. Breitweg et al.
C  Eq. 15 of Eur. Phys. J. C2, 247 (1998).  SRK 11/10/2000
C      nrbw = (ANORM*mass*GammaPrim/(bb*bb+cc*cc))**2

      nrbw = C*nrbw

      RETURN
      END


