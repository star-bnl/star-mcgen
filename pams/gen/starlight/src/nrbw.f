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
      DOUBLE PRECISION aa,bb,cc,dd



C  width depends on energy - Jackson Eq. A.2

      ppi=DSQRT( (W/2.)**2 - mpi*mpi )
      ppi0=0.358


C handle phi-->K+K- properly
      if (ip .eq. 333) then
         ppi=DSQRT( (W/2.)**2- mK*mK)
         ppi0=DSQRT( (mass/2)**2-mK*mK)
      endif

      rat=ppi/ppi0
      GammaPrim=width*(mass/W)*rat*rat*rat

      aa=ANORM*DSQRT(GammaPrim*mass*W)
      bb=W*W-mass*mass
      cc=mass*GammaPrim
      dd=B
C      nrbw = ( (aa*bb)/(bb*bb+cc*cc) + dd )**2
C      nrbw = nrbw + ( (aa*cc)/(bb*bb+cc*cc) )**2
C  Put in a simple, no-background BW, following J. Breitweg et al.
C  Eq. 15 of Eur. Phys. J. C2, 247 (1998).  SRK 11/10/2000

      nrbw = (ANORM*mass*GammaPrim/(bb*bb+cc*cc))**2
      nrbw = C*nrbw

      RETURN
      END


