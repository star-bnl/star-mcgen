      DOUBLE PRECISION FUNCTION nrbw(W,ANORM,B,C)

C     Relativistic Breit-Wigner according to J.D. Jackson,
C     Nuovo Cimento, 1964, with nonresonant term. A is the strength
C     of the resnonant term and b the strength of the non-resonant
C     term. C is an overall normalization. 

      IMPLICIT NONE

      include 'global.inc'
      include 'const.inc'
      DOUBLE PRECISION W,ANORM,B,C
      DOUBLE PRECISION ppi,ppi0,GammaPrim,rat
      DOUBLE PRECISION aa,bb,cc,dd

      ppi=DSQRT( (W/2.)**2 - mpi*mpi )
      ppi0=0.358
      rat=ppi/ppi0
      GammaPrim=width*(mass/W)*rat*rat*rat
      aa=ANORM*DSQRT(GammaPrim*mass*W)
      bb=W*W-mass*mass
      cc=mass*GammaPrim
      dd=B
      nrbw = ( (aa*bb)/(bb*bb+cc*cc) + dd )**2
      nrbw = nrbw + ( (aa*cc)/(bb*bb+cc*cc) )**2
      nrbw = C*nrbw

      RETURN
      END
