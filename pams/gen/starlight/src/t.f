c	JS	This code calculates the nuclear thickness function as per Eq. 4 in
c	Klein and Nystrand, PRC 60.

      DOUBLE PRECISION FUNCTION T(b)

      IMPLICIT NONE
      DOUBLE PRECISION b
      DOUBLE PRECISION zsp,zmin,zmax,r,sum,rws
      DOUBLE PRECISION xg(1:5),ag(1:5)
      INTEGER I,NGAUSS
      PARAMETER(zmin=0.0)
      PARAMETER(zmax=15.0)

      PARAMETER(NGAUSS=5)
C     >> DATA FOR GAUSS INTEGRATION
      DATA xg/0.1488743390,0.4333953941,0.6794095683,0.8650633667,
     $        0.9739065285/
      DATA ag/0.2955242247,0.2692667193,0.2190863625,0.1494513492,
     $        0.0666713443/

      sum=0.
      DO 100 I=1,NGAUSS

        zsp = 0.5*(zmax-zmin)*xg(I) + 0.5*(zmax+zmin)
        r   = DSQRT(b*b+zsp*zsp)
        sum = sum + ag(I)*rws(r)
        zsp = 0.5*(zmax-zmin)*(-xg(I)) + 0.5*(zmax+zmin)
        r   = DSQRT(b*b+zsp*zsp)
        sum = sum + ag(I)*rws(r)

 100  CONTINUE
      sum = 0.5*(zmax-zmin)*2.*sum

      T=sum

      RETURN
      END
