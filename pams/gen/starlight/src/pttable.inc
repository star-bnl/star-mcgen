        DOUBLE PRECISION ptarray,fptarray,ptmax,dpt
	INTEGER NPT
      common/pttable/ptarray(1000),fptarray(1000,1000),ptmax,dpt,NPT
C  This is a 2-dimensional (y,pt) lookup table for pt
C    for use when calculating vector meson pt with interference
C  ptarray contains the pt values
C  the y values are in array yarray in common Ftable.inc
C  fptarray(y,pt) is the integral probablity to find a vector meson
C  at rapidity y and transverse momentum pt

