C     >> Function for the gamma-proton --> VectorMeson 
C     >> cross section. Wgp is the gamma-proton CM energy.
C     >> Unit for cross section: fm**2
      DOUBLE PRECISION FUNCTION sigmagp(Wgp)
 
      IMPLICIT NONE

      include 'inputp.inc'
      DOUBLE PRECISION Wgp

      IF(ip.eq.113)THEN
c	these are rhos
       sigmagp=1.E-4*(5.0*DEXP(0.22*DLOG(Wgp))+26.0*
     &		DEXP(-1.23*DLOG(Wgp)))
      ELSEIF(ip.eq.223)THEN
c	these are omegas
       sigmagp=1.E-4*(0.55*DEXP(0.22*DLOG(Wgp))+18.0*
     &		DEXP(-1.92*DLOG(Wgp)))
      ELSEIF(ip.eq.333)THEN
c	these are phis
        sigmagp=1.E-4*0.34*DEXP(0.22*DLOG(Wgp))
      ELSEIF(ip.eq.443)THEN
c	these are J/psis
        sigmagp=1.E-4*0.0015*DEXP(0.80*DLOG(Wgp))
      ELSE
        WRITE(*,*) 'ERROR: Unidentified Meson:',ip
      ENDIF

      RETURN
      END
