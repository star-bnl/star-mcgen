* $Id: hijjet.F,v 1.1 1998/04/20 18:19:54 didenko Exp $
* $Log: hijjet.F,v $
* Revision 1.1  1998/04/20 18:19:54  didenko
* version hijing1.34 released
*
*
      INTEGER*4 FUNCTION HIJJET(par_h, par) 
      IMPLICIT NONE
*:>--------------------------------------------------------------------
*: ROUTINE:    HIJING
*: generate one hijing event and put into table
*:>--------------------------------------------------------------------
      EXTERNAL  HIJEV
#include "hijjet.inc"
#include "hicdes/headpss.inc"
#include "hicdes/hepevt.inc"
#include "hicdes/himevt2.inc"
*--
      INTEGER NATT, JATT, NT, NP, N0, N01, N10, N11
      REAL  EATT
      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/
*
      INTEGER I, J, IK, IP
      DOUBLE PRECISION pp, pmass
*
       CALL HIJEV     
*
*-- Sanity check of header for table "par"...
*
      par_h.nok = 0
      IF(NATT .GT. par_h.maxlen) THEN
         HIJJET = STAFCV_BAD
         GO TO 999
      ENDIF
*
*-- Start event loop
*
        ik = 1
*
        CALL DATIME(I,J)
        VSSHEP(2) = VSSHEP(2) + 1
        VSSHEP(3) = FLOAT(I)
        VSSHEP(4) = FLOAT(J)
        par(1).isthep = 10
        par(1).idhep = 9999999
        par(1).jmohep(1) = 0
        par(1).jmohep(2) = 0
        par(1).jdahep(1) = 0
        par(1).jdahep(2) = 0
        par(1).phep(1) = 0.0
        par(1).phep(2) = 0.0
        par(1).phep(3) = psshep(3)
        par(1).phep(4) = psshep(4)
        par(1).phep(5) = psshep(5)
        par(1).vhep(1) = vsshep(1)
        par(1).vhep(2) = vsshep(2)
        par(1).vhep(3) = vsshep(3)
        par(1).vhep(4) = vsshep(4)
*--
      DO IP = 1, NATT
        if(katt(ip,4) .eq. 1) then
        ik = ik + 1
        par(ik).isthep     = katt(ip,4)    !  status code
        par(ik).idhep      = katt(ip,1)    !  particle identity
        par(ik).jmohep(1)  = katt(ip,3)    !  pointer(s) to position 
        par(ik).jmohep(2)  = 0             !  where the mothers(s) stored   
        par(ik).jdahep(1)  = 0             !  pointers to position of 
        par(ik).jdahep(2)  = 0             !  the first/last daughter  
        DO J = 1, 4
          par(ik).phep(j)  = patt(ip,j)    !  p4 
          par(ik).vhep(j)  = 0.            !  production vertex (mm) and time (mm/c) 
        END DO
        pp = patt(ip,4)**2 - patt(ip,1)**2 - patt(ip,2)**2 - patt(ip,3)**2
        pmass = dsqrt(pp)   
        par(ik).phep(5)    =  pmass        !  mass (GeV)
      endif  
      END DO
      par_h.nok = ik
*
      HIJJET=STAFCV_OK
  999 END
