* $Id: hijjet.F,v 1.2 1998/05/13 18:55:50 didenko Exp $
* $Log: hijjet.F,v $
* Revision 1.2  1998/05/13 18:55:50  didenko
* modified version with full particles history
*
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
      COMMON/WEVENT/IEVEF
      INTEGER IEVEF
      SAVE /WEVENT/
*
      REAL px, py, pz, pt, eta, theta, phi, pe, yyr, HMEMOR 
      INTEGER I, J, IK, IP
*
      COMMON/PAWC/HMEMOR(20000)
      CHARACTER*6 HTAGP(2)
      CHARACTER*6 HTAGE(2)
      CHARACTER*6 ITAGN(2)
      DIMENSION PIEV(2)
      DIMENSION PREV(2)
      DIMENSION NMEV(2)
      REAL PIEV, PREV, NMEV, NMPR, NMPI
      EQUIVALENCE (PIEV(1), PTPI), (PIEV(2), ETHAPI)
      EQUIVALENCE (PREV(1), PTPR), (PREV(2), ETHAPR)
      EQUIVALENCE (NMEV(1), NMPI), (NMEV(2), NMPR)
      REAL PTPR, PTPI, ETHAPR, ETHAPI, PIEV, PREV
      DATA HTAGP/'PTPI','ETHAPI'/
      DATA HTAGE/'PTPR','ETHAPR'/
      DATA ITAGN/'NMPI','NMPR'/ 
      LOGICAL FIRST/.TRUE./
**
      IF (FIRST) THEN
       FIRST=.FALSE.
       CALL HBOOKN(101,'PARTICLE PT',2,' ',1000, htagp)
       CALL HBOOKN(102,'PARTICLE ETHA',2,' ',1000, htage)
       CALL HBOOKN(103,'PARTICLE MULPYPLICITY',2,' ',1000,itagn)
      ENDIF
C
       IEVEF = 0
       CALL HIJEV
       IF(IEVEF .eq. 0)then     
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
        nmev(1) = 0.
        nmev(2) = 0.
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
        par(1).phep(1) = psshep(1) 
        par(1).phep(2) = psshep(2)
        par(1).phep(3) = psshep(3)
        par(1).phep(4) = psshep(4)
        par(1).phep(5) = psshep(5)
        par(1).vhep(1) = vsshep(1)
        par(1).vhep(2) = vsshep(2)
        par(1).vhep(3) = vsshep(3)
        par(1).vhep(4) = vsshep(4)
*--
      DO IP = 1, NATT
        ik = ik + 1
        par(ik).isthep     = katt(ip,4)    !  status code
        par(ik).idhep      = katt(ip,1)    !  particle identity
        par(ik).jmohep(1)  = katt(ip,3)    !  pointer(s) to position 
        par(ik).jmohep(2)  = 0             !  where the mothers(s) stored   
        par(ik).jdahep(1)  = katt(ip,5)    !  pointers to position of 
        par(ik).jdahep(2)  = katt(ip,6)    !  the first/last daughter  
        DO J = 1, 4
          par(ik).phep(j)  = patt(ip,j)    !  p4 
          par(ik).vhep(j)  = vatt(ip,j)    !  decay vertex (mm) and time(mm/c) 
        END DO
          par(ik).phep(5) = patt(ip,5)
C
             px = patt(ip,1)
             py = patt(ip,2)
             pz = patt(ip,3)
             pe = patt(ip,4)
C
             pt = sqrt(px*px+py*py)
             theta = atan2 ( pt, pz )
             eta   = -log(tan(theta/2.))
             yyr = 0.5*log((pe + pz)/(pe - pz))                   
             phi  = atan2 ( py, px )
C                 
             if(abs(katt(ip,1)) .eq. 211) then
                  nmev(1) = nmev(1) + 1
                  piev(1) = pt
                  piev(2) = eta
                  call hfn(101, piev)
             elseif(abs(katt(ip,1)) .eq. 2212) then
                  nmev(2) = nmev(2) + 1
                  prev(1) = pt
                  prev(2) = eta
                  if(abs(eta).lt.20.)then
                  call hfn(102, prev)
                  endif
             endif
C  
      END DO
                 call hfn(103, nmev)
      par_h.nok = ik
*
      HIJJET=STAFCV_OK
      ENDIF
C
  999 END


