* $Id: hijjet.f,v 1.1 2000/06/16 14:58:04 longacre Exp $
* $Log: hijjet.f,v $
* Revision 1.1  2000/06/16 14:58:04  longacre
* hijing ntuple maker
*
* Revision 1.6  2000/01/22 17:37:41  nevski
* clean up unused variables
*
* Revision 1.5  1999/08/14 21:49:00  fisyak
* eliminate include and source level
*
* Revision 1.4  1999/08/13 23:55:08  fisyak
* remove hepevt
*
* Revision 1.3  1999/01/16 20:53:22  didenko
* add some usefull histogram for jets
*
* Revision 1.2  1998/05/13 18:55:50  didenko
* modified version with full particles history
*
* Revision 1.1  1998/04/20 18:19:54  didenko
* version hijing1.34 released
*
*
      PROGRAM HIJJET
      IMPLICIT NONE
*:>--------------------------------------------------------------------
*: ROUTINE:    HIJING
*: generate one hijing event and put into table
*:>--------------------------------------------------------------------
      EXTERNAL  HIJEV
CCC#include "headpss.inc"
      COMMON/HEADPSS/PSSHEP(5),VSSHEP(4),IFIRST,IRUN
      REAL PSSHEP,VSSHEP
      INTEGER IFIRST, IRUN
      SAVE /HEADPSS/
CCCC#include "himevt2.inc"
      INTEGER    KATT
      REAL  PATT, VATT
      COMMON/HIMAIN2/KATT(100000,6), PATT(100000,5), VATT(100000,5)
      SAVE  /HIMAIN2/
*--
      INTEGER NATT, JATT, NT, NP, N0, N01, N10, N11
      REAL  EATT
      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/
*
      COMMON/WEVENT/EVEF
      REAL EVEF
      SAVE /WEVENT/
*
      REAL PT, ETA, ptjet, px1, py1, P(3), V(3), E, AMASS, X4
      REAL px, py, pz, theta, HMEMOR 
      INTEGER I, J, IK, IP, imo(2), idau(2)
*
      INTEGER ii, nd1, nd2, kk, mm, ll, km, ISTAT
      INTEGER nd3, nd4, nd5, nd6, nd7, nd8
      INTEGER nd9, nd10, nd11, nd12, njets, njetp
C
       EVEF = 0.
       X4=0.0
       CALL HIJEV
 999   EVEF=EVEF+1.
        
       IF(EVEF .GT. VSSHEP(4)) then     
         CALL HEPEnd(' ')
         STOP
       ENDIF
       IF(EVEF.GT.1) CALL HIJEV
*
*
*
*-- Start event loop
*
*
C
          njets = 0
C 
          do ii = 1, natt
        if(abs(katt(ii,1)).ge.91.and.abs(katt(ii,1)).le.93) then
          pt=SQRT(patt(ii,1)**2+patt(ii,2)**2)
          if(pt .ge. 10.) njets = njets + 1 
          endif
         enddo
C
         WRITE(6,888) njets
 888     FORMAT(' njets ',I10)
         vsshep(3)=float(njets)
C--
       CALL HEPEvent('hijing',IRUN,NATT,vsshep(1),vsshep(2),psshep(5)
      1,vsshep(3),psshep(1),psshep(2),psshep(3),psshep(4))
        ik = 0
        DO IP = 1, NATT
        ik = ik + 1
        ISTAT=katt(ip,4)
        imo(1)=katt(ip,3)
        imo(2)=0
        idau(1)=katt(ip,5)
        idau(2)=katt(ip,6)
        IF(ik.eq.NATT)  ik = -1
        DO J = 1, 3
          P(j)  = patt(ip,j)    
          V(j)  = 0.        
        END DO
          E     = patt(ip,4)    
          AMASS = patt(ip,5)
      CALL HEPPart(ik,ISTAT,katt(ip,1),imo,idau,P,E,AMASS,V,X4) 
CCCC        WRITE(61,6116)ik,ISTAT,katt(ip,1),imo,idau
CCCC     1,P,E,AMASS
 6116   FORMAT(7I6,5G12.5)
      END DO
      GO TO 999
      END
      FUNCTION RNDM(IDUM)
        INTEGER LEN,IDUM 
        REAL*4 RNDM, RVEC(1)         
        SAVE LEN
        DATA LEN/1/
        CALL RANLUX(RVEC,LEN)
        RNDM=RVEC(1)
        RETURN
        END


