* Id: gst2xdf.F, v 1.1 1998/03/17 14:35:30 longacre Exp $
*$Log: gst2xdf.F,v $
*Revision 1.5  1998/08/02 21:38:15  longacre
*added vertex cm -> mm
*
*Revision 1.4  1998/07/01 19:58:37  longacre
*make g2x read Lanny's txt files
*
*Revision 1.3  1998/05/29 14:31:18  longacre
*now makes correct particle for new txtfile
*
*Revision 1.3  1998 /5/25  (M. Kramer)
*Add new formatting  (with keywords) 
*Mostly start new event when the keyword 'EVENT' is found.
*Add IEVENT to keep track of the number of events found,
*and IEVTREAD to keep track if there is unprocessed data in LINE.
*Revision 1.2  1998/04/03 20:33:05  longacre
*correct exit for too many particles
*
*Revision 1.1  1998/03/25 16:49:16  longacre
*text2xdf tables
* 
      INTEGER FUNCTION GST2XDF(par_h, par) 
  	Implicit None
C  Still need to do end ofg event for new format.
C
*:>--------------------------------------------------------------------
*: ROUTINE:    GST2XDF
*: GET one text ASCII FILE event and put into table
*:>--------------------------------------------------------------------
#include "gst2xdf.inc"
#include "g2xcdes/hepevt.inc"
C---------------------------------------------------------------------
C      IVENT records the number of times an 'EVENT' flag was found.
C      Actually it gets incremented TWICE per event after the first one.
C  IEVTREAD = 0   Need new data
C           = 1   Event in :  not processed
C           = 2   Event read and processed

      CHARACTER*120 LINE 
      INTEGER LENOCC,INDEX,LI,IEVEN,NTRAC,NVERT,ITR,IVT,
     *LABELTR,LABELVX,GE_PID,EG_PID,STARTVX,STOPVX,I, EG_PROC,PARENT   
      REAL VERSION,EAST_Z,EAST_A,WEST_Z,WEST_A,SQRTS,B_MAX,PP(3),VERT(4)
     *,A,B 
      INTEGER IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN,IEOTRI,
     *IEVENT,ISWIT,IFINIT,NEVENT,NRNDM
      COMMON/GCFLAG/ IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN,
     *IEOTRI, IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      SAVE IFIRST,IEVTREAD,ITR
      DOUBLE PRECISION hepmass
      DOUBLE PRECISION ionmass
      DOUBLE PRECISION masss
      INTEGER I,J
      INTEGER IFIRST
      INTEGER IDTRAN
      INTEGER ICODE
      INTEGER IEVTREAD, IVENT, ITRSAV
      DATA IFIRST/0/, IVENT/0/, IEVTREAD/0/, ITRSAV/0/
      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         Open(UNIT=22,FILE='evgen.evt',READONLY)
      ENDIF
      par_h.nok = 0
C      ITRSAV=ITR
      IF(NMXHEP .GT. par_h.maxlen) THEN
         GO TO 999
      ENDIF
      LI=22
      NTRAC=999
      NVERT=999
      ITR=0
      IVT=0
C 
C    Loop here 
1011     CONTINUE
C********IF((ITR.GE.NTRAC).AND.(IVT.GE.NVERT)) GO TO 1012
         IF(IEVTREAD.EQ.10) STOP  ! Quit if EOF had been found
C  Only read from disk if a new record is needed.
         IF(IEVTREAD.EQ.0) THEN 
             READ(LI,'(a)',ERR=1020,END=2000)LINE
             IEVTREAD=1
         ENDIF
         WRITE(*,1030)LINE
1030     FORMAT(1X,A)
        IF(LINE(1:1).EQ.'*') THEN
          IEVTREAD=0
          GO TO 1011
        ENDIF
        IF(LINE(1:1).EQ.' ') THEN
         IF(INDEX(LINE,'event').GT.0 .AND. ITR+IVT.EQ.0) THEN
C  No keywords.  Must be "OLD" format.
C  Read in the event and output it.  
            I=INDEX(LINE,'event')
            IEVTREAD=0
            LINE(I:I+6)='  ' 
            READ(LINE,*,ERR=1020)NTRAC,IEVEN 
            WRITE(*,1090)NTRAC,IEVEN 
1090        FORMAT(' gstar_ReadOld: ',I8,' event# ',I6)
            DO 1101 ITR=1,NTRAC 
               READ(LI,*,ERR=1020)GE_PID,phep(1), phep(2), phep(3) 
               idhep=idtran(GE_PID)
CCCCCCC               WRITE(*,1110)idhep, phep(1), phep(2), phep(3)
1110           FORMAT(16X,I7,3G10.3)
               ICODE = 1
               IF(IABS(idhep).LT.7) ICODE=3
               IF(idhep.EQ.21) ICODE=3
               IF(idhep.LT.6999999) THEN
                  masss    = hepmass(idhep)      !  mass (GeV)  
               ELSE
                  masss    = ionmass(GE_PID)        !  mass (GeV)  
               ENDIF
               phep(4)=dsqrt(phep(1)**2+phep(2)**2+phep(3)**2+masss**2)
               par(ITR).isthep     = ICODE    !  status code
               par(ITR).idhep      = idhep    !  particle identity
               par(ITR).jmohep(1)  = 0        !  pointer(s) to position 
               par(ITR).jmohep(2)  = 0        !  where the mothers(s) stored   
               par(ITR).jdahep(1)  = 0        !  pointers to position of 
               par(ITR).jdahep(2)  = 0        !  the first/last daughter  
               DO J = 1, 4
                 par(ITR).phep(j)  = phep(j)  !  p4 
                 par(ITR).vhep(j)  = 0.0 !production vertex(mm) and time(mm/c) 
               END DO
               par(ITR).phep(5)    = masss    !  mass (GeV)  
1101        CONTINUE 
            par_h.nok = NTRAC 
            GST2XDF=STAFCV_OK
1102        CONTINUE
            GO TO 1012
         ENDIF
          IEVTREAD=0
          GO TO 1011
        ENDIF
C    * 
C       Check LINE(1:5) for FORMAT key word 
        IF (LINE(1:5).EQ.'GENER') THEN
            READ(LINE(11:),*,ERR=1020)VERSION,EAST_Z,EAST_A,WEST_Z,
     *       WEST_A,SQRTS,B_MAX
            WRITE(*,1040)VERSION,EAST_Z,EAST_A,WEST_Z,WEST_A,SQRTS,B_MAX 
1040        FORMAT(' gstar_ReadNew: GENER :',F8.3,4F8.0,F8.1,F8.3)
            IEVTREAD=0
         ELSE IF(LINE(1:5).EQ.'EVENT') THEN 
C
C  We assume this is a new event.
C  So first finish off previous event.
C  However if this is the first event there is nothing to write.
C  So test for first event.
C  Set flag to know an EVENT record was read and is stored in 
c   'LINE' and not to reread it from disk for the next event.
            IVENT=IVENT+1
            WRITE(*,1045) IEVTREAD,IVENT, ITRSAV
 1045       FORMAT(' IEVTREAD,IVENT,ITRSAV = ', 3i7)
            IF(IEVTREAD.EQ.1) THEN ! Just read new EVENT from disk
                IEVTREAD=2
                IF(IVENT.GT.1) THEN ! Write previous event each time 
*                                   !an EVENT flag is seen. (Skip for
*                                   ! first flag because there is no
*                                   !previous event).
                   par_h.nok=ITRSAV
                   GST2XDF= STAFCV_OK
                   RETURN
                 ENDIF
            ENDIF
10         READ(LINE(7:),*,END=  1022,ERR=1022)IEVEN,NTRAC,NVERT,A,B
           GO TO 1023
1022       READ(LINE(7:),*,END=  1050,ERR=1020)IEVEN,NTRAC,NVERT
1023       WRITE(*,1060)IEVEN,NTRAC,NVERT,A,B
1060        FORMAT(' gstar_ReadNew: EVENT :',3I8,2F8.3, ' changed ')
            IEVTREAD=0
1050        IF (IEVEN.LT.-999) GOTO 1020    ! Error on read
            IF (IEVEN.EQ.-999) THEN
C              end of data       
               WRITE( *, 1055)
 1055          FORMAT(' Found EVENT -999. Terminating. -999 not written')
               STOP
            ENDIF
         ELSE IF(LINE(1:5).EQ.'TRACK') THEN 
            IEVTREAD=0
            READ(LINE(7:),*,ERR=1020)GE_PID,PP,LABELTR,STARTVX,STOPVX,
     *      EG_PID  
            WRITE(*,1070)GE_PID,PP,LABELTR,STARTVX,STOPVX,EG_PID 
1070        FORMAT(16X,'TRACK :',I6,3F9.3,4I6)
            ITR=ITR+1 
            IF(LABELVX.NE.0) THEN
            IF(STARTVX.NE.LABELVX) THEN
               WRITE(*,1072) STARTVX, LABELVX
1072           FORMAT(' VERTEX GOES BEFORE TRACK AND MUST MATCH',/
     *         ,' TRACK VERTEX IS ',I10,'VERTEX IS ',I10)
               STOP
            ENDIF
            ENDIF
            ITRSAV=ITR
            phep(1)=pp(1)
            phep(2)=pp(2)
            phep(3)=pp(3)
            idhep=idtran(GE_PID)
c
C  A test write commented out.
CCCCCCC               WRITE(*,1110)idhep, phep(1), phep(2), phep(3)
            ICODE = 1
            IF(IABS(idhep).LT.7) ICODE=3
            IF(idhep.EQ.21) ICODE=3
            IF(idhep.LT.6999999) THEN
                  masss    = hepmass(idhep)      !  mass (GeV)  
            ELSE
                  masss    = ionmass(GE_PID)        !  mass (GeV)  
            ENDIF
            phep(4)=dsqrt(phep(1)**2+phep(2)**2+phep(3)**2+masss**2)
            par(ITR).isthep     = ICODE    !  status code
            par(ITR).idhep      = idhep    !  particle identity
            par(ITR).jmohep(1)  = 0        !  pointer(s) to position 
            par(ITR).jmohep(2)  = 0        !  where the mothers(s) stored   
            par(ITR).jdahep(1)  = 0        !  pointers to position of 
            par(ITR).jdahep(2)  = 0        !  the first/last daughter  
            DO J = 1, 4
              par(ITR).phep(j)  = phep(j)  !  p4 
              par(ITR).vhep(j)  = VERT(j)*10. ! vertex(mm) and time(mm/c)
            END DO
            par(ITR).phep(5)    = masss    !  mass (GeV)   
         ELSE IF(LINE(1:6).EQ.'VERTEX') THEN 
            IEVTREAD=0
            READ(LINE(8:),*,ERR=1020)VERT,LABELVX,EG_PROC,PARENT 
            WRITE(*,1080)VERT,LABELVX,EG_PROC,PARENT
1080        FORMAT(16X,'VERTEX:',4F10.6,3I6)
            IVT=IVT+1
         ELSE IF(LENOCC(LINE).GT.0) THEN  
            WRITE(*,1120)LINE(1:LENOCC(LINE))
1120        FORMAT(' unknown line : ',A)
            IEVTREAD=0
         END IF 
      GO TO 1011
1012  CONTINUE 
      GO TO 1 

C   Error Returns
1020  WRITE(*,1130)LINE 
1130  FORMAT(' gstar_readTXT error in line '/1X,A) 
      GST2XDF  = STAFCV_BAD
1     CONTINUE 
      RETURN
999   WRITE(*,1131)
1131  FORMAT('gst2xdf: SIZE problem  ')
      GST2XDF  = STAFCV_BAD
      RETURN
C  Come here at end of tape.  Set up to write out the last event.
2000  par_h.nok=ITRSAV
      GST2XDF=STAFCV_OK
 2005 WRITE(*,2006)
 2006 FORMAT( ' Found EOF ')
C  Is this legitimate end?   i.e. after event -999.
      IF (IEVEN.NE.-999)  Write (*,2010)
 2010 FORMAT( ' Found EOF before event -999  ')
C   Tell program EOF was found
      IEVTREAD=10
      CONTINUE
      RETURN
      End











