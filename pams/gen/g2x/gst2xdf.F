* Id: gst2xdf.F, v 1.1 1998/03/17 14:35:30 longacre Exp $
*$Log: gst2xdf.F,v $
*Revision 1.9  1999/10/08 22:24:11  nevski
*pgf->g77
*
*Revision 1.8  1999/10/08 16:31:16  nevski
*no comments anymore
*
*Revision 1.7  1999/10/08 15:50:10  nevski
*pgf->g77 cleanup
*
*Revision 1.6  1999/08/13 21:54:15  fisyak
*move hepevt
*
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
#include "hepevt.inc"
C---------------------------------------------------------------------
C      IVENT records the number of times an 'EVENT' flag was found.
C      Actually it gets incremented TWICE per event after the first one.
C  IEVTREAD = 0   Need new data
C           = 1   Event in :  not processed
C           = 2   Event read and processed

      CHARACTER*120 LINE 
      INTEGER LENOCC,INDEX,LI,IEVEN,NTRAC,NVERT,ITR,IVT,
     >LABELTR,LABELVX,GE_PID,EG_PID,STARTVX,STOPVX,EG_PROC,PARENT   
      REAL VERSION,EAST_Z,EAST_A,WEST_Z,WEST_A,SQRTS,B_MAX,PP(3),
     >     VERT(4),A,B 
      INTEGER IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN,IEOTRI,
     >IEVENT,ISWIT,IFINIT,NEVENT,NRNDM
      COMMON/GCFLAG/ IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN,
     >IEOTRI, IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
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
C -------------------------------------------------------------------

      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         Open(UNIT=22,FILE='evgen.evt',status='old')
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
111     CONTINUE
C********IF((ITR.GE.NTRAC).AND.(IVT.GE.NVERT)) GO TO 112
         IF(IEVTREAD.EQ.10) STOP  ! Quit if EOF had been found
C  Only read from disk if a new record is needed.
         IF(IEVTREAD.EQ.0) THEN 
             READ(LI,'(a)',ERR=120,END=2000)LINE
             IEVTREAD=1
         ENDIF
         WRITE(*,130)LINE
130      FORMAT(1X,A)
        IF(LINE(1:1).EQ.'*') THEN
          IEVTREAD=0
          GO TO 111
        ENDIF
        IF(LINE(1:1).EQ.' ') THEN
         IF(INDEX(LINE,'event').GT.0 .AND. ITR+IVT.EQ.0) THEN
C  No keywords.  Must be "OLD" format.
C  Read in the event and output it.  
            I=INDEX(LINE,'event')
            IEVTREAD=0
            LINE(I:I+6)='  ' 
            READ(LINE,*,ERR=120)NTRAC,IEVEN 
            WRITE(*,190)NTRAC,IEVEN 
190         FORMAT(' gstar_ReadOld: ',I8,' event# ',I6)
            DO 201 ITR=1,NTRAC 
               READ(LI,*,ERR=120)GE_PID,phep(1), phep(2), phep(3) 
               idhep=idtran(GE_PID)
CCCCCCC               WRITE(*,210)idhep, phep(1), phep(2), phep(3)
210            FORMAT(16X,I7,3G10.3)
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
201         CONTINUE 
            par_h.nok = NTRAC 
            GST2XDF=STAFCV_OK
202         CONTINUE
            GO TO 112
         ENDIF
          IEVTREAD=0
          GO TO 111
        ENDIF
C    * 
C       Check LINE(1:5) for FORMAT key word 
        IF (LINE(1:5).EQ.'GENER') THEN
            READ(LINE(11:),*,ERR=120)VERSION,EAST_Z,EAST_A,WEST_Z,
     *       WEST_A,SQRTS,B_MAX
            WRITE(*,140)VERSION,EAST_Z,EAST_A,WEST_Z,WEST_A,SQRTS,B_MAX 
140         FORMAT(' gstar_ReadNew: GENER :',F8.3,4F8.0,F8.1,F8.3)
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
            WRITE(*,145) IEVTREAD,IVENT, ITRSAV
 145        FORMAT(' IEVTREAD,IVENT,ITRSAV = ', 3i7)
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
10         READ(LINE(7:),*,END=  122,ERR=122)IEVEN,NTRAC,NVERT,A,B
           GO TO 123
122        READ(LINE(7:),*,END=  150,ERR=120)IEVEN,NTRAC,NVERT
123        WRITE(*,160) IEVEN,NTRAC,NVERT,A,B
160        FORMAT(' gstar_ReadNew: EVENT :',3I8,2F8.3, ' changed ')
            IEVTREAD=0
150        IF (IEVEN.LT.-999) GOTO 120    ! Error on read
            IF (IEVEN.EQ.-999) THEN
C              end of data       
               WRITE( *, 155)
 155           FORMAT(' Found EVENT -999. Terminating. -999 not written')
               STOP
            ENDIF
         ELSE IF(LINE(1:5).EQ.'TRACK') THEN 
            IEVTREAD=0
            READ(LINE(7:),*,ERR=120)GE_PID,PP,LABELTR,STARTVX,STOPVX,
     *      EG_PID  
            WRITE(*,170) GE_PID,PP,LABELTR,STARTVX,STOPVX,EG_PID 
170         FORMAT(16X,'TRACK :',I6,3F9.3,4I6)
            ITR=ITR+1 
            IF(LABELVX.NE.0) THEN
            IF(STARTVX.NE.LABELVX) THEN
               WRITE(*,172) STARTVX, LABELVX
172            FORMAT(' VERTEX GOES BEFORE TRACK AND MUST MATCH',/
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
CCCCCCC               WRITE(*,210)idhep, phep(1), phep(2), phep(3)
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
            READ(LINE(8:),*,ERR=120)VERT,LABELVX,EG_PROC,PARENT 
            WRITE(*,180) VERT,LABELVX,EG_PROC,PARENT
180         FORMAT(16X,'VERTEX:',4F10.6,3I6)
            IVT=IVT+1
         ELSE IF(LENOCC(LINE).GT.0) THEN  
            WRITE(*,220) LINE(1:LENOCC(LINE))
220         FORMAT(' unknown line : ',A)
            IEVTREAD=0
         END IF 
      GO TO 111
112   CONTINUE 
      return

C   Error Returns
120   WRITE(*,230)LINE 
230   FORMAT(' gstar_readTXT error in line '/1X,A) 
      GST2XDF  = STAFCV_BAD
      RETURN

999   WRITE(*,231)
231  FORMAT('gst2xdf: SIZE problem  ')
      GST2XDF  = STAFCV_BAD
      RETURN

C  Come here at end of tape.  Set up to write out the last event.
 2000 par_h.nok=ITRSAV
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











