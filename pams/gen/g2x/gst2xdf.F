* Id: gst2xdf.F, v 1.1 1998/03/17 14:35:30 longacre Exp $
*$Log: gst2xdf.F,v $
*Revision 1.2  1998/04/03 20:33:05  longacre
*correct exit for too many particles
*
*Revision 1.1  1998/03/25 16:49:16  longacre
*text2xdf tables
* 
      INTEGER FUNCTION GST2XDF(par_h, par) 
  	Implicit None
*:>--------------------------------------------------------------------
*: ROUTINE:    GST2XDF
*: GET one text ASCII FILE event and put into table
*:>--------------------------------------------------------------------
#include "gst2xdf.inc"
#include "g2xcdes/hepevt.inc"
C---------------------------------------------------------------------
      CHARACTER*120 LINE 
      INTEGER LENOCC,INDEX,LI,IEVEN,NTRAC,NVERT,ITR,IVT,
     *LABELTR,LABELVX,GE_PID,EG_PID,STARTVX,STOPVX,I, EG_PROC,PARENT
      REAL VERSION,EAST_Z,EAST_A,WEST_Z,WEST_A,SQRTS,B_MAX,PP(3),VERT(4)
     *,A,B 
      INTEGER IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN,IEOTRI,
     *IEVENT,ISWIT,IFINIT,NEVENT,NRNDM
      COMMON/GCFLAG/ IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN,
     *IEOTRI, IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      DOUBLE PRECISION hepmass
      DOUBLE PRECISION ionmass
      DOUBLE PRECISION masss
      INTEGER I
      INTEGER J
      INTEGER IFIRST
      INTEGER IDTRAN
      INTEGER ICODE
      DATA IFIRST/0/
      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         Open(UNIT=22,FILE='evgen.evt',READONLY)
      ENDIF
      par_h.nok = 0
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
1011     IF((ITR.GE.NTRAC).AND.(IVT.GE.NVERT)) GO TO 1012
         LINE='end-of-file' 
         READ(LI,'(a)',ERR=  1020)LINE
         WRITE(*,1030)LINE
1030     FORMAT(1X,A)
C    * 
C       Check LINE(1:5)=='GENER' 
            IF (LINE(1:5).EQ.'GENER') THEN
            READ(LINE(11:),*,ERR=1020)VERSION,EAST_Z,EAST_A,WEST_Z,
     *      WEST_A,SQRTS,B_MAX
         WRITE(*,1040)VERSION,EAST_Z,EAST_A,WEST_Z,WEST_A,SQRTS,B_MAX 
1040        FORMAT(' gstar_ReadNew: GENER :',F8.3,4F8.0,F8.1,F8.3)
         ELSE IF(LINE(1:5).EQ.'EVENT') THEN 
            READ(LINE(7:),*,END=  1050,ERR=1020)IEVEN,NTRAC,NVERT,A,B
            WRITE(*,1060)IEVEN,NTRAC,NVERT,A,B
1060        FORMAT(' gstar_ReadNew: EVENT :',3I8,2F8.3)
C        end of data
1050        IF (IEVEN.LE.-999) GOTO 1020 
         ELSE IF(LINE(1:5).EQ.'TRACK') THEN 
            READ(LINE(7:),*,ERR=1020)GE_PID,PP,LABELTR,STARTVX,STOPVX,
     *      EG_PID  
            WRITE(*,1070)GE_PID,PP,LABELTR,STARTVX,STOPVX,EG_PID 
1070        FORMAT(16X,'TRACK :',I6,3F9.3,4I6)
            ITR=ITR+1 
         ELSE IF(LINE(1:6).EQ.'VERTEX') THEN 
            READ(LINE(8:),*,ERR=1020)VERT,LABELVX,EG_PROC,PARENT 
            WRITE(*,1080)VERT,LABELVX,EG_PROC,PARENT
1080        FORMAT(16X,'VERTEX:',4F10.6,3I6)
            IVT=IVT+1
         ELSE IF(INDEX(LINE,'event').GT.0 .AND. ITR+IVT.EQ.0) THEN
            I=INDEX(LINE,'event')
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
        masss    = hepmass(idhep)         !  mass (GeV)  
        ELSE
        masss    = ionmass(GE_PID)        !  mass (GeV)  
        ENDIF
        phep(4)=dsqrt(phep(1)**2+phep(2)**2+phep(3)**2+masss**2)
        par(ITR).isthep     = ICODE       !  status code
        par(ITR).idhep      = idhep       !  particle identity
        par(ITR).jmohep(1)  = 0           !  pointer(s) to position 
        par(ITR).jmohep(2)  = 0           !  where the mothers(s) stored   
        par(ITR).jdahep(1)  = 0           !  pointers to position of 
        par(ITR).jdahep(2)  = 0           !  the first/last daughter  
        DO J = 1, 4
          par(ITR).phep(j)  = phep(j)     !  p4 
          par(ITR).vhep(j)  = 0.0         !  production vertex (mm) and time (mm/c) 
        END DO
        par(ITR).phep(5)    = masss       !  mass (GeV)  
1101        CONTINUE 
      par_h.nok = NTRAC 
      GST2XDF=STAFCV_OK
1102        CONTINUE
            GO TO 1012
         ELSE IF(LENOCC(LINE).GT.0) THEN  
            WRITE(*,1120)LINE(1:LENOCC(LINE))
1120        FORMAT(' unknown line : ',A)
         END IF 
      GO TO 1011
1012  CONTINUE 
      GO TO 1 
1020  WRITE(*,1130)LINE 
1130  FORMAT(' gstar_readTXT error in line '/1X,A) 
      GST2XDF  = STAFCV_BAD
1     CONTINUE 
      RETURN
999   WRITE(*,1131)
1131  FORMAT('gst2xdf: SIZE problem  ')
      GST2XDF  = STAFCV_BAD
      End



