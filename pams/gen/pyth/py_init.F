      SUBROUTINE PY_INIT (LUN)
************************************************************************
*DOC PY_INIT                                                           *
*DOC  Operation   :                                                    *
*DOC  Description :                                                    *
*DOC  Reads :        PYGIVE cards                                      *
*DOC  Changes :                                                        *
*DOC  Sends :                                                          *
*DOC  Assumes :                                                        *
*DOC  Results :                                                        *
*DOC  Errors      : None                                               *
*                                                                      *
*DOC  References  : None                                               *
*DOC  Created     : 15-MAR-1997   Author : Yuri Fisyak                 *
************************************************************************
*                                                                      *
      IMPLICIT NONE
*
* Formal parameters
*____________________
*
      INTEGER   LUN
*
* Global sections
*__________________
*
      EXTERNAL PYDATA
#include "geant321/gcflag.inc"
      INTEGER         MSTU, MSTJ
#include "pydat1.inc"
      INTEGER         MDCY ,MDME, KFDP
#include "pydat3.inc"
      INTEGER         PYCOMP
*
* Local variables
*__________________
*
      CHARACTER*4 FRAME, BEAM, TARG, TAG
      CHARACTER*80 FNAME, C_FNAME, CERN_ROOT
      CHARACTER*4 KEYSTR
      INTEGER   LINP, NCHW, LENNOC, LENOCC
      INTEGER   JTDKY, JTEVT, JTCOM, JTLIS, IFL
      CHARACTER       CBUF*80
      INTEGER         IER, KC, IFCH, IDC, I, KF, KF1, KF2, KF3
      INTEGER         KF4, KF5, M1S, M2S 
      INTEGER         LINP
*
*     General data card variables
*     ---------------------------
*
*
      INTEGER       IFRAME, IBEAM, ITARG
      LOGICAL       PYTHIA, HERWIG, ISAJET, DTUJET
      INTEGER       KBEAM, KTARG, KPROC
      REAL          PCMS, PBEAM, PTARG, PTMINS
*
* Executable statements
*________________________
*
      LINP = LUN
      IF (LINP .LE. 0) LINP = 5
      WRITE (*, *) ' PY_INIT called with LINP =', LINP
      WRITE (*, *)
#if defined(NEVER)
      IF (KEYSTR .EQ. 'PYTH') THEN        
#endif
*
*
      PYTHIA = .FALSE.
      FRAME  = 'cms' !
      BEAM   = 'p'   !
      TARG   = 'p'   !
      PCMS   = 200.

*
*  **   Read in data cards
*
   10   READ (LINP, '(A)', END=199) CBUF
        TAG = CBUF(1:4)
        CALL CLTOU (TAG)
        IF (TAG(1:3) .EQ. 'END') THEN
*           End of block
          GOTO 199
        END IF
        IF (TAG(1:4) .EQ. 'PYTH') THEN
          PYTHIA = .TRUE.
          READ (CBUF(5:80),*) FRAME, BEAM, TARG, PCMS 
          GO TO 10
        END IF
        IF (.NOT. PYTHIA) GO TO 10
        IF (CBUF(1:1) .EQ. '*' .OR. CBUF(1:6) .EQ. '      ' .OR.
     &      CBUF(2:5) .EQ.'DECK')                           THEN
*           Ignore
        ELSEIF (TAG(1:4) .EQ. 'DECA')                       THEN
*           Decay blocks
          READ (CBUF, '(6X,6I5)') KF, KF1, KF2, KF3, KF4, KF5
          KC     = PYCOMP (KF)
          IFCH   = 0
          IF (KC.GT.0) THEN
            DO IDC = MDCY(KC,2), MDCY(KC,2)+MDCY(KC,3)-1
              IF (KFDP(IDC,1) .NE. KF1 .OR. KFDP(IDC,2) .NE. KF2 .OR.
     &            KFDP(IDC,3) .NE. KF3 .OR. KFDP(IDC,4) .NE. KF4 .OR.
     &            KFDP(IDC,5) .NE. KF5)                     THEN
                IF (MDME(IDC,1) .GT. 0) THEN
                  MDME(IDC,1) = -1
                  BRAT(IDC)   = 0.0
                  IFCH        = IFCH + 1
                ENDIF
              ELSE
                IF (MDME(IDC,1) .NE. 1) THEN
                  MDME(IDC,1) = 1
                  BRAT(IDC)   = 1.0
                  IFCH        = IFCH + 1
                ENDIF
              ENDIF
            ENDDO
            IF (IFCH .EQ. 0) THEN
*
* *** force decay even if there is this mode in JETSET tables
*
              IDC = MDCY(KC,2)
              MDME(IDC,1) = 1
              BRAT(IDC)   = 1.0
              IFCH        = IFCH + 1
              KFDP(IDC,1) = KF1
              KFDP(IDC,2) = KF2
              KFDP(IDC,3) = KF3
              KFDP(IDC,4) = KF4
              KFDP(IDC,5) = KF5
            END IF
            IF (IFCH .GT. 0 .AND. IDEBUG .GT. 0) THEN
              M1S     = MSTU(1)
              M2S     = MSTU(2)
              MSTU(1) = KF
              MSTU(2) = KF
              CALL PYLIST (12)
              MSTU(1) = M1S
              MSTU(2) = M2S
            ENDIF
          ENDIF
        ELSE
*           Other control cards
          I       = INDEX (CBUF, '!')
          IF (I.LE.0) I = 81
          CALL PYGIVE (CBUF(1:I-1))
        ENDIF
        GOTO 10
 199    IF (PYTHIA) THEN
*-----------------------------------------------------------------------
*
*         CALL UHTOC (IFRAME, 4, FRAME, 4)
*         CALL UHTOC (IBEAM , 4, BEAM , 4)
*         CALL UHTOC (ITARG , 4, TARG , 4)
          IF (FRAME .NE. 'none') THEN
            CALL PYINIT(FRAME,BEAM,TARG,DBLE(PCMS)) 
            IFINIT(3) = 1
          END IF
        END IF
#if defined(NEVER)
      ELSEIF (KEYSTR .EQ. 'ISAJ') THEN        
*
* *** Initialize ISAJET
*
        IF (ISAJET) THEN
          CALL GETENVF ('CERN_ROOT',CERN_ROOT)        ! Find value of $CERN_ROOT 
          LENNOC = LENOCC(CERN_ROOT)
          FNAME = CERN_ROOT(1:LENNOC)//'/lib/isajet.dat'
#if defined(D0_VAX)
*
* *** Convert file-name syntax UNIX -> VAX :
*    UNIX:   //node/=log/dir/a/b/c/file.ext;cy
*       or    //node/log/dir/a/b/c/file.ext;cy
*    gives:   node::log:[dir.a.b.c]file.ext;cy
*
*    UNIX:   /=log/dir/a/b/c/file.ext;cy
*       or    /log/dir/a/b/c/file.ext;cy
*    gives:   log:[dir.a.b.c]file.ext;cy
*                 
*    UNIX:   /(dir/a/b/c/file.ext;cy
*    gives:   [dir.a.b.c]file.ext;cy
*
*    UNIX:      a/b/c/file.ext;cy
*    gives:   [.a.b.c]file.ext;cy
*
*    UNIX:      ../a/b/c/file.ext;cy
*    gives:     [-.a.b.c]file.ext;cy
*
*    UNIX:           ~/a/b/c/f.e    or          ~/f.e
*    gives   log:[name.a.b.c]f.e        log:[name]f.e
*
          LENNOC = LENOCC (FNAME)
          CALL FTOVAX (FNAME, LENNOC)
          IF (LENNOC .LT. 1) THEN
            WRITE(*,*)
     &        '**IO_FILE** Error: file name =', FNAME, ' failed'
            GOTO 299
          ENDIF
*-----------------------------------------------------------
#endif
*
          LENNOC = LENOCC (FNAME)
          C_FNAME = FNAME(1:LENNOC)
          JTDKY = -1
          JTEVT = -2
          JTCOM =  LINP
          JTLIS =  6
          
          OPEN(-JTDKY,FILE=C_FNAME,STATUS='OLD')
*
          CALL ISAINI(JTDKY,JTEVT,JTCOM,JTLIS)
*
*     Read Isajet control cards; terminate for STOP command or error.
*
          IFL=0
          CALL ISABEG(IFL)
          IF (IFL .NE. 0) IEOTRI = 99
        ELSE                       ! skip upto stop
   30     READ (LINP, '(A)', END=299) CBUF
          IF (TAG(1:3) .EQ. 'END') THEN
*           End of block
            GOTO 299
          END IF
          GO TO 30
  299     CONTINUE
        END IF
      ELSEIF (KEYSTR .EQ. 'HERW') THEN        
*
* *** Initialize HERWIG
*
* --- MAX NUMBER OF EVENTS THIS RUN
        MAXEV=10
* --- BEAM PARTICLES
        PART1 = ' '
        PART2 = ' '
        CALL UHTOC (KBEAM , 4, PART1, 4)
        CALL UHTOC (KTARG , 4, PART2, 4)
* --- BEAM MOMENTA
        PBEAM1=PBEAM
        PBEAM2=PTARG
* --- PROCESS
        IPROC=KPROC
* --- INITIALISE OTHER COMMON BLOCKS
        CALL HWIGIN
* --- USER CAN RESET PARAMETERS AT
* --- THIS POINT, OTHERWISE DEFAULT
* --- VALUES IN HWIGIN WILL BE USED.
        PTMIN=PTMINS
* --- COMPUTE PARAMETER-DEPENDENT CONSTANTS
        CALL HWUINC
* --- CALL HWUSTA TO MAKE ANY PARTICLE STABLE
        CALL HWUSTA('PI0     ')
* --- USER'S INITIAL CALCULATIONS
        CALL HWABEG
* --- INITIALISE ELEMENTARY PROCESS
        CALL HWEINI
      ENDIF
#endif
*                                                        end of PY_INIT *
 999  END
