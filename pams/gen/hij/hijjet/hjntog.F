C=======================================================================
      SUBROUTINE HJNTOG(LUID,GEID)
C=======================================================================
C
C     Input Parameters:
C     =================
C
C     Output Parameters:
C     ==================
C
C     Calling:
C     ========
C
C     Functional Description :
C     ========================
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Subroutine to convert HIJING (PDG) id's to GEANT ID's         *
C.    *                                                                *
C.    ******************************************************************
C.
C     Author:-
C     ========
C     Richard Morse  LBL (510)-486-6349
C
C     Creation Date:
C     Modified     : 
C	3 April 1995 Matt Bloomer - extended GLTBLE to include Delta
C		     resonances.  STAR GEANT will ignore extensions.
C
C     Global specifications:
C     ======================
C
C     External Specifications:
C     ========================
      INTEGER IUCOMP
C
C     Local Specifications:
C     =====================
      INTEGER LUID, GEID, SAVID
C
      INTEGER NUMID
      PARAMETER (NUMID=75)
      INTEGER GLTBLE(1:75)
      SAVE GLTBLE
C
C     Data statements:
C     ================
      DATA GLTBLE/   22,  -11,   11,   12,  -13
     +,              13,  111,  211, -211,  130
     +,             321, -321, 2112, 2212,-2212
     +,             310,  221, 3122, 3222, 3212
     +,            3112, 3322, 3312, 3334,-2112
     +,           -3122,-3112,-3212,-3222,-3322
     +,           -3312,-3334,  -15,   15,  411
     +,            -411,  421, -421,  431, -431
     +,            4122,   24,  -24,   23,    0
     +,               0,    0,    0,  223,  333
     +,             113,  213, -213,  411,  421
     +,            1114, 2114, 2214, 2224,    0
     +,               1,   -1,    2,   -2,    3
     +,              -3,    4,   -4,    5,   -5
     +,               6,   -6,    21, -21,    0/
C
C     Executable Statements:
C=======================================================================
      SAVID=IUCOMP(LUID,GLTBLE,75)
c
c     We must map all neutrinos onto a single ID...
c
      IF(ABS(LUID).EQ.12 .OR.
     +   ABS(LUID).EQ.14 .OR.
     +   ABS(LUID).EQ.16)THEN
         GEID=4
      ELSE
         GEID = SAVID
      ENDIF
c
c     The rest become GEANTINOS...
c
      IF(GEID.EQ.0)GEID=48

999   Continue
      Return
      End
