C*******************************************************************************
        Subroutine HIJ_WRTXT( outlun, i_event, istat )
C*******************************************************************************
c
c       Description:-
c       =============
c	Writes ASCII text file output of HIJING event.  This format is not
c	standard (but is compatible with current GEANT text format) and 
c	might be replaced in the future. 
c
c
c       Author:  Matthew Bloomer
c       ======   
c
c       Last Modified:  17-Mar-1995 (written)
c       =============   23-Apr-1995 (Matt Bloomer): hjntog -> 
c                                   hijing_to_geant
c			13-May-1995 (M. Bloomer): status code 0 =success
c       Implicit inputs, outputs, side effects:-
c       ========================================
c	Assumes output file has already been open for formatted write.
c
c	Return Values:
c	=============
c	 0 = success
c	 1 = error writing event
c
c       Global Variables:-
c       ================
c
        Implicit None

c
c  Control parameters
c
      CHARACTER*4     FRAME, PROJ, TARG
      REAL            EFRM, BMIN, BMAX
      INTEGER         IAP, IZP, IAT, IZT, MAXEVT, N_EVENT_MSG
      INTEGER		      IFORMAT, IVERSION
      REAL		ETAMIN, ETAMAX

      COMMON/HI_CONTROL/ EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT,
     1                   BMIN, BMAX, MAXEVT, N_EVENT_MSG,
     1                   IFORMAT, IVERSION, ETAMIN, ETAMAX

c
c HIJING main common blocks:  put data type statements BEFORE COMMON
c statement, and include dimension.
c
        INTEGER  NATT, JATT, NT, NP, N0, N01, N10, N11
	INTEGER  KATT
        REAL     EATT, PATT

	COMMON /HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
	COMMON /HIMAIN2/KATT(130000,4),PATT(130000,4)

c
c HIJING parameter arrays
c
	INTEGER IHPR2, IHNT2
	REAL    HIPR1, HINT1

	COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
	
c
C
C       Local Specifications:-
C       ======================
c
        Integer         istat, i_event, size, outlun, luid, geid
	Integer	        j_particle, k, n_particles

c
c       Executable statements:-
c       ======================
c
c     Write out event in GEANT text file format.  First do
c     translation to GEANT particle id.  Translate and write
c     out ALL particles in KATT.  It is up to the user to decide
c     what subset of parents/daughters are to be placed in KATT.
c     Note that GEANT doesnt do particle filtering on TXT 
c     format files (i.e. it is up to you).
c
c     If hard scattering has been turned on by user, write out parton
c     momenta with rest of particle list.
c 
      If (ihpr2(3) .ge. 1 ) Then
	 n_particles = natt + 2
      Else
	 n_particles = natt
      Endif 

      Write( outlun, 250, Err=920) n_particles, 'event #', i_event,
     >         ' ID     Px(GeV)','  Py(GeV)','  Pz(GeV)'
 250  Format(1x,i6,a,i6,3x,a,5x,a,3x,a)

c
c  Write out particle info into text file.  Note that the parent
c  index is appended to particle momenta.
c
  
      Do 900 j_particle = 1, natt

         luid = katt(j_particle,1)
         Call HIJING_TO_GEANT(luid, geid)

         Write( outlun, 450, Err=920) geid, (patt(j_particle, k), k=1,3),
     >                                katt(j_particle,3) 
 450     Format(18x,i2,2x,G12.5,1x,G12.5,1x,G12.5,1x,i7)

 900  Continue

C
C  Add pre-radiation partons to particle list for triggered
C  jet production.  All partons are assumed to be up quarks
C  (I don't know what their real identity is).
C
      If (ihpr2(3) .ge. 1) Then

	luid = 2
        Call HIJING_TO_GEANT(luid, geid)
        Write( outlun, 450, Err=920) geid, hint1(21), hint1(22), 
     >       hint1(23), 0
	Write( outlun, 450, Err=920) geid, hint1(31), hint1(32), 
     >       hint1(33), 0

      EndIf

      istat = 0
      goto 999
 
 920  Write  ( 6, 930) i_event
 930  Format (1X, '%HIJ_WRTXT: Error while writing event ',I11)
      istat = 1

 999  Return
      End

