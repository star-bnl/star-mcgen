C
C
C*******************************************************************************
        PROGRAM HIJEVT
C*******************************************************************************
C
C	********************************************************************
C	*			HIJING Monte Carlo Program		   *
C	*  		For Relativisitic Heavy Ion Collision		   *
C	*								   *
C	*		  Miklos Gyulassy and Xin-Nian Wang		   *
C	*		Nuclear Science Division, MS 70A-3307		   *
C	*	Lawrence Berkeley Laboratory and University of California  *
C	*	    		Berkeley, CA 94720 USA			   *
C	*		       Bitnet: XNWANG@LBL.bitnet		   *
C	********************************************************************
C
C*******************************************************************************
C
c       Description:-
c       =============
c           General routine to run Hijing and write events to disk for future
c       analysis.  The controlling parameters are read from a user-specified
c       file in HIJ_READ_PARAMETERS.  The output file header is written in
c       HIJ_WRITE_HEADER, and event by event data are written to file in
c       HIJ_WRITE_EVENT.
c
c
c       Author:  Chuck Naudet (a long time ago)
c       ======   Matthew Bloomer
c		 Richard Morse
c
c       Last Modified:  1-Feb-1991
c       =============
c       14-MAY-1991 R. Morse MIT (added BNL RHIC MC interface)
c       24-JAN-1992 R. Morse LBL (updated to HIJING 2.0 from CMS
c	16-Mar-1995 M.A. Bloomer (reorganization of hijing in library and 
c                     mods needed for new OS's).
c       24-Apr-1995 M.A. Bloomer (new EGZ implemented:  
c                     - different include files
c                     - changes to egzinit and egzout
c                     - isgeant variable defined
c		      - different return codes for egz routines: 0 = success)
c       13-May-1995 M.A. Bloomer (prompt user for file name, Z_TAPE disappears)
c
c       Implicit inputs, outputs, side effects:-
c       ========================================
c
c
c       Global Variables:-
c       ================
c
        Implicit None
        Integer         INDEX

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

c->rm
c     
c     include PAWC common in order to initialise KUIP...
c     
      Integer hmemor
      common/pawc/hmemor(40000)

#include "evt_zebra.inc"
#include "evt_ctrl.inc"

      integer*4 egzinit
      integer*4 egzout_hijing
      integer*4 egzend
c->rm
	
c
C
C       Local Specifications:-
C       ======================
c
        Character*80    outfile
        Character*50    routine
        Integer         istat, i_event, outlun, size, isgeant, length

        data evtcod/hijng/
        data outlun/22/
	data isgeant/0/

c
c       Executable statements:-
c       ======================
c

      Write( *, 100)
 100  Format('%HIJEVT:  Main program to generate HIJING data files',//)

c
c   Read in control parameters
c
      CALL HIJ_READ_PARAMETERS( istat )
      If( istat .NE. -1 ) goto 7000

c ->rm
c     Initialise KUIP in order to use kupryn...
c
      Call HLIMIT(40000)
      Call KUINIT(10000)
c ->rm

c Prompt for file name 

      Call KUPROS( 'Enter filename for output', outfile, length )
      If ( length .eq. 0 ) Goto 7500

c
c  Initialization for output:
c     text: Open input file for text output format 
c     egz:  startup egz
c
c      call getenv('Z_TAPE',outfile)

      If ( iformat .eq. 1 ) Then

c
c     Initialize the ZEBRA store
c     ==========================
c
         call mzebra(-3)
         z_input = .false.
c
c init egz package, isgeant=0 for non-GEANT applications
c MAB note that egzinit returns 0 for sucess (GDMIT!)
c
         istat = egzinit( isgeant, outfile )
         if( istat .ne. 0 ) goto 8000

      ElseIf ( iformat .eq. 2 ) Then

           Open( outlun, file=outfile, status='UNKNOWN', 
     1        access='SEQUENTIAL', form='FORMATTED' )

      Else

         goto 8100

      Endif

c
c  Setup HIJING for running
c
      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)

c
c  START of EVENT Loop
c
      Do 1000 i_event = 1, maxevt

         CALL HIJING( FRAME, BMIN, BMAX) ! get a HIJING event

         If (iformat .eq. 1 ) Then

            istat = egzout_hijing()                 ! write out ZEBRA event
            If( istat .ne. 0 ) goto  9000

         ElseIf (iformat .eq. 2 ) Then

            CALL HIJ_WRTXT( outlun, i_event, istat) ! write out TEXT event
            If( istat .ne. 0 ) goto  9100

         Endif

         If( mod( i_event, n_event_msg) .eq. 0 )Then
             Write (*, 400 ) i_event
 400         Format( '%HIJEVT: Number of events generated = ',i10)
         Endif

 1000 CONTINUE

      If ( iformat .eq. 1 ) Then
c
c     Finish up Zebra
c     ===============
         istat = egzend(isgeant)

      ElseIf ( iformat .eq. 2 ) Then

         CLOSE ( outlun )
       
      Endif
      STOP

c
c  Error messages
c
 7000   size = 19
        routine(1:size) = 'HIJ_READ_PARAMETERS'
        Write( *, 9900 ) routine(1:size), istat
        STOP

 7500   Write( *, 7600 ) 
 7600   Format( '%HIJEVT:  Invalid file name for output file' )
        STOP

 8000   size = 16
        routine(1:size) = 'egzinit'
        Write( *, 9900 ) routine(1:size), istat
        STOP

 8100   Write( *, 8150 ) iformat
 8150   Format( '%HIJEVT: format type ',i5,' for output file not understood')
        STOP

 9000   size = 6
        routine(1:size) = 'egzout'
        Write( *, 9900 ) routine(1:size), istat
        STOP

 9100   size = 10
        routine(1:size) = 'HIJ_WRTXT'
        Write( *, 9900 ) routine(1:size), istat
        STOP

 9900   Format( '%HIJEVT:  Status error after ',a,', istat=', i6 )

	END
C  DEC/CMS REPLACEMENT HISTORY, Element HIJEVT.FOR
C  *4    23-JUL-1991 12:06:56 BLOOMER "contains wrong version number"
C  *3    23-JUL-1991 12:03:58 BLOOMER "cosmetic changes"
C  *2     5-APR-1991 12:22:35 BLOOMER "Change HI_CONTROL common block"
C  *1     5-FEB-1991 14:26:09 BLOOMER "Main program to generate disk files of HI
C  DEC/CMS REPLACEMENT HISTORY, Element HIJEVT.FOR
C  DEC/CMS REPLACEMENT HISTORY, Element HIJEVT.FOR
C  *1    24-JAN-1992 14:53:21 RMORSE "EVZ:- HIJING 2.0 with ZEBRA I/O MAIN"
C  DEC/CMS REPLACEMENT HISTORY, Element HIJEVT.FOR

