C*******************************************************************************
        Subroutine HIJ_READ_PARAMETERS( istat )
C*******************************************************************************
c
c       Description:-
c       =============
c           Reads user-specified input file of control parameters for running
c       HIJING.
c
c       Author:  Chuck Naudet
c       ======   Matthew Bloomer
c
c       Last Modified:  1-Feb-1991
c       =============
c	17-Mar-1995 Matthew Bloomer : modify format of parameter 
c		    statement.  Etamin and Etamax discontinued.
c
c       Output parameters:-
c       ==================
c       istat       (I*4) : Execution status code  -1 = succesful
c                                                   0 = error
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
	INTEGER		IFORMAT, IVERSION
	REAL		ETAMIN, ETAMAX

        COMMON/HI_CONTROL/ EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT,
     1                     BMIN, BMAX, MAXEVT, N_EVENT_MSG,
     1			   IFORMAT, IVERSION, ETAMIN, ETAMAX

c
c iformat:  output format (1=egz, 2=text)
c

	DATA 		IFORMAT, IVERSION / 1, 2 /


        CHARACTER*72    HEADER(100)
        INTEGER         NLINE
        COMMON/TEXT/ HEADER, NLINE

c
c  Event options and parameters (internal to HIJING).
c
        COMMON/HIPARNT/HIPR1, IHPR2, HINT1, IHNT2

        INTEGER         IHPR2(50), IHNT2(50)
        REAL            HIPR1(100), HINT1(100)

c  Hallman improvement:  seed values entered by user

	INTEGER		NSEED
        COMMON/RANSEED/NSEED

C
C       Local Specifications:-
C       ======================
C
        Character*120   infile
        Character*72    line
        Integer         istat, inlun, index, iline
	Real		value
	Character*5	array

        Data inlun/20/

c
c       Executable statements:-
c       ======================
c

        Write( *, 200 )
 200    Format( 'Enter filename with input parameters (C*120)',/,
     1          'File: ',$)
        Read( *, 220 ) infile
 220    Format( A )
        Open( inlun, file=infile, status='OLD', form='FORMATTED',
     1        access='SEQUENTIAL' )

c
c  Read text header of input file.  First 100 lines of text comment will
c  be outputted as header of file.
c
        nline = 0
 300    Read( inlun, 310, end=900 ) line
 310    Format( A )
        If( nline .ge. 100 )Then
            nline = 100
            header(nline) = '$$$$'
            goto 400
        Endif
        nline = nline + 1
        header(nline) = line
        If( line(1:4) .NE. '$$$$' )Goto 300

c
c  Read in control parameters: etamin and etamax are obsoleted
c
 400    Read( inlun, *, end=900 )   efrm
        Read( inlun, *, end=900 )   frame
	Read( inlun, *, end=900 )   proj, targ
        Read( inlun, *, end=900 )   iap, izp, iat, izt
        Read( inlun, *, end=900 )   maxevt
        Read( inlun, *, end=900 )   bmin, bmax
c	Read( inlun, *, end=900 )   etamin, etamax
	Read( inlun, *, end=900 )   nseed
        Read( inlun, *, end=900 )   n_event_msg
        Read( inlun, *, end=900 )   iformat
c
c  Force user to specify commonly used parameters (Yepes-inspired)
c
c  quench parameter (0=no quench, 1=quench)
c
        Read( inlun, *, end=900 )   ihpr2(4)
c
c  Decays: keep parents with daughters in particle list (=1)
c          or stable daughters only (=0)
c
        Read( inlun, *, end=900 )   ihpr2(21)
c
c  Sample jets:  one hard jet in each event (=1), else min-bias (=0)
c                pt of sample jet (<0 means jets of pt > abs(value))
c
        Read( inlun, *, end=900 )   ihpr2(3), hipr1(10)

c
c  Any remaining lines are for changes in the parameter arrays of
c  common/HIPARNT/

 500	Read( inlun, *, end=599) array, index, value

	If( 	array .eq. 'HIPR1' .and.
     +   	index .ge. 1 .and. index .le. 100 )then
		hipr1(index) = value
	Elseif( array .eq. 'IHPR2'  .and.
     +	        index .ge. 1 .and. index .le. 50 )then
		ihpr2(index) = value
	Elseif( array .eq. 'HINT1' .and.
     +   	index .ge. 1 .and. index .le. 100 )then
		hint1(index) = value
	Elseif( array .eq. 'IHNT2'  .and.
     +  	index .ge. 1 .and. index .le. 50 )then
		ihnt2(index) = value
	Else
		Write(*, 540) array, index, value
 540		Format( 'Error reading line: ',A,2x,'index= ',i8,2x,
     +		'value= ',g12.5)
	Endif
	Goto 500

c
c  Echo values out back to user
c

 599	Continue
  	Write( *, 615 )
 	Do iline = 1, nline-1
	   Write( *, 610 ) header(iline)
 610	   Format( A )
	Enddo

	Write( *, 620 )  efrm, frame
	Write( *, 630 )  proj, iap, izp
	Write( *, 640 )  targ, iat, izt
	Write( *, 650 )  maxevt, bmin, bmax
c	Write( *, 660 )  etamin, etamax
	Write( *, 670 )  nseed
        Write( *, 680 )  n_event_msg, iformat
        Write( *, 690 )  ihpr2(4)
        Write( *, 700 )  ihpr2(21)
        Write( *, 710 )  ihpr2(3), hipr1(10)

#ifdef AIX
C
C       AIX uses the C-Library random number generator...
C	call SRAND to initialize RAND...
C
        if(nseed .ne. 0)then
           call SRAND(NSEED)
        endif
#endif

 615	Format(/,'Summary of parameter file: ',/)
 620	Format('Energy of collision = ',g12.5,' GeV in ',A,' frame')
 630	Format('Projectile : ',A,3x,'(A= ',i4,', Z=',i4,')')
 640	Format('Target     : ',A,3x,'(A= ',i4,', Z=',i4,')')
 650	Format('Nevents to be generated = ',i10,' between ',
     +  	 f6.2,' fm and ',f6.2,' fm')
 660	Format('Eta limits for outputted particles : ',g12.5,' to ',
     +		 g12.5 )
 670	Format('Seed value of random number generator : ',i15)
 680    Format('Event message frequency:  every ', i7, ' events',/,
     +         'Output file format (1=egz, 2=text):', i5)
 690    Format('Options: quenching     -> ihpr2( 4) = ',i5)
 700    Format('Options: decay history -> ihpr2(21) = ',i5)
 710    Format('Options: sample jets   -> ihpr2( 3) = ',i5,/,
     +         'Options: jet pt        -> hipr1(10) = ',g12.5)
 
c
c  Sucessful read of file.  Prompt user for frequency of event messages
c

 799    Close (inlun)

        istat = -1
        RETURN

 900    istat = 0
        RETURN

        END

#ifdef AIX
        REAL FUNCTION RAN(ISEED)
C
C       AIX uses the C-Library random number generator...
C
        INTEGER ISEED
        RAN = RAND()
        RETURN
        END
#endif



