C*******************************************************************************
      Program HIJ_TEST
C*******************************************************************************
C
C     Description:
C
C     Very simple shell to read in HIJING events from a text file and sort
C     results into reference histograms to be stored along with the
C     package.
C
C     Author:         Matthew Bloomer
C     Created:        2 April, 1994
C     Last modified:
C
C
      Implicit None

      Integer  istat, ierr, maxp, inlun
      Parameter (maxp=25000, inlun=22)

      Integer  HIJ_TEST_BOOK, HIJ_TEST_SORT, HIJ_TEST_END
      Integer  HIJ_RDTXT
      Integer  n_event, i_event, j_event, max_events

      Integer  n_particles, geid(maxp), iparent(maxp)
      Real     px(maxp), py(maxp), pz(maxp)

      Character*132  infile

C*******************************************************************************

C
C Prompt user for number of events (necessary for fifo running)
C
 40   Write( *, 50 )
 50   Format( 1X, 'Enter number of events to analyse: ',$)
      Read( *, *, err=40 ) max_events
      
C
C Define histograms
C
      istat = HIJ_TEST_BOOK()
      If ( istat .ne. -1 ) Then
         Write (*, 110 ) istat
 110     Format (1X,'Error in HIJ_TEST_BOOK, status = ', i10)
         Goto 990
      EndIf

C
C Open file, which could be a fifo file (hence default name of
C /tmp/hijing.fifo).  Use filename from Z_TAPE environment variable 
C if set.
C 
      Call getenv('Z_TAPE',infile)
      If (infile .eq. '') infile = '/tmp/hijing.fifo'

      Open( inlun, file=infile, status='OLD', access='SEQUENTIAL',
     >     form='FORMATTED', err=910, iostat=ierr )

C
C Initialize event counter.  Loop until max_events or until
C EOF reached in HIJ_RDTXT.
C
      n_event = 0

      Do j_event = 1, max_events
C 
C Read one event
C

 300     istat = HIJ_RDTXT( inlun, i_event, n_particles, geid, 
     >        px, py, pz, iparent )

         If (istat .eq. 1 ) Then
C
C End of file, go on to edit histograms
C
            Goto 400

         ElseIf (istat .ne. -1 ) Then
C
C Error reading file
C
            Write (*, 310 ) istat
 310        Format (1X,'Error in HIJ_RDTXT, status = ', i10)
            Goto 990
         
         Else
C
C  Successful read, sort histograms
C
            istat = HIJ_TEST_SORT( i_event, n_particles, geid, 
     >           px, py, pz, iparent )
            If (istat .ne. -1 ) Then
               Write (*, 340 ) istat
 340           Format (1X,'Error in HIJ_TEST_SORT, status = ', i10)
               Goto 990
            EndIf

            n_event = n_event + 1 

         EndIf
      EndDo

C
C End of file reached.  Edit histograms, store and exit
C

 400  Close (inlun)

      istat = HIJ_TEST_END(n_event)
      If (istat .ne. -1 ) Then
            Write (*, 440 ) istat
 440        Format (1X,'Error in HIJ_TEST_END, status = ', i10)
            Goto 990
      EndIf
      Goto 990

 910  Write (*, 911) infile, ierr 
 911  Format (1X,'Error reading file ',A,/,'ierr = ',i10)

 990  STOP
      End


