c     This subroutine takes the event stored in the jetset
c     common block and writes it to a file that can be
c     read by the old version of GEANT

      subroutine writejetsetText(ievent)

      implicit NONE

      include 'lujets.inc'
      integer i,ipart,ievent,jtog

c     figure out how many final state particles exist in an event
      ipart = 0
      do 100 i = 1,n
c     1 means that this particle is considered stable by JETSET
c     4 means that JETSET could have decayed the particle but didn't
        if((k(i,1).eq.1).or.(k(i,1).eq.4)) ipart = ipart + 1
 100  continue

C     write the event to the file
      write (25,4000) ipart,ievent

      do 200 i = 1,n
         if((k(i,1).eq.1).or.(k(i,1).eq.4)) write (25,4001)
     &   jtog(k(i,2)),p(i,1),p(i,2),p(i,3)
  200  continue

 4000 format(1X,I6,7X,I6)
 4001 format(18X,I2,2X,G12.5,1X,G12.5,1X,G12.5)

      return
      end
