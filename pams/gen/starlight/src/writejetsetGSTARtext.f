c     This subroutine takes the event stored in the jetset common
c     block and writes it to a readable gstar file
c     for a description of the various parameters consult the
c     documentation of gstar and jetset.
      subroutine writejetsetGSTARtext(ievent)

      implicit NONE
        
      include 'lujets.inc'
C     variables for vertex
      real x(500),y(500),zz(500),t(500)
      integer iparent_track(500),ndaughter(500)
c     variables for track
      integer start_vertex(4000),stop_vertex(4000)
      integer i,n_tracks,n_vertices,iver,iorg,ievent,jtog

      do 100 i = 1,n
         start_vertex(i)=1
         stop_vertex(i) = 0
 100   continue

      n_tracks = n
      n_vertices = 1
      iver = 1
      x(iver) = 0.
      y(iver) = 0.
      zz(iver) = 0.
      t(iver) = 0.
      iparent_track(iver) = 0
      ndaughter(iver) = 0
      iorg = 0

      do 200 i = 1,n_tracks
         if(k(i,3).eq.iorg) then
           ndaughter(iver) = ndaughter(iver) + 1
         else
           iver = iver + 1
           stop_vertex(k(i,3)) = iver
           ndaughter(iver) = 1
           iorg = k(i,3)
c     jetset works in mm and GEANT in cm
           x(iver) = v(i,1)/10.
           y(iver) = v(i,2)/10.
           zz(iver) = v(i,3)/10.
c     jetset works in mm/c and GEANT in sec
           t(iver) = v(i,4)*3.336*10**(-12)
           iparent_track(iver) = iorg
           n_vertices = n_vertices + 1
         endif

         start_vertex(i) = iver
 200   continue

      write (25,5000) ievent,n_tracks,n_vertices

      do 300 i = 1,n_vertices
        write (25,5002) x(i),y(i),zz(i),t(i),i,0,
     &        iparent_track(i),ndaughter(i)
 300  continue

      do 400 i = 1, n_tracks
         write (25,5001) jtog(k(i,2)),p(i,1),p(i,2),p(i,3),i,
     &         start_vertex(i),stop_vertex(i),k(i,2)
 400  continue

 5000  format ('EVENT:',3x,3(1x,i6))
 5001  format ('TRACK:',1x,i6,3(1x,g12.5),4(1x,i6))
 5002  format ('VERTEX:',4(1x,g12.5),4(1x,i6))

      return
      end
