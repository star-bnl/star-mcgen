c     This subroutine writes two or four-particle decays to a text file
c     that can be read by GSTAR

      subroutine writeGSTARtext (ievent,n,ipid,px,py,pz)

      implicit NONE

      include 'inputp.inc'
      integer ievent, n, ipid,i,jtog
      real px(4),py(4),pz(4),x,ran

c     write header line
      write (25,5000) ievent,n,1
      write (25,5002) 0.,0.,0.,0.,1,0,0,n

c     randomly assign particle/antiparticle designations for each pair 
c     and write out info
      do 100 i = 1,n,2
c      call ranmar(x,1)
      x = ran(iseed)
      if (x.lt.0.5) then
        write (25,5001)  jtog(ipid), px(i),py(i),pz(i),i,1,0,ipid
        write (25,5001)  jtog(-ipid),px(i+1),py(i+1),pz(i+1),i+1,
     &		1,0,-ipid
      else
        write (25,5001)  jtog(-ipid), px(i),py(i),pz(i),i,1,0,-ipid
       write (25,5001)  jtog(ipid),px(i+1),py(i+1),pz(i+1),(i+1),
     &		1,0,ipid
      endif
 100  continue


 5000 format ('EVENT:',3x,3(1x,i6))
 5001 format ('TRACK:',1x,i6,3(1x,g12.5),4(1x,i6))
 5002 format ('VERTEX:',4(1x,g12.5),4(1x,i6))

      return
      end


