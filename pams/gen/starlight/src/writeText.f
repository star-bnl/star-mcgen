c     This subroutine writes two- and four-body decays to a file
c     that can be read by the old version of GEANT
      subroutine writeText (ievent,n,ipid,pxdec,pydec,pzdec)

      implicit NONE

      include 'inputp.inc'
      integer ievent,ipid,n,i,jtog
      real pxdec(4),pydec(4),pzdec(4),x,ran

c     randomly tag particles in each pair as particle/antiparticle

	do 100 i = 1,n,2
c          call ranmar(x,1)
	  x = ran(iseed)
          if (x.lt.0.5) then
            write (25,4001)  ievent,i,jtog(ipid), pxdec(i),pydec(i),
     &		pzdec(i)
            write (25,4001)  ievent,i+1,jtog(-ipid),pxdec(i+1),
     &		pydec(i+1),pzdec(i+1)
          else
            write (25,4001)  ievent,i,jtog(-ipid),pxdec(i),pydec(i),
     &		pzdec(i)
            write (25,4001)  ievent,i+1,jtog(ipid),pxdec(i+1),
     &		pydec(i+1),pzdec(i+1)
          endif
 100	continue

 4000 format(1X,I6,7X,I6)
 4001 format(1X,I6,1X,I6,1X,I2,1X,G12.5,1X,G12.5,1X,G12.5)

      return
      end
