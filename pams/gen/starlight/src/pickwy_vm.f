c     This subroutine picks a w and y value for the vector-meson
c     calculation.
      subroutine pickwy_vm(W,Y)

      IMPLICIT NONE

      include 'D2LParam.inc'
      include 'global.inc'
      include 'Ftable.inc'
      include 'range.inc'
      include 'inputp.inc'
      real W,Y,IW,IY,xw,xy,xtest,ran,dW,dY      
      integer ISEED

      dW = (Wtop-Wmin)/DFLOAT(numw)
      dY = (Ytop-Ymin)/DFLOAT(numy)

C       >> DRAW xw,xy
 201    xw = ran(ISEED)
        W = Wmin + xw*(Wtop-Wmin)
        IW = INT((W-Wmin)/dW) + 1
        xy = ran(ISEED)
        Y = Ymin + xy*(Ytop-Ymin)
        IY = INT((Y-Ymin)/dY) + 1

C       >> Check if this matches the map
        xtest = ran(ISEED)
        IF(xtest.gt.farray(IW,IY))GOTO 201

C       >> W,Y should now have the right distribution

 	return
        END
