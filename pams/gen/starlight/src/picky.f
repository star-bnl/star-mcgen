c     given a w, this function picks a y
      subroutine picky (y)

      implicit NONE

      include 'inputp.inc'
      include 'Ftable.inc'
      include 'sig.inc'

      real sigofy(1000),sgfint(1000),remainw,x,remainarea,remainy
      real a,b,c,y,sgf,signorm,ran
      integer ivalw,j,ivaly,i

      ivalw = ivalwd
      remainw = remainwd

c     avarage over two columns to get y array
      do 100 j = 1,numy
        sigofy(j) = sigma(ivalw,j) +
     &           (sigma(ivalw+1,j)-sigma(ivalw,j)) * remainw
 100  continue

c     calculate the unormalized sgfint
      do 150 j = 1,numy-1
 150  continue

      sgfint(1) = 0.
      do 200 j = 1,numy-1
        sgf = (sigofy(j+1)+sigofy(j))/2.
C        print *,yarray(j),sgf
        sgfint(j+1) = sgfint(j) + sgf*(yarray(j+1)-yarray(j))
 200  continue

c     normalize the sgfint array
      signorm = sgfint(numy)
      do 300 j = 1,numy
        sgfint(j) = sgfint(j)/signorm
 300  continue

c     pick a random number
      x = ran(iseed)
c      call ranmar(x,1)

c     compare x and sgfint to find the ivalue which is just less than
c     the random number x
      do 500 i = 1,numy
        if(x.gt.sgfint(i)) ivaly = i
 500  continue

c     remainder above ivaly
      remainarea =x-sgfint(ivaly)

c     figure what point corresponds to the leftover area in
c     remainarea
      c = - remainarea * signorm / (yarray(ivaly+1)-yarray(ivaly))
      b = sigofy(ivaly)
      a = (sigofy(ivaly+1) - sigofy(ivaly) ) /2.
      if (a.eq.0) then
        remainy = - c/b
      else
        remainy =  (-b + sqrt(b**2 - 4.*a*c)) / (2. * a)
      endif

c     calculate the y value
      y = yarray(ivaly) + (yarray(ivaly+1)-yarray(ivaly))*remainy

      return
      end

