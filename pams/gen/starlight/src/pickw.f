c     This function picks a w for the 2-photon calculation.
      subroutine pickw(w)

      implicit NONE

      include 'inputp.inc'
      include 'Ftable.inc'
      include 'sig.inc'

      real sigofw(1000),sgfint(1000),ivalw,remainarea
      real remainw,w,sgf,signorm,x,a,b,c,ran
      integer i,j

c     DEAL WITH THE DELTA FUNCTION TYPE SIGMA CASE
      if (wdelta.ne.0.) then
        w=wdelta
        ivalw = ivalwd
        remainw = remainwd
      else

c     DEAL WITH THE CASE WHERE SIGMA IS AN ARRAY
c     sigofw is sigma integrated over y using a linear interpolation
c     sgfint is the integral of sgfint, normalized

c     integrate sigma down to a function of just w
        do 200 i = 1,numw
        sigofw(i) = 0
        do 100 j = 1,numy-1
          sigofw(i) =  sigofw(i) + (yarray(j+1)-yarray(j))  *
     &               (sigma(i,j+1)+sigma(i,j))/2.
 100    continue
 200    continue

c     calculate the unormalized sgfint array
        sgfint(1) = 0.
        do 300 i = 1,numw-1
          sgf = (sigofw(i+1)+sigofw(i))*(warray(i+1)-warray(i))/2.
          sgfint(i+1) = sgfint(i) + sgf
 300    continue

c     normalize the sgfint array
        signorm = sgfint(numw)
        do 400 i = 1,numw
          sgfint(i) = sgfint(i)/signorm
 400    continue

c     pick a random number
       x = ran(iseed)
c      call ranmar (x,1)

c     compare x and sgfint to find the ivalue which is just less than
c     the random number x
        do 500 i = 1,numw
          if(x.gt.sgfint(i)) ivalw = i
 500    continue

c     remainder above ivalw
        remainarea =x-sgfint(ivalw)

c     figure out what point corresponds to the excess area in
c     remainarea
        c = - remainarea * signorm / (warray(ivalw+1)-warray(ivalw))
        b = sigofw(ivalw)
        a = (sigofw(ivalw+1) - sigofw(ivalw) ) /2.
        if (a.eq.0) then
          remainw = - c/b
        else
          remainw =  (-b + sqrt(b**2 - 4.*a*c)) / (2. * a)
        endif

      ivalwd = ivalw
      remainwd = remainw

c     calculate the w value
        w = warray(ivalw) + (warray(ivalw+1)-warray(ivalw))*remainw
      endif

      return
      end

