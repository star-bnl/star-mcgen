c     This subroutine reads in a table of F values
c     parameterized by rapidity and mass

      subroutine readDiffLum

        include 'Ftable.inc'
        include 'inputp.inc'

c     numw is the number of w values
c     numy is the number of y values
c     warray is array of w values
c     y is the array of y values
c     f is the twoD array of f values

      open (unit=20,file='starlight.dat',status='unknown')

c     skip first eight entries, Z, A,Gamma,wmax,numw,ymax,numw,gg_or_gP
      read (20,*) dummy
      read (20,*) dummy
      read (20,*) dummy
      read (20,*) dummy
      read (20,*) dummy
      read (20,*) dummy
      read (20,*) dummy
      read (20,*) dummy

      do 100 i = 1,numw
        read(20,*) warray(i)
 100  continue
      do 200 i = 1,numy
        read(20,*) yarray(i)
 200  continue

      do 300 i = 1,numw
        do 400 j = 1,numy
          read(20,*) farray(i,j)
 400    continue
 300  continue

      close (unit=20)

      return
      end
