c     This subroutine reads in a table of F values
c     parameterized by rapidity and mass
c     comment added for newer modify date
C    when VM interference is present, it also reads a pt table

      subroutine readDiffLum

CVM: this is now in Ftable.inc        REAL f_max

        include 'Ftable.inc'
        include 'inputp.inc'
        include 'pttable.inc'
        DOUBLE PRECISION fpart
        DOUBLE PRECISION finterm(1000,1000)
c     numw is the number of w values
c     numy is the number of y values
c     warray is array of w values
c     y is the array of y values
c     f is the twoD array of f values

      f_max=0.0

      open (unit=20,file='starlight.dat',status='unknown')

c  Skip first 14 entries, 
C  Z, A,Gamma,Wmax,wmin,
C  numw,Ymax,numw,gg_or_gP,ibreakup
C  iinterfere,xinterfere,ptmax,NPT

      do 10 j=1,14
 10      read(20,*)dummy

         WRITE(6,11)
 11      FORMAT(' Reading 14 parameters from starlight.dat')


      do 100 i = 1,numw
        read(20,*) warray(i)
 100  continue
      do 200 i = 1,numy
        read(20,*) yarray(i)
 200  continue

      do 300 i = 1,numw
        do 400 j = 1,numy
          read(20,*) farray(i,j)
          IF( farray(i,j) .gt. f_max )f_max=farray(i,j)
 400    continue
 300  continue

C     Normalize farray (JN 010402)
      do 500 i = 1,numw
        do 600 j = 1,numy
          farray(i,j) = farray(i,j)/f_max
 600    continue
 500  continue

C  Is there a pt table (for vector mesons with interference)?

      if (gg_or_gP .eq. 1) goto 1000
      if (iinterfere .eq. 0) goto 1000

C  only numy/2 y bins here, from 0 (not -ymax) to ymax

      do 800 j=1,numy/2

         fmax=0.

C  we want to convert fptarray to an integral array where
C  fpt(j,1) is near 0, and fpt(j,NPT) ~ 1.  This will facilitate
C  a simple table lookup 


         fptsum=0.
         do 700 k=1,NPT
            read(20,*)fpart
            finterm(j,k)=fpart
            fptarray(j,k)=0.
            fptsum=fptsum+fpart
 700     continue

C  convert array to integral

         fptarray(j,1)=finterm(j,1)/fptsum

         do 750 k=2,NPT
            do 730 kk=1,k
               fptarray(j,k)=fptarray(j,k)+finterm(j,kk)
 730        continue
            fptarray(j,k)=fptarray(j,k)/fptsum

C            write(6,740)j,k,fptarray(j,k),fmax
C 740        format(' ',I3,2X,I3,E11.3,2X,E11.3)

 750        continue

 800  continue

 1000 close (unit=20)
      return
      end




