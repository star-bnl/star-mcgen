c     This subroutine gets the input parameters from starlight.in
c     See starlight.doc for an explanation of input parameters

      subroutine input

      implicit NONE

      include 'inputp.inc'
      include 'D2LParam.inc'
      include 'pttable.inc'

      open (unit=15,file='starlight.in',status='unknown')
      read (15,*) Z
      read (15,*) A
      read (15,*) gamma_em
      read (15,*) Wmax
      read (15,*) Wmin
      read (15,*) numw
      read (15,*) Ymax
      read (15,*) numy
      read (15,*) gg_or_gP
      read (15,*) ievents
      read (15,*) ip
      read (15,*) iseed
      read (15,*) iout

C added 4/12/2001 by SRK ibreakup: 1=don't care; 2= mutual excitation
C 3= no excitation 4=exactly 1 excitation 5=1+ excitation

      read (15,*) ibreakup

C  added 8/2002 by SRK for interference

      read(15,*) iinterfere
      read(15,*) xinterfere
      read(15,*) ptmax
      read(15,*) NPT
      dpt = ptmax/DFLOAT(NPT)
      close (unit=15)

      return
      end
