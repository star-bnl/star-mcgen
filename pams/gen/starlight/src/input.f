c     This subroutine gets the input parameters from starlight.in
c     See starlight.doc for an explanation of input parameters

      subroutine input

      implicit NONE

      include 'inputp.inc'
      include 'D2LParam.inc'

      open (unit=15,file='starlight.in',status='unknown')
      read (15,*) Z
      read (15,*) A
      read (15,*) gamma_em
      read (15,*) wmax
      read (15,*) numw
      read (15,*) ymax
      read (15,*) numy
      read (15,*) gg_or_gP
      read (15,*) ievents
      read (15,*) ip
      read (15,*) iseed
      read (15,*) iout
      close (unit=15)

      return
      end
