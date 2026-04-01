!c $Id: boxinc.f,v 1.3 1999/01/18 09:56:53 ernst Exp $
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     Unit     : all modules using the box
!c     Version  : 1.0
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c       max count of different species
        integer:: bptmax
        parameter(bptmax=20)

!c      counter
        integer,target:: cbox
!c      Flags
        integer,target:: boxflag
        integer,target:: edensflag, mbflag
        integer,target:: para,solid,mtest
!c      number of different species
        integer,target:: mbox
!c      particle counters
        integer,target:: bptpart(bptmax),bptityp(bptmax),bptiso3(bptmax)
        real*8,target:: bptpmax(bptmax)
        real*8,target:: edens
!c      edge length of a cube in fm
        real*8,target::  lbox
!c      half edge lengt of a cube 
        real*8,target::  lboxhalbe
!c      double edge length of a cube
        real*8,target::  lboxd
!c      momenta
        real*8,target:: mbp0, mbpx, mbpy, mbpz
        
        common /boxic/ cbox,boxflag, mbox, bptityp, bptiso3, bptpart
        common /boxic/ edensflag, para, solid, mbflag, mtest
        common /boxrc/ lbox, lboxhalbe, lboxd, bptpmax, edens
        common /boxrc/ mbp0, mbpx, mbpy, mbpz
        
