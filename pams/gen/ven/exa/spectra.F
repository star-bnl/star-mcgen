c-*-Fortran-*-
C
C File: spectra.F -- fills dN/dy and 1/pt dN/dpt histograms
C                    for specified particles.
C
C Created:  9 May 1995 Matthew Bloomer
C Based on:  hij/src/hij_test_routines.F
C
C This file contains 3 subroutines:
C
C aninit
C anuser
C anfin
C
C These routines are called by read_zebra_evt.F which can be used
C as an EGZ file analysis shell.  This file provides a (simple) 
C example of performing analysis of EGZ event files generated 
C by any event generator. 
C
C*******************************************************************************
      Subroutine aninit
C*******************************************************************************
C
C     Description:
C
C     Books histograms.
C
C     Author:         Matthew Bloomer
C     Created:        9 May, 1995
C     Last modified:
C
C*******************************************************************************

      Implicit None

      Integer   nwpaw, hmemor
      Parameter (nwpaw=200000)

      common/pawc/hmemor(nwpaw)

C  Counting number of events analysis

      integer n_event
      common/spectra/n_event

C*******************************************************************************

C
C  dN/deta charged = stable charged particles
C
C  Initialize Hbook (I don't know if the mzebra call is really necessary, but
C  the - sign in the hlimit call sure is...)
C
      Call mzebra(-3)
      Call hlimit(-nwpaw)

      Call HBOOK1( 1,'Geant ID',        70, 0., 70., 0. )
      Call HBOOK1( 2,'dN/dy pi+',      140, -7., 7., 0. )
      Call HBOOK1( 3,'dN/dy pi-',      140, -7., 7., 0. )
      Call HBOOK1( 4,'dN/dy pi0',      140, -7., 7., 0. )
      Call HBOOK1( 5,'dN/dy K+',       140, -7., 7., 0. )
      Call HBOOK1( 6,'dN/dy K-',       140, -7., 7., 0. )
      Call HBOOK1( 7,'dN/dy K0',       140, -7., 7., 0. )
      Call HBOOK1( 8,'dN/dy quarks',   140, -7., 7., 0. )
      Call HBOOK1( 9,'dN/dy gluons',   140, -7., 7., 0. )
      Call HBOOK1(10,'dN/eta charged', 140, -7., 7., 0. )

      Call HBOOK1(20,'dN/dy p',        140, -7., 7., 0. )
      Call HBOOK1(21,'dN/dy p-bar',    140, -7., 7., 0. )
      Call HBOOK1(22,'dN/dy Lambda',   140, -7., 7., 0. )
      Call HBOOK1(23,'dN/dy Lambda-bar', 140, -7., 7., 0. )
      Call HBOOK1(24,'dN/dy Sigma+',   140, -7., 7., 0. )
      Call HBOOK1(25,'dN/dy Sigma-',   140, -7., 7., 0. )
      Call HBOOK1(26,'dN/dy Sigma0',   140, -7., 7., 0. )

      Call HBOOK1(40,'1/pt dN/dpt pi+ (ycm)',    50, 0., 5., 0. )
      Call HBOOK1(41,'1/pt dN/dpt pi- (ycm)',    50, 0., 5., 0. )
      Call HBOOK1(42,'1/pt dN/dpt pi0 (ycm)',    50, 0., 5., 0. )
      Call HBOOK1(43,'1/pt dN/dpt K+ (ycm)',     50, 0., 5., 0. )
      Call HBOOK1(44,'1/pt dN/dpt K- (ycm)',     50, 0., 5., 0. )
      Call HBOOK1(45,'1/pt dN/dpt K0  (ycm)',    50, 0., 5., 0. )
      Call HBOOK1(46,'1/pt dN/dpt quarks',       50, 0., 5., 0. )
      Call HBOOK1(47,'1/pt dN/dpt gluons',       50, 0., 5., 0. )

      Call HBOOK1(50,'1/pt dN/dpt p (ycm)',      50, 0., 5., 0. )
      Call HBOOK1(51,'1/pt dN/dpt p-bar (ycm)',  50, 0., 5., 0. )
      Call HBOOK1(52,'1/pt dN/dpt Lambda (ycm)', 50, 0., 5., 0. )
      Call HBOOK1(53,'1/pt dN/dpt Lambda-bar(ycm)', 50, 0., 5., 0. )
      Call HBOOK1(54,'1/pt dN/dpt Sigma+ (ycm)', 50, 0., 5., 0. )
      Call HBOOK1(55,'1/pt dN/dpt Sigma- (ycm)', 50, 0., 5., 0. )
      Call HBOOK1(56,'1/pt dN/dpt Sigma0 (ycm)', 50, 0., 5., 0. )

      n_event = 0

      Return
      End

C*******************************************************************************
      Subroutine anuser(istat)
C*******************************************************************************
C
C     Description:
C
C     Sorts histograms.
C
C     Author:         Matthew Bloomer
C     Created:        2 April, 1994
C     Last modified:
C
C     Return values:
C         0  Success
C        -1  Not Success (never happens, of course)
C*******************************************************************************

      Implicit None

#include "event.inc"

C  Counting number of events analysis

      integer n_event
      common/spectra/n_event


      Integer  j_part, gid, istat
      Real     pperp, eta, rap

c     
c List of assumed GEANT ids (taken from hijing_to_geant.F) is
c same as for GEANT except for extensions described below:
c     
c     1       : gamma
c     2       : positron
c     3       : electron
c     4       : neutrino
c     5       : Muon+
c     6       : Muon-
c     7       : Pion 0
c     8       : Pion +
c     9       : Pion -
c     10      : Kaon 0 long
c     11      : kaon +
c     12      : kaon -
c     13      : neutron
c     14      : proton
c     15      : Antiproton
c     16      : Kaon 0 short
c     17      : Eta
c     18      : lambda
c     25      : Antineutron
c     26      : Antilambda
c
c extensions: 
c
c     61      : u 
c     62      : u-bar
c     63      : d
c     64      : d-bar
c     65      : c
c     66      : c-bar
c     67      : s
c     68      : s-bar
c     69      : t
c     70      : t-bar
c     71      : b
c     72      : b-bar
c     73      : g
c     74      : g-bar
 
      Logical  midrap

C*******************************************************************************

      Do j_part = 1, nptls

C
C Calculate eta for each particle.
C

         gid = gtype(j_part)

         eta = tan(theta(j_part)/2.)
         If( eta .le. 0 ) Then
            eta = -999.
         Else
            eta = -log(eta)
         Endif

         rap = y(j_part)
         pperp = pt(j_part)

C
C Definition of midrapidity particle
C

         If ( abs(rap) .le. 0.5 ) Then
            midrap = .true.
         Else
            midrap = .false.
         EndIf

C
C Now loop through particles of interest.
C
         Call HFILL( 1, float(gid), 0., 1. )

         If( gid .eq. 2 .or. gid .eq. 3) Then
            Call HFILL(10, eta, 0., 1. )

         ElseIf( gid .eq. 5 .or. gid .eq. 6 )Then
            Call HFILL(10, eta, 0., 1. )

         ElseIf( gid .eq. 7 )Then
            Call HFILL(4, rap, 0., 1. )
            If (midrap) Call HFILL(42, pperp, 0., 1. )

         ElseIf( gid .eq. 8 )Then
            Call HFILL(2, rap, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(40, pperp, 0., 1. )
            
         ElseIf( gid .eq. 9 )Then
            Call HFILL(3, rap, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(41, pperp, 0., 1. )

         ElseIf( gid .eq. 10 .or. gid .eq. 16 )Then
            Call HFILL(7, rap, 0., 1. )
            If (midrap) Call HFILL(45, pperp, 0., 1. )

         ElseIf( gid .eq. 11 )Then
            Call HFILL(5, rap, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(43, pperp, 0., 1. )

         ElseIf( gid .eq. 12 )Then
            Call HFILL(6, rap, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(44, pperp, 0., 1. )

**         ElseIf( gid .eq. 14 )Then
         ElseIf( gid .eq. 14 .or. gid .eq. 13)Then
            Call HFILL(20, rap, 0., 1. )
**            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(50, pperp, 0., 1. )

**         ElseIf( gid .eq. 15 )Then
         ElseIf( gid .eq. 15 .or. gid .eq. 25)Then
            Call HFILL(21, rap, 0., 1. )
**            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(51, pperp, 0., 1. )

         ElseIf( gid .eq. 18 )Then
            Call HFILL(22, rap, 0., 1. )
            If (midrap) Call HFILL(52, pperp, 0., 1. )

         ElseIf( gid .eq. 19 )Then
            Call HFILL(24, rap, 0., 1. )
            If (midrap) Call HFILL(54, pperp, 0., 1. )

         ElseIf( gid .eq. 20 )Then
            Call HFILL(26, rap, 0., 1. )
            If (midrap) Call HFILL(56, pperp, 0., 1. )

         ElseIf( gid .eq. 21 )Then
            Call HFILL(25, rap, 0., 1. )
            If (midrap) Call HFILL(55, pperp, 0., 1. )

         ElseIf( gid .eq. 26 )Then
            Call HFILL(23, rap, 0., 1. )
            If (midrap) Call HFILL(53, pperp, 0., 1. )

         ElseIf( gid .ge. 61 .and. gid .le. 72 )Then
            Call HFILL(9, rap, 0., 1. )
            Call HFILL(46, pperp, 0., 1. )

         ElseIf( gid .eq. 73 .or. gid .eq. 74 )Then
            Call HFILL(10, rap, 0., 1. )
            Call HFILL(47, pperp, 0., 1. )
 
         EndIf

      EndDo

      istat = 0
      n_event = n_event + 1

      Write( *, *) 'SPECTRA: analysed event ',n_event

      Return
      End

C*******************************************************************************
      Subroutine anfin
C*******************************************************************************
C
C     Description:
C
C     Edits, stores histograms
C
C     Author:         Matthew Bloomer
C     Created:        9 May, 1995
C     Last modified:
C
C*******************************************************************************

      Implicit None

C  Counting number of events analysis

      integer n_event
      common/spectra/n_event


      Integer  i_hist, nchan
      Parameter (nchan=50)
      Real     h_norm
      Logical  HEXIST

      Real     x(nchan), y(nchan), ex(nchan), ey(nchan), x0(nchan)
      Data     x0 / nchan*0. /

C*******************************************************************************

C
C Normalize spectra and rapidity distributions
C
      Do i_hist = 2, 39
         If ( HEXIST(i_hist) .and. n_event .gt. 0 )Then

            h_norm = 1./(float(n_event)*0.1)
            Call HOPERA( i_hist, '+', i_hist, i_hist, h_norm, 0. )

         EndIf
      EndDo

      Do i_hist = 40, 60
         If ( HEXIST(i_hist) .and. n_event .gt. 0 )Then

            Call HCOPY( i_hist, 99, 'temp')
            Call HREBIN( 99, x, y, ex, ey, nchan, 1, nchan )
            Call HPAK( 99, x )
            Call HPAKE( 99, x0 )
            
            h_norm = float(n_event)*0.1
            Call HOPERA( i_hist, '/', 99, i_hist, 1., h_norm )
            Call HDELET(99)

         EndIf
      EndDo

C
C Store results in local directory
C

      Call HRPUT( 0, 'spectra.hbook', 'n' )

      Return
      End
