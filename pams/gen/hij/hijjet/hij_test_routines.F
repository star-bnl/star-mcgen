C
C HIJ_TEST_ROUTINES: 	hij_test_book
C			hij_test_sort
C			hij_test_end
C
C*******************************************************************************
      Integer Function HIJ_TEST_BOOK()
C*******************************************************************************
C
C     Description:
C
C     Books histograms.
C
C     Author:         Matthew Bloomer
C     Created:        2 April, 1994
C     Last modified:
C
C     Return values:
C        -1  Success
C         0  Not Success
C*******************************************************************************

      Implicit None

C
C Make size of pawc same as in hijevt.F, for compatibility of symbol
C definitions within same .a file (this seems to be an issue
C with AIX)
C
      Integer   nwpaw, hmemor
      Parameter (nwpaw=40000)

      common/pawc/hmemor(nwpaw)

C*******************************************************************************

C
C  dN/deta charged = stable charged particles
C

      Call hlimit(nwpaw)

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

      hij_test_book = -1

      Return
      End

C*******************************************************************************
      Integer Function HIJ_TEST_SORT(i_event, n_particles, geid, px, py, pz,
     >     iparent )
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
C        -1  Success
C         0  Not Success
C*******************************************************************************

      Implicit None

      Integer  i_event, j_part, gid, maxg
      Parameter (maxg=75)
      Real     mass, pperp, energy, theta, eta, y, gmass(maxg)

C
C  The GEANT mass table, according to GEANT id codes assumed in
C  HJNTOG.  All quarks/gluons are massless.  The Delta masses 
C  are taken from the JETSET parameters.
C
      Data  gmass /    0.0, 0.000511, 0.000511, 0.0,     0.10566, 
     +             0.10566, 0.13497,  0.13957,  0.13957, 0.49767, 
     +             0.49365, 0.49365,  0.93957,  0.93827, 0.93827, 
     +             0.49767, 0.54745,  1.11563,  1.18937, 1.19255,
     +             1.19747, 1.31485,  1.32133,  1.67243, 0.93957, 
     +             1.11563, 1.18937,  1.19255,  1.19743, 1.31485,
     +             1.32133, 1.67243,  0.0,      0.0,     0.0, 
     +             0.0,     0.0,      0.0,      0.0,     0.0,
     +             0.0,     0.0,      0.0,      0.0,     1.87561,  
     +             2.80925, 3.72742,  0.0,      0.78194, 1.01941,
     +             0.7699,  0.7699,   0.7699,   1.8694,  1.8646,
     +             1.234,   1.233,    1.232,    1.231,   0.0,
     +             0.0,     0.0,      0.0,      0.0,     0.0,
     +             0.0,     0.0,      0.0,      0.0,     0.0,
     +		   0.0,     0.0,      0.0,      0.0,     0.0 / 
 
      Integer  n_particles, geid(*), iparent(*)
      Real     px(*), py(*), pz(*)

      Logical  midrap

C*******************************************************************************

      Do j_part = 1, n_particles

C
C For each particle, calculate kinematics based on Geant id
C and momenta.
C

         gid = geid(j_part)
         mass = gmass(gid)

         pperp = sqrt ( px(j_part)**2 + py(j_part)**2 )
         energy = sqrt ( px(j_part)**2 + py(j_part)**2 + 
     >        pz(j_part)**2 + mass**2 )

         theta = atan2( pperp, pz(j_part) )
         eta = -log( tan (theta/2.) )

         If ( (energy-pz(j_part)) .gt. 0. )Then
            y = 0.5 * log ( (energy+pz(j_part))/(energy-pz(j_part)) )
         Else
            y = -999.
         EndIf

C
C Definition of midrapidity particle
C

         If ( abs(y) .le. 0.5 ) Then
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
            Call HFILL(4, y, 0., 1. )
            If (midrap) Call HFILL(42, pperp, 0., 1. )

         ElseIf( gid .eq. 8 )Then
            Call HFILL(2, y, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(40, pperp, 0., 1. )
            
         ElseIf( gid .eq. 9 )Then
            Call HFILL(3, y, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(41, pperp, 0., 1. )

         ElseIf( gid .eq. 10 .or. gid .eq. 16 )Then
            Call HFILL(7, y, 0., 1. )
            If (midrap) Call HFILL(45, pperp, 0., 1. )

         ElseIf( gid .eq. 11 )Then
            Call HFILL(5, y, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(43, pperp, 0., 1. )

         ElseIf( gid .eq. 12 )Then
            Call HFILL(6, y, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(44, pperp, 0., 1. )

         ElseIf( gid .eq. 14 )Then
            Call HFILL(20, y, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(50, pperp, 0., 1. )

         ElseIf( gid .eq. 15 )Then
            Call HFILL(21, y, 0., 1. )
            Call HFILL(10, eta, 0., 1. )
            If (midrap) Call HFILL(51, pperp, 0., 1. )

         ElseIf( gid .eq. 18 )Then
            Call HFILL(22, y, 0., 1. )
            If (midrap) Call HFILL(52, pperp, 0., 1. )

         ElseIf( gid .eq. 19 )Then
            Call HFILL(24, y, 0., 1. )
            If (midrap) Call HFILL(54, pperp, 0., 1. )

         ElseIf( gid .eq. 20 )Then
            Call HFILL(26, y, 0., 1. )
            If (midrap) Call HFILL(56, pperp, 0., 1. )

         ElseIf( gid .eq. 21 )Then
            Call HFILL(25, y, 0., 1. )
            If (midrap) Call HFILL(55, pperp, 0., 1. )

         ElseIf( gid .eq. 26 )Then
            Call HFILL(23, y, 0., 1. )
            If (midrap) Call HFILL(53, pperp, 0., 1. )

         ElseIf( gid .ge. 61 .and. gid .le. 72 )Then
            Call HFILL(9, y, 0., 1. )
            Call HFILL(46, pperp, 0., 1. )

         ElseIf( gid .eq. 73 .or. gid .eq. 74 )Then
            Call HFILL(10, y, 0., 1. )
            Call HFILL(47, pperp, 0., 1. )
 
         EndIf

      EndDo

      hij_test_sort = -1

      Return
      End

C*******************************************************************************
      Integer Function HIJ_TEST_END( n_event )
C*******************************************************************************
C
C     Description:
C
C     Edits, stores histograms
C
C     Author:         Matthew Bloomer
C     Created:        2 April, 1994
C     Last modified:
C
C     Return values:
C        -1  Success
C         0  Not Success
C*******************************************************************************

      Implicit None

      Integer  n_event, i_hist, nchan
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
C Store results in local directory (probably will usually be
C hij/exa).
C

      Call HRPUT( 0, 'hij_test.hbook', 'n' )
      hij_test_end = -1

      Return
      End
