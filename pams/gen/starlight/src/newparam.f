c     this subroutine determines if the luminosity function
c     needs to be recalculated or not-- will be calculated if certain 
c     input parameters (Z,A,gamma_em,numw,Wmax,numy,Ymax or gg_or_gP)
c     are different for this run than for the previous one
      subroutine newparam(new)

      implicit NONE

      include 'inputp.inc'
      include 'D2LParam.inc'
      include 'pttable.inc'
      logical new
      integer Ztest,Atest,gg_or_gPtest,numwtest,numytest,ibreakuptest
      double precision Gammatest,wmintest,wmaxtest,ymaxtest

      integer NPTtest,iinterferetest
      double precision xinterferetest,ptmaxtest

c     initialize test variables to ensure that if the file starlight.dat
c     does not already exist, the parameter new will be set to TRUE
      Ztest = 0
      Atest = 0
      Gammatest = 0.
      wmaxtest = 0.
      numwtest = 0
      ymaxtest = 0.
      numytest = 0
      gg_or_gPtest = 0
      ibreakuptest=0

c     open file and read in parameters from previous run
      open (unit=20,file='starlight.dat',status='unknown')
      read (20,*, END=50) Ztest
      read (20,*) Atest
      read (20,*) Gammatest
      read (20,*) wmaxtest
      read (20,*) wmintest
      read (20,*) numwtest
      read (20,*) ymaxtest
      read (20,*) numytest
      read (20,*) gg_or_gPtest
      read (20,*) ibreakuptest
      read (20,*) iinterferetest
      read (20,*) xinterferetest
      read (20,*) ptmaxtest
      read (20,*) NPTtest

c     if parameters are not the same as before, set new to TRUE
 50      new = .not.(  (Z.eq.Ztest).and.(A.eq.Atest).and.
     &    (gamma_em.eq.Gammatest).and.(numw.eq.numwtest).and.
     &    (numy.eq.numytest).and.(wmaxtest.eq. Wmax)
     &    .and.(wmintest.eq.wmin).and.
     &    (ymaxtest.eq.Ymax).and.(gg_or_gPtest.eq.gg_or_gP)  
     &    .and. (ibreakup .eq. ibreakuptest) .and.
     &    (iinterferetest.eq.iinterfere).and.
     &    (xinterferetest.eq.xinterfere).and.
     &    (ptmaxtest.eq.ptmax) .and. (NPTtest.eq.NPT))
      close (unit=20)
      return
      end


