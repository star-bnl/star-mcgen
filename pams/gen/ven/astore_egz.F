c     -*-Fortran-*-
c-----------------------------------------------------------------------
      subroutine astore_egz
c     writes the results of a simulation on tape.
c     contains a description of the stored variables.
c-----------------------------------------------------------------------
c     
c     Common and variables relating to VENUS
c     ======================================
#include "venus.inc"
c     
c     include PAWC common in order to initialise KUIP...
c     
      common/nevnt/nevnt
      common/outfgdt/ifgdt
c  outfgdt directs output
c	flow.ifgdt=1 is standard ZEBRA output(matt bloomer) 
c	     ifgdt=2 is standard text output(richard morris)
c	     ifgdt=4 is gstar text output (peter jacobs)
c	     ifgdt=8 is gstar xdf output (craig tull)
      common/pawc/hmemor(40000)
c     
      integer egzinit
      integer egzout_venus
      real scfermi
c     
c     Global Declarations
c     ===================
c     
#include "evt_ctrl.inc"
#include "evt_zebra.inc"
c     
C- XDF OUTPUT --------------------------------------------------------
#include "table_header.inc"
#include "eg_gener.inc"
#include "eg_event.inc"
#include "eg_track.inc"
#include "eg_vertex.inc"
CC
      INTEGER*4 GMAXLEN
      PARAMETER (GMAXLEN=1)
      RECORD /EG_GENER_ST/ gtable(GMAXLEN)
      INTEGER*4 GNOK
      RECORD /TABLE_HEAD_ST/ gtable_h
CC
      INTEGER*4 EMAXLEN
      PARAMETER (EMAXLEN=1)
      RECORD /EG_EVENT_ST/ etable(EMAXLEN)
      INTEGER*4 ENOK
      RECORD /TABLE_HEAD_ST/ etable_h
CC
      INTEGER*4 TMAXLEN
      PARAMETER (TMAXLEN=100000)
      RECORD /EG_TRACK_ST/ ttable(TMAXLEN)
      INTEGER*4 TNOK
      RECORD /TABLE_HEAD_ST/ ttable_h
CC
      INTEGER*4 VMAXLEN
      PARAMETER (VMAXLEN=50000)
      RECORD /EG_VERTEX_ST/ vtable(VMAXLEN)
      INTEGER*4 VNOK
      RECORD /TABLE_HEAD_ST/ vtable_h
C---------------------------------------------------------------------


      character*80 chkine
      logical first
c     
      data first  /.true./
      data scfermi/0.100E-12/
c
c MAB

      integer  i, length 
      Real*8   pz, x, u, y, mperp
      logical  ynn
      character*4 chfram
      data chfram /'CMS '/

c MAB end

      evtcod = venus
c     
      ish0=ish
      ishs0=ishsub
      if(ishevt.ne.0.and.nrevt.ne.ishevt)ish=0
      if(ishevt.ne.0.and.nrevt.ne.ishevt)ishsub=0
      if(ishsub/100.eq.8)ish=mod(ishsub,100)
      if(ish.eq.14)call uttima('*** astore *** ')
      if(ish.eq.17)write(ifch,*)' '
      if(ish.eq.17)write(ifch,*)'astore (entry)'
      if(ichoic.eq.4.or.ichoic.eq.5)goto1000
      if(nevt.ne.1)goto1000
      nptlu=max0(nptl,nptlu)
      if(iwzzzz.eq.1)goto9999
      if(imihis.eq.1)goto9999
      if(isphis.eq.1)goto9999
      if(iclhis.eq.1)goto9999
      if(ijphis.eq.1)goto9999
      levt=1
      lptl=3
      
      nptls=0
      ivtxs=0
      xstorx=100000.
      xstory=100000.
      xstorz=100000.
      do 1 i=1,nptl
          nrptla(i)=-999
          if(istptl(i).gt.istmax)goto1
        if(xorptl(1,i).eq.xstorx.and.xorptl(2,i).eq.xstory.and.xorptl(3,i).eq.xstorz) then
        else
        ivtxs=ivtxs+1
        endif
          nptls=nptls+1
          nrptla(i)=nptls
      xstorx=xorptl(1,i)
      xstory=xorptl(2,i)
      xstorz=xorptl(3,i)
 1    continue
      
c     
c     if first time write reaction information
c     

      if (first) then
          first = .false.
CCCCC STRATOS
      if(ifgdt.eq.2) OPEN(UNIT=32, FILE='venus.txt', status='unknown')
      if(IBITS(ifgdt,2,1).NE.0)
     &		OPEN(UNIT=32, FILE='venus.gst', status='unknown')
      if(IBITS(ifgdt,2,1).NE.0) then
        write(32,3) laproj,maproj,latarg,matarg,engy,bmaxim
 3      format('GENER: ven 4.12 ',4I12,2G12.5,' CMS ')
      endif
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
	gnok = 1
CRAIG IS DONE EDITING HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	gtable(gnok).east_z			= laproj
        gtable(gnok).east_a			= maproj
        gtable(gnok).west_z			= latarg
        gtable(gnok).west_a			= matarg
        gtable(gnok).sqrts			= engy
        gtable(gnok).b_max			= bmaxim
      Endif
C---------------------------------------------------------------------
CCCCC STRATOS
c     
c     Initialise KUIP in order to use kupryn...
c     
          call hlimit(40000)
          call kuinit(10000)
c     Initialize ZEBRA
          call mzebra(-3)
c     Initialize EGZ package
          z_input = .false.
c     init egz package, 0 for not GEANT
          chkine = 'venus_egz.data'
          istat = egzinit(0, chkine)
          if (istat.ne.0) then
              write(6,*)'Failed zebra file open'
              stop
          endif

c MAB specify frame of reference of output first time round

 220      call kuproc('Which frame of reference:  projectile (PROJ)'//
     >         ' or center of mass (CMS)', chfram, length)

          If ( length .eq. 0 .or. 
     >         ( chfram .ne. 'PROJ' .and. chfram .ne. 'CMS' )) Then
             Write (*, * ) 'Invalid entry for reference frame!'
             Stop
          Else
             If( chfram .eq. 'CMS' ) Then
                ynn = .true.
             Else
                ynn = .false.
             Endif
          Endif
c MAB end

      endif
      
      if(ichoic.eq.3)goto1000
      

c     store particle variables:
c     
c MAB
c Transform to yhaha frame if 'CMS' has been selected as
c frame of reference.  
c
c Usual Lorentz transformation does not work well (ie, 
c p* = gamma*p - gamma*beta*E).  The problem is that sometimes 
c pz>E especially at rapidities of 2*yhaha, for yhaha=5.3.  
c This is due to imprecision of real*4 array pptl.
c
c   
c  

      If ( ynn ) Then

         Do i = 1, nptl

            pz = pptl(3,i)
            mperp = sqrt ( pptl(1,i)**2 + pptl(2,i)**2 + pptl(5,i)**2 )

c Use relation pz = mperp * sinh(y) to get rapidity of particle.
c Note the use of a log expression to get asinh.

            If ( mperp .gt. 0. ) Then
               x = pz/mperp
               u = sqrt( x**2 + 1. )
               y = log ( x + u )    ! asinh
            Else
               Write (*, *) 'astore_egz: strange mperp value of ', 
     >              mperp
               y = 999.
            EndIf

            pptl(3,i) = mperp*sinh(y-yhaha)
            pptl(4,i) = mperp*cosh(y-yhaha)

         EndDo

      Endif              ! CMS ref frame

c MAB end

      if(ifgdt.eq.1) istat = egzout_venus ()
      if(istat.ne.0) then
          write(6,*)'Output error'
          stop
      endif

CCCC STRATOS
CCCC Try to write text output.
      if(ifgdt.eq.2) then
      write(32,98) nptls,'event #', nrevt, ' ID     Px(GeV)','  Py(GeV)', '  Pz(GeV)'
98    format(1x,I6,a,I6,3x,a,5x,a,3x,a)
      Do 9918 i=1, nptl
        if(nrptla(i).lt.0) goto 9918
        Call isa_to_geant(idptl(i),idtype)
c        write(32, 9910) idtype, pptl(1,i), pptl(2,i), pptl(3,i), xorptl(1,i),xorptl(2,i), xorptl(3,i), xorptl(4,i)
        write(32, 9911) idtype, pptl(1,i), pptl(2,i), pptl(3,i)
 9911   format(17x,I3,2x,G12.5,1x,G12.5,1x,G12.5)
 9918 Continue
      endif
      if(IBITS(ifgdt,2,1).NE.0				!gst text
     &.OR.IBITS(ifgdt,3,1).NE.0) then			!gst xdf
      if(IBITS(ifgdt,2,1).NE.0)
     &	write(32,91) nrevt,nptls,ivtxs,bimevt,phievt
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
      enok = 1
      etable(enok).n_event		= nrevt
      etable(enok).n_track		= nptls
      etable(enok).n_vertex		= ivtxs
      etable(enok).b_impact		= bimevt
      etable(enok).phi_impact		= phievt
      tnok = 0
      vnok = 0
      Endif
C---------------------------------------------------------------------
 91   format('EVENT:',3x,3(1x,i6),1x,G12.5,1x,G12.5)
      it=0
      iv=0
      np=0
      iz=0
      xstorx=100000.
      xstory=100000.
      xstorz=100000.
      Do 9910 i=1, nptl
        if(nrptla(i).lt.0) goto 9910
        Call isa_to_geant(idptl(i),idtype)
        it=it+1
        np=np+1
        if(xorptl(1,i).eq.xstorx.and.xorptl(2,i).eq.xstory.and.xorptl(3,i).eq.xstorz) then
      if(IBITS(ifgdt,2,1).NE.0)
     &	write(32,92) idtype, pptl(1,i), pptl(2,i), pptl(3,i), it, iv
     &	, iz, idtype
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
      if(tnok.LT.tmaxlen)then
	 tnok = tnok+1
	 ttable(tnok).ge_pid		= idtype
	 ttable(tnok).p(1)		= pptl(1,i)
	 ttable(tnok).p(2)		= pptl(2,i)
	 ttable(tnok).p(3)		= pptl(3,i)
	 ttable(tnok).itrack		= it
	 ttable(tnok).ivertex		= iv
	 ttable(tnok).iz		= iz
      endif
      Endif
C---------------------------------------------------------------------
        else
        np=1
        iv=iv+1
        if(iv.gt.1) then
      if(IBITS(ifgdt,2,1).NE.0)
     &	write(32,93) xstrx,xstry,xstrz,xstrt,ivstor,iz,iz,npstor
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
      if(vnok.LT.vmaxlen)then
	 vnok = vnok+1
	 vtable(vnok).x(1)		= xstrx
	 vtable(vnok).x(2)		= xstry
	 vtable(vnok).x(3)		= xstrz
	 vtable(vnok).t		= xstrt
	 vtable(vnok).ivstor		= ivstor
	 vtable(vnok).iz		= iz
	 vtable(vnok).npstor		= npstor
      endif
      Endif
C---------------------------------------------------------------------
        endif
      if(IBITS(ifgdt,2,1).NE.0)
     &	write(32,92) idtype, pptl(1,i), pptl(2,i), pptl(3,i), it, iv
     &		, iz, idtype
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
      if(tnok.LT.tmaxlen)then
	 tnok = tnok+1
	 ttable(tnok).ge_pid		= idtype
	 ttable(tnok).p(1)		= pptl(1,i)
	 ttable(tnok).p(2)		= pptl(2,i)
	 ttable(tnok).p(3)		= pptl(3,i)
	 ttable(tnok).itrack		= it
	 ttable(tnok).ivertex		= iv
	 ttable(tnok).iz		= iz
      endif
      Endif
C---------------------------------------------------------------------
        endif
 92   format('TRACK:',1x,i6,3(1x,g12.5),4(1x,i6))
 93   format('VERTEX:',4(1x,g12.5),4(1x,i6))
      xstorx=xorptl(1,i)
      xstory=xorptl(2,i)
      xstorz=xorptl(3,i)
      xstrx=xorptl(1,i)*scfermi
      xstry=xorptl(2,i)*scfermi
      xstrz=xorptl(3,i)*scfermi
      xstrt=xorptl(4,i)*scfermi
      ivstor=iv
      npstor=np
9910   Continue
      if(IBITS(ifgdt,2,1).NE.0)
     &	write(32,93) xstrx,xstry,xstrz,xstrt,ivstor,iz,iz,npstor
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
      if(vnok.LT.vmaxlen)then
	 vnok = vnok+1
	 vtable(vnok).x(1)		= xstrx
	 vtable(vnok).x(2)		= xstry
	 vtable(vnok).x(3)		= xstrz
	 vtable(vnok).t		= xstrt
	 vtable(vnok).ivstor		= ivstor
	 vtable(vnok).iz		= iz
	 vtable(vnok).npstor		= npstor
      endif
      Endif
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      If(IBITS(ifgdt,3,1).NE.0)Then
      gtable_h.nok = gnok
      gtable_h.maxlen = gmaxlen
      etable_h.nok = enok
      etable_h.maxlen = emaxlen
      ttable_h.nok = tnok
      ttable_h.maxlen = tmaxlen
      vtable_h.nok = vnok
      vtable_h.maxlen = vmaxlen
      call astore_xdf(gtable_h,gtable
     &               ,etable_h,etable
     &               ,ttable_h,ttable
     &               ,vtable_h,vtable)
      Endif
C---------------------------------------------------------------------
      ifinal=-999
      iz=0
      aza=0.0
      if(IBITS(ifgdt,2,1).NE.0)then
     		if(nrevt.ge.nevnt) write(32,91) ifinal,iz,iz,aza,aza
      endif
      endif
CCCC STRATOS

c     lptl ......... record label (lptl=3)
c     nrevt ........ event number
c     nrptl ........ particle number
c     i ............ original ptl number
c     idptl ........ particle id
c     pptl ......... 5-momentum (px,py,pz,en,mass) in lab
c     ioptl ........ origin (-999:parent not stored, -1,0:no parent)
c     joptl ........ origin (second parent)
c     istptl ....... stable (=0) or not (=1)
c     xorptl ....... space-time point (x,y,z,t) on ptl track (pp-cm)
c     tivptl ....... time interval of existence
c     nqjptl ....... quark numbers of jets
      
 5    continue
 4    continue
      goto1000
      
c     fill histograms
c     ---------------
 9999 continue
      if(iwzzzz.eq.1)call wzzzz
      if(imihis.eq.1)call wimi
      if(isphis.eq.1)call wspa
      if(iclhis.eq.1)call wclu
      if(ijphis.eq.1)call wjpsi
      
 1000 continue
      call ranfgt(seedf)
      call uttimx(timaa)
      if(isup.ne.1.and.mod(nrevt,modsho).eq.0)
     *write(ifmt,*)nrevt,seedf,timaa
      seedc=seedf
      if(ish.eq.14)call uttima('    astore f   ')
      if(ish.eq.17)write(ifch,*)'astore (exit)'
      ish=ish0
      ishsub=ishs0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine add_gener_row(laproj,maproj,latarg,matarg,engy,bmaxim,
c    &		maxlen,g,nok)
c     integer*4 laproj, maproj, latarg, matarg
c     real*4 engy, bmaxim
c     integer*4 maxlen,nok
c     record /GENER_ST/ g(*)
c
c     if( nok < maxlen )then
c        nok = nok + 1
c        g(nok).laproj			= laproj
c        g(nok).maproj			= maproj
c        g(nok).latarg			= latarg
c        g(nok).matarg			= matarg
c        g(nok).engy			= engy
c        g(nok).bmaxim			= bmaxim
c     endif

c     return
c     end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine add_event_row(nrevt,nptls,ivtxs,bimevt,phievt,
c    &		maxlen,e,nok)
c     integer*4 nrevt, nptls, ivtxs
c     real*4 bimevt, phievt
c     integer*4 maxlen,nok
c     record /EVENT_ST/ e(*)
c
c     if(nok.LT.maxlen)then
c        nok = nok+1
c        e(nok).n_event		= nrevt
c        e(nok).n_track		= nptls
c        e(nok).n_vertex		= ivtxs
c        e(nok).b_impact		= bimevt
c        e(nok).phi_impact		= phievt
c     endif
c     return
c     end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine add_track_row(idtype, pptl1, pptl2, pptl3, it, iv, iz,
c    &		maxlen,t,nok)
c     integer*4 idtype
c     real*4 pptl1, pptl2, pptl3
c     integer*4 it, iv, iz,
c     integer*4 maxlen,nok
c     record /TRACK_ST/ t(*)
c
c     if(nok.LT.maxlen)then
c        nok = nok+1
c        t(nok).ge_pid		= idtype
c        t(nok).p(1)		= pptl1
c        t(nok).p(2)		= pptl2
c        t(nok).p(3)		= pptl3
c        t(nok).itrack		= it
c        t(nok).ivertex		= iv
c        t(nok).iz		= iz
c     endif
c     return
c     end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine add_vertex_row(xstrx,xstry,xstrz,xstrt,ivstor,iz,npstor,
c    &		maxlen,v,nok)
c     integer*4 maxlen,nok
c     record /VERTEX_ST/ v(*)
c
c     if(nok.LT.maxlen)then
c        nok = nok+1
c        v(nok).x(1)		= xstrx
c        v(nok).x(2)		= xstry
c        v(nok).x(3)		= xstrz
c        v(nok).t		= xstrt
c        v(nok).ivstor		= ivstor
c        v(nok).iz		= iz
c        v(nok).npstor		= npstor
c     endif
c     return
c     end
