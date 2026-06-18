c $Id: init.f,v 1.22 2007/01/30 14:50:25 bleicher Exp $
C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine uinit(io)
c
c Revision : 1.0
c
cinput io : flag for call to {\tt input(io)} 
c
c This subroutine calls initialization procedures for uqmd
c i.e. random generater, etc. 
c Routines called before the first (physical) event of uqmd should
c enter here. an exception is the subroutine {\tt init}. {\tt init} should NOT
c be included in uinit. 
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      include 'coms.inc'
      include 'options.inc'
      include 'boxinc.inc'
      include 'inputs.inc'
      integer io

      call set0
      call params

c.. display logo
      call urqmdlogo

c read input file
      call input(io)

      if(CTOption(28).gt.0)then
        write(6,*)'CTOption 28 >0 not used anymore. Set to 0.'
        CTOption(28)=0
      endif     

      call strini

      call output(19)

      firstseed=.true.
      fixedseed=ranseed.gt.0
      if(fixedseed)write(6,*)'fixed random number:',ranseed
      call sseed(ranseed)
      call loginit
      if(CTOption(33).eq.0.or.CTOption(9).eq.0) call loadwtab(io)
c in case of CASCADE mode, the potentials need not be calculated
      if(EoS.ne.0) then
cJS      if tabulated EoS is used, read the table
         if (CTOption(63).ge.1)then
           call readCMF(CTOption(63))
         endif
         if(logSky) call potdww
         if(logYuk) call potYuk
         if(logPau) call potPau
         if(logCb)  call potCb
      endif
c
c  calculate the normalization of resonances distribution...
c
      call norm_init

      if(io.ne.0) return

c     do not initialize projectile and target if old event is read in
      if(CTOption(40).ne.0) return

      if(boxflag.eq.0) then
c
c initialize nuclei (projectile and target) and store them
c
c initialize normal projectile
         if(prspflg.eq.0) then
c         if(eos.eq.0) then
            call cascinit(Zp,Ap,1)
c         else
c            write(6,*)'illegal EOS in init.'
c            stop 137
c         endif
         endif
c initialize normal target
         if(At.ne.0) then
            if(trspflg.eq.0) then
c            if(eos.eq.0) then
               call cascinit(Zt,At,2)
c            else
c               write(6,*)'illegal EOS in init.'
c               stop 137
c            endif
            endif
         endif
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine init
c
c     Revision : 1.0
c     This subroutine calls initialization procedures for different
c     equations of state and calculation modi.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.inc'
      include 'options.inc'
      include 'comres.inc'
      include 'inputs.inc'
      include 'freezeout.inc'
      include 'newpart.inc'
      include 'boxinc.inc'
      include 'colltab.inc'

      integer j,k,icount,npold,getspin,fchg,indsp(2),isrt,ipbm
      real*8 nucrad,dstt,dstp,pcm,eb,embeam,emtarget
      real*8 massit,ranf,pcms,dectim
      real*8 gaeq,beeq,galab,belab,ppeq,pteq
      real*8 pboost
      real*8 ratio
      integer AAp, AAt
        integer nnuc
        parameter (nnuc=11)
        save isrt,ipbm
        logical bcorr
        common /ini/ bcorr

        integer i,flagy
        real*8 alf,regula
c momenta
       real*8 mx,my,mz  
c important: never change!
      npart = 0
      npold = 0
      nbar=0
      nmes=0
      apt=0
      uid_cnt=0
c reset counters
c     all collisions/decays
      ctag  = 0
c     all decays
      dectag = 0
c     number of prod. hard resonances
      NHardRes=0
c     number of prod. soft resonances
      NSoftRes=0
c     number of prod. resonances via decay
      NDecRes=0
c     number of blocked collisions
      NBlColl=0
c     number of elastic collisions
      NElColl=0
c
      eb=0D0
c icount is the number of EXTRAordinary pro/tar combinations (i.e. pion ...)
      icount = 0
c reset particle vectors
      do 20 j=1,nmax
        spin(j)  = 0
        ncoll(j) = 0
        lstcoll(j)=0
        r0(j) = 0.0
        rx(j)    = 0.0
        ry(j)    = 0.0
        rz(j)    = 0.0
        p0(j)    = 0.0
        px(j)    = 0.0
        py(j)    = 0.0
        pz(j)    = 0.0
        frr0(j) = 0.0
        frrx(j)    = 0.0
        frry(j)    = 0.0
        frrz(j)    = 0.0
        frp0(j)    = 0.0
        frpx(j)    = 0.0
        frpy(j)    = 0.0
        frpz(j)    = 0.0
        fmass(j) = 0.0
        charge(j)= 0
        iso3(j)  = 0
        ityp(j)  = 0
        dectime(j)= 0.0
        origin(j)=0
        tform(j)=0.0
        xtotfac(j)=1.0
        uid(j)=0
         ffermpx(j) = 0.0
         ffermpy(j) = 0.0
         ffermpz(j) = 0.0
ctd
         do 21 k=1,2
            p0td(k,j)=0.d0
            pxtd(k,j)=0.d0
            pytd(k,j)=0.d0
            pztd(k,j)=0.d0
            fmasstd(k,j)=0.d0
            ityptd(k,j)=0
            iso3td(k,j)=0
 21      continue
 20   continue


      if(CTOption(40).ne.0) then
         call getoldevent
         return
      endif


      if (boxflag.eq.1) then 
         mbpx=0
         mbpy=0
         mbpz=0
         mx=0
         my=0
         mz=0

         nbar=0
         nmes=0
         flagy=edensflag
         ctoption(30)=0         ! no frozen fermi for box
         do 100 cbox=1,mbox
            if (flagy.ge.1) then
               bptpmax(cbox)=edens/mbox
            endif       
            call bptinit(cbox)
 100     Continue

c prevents a collective motion in the box
         do 143 i=1,npart
            mbpx=mbpx+px(i)
            mbpy=mbpy+py(i)
            mbpz=mbpz+pz(i)
 143     continue
         do 142 i=1,npart
            px(i)=px(i)-mbpx/npart
            py(i)=py(i)-mbpy/npart
            pz(i)=pz(i)-mbpz/npart
            call setonshell(i)
 142     continue
         
         if (flagy.ge.1) then
            alf=Regula(edens)
            do 42  i=1,npart   
               px(i) = alf*px(i)
               py(i) = alf*py(i)
               pz(i) = alf*pz(i)
               call setonshell(i)
 42         continue
         endif
         Write(*,*) 'walls selected'
         mbflag=2
         if (solid.gt.0) Write(*,*)'solid walls selected'
         
         Return
      EndIf

      if(At.ne.0) then
         if (CTParam(21).eq.0.0) then
            dstp = nucrad(Ap)+CTParam(41)
            dstt = nucrad(At)+CTParam(41)
         else
            ratio = sqrt((1 + 4.0*CTParam(21)/3.0) / 
     $           (1 - 2.0*CTParam(21)/3.0) )
            dstp = nucrad(Ap)*ratio**(2.0/3.0)+CTParam(41)
            dstt = nucrad(At)*ratio**(2.0/3.0)+CTParam(41)
         endif
c add radius offset
         dstp=dstp+CTParam(30)
         dstt=dstt+CTParam(30)      
c         
c         dst=0.5d0*(dstt+dstp)
      else
c            dst=0.d0
            dstp=0d0
            dstt=0d0
      endif

ce For anti nuclei:
      AAp = abs(Ap)
      AAt = abs(At)

c
c fix masses of projectile and target for calculation of pbeam,ecm,pcm
      if(prspflg.eq.0) then
         embeam=AAp*EMNUC
      elseif(prspflg.eq.1) then
         icount=icount+1
c!!!  sofar only pro/tar with fixed masses allowed
         embeam=massit(spityp(1))
      endif
      if(trspflg.eq.0) then
         emtarget=AAt*EMNUC
      elseif(trspflg.eq.1) then
         icount=icount+1
         emtarget=massit(spityp(2))
      endif
c         
c p(equal_speed) with given elab  cccccccccccccccccccccccccccccccccccccc

      if(srtflag.eq.0) then
         
         ebeam=AAp*ebeam
         eb=ebeam+embeam
         pbeam=sqrt(ebeam*(ebeam+2.0d0*embeam))
         
c     p(equal_speed) with given sqrt(s) ccccccccccccccccccccccccccccccccccccc
         
      elseif(srtflag.eq.1) then
         
c     in this mode, everything has to calculated on a per particle basis
         embeam=embeam/dble(AAp)
         emtarget=emtarget/dble(AAt)
         
         if(emtarget+embeam.gt.srtmin)then
            srtmin=emtarget+embeam+1d-2
            write(6,*)' *** error:initial energy below treshold'
            write(6,*)'     c.m. energy will be increased to:',
     &           srtmin
            srtmax=max(srtmax,srtmin)
         end if
         if(efuncflag.eq.0) then
            ecm=srtmin
         else                   ! if(efuncflag.eq.1) then
c     excitation function
            if(mod((event-1)*nsrt,nevents).eq.0
     &           .and.firstev.gt.0)then
               if(efuncflag.eq.1)then 
                  ecm=ecm+(srtmax-srtmin)/dble(nsrt-1)
               else if(efuncflag.eq.2)then
                  isrt=isrt+1
                  ecm=srtmin*exp( 
     &                 (dlog(srtmax/srtmin)) 
     &                 *isrt/(nsrt-1))
               end if
            elseif(firstev.eq.0) then
               isrt=0
               firstev=1
               ecm=srtmin
            endif
         endif
c     
c     this is all on a per particle basis
         pcm=pcms(ecm,embeam,emtarget)
         ebeam=sqrt(embeam**2 + (pcm*ecm/emtarget)**2) - embeam
         pbeam=pcm*ecm/emtarget
c     now revert to full quantities
         ebeam=AAp*ebeam
         pbeam=AAp*pbeam
         embeam=embeam*AAp
         emtarget=emtarget*AAt
         eb=sqrt(embeam**2+pbeam**2)
         

         
c     p(equal_speed) with given plab ccccccccccccccccccccccccccccccccccc
      elseif(srtflag.eq.2) then
         if(efuncflag.gt.0) then
c     excitation function
c     calculate the next pbeam if number of events at current pbeam exceeds nevents/npb 
            if (mod((event-1)*npb,nevents).eq.0
     &           .and.firstev.gt.0) then
c     excitation function linear in pbeam
               if (efuncflag.eq.1) then 
                  pbeam=pbeam+(pbmax-pbmin)/dble(npb-1)
c     else use a logaritmic excitation function
               else if(efuncflag.eq.2) then
                  ipbm=ipbm+1
                  pbeam=pbmin*exp( 
     &                 (dlog(pbmax/pbmin)) 
     &                 *ipbm/(npb-1))
               end if
            else if (firstev.eq.0) then
               ipbm=0
               firstev=1
               pbeam=pbmin
            endif
         endif
c     input was pbeam per particle
         pbeam=AAp*pbeam
         eb=sqrt(embeam**2+pbeam**2)
         ebeam=eb-embeam
      endif   
      
c     now do the calculation of equal_speed quantities
      galab=eb/embeam           ! gamma_lab
      belab=pbeam/eb            ! beta_lab
      gaeq=sqrt(0.5*(1+galab))  ! gamma_equal_speed
      beeq=belab*galab/(1+galab) ! beta_equal_speed
      ppeq=gaeq*beeq*embeam     ! p_projectile(eq)
      pteq=-(gaeq*beeq*emtarget) ! p_target(eq)
      
c     reduce to per particle quantities
      ppeq=ppeq/dble(AAp)
      if(AAt.ne.0) then
         pteq=pteq/dble(AAt)
         emtarget=emtarget/dble(AAt)
      endif
      embeam=embeam/dble(AAp)
      pbeam=pbeam/dble(AAp)
      ebeam=ebeam/dble(AAp)
c     the following is the NN sqrt(s)
      ecm=sqrt(embeam**2+2*eb/dble(AAp)*emtarget+emtarget**2)  
ccccccccccccccccccccccccccccccccccccccccccccc
c compute transformation betas for output

         pcm=max(1d-10,pbeam*emtarget/ecm)

         if(CTOption(27).eq.0) then
            betann=0.d0
            betatar=pcm/sqrt(pcm**2+emtarget**2)
            betapro=-(1.*pcm/sqrt(pcm**2+embeam**2))
         elseif(CTOption(27).eq.1) then
            betann=-(1*pcm/sqrt(pcm**2+emtarget**2))
            betatar=0.d0
            betapro=-(1*pbeam/sqrt(pbeam**2+embeam**2))
         elseif(CTOption(27).eq.2) then
            betann=pcm/sqrt(pcm**2+emtarget**2)
            betatar=pbeam/sqrt(pbeam**2+embeam**2)
            betapro=0.d0
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c determine impact parameter
         if(CTOption(5).eq.0) then
            bimp=bdist
         elseif(CTOption(5).eq.1) then
            if(bdist.gt.(nucrad(Ap)+nucrad(At)+2*CTParam(30))) 
     &           bdist=nucrad(Ap)+nucrad(At)+2*CTParam(30)
            bimp=sqrt(bmin**2 + ranf(0) * (bdist**2 - bmin**2))
         elseif(CTOption(5).eq.2) then
            if(bdist.gt.(nucrad(Ap)+nucrad(At)+2*CTParam(30))) 
     &           bdist=nucrad(Ap)+nucrad(At)+2*CTParam(30)
            bimp=bmin+ranf(0)*(bdist-bmin)
         else
            write(6,*)'illegal CTOption(5) :',CTOption(5)
            stop 137
         endif

         if(At.eq.0) bimp=0.d0

c initialize normal projectile
         if(prspflg.eq.0) then
            if(mod(event,nnuc).eq.0)then
               call cascinit(Zp,Ap,1)
            endif
            call getnucleus(1,npart)
            npart=npart+AAp
c change reference frame
            if (CTOption(27).eq.1) then   
               pboost = -pbeam
            elseif (CTOption(27).eq.2) then
               pboost = 0.d0
            else
               pboost = -ppeq
            endif  
            call boostnuc(npold+1,npold+AAp,
     &                    pboost,0.5*bimp,-dstp)
c save fermi motion
            if (CTOption(30).eq.1) then
               call savefermi(npold+1,npold+AAp,-pboost)
            endif
            npold=npart
            nbar=nbar+AAp
            if (CTParam(20).ne.0) then
               call getnucleus(1,npart)
               npart=npart+AAp
               call boostnuc(npold+1,npold+AAp,
     &              pboost,0.5*bimp,-dstp+CTParam(20))
               if (CTOption(30).eq.1) then
                  call savefermi(npold+1,npold+AAp,-pboost)
               endif
               npold=npart
               nbar=nbar+AAp
            endif
         endif


c initialize normal target
         if(At.ne.0) then
            if(trspflg.eq.0) then
               if(mod(event,nnuc).eq.0)then
                  call cascinit(Zt,At,2)
               endif
               call getnucleus(2,npart)
               npart=npart+AAt
c change ref. frame
               if(CTOption(27).eq.1) then  
                  pboost = 0.d0
               elseif(CTOption(27).eq.2) then
                  pboost = pbeam
               else
                  pboost = -pteq
               endif
               call boostnuc(npold+1,npold+AAt,
     &                       pboost,-(0.5*bimp),dstt)
c save fermi motion
               if (CTOption(30).eq.1) then
                  call savefermi(npold+1,npold+AAt,-pboost)
               endif
               npold=npart
               nbar = nbar + AAt
            endif
         endif

c set unique ID-tag counter (is not initialized with getnucleus calls)
         uid_cnt=npart

         if(icount.eq.0) then
c set counter for collupd 
            apt=Ap
            return
c initialize special PRO/TAR combinations
         elseif(icount.eq.1) then
            if(prspflg.eq.1) then
               indsp(1)=1
c the "regular" target sits first in the arrays
               apt=At
            else
               indsp(1)=2
               apt=Ap
            endif
         elseif(icount.eq.2) then
            if(abs(spityp(1)).le.abs(spityp(2))) then
               indsp(1)=1
               indsp(2)=2
               apt=Ap
            else
               indsp(1)=2
               indsp(2)=1
               apt=At
            endif
         endif
         do 40 j=1,icount
            npart=npart+1
            if(abs(spityp(indsp(j))).lt.minmes) then
               nbar=nbar+1
            else
               nmes=nmes+1
            endif
            iso3(npart) = spiso3(indsp(j))
            ityp(npart) = spityp(indsp(j))
            uid_cnt=uid_cnt+1
            uid(npart)=uid_cnt
            spin(npart) = getspin(ityp(npart),-1)
            charge(npart)=fchg(iso3(npart),ityp(npart))
            fmass(npart) = massit(ityp(npart))
            rx(npart) = 0.d0
            ry(npart) = 0.d0
            rz(npart) = 0.d0
            px(npart) = 0.d0 
            py(npart) = 0.d0
c       pz ist stored in pbeam,p?eq! 
            pz(npart) = 0.d0
            p0(npart)=sqrt(px(npart)**2+py(npart)**2+pz(npart)**2
     &           +fmass(npart)**2)
            if(indsp(j).eq.1) then
               if(CTOption(27).eq.1) then   
                  call boostnuc(npart,npart,-pbeam,0.5*bimp,-dstp)
               elseif(CTOption(27).eq.2) then
                  call boostnuc(npart,npart,0.d0,0.5*bimp,-dstp)
               else
                  call boostnuc(npart,npart,-ppeq,0.5*bimp,-dstp)
               endif
            elseif(indsp(j).eq.2) then
               if(CTOption(27).eq.1) then
                  call boostnuc(npart,npart,0.d0,-(0.5*bimp),dstt)
               elseif(CTOption(27).eq.2) then
                  call boostnuc(npart,npart,pbeam,-(0.5*bimp),dstt)
               else
                  call boostnuc(npart,npart,-pteq,-(0.5*bimp),dstt)
               endif
            endif
            dectime(npart) = dectim(npart,1)

 40      continue
      return
      end

      subroutine savefermi(i1,i2,p)
c
c     Revision : 1.0
c     Store fermi momentum in fferm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.inc'
      
      integer i,i1,i2
      real*8 p

      if (i1.eq.0) return

      if (ncoll(i1).gt.0) return        

      do i=i1,i2
         ffermpx(i)=px(i)
         ffermpy(i)=py(i)
         ffermpz(i)=pz(i)-p
         px(i)=0.0
         py(i)=0.0
         pz(i)=p
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine addfermi(ind,p)
c
c     Revision : 1.0
c     Restore fermi momentum from fferm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.inc'
      
      integer ind
      real*8 p

      if (ind.eq.0) return
      
      p = pz(ind)
      px(ind) = px(ind)+ffermpx(ind)
      py(ind) = py(ind)+ffermpy(ind)
      pz(ind) = pz(ind)+ffermpz(ind)
      ffermpx(ind) = 0.0
      ffermpy(ind) = 0.0
      ffermpz(ind) = 0.0

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c T. Reichert 2024      
      subroutine init_eccentricity()
c
c     Revision : 1.0
c     Calculate initial eccentricity and Glauber-like participants
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'coms.inc'

      logical Nparticipant(Ap+At)
      integer i,j
      real*8 distance,d0
      real*8 Nrx,Nry
      real*8 rx_event,ry_event,rx2_event,ry2_event
      real*8 rxy_event,rxy2_event
      real*8 sig_rx2,sig_ry2,sig_rxy
      real*8 phi
      real*8 rcos1,rsin1,rcos2,rsin2
      real*8 rcos3,rsin3,rcos4,rsin4
      real*8 r1,r2,r3,r4,Nr
      integer Ncollisions(Ap+At)

      d0=1.13d0 ! SQRT(sigtot_el/PI)

      num_part=0d0
      proton_part=0d0
      neutron_part=0d0
      Nrx=0d0
      Nry=0d0
      rx_event=0d0
      rx2_event=0d0
      ry_event=0d0
      ry2_event=0d0
      rxy_event=0d0
      rxy2_event=0d0
      sig_rx2=0d0
      sig_ry2=0d0
      sig_rxy=0d0
      eccentricity=0d0
      Nr=0d0
      r1=0d0
      r2=0d0
      r3=0d0
      r4=0d0      
      rcos1=0d0
      rsin1=0d0
      rcos2=0d0
      rsin2=0d0
      rcos3=0d0
      rsin3=0d0
      rcos4=0d0
      rsin4=0d0
      phi=0d0
      eps1=0d0
      eps2=0d0
      eps3=0d0
      eps4=0d0

      do i=1,Ap+At
       Nparticipant(i)=.false.
       Ncollisions(i)=0
      enddo

      do i=1,Ap
       do j=Ap+1,Ap+At
        distance=SQRT((rx(i)-rx(j))**2+(ry(i)-ry(j))**2)
        if(distance.le.d0)then
         Nparticipant(i)=.true.
         Nparticipant(j)=.true.
         Ncollisions(i)=Ncollisions(i)+1
         Ncollisions(j)=Ncollisions(j)+1
        endif
       enddo
      enddo

      do i=1,Ap+At
       if(Nparticipant(i))then
        num_part=num_part+1d0
        if(charge(i).eq.1)then
         proton_part=proton_part+1d0
        elseif(charge(i).eq.0)then
         neutron_part=neutron_part+1d0
        endif
        Nrx=rx(i)
        Nry=ry(i)
        Nr=SQRT(Nrx**2+Nry**2)
        rx_event = rx_event + Nrx
        ry_event = ry_event + Nry
        rx2_event = rx2_event + Nrx**2
        ry2_event = ry2_event + Nry**2
        rxy_event = rxy_event + Nrx * Nry
        rxy2_event = rxy2_event + (Nrx * Nry)**2
       else
ctr Uncomment to shift Glauber-like participants to infinity
ctr     rz(i)=pz(i)/dabs(pz(i))*1000d0+rz(i)
       if (.not.(Nry.eq.0d0.and.Nrx.eq.0d0)) then ! exclude particles at (0,0) where atan is undefined
        phi=ATAN2(Nry,Nrx)
        rcos1=rcos1+(Nr*COS(phi))
        rsin1=rsin1+(Nr*SIN(phi))
        rcos2=rcos2+(Nr**2*COS(2*phi))
        rsin2=rsin2+(Nr**2*SIN(2*phi))
        rcos3=rcos3+(Nr**3*COS(3*phi))
        rsin3=rsin3+(Nr**3*SIN(3*phi))
        rcos4=rcos4+(Nr**4*COS(4*phi))
        rsin4=rsin4+(Nr**4*SIN(4*phi))
        r1=r1+Nr
        r2=r2+Nr**2
        r3=r3+Nr**3
        r4=r4+Nr**4
       endif 
       endif
      enddo

      if(num_part.ne.0d0)then
       rx_event=rx_event/num_part
       rx2_event=rx2_event/num_part
       ry_event=ry_event/num_part
       ry2_event=ry2_event/num_part
       rxy_event=rxy_event/num_part
       rxy2_event=rxy2_event/num_part

       rcos1=rcos1/num_part
       rsin1=rsin1/num_part
       rcos2=rcos2/num_part
       rsin2=rsin2/num_part
       rcos3=rcos3/num_part
       rsin3=rsin3/num_part
       rcos4=rcos4/num_part
       rsin4=rsin4/num_part
       r1=r1/num_part
       r2=r2/num_part
       r3=r3/num_part
       r4=r4/num_part
      endif

      sig_rx2 = rx2_event - rx_event**2
      sig_ry2 = ry2_event - ry_event**2
      sig_rxy = rxy_event - rx_event * ry_event
      if (sig_rx2+sig_ry2.ne.0d0) then
       eccentricity=SQRT((sig_rx2-sig_ry2)**2+4*sig_rxy**2)
     &              /(sig_rx2+sig_ry2)
      endif

      if (r1.ne.0d0) then
       eps1=SQRT(rcos1**2+rsin1**2)/r1
      endif
      if (r2.ne.0d0) then
       eps2=SQRT(rcos2**2+rsin2**2)/r2
      endif
      if (r3.ne.0d0) then
       eps3=SQRT(rcos3**2+rsin3**2)/r3
      endif
      if (r4.ne.0d0) then
       eps4=SQRT(rcos4**2+rsin4**2)/r4
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c T. Reichert
c
c This subroutine initializes all arrays for the resonance reconstruction
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_resrec

      implicit none
      include 'comhis.inc'

      integer iww,i1

ctr set arrays for resonance reconstruction to zero
      do iww=1,wwmax
       hisnin(iww)=0
       hisnexit(iww)=0
       hisiline(iww)=0
       hisctag(iww)=0
       hisacttime(iww)=0
       hissqrts(iww)=0
       hisstot(iww)=0
       hissigpart(iww)=0
       hiscdens(iww)=0
       do i1=1,inmax
        INind(iww,i1)=0
        INr0(iww,i1)=0d0
        INrx(iww,i1)=0d0
        INry(iww,i1)=0d0
        INrz(iww,i1)=0d0
        INp0(iww,i1)=0d0
        INpx(iww,i1)=0d0
        INpy(iww,i1)=0d0
        INpz(iww,i1)=0d0
        INmass(iww,i1)=0d0
        INityp(iww,i1)=0
        INiso3(iww,i1)=0
        INch(iww,i1)=0
        INlcoll(iww,i1)=0
        INcoll(iww,i1)=0
        INistr(iww,i1)=0
        INorigin(iww,i1)=0
       enddo
       do i1=1,outmax
        OUTind(iww,i1)=0
        OUTr0(iww,i1)=0d0
        OUTrx(iww,i1)=0d0
        OUTry(iww,i1)=0d0
        OUTrz(iww,i1)=0d0
        OUTp0(iww,i1)=0d0
        OUTpx(iww,i1)=0d0
        OUTpy(iww,i1)=0d0
        OUTpz(iww,i1)=0d0
        OUTmass(iww,i1)=0d0
        OUTityp(iww,i1)=0
        OUTiso3(iww,i1)=0
        OUTch(iww,i1)=0
        OUTlcoll(iww,i1)=0
        OUTcoll(iww,i1)=0
        OUTistr(iww,i1)=0
        OUTorigin(iww,i1)=0
       enddo
      enddo
      do i1=1,resomax
       RECr0(i1)=0d0
       RECrx(i1)=0d0
       RECry(i1)=0d0
       RECrz(i1)=0d0
       RECp0(i1)=0d0
       RECpx(i1)=0d0
       RECpy(i1)=0d0
       RECpz(i1)=0d0
       RECmass(i1)=0d0
       RECityp(i1)=0
       RECiso3(i1)=0
       RECch(i1)=0
       REClcoll(i1)=0
       RECcoll(i1)=0
       RECorigin(i1)=0
       do iww=0,4
        RECind(i1,iww)=0
       enddo
       ORIr0(i1)=0d0
       ORIrx(i1)=0d0
       ORIry(i1)=0d0
       ORIrz(i1)=0d0
       ORIp0(i1)=0d0
       ORIpx(i1)=0d0
       ORIpy(i1)=0d0
       ORIpz(i1)=0d0
       ORImass(i1)=0d0
       ORIityp(i1)=0
       ORIiso3(i1)=0
       ORIch(i1)=0
       ORIlcoll(i1)=0
       ORIcoll(i1)=0
       ORIorigin(i1)=0
      enddo
      resctr=0

      return
      end
