c $Id: output.f,v 1.24 2007/05/23 14:28:50 bleicher Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine output(iunit)
c
c     Revision : 1.0
c
c     This subroutine writes the event-header to file(iunit)
C
c
cinput iunit  : output-unit
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'
      include 'comhis.inc'

c
      integer iunit,i,ttime,iu,app,att,zpp,ztt,l
      integer iiunit,isunit, id, pdgid,iout
      integer timestep,itotcoll,iinelcoll
      real*8  sigmatot,ptsigtot,stot,otime

      common /outco2/sigmatot


      character*3 abox3, abox4, abox5
      character*5 abox6
      character*4 reffram
      character*20 aa,ah,ai,ak,abox1,abox2
      character*36 ae,abt
      character*31 aee
      character*15 ab,aj,al,am
      character*13 ac,ag,pds,tds
      character*12 ad
      character*7 af
      character*9 ag2
      character*1 add
      character*171 apa14,apa15,apav,line
      character*2 apa,aop
      character*35 aboxhead
      character*54 eccAp1,eccAp2

      real*8 ecc,ecc1,ecc2,ecc3,ecc4
      real*8 number_part,numprot_part,numneut_part

c file15out
      integer ind,ind1,ind2,nin
      integer istr,ich,ii,iid

      real*8 sqrts, sigpart, colldens, cdens,cdens_
      logical bdum,paulibl

      include 'outcom.inc'

      integer fchg,strit
      character*1 echar
c     temporary arrays for CTO and CTP when read in from old event
c     (before they are overwritten)
      integer CTOtmp(numcto)
      real*8 CTPtmp(numctp)
c     CTOtc and CTPtc are the temporary fields for CTOdc and CTPdc.
      character CTOtc(numcto)*2
      character CTPtc(numctp)*2
      integer ctoforty, ctofoone
      integer ctolines,ctplines

      integer iou(13:20)

chp variables necessary for hydro evolution output in f15
      real*8 thydro_start,thydro
chp hydro flag for visualization output
      logical hydro_flag
chp new variable for vis-out to count correct npart
      integer npart_form


      save

      data iou/13,14,15,16,17,18,19,20/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              output formats
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c fileheader
 101  format(a20,3i7,a15,i2)
 301  format(a13,a13,i4,i4,a12,a13,i4,i4,a1)
c 305  format(a36,3f10.7)
 304  format(a36,3f6.2,a31,1f9.2)
 302  format (a7,i9,a13,i12,a9,a20,i7,a20,f11.3)
c 303  format(a20,i3,a15,e10.4,a15,e10.4,a15,e10.4)
 102  format(a2,15(i3,a2))
 103  format(a2,12(e11.4,a2))
 306  format(a171)

 305  format(a36,3f11.7)
 303  format(a20,i3,a15,e11.4,a15,e11.4,a15,e11.4)

c standard particle information vector
 201  format(9e16.8,i11,2i3,i9,i5,i4)
cLHC 201  format(9e24.16,i11,2i3,i9,i5,i4)

c special output for cto40 (restart of old event)
 210  format(9e16.8,i11,2i3,i9,i5,i10,3e16.8,i8)
cLHC 210  format(9e24.16,i11,2i3,i9,i5,i10,3e24.16,i8)

c special output for mmaker
 203  format(9e16.8,i5,2i3,i6,i5,i4,i5,2e16.8)
cLHC 203  format(9e24.16,i5,2i3,i6,i5,i4,i5,2e24.16)

c same with index for file15
 501  format(i5,9e16.8,i11,2i3,i9,i5,i3,i15)
cLHC 501  format(i5,9e24.16,i11,2i3,i9,i5,i3,i15)

c enhanced file16
 503  format(9e15.7,i11,2i3,i9,i5,i4,2i4)
cLHC 503  format(9e24.16,i11,2i3,i9,i5,i4,2i4)

c same including freeze-out coordinates
 213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8)
cLHC 213  format(9e24.16,i11,2i3,i9,i5,i4,8e24.16)

c collsision stats for file14
 202  format(8i8)
c same with EndOfEvent tag for file16
 602  format(a1,8i8)

c header-line for each collision in file15
 502  format(i8,i8,i4,i7,f8.3,4e12.4)

c eccentricty
 1010 format(a54,3f6.1,2x,a54,5f10.6)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c header-line for box-mode
 504  format(a20,e14.6,a20,e14.6,a3,i1,a3,i1,a3,i3)
 505  format(a35)
 506  format(a5,2i4,i8,e14.6)

c
      if(iunit.eq.17)return
      if(bf13.and.(iunit.eq.13)) return
      if(bf14.and.(iunit.eq.14)) return
      if(bf15.and.(iunit.eq.15)) return
      if(bf16.and.(iunit.eq.16)) return

c     copy projectile/target info to local vars
      app=ap
      zpp=zp
      att=at
      ztt=zt
c     copy eccentricity/Glauber-participants to local vars      
      ecc=eccentricity
      number_part=num_part
      numprot_part=proton_part
      numneut_part=neutron_part
      ecc1=eps1
      ecc2=eps2
      ecc3=eps3
      ecc4=eps4

      if(iunit.eq.19) return
c
      aa='UQMD   version:     '
      ab='  output_file '
      abt='transformation betas (NN,lab,pro) '
      ac='projectile:  '
      ad='   target: '
      add=' '
      ae='impact_parameter_real/min/max(fm):  '
      aee='  total_cross_section(mbarn):  '
      af='event# '
      ag=' random seed:'
      ah='equation_of_state: '
      ai=' total_time(fm/c): '
      aj='  E_lab(GeV/u):'
      ak=' Delta(t)_O(fm/c): '
      al='  sqrt(s)(GeV):'
      am='  p_lab(GeV/u):'
      apa='pa'
      aop='op'
      line=''
      eccAp1='Participants_Glauber (Nuc_part,prot_part,neut_part):  '
      eccAp2='Eccentricities_Glauber (ecc{2},eps1,eps2,eps3,eps4):  '


      abox1='boxmode length(fm): '
      abox2=' tot. energy (GeV): '
      abox3=' s:'
      abox4=' p:'
      abox5=' #:'
      abox6='box: '
      aboxhead='boxh ityp 2i3       N     pmax(GeV)'
      apa14='pvec: '//
     & 'r0              rx              ry              rz          '//
     & '    p0              px              py              pz      '//
     & '        m          ityp 2i3 chg lcl#  ncl or'
      apa15='pvec:ind   '//
     & 'r0              rx              ry              rz          '//
     & '    p0              px              py              pz      '//
     & '        m          ityp 2i3 chg lcl#  ncl st'
      if(iunit.eq.15) then
         apav=apa15
      else
         apav=apa14
      endif

      if(fixedseed) then
         ag2=' (fixed) '
      else
         ag2=' (auto)  '
      endif
      if(prspflg.eq.1) then
         pds='(ityp, char) '
         app=spityp(1)
         zpp=fchg(spiso3(1),app)
      else
         pds='(mass, char) '
      endif
      if(trspflg.eq.1) then
         tds='(ityp, char) '
         att=spityp(2)
         ztt=fchg(spiso3(2),att)
      else
         tds='(mass, char) '
      endif

c determine cross section of the projectile-target system
      sigmatot = ptsigtot()
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      otime=outsteps*dtimestep
      ttime=int(nsteps*dtimestep+0.01)

cbb File 15 has the abbreviated event header unless CTO(58) is 1
      if(iunit.eq.15.and.CTOption(58).ne.1)then
       write(iou(15),502)-1,event,Ap,At,bimp,ecm
     ,     ,sigmatot,ebeam,pbeam
      else
cbb How many cto-lines and ctp-lines are to be written:
cbb write all if extended output is on OR if legacy mode CTO(57) is NOT
cbb set.
      if (CTOption(41).ne.0.or.CTOption(57).ne.0) then
        ctolines = 4  ! * 15 = 60 CTOs
        ctplines = 7  ! * 12 = 84 CTPs
      else
        ctolines = 3  ! * 15 = 45 CTOs
        ctplines = 4  ! * 12 = 48 CTPs
      endif
      write(iou(iunit),101) aa,version, sigver, laires, ab,iunit
      write(iou(iunit),301) ac,pds, App, Zpp, ad,tds, Att, Ztt,add
      write(iou(iunit),305) abt,betann,betatar,betapro
      write(iou(iunit),304) ae,bimp,bmin,bdist,aee,sigmatot
      write(iou(iunit),1010) 
     &                  eccAp1,number_part,numprot_part,numneut_part,
     &                  eccAp2,ecc,ecc1,ecc2,ecc3,ecc4
      write(iou(iunit),303) ah,eos,aj,ebeam,al,ecm,am,pbeam
      write(iou(iunit),302) af,event,ag,ranseed,ag2,ai,ttime,ak,otime
      do ii = 0,ctolines-1
       write(iou(iunit),102) aop,(CTOption(i),CTOdc(i)
     &                                       ,i=ii*15+1,ii*15+15)
      enddo
      do ii = 0,ctplines-1
       write(iou(iunit),103) apa,(CTParam(i),CTPdc(i)
     &                                      ,i=ii*12+1,ii*12+12)
      enddo
      if(boxflag.eq.1) then
         write(iou(iunit),504) abox1, lbox, abox2, edens, abox3, solid,
     1                         abox4, para, abox5, mbox
         write(iou(iunit),505) aboxhead
         do 507 l=1,mbox
            write(iou(iunit),506) abox6, bptityp(l), bptiso3(l),
     1                            bptpart(l), bptpmax(l)
 507     continue
      end if
      write(iou(iunit),306) apav
      end if

c
      return
c.....
      entry uounit(iiunit,isunit)
      iou(iiunit)=isunit
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file14out(timestep)
c
c     Revision : 1.0
c
c     This subroutine writes the standard output-file (unit 14)
c
cinput timestep  : timestep of output
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

c
      if(bf14)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      write(iou(14),*) npart,ttime
      write(iou(14),202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes

c now write particle-output

      if(CTOption(41).eq.0) then
         do 13 i=1,npart
            write(iou(14),201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100)
 13      continue
      else
         do 31 i=1,npart
            write(iou(14),210) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),origin(i),
     @           dectime(i),tform(i),xtotfac(i),uid(i)
 31      continue
      endif
c
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file13out(timestep)
c
c     Revision : 1.0
c
c     This subroutine writes the standard output-file (unit 13),
c     including the freeze-out configuration of the particles
c
cinput timestep  : timestep of output
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

c
      if(bf13)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl

      if(CTOption(64).eq.0)write(iou(13),*) npart,ttime
      if(CTOption(64).eq.1)write(iou(13),*) npartcoal,ttime

      write(iou(13),202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes

c Reconstruct resonances if cto(68) is 1
      if (CTOption(68).eq.1) then
         call reconstruct_resonances()
      endif      

c now write particle-output


      do 90 i=1,npart
        if(lstcoll(i).ge.-1)then
         if(ncoll(i).eq.0) then
          if(CTOption(68).eq.0)then
            write(iou(13),213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           r0(i),rx(i),ry(i),rz(i),p0(i),px(i)+ffermpx(i),
     @           py(i)+ffermpy(i),pz(i)+ffermpz(i)
           else
            write(iou(13),22213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           r0(i),rx(i),ry(i),rz(i),p0(i),px(i)+ffermpx(i),
     @           py(i)+ffermpy(i),pz(i)+ffermpz(i),uid(i),0,0,0,0
           endif
         else
           if(CTOption(68).eq.0)then
            write(iou(13),213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           frr0(i),frrx(i),frry(i),frrz(i),frp0(i),frpx(i),
     @           frpy(i),frpz(i)
           else
            write(iou(13),22213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           frr0(i),frrx(i),frry(i),frrz(i),frp0(i),frpx(i),
     @           frpy(i),frpz(i),uid(i),0,0,0,0
           endif
         endif
         endif
 90   continue

ctr write reconstructable resonances after all stable particles       
       if (CTOption(68).eq.1) then

        do i=1,resctr
         write(iou(13),22213)ORIr0(i),ORIrx(i),ORIry(i),ORIrz(i),
     @           ORIp0(i),ORIpx(i),ORIpy(i),ORIpz(i),
     @           RECmass(i),
     @           RECityp(i),RECiso3(i),RECch(i),
     @           REClcoll(i),RECcoll(i),mod(RECorigin(i),100),
     @           RECr0(i),RECrx(i),RECry(i),RECrz(i),RECp0(i),RECpx(i),
     @           RECpy(i),RECpz(i),RECind(i,0),
     @           RECind(i,1),RECind(i,2),RECind(i,3),RECind(i,4)
        enddo
       endif


cJSF write nuclear clusters from Coalescence

       if(CTOption(64).eq.1)then
          if(CTOption(68).eq.0)then
           do i=1,nclus
            write(iou(13),213) 0d0,0d0,0d0,0d0,0d0,
     @           0d0,0d0,0d0,mcl(i),
     @           itypcl(i),isocl(i),chgcl(i),
     @           lstcollcl(i),ncollcl(i),0,
     @           r0cl(i),rxcl(i),rycl(i),rzcl(i),p0cl(i),pxcl(i),
     @           pycl(i),pzcl(i)
           enddo
          else
           do i=1,nclus
            write(iou(13),22213) 0d0,0d0,0d0,0d0,0d0,
     @           0d0,0d0,0d0,mcl(i),
     @           itypcl(i),isocl(i),chgcl(i),
     @           lstcollcl(i),ncollcl(i),0,
     @           r0cl(i),rxcl(i),rycl(i),rzcl(i),p0cl(i),pxcl(i),
     @           pycl(i),pzcl(i),0,0,0,0,0
           enddo

       endif
       endif       

22213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8,i6,4i6)         

c
      return

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      entry file15out(ind1,ind2,sqrts,stot,sigpart)
c
c     Revision : 1.0
c
c     This subroutine writes information about the in-channel to file15
c     (the collision statistics file)
c
cinput        ind1    : index of particle 1
cinput        ind2    : index of particle 2 (=0 for decay of {\tt ind1})
cinput        sqrts   : $\sqrt{s}$ of collision
cinput        stot      : total cross section
cinput        sigpart   : partial cross section
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c determine tag for scatter-input or decay-input
c and store entry channel in temporary observables
      bdum=paulibl(ind1,cdens,-1000)
      tsqrts=sqrts
      tstot=stot
      tsigpart=sigpart

      tind(1)=ind1
      tr0(1)=r0(ind1)
      trx(1)=rx(ind1)
      try(1)=ry(ind1)
      trz(1)=rz(ind1)
      tp0(1)=p0(ind1)
      tpx(1)=px(ind1)
      tpy(1)=py(ind1)
      tpz(1)=pz(ind1)
      tm(1)=fmass(ind1)
      tityp(1)=ityp(ind1)
      tiso3(1)=iso3(ind1)
      tstrange(1) = strit(tityp(1))
      tcoll(1) = ncoll(ind1)
      tlcoll(1)=lstcoll(ind1)
      torigin(1)=origin(ind1)
      tuid(1)=uid(ind1)
      if(ind2.le.0) then
         nin=1
      elseif(ind2.gt.0) then
         bdum=paulibl(ind2,cdens_,-1000)
         cdens=5d-1*(cdens+cdens_)
         nin=2
         tind(2)=ind2
         tr0(2)=r0(ind2)
         trx(2)=rx(ind2)
         try(2)=ry(ind2)
         trz(2)=rz(ind2)
         tp0(2)=p0(ind2)
         tpx(2)=px(ind2)
         tpy(2)=py(ind2)
         tpz(2)=pz(ind2)
         tm(2)=fmass(ind2)
         tityp(2)=ityp(ind2)
         tiso3(2)=iso3(ind2)
         tstrange(2)=strit(tityp(2))
         tcoll(2) = ncoll(ind2)
         tlcoll(2)=lstcoll(ind2)
         torigin(2)=origin(ind2)
         tuid(2)=uid(ind2)
      endif

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15outch(colldens)

c Memorize history to reconstruct resonances      
      if (CTOption(68).eq.1.and..not.bf13) then
       www=ctag
       if(www.gt.wwmax)then
        write(*,*)'WARNING! ','Increase wwmax!'
       endif
       hisnin(www)=nin
       hisnexit(www)=nexit
       hisiline(www)=iline
       hisctag(www)=ctag
       hisacttime(www)=acttime
       hissqrts(www)=tsqrts
       hisstot(www)=tstot
       hissigpart(www)=tsigpart
       hiscdens(www)=cdens
       do i=1,nin
        istr=strit(tityp(i))
        ich=fchg(tiso3(i),tityp(i))
        iid=tuid(i)
        INind(www,i)=iid
        INr0(www,i)=tr0(i)
        INrx(www,i)=trx(i)
        INry(www,i)=try(i)
        INrz(www,i)=trz(i)
        INp0(www,i)=tp0(i)
        INpx(www,i)=tpx(i)
        INpy(www,i)=tpy(i)
        INpz(www,i)=tpz(i)
        INmass(www,i)=tm(i)
        INityp(www,i)=tityp(i)
        INiso3(www,i)=tiso3(i)
        INch(www,i)=ich
        INlcoll(www,i)=tlcoll(i)
        INcoll(www,i)=tcoll(i)
        INistr(www,i)=istr
        INorigin(www,i)=torigin(i)
       enddo
       if(nexit.gt.outmax)then
        write(*,*)'WARNING! ','Increase outmax!'
       endif
       do ii=1,nexit
        i=inew(ii)
        istr=strit(ityp(i))
        iid=uid(i)
        OUTind(www,ii)=iid
        OUTr0(www,ii)=r0(i)
        OUTrx(www,ii)=rx(i)
        OUTry(www,ii)=ry(i)
        OUTrz(www,ii)=rz(i)
        OUTp0(www,ii)=p0(i)
        OUTpx(www,ii)=px(i)
        OUTpy(www,ii)=py(i)
        OUTpz(www,ii)=pz(i)
        OUTmass(www,ii)=fmass(i)
        OUTityp(www,ii)=ityp(i)
        OUTiso3(www,ii)=iso3(i)
        OUTch(www,ii)=charge(i)
        OUTlcoll(www,ii)=lstcoll(i)
        OUTcoll(www,ii)=ncoll(i)
        OUTistr(www,ii)=istr
        OUTorigin(www,ii)=origin(i)
       enddo
      endif
       

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...

      write(iou(15),502) nin,nexit,iline,ctag,acttime,tsqrts
     ,     ,tstot,tsigpart,cdens
      do 11 i=1,nin
         istr=strit(tityp(i))
         ich = fchg(tiso3(i),tityp(i))
cbb iid is the particle slot id (standard) or the unique particle id (if
c   CTO(56) == 1).
         iid = tind(i)
         if( CTOption(56).eq.1 ) then
             iid = tuid(i)
         endif

         write(iou(15),501) iid,tr0(i),trx(i),try(i),trz(i),
     @                   tp0(i),tpx(i),tpy(i),tpz(i),tm(i),
     @                   tityp(i),tiso3(i),ich,tlcoll(i),
     @                   tcoll(i),istr,torigin(i)
 11   continue
      do 20 ii=1,nexit
         i=inew(ii)
         istr=strit(ityp(i))
         iid = i
         if( CTOption(56).eq.1 ) then
             iid = uid(i)
         endif
         write(iou(15),501) iid,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i)
 20   continue

      return



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15hyin(thydro_start)

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...
      nin=npart
      nexit=0

      write(iou(15),502) nin,nexit,91,0,thydro_start,0.d0
     ,     ,0.d0,0.d0,0.d0
      do 250 i=1,nin
         write(iou(15),501) i,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i)
 250  continue

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15hyout(thydro_start,thydro)

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...
      nin=0
      nexit=npart

      write(iou(15),502) nin,nexit,96,0,thydro_start+thydro,0.d0
     ,     ,0.d0,0.d0,0.d0
      do 251 i=1,nexit
       if(tform(i).lt.1.d-8)then
         write(iou(15),501) i,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i)
       else
         write(iou(15),501) i,tform(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i)
       end if
 251  continue

      return





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15outhy(thydro_start,thydro)

cbb Not sure what hydro_flag does and if it cannot be set somewhere
c   else. Nor what it really does. If it can be set somewhere else, I'd
c   prefer doing that; we could get rid of this entry altogether. For
c   the time being, though, we'll leave it in.
      if(thydro.gt.0.d0)then
       hydro_flag=.true.
      end if
      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f16outch

      if (bf16.or.(CTOption(13).eq.0)) return
      if (nin.eq.1) then
         tityp(2)=0
      endif

      do 22 ii=1,nexit
         i=inew(ii)
         write(iou(16),503) r0(i),rx(i),ry(i),rz(i),
     @        p0(i),px(i)+ffermpx(i),py(i)+ffermpy(i),
     @        pz(i)+ffermpz(i),fmass(i),
     @        ityp(i),iso3(i),charge(i),lstcoll(i),
     @        ncoll(i),mod(origin(i),100),tityp(1),tityp(2)
 22   continue

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16out

      echar='E'
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl

c
      if(bf16) return

c now write particle-output
      if (CToption(13).eq.0) then

         do 12 i=1,npart
            write(iou(16),201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           dectag+lstcoll(i),ncoll(i),mod(origin(i),100)

 12      continue
      else
         do 14 i=1,npart
            write(iou(16),503) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           dectag+lstcoll(i),ncoll(i),mod(origin(i),100),-99,-99
 14      continue
      endif
c
c write collision counters etc.
       write(iou(16),602) echar,itotcoll,NElColl,iinelcoll,NBlColl,
     @     dectag,NHardRes,NSoftRes,NDecRes
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16entry(ind)
c
c     This entry stores one decay for later output (must be done, in case
c     of pauli-blocked decay)
c
      tr0(3)=r0(ind)
      trx(3)=rx(ind)
      try(3)=ry(ind)
      trz(3)=rz(ind)
      tp0(3)=p0(ind)
      tpx(3)=px(ind)
      tpy(3)=py(ind)
      tpz(3)=pz(ind)
      tm(3)=fmass(ind)
      tityp(3)=ityp(ind)
      tind(3)=ind
      tiso3(3)=iso3(ind)
      tcharge(3)=charge(ind)

c     lstcoll is negative to identify decayed particles


      tlcoll(3)=-(1*lstcoll(ind))
      tcoll(3)=ncoll(ind)
      torigin(3)=origin(ind)
      tuid(3)=uid(ind)

      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16write

c     This entry writes the decay to file
      i=3

      if(bf16)return
      if (CTOption(13).eq.0) then

      write(iou(16),201) tr0(i),trx(i),try(i),trz(i),tp0(i),tpx(i),
     @        tpy(i),tpz(i),tm(i),tityp(i),tiso3(i),tcharge(i),
     @        tlcoll(i),tcoll(i),mod(torigin(i),100)
      else
      write(iou(16),503) tr0(i),trx(i),try(i),trz(i),tp0(i),tpx(i),
     @        tpy(i),tpz(i),tm(i),tityp(i),tiso3(i),tcharge(i),
     @        tlcoll(i),tcoll(i),mod(torigin(i),100),-98,-98
      endif

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_header

      if (bf19) return

chp initialize flag
      hydro_flag=.false.

      write (19,901) 'OSC1997A    '
      write (19,901) 'final_id_p_x'

 901  format (a12)

      if (CTOption(27).eq.0) then
         reffram='eqsp'
      elseif (CTOption(27).eq.1) then
         reffram='tar'
      elseif (CTOption(27).eq.2) then
         reffram='pro'
      else
         call error ('osc_header','Unknown Ref-Frame',
     .        dble(CTOption(27)),2)
         reffram='----'
      endif

      write (19,902) 'UrQMD', versiontxt, app, zpp, att, ztt,
     .     reffram, ebeam, 1

 902  format (2(a8,2x),'(',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,
     &     e10.4,2x,i8)

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_header

c header for OSCAR 99A output format

      if (bf20) return

      write (20,991)
      write (20,992)

 991  format ('# OSC1999A')
 992  format ('# full_event_history')

      if (CTOption(27).eq.0) then
         reffram='nncm'
      elseif (CTOption(27).eq.1) then
         reffram='tar'
      elseif (CTOption(27).eq.2) then
         reffram='pro'
      else
         call error ('osc_header','Unknown Ref-Frame',
     .        dble(CTOption(27)),2)
         reffram='----'
      endif

      write (20,993) versiontxt
 993  format ('# UrQMD ',a8)

      write (20,994) app, zpp, att, ztt,reffram, ebeam, 1

 994  format ('# (',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,
     &     e10.4,2x,i8)

      return



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_event

c body for OSCAR 97A format

      if (bf19) return


      if(CTOption(64).eq.0)write (19,903) event, npart, bimp, 0D0

      if(CTOption(64).eq.1)write (19,903) event,npartcoal, bimp, 0D0

 903  format (i10,2x,i10,2x,f8.3,2x,f8.3)


 904  format (i10,2x,i10,2x,9(e12.6,2x))


c particles, original
      iout=0
      do 99 i=1,npart
         if(lstcoll(i).ge.-1)then
         id = pdgid(ityp(i), iso3(i))
         iout=iout+1
         write(19,904) iout, id,
     .        px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i),
     .        p0(i), fmass(i),
     .        frrx(i), frry(i), frrz(i), frr0(i)
         endif
 99   continue

      if(CTOption(64).eq.1)then
       do i=1,nclus
         iout=iout+1
         id=itypcl(i)+99000
         write(19,904) iout,id,
     &        pxcl(i),pycl(i),pzcl(i),
     &        p0cl(i),mcl(i),
     &        rxcl(i),rycl(i),rzcl(i),r0cl(i)
       end do
      endif     


      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_vis(timestep)

c body for OSCAR 97A format adjusted for visualization

      if (bf19) return

      npart_form=0

       if(hydro_flag.eqv..true.)then
        do i=1,npart
         if(tform(i).le.acttime) then
           npart_form=npart_form+1
         end if
        end do
       end if

       if(hydro_flag.eqv..true.) then
           write (19,2301) event, npart_form, bimp, 0D0,nsteps,timestep
       else
           write(19,2301) event, npart, bimp,0D0,nsteps,timestep
       end if

 2301       FORMAT(I10,2X,I10,2X,F8.3,2X,F8.3,2x,i4,2x,i4)


chp modification for visualization output
      if(hydro_flag.eqv..true.)then
       do 990 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         if(tform(i).le.acttime)then
         write(19,2302) i, id,
     .       px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i),
     .       p0(i), fmass(i),
     .       rx(i), ry(i), rz(i), r0(i),tform(i),
     .       frrx(i), frry(i), frrz(i), frr0(i),ncoll(i)
         end if
 990    continue
      else
       do 980 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         write(19,2302) i, id,
     .       px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i),
     .       p0(i), fmass(i),
     .       rx(i), ry(i), rz(i), r0(i),tform(i),
     .       frrx(i), frry(i), frrz(i), frr0(i),ncoll(i)
 980    continue
      end if

 2302  FORMAT(I10,2X,I10,14(2X,E12.6),I10)



      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_event(ind)

c full event info for OSCAR 99A format

      if (bf20) return

      if(ind.eq.-1) then
         write (20,995) 0, npart, event, bimp, 0D0
      elseif(ind.eq.1) then
         write (20,996) npart, 0
      else
         write(6,*) 'fatal error in osc_99_event: wrong tag'
         stop 137
      endif

 995  format (3(i7,2x),2(f8.3,2x))
 996  format (2(i7,2x))

c particles

      do 88 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         write(20,997) uid(i), id, 0,
     .        px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i),
     .        p0(i), fmass(i),
     .        frrx(i), frry(i), frrz(i), frr0(i)
 88   continue

 997  format (3(i10,2x),9(e12.6,2x))

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_coll

      if (bf20) return
c     This entry writes information about the collision to file20
c     one line for each particle:
c     format: x y z px py pz ityp iso3...

      write(iou(20),999) nin,nexit,iline,ctag,acttime,tsqrts
     ,     ,tstot,tsigpart,cdens

      do 911 i=1,nin
         id = pdgid(tityp(i), tiso3(i))
         write(20,997) tuid(i), id, 0,
     .        tpx(i), tpy(i), tpz(i),tp0(i),tm(i),
     .        trx(i), try(i), trz(i), tr0(i)
 911   continue
      do 912 ii=1,nexit
         i=inew(ii)
         id = pdgid(ityp(i), iso3(i))
         write(20,997) uid(i), id, 0,
     .        px(i), py(i), pz(i),p0(i),fmass(i),
     .        rx(i), ry(i), rz(i), r0(i)
 912  continue


 999  format(3(i7,2x),i7,f8.3,4e12.4)
      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_eoe

c end of event tag for OSCAR 99A format
      if (bf20) return

      write(20,996) 0,0

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry getoldevent

c     read event header
      read(10,*,end=666)
     @ aa

      read(10,301) ac,pds, App, Zpp, ad,tds, Att, Ztt,add
      read(10,305) abt,betann,betatar,betapro
      read(10,304) ae,bimp,bmin,bdist,aee,sigmatot
      read(10,303) ah,eos,aj,ebeam,al,ecm,am,pbeam
      read(10,302) af,event,ag,ranseed,ag2,ai,ttime,ak,otime
C read unspecified number of lines with CTOptions (no more lines than
C possible; since 15 CTOs per line and no more than numcto options,
C there cannot be more than numcto/15+1 lines (note that numcto/15 is a
C truncated integer):
      do ii=0, numcto/15
C get next line:
        read(10,'(A)') line
C if that doesn't start with "op", we have reached the first
C non-CTO-line
        if(line(1:2).ne.'op') then
C ... and need to exit the loop.
          exit
        endif
        ! (else we are where we want to be and read Options)
        read(line,102) aop,(CTOtmp(i),CTOtc(i),i=ii*15+1,ii*15+15)
      enddo
C Same game for CTParams. Here, we inherit one line from the CTO code
C above.
      do ii=0, numctp/12
        if(line(1:2).ne.'pa') then
C doesn't start with 'pa'? Get outta here!
          exit
        endif
        ! (else we are where we want to be and read Parameters)
        read(line,103) apa,(CTPtmp(i),CTPtc(i),i=ii*12+1,ii*12+12)
C read next line for next iteration of loop.
        read(10,'(A)') line
      enddo
C if we got out of the loop, we might have a box header at hand.
      if (line(1:3).eq.'box') then
        boxflag = 1
        ! box header with length, energy, flags and number of particle
        ! specification lines to follow.
        read(line,504) abox1, lbox, abox2, edens, abox3, solid, abox4,
     1                 para, abox5, mbox
        read(10,505) aboxhead
        do l=1, mbox
         ! read from line
         read(10,506) abox6, bptityp(l), bptiso3(l),
     1                       bptpart(l), bptpmax(l)
        enddo
        ! read next line (which we would have if we didn't use up all
        ! our mojo for the box headers):
        read(10,'(A)') line
      endif
C Now, we expect the next line to be the "Particle Vector" line.
      read(line,306) apav

C set CTO only if CTOption(40) is 1 (standard behaviour) or if the
C asterisk '*' has been set (and, implicitly, CTOption(40) > 1).
C CTOption(40) will be overwritten in this loop, so we save the original
C value.
C Also, we set CTOdc(i) only if the above is true: We want to keep
C CTOdc(i) set from the inputfile.
      ctoforty=CTOption(40)
      ctofoone=CTOption(41)
      do i=1, numcto
        if( (ctoforty.eq.1.and.CTOption(i).ne.CTOtmp(i))
     &   .or.CTOtc(i).eq.' *') then
          write(*,*) 'Setting Option ',i,' (default: ',CTOption(i),
     &               ') to ',CTOtmp(i)
          CTOption(i) = CTOtmp(i)
          CTOdc(i) = CTOtc(i)
        endif
      enddo
C same as above, but for CTParam and CTPdc.
      do i=1, numctp
        if( (ctoforty.eq.1.and.CTParam(i).ne.CTPtmp(i))
     &   .or.CTPtc(i).eq.' *') then
          write(*,*) 'Setting Parameter ',i,' (default: ',CTParam(i),
     &               ') to ',CTPtmp(i)
          CTParam(i) = CTPtmp(i)
          CTPdc(i) = CTPtc(i)
        endif
      enddo
c check that CTO(41) is turned on. If not, issue a warning.
      if (CTOption(41).eq.0) then
        write(6,*) 'Warning: Reading in old event which seems to be ',
     +             'generated with CTO(41) = 0!'
      endif
c reset option 40
      CTOption(40)=ctoforty
c also reset option 41: Of course, it is turned on in the output we pipe
c in (or, really, we might try to take a look at that and warn if it
c isn't), but that shouldn't determine if we have the output of the
c continuation still in extended format or not; the current input file
c should do that instead.
      CTOption(41)=ctofoone

c read event body
      read(10,*) npart,ttime
      read(10,202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes
c      timestep=dble(ttime)/dtimestep
      ctag=itotcoll+dectag
c now read particle-output
      nbar=0
      do 39 i=1,npart
         read(10,210) r0(i),rx(i),ry(i),rz(i),p0(i),
     @        px(i),py(i),pz(i),fmass(i),
     @        ityp(i),iso3(i),charge(i),
     @        lstcoll(i),ncoll(i),origin(i),
     @        dectime(i),tform(i),xtotfac(i)
      if(abs(ityp(i)).le.maxbar)nbar=nbar+1
 39   continue
      nmes=npart-nbar
      acttime=r0(1)
c     read options-file
      call getparams
      success = .true.
      return
c handle EOF. This is not necessarily an error.
 666  write(*,*) 'No more events to read'
      success = .false.
      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry specout(ind,iu)
      i=ind
      if (CTOption(28).lt.0) return
      if (iu.eq.16.and.bf16) return
      if (iu.eq.14.and.bf14) return
      write(iu,201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @     px(i)+ffermpx(i),py(i)+ffermpy(i),
     @     pz(i)+ffermpz(i),fmass(i),
     @     ityp(i),iso3(i),charge(i),
     @     dectag+lstcoll(i),ncoll(i),mod(origin(i),100)

      return

      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spectrans(tstep)
c
c  (when cto 28 is set to 2 this subroutine is called
c  to propagate the spectators along straight lines)
c
cinput tstep : timestep
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 dtime,energ,tstep
      integer j
      include 'coms.inc'

      dtime=tstep

      do 1 j=1,nspec
         energ = p0s(j)
         r0s(j) = r0s(j) + dtime
         rxs(j) = rxs(j) + pxs(j)/energ*dtime
         rys(j) = rys(j) + pys(j)/energ*dtime
         rzs(j) = rzs(j) + pzs(j)/energ*dtime
1     continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function ptsigtot()
c
c     Revision : 1.0
c
c     This function caculates the total cross section of the reaction.
c     (Projectile - target total cross section)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'

      integer indmn,indmx,itypmn,iso3mn,itypmx,iso3mx
      integer isigline,iline,collclass
      real*8 stot,sigel
      real*8 sigtot

c determine total cross section for reaction:
      if(abs(Ap)+abs(At).gt.2) then
         stot=10.d0*pi*(bdist**2-bmin**2)
      elseif(abs(Ap)+abs(At).eq.2) then
         stot=sigtot(1,2,ecm)
cccccccc for CTOption(7)=1 no elastic cross section:
         if(CTOption(7).eq.1) then
c first sort the two itypes for call to collclass and anndec
            if(abs(ityp(1)).lt.abs(ityp(2))) then
               indmn=1
               indmx=2
            else
               indmn=2
               indmx=1
            endif

            itypmn=ityp(indmn)
            iso3mn=iso3(indmn)
            itypmx=ityp(indmx)
            iso3mx=iso3(indmx)
            isigline=collclass(itypmx,iso3mx,itypmn,iso3mn)
c     the elastic cross section is always the first entry (#3)
            iline=SigmaLn(3,1,isigline)
c!!!!DANGER: does not work for unstable particles (-> detailed balance)
            call crossx(iline,ecm,ityp(1),iso3(1),
     &              fmass(1),ityp(2),iso3(2),fmass(2),sigel)
c
            if(stot-sigel.gt.0) then
               stot=stot-sigel
            else
               stot=sigel
            endif
         endif
      else
         stot=0.d0
      endif
c
      ptsigtot=stot
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine urqmdlogo
c
c Displays the UrQMD Logo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'options.inc'
c     we need coms.f for versiontxt
      include 'coms.inc'

      integer firsttime
      save firsttime
      if (firsttime.eq.1)return
      firsttime=1

      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "## UrQMD ",versiontxt,
     $                    "   University of Frankfurt                ##"
      write (*,*)
     $   "##                  http://urqmd.org                       ##"
      write (*,*)
     $   "##                  bleicher@itp.uni-frankfurt.de          ##"
      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "##     Please cite when using this model:                  ##"
      write (*,*)
     $   "##     S.A.Bass et al., Prog.Part.Nucl.Phys. 41 (1998) 225 ##"
      write (*,*)
     $   "##     M.Bleicher et al., J.Phys. G25  (1999) 1859         ##"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##     UrQMD uses Pythia6.409 by T. Sjorstrand             ##"
      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "##     If hydrodynamic evolution is switched on (CTO 45 1) ##"
      write (*,*)
     $   "##     UrQMD uses the SHASTA algorithm by D. Rischke       ##"
      write (*,*)
     $   "##     Please cite when using the hybrid code:             ##"
      write (*,*)
     $   "##     D. Rischke et al., Nucl.Phys. A 595 (1995) 346      ##"
      write (*,*)
     $   "##     D. Rischke et al., Nucl.Phys. A 595 (1995) 383      ##"
      write (*,*)
     $   "##     H. Petersen et al., Phys.Rev. C78 (2008) 044901     ##"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "#############################################################"


      return
      end
c-------------------------------------------------------------------------

      subroutine makedeutfo(chg)

      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'

      integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
      integer lamb(1:1000),c1,c2,c3,casc(1:1000),cascmax,c4
      integer protmax,neutmax,lambmax,deutcl,deutn,nuclmax
      real*8 vcm(1:3)
      real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
      real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
      real*8 dP,dR,dR_max,dP_max
      real*8 ranf
      real*8 ndeut
      integer chg

      real*8 E9,pz9,y9,delt

      integer iou(13:20)
      data iou/13,14,15,16,17,18,19,20/
      common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax

      deutn=0


      do i2=1,protmax
         do i3=1,neutmax

c     boost momenta to local restframe

            pp1(0)=frp0(prot(i2))
            pp1(1)=frpx(prot(i2))
            pp1(2)=frpy(prot(i2))
            pp1(3)=frpz(prot(i2))

            pp2(0)=frp0(neut(i3))
            pp2(1)=frpx(neut(i3))
            pp2(2)=frpy(neut(i3))
            pp2(3)=frpz(neut(i3))


            xx1(0)=frr0(prot(i2))
            xx1(1)=frrx(prot(i2))
            xx1(2)=frry(prot(i2))
            xx1(3)=frrz(prot(i2))

            xx2(0)=frr0(neut(i3))
            xx2(1)=frrx(neut(i3))
            xx2(2)=frry(neut(i3))
            xx2(3)=frrz(neut(i3))

c     protons
            if (ncoll(prot(i2)).eq.0) then
               pp1(0)=p0(prot(i2))
               pp1(1)=px(prot(i2))+ffermpx(prot(i2))
               pp1(2)=py(prot(i2))+ffermpy(prot(i2))
               pp1(3)=pz(prot(i2))+ffermpz(prot(i2))
               xx1(0)=r0(prot(i2))
               xx1(1)=rx(prot(i2))
               xx1(2)=ry(prot(i2))
               xx1(3)=rz(prot(i2))
            end if
c     neutrons
            if (ncoll(neut(i3)).eq.0) then
               pp2(0)=p0(neut(i3))
               pp2(1)=px(neut(i3))+ffermpx(neut(i3))
               pp2(2)=py(neut(i3))+ffermpy(neut(i3))
               pp2(3)=pz(neut(i3))+ffermpz(neut(i3))
               xx2(0)=r0(neut(i3))
               xx2(1)=rx(neut(i3))
               xx2(2)=ry(neut(i3))
               xx2(3)=rz(neut(i3))
            end if

            vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
            vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
            vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

            call Lorentz(pp1,pp2,vcm,p1n,p2n)
            call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
            if (x2n(0).gt.x1n(0)) then
               delt=x2n(0)-x1n(0)
               do i1=0,3
                  x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
               end do
            else
c..   x1n .gt. x2n, propagate x2n
               delt=x1n(0)-x2n(0)
               do i1=0,3
                  x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
               end do
            end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

            call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

            dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
            dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

            dR_max=4.0
            dP_max=0.25
            if(eos.eq.1) dP_max=0.33


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

            if (ncoll(prot(i2)).eq.0)then
               c1=0
            else
               c1=1
            endif

            if (ncoll(neut(i3)).eq.0)then
               c2=0
            else
               c2=1
            endif

            deutcl=c1+c2

            if (lstcoll(prot(i2)).ge.-1.and.lstcoll(neut(i3)).ge.-1
     &           .and.ranf(0).le.0.375d0
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then
               lstcoll(prot(i2))=-2
               lstcoll(neut(i3))=-2

               if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
               xx12(0)=xx1(0)
               xx12(1)=(xx1(1)+xx2(1))/2d0
               xx12(2)=(xx1(2)+xx2(2))/2d0
               xx12(3)=(xx1(3)+xx2(3))/2d0
c     two-particle momentum
               pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &              (pp1(3)+pp2(3))**2
     &              +  (fmass(prot(i2))+fmass(neut(i3)))**2)
               pp12(1)=pp1(1)+pp2(1)
               pp12(2)=pp1(2)+pp2(2)
               pp12(3)=pp1(3)+pp2(3)

               deutn=deutn+1


                  nclus=nclus+1

                  mcl(nclus)=(fmass(prot(i2))+fmass(neut(i3)))
                  itypcl(nclus)=200*chg
                  chgcl(nclus)=1*chg
                  isocl(nclus)=0*chg
                  ncollcl(nclus)=deutcl
                  lstcollcl(nclus)=-2
                  r0cl(nclus)=xx12(0)
                  rxcl(nclus)=xx12(1)
                  rycl(nclus)=xx12(2)
                  rzcl(nclus)=xx12(3)
                  p0cl(nclus)=pp12(0)
                  pxcl(nclus)=pp12(1)
                  pycl(nclus)=pp12(2)
                  pzcl(nclus)=pp12(3)

            end if
         end do
      end do

      return

      end

c-------------------------------------------------------------------------

      subroutine maketritfo(chg)

      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'

      integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
      integer lamb(1:1000),c1,c2,c3,casc(1:1000),cascmax,c4
      integer protmax,neutmax,lambmax,deutcl,deutn,tricl,nuclmax
      real*8 vcm(1:3),vcm2(1:3)
      real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
      real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
      real*8 p1n2(0:3),p2n2(0:3),pp22(0:3)
      real*8 x1n2(0:3),x2n2(0:3),xx22(0:3)
      real*8 dP,dR,dR_max,dP_max,dP2,dR2,dR_max2,dP_max2
      real*8 ranf
      real*8 ndeut
      integer chg

      real*8 E9,pz9,y9,delt

      integer iou(13:20)
      data iou/13,14,15,16,17,18,19,20/
      common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax

      deutn=0

      do i2=1,protmax
         do i3=1,neutmax

c     boost momenta to local restframe

            pp1(0)=frp0(prot(i2))
            pp1(1)=frpx(prot(i2))
            pp1(2)=frpy(prot(i2))
            pp1(3)=frpz(prot(i2))

            pp2(0)=frp0(neut(i3))
            pp2(1)=frpx(neut(i3))
            pp2(2)=frpy(neut(i3))
            pp2(3)=frpz(neut(i3))


            xx1(0)=frr0(prot(i2))
            xx1(1)=frrx(prot(i2))
            xx1(2)=frry(prot(i2))
            xx1(3)=frrz(prot(i2))

            xx2(0)=frr0(neut(i3))
            xx2(1)=frrx(neut(i3))
            xx2(2)=frry(neut(i3))
            xx2(3)=frrz(neut(i3))

c     protons
            if (ncoll(prot(i2)).eq.0) then
               pp1(0)=p0(prot(i2))
               pp1(1)=px(prot(i2))+ffermpx(prot(i2))
               pp1(2)=py(prot(i2))+ffermpy(prot(i2))
               pp1(3)=pz(prot(i2))+ffermpz(prot(i2))
               xx1(0)=r0(prot(i2))
               xx1(1)=rx(prot(i2))
               xx1(2)=ry(prot(i2))
               xx1(3)=rz(prot(i2))
            end if
c     neutrons
            if (ncoll(neut(i3)).eq.0) then
               pp2(0)=p0(neut(i3))
               pp2(1)=px(neut(i3))+ffermpx(neut(i3))
               pp2(2)=py(neut(i3))+ffermpy(neut(i3))
               pp2(3)=pz(neut(i3))+ffermpz(neut(i3))
               xx2(0)=r0(neut(i3))
               xx2(1)=rx(neut(i3))
               xx2(2)=ry(neut(i3))
               xx2(3)=rz(neut(i3))
            end if

            vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
            vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
            vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

            call Lorentz(pp1,pp2,vcm,p1n,p2n)
            call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
            if (x2n(0).gt.x1n(0)) then
               delt=x2n(0)-x1n(0)
               do i1=0,3
                  x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
               end do
            else
c..   x1n .gt. x2n, propagate x2n
               delt=x1n(0)-x2n(0)
               do i1=0,3
                  x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
               end do
            end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

            call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

            dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
            dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

            dR_max=3.5
c            dP_max=0.32
            dP_max=0.4
            if(eos.eq.1) dP_max=0.45


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

            if (lstcoll(prot(i2)).ge.-1.and.lstcoll(neut(i3)).ge.-1
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then
c               write(*,*)nuclmax
               do i4=1,nuclmax
              if((nucl(i4).ne.neut(i3)).and.(nucl(i4).ne.prot(i2)))then
c     boost momenta to local restframe
                     if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
c     two-particle momentum
                     pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &                    (pp1(3)+pp2(3))**2
     &                    +  (fmass(prot(i2))+fmass(neut(i3)))**2)
                     pp12(1)=pp1(1)+pp2(1)
                     pp12(2)=pp1(2)+pp2(2)
                     pp12(3)=pp1(3)+pp2(3)

                     pp22(0)=frp0(nucl(i4))
                     pp22(1)=frpx(nucl(i4))
                     pp22(2)=frpy(nucl(i4))
                     pp22(3)=frpz(nucl(i4))


                     xx12(0)=xx1(0)
                     xx12(1)=(xx1(1)+xx2(1))/2d0
                     xx12(2)=(xx1(2)+xx2(2))/2d0
                     xx12(3)=(xx1(3)+xx2(3))/2d0

                     xx22(0)=frr0(nucl(i4))
                     xx22(1)=frrx(nucl(i4))
                     xx22(2)=frry(nucl(i4))
                     xx22(3)=frrz(nucl(i4))

                     if (ncoll(nucl(i4)).eq.0) then
                        pp22(0)=p0(nucl(i4))
                        pp22(1)=px(nucl(i4))+ffermpx(nucl(i4))
                        pp22(2)=py(nucl(i4))+ffermpy(nucl(i4))
                        pp22(3)=pz(nucl(i4))+ffermpz(nucl(i4))
                        xx22(0)=r0(nucl(i4))
                        xx22(1)=rx(nucl(i4))
                        xx22(2)=ry(nucl(i4))
                        xx22(3)=rz(nucl(i4))
                     end if

                     vcm2(1)=(pp12(1)+pp22(1))/(pp12(0)+pp22(0))
                     vcm2(2)=(pp12(2)+pp22(2))/(pp12(0)+pp22(0))
                     vcm2(3)=(pp12(3)+pp22(3))/(pp12(0)+pp22(0))

c..   Lorentzboost to 2-particle restframe
                     call Lorentz(pp12,pp22,vcm2,p1n2,p2n2)
                     call Lorentz(xx12,xx22,vcm2,x1n2,x2n2)
c     propagate to the later particle
                     if (x2n2(0).gt.x1n2(0)) then
                        delt=x2n2(0)-x1n2(0)
                        do i1=0,3
                           x1n2(i1)=x1n2(i1)+p1n2(i1)/p1n2(0)*delt
                        end do
                     else
c..   x1n .gt. x2n, propagate x2n
                        delt=x1n2(0)-x2n2(0)
                        do i1=0,3
                           x2n2(i1)=x2n2(i1)+p2n2(i1)/p2n2(0)*delt
                        end do
                     end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
                     call Lorentz(x1n2,x2n2,-vcm2,xx12,xx22)

c..
             dP2=sqrt((p1n2(1)-p2n2(1))**2+(p1n2(2)-p2n2(2))**2+(p1n2(3)
     &                    -p2n2(3))**2)
             dR2=sqrt((x1n2(1)-x2n2(1))**2+(x1n2(2)-x2n2(2))**2+(x1n2(3)
     &                    -x2n2(3))**2)

                     dR_max2=3.5
c                     dP_max2=0.32
                     dP_max2=0.4
                     if(eos.eq.1) dP_max2=0.45


                     if (ncoll(prot(i2)).eq.0)then
                        c1=0
                     else
                        c1=1
                     endif

                     if (ncoll(neut(i3)).eq.0)then
                        c2=0
                     else
                        c2=1
                     endif

                     if (ncoll(nucl(i4)).eq.0)then
                        c3=0
                     else
                        c3=1
                     endif

                     tricl=c1+c2+c3

                     if (lstcoll(prot(i2)).ge.-1
     $                    .and.lstcoll(neut(i3)).ge.-1
     &                    .and.dR2.lt.dR_max2.and.dP2.lt.dP_max2
     &                    .and.lstcoll(nucl(i4)).ge.-1
     $                    .and.ranf(0).le.0.0417d0
     $                    )then
                        lstcoll(prot(i2))=-3
                        lstcoll(neut(i3))=-3
                        lstcoll(nucl(i4))=-3

c...  production time of the 3He is the later time of neutron,deuteron
                        if (xx22(0).gt.xx12(0)) xx12(0)=xx22(0)


                  nclus=nclus+1

                  mcl(nclus)=(fmass(prot(i2))+fmass(neut(i3))+
     &     fmass(nucl(i4)))
                  ncollcl(nclus)=tricl
                  lstcollcl(nclus)=-3
                  r0cl(nclus)=xx12(0)
                  rxcl(nclus)=(2.*xx12(1)+xx22(1))/3d0
                  rycl(nclus)=(2.*xx12(2)+xx22(2))/3d0
                  rzcl(nclus)=(2.*xx12(3)+xx22(3))/3d0
                  p0cl(nclus)=sqrt((pp12(1)+pp22(1))**2
     &         +(pp12(2)+pp22(2))**2+(pp12(3)+pp22(3))**2
     &         +(fmass(prot(i2))+fmass(neut(i3))+fmass(nucl(i4)))**2)
                  pxcl(nclus)=(pp12(1)+pp22(1))
                  pycl(nclus)=(pp12(2)+pp22(2))
                  pzcl(nclus)=(pp12(3)+pp22(3))

                  if(charge(nucl(i4)).eq.0)then
                  itypcl(nclus)=201*chg
                  isocl(nclus)=-1*chg
                  chgcl(nclus)=1*chg
                  endif

                  if(charge(nucl(i4)).eq.1*chg)then
                  itypcl(nclus)=202*chg
                  isocl(nclus)=1*chg
                  chgcl(nclus)=2*chg
                  endif

                     end if
                  end if
               end do
            end if
         end do
      end do


      return


      end

c-------------------------------------------------------------------------

      subroutine makecasctrifo(chg)

      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'

      integer i1,i2,i3,i4,prot(1:1000),neut(1:1000)
      integer lamb(1:1000),c1,c2,c3,casc(1:1000),cascmax,c4
      integer protmax,neutmax,lambmax,deutcl,deutn,tricl
      integer nucl(1:1000),nuclmax,i5
      real*8 vcm(1:3),vcm2(1:3)
      real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
      real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
      real*8 p1n2(0:3),p2n2(0:3),pp22(0:3)
      real*8 x1n2(0:3),x2n2(0:3),xx22(0:3)
      real*8 dP,dR,dR_max,dP_max,dP2,dR2,dR_max2,dP_max2
      real*8 ranf
      real*8 ndeut
      integer chg

      real*8 E9,pz9,y9,delt

      integer iou(13:20)
      data iou/13,14,15,16,17,18,19,20/
      common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax

      deutn=0

      do i2=1,nuclmax
         do i3=1,nuclmax
            if(nucl(i2).ne.nucl(i3))then
c     boost momenta to local restframe

               pp1(0)=frp0(nucl(i2))
               pp1(1)=frpx(nucl(i2))
               pp1(2)=frpy(nucl(i2))
               pp1(3)=frpz(nucl(i2))

               pp2(0)=frp0(nucl(i3))
               pp2(1)=frpx(nucl(i3))
               pp2(2)=frpy(nucl(i3))
               pp2(3)=frpz(nucl(i3))


               xx1(0)=frr0(nucl(i2))
               xx1(1)=frrx(nucl(i2))
               xx1(2)=frry(nucl(i2))
               xx1(3)=frrz(nucl(i2))

               xx2(0)=frr0(nucl(i3))
               xx2(1)=frrx(nucl(i3))
               xx2(2)=frry(nucl(i3))
               xx2(3)=frrz(nucl(i3))

c     protons
               if (ncoll(nucl(i2)).eq.0) then
                  pp1(0)=p0(nucl(i2))
                  pp1(1)=px(nucl(i2))+ffermpx(nucl(i2))
                  pp1(2)=py(nucl(i2))+ffermpy(nucl(i2))
                  pp1(3)=pz(nucl(i2))+ffermpz(nucl(i2))
                  xx1(0)=r0(nucl(i2))
                  xx1(1)=rx(nucl(i2))
                  xx1(2)=ry(nucl(i2))
                  xx1(3)=rz(nucl(i2))
               end if
c     neutrons
               if (ncoll(nucl(i3)).eq.0) then
                  pp2(0)=p0(nucl(i3))
                  pp2(1)=px(nucl(i3))+ffermpx(nucl(i3))
                  pp2(2)=py(nucl(i3))+ffermpy(nucl(i3))
                  pp2(3)=pz(nucl(i3))+ffermpz(nucl(i3))
                  xx2(0)=r0(nucl(i3))
                  xx2(1)=rx(nucl(i3))
                  xx2(2)=ry(nucl(i3))
                  xx2(3)=rz(nucl(i3))
               end if

               vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
               vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
               vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

               call Lorentz(pp1,pp2,vcm,p1n,p2n)
               call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
               if (x2n(0).gt.x1n(0)) then
                  delt=x2n(0)-x1n(0)
                  do i1=0,3
                     x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
                  end do
               else
c..   x1n .gt. x2n, propagate x2n
                  delt=x1n(0)-x2n(0)
                  do i1=0,3
                     x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
                  end do
               end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

               call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

               dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &              -p2n(3))**2)
               dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &              -x2n(3))**2)

               dR_max=4.0
               dP_max=0.25
               if(eos.eq.1) dP_max=0.33


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

               if (lstcoll(nucl(i2)).ge.-1.and.lstcoll(nucl(i3)).ge.-1
     &              .and.dR.lt.dR_max.and.dP.lt.dP_max)then

                  do i4=1,cascmax
c     boost momenta to local restframe
                     if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
c     two-particle momentum
                     pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &                    (pp1(3)+pp2(3))**2
     &                    +  (fmass(nucl(i2))+fmass(nucl(i3)))**2)
                     pp12(1)=pp1(1)+pp2(1)
                     pp12(2)=pp1(2)+pp2(2)
                     pp12(3)=pp1(3)+pp2(3)

                     pp22(0)=frp0(casc(i4))
                     pp22(1)=frpx(casc(i4))
                     pp22(2)=frpy(casc(i4))
                     pp22(3)=frpz(casc(i4))


                     xx12(0)=xx1(0)
                     xx12(1)=(xx1(1)+xx2(1))/2d0
                     xx12(2)=(xx1(2)+xx2(2))/2d0
                     xx12(3)=(xx1(3)+xx2(3))/2d0

                     xx22(0)=frr0(casc(i4))
                     xx22(1)=frrx(casc(i4))
                     xx22(2)=frry(casc(i4))
                     xx22(3)=frrz(casc(i4))

                     if (ncoll(casc(i4)).eq.0) then
                        pp22(0)=p0(casc(i4))
                        pp22(1)=px(casc(i4))+ffermpx(casc(i4))
                        pp22(2)=py(casc(i4))+ffermpy(casc(i4))
                        pp22(3)=pz(casc(i4))+ffermpz(casc(i4))
                        xx22(0)=r0(casc(i4))
                        xx22(1)=rx(casc(i4))
                        xx22(2)=ry(casc(i4))
                        xx22(3)=rz(casc(i4))
                     end if

                     vcm2(1)=(pp12(1)+pp22(1))/(pp12(0)+pp22(0))
                     vcm2(2)=(pp12(2)+pp22(2))/(pp12(0)+pp22(0))
                     vcm2(3)=(pp12(3)+pp22(3))/(pp12(0)+pp22(0))

c..   Lorentzboost to 2-particle restframe
                     call Lorentz(pp12,pp22,vcm2,p1n2,p2n2)
                     call Lorentz(xx12,xx22,vcm2,x1n2,x2n2)
c     propagate to the later particle
                     if (x2n2(0).gt.x1n2(0)) then
                        delt=x2n2(0)-x1n2(0)
                        do i1=0,3
                           x1n2(i1)=x1n2(i1)+p1n2(i1)/p1n2(0)*delt
                        end do
                     else
c..   x1n .gt. x2n, propagate x2n
                        delt=x1n2(0)-x2n2(0)
                        do i1=0,3
                           x2n2(i1)=x2n2(i1)+p2n2(i1)/p2n2(0)*delt
                        end do
                     end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
                     call Lorentz(x1n2,x2n2,-vcm2,xx12,xx22)

c..
             dP2=sqrt((p1n2(1)-p2n2(1))**2+(p1n2(2)-p2n2(2))**2+(p1n2(3)
     &                    -p2n2(3))**2)
             dR2=sqrt((x1n2(1)-x2n2(1))**2+(x1n2(2)-x2n2(2))**2+(x1n2(3)
     &                    -x2n2(3))**2)

                     dR_max2=9.5
                     dP_max2=0.15
c
c
                     if (ncoll(nucl(i2)).eq.0)then
                        c1=0
                     else
                        c1=1
                     endif

                     if (ncoll(nucl(i3)).eq.0)then
                        c2=0
                     else
                        c2=1
                     endif

                     if (ncoll(casc(i4)).eq.0)then
                        c3=0
                     else
                        c3=1
                     endif

                     tricl=c1+c2+c3

                     if (lstcoll(nucl(i2)).ge.-1
     $                    .and.lstcoll(nucl(i3)).ge.-1
     &                    .and.dR2.lt.dR_max2.and.dP2.lt.dP_max2
     &                    .and.lstcoll(casc(i4)).ge.-1
     $                    .and.ranf(0).le.0.0417d0
     $                     )then
                        lstcoll(nucl(i2))=-3
                        lstcoll(nucl(i3))=-3
                        lstcoll(casc(i4))=-3



c...  production time of the triton is the later time of neutron,deuteron
                        if (xx22(0).gt.xx12(0)) xx12(0)=xx22(0)


                  nclus=nclus+1

                  mcl(nclus)=(fmass(nucl(i2))+fmass(nucl(i3))+
     &     fmass(casc(i4)))
                  itypcl(nclus)=222*chg
                  isocl(nclus)=0*chg
                  chgcl(nclus)=0*chg
                  ncollcl(nclus)=tricl
                  lstcollcl(nclus)=-3
                  r0cl(nclus)=xx12(0)
                  rxcl(nclus)=(2.*xx12(1)+xx22(1))/3d0
                  rycl(nclus)=(2.*xx12(2)+xx22(2))/3d0
                  rzcl(nclus)=(2.*xx12(3)+xx22(3))/3d0
                  p0cl(nclus)=sqrt((pp12(1)+pp22(1))**2
     &         +(pp12(2)+pp22(2))**2+(pp12(3)+pp22(3))**2
     &         +(fmass(nucl(i2))+fmass(nucl(i3))+fmass(casc(i4)))**2)
                  pxcl(nclus)=(pp12(1)+pp22(1))
                  pycl(nclus)=(pp12(2)+pp22(2))
                  pzcl(nclus)=(pp12(3)+pp22(3))



                     end if
                  end do
               end if
            endif
         end do
      end do


         return

      end


c..............................................................
c     boost
      subroutine Lorentz(p1,p2,vcm,p1n,p2n)
      implicit none
      real*8 p1(0:3),p2(0:3),p1n(0:3),p2n(0:3)
      real*8 vcm(1:3), g
      real*8 L(0:3,0:3)
      real*8 v(0:3)
      integer i,j

      v(0)=1./sqrt(1-(vcm(1)**2+vcm(2)**2+vcm(3)**2) )
      do i=1,3
         v(i)=vcm(i)
      end do

      L(0,0)=v(0)
      do i=1,3
         L(0,i)=-v(0)*v(i)
         L(i,0)=-v(0)*v(i)
      end do
      do i=1,3
         do j=1,3
            if (v(0).gt.1d-6) then
               L(i,j) =(v(0)-1.)*v(i)*v(j)/(v(1)**2+v(2)**2+v(3)**2)
            else
               L(i,j) = 0.5*v(i)*v(j)
            end if
         end do
      end do
      L(1,1)=L(1,1)+1.
      L(2,2)=L(2,2)+1.
      L(3,3)=L(3,3)+1.

      do i=0,3
         p1n(i)=0
         p2n(i)=0
      end do

      do i = 0,3
         do j = 0,3
            p1n(i)=p1n(i)+L(i,j)*p1(j)
            p2n(i)=p2n(i)+L(i,j)*p2(j)
         end do
      end do

      return

      end

c-------------------------------------------------------------------------

      subroutine makehyptritfo(chg)

      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'

      integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
      integer lamb(1:1000),c1,c2,c3,casc(1:1000),cascmax,c4
      integer protmax,neutmax,lambmax,deutcl,deutn,tricl,nuclmax
      real*8 vcm(1:3),vcm2(1:3)
      real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
      real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
      real*8 p1n2(0:3),p2n2(0:3),pp22(0:3)
      real*8 x1n2(0:3),x2n2(0:3),xx22(0:3)
      real*8 dP,dR,dR_max,dP_max,dP2,dR2,dR_max2,dP_max2
      real*8 ranf
      real*8 ndeut
      integer chg

      real*8 E9,pz9,y9,delt

      integer iou(13:20)
      data iou/13,14,15,16,17,18,19,20/
      common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax

      deutn=0


      do i2=1,protmax
         do i3=1,neutmax

c     boost momenta to local restframe

            pp1(0)=frp0(prot(i2))
            pp1(1)=frpx(prot(i2))
            pp1(2)=frpy(prot(i2))
            pp1(3)=frpz(prot(i2))

            pp2(0)=frp0(neut(i3))
            pp2(1)=frpx(neut(i3))
            pp2(2)=frpy(neut(i3))
            pp2(3)=frpz(neut(i3))


            xx1(0)=frr0(prot(i2))
            xx1(1)=frrx(prot(i2))
            xx1(2)=frry(prot(i2))
            xx1(3)=frrz(prot(i2))

            xx2(0)=frr0(neut(i3))
            xx2(1)=frrx(neut(i3))
            xx2(2)=frry(neut(i3))
            xx2(3)=frrz(neut(i3))

c     protons
            if (ncoll(prot(i2)).eq.0) then
               pp1(0)=p0(prot(i2))
               pp1(1)=px(prot(i2))+ffermpx(prot(i2))
               pp1(2)=py(prot(i2))+ffermpy(prot(i2))
               pp1(3)=pz(prot(i2))+ffermpz(prot(i2))
               xx1(0)=r0(prot(i2))
               xx1(1)=rx(prot(i2))
               xx1(2)=ry(prot(i2))
               xx1(3)=rz(prot(i2))
            end if
c     neutrons
            if (ncoll(neut(i3)).eq.0) then
               pp2(0)=p0(neut(i3))
               pp2(1)=px(neut(i3))+ffermpx(neut(i3))
               pp2(2)=py(neut(i3))+ffermpy(neut(i3))
               pp2(3)=pz(neut(i3))+ffermpz(neut(i3))
               xx2(0)=r0(neut(i3))
               xx2(1)=rx(neut(i3))
               xx2(2)=ry(neut(i3))
               xx2(3)=rz(neut(i3))
            end if

            vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
            vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
            vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

            call Lorentz(pp1,pp2,vcm,p1n,p2n)
            call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
            if (x2n(0).gt.x1n(0)) then
               delt=x2n(0)-x1n(0)
               do i1=0,3
                  x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
               end do
            else
c..   x1n .gt. x2n, propagate x2n
               delt=x1n(0)-x2n(0)
               do i1=0,3
                  x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
               end do
            end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

            call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

            dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
            dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

 
            dR_max=4.0
            dP_max=0.25
            if(eos.eq.1) dP_max=0.33


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

            if (lstcoll(prot(i2)).ge.-1.and.lstcoll(neut(i3)).ge.-1
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then

               do i4=1,lambmax
                  if (i4.ne.i3) then
c     boost momenta to local restframe
                     if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
c     two-particle momentum
                     pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &                    (pp1(3)+pp2(3))**2
     &                    +  (fmass(prot(i2))+fmass(neut(i3)))**2)
                     pp12(1)=pp1(1)+pp2(1)
                     pp12(2)=pp1(2)+pp2(2)
                     pp12(3)=pp1(3)+pp2(3)

                     pp22(0)=frp0(lamb(i4))
                     pp22(1)=frpx(lamb(i4))
                     pp22(2)=frpy(lamb(i4))
                     pp22(3)=frpz(lamb(i4))


                     xx12(0)=xx1(0)
                     xx12(1)=(xx1(1)+xx2(1))/2d0
                     xx12(2)=(xx1(2)+xx2(2))/2d0
                     xx12(3)=(xx1(3)+xx2(3))/2d0

                     xx22(0)=frr0(lamb(i4))
                     xx22(1)=frrx(lamb(i4))
                     xx22(2)=frry(lamb(i4))
                     xx22(3)=frrz(lamb(i4))

                     vcm2(1)=(pp12(1)+pp22(1))/(pp12(0)+pp22(0))
                     vcm2(2)=(pp12(2)+pp22(2))/(pp12(0)+pp22(0))
                     vcm2(3)=(pp12(3)+pp22(3))/(pp12(0)+pp22(0))

c..   Lorentzboost to 2-particle restframe
                     call Lorentz(pp12,pp22,vcm2,p1n2,p2n2)
                     call Lorentz(xx12,xx22,vcm2,x1n2,x2n2)
c     propagate to the later particle
                     if (x2n2(0).gt.x1n2(0)) then
                        delt=x2n2(0)-x1n2(0)
                        do i1=0,3
                           x1n2(i1)=x1n2(i1)+p1n2(i1)/p1n2(0)*delt
                        end do
                     else
c..   x1n .gt. x2n, propagate x2n
                        delt=x1n2(0)-x2n2(0)
                        do i1=0,3
                           x2n2(i1)=x2n2(i1)+p2n2(i1)/p2n2(0)*delt
                        end do
                     end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
                     call Lorentz(x1n2,x2n2,-vcm2,xx12,xx22)

c..
             dP2=sqrt((p1n2(1)-p2n2(1))**2+(p1n2(2)-p2n2(2))**2+(p1n2(3)
     &                    -p2n2(3))**2)
             dR2=sqrt((x1n2(1)-x2n2(1))**2+(x1n2(2)-x2n2(2))**2+(x1n2(3)
     &                    -x2n2(3))**2)

                     dR_max2=9.5
                     dP_max2=0.15

                     if(eos.eq.1)dP_max2=0.18

                     if (ncoll(prot(i2)).eq.0)then
                        c1=0
                     else
                        c1=1
                     endif

                     if (ncoll(neut(i3)).eq.0)then
                        c2=0
                     else
                        c2=1
                     endif

                     if (ncoll(lamb(i4)).eq.0)then
                        c3=0
                     else
                        c3=1
                     endif

                     tricl=c1+c2+c3
c                     if(dR2.lt.dR_max2.and.dP2.lt.dP_max2)then
c                     write(*,*)i2,i3,i4,dR2,dP2,tricl
c                     endif

                     if (lstcoll(prot(i2)).ge.-1
     $                    .and.lstcoll(neut(i3)).ge.-1
     &                    .and.dR2.lt.dR_max2.and.dP2.lt.dP_max2
     &                    .and.lstcoll(lamb(i4)).ge.-1
     $                    .and.ranf(0).le.0.0833d0
     $                    )then
                        lstcoll(prot(i2))=-3
                        lstcoll(neut(i3))=-3
                        lstcoll(lamb(i4))=-3

c                       write(*,*)"hypertriton!"


c...  production time of the triton is the later time of neutron,deuteron
                        if (xx22(0).gt.xx12(0)) xx12(0)=xx22(0)


                  nclus=nclus+1

                  mcl(nclus)=(fmass(prot(i2))+fmass(neut(i3))+
     &     fmass(lamb(i4)))
                  itypcl(nclus)=211*chg
                  chgcl(nclus)=1*chg
                  isocl(nclus)=0*chg
                  ncollcl(nclus)=tricl
                  lstcollcl(nclus)=-3
                  r0cl(nclus)=xx12(0)
                  rxcl(nclus)=(2.*xx12(1)+xx22(1))/3d0
                  rycl(nclus)=(2.*xx12(2)+xx22(2))/3d0
                  rzcl(nclus)=(2.*xx12(3)+xx22(3))/3d0
                  p0cl(nclus)=sqrt((pp12(1)+pp22(1))**2
     &         +(pp12(2)+pp22(2))**2+(pp12(3)+pp22(3))**2
     &         +(fmass(prot(i2))+fmass(neut(i3))+fmass(lamb(i4)))**2)
                  pxcl(nclus)=(pp12(1)+pp22(1))
                  pycl(nclus)=(pp12(2)+pp22(2))
                  pzcl(nclus)=(pp12(3)+pp22(3))

                     end if
                  end if
               end do
            end if
         end do
      end do


      return

      end
c-------------------------------------------------------------------------

      subroutine makecascdifo(chg)

      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'

      integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
      integer lamb(1:1000),c1,c2,c3,casc(1:1000),cascmax,c4
      integer protmax,neutmax,lambmax,deutcl,deutn,nuclmax
      real*8 vcm(1:3)
      real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
      real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
      real*8 dP,dR,dR_max,dP_max
      real*8 ranf
      real*8 ndeut
      integer chg

      real*8 E9,pz9,y9,delt

      integer iou(13:20)
      data iou/13,14,15,16,17,18,19,20/
      common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax

      deutn=0

      do i2=1,cascmax
         do i3=1,nuclmax

c     boost momenta to local restframe

            pp1(0)=frp0(casc(i2))
            pp1(1)=frpx(casc(i2))
            pp1(2)=frpy(casc(i2))
            pp1(3)=frpz(casc(i2))

            pp2(0)=frp0(nucl(i3))
            pp2(1)=frpx(nucl(i3))
            pp2(2)=frpy(nucl(i3))
            pp2(3)=frpz(nucl(i3))


            xx1(0)=frr0(casc(i2))
            xx1(1)=frrx(casc(i2))
            xx1(2)=frry(casc(i2))
            xx1(3)=frrz(casc(i2))

            xx2(0)=frr0(nucl(i3))
            xx2(1)=frrx(nucl(i3))
            xx2(2)=frry(nucl(i3))
            xx2(3)=frrz(nucl(i3))

c     protons
            if (ncoll(casc(i2)).eq.0) then
               pp1(0)=p0(casc(i2))
               pp1(1)=px(casc(i2))+ffermpx(casc(i2))
               pp1(2)=py(casc(i2))+ffermpy(casc(i2))
               pp1(3)=pz(casc(i2))+ffermpz(casc(i2))
               xx1(0)=r0(casc(i2))
               xx1(1)=rx(casc(i2))
               xx1(2)=ry(casc(i2))
               xx1(3)=rz(casc(i2))
            end if
c     neutrons
            if (ncoll(nucl(i3)).eq.0) then
               pp2(0)=p0(nucl(i3))
               pp2(1)=px(nucl(i3))+ffermpx(nucl(i3))
               pp2(2)=py(nucl(i3))+ffermpy(nucl(i3))
               pp2(3)=pz(nucl(i3))+ffermpz(nucl(i3))
               xx2(0)=r0(nucl(i3))
               xx2(1)=rx(nucl(i3))
               xx2(2)=ry(nucl(i3))
               xx2(3)=rz(nucl(i3))
            end if

            vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
            vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
            vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

            call Lorentz(pp1,pp2,vcm,p1n,p2n)
            call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
            if (x2n(0).gt.x1n(0)) then
               delt=x2n(0)-x1n(0)
               do i1=0,3
                  x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
               end do
            else
c..   x1n .gt. x2n, propagate x2n
               delt=x1n(0)-x2n(0)
               do i1=0,3
                  x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
               end do
            end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

            call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

            dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
            dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

            dR_max=9.5
            dP_max=0.15


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

            if (ncoll(casc(i2)).eq.0)then
               c1=0
            else
               c1=1
            endif

            if (ncoll(nucl(i3)).eq.0)then
               c2=0
            else
               c2=1
            endif

            deutcl=c1+c2

            if (lstcoll(casc(i2)).ge.-1.and.lstcoll(nucl(i3)).ge.-1
     &           .and.ranf(0).le.0.375d0
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then
               lstcoll(casc(i2))=-2
               lstcoll(nucl(i3))=-2

               if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
               xx12(0)=xx1(0)
               xx12(1)=(xx1(1)+xx2(1))/2d0
               xx12(2)=(xx1(2)+xx2(2))/2d0
               xx12(3)=(xx1(3)+xx2(3))/2d0
c     two-particle momentum
               pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &              (pp1(3)+pp2(3))**2
     &              +  (fmass(casc(i2))+fmass(nucl(i3)))**2)
               pp12(1)=pp1(1)+pp2(1)
               pp12(2)=pp1(2)+pp2(2)
               pp12(3)=pp1(3)+pp2(3)

               deutn=deutn+1

                  nclus=nclus+1

                  mcl(nclus)=(fmass(casc(i2))+fmass(nucl(i3)))
                  itypcl(nclus)=221*chg
                  chgcl(nclus)=0*chg
                  isocl(nclus)=0*chg
                  ncollcl(nclus)=deutcl
                  lstcollcl(nclus)=-2
                  r0cl(nclus)=xx12(0)
                  rxcl(nclus)=xx12(1)
                  rycl(nclus)=xx12(2)
                  rzcl(nclus)=xx12(3)
                  p0cl(nclus)=pp12(0)
                  pxcl(nclus)=pp12(1)
                  pycl(nclus)=pp12(2)
                  pzcl(nclus)=pp12(3)


            

            end if
         end do
      end do

      return

      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine makehdibfo(chg)

      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'

      integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
      integer lamb(1:1000),c1,c2,c3,casc(1:1000),cascmax,c4
      integer protmax,neutmax,lambmax,deutcl,deutn,nuclmax
      real*8 vcm(1:3)
      real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
      real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
      real*8 dP,dR,dR_max,dP_max
      real*8 ranf
      real*8 ndeut
      integer chg

      real*8 E9,pz9,y9,delt

      integer iou(13:20)
      data iou/13,14,15,16,17,18,19,20/
      common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax


      deutn=0

      do i2=1,lambmax
         do i3=1,lambmax
            if(i2.ne.i3)then
c     boost momenta to local restframe

            pp1(0)=frp0(lamb(i2))
            pp1(1)=frpx(lamb(i2))
            pp1(2)=frpy(lamb(i2))
            pp1(3)=frpz(lamb(i2))

            pp2(0)=frp0(lamb(i3))
            pp2(1)=frpx(lamb(i3))
            pp2(2)=frpy(lamb(i3))
            pp2(3)=frpz(lamb(i3))


            xx1(0)=frr0(lamb(i2))
            xx1(1)=frrx(lamb(i2))
            xx1(2)=frry(lamb(i2))
            xx1(3)=frrz(lamb(i2))

            xx2(0)=frr0(lamb(i3))
            xx2(1)=frrx(lamb(i3))
            xx2(2)=frry(lamb(i3))
            xx2(3)=frrz(lamb(i3))

c     protons
            if (ncoll(lamb(i2)).eq.0) then
               pp1(0)=p0(lamb(i2))
               pp1(1)=px(lamb(i2))+ffermpx(lamb(i2))
               pp1(2)=py(lamb(i2))+ffermpy(lamb(i2))
               pp1(3)=pz(lamb(i2))+ffermpz(lamb(i2))
               xx1(0)=r0(lamb(i2))
               xx1(1)=rx(lamb(i2))
               xx1(2)=ry(lamb(i2))
               xx1(3)=rz(lamb(i2))
            end if
c     neutrons
            if (ncoll(nucl(i3)).eq.0) then
               pp2(0)=p0(lamb(i3))
               pp2(1)=px(lamb(i3))+ffermpx(lamb(i3))
               pp2(2)=py(lamb(i3))+ffermpy(lamb(i3))
               pp2(3)=pz(lamb(i3))+ffermpz(lamb(i3))
               xx2(0)=r0(lamb(i3))
               xx2(1)=rx(lamb(i3))
               xx2(2)=ry(lamb(i3))
               xx2(3)=rz(lamb(i3))
            end if

            vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
            vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
            vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

            call Lorentz(pp1,pp2,vcm,p1n,p2n)
            call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
            if (x2n(0).gt.x1n(0)) then
               delt=x2n(0)-x1n(0)
               do i1=0,3
                  x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
               end do
            else
c..   x1n .gt. x2n, propagate x2n
               delt=x1n(0)-x2n(0)
               do i1=0,3
                  x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
               end do
            end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

            call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

            dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
            dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

            dR_max=9.5
            dP_max=0.15


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

            if (ncoll(lamb(i2)).eq.0)then
               c1=0
            else
               c1=1
            endif

            if (ncoll(lamb(i3)).eq.0)then
               c2=0
            else
               c2=1
            endif

            deutcl=c1+c2

            if (lstcoll(lamb(i2)).ge.-1.and.lstcoll(lamb(i3)).ge.-1
     &           .and.ranf(0).le.0.1875d0
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then
               lstcoll(lamb(i2))=-2
               lstcoll(lamb(i3))=-2

               if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
               xx12(0)=xx1(0)
               xx12(1)=(xx1(1)+xx2(1))/2d0
               xx12(2)=(xx1(2)+xx2(2))/2d0
               xx12(3)=(xx1(3)+xx2(3))/2d0
c     two-particle momentum
               pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &              (pp1(3)+pp2(3))**2
     &              +  (fmass(lamb(i2))+fmass(lamb(i3)))**2)
               pp12(1)=pp1(1)+pp2(1)
               pp12(2)=pp1(2)+pp2(2)
               pp12(3)=pp1(3)+pp2(3)

               deutn=deutn+1

                  nclus=nclus+1

                  mcl(nclus)=(fmass(lamb(i2))+fmass(lamb(i3)))
                  itypcl(nclus)=220*chg
                  chgcl(nclus)=0*chg
                  isocl(nclus)=0*chg
                  ncollcl(nclus)=deutcl
                  lstcollcl(nclus)=-2
                  r0cl(nclus)=xx12(0)
                  rxcl(nclus)=xx12(1)
                  rycl(nclus)=xx12(2)
                  rzcl(nclus)=xx12(3)
                  p0cl(nclus)=pp12(0)
                  pxcl(nclus)=pp12(1)
                  pycl(nclus)=pp12(2)
                  pzcl(nclus)=pp12(3)

            end if
         endif
        end do
      end do

      return

      end


c-------------------------------------------------------------------------

            subroutine makehefo(chg)

            implicit none

            include 'comres.inc'
            include 'coms.inc'
            include 'options.inc'
            include 'inputs.inc'
            include 'newpart.inc'
            include 'freezeout.inc'
            include 'boxinc.inc'

          integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
            integer lamb(1:1000),c1,c2,c3,c4,casc(1:1000),cascmax
            integer protmax,neutmax,lambmax,deutcl,hen,tricl,nuclmax
            integer hecl
            real*8 vcm(1:3),vcm2(1:3)
            real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
            real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
            real*8 p1n2(0:3),p2n2(0:3),pp22(0:3)
            real*8 x1n2(0:3),x2n2(0:3),xx22(0:3)
            real*8 p1n3(0:3),p2n3(0:3),pp23(0:3),pp13(0:3)
            real*8 x1n3(0:3),x2n3(0:3),xx23(0:3),xx13(0:3)
            real*8 dP,dR,dR_max,dP_max,dP2,dR2,dR_max2,dP_max2
            real*8 dR_max3,dP_max3
            real*8 ranf
            real*8 nhe
            integer chg

            real*8 E9,pz9,y9,delt

            integer iou(13:20)
            data iou/13,14,15,16,17,18,19,20/
            common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax


c            hen=0

            do i2=1,protmax
               do i3=1,neutmax

c     boost momenta to local restframe

                  pp1(0)=frp0(prot(i2))
                  pp1(1)=frpx(prot(i2))
                  pp1(2)=frpy(prot(i2))
                  pp1(3)=frpz(prot(i2))

                  pp2(0)=frp0(neut(i3))
                  pp2(1)=frpx(neut(i3))
                  pp2(2)=frpy(neut(i3))
                  pp2(3)=frpz(neut(i3))


                  xx1(0)=frr0(prot(i2))
                  xx1(1)=frrx(prot(i2))
                  xx1(2)=frry(prot(i2))
                  xx1(3)=frrz(prot(i2))

                  xx2(0)=frr0(neut(i3))
                  xx2(1)=frrx(neut(i3))
                  xx2(2)=frry(neut(i3))
                  xx2(3)=frrz(neut(i3))

c     protons
                  if (ncoll(prot(i2)).eq.0) then
                     pp1(0)=p0(prot(i2))
                     pp1(1)=px(prot(i2))+ffermpx(prot(i2))
                     pp1(2)=py(prot(i2))+ffermpy(prot(i2))
                     pp1(3)=pz(prot(i2))+ffermpz(prot(i2))
                     xx1(0)=r0(prot(i2))
                     xx1(1)=rx(prot(i2))
                     xx1(2)=ry(prot(i2))
                     xx1(3)=rz(prot(i2))
                  end if
c     neutrons
                  if (ncoll(neut(i3)).eq.0) then
                     pp2(0)=p0(neut(i3))
                     pp2(1)=px(neut(i3))+ffermpx(neut(i3))
                     pp2(2)=py(neut(i3))+ffermpy(neut(i3))
                     pp2(3)=pz(neut(i3))+ffermpz(neut(i3))
                     xx2(0)=r0(neut(i3))
                     xx2(1)=rx(neut(i3))
                     xx2(2)=ry(neut(i3))
                     xx2(3)=rz(neut(i3))
                  end if

                  vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
                  vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
                  vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

                  call Lorentz(pp1,pp2,vcm,p1n,p2n)
                  call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
                  if (x2n(0).gt.x1n(0)) then
                     delt=x2n(0)-x1n(0)
                     do i1=0,3
                        x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
                     end do
                  else
c..   x1n .gt. x2n, propagate x2n
                     delt=x1n(0)-x2n(0)
                     do i1=0,3
                        x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
                     end do
                  end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

                  call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

          dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
          dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

                  dR_max=3.5
                  dP_max=0.41

                  if(eos.eq.1) dP_max=0.55


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

       if (lstcoll(prot(i2)).ge.-1.and.lstcoll(neut(i3)).ge.-1
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then
c               write(*,*)nuclmax
             do i4=1,protmax
            if(prot(i4).ne.prot(i2))then
c     boost momenta to local restframe
                           if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
c     two-particle momentum
             pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &                    (pp1(3)+pp2(3))**2
     &                    +  (fmass(prot(i2))+fmass(neut(i3)))**2)
                         pp12(1)=pp1(1)+pp2(1)
                         pp12(2)=pp1(2)+pp2(2)
                         pp12(3)=pp1(3)+pp2(3)

                         pp22(0)=frp0(prot(i4))
                         pp22(1)=frpx(prot(i4))
                         pp22(2)=frpy(prot(i4))
                         pp22(3)=frpz(prot(i4))


                         xx12(0)=xx1(0)
                         xx12(1)=(xx1(1)+xx2(1))/2d0
                         xx12(2)=(xx1(2)+xx2(2))/2d0
                         xx12(3)=(xx1(3)+xx2(3))/2d0

                         xx22(0)=frr0(prot(i4))
                         xx22(1)=frrx(prot(i4))
                         xx22(2)=frry(prot(i4))
                         xx22(3)=frrz(prot(i4))

                       if (ncoll(prot(i4)).eq.0) then
                          pp22(0)=p0(prot(i4))
                          pp22(1)=px(prot(i4))+ffermpx(prot(i4))
                          pp22(2)=py(prot(i4))+ffermpy(prot(i4))
                          pp22(3)=pz(prot(i4))+ffermpz(prot(i4))
                          xx22(0)=r0(prot(i4))
                          xx22(1)=rx(prot(i4))
                          xx22(2)=ry(prot(i4))
                          xx22(3)=rz(prot(i4))
                       end if

                       vcm2(1)=(pp12(1)+pp22(1))/(pp12(0)+pp22(0))
                       vcm2(2)=(pp12(2)+pp22(2))/(pp12(0)+pp22(0))
                       vcm2(3)=(pp12(3)+pp22(3))/(pp12(0)+pp22(0))

c..   Lorentzboost to 2-particle restframe
                       call Lorentz(pp12,pp22,vcm2,p1n2,p2n2)
                       call Lorentz(xx12,xx22,vcm2,x1n2,x2n2)
c     propagate to the later particle
                       if (x2n2(0).gt.x1n2(0)) then
                          delt=x2n2(0)-x1n2(0)
                          do i1=0,3
                             x1n2(i1)=x1n2(i1)+p1n2(i1)/p1n2(0)*delt
                          end do
                       else
c..   x1n .gt. x2n, propagate x2n
                          delt=x1n2(0)-x2n2(0)
                          do i1=0,3
                             x2n2(i1)=x2n2(i1)+p2n2(i1)/p2n2(0)*delt
                          end do
                       end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
                       call Lorentz(x1n2,x2n2,-vcm2,xx12,xx22)

c..
             dP2=sqrt((p1n2(1)-p2n2(1))**2+(p1n2(2)-p2n2(2))**2+(p1n2(3)
     &                    -p2n2(3))**2)
             dR2=sqrt((x1n2(1)-x2n2(1))**2+(x1n2(2)-x2n2(2))**2+(x1n2(3)
     &                    -x2n2(3))**2)

                     dR_max2=3.5
                     dP_max2=0.41
                     if(eos.eq.1)dP_max2=0.55

                     if (lstcoll(prot(i2)).ge.-1
     $                    .and.lstcoll(neut(i3)).ge.-1
     &                    .and.dR2.lt.dR_max2.and.dP2.lt.dP_max2
     &                    .and.lstcoll(prot(i4)).ge.-1
     $                    )then


cJS add 4th particle
c               write(*,*)nuclmax
       do i5=1,neutmax
       if(neut(i5).ne.neut(i3))then
c     boost momenta to local restframe
         if (xx22(0).gt.xx12(0)) xx12(0)=xx22(0)
c     two-particle momentum
        pp13(0)=sqrt((pp12(1)+pp22(1))**2+(pp12(2)+pp22(2))**2+
     &                    (pp12(3)+pp22(3))**2
     &   +  (fmass(prot(i2))+fmass(neut(i3))+fmass(prot(i4)))**2)
         pp13(1)=pp12(1)+pp22(1)
         pp13(2)=pp12(2)+pp22(2)
         pp13(3)=pp12(3)+pp22(3)

         pp23(0)=frp0(neut(i5))
         pp23(1)=frpx(neut(i5))
         pp23(2)=frpy(neut(i5))
         pp23(3)=frpz(neut(i5))


         xx13(0)=xx12(0)
         xx13(1)=(xx12(1)+xx22(1))/2d0
         xx13(2)=(xx12(2)+xx22(2))/2d0
         xx13(3)=(xx12(3)+xx22(3))/2d0

         xx23(0)=frr0(neut(i5))
         xx23(1)=frrx(neut(i5))
         xx23(2)=frry(neut(i5))
         xx23(3)=frrz(neut(i5))

         if (ncoll(neut(i5)).eq.0) then
            pp23(0)=p0(neut(i5))
            pp23(1)=px(neut(i5))+ffermpx(neut(i5))
            pp23(2)=py(neut(i5))+ffermpy(neut(i5))
            pp23(3)=pz(neut(i5))+ffermpz(neut(i5))
            xx23(0)=r0(neut(i5))
            xx23(1)=rx(neut(i5))
            xx23(2)=ry(neut(i5))
            xx23(3)=rz(neut(i5))
         end if

         vcm2(1)=(pp13(1)+pp23(1))/(pp13(0)+pp23(0))
         vcm2(2)=(pp13(2)+pp23(2))/(pp13(0)+pp23(0))
         vcm2(3)=(pp13(3)+pp23(3))/(pp13(0)+pp23(0))

c..   Lorentzboost to 2-particle restframe
         call Lorentz(pp13,pp23,vcm2,p1n3,p2n3)
         call Lorentz(xx13,xx23,vcm2,x1n3,x2n3)
c     propagate to the later particle
         if (x2n3(0).gt.x1n3(0)) then
            delt=x2n3(0)-x1n3(0)
            do i1=0,3
               x1n3(i1)=x1n3(i1)+p1n3(i1)/p1n3(0)*delt
            end do
         else
c..   x1n .gt. x2n, propagate x2n
            delt=x1n3(0)-x2n3(0)
            do i1=0,3
               x2n3(i1)=x2n3(i1)+p2n3(i1)/p2n3(0)*delt
            end do
         end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
         call Lorentz(x1n3,x2n3,-vcm2,xx13,xx23)

c..
       dP2=sqrt((p1n3(1)-p2n3(1))**2+(p1n3(2)-p2n3(2))**2+(p1n3(3)
     &                   -p2n3(3))**2)
       dR2=sqrt((x1n3(1)-x2n3(1))**2+(x1n3(2)-x2n3(2))**2+(x1n3(3)
     &                    -x2n3(3))**2)

         dR_max3=3.5
         dP_max3=0.41
         if(eos.eq.1) dP_max3=0.55

         if (lstcoll(prot(i2)).ge.-1
     $                    .and.lstcoll(neut(i3)).ge.-1
     &                .and.dR2.lt.dR_max3.and.dP2.lt.dP_max3
     &                    .and.lstcoll(prot(i4)).ge.-1
     &                    .and.lstcoll(neut(i5)).ge.-1
     $                    .and.ranf(0).le.0.002604166d0
     $                    )then


        if (ncoll(prot(i2)).eq.0)then
          c1=0
        else
          c1=1
        endif

         if (ncoll(neut(i3)).eq.0)then
          c2=0
          else
          c2=1
          endif

          if (ncoll(prot(i4)).eq.0)then
              c3=0
          else
             c3=1
           endif

          if (ncoll(neut(i5)).eq.0)then
              c4=0
           else
           c4=1
          endif

          hecl=c1+c2+c3+c4


              lstcoll(prot(i2))=-4
              lstcoll(neut(i3))=-4
              lstcoll(prot(i4))=-4
              lstcoll(neut(i5))=-4


c...  production time of the helium is the later time of neutron, triton
              if (xx23(0).gt.xx13(0)) xx13(0)=xx23(0)

                  nclus=nclus+1

                  mcl(nclus)=(fmass(prot(i2))+fmass(neut(i3))+
     &     fmass(prot(i4))+fmass(neut(i5)))
                  itypcl(nclus)=203*chg
                  chgcl(nclus)=2*chg
                  isocl(nclus)=0*chg
                  ncollcl(nclus)=hecl
                  lstcollcl(nclus)=-4
                  r0cl(nclus)=xx13(0)
                  rxcl(nclus)=(3.*xx13(1)+xx23(1))/4d0
                  rycl(nclus)=(3.*xx13(2)+xx23(2))/4d0
                  rzcl(nclus)=(3.*xx13(3)+xx23(3))/4d0
                  p0cl(nclus)=sqrt((pp13(1)+pp23(1))**2
     &         +(pp13(2)+pp23(2))**2+(pp13(3)+pp23(3))**2
     &         +(fmass(prot(i2))+fmass(neut(i3))+fmass(prot(i4))+
     &         fmass(neut(i5)))**2)
                  pxcl(nclus)=(pp13(1)+pp23(1))
                  pycl(nclus)=(pp13(2)+pp23(2))
                  pzcl(nclus)=(pp13(3)+pp23(3))


                       end if
                       end if
                       end do
                      end if
                      end if
                     end do
                  end if
               end do
            end do


            return

c     same including freeze-out coordinates
  213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8)
22213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8,i6,4i6)

            end
c-------------------------------------------------------------------------

c-------------------------------------------------------------------------

            subroutine makehyphyfo(chg)

            implicit none

            include 'comres.inc'
            include 'coms.inc'
            include 'options.inc'
            include 'inputs.inc'
            include 'newpart.inc'
            include 'freezeout.inc'
            include 'boxinc.inc'

          integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
            integer lamb(1:1000),c1,c2,c3,c4,casc(1:1000),cascmax
            integer protmax,neutmax,lambmax,deutcl,hen,tricl,nuclmax
            integer hyphycl
            real*8 vcm(1:3),vcm2(1:3)
            real*8 p1n(0:3),p2n(0:3),pp1(0:3),pp2(0:3),pp12(0:3)
            real*8 x1n(0:3),x2n(0:3),xx1(0:3),xx2(0:3),xx12(0:3)
            real*8 p1n2(0:3),p2n2(0:3),pp22(0:3)
            real*8 x1n2(0:3),x2n2(0:3),xx22(0:3)
            real*8 p1n3(0:3),p2n3(0:3),pp23(0:3),pp13(0:3)
            real*8 x1n3(0:3),x2n3(0:3),xx23(0:3),xx13(0:3)
            real*8 dP,dR,dR_max,dP_max,dP2,dR2,dR_max2,dP_max2
            real*8 dR_max3,dP_max3
            real*8 ranf
            real*8 nhyphy
            integer chg

            real*8 E9,pz9,y9,delt

            integer iou(13:20)
            data iou/13,14,15,16,17,18,19,20/
            common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     $     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax


            nhyphy=0

            do i2=1,protmax
               do i3=1,neutmax

c     boost momenta to local restframe

                  pp1(0)=frp0(prot(i2))
                  pp1(1)=frpx(prot(i2))
                  pp1(2)=frpy(prot(i2))
                  pp1(3)=frpz(prot(i2))

                  pp2(0)=frp0(neut(i3))
                  pp2(1)=frpx(neut(i3))
                  pp2(2)=frpy(neut(i3))
                  pp2(3)=frpz(neut(i3))


                  xx1(0)=frr0(prot(i2))
                  xx1(1)=frrx(prot(i2))
                  xx1(2)=frry(prot(i2))
                  xx1(3)=frrz(prot(i2))

                  xx2(0)=frr0(neut(i3))
                  xx2(1)=frrx(neut(i3))
                  xx2(2)=frry(neut(i3))
                  xx2(3)=frrz(neut(i3))

c     protons
                  if (ncoll(prot(i2)).eq.0) then
                     pp1(0)=p0(prot(i2))
                     pp1(1)=px(prot(i2))+ffermpx(prot(i2))
                     pp1(2)=py(prot(i2))+ffermpy(prot(i2))
                     pp1(3)=pz(prot(i2))+ffermpz(prot(i2))
                     xx1(0)=r0(prot(i2))
                     xx1(1)=rx(prot(i2))
                     xx1(2)=ry(prot(i2))
                     xx1(3)=rz(prot(i2))
                  end if
c     neutrons
                  if (ncoll(neut(i3)).eq.0) then
                     pp2(0)=p0(neut(i3))
                     pp2(1)=px(neut(i3))+ffermpx(neut(i3))
                     pp2(2)=py(neut(i3))+ffermpy(neut(i3))
                     pp2(3)=pz(neut(i3))+ffermpz(neut(i3))
                     xx2(0)=r0(neut(i3))
                     xx2(1)=rx(neut(i3))
                     xx2(2)=ry(neut(i3))
                     xx2(3)=rz(neut(i3))
                  end if

                  vcm(1)=(pp1(1)+pp2(1))/(pp1(0)+pp2(0))
                  vcm(2)=(pp1(2)+pp2(2))/(pp1(0)+pp2(0))
                  vcm(3)=(pp1(3)+pp2(3))/(pp1(0)+pp2(0))

                  call Lorentz(pp1,pp2,vcm,p1n,p2n)
                  call Lorentz(xx1,xx2,vcm,x1n,x2n)

c     propagate to the later particle
                  if (x2n(0).gt.x1n(0)) then
                     delt=x2n(0)-x1n(0)
                     do i1=0,3
                        x1n(i1)=x1n(i1)+p1n(i1)/p1n(0)*delt
                     end do
                  else
c..   x1n .gt. x2n, propagate x2n
                     delt=x1n(0)-x2n(0)
                     do i1=0,3
                        x2n(i1)=x2n(i1)+p2n(i1)/p2n(0)*delt
                     end do
                  end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)

                  call Lorentz(x1n,x2n,-vcm,xx1,xx2)
c..

          dP=sqrt((p1n(1)-p2n(1))**2+(p1n(2)-p2n(2))**2+(p1n(3)
     &           -p2n(3))**2)
          dR=sqrt((x1n(1)-x2n(1))**2+(x1n(2)-x2n(2))**2+(x1n(3)
     &           -x2n(3))**2)

                  dR_max=3.5
                  dP_max=0.41

                  if(eos.eq.1) dP_max=0.55


c..   include spin and isospin projection stochastically (factor 1/2*3/4)
c...  lstcoll is set to -2 for protons,neutrons that are now in a deuteron
c..   deuteron ityp is set to ityp=200

       if (lstcoll(prot(i2)).ge.-1.and.lstcoll(neut(i3)).ge.-1
     &           .and.dR.lt.dR_max.and.dP.lt.dP_max)then
c               write(*,*)nuclmax
             do i4=1,neutmax
            if(neut(i4).ne.neut(i3))then
c     boost momenta to local restframe
                           if (xx2(0).gt.xx1(0)) xx1(0)=xx2(0)
c     two-particle momentum
             pp12(0)=sqrt((pp1(1)+pp2(1))**2+(pp1(2)+pp2(2))**2+
     &                    (pp1(3)+pp2(3))**2
     &                    +  (fmass(prot(i2))+fmass(neut(i3)))**2)
                         pp12(1)=pp1(1)+pp2(1)
                         pp12(2)=pp1(2)+pp2(2)
                         pp12(3)=pp1(3)+pp2(3)

                         pp22(0)=frp0(neut(i4))
                         pp22(1)=frpx(neut(i4))
                         pp22(2)=frpy(neut(i4))
                         pp22(3)=frpz(neut(i4))


                         xx12(0)=xx1(0)
                         xx12(1)=(xx1(1)+xx2(1))/2d0
                         xx12(2)=(xx1(2)+xx2(2))/2d0
                         xx12(3)=(xx1(3)+xx2(3))/2d0

                         xx22(0)=frr0(neut(i4))
                         xx22(1)=frrx(neut(i4))
                         xx22(2)=frry(neut(i4))
                         xx22(3)=frrz(neut(i4))

                       if (ncoll(neut(i4)).eq.0) then
                          pp22(0)=p0(neut(i4))
                          pp22(1)=px(neut(i4))+ffermpx(neut(i4))
                          pp22(2)=py(neut(i4))+ffermpy(neut(i4))
                          pp22(3)=pz(neut(i4))+ffermpz(neut(i4))
                          xx22(0)=r0(neut(i4))
                          xx22(1)=rx(neut(i4))
                          xx22(2)=ry(neut(i4))
                          xx22(3)=rz(neut(i4))
                       end if

                       vcm2(1)=(pp12(1)+pp22(1))/(pp12(0)+pp22(0))
                       vcm2(2)=(pp12(2)+pp22(2))/(pp12(0)+pp22(0))
                       vcm2(3)=(pp12(3)+pp22(3))/(pp12(0)+pp22(0))

c..   Lorentzboost to 2-particle restframe
                       call Lorentz(pp12,pp22,vcm2,p1n2,p2n2)
                       call Lorentz(xx12,xx22,vcm2,x1n2,x2n2)
c     propagate to the later particle
                       if (x2n2(0).gt.x1n2(0)) then
                          delt=x2n2(0)-x1n2(0)
                          do i1=0,3
                             x1n2(i1)=x1n2(i1)+p1n2(i1)/p1n2(0)*delt
                          end do
                       else
c..   x1n .gt. x2n, propagate x2n
                          delt=x1n2(0)-x2n2(0)
                          do i1=0,3
                             x2n2(i1)=x2n2(i1)+p2n2(i1)/p2n2(0)*delt
                          end do
                       end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
                       call Lorentz(x1n2,x2n2,-vcm2,xx12,xx22)

c..
             dP2=sqrt((p1n2(1)-p2n2(1))**2+(p1n2(2)-p2n2(2))**2+(p1n2(3)
     &                    -p2n2(3))**2)
             dR2=sqrt((x1n2(1)-x2n2(1))**2+(x1n2(2)-x2n2(2))**2+(x1n2(3)
     &                    -x2n2(3))**2)

                     dR_max2=3.5
                     dP_max2=0.41
                     if(eos.eq.1)dP_max2=0.55

                     if (lstcoll(prot(i2)).ge.-1
     $                    .and.lstcoll(neut(i3)).ge.-1
     &                    .and.dR2.lt.dR_max2.and.dP2.lt.dP_max2
     &                    .and.lstcoll(neut(i4)).ge.-1
     $                    )then


cJS add 4th particle
c               write(*,*)nuclmax
       do i5=1,lambmax
c       if(neut(i5).ne.neut(i3))then
c     boost momenta to local restframe
         if (xx22(0).gt.xx12(0)) xx12(0)=xx22(0)
c     two-particle momentum
        pp13(0)=sqrt((pp12(1)+pp22(1))**2+(pp12(2)+pp22(2))**2+
     &                    (pp12(3)+pp22(3))**2
     &   +  (fmass(prot(i2))+fmass(neut(i3))+fmass(neut(i4)))**2)
         pp13(1)=pp12(1)+pp22(1)
         pp13(2)=pp12(2)+pp22(2)
         pp13(3)=pp12(3)+pp22(3)

         pp23(0)=frp0(lamb(i5))
         pp23(1)=frpx(lamb(i5))
         pp23(2)=frpy(lamb(i5))
         pp23(3)=frpz(lamb(i5))


         xx13(0)=xx12(0)
         xx13(1)=(xx12(1)+xx22(1))/2d0
         xx13(2)=(xx12(2)+xx22(2))/2d0
         xx13(3)=(xx12(3)+xx22(3))/2d0

         xx23(0)=frr0(lamb(i5))
         xx23(1)=frrx(lamb(i5))
         xx23(2)=frry(lamb(i5))
         xx23(3)=frrz(lamb(i5))

         if (ncoll(lamb(i5)).eq.0) then
            pp23(0)=p0(lamb(i5))
            pp23(1)=px(lamb(i5))+ffermpx(lamb(i5))
            pp23(2)=py(lamb(i5))+ffermpy(lamb(i5))
            pp23(3)=pz(lamb(i5))+ffermpz(lamb(i5))
            xx23(0)=r0(lamb(i5))
            xx23(1)=rx(lamb(i5))
            xx23(2)=ry(lamb(i5))
            xx23(3)=rz(lamb(i5))
         end if

         vcm2(1)=(pp13(1)+pp23(1))/(pp13(0)+pp23(0))
         vcm2(2)=(pp13(2)+pp23(2))/(pp13(0)+pp23(0))
         vcm2(3)=(pp13(3)+pp23(3))/(pp13(0)+pp23(0))

c..   Lorentzboost to 2-particle restframe
         call Lorentz(pp13,pp23,vcm2,p1n3,p2n3)
         call Lorentz(xx13,xx23,vcm2,x1n3,x2n3)
c     propagate to the later particle
         if (x2n3(0).gt.x1n3(0)) then
            delt=x2n3(0)-x1n3(0)
            do i1=0,3
               x1n3(i1)=x1n3(i1)+p1n3(i1)/p1n3(0)*delt
            end do
         else
c..   x1n .gt. x2n, propagate x2n
            delt=x1n3(0)-x2n3(0)
            do i1=0,3
               x2n3(i1)=x2n3(i1)+p2n3(i1)/p2n3(0)*delt
            end do
         end if
c..   boost the equal-time proton,neutron freeze-out 4-coordinates back
c..   to the computational frame (use old variables xx1,xx2 again)
         call Lorentz(x1n3,x2n3,-vcm2,xx13,xx23)

c..
       dP2=sqrt((p1n3(1)-p2n3(1))**2+(p1n3(2)-p2n3(2))**2+(p1n3(3)
     &                   -p2n3(3))**2)
       dR2=sqrt((x1n3(1)-x2n3(1))**2+(x1n3(2)-x2n3(2))**2+(x1n3(3)
     &                    -x2n3(3))**2)

            dR_max3=9.5
            dP_max3=0.2
            if(eos.eq.1)dP_max3=0.25



         if (lstcoll(prot(i2)).ge.-1
     $                    .and.lstcoll(neut(i3)).ge.-1
     &                .and.dR2.lt.dR_max3.and.dP2.lt.dP_max3
     &                    .and.lstcoll(neut(i4)).ge.-1
     &                    .and.lstcoll(lamb(i5)).ge.-1
     $                    .and.ranf(0).le.0.00520833d0
     $                    )then


        if (ncoll(prot(i2)).eq.0)then
          c1=0
        else
          c1=1
        endif

         if (ncoll(neut(i3)).eq.0)then
          c2=0
          else
          c2=1
          endif

          if (ncoll(neut(i4)).eq.0)then
              c3=0
          else
             c3=1
           endif

          if (ncoll(lamb(i5)).eq.0)then
              c4=0
           else
           c4=1
          endif

          hyphycl=c1+c2+c3+c4


              lstcoll(prot(i2))=-4
              lstcoll(neut(i3))=-4
              lstcoll(neut(i4))=-4
              lstcoll(lamb(i5))=-4


c...  production time of the hyperhydro4 is the later time of Lambda, triton
              if (xx23(0).gt.xx13(0)) xx13(0)=xx23(0)

                  nclus=nclus+1

                  mcl(nclus)=(fmass(prot(i2))+fmass(neut(i3))+
     &     fmass(neut(i4))+fmass(lamb(i5)))
                  itypcl(nclus)=213*chg
                  chgcl(nclus)=1*chg
                  isocl(nclus)=-1*chg
                  ncollcl(nclus)=hyphycl
                  lstcollcl(nclus)=-4
                  r0cl(nclus)=xx13(0)
                  rxcl(nclus)=(3.*xx13(1)+xx23(1))/4d0
                  rycl(nclus)=(3.*xx13(2)+xx23(2))/4d0
                  rzcl(nclus)=(3.*xx13(3)+xx23(3))/4d0
                  p0cl(nclus)=sqrt((pp13(1)+pp23(1))**2
     &         +(pp13(2)+pp23(2))**2+(pp13(3)+pp23(3))**2
     &         +(fmass(prot(i2))+fmass(neut(i3))+fmass(neut(i4))+
     &         fmass(lamb(i5)))**2)
                  pxcl(nclus)=(pp13(1)+pp23(1))
                  pycl(nclus)=(pp13(2)+pp23(2))
                  pzcl(nclus)=(pp13(3)+pp23(3))
                 
                 end if
                
                       end do
                      end if
                      end if
                     end do
                  end if
               end do
            end do


            return

c     same including freeze-out coordinates
  213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8)
22213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8,i6,4i6)

            end
c-------------------------------------------------------------------------


c-------------------------------------------------------------------------
ctr This subroutine determines the reconstructable resonances by
ctr following each daughter particles of every decaying resonance
ctr throug the whole collision history. A resonance is treated as
ctr reconstructable by an experiment if both (all) daughter particles
ctr escape from the system without further rescattering.
c
      subroutine reconstruct_resonances()

      implicit none
      include 'comhis.inc'

      integer j,k,l,m,n
      integer j0,k0
      integer maxchain
      parameter (maxchain=20)
      integer j00(maxchain)
      integer rtot,r,q

      rtot=0
      r=0
      do q=1,maxchain
       j00(q)=0
      enddo

c.. loop over all interactions www
      do j=1,www
       j0=j
       rtot=0
       r=0
       do q=1,maxchain
        j00(q)=0
       enddo
 91    continue
c.. resonance decay, NOT pauli blocked
       if (hisiline(j0).eq.20.and.hisnexit(j0).ne.0) then
c.. loop over daughter particles
        do k=1,hisnexit(j0)
c.. loop over future from resonances perspective
         do l=j0+1,www
c.. loop over incoming particles
          do m=1,hisnin(l)
c.. check if daughters interact again, NOT pauli blocked
           if(OUTind(j0,k).eq.INind(l,m)
     &        .and.hisnexit(l).ne.0) then
c.. if daughter is resonance, it can still be reconstructed
c.. if the daughter resonance can also be reconstructed
c.. go to beginning of algorith with shifted starting point
            if(hisiline(l).eq.20)then
             rtot=rtot+1
             j00(rtot)=l
             goto 80 ! skip goto 77
            endif
            goto 77
 80         continue
           endif
          enddo
         enddo
        enddo
        if(SUM(j00).ne.0)then
         r=r+1
         if(r.gt.rtot)then
          goto 92
         endif
         j0=j00(r)
         goto 91
        endif
 92     continue
c.. if this point is reached, all daughters have escaped
c.. without interaction, store information
        resctr=resctr+1
        RECr0(resctr)=INr0(j,1)
        RECrx(resctr)=INrx(j,1)
        RECry(resctr)=INry(j,1)
        RECrz(resctr)=INrz(j,1)
        RECp0(resctr)=INp0(j,1)
        RECpx(resctr)=INpx(j,1)
        RECpy(resctr)=INpy(j,1)
        RECpz(resctr)=INpz(j,1)
        RECmass(resctr)=INmass(j,1)
        RECityp(resctr)=INityp(j,1)
        RECiso3(resctr)=INiso3(j,1)
        RECch(resctr)=INch(j,1)
        REClcoll(resctr)=INlcoll(j,1)
        RECcoll(resctr)=INcoll(j,1)
        RECorigin(resctr)=INorigin(j,1)
        do k=0,4
         RECind(resctr,k)=0
        enddo
        RECind(resctr,0)=INind(j,1)
        do k=1,hisnexit(j)
         RECind(resctr,k)=OUTind(j,k)
        enddo
        do l=1,j-1
         do m=1,hisnexit(j-l)
          if(OUTind(j-l,m).eq.INind(j,1))then
           ORIr0(resctr)=OUTr0(j-l,m)
           ORIrx(resctr)=OUTrx(j-l,m)
           ORIry(resctr)=OUTry(j-l,m)
           ORIrz(resctr)=OUTrz(j-l,m)
           ORIp0(resctr)=OUTp0(j-l,m)
           ORIpx(resctr)=OUTpx(j-l,m)
           ORIpy(resctr)=OUTpy(j-l,m)
           ORIpz(resctr)=OUTpz(j-l,m)
           ORImass(resctr)=OUTmass(j-l,m)
           ORIityp(resctr)=OUTityp(j-l,m)
           ORIiso3(resctr)=OUTiso3(j-l,m)
           ORIch(resctr)=OUTch(j-l,m)
           ORIlcoll(resctr)=OUTlcoll(j-l,m)
           ORIcoll(resctr)=OUTcoll(j-l,m)
           ORIorigin(resctr)=OUTorigin(j-l,m)
           goto 77
          endif
         enddo
        enddo
       endif
 77    continue
      enddo

      return
      end
c-------------------------------------------------------------------------


      subroutine coalescence


      implicit none

      include 'comres.inc'
      include 'coms.inc'
      include 'options.inc'
      include 'inputs.inc'
      include 'newpart.inc'
      include 'freezeout.inc'
      include 'boxinc.inc'
      include 'comhis.inc'

      integer i1,i2,i3,i4,i5,prot(1:1000),neut(1:1000),nucl(1:1000)
      integer lamb(1:1000),c1,c2,c3,i6,casc(1:1000),c4
      integer protmax,neutmax,lambmax,deutcl,deutn,nuclmax,cascmax


       common /barlistn/i1,i2,i3,i4,prot,neut,lamb,c1,c2,c3,c4
     &     ,protmax,neutmax,lambmax,nuclmax,i5,nucl,casc,cascmax

      nclus=0

c     JS create list of baryons
      i2=1
      i3=1
      i4=1
      i5=1
      i6=1
      do i1 = 1,npart
c     generate list of protons and neutrons and lambda
         if(lstcoll(i1).lt.0)then
            lstcoll(i1)=0
         endif
         if (ityp(i1).eq.1.and.charge(i1).eq.1) then
            prot(i2)=i1
            i2=i2+1
         end if
         if (ityp(i1).eq.1.and.charge(i1).eq.0) then
            neut(i3)=i1
            i3=i3+1
         end if
         if (ityp(i1).eq.41.and.charge(i1).eq.0) then
            lamb(i4)=i1
            i4=i4+1
         endif
         if (ityp(i1).eq.1) then
            nucl(i5)=i1
            i5=i5+1
         end if
         if (ityp(i1).eq.63)then
            casc(i6)=i1
            i6=i6+1
         endif
      end do
      protmax=i2-1
      neutmax=i3-1
      lambmax=i4-1
      nuclmax=i5-1
      cascmax=i6-1
c     JS first start with heavy nuclei
      
      if(protmax.ge.1.and.neutmax.ge.2.and.lambmax.ge.1)then
         call makehyphyfo(1)
      endif
      if(protmax.ge.2.and.neutmax.ge.2)then
         call makehefo(1)
      endif
      if(cascmax.ge.1)then
         call makecasctrifo(1)
         call makecascdifo(1)
      endif
      if(lambmax.ge.2) then
         call makehdibfo(1)
      endif
      if(lambmax.ge.1) then
        call makehyptritfo(1)
      endif
      if(protmax.ge.1.and.neutmax.ge.1.and.nuclmax.ge.3)then
        call maketritfo(1)
      endif
      if(protmax.ge.1.and.neutmax.ge.1)then
        call makedeutfo(1)
      endif
cJS Now do anti-nuclei
c     JS create list of anti-baryons
      i2=1
      i3=1
      i4=1
      i5=1
      i6=1
      do i1 = 1,npart
c     generate list of anti-protons and -neutrons and -lambda
         if (ityp(i1).eq.-1.and.charge(i1).eq.-1) then
            prot(i2)=i1
            i2=i2+1
         end if
         if (ityp(i1).eq.-1.and.charge(i1).eq.0) then
            neut(i3)=i1
            i3=i3+1
         end if
         if (ityp(i1).eq.-41.and.charge(i1).eq.0) then
            lamb(i4)=i1
            i4=i4+1
         endif
         if (ityp(i1).eq.-1) then
            nucl(i5)=i1
            i5=i5+1
         end if
         if (ityp(i1).eq.-63)then
            casc(i6)=i1
            i6=i6+1
         endif
      end do
      protmax=i2-1
      neutmax=i3-1
      lambmax=i4-1
      nuclmax=i5-1
      cascmax=i6-1
c     JS first start with heavy anti-nuclei


      if(protmax.ge.1.and.neutmax.ge.2.and.lambmax.ge.1)then
         call makehyphyfo(-1)
      endif
      if(protmax.ge.2.and.neutmax.ge.2)then
         call makehefo(-1)
      endif
      if(cascmax.ge.1)then
         call makecasctrifo(-1)
         call makecascdifo(-1)
      endif
      if(lambmax.ge.2) then
         call makehdibfo(-1)
      endif
      if(lambmax.ge.1) then
        call makehyptritfo(-1)
      endif
      if(protmax.ge.1.and.neutmax.ge.1.and.nuclmax.ge.3)then
        call maketritfo(-1)
      endif
      if(protmax.ge.1.and.neutmax.ge.1)then
        call makedeutfo(-1)
      endif

      npartcoal=0

      do i1 =1,npart
         if(lstcoll(i1).ge.-1)then
                 npartcoal=npartcoal+1
         endif 
      end do

      npartcoal=npartcoal+nclus

      return
      end


