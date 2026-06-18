
c Setting of global paramters
c
      subroutine params
      implicit none
      real*8 A0, chi
      include 'coms.inc'
      include 'options.inc'
c     gw = 0.25 fm^-2 width of the gaussian

      logSky = .true.
      logYuk = .true.
      logCb  = .true.
      logPau = .false.

      gw     = 0.25
      sgw    = sqrt(gw)
      Cb0    = 1.44
      Yuk0   = 0.0 !-85.0
      gamYuk = 1.4
      drPau  = 9.0
      dpPau  = 0.0144
      Pau0   = 0.0 !99.5*(hqc/sqrt(drPau*dpPau))**3

C
C hard Skyrme EOS (usual stuff)
C
c      Sky30  = 70.5 
c      gamSky = 2.0 
c      A0     = -124.2 * 0.5
C

C hard Skyrme EOS (JK parametrisation corrected for Gausswidth gw=0.25)
C
       Sky30  = 125.93
       gamSky = 1.676
       A0     = -87.67
       chi    = 0.93 

C soft Skyrme EoS

c       Sky30  = 303.0
c       gamSky = 7.0/6.0
c       A0     = -356.0*0.5
c       chi    = 0.93

      Sky20 = chi*2.0*A0
      Yuk0 = (1.0-chi)/(2.0*pi*gamYuk**2)*A0

      delr = 0.2
      fdel = delr*delr/6.0
      da = -(1.0/delr)
      db = -da

      cutdww = 20.0
      cutPau = 20.0
      cutYuk = 20.0
      cutCb  = 20.0

      dtimestep=0.2     
      dt  = 0.02
c      dt2 = 0.5*dt
c      dt6 = dt/6.0

      return
      end


c Reset of all indexed variables
c
      subroutine set0
      implicit none
      integer i, j
      include 'coms.inc'

      do 10 i=1,nspl
        spPauy(i) = 0.0 
        outPau(i) = 0.0 
        spCby(i)  = 0.0  
        outCb(i)  = 0.0
        spYuky(i) = 0.0 
        outYuk(i) = 0.0
        spSkyy(i) = 0.0 
        outSky(i) = 0.0
        spdwwy(i) = 0.0 
        outdww(i) = 0.0
  10  continue

      do 20 j=1,nmax
        spin(j)  = 0
        iso3(j)   = 0
        ncoll(j) = 0
        rx(j)    = 0.0
        ry(j)    = 0.0
        rz(j)    = 0.0
        px(j)    = 0.0
        py(j)    = 0.0
        pz(j)    = 0.0
        fmass(j) = 0.0
  20  continue
      return
      end



      subroutine derivs(row)
      implicit none
      integer j, k, index, row
      real*8 spu, spo, outu, outo, tmp, a, b, dy, dp, drdp, dpj
      real*8 rxjku, ryjku, rzjku, rjku, pxjku, pyjku, pzjku, pjku
      real*8 dpdr
      real*8 rmom,delpx,delpy,delpz
      integer outdum
      logical iPau
      include 'coms.inc'
      include 'options.inc'
C aopx(j,?) = -dH/dx_j
C aorx(j,?) = dH/dpx_j


      do 10 j=1,nbar
        rww(j) = 0.0
   10 continue  
      if(logSky) then 
      do 20 j=1,nbar
        do 30 k=j+1,nbar
          rxjku = (airx(j)-airx(k))
          ryjku = (airy(j)-airy(k))
          rzjku = (airz(j)-airz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          if(rjku.lt.cutdww) then
            index = int(rjku/delr)+1
            a = dble(index) - rjku/delr
            b = 1.0 - a 
            tmp = a*spdwwy(index)+b*spdwwy(index+1)+
     +       ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
c            rww(j) = rww(j) + tmp
c            rww(k) = rww(k) + tmp
            rww(j) = rww(j) + ityp(k)/iabs(ityp(k))*tmp
            rww(k) = rww(k) + ityp(j)/iabs(ityp(j))*tmp
c            write(*,*)j,rww(j),ityp(j)
          end if
  30    continue
  20  continue
      end if

      do 40 j=1,nbar
        aopx(j,row) = 0.0
        aopy(j,row) = 0.0
        aopz(j,row) = 0.0
        dpj = 1.0/sqrt(aipx(j)*aipx(j)+aipy(j)*aipy(j)+aipz(j)*aipz(j)+
     +             fmass(j)*fmass(j))
        aorx(j,row) = aipx(j)*dpj
        aory(j,row) = aipy(j)*dpj
        aorz(j,row) = aipz(j)*dpj
  40  continue

      do 50 j=1,nbar
        do 60 k=j+1,nbar
          rxjku = (airx(j)-airx(k))
          ryjku = (airy(j)-airy(k))
          rzjku = (airz(j)-airz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          if (rjku.ge.1.0E-8) then
            rxjku = rxjku/rjku
            ryjku = ryjku/rjku
            rzjku = rzjku/rjku
          else
            rxjku = 0.0
            ryjku = 0.0
            rzjku = 0.0
          end if
          spu  = 0.0
          spo  = 0.0
          outu = 0.0
          outo = 0.0
          dy = 0.0
          index = int(rjku/delr)+1
          a = dble(index)-rjku/delr
          b = 1.0-a
          if(logYuk.and.rjku.lt.cutYuk) then
            spu  = spu  + spYuky(index)
            spo  = spo  + spYuky(index+1)
            outu = outu + outYuk(index)
            outo = outo + outYuk(index+1)
          end if
c            write(*,*)j,k,rjku 
          if(logSky.and.rjku.lt.cutdww) then
            if(ctoption(63).eq.1)then
            tmp=dpdr(rww(j),0.0d0,0)+dpdr(rww(k),0.0d0,0)
c            if(rww(j).lt.0.1.or.rww(k).lt.0.1)then
c            write(*,*)"cto1 ",j,k,tmp,
c     &       rww(j),dpdr(rww(j),0.0d0,0),rww(k),dpdr(rww(k),0.0d0,0)
c            endif
            else if(ctoption(63).ge.2)then
            delpx=aipx(j)-aipx(k)
            delpy=aipy(j)-aipy(k)
            delpz=aipz(j)-aipz(k)
c            outdum=0
c            if(ityp(j).le.0.or.ityp(k).le.0)then
c             outdum=1
c             write(*,*)"before: ",ityp(j),ityp(k),rww(j),rww(k)
c     &      ,dpdr(rww(j),rmom,ityp(j)),dpdr(rww(k),rmom,ityp(k))
c     &      ,spdwwy(index)
c            endif
            rmom=dsqrt(delpx**2+delpy**2+delpz**2)
            tmp=dpdr(rww(j),rmom,ityp(j))+dpdr(rww(k),rmom,ityp(k))

c            if(outdum.eq.1)then
c             write(*,*)"after: ",ityp(j),ityp(k),rww(j),rww(k)
c            endif

c             if(rww(j).lt.0.1.or.rww(k).lt.0.1)then
c            write(*,*)"cto2", j,k,tmp,
c     &       rww(j),dpdr(rww(j),rmom,ityp(j))
c     &      ,rww(k),dpdr(rww(k),rmom,ityp(k)),rmom
c              endif
            else
            tmp = Sky20 + Sky30*gamSky/(gamSky+1.0)*
     *            (rww(j)**(gamSky-1.0)+rww(k)**(gamSky-1.0))
c            write(*,*)"cto0 ",tmp
            endif
            spu  = spu  + spdwwy(index)*tmp
            spo  = spo  + spdwwy(index+1)*tmp
            outu = outu + outdww(index)*tmp
            outo = outo + outdww(index+1)*tmp
          end if    
          dy = da*spu+db*spo+
     +      ((3.0*a**2-1.0)*da*outu+(3.0*b**2-1.0)*db*outo)*fdel

          if(logCb) then
          if(rjku.lt.cutCb) then
            dy = dy + (da*spCby(index)+db*spCby(index+1)+
     +                ((3.0*a**2-1.0)*da*outCb(index)+
     +                (3.0*b**2-1.0)*db*outCb(index+1))*fdel)*
     *                dble(charge(j)*charge(k))
          else 
            dy = dy - Cb0/rjku/rjku*dble(charge(j)*charge(k))
          end if
          end if

          if(logPau.and.iPau(j,k)) then
            pxjku = (aipx(j)-aipx(k))
            pyjku = (aipy(j)-aipy(k))
            pzjku = (aipz(j)-aipz(k))
            pjku = sqrt(pxjku**2+pyjku**2+pzjku**2)
            if (pjku.ge.1.0E-8) then
              pxjku = pxjku/pjku
              pyjku = pyjku/pjku
              pzjku = pzjku/pjku
            else
              pxjku = 0.0
              pyjku = 0.0
              pzjku = 0.0
            end if
            drdp = 0.5*(pjku*pjku/dpPau+rjku*rjku/drPau)
            if(drdp.lt.cutPau) then
              index = int(drdp/delr)+1   
              a = dble(index)-drdp/delr
              b = 1.0-a
              tmp = da*spPauy(index)+db*spPauy(index+1)+
     +              ((3.0*a**2-1.0)*da*outPau(index)+
     +               (3.0*b**2-1.0)*db*outPau(index+1))*fdel
              dy = dy+tmp*rjku/drPau
              dp = tmp*pjku/dpPau*0.001
              aorx(j,row) = aorx(j,row)+dp*pxjku
              aory(j,row) = aory(j,row)+dp*pyjku
              aorz(j,row) = aorz(j,row)+dp*pzjku
              aorx(k,row) = aorx(k,row)-dp*pxjku
              aory(k,row) = aory(k,row)-dp*pyjku
              aorz(k,row) = aorz(k,row)-dp*pzjku
            end if
          end if 
          dy = -(0.001*dy)
          aopx(j,row) = aopx(j,row)+dy*rxjku
          aopy(j,row) = aopy(j,row)+dy*ryjku
          aopz(j,row) = aopz(j,row)+dy*rzjku
          aopx(k,row) = aopx(k,row)-dy*rxjku
          aopy(k,row) = aopy(k,row)-dy*ryjku
          aopz(k,row) = aopz(k,row)-dy*rzjku
  60    continue
  50  continue

c        if(logCb.and.CTOption(64).eq.1) then
c        do j = 1,npart
c          do k=j+1,npart
c          if(logCb.and.CTOption(64).eq.0) then
c          if(rjku.lt.cutCb) then
c            dy = dy + (da*spCby(index)+db*spCby(index+1)+
c     +                ((3.0*a**2-1.0)*da*outCb(index)+
c     +                (3.0*b**2-1.0)*db*outCb(index+1))*fdel)*
c     *                dble(charge(j)*charge(k))
c          else
c            dy = dy - Cb0/rjku/rjku*dble(charge(j)*charge(k))
c          end if
c          end if
c        enddo
c        enddo
c        end if


      return
      end

      real*8 function Ekintot() 
      implicit none
      integer j
      real*8 Ekin
      include 'coms.inc'

      Ekintot = 0.0
      do 3 j=1,npart
         Ekintot= Ekintot+Ekin(j)
 3    continue
      return
      end

      real*8 function EtotJK() 
      implicit none
      real*8 Etot
      integer j
      include 'coms.inc'

      EtotJK = Etot()
      do 3 j=1,npart
         EtotJK= EtotJK-fmass(j)
 3    continue
      EtotJK = EtotJK/npart
      return
      end


      real*8 function Etot() 
      implicit none
      integer j, k, index
      real*8 a, b, y, drdp, tp, tr, tmp, Ekintot
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau, Ekin
      real*8 rxjku, ryjku, rzjku, rjku, pxjku, pyjku, pzjku, pjku
      logical iPau
      real*8 ecmf
      include 'coms.inc'
      include 'options.inc'
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau

      Etot = 0.0
      Ekinbar = 0.0
      Ekinmes = 0.0
      ESky2 = 0.0
      ESky3 = 0.0
      EYuk = 0.0
      ECb = 0.0
      EPau = 0.0

      if(EoS.eq.0) then

c CASCADE mode
         Etot=Ekintot()
         return
      else
c with potentials
c kinetic energies of mesons first
         do 4 j=nbar+1,npart
            Etot= Etot+Ekin(j)
            Ekinmes = Ekinmes+Ekin(j)
 4       continue

      do 10 j=1,nbar
        rww(j) = 0.0
   10 continue  
      if(logSky) then 
      do 20 j=1,nbar
        do 30 k=j+1,nbar
          rxjku = (rx(j)-rx(k))
          ryjku = (ry(j)-ry(k))
          rzjku = (rz(j)-rz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          if(rjku.lt.cutdww) then
            index = int(rjku/delr)+1
            a = dble(index) - rjku/delr
            b = 1.0 - a 
            tmp = a*spdwwy(index)+b*spdwwy(index+1)+
     +       ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
            rww(j) = rww(j) + tmp
            rww(k) = rww(k) + tmp
          end if
  30    continue
  20  continue
      end if

      do 40 j=1,nbar
       if(ctoption(63).ge.1)then
        Etot = Etot + Ekin(j)+ecmf(rww(j),0.0d0,ityp(j))/1000.0d0
       else
 
       Etot = Etot + Ekin(j) + 0.0005*Sky20*rww(j) +
     +         0.001*Sky30/(gamSky+1.0)*rww(j)**gamSky
       endif
        Ekinbar = Ekinbar + Ekin(j)
        ESky2 = ESky2 + 0.0005*Sky20*rww(j)
        ESky3 = ESky3 + 0.001*Sky30/(gamSky+1.0)*rww(j)**gamSky
        do 50 k=j+1,nbar
          rxjku = (rx(j)-rx(k))
          ryjku = (ry(j)-ry(k))
          rzjku = (rz(j)-rz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          index = int(rjku/delr)+1
          a = dble(index)-rjku/delr
          b = 1.0-a
          if(logYuk.and.rjku.lt.cutYuk) then
            y = a*spYuky(index)+b*spYuky(index+1)+
     +          ((a**3-a)*outYuk(index)+(b**3-b)*outYuk(index+1))*fdel
            Etot = Etot + 0.001*y
            EYuk = EYuk + 0.001*y
          end if
          if(logCb) then
          if(rjku.lt.cutCb) then
            y = (a*spCby(index)+b*spCby(index+1)+
     +          ((a**3-a)*outCb(index)+(b**3-b)*outCb(index+1))*fdel)*
     *          dble(charge(j)*charge(k))
          else
            y = Cb0/rjku*dble(charge(j)*charge(k))
          end if
          Etot = Etot + 0.001*y
          ECb = ECb + 0.001*y
          end if
          if(logPau.and.iPau(j,k)) then
            pxjku = (px(j)-px(k))
            pyjku = (py(j)-py(k))
            pzjku = (pz(j)-pz(k))
            pjku = sqrt(pxjku**2+pyjku**2+pzjku**2)
            tp = pjku
            tr = rjku
            drdp = 0.5*(pjku*pjku/dpPau+rjku*rjku/drPau)
            if(drdp.lt.cutPau) then
              index = int(drdp/delr)+1
              a = dble(index)-drdp/delr
              b = 1.0-a
              y = a*spPauy(index)+b*spPauy(index+1)+
     +            ((a**3-a)*outPau(index)+(b**3-b)*outPau(index+1))*fdel
              Etot = Etot + 0.001*y
              EPau = EPau + 0.001*y
            end if
          end if
  50    continue
  40  continue
      Ekinbar = Ekinbar/dble(nbar)
      Ekinmes = Ekinmes/dble(max(1,npart-nbar))
      ESky2 = ESky2/dble(nbar)
      ESky3 = ESky3/dble(nbar)
      EYuk = EYuk/dble(nbar)
      ECb = ECb/dble(nbar)
      EPau = EPau/dble(nbar)
      end if
      return
      end
          

      subroutine cascstep(tim,dtime)
      implicit none
      real*8 tim,dtime,energ
      integer j
      include 'coms.inc'
      include 'boxinc.inc'
      include 'options.inc'

      do 1 j=1,npart
         energ = sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
         r0(j) = r0(j) + dtime
         rx(j) = rx(j) + px(j)/energ*dtime
         rz(j) = rz(j) + pz(j)/energ*dtime
         ry(j) = ry(j) + py(j)/energ*dtime
1     continue          
      return
      end
 
      subroutine proprk(tim,dtime)

      implicit none
      real*8 tim,dtime,energ,dt2, dt6
      integer j

      include 'coms.inc'
      include 'boxinc.inc'
      include 'options.inc'

      if (EoS.eq.0) then
c  cascade mode
         do 1 j=1,npart
               energ = p0(j)    ! sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
               r0(j) = r0(j) + dtime
               rx(j) = rx(j) + px(j)/energ*dtime
               ry(j) = ry(j) + py(j)/energ*dtime
               rz(j) = rz(j) + pz(j)/energ*dtime
1        continue          
         return
      else
c propagation with potentials
c propagate mesons on straight lines
         do 2 j=nbar+1,npart
            energ = p0(j) ! sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
            r0(j) = r0(j) + dtime
            rx(j) = rx(j) + px(j)/energ*dtime
            ry(j) = ry(j) + py(j)/energ*dtime
            rz(j) = rz(j) + pz(j)/energ*dtime
2        continue          

c propagate baryons
c adjust time-step parameters
      dt = dtime
      dt2 = dtime/2.0d0
      dt6 = dtime/6.0d0

      do 10 j=1,nbar
        airx(j) = rx(j) 
        airy(j) = ry(j) 
        airz(j) = rz(j) 
        aipx(j) = px(j)
        aipy(j) = py(j)
        aipz(j) = pz(j)
  10  continue
        
      call derivs(1)
      do 20 j=1,nbar
        airx(j) = rx(j) + dt2*aorx(j,1) 
        airy(j) = ry(j) + dt2*aory(j,1) 
        airz(j) = rz(j) + dt2*aorz(j,1) 
        aipx(j) = px(j) + dt2*aopx(j,1)
        aipy(j) = py(j) + dt2*aopy(j,1)
        aipz(j) = pz(j) + dt2*aopz(j,1)
  20  continue

      call derivs(2)
        
      do 30 j=1,nbar
        airx(j) = rx(j) + dt2*aorx(j,2) 
        airy(j) = ry(j) + dt2*aory(j,2) 
        airz(j) = rz(j) + dt2*aorz(j,2) 
        aipx(j) = px(j) + dt2*aopx(j,2)
        aipy(j) = py(j) + dt2*aopy(j,2)
        aipz(j) = pz(j) + dt2*aopz(j,2)
  30  continue
      
      call derivs(3)

      do 40 j=1,nbar
        airx(j) = rx(j) + dt*aorx(j,3)
        airy(j) = ry(j) + dt*aory(j,3)
        airz(j) = rz(j) + dt*aorz(j,3)
        aipx(j) = px(j) + dt*aopx(j,3)
        aipy(j) = py(j) + dt*aopy(j,3)
        aipz(j) = pz(j) + dt*aopz(j,3)
  40  continue
   
      call derivs(4)

      do 50 j=1,nbar
         r0(j) = r0(j) + dtime
        rx(j)=rx(j)+dt6*(aorx(j,1)+2.0*(aorx(j,2)+aorx(j,3))+aorx(j,4))
        ry(j)=ry(j)+dt6*(aory(j,1)+2.0*(aory(j,2)+aory(j,3))+aory(j,4))
        rz(j)=rz(j)+dt6*(aorz(j,1)+2.0*(aorz(j,2)+aorz(j,3))+aorz(j,4))
        px(j)=px(j)+dt6*(aopx(j,1)+2.0*(aopx(j,2)+aopx(j,3))+aopx(j,4))
        py(j)=py(j)+dt6*(aopy(j,1)+2.0*(aopy(j,2)+aopy(j,3))+aopy(j,4))
        pz(j)=pz(j)+dt6*(aopz(j,1)+2.0*(aopz(j,2)+aopz(j,3))+aopz(j,4))
        p0(j)=sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
 50   continue
      end if

      return
      end


      subroutine potPau
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy, Pau
      include 'coms.inc'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      px(1) = 0.0d0
      py(1) = 0.0d0
      pz(1) = 0.0d0
      px(1) = 0.0d0
      py(2) = 0.0d0
      pz(2) = 0.0d0
      Ecut = 1.0E-5
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = sqrt(2.0*dr*drPau)
      spx(i) = dr
      spPauy(i) = Pau(1,2)
      if(spPauy(i).lt.Ecut) then
        spPauy(i) = 0.0
        cutPau = dr
        abl0 = -Pau0 
        abln = 0.0
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spPauy,ncut,abl0,abln,outPau)

      write(6,'(''Pauli-Potential    '',e10.3,i5,f7.1)') 
     +      Ecut, ncut, cutPau

      do 10 i=0,20
        dr = 0.323*dble(i)
        if(dr.gt.cutPau) then
          y = 0.0
          dy = 0.0
        else
        rx(2) = dr
        dr = 0.5*dr*dr/drPau
        index = int(dr/delr)+1
        a = dble(index) - dr/delr
        b = 1.0 - a
        y = a*spPauy(index)+b*spPauy(index+1)+
     +     ((a**3-a)*outPau(index)+
     +      (b**3-b)*outPau(index+1))*fdel
        dy = da*spPauy(index)+db*spPauy(index+1)+
     +        ((3.0*a**2-1.0)*da*outPau(index)+
     +         (3.0*b**2-1.0)*db*outPau(index+1))*fdel
        dy = dy*sqrt(2.0*dr*drPau)/drPau
        end if
  10  continue
      return
      end

      subroutine potCb
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy, dCb, Cb
      include 'coms.inc'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      Ecut = 1.0E-5
      iso3(1) = 1
      iso3(2) = 1
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = dr
      spx(i) = dr
      spCby(i) = Cb(1,2)
      if(abs(spCby(i)*dr-Cb0)/max(dr,1.0d-5).lt.Ecut) then
        spCby(i) = Cb0/dr
        cutCb  = dr
        abln = dCb(1,2)
        abl0 = 0.0 
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spCby,ncut,abl0,abln,outCb)

      write(6,'(''Coulomb-Potential  '',e10.3,i5,f7.1)')
     +      Ecut, ncut, cutCb

      do 10 i=0,20
        dr = 0.2*dble(i)+0.01212
        rx(2) = dr
        if(dr.ge.cutCb) then
          y = Cb0/dr
          dy = -(Cb0/dr/dr)
        else
          index = int(dr/delr)+1
          a = dble(index) - dr/delr
          b = 1.0 - a
          y = a*spCby(index)+b*spCby(index+1)+
     +       ((a**3-a)*outCb(index)+(b**3-b)*outCb(index+1))*fdel
          dy = da*spCby(index)+db*spCby(index+1)+
     +         ((3.0*a**2-1.0)*da*outCb(index)+
     +           (3.0*b**2-1.0)*db*outCb(index+1))*fdel
        end if
  10  continue
      return
      end

      subroutine potYuk
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy
      real*8 Yuk
      include 'coms.inc'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      Ecut = 1.0E-5
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = dr 
      spx(i) = dr
      spYuky(i) = Yuk(1,2)
      if(abs(spYuky(i)).lt.Ecut) then
        spYuky(i) = 0.0
        cutYuk = dr
        abl0 = 0.0 
        abln = 0.0
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spYuky,ncut,abl0,abln,outYuk)

      write(6,'(''Yukawa-Potential   '',e10.3,i5,f7.1)') 
     +      Ecut, ncut, cutYuk

      do 10 i=0,40
        dr = 0.2*dble(i)
        rx(2) = dr
        if(dr.gt.cutYuk) then
          y = 0.0
          dy = 0.0
        else
        index = int(dr/delr)+1
        a = dble(index) - dr/delr
        b = 1.0 - a
        y = a*spYuky(index)+b*spYuky(index+1)+
     +     ((a**3-a)*outYuk(index)+(b**3-b)*outYuk(index+1))*fdel
        dy = da*spYuky(index)+db*spYuky(index+1)+
     +        ((3.0*a**2-1.0)*da*outYuk(index)+
     +         (3.0*b**2-1.0)*db*outYuk(index+1))*fdel
        end if
  10  continue
      return
      end


      subroutine potdww
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy
      real*8 dww
      include 'coms.inc'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      Ecut = 1.0E-8
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = dr 
      spx(i) = dr
      spdwwy(i) = dww(1,2)
      if(abs(spdwwy(i)).lt.Ecut) then
        spdwwy(i) = 0.0
        cutdww = dr
        abl0 = 0.0 
        abln = 0.0
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spdwwy,ncut,abl0,abln,outdww)

      write(6,'(''Interaction-Density'',e10.3,i5,f7.1)') 
     +      Ecut, ncut, cutdww

      do 10 i=0,20
        dr = 0.295*dble(i)
        rx(2) = dr
        if(dr.gt.cutdww) then
          y = 0.0
          dy = 0.0
        else
        index = int(dr/delr)+1
        a = dble(index) - dr/delr
        b = 1.0 - a
        y = a*spdwwy(index)+b*spdwwy(index+1)+
     +     ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
        dy = da*spdwwy(index)+db*spdwwy(index+1)+
     +        ((3.0*a**2-1.0)*da*outdww(index)+
     +         (3.0*b**2-1.0)*db*outdww(index+1))*fdel
        end if
  10  continue
      return
      end

c Kinetic Energy
c
      function Ekin(j)
      implicit none
      integer j
      real*8 Ekin
      include 'coms.inc'
 
      Ekin = sqrt((px(j)+ffermpx(j))*(px(j)+ffermpx(j))+
     +            (py(j)+ffermpy(j))*(py(j)+ffermpy(j))+
     +            (pz(j)+ffermpz(j))*(pz(j)+ffermpz(j))+
     +            fmass(j)*fmass(j)) 

      return
      end

c Derivative for Kinetic Energy
c
      function dEkin(j)
      implicit none
      integer j
      real*8 dEkin
      include 'coms.inc'
      
      dEkin = 1.0/sqrt(px(j)*px(j)+py(j)*py(j)+pz(j)*pz(j)+
     +                 fmass(j)*fmass(j)) 
      return
      end

c Skyrme Potential (3-body) rwwterm
c
      function dww(j,k)
      implicit none
      integer j, k
      real*8 dww, rjk
      include 'coms.inc'
      
      dww = gw/pi*sqrt(gw/pi)*exp(-(gw*rjk(j,k)*rjk(j,k)))/
     /      rho0
      return
      end

cc Skyrme Potential
cc
      function Sky(j,k)
      implicit none
      integer j, k
      real*8 Sky, rjk
      include 'coms.inc'
      
      Sky = Sky20*gw/pi*sqrt(gw/pi)*exp(-(gw*rjk(j,k)*rjk(j,k))) 
      return
      end

c Coulomb Potential
c
      function Cb(j,k)
      implicit none
      integer j, k
      real*8 Cb, rjk
      real*8 erf
      include 'coms.inc'
      
      if (iso3(j).eq.1.and.iso3(k).eq.1) then
        if (rjk(j,k).lt.eps) then
          Cb = Cb0*er0*sgw
        else
          Cb = Cb0/rjk(j,k)*erf(sgw*rjk(j,k))
        end if
      else 
        Cb = 0.0
      end if 
      return
      end

c Derivative for Coulomb Potential
c
      function dCb(j,k)
      implicit none
      integer j, k
      real*8 dCb, rjk
      real*8 erf
      include 'coms.inc'
      
      if (iso3(j).eq.1.and.iso3(k).eq.1) then
        if (rjk(j,k).lt.eps) then
          dCb = 0.0
        else
          dCb = Cb0*(er0*exp(-(gw*rjk(j,k)*rjk(j,k)))*sgw*rjk(j,k)-
     +               erf(sgw*rjk(j,k)))/rjk(j,k)/rjk(j,k)
        end if
      else 
        dCb = 0.0
      end if 
      return
      end

c Yukawa Potential
c
      function Yuk(j,k)
      implicit none
      integer j, k
      real*8 Yuk, rjk
      real*8 erf
      include 'coms.inc'
      
      if(rjk(j,k).lt.eps) then
        Yuk = Yuk0*(er0*sgw-exp(0.25/gamYuk/gamYuk/gw)/gamYuk*
     *              (1.0-erf(0.5/gamYuk/sgw)))
      else
        Yuk = Yuk0*0.5/rjk(j,k)*exp(0.25/gamYuk/gamYuk/gw)*
     *           (exp(-(rjk(j,k)/gamYuk))*
     +            (1.0-erf(0.5/gamYuk/sgw-sgw*rjk(j,k)))-
     -            exp(rjk(j,k)/gamYuk)*
     +            (1.0-erf(0.5/gamYuk/sgw+sgw*rjk(j,k))))
      end if
      return
      end

c Derivative for Yukawa Potential
c
      function dYuk(j,k)
      implicit none
      integer j, k
      real*8 dYuk, rjk
      real*8 erf
      include 'coms.inc'
      
      if(rjk(j,k).lt.eps) then
        dYuk = 0.0
      else
        dYuk = 0.5*Yuk0/rjk(j,k)*( exp(0.25/gamYuk/gamYuk/gw)*(
     *          (-(1.0/rjk(j,k))-1.0/gamYuk)*exp(-(rjk(j,k)/gamYuk))*
     *             (1.0-erf(0.5/gamYuk/sgw-sgw*rjk(j,k))) +
     *          (1.0/rjk(j,k)-1.0/gamYuk)*exp(rjk(j,k)/gamYuk)*
     *             (1.0-erf(0.5/gamYuk/sgw+sgw*rjk(j,k))) ) +
     +          sgw*er0*2.0*exp(-(gw*rjk(j,k)*rjk(j,k))) )
      end if
      return
      end


c Pauli Potential
c
      function Pau(j,k)
      implicit none
      integer j, k
      real*8 Pau, pjk, rjk
      include 'coms.inc'
      
      Pau = Pau0*exp(-(0.5*rjk(j,k)*rjk(j,k)/drPau))*
     *           exp(-(0.5*pjk(j,k)*pjk(j,k)/dpPau))
      return
      end

c Derivative (p) for Pauli Potential
c
      function dPaup(j,k)
      implicit none
      integer j, k
      real*8 dPaup, pjk, rjk
      include 'coms.inc'
      
      dPaup = -(Pau0/dpPau*pjk(j,k)*
     *                    exp(-(0.5*rjk(j,k)*rjk(j,k)/drPau))*
     *                    exp(-(0.5*pjk(j,k)*pjk(j,k)/dpPau)))
      return
      end
  
c Derivative (r) for Pauli Potential
c
      function dPaur(j,k)
      implicit none
      integer j, k
      real*8 dPaur, pjk, rjk
      include 'coms.inc'
      
      dPaur = -(Pau0/drPau*rjk(j,k)*
     *                   exp(-(0.5*rjk(j,k)*rjk(j,k)/drPau))*
     *                   exp(-(0.5*pjk(j,k)*pjk(j,k)/dpPau)))
      return
      end

      function rjk(j,k)
      implicit none 
      integer j, k
      real*8 rjk
      include 'coms.inc'

      rjk = sqrt((rx(j)-rx(k))**2+(ry(j)-ry(k))**2+(rz(j)-rz(k))**2)           
      return
      end

      function pjk(j,k)
      implicit none 
      integer j, k
      real*8 pjk
      include 'coms.inc'

      pjk = sqrt((px(j)-px(k))**2+(py(j)-py(k))**2+(pz(j)-pz(k))**2)           
      return
      end

      function iPau(j,k)
      implicit none 
      integer j, k
      logical iPau
      include 'coms.inc'

      iPau = .false.
      if (iso3(j).eq.iso3(k).and.ityp(j).eq.ityp(k)) iPau = .true.
      return
      end


c---------------------------
      subroutine readCMF(nmom)
      implicit none
      
      real*8 rnb,rps,rder,udum
      real*8 mom,vnst,dvnst,vdel,dvdel,vdelst,dvdelst
      real*8 vlam,dvlam,vlamst,dvlamst,vsig,dvsig
      real*8 vsigst,dvsigst,vxi,dvxi,vxist,dvxist
      real*8 CMFskyrmetab(0:2,0:1000)
      real*8 CMFaskyrmetab(0:2,0:1000)
      real*8 CMFmomSkyrme(0:601,0:300,0:19)
      real*8 CMFmomaSkyrme(0:601,0:300,0:19)

      integer ir,im,nmom,ire
      common/CMF/ CMFskyrmetab,CMFmomSkyrme,
     & CMFaskyrmetab,CMFmomaSkyrme

         
        if(nmom.eq.1)then
c       read non-mom dep. table
            open(unit=76,
     $     file='eosfiles/CMF_skyrme.dat')


         do ir = 0,600
            read(76,*) rnb,rps,rder
            CMFskyrmetab(0,ir)=rnb
            CMFskyrmetab(1,ir)=rps
            CMFskyrmetab(2,ir)=rder
        enddo
        close(76)
         
        else
c       read mom dep. table
        open(unit=76,
     $     file='eosfiles/CMF_mom_dep.dat')
        open(unit=77,
     $     file='eosfiles/CMF_mom_dep_a.dat')

        do ire=0,180900
          read(76,*) rnb,mom,Udum,rps,rder,Udum,vnst,dvnst,
     &      Udum,vdel,dvdel,Udum,vdelst,dvdelst,Udum,vlam,dvlam,
     &      Udum,vlamst,dvlamst,Udum,vsig,dvsig,
     &      Udum,vsigst,dvsigst,udum,vxi,dvxi,
     &      Udum,vxist,dvxist

           ir=(rnb+0.001)/0.05
           im=(mom+0.001)/0.05

           CMFmomSkyrme(ir,im,0)=rps
           CMFmomSkyrme(ir,im,1)=rder
           CMFmomSkyrme(ir,im,2)=vnst
           CMFmomSkyrme(ir,im,3)=dvnst
           CMFmomSkyrme(ir,im,4)=vdel
           CMFmomSkyrme(ir,im,5)=dvdel
           CMFmomSkyrme(ir,im,6)=vdelst
           CMFmomSkyrme(ir,im,7)=dvdelst
           CMFmomSkyrme(ir,im,8)=vlam
           CMFmomSkyrme(ir,im,9)=dvlam
           CMFmomSkyrme(ir,im,10)=vlamst
           CMFmomSkyrme(ir,im,11)=dvlamst
           CMFmomSkyrme(ir,im,12)=vsig
           CMFmomSkyrme(ir,im,13)=dvsig
           CMFmomSkyrme(ir,im,14)=vsigst
           CMFmomSkyrme(ir,im,15)=dvsigst
           CMFmomSkyrme(ir,im,16)=vxi
           CMFmomSkyrme(ir,im,17)=dvxi
           CMFmomSkyrme(ir,im,18)=vxist
           CMFmomSkyrme(ir,im,19)=dvxist

          read(77,*) rnb,mom,Udum,rps,rder,Udum,vnst,dvnst,
     &      Udum,vdel,dvdel,Udum,vdelst,dvdelst,Udum,vlam,dvlam,
     &      Udum,vlamst,dvlamst,Udum,vsig,dvsig,
     &      Udum,vsigst,dvsigst,udum,vxi,dvxi,
     &      Udum,vxist,dvxist

           ir=(rnb+0.001)/0.05
           im=(mom+0.001)/0.05

           CMFmomaSkyrme(ir,im,0)=rps
           CMFmomaSkyrme(ir,im,1)=rder
           CMFmomaSkyrme(ir,im,2)=vnst
           CMFmomaSkyrme(ir,im,3)=dvnst
           CMFmomaSkyrme(ir,im,4)=vdel
           CMFmomaSkyrme(ir,im,5)=dvdel
           CMFmomaSkyrme(ir,im,6)=vdelst
           CMFmomaSkyrme(ir,im,7)=dvdelst
           CMFmomaSkyrme(ir,im,8)=vlam
           CMFmomaSkyrme(ir,im,9)=dvlam
           CMFmomaSkyrme(ir,im,10)=vlamst
           CMFmomaSkyrme(ir,im,11)=dvlamst
           CMFmomaSkyrme(ir,im,12)=vsig
           CMFmomaSkyrme(ir,im,13)=dvsig
           CMFmomaSkyrme(ir,im,14)=vsigst
           CMFmomaSkyrme(ir,im,15)=dvsigst
           CMFmomaSkyrme(ir,im,16)=vxi
           CMFmomaSkyrme(ir,im,17)=dvxi
           CMFmomaSkyrme(ir,im,18)=vxist
           CMFmomaSkyrme(ir,im,19)=dvxist


        enddo
        endif


      return
      end


c--------------------------
      function ecmf(densb,mom,id)
      implicit none


      include 'options.inc'       


      integer ir,im,id,line
      real*8 ecmf,densb,drho,step,decmf
      real*8 pa1,pa2,pa3,pa4,pa5,pa6,pa7
      real*8 CMFskyrmetab(0:2,0:1000)
      real*8 CMFaskyrmetab(0:2,0:1000)
      real*8 mom
      real*8 CMFmomSkyrme(0:601,0:300,0:19)
      real*8 CMFmomaSkyrme(0:601,0:300,0:19)
      real*8 p1,p2,p3,p4,p12,p13,p24,p34,dn,de
      common/CMF/ CMFskyrmetab,CMFmomSkyrme,
     & CMFaskyrmetab,CMFmomaSkyrme

      if (densb.gt.30.00d0)densb=30.00d0
      if (mom.gt.14.999d0)mom=14.999d0
      if (densb.lt.0.0d0)densb=0.0d0

      if(ctoption(65).eq.1)then
         pa1=ctparam(76)
         pa2=ctparam(77)
         pa3=ctparam(78)
         pa4=ctparam(79)
         pa5=ctparam(80)
         pa6=ctparam(81)
         pa7=ctparam(82)
         
         if(densb.gt.2.d0)then
            ecmf=-22.07
     $           +pa1*(densb-2.0d0)
     $           +pa2*(densb-2.0d0)**2
     $           +pa3*(densb-2.0d0)**3
     $           +pa4*(densb-2.0d0)**4
     $           +pa5*(densb-2.0d0)**5
     $           +pa6*(densb-2.0d0)**6
     $           +pa7*(densb-2.0d0)**7
            return
         endif     
      endif    

      if(ctoption(63).eq.1)then       
      do ir = 0,599
         if(CMFskyrmetab(0,ir).lt.densb
     $   .and.CMFskyrmetab(0,ir+1).gt.densb)then
         drho=CMFskyrmetab(0,ir+1)-CMFskyrmetab(0,ir)
         step=densb-CMFskyrmetab(0,ir)
         decmf=CMFskyrmetab(1,ir+1)-CMFskyrmetab(1,ir)
         ecmf=CMFskyrmetab(1,ir)+step/drho*decmf
         return
         endif
      end do
      endif  

      if(ctoption(63).eq.2)then

cJS     First get correct potential from itype
c N
      line=0
c N*
      if(id.ge.4.and.id.le.23)line=2
c Delta
      if(id.ge.24.and.id.le.26)line=4
c Delta*
      if(id.ge.27.and.id.le.40)line=6
c Lambda
      if(id.ge.41.and.id.le.45)line=8
c Lambda*
      if(id.ge.46.and.id.le.53)line=10
c Sigma
      if(id.ge.54.and.id.le.57)line=12
c Sigma*
      if(id.ge.58.and.id.le.62)line=14
c Xi
      if(id.ge.63.and.id.le.65)line=16
c Xi*
      if(id.ge.66.and.id.le.68)line=18
c Omega (no Omega yet)
      if(id.ge.69)then
      ecmf=0.0d0   
      return
      endif


      if(ctoption(63).eq.2)line=0


      dn=0.05d0
      de=0.05d0
      ir=densb/dn
      im=mom/dn
      p1 = CMFmomSkyrme(ir,im,line)
      p2 = CMFmomSkyrme(ir+1,im,line)
      p3 = CMFmomSkyrme(ir,im+1,line)
      p4 = CMFmomSkyrme(ir+1,im+1,line)
      p13 = (p3-p1)/dn*(densb-dfloat(idint(densb/dn))*dn)+p1
      p24 = (p4-p2)/dn*(densb-dfloat(idint(densb/dn))*dn)+p2
      p12 = (p2-p1)/de*(mom-dfloat(idint(mom/de))*de) + p1
      p34 = (p4-p3)/de*(mom-dfloat(idint(mom/de))*de) + p3

      ecmf = (p24-p13)/de*(mom-dfloat(idint(mom/de))*de) + p13
           
      return

      endif


      return
      end

c--------------------------
      function dpdr(densb,mom,id)
      implicit none

      include 'options.inc'

      integer ir,im,id,line,aid,typ
      real*8 dpdr,densb,drho,step,decmf
      real*8 rhob
      real*8 CMFskyrmetab(0:2,0:1000)
      real*8 CMFaskyrmetab(0:2,0:1000)
      real*8 pa1,pa2,pa3,pa4,pa5,pa6,pa7
      real*8 mom,apf
      real*8 CMFmomSkyrme(0:601,0:300,0:19)
      real*8 CMFmomaSkyrme(0:601,0:300,0:19)
      real*8 p1,p2,p3,p4,p12,p13,p24,p34,dn,de
      common/CMF/ CMFskyrmetab,CMFmomSkyrme,
     & CMFaskyrmetab,CMFmomaSkyrme

      apf=1.0d0
      rhob=densb
      typ=id

      if (rhob.lt.0.0d0.and.ctoption(63).lt.4)rhob=0.0d0
      if (rhob.lt.0.0d0.and.ctoption(63).eq.4) then
        apf=-1.0d0
        rhob=-rhob
        typ=-typ
      endif     
      if (rhob.gt.30.00d0)rhob=30.00d0
      if (mom.gt.14.999d0)mom=14.999d0

      if(ctoption(65).eq.1)then
         pa1=ctparam(76)
         pa2=ctparam(77)
         pa3=ctparam(78)
         pa4=ctparam(79)
         pa5=ctparam(80)
         pa6=ctparam(81)
         pa7=ctparam(82)
         
         if(densb.gt.2.d0)then
            dpdr=pa1
     $           +pa2*2.0d0*(rhob-2.0d0)
     $           +pa3*3.0d0*(rhob-2.0d0)**2
     $           +pa4*4.0d0*(rhob-2.0d0)**3
     $           +pa5*5.0d0*(rhob-2.0d0)**4
     $           +pa6*6.0d0*(rhob-2.0d0)**5
     $           +pa7*7.0d0*(rhob-2.0d0)**6
            return
         endif         
      endif

      if(ctoption(63).eq.1)then

c        if(densb.gt.1.60d0)then
c        dpdr=0.0d0
c        return
c       endif      
        do ir = 0,599
         if(CMFskyrmetab(0,ir).lt.rhob
     $   .and.CMFskyrmetab(0,ir+1).gt.rhob)then
         drho=CMFskyrmetab(0,ir+1)-CMFskyrmetab(0,ir)
         step=rhob-CMFskyrmetab(0,ir)
         decmf=CMFskyrmetab(2,ir+1)-CMFskyrmetab(2,ir)
         dpdr=CMFskyrmetab(2,ir)+step/drho*decmf
c         write(*,*)rhob,dpdr
         return
         endif
      end do
      endif

      if(ctoption(63).ge.2)then

      if(ctoption(63).le.3.and.id.lt.0)then
c     No anti-baryon potentials
        dpdr=0.0
        return
      endif

cJS   First get correct potential from itype

      aid=iabs(typ)
c N
      if(aid.ge.1)line=1
c N*
      if(aid.ge.4.and.aid.le.23)line=3
c Delta
      if(aid.ge.24.and.aid.le.26)line=5
c Delta*
      if(aid.ge.27.and.aid.le.40)line=7
c Lambda
      if(aid.ge.41.and.aid.le.45)line=9
c Lambda*
      if(aid.ge.46.and.aid.le.53)line=11
c Sigma
      if(aid.ge.54.and.aid.le.57)line=13
c Sigma*
      if(aid.ge.58.and.aid.le.62)line=15
c Xi
      if(aid.ge.63.and.aid.le.65)line=17
c Xi*
      if(aid.ge.66.and.aid.le.68)line=19
c Omega (no Omega yet)
      if(aid.ge.69)then
      dpdr=0.0d0
      return      
      endif


      if(ctoption(63).eq.2)line=1
  
       
      if(typ.gt.0)then
      dn=0.05d0
      de=0.05d0
      ir=rhob/dn
      im=mom/de
      p1 = CMFmomSkyrme(ir,im,line)
      p2 = CMFmomSkyrme(ir+1,im,line)
      p3 = CMFmomSkyrme(ir,im+1,line)
      p4 = CMFmomSkyrme(ir+1,im+1,line)
      p13 = (p3-p1)/dn*(rhob-dfloat(idint(rhob/dn))*dn)+p1
      p24 = (p4-p2)/dn*(rhob-dfloat(idint(rhob/dn))*dn)+p2
      p12 = (p2-p1)/de*(mom-dfloat(idint(mom/de))*de) + p1
      p34 = (p4-p3)/de*(mom-dfloat(idint(mom/de))*de) + p3

      dpdr = (p24-p13)/de*(mom-dfloat(idint(mom/de))*de) + p13
      dpdr = dpdr
      return
      endif

      if(typ.lt.0.and.ctoption(63).ge.4)then
      dn=0.05d0
      de=0.05d0
      ir=rhob/dn
      im=mom/de
      p1 = CMFmomaSkyrme(ir,im,line)
      p2 = CMFmomaSkyrme(ir+1,im,line)
      p3 = CMFmomaSkyrme(ir,im+1,line)
      p4 = CMFmomaSkyrme(ir+1,im+1,line)
      p13 = (p3-p1)/dn*(rhob-dfloat(idint(rhob/dn))*dn)+p1
      p24 = (p4-p2)/dn*(rhob-dfloat(idint(rhob/dn))*dn)+p2
      p12 = (p2-p1)/de*(mom-dfloat(idint(mom/de))*de) + p1
      p34 = (p4-p3)/de*(mom-dfloat(idint(mom/de))*de) + p3

      dpdr = (p24-p13)/de*(mom-dfloat(idint(mom/de))*de) + p13
      dpdr = dpdr
      return
      endif

      
      dpdr=0.0d0

      endif



      return
      end

