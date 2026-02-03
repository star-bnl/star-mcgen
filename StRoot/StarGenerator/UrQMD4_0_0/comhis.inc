c************************************************************8
c.. Author: Marcus Bleicher, LBNL, Dec. 5. 1999
c.. file 15 analysis...
c.. Author: Tom Reichert, GU, May 13. 2020
c.. kinetic freeze-out analysis
c.. Author: Tom Reichert, FIAS, Dec 29. 2023
c.. integrate analysis for all reconstructable resonances 
c.. into urqmd natively, write resonances to f13 output      
c.. 
c.. This file contains array declarations for the resonance reconstruction
c.. Therefore a single event has to be remembered until kinetic fo      
c************************************************************8

c. array sizes (original: 20000, 10, 50)
      integer wwmax,inmax,outmax,resomax
      parameter (wwmax=20000)
      parameter (inmax=10)
      parameter (outmax=50)
      parameter (resomax=40000)

c. interaction header, one for every interaction ww
      integer hisnin(wwmax),hisnexit(wwmax),hisiline(wwmax)
      integer hisctag(wwmax)
      real*8 hisacttime(wwmax),hissqrts(wwmax),hisstot(wwmax)
      real*8 hissigpart(wwmax),hiscdens(wwmax)

c.. INgoing particles
      integer INind(wwmax,inmax)
      real*8 INr0(wwmax,inmax),INrx(wwmax,inmax)
      real*8 INry(wwmax,inmax),INrz(wwmax,inmax)
      real*8 INp0(wwmax,inmax),INpx(wwmax,inmax)
      real*8 INpy(wwmax,inmax),INpz(wwmax,inmax)
      real*8 INmass(wwmax,inmax)
      integer INityp(wwmax,inmax),INiso3(wwmax,inmax)
      integer INch(wwmax,inmax),INlcoll(wwmax,inmax)
      integer INcoll(wwmax,inmax),INistr(wwmax,inmax)
      integer INorigin(wwmax,inmax)

c.. OUTgoing particles
      integer OUTind(wwmax,outmax)
      real*8 OUTr0(wwmax,outmax),OUTrx(wwmax,outmax)
      real*8 OUTry(wwmax,outmax),OUTrz(wwmax,outmax)
      real*8 OUTp0(wwmax,outmax),OUTpx(wwmax,outmax)
      real*8 OUTpy(wwmax,outmax),OUTpz(wwmax,outmax)
      real*8 OUTmass(wwmax,outmax)
      integer OUTityp(wwmax,outmax),OUTiso3(wwmax,outmax)
      integer OUTch(wwmax,outmax),OUTlcoll(wwmax,outmax)
      integer OUTcoll(wwmax,outmax),OUTistr(wwmax,outmax)
      integer OUTorigin(wwmax,outmax)

      integer www

c.. REConstructable resonances, store particle vector
      real*8 RECr0(resomax),RECrx(resomax)
      real*8 RECry(resomax),RECrz(resomax)
      real*8 RECp0(resomax),RECpx(resomax)
      real*8 RECpy(resomax),RECpz(resomax)
      real*8 RECmass(resomax)
      integer RECityp(resomax),RECiso3(resomax)
      integer RECch(resomax),REClcoll(resomax)
      integer RECcoll(resomax),RECorigin(resomax)
      integer RECind(resomax,0:4)

      real*8 ORIr0(resomax),ORIrx(resomax)
      real*8 ORIry(resomax),ORIrz(resomax)
      real*8 ORIp0(resomax),ORIpx(resomax)
      real*8 ORIpy(resomax),ORIpz(resomax)
      real*8 ORImass(resomax)
      integer ORIityp(resomax),ORIiso3(resomax)
      integer ORIch(resomax),ORIlcoll(resomax)
      integer ORIcoll(resomax),ORIorigin(resomax)

      integer resctr

      common /historyhead/hisnin,hisnexit,hisiline,hisctag,
     &        hisacttime,hissqrts,hisstot,hissigpart,hiscdens
      common /historyin/INind,INr0,INrx,INry,INrz,
     &        INp0,INpx,INpy,INpz,INmass,INityp,
     &        INiso3,INch,INlcoll,INcoll,INistr,INorigin
      common /historyout/OUTind,OUTr0,OUTrx,OUTry,OUTrz,
     &        OUTp0,OUTpx,OUTpy,OUTpz,OUTmass,OUTityp,
     &        OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTistr,OUTorigin
      common /historyrec/resctr,www,RECr0,RECrx,RECry,RECrz,
     &        RECp0,RECpx,RECpy,RECpz,RECmass,RECityp,
     &        RECiso3,RECch,REClcoll,RECcoll,RECorigin,
     &        RECind,
     &        ORIr0,ORIrx,ORIry,ORIrz,
     &        ORIp0,ORIpx,ORIpy,ORIpz,ORImass,ORIityp,
     &        ORIiso3,ORIch,ORIlcoll,ORIcoll,ORIorigin
