c-*-Fortran-*-
c
c     This file contains the common blocks from VENUS_4.12
c
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      dimension nrptla(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      double precision seedf,seedc
      common /cseed/seedc /cptlu/nptlu
