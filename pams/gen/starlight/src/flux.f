      DOUBLE PRECISION FUNCTION flux(Egamma)


C  This routine gives the photon flux as a function of energy Egamma
C  It works for arbitrary nuclei and gamma; the first time it is
C  called, it calculates a lookup table which is used on 
C  subsequent calls

C it returns dN_gamma/dE (dimensions 1/E), not dI/dE
C energies are in GeV, in the lab frame

C  rewritten 4/25/2001 by SRK

      IMPLICIT NONE

      include 'D2LParam.inc'
      include 'const.inc'
      include 'inputp.inc'

      DOUBLE PRECISION Egamma,lEgamma,Emin,Emax,lnEmin,lnEmax
      DOUBLE PRECISION stepmult,energy,rZ,rmult,beta
      INTEGER j,nbstep,jb,nrstep,jr,nphistep,jphi,nstep
      DOUBLE PRECISION bmin,bmax,bmult,biter,bold,integratedflux
      DOUBLE PRECISION fluxelement,rmin,deltar,riter
      DOUBLE PRECISION deltaphi,phiiter,dist
      DOUBLE PRECISION dide(1:400),dlnE,lnElt

      DOUBLE PRECISION Xvar,Dbesk1,Dbesk0

      DOUBLE PRECISION binc,ra,phad,P1N,PXN,omax

      INTEGER Icheck,I,Ilt

      INTEGER ILOW,IHIGH
      DOUBLE PRECISION deltab,brange,singlep,bin

      DOUBLE PRECISION gammatarg,rhad,pnohad

C  This holds the 'chosen' table
      DOUBLE PRECISION prob(1070),b1(1070)

      SAVE Icheck,dide,lnEMax,lnEmin,dlnE

      DATA Icheck/0/

C   first call?  - initialize - calculate photon flux

      Icheck=Icheck+1
      IF(Icheck.gt.1)GOTO 1000

      rz=float(Z)
      ra=float(A)

C  is nuclear breakup considered?

      if (ibreakup .eq. 1) then
         write(6,111)
 111     format(' no nuclear breakup criteria')
         goto 220
      endif


C  !!!!Breakup is considered!!!!!


C  determine hadronic breakup probability from hadronbreakup.f
C  determine EM breakup probability from photonbreakup.f

C  start lookup table at bmin=2R_A; multiplicative steps of 1.1%
      bmin=2.*RNuc
      rmult=0.01
      binc=10.**rmult

C  If we allow for hadronic breakup, try a slightly smaller
C  bmin

      IF (ibreakup .NE. 1) bmin=0.95*bmin


C  different photon energy cutoffs for RHIC & LHC
C  in MeV

      omax=1.E7
      if (gamma_em .gt. 500) omax=1.E10


C  find gamma in target system

      gammatarg=2.*gamma_em*gamma_em-1.

C  both nuclei break up (XnXn)

           if (ibreakup .eq. 2) THEN
              write(6,213)bmin
 213      FORMAT(' Requiring both nuclei to break up (XnXn). ',
     *'Coulomb Only. bmin= ',F7.3)

              do 203 j=1,996
                 bmin=bmin*binc
                 call hadronbreakup(rhad,bmin,rz,ra,gamma_em)
                 call photonbreakup(P1N,PXN,bmin,rz,ra,gammatarg,omax)
                 b1(j)=bmin
                 pnohad=dexp(-rhad)
 203             prob(j)= (1.-dexp(-PXN))*(1-dexp(-PXN))*pnohad
              goto 220

           ENDIF

C ibreakup=3  --> both nuclear emit 1 neutron  (1n1n)

           if (ibreakup .eq. 3) THEN

              WRITE(6,214)bmin
 214   FORMAT(' Req. both nuclei to emit 1 neutron (1n1n) . bmin = '
     *,F7.3)

              do 204 j=1,996
                 bmin=bmin*binc
                 call hadronbreakup(rhad,bmin,rz,ra,gamma_em)
                 call photonbreakup(P1N,PXN,bmin,rz,ra,gammatarg,omax)
                 b1(j)=bmin
                 pnohad=dexp(-rhad)                 

C  old version
C 204            prob(j)= (1.-exp(-P1N))*(1-exp(-P1N))*pnohad

C following Eq. 7 of Baltz, Klein & Nystrand
C the 1n excitation cannot be accompanied by Xn
C 204             prob(j)=(exp(P1N-PXN)-exp(-PXN))**2*pnohad

C  Version following the 'old' version of Joakim's Erice talk

C  204             prob(j)=(dexp(P1N)-1.)*dexp(-2.*PXN)*pnohad

C Following the 'new' version of Joakim's Erice talk

  204             prob(j)=(P1N**2)*dexp(-2.*PXN)*pnohad

                 goto 220
           ENDIF

C ibreakup=4 --> neither nucleus is excited

           if (ibreakup .eq. 4) THEN
              do 205 j=1,996
                 bmin=bmin*binc
                 call hadronbreakup(rhad,bmin,rz,ra,gamma_em)
                 call photonbreakup(P1N,PXN,bmin,rz,ra,gammatarg,omax)
                 b1(j)=bmin
                 pnohad=dexp(-rhad)
 205             prob(j)= pnohad*dexp(-PXN)*dexp(-PXN)
                WRITE(6,215)
 215            FORMAT(' Requiring that neither nucleus break up.')
                GOTO 220
           ENDIF

C  ibreakup=5 - no hadronic breakup

           if (ibreakup .eq. 5) THEN

              WRITE(6,224)bmin
 224   FORMAT(' Requiring no hadronic breakup. bmin = '
     *,F7.3)

              do 225 j=1,996
                 bmin=bmin*binc
                 call hadronbreakup(rhad,bmin,rz,ra,gamma_em)
                 b1(j)=bmin
                 pnohad=dexp(-rhad)                 
 225             prob(j)=pnohad
                goto 220
           ENDIF

           WRITE(6,227)ibreakup
 227       FORMAT(' Ibreakup =',I7,'. Not understood')
           STOP

 220       continue

C  collect number of integration steps here, in one place

          nbstep=400
          nrstep=60
          nphistep=40

C  this last one is the number of energy steps
         nstep=100


C  following previous choices, take Emin=10 keV at LHC, Emin = 1 MeV at RHIC

        Emin=1.E-5
        if (gamma_em .lt. 500) Emin=1.E-3

C  maximum energy is 12 times the cutoff
C  25 GeV for gold at RHIC, 650 GeV for lead at LHC

        Emax=12.*hbarc*gamma_em/RNuc

C     >> lnEmin <-> ln(Egamma) for the 0th bin
C     >> lnEmax <-> ln(Egamma) for the last bin

        lnEmin=DLOG(Emin)
        lnEmax=DLOG(Emax)
        dlnE=(lnEmax-lnEmin)/nstep

        write(6,5)Emin,Emax
 5      format(' flux.f.  Calculating flux for photon energies', 
     &'from E=',E11.3,' to ',E11.3,' GeV (lab frame).')


        stepmult= dexp(dlog(Emax/Emin)/float(nstep))
        energy=Emin

        do 100 j=1,nstep
           energy=energy*stepmult

C  integrate flux over 2R_A < b < 2R_A+ 6* gamma hbar/energy 
C  use exponential steps

          bmin=2.*RNuc
          bmax=bmin + 6.*hbarc*gamma_em/energy

          bmult=dexp(dlog(bmax/bmin)/float(nbstep))
          biter=bmin
          integratedflux=0.

          do 90 jb=1,nbstep

             bold=biter
             biter=biter*bmult


C  When we get to b>20R_A change methods - just take the photon flux
C  at the center of the nucleus.


             if (biter .gt. 10.*RNuc) THEN

C if there is no nuclear breakup or only hadronic breakup, which only
C occurs at smaller b, we can analytically integrate the flux from b~20R_A
C to infinity, following Jackson (2nd edition), Eq. 15.54

                   Xvar=energy*biter/(hbarc*gamma_em)

C                if (ibreakup .eq. 1 .or. ibreakup .eq. 5) then

C                   beta=dsqrt(1.-1./gamma_em**2)
           
C                   fluxelement=(2.0/pi)*rZ*rZ*alpha/energy*
C     * (Xvar*Dbesk0(Xvar)*Dbesk1(Xvar) -
C     *beta*beta/2*Xvar*Xvar*(Dbesk1(Xvar)**2-Dbesk0(Xvar)**2))
           
C                   integratedflux=integratedflux+fluxelement

C branch out of the b loop.  We're done at this photon energy

C                   GOTO 95
C                ENDIF

C  Here, there is nuclear breakup.  So, we can't use the integrated flux
C  However, we can do a single flux calculation, at the center of the
C  nucleus

C  Eq. 41 of Vidovic, Greiner and Soff, Phys. Rev. C47, 2308 (1993), among other places
C  this is the flux per unit area

                fluxelement  = (rZ*rZ*alpha*energy)*(Dbesk1(Xvar))**2/
     *((pi*gamma_em*hbarc)**2)
           
             ELSE

C integrate over nuclear surface. n.b. this assumes total shadowing - 
C treat photons hitting the nucleus the same no matter where they strike

                fluxelement=0.
                deltar=RNuc/float(nrstep)
                riter=-deltar/2.

                do 80 jr=1,nrstep
                   riter=riter+deltar

C use symmetry;  only integrate from 0 to pi (half circle)

                   deltaphi=pi/float(nphistep)
                   phiiter=0.

                   do 70 jphi=1,nphistep
                      phiiter=(float(jphi)-0.5)*deltaphi
                   
C  dist is the distance from the center of the emitting nucleus to the point in question

                      dist=sqrt((biter+riter*cos(phiiter))**2+
     *(riter*sin(phiiter))**2)

                      Xvar=energy*dist/(hbarc*gamma_em)

                      flux = (rZ*rZ*alpha*energy)*(Dbesk1(Xvar))**2/
     *((pi*gamma_em*hbarc)**2)
                      
C  The surface  element is 2.* delta phi* r * delta r
C  The '2' is because the phi integral only goes from 0 to pi

         fluxelement=fluxelement+flux*2.*deltaphi*riter*deltar

C  end phi and r integrations
                     
 70                continue
 80             continue

C  average fluxelement over the nuclear surface

                fluxelement=fluxelement/(pi*RNuc**2)

C                 write(6,85)biter,fluxelement,Xvar
C  85    format(' b<2R_A = ',F7.3,' fluxelement= ',E11.3,' Xvar= ',E11.3)

             ENDIF

C  multiply by volume element to get total flux in the volume element

             fluxelement=fluxelement*2.*pi*biter*(biter-bold)

C  modulate by the probability of nuclear breakup as f(biter)

             IF (ibreakup .gt. 1) THEN

C This was the original line - I don't think the 0.5 belongs here
C since we do the linear interpolation
C                 bin = 0.5 + DLOG(biter/10.)/(rmult*LOG(10.))
                  bin = DLOG(biter/bmin)/(rmult*LOG(10.))
                 ILOW = INT(bin)
                 IHIGH = ILOW + 1
                 brange = b1(IHIGH)-b1(ILOW)
                 deltab = biter - b1(ILOW)        
                 singlep = (1-(deltab/brange))*prob(ILOW) + 
     &(deltab/brange)*prob(IHIGH)
                 fluxelement=fluxelement*singlep
              ENDIF

             integratedflux=integratedflux+fluxelement

 90          continue
C  end energy integration

 95          continue

C  In lookup table, store k*dN/dk because it changes less
C  so the interpolation should be better


             dide(j)=integratedflux*energy

C             write(6,97)energy,dide(j)
C 97          format(' Energy,flux= ',E11.3,'  ',E11.3)


 100         continue
                   
C  for 2nd and subsequent calls, use lookup table immediately



 1000   lEgamma=DLOG(Egamma)
        IF (lEgamma.lt.(lnEmin+dlnE).or. lEgamma .gt.lnEmax)THEN
            WRITE(6,1005)Egamma
 1005   FORMAT(' ERROR: Egamma outside defined range. Egamma= ',E11.3)
            flux=0.0
        ELSE
C       >> Egamma between Ilt and Ilt+1

            Ilt = INT((lEgamma-lnEmin)/dlnE)

C       >> ln(Egamma) for first point

            lnElt = lnEmin + Ilt*dlnE

C       >> Interpolate

            flux = dide(Ilt) + ((lEgamma-lnElt)/dlnE)*
     &(dide(Ilt+1)- dide(Ilt))
            flux = flux/Egamma
      ENDIF
      RETURN
      END

