C      OPTIONS/ EXTEND_SOURCE

      SUBROUTINE PTTABLEGEN

C  Calculates the pt spectra for VM production with interference
C  Follows S. Klein and J. Nystrand, Phys. Rev Lett. 84, 2330 (2000).
C  Written by S. Klein, 8/2002

C  fill in table pttable in one call
C  Integrate over all y (using the same y values as in table yarray
C  note that the cross section goes from ymin (<0) to ymax (>0), in numy points
C  here,  we go from 0 to ymax in (numy/2)+1 points
C  numy must be even.

C  At each y, calculate the photon energies Egamma1 and Egamma2
C  and the two photon-A cross sections

C  loop over each p_t entry in the table.

C  Then, loop over b and phi (the angle between the VM \vec{p_t} and \vec{b} 
C  and calculate the cross section at each step.
C  Put the results in pttable

      IMPLICIT NONE

      include 'Ftable.inc'
      include 'pttable.inc'
      include 'inputp.inc'
      include 'const.inc'
      include 'global.inc'
      include 'D2LParam.inc'

      DOUBLE PRECISION nofe
      DOUBLE PRECISION b,bmin,bmax,db,sum1,sumg, sumint
      DOUBLE PRECISION pt,Egamma1,Egamma2,theta,dtheta
      DOUBLE PRECISION int1,sig_ga_1,sig_ga_2
      DOUBLE PRECISION A1,A2,amp_i_2
      DOUBLE PRECISION Wgp,cs,cvma,Av,csgA,t,tmin,tmax,formf,ax,bx
      DOUBLE PRECISION norm2,dY
      DOUBLE PRECISION sigmagp,sigma_A
      DOUBLE PRECISION xg(1:16),ag(1:16)
      INTEGER          I,J,K,NBIN,NThetaBIN,NGAUSS,jy
      DOUBLE PRECISION ptparam1(300),ptparam2(300),Yp
      REAL PofB


      DATA xg/.0483076656877383162D0,.144471961582796493D0,
     $        .239287362252137075D0, .331868602282127650D0,
     $        .421351276130635345D0, .506899908932229390D0,
     $        .587715757240762329D0, .663044266930215201D0,
     $        .732182118740289680D0, .794483795967942407D0,
     $        .849367613732569970D0, .896321155766052124D0,
     $        .934906075937739689D0, .964762255587506430D0,
     $        .985611511545268335D0, .997263861849481564D0/
      DATA ag/.0965400885147278006D0, .0956387200792748594D0,
     $        .0938443990808045654D0, .0911738786957638847D0,
     $        .0876520930044038111D0, .0833119242269467552D0,
     $        .0781938957870703065D0, .0723457941088485062D0,
     $        .0658222227763618468D0, .0586840934785355471D0,
     $        .0509980592623761762D0, .0428358980222266807D0,
     $        .0342738629130214331D0, .0253920653092620595D0,
     $        .0162743947309056706D0, .00701861000947009660D0/
      NGAUSS=16

C  loop over y from 0 (not -ymax) to ymax
 
      WRITE(6,11)Ymax,numy/2
 11   FORMAT(' PTTABLEGEN. Integrating Y from 0 --> ',F7.4,
     *' in ', I6,' steps.')

      dpt=ptmax/DFLOAT(NPT)

      write(6,12)ptmax,NPT
 12   format(' Scanning pt from 0-> ',F8.4,
     *' in ',I5,' steps.')


      dY=(2.*Ymax)/numy
      DO 1000 jy=1,numy/2
         Yp=(DFLOAT(jy)-0.5)*dY

C  Find the photon energies.  Yp >= 0, so Egamma2 is smaller

C  Use the vector meson mass for W here - neglect the width

         Egamma1 = 0.5*mass*DEXP(Yp)
         Egamma2 = 0.5*mass*DEXP(-Yp)

C  Find the sigma(gammaA) for the two directions

C  Photonuclear Cross Section 1

C  Gamma-proton CM energy

        Wgp=DSQRT(2.*Egamma1*(Ep+DSQRT(Ep*Ep-mp*mp))+mp*mp)

C Calculate V.M.+proton cross section

        cs=DSQRT(16.*pi*f2o4pi*bslope*hbarc*hbarc*sigmagp(Wgp)/alpha)

C Calculate V.M.+Nucleus cross section

        cvma=sigma_A(cs)

C Calculate Av = dsigma/dt(t=0) Note Units: fm**2/Gev**2

        Av=(alpha*cvma*cvma)/(16.*pi*f2o4pi*hbarc*hbarc)

        tmin  = ( (mass**2)/(4.*Egamma1*gamma_em) )**2
        tmax  = tmin + 0.25
        ax    = 0.5*(tmax-tmin)
        bx    = 0.5*(tmax+tmin)
        csgA  = 0.
        DO 112 K=1,NGAUSS
          t     = DSQRT(ax*xg(K)+bx)
          csgA  = csgA + ag(K)*formf(t)*formf(t)
          t     = DSQRT(ax*(-xg(K))+bx)
          csgA  = csgA + ag(K)*formf(t)*formf(t)
 112    CONTINUE
        csgA = 0.5*(tmax-tmin)*csgA
        csgA = Av*csgA
      sig_ga_1 = csgA

C Photonuclear Cross Section 2

        Wgp=DSQRT(2.*Egamma2*(Ep+DSQRT(Ep*Ep-mp*mp))+mp*mp)

        cs=DSQRT(16.*pi*f2o4pi*bslope*hbarc*hbarc*sigmagp(Wgp)/alpha)

        cvma=sigma_A(cs)

        Av=(alpha*cvma*cvma)/(16.*pi*f2o4pi*hbarc*hbarc)

        tmin  = ( (mass**2)/(4.*Egamma2*gamma_em) )**2
        tmax  = tmin + 0.25
        ax    = 0.5*(tmax-tmin)
        bx    = 0.5*(tmax+tmin)
        csgA  = 0.
        DO 122 K=1,NGAUSS
          t     = DSQRT(ax*xg(K)+bx)
          csgA  = csgA + ag(K)*formf(t)*formf(t)
          t     = DSQRT(ax*(-xg(K))+bx)
          csgA  = csgA + ag(K)*formf(t)*formf(t)
 122    CONTINUE
        csgA = 0.5*(tmax-tmin)*csgA
        csgA = Av*csgA
      sig_ga_2 = csgA

C       WRITE(*,*) 'Y,Egamma1,Egamma2: ',Yp,Egamma1,Egamma2
C       WRITE(*,*) 'sig1,sig2        : ',sig_ga_1,sig_ga_2

C  Set up pttables - they find the reduction in sigma(pt)
C  due to the nuclear form factors.

C  Use the vector meson mass for W here - neglect width in
C  interference calculation

      CALL VMSIGMAPT(mass,Egamma1,ptparam1)
      CALL VMSIGMAPT(mass,Egamma2,ptparam2)

C  set  bmax according to the smaller photon energy, following flux.f

         bmax=bmin+6.*hbarc*gamma_em/Egamma2

         bmin = 2.*RNuc
C  if we allow for nuclear breakup, use a slightly smaller bmin

         if (ibreakup .ne. 1) bmin=0.95*bmin


C  set number of bins to a reasonable number to start

      NBIN = 2000
      NThetaBIN = 1000

      db   = (bmax-bmin)/DFLOAT(NBIN)

C loop over pt

         DO 101 I=1,NPT

            pt = (DFLOAT(I)-0.5)*dpt

           sum1=0.0
 
C loop over b

           DO 102 J=1,NBIN
  
                b = bmin + (DFLOAT(J)-0.5)*db

C  nofe is the photon flux function

               A1 = Egamma1*nofe(Egamma1,b)*sig_ga_1*ptparam1(I)
               A2 = Egamma2*nofe(Egamma2,b)*sig_ga_2*ptparam2(I)

               sumg=0.0

C  do this as a Gaussian integral, from 0 to pi

               DO 103 K=1,NGAUSS

                    theta=xg(K)*pi

C  allow for a linear sum of interfering and non-interfering amplitudes

                    amp_i_2 = A1 + A2 
     *- 2.*xinterfere*DSQRT(A1*A2)*DCOS(pt*b*DCOS(theta)/hbarc)

                    sumg  = sumg+ag(K)*amp_i_2

 103      CONTINUE

C  this is dn/dpt^2

C  The factor of 2 is because the theta integral is only from 0 to pi

          sumint=2.*sumg*b*db
          if (ibreakup .gt. 1) sumint=sumint*PofB(b)
          sum1 = sum1 + sumint
 
C end of b loop
 102    CONTINUE

C  normalization is done in readDiffLum.f
C  This is d^2sigma/dpt^2; convert to dsigma/dpt

      write(20,*)sum1*pt*dpt

C  end of pt loop
 101  CONTINUE

C  end of y loop
1000  CONTINUE

      RETURN 
      END 






