C      OPTIONS/ EXTEND_SOURCE
C
      SUBROUTINE VMSIGMAPT(W,Egamma,SIGMAPT)

C  This subroutine calculates the effect of the nuclear form factor
C  on the pt spectrum, for use in interference calculations
C  For an interaction with mass W and photon energy Egamma,
C  it calculates the cross section suppression SIGMAPT(PT)
C  as a function of pt.

C  The input pt values come from pttable.inc

      IMPLICIT NONE

      INCLUDE 'const.inc'
      INCLUDE 'inputp.inc'
      INCLUDE 'D2LParam.inc'
      include 'pttable.inc'

      DOUBLE PRECISION SIGMAPT(500)

      DOUBLE PRECISION px,py,pymax,pxmax,dx,px0,py0
      DOUBLE PRECISION pt,pt1,pt2,q1,q2
      DOUBLE PRECISION W,Egamma,Epom,sum,sumy
      DOUBLE PRECISION formf,f1,f2,norm
      INTEGER          K,I,J,Nxbin,ybin1,ybin2
      DOUBLE PRECISION xg(1:16),ag(1:16)
      INTEGER NGAUSS

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

C     >> Initialize 
      pxmax = 10.*(hbarc/Rnuc) 
      pymax = 10.*(hbarc/RNuc)

      Nxbin = 500

      dx = 2.*pxmax/DFLOAT(Nxbin)

C      WRITE(*,*) 'Integration Limit, Nbins   (X): ',pxmax,Nxbin
C      WRITE(*,*) 'Integration Limit, Nbins   (Y): ',pymax,Nybin

      Epom   = W*W/(4.*Egamma)
C     WRITE(*,*) 'Egamma,Epom=',Egamma,Epom

C     >> Loop over total Pt to find distribution

      DO 100 K=1,NPT

        pt=dpt*(DFLOAT(K)-0.5)

        px0 = pt
        py0 = 0.0
C       WRITE(*,*) 'K,pt: ',K,pt

C  For each total Pt, integrate over Pt1, , the photon pt
C  The pt of the Pomeron  is the difference
C  pt1 is  
        sum=0.
        DO 101 I=1,Nxbin
          px = -pxmax + (DFLOAT(I)-0.5)*dx 

          sumy=0.0
          DO 102 J=1,NGAUSS

            py = 0.5*pymax*xg(J)+0.5*pymax
C  photon pt
            pt1 = DSQRT( px*px + py*py )
C  pomeron pt
            pt2 = DSQRT( (px-px0)*(px-px0) + (py-py0)*(py-py0) )
            q1  = DSQRT( (Egamma/gamma_em)**2 + pt1**2 )
            q2  = DSQRT( (Epom/gamma_em)**2   + pt2**2 )

C  photon form factor

C add in phase space factor?

            f1  = (formf(q1*q1)*formf(q1*q1)*pt1*pt1)/(q1*q1*q1*q1)

C  Pomeron form factor

            f2  = formf(q2*q2)*formf(q2*q2)
            sumy= sumy + ag(J)*f1*f2

C  now consider other half of py phase space - why is this split?

            py = 0.5*pymax*(-xg(J))+0.5*pymax
            pt1 = DSQRT( px*px + py*py )
            pt2 = DSQRT( (px-px0)*(px-px0) + (py-py0)*(py-py0) )
            q1  = DSQRT( (Egamma/gamma_em)**2 + pt1**2 )
            q2  = DSQRT( (Epom/gamma_em)**2   + pt2**2 )
C  add in phase space factor?
            f1  = (formf(q1*q1)*formf(q1*q1)*pt1*pt1)/(q1*q1*q1*q1)
            f2  = formf(q2*q2)*formf(q2*q2)
            sumy= sumy + ag(J)*f1*f2

 102      CONTINUE
C         >> This is to normalize the gaussian integration
          sumy = 0.5*pymax*sumy
C         >> The 2 is to account for py: 0 -- pymax
          sum  = sum + 2.*sumy*dx

 101    CONTINUE

        IF(K.eq.1)norm=1./sum
        SIGMAPT(K)=sum*norm

 100  CONTINUE

      RETURN
      END



