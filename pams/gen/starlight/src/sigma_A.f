C     >> Nuclear Cross Section
C     >> sig_N,sigma_A in (fm**2)

      DOUBLE PRECISION FUNCTION sigma_A(sig_N)

      IMPLICIT NONE

      include 'global.inc'
      include 'const.inc'
      DOUBLE PRECISION sig_N
      DOUBLE PRECISION rho0,sum
      DOUBLE PRECISION b,bmax,Pint,t,arg
      DOUBLE PRECISION xg(1:16),ag(1:16)
      INTEGER Ib,NGAUSS

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

      bmax = 25.0
      sum  = 0.
C     >> CALCULATE P(int) FOR b=0.0 - bmax (fm)
      DO 100 IB=1,NGAUSS

        b = 0.5*bmax*xg(IB)+0.5*bmax

        arg=-sig_n*rho0*T(b)
        Pint=1.0-DEXP(arg)
        sum=sum+2.*pi*b*Pint*ag(IB)

        b = 0.5*bmax*(-xg(IB))+0.5*bmax

        arg=-sig_n*rho0*T(b)
        Pint=1.0-DEXP(arg)
        sum=sum+2.*pi*b*Pint*ag(IB)
c	write(*,*) 'T,xg,ag',T(b),xg(IB),ag(IB)

 100  CONTINUE
      sum=0.5*bmax*sum

      sigma_A=sum
C      WRITE(*,*) 'Cross Section (mb):',10.*sum

      RETURN
      END
