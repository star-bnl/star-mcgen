      OPTIONS/ EXTEND_SOURCE
C
      SUBROUTINE ptconv(W,Y)
C
      IMPLICIT NONE
C     >> Declare Global Variables
      INCLUDE 'global_var.inc'
C     >> Declare Arguments
      DOUBLE PRECISION W,Y
C     >> Declare Local Variables
      DOUBLE PRECISION px,py,pymax,pxmax,dx,dy,px0,py0
      DOUBLE PRECISION pt,pt1,pt2,q1,q2
      DOUBLE PRECISION W,Y,Egamma,Epom,sum,sumy
      DOUBLE PRECISION formf,f1,f2,norm
      INTEGER          K,I,J,Nxbin,Nybin,ybin1,ybin2
      REAL             xfill,yfill,xlo,xhi 
      DOUBLE PRECISION xg(1:8),ag(1:8)
      INTEGER NGAUSS
      DATA xg/.0950125098376374402D0,.281603550779258913D0,
     $        .458016777657227386D0, .617876244402643748D0,
     $        .755404408355003034D0, .865631202387831744D0,
     $        .944575023073232576D0, .989400934991649933D0/
      DATA ag/.189450610455068496D0, .182603415044923589D0,
     $        .169156519395002538D0, .149595988816576732D0,
     $        .124628971255533872D0, .0951585116824927848D0,
     $        .0622535239386478929D0,.0271524594117540949D0/
      NGAUSS=8
C      DATA xg/.0483076656877383162D0,.144471961582796493D0,
C     $        .239287362252137075D0, .331868602282127650D0,
C     $        .421351276130635345D0, .506899908932229390D0,
C     $        .587715757240762329D0, .663044266930215201D0,
C     $        .732182118740289680D0, .794483795967942407D0,
C     $        .849367613732569970D0, .896321155766052124D0,
C     $        .934906075937739689D0, .964762255587506430D0,
C     $        .985611511545268335D0, .997263861849481564D0/
C      DATA ag/.0965400885147278006D0, .0956387200792748594D0,
C     $        .0938443990808045654D0, .0911738786957638847D0,
C     $        .0876520930044038111D0, .0833119242269467552D0,
C     $        .0781938957870703065D0, .0723457941088485062D0,
C     $        .0658222227763618468D0, .0586840934785355471D0,
C     $        .0509980592623761762D0, .0428358980222266807D0,
C     $        .0342738629130214331D0, .0253920653092620595D0,
C     $        .0162743947309056706D0, .00701861000947009660D0/
C      NGAUSS=16

C     >> Initialize 
      pxmax = 10.*(hbarc/Rnuc) 
      pymax = 10.*(hbarc/RNuc)
      Nxbin = 1000
      Nybin = 500
      dx = 2.*pxmax/DFLOAT(Nxbin)
      dy = pymax/DFLOAT(Nybin)
C      WRITE(*,*) 'Integration Limit, Nbins   (X): ',pxmax,Nxbin
C      WRITE(*,*) 'Integration Limit, Nbins   (Y): ',pymax,Nybin

C     >> FIRST DISTRIBUTION
      Egamma = 0.5*W*DEXP(Y)
      Epom   = 0.5*W*DEXP(-Y)
C      WRITE(*,*) 'Egamma,Epom=',Egamma,Epom

C     >> Loop over total Pt to find distribution
      DO 100 K=1,100
C        pt = 0.000625*DFLOAT(K) - 0.0003125
        pt = 0.0025*DFLOAT(K) - 0.00125
        px0 = pt
        py0 = 0.0
C        WRITE(*,*) 'K,pt: ',K,pt

C       >> For each total Pt, integrate over Pt of one of the sources
        sum=0.
        DO 101 I=1,Nxbin
          px = -pxmax + (DFLOAT(I)-0.5)*dx 

          sumy=0.0
          DO 102 J=1,NGAUSS

            py = 0.5*pymax*xg(J)+0.5*pymax
            pt1 = DSQRT( px*px + py*py )
            pt2 = DSQRT( (px-px0)*(px-px0) + (py-py0)*(py-py0) )
            q1  = DSQRT( (Egamma/gamma_ta)**2 + pt1**2 )
            q2  = DSQRT( (Epom/gamma_ta)**2   + pt2**2 )
            f1  = (formf(q1)*formf(q1)*pt1*pt1)/(q1*q1*q1*q1)
            f2  = formf(q2)*formf(q2)
            sumy= sumy + ag(J)*f1*f2

            py = 0.5*pymax*(-xg(J))+0.5*pymax
            pt1 = DSQRT( px*px + py*py )
            pt2 = DSQRT( (px-px0)*(px-px0) + (py-py0)*(py-py0) )
            q1  = DSQRT( (Egamma/gamma_ta)**2 + pt1**2 )
            q2  = DSQRT( (Epom/gamma_ta)**2   + pt2**2 )
            f1  = (formf(q1)*formf(q1)*pt1*pt1)/(q1*q1*q1*q1)
            f2  = formf(q2)*formf(q2)
            sumy= sumy + ag(J)*f1*f2

 102      CONTINUE
C         >> This is to normalize the gaussian integration
          sumy = 0.5*pymax*sumy
C         >> The 2 is to account for py: 0 -- pymax
          sum  = sum + 2.*sumy*dx

 101    CONTINUE

        IF(K.eq.1)norm=1./sum
        ptparamy1(K)=sum*norm

 100  CONTINUE

C     >> Second Distribution
      Egamma = 0.5*W*DEXP(-Y)
      Epom   = 0.5*W*DEXP(Y)
C      WRITE(*,*) 'Egamma,Epom=',Egamma,Epom

C     >> Loop over total Pt to find distribution
      DO 200 K=1,100
C        pt = 0.000625*DFLOAT(K) - 0.0003125
        pt = 0.0025*DFLOAT(K) - 0.00125
        px0 = pt
        py0 = 0.0
C        WRITE(*,*) 'K,pt: ',K,pt

C       >> For each total Pt, integrate over Pt of one of the sources
        sum=0.
        DO 201 I=1,Nxbin
          px = -pxmax + (DFLOAT(I)-0.5)*dx 

          sumy=0.0
          DO 202 J=1,NGAUSS

            py = 0.5*pymax*xg(J)+0.5*pymax
            pt1 = DSQRT( px*px + py*py )
            pt2 = DSQRT( (px-px0)*(px-px0) + (py-py0)*(py-py0) )
            q1  = DSQRT( (Egamma/gamma_ta)**2 + pt1**2 )
            q2  = DSQRT( (Epom/gamma_ta)**2   + pt2**2 )
            f1  = (formf(q1)*formf(q1)*pt1*pt1)/(q1*q1*q1*q1)
            f2  = formf(q2)*formf(q2)
            sumy= sumy + ag(J)*f1*f2

            py = 0.5*pymax*(-xg(J))+0.5*pymax
            pt1 = DSQRT( px*px + py*py )
            pt2 = DSQRT( (px-px0)*(px-px0) + (py-py0)*(py-py0) )
            q1  = DSQRT( (Egamma/gamma_ta)**2 + pt1**2 )
            q2  = DSQRT( (Epom/gamma_ta)**2   + pt2**2 )
            f1  = (formf(q1)*formf(q1)*pt1*pt1)/(q1*q1*q1*q1)
            f2  = formf(q2)*formf(q2)
            sumy= sumy + ag(J)*f1*f2

 202      CONTINUE
C         >> This is to normalize the gaussian integration
          sumy = 0.5*pymax*sumy
C         >> The 2 is to account for py: 0 -- pymax
          sum  = sum + 2.*sumy*dx

 201    CONTINUE

        IF(K.eq.1)norm=1./sum
        ptparam2(K)=sum*norm

 200  CONTINUE

      ybin1 = INT(5.*(9.5+y))+1
      ybin2 = INT(5.*(9.5-y))+1
      WRITE(*,*) 'y,ybin1,ybin2:',y,ybin1,ybin2
      WRITE(*,*) 'dn/dy(1),dn/dy(2):',dsigdy(ybin1),dsigdy(ybin2)

      DO 300 K=1,100

        pt = 0.000625*DFLOAT(K) - 0.0003125
C        pt = 0.0025*DFLOAT(K) - 0.00125
        xfill = pt
        yfill = ptparam1(K)
        CALL HFILL(901,xfill,0.0,yfill)
        yfill = ptparam2(K)
        CALL HFILL(902,xfill,0.0,yfill)
        yfill = 0.5*(ptparam1(K)+ptparam2(K))
        CALL HFILL(903,xfill,0.0,yfill)
        yfill = ptparam1(K)*dsigdy(ybin1) + ptparam2(K)*dsigdy(ybin2)
        yfill = yfill/(dsigdy(ybin1)+dsigdy(ybin2))
        CALL HFILL(904,xfill,0.0,yfill)

 300  CONTINUE

      RETURN
      END
