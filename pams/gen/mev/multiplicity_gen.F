      Subroutine track_gen(gpid,N,model_type,T,sigma,expvel,
     1     pt_cut_min,pt_cut_max,eta_cut_min,eta_cut_max,
     2     phi_cut_min,phi_cut_max,n_scan_pts,rad,pid,status)
      implicit none
#include "mevcdes/multi_gen.inc"
      Common/track/ pout(npid,3,factorial_max)
      real*4 pout
      Common/Geant/geant_mass(60)
      real*4 geant_mass

      real*4 T,sigma,expvel,pt_cut_min,pt_cut_max,eta_cut_min
      real*4 eta_cut_max,phi_cut_min,phi_cut_max,mass
      real*4 dpt,deta,facmax,pt,eta,y,rapidity,dNdpty,dNdp
      real*4 pt_trial,eta_trial,y_trial,rndm,rad,phi

      integer gpid,N,model_type,npt,neta,n_scan_pts,ipt,ieta
      integer Track_count,i,j,pid,status

      do i = 1,factorial_max
      do j = 1,3
         pout(pid,j,i) = 0.0
      end do
      end do

      mass = geant_mass(gpid)
      npt  = n_scan_pts
      neta = n_scan_pts

CCC  Determine maximum value of model distribution in (pt,eta) range:
    
      dpt = (pt_cut_max - pt_cut_min)/float(npt - 1)
      deta = (eta_cut_max - eta_cut_min)/float(neta - 1)
      facmax = 0.0
      do ipt = 1,npt
         pt = pt_cut_min + dpt*float(ipt - 1)
         do ieta = 1,neta
            eta = eta_cut_min + deta*float(ieta - 1)
            y   = rapidity(mass,pt,eta)
            dNdp = dNdpty(1.0,pt,eta,y,mass,T,sigma,expvel,
     1                    model_type,1,pid)
            if(dNdp .gt. facmax) facmax = dNdp
         end do
      end do

CCC   If dNdp always underflows exp() range, then facmax will stay 
CCC   equal to 0.0, thus causing a divide by 0.0 error below. 
CCC   Check for this and Return if this is the case.  This event will
CCC   be aborted in this instance.

      if(facmax .eq. 0.0) then
         status = -86
         Return
      else
         status = 1
      end if

CCC Start Random Track Selection:

      Track_count = 0

100   pt_trial = rndm(-1)*(pt_cut_max - pt_cut_min)+pt_cut_min
      eta_trial=rndm(-1)*(eta_cut_max-eta_cut_min)+eta_cut_min
      y_trial = rapidity(mass,pt_trial,eta_trial)
      dNdp = dNdpty(1.0,pt_trial,eta_trial,y_trial,mass,T,sigma,
     1       expvel,model_type,1,pid)/facmax 
      if(rndm(-1) .le. dNdp) then
         Track_count = Track_count + 1
         phi = (rndm(-1)*(phi_cut_max - phi_cut_min) + 
     1         phi_cut_min)/rad
         Call kinematics(mass,pt_trial,y_trial,phi,Track_count,pid)
         if(Track_count .lt. N) then
            go to 100
         else if(Track_count .ge. N) then
            Return
         end if

      else
         go to 100
      end if

      END

      Real*4 Function rapidity(m,pt,eta)
      implicit none
      real*4 m,pt,eta,theta,pz,E,pi,expR

      pi = 3.141592654
      theta = 2.0*ATAN(expR(-eta))
      if(abs(theta - pi/2.0) .lt. 0.000001) then
         pz = 0.0
      else
         pz = pt/tan(theta)
      end if
      E = sqrt(pt*pt + pz*pz + m*m)
      rapidity = 0.5*log((E+pz)/(E-pz))
      Return
      END

      Subroutine kinematics(m,pt,y,phi,index,pid)
      implicit none

CCC  This SUBR takes the input particle mass (m), pt, y and phi and
CCC  computes E, px, py, pz and loads the momenta into the index-th
CCC  row of array pout(,,) in Common/track/ .

      integer index,pid
#include "mevcdes/multi_gen.inc"
      Common/track/ pout(npid,3,factorial_max)
      real*4 pout

      real*4 m,pt,y,phi,mt,coshy,E,pzmag,px,py,pz,expR

      mt = sqrt(m*m + pt*pt)
      coshy = 0.5*(expR(y) + expR(-y))
      E = mt*coshy
      pzmag = sqrt(abs(E*E - mt*mt))
      if(y.eq.0.0) then
         pz = 0.0
      else
         pz = (y/abs(y))*pzmag
      end if
      px = pt*cos(phi)
      py = pt*sin(phi)

      pout(pid,1,index) = px
      pout(pid,2,index) = py
      pout(pid,3,index) = pz

      Return
      END

      Real*4 Function dNdpty(A,pt,eta,y,m,T,sigma,vel,control,ktl,pid)
      implicit none

      real*4 A,pt,eta,y,m,T,sigma,vel
      real*4 pi,mt,coshy,E,ptot,FAC,gamma,yp,sinhyp,coshyp
      real*4 FAC2,FAC3,expR
      integer control, ktl,pid,pt_index,eta_index,index_locate
#include "mevcdes/multi_gen.inc"
CCC   Calculates dN/dp^3 using several models:
CCC
CCC      control = 1,  Humanic factorized model
CCC              = 2,  Pratt non-expanding spherical thermal source
CCC              = 3,  Bertsch non-expanding spherical thermal source
CCC              = 4,  Pratt spherical expanding thermally equilibrated
CCC                    source.
CCC              = 5,  Factorized pt, eta bin-by-bin distribution.
CCC              = 6,  Full 2D pt, eta bin-by-bin distribution.
CCC
CCC      ktl     = 0,  to return value of dN/dp^3
CCC      ktl     = 1,  to return value of dN/dpt*dy

      pi = 3.141592654
      mt = sqrt(pt*pt + m*m)
      coshy = 0.5*(expR(y) + expR(-y))
      E = mt*coshy
      ptot = sqrt(E*E - m*m)
      if(ktl .eq. 0) then
         FAC = 2.0*pi*pt*E
      else if(ktl .eq. 1) then
         FAC = 1.0
      end if

      if(control .eq. 1) then
         dNdpty = A*pt*expR(-mt/T)*expR(-y*y/(2.0*sigma*sigma))
         dNdpty = dNdpty/FAC

      else if(control .eq. 2) then
         dNdpty = A*pt*E*expR(-E/T)
         dNdpty = dNdpty/FAC

      else if(control .eq. 3) then
         dNdpty = A*pt*E/(expR(E/T) - 1.0)
         dNdpty = dNdpty/FAC

      else if(control .eq. 4) then
         gamma = 1.0/sqrt(1.0 - vel*vel)
         yp = gamma*vel*ptot/T
         sinhyp = 0.5*(expR(yp) - expR(-yp))
         coshyp = 0.5*(expR(yp) + expR(-yp))
         if(yp .ne. 0.0) then
            FAC2 = sinhyp/yp
         else
            FAC2 = 1.0
         end if
         FAC3 = FAC2+ (T/(gamma*E))*(FAC2 - coshyp)
         dNdpty = A*pt*E*expR(-gamma*E/T)*FAC3
         dNdpty = dNdpty/FAC

      else if(control .eq. 5) then
         pt_index = index_locate(pid,pt,1)
         eta_index = index_locate(pid,eta,2)
         dNdpty = A*pt_bin(pid,pt_index)*eta_bin(pid,eta_index)
         dNdpty = dNdpty/FAC

      else if(control .eq. 6) then
         pt_index = index_locate(pid,pt,1)
         eta_index = index_locate(pid,eta,2)
         dNdpty = A*pt_eta_bin(pid,pt_index,eta_index)
         dNdpty = dNdpty/FAC

      else
         dNdpty = -86.0

      end if

      return
      END

      Integer Function index_locate(pid,arg,kind)
      implicit none

#include "mevcdes/multi_gen.inc"
      integer pid,kind,ibin
      real*4 arg

CCC   This Function locates the pt or eta bin number corresponding to the
CCC   input bin mesh, the requested value of pt or eta, for the current
CCC   value of particle type.
CCC
CCC   If kind = 1, then pt bin number is located.
CCC   If kind = 2, then eta bin number is located.

      if(kind .eq. 1) then
         do ibin = 1,n_pt_bins(pid)
            if(arg.le.pt_bin_mesh(pid,ibin)) then
            index_locate = ibin
            Return
            end if
         end do
         index_locate = n_pt_bins(pid)
         write(8,10) pid,arg
10       Format(//10x,'In Function index_locate, for pid = ',I5,
     1   'pt  =',E15.6,' is out of range - use last bin #')
         Return

      else if(kind .eq. 2) then

         do ibin = 1,n_eta_bins(pid)
            if(arg.le.eta_bin_mesh(pid,ibin)) then
            index_locate = ibin
            Return
            end if
         end do
         index_locate = n_eta_bins(pid)
         write(8,11) pid,arg
11       Format(//10x,'In Function index_locate, for pid = ',I5,
     1   'eta =',E15.6,' is out of range - use last bin #')
         Return

      end if
      END

      Real*4 Function expR(x)
      implicit none
      real*4 x
      if(x .gt. 69.0) then
         write(8,10) x
         STOP
      else if(x .lt. -69.0) then
         expR = 0.0
      else
         expR = exp(x)
      end if
10    Format(///10x,'Func. expR(x) called with x = ',E15.6,
     1    ', gt 69.0 - STOP')
      Return
      END

      SUBROUTINE LAGRNG (X,ARG,Y,VAL,NDIM,NFS,NPTS,MAXARG,MAXFS)
	IMPLICIT REAL*4(A-H,O-Z)
C
C     LAGRANGE INTERPOLATION,UNEQUALLY SPACED POINTS
C     ROUTINE OBTAINED FROM R. LANDAU, UNIV. OF OREGON.
C     ARG=VECTOR OF INDEPENDENT VARIABLE CONTAINING MAXARG VALUES.
C     VAL=MATRIX OF FUNCTION VALUES CORRESPONDING TO ARG. (MAXFS
C         FUNCTIONS AT MAXARG VALUES.)
C     X  =VALUE OF INDEP. VARIABLE FOR WHICH INTERPOLATION IS DESIRED.
C     Y  =VECTOR OF MAXFS FUNCTION VALUES RESULTING FROM SIMUL. INTERP.
C     NDIM=NUMBER OF ARG VALUES TO BE USED. (NDIM.LE.MAXARG)
C     NFS=NUMBER OF FUNCTIONS SIMUL. INTERP (NFS.LE.MAXFS)
C     NPTS=NUMBER OF POINTS USED IN INTERPOLATION. (NPTS=2,3,4,5,6)
C
      DIMENSION ARG(MAXARG), VAL(MAXFS,MAXARG), Y(MAXFS)
C
C     -----FIND X0, THE CLOSEST POINT TO X.
C
      NI=1
      NF=NDIM
   10 IF ((X.LE.ARG(NI)).OR.(X.GE.ARG(NF))) GO TO 30
      IF ((NF-NI+1).EQ.2) GO TO 70
      NMID=(NF+NI)/2
      IF (X.GT.ARG(NMID)) GO TO 20
      NF=NMID
      GO TO 10
   20 NI=NMID
      GO TO 10
C
C     ------ X IS ONE OF THE TABLULATED VALUES.
C
   30 IF (X.LE.ARG(NI)) GO TO 60
      NN=NF
   40 NUSED=0
      DO 50 N=1,NFS
   50 Y(N)=VAL(N,NN)
      RETURN
   60 NN=NI
      GO TO 40
C
C     ------- 2 PTS LEFT, CHOOSE SMALLER ONE.
C
   70 N0=NI
      NN=NPTS-2
      GO TO (110,100,90,80), NN
   80 CONTINUE
      IF (((N0+3).GT.NDIM).OR.((N0-2).LT.1)) GO TO 90
      NUSED=6
      GO TO 130
   90 CONTINUE
      IF ((N0+2).GT.NDIM) GO TO 110
      IF ((N0-2).LT.1) GO TO 100
      NUSED=5
      GO TO 130
  100 CONTINUE
      IF (((N0+2).GT.NDIM).OR.((N0-1).LT.1)) GO TO 110
      NUSED=4
      GO TO 130
  110 IF ((N0+1).LT.NDIM) GO TO 120
C
C     ------N0=NDIM, SPECIAL CASE.
C
      NN=NDIM
      GO TO 40
  120 NUSED=3
      IF ((N0-1).LT.1) NUSED=2
  130 CONTINUE
C
C     ------AT LEAST 2 PTS LEFT.
C
      Y0=X-ARG(N0)
      Y1=X-ARG(N0+1)
      Y01=Y1-Y0
      C0=Y1/Y01
      C1=-Y0/Y01
      IF (NUSED.EQ.2) GO TO 140
C
C     ------AT LEAST 3 PTS.
C
      YM1=X-ARG(N0-1)
      Y0M1=YM1-Y0
      YM11=Y1-YM1
      CM1=-Y0*Y1/Y0M1/YM11
      C0=C0*YM1/Y0M1
      C1=-C1*YM1/YM11
      IF (NUSED.EQ.3) GO TO 160
C
C     ------AT LEAST 4 PTS
C
      Y2=X-ARG(N0+2)
      YM12=Y2-YM1
      Y02=Y2-Y0
      Y12=Y2-Y1
      CM1=CM1*Y2/YM12
      C0=C0*Y2/Y02
      C1=C1*Y2/Y12
      C2=-YM1*Y0*Y1/YM12/Y02/Y12
      IF (NUSED.EQ.4) GO TO 180
C
C     ------AT LEAST 5 PTS.
C
      YM2=X-ARG(N0-2)
      YM2M1=YM1-YM2
      YM20=Y0-YM2
      YM21=Y1-YM2
      YM22=Y2-YM2
      CM2=YM1*Y0*Y1*Y2/YM2M1/YM20/YM21/YM22
      CM1=-CM1*YM2/YM2M1
      C0=-C0*YM2/YM20
      C1=-C1*YM2/YM21
      C2=-C2*YM2/YM22
      IF (NUSED.EQ.5) GO TO 200
C
C     ------AT LEAST 6 PTS.
C
      Y3=X-ARG(N0+3)
      YM23=Y3-YM2
      YM13=Y3-YM1
      Y03=Y3-Y0
      Y13=Y3-Y1
      Y23=Y3-Y2
      CM2=CM2*Y3/YM23
      CM1=CM1*Y3/YM13
      C0=C0*Y3/Y03
      C1=C1*Y3/Y13
      C2=C2*Y3/Y23
      C3=YM2*YM1*Y0*Y1*Y2/YM23/YM13/Y03/Y13/Y23
      GO TO 220
  140 CONTINUE
      DO 150 N=1,NFS
  150 Y(N)=C0*VAL(N,N0)+C1*VAL(N,N0+1)
      GO TO 240
  160 CONTINUE
      DO 170 N=1,NFS
  170 Y(N)=CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)
      GO TO 240
  180 CONTINUE
      DO 190 N=1,NFS
  190 Y(N)=CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)+C2*VAL(N,N0+2)
      GO TO 240
  200 CONTINUE
      DO 210 N=1,NFS
  210 Y(N)=CM2*VAL(N,N0-2)+CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)+C
     12*VAL(N,N0+2)
      GO TO 240
  220 CONTINUE
      DO 230 N=1,NFS
  230 Y(N)=CM2*VAL(N,N0-2)+CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)+C
     12*VAL(N,N0+2)+C3*VAL(N,N0+3)
  240 RETURN
C
      END

      Subroutine FACTORIAL
      implicit none

      integer n
#include "mevcdes/multi_gen.inc"
      real*4 FN

CCC   Computes the natural log of n! for n = 0,1,2...(factorial_max -1)
CCC   and puts the result in array FACLOG().
C
CCC   FACLOG(1) = log(0!) = 0.0
CCC   FACLOG(2) = log(1!) = 0.0
CCC   FACLOG(n+1) = log(n!)

      FACLOG(1) = 0.0
      FACLOG(2) = 0.0
      FN = 1.0
      do n = 3,factorial_max
      FN = FN + 1.0
      FACLOG(n) = FACLOG(n-1) + log(FN)
      end do
      Return
      END

      Subroutine MinMax(mean,stdev,n_stdev,min,max)
      implicit none

CCC   Computes range of integration for random number selections

      real*4 mean,stdev,n_stdev,min,max

      min = mean - n_stdev*stdev
      if(min .lt. 0.0) min = 0.0
      max = mean + n_stdev*stdev
      Return
      END

      Subroutine Poisson(min,max,mean,nsteps,integ,xfunc,ndim)
      implicit none

CCC   Computes Poisson distribution from n = min to max;
CCC   Integrates this distribution and records result at each step in
CCC      array integ();
CCC   Records the coordinates in array xfunc().

      integer min,max,mean,nsteps,ndim,i,n
#include "mevcdes/multi_gen.inc"
      real*4 mean_real,mean_real_log,expR
      real*4 integ(ndim)
      real*4 xfunc(ndim)
      real*4 Poisson_dist(n_mult_max_steps)
CCC Initialize arrays to zero:

      do i = 1,ndim
         integ(i) = 0.0
         xfunc(i) = 0.0
         Poisson_dist(i) = 0.0
      end do

      mean_real = float(mean)
      mean_real_log = log(mean_real) 

CCC   Compute Poisson distribution from n = min to max:
      do i = 1,nsteps
      n = (i - 1) + min
      Poisson_dist(i) = expR(-mean_real + float(n)*mean_real_log
     1      - FACLOG(n+1))
      end do

CCC   Integrate the Poisson distribution:
      integ(1) = 0.0
      do i = 2,nsteps
      integ(i) = 0.5*(Poisson_dist(i-1) + Poisson_dist(i)) + integ(i-1)
      end do

CCC   Normalize the integral to unity:
      do i = 1,nsteps
      integ(i) = integ(i)/integ(nsteps)
      end do

CCC   Fill xfunc array:
      do i = 1,nsteps
      xfunc(i) = float(i - 1 + min)
      end do

CCC  Extend integ() and xfunc() by one additional mesh point past the
CCC  end point in order to avoid a bug in the Lagrange interpolation
CCC  subroutine that gives erroneous interpolation results within the
CCC  the last mesh bin.

      integ(nsteps + 1) = integ(nsteps) + 0.01
      xfunc(nsteps + 1) = xfunc(nsteps)

      Return
      END

      Subroutine Gaussian(min,max,mean,stdev,npts,integ,xfunc,ndim)
      implicit none

CCC   Compute Gaussian distribution from x = min to max at npts;
CCC   Integrate this distribution and record result at each mesh in
CCC      array integ();
CCC   Record the coordinates in array xfunc().

#include "mevcdes/multi_gen.inc"
      integer npts,ndim,i
      real*4 min,max,mean,stdev,integ(ndim),xfunc(ndim)
      real*4 dm,x,Gauss_dist(nmax_integ),FAC1,FAC2,pi,expR

CCC   Initialize arrays to zero:
      do i = 1,ndim
         integ(i) = 0.0
         xfunc(i) = 0.0
         Gauss_dist(i) = 0.0
      end do

      pi = 3.141592654
      FAC1 = 1.0/(sqrt(2.0*pi)*stdev)
      FAC2 = 2.0*stdev*stdev
      dm = (max - min)/float(npts - 1)

CCC   Compute normalized Gaussian distribution:
      do i = 1,npts
      x = min + dm*float(i-1)
      xfunc(i) = x
      Gauss_dist(i) = FAC1*expR(-((x - mean)**2)/FAC2)
      end do

CCC   Integrate Gaussian distribution over specified range
      integ(1) = 0.0
      do i = 2,npts
      integ(i) = 0.5*(Gauss_dist(i-1) + Gauss_dist(i))*dm + integ(i-1)
      end do

CCC   Normalize integral to unity:
      do i = 1,npts
      integ(i) = integ(i)/integ(npts)
      end do

CCC  Extend integ() and xfunc() by one mesh point to avoid Lagrange
CCC  interpolation subroutine bug:
      integ(npts + 1) = integ(npts) + 0.01
      xfunc(npts + 1) = xfunc(npts)

      Return
      END

      Subroutine Particle_prop
      implicit none

      Common/Geant/geant_mass(60)
      real*4 geant_mass

CCC   Fills array geant_mass() with particle mass in GeV,
CCC   Row # = Geant Particle ID # code (see Geant Manual 3.10,
CCC   User's Guide, CONS 300-1 and 300-2).

      geant_mass(1)  = 0.0           ! Gamma
      geant_mass(2)  = 0.000511      ! Positron
      geant_mass(3)  = 0.000511      ! Electron
      geant_mass(4)  = 0.0           ! Neutrino
      geant_mass(5)  = 0.105659      ! Muon +
      geant_mass(6)  = 0.105659      ! Muon -
      geant_mass(7)  = 0.134693      ! Pion 0
      geant_mass(8)  = 0.139567      ! Pion +
      geant_mass(9)  = 0.139567      ! Pion -
      geant_mass(10) = 0.49767       ! Kaon 0 Long
      geant_mass(11) = 0.493667      ! Kaon +
      geant_mass(12) = 0.493667      ! Kaon -
      geant_mass(13) = 0.939573      ! Neutron
      geant_mass(14) = 0.93828       ! Proton
      geant_mass(15) = 0.93828       ! Antiproton
      geant_mass(16) = 0.49767       ! Kaon 0 Short
      geant_mass(17) = 0.5488        ! Eta
      geant_mass(18) = 1.11560       ! Lambda
      geant_mass(19) = 1.18936       ! Sigma +
      geant_mass(20) = 1.19246       ! Sigma 0
      geant_mass(21) = 1.19734       ! Sigma -
      geant_mass(22) = 1.31490       ! Xi 0
      geant_mass(23) = 1.32132       ! Xi -
      geant_mass(24) = 1.67245       ! Omega
      geant_mass(25) = 0.939573      ! Antineutron
      geant_mass(26) = 1.11560       ! Antilambda
      geant_mass(27) = 1.18936       ! Antisigma -
      geant_mass(28) = 1.19246       ! Antisigma 0
      geant_mass(29) = 1.19734       ! Antisigma +
      geant_mass(30) = 1.3149        ! Antixi 0
      geant_mass(31) = 1.32132       ! Antixi +
      geant_mass(32) = 1.67245       ! Antiomega +
      geant_mass(33) = 1.7842        ! Tau +
      geant_mass(34) = 1.7842        ! Tau -
      geant_mass(35) = 1.8694        ! D+
      geant_mass(36) = 1.8694        ! D-
      geant_mass(37) = 1.8647        ! D0
      geant_mass(38) = 1.8647        ! Anti D0
      geant_mass(39) = 1.9710        ! F+
      geant_mass(40) = 1.9710        ! F-
      geant_mass(41) = 2.2822        ! Lambda C+
      geant_mass(42) = 80.8000       ! W+
      geant_mass(43) = 80.8000       ! W-
      geant_mass(44) = 92.9000       ! Z0
      geant_mass(45) = 1.877         ! Deuteron
      geant_mass(46) = 2.817         ! Tritium
      geant_mass(47) = 3.755         ! Alpha
      geant_mass(48) = 0.0           ! Geantino

      Return
      END















