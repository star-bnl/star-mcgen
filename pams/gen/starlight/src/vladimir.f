C Calculate the Luminosity-function d^2L/dMdY where M is the invariant
C mass and Y the rapidity C of the photon-photon-system
C
C Calculations are performed as brute-force implementation of 
C Eq. 3 in in STAR Note 380, except that a realistic nucleus-nucleus
C interactin probability is substituted for the theta-function in that equation

       REAL FUNCTION  D2LDMDY(M,Y)
       IMPLICIT NONE
       REAL M,Y

C M : invariant mass (MeV)
C Y : rapidity
C D2LDMDY : diff. Lum. (MeV**-1)

       include 'D2LParam.inc'
       include 'const.inc'
       Double Precision w1,w2,gamma
       double precision Integral
       External Integral
       common/PhotonEnergies/W1,W2,gamma

C calculate photon frequencies from M,Y
     
       W1 = M/2.0 * EXP (Y)
       W2 = M/2.0 * EXP (-Y)
       gamma = gamma_em

C calculate Rho in MeV^-1 instead of fm
C use A or R, if A<= 0
  
       D2LDMDY = 2.0 / M * Z**4 * alpha**2 * Integral()
       Normalize = D2LDMDY*M/(2.0 * Z**4 * alpha**2)
       RETURN
       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       DOUBLE PRECISION FUNCTION Integral()
       implicit none
       external Integrand
       external RADMUL
       DOUBLE PRECISION RM
       include 'const.inc'
       include 'D2LParam.inc'
       double precision w1,w2,gamma
       common/PhotonEnergies/W1,W2,gamma
       Integer NFNEVL,Summary ! NFNEVL - number of iterative evaluations done during the integration, Summary - summary report of integration
       REAL EPS,WK(2799000)! eps for integration, wk - working space
       REAL Upper(3),Lower(3),Result, ResErr !limits of integration, result of integration and integration error
       REAL u1,u2,b1,b2
       integer totsummary,NIter,NEval,NIterMin

C Niter is maximal num evaluations of integrand, NEval - actual number, NIterMin - minimal number of allowed iterations
C totsummary is useful for debuggin this code. This variable is in format 4xxxx or 20x0x or 2x0x0 or 1000x
C 4,2 or 1 stand for number of sub-integrals that a big integral is split into
C x is "summary" variable from each of subintegrals. If any of x is not 1 or 0 you're in trouble, and have to investigate further

       EPS = .01*Normalize  ! this is EPS for integration, 1% of previous integral value 
       RM = RNuc / hbarcmev
   
CVM  CERNLIB subroutine RADMUL has the following syntax:

CVM  CALL RADMUL(Integrand,NumDimensions,LowerIntLimits,UpperIntLimits,MinNumIterations,MaxNumIterations,EPS,WorkingSpaceVector,LengthWorkingSpaceVector,Result,ResultError,NumIterationsDone,Summary)
CVM  the last four parameters are outputs of the subroutine

CVM    Set up limits of integration:

        NIter = 10000 + 1000000*Normalize ! if integral value is very small, we don't do too many iterations to get precision down to 1%
        NIterMin = 600
        u1 = 9.*gamma/w1 ! upper boundary in B1
        u2 = 9.*gamma/w2 ! upper boundary in B2
        B1 = .4*gamma/w1 ! intermediate boundary in B1
        B2 = .4*gamma/w2 ! intermediate boundary in B1

C the trick is that u1,2 and b1,2 could be less than RM - the lower integration boundary, thus integration area splits into 4, 2 or 1 pieces
       if (u1.lt.RM) then
          integral = 0
          totsummary = 0
          NEval = 0
       elseif (B1.gt.RM) then
          if (u2.lt.RM) then
             integral = 0
             totsummary = 0
             NEval = 0
          elseif (B2.gt.RM) then
CVM             integral has 4 parts
             Integral = 0 
	     totsummary = 40000
             NEval = 0
             Lower(1)=RM
             Lower(2)=RM
             Lower(3)=0.
             Upper(3)=2.*pi
             Upper(1)=B1
             Upper(2)=B2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL

             Integral = Integral + Result
 	     totsummary = totsummary + 1000*summary
             NEval = NEval + NFNEVL
             Upper(1)= u1
             Upper(2)= B2
             lower(1)= b1
             lower(2)= RM

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL

             Integral = Integral + Result
             totsummary = totsummary + 100*summary
             NEval = NEval + NFNEVL
             Upper(1)= B1
             Upper(2)= u2
             lower(1)= RM
             lower(2)= B2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL

             Integral = Integral + Result
             totsummary = totsummary + 10*summary
             NEval = NEval + NFNEVL             
             Upper(1)= u1
             Upper(2)= u2
             lower(1)= b1
             lower(2)= b2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL

             Integral = Integral + Result
             totsummary = totsummary + summary
             NEval = NEval + NFNEVL
          else 
CVM             integral has 2 parts, b2 integral has only 1 component
             Integral = 0
	     totsummary = 20000
             NEval = 0
             Lower(1)=RM
             Lower(2)=RM
             Lower(3)=0.
             Upper(3)=2.*pi
             Upper(1)=B1
             Upper(2)=u2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL
 
             Integral = Integral + Result
             totsummary = totsummary + 100*summary
             NEval = NEval + NFNEVL
             Upper(1)= u1
             lower(1)= b1

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL

             Integral = Integral + Result         
             totsummary = totsummary + summary
             NEval = NEval + NFNEVL            
         
          endif
       else
          if (u2.lt.RM) then
             integral = 0
             totsummary = 0
             NEval = 0
          elseif (B2.gt.RM) then
CVM             integral has 2 parts, b1 integral has only 1 component
             Integral = 0
	     totsummary = 20000
             NEval = 0
             Lower(1)=RM
             Lower(2)=RM
             Lower(3)=0.
             Upper(3)=2.*pi
             Upper(1)=u1
             Upper(2)=b2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL
 
             Integral = Integral + Result
             totsummary = totsummary + 1000*summary
             NEval = NEval + NFNEVL 
             Upper(2)= u2
             lower(2)= b2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL
 
             Integral = Integral + Result         
             totsummary = totsummary + 10*summary
             NEval = NEval + NFNEVL            
          else 
C             integral has 1 parts
             Integral = 0
	     totsummary = 10000 
             NEval = 0
             Lower(1)=RM
             Lower(2)=RM
             Lower(3)=0.
             Upper(3)=2.*pi
             Upper(1)=u1
             Upper(2)=u2

             call radmul(Integrand,3,Lower,Upper,NIterMin,NIter,
     *       Eps,WK,NIter,Result,ResErr,NFNEVL,Summary)
C            print*,'Int = ',Result,' Err = ',ResErr,' NFNEVL= ',NFNEVL

             Integral = Integral + Result         
             totsummary = totsummary + summary             
             NEval = NEval + NFNEVL
          endif
       endif

CVM    multiply your integral by 2*pi to account for one ommitted phi integral
       Integral = 2*pi*Integral

CVM  write integral, M,Y,summary, numiterations,Eps  - this line is for debugging only
C       write(77,88)Integral,w1,w2,2*sqrt(w1*w2),.5*log(w1/w2),
C    * totsummary,NEval,Eps
 88    format(F10.8,3x,F10.3,3x,F10.3,3x,F10.5,3x,F10.8,3x,I5,3x,
     * I8,3X,F8.6)

       return
       end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       REAL FUNCTION Nphoton (W,gamma,Rho)

       DOUBLE PRECISION gamma,w, pi
       parameter (pi = 3.141592654 )
       REAL Rho, WGamma, WPhib
       REAL WGR
C       print*,W,Rho,gamma
          WGamma= W/gamma
          WGR= 1.0*Wgamma*Rho

CVM  factor of (Z^2*alpha) is omitted

       Wphib = WGamma*besk1(WGR)
       Nphoton = 1.0/PI**2 * Wphib**2
C       print*,'Nphoton= ',Nphoton
       RETURN
       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
       REAL FUNCTION integrand(N,X)
       Implicit none

C  b1,b2 are in units of MeV^-1
       Integer N,i ! N is dimension of our integral
       REAL X
       Dimension X(*)
       real b1,b2,theta
       real nphoton,PofB
       double precision W1,W2,gamma,hbarcmev,D
       parameter( hbarcmev = 197.327053)
       common/PhotonEnergies/W1,W2,gamma
  
       b1=x(1)
       b2=x(2)
       theta=x(3)
       
C  breakup expects distances in fermis, 
C  so convert to fermis (factor of hbarcmev)
       D = sqrt(b1**2 + b2**2 - 2*b1*b2*cos(theta))*hbarcmev

       integrand = Nphoton(w1,gamma,b1)*Nphoton(w2,gamma,b2)
     * *b1*b2*PofB(D)

C       write(*,75)D,w1,w2,i,PofB,integrand        ! for debugging purposes
 75    FORMAT(F7.2,1x,F9.2,1x,F10.2,1x,I4,1x,F10.8,1x,F12.4)
       return 
       end
 
