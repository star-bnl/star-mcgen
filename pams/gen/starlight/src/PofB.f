      real function PofB(D)
C calculate photon and hadron breakup probability as function of d - distance b/w niclei
C step in d is .1 fm
C all d must be later converted to units of Mev^-1, which is d(fm)/(hbarc*1000)

      IMPLICIT NONE
      include 'Breakup.inc'
      include 'D2LParam.inc'
      include 'inputp.inc'
      real*8 Step,BIter,Zz,Az
      real*8 P1n,PXn,PHadr,PPhoton,kmax,gammatarg,D
      integer i,k,NStep,ifirst
      real Dlow,DeltaP,DeltaD

      save ifirst,Nstep
      data ifirst /0/

      if (ifirst.eq.0) then  

C  Nstep is the number of entries in the table
    
       NStep = 1000
       kmax= 1.D7
       if (gamma_em.gt.500) kmax = 1.E10
       Az = A
       Zz = Z
       Bmin = 1.75*RNuc  ! for hard sphere ions BMin = 2*RNuc
       Step = 1.01  ! we will multiplicatively increase BIter by 1% 


        if (ibreakup .eq. 1)   write(6,31)2.*RNuc
 31     format(' Hard Sphere Breakup criteria.  b > ',F7.3)

        if (ibreakup .eq. 2) write(6,32)bmin
 32     format('Requiring XnXn [Coulomb] breakup. bmin=',F7.4)

        if (ibreakup .eq. 3) write(6,33)bmin
 33     format(' Requiring 1n1n [Coulomb only] breakup. bmin =',F7.4)

        if (ibreakup .eq. 4) write(6,34)bmin
 34     format(' Requiring both nuclei remain intact. bmin = ',F7.4)

        if (ibreakup .eq. 5) write(6,35)bmin
 35     format(' Requiring no hadronic interactions bmin = ',F7.4)
     
       BIter = 1.0*Bmin
       do 70 k = 1,NStep
       	PHadr = 0.
       	PPhoton = 0.
      	call hadronbreakup(PHadr,BIter,Zz,Az,gamma_em)
        gammatarg=2.*gamma_em*gamma_em-1
      	call photonbreakup(P1n,Pxn,BIter,Zz,Az,gammatarg,kmax)

        if (ibreakup.eq.1) then
           PPhoton = 1
           if (BIter.gt.2*RNuc) then
              Phadr = 0.0
           else 
              Phadr = 999.
           endif
        endif 
  
         if (ibreakup.eq.2) PPhoton = (1-exp(-1*Pxn))**2

         if (ibreakup.eq.3) PPhoton=(p1n*exp(-1*pxn))**2

         if (ibreakup.eq.4) PPhoton=exp(-2*pxn)

         if (ibreakup.eq.5) PPhoton=1

C      	write(77,78) BIter,Phadr,P1n,Pxn,PPhoton  ! useful for debugging
        BIter = BIter*Step
         ProbTot(k)=exp(-1*Phadr)*PPhoton
   70   continue
      ifirst = 1
      endif

C  use the lookup table and return

        PofB = 1.
        if (D.gt.0.0) then
C  Now we must determine which step number in d corresponds 
C  to this D, and use appropriate Ptot(D_i)
        i = log(D/Bmin)/log(1.01)
        if (i.le.0) then
           PofB = ProbTot(1) 
        else
           if(i.ge.Nstep) then
              PofB = ProbTot(Nstep)
           else                      ! interpolate
              DLow = BMin*(1.01)**i
              DeltaD = .01*DLow
              DeltaP = ProbTot(i+1) - ProbTot(i)
              PofB = ProbTot(i) + DeltaP*(D-DLow)/DeltaD
           endif
         endif
       endif
       return
       write(6,74)D
 74    format(' ERROR in function PofB.  D= ',E11.3)
       STOP

 78   FORMAT(F10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5)     
      end


