C     This routine decays a particle into two particles of mass mdec,
c     taking spin into account

      subroutine twodecay (ipid,E,px0,py0,pz0,mdec,px1,py1,pz1,E1,
     &	px2,py2,pz2,E2,iFbadevent)

	implicit NONE

        include 'const.inc'
	include 'global.inc'
	include 'inputp.inc'
	include 'tables.inc'
	real mdec,px0,py0,pz0,px1,py1,pz1,px2,py2,pz2,ran,E
	real pmag, anglelep(0:100),ytest
	real phi,theta,xtest,dndtheta,thetalep,Ecm,E1,E2
	integer ipid,iFbadevent,i
	double precision betax,betay,betaz

c	set the mass of the daughter particles
	if((ip.eq.11).or.(ip.eq.13).or.(ip.eq.15)) mdec = mass
        if((ip.eq.113).or.(ip.eq.223).or.(ip.eq.33).or.(ip.eq.225)) 
     &		mdec = mpi
        if(ip.eq.333) mdec = mK
        if(ip.eq.335) then
c	decays 50% to K+/K-, 50% to K_0's
           ytest = RAN(ISEED)
           if(ytest.ge.0.5) then
                mdec = mK
           else
                mdec = 0.493677
           endif
	endif
        if(ip.eq.443) then
c	decays 50% to e+/e-, 50% to mu+/mu-
           ytest = RAN(ISEED)
           if(ytest.ge.0.5) then
                mdec = mel
           else
                mdec = mmu
           endif
        endif


c     calculate the magnitude of the momentum
	if(ip.eq.33) then
c	the rho pairs are produced at threshold
	      pmag = sqrt(mass*mass/4. - mdec*mdec)
	else
	      IF(E.lt.2*mdec) then
		WRITE(*,*) 'ERROR: E=',E
		iFbadevent = 1
	        return
	      endif
	      pmag = sqrt(E*E/4. - mdec*mdec)
	endif

c     pick an orientation, based on the spin
c	phi has a flat distribution in 2*pi
      phi = ran(ISEED)* 2.*pi

c     find theta, the angle between one of the outgoing particles and
c     the beamline, in the frame of the two photons

      if(spin.eq.0.) then
 100	theta = pi*ran(ISEED)
	xtest = ran(ISEED)
	dndtheta = sin(theta)	
	if(xtest.gt.dndtheta) goto 100

      elseif(spin.eq.0.5) then
C     calculate a table of integrated angle values for leptons
        anglelep(0) = 0
        do 125 i = 1,100
          theta = pi * float(i) /100.
          anglelep(i) = anglelep(i-1) + thetalep(E,theta)
 125    continue
        theta = 0.
	xtest = ran(ISEED)
        do 150 i = 1,100
          if(xtest.gt.(anglelep(i)/anglelep(100)))
     &         theta = pi * float(i) / 100.
 150     continue

      elseif(spin.eq.1.) then
 200	theta = pi*ran(ISEED)
	xtest = ran(ISEED)
	dndtheta = sin(theta)*cos(theta)*cos(theta)	
	if(xtest.gt.dndtheta) goto 200

      elseif(spin.eq.2.) then
 300	theta = pi*ran(ISEED)
	xtest = ran(ISEED)
	dndtheta = sin(theta)**5	
	if(xtest.gt.dndtheta) goto 300
      else
	write(*,*) 'This model cannot yet handle this spin 
     &		value: ',spin
      endif

c     compute unboosted momenta
      px1 = sin(theta)*cos(phi)*pmag
      py1 = sin(theta)*sin(phi)*pmag
      pz1 = cos(theta)*pmag
      px2 = -px1
      py2 = -py1
      pz2 = -pz1

c	compute energies
      if(ip.eq.33) then
      	Ecm = sqrt(mass**2+px0**2+py0**2+pz0**2)
      else
	Ecm = sqrt(E**2 + px0**2 + py0**2 + pz0**2)
      endif
      E1 = sqrt(mdec**2+px1**2+py1**2+pz1**2)
      E2 = sqrt(mdec**2+px2**2+py2**2+pz2**2)

c	decay tauons to electrons
c	note that after this routine px1, etc., refer to the electrons
      if(ip.eq.15) call taudecay(px1,py1,pz1,E1,px2,py2,pz2,E2)

c     lorentz transform into the lab frame
      betax = -(px0/Ecm)
      betay = -(py0/Ecm)
      betaz = -(pz0/Ecm)
      call transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent)
      call transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent)
	if(iFbadevent.eq.1) return

c       change particle id from that of parent to that of daughters
c       rhos and omegas, f2(1270) - each goes to pi+/pi-
          if((ip.eq.113).or.(ip.eq.223).or.(ip.eq.33).or.
     &		(ip.eq.225)) ipid=211
c	f2(1525)-- 50% to K+/K-, 50% to K0_S,K0_L
	  if(ip.eq.335) then 
		if(ytest.ge.0.5) then
		   ipid = 310
		else
		   ipid = 321
		endif
	  endif
c       phis-- this is decay to K+/K-
          if(ip.eq.333) ipid=321
c       J/psi -- 50% to e+/e-, 50% to mu+/mu-
          if(ip.eq.443) then
	     if(ytest.ge.0.5) then
                ipid = 11
             else
                ipid = 13
             endif
	  endif
c	taus decay into electrons
	  if(ip.eq.15) ipid = 11
c	electrons remain electrons, muons, remain muons
	  if ((ip.eq.11).or.(ip.eq.13)) ipid = ip

      return
      end

