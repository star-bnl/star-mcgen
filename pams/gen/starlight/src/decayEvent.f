c     this subroutine decays particles and writes events to a file

      subroutine decayEvent(ievent,W,E,px,py,pz)

      implicit NONE

      include 'sig.inc'
      include 'lujets.inc'
      include 'inputp.inc'
      integer ievent
      real w,px,py,pz,W,theta,phi,px0,py0,pz0,mdec
      real pxdec(4),pydec(4),pzdec(4),Edec(4)
      integer ipid,iFbadevent,i
      double precision E

c     zero out the decay momenta
      do 100 i = 1,4
	  pxdec(i) = 0
	  pydec(i) = 0
	  pzdec(i) = 0
	  Edec(i) = 0
 100  continue

c     deal with single particle resonances  of the eta, eta',
c     and f0, which can be handled by jetset
      if ( (ip.eq.221).or.(ip.eq.331).or.
     & (ip.eq.10221).or.(ip.eq.441).or.(ip.eq.111).or.(ip.eq.115)) then
c       calculate theta and phi of the particle
        call thephi(W,px,py,pz,E,theta,phi)
c       let jetset decay the particle
        call lu1ent(0,ip,E,theta,phi)
C       call lulist(1)
c       write the decay to a file to be read by GEANT or gstar
        if (iout.eq.1) call writejetsetText(ievent)
        if (iout.eq.2) call writejetsetGSTARtext(ievent)
        if (iout.eq.3) call writejetsetNtuple(ievent)
c       increment the decay counter
        ievent = ievent + 1

c     deal with rho pairs at threshold
       elseif (ip.eq.33) then
c      decay each of the rho0's  to pi+pi-
         px0 = px /2.
         py0 = py /2.
         pz0 = pz /2.
         call twodecay (ipid,E,px0,py0,pz0,mdec,pxdec(1),pydec(1),
     &	pzdec(1),Edec(1),pxdec(2),pydec(2),pzdec(2),Edec(2),iFbadevent)
         call twodecay (ipid,E,px0,py0,pz0,mdec,pxdec(3),pydec(3),
     &	pzdec(3),Edec(3),pxdec(4),pydec(4),pzdec(4),Edec(4),iFbadevent)

        if (iout.eq.1)
     &     call writeText(ievent,4,ipid,pxdec,pydec,pzdec)
        if (iout.eq.2)
     &     call writeGSTARtext(ievent,4,ipid,pxdec,pydec,pzdec)
        if (iout.eq.3)
     &     call writeNtuple(ievent,4,ipid,pxdec,pydec,pzdec,Edec)
c       increment the decay counter
        ievent = ievent + 1

c     deal with electrons, muons tauons and f0(1525) to pi+ pi-
        elseif((ip.eq.11).or.(ip.eq.13).or.(ip.eq.15).or.(ip.eq.225)
     &  .or.(ip.eq.335)) then
          iFbadevent = 0
         call twodecay (ipid,E,px,py,pz,mdec,pxdec(1),pydec(1),
     &	pzdec(1),Edec(1),pxdec(2),pydec(2),pzdec(2),Edec(2),iFbadevent)
          if (iFbadevent.eq.0) then
             if (iout.eq.1)
     &         call writeText(ievent,2,ipid,pxdec,pydec,pzdec)
             if (iout.eq.2)
     &         call writeGSTARtext(ievent,2,ipid,pxdec,pydec,pzdec)
	     if (iout.eq.3)
     & 		call writeNtuple(ievent,2,ipid,mdec,pxdec,pydec,
     &		pzdec,Edec)
  
c       increment the decay counter
             ievent = ievent + 1
           endif

c       deal with vector meson decays
        elseif((ip.eq.113).or.(ip.eq.223).or.(ip.eq.333).or.
     &          (ip.eq.443)) then
         call twodecay (ipid,E,px,py,pz,mdec,pxdec(1),pydec(1),
     &	pzdec(1),Edec(1),pxdec(2),pydec(2),pzdec(2),Edec(2),iFbadevent)
             if (iout.eq.1)
     &         call writeText(ievent,2,ipid,pxdec,pydec,pzdec)
             if (iout.eq.2)
     &         call writeGSTARtext(ievent,2,ipid,pxdec,pydec,pzdec)
	     if (iout.eq.3)
     &		call writeNtuple(ievent,2,ipid,mdec,pxdec,pydec,
     &		pzdec,Edec)

c       increment the decay counter
             ievent = ievent + 1

        else
c       deal with invalid particle ID's
        write (25,*)
     &  'THIS CODE IS NOT YET ABLE TO HANDLE THAT DECAY CHANNEL'
c       increment the decay counter
             ievent = ievent + 1
        endif

      return
      end
