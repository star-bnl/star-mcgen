c     this routine assumes that the tauons decay to electrons and
c     calculates the directons of the decays
      subroutine taudecay(px1,py1,pz1,E1,px2,py2,pz2,E2)

      implicit NONE

      include 'const.inc'
      include 'tables.inc'
      include 'global.inc'
      include 'inputp.inc'
      real ran1,ran2,px1,py1,pz1,E1,px2,py2,pz2,E2
      real Ee1,Ee2,theta1,theta2,phi1,phi2
      real pmag1,pex1,pey1,pez1,pex2,pey2,pez2,pmag2
      real betax,betay,betaz,dir,ran
      integer i

c     the highest energy attainable by the electron in this system is
c     .5 * mass of the tau

c     get two random numbers to compare with
c      call ranmar(ran1,1)
c      call ranmar(ran2,2)
      ran1 = ran(iseed) * dgammade(100)
      ran2 = ran(iseed) * dgammade(100)

c     compute the energies that correspond to those numbers
      Ee1 = 0.
      Ee2 = 0.
      do 200 i = 1,100
        if (ran1.gt.dgammade(i)) Ee1 = float(i) /100. * .5 * mass
        if (ran2.gt.dgammade(i)) Ee2 = float(i) /100. * .5 * mass
 200  continue

c     to find the direction of emmission, first
c     we determine if the tauons have spin of +1 or -1 along the
c     direction of the beam line
c      call ranmar(ran1,1)
      dir = 1.
      if (ran(iseed).lt.0.5) dir = -1.

c     get two random numbers to compare with
c      call ranmar(ran1,1)
c      call ranmar(ran2,2)
      ran1 = ran(iseed) * tautolangle(100)
      ran2 = ran(iseed) * tautolangle(100)

c     find the angles corrsponding to those numbers
      theta1 = 0.
      theta2 = 0.
      do 400 i = 1,100
        if (ran1.gt.tautolangle(i)) theta1 = pi * float(i) /100.
        if (ran2.gt.tautolangle(i)) theta2 = pi * float(i) /100.
 400  continue

c     grab another two random numbers to determine phi's
c      call ranmar(ran1,1)
c      call ranmar(ran2,2)
      phi1 = ran(iseed) * 2. * pi
      phi2 = ran(iseed) * 2. * pi
c     figure out the momenta of the electron in the frames of the
c     tauons from which they decayed, that is electron1 is in the
c     rest frame of tauon1 and e2 is in the rest fram of tau2
c     again the electrons are assumed to be massless
      pmag1 = Ee1
      pex1 = cos(phi1)*sin(theta1)*pmag1
      pey1 = sin(phi1)*sin(theta1)*pmag1
      pez1 = cos(theta1)*pmag1*dir
      pmag2 = Ee2
      pex2 = cos(phi2)*sin(theta2)*pmag2
      pey2 = sin(phi2)*sin(theta2)*pmag2
      pez2 = cos(theta2)*pmag2*(-1.*dir)
c     now Lorentz transform into the frame of each of the particles
c     do particle one first
      betax = -(px1/E1)
      betay = -(py1/E1)
      betaz = -(pz1/E1)
      call transform (betax,betay,betaz,Ee1,pex1,pey1,pez1)
c     then do particle two
      betax = -(px1/E2)
      betay = -(py1/E2)
      betaz = -(pz1/E2)
      call transform (betax,betay,betaz,Ee2,pex2,pey2,pez2)
c     finally dump the electron values into the approriate
c     variables
      E1 = Ee1
      E2 = Ee2
      px1 = pex1
      px2 = pex2
      py1 = pey1
      py2 = pey2
      pz1 = pez1
      pz2 = pez2

      return
      end

