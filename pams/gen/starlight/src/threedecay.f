c     This routine decays a particle of mass m0 into three
c     particles of mass m1,m2,m3.  Based ONLY on phase space--
c     no matrix element
c     resultinig momenta are given in the lab frame
      subroutine threedecay (m0,m1,m2,m3,px0,py0,pz0,px1,py1,pz1,
     &   px2,py2,pz2,px3,py3,pz3)

        include 'const.inc'
      real m0,m1,m2,m3,px1,py1,pz1,px2,py2,pz2,px3,py3,pz3
      real emax1,emax2,p1,p2,theta12,phi,theta,lambda,ran1,ran2,ran
      real px01,py01,px02,py02,px03,py03

c     make sure that the decay is kinematically alowed
      if ((m0-m1-m2-m3).le.0) then
            print *,'kinematically unallowed decay'
            stop
      endif

c     figure out emax
c     kinematic emax

      emax1 = sqrt(m1**2 + (4./9.)*(m0**2-m1**2-m2*2-m3**2))
      emax2 = sqrt(m2**2 + (4./9.)*(m0**2-m1**2-m2*2-m3**2))
c     also just mass differences, take smaller
      em1 = m0 - m2 - m3
      em2 = m0 - m2 - m3
      if (em1.lt.emax1) emax1 = em1
      if (em2.lt.emax2) emax2 = em2

c     pick an e1 and e2 that fall within the allowed kinematic
c     region.

 10   call ranmar(ran1,1)
      call ranmar(ran2,2)
      call ranmar(ran3,3)

      e1 = m1 + (emax1 - m1) * ran1
      e2 = m2 + (emax2 - m2) * ran2
      theta12 = pi * ran3

      p1 = sqrt(e1**2 - m1**2)
      p2 = sqrt(e2**2 - m2**2)
      p3 = sqrt((m0 - e1 - e2)**2  - m3**2)
      e3min = sqrt((p1-p2)**2 + m3**2)

      if ( (e1+e2+e3min-m0).gt.0) go to 10
      if ( (p3-p1-p2).gt.0) go to 10

c     find the angle between p1 and p2
      theta12 = acos ((p3**2 - p1**2 - p2**2) / (2.* p1 * p2))

      px01 = p1
      py01 = 0
      px02 = p2 * cos(theta12)
      py02 = p2 * sin(theta12)
      px03 = - p1 - p2 * cos(theta12)
      py03 = - p2 * sin(theta12)

c     phi and theta pick the direction of p1, lambda is
c     the oreintation angle about that axis.
      call ranmar(ran,1)

      phi = ran * 2. * pi
      call ranmar(ran,1)
      theta = ran * pi
      call ranmar(ran,1)
      lambda = ran * 2. * pi


c     transform to unboosted lab coordinates
      px1 = px01*sin(theta)*cos(phi) +
     &      py01*cos(lambda)*cos(theta)*cos(phi) -
     &      py01*sin(lambda)*sin(phi)
      py1 = px01*sin(theta)*sin(phi) +
     &      py01*cos(lambda)*cos(theta)*sin(phi) +
     &      py01*sin(lambda)*cos(phi)
      pz1 = px01*cos(theta) -
     &      py01*cos(lambda)*sin(theta)
      px2 = px02*sin(theta)*cos(phi) +
     &      py02*cos(lambda)*cos(theta)*cos(phi) -
     &      py02*sin(lambda)*sin(phi)
      py2 = px02*sin(theta)*sin(phi) +
     &      py02*cos(lambda)*cos(theta)*sin(phi) +
     &      py02*sin(lambda)*cos(phi)
      pz2 = px02*cos(theta) -
     &      py02*cos(lambda)*sin(theta)
      px3 = px03*sin(theta)*cos(phi) +
     &      py03*cos(lambda)*cos(theta)*cos(phi) -
     &      py03*sin(lambda)*sin(phi)
      py3 = px03*sin(theta)*sin(phi) +
     &      py03*cos(lambda)*cos(theta)*sin(phi) +
     &      py03*sin(lambda)*cos(phi)
      pz3 = px03*cos(theta) -
     &      py03*cos(lambda)*sin(theta)


c       lorentz transform into the lab frame
         Ecm = sqrt(m0**2+px0**2+py0**2+pz0**2)
         E1 = sqrt(m1**2+px1**2+py1**2+pz1**2)
         E2 = sqrt(m2**2+px2**2+py2**2+pz2**2)
         E3 = sqrt(m3**2+px3**2+py3**2+pz3**2)

         betax = -(px0/Ecm)
         betay = -(py0/Ecm)
         betaz = -(pz0/Ecm)
         call transform (betax,betay,betaz,E1,px1,py1,pz1)
         call transform (betax,betay,betaz,E2,px2,py2,pz2)
         call transform (betax,betay,betaz,E3,px3,py3,pz3)


      return
      end
