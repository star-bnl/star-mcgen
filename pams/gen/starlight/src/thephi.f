c     This subroutine calculates angles for channels decayed by jetset.
      subroutine thephi(W,px,py,pz,E,theta,phi)

      implicit NONE

      include 'const.inc'
      real px,py,pz,E,theta,phi,W

      E = sqrt (W**2+px**2+py**2+pz**2)

      theta = acos(pz/sqrt(px**2+py**2+pz**2))
      phi = acos(px/sqrt(px**2+py**2))
      if ((px.eq.0).and.(py.eq.0)) phi = 0.
      if (py.lt.0) phi = 2*pi - phi

      return
      end
