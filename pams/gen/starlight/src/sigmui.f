c     This function gives the two muon cross section as a
c     function of Y and W.  Using the formula given in G. Soff
c     et. al Nuclear Equation of State, part B, 579

      real function sigmui(w)

      include 'const.inc'
      include 'global.inc'
      real s,Etest
      real w,deltat,xL

      s = w**2
      Etest = 4.*mass**2/s
      deltat = s * sqrt(1.-Etest)
      xL = 2.*log(sqrt(s)/(2.*mass) + sqrt(1/Etest - 1.))

      sigmui = 4.*pi*alpha**2/s * hbarc**2 *
     & ((1 + Etest -0.5*Etest**2) * xL -
     &            (1/s + Etest/s) * deltat)

c	write(*,*) 'w,deltat,xL,sigmui',w,deltat,xL,sigmui
      if (Etest.gt.1.) sigmui = 0.

      continue

      end

