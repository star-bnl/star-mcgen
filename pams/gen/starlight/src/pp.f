      double precision function pp(E)
C     returns one random draw from pp(E) distribution
      implicit none
      real E,satisfy,x,ereds
      double precision u,test,formf,coef,Cm

      include 'const.inc'
      include 'D2LParam.inc' 
      include 'inputp.inc'

      satisfy = 0.

C
      ereds=(E/gamma_em)**2

C sqrt(3)*E/gamma_em is p_t where the distribution is a maximum

      Cm=sqrt(3.)*E/gamma_em

C the amplitude of the p_t spectrum at the maximum

      Coef = 3.0*(formf(Cm**2+ereds)**2)*Cm**3/
     * (2*pi*(ereds+ Cm**2))**2


C pick a test value pp, and find the amplitude there

      x=ran(iseed)
      pp = x*5.*hbarc/RNuc

      test = (formf(pp**2+ereds)**2)*pp**3/
     * (2*pi*(ereds + pp**2))**2


      do while (satisfy.eq.0.)
        u = ran(iseed)
C       write(76,*)u,Coef,test,Cmax,E
        if (u*Coef.le.test) then 
          satisfy = 1.
        else 
          x=ran(iseed)
          pp = 5*hbarc/RNuc*x
          test = (formf(pp**2+ereds)**2)*pp**3/
     *    (2*pi*(ereds + pp**2))**2
        endif
      enddo
      return
      end
