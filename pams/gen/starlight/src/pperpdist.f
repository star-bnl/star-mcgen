c     this is used so much that we can't really do it with a table
C     unfortunately the expression for integrated pperp involves a nasty
c     taylor expansion.  This expression is good to around u1 u0 lt 10
c     which is twice the nuclear radius

      real function pperpdist (u0,u1)

      implicit NONE

      real u0,u1
      real u2,lastterm

      u2 = u0+u1
      lastterm = log (1. + u1/u0) - u1 +
     &  (u2**2 - u0**2)/    4. -       (u2**3 - u0**3) /  18. +
     &  (u2**4 - u0**4) /   96. -      (u2**5 - u0**5) /  600. +
     &  (u2**6 - u0**6) /   4320. -    (u2**7 - u0**7) /  35480. +
     &  (u2**8 - u0**8) /   322460. -  (u2**9 - u0**9) /  3265920. +
     &  (u2**10 - u0**10) / 36288000. - (u2**11-u0**11)/439084800. +
     &  (u2**12 - u0**12)/5748019000. - (u2**13-u0**13)/80951270000. +
     &  (u2**14 - u0**14)/1.220496E12  - (u2**15-u0**15)/1.961511E13  +
     &  (u2**16 - u0**16) / (16.*2.092279E13) -
     &  (u2**17 - u0**17) / (17.*3.556874E14) +
     &  (u2**18 - u0**18) / (18.*6.402374E15) -
     &  (u2**19 - u0**19) / 2.311257E18  +
     &  (u2**20 - u0**20) / 4.865804E19
      pperpdist = -exp(-u0) + u0*exp(-u1-u0)/(u1+u0) + (u0+1.)*lastterm

      return
      end
