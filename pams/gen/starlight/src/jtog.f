C     this function takes a jetset particle number and converts it
c      to GEANT

      integer function jtog(ip)

      implicit NONE

      integer ip

      jtog = 0
      if(ip.eq.22) jtog = 1
      if(ip.eq.-11) jtog = 2
      if(ip.eq.11) jtog = 3
      if((ip.eq.12).or.(ip.eq.-12).or.(ip.eq.14).or.(ip.eq.-14).or.
     &        (ip.eq.16).or.(ip.eq.-16)) jtog = 4
      if(ip.eq.13) jtog = 5
      if(ip.eq.-13) jtog = 6
      if(ip.eq.111) jtog = 7
      if(ip.eq.211) jtog = 8
      if(ip.eq.-211) jtog = 9
      if(ip.eq.130) jtog = 10
      if(ip.eq.321) jtog = 11
      if(ip.eq.-321) jtog = 12
      if(ip.eq.2112) jtog = 13
      if(ip.eq.2212) jtog = 14
      if(ip.eq.-2212) jtog = 15
      if(ip.eq.310) jtog = 16
      if(ip.eq.221) jtog = 17
      if(ip.eq.3122) jtog = 18
      if(ip.eq.3222) jtog = 19
      if(ip.eq.3212) jtog = 20
      if(ip.eq.3112) jtog = 21
      if(ip.eq.3322) jtog = 22
      if(ip.eq.3312) jtog = 23
      if(ip.eq.3334) jtog = 24
      if(ip.eq.-2112) jtog = 25
      if(ip.eq.-3122) jtog = 26
      if(ip.eq.-3222) jtog = 27
      if(ip.eq.-3212) jtog = 28
      if(ip.eq.-3112) jtog = 29
      if(ip.eq.-3322) jtog = 30
      if(ip.eq.-3312) jtog = 31
      if(ip.eq.-3334) jtog = 32

      return
      end

