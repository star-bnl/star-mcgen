c	This subroutine chooses the appropriate cross section code
c	to call, depending on which channel is being studied.

	subroutine sigmacalc

	implicit NONE

        include 'inputp.inc'
        include 'Ftable.inc'
        include 'sig.inc'


C       take care of one particle states
        if ( (ip.eq.221).or.(ip.eq.331).or.(ip.eq.441).or.
     &       (ip.eq.10221).or.(ip.eq.225).or.(ip.eq.115).or.
     &    (ip.eq.335).or.(ip.eq.33) ) then
           write(*,*) 'Calling sigmadelta...'
           call sigmadelta

c     	take care of two particle decays
          elseif((ip.eq.11).or.(ip.eq.13).or.(ip.eq.15)) then
            write(*,*) 'Calling sigma2...'
            call sigma2

c    	take care of vector mesons
          elseif((ip.eq.113).or.(ip.eq.223).or.(ip.eq.333).or.
     &    (ip.eq.443)) then
            if(gg_or_gP.eq.2) then
              write(*,*) 'Calling sigmavm...'
              call sigmavm
            elseif(gg_or_gP.eq.3) then
              write(*,*) 'Calling sigmavmw...'
              call sigmavmw
            endif
          else
            print *,'invalid particle'
            stop
          endif

      return
      end
