c	this subroutine writes out the ntuple format currently preferred
c	by GSTAR, which fills the particle table, enables use of StMcEvent, etc.

        subroutine writejetsetNtuple(ievent)

     	implicit NONE

        include 'lujets.inc'
	integer MM(2),DD(2),i,ievent
	real PP(3),VV(3),vt

c	initiate the event-- put header line in ntuple
	call hepevent ('starlight',1,n,99.,0.5,100.,0.1,197.,97.,197.,97.)

c	initialize variables-- vertex at origin, no mothers, no daughters
 	vv(1) = 0.
      	vv(2) = 0.
      	vv(3) = 0.
      	vt = 0.
	MM(1) = 0
	MM(2) = 0
	DD(1) = 0
	DD(2) = 0

c	do the particle loop
      	do 100 i = 1,n

c	  define mother and daughter tracks, if existing, from the jetset
c	  event record; MM(1) contains the track number of the mother;
c	  DD(1) contains the track number of the first daughter; DD(2)
c	  contains the track number of the last daughter
	  MM(1) = k(i,3)
	  DD(1) = k(i,4)
	  DD(2) = k(i,5)

c	  define momentum and vertex info from jetset event record
	  pp(1) = p(i,1)
	  pp(2) = p(i,2)
	  pp(3) = p(i,2)

c      	  jetset works in mm and GEANT in cm
          vv(1) = v(i,1)/10.
          vv(2) = v(i,2)/10.
          vv(3) = v(i,3)/10.
c     	  jetset works in mm/c and GEANT in sec
          vt = v(i,4)*3.336*10**(-12)

c	  write the particle info to the ntuple
	  call heppart (i,k(i,1),k(i,2),MM,DD,pp,p(i,4),p(i,5),vv,vt)
 100   	continue

      	return
      	end
