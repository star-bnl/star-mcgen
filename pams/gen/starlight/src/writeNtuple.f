c	this subroutine writes out the ntuple format currently preferred
c	by GSTAR, which fills the particle table, enables use of StMcEvent, etc.

      	subroutine writeNtuple(ievent,n,ipid,mdec,pxdec,pydec,pzdec,Edec)

	implicit NONE

	include 'inputp.inc'
	integer MM(2),DD(2),ievent,n,ipid,i
	real PP(3),VV(3),x,mdec,pxdec(4),pydec(4),pzdec(4),Edec(4)
    	real E,vt,ran
c----------------------------------------------------------------

c	initiate the event-- put header line in ntuple
	call hepevent ('starlight',1,2,99.,0.5,100.,0.1,197.,97.,197.,97.)

c----------------------------------------------------------------

c	initialize variables-- vertex at origin, no mothers, no daughters
 	vv(1) = 0.
      	vv(2) = 0.
      	vv(3) = 0.
      	vt = 0.
	MM(1) = 0
	MM(2) = 0
	DD(1) = 0
	DD(2) = 0

c----------------------------------------------------------------
c      loop over pairs
       do 100 i = 1,n,2

c	write the first particle of the pair

	pp(1) = pxdec(i)
	pp(2) = pydec(i)
	pp(3) = pzdec(i)
	E = Edec(i)
c	randomly assign one particle to be the anti-particle, and write
c	the particle info to the ntuple
c      	call ranmar(x,1)
   	x = ran(iseed)
      	if (x.lt.0.5) then
	  call heppart (i,1,ipid,MM,DD,pp,E,mdec,vv,vt)
      	else
	  call heppart (i,1,-ipid,MM,DD,pp,E,mdec,vv,vt)
	endif


c	write the second particle in the pair

	  pp(1) = pxdec(i+1)
	  pp(2) = pydec(i+1)
	  pp(3) = pzdec(i+1)
	  E = Edec(i+1)
      if (x.lt.0.5) then
c	  write the particle info to the ntuple
	  call heppart (i+1,1,-ipid,MM,DD,pp,E,mdec,vv,vt)
      	else
	  call heppart (i+1,1,ipid,MM,DD,pp,E,mdec,vv,vt)
	endif


 100   continue
c----------------------------------------------------------------

      	return
      	end
