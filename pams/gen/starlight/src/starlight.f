      Program starlight
c     Version 2.0
c     Evan Scannapieco, Joakim Nystrand and Janet Seger
c     April 28, 2000
c     This program is a monte carlo for two photon and
c     photon-Pomeron interactions in peripheral heavy ion collisions.
c	___________
c	How to use:
c	  1)  copy contents of src directory to your own src directory
c	  2)  copy jet.dat and starlight.in from the bin directory into
c	your own bin directory
c	  3)  edit the Makefile to use the appropriate lines depending on
c	whether you are compiling on Linux or Solaris; then type "make" to
c	compile-- the executable "starlight" will be placed in the bin
c	directory (see starlight.doc in the doc directory for an
c	explanation of which lines in the Makefile to edit)
c	  4)  edit starlight.in to reflect your preferred choice of input
c	parameters-- see starlight.doc in the doc directory for an
c	explanation of input parameters
c	  5)  type "starlight" to run; if you selected either text output
c	format, the output will be written to the file starlight.out in
c	the bin directory; if you selected ntuple output, the output will
c	be written to the file evgen.1.nt.gz, and will need to be unzipped
c	(type "gunzip evgen.1.nt.gz")
c
c	___________
c	Subroutines called:
c		luupda-- reads jetset branching ratios from jet.dat
c		input-- reads input parameters from starlight.in
c		setConst-- sets constants like pi, mass, width, etc.
c		newparam-- determines whether or not differential
c			luminosity tables need to be re-calculated
c		diffLum_2gamma OR diffLum_vm-- calculates differntial
c			luminosity tables, if necessary, and writes to
c			starlight.dat file; separate routines for 2-photon
c			and vector meson channels
c		readDiffLum-- read in luminosity tables from starlight.dat
c		sigmacalc-- calculates cross section
c		  calls:  sigmadelta OR sigma2 OR sigmavm OR sigmavmw, 
c		  depending on the type of channel; sigmadelta is called
c		  for one-particle 2-photon channels, sigma2 is called
c		  for 2-particle 2-photon channels, sigmavm is called for
c		  vector meson channels using narrow resonance, and
c		  sigmavmw is called for vector meson channels using 
c		  Breit-Wigner resonance
c		tablecalc-- calculates tables used in tau decay 
c		pickw and picky OR pickwy_vm-- picks values for invariant
c			mass and rapidity of the particle, based on 
c			differential luminosity tables; there are
c			currently two separate routines for the 2-photon
c			channels, and one for the vector meson channels
c		momenta OR vmomenta-- given the w and y chosen, calculates
c			the momentum components and energy of the
c			particle; there are currently separate routines 
c			for the 2-photon channels and the vector meson
c			channels
c		decayEvent-- decays the particle, then writes results in
c			one of three formats-- text, GSTAR text, or ntuple
c		  calls: thephi and lu1ent if channel is to be decayed
c		  using jetset; twodecay if channel decays into a
c		  particle/antiparticle pair; taudecay for the unique case
c		  of tau pairs; some write-out routine-- if the channel is
c		  decayed through jetset, will call writejetsetText Or
c		  writejetsetGSTARtext OR writejetsetNtuple-- otherwise,
c		  will call writeText OR writeGSTARtext or writeNtuple 
c	___________
c
c     All units are in fm, sec, and GeV
c  ************************************************************

	implicit NONE

	include 'ludat1.inc'
	include 'const.inc'
	include 'D2LParam.inc'
	include 'inputp.inc'
	include 'global.inc'

	integer tcheck,i
	real w,y,E,px,py,pz
       	logical new

c  ************************************************************

c     Give user some indication that the program is working
      print *,'                     ',
     &  '************** Starlight ***************'
      print *,'                     ',
     &  '*******The two-photon Monte Carlo*******'
      print *,'                     ',
     &  'Simulates two photon collisions at RHIC'
      print *,'                                     and uses:'

c  ************************************************************
c     setting up jetset

c     setup jetset such that it does not decay particles with
c     lifetimes longer than 1 mm/c
      mstj(22) = 2
      parj(71) = 1.

c     read in modified branching ratios for jetset
      open (unit=30,file='jet.dat',status='unknown')
      call luupda(2,30)
      close(unit=30)

c  ************************************************************
c     read input parameters

      call input

c  ************************************************************
c     define constants

      call setConst

c  ************************************************************
c     Check to see if the luminosity function needs to be 
c     re-calculated, and write it out if necessary.  Read the values 
c     in any case.

      call newparam(new)
  
C     write values to a file if necessary
      if (new) then
c	gg_or_gP is specified in starlight.in-- a 1 represents 2-photon
c	channels, a 2 represents vector meson channels with narrow
c	resonance, and a 3 represents vector meson channels with a wide
c	(Breit-Wigner) resonance
	if (gg_or_gP.eq.1) then
	  write(*,*) 'Calling diffLum_2gamma...'
	  call diffLum_2gamma
	elseif ((gg_or_gP.eq.2).or.(gg_or_gP.eq.3)) then 
	  write(*,*) 'Calling diffLum_vm...'
 	  call diffLum_vm
	else
	  write(*,*) 'ERROR:  Invalid entry for gg_or_gP'
	endif
      endif

c     read in the table of differential luminosity values
      call readDiffLum

c  ************************************************************
c     calculate the tables of cross section times F for
c     the various decays and calculates cross section

      call sigmacalc

c  ************************************************************
c     calculate some of the tables that will be used to decay the
c     particles

      call tablecalc

c  ************************************************************
C     open output file for text output, if desired
      if ((iout.eq.1).or.(iout.eq.2)) then
	open (unit=25,file='starlight.out',status='unknown')
      endif

      if (iout.eq.2) then
c	write header line to output file
        WRITE(25,4999) 'GENER:','STL','1.0',Z,A,Z,A,200.0,999.999,'CMS'
      endif

c  ************************************************************
c     loop through events
	write(*,*) 'Starting event loop...'
      i = 1
  
 50   if(i.gt.ievents) go to 100 

c     choose a w and a y based on the cross-section tables;
c     w is the center-of-mass energy of the 2-photon or photon-Pomeron
c     system; y is the rapidity
      if(gg_or_gp.eq.1) then
        call pickw(w)
        call picky(y)
      else
	call pickwy_vm(w,y)
      endif

C     use these to calculate momentum components and energy
	if(gg_or_gP.eq.1) call momenta(w,y,E,px,py,pz)
	if((gg_or_gP.eq.2).or.(gg_or_gP.eq.3)) then
	  tcheck = 0
	  call vmomenta(w,y,E,px,py,pz,tcheck)
	  if(tcheck.eq.1) go to 50
	endif

c     decay this channel and write out the decay products in the chosen
c     format
      call decayEvent(i,w,E,px,py,pz)
      
      go to 50 

 100  continue

      write(*,*),'Number of events processed:',i-1

c  ************************************************************

c     close the output file 
      if (iout.eq.2) write(25,5000) -999,0,0
      if ((iout.eq.1).or.(iout.eq.2)) close(unit=25)
      if (iout.eq.3) call hepend('z')

c  ************************************************************

 4999 format (A6,1x,A3,1x,A4,1x,4(I12),1x,F7.2,1x,F11.3,5x,A3) 
 5000 format ('EVENT:',3x,3(1x,i6))

      stop
      end

