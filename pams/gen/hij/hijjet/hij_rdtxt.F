C*******************************************************************************
      Integer Function HIJ_RDTXT( inlun, i_event, n_particles,
     >     geid, px, py, pz, iparent )
C*******************************************************************************
C
C     Description:
C
C     Reads one event from HIJING text file.  Assumes file has already been
C     opened.
C
C     Author:         Matthew Bloomer
C     Created:        2 April, 1994
C     Last modified:
C
C     Return values:
C        -1  Success
C         1  End of File
C         0  Not Success
C*******************************************************************************

      Implicit None

      Integer  inlun, i_event, j_particle

      Integer  n_particles, geid(*), iparent(*)
      Real     px(*), py(*), pz(*)

C*******************************************************************************

C
C Read event header for number of particles to read
C

      Read( inlun, 250, Err=920, End=950 ) n_particles, i_event
 250  Format(1x,i6,7X,i6)

C
C Read particle list
C

      Do 400 j_particle = 1, n_particles

         Read( inlun, 350, Err=920, End=950 ) geid(j_particle), 
     >        px(j_particle), 
     >        py(j_particle),
     >        pz(j_particle), 
     >        iparent(j_particle)
 350     Format(18x,i2,2x,G12.5,1x,G12.5,1x,G12.5,1x,i7)

 400  Continue

C
C Successful read
C
      hij_rdtxt = -1
      Return

C
C Error while reading event, abort
C
 920  Write  ( 6, 930) i_event
 930  Format (1X, '%HIJ_RDTXT: Error while reading event ',I11)
      hij_rdtxt = 0
      Return

C
C End of file reached
C
 950  hij_rdtxt = 1
      Return

 999  End
