* Id: mevsim.F, v 1.1 1998/03/30 14:20:50 longacre Exp $
*$Log: mevsim.F,v $
*Revision 1.1  1998/04/01 15:28:18  longacre
*first
* 
      INTEGER FUNCTION MEVSIM(par_h, par) 
      implicit none

CCC   Documentation and Description:
CCC
C     This code is intended to provide a quick means of producing
C     uncorrelated simulated events for event-by-event studies,
C     detector acceptance and efficiency studies, etc.  The
C     user selects the number of events, the one-particle distribution
C     model, the particles to include, the ranges in transverse
C     momentum, pseudorapidity and azimuthal angle, the mean
C     multiplicity for each particle type for the event run, the
C     mean temperature, Rapidity width, etc., and the standard deviations
C     for the event-to-event variation in the model parameters.
C    
C     The user may select either to have the same multiplicity per
C     particle type for each event or to let the multiplicity vary
C     randomly according to a Poisson distribution.  Similarly the
C     parameters of the one-particle distribution models may either
C     be fixed to the same value for each event or allowed to randomly
C     vary about a specified mean with a specified standard deviation
C     and assuming a Gaussian distribution.
C
C     The input is from a file, named 'mult_gen.in'.  The output is
C     particle table for a given event.
C     Header information is store in the first particle. Event header 
C     information store in the first particle includes: Impact parameter
C     and phi angle of collision; Generator number 62 mev giving VENUS
C     like central collision and number 80 giving VNI like central; 
C     Energy per nucleon-nucleon collision; Colliding system A.A;
C     Run number; Event number; Day; and Time. 
C     A log file, 'mult_gen.log' is also written which contain config. 
C     messages.
C
C     The method for generating random multiplicities and model parameter
C     values involves the following steps:
C        (1) The Poisson or Gaussian distributions are computed and
C            loaded into function f().
C        (2) The distribution f(x') is integrated from xmin to x
C            and saved from x = xmin to x = xmax.  The range and mesh
C            spaces are specified by the user.
C        (3) The integral of f is normalized to unity where 
C            integral[f(x')](at x = xmin) = 0.0
C            integral[f(x')](at x = xmax) = 1.0
C        (4) A random number generator is called which delivers values
C            between 0.0 and 1.0.  
C        (5) We consider the coordinate x (from xmin to xmax) to be
C            dependent on the integral[f].  Using the random number
C            for the selected value of integral[f] the value of x
C            is obtained by interpolation.
C
C     An interpolation subroutine from Rubin Landau, Oregon State Univ.,
C     is used to do this interpolation; it involves uneven mesh point 
C     spacing.
C
C     The method for generating the particle momenta uses the
C     standard random elimination method and involves the following
C     steps:
C        (1) The transverse momentum (pt) and pseudorapidity (eta) are 
C            randomly chosen within the specified ranges.
C        (2) The one-particle model distribution is calculated at this
C            point and its ratio to the maximum value throughout the
C            (pt,eta) region is calculated.
C        (3) Another random number is called and if less than the ratio
C            from step(2) the particle momentum is used; if not then
C            another trial value of (pt,eta) is obtained.
C        (4) This continues until the required multiplicity for the 
C            specific event and particle type has been satisfied.
C        (5) This process is repeated for the requested number of particle
C            types and events. 
C
C     Problematic parameter values are tested, bad input values are checked
C     and in some cases may be changed so that the program will not crash.
C     In some cases the code execution is stopped.
C     Some distributions and/or unusual model parameter values may cause the
C     code to hang up due to the poor performance of the "elimination"
C     method for very strongly peaked distributions.  These are tested for
C     certain problematic values and if necessary these events are aborted.
C     A message, "*** Event No.    2903 ABORTED:" for example is printed
C     in the 'mult_gen.out' file.  Temperatures .le. 0.01 GeV and rapidity
C     width parameters .le. 0.01 will cause the event to abort.
C
C     The input is described below in the 'read' statements and also in
C     the sample input file.  Some additional comments are as follows:
C
C     (1) n_events - Selected number of events in run. Can be anything
C                    .ge. 1.
C     (2) n_pid_type - Number of particle ID types to include in the
C                      particle list. e.g. pi(+) and pi(-) are counted
C                      separately.  The limit is set by parameter npid
C                      in the accompanying include file 'Parameter_values.inc'
C                      and is presently set at 10.
C     (3) model_type - equals 1,2,3,4,5 or 6 so far.  See comments in
C                      Function dNdpty to see what is calculated.
C                      The models included are:
C                    = 1, Factorized mt exponential, Gaussian rapidity model
C                    = 2, Pratt non-expanding, spherical thermal source model
C                    = 3, Bertsch non-expanding spherical thermal source model
C                    = 4, Pratt spherically expanding, thermally equilibrated
C                         source model.
C                    = 5, Factorized pt and eta distributions input bin-by-bin.
C                    = 6, Fully 2D pt,eta distributions input bin-by-bin.
C     (4) pt_cut_min,pt_cut_max - Range of transverse momentum in GeV/c.
C     (5) eta_cut_min,eta_cut_max - Pseudorapidity range
C     (6) phi_cut_min,phi_cut_max - Azimuthal angular range in degrees.
C     (7) n_stdev_mult - Number of standard deviations about the mean value
C                        of multiplicity to include in the random event-to-
C                        event selection process.  The maximum number of
C                        steps that can be covered is determined by
C                        parameter n_mult_max_steps in the accompanying
C                        include file 'Parameter_values.inc' which is
C                        presently set at 1000, but the true upper limit for
C                        this is n_mult_max_steps - 1 = 999.
C     (8) n_stdev_temp - Same, except for the "Temperature" parameter.
C     (9) n_stdev_sigma- Same, except for the rapidity width parameter.
C    (10) n_stdev_expvel - Same, except for the expansion velocity parameter.
C    (11) n_integ_pts - Number of mesh points to use in the random model
C                       parameter selection process.  The upper limit is
C                       set by parameter nmax_integ in the accompanying
C                       include file 'Parameter_values.inc' which is presently
C                       set at 100, but the true upper limit for n_integ_pts
C                       is nmax_integ - 1 = 99. 
C    (12) n_scan_pts  - Number of mesh points to use to scan the (pt,y)
C                       dependence of the model distributions looking for
C                       the maximum value.  The 2-D grid has
C                       n_scan_pts * n_scan_pts points; no limit to size of
C                       n_scan_pts.
C    (13) irand       - Starting random number seed.
C
C**************************************************************************
C    FOR MODEL_TYPE = 1,2,3 or 4:
C    Input the following 5 lines for each particle type; repeat these
C    set of lines n_pid_type times:
C
C         (a) gpid - Geant Particle ID code number
C         (b) mult_mean,mult_variance_control - Mean multiplicity and
C                                               variance control where:
C             mult_variance_control = 0 for no variance in multiplicity 
C             mult_variance_control = 1 to allow Poisson distribution for
C                                       particle multiplicities for all events.
C             Note that a hard limit exists for the maximum possible
C             multiplicity for a given particle type per event.  This is
C             determined by parameter factorial_max in accompanying include
C             file 'common_facfac.inc' and is presently set at 5000.
C         (c) Temp_mean, Temp_stdev - Temperature parameter mean (in GeV)
C             and standard deviation (Gaussian distribution assumed).
C         (d) sigma_mean, sigma_stdev - Rapidity distribution width (sigma)
C             parameter mean and standard deviation (Gaussian distribution
C             assumed).
C         (e) expvel_mean, expvel_stdev - S. Pratt expansion velocity
C             (in units of c) mean and standard deviation (Gaussian 
C             distribution assumed).
C
C**************************************************************************
C    FOR MODEL_TYPE = 5 input the following set of lines for each particle
C                       type; repeat these n_pid_type times.
C
C         (a) gpid - Geant Particle ID code number
C         (b) mult_mean,mult_variance_control - Mean multiplicity and
C                                               variance control where:
C             mult_variance_control = 0 for no variance in multiplicity
C             mult_variance_control = 1 to allow Poisson distribution for
C                                       particle multiplicities for all events.
C         (c) pt_start, eta_start - minimum starting values for pt, eta 
C                                   input for the bin-by-bin distributions.
C         (d) n_pt_bins, n_eta_bins - # input pt and eta bins.
C         (e) delta_pt, pt_bin - pt bin size and function value, repeat for
C                                each pt bin.
C         (f) delta_eta, eta_bin - eta bin size and function value, repeat
C                                  for each eta bin.
C
C         NOTE: The pt, eta ranges must fully include the requested ranges
C               in input #4 and 5 above; else the code execution will stop.
C
C         Also, variable bin sizes are permitted for the input distributions.
C
C         Also, this input distribution is used for all events in the run;
C         no fluctuations in this "parent" distribution are allowed from 
C         event-to-event.
C
C**************************************************************************
C    FOR MODEL_TYPE = 6 input the following set of lines for each particle
C                       type; repeat these n_pid_type times.
C
C         (a) gpid - Geant Particle ID code number
C         (b) mult_mean,mult_variance_control - Mean multiplicity and
C                                               variance control where:
C             mult_variance_control = 0 for no variance in multiplicity
C             mult_variance_control = 1 to allow Poisson distribution for
C                                       particle multiplicities for all events.
C         (c) pt_start, eta_start - minimum starting values for pt, eta
C                                   input for the bin-by-bin distributions.
C         (d) n_pt_bins, n_eta_bins - # input pt and eta bins.
C         (e) delta_pt - pt bin size, repeat for each pt bin. 
C         (f) delta_eta - eta bin size, repeat for each eta bin.
C         (g) i,j,pt_eta_bin(i,j) - read pt (index = i) and eta (index = j)
C                                   bin numbers and bin value for full 2D space.
C
C         NOTE: The pt, eta ranges must fully include the requested ranges
C               in input #4 and 5 above; else the code execution will stop.
C
C         Also, variable bin sizes are permitted for the input distributions.
C
C         Also, this input distribution is used for all events in the run;
C         no fluctuations in this "parent" distribution are allowed from
C         event-to-event.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*:>--------------------------------------------------------------------
*: ROUTINE:    MEVSIM
*:>--------------------------------------------------------------------
#include "mevsim.inc"
#include "mevcdes/hepevt.inc"
#include "mevcdes/headpss.inc"
C---------------------------------------------------------------------
      DOUBLE PRECISION HEPMASS
      INTEGER IFIRST
      INTEGER IDTRAN
      INTEGER I
      INTEGER II
      INTEGER J
#include "mevcdes/multi_gen.inc"
      Common/track/ pout(npid,3,factorial_max)
      real*4 pout
      Common/Geant/geant_mass(60)
      real*4 geant_mass

      integer n_events, n_pid_type, model_type, n_integ_pts, irand
      integer gpid(npid),mult_mean(npid),mult_variance_control(npid)
      integer i,j,k,pid,mult_min,mult_max,n_mult_steps(npid),ievent
      integer mult_event(npid),n_scan_pts,total_mult,n_vertices
      integer event_abort,status_abort, status
      integer iptbin, ietabin

      real*4 pt_cut_min,pt_cut_max,eta_cut_min,eta_cut_max
      real*4 phi_cut_min,phi_cut_max,n_stdev_mult,n_stdev_temp
      real*4 n_stdev_sigma,n_stdev_expvel,Temp_mean(npid)
      real*4 Temp_stdev(npid),pi,rad,mult_mean_real,mult_stdev
      real*4 mult_min_real,mult_max_real,ran
      real*4 Temp_abort, sigma_abort, bin_value

      real*4 mult_integ(n_mult_max_steps),mult_xfunc(n_mult_max_steps)
      real*4 mult_integ_save(npid,n_mult_max_steps)
      real*4 mult_xfunc_save(npid,n_mult_max_steps)
      real*4 mult_event_real

      real*4 Temp_min,Temp_max,integ(nmax_integ),xfunc(nmax_integ)
      real*4 Temp_integ_save(npid,nmax_integ)
      real*4 Temp_xfunc_save(npid,nmax_integ)
      real*4 Temp_event(npid)

      real*4 sigma_stdev(npid),sigma_mean(npid),sigma_min,sigma_max
      real*4 sigma_integ_save(npid,nmax_integ)
      real*4 sigma_xfunc_save(npid,nmax_integ)
      real*4 sigma_event(npid)

      real*4 expvel_stdev(npid), expvel_mean(npid)
      real*4 expvel_min, expvel_max
      real*4 expvel_integ_save(npid,nmax_integ)
      real*4 expvel_xfunc_save(npid,nmax_integ)
      real*4 expvel_event(npid)
      DATA IFIRST/0/
      IF(IFIRST.EQ.0) THEN

CCC  Open I/O Files:

      open(unit=4,type='old',access='sequential',name='mult_gen.in')
      open(unit=8,type='new',access='sequential',name='mult_gen.log')

CCC   File 'mult_gen.in' is the input file for the run.
CCC   File 'mult_gen.log' is a log file for the run.

CCC   Initialize Arrays to Zero:

      do i = 1,npid
         gpid(i) = 0
         mult_mean(i) = 0
         mult_variance_control(i) = 0
         n_mult_steps(i) = 0
         Temp_mean(i) = 0.0
         Temp_stdev(i) = 0.0
         sigma_mean(i) = 0.0
         sigma_stdev(i) = 0.0
         expvel_mean(i) = 0.0
         expvel_stdev(i) = 0.0
         mult_event(i) = 0
         Temp_event(i) = 0.0
         sigma_event(i) = 0.0
         expvel_event(i) = 0.0

         do j = 1,n_mult_max_steps
            mult_integ_save(i,j) = 0.0
            mult_xfunc_save(i,j) = 0.0
         end do

         do j = 1,nmax_integ
            Temp_integ_save(i,j) = 0.0
            Temp_xfunc_save(i,j) = 0.0
            sigma_integ_save(i,j) = 0.0
            sigma_xfunc_save(i,j) = 0.0
            expvel_integ_save(i,j) = 0.0
            expvel_xfunc_save(i,j) = 0.0
         end do
      end do

      do i = 1,n_mult_max_steps
         mult_integ(i) = 0.0
         mult_xfunc(i) = 0.0
      end do

      do i = 1,nmax_integ
         integ(i) = 0.0
         xfunc(i) = 0.0
      end do

      do i = 1,factorial_max
         FACLOG(i) = 0.0
      end do

      do i = 1,60
         geant_mass(i) = 0.0
      end do

      do i = 1,npid
         pt_start(i) = 0.0
         pt_stop(i)  = 0.0
         eta_start(i) = 0.0
         eta_stop(i)  = 0.0
         n_pt_bins(i) = 0
         n_eta_bins(i) = 0
         do j = 1,n_bins_max
            delta_pt(i,j)   = 0.0
            delta_eta(i,j)  = 0.0
            pt_bin(i,j)     = 0.0
            eta_bin(i,j)    = 0.0
            do k = 1,n_bins_max
               pt_eta_bin(i,j,k) = 0.0
            end do
         end do
      end do

CCC  Read Input:

      read(4,*) VSSHEP(1)              ! Run number
      read(4,*) VSSHEP(2)              ! Starting No. for the  events to generate
      read(4,*) PSSHEP(3)              ! Generator Code
      read(4,*) PSSHEP(4)              ! Energy of the Nucleon-Nucleon CM
      read(4,*) PSSHEP(5)              ! Awest.Aeast for the colliding systems
      read(4,*) n_pid_type             ! No. of Geant PID types to include
      read(4,*) model_type             ! Distribution model type (see
CCC                                    ! Function dNdpty for explanation).
      read(4,*) pt_cut_min,pt_cut_max  ! Min/Max pt range in GeV/c
      read(4,*) eta_cut_min,eta_cut_max ! Min/Max pseudorapidity range
      read(4,*) phi_cut_min,phi_cut_max ! Min/Max azimuthal angular range (deg)
      read(4,*) n_stdev_mult            ! No.(+/-) standard deviation range
CCC                                     ! for multiplicity
      read(4,*) n_stdev_temp            ! No.(+/-) st.dev. range for Temp.
      read(4,*) n_stdev_sigma           ! No.(+/-) st.dev. range for rapidity
CCC                                     ! width, sigma.
      read(4,*) n_stdev_expvel          ! No.(+/-) st.dev. range for expansion
CCC                                     ! velocity.
      read(4,*) n_integ_pts             ! No. of integration mesh points to use
CCC                                     ! for random parameter fluctuations.
      read(4,*) n_scan_pts              ! No. of pt and eta mesh points to use
CCC                                     ! in scan for maximum value of dN/dpt*dy
      read(4,*) irand                   ! Random number seed; default=12345.

CCC   Check Validity and Consistency of Input Parameters so far:

      if(n_pid_type .le. 0) n_pid_type = 1
      if(pt_cut_min .gt. pt_cut_max) then
         write(8,40) pt_cut_min,pt_cut_max
         STOP
      end if
      if(eta_cut_min .gt. eta_cut_max) then
         write(8,41) eta_cut_min,eta_cut_max
         STOP
      end if     
      if(phi_cut_min .gt. phi_cut_max) then
         write(8,42) phi_cut_min,phi_cut_max
         STOP
      end if
40    Format(//10x,'pt_cut_min = ',F7.4,' gt max = ',F7.4,' -STOP')
41    Format(//10x,'eta_cut_min = ',F7.4,' gt max = ',F7.4,' -STOP')
42    Format(//10x,'phi_cut_min = ',F7.4,' gt max = ',F7.4,' -STOP')
      if(n_stdev_mult   .lt. 0.0) n_stdev_mult   = 1.0
      if(n_stdev_temp   .lt. 0.0) n_stdev_temp   = 1.0
      if(n_stdev_sigma  .lt. 0.0) n_stdev_sigma  = 1.0
      if(n_stdev_expvel .lt. 0.0) n_stdev_expvel = 1.0
      if(n_integ_pts .le. 0) n_integ_pts = 10
      if(n_scan_pts  .le. 0) n_scan_pts  = 10

      if(irand .le. 0) irand = 12345
      if(n_pid_type .gt. npid) then
         write(8,10) n_pid_type, npid
10       Format(//10x,'No. requested PID types = ',I7,
     1   ', exceeds maximum of ',I7,'; reset')
         n_pid_type = npid
      end if
      if(model_type .lt. 0 .or. model_type .gt. 6) then
         write(8,11) model_type
11       Format(/10x,'model_type = ',I5,' is not allowed; STOP')
         STOP
      end if
      if(n_integ_pts .gt. nmax_integ) then
         write(8,12) n_integ_pts, nmax_integ
12       Format(/10x,'No. integ. pts = ',I7,
     1   ', exceeds maximum of ',I7,'; reset')
         n_integ_pts = nmax_integ
      end if

CCC   FOR MODEL_TYPE = 1,2,3 or 4; 
CCC   Repeat the following lines of input for each particle ID type:
     
      IF(model_type.ge.1 .and. model_type.le.4) then

      do pid = 1,n_pid_type

         read(4,*) gpid(pid)            ! Geant Particle ID code number

         read(4,*) mult_mean(pid), mult_variance_control(pid)
CCC      Mean multiplicity and variance control where:
CCC           mult_variance_control = 0 for no variance in multiplicity
CCC           mult_variance_control = 1 to allow Poisson distribution for 
CCC                                     particle multiplicities for all events.

         read(4,*) Temp_mean(pid), Temp_stdev(pid)
CCC      Temperature parameter mean (in GeV) and standard deviation (Gaussian
CCC      distribution assumed).

         read(4,*) sigma_mean(pid), sigma_stdev(pid)
CCC      Rapidity distribution width (sigma) parameter mean and standard
CCC      deviation (Gaussian distribution assumed).

         read(4,*) expvel_mean(pid), expvel_stdev(pid)
CCC      S. Pratt expansion velocity (in units of c) mean and standard
CCC      deviation (Gaussian distribution assumed).

CCC      Check Validity and Consistency of Input Parameters

         if(Temp_mean(pid).le.0.0 .or. Temp_stdev(pid).lt.0.0) then
         write(8,45) pid,Temp_mean(pid),Temp_stdev(pid)
45       Format(//10x,'For pid # ',I7,', Temp_mean or Temp_stdev= ',
     1   2F7.4,' - not valid -STOP')
         STOP
         end if
         if(sigma_mean(pid).le.0.0 .or. sigma_stdev(pid).lt.0.0) then
         write(8,46) pid,sigma_mean(pid),sigma_stdev(pid)
46       Format(//10x,'For pid # ',I7,', sigma_mean or sigma_stdev= ',
     1   2F7.4,' - not valid -STOP')
         STOP
         end if
         if(expvel_mean(pid).lt.0.0.or.expvel_stdev(pid).lt.0.0) then
         write(8,47) pid,expvel_mean(pid),expvel_stdev(pid)
47       Format(//10x,'For pid #',I7,',expvel_mean or expvel_stdev=',
     1   2F7.4,' - not valid -STOP')
         STOP
         end if

      end do  !  End of loop over of Particle ID types.

      ELSE IF (model_type .eq. 5) then

CCC      Input for Factorized pt, eta bin-by-bin distribution:

         do pid = 1,n_pid_type
            read(4,*) gpid(pid)
            read(4,*) mult_mean(pid), mult_variance_control(pid)
            read(4,*) pt_start(pid), eta_start(pid)
            read(4,*) n_pt_bins(pid), n_eta_bins(pid)
            do i = 1,n_pt_bins(pid)
               read(4,*) delta_pt(pid,i), pt_bin(pid,i)
            end do
            do i = 1,n_eta_bins(pid)
               read(4,*) delta_eta(pid,i), eta_bin(pid,i)
            end do

CCC      Evaluate grid and find maximum values of pt and eta for input bins:

            pt_stop(pid) = pt_start(pid)
            do i = 1,n_pt_bins(pid)
            pt_stop(pid) = pt_stop(pid) + delta_pt(pid,i)
            pt_bin_mesh(pid,i) = pt_stop(pid)
            end do
            eta_stop(pid) = eta_start(pid)
            do i = 1,n_eta_bins(pid)
            eta_stop(pid) = eta_stop(pid) + delta_eta(pid,i)
            eta_bin_mesh(pid,i) = eta_stop(pid)
            end do

CCC      Check ranges of pt and eta coverage:

            if(pt_cut_min .lt. pt_start(pid)) then
               write(8,50) pt_cut_min,pt_start(pid)
50             Format(//10x,'pt_cut_min = ',F10.7,' .lt. pt_start =',
     1         F10.7,' - STOP')
               STOP
            end if
            if(pt_cut_max .gt. pt_stop(pid)) then
               write(8,51) pt_cut_max,pt_stop(pid)
51             Format(//10x,'pt_cut_max = ',F10.7,' .gt. pt_stop =',
     1         F10.7,' - STOP')
               STOP
            end if
            if(eta_cut_min .lt. eta_start(pid)) then
               write(8,52) eta_cut_min,eta_start(pid)
52             Format(//10x,'eta_cut_min = ',F10.7,'.lt.eta_start =',
     1         F10.7,' - STOP')
               STOP
            end if
            if(eta_cut_max .gt. eta_stop(pid)) then
               write(8,53) eta_cut_max,eta_stop(pid)
53             Format(//10x,'eta_cut_max = ',F10.7,'.gt.eta_stop =',
     1         F10.7,' - STOP')
               STOP
            end if
         end do ! End loop over particle ID types.

      ELSE IF (model_type .eq. 6) then

CCC      Input for Full 2D (pt,eta) bin-by-bin distribution:

         do pid = 1,n_pid_type
            read(4,*) gpid(pid)
            read(4,*) mult_mean(pid), mult_variance_control(pid)
            read(4,*) pt_start(pid), eta_start(pid)
            read(4,*) n_pt_bins(pid), n_eta_bins(pid)
            do i = 1,n_pt_bins(pid)
               read(4,*) delta_pt(pid,i)
            end do
            do i = 1,n_eta_bins(pid)
               read(4,*) delta_eta(pid,i)
            end do

CCC      Evaluate grid and find maximum values of pt and eta for input bins:

            pt_stop(pid) = pt_start(pid)
            do i = 1,n_pt_bins(pid)
            pt_stop(pid) = pt_stop(pid) + delta_pt(pid,i)
            pt_bin_mesh(pid,i) = pt_stop(pid)
            end do
            eta_stop(pid) = eta_start(pid)
            do i = 1,n_eta_bins(pid)
            eta_stop(pid) = eta_stop(pid) + delta_eta(pid,i)
            eta_bin_mesh(pid,i) = eta_stop(pid)
            end do

CCC   Load 2D bin-by-bin distribution:

            do i = 1,n_pt_bins(pid)*n_eta_bins(pid)
               read(4,*) iptbin,ietabin,bin_value
               pt_eta_bin(pid,iptbin,ietabin) = bin_value
            end do

CCC      Check ranges of pt and eta coverage:

            if(pt_cut_min .lt. pt_start(pid)) then
               write(8,50) pt_cut_min,pt_start(pid)
               STOP
            end if
            if(pt_cut_max .gt. pt_stop(pid)) then
               write(8,51) pt_cut_max,pt_stop(pid)
               STOP
            end if
            if(eta_cut_min .lt. eta_start(pid)) then
               write(8,52) eta_cut_min,eta_start(pid)
               STOP
            end if
            if(eta_cut_max .gt. eta_stop(pid)) then
               write(8,53) eta_cut_max,eta_stop(pid)
               STOP
            end if
         end do ! End loop over particle ID types.

      END IF !  End of MODEL_TYPE Options input:
      
CCC   Check some more input parameters; Set constants, and load data tables:

      do pid = 1,n_pid_type
      if(gpid(pid).le.0 .or. gpid(pid).gt.48) then
      write(8,43) pid,gpid(pid)
43    Format(//10x,'For pid # ',I7,', gpid = ',I7,' - not valid -STOP')
      STOP
      end if
      if(mult_variance_control(pid) .lt. 0  .or.
     1   mult_variance_control(pid) .gt. 1)
     2   mult_variance_control(pid) = 0
      if(mult_mean(pid) .le. 0) then
      write(8,44) pid,mult_mean(pid)
44    Format(//10x,'For pid # ',I7,', mult_mean= ',I7,
     1   ' - not valid -STOP')
      STOP
      end if
      end do ! End Particle ID input parameter check and verification loop.

      pi = 3.141592654
      rad = 180.0/pi
      Temp_abort  = 0.01
      sigma_abort = 0.01

CCC   Load particle properties array; mass in GeV:

      Call Particle_prop

CCC   Load log(n!) values into Common/FACFAC/

      Call FACTORIAL

CCC   Obtain integrals of 1D distribution functions of multiplicity
CCC   (Poisson, integer) and parameters (Gaussian, real*4) for the
CCC   particle model distributions, for each particle ID type.
CCC   These integrated quantities are then used with random number
CCC   selection to generate a distribution of multiplicities and
CCC   parameter values for each event in the run.

      do pid = 1,n_pid_type

         if(mult_variance_control(pid) .ne. 0) then
            mult_mean_real = float(mult_mean(pid))
            mult_stdev = sqrt(mult_mean_real)
            Call MinMax(mult_mean_real,mult_stdev,n_stdev_mult,
     1                  mult_min_real,mult_max_real)
            mult_min = int(mult_min_real)
            mult_max = int(mult_max_real + 1)
            n_mult_steps(pid) = mult_max - mult_min + 1
            if((n_mult_steps(pid) + 1) .gt. n_mult_max_steps) then
               write(8,20) n_mult_steps(pid),n_mult_max_steps
20             Format(//10x,'No. steps in multiplicity integral = ',
     1         I7,' + 1, exceeds max no. of ',I7,'; reset')
            mult_min = mult_mean(pid) - (n_mult_max_steps - 1)/2
            mult_max = mult_mean(pid) + (n_mult_max_steps - 1)/2
            n_mult_steps(pid) = mult_max - mult_min + 1
            end if
            if((mult_max + 1) .gt. factorial_max) then
               write(8,30) mult_max, factorial_max
30             Format(//10x,'In Poisson multiplicity distribution,',
     1         ' max = ',I7,', exceeds limit of ',I7,' - STOP')
               STOP
            end if

            Call Poisson(mult_min,mult_max,mult_mean(pid),
     1      n_mult_steps(pid),mult_integ,mult_xfunc,n_mult_max_steps)
            do i = 1,n_mult_steps(pid) + 1
               mult_integ_save(pid,i) = mult_integ(i)
               mult_xfunc_save(pid,i) = mult_xfunc(i)
            end do
         else if (mult_variance_control(pid) .eq. 0) then
            mult_event(pid) = mult_mean(pid)
         end if

         if(model_type .le. 4) then
         if(Temp_stdev(pid) .ne. 0.0) then
            Call MinMax(Temp_mean(pid),Temp_stdev(pid),n_stdev_temp,
     1                  Temp_min, Temp_max)
            Call Gaussian(Temp_min,Temp_max,Temp_mean(pid),
     1         Temp_stdev(pid),n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               Temp_integ_save(pid,i) = integ(i)
               Temp_xfunc_save(pid,i) = xfunc(i)
            end do
         else if(Temp_stdev(pid) .eq. 0.0) then
            Temp_event(pid) = Temp_mean(pid)
         end if

         if(sigma_stdev(pid) .ne. 0.0) then
            Call MinMax(sigma_mean(pid),sigma_stdev(pid),n_stdev_sigma,
     1                  sigma_min, sigma_max)
            Call Gaussian(sigma_min,sigma_max,sigma_mean(pid),
     1         sigma_stdev(pid),n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               sigma_integ_save(pid,i) = integ(i)
               sigma_xfunc_save(pid,i) = xfunc(i)
            end do
         else if(sigma_stdev(pid) .eq. 0.0) then
            sigma_event(pid) = sigma_mean(pid)
         end if

         if(expvel_stdev(pid) .ne. 0.0) then
          Call MinMax(expvel_mean(pid),expvel_stdev(pid),n_stdev_expvel,
     1                  expvel_min, expvel_max)
            Call Gaussian(expvel_min,expvel_max,expvel_mean(pid),
     1         expvel_stdev(pid),n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               expvel_integ_save(pid,i) = integ(i)
               expvel_xfunc_save(pid,i) = xfunc(i)
            end do
         else if(expvel_stdev(pid) .eq. 0.0) then
            expvel_event(pid) = expvel_mean(pid)
         end if 
         end if ! End model_type .le. 4 options.

      end do  !  End of PID Loop:

CCC   Write Run Header Output:

      write(8,200)
      write(8,201) n_events,n_pid_type
      if(model_type .eq. 1) write(8,202)
      if(model_type .eq. 2) write(8,203)
      if(model_type .eq. 3) write(8,204)
      if(model_type .eq. 4) write(8,205)
      if(model_type .eq. 5) write(8,2051)
      if(model_type .eq. 6) write(8,2052)
      write(8,206) pt_cut_min, pt_cut_max
      write(8,207) eta_cut_min, eta_cut_max
      write(8,208) phi_cut_min, phi_cut_max
      write(8,209) n_stdev_mult,n_stdev_temp,n_stdev_sigma,
     1             n_stdev_expvel
      write(8,210) n_integ_pts
      write(8,211) n_scan_pts
      write(8,212) irand
      write(8,213) (gpid(pid),pid = 1,n_pid_type)
      write(8,214) (mult_mean(pid),pid = 1,n_pid_type)
      write(8,215) (mult_variance_control(pid),pid = 1,n_pid_type)
      if(model_type .le. 4) then
         write(8,216) (Temp_mean(pid),pid = 1,n_pid_type)
         write(8,217) (Temp_stdev(pid),pid = 1,n_pid_type)
         write(8,218) (sigma_mean(pid),pid = 1,n_pid_type)
         write(8,219) (sigma_stdev(pid),pid = 1,n_pid_type)
         write(8,220) (expvel_mean(pid),pid = 1,n_pid_type)
         write(8,221) (expvel_stdev(pid),pid = 1,n_pid_type)
      else if(model_type .ge. 5) then
         write(8,222) (n_pt_bins(pid),pid = 1,n_pid_type)
         write(8,223) (n_eta_bins(pid),pid = 1,n_pid_type)
         write(8,224) (pt_start(pid),pid = 1,n_pid_type)
         write(8,225) (pt_stop(pid),pid = 1,n_pid_type)
         write(8,226) (eta_start(pid),pid = 1,n_pid_type)
         write(8,227) (eta_stop(pid),pid = 1,n_pid_type)
      end if

200   Format('***  Multiplicity Generator Run Header ***')
201   Format('*    Number of events = ',I7,'; number of Particle ID',
     1       ' types = ',I5)
202   Format('* Factorized mt exponential, Gaussian rapidity model')
203   Format('* Pratt non-expanding, spherical thermal source model')
204   Format('* Bertsch non-expanding spherical thermal source model')
205   Format('* Pratt spherically expanding, thermally equilibrated ',
     1        'source model')
2051  Format('* Factorized pt,eta bin-by-bin Distribution')
2052  Format('* Full 2D pt,eta bin-by-bin Distribution')
206   Format('* Min, Max range in pt              = ', 2G12.5)
207   Format('* Min, Max range in pseudorapidity  = ', 2G12.5)
208   Format('* Min, Max range in azimuthal angle = ', 2G12.5)
209   Format('* No. std. dev. range used for mult and parameters = ',
     1       4F5.2)
210   Format('* No. integration points for parameter variance = ',
     1       I6)
211   Format('* No. of dN/dp(pt,y) scanning grid points to find ',
     1   'maximum = ', I6)
212   Format('* Starting Random Number Seed = ',I10)
213   Format('* Geant PID:          ',10I7)
214   Format('* Mean Multiplicity:  ',10I7)
215   Format('* Mult. Var. Control: ',10I7)
216   Format('* Mean Temperature:   ',10F7.4)
217   Format('* Std. Dev. Temp:     ',10F7.4)
218   Format('* Mean Rap. sigma:    ',10F7.4)
219   Format('* Std. Dev. y-sigma:  ',10F7.4)
220   Format('* Mean expansion vel: ',10F7.4)
221   Format('* Std. Dev. exp. vel: ',10F7.4)
222   Format('* No. input pt bins:  ',10I7)
223   Format('* No. input eta bins: ',10I7)
224   Format('* Input pt start:     ',10F7.4)
225   Format('* Input pt stop:      ',10F7.4)
226   Format('* Input eta start:    ',10F7.4)
227   Format('* Input eta stop:     ',10F7.4)

CCC  END Run Header Output
         IFIRST=1
      ENDIF

C**************************************************************
C                                                            **
C                     START EVENT LOOP                       **
C                                                            **
C**************************************************************
      CALL DATIME(I,J)
      VSSHEP(2) = VSSHEP(2) + 1
      VSSHEP(3) = FLOAT(I)
      VSSHEP(4) = FLOAT(J)
      par(1).isthep = 10
      par(1).idhep = 9999999
      par(1).jmohep(1) = 0
      par(1).jmohep(2) = 0
      par(1).jdahep(1) = 0
      par(1).jdahep(2) = 0
      par(1).phep(1) = 0.0
      par(1).phep(2) = 0.0
      par(1).phep(3) = psshep(3)
      par(1).phep(4) = psshep(4)
      par(1).phep(5) = psshep(5)
      par(1).vhep(1) = vsshep(1)
      par(1).vhep(1) = vsshep(1)
      par(1).vhep(2) = vsshep(2)
      par(1).vhep(3) = vsshep(3)
      par(1).vhep(4) = vsshep(4)
         do pid = 1,n_pid_type
            if(mult_variance_control(pid) .ne. 0) then
               do i = 1,n_mult_steps(pid) + 1
               mult_integ(i) = mult_integ_save(pid,i)
               mult_xfunc(i) = mult_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),mult_integ,mult_event_real,
     1         mult_xfunc,n_mult_steps(pid)+1,1,5,
     2         n_mult_steps(pid)+1,1)
               mult_event(pid) = mult_event_real
            end if

         if(model_type .le. 4) then

            if(Temp_stdev(pid) .ne. 0.0) then
               do i = 1,n_integ_pts + 1
               integ(i) = Temp_integ_save(pid,i)
               xfunc(i) = Temp_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),integ,Temp_event(pid),xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
            end if

            if(sigma_stdev(pid) .ne. 0.0) then
               do i = 1,n_integ_pts + 1
               integ(i) = sigma_integ_save(pid,i)
               xfunc(i) = sigma_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),integ,sigma_event(pid),xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
            end if

            if(expvel_stdev(pid) .ne. 0.0) then
               do i = 1,n_integ_pts + 1
               integ(i) = expvel_integ_save(pid,i)
               xfunc(i) = expvel_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),integ,expvel_event(pid),xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
            end if
         end if
         end do !  End Particle ID setup Loop

         event_abort = 1

         if(model_type .le. 4) then
CCC      Check validity of Temperature and sigma parameters:
         do pid = 1,n_pid_type
            if(Temp_event(pid) .le. Temp_abort) event_abort = -86
            if(sigma_event(pid).le.sigma_abort) event_abort = -86
         end do
         end if

         if(event_abort .eq. 1) then

            total_mult = 0
            do pid = 1,n_pid_type
            total_mult = total_mult + mult_event(pid)
            end do
            n_vertices = 0
            status_abort = 1
            do pid = 1,n_pid_type
CCC   For the specified Geant PID, multiplicity, model type, temperature,
CCC   rapidity width (sigma), and expansion velocity parameter, generate
CCC   random track list:

            Call track_gen(gpid(pid),mult_event(pid),model_type,
     1      Temp_event(pid),sigma_event(pid),expvel_event(pid),
     2      pt_cut_min,pt_cut_max,eta_cut_min,eta_cut_max,
     3      phi_cut_min,phi_cut_max,n_scan_pts,irand,rad,pid,status)
            if(status .eq. -86) status_abort = -86
            end do
         end if

         if(event_abort.eq.1 .and. status_abort.eq.1) then
CCC Event Header Output:
           IF(total_mult.GT.par_h.maxlen) STOP
           write(8,213) (gpid(pid),pid = 1,n_pid_type)
           write(8,231) (mult_event(pid),pid = 1,n_pid_type)
           if(model_type .le. 4) then
              write(8,232) (Temp_event(pid),pid = 1,n_pid_type)
              write(8,233) (sigma_event(pid),pid = 1,n_pid_type)
              write(8,234) (expvel_event(pid),pid = 1,n_pid_type)
           end if

230   Format(/'***  Next Event:  No. ',I7,'; Total # tracks = ',I10)
231   Format('* Multiplicity:         ',10I7)
232   Format('* Temperature:          ',10F7.4)
233   Format('* Rapidity Dist. sigma: ',10F7.4)
234   Format('* Expansion Velocity:   ',10F7.4)

           write(8,235) ievent,total_mult,n_vertices
235        Format('EVENT:',3x,3(1x,I6))
           write(8,236)
           write(8,237)
236        Format('*** Tracks for this event ***')
237        Format('* Geant PID #      px          py          pz')
CCC   End Event Header Output:

CCC   Output track kinematics for ievent and pid:
        II = 1
           do pid = 1,n_pid_type
           do i = 1,mult_event(pid)
        II = II + 1
        idhep=idtran(gpid(pid))
        phep(1) = pout(pid,1,i)
        phep(2) = pout(pid,2,i)
        phep(3) = pout(pid,3,i)
        phep(4)=dsqrt(phep(1)**2+phep(2)**2+phep(3)**2+hepmass(idhep)**2)
        par(ii).isthep     = 1           !  status code
        par(ii).idhep      = idhep       !  particle identity
        par(ii).jmohep(1)  = 0           !  pointer(s) to position 
        par(ii).jmohep(2)  = 0           !  where the mothers(s) stored   
        par(ii).jdahep(1)  = 0           !  pointers to position of 
        par(ii).jdahep(2)  = 0           !  the first/last daughter  
        DO J = 1, 4
          par(ii).phep(j)  = phep(j)     !  p4 
          par(ii).vhep(j)  = 0.0         !  production vertex (mm) and time (mm/c) 
        END DO
        par(ii).phep(5)    = hepmass(idhep)!  mass (GeV)  
           end do
           end do
CCC   End track kinematics output.
      par_h.nok = II
      MEVSIM=STAFCV_OK

         else
           write(8,250) ievent, event_abort, status_abort
250        Format(/'*** Event No. ',I7,' ABORTED: event_abort,',
     1     'status_abort = ',2I7)
         end if
      RETURN
      END


