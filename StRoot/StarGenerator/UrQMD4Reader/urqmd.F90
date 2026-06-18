MODULE UrQMDmod
  use, intrinsic :: iso_c_binding
  include 'coms90.inc'
  include 'comres90.inc'
  include 'options90.inc'
  include 'colltab90.inc'
  include 'inputs90.inc'
  include 'newpart90.inc'
  include 'boxinc90.inc'

  integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, noc, it1,it2
  real*8 sqrts,otime,xdummy,st
  external sqrts
  logical isstable
  integer stidx,CTOsave
  real*8,target:: Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
  common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
  integer cti1sav,cti2sav

  real*8 thydro_start,thydro,nucrad
  logical lhydro    

  external nucrad

  contains
 
  subroutine urqmd_init()  bind(c, name="urqmd_init_")

!c
!c     numerical/technical initialisation
!c
      call uinit(0)

!c
!c  Main program
!c

      mc=0
      mp=0
      noc=0

  end subroutine urqmd_init

  subroutine urqmd_make() bind(c, name="urqmd_make_")

!     start event here
!    

!     time is the system time at the BEGINNING of every timestep
      time = 0.0
!hp hydro flag, hydro should be called only once
      lhydro=.true.

!     initialize random number generator
!     call auto-seed generator only for first event and if no seed was fixed
      if(.not.firstseed.and.(.not.fixedseed)) then
         ranseed=-(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed=.false.
      endif

      write(6,*)'event# ',event,ranseed

!tr Init resonance reconstruction (only f13)      
      if (.not.bf13.and.CTOption(68).eq.1) then
       call init_resrec
      endif

!
!     initialisation of physics quantities 
!
      call init
      call init_eccentricity
!bb if we are reading old events, check the success of the read-in:
      if (CTOption(40).ne.0.and.(.not.success)) then
         return !exit  (used to be in a loop)
      endif

!hp hydro switch      
      if (CToption(45).eq.1) then
!hp hydro start time (nuclei have passed each other)
!hp ebeam is only the kinetic energy 
!hp CTParam(65) is useful for the variation of the start time 
!hp default value is one

         thydro_start=CTParam(65)*2.d0 * nucrad(ap) !$$$ *nucrad(Ap)*sqrt(2.d0*emnuc/ebeam)
         
         thydro_start=CTParam(65)*2.d0*nucrad(Ap)*sqrt(2.d0*emnuc/ebeam)

       write(*,300) 'hydro starts after',thydro_start
!hp lower limit for hydro start time
       if(thydro_start.lt.CTParam(63)) then
        thydro_start=CTParam(63)
        write(6,300) '... extended to',CTParam(63)
       end if
      end if
 300  format(a18,x,f5.2,' fm/c')

! old time if an old fort.14 is used 
      if(CTOption(40).ne.0)time=acttime

! output preparation

! write headers to file
      call output(13)
      call output(14)
      call output(15)
      call output(16)
      if(event.eq.1) then
        call output(17)
        call osc_header
        call osc99_header
      endif
      call osc99_event(-1)

! for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)

!     participant/spectator model:
      if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

!     compute time of output
      otime = outsteps*dtimestep
      
!     reset time step counter
      steps = 0
      
!     loop over all timesteps

      do 20  steps=1,nsteps
         
!     store coordinates in arrays with *_t
!     this is needed for MD type propagation
         if (eos.ne.0) then
            do 23 j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
 23         continue
         end if
         
!     we are at the beginning of the timestep, set current time (acttime) 
         acttime = time
!     option for MD without collision term
         if(CTOption(16).ne.0) goto 103

!     Load collision table with next collisions in current timestep
         call colload

!     check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then
            
!     entry-point for collision loop in case of full colload after every coll. 
 101        continue
            k = 0
            
!     normal entry-point for collision loop 
 100        continue
            
!     get next collision
            call getnext(k)
            

!     exit collision loop if no collisions are left
            if (k.eq.0) goto 102
!     hp call hydro if start time is reached
            if(CTOption(45).eq.1)then         
               if(cttime(k).gt.thydro_start.and.lhydro)then

              if(CTOption(62).eq.1)then 
               call prepout
               call file14out(0)
               call restore
              endif
             
                  st=thydro_start-acttime
                  call cascstep(acttime,st)
!     hp all particle arrays will be modified by hydro
                  write(*,*)'starting hydro'
                  call hydro(thydro_start,thydro)
                  acttime=thydro_start
                  lhydro=.false.
                  if(CTOption(50).eq.1) goto 10
                  if(thydro.gt.1.d-8.or.CTOption(48).eq.1)then
!     hp full update of collision table
                     call colload
                     
                     go to 101
                  end if
               end if    
            end if 
        

!     propagate all particles to next collision time
!     store actual time in acttime, propagation time st=cttime(k)-acttime
            st=cttime(k)-acttime
            call cascstep(acttime,st)
!     new actual time (for upcoming collision)
            acttime = cttime(k)
            
!     perform collision 
            
            if(cti2(k).gt.0.and. &
               abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(6,*)' ***(E) wrong collision update (col) ***'
               write(6,*)k,cti1(k),cti2(k), &
                   ctsqrts(k),sqrts(cti1(k),cti2(k))
            else if(cti2(k).eq.0.and. &
                 abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
               write(6,*)' *** main(W) wrong collision update (decay)'
               write(6,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)), &
                  fmass(cti1(k)),ctsqrts(k)
            endif
            
            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))

!     store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))

!     increment "dirty" collision counter
            if(cti2(k).gt.0)then !scatter
               mc=mc+1
            endif



!     perform scattering/decay
            cti1sav = cti1(k)               
            cti2sav = cti2(k)    
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k), &
                ctcolfluc(k))
            


!     
!     update collision table 
!     
!     normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
!     new collision partners for pauli-blocked states (nexit=0)
                  if (cti1(k).ne.cti1sav.or.cti2(k).ne.cti2sav) then
                     cti1(k) = cti1sav 
                     cti2(k) = cti2sav 
                  endif
                  call collupd(cti1(k),1)
                  if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                  ncharge=0
!     new collision partners for scattered/produced particles (nexit><0)
                  do 30 i=1,nexit
!     ncharge is used for charge conservation check
                     ncharge=ncharge+charge(inew(i))
                     call collupd(inew(i),1)
 30               continue
                  
!     charge conservation check
                  if(ocharge.ne.ncharge) then
                     write(6,*)'ch-conservation error coll/dec ',ctag
                     write(6,*)'   it1:',it1,'   it2:',it2
                     write(6,*)'   ch:',ocharge,ncharge
                     write(6,*)'cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)'
                     write(6,*)cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)
                  endif
               endif
               
!     update collisions for partners of annihilated particles
               do 55 ii=1,nsav
                  call collupd(ctsav(ii),1)
 55            continue
               nsav=0
            else                ! (CTOption(17).ne.0)
!     full collision load
               call colload
            endif
            
            if (CTOption(17).eq.0) goto 100
            goto 101

!     this is the point to jump to after all collisions in the timestep
!     have been taken care of
 102        continue

         endif                  ! (nct.gt.0)
         
!     After all collisions in the timestep are done, propagate to end of 
!     the timestep.
         
!     point to jump to in case of MD without collision term
 103     continue
        
!     increment timestep
         time = time+dtimestep
         
!     After all collisions in the timestep are done, propagate to end of 
!     the timestep.
         call cascstep(acttime,time-acttime)

!     in case of potential interaction, do MD propagation step
         if (eos.ne.0) then
            
! set initial conditions for MD propagation-step
            do 24 j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
 24         continue

! now molecular dynamics trajectories
            call proprk(time,dtimestep)

         end if                 ! (eos.ne.0)
         
!     perform output if desired
         if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then 
            if(CTOption(28).eq.2)call spectrans(otime)
            if(CTOption(62).eq.1)call prepout
            call file14out(steps)
            if(CTOption(64).eq.1)then
            call file13out(steps)
            endif
            if(CTOption(62).eq.1)then
               call restore
               call colload
            endif
            if(CTOption(55).eq.1)then
               call osc_vis(steps)
            endif  
         endif                  ! output handling
         
 20   continue                  ! time step loop
      
      
!     e
      acttime=time
!     optional decay of all unstable particles before final output
!     DANGER: pauli-blocked decays are not performed !!!
      if(CTOption(18).eq.0) then
!     no do-loop is used because npart changes in loop-structure
         i=0
         nct=0
         actcol=0
!     disable Pauli-Blocker for final decays
         CTOsave=CTOption(10)
         CTOption(10)=1
!     decay loop structure starts here
 40      continue
         i=i+1
         
!     is particle unstable
         if(dectime(i).lt.1.d30) then
 41         continue
            isstable = .false.
               do 44 stidx=1,nstable
                  if (ityp(i).eq.stabvec(stidx)) then
!                     write (6,*) 'no decay of particle ',ityp(i)
                     isstable = .true.
                  endif
 44            enddo
               if (.not.isstable) then
!     perform decay
                  call scatter(i,0,0.d0,fmass(i),xdummy)
!     backtracing if decay-product is unstable itself
                  if(dectime(i).lt.1.d30) goto 41
               endif
            endif
!     check next particle
            if(i.lt.npart) goto 40
         endif ! final decay
         CTOption(10)=CTOsave
! final output

         if(CTOption(64).eq.1)call coalescence

         call file13out(nsteps)
         if(CTOption(50).eq.0)then
          call file14out(nsteps)
         end if
         call file16out
         if(CTOption(50).eq.0.and.CTOption(55).eq.0)then
          call osc_event
         end if
         if(CTOption(50).eq.0.and.CTOption(55).eq.1)then
          call osc_vis(nsteps)
         end if
         call osc99_event(1)
         call osc99_eoe
      
         mp=mp+npart
         if(ctag.eq.0)then
            write(*,*)'(W) No collision in event ',event
            noc=noc+1
         endif

!     end of event loop
         10 continue

      write(6,301)'no. of collisions = ',mc/dble(nevents),' (per event)'
      write(6,301)'final particles   = ',mp/dble(nevents),' (per event)'
      write(6,302)'empty events      : ', noc,noc*1d2/dble(nevents)
 301  format(a19,f8.1,a12)
 302  format(a19,i8,' = ',f5.1,'%')


  end subroutine urqmd_make

  subroutine urqmd_finish
  end subroutine urqmd_finish



  Function address_of_energies() bind(c,name="address_of_energies")
      use iso_c_binding
      Type(C_Ptr) :: address_of_energies
      address_of_energies = c_loc( Ekinbar )
      return
    End Function address_of_energies

    Function address_of_sys() bind(c,name="address_of_sys")
      use iso_c_binding
      Type(C_Ptr) :: address_of_sys
      address_of_sys = c_loc ( npart )
      return
    End Function address_of_sys
      
    Function address_of_rsys() bind(c,name="address_of_rsys")
      use iso_c_binding
      Type(C_Ptr) :: address_of_rsys
      address_of_rsys = c_loc ( time )
      return
    End Function address_of_rsys

    Function address_of_cuts() bind(c,name="address_of_cuts")
      use iso_c_binding
      Type(C_Ptr) :: address_of_cuts
      address_of_cuts = c_loc ( cutmax )
      return
    End Function address_of_cuts

    Function address_of_spdata() bind(c,name="address_of_spdata")
      use iso_c_binding
      type(c_ptr):: address_of_spdata
      address_of_spdata = c_loc ( spx(1) )
      return
    End Function address_of_spdata
    
    Function address_of_isys() bind(c,name="address_of_isys")
      use iso_c_binding
      type(c_ptr):: address_of_isys
      address_of_isys = c_loc ( spin(1) )
      return
    End Function address_of_isys

    Function address_of_coor() bind(c,name="address_of_coor")
      use iso_c_binding
      type(c_ptr)::address_of_coor
      address_of_coor = c_loc ( r0(1) )
      return
    End Function address_of_coor
      
    Function address_of_frag() bind(c,name="address_of_frag")
      use iso_c_binding
      type(c_ptr)::address_of_frag
      address_of_frag = c_loc ( tform(1) )
      return
    End Function address_of_frag
      
    Function address_of_aios() bind(c,name="address_of_aios")
      use iso_c_binding
      type(c_ptr):: address_of_aios
      address_of_aios = c_loc ( airx(1) )
      return
    End Function address_of_aios

    Function address_of_pots() bind(c,name="address_of_pots")
      use iso_c_binding
      type(c_ptr):: address_of_pots
      address_of_pots = c_loc ( Cb0 )
      return
    End Function address_of_pots

    Function address_of_scoor() bind(c,name="address_of_scoor")
      use iso_c_binding
      type(c_ptr):: address_of_scoor
      address_of_scoor = c_loc ( r0s(1) )
      return
    End Function address_of_scoor

    Function address_of_sisys() bind(c,name="address_of_sisys")
      use iso_c_binding
      type(c_ptr)::address_of_sisys
      address_of_sisys = c_loc(sspin(1))
      return
    End Function address_of_sisys

    Function address_of_ssys() bind(c,name="address_of_ssys")
      use iso_c_binding
      type(c_ptr):: address_of_ssys
      address_of_ssys = c_loc ( nspec )
      return
    End Function address_of_ssys

    Function address_of_rtdelay() bind(c,name="address_of_rtdelay")
      use iso_c_binding
      type(c_ptr):: address_of_rtdelay
      address_of_rtdelay = c_loc ( p0td(1,1) )
      return
    End Function address_of_rtdelay

    Function address_of_itdelay() bind(c,name="address_of_itdelay")
      use iso_c_binding
      type(c_ptr):: address_of_itdelay
      address_of_itdelay = c_loc ( ityptd(1,1) )
      return
    End Function address_of_itdelay

    Function address_of_svinfo()  bind(c,name="address_of_svinfo")
      use iso_c_binding
      type(c_ptr):: address_of_svinfo
      address_of_svinfo = c_loc ( itypt(1) )
      return
    End Function address_of_svinfo

    Function address_of_ffermi() bind(c,name="address_of_ffermi")
      use iso_c_binding
      type(c_ptr):: address_of_ffermi
      address_of_ffermi = c_loc ( ffermpx(1) )
      return
    End Function address_of_ffermi

    Function address_of_peq() bind(c,name="address_of_peq")
      use iso_c_binding
      type(c_ptr):: address_of_peq
      address_of_peq = c_loc ( peq1 )
      return
    End Function address_of_peq


END MODULE UrQMDmod


