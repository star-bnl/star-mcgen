c-*-Fortran-*-
      program read_zebra_evt
c     =================================================================
c     
c     Program to read events from an Exchange format ZEBRA tape.
c     The user must supply information about the Program which created
c     the tape since the bank format for the run_divison is dependent
c     on the tape format. The data format for events will be identical
c     for all partitions.
c     
c     The user analysis can be done thru implementing the following 3
c     routines.
c     aninit  : Analysis init
c     anuser  : Analyze one event
c     anfin   : Terminate analysis i.e. output results
c     
c     It should only be necessary for the user to have included
c     the common blocks defined by the include file
c     event.inc       to have access to the necessary event information
c     
c     
      implicit none
c     
c     include files
#include "event.inc"
#include "evt_ctrl.inc"
#include "evt_zebra.inc"
c     
      integer*4 egzin, egzinit, egzend
      integer*4 read_run_bank_hijet,
     1          read_run_bank_venus,
     2          read_run_bank_fritiof,
     3          eg_read_code_bank
      integer*4 istat, jstat
      integer*4 idummy
      character*1 chfile*80,chstat*80
      logical done
c     
c     Executable code
c     ===============
c     
c     This program is for reading
c     ===========================
      z_input = .true.
c     Call the user initialization
c     
      call aninit
c     

c     Ask for filename
      call kupros('Input ZEBRA/FZ filename',chfile,idummy)

c     init egz package, 0 for not GEANT
      istat = egzinit(0, chfile)
      if (istat.ne.0) goto 9000
c     
c     File loop
c     =========
      jstat = 0
      done = .false.
      do while (.not.done)
c     
c     Event loop
c     ==========
          istat = 0
          do while (istat.eq.0)
              istat = egzin()
              if (istat.eq.0) call anuser (istat)
          end do
c
c     check if additional files
c     
          call kupros('Input additional files? (y/n)',chstat,
     $    idummy)
          if(chstat(1:1).eq.'Y' .or. chstat(1:1).eq.'y')then
c     
              call fzendi(z_lun_in,'T')
              close(z_lun_in)
              call kupros('Input ZEBRA/FZ filename',chfile,idummy)

c     
c     open file using exchange format
c     
              open (z_lun_in, file=chfile, status='OLD',
     $        form='UNFORMATTED', access='DIRECT',
     $        recl=z_blocksize)
              
c
c     set ZEBRA file
c
              call fzfile(z_lun_in,z_blocksize/4,z_fz_option_in)
              call fzlogl(z_lun_in, z_fzlogl) ! log = 2 lots of messages
c     
              istat = eg_read_code_bank(6)
              if (istat.ne.0) then ! error
                  goto 9000
              end if
c     
c     Generator dependant RUN bank
c     
              if (evtcod.eq.hijet) then
                  istat = read_run_bank_hijet(6)
              else if (evtcod.eq.venus) then
                  istat = read_run_bank_venus(6)
              else if (evtcod.eq.frtiof) then
                  istat = read_run_bank_fritiof(6)
              else if (evtcod.eq.hijng) then
                  istat = read_run_bank_fritiof(6)
              else
                  write(6,*)'  <E> Unknown event generator:',
     $            evtcod
                  goto 9000
              end if
c     
              call def_event_form
          else
              done  =  .true.
              jstat = 1
          endif
      enddo
c     
      call anfin

c     0 for not GEANT
      istat = egzend (0)
c     
 9000 continue
c     call exit

      end
