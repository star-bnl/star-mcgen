* $Id: pythia.F,v 1.3 1999/01/11 22:35:57 nevski Exp $
* $Log: pythia.F,v $
* Revision 1.3  1999/01/11 22:35:57  nevski
* force loading of rluset, pystat
*
* Revision 1.2  1998/03/27 14:16:16  fisyak
* force load of py_init
*
      INTEGER*4 FUNCTION PYTHIA(par_h, par) 
      IMPLICIT NONE
*:>--------------------------------------------------------------------
*: ROUTINE:    PYTHIA
*: generate one pythia event and put into table
*:>--------------------------------------------------------------------
      EXTERNAL PYDATA, PYEVNT, PYHEPC, PY_INIT
#include "pythia.inc"
#include "pycdes/hepevt.inc"
      INTEGER I, J, LOAD
      DATA LOAD /-1/
      SAVE LOAD
      IF (LOAD .EQ. 0) then
*         just make it loaded
          CALL PY_INIT(4)
          CALL PYSTAT(0)
          CALL RLUSET(0)
      endif
*---     
      CALL PYEVNT
      CALL PYHEPC (1)
*
*-- Sanity check of header for table "par"...
*
      par_h.nok = 0
      IF(NHEP .GT. par_h.maxlen) THEN
         PYTHIA = STAFCV_BAD
         GO TO 999
      ENDIF
      DO I = 1, NHEP
        par(i).isthep     = isthep(i)   !  status code
        par(i).idhep      = idhep(i)    !  particle identity
        par(i).jmohep(1)  = jmohep(1,i) !  pointer(s) to position 
        par(i).jmohep(2)  = jmohep(2,i) !  where the mothers(s) stored   
        par(i).jdahep(1)  = jdahep(1,i) !  pointers to position of 
        par(i).jdahep(2)  = jdahep(2,i) !  the first/last daughter  
        DO J = 1, 4
          par(i).phep(j)  = phep(j,i)   !  p4 
          par(i).vhep(j)  = vhep(j,i)   !  production vertex (mm) and time (mm/c) 
        END DO
        par(i).phep(5)    = phep(5,i)   !  mass (GeV)  
      END DO
      par_h.nok = nhep
*
      PYTHIA=STAFCV_OK
  999 END
