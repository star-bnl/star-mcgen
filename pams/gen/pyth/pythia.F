* $Id: pythia.F,v 1.1 1998/02/10 15:13:46 fisyak Exp $
* $Log: pythia.F,v $
* Revision 1.1  1998/02/10 15:13:46  fisyak
* Add interface with Staf
*
      INTEGER*4 FUNCTION PYTHIA(par_h, par) 
      IMPLICIT NONE
*:>--------------------------------------------------------------------
*: ROUTINE:    PYTHIA
*: generate one pythia event and put into table
*:>--------------------------------------------------------------------
      EXTERNAL PYDATA, PYEVNT, PYHEPC
#include "pythia.inc"
#include "pycdes/hepevt.inc"
      INTEGER I, J
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
