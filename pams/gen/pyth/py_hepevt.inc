*
* $Id: py_hepevt.inc,v 1.1 1999/08/13 23:45:11 fisyak Exp $
*
* $Log: py_hepevt.inc,v $
* Revision 1.1  1999/08/13 23:45:11  fisyak
* introduce py_hepevt
*
* Revision 1.2  1998/02/10 15:14:38  fisyak
* Define variables
*
* Revision 1.1  1998/02/05 22:56:30  fisyak
* Import from d0gstar
*
*Revision 1.1  1997/09/02  17:35:38  fisyak
*Add pythia version 6.103
*
*
*
* hepevt.inc
*
      INTEGER    NMXHEP
      PARAMETER (NMXHEP=4000)
      INTEGER       NEVHEP,NHEP,ISTHEP,        IDHEP,
     &JMOHEP          ,JDAHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/
