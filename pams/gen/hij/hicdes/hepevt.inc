*
* $Id: hepevt.inc,v 1.2 1998/05/15 20:03:54 didenko Exp $
*
* $Log: hepevt.inc,v $
* Revision 1.2  1998/05/15 20:03:54  didenko
* updating hepevt
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
      PARAMETER (NMXHEP=9000)
      INTEGER       NEVHEP,NHEP,ISTHEP,        IDHEP,
     &JMOHEP          ,JDAHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      REAL PHEP,VHEP
      SAVE /HEPEVT/
