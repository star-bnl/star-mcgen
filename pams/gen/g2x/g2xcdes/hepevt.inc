*
* $Id: hepevt.inc,v 1.1 1998/03/25 16:49:14 longacre Exp $
*

*
* hepevt.inc
*
      INTEGER    NMXHEP
      PARAMETER (NMXHEP=30000)
      INTEGER       NEVHEP,NHEP,ISTHEP,        IDHEP,
     &JMOHEP          ,JDAHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP,IDHEP,
     &JMOHEP(2),JDAHEP(2),PHEP(5),VHEP(4)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/