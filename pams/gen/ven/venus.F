* Id: venus.F, v 1.1 1998/04/02 11:10:40 longacre Exp $
*$Log: venus.F,v $
*Revision 1.8  1999/10/09 01:22:03  nevski
*pgf->g77 cleanup
*
*Revision 1.7  1999/08/14 21:49:57  fisyak
*eliminate include and source level
*
*Revision 1.6  1999/08/13 23:35:26  fisyak
*unify hepevt
*
*Revision 1.5  1999/06/25 16:18:14  longacre
*add flag icmtst keep beamgas in lab fram
*
*Revision 1.4  1998/08/02 21:38:16  longacre
*added vertex cm -> mm
*
*Revision 1.3  1998/07/16 13:49:32  longacre
*fixed bugs with ks and kl ids
*
*Revision 1.2  1998/04/20 21:01:36  longacre
*getting ZEBRA out
*
*Revision 1.1  1998/04/16 17:26:44  longacre
*checking in new venus staf module
*
c     -*-Fortran-*-
c-----------------------------------------------------------------------
      INTEGER FUNCTION VENUS(par_h, par)
      implicit none
c     writes the results of a simulation into particle table.
c-----------------------------------------------------------------------
c     
c     Common and variables relating to VENUS
c     ======================================
#include "venus.inc"
#include "venpack.inc"
#include "hepevt.inc"
#include "headpss.inc"
#include "outfgdt.inc"
c     
      real scfermi
C---------------------------------------------------------------------
      INTEGER isahep
      INTEGER II
      INTEGER I
      INTEGER J
      INTEGER nevout
      INTEGER idtype
      INTEGER nptls
      Real*8   pz, x, u, y, mperp
c     
c
c MAB

      data scfermi/0.100E-12/
      IF(IFIRST.EQ.0) THEN
CCCCC STRATOS
        call afiles
        call ainitl(1)
      if(ifgdt.eq.1) OPEN(UNIT=32, FILE='venus.txt', status='unknown')
      IFIRST=1
      ENDIF
      CALL DATIME(I,J)
      VSSHEP(2) = VSSHEP(2) + 1
      nevout = VSSHEP(2)
      VSSHEP(3) = FLOAT(I)     
      VSSHEP(4) = FLOAT(J)
      par(1).isthep = 10
      par(1).idhep = 9999999
      par(1).jmohep(1) = 0
      par(1).jmohep(2) = 0
      par(1).jdahep(1) = 0
      par(1).jdahep(2) = 0
      par(1).phep(1) = psshep(1)
      par(1).phep(2) = psshep(2)
      par(1).phep(3) = psshep(3)
      par(1).phep(4) = psshep(4)
      par(1).phep(5) = psshep(5)
      par(1).vhep(1) = vsshep(1)
      par(1).vhep(1) = vsshep(1)
      par(1).vhep(2) = vsshep(2)
      par(1).vhep(3) = vsshep(3)
      par(1).vhep(4) = vsshep(4)
      call avenus
      par(1).phep(1) = bimevt
      par(1).phep(2) = phievt
      nptls=0
      do 10 i=1,nptl
          nrptla(i)=-999
          if(istptl(i).gt.istmax)goto10
          nptls=nptls+1
          nrptla(i)=nptls
 10       continue
         Do i = 1, nptl

            pz = pptl(3,i)
            mperp = sqrt ( pptl(1,i)**2 + pptl(2,i)**2 + pptl(5,i)**2 )

c Use relation pz = mperp * sinh(y) to get rapidity of particle.
c Note the use of a log expression to get asinh.

            If ( mperp .gt. 0. ) Then
               x = pz/mperp
               u = sqrt( x**2 + 1. )
               y = log ( x + u )    ! asinh
            Else
               Write (*, *) 'astore_egz: strange mperp value of ', 
     >              mperp
               y = 999.
            EndIf

            if(icmtst.lt.0) then
            pptl(3,i) = mperp*sinh(y-yhaha)
            pptl(4,i) = mperp*cosh(y-yhaha)
            endif

         EndDo
CCCC STRATOS
CCCC Try to write text output.
      if(ifgdt.eq.1) then
      write(32,98) nptls, nevout, bimevt, phievt
98    FORMAT(1X,I6,7X,I6,8X,' BTOT  =',G12.5,' PHI =',G12.5)
      Do 9918 i=1, nptl
        if(nrptla(i).lt.0) goto 9918
        Call isa_to_geant(idptl(i),idtype)
        write(32, 9911) idtype, pptl(1,i), pptl(2,i), 
     >  pptl(3,i), pptl(4,i), xorptl(1,i), xorptl(2,i),
     >  xorptl(3,i), xorptl(4,i)
 9911   FORMAT(I9,9X,8F12.5)
 9918 Continue
      endif
      II = 1
      Do 9910 i=1, nptl
        if(nrptla(i).lt.0) goto 9910
        II = II + 1
        par(ii).isthep     = 1           !  status code
        par(ii).idhep      = isahep(idptl(i)) !  particle identity
        par(ii).jmohep(1)  = 0           !  pointer(s) to position 
        par(ii).jmohep(2)  = 0           !  where the mothers(s) stored   
        par(ii).jdahep(1)  = 0           !  pointers to position of 
        par(ii).jdahep(2)  = 0           !  the first/last daughter  
        DO J = 1, 4
          par(ii).phep(j)  = pptl(j,i)   !  p4 
          par(ii).vhep(j)  = xorptl(j,i)*scfermi*10.! vertex (mm) and time (mm)
        END DO
        par(ii).phep(5)    = pptl(5,i)   !  mass (GeV)  
9910    CONTINUE
      par_h.nok = II
      VENUS=STAFCV_OK
C---------------------------------------------------------------------
c     lptl ......... record label (lptl=3)
c     nrevt ........ event number
c     nrptl ........ particle number
c     i ............ original ptl number
c     idptl ........ particle id
c     pptl ......... 5-momentum (px,py,pz,en,mass) in lab
c     ioptl ........ origin (-999:parent not stored, -1,0:no parent)
c     joptl ........ origin (second parent)
c     istptl ....... stable (=0) or not (=1)
c     xorptl ....... space-time point (x,y,z,t) on ptl track (pp-cm)
c     tivptl ....... time interval of existence
c     nqjptl ....... quark numbers of jets
      return
      end
