         FUNCTION PAR_ANAL(par_h, par) 
         IMPLICIT NONE
*
*********************************************************************
*                                                                   *
*         This module gives an example how to create ntuples        *
*         or histograms from particle table                         *
*                                                                   *
*********************************************************************
*
#include "table_header.inc"
#include "particle.inc"
         RECORD/table_head_st/ par_h
         RECORD/particle_st/   par(*)
C
         INTEGER PAR_ANAL
         INTEGER NATT, IP, I
C
         DIMENSION PIEV(2)
         DIMENSION PREV(2)
         DIMENSION NMEV(2)
         REAL PIEV, PREV, NMEV
         REAL ist, idh, jmo, jda1, jda2, px, py, pz, pe, pm, vt
         REAL vx, vy, vz, pt, theta, eta, phi
C
          ist(i)  = par(i).isthep
          idh(i)  = par(i).idhep
          jmo(i)  = par(i).jmohep(1)
          jda1(i) = par(i).jdahep(1)
          jda2(i) = par(i).jdahep(2)
          px(i)   = par(i).phep(1)
          py(i)   = par(i).phep(2)
          pz(i)   = par(i).phep(3)
          pe(i)   = par(i).phep(4)
          pm(i)   = par(i).phep(5)
          vx(i)   = par(i).vhep(1)
          vy(i)   = par(i).vhep(2)
          vz(i)   = par(i).vhep(3)
          vt(i)   = par(i).vhep(4) 
C
         NATT = par_h.nok
C
         nmev(1) = 0.
         nmev(2) = 0.
C       
         DO IP = 1 ,NATT
           if(ist(ip).ne.1) goto 99
C
             pt = sqrt(px(ip)**2+py(ip)**2)
             theta = atan2 ( pt, pz(ip) )
             eta   = -log(tan(theta/2.) )                   
             phi  = atan2 ( py(ip), px(ip) )
C                 
             if(abs(idh(ip)) .eq. 211) then
                  nmev(1) = nmev(1) + 1
                  piev(1) = pt
                  piev(2) = eta
                  call hfn(101, piev)
             elseif(abs(idh(ip)) .eq. 2212) then
                  nmev(2) = nmev(2) + 1
                  prev(1) = pt
                  prev(2) = eta
                  if(abs(eta).lt.20.)then
                  call hfn(102, prev)
                  endif
             endif
C  
 99       END DO
                 call hfn(103, nmev)
C
         PAR_ANAL=1
C
        END
******************************************************************************
        SUBROUTINE PAR_START
*******************************************************************
*                                                                 *
*        This subroutine initialize ntuples                       *
*                                                                 *
*        The folowing text is an example                          *
*                                                                 *
*******************************************************************
*
        CHARACTER*6 HTAGP(2)/'PTPI','ETHAPI'/
        CHARACTER*6 HTAGE(2)/'PTPR','ETHAPR'/
        CHARACTER*6 ITAGN(2)/'NMPI','NMPR'/
        REAL PIEV(2),PREV(2),NMEV(2)

        EQUIVALENCE (PIEV(1), PTPI), (PIEV(2), ETHAPI)
        EQUIVALENCE (PREV(1), PTPR), (PREV(2), ETHAPR)
        EQUIVALENCE (NMEV(1), NMPI), (NMEV(2), NMPR)
        REAL PTPR, PTPI, ETHAPR, ETHAPI, NMPR, NMPI
C
        CALL HBOOKN(101,'PARTICLE PT',2,' ',1000, htagp)
        CALL HBOOKN(102,'PARTICLE ETHA',2,' ',1000, htage)
        CALL HBOOKN(103,'PARTICLE MULTIPLICITY',2,' ',1000,itagn)
C
        RETURN
        END
