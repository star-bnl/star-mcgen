      SUBROUTINE PYINIT(FRAME,BEAM,TARGET,WIN)  
    
C...Initializes the generation procedure; finds maxima of the   
C...differential cross-sections to be used for weighting.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/LUDAT4/CHAF(500)   
      CHARACTER CHAF*8  
      SAVE /LUDAT4/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      CHARACTER*(*) FRAME,BEAM,TARGET   
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHMO(12)*3,CHLH(2)*6 
      DATA CHMO/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',  
     &'Oct','Nov','Dec'/, CHLH/'lepton','hadron'/   
    
C...Write headers.  
C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1000) MSTP(181),MSTP(182),  
C     &MSTP(185),CHMO(MSTP(184)),MSTP(183)   
      CALL LULIST(0)
C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1100)  
    
C...Identify beam and target particles and initialize kinematics.   
      CHFRAM=FRAME//' ' 
      CHBEAM=BEAM//' '  
      CHTARG=TARGET//' '    
      CALL PYINKI(CHFRAM,CHBEAM,CHTARG,WIN) 
    
C...Select partonic subprocesses to be included in the simulation.  
      IF(MSEL.NE.0) THEN    
        DO 100 I=1,200  
  100   MSUB(I)=0   
      ENDIF 
      IF(MINT(43).EQ.1.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN   
C...Lepton+lepton -> gamma/Z0 or W. 
        IF(MINT(11)+MINT(12).EQ.0) MSUB(1)=1    
        IF(MINT(11)+MINT(12).NE.0) MSUB(2)=1    
      ELSEIF(MSEL.EQ.1) THEN    
C...High-pT QCD processes:  
        MSUB(11)=1  
        MSUB(12)=1  
        MSUB(13)=1  
        MSUB(28)=1  
        MSUB(53)=1  
        MSUB(68)=1  
        IF(MSTP(82).LE.1.AND.CKIN(3).LT.PARP(81)) MSUB(95)=1    
        IF(MSTP(82).GE.2.AND.CKIN(3).LT.PARP(82)) MSUB(95)=1    
      ELSEIF(MSEL.EQ.2) THEN    
C...All QCD processes:  
        MSUB(11)=1  
        MSUB(12)=1  
        MSUB(13)=1  
        MSUB(28)=1  
        MSUB(53)=1  
        MSUB(68)=1  
        MSUB(91)=1  
        MSUB(92)=1  
        MSUB(93)=1  
        MSUB(95)=1  
      ELSEIF(MSEL.GE.4.AND.MSEL.LE.8) THEN  
C...Heavy quark production. 
        MSUB(81)=1  
        MSUB(82)=1  
        DO 110 J=1,MIN(8,MDCY(21,3))    
  110   MDME(MDCY(21,2)+J-1,1)=0    
        MDME(MDCY(21,2)+MSEL-1,1)=1 
      ELSEIF(MSEL.EQ.10) THEN   
C...Prompt photon production:   
        MSUB(14)=1  
        MSUB(18)=1  
        MSUB(29)=1  
      ELSEIF(MSEL.EQ.11) THEN   
C...Z0/gamma* production:   
        MSUB(1)=1   
      ELSEIF(MSEL.EQ.12) THEN   
C...W+/- production:    
        MSUB(2)=1   
      ELSEIF(MSEL.EQ.13) THEN   
C...Z0 + jet:   
        MSUB(15)=1  
        MSUB(30)=1  
      ELSEIF(MSEL.EQ.14) THEN   
C...W+/- + jet: 
        MSUB(16)=1  
        MSUB(31)=1  
      ELSEIF(MSEL.EQ.15) THEN   
C...Z0 & W+/- pair production:  
        MSUB(19)=1  
        MSUB(20)=1  
        MSUB(22)=1  
        MSUB(23)=1  
        MSUB(25)=1  
      ELSEIF(MSEL.EQ.16) THEN   
C...H0 production:  
        MSUB(3)=1   
        MSUB(5)=1   
        MSUB(8)=1   
        MSUB(102)=1 
      ELSEIF(MSEL.EQ.17) THEN   
C...H0 & Z0 or W+/- pair production:    
        MSUB(24)=1  
        MSUB(26)=1  
      ELSEIF(MSEL.EQ.21) THEN   
C...Z'0 production: 
        MSUB(141)=1 
      ELSEIF(MSEL.EQ.22) THEN   
C...H+/- production:    
        MSUB(142)=1 
      ELSEIF(MSEL.EQ.23) THEN   
C...R production:   
        MSUB(143)=1 
      ENDIF 
    
C...Count number of subprocesses on.    
      MINT(44)=0    
      DO 120 ISUB=1,200 
      IF(MINT(43).LT.4.AND.ISUB.GE.91.AND.ISUB.LE.96.AND.   
     &MSUB(ISUB).EQ.1) THEN 
        WRITE(MSTU(11),1200) ISUB,CHLH(MINT(41)),CHLH(MINT(42)) 
        STOP    
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).EQ.-1) THEN 
        WRITE(MSTU(11),1300) ISUB   
        STOP    
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).LE.-2) THEN 
        WRITE(MSTU(11),1400) ISUB   
        STOP    
      ELSEIF(MSUB(ISUB).EQ.1) THEN  
        MINT(44)=MINT(44)+1 
      ENDIF 
  120 CONTINUE  
      IF(MINT(44).EQ.0) THEN    
        WRITE(MSTU(11),1500)    
        STOP    
      ENDIF 
      MINT(45)=MINT(44)-MSUB(91)-MSUB(92)-MSUB(93)-MSUB(94) 
    
C...Maximum 4 generations; set maximum number of allowed flavours.  
      MSTP(1)=MIN(4,MSTP(1))    
      MSTU(114)=MIN(MSTU(114),2*MSTP(1))    
      MSTP(54)=MIN(MSTP(54),2*MSTP(1))  
    
C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton. 
      DO 140 I=-20,20   
      VINT(180+I)=0.    
      IA=IABS(I)    
      IF(IA.GE.1.AND.IA.LE.2*MSTP(1)) THEN  
        DO 130 J=1,MSTP(1)  
        IB=2*J-1+MOD(IA,2)  
        IPM=(5-ISIGN(1,I))/2    
        IDC=J+MDCY(IA,2)+2  
  130   IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) VINT(180+I)= 
     &  VINT(180+I)+VCKM((IA+1)/2,(IB+1)/2) 
      ELSEIF(IA.GE.11.AND.IA.LE.10+2*MSTP(1)) THEN  
        VINT(180+I)=1.  
      ENDIF 
  140 CONTINUE  
    
C...Choose Lambda value to use in alpha-strong. 
      MSTU(111)=MSTP(2) 
      IF(MSTP(3).GE.1) THEN 
        ALAM=PARP(1)    
        IF(MSTP(51).EQ.1) ALAM=0.2  
        IF(MSTP(51).EQ.2) ALAM=0.29 
        IF(MSTP(51).EQ.3) ALAM=0.2  
        IF(MSTP(51).EQ.4) ALAM=0.4  
        IF(MSTP(51).EQ.11) ALAM=0.16    
        IF(MSTP(51).EQ.12) ALAM=0.26    
        IF(MSTP(51).EQ.13) ALAM=0.36    
        PARP(1)=ALAM    
        PARP(61)=ALAM   
        PARU(112)=ALAM  
        PARJ(81)=ALAM   
      ENDIF 
    
C...Initialize widths and partial widths for resonances.    
      CALL PYINRE   
    
C...Reset variables for cross-section calculation.  
      DO 150 I=0,200    
      DO 150 J=1,3  
      NGEN(I,J)=0   
  150 XSEC(I,J)=0.  
      VINT(108)=0.  
    
C...Find parametrized total cross-sections. 
      IF(MINT(43).EQ.4) CALL PYXTOT 
    
C...Maxima of differential cross-sections.  
      IF(MSTP(121).LE.0) CALL PYMAXI    
    
C...Initialize possibility of overlayed events. 
      IF(MSTP(131).NE.0) CALL PYOVLY(1) 
    
C...Initialize multiple interactions with variable impact parameter.    
      IF(MINT(43).EQ.4.AND.(MINT(45).NE.0.OR.MSTP(131).NE.0).AND.   
     &MSTP(82).GE.2) CALL PYMULT(1) 
C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1600)  
    
C...Formats for initialization information. 
 1000 FORMAT(///20X,'The Lund Monte Carlo - PYTHIA version ',I1,'.',I1/ 
     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)  
 1100 FORMAT('1',18('*'),1X,'PYINIT: initialization of PYTHIA ',    
     &'routines',1X,17('*'))    
 1200 FORMAT(1X,'Error: process number ',I3,' not meaningful for ',A6,  
     &'-',A6,' interactions.'/1X,'Execution stopped!')  
 1300 FORMAT(1X,'Error: requested subprocess',I4,' not implemented.'/   
     &1X,'Execution stopped!')  
 1400 FORMAT(1X,'Error: requested subprocess',I4,' not existing.'/  
     &1X,'Execution stopped!')  
 1500 FORMAT(1X,'Error: no subprocess switched on.'/    
     &1X,'Execution stopped.')  
 1600 FORMAT(/1X,22('*'),1X,'PYINIT: initialization completed',1X,  
     &22('*'))  
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYTHIA 
    
C...Administers the generation of a high-pt event via calls to a number 
C...of subroutines; also computes cross-sections.   
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
    
C...Loop over desired number of overlayed events (normally 1).  
      MINT(7)=0 
      MINT(8)=0 
      NOVL=1    
      IF(MSTP(131).NE.0) CALL PYOVLY(2) 
      IF(MSTP(131).NE.0) NOVL=MINT(81)  
      MINT(83)=0    
      MINT(84)=MSTP(126)    
      MSTU(70)=0    
      DO 190 IOVL=1,NOVL    
      IF(MINT(84)+100.GE.MSTU(4)) THEN  
        CALL LUERRM(11, 
     &  '(PYTHIA:) no more space in LUJETS for overlayed events')   
        IF(MSTU(21).GE.1) GOTO 200  
      ENDIF 
      MINT(82)=IOVL 
    
C...Generate variables of hard scattering.  
  100 CONTINUE  
      IF(IOVL.EQ.1) NGEN(0,2)=NGEN(0,2)+1   
      MINT(31)=0    
      MINT(51)=0    
      CALL PYRAND   
      ISUB=MINT(1)  
      IF(IOVL.EQ.1) THEN    
        NGEN(ISUB,2)=NGEN(ISUB,2)+1 
    
C...Store information on hard interaction.  
        DO 110 J=1,200  
        MSTI(J)=0   
  110   PARI(J)=0.  
        MSTI(1)=MINT(1) 
        MSTI(2)=MINT(2) 
        MSTI(11)=MINT(11)   
        MSTI(12)=MINT(12)   
        MSTI(15)=MINT(15)   
        MSTI(16)=MINT(16)   
        MSTI(17)=MINT(17)   
        MSTI(18)=MINT(18)   
        PARI(11)=VINT(1)    
        PARI(12)=VINT(2)    
        IF(ISUB.NE.95) THEN 
          DO 120 J=13,22    
  120     PARI(J)=VINT(30+J)    
          PARI(33)=VINT(41) 
          PARI(34)=VINT(42) 
          PARI(35)=PARI(33)-PARI(34)    
          PARI(36)=VINT(21) 
          PARI(37)=VINT(22) 
          PARI(38)=VINT(26) 
          PARI(41)=VINT(23) 
        ENDIF   
      ENDIF 
    
      IF(MSTP(111).EQ.-1) GOTO 160  
      IF(ISUB.LE.90.OR.ISUB.GE.95) THEN 
C...Hard scattering (including low-pT): 
C...reconstruct kinematics and colour flow of hard scattering.  
        CALL PYSCAT 
        IF(MINT(51).EQ.1) GOTO 100  
    
C...Showering of initial state partons (optional).  
        IPU1=MINT(84)+1 
        IPU2=MINT(84)+2 
        IF(MSTP(61).GE.1.AND.MINT(43).NE.1.AND.ISUB.NE.95)  
     &  CALL PYSSPA(IPU1,IPU2)  
        NSAV1=N 
    
C...Multiple interactions.  
        IF(MSTP(81).GE.1.AND.MINT(43).EQ.4.AND.ISUB.NE.95)  
     &  CALL PYMULT(6)  
        MINT(1)=ISUB    
        NSAV2=N 
    
C...Hadron remnants and primordial kT.  
        CALL PYREMN(IPU1,IPU2)  
        IF(MINT(51).EQ.1) GOTO 100  
        NSAV3=N 
    
C...Showering of final state partons (optional).    
        IPU3=MINT(84)+3 
        IPU4=MINT(84)+4 
        IF(MSTP(71).GE.1.AND.ISUB.NE.95.AND.K(IPU3,1).GT.0.AND. 
     &  K(IPU3,1).LE.10.AND.K(IPU4,1).GT.0.AND.K(IPU4,1).LE.10) THEN    
          QMAX=SQRT(PARP(71)*VINT(52))  
          IF(ISUB.EQ.5) QMAX=SQRT(PMAS(23,1)**2)    
          IF(ISUB.EQ.8) QMAX=SQRT(PMAS(24,1)**2)    
          CALL LUSHOW(IPU3,IPU4,QMAX)   
        ENDIF   
    
C...Sum up transverse and longitudinal momenta. 
        IF(IOVL.EQ.1) THEN  
          PARI(65)=2.*PARI(17)  
          DO 130 I=MSTP(126)+1,N    
          IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 130  
          PT=SQRT(P(I,1)**2+P(I,2)**2)  
          PARI(69)=PARI(69)+PT  
          IF(I.LE.NSAV1.OR.I.GT.NSAV3) PARI(66)=PARI(66)+PT 
          IF(I.GT.NSAV1.AND.I.LE.NSAV2) PARI(68)=PARI(68)+PT    
  130     CONTINUE  
          PARI(67)=PARI(68) 
          PARI(71)=VINT(151)    
          PARI(72)=VINT(152)    
          PARI(73)=VINT(151)    
          PARI(74)=VINT(152)    
        ENDIF   
    
C...Decay of final state resonances.    
        IF(MSTP(41).GE.1.AND.ISUB.NE.95) CALL PYRESD    
    
      ELSE  
C...Diffractive and elastic scattering. 
        CALL PYDIFF 
        IF(IOVL.EQ.1) THEN  
          PARI(65)=2.*PARI(17)  
          PARI(66)=PARI(65) 
          PARI(69)=PARI(65) 
        ENDIF   
      ENDIF 
    
C...Recalculate energies from momenta and masses (if desired).  
      IF(MSTP(113).GE.1) THEN   
        DO 140 I=MINT(83)+1,N   
  140   IF(K(I,1).GT.0.AND.K(I,1).LE.10) P(I,4)=SQRT(P(I,1)**2+ 
     &  P(I,2)**2+P(I,3)**2+P(I,5)**2)  
      ENDIF 
    
C...Rearrange partons along strings, check invariant mass cuts. 
      MSTU(28)=0    
      CALL LUPREP(MINT(84)+1)   
      IF(MSTP(112).EQ.1.AND.MSTU(28).EQ.3) GOTO 100 
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) THEN 
        DO 150 I=MINT(84)+1,N   
        IF(K(I,2).NE.94) GOTO 150   
        K(I+1,3)=MOD(K(I+1,4)/MSTU(5),MSTU(5))  
        K(I+2,3)=MOD(K(I+2,4)/MSTU(5),MSTU(5))  
  150   CONTINUE    
        CALL LUEDIT(12) 
        CALL LUEDIT(14) 
        IF(MSTP(125).EQ.0) CALL LUEDIT(15)  
        IF(MSTP(125).EQ.0) MINT(4)=0    
      ENDIF 
    
C...Introduce separators between sections in LULIST event listing.  
      IF(IOVL.EQ.1.AND.MSTP(125).LE.0) THEN 
        MSTU(70)=1  
        MSTU(71)=N  
      ELSEIF(IOVL.EQ.1) THEN    
        MSTU(70)=3  
        MSTU(71)=2  
        MSTU(72)=MINT(4)    
        MSTU(73)=N  
      ENDIF 
    
C...Perform hadronization (if desired). 
      IF(MSTP(111).GE.1) CALL LUEXEC    
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) CALL LUEDIT(14)  
    
C...Calculate Monte Carlo estimates of cross-sections.  
  160 IF(IOVL.EQ.1) THEN    
        IF(MSTP(111).NE.-1) NGEN(ISUB,3)=NGEN(ISUB,3)+1 
        NGEN(0,3)=NGEN(0,3)+1   
        XSEC(0,3)=0.    
        DO 170 I=1,200  
        IF(I.EQ.96) THEN    
          XSEC(I,3)=0.  
        ELSEIF(MSUB(95).EQ.1.AND.(I.EQ.11.OR.I.EQ.12.OR.I.EQ.13.OR. 
     &  I.EQ.28.OR.I.EQ.53.OR.I.EQ.68)) THEN    
          XSEC(I,3)=XSEC(96,2)*NGEN(I,3)/MAX(1.,FLOAT(NGEN(96,1))*  
     &    FLOAT(NGEN(96,2)))    
        ELSEIF(NGEN(I,1).EQ.0) THEN 
          XSEC(I,3)=0.  
        ELSEIF(NGEN(I,2).EQ.0) THEN 
          XSEC(I,3)=XSEC(I,2)*NGEN(0,3)/(FLOAT(NGEN(I,1))*  
     &    FLOAT(NGEN(0,2))) 
        ELSE    
          XSEC(I,3)=XSEC(I,2)*NGEN(I,3)/(FLOAT(NGEN(I,1))*  
     &    FLOAT(NGEN(I,2))) 
        ENDIF   
  170   XSEC(0,3)=XSEC(0,3)+XSEC(I,3)   
        IF(MSUB(95).EQ.1) THEN  
          NGENS=NGEN(91,3)+NGEN(92,3)+NGEN(93,3)+NGEN(94,3)+NGEN(95,3)  
          XSECS=XSEC(91,3)+XSEC(92,3)+XSEC(93,3)+XSEC(94,3)+XSEC(95,3)  
          XMAXS=XSEC(95,1)  
          IF(MSUB(91).EQ.1) XMAXS=XMAXS+XSEC(91,1)  
          IF(MSUB(92).EQ.1) XMAXS=XMAXS+XSEC(92,1)  
          IF(MSUB(93).EQ.1) XMAXS=XMAXS+XSEC(93,1)  
          IF(MSUB(94).EQ.1) XMAXS=XMAXS+XSEC(94,1)  
          FAC=1.    
          IF(NGENS.LT.NGEN(0,3)) FAC=(XMAXS-XSECS)/(XSEC(0,3)-XSECS)    
          XSEC(11,3)=FAC*XSEC(11,3) 
          XSEC(12,3)=FAC*XSEC(12,3) 
          XSEC(13,3)=FAC*XSEC(13,3) 
          XSEC(28,3)=FAC*XSEC(28,3) 
          XSEC(53,3)=FAC*XSEC(53,3) 
          XSEC(68,3)=FAC*XSEC(68,3) 
          XSEC(0,3)=XSEC(91,3)+XSEC(92,3)+XSEC(93,3)+XSEC(94,3)+    
     &    XSEC(95,1)    
        ENDIF   
    
C...Store final information.    
        MINT(5)=MINT(5)+1   
        MSTI(3)=MINT(3) 
        MSTI(4)=MINT(4) 
        MSTI(5)=MINT(5) 
        MSTI(6)=MINT(6) 
        MSTI(7)=MINT(7) 
        MSTI(8)=MINT(8) 
        MSTI(13)=MINT(13)   
        MSTI(14)=MINT(14)   
        MSTI(21)=MINT(21)   
        MSTI(22)=MINT(22)   
        MSTI(23)=MINT(23)   
        MSTI(24)=MINT(24)   
        MSTI(25)=MINT(25)   
        MSTI(26)=MINT(26)   
        MSTI(31)=MINT(31)   
        PARI(1)=XSEC(0,3)   
        PARI(2)=XSEC(0,3)/MINT(5)   
        PARI(31)=VINT(141)  
        PARI(32)=VINT(142)  
        IF(ISUB.NE.95.AND.MINT(7)*MINT(8).NE.0) THEN    
          PARI(42)=2.*VINT(47)/VINT(1)  
          DO 180 IS=7,8 
          PARI(36+IS)=P(MINT(IS),3)/VINT(1) 
          PARI(38+IS)=P(MINT(IS),4)/VINT(1) 
          I=MINT(IS)    
          PR=MAX(1E-20,P(I,5)**2+P(I,1)**2+P(I,2)**2)   
          PARI(40+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/    
     &    SQRT(PR),1E20)),P(I,3))   
          PR=MAX(1E-20,P(I,1)**2+P(I,2)**2) 
          PARI(42+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/    
     &    SQRT(PR),1E20)),P(I,3))   
          PARI(44+IS)=P(I,3)/SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)    
          PARI(46+IS)=ULANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))  
          PARI(48+IS)=ULANGL(P(I,1),P(I,2)) 
  180     CONTINUE  
        ENDIF   
        PARI(61)=VINT(148)  
        IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN 
          MSTU(161)=MINT(21)    
          MSTU(162)=0   
        ELSE    
          MSTU(161)=MINT(21)    
          MSTU(162)=MINT(22)    
        ENDIF   
      ENDIF 
    
C...Prepare to go to next overlayed event.  
      MSTI(41)=IOVL 
      IF(IOVL.GE.2.AND.IOVL.LE.10) MSTI(40+IOVL)=ISUB   
      IF(MSTU(70).LT.10) THEN   
        MSTU(70)=MSTU(70)+1 
        MSTU(70+MSTU(70))=N 
      ENDIF 
      MINT(83)=N    
      MINT(84)=N+MSTP(126)  
  190 CONTINUE  
    
C...Information on overlayed events.    
      IF(MSTP(131).EQ.1.AND.MSTP(133).GE.1) THEN    
        PARI(91)=VINT(132)  
        PARI(92)=VINT(133)  
        PARI(93)=VINT(134)  
        IF(MSTP(133).EQ.2) PARI(93)=PARI(93)*XSEC(0,3)/VINT(131)    
      ENDIF 
    
C...Transform to the desired coordinate frame.  
  200 CALL PYFRAM(MSTP(124))    
    
      RETURN    
      END   
    
C***********************************************************************    
    
      SUBROUTINE PYSTAT(MSTAT)  
    
C...Prints out information about cross-sections, decay widths, branching    
C...ratios, kinematical limits, status codes and parameter values.  
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      COMMON/PYINT6/PROC(0:200) 
      CHARACTER PROC*28 
      SAVE /PYINT6/ 
      CHARACTER CHAU*16,CHPA(-40:40)*12,CHIN(2)*12, 
     &STATE(-1:5)*4,CHKIN(21)*18    
      DATA STATE/'----','off ','on  ','on/+','on/-','on/1','on/2'/, 
     &CHKIN/' m_hard (GeV/c^2) ',' p_T_hard (GeV/c) ',  
     &'m_finite (GeV/c^2)','   y*_subsystem   ','     y*_large     ',   
     &'     y*_small     ','    eta*_large    ','    eta*_small    ',   
     &'cos(theta*)_large ','cos(theta*)_small ','       x_1        ',   
     &'       x_2        ','       x_F        ',' cos(theta_hard)  ',   
     &'m''_hard (GeV/c^2) ','       tau        ','        y*        ',  
     &'cos(theta_hard^-) ','cos(theta_hard^+) ','      x_T^2       ',   
     &'       tau''       '/    
    
C...Cross-sections. 
      IF(MSTAT.LE.1) THEN   
        WRITE(MSTU(11),1000)    
        WRITE(MSTU(11),1100)    
        WRITE(MSTU(11),1200) 0,PROC(0),NGEN(0,3),NGEN(0,1),XSEC(0,3)    
        DO 100 I=1,200  
        IF(MSUB(I).NE.1) GOTO 100   
        WRITE(MSTU(11),1200) I,PROC(I),NGEN(I,3),NGEN(I,1),XSEC(I,3)    
  100   CONTINUE    
        WRITE(MSTU(11),1300) 1.-FLOAT(NGEN(0,3))/   
     &  MAX(1.,FLOAT(NGEN(0,2)))    
    
C...Decay widths and branching ratios.  
      ELSEIF(MSTAT.EQ.2) THEN   
        DO 110 KF=-40,40    
        CALL LUNAME(KF,CHAU)    
  110   CHPA(KF)=CHAU(1:12) 
        WRITE(MSTU(11),1400)    
        WRITE(MSTU(11),1500)    
C...Off-shell branchings.   
        DO 130 I=1,17   
        KC=I    
        IF(I.GE.9) KC=I+2   
        IF(I.EQ.17) KC=21   
        WRITE(MSTU(11),1600) CHPA(KC),0.,0.,STATE(MDCY(KC,1)),0.    
        DO 120 J=1,MDCY(KC,3)   
        IDC=J+MDCY(KC,2)-1  
  120   IF(MDME(IDC,2).EQ.102) WRITE(MSTU(11),1700) CHPA(KFDP(IDC,1)),  
     &  CHPA(KFDP(IDC,2)),0.,0.,STATE(MDME(IDC,1)),0.   
  130   CONTINUE    
C...On-shell decays.    
        DO 150 I=1,6    
        KC=I+22 
        IF(I.EQ.4) KC=32    
        IF(I.EQ.5) KC=37    
        IF(I.EQ.6) KC=40    
        IF(WIDE(KC,0).GT.0.) THEN   
          WRITE(MSTU(11),1600) CHPA(KC),WIDP(KC,0),1.,  
     &    STATE(MDCY(KC,1)),1.  
          DO 140 J=1,MDCY(KC,3) 
          IDC=J+MDCY(KC,2)-1    
  140     WRITE(MSTU(11),1700) CHPA(KFDP(IDC,1)),CHPA(KFDP(IDC,2)), 
     &    WIDP(KC,J),WIDP(KC,J)/WIDP(KC,0),STATE(MDME(IDC,1)),  
     &    WIDE(KC,J)/WIDE(KC,0) 
        ELSE    
          WRITE(MSTU(11),1600) CHPA(KC),WIDP(KC,0),1.,  
     &    STATE(MDCY(KC,1)),0.  
        ENDIF   
  150   CONTINUE    
        WRITE(MSTU(11),1800)    
    
C...Allowed incoming partons/particles at hard interaction. 
      ELSEIF(MSTAT.EQ.3) THEN   
        WRITE(MSTU(11),1900)    
        CALL LUNAME(MINT(11),CHAU)  
        CHIN(1)=CHAU(1:12)  
        CALL LUNAME(MINT(12),CHAU)  
        CHIN(2)=CHAU(1:12)  
        WRITE(MSTU(11),2000) CHIN(1),CHIN(2)    
        DO 160 KF=-40,40    
        CALL LUNAME(KF,CHAU)    
  160   CHPA(KF)=CHAU(1:12) 
        IF(MINT(43).EQ.1) THEN  
          WRITE(MSTU(11),2100) CHPA(MINT(11)),STATE(KFIN(1,MINT(11))),  
     &    CHPA(MINT(12)),STATE(KFIN(2,MINT(12)))    
        ELSEIF(MINT(43).EQ.2) THEN  
          WRITE(MSTU(11),2100) CHPA(MINT(11)),STATE(KFIN(1,MINT(11))),  
     &    CHPA(-MSTP(54)),STATE(KFIN(2,-MSTP(54)))  
          DO 170 I=-MSTP(54)+1,-1   
  170     WRITE(MSTU(11),2200) CHPA(I),STATE(KFIN(2,I)) 
          DO 180 I=1,MSTP(54)   
  180     WRITE(MSTU(11),2200) CHPA(I),STATE(KFIN(2,I)) 
          WRITE(MSTU(11),2200) CHPA(21),STATE(KFIN(2,21))   
        ELSEIF(MINT(43).EQ.3) THEN  
          WRITE(MSTU(11),2100) CHPA(-MSTP(54)),STATE(KFIN(1,-MSTP(54))),    
     &    CHPA(MINT(12)),STATE(KFIN(2,MINT(12)))    
          DO 190 I=-MSTP(54)+1,-1   
  190     WRITE(MSTU(11),2300) CHPA(I),STATE(KFIN(1,I)) 
          DO 200 I=1,MSTP(54)   
  200     WRITE(MSTU(11),2300) CHPA(I),STATE(KFIN(1,I)) 
          WRITE(MSTU(11),2300) CHPA(21),STATE(KFIN(1,21))   
        ELSEIF(MINT(43).EQ.4) THEN  
          DO 210 I=-MSTP(54),-1 
  210     WRITE(MSTU(11),2100) CHPA(I),STATE(KFIN(1,I)),CHPA(I),    
     &    STATE(KFIN(2,I))  
          DO 220 I=1,MSTP(54)   
  220     WRITE(MSTU(11),2100) CHPA(I),STATE(KFIN(1,I)),CHPA(I),    
     &    STATE(KFIN(2,I))  
          WRITE(MSTU(11),2100) CHPA(21),STATE(KFIN(1,21)),CHPA(21), 
     &    STATE(KFIN(2,21)) 
        ENDIF   
        WRITE(MSTU(11),2400)    
    
C...User-defined and derived limits on kinematical variables.   
      ELSEIF(MSTAT.EQ.4) THEN   
        WRITE(MSTU(11),2500)    
        WRITE(MSTU(11),2600)    
        SHRMAX=CKIN(2)  
        IF(SHRMAX.LT.0.) SHRMAX=VINT(1) 
        WRITE(MSTU(11),2700) CKIN(1),CHKIN(1),SHRMAX    
        PTHMIN=MAX(CKIN(3),CKIN(5)) 
        PTHMAX=CKIN(4)  
        IF(PTHMAX.LT.0.) PTHMAX=0.5*SHRMAX  
        WRITE(MSTU(11),2800) CKIN(3),PTHMIN,CHKIN(2),PTHMAX 
        WRITE(MSTU(11),2900) CHKIN(3),CKIN(6)   
        DO 230 I=4,14   
  230   WRITE(MSTU(11),2700) CKIN(2*I-1),CHKIN(I),CKIN(2*I) 
        SPRMAX=CKIN(32) 
        IF(SPRMAX.LT.0.) SPRMAX=VINT(1) 
        WRITE(MSTU(11),2700) CKIN(31),CHKIN(13),SPRMAX  
        WRITE(MSTU(11),3000)    
        WRITE(MSTU(11),3100)    
        WRITE(MSTU(11),2600)    
        DO 240 I=16,21  
  240   WRITE(MSTU(11),2700) VINT(I-5),CHKIN(I),VINT(I+15)  
        WRITE(MSTU(11),3000)    
    
C...Status codes and parameter values.  
      ELSEIF(MSTAT.EQ.5) THEN   
        WRITE(MSTU(11),3200)    
        WRITE(MSTU(11),3300)    
        DO 250 I=1,100  
  250   WRITE(MSTU(11),3400) I,MSTP(I),PARP(I),100+I,MSTP(100+I),   
     &  PARP(100+I) 
      ENDIF 
    
C...Formats for printouts.  
 1000 FORMAT('1',9('*'),1X,'PYSTAT:  Statistics on Number of ', 
     &'Events and Cross-sections',1X,9('*'))    
 1100 FORMAT(/1X,78('=')/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',12X, 
     &'Subprocess',12X,'I',6X,'Number of points',6X,'I',4X,'Sigma',3X,  
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',34('-'),'I',28('-'),    
     &'I',4X,'(mb)',4X,'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',1X,    
     &'N:o',1X,'Type',25X,'I',4X,'Generated',9X,'Tried',1X,'I',12X, 
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/1X,'I',34X,'I',28X, 
     &'I',12X,'I')  
 1200 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I12,1X,I13,1X,'I',1X,1P, 
     &E10.3,1X,'I') 
 1300 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')//    
     &1X,'********* Fraction of events that fail fragmentation ',   
     &'cuts =',1X,F8.5,' *********'/)   
 1400 FORMAT('1',17('*'),1X,'PYSTAT:  Decay Widths and Branching ', 
     &'Ratios',1X,17('*'))  
 1500 FORMAT(/1X,78('=')/1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/ 
     &1X,'I',1X,'Branching/Decay Channel',5X,'I',1X,'Width (GeV)',1X,   
     &'I',7X,'B.R.',1X,'I',1X,'Stat',1X,'I',2X,'Eff. B.R.',1X,'I'/1X,   
     &'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,78('='))    
 1600 FORMAT(1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,'I',1X,   
     &A12,1X,'->',13X,'I',2X,1P,E10.3,0P,1X,'I',1X,1P,E10.3,0P,1X,'I',  
     &1X,A4,1X,'I',1X,1P,E10.3,0P,1X,'I')   
 1700 FORMAT(1X,'I',1X,A12,1X,'+',1X,A12,1X,'I',2X,1P,E10.3,0P,1X,'I',  
     &1X,1P,E10.3,0P,1X,'I',1X,A4,1X,'I',1X,1P,E10.3,0P,1X,'I') 
 1800 FORMAT(1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,78('='))  
 1900 FORMAT('1',7('*'),1X,'PYSTAT: Allowed Incoming Partons/', 
     &'Particles at Hard Interaction',1X,7('*'))    
 2000 FORMAT(/1X,78('=')/1X,'I',38X,'I',37X,'I'/1X,'I',1X,  
     &'Beam particle:',1X,A,10X,'I',1X,'Target particle:',1X,A,7X,  
     &'I'/1X,'I',38X,'I',37X,'I'/1X,'I',1X,'Content',9X,'State',16X,    
     &'I',1X,'Content',9X,'State',15X,'I'/1X,'I',38X,'I',37X,'I'/1X,    
     &78('=')/1X,'I',38X,'I',37X,'I')   
 2100 FORMAT(1X,'I',1X,A,5X,A,16X,'I',1X,A,5X,A,15X,'I')    
 2200 FORMAT(1X,'I',38X,'I',1X,A,5X,A,15X,'I')  
 2300 FORMAT(1X,'I',1X,A,5X,A,16X,'I',37X,'I')  
 2400 FORMAT(1X,'I',38X,'I',37X,'I'/1X,78('=')) 
 2500 FORMAT('1',12('*'),1X,'PYSTAT: User-Defined Limits on ',  
     &'Kinematical Variables',1X,12('*'))   
 2600 FORMAT(/1X,78('=')/1X,'I',76X,'I')    
 2700 FORMAT(1X,'I',16X,1P,E10.3,0P,1X,'<',1X,A,1X,'<',1X,1P,E10.3,0P,  
     &16X,'I')  
 2800 FORMAT(1X,'I',3X,1P,E10.3,0P,1X,'(',1P,E10.3,0P,')',1X,'<',1X,A,  
     &1X,'<',1X,1P,E10.3,0P,16X,'I')    
 2900 FORMAT(1X,'I',29X,A,1X,'=',1X,1P,E10.3,0P,16X,'I')    
 3000 FORMAT(1X,'I',76X,'I'/1X,78('=')) 
 3100 FORMAT(////1X,5('*'),1X,'PYSTAT: Derived Limits on Kinematical ', 
     &'Variables Used in Generation',1X,5('*')) 
 3200 FORMAT('1',12('*'),1X,'PYSTAT: Summary of Status Codes and ', 
     &'Parameter Values',1X,12('*'))    
 3300 FORMAT(/3X,'I',4X,'MSTP(I)',9X,'PARP(I)',20X,'I',4X,'MSTP(I)',9X, 
     &'PARP(I)'/)   
 3400 FORMAT(1X,I3,5X,I6,6X,1P,E10.3,0P,18X,I3,5X,I6,6X,1P,E10.3)   
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYINKI(CHFRAM,CHBEAM,CHTARG,WIN)   
    
C...Identifies the two incoming particles and sets up kinematics,   
C...including rotations and boosts to/from CM frame.    
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHCOM(3)*8,CHALP(2)*26,  
     &CHIDNT(3)*8,CHTEMP*8,CHCDE(18)*8,CHINIT*76    
      DIMENSION LEN(3),KCDE(18) 
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',  
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/ 
      DATA CHCDE/'e-      ','e+      ','nue     ','nue~    ',   
     &'mu-     ','mu+     ','numu    ','numu~   ','tau-    ',   
     &'tau+    ','nutau   ','nutau~  ','pi+     ','pi-     ',   
     &'n       ','n~      ','p       ','p~      '/  
      DATA KCDE/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,  
     &211,-211,2112,-2112,2212,-2212/   
    
C...Convert character variables to lowercase and find their length. 
      CHCOM(1)=CHFRAM   
      CHCOM(2)=CHBEAM   
      CHCOM(3)=CHTARG   
      DO 120 I=1,3  
      LEN(I)=8  
      DO 100 LL=8,1,-1  
      IF(LEN(I).EQ.LL.AND.CHCOM(I)(LL:LL).EQ.' ') LEN(I)=LL-1   
      DO 100 LA=1,26    
  100 IF(CHCOM(I)(LL:LL).EQ.CHALP(2)(LA:LA)) CHCOM(I)(LL:LL)=   
     &CHALP(1)(LA:LA)   
      CHIDNT(I)=CHCOM(I)    
      DO 110 LL=1,6 
      IF(CHIDNT(I)(LL:LL+2).EQ.'bar') THEN  
        CHTEMP=CHIDNT(I)    
        CHIDNT(I)=CHTEMP(1:LL-1)//'~'//CHTEMP(LL+3:8)//'  ' 
      ENDIF 
  110 CONTINUE  
      DO 120 LL=1,8 
      IF(CHIDNT(I)(LL:LL).EQ.'_') THEN  
        CHTEMP=CHIDNT(I)    
        CHIDNT(I)=CHTEMP(1:LL-1)//CHTEMP(LL+1:8)//' '   
      ENDIF 
  120 CONTINUE  
    
C...Set initial state. Error for unknown codes. Reset variables.    
      N=2   
      DO 140 I=1,2  
      K(I,2)=0  
      DO 130 J=1,18 
  130 IF(CHIDNT(I+1).EQ.CHCDE(J)) K(I,2)=KCDE(J)    
      P(I,5)=ULMASS(K(I,2)) 
      MINT(40+I)=1  
      IF(IABS(K(I,2)).GT.100) MINT(40+I)=2  
      DO 140 J=1,5  
  140 V(I,J)=0. 
      IF(K(1,2).EQ.0) WRITE(MSTU(11),1000) CHBEAM(1:LEN(2)) 
      IF(K(2,2).EQ.0) WRITE(MSTU(11),1100) CHTARG(1:LEN(3)) 
      IF(K(1,2).EQ.0.OR.K(2,2).EQ.0) STOP   
      DO 150 J=6,10 
  150 VINT(J)=0.    
      CHINIT=' '    
    
C...Set up kinematics for events defined in CM frame.   
      IF(CHCOM(1)(1:2).EQ.'cm') THEN    
        IF(CHCOM(2)(1:1).NE.'e') THEN   
          LOFFS=(34-(LEN(2)+LEN(3)))/2  
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for a '//  
     &    CHCOM(2)(1:LEN(2))//'-'//CHCOM(3)(1:LEN(3))//' collider'//' ' 
        ELSE    
          LOFFS=(33-(LEN(2)+LEN(3)))/2  
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for an '// 
     &    CHCOM(2)(1:LEN(2))//'-'//CHCOM(3)(1:LEN(3))//' collider'//' ' 
        ENDIF   
C        WRITE(MSTU(11),1200) CHINIT 
C        WRITE(MSTU(11),1300) WIN    
        S=WIN**2    
        P(1,1)=0.   
        P(1,2)=0.   
        P(2,1)=0.   
        P(2,2)=0.   
        P(1,3)=SQRT(((S-P(1,5)**2-P(2,5)**2)**2-(2.*P(1,5)*P(2,5))**2)/ 
     &  (4.*S)) 
        P(2,3)=-P(1,3)  
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)    
        P(2,4)=SQRT(P(2,3)**2+P(2,5)**2)    
    
C...Set up kinematics for fixed target events.  
      ELSEIF(CHCOM(1)(1:3).EQ.'fix') THEN   
        LOFFS=(29-(LEN(2)+LEN(3)))/2    
        CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//  
     &  CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//    
     &  ' fixed target'//' '    
C        WRITE(MSTU(11),1200) CHINIT 
C        WRITE(MSTU(11),1400) WIN    
        P(1,1)=0.   
        P(1,2)=0.   
        P(2,1)=0.   
        P(2,2)=0.   
        P(1,3)=WIN  
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)    
        P(2,3)=0.   
        P(2,4)=P(2,5)   
        S=P(1,5)**2+P(2,5)**2+2.*P(2,4)*P(1,4)  
        VINT(10)=P(1,3)/(P(1,4)+P(2,4)) 
        CALL LUROBO(0.,0.,0.,0.,-VINT(10))  
C        WRITE(MSTU(11),1500) SQRT(S)    
    
C...Set up kinematics for events in user-defined frame. 
      ELSEIF(CHCOM(1)(1:3).EQ.'use') THEN   
        LOFFS=(13-(LEN(1)+LEN(2)))/2    
        CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//  
     &  CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//    
     &  'user-specified configuration'//' ' 
C        WRITE(MSTU(11),1200) CHINIT 
C        WRITE(MSTU(11),1600)    
C        WRITE(MSTU(11),1700) CHCOM(2),P(1,1),P(1,2),P(1,3)  
C        WRITE(MSTU(11),1700) CHCOM(3),P(2,1),P(2,2),P(2,3)  
        P(1,4)=SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)    
        P(2,4)=SQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)    
        DO 160 J=1,3    
  160   VINT(7+J)=(DBLE(P(1,J))+DBLE(P(2,J)))/DBLE(P(1,4)+P(2,4))   
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))  
        VINT(7)=ULANGL(P(1,1),P(1,2))   
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)   
        VINT(6)=ULANGL(P(1,3),P(1,1))   
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)   
        S=P(1,5)**2+P(2,5)**2+2.*(P(1,4)*P(2,4)-P(1,3)*P(2,3))  
C        WRITE(MSTU(11),1500) SQRT(S)    
    
C...Unknown frame. Error for too low CM energy. 
      ELSE  
        WRITE(MSTU(11),1800) CHFRAM(1:LEN(1))   
        STOP    
      ENDIF 
      IF(S.LT.PARP(2)**2) THEN  
        WRITE(MSTU(11),1900) SQRT(S)    
        STOP    
      ENDIF 
    
C...Save information on incoming particles. 
      MINT(11)=K(1,2)   
      MINT(12)=K(2,2)   
      MINT(43)=2*MINT(41)+MINT(42)-2    
      VINT(1)=SQRT(S)   
      VINT(2)=S 
      VINT(3)=P(1,5)    
      VINT(4)=P(2,5)    
      VINT(5)=P(1,3)    
    
C...Store constants to be used in generation.   
      IF(MSTP(82).LE.1) VINT(149)=4.*PARP(81)**2/S  
      IF(MSTP(82).GE.2) VINT(149)=4.*PARP(82)**2/S  
    
C...Formats for initialization and error information.   
 1000 FORMAT(1X,'Error: unrecognized beam particle ''',A,'''.'/ 
     &1X,'Execution stopped!')  
 1100 FORMAT(1X,'Error: unrecognized target particle ''',A,'''.'/   
     &1X,'Execution stopped!')  
 1200 FORMAT(/1X,78('=')/1X,'I',76X,'I'/1X,'I',A76,'I') 
 1300 FORMAT(1X,'I',18X,'at',1X,F10.3,1X,'GeV center-of-mass energy',   
     &19X,'I'/1X,'I',76X,'I'/1X,78('='))    
 1400 FORMAT(1X,'I',22X,'at',1X,F10.3,1X,'GeV/c lab-momentum',22X,'I')  
 1500 FORMAT(1X,'I',76X,'I'/1X,'I',11X,'corresponding to',1X,F10.3,1X,  
     &'GeV center-of-mass energy',12X,'I'/1X,'I',76X,'I'/1X,78('='))    
 1600 FORMAT(1X,'I',76X,'I'/1X,'I',24X,'px (GeV/c)',3X,'py (GeV/c)',3X, 
     &'pz (GeV/c)',16X,'I') 
 1700 FORMAT(1X,'I',15X,A8,3(2X,F10.3,1X),15X,'I')  
 1800 FORMAT(1X,'Error: unrecognized coordinate frame ''',A,'''.'/  
     &1X,'Execution stopped!')  
 1900 FORMAT(1X,'Error: too low CM energy,',F8.3,' GeV for event ', 
     &'generation.'/1X,'Execution stopped!')    
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYINRE 
    
C...Calculates full and effective widths of guage bosons, stores masses 
C...and widths, rescales coefficients to be used for resonance  
C...production generation.  
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT6/PROC(0:200) 
      CHARACTER PROC*28 
      SAVE /PYINT6/ 
      DIMENSION WDTP(0:40),WDTE(0:40,0:5)   
    
C...Calculate full and effective widths of gauge bosons.    
      AEM=PARU(101) 
      XW=PARU(102)  
      DO 100 I=21,40    
      DO 100 J=0,40 
      WIDP(I,J)=0.  
  100 WIDE(I,J)=0.  
    
C...W+/-:   
      WMAS=PMAS(24,1)   
      WFAC=AEM/(24.*XW)*WMAS    
      CALL PYWIDT(24,WMAS,WDTP,WDTE)    
      WIDS(24,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(24,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(24,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 110 I=0,40 
      WIDP(24,I)=WFAC*WDTP(I)   
  110 WIDE(24,I)=WFAC*WDTE(I,0) 
    
C...H+/-:   
      HCMAS=PMAS(37,1)  
      HCFAC=AEM/(8.*XW)*(HCMAS/WMAS)**2*HCMAS   
      CALL PYWIDT(37,HCMAS,WDTP,WDTE)   
      WIDS(37,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(37,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(37,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 120 I=0,40 
      WIDP(37,I)=HCFAC*WDTP(I)  
  120 WIDE(37,I)=HCFAC*WDTE(I,0)    
    
C...Z0: 
      ZMAS=PMAS(23,1)   
      ZFAC=AEM/(48.*XW*(1.-XW))*ZMAS    
      CALL PYWIDT(23,ZMAS,WDTP,WDTE)    
      WIDS(23,1)=((WDTE(0,1)+WDTE(0,2))**2+ 
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(23,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(23,3)=0. 
      DO 130 I=0,40 
      WIDP(23,I)=ZFAC*WDTP(I)   
  130 WIDE(23,I)=ZFAC*WDTE(I,0) 
    
C...H0: 
      HMAS=PMAS(25,1)   
      HFAC=AEM/(8.*XW)*(HMAS/WMAS)**2*HMAS  
      CALL PYWIDT(25,HMAS,WDTP,WDTE)    
      WIDS(25,1)=((WDTE(0,1)+WDTE(0,2))**2+ 
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(25,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(25,3)=0. 
      DO 140 I=0,40 
      WIDP(25,I)=HFAC*WDTP(I)   
  140 WIDE(25,I)=HFAC*WDTE(I,0) 
    
C...Z'0:    
      ZPMAS=PMAS(32,1)  
      ZPFAC=AEM/(48.*XW*(1.-XW))*ZPMAS  
      CALL PYWIDT(32,ZPMAS,WDTP,WDTE)   
      WIDS(32,1)=((WDTE(0,1)+WDTE(0,2)+WDTE(0,3))**2+   
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(32,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(32,3)=0. 
      DO 150 I=0,40 
      WIDP(32,I)=ZPFAC*WDTP(I)  
  150 WIDE(32,I)=ZPFAC*WDTE(I,0)    
    
C...R:  
      RMAS=PMAS(40,1)   
      RFAC=0.08*RMAS/((MSTP(1)-1)*(1.+6.*(1.+ULALPS(RMAS**2)/PARU(1)))) 
      CALL PYWIDT(40,RMAS,WDTP,WDTE)    
      WIDS(40,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(40,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(40,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 160 I=0,40 
      WIDP(40,I)=WFAC*WDTP(I)   
  160 WIDE(40,I)=WFAC*WDTE(I,0) 
    
C...Q:  
      KFLQM=1   
      DO 170 I=1,MIN(8,MDCY(21,3))  
      IDC=I+MDCY(21,2)-1    
      IF(MDME(IDC,1).LE.0) GOTO 170 
      KFLQM=I   
  170 CONTINUE  
      MINT(46)=KFLQM    
      KFPR(81,1)=KFLQM  
      KFPR(81,2)=KFLQM  
      KFPR(82,1)=KFLQM  
      KFPR(82,2)=KFLQM  
    
C...Set resonance widths and branching ratios in JETSET.    
      DO 180 I=1,6  
      IF(I.LE.3) KC=I+22    
      IF(I.EQ.4) KC=32  
      IF(I.EQ.5) KC=37  
      IF(I.EQ.6) KC=40  
      PMAS(KC,2)=WIDP(KC,0) 
      PMAS(KC,3)=MIN(0.9*PMAS(KC,1),10.*PMAS(KC,2)) 
      DO 180 J=1,MDCY(KC,3) 
      IDC=J+MDCY(KC,2)-1    
      BRAT(IDC)=WIDE(KC,J)/WIDE(KC,0)   
  180 CONTINUE  
    
C...Special cases in treatment of gamma*/Z0: redefine process name. 
      IF(MSTP(43).EQ.1) THEN    
        PROC(1)='f + fb -> gamma*'  
      ELSEIF(MSTP(43).EQ.2) THEN    
        PROC(1)='f + fb -> Z0'  
      ELSEIF(MSTP(43).EQ.3) THEN    
        PROC(1)='f + fb -> gamma*/Z0'   
      ENDIF 
    
C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name. 
      IF(MSTP(44).EQ.1) THEN    
        PROC(141)='f + fb -> gamma*'    
      ELSEIF(MSTP(44).EQ.2) THEN    
        PROC(141)='f + fb -> Z0'    
      ELSEIF(MSTP(44).EQ.3) THEN    
        PROC(141)='f + fb -> Z''0'  
      ELSEIF(MSTP(44).EQ.4) THEN    
        PROC(141)='f + fb -> gamma*/Z0' 
      ELSEIF(MSTP(44).EQ.5) THEN    
        PROC(141)='f + fb -> gamma*/Z''0'   
      ELSEIF(MSTP(44).EQ.6) THEN    
        PROC(141)='f + fb -> Z0/Z''0'   
      ELSEIF(MSTP(44).EQ.7) THEN    
        PROC(141)='f + fb -> gamma*/Z0/Z''0'    
      ENDIF 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYXTOT 
    
C...Parametrizes total, double diffractive, single diffractive and  
C...elastic cross-sections for different energies and beams.    
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      DIMENSION BCS(5,8),BCB(2,5),BCC(3)    
    
C...The following data lines are coefficients needed in the 
C...Block, Cahn parametrization of total cross-section and nuclear  
C...slope parameter; see below. 
      DATA ((BCS(I,J),J=1,8),I=1,5)/    
     1 41.74, 0.66, 0.0000, 337.,  0.0, 0.0, -39.3, 0.48,   
     2 41.66, 0.60, 0.0000, 306.,  0.0, 0.0, -34.6, 0.51,   
     3 41.36, 0.63, 0.0000, 299.,  7.3, 0.5, -40.4, 0.47,   
     4 41.68, 0.63, 0.0083, 330.,  0.0, 0.0, -39.0, 0.48,   
     5 41.13, 0.59, 0.0074, 278., 10.5, 0.5, -41.2, 0.46/   
      DATA ((BCB(I,J),J=1,5),I=1,2)/    
     1 10.79, -0.049, 0.040, 21.5, 1.23,    
     2  9.92, -0.027, 0.013, 18.9, 1.07/    
      DATA BCC/2.0164346,-0.5590311,0.0376279/  
    
C...Total cross-section and nuclear slope parameter for pp and p-pbar   
      NFIT=MIN(5,MAX(1,MSTP(31)))   
      SIGP=BCS(NFIT,1)+BCS(NFIT,2)*(-0.25*PARU(1)**2*   
     &(1.-0.25*BCS(NFIT,3)*PARU(1)**2)+(1.+0.5*BCS(NFIT,3)*PARU(1)**2)* 
     &(LOG(VINT(2)/BCS(NFIT,4)))**2+BCS(NFIT,3)*    
     &(LOG(VINT(2)/BCS(NFIT,4)))**4)/   
     &((1.-0.25*BCS(NFIT,3)*PARU(1)**2)**2+2.*BCS(NFIT,3)*  
     &(1.+0.25*BCS(NFIT,3)*PARU(1)**2)*(LOG(VINT(2)/BCS(NFIT,4)))**2+   
     &BCS(NFIT,3)**2*(LOG(VINT(2)/BCS(NFIT,4)))**4)+BCS(NFIT,5)*    
     &VINT(2)**(BCS(NFIT,6)-1.)*SIN(0.5*PARU(1)*BCS(NFIT,6))    
      SIGM=-BCS(NFIT,7)*VINT(2)**(BCS(NFIT,8)-1.)*  
     &COS(0.5*PARU(1)*BCS(NFIT,8))  
      REFP=BCS(NFIT,2)*PARU(1)*LOG(VINT(2)/BCS(NFIT,4))/    
     &((1.-0.25*BCS(NFIT,3)*PARU(1)**2)**2+2.*BCS(NFIT,3)*  
     &(1.+0.25*BCS(NFIT,3)*PARU(1)**2)+(LOG(VINT(2)/BCS(NFIT,4)))**2+   
     &BCS(NFIT,3)**2*(LOG(VINT(2)/BCS(NFIT,4)))**4)-BCS(NFIT,5)*    
     &VINT(2)**(BCS(NFIT,6)-1.)*COS(0.5*PARU(1)*BCS(NFIT,6))    
      REFM=-BCS(NFIT,7)*VINT(2)**(BCS(NFIT,8)-1.)*  
     &SIN(0.5*PARU(1)*BCS(NFIT,8))  
      SIGMA=SIGP-ISIGN(1,MINT(11)*MINT(12))*SIGM    
      RHO=(REFP-ISIGN(1,MINT(11)*MINT(12))*REFM)/SIGMA  
    
C...Nuclear slope parameter B, curvature C: 
      NFIT=1    
      IF(MSTP(31).GE.4) NFIT=2  
      BP=BCB(NFIT,1)+BCB(NFIT,2)*LOG(VINT(2))+  
     &BCB(NFIT,3)*(LOG(VINT(2)))**2 
      BM=BCB(NFIT,4)+BCB(NFIT,5)*LOG(VINT(2))   
      B=BP-ISIGN(1,MINT(11)*MINT(12))*SIGM/SIGP*(BM-BP) 
      VINT(121)=B   
      C=-0.5*BCC(2)/BCC(3)*(1.-SQRT(MAX(0.,1.+4.*BCC(3)/BCC(2)**2*  
     &(1.E-03*VINT(1)-BCC(1)))))    
      VINT(122)=C   
    
C...Elastic scattering cross-section (fixed by sigma-tot, rho and B).   
      SIGEL=SIGMA**2*(1.+RHO**2)/(16.*PARU(1)*PARU(5)*B)    
    
C...Single diffractive scattering cross-section from Goulianos: 
      SIGSD=2.*0.68*(1.+36./VINT(2))*LOG(0.6+0.1*VINT(2))   
    
C...Double diffractive scattering cross-section (essentially fixed by   
C...sigma-sd and sigma-el). 
      SIGDD=SIGSD**2/(3.*SIGEL) 
    
C...Total non-elastic, non-diffractive cross-section.   
      SIGND=SIGMA-SIGDD-SIGSD-SIGEL 
    
C...Rescale for pions.  
      IF(IABS(MINT(11)).EQ.211.AND.IABS(MINT(12)).EQ.211) THEN  
        SIGMA=4./9.*SIGMA   
        SIGDD=4./9.*SIGDD   
        SIGSD=4./9.*SIGSD   
        SIGEL=4./9.*SIGEL   
        SIGND=4./9.*SIGND   
      ELSEIF(IABS(MINT(11)).EQ.211.OR.IABS(MINT(12)).EQ.211) THEN   
        SIGMA=2./3.*SIGMA   
        SIGDD=2./3.*SIGDD   
        SIGSD=2./3.*SIGSD   
        SIGEL=2./3.*SIGEL   
        SIGND=2./3.*SIGND   
      ENDIF 
    
C...Save cross-sections in common block PYPARA. 
      VINT(101)=SIGMA   
      VINT(102)=SIGEL   
      VINT(103)=SIGSD   
      VINT(104)=SIGDD   
      VINT(106)=SIGND   
      XSEC(95,1)=SIGND  
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYMAXI 
    
C...Finds optimal set of coefficients for kinematical variable selection    
C...and the maximum of the part of the differential cross-section used  
C...in the event weighting. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      COMMON/PYINT6/PROC(0:200) 
      CHARACTER PROC*28 
      SAVE /PYINT6/ 
      CHARACTER CVAR(4)*4   
      DIMENSION NPTS(4),MVARPT(200,4),VINTPT(200,30),SIGSPT(200),   
     &NAREL(6),WTREL(6),WTMAT(6,6),COEFU(6),IACCMX(4),SIGSMX(4),    
     &SIGSSM(3) 
      DATA CVAR/'tau ','tau''','y*  ','cth '/   
    
C...Select subprocess to study: skip cases not applicable.  
      VINT(143)=1.  
      VINT(144)=1.  
      XSEC(0,1)=0.  
      DO 350 ISUB=1,200 
      IF(ISUB.GE.91.AND.ISUB.LE.95) THEN    
        XSEC(ISUB,1)=VINT(ISUB+11)  
        IF(MSUB(ISUB).NE.1) GOTO 350    
        GOTO 340    
      ELSEIF(ISUB.EQ.96) THEN   
        IF(MINT(43).NE.4) GOTO 350  
        IF(MSUB(95).NE.1.AND.MSTP(81).LE.0.AND.MSTP(131).LE.0) GOTO 350 
      ELSEIF(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.ISUB.EQ.28.OR.   
     &ISUB.EQ.53.OR.ISUB.EQ.68) THEN    
        IF(MSUB(ISUB).NE.1.OR.MSUB(95).EQ.1) GOTO 350   
      ELSE  
        IF(MSUB(ISUB).NE.1) GOTO 350    
      ENDIF 
      MINT(1)=ISUB  
      ISTSB=ISET(ISUB)  
      IF(ISUB.EQ.96) ISTSB=2    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1000) ISUB  
    
C...Find resonances (explicit or implicit in cross-section).    
      MINT(72)=0    
      KFR1=0    
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3) THEN 
        KFR1=KFPR(ISUB,1)   
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN    
        KFR1=25 
      ENDIF 
      IF(KFR1.NE.0) THEN    
        TAUR1=PMAS(KFR1,1)**2/VINT(2)   
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2) 
        MINT(72)=1  
        MINT(73)=KFR1   
        VINT(73)=TAUR1  
        VINT(74)=GAMR1  
      ENDIF 
      IF(ISUB.EQ.141) THEN  
        KFR2=23 
        TAUR2=PMAS(KFR2,1)**2/VINT(2)   
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2) 
        MINT(72)=2  
        MINT(74)=KFR2   
        VINT(75)=TAUR2  
        VINT(76)=GAMR2  
      ENDIF 
    
C...Find product masses and minimum pT of process.  
      SQM3=0.   
      SQM4=0.   
      MINT(71)=0    
      VINT(71)=CKIN(3)  
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN 
        IF(KFPR(ISUB,1).NE.0) SQM3=PMAS(KFPR(ISUB,1),1)**2  
        IF(KFPR(ISUB,2).NE.0) SQM4=PMAS(KFPR(ISUB,2),1)**2  
        IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) MINT(71)=1 
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5)) 
        IF(ISUB.EQ.96.AND.MSTP(82).LE.1) VINT(71)=PARP(81)  
        IF(ISUB.EQ.96.AND.MSTP(82).GE.2) VINT(71)=0.08*PARP(82) 
      ENDIF 
      VINT(63)=SQM3 
      VINT(64)=SQM4 
    
C...Number of points for each variable: tau, tau', y*, cos(theta-hat).  
      NPTS(1)=2+2*MINT(72)  
      IF(MINT(43).EQ.1.AND.(ISTSB.EQ.1.OR.ISTSB.EQ.2)) NPTS(1)=1    
      NPTS(2)=1 
      IF(MINT(43).GE.2.AND.(ISTSB.EQ.3.OR.ISTSB.EQ.4)) NPTS(2)=2    
      NPTS(3)=1 
      IF(MINT(43).EQ.4) NPTS(3)=3   
      NPTS(4)=1 
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) NPTS(4)=5    
      NTRY=NPTS(1)*NPTS(2)*NPTS(3)*NPTS(4)  
    
C...Reset coefficients of cross-section weighting.  
      DO 100 J=1,20 
  100 COEF(ISUB,J)=0.   
      COEF(ISUB,1)=1.   
      COEF(ISUB,7)=0.5  
      COEF(ISUB,8)=0.5  
      COEF(ISUB,10)=1.  
      COEF(ISUB,15)=1.  
      MCTH=0    
      MTAUP=0   
      CTH=0.    
      TAUP=0.   
      SIGSAM=0. 
    
C...Find limits and select tau, y*, cos(theta-hat) and tau' values, 
C...in grid of phase space points.  
      CALL PYKLIM(1)    
      NACC=0    
      DO 120 ITRY=1,NTRY    
      IF(MOD(ITRY-1,NPTS(2)*NPTS(3)*NPTS(4)).EQ.0) THEN 
        MTAU=1+(ITRY-1)/(NPTS(2)*NPTS(3)*NPTS(4))   
        CALL PYKMAP(1,MTAU,0.5) 
        IF(ISTSB.EQ.3.OR.ISTSB.EQ.4) CALL PYKLIM(4) 
      ENDIF 
      IF((ISTSB.EQ.3.OR.ISTSB.EQ.4).AND.MOD(ITRY-1,NPTS(3)*NPTS(4)).    
     &EQ.0) THEN    
        MTAUP=1+MOD((ITRY-1)/(NPTS(3)*NPTS(4)),NPTS(2)) 
        CALL PYKMAP(4,MTAUP,0.5)    
      ENDIF 
      IF(MOD(ITRY-1,NPTS(3)*NPTS(4)).EQ.0) CALL PYKLIM(2)   
      IF(MOD(ITRY-1,NPTS(4)).EQ.0) THEN 
        MYST=1+MOD((ITRY-1)/NPTS(4),NPTS(3))    
        CALL PYKMAP(2,MYST,0.5) 
        CALL PYKLIM(3)  
      ENDIF 
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN 
        MCTH=1+MOD(ITRY-1,NPTS(4))  
        CALL PYKMAP(3,MCTH,0.5) 
      ENDIF 
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2) 
    
C...Calculate and store cross-section.  
      MINT(51)=0    
      CALL PYKLIM(0)    
      IF(MINT(51).EQ.1) GOTO 120    
      NACC=NACC+1   
      MVARPT(NACC,1)=MTAU   
      MVARPT(NACC,2)=MTAUP  
      MVARPT(NACC,3)=MYST   
      MVARPT(NACC,4)=MCTH   
      DO 110 J=1,30 
  110 VINTPT(NACC,J)=VINT(10+J) 
      CALL PYSIGH(NCHN,SIGS)    
      SIGSPT(NACC)=SIGS 
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1100) MTAU,MTAUP,MYST,MCTH, 
     &VINT(21),VINT(22),VINT(23),VINT(26),SIGS  
  120 CONTINUE  
      IF(SIGSAM.EQ.0.) THEN 
        WRITE(MSTU(11),1200) ISUB   
        STOP    
      ENDIF 
    
C...Calculate integrals in tau and y* over maximal phase space limits.  
      TAUMIN=VINT(11)   
      TAUMAX=VINT(31)   
      ATAU1=LOG(TAUMAX/TAUMIN)  
      ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN) 
      IF(NPTS(1).GE.3) THEN 
        ATAU3=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1    
        ATAU4=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/  
     &  GAMR1   
      ENDIF 
      IF(NPTS(1).GE.5) THEN 
        ATAU5=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2    
        ATAU6=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/  
     &  GAMR2   
      ENDIF 
      YSTMIN=0.5*LOG(TAUMIN)    
      YSTMAX=-YSTMIN    
      AYST0=YSTMAX-YSTMIN   
      AYST1=0.5*(YSTMAX-YSTMIN)**2  
      AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))    
    
C...Reset. Sum up cross-sections in points calculated.  
      DO 230 IVAR=1,4   
      IF(NPTS(IVAR).EQ.1) GOTO 230  
      IF(ISUB.EQ.96.AND.IVAR.EQ.4) GOTO 230 
      NBIN=NPTS(IVAR)   
      DO 130 J1=1,NBIN  
      NAREL(J1)=0   
      WTREL(J1)=0.  
      COEFU(J1)=0.  
      DO 130 J2=1,NBIN  
  130 WTMAT(J1,J2)=0.   
      DO 140 IACC=1,NACC    
      IBIN=MVARPT(IACC,IVAR)    
      NAREL(IBIN)=NAREL(IBIN)+1 
      WTREL(IBIN)=WTREL(IBIN)+SIGSPT(IACC)  
    
C...Sum up tau cross-section pieces in points used. 
      IF(IVAR.EQ.1) THEN    
        TAU=VINTPT(IACC,11) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAU1/ATAU2)/TAU   
        IF(NBIN.GE.3) THEN  
          WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAU1/ATAU3)/(TAU+TAUR1) 
          WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ATAU1/ATAU4)*TAU/    
     &    ((TAU-TAUR1)**2+GAMR1**2) 
        ENDIF   
        IF(NBIN.GE.5) THEN  
          WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ATAU1/ATAU5)/(TAU+TAUR2) 
          WTMAT(IBIN,6)=WTMAT(IBIN,6)+(ATAU1/ATAU6)*TAU/    
     &    ((TAU-TAUR2)**2+GAMR2**2) 
        ENDIF   
    
C...Sum up tau' cross-section pieces in points used.    
      ELSEIF(IVAR.EQ.2) THEN    
        TAU=VINTPT(IACC,11) 
        TAUP=VINTPT(IACC,16)    
        TAUPMN=VINTPT(IACC,6)   
        TAUPMX=VINTPT(IACC,26)  
        ATAUP1=LOG(TAUPMX/TAUPMN)   
        ATAUP2=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAUP1/ATAUP2)*(1.-TAU/TAUP)**3/   
     &  TAUP    
    
C...Sum up y* and cos(theta-hat) cross-section pieces in points used.   
      ELSEIF(IVAR.EQ.3) THEN    
        YST=VINTPT(IACC,12) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+(AYST0/AYST1)*(YST-YSTMIN)  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(AYST0/AYST1)*(YSTMAX-YST)  
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(AYST0/AYST3)/COSH(YST) 
      ELSE  
        RM34=2.*SQM3*SQM4/(VINTPT(IACC,11)*VINT(2))**2  
        RSQM=1.+RM34    
        CTHMAX=SQRT(1.-4.*VINT(71)**2/(TAUMAX*VINT(2))) 
        CTHMIN=-CTHMAX  
        IF(CTHMAX.GT.0.9999) RM34=MAX(RM34,2.*VINT(71)**2/  
     &  (TAUMAX*VINT(2)))   
        ACTH1=CTHMAX-CTHMIN 
        ACTH2=LOG(MAX(RM34,RSQM-CTHMIN)/MAX(RM34,RSQM-CTHMAX))  
        ACTH3=LOG(MAX(RM34,RSQM+CTHMAX)/MAX(RM34,RSQM+CTHMIN))  
        ACTH4=1./MAX(RM34,RSQM-CTHMAX)-1./MAX(RM34,RSQM-CTHMIN) 
        ACTH5=1./MAX(RM34,RSQM+CTHMIN)-1./MAX(RM34,RSQM+CTHMAX) 
        CTH=VINTPT(IACC,13) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ACTH1/ACTH2)/MAX(RM34,RSQM-CTH)    
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ACTH1/ACTH3)/MAX(RM34,RSQM+CTH)    
        WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ACTH1/ACTH4)/MAX(RM34,RSQM-CTH)**2 
        WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ACTH1/ACTH5)/MAX(RM34,RSQM+CTH)**2 
      ENDIF 
  140 CONTINUE  
    
C...Check that equation system solvable; else trivial way out.  
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1300) CVAR(IVAR)    
      MSOLV=1   
      DO 150 IBIN=1,NBIN    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1400) (WTMAT(IBIN,IRED),    
     &IRED=1,NBIN),WTREL(IBIN)  
  150 IF(NAREL(IBIN).EQ.0) MSOLV=0  
      IF(MSOLV.EQ.0) THEN   
        DO 160 IBIN=1,NBIN  
  160   COEFU(IBIN)=1.  
    
C...Solve to find relative importance of cross-section pieces.  
      ELSE  
        DO 170 IRED=1,NBIN-1    
        DO 170 IBIN=IRED+1,NBIN 
        RQT=WTMAT(IBIN,IRED)/WTMAT(IRED,IRED)   
        WTREL(IBIN)=WTREL(IBIN)-RQT*WTREL(IRED) 
        DO 170 ICOE=IRED,NBIN   
  170   WTMAT(IBIN,ICOE)=WTMAT(IBIN,ICOE)-RQT*WTMAT(IRED,ICOE)  
        DO 190 IRED=NBIN,1,-1   
        DO 180 ICOE=IRED+1,NBIN 
  180   WTREL(IRED)=WTREL(IRED)-WTMAT(IRED,ICOE)*COEFU(ICOE)    
  190   COEFU(IRED)=WTREL(IRED)/WTMAT(IRED,IRED)    
      ENDIF 
    
C...Normalize coefficients, with piece shared democratically.   
      COEFSU=0. 
      DO 200 IBIN=1,NBIN    
      COEFU(IBIN)=MAX(0.,COEFU(IBIN))   
  200 COEFSU=COEFSU+COEFU(IBIN) 
      IF(IVAR.EQ.1) IOFF=0  
      IF(IVAR.EQ.2) IOFF=14 
      IF(IVAR.EQ.3) IOFF=6  
      IF(IVAR.EQ.4) IOFF=9  
      IF(COEFSU.GT.0.) THEN 
        DO 210 IBIN=1,NBIN  
  210   COEF(ISUB,IOFF+IBIN)=PARP(121)/NBIN+(1.-PARP(121))*COEFU(IBIN)/ 
     &  COEFSU  
      ELSE  
        DO 220 IBIN=1,NBIN  
  220   COEF(ISUB,IOFF+IBIN)=1./NBIN    
      ENDIF 
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1500) CVAR(IVAR),   
     &(COEF(ISUB,IOFF+IBIN),IBIN=1,NBIN)    
  230 CONTINUE  
    
C...Find two most promising maxima among points previously determined.  
      DO 240 J=1,4  
      IACCMX(J)=0   
  240 SIGSMX(J)=0.  
      NMAX=0    
      DO 290 IACC=1,NACC    
      DO 250 J=1,30 
  250 VINT(10+J)=VINTPT(IACC,J) 
      CALL PYSIGH(NCHN,SIGS)    
      IEQ=0 
      DO 260 IMV=1,NMAX 
  260 IF(ABS(SIGS-SIGSMX(IMV)).LT.1E-4*(SIGS+SIGSMX(IMV))) IEQ=IMV  
      IF(IEQ.EQ.0) THEN 
        DO 270 IMV=NMAX,1,-1    
        IIN=IMV+1   
        IF(SIGS.LE.SIGSMX(IMV)) GOTO 280    
        IACCMX(IMV+1)=IACCMX(IMV)   
  270   SIGSMX(IMV+1)=SIGSMX(IMV)   
        IIN=1   
  280   IACCMX(IIN)=IACC    
        SIGSMX(IIN)=SIGS    
        IF(NMAX.LE.1) NMAX=NMAX+1   
      ENDIF 
  290 CONTINUE  
    
C...Read out starting position for search.  
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1600)   
      SIGSAM=SIGSMX(1)  
      DO 330 IMAX=1,NMAX    
      IACC=IACCMX(IMAX) 
      MTAU=MVARPT(IACC,1)   
      MTAUP=MVARPT(IACC,2)  
      MYST=MVARPT(IACC,3)   
      MCTH=MVARPT(IACC,4)   
      VTAU=0.5  
      VYST=0.5  
      VCTH=0.5  
      VTAUP=0.5 
    
C...Starting point and step size in parameter space.    
      DO 320 IRPT=1,2   
      DO 310 IVAR=1,4   
      IF(NPTS(IVAR).EQ.1) GOTO 310  
      IF(IVAR.EQ.1) VVAR=VTAU   
      IF(IVAR.EQ.2) VVAR=VTAUP  
      IF(IVAR.EQ.3) VVAR=VYST   
      IF(IVAR.EQ.4) VVAR=VCTH   
      IF(IVAR.EQ.1) MVAR=MTAU   
      IF(IVAR.EQ.2) MVAR=MTAUP  
      IF(IVAR.EQ.3) MVAR=MYST   
      IF(IVAR.EQ.4) MVAR=MCTH   
      IF(IRPT.EQ.1) VDEL=0.1    
      IF(IRPT.EQ.2) VDEL=MAX(0.01,MIN(0.05,VVAR-0.02,0.98-VVAR))    
      IF(IRPT.EQ.1) VMAR=0.02   
      IF(IRPT.EQ.2) VMAR=0.002  
      IMOV0=1   
      IF(IRPT.EQ.1.AND.IVAR.EQ.1) IMOV0=0   
      DO 300 IMOV=IMOV0,8   
    
C...Define new point in parameter space.    
      IF(IMOV.EQ.0) THEN    
        INEW=2  
        VNEW=VVAR   
      ELSEIF(IMOV.EQ.1) THEN    
        INEW=3  
        VNEW=VVAR+VDEL  
      ELSEIF(IMOV.EQ.2) THEN    
        INEW=1  
        VNEW=VVAR-VDEL  
      ELSEIF(SIGSSM(3).GE.MAX(SIGSSM(1),SIGSSM(2)).AND. 
     &VVAR+2.*VDEL.LT.1.-VMAR) THEN 
        VVAR=VVAR+VDEL  
        SIGSSM(1)=SIGSSM(2) 
        SIGSSM(2)=SIGSSM(3) 
        INEW=3  
        VNEW=VVAR+VDEL  
      ELSEIF(SIGSSM(1).GE.MAX(SIGSSM(2),SIGSSM(3)).AND. 
     &VVAR-2.*VDEL.GT.VMAR) THEN    
        VVAR=VVAR-VDEL  
        SIGSSM(3)=SIGSSM(2) 
        SIGSSM(2)=SIGSSM(1) 
        INEW=1  
        VNEW=VVAR-VDEL  
      ELSEIF(SIGSSM(3).GE.SIGSSM(1)) THEN   
        VDEL=0.5*VDEL   
        VVAR=VVAR+VDEL  
        SIGSSM(1)=SIGSSM(2) 
        INEW=2  
        VNEW=VVAR   
      ELSE  
        VDEL=0.5*VDEL   
        VVAR=VVAR-VDEL  
        SIGSSM(3)=SIGSSM(2) 
        INEW=2  
        VNEW=VVAR   
      ENDIF 
    
C...Convert to relevant variables and find derived new limits.  
      IF(IVAR.EQ.1) THEN    
        VTAU=VNEW   
        CALL PYKMAP(1,MTAU,VTAU)    
        IF(ISTSB.EQ.3.OR.ISTSB.EQ.4) CALL PYKLIM(4) 
      ENDIF 
      IF(IVAR.LE.2.AND.(ISTSB.EQ.3.OR.ISTSB.EQ.4)) THEN 
        IF(IVAR.EQ.2) VTAUP=VNEW    
        CALL PYKMAP(4,MTAUP,VTAUP)  
      ENDIF 
      IF(IVAR.LE.2) CALL PYKLIM(2)  
      IF(IVAR.LE.3) THEN    
        IF(IVAR.EQ.3) VYST=VNEW 
        CALL PYKMAP(2,MYST,VYST)    
        CALL PYKLIM(3)  
      ENDIF 
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN 
        IF(IVAR.EQ.4) VCTH=VNEW 
        CALL PYKMAP(3,MCTH,VCTH)    
      ENDIF 
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2) 
    
C...Evaluate cross-section. Save new maximum. Final maximum.    
      CALL PYSIGH(NCHN,SIGS)    
      SIGSSM(INEW)=SIGS 
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1700) IMAX,IVAR,MVAR,IMOV,  
     &VNEW,VINT(21),VINT(22),VINT(23),VINT(26),SIGS 
  300 CONTINUE  
  310 CONTINUE  
  320 CONTINUE  
      IF(IMAX.EQ.1) SIGS11=SIGSAM   
  330 CONTINUE  
      XSEC(ISUB,1)=1.05*SIGSAM  
  340 IF(ISUB.NE.96) XSEC(0,1)=XSEC(0,1)+XSEC(ISUB,1)   
  350 CONTINUE  
    
C...Print summary table.    
      IF(MSTP(122).GE.1) THEN   
        WRITE(MSTU(11),1800)    
        WRITE(MSTU(11),1900)    
        DO 360 ISUB=1,200   
        IF(MSUB(ISUB).NE.1.AND.ISUB.NE.96) GOTO 360 
        IF(ISUB.EQ.96.AND.MINT(43).NE.4) GOTO 360   
        IF(ISUB.EQ.96.AND.MSUB(95).NE.1.AND.MSTP(81).LE.0) GOTO 360 
        IF(MSUB(95).EQ.1.AND.(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.    
     &  ISUB.EQ.28.OR.ISUB.EQ.53.OR.ISUB.EQ.68)) GOTO 360   
        WRITE(MSTU(11),2000) ISUB,PROC(ISUB),XSEC(ISUB,1)   
  360   CONTINUE    
        WRITE(MSTU(11),2100)    
      ENDIF 
    
C...Format statements for maximization results. 
 1000 FORMAT(/1X,'Coefficient optimization and maximum search for ',    
     &'subprocess no',I4/1X,'Coefficient modes     tau',10X,'y*',9X,    
     &'cth',9X,'tau''',7X,'sigma')  
 1100 FORMAT(1X,4I4,F12.8,F12.6,F12.7,F12.8,1P,E12.4)   
 1200 FORMAT(1X,'Error: requested subprocess ',I3,' has vanishing ',    
     &'cross-section.'/1X,'Execution stopped!')
 1300 FORMAT(1X,'Coefficients of equation system to be solved for ',A4) 
 1400 FORMAT(1X,1P,7E11.3)  
 1500 FORMAT(1X,'Result for ',A4,':',6F9.4) 
 1600 FORMAT(1X,'Maximum search for given coefficients'/2X,'MAX VAR ',  
     &'MOD MOV   VNEW',7X,'tau',7X,'y*',8X,'cth',7X,'tau''',7X,'sigma') 
 1700 FORMAT(1X,4I4,F8.4,F11.7,F9.3,F11.6,F11.7,1P,E12.4)   
 1800 FORMAT(/1X,8('*'),1X,'PYMAXI: summary of differential ',  
     &'cross-section maximum search',1X,8('*')) 
 1900 FORMAT(/11X,58('=')/11X,'I',38X,'I',17X,'I'/11X,'I  ISUB  ',  
     &'Subprocess name',15X,'I  Maximum value  I'/11X,'I',38X,'I',  
     &17X,'I'/11X,58('=')/11X,'I',38X,'I',17X,'I')  
 2000 FORMAT(11X,'I',2X,I3,3X,A28,2X,'I',2X,1P,E12.4,3X,'I')    
 2100 FORMAT(11X,'I',38X,'I',17X,'I'/11X,58('='))   
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYOVLY(MOVLY)  
    
C...Initializes multiplicity distribution and selects mutliplicity  
C...of overlayed events, i.e. several events occuring at the same   
C...beam crossing.  
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      DIMENSION WTI(0:100)  
      SAVE IMAX,WTI,WTS 
    
C...Sum of allowed cross-sections for overlayed events. 
      IF(MOVLY.EQ.1) THEN   
        VINT(131)=VINT(106) 
        IF(MSTP(132).GE.2) VINT(131)=VINT(131)+VINT(104)    
        IF(MSTP(132).GE.3) VINT(131)=VINT(131)+VINT(103)    
        IF(MSTP(132).GE.4) VINT(131)=VINT(131)+VINT(102)    
    
C...Initialize multiplicity distribution for unbiased events.   
        IF(MSTP(133).EQ.1) THEN 
          XNAVE=VINT(131)*PARP(131) 
          IF(XNAVE.GT.40.) WRITE(MSTU(11),1000) XNAVE   
          WTI(0)=EXP(-MIN(50.,XNAVE))   
          WTS=0.    
          WTN=0.    
          DO 100 I=1,100    
          WTI(I)=WTI(I-1)*XNAVE/I   
          IF(I-2.5.GT.XNAVE.AND.WTI(I).LT.1E-6) GOTO 110    
          WTS=WTS+WTI(I)    
          WTN=WTN+WTI(I)*I  
  100     IMAX=I    
  110     VINT(132)=XNAVE   
          VINT(133)=WTN/WTS 
          VINT(134)=WTS 
    
C...Initialize mutiplicity distribution for biased events.  
        ELSEIF(MSTP(133).EQ.2) THEN 
          XNAVE=VINT(131)*PARP(131) 
          IF(XNAVE.GT.40.) WRITE(MSTU(11),1000) XNAVE   
          WTI(1)=EXP(-MIN(50.,XNAVE))*XNAVE 
          WTS=WTI(1)    
          WTN=WTI(1)    
          DO 120 I=2,100    
          WTI(I)=WTI(I-1)*XNAVE/(I-1)   
          IF(I-2.5.GT.XNAVE.AND.WTI(I).LT.1E-6) GOTO 130    
          WTS=WTS+WTI(I)    
          WTN=WTN+WTI(I)*I  
  120     IMAX=I    
  130     VINT(132)=XNAVE   
          VINT(133)=WTN/WTS 
          VINT(134)=WTS 
        ENDIF   
    
C...Pick multiplicity of overlayed events.  
      ELSE  
        IF(MSTP(133).EQ.0) THEN 
          MINT(81)=MAX(1,MSTP(134)) 
        ELSE    
          WTR=WTS*RLU(0)    
          DO 140 I=1,IMAX   
          MINT(81)=I    
          WTR=WTR-WTI(I)    
          IF(WTR.LE.0.) GOTO 150    
  140     CONTINUE  
  150     CONTINUE  
        ENDIF   
      ENDIF 
    
C...Format statement for error message. 
 1000 FORMAT(1X,'Warning: requested average number of events per bunch',    
     &'crossing too large, ',1P,E12.4)  
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYRAND 
    
C...Generates quantities characterizing the high-pT scattering at the   
C...parton level according to the matrix elements. Chooses incoming,    
C...reacting partons, their momentum fractions and one of the possible  
C...subprocesses.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
    
C...Initial values, specifically for (first) semihard interaction.  
      MINT(17)=0    
      MINT(18)=0    
      VINT(143)=1.  
      VINT(144)=1.  
      IF(MSUB(95).EQ.1.OR.MINT(82).GE.2) CALL PYMULT(2) 
      ISUB=0    
  100 MINT(51)=0    
    
C...Choice of process type - first event of overlay.    
      IF(MINT(82).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96)) THEN 
        RSUB=XSEC(0,1)*RLU(0)   
        DO 110 I=1,200  
        IF(MSUB(I).NE.1) GOTO 110   
        ISUB=I  
        RSUB=RSUB-XSEC(I,1) 
        IF(RSUB.LE.0.) GOTO 120 
  110   CONTINUE    
  120   IF(ISUB.EQ.95) ISUB=96  
    
C...Choice of inclusive process type - overlayed events.    
      ELSEIF(MINT(82).GE.2.AND.ISUB.EQ.0) THEN  
        RSUB=VINT(131)*RLU(0)   
        ISUB=96 
        IF(RSUB.GT.VINT(106)) ISUB=93   
        IF(RSUB.GT.VINT(106)+VINT(104)) ISUB=92 
        IF(RSUB.GT.VINT(106)+VINT(104)+VINT(103)) ISUB=91   
      ENDIF 
      IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)+1   
      IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)+1 
      MINT(1)=ISUB  
    
C...Find resonances (explicit or implicit in cross-section).    
      MINT(72)=0    
      KFR1=0    
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        KFR1=KFPR(ISUB,1)   
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN    
        KFR1=25 
      ENDIF 
      IF(KFR1.NE.0) THEN    
        TAUR1=PMAS(KFR1,1)**2/VINT(2)   
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2) 
        MINT(72)=1  
        MINT(73)=KFR1   
        VINT(73)=TAUR1  
        VINT(74)=GAMR1  
      ENDIF 
      IF(ISUB.EQ.141) THEN  
        KFR2=23 
        TAUR2=PMAS(KFR2,1)**2/VINT(2)   
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2) 
        MINT(72)=2  
        MINT(74)=KFR2   
        VINT(75)=TAUR2  
        VINT(76)=GAMR2  
      ENDIF 
    
C...Find product masses and minimum pT of process,  
C...optionally with broadening according to a truncated Breit-Wigner.   
      VINT(63)=0.   
      VINT(64)=0.   
      MINT(71)=0    
      VINT(71)=CKIN(3)  
      IF(MINT(82).GE.2) VINT(71)=0. 
      IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        DO 130 I=1,2    
        IF(KFPR(ISUB,I).EQ.0) THEN  
        ELSEIF(MSTP(42).LE.0) THEN  
          VINT(62+I)=PMAS(KFPR(ISUB,I),1)**2    
        ELSE    
          VINT(62+I)=ULMASS(KFPR(ISUB,I))**2    
        ENDIF   
  130   CONTINUE    
        IF(MIN(VINT(63),VINT(64)).LT.CKIN(6)**2) MINT(71)=1 
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5)) 
      ENDIF 
    
      IF(ISET(ISUB).EQ.0) THEN  
C...Double or single diffractive, or elastic scattering:    
C...choose m^2 according to 1/m^2 (diffractive), constant (elastic) 
        IS=INT(1.5+RLU(0))  
        VINT(63)=VINT(3)**2 
        VINT(64)=VINT(4)**2 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) VINT(62+IS)=PARP(111)**2   
        IF(ISUB.EQ.93) VINT(65-IS)=PARP(111)**2 
        SH=VINT(2)  
        SQM1=VINT(3)**2 
        SQM2=VINT(4)**2 
        SQM3=VINT(63)   
        SQM4=VINT(64)   
        SQLA12=(SH-SQM1-SQM2)**2-4.*SQM1*SQM2   
        SQLA34=(SH-SQM3-SQM4)**2-4.*SQM3*SQM4   
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH    
        THTER2=SQRT(MAX(0.,SQLA12))*SQRT(MAX(0.,SQLA34))/SH 
        THL=0.5*(THTER1-THTER2) 
        THU=0.5*(THTER1+THTER2) 
        THM=MIN(MAX(THL,PARP(101)),THU) 
        JTMAX=0 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) JTMAX=ISUB-91  
        DO 140 JT=1,JTMAX   
        MINT(13+3*JT-IS*(2*JT-3))=1 
        SQMMIN=VINT(59+3*JT-IS*(2*JT-3))    
        SQMI=VINT(8-3*JT+IS*(2*JT-3))**2    
        SQMJ=VINT(3*JT-1-IS*(2*JT-3))**2    
        SQMF=VINT(68-3*JT+IS*(2*JT-3))  
        SQUA=0.5*SH/SQMI*((1.+(SQMI-SQMJ)/SH)*THM+SQMI-SQMF-    
     &  SQMJ**2/SH+(SQMI+SQMJ)*SQMF/SH+(SQMI-SQMJ)**2/SH**2*SQMF)   
        QUAR=SH/SQMI*(THM*(THM+SH-SQMI-SQMJ-SQMF*(1.-(SQMI-SQMJ)/SH))+  
     &  SQMI*SQMJ-SQMJ*SQMF*(1.+(SQMI-SQMJ-SQMF)/SH))   
        SQMMAX=SQUA+SQRT(MAX(0.,SQUA**2-QUAR))  
        IF(ABS(QUAR/SQUA**2).LT.1.E-06) SQMMAX=0.5*QUAR/SQUA    
        SQMMAX=MIN(SQMMAX,(VINT(1)-SQRT(SQMF))**2)  
        VINT(59+3*JT-IS*(2*JT-3))=SQMMIN*(SQMMAX/SQMMIN)**RLU(0)    
  140   CONTINUE    
C...Choose t-hat according to exp(B*t-hat+C*t-hat^2).   
        SQM3=VINT(63)   
        SQM4=VINT(64)   
        SQLA34=(SH-SQM3-SQM4)**2-4.*SQM3*SQM4   
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH    
        THTER2=SQRT(MAX(0.,SQLA12))*SQRT(MAX(0.,SQLA34))/SH 
        THL=0.5*(THTER1-THTER2) 
        THU=0.5*(THTER1+THTER2) 
        B=VINT(121) 
        C=VINT(122) 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) THEN   
          B=0.5*B   
          C=0.5*C   
        ENDIF   
        THM=MIN(MAX(THL,PARP(101)),THU) 
        EXPTH=0.    
        THARG=B*(THM-THU)   
        IF(THARG.GT.-20.) EXPTH=EXP(THARG)  
  150   TH=THU+LOG(EXPTH+(1.-EXPTH)*RLU(0))/B   
        TH=MAX(THM,MIN(THU,TH)) 
        RATLOG=MIN((B+C*(TH+THM))*(TH-THM),(B+C*(TH+THU))*(TH-THU)) 
        IF(RATLOG.LT.LOG(RLU(0))) GOTO 150  
        VINT(21)=1. 
        VINT(22)=0. 
        VINT(23)=MIN(1.,MAX(-1.,(2.*TH-THTER1)/THTER2)) 
    
C...Note: in the following, by In is meant the integral over the    
C...quantity multiplying coefficient cn.    
C...Choose tau according to h1(tau)/tau, where  
C...h1(tau) = c0 + I0/I1*c1*1/tau + I0/I2*c2*1/(tau+tau_R) +    
C...I0/I3*c3*tau/((s*tau-m^2)^2+(m*Gamma)^2) +  
C...I0/I4*c4*1/(tau+tau_R') +   
C...I0/I5*c5*tau/((s*tau-m'^2)^2+(m'*Gamma')^2), and    
C...c0 + c1 + c2 + c3 + c4 + c5 = 1 
      ELSEIF(ISET(ISUB).GE.1.AND.ISET(ISUB).LE.4) THEN  
        CALL PYKLIM(1)  
        IF(MINT(51).NE.0) GOTO 100  
        RTAU=RLU(0) 
        MTAU=1  
        IF(RTAU.GT.COEF(ISUB,1)) MTAU=2 
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)) MTAU=3    
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)) MTAU=4   
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)) 
     &  MTAU=5  
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+ 
     &  COEF(ISUB,5)) MTAU=6    
        CALL PYKMAP(1,MTAU,RLU(0))  
    
C...2 -> 3, 4 processes:    
C...Choose tau' according to h4(tau,tau')/tau', where   
C...h4(tau,tau') = c0 + I0/I1*c1*(1 - tau/tau')^3/tau', and 
C...c0 + c1 = 1.    
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          CALL PYKLIM(4)    
          IF(MINT(51).NE.0) GOTO 100    
          RTAUP=RLU(0)  
          MTAUP=1   
          IF(RTAUP.GT.COEF(ISUB,15)) MTAUP=2    
          CALL PYKMAP(4,MTAUP,RLU(0))   
        ENDIF   
    
C...Choose y* according to h2(y*), where    
C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +    
C...I0/I3*c3*1/cosh(y*), I0 = y*max-y*min, and c1 + c2 + c3 = 1.    
        CALL PYKLIM(2)  
        IF(MINT(51).NE.0) GOTO 100  
        RYST=RLU(0) 
        MYST=1  
        IF(RYST.GT.COEF(ISUB,7)) MYST=2 
        IF(RYST.GT.COEF(ISUB,7)+COEF(ISUB,8)) MYST=3    
        CALL PYKMAP(2,MYST,RLU(0))  
    
C...2 -> 2 processes:   
C...Choose cos(theta-hat) (cth) according to h3(cth), where 
C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +    
C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,    
C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products), 
C...and c0 + c1 + c2 + c3 + c4 = 1. 
        CALL PYKLIM(3)  
        IF(MINT(51).NE.0) GOTO 100  
        IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN 
          RCTH=RLU(0)   
          MCTH=1    
          IF(RCTH.GT.COEF(ISUB,10)) MCTH=2  
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)) MCTH=3    
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)+COEF(ISUB,12)) MCTH=4  
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)+COEF(ISUB,12)+ 
     &    COEF(ISUB,13)) MCTH=5 
          CALL PYKMAP(3,MCTH,RLU(0))    
        ENDIF   
    
C...Low-pT or multiple interactions (first semihard interaction).   
      ELSEIF(ISET(ISUB).EQ.5) THEN  
        CALL PYMULT(3)  
        ISUB=MINT(1)    
      ENDIF 
    
C...Choose azimuthal angle. 
      VINT(24)=PARU(2)*RLU(0)   
    
C...Check against user cuts on kinematics at parton level.  
      MINT(51)=0    
      IF(ISUB.LE.90.OR.ISUB.GT.100) CALL PYKLIM(0)  
      IF(MINT(51).NE.0) GOTO 100    
      IF(MINT(82).EQ.1.AND.MSTP(141).GE.1) THEN 
        MCUT=0  
        IF(MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+MSUB(95).EQ.0)   
     &  CALL PYKCUT(MCUT)   
        IF(MCUT.NE.0) GOTO 100  
      ENDIF 
    
C...Calculate differential cross-section for different subprocesses.    
      CALL PYSIGH(NCHN,SIGS)    
    
C...Calculations for Monte Carlo estimate of all cross-sections.    
      IF(MINT(82).EQ.1.AND.ISUB.LE.90.OR.ISUB.GE.96) THEN   
        XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS  
      ELSEIF(MINT(82).EQ.1) THEN    
        XSEC(ISUB,2)=XSEC(ISUB,2)+XSEC(ISUB,1)  
      ENDIF 
    
C...Multiple interactions: store results of cross-section calculation.  
      IF(MINT(43).EQ.4.AND.MSTP(82).GE.3) THEN  
        VINT(153)=SIGS  
        CALL PYMULT(4)  
      ENDIF 
    
C...Weighting using estimate of maximum of differential cross-section.  
      VIOL=SIGS/XSEC(ISUB,1)    
      IF(VIOL.LT.RLU(0)) GOTO 100   
    
C...Check for possible violation of estimated maximum of differential   
C...cross-section used in weighting.    
      IF(MSTP(123).LE.0) THEN   
        IF(VIOL.GT.1.) THEN 
          WRITE(MSTU(11),1000) VIOL,NGEN(0,3)+1 
          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26) 
          STOP  
        ENDIF   
      ELSEIF(MSTP(123).EQ.1) THEN   
        IF(VIOL.GT.VINT(108)) THEN  
          VINT(108)=VIOL    
C          IF(VIOL.GT.1.) THEN   
C            WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1   
C            WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),   
C     &      VINT(26)    
C          ENDIF 
        ENDIF   
      ELSEIF(VIOL.GT.VINT(108)) THEN    
        VINT(108)=VIOL  
        IF(VIOL.GT.1.) THEN 
          XDIF=XSEC(ISUB,1)*(VIOL-1.)   
          XSEC(ISUB,1)=XSEC(ISUB,1)+XDIF    
          IF(MSUB(ISUB).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96))    
     &    XSEC(0,1)=XSEC(0,1)+XDIF  
C          WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1 
C          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26) 
C          IF(ISUB.LE.9) THEN    
C            WRITE(MSTU(11),1300) ISUB,XSEC(ISUB,1)  
C          ELSEIF(ISUB.LE.99) THEN   
C            WRITE(MSTU(11),1400) ISUB,XSEC(ISUB,1)  
C          ELSE  
C            WRITE(MSTU(11),1500) ISUB,XSEC(ISUB,1)  
C          ENDIF 
          VINT(108)=1.  
        ENDIF   
      ENDIF 
    
C...Multiple interactions: choose impact parameter. 
      VINT(148)=1.  
      IF(MINT(43).EQ.4.AND.(ISUB.LE.90.OR.ISUB.GE.96).AND.MSTP(82).GE.3)    
     &THEN  
        CALL PYMULT(5)  
        IF(VINT(150).LT.RLU(0)) GOTO 100    
      ENDIF 
      IF(MINT(82).EQ.1.AND.MSUB(95).EQ.1) THEN  
        IF(ISUB.LE.90.OR.ISUB.GE.95) NGEN(95,1)=NGEN(95,1)+1    
        IF(ISUB.LE.90.OR.ISUB.GE.96) NGEN(96,2)=NGEN(96,2)+1    
      ENDIF 
      IF(ISUB.LE.90.OR.ISUB.GE.96) MINT(31)=MINT(31)+1  
    
C...Choose flavour of reacting partons (and subprocess).    
      RSIGS=SIGS*RLU(0) 
      QT2=VINT(48)  
      RQQBAR=PARP(87)*(1.-(QT2/(QT2+(PARP(88)*PARP(82))**2))**2)    
      IF(ISUB.NE.95.AND.(ISUB.NE.96.OR.MSTP(82).LE.1.OR.    
     &RLU(0).GT.RQQBAR)) THEN   
        DO 190 ICHN=1,NCHN  
        KFL1=ISIG(ICHN,1)   
        KFL2=ISIG(ICHN,2)   
        MINT(2)=ISIG(ICHN,3)    
        RSIGS=RSIGS-SIGH(ICHN)  
        IF(RSIGS.LE.0.) GOTO 210    
  190   CONTINUE    
    
C...Multiple interactions: choose qqbar preferentially at small pT. 
      ELSEIF(ISUB.EQ.96) THEN   
        CALL PYSPLI(MINT(11),21,KFL1,KFLDUM)    
        CALL PYSPLI(MINT(12),21,KFL2,KFLDUM)    
        MINT(1)=11  
        MINT(2)=1   
        IF(KFL1.EQ.KFL2.AND.RLU(0).LT.0.5) MINT(2)=2    
    
C...Low-pT: choose string drawing configuration.    
      ELSE  
        KFL1=21 
        KFL2=21 
        RSIGS=6.*RLU(0) 
        MINT(2)=1   
        IF(RSIGS.GT.1.) MINT(2)=2   
        IF(RSIGS.GT.2.) MINT(2)=3   
      ENDIF 
    
C...Reassign QCD process. Partons before initial state radiation.   
  210 IF(MINT(2).GT.10) THEN    
        MINT(1)=MINT(2)/10  
        MINT(2)=MOD(MINT(2),10) 
      ENDIF 
      MINT(15)=KFL1 
      MINT(16)=KFL2 
      MINT(13)=MINT(15) 
      MINT(14)=MINT(16) 
      VINT(141)=VINT(41)    
      VINT(142)=VINT(42)    
    
C...Format statements for differential cross-section maximum violations.    
 1000 FORMAT(1X,'Error: maximum violated by',1P,E11.3,1X,   
     &'in event',1X,I7,'.'/1X,'Execution stopped!') 
 1100 FORMAT(1X,'ISUB = ',I3,'; Point of violation:'/1X,'tau =',1P, 
     &E11.3,', y* =',E11.3,', cthe = ',0P,F11.7,', tau'' =',1P,E11.3)   
 1200 FORMAT(1X,'Warning: maximum violated by',1P,E11.3,1X, 
     &'in event',1X,I7) 
 1300 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,E11.3) 
 1400 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,E11.3) 
 1500 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,E11.3) 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYSCAT 
    
C...Finds outgoing flavours and event type; sets up the kinematics  
C...and colour flow of the hard scattering. 
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      DIMENSION WDTP(0:40),WDTE(0:40,0:5),PMQ(2),Z(2),CTHE(2),PHI(2)    
    
C...Choice of subprocess, number of documentation lines.    
      ISUB=MINT(1)  
      IDOC=6+ISET(ISUB) 
      IF(ISUB.EQ.95) IDOC=8 
      MINT(3)=IDOC-6    
      IF(IDOC.GE.9) IDOC=IDOC+2 
      MINT(4)=IDOC  
      IPU1=MINT(84)+1   
      IPU2=MINT(84)+2   
      IPU3=MINT(84)+3   
      IPU4=MINT(84)+4   
      IPU5=MINT(84)+5   
      IPU6=MINT(84)+6   
    
C...Reset K, P and V vectors. Store incoming particles. 
      DO 100 JT=1,MSTP(126)+10  
      I=MINT(83)+JT 
      DO 100 J=1,5  
      K(I,J)=0  
      P(I,J)=0. 
  100 V(I,J)=0. 
      DO 110 JT=1,2 
      I=MINT(83)+JT 
      K(I,1)=21 
      K(I,2)=MINT(10+JT)    
      P(I,1)=0. 
      P(I,2)=0. 
      P(I,5)=VINT(2+JT) 
      P(I,3)=VINT(5)*(-1)**(JT+1)   
  110 P(I,4)=SQRT(P(I,3)**2+P(I,5)**2)  
      MINT(6)=2 
      KFRES=0   
    
C...Store incoming partons in their CM-frame.   
      SH=VINT(44)   
      SHR=SQRT(SH)  
      SHP=VINT(26)*VINT(2)  
      SHPR=SQRT(SHP)    
      SHUSER=SHR    
      IF(ISET(ISUB).GE.3) SHUSER=SHPR   
      DO 120 JT=1,2 
      I=MINT(84)+JT 
      K(I,1)=14 
      K(I,2)=MINT(14+JT)    
      K(I,3)=MINT(83)+2+JT  
  120 P(I,5)=ULMASS(K(I,2)) 
      IF(P(IPU1,5)+P(IPU2,5).GE.SHUSER) THEN    
        P(IPU1,5)=0.    
        P(IPU2,5)=0.    
      ENDIF 
      P(IPU1,4)=0.5*(SHUSER+(P(IPU1,5)**2-P(IPU2,5)**2)/SHUSER) 
      P(IPU1,3)=SQRT(MAX(0.,P(IPU1,4)**2-P(IPU1,5)**2)) 
      P(IPU2,4)=SHUSER-P(IPU1,4)    
      P(IPU2,3)=-P(IPU1,3)  
    
C...Copy incoming partons to documentation lines.   
      DO 130 JT=1,2 
      I1=MINT(83)+4+JT  
      I2=MINT(84)+JT    
      K(I1,1)=21    
      K(I1,2)=K(I2,2)   
      K(I1,3)=I1-2  
      DO 130 J=1,5  
  130 P(I1,J)=P(I2,J)   
    
C...Choose new quark flavour for relevant annihilation graphs.  
      IF(ISUB.EQ.12.OR.ISUB.EQ.53) THEN 
        CALL PYWIDT(21,SHR,WDTP,WDTE)   
        RKFL=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*RLU(0) 
        DO 140 I=1,2*MSTP(1)    
        KFLQ=I  
        RKFL=RKFL-(WDTE(I,1)+WDTE(I,2)+WDTE(I,4))   
        IF(RKFL.LE.0.) GOTO 150 
  140   CONTINUE    
  150   CONTINUE    
      ENDIF 
    
C...Final state flavours and colour flow: default values.   
      JS=1  
      MINT(21)=MINT(15) 
      MINT(22)=MINT(16) 
      MINT(23)=0    
      MINT(24)=0    
      KCC=20    
      KCS=ISIGN(1,MINT(15)) 
    
      IF(ISUB.LE.10) THEN   
      IF(ISUB.EQ.1) THEN    
C...f + fb -> gamma*/Z0.    
        KFRES=23    
    
      ELSEIF(ISUB.EQ.2) THEN    
C...f + fb' -> W+/- .   
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))   
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))   
        KFRES=ISIGN(24,KCH1+KCH2)   
    
      ELSEIF(ISUB.EQ.3) THEN    
C...f + fb -> H0.   
        KFRES=25    
    
      ELSEIF(ISUB.EQ.4) THEN    
C...gamma + W+/- -> W+/-.   
    
      ELSEIF(ISUB.EQ.5) THEN    
C...Z0 + Z0 -> H0.  
        XH=SH/SHP   
        MINT(21)=MINT(15)   
        MINT(22)=MINT(16)   
        PMQ(1)=ULMASS(MINT(21)) 
        PMQ(2)=ULMASS(MINT(22)) 
  240   JT=INT(1.5+RLU(0))  
        ZMIN=2.*PMQ(JT)/SHPR    
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))  
        ZMAX=MIN(1.-XH,ZMAX)    
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)   
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.  
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 240 
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 240  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP) 
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1   
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))  
        Z(3-JT)=1.-XH/(1.-Z(JT))    
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 240  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP) 
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1 
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))  
        PHIR=PARU(2)*RLU(0) 
        CPHI=COS(PHIR)  
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI    
        Z1=2.-Z(JT) 
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP) 
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP    
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*    
     &  PMQ(3-JT)**2/SHP))  
        ZMIN=2.*PMQ(3-JT)/SHPR  
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))    
        ZMAX=MIN(1.-XH,ZMAX)    
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 240 
        KCC=22  
        KFRES=25    
    
      ELSEIF(ISUB.EQ.6) THEN    
C...Z0 + W+/- -> W+/-.  
    
      ELSEIF(ISUB.EQ.7) THEN    
C...W+ + W- -> Z0.  
    
      ELSEIF(ISUB.EQ.8) THEN    
C...W+ + W- -> H0.  
        XH=SH/SHP   
  250   DO 280 JT=1,2   
        I=MINT(14+JT)   
        IA=IABS(I)  
        IF(IA.LE.10) THEN   
          RVCKM=VINT(180+I)*RLU(0)  
          DO 270 J=1,MSTP(1)    
          IB=2*J-1+MOD(IA,2)    
          IPM=(5-ISIGN(1,I))/2  
          IDC=J+MDCY(IA,2)+2    
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 270  
          MINT(20+JT)=ISIGN(IB,I)   
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)   
          IF(RVCKM.LE.0.) GOTO 280  
  270     CONTINUE  
        ELSE    
          IB=2*((IA+1)/2)-1+MOD(IA,2)   
          MINT(20+JT)=ISIGN(IB,I)   
        ENDIF   
  280   PMQ(JT)=ULMASS(MINT(20+JT)) 
        JT=INT(1.5+RLU(0))  
        ZMIN=2.*PMQ(JT)/SHPR    
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))  
        ZMAX=MIN(1.-XH,ZMAX)    
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)   
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.  
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 250 
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 250  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(24,1)**2-PMQ(JT)**2)/(Z(JT)*SHP) 
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1   
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))  
        Z(3-JT)=1.-XH/(1.-Z(JT))    
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 250  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(24,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP) 
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1 
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))  
        PHIR=PARU(2)*RLU(0) 
        CPHI=COS(PHIR)  
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI    
        Z1=2.-Z(JT) 
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP) 
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP    
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*    
     &  PMQ(3-JT)**2/SHP))  
        ZMIN=2.*PMQ(3-JT)/SHPR  
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))    
        ZMAX=MIN(1.-XH,ZMAX)    
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 250 
        KCC=22  
        KFRES=25    
      ENDIF 
    
      ELSEIF(ISUB.LE.20) THEN   
      IF(ISUB.EQ.11) THEN   
C...f + f' -> f + f'; th = (p(f)-p(f))**2.  
        KCC=MINT(2) 
        IF(MINT(15)*MINT(16).LT.0) KCC=KCC+2    
    
      ELSEIF(ISUB.EQ.12) THEN   
C...f + fb -> f' + fb'; th = (p(f)-p(f'))**2.   
        MINT(21)=ISIGN(KFLQ,MINT(15))   
        MINT(22)=-MINT(21)  
        KCC=4   
    
      ELSEIF(ISUB.EQ.13) THEN   
C...f + fb -> g + g; th arbitrary.  
        MINT(21)=21 
        MINT(22)=21 
        KCC=MINT(2)+4   
    
      ELSEIF(ISUB.EQ.14) THEN   
C...f + fb -> g + gam; th arbitrary.    
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=21  
        MINT(23-JS)=22  
        KCC=17+JS   
    
      ELSEIF(ISUB.EQ.15) THEN   
C...f + fb -> g + Z0; th arbitrary. 
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=21  
        MINT(23-JS)=23  
        KCC=17+JS   
    
      ELSEIF(ISUB.EQ.16) THEN   
C...f + fb' -> g + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2. 
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))   
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))   
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2  
        MINT(20+JS)=21  
        MINT(23-JS)=ISIGN(24,KCH1+KCH2) 
        KCC=17+JS   
    
      ELSEIF(ISUB.EQ.17) THEN   
C...f + fb -> g + H0; th arbitrary. 
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=21  
        MINT(23-JS)=25  
        KCC=17+JS   
    
      ELSEIF(ISUB.EQ.18) THEN   
C...f + fb -> gamma + gamma; th arbitrary.  
        MINT(21)=22 
        MINT(22)=22 
    
      ELSEIF(ISUB.EQ.19) THEN   
C...f + fb -> gamma + Z0; th arbitrary. 
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=22  
        MINT(23-JS)=23  
    
      ELSEIF(ISUB.EQ.20) THEN   
C...f + fb' -> gamma + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2. 
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))   
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))   
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2  
        MINT(20+JS)=22  
        MINT(23-JS)=ISIGN(24,KCH1+KCH2) 
      ENDIF 
    
      ELSEIF(ISUB.LE.30) THEN   
      IF(ISUB.EQ.21) THEN   
C...f + fb -> gamma + H0; th arbitrary. 
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=22  
        MINT(23-JS)=25  
    
      ELSEIF(ISUB.EQ.22) THEN   
C...f + fb -> Z0 + Z0; th arbitrary.    
        MINT(21)=23 
        MINT(22)=23 
    
      ELSEIF(ISUB.EQ.23) THEN   
C...f + fb' -> Z0 + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.    
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))   
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))   
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2  
        MINT(20+JS)=23  
        MINT(23-JS)=ISIGN(24,KCH1+KCH2) 
    
      ELSEIF(ISUB.EQ.24) THEN   
C...f + fb -> Z0 + H0; th arbitrary.    
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=23  
        MINT(23-JS)=25  
    
      ELSEIF(ISUB.EQ.25) THEN   
C...f + fb -> W+ + W-; th = (p(f)-p(W-))**2.    
        MINT(21)=-ISIGN(24,MINT(15))    
        MINT(22)=-MINT(21)  
    
      ELSEIF(ISUB.EQ.26) THEN   
C...f + fb' -> W+/- + H0; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.    
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))   
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))   
        IF(MINT(15)*(KCH1+KCH2).GT.0) JS=2  
        MINT(20+JS)=ISIGN(24,KCH1+KCH2) 
        MINT(23-JS)=25  
    
      ELSEIF(ISUB.EQ.27) THEN   
C...f + fb -> H0 + H0.  
    
      ELSEIF(ISUB.EQ.28) THEN   
C...f + g -> f + g; th = (p(f)-p(f))**2.    
        KCC=MINT(2)+6   
        IF(MINT(15).EQ.21) KCC=KCC+2    
        IF(MINT(15).NE.21) KCS=ISIGN(1,MINT(15))    
        IF(MINT(16).NE.21) KCS=ISIGN(1,MINT(16))    
    
      ELSEIF(ISUB.EQ.29) THEN   
C...f + g -> f + gamma; th = (p(f)-p(f))**2.    
        IF(MINT(15).EQ.21) JS=2 
        MINT(23-JS)=22  
        KCC=15+JS   
        KCS=ISIGN(1,MINT(14+JS))    
    
      ELSEIF(ISUB.EQ.30) THEN   
C...f + g -> f + Z0; th = (p(f)-p(f))**2.   
        IF(MINT(15).EQ.21) JS=2 
        MINT(23-JS)=23  
        KCC=15+JS   
        KCS=ISIGN(1,MINT(14+JS))    
      ENDIF 
    
      ELSEIF(ISUB.LE.40) THEN   
      IF(ISUB.EQ.31) THEN   
C...f + g -> f' + W+/-; th = (p(f)-p(f'))**2; choose flavour f'.    
        IF(MINT(15).EQ.21) JS=2 
        I=MINT(14+JS)   
        IA=IABS(I)  
        MINT(23-JS)=ISIGN(24,KCHG(IA,1)*I)  
        RVCKM=VINT(180+I)*RLU(0)    
        DO 220 J=1,MSTP(1)  
        IB=2*J-1+MOD(IA,2)  
        IPM=(5-ISIGN(1,I))/2    
        IDC=J+MDCY(IA,2)+2  
        IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 220    
        MINT(20+JS)=ISIGN(IB,I) 
        RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2) 
        IF(RVCKM.LE.0.) GOTO 230    
  220   CONTINUE    
  230   KCC=15+JS   
        KCS=ISIGN(1,MINT(14+JS))    
    
      ELSEIF(ISUB.EQ.32) THEN   
C...f + g -> f + H0; th = (p(f)-p(f))**2.   
        IF(MINT(15).EQ.21) JS=2 
        MINT(23-JS)=25  
        KCC=15+JS   
        KCS=ISIGN(1,MINT(14+JS))    
    
      ELSEIF(ISUB.EQ.33) THEN   
C...f + gamma -> f + g. 
    
      ELSEIF(ISUB.EQ.34) THEN   
C...f + gamma -> f + gamma. 
    
      ELSEIF(ISUB.EQ.35) THEN   
C...f + gamma -> f + Z0.    
    
      ELSEIF(ISUB.EQ.36) THEN   
C...f + gamma -> f' + W+/-. 
    
      ELSEIF(ISUB.EQ.37) THEN   
C...f + gamma -> f + H0.    
    
      ELSEIF(ISUB.EQ.38) THEN   
C...f + Z0 -> f + g.    
    
      ELSEIF(ISUB.EQ.39) THEN   
C...f + Z0 -> f + gamma.    
    
      ELSEIF(ISUB.EQ.40) THEN   
C...f + Z0 -> f + Z0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.50) THEN   
      IF(ISUB.EQ.41) THEN   
C...f + Z0 -> f' + W+/-.    
    
      ELSEIF(ISUB.EQ.42) THEN   
C...f + Z0 -> f + H0.   
    
      ELSEIF(ISUB.EQ.43) THEN   
C...f + W+/- -> f' + g. 
    
      ELSEIF(ISUB.EQ.44) THEN   
C...f + W+/- -> f' + gamma. 
    
      ELSEIF(ISUB.EQ.45) THEN   
C...f + W+/- -> f' + Z0.    
    
      ELSEIF(ISUB.EQ.46) THEN   
C...f + W+/- -> f' + W+/-.  
    
      ELSEIF(ISUB.EQ.47) THEN   
C...f + W+/- -> f' + H0.    
    
      ELSEIF(ISUB.EQ.48) THEN   
C...f + H0 -> f + g.    
    
      ELSEIF(ISUB.EQ.49) THEN   
C...f + H0 -> f + gamma.    
    
      ELSEIF(ISUB.EQ.50) THEN   
C...f + H0 -> f + Z0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.60) THEN   
      IF(ISUB.EQ.51) THEN   
C...f + H0 -> f' + W+/-.    
    
      ELSEIF(ISUB.EQ.52) THEN   
C...f + H0 -> f + H0.   
    
      ELSEIF(ISUB.EQ.53) THEN   
C...g + g -> f + fb; th arbitrary.  
        KCS=(-1)**INT(1.5+RLU(0))   
        MINT(21)=ISIGN(KFLQ,KCS)    
        MINT(22)=-MINT(21)  
        KCC=MINT(2)+10  
    
      ELSEIF(ISUB.EQ.54) THEN   
C...g + gamma -> f + fb.    
    
      ELSEIF(ISUB.EQ.55) THEN   
C...g + Z0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.56) THEN   
C...g + W+/- -> f + fb'.    
    
      ELSEIF(ISUB.EQ.57) THEN   
C...g + H0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.58) THEN   
C...gamma + gamma -> f + fb.    
    
      ELSEIF(ISUB.EQ.59) THEN   
C...gamma + Z0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.60) THEN   
C...gamma + W+/- -> f + fb'.    
      ENDIF 
    
      ELSEIF(ISUB.LE.70) THEN   
      IF(ISUB.EQ.61) THEN   
C...gamma + H0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.62) THEN   
C...Z0 + Z0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.63) THEN   
C...Z0 + W+/- -> f + fb'.   
    
      ELSEIF(ISUB.EQ.64) THEN   
C...Z0 + H0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.65) THEN   
C...W+ + W- -> f + fb.  
    
      ELSEIF(ISUB.EQ.66) THEN   
C...W+/- + H0 -> f + fb'.   
    
      ELSEIF(ISUB.EQ.67) THEN   
C...H0 + H0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.68) THEN   
C...g + g -> g + g; th arbitrary.   
        KCC=MINT(2)+12  
        KCS=(-1)**INT(1.5+RLU(0))   
    
      ELSEIF(ISUB.EQ.69) THEN   
C...gamma + gamma -> W+ + W-.   
    
      ELSEIF(ISUB.EQ.70) THEN   
C...gamma + W+/- -> gamma + W+/-    
      ENDIF 
    
      ELSEIF(ISUB.LE.80) THEN   
      IF(ISUB.EQ.71.OR.ISUB.EQ.72) THEN 
C...Z0 + Z0 -> Z0 + Z0; Z0 + Z0 -> W+ + W-. 
        XH=SH/SHP   
        MINT(21)=MINT(15)   
        MINT(22)=MINT(16)   
        PMQ(1)=ULMASS(MINT(21)) 
        PMQ(2)=ULMASS(MINT(22)) 
  290   JT=INT(1.5+RLU(0))  
        ZMIN=2.*PMQ(JT)/SHPR    
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))  
        ZMAX=MIN(1.-XH,ZMAX)    
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)   
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.  
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 290 
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 290  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP) 
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1   
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))  
        Z(3-JT)=1.-XH/(1.-Z(JT))    
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 290  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP) 
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1 
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))  
        PHIR=PARU(2)*RLU(0) 
        CPHI=COS(PHIR)  
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI    
        Z1=2.-Z(JT) 
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP) 
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP    
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*    
     &  PMQ(3-JT)**2/SHP))  
        ZMIN=2.*PMQ(3-JT)/SHPR  
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))    
        ZMAX=MIN(1.-XH,ZMAX)    
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 290 
        KCC=22  
    
      ELSEIF(ISUB.EQ.73) THEN   
C...Z0 + W+/- -> Z0 + W+/-. 
        XH=SH/SHP   
  300   JT=INT(1.5+RLU(0))  
        I=MINT(14+JT)   
        IA=IABS(I)  
        IF(IA.LE.10) THEN   
          RVCKM=VINT(180+I)*RLU(0)  
          DO 320 J=1,MSTP(1)    
          IB=2*J-1+MOD(IA,2)    
          IPM=(5-ISIGN(1,I))/2  
          IDC=J+MDCY(IA,2)+2    
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 320  
          MINT(20+JT)=ISIGN(IB,I)   
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)   
          IF(RVCKM.LE.0.) GOTO 330  
  320     CONTINUE  
        ELSE    
          IB=2*((IA+1)/2)-1+MOD(IA,2)   
          MINT(20+JT)=ISIGN(IB,I)   
        ENDIF   
  330   PMQ(JT)=ULMASS(MINT(20+JT)) 
        MINT(23-JT)=MINT(17-JT) 
        PMQ(3-JT)=ULMASS(MINT(23-JT))   
        JT=INT(1.5+RLU(0))  
        ZMIN=2.*PMQ(JT)/SHPR    
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))  
        ZMAX=MIN(1.-XH,ZMAX)    
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)   
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.  
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 300 
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 300  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP) 
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1   
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))  
        Z(3-JT)=1.-XH/(1.-Z(JT))    
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 300  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP) 
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1 
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))  
        PHIR=PARU(2)*RLU(0) 
        CPHI=COS(PHIR)  
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI    
        Z1=2.-Z(JT) 
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP) 
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP    
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*    
     &  PMQ(3-JT)**2/SHP))  
        ZMIN=2.*PMQ(3-JT)/SHPR  
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))    
        ZMAX=MIN(1.-XH,ZMAX)    
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 300 
        KCC=22  
    
      ELSEIF(ISUB.EQ.74) THEN   
C...Z0 + H0 -> Z0 + H0. 
    
      ELSEIF(ISUB.EQ.75) THEN   
C...W+ + W- -> gamma + gamma.   
    
      ELSEIF(ISUB.EQ.76.OR.ISUB.EQ.77) THEN 
C...W+ + W- -> Z0 + Z0; W+ + W- -> W+ + W-. 
        XH=SH/SHP   
  340   DO 370 JT=1,2   
        I=MINT(14+JT)   
        IA=IABS(I)  
        IF(IA.LE.10) THEN   
          RVCKM=VINT(180+I)*RLU(0)  
          DO 360 J=1,MSTP(1)    
          IB=2*J-1+MOD(IA,2)    
          IPM=(5-ISIGN(1,I))/2  
          IDC=J+MDCY(IA,2)+2    
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 360  
          MINT(20+JT)=ISIGN(IB,I)   
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)   
          IF(RVCKM.LE.0.) GOTO 370  
  360     CONTINUE  
        ELSE    
          IB=2*((IA+1)/2)-1+MOD(IA,2)   
          MINT(20+JT)=ISIGN(IB,I)   
        ENDIF   
  370   PMQ(JT)=ULMASS(MINT(20+JT)) 
        JT=INT(1.5+RLU(0))  
        ZMIN=2.*PMQ(JT)/SHPR    
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))  
        ZMAX=MIN(1.-XH,ZMAX)    
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)   
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.  
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 340 
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 340  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(24,1)**2-PMQ(JT)**2)/(Z(JT)*SHP) 
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1   
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))  
        Z(3-JT)=1.-XH/(1.-Z(JT))    
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)    
        IF(SQC1.LT.1.E-8) GOTO 340  
        C1=SQRT(SQC1)   
        C2=1.+2.*(PMAS(24,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP) 
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1 
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))  
        PHIR=PARU(2)*RLU(0) 
        CPHI=COS(PHIR)  
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI    
        Z1=2.-Z(JT) 
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP) 
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP    
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*    
     &  PMQ(3-JT)**2/SHP))  
        ZMIN=2.*PMQ(3-JT)/SHPR  
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))    
        ZMAX=MIN(1.-XH,ZMAX)    
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 340 
        KCC=22  
    
      ELSEIF(ISUB.EQ.78) THEN   
C...W+/- + H0 -> W+/- + H0. 
    
      ELSEIF(ISUB.EQ.79) THEN   
C...H0 + H0 -> H0 + H0. 
      ENDIF 
    
      ELSEIF(ISUB.LE.90) THEN   
      IF(ISUB.EQ.81) THEN   
C...q + qb -> Q' + Qb'; th = (p(q)-p(q'))**2.   
        MINT(21)=ISIGN(MINT(46),MINT(15))   
        MINT(22)=-MINT(21)  
        KCC=4   
    
      ELSEIF(ISUB.EQ.82) THEN   
C...g + g -> Q + Qb; th arbitrary.  
        KCS=(-1)**INT(1.5+RLU(0))   
        MINT(21)=ISIGN(MINT(46),KCS)    
        MINT(22)=-MINT(21)  
        KCC=MINT(2)+10  
      ENDIF 
    
      ELSEIF(ISUB.LE.100) THEN  
      IF(ISUB.EQ.95) THEN   
C...Low-pT ( = energyless g + g -> g + g).  
        KCC=MINT(2)+12  
        KCS=(-1)**INT(1.5+RLU(0))   
    
      ELSEIF(ISUB.EQ.96) THEN   
C...Multiple interactions (should be reassigned to QCD process).    
      ENDIF 
    
      ELSEIF(ISUB.LE.110) THEN  
      IF(ISUB.EQ.101) THEN  
C...g + g -> gamma*/Z0. 
        KCC=21  
        KFRES=22    
    
      ELSEIF(ISUB.EQ.102) THEN  
C...g + g -> H0.    
        KCC=21  
        KFRES=25    
      ENDIF 
    
      ELSEIF(ISUB.LE.120) THEN  
      IF(ISUB.EQ.111) THEN  
C...f + fb -> g + H0; th arbitrary. 
        IF(RLU(0).GT.0.5) JS=2  
        MINT(20+JS)=21  
        MINT(23-JS)=25  
        KCC=17+JS   
    
      ELSEIF(ISUB.EQ.112) THEN  
C...f + g -> f + H0; th = (p(f) - p(f))**2. 
        IF(MINT(15).EQ.21) JS=2 
        MINT(23-JS)=25  
        KCC=15+JS   
        KCS=ISIGN(1,MINT(14+JS))    
    
      ELSEIF(ISUB.EQ.113) THEN  
C...g + g -> g + H0; th arbitrary.  
        IF(RLU(0).GT.0.5) JS=2  
        MINT(23-JS)=25  
        KCC=22+JS   
        KCS=(-1)**INT(1.5+RLU(0))   
    
      ELSEIF(ISUB.EQ.114) THEN  
C...g + g -> gamma + gamma; th arbitrary.   
        IF(RLU(0).GT.0.5) JS=2  
        MINT(21)=22 
        MINT(22)=22 
        KCC=21  
    
      ELSEIF(ISUB.EQ.115) THEN  
C...g + g -> gamma + Z0.    
    
      ELSEIF(ISUB.EQ.116) THEN  
C...g + g -> Z0 + Z0.   
    
      ELSEIF(ISUB.EQ.117) THEN  
C...g + g -> W+ + W-.   
      ENDIF 
    
      ELSEIF(ISUB.LE.140) THEN  
      IF(ISUB.EQ.121) THEN  
C...g + g -> f + fb + H0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.160) THEN  
      IF(ISUB.EQ.141) THEN  
C...f + fb -> gamma*/Z0/Z'0.    
        KFRES=32    
    
      ELSEIF(ISUB.EQ.142) THEN  
C...f + fb' -> H+/-.    
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))   
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))   
        KFRES=ISIGN(37,KCH1+KCH2)   
    
      ELSEIF(ISUB.EQ.143) THEN  
C...f + fb' -> R.   
        KFRES=ISIGN(40,MINT(15)+MINT(16))   
      ENDIF 
    
      ELSE  
      IF(ISUB.EQ.161) THEN  
C...g + f -> H+/- + f'; th = (p(f)-p(f))**2.    
        IF(MINT(16).EQ.21) JS=2 
        IA=IABS(MINT(17-JS))    
        MINT(20+JS)=ISIGN(37,KCHG(IA,1)*MINT(17-JS))    
        JA=IA+MOD(IA,2)-MOD(IA+1,2) 
        MINT(23-JS)=ISIGN(JA,MINT(17-JS))   
        KCC=18-JS   
        IF(MINT(15).NE.21) KCS=ISIGN(1,MINT(15))    
        IF(MINT(16).NE.21) KCS=ISIGN(1,MINT(16))    
      ENDIF 
      ENDIF 
    
      IF(IDOC.EQ.7) THEN    
C...Resonance not decaying: store colour connection indices.    
        I=MINT(83)+7    
        K(IPU3,1)=1 
        K(IPU3,2)=KFRES 
        K(IPU3,3)=I 
        P(IPU3,4)=SHUSER    
        P(IPU3,5)=SHUSER    
        K(IPU1,4)=IPU2  
        K(IPU1,5)=IPU2  
        K(IPU2,4)=IPU1  
        K(IPU2,5)=IPU1  
        K(I,1)=21   
        K(I,2)=KFRES    
        P(I,4)=SHUSER   
        P(I,5)=SHUSER   
        N=IPU3  
        MINT(21)=KFRES  
        MINT(22)=0  
    
      ELSEIF(IDOC.EQ.8) THEN    
C...2 -> 2 processes: store outgoing partons in their CM-frame. 
        DO 390 JT=1,2   
        I=MINT(84)+2+JT 
        K(I,1)=1    
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3   
        K(I,2)=MINT(20+JT)  
        K(I,3)=MINT(83)+IDOC+JT-2   
        IF(IABS(K(I,2)).LE.10.OR.K(I,2).EQ.21) THEN 
          P(I,5)=ULMASS(K(I,2)) 
        ELSE    
          P(I,5)=SQRT(VINT(63+MOD(JS+JT,2)))    
        ENDIF   
  390   CONTINUE    
        IF(P(IPU3,5)+P(IPU4,5).GE.SHR) THEN 
          KFA1=IABS(MINT(21))   
          KFA2=IABS(MINT(22))   
          IF((KFA1.GT.3.AND.KFA1.NE.21).OR.(KFA2.GT.3.AND.KFA2.NE.21))  
     &    THEN  
            MINT(51)=1  
            RETURN  
          ENDIF 
          P(IPU3,5)=0.  
          P(IPU4,5)=0.  
        ENDIF   
        P(IPU3,4)=0.5*(SHR+(P(IPU3,5)**2-P(IPU4,5)**2)/SHR) 
        P(IPU3,3)=SQRT(MAX(0.,P(IPU3,4)**2-P(IPU3,5)**2))   
        P(IPU4,4)=SHR-P(IPU3,4) 
        P(IPU4,3)=-P(IPU3,3)    
        N=IPU4  
        MINT(7)=MINT(83)+7  
        MINT(8)=MINT(83)+8  
    
C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4). 
        CALL LUDBRB(IPU3,IPU4,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)  
    
      ELSEIF(IDOC.EQ.9) THEN    
C'''2 -> 3 processes:   
    
      ELSEIF(IDOC.EQ.11) THEN   
C...Z0 + Z0 -> H0, W+ + W- -> H0: store Higgs and outgoing partons. 
        PHI(1)=PARU(2)*RLU(0)   
        PHI(2)=PHI(1)-PHIR  
        DO 400 JT=1,2   
        I=MINT(84)+2+JT 
        K(I,1)=1    
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3   
        K(I,2)=MINT(20+JT)  
        K(I,3)=MINT(83)+IDOC+JT-2   
        P(I,5)=ULMASS(K(I,2))   
        IF(0.5*SHPR*Z(JT).LE.P(I,5)) P(I,5)=0.  
        PABS=SQRT(MAX(0.,(0.5*SHPR*Z(JT))**2-P(I,5)**2))    
        PTABS=PABS*SQRT(MAX(0.,1.-CTHE(JT)**2)) 
        P(I,1)=PTABS*COS(PHI(JT))   
        P(I,2)=PTABS*SIN(PHI(JT))   
        P(I,3)=PABS*CTHE(JT)*(-1)**(JT+1)   
        P(I,4)=0.5*SHPR*Z(JT)   
        IZW=MINT(83)+6+JT   
        K(IZW,1)=21 
        K(IZW,2)=23 
        IF(ISUB.EQ.8) K(IZW,2)=ISIGN(24,LUCHGE(MINT(14+JT)))    
        K(IZW,3)=IZW-2  
        P(IZW,1)=-P(I,1)    
        P(IZW,2)=-P(I,2)    
        P(IZW,3)=(0.5*SHPR-PABS*CTHE(JT))*(-1)**(JT+1)  
        P(IZW,4)=0.5*SHPR*(1.-Z(JT))    
  400   P(IZW,5)=-SQRT(MAX(0.,P(IZW,3)**2+PTABS**2-P(IZW,4)**2))    
        I=MINT(83)+9    
        K(IPU5,1)=1 
        K(IPU5,2)=KFRES 
        K(IPU5,3)=I 
        P(IPU5,5)=SHR   
        P(IPU5,1)=-P(IPU3,1)-P(IPU4,1)  
        P(IPU5,2)=-P(IPU3,2)-P(IPU4,2)  
        P(IPU5,3)=-P(IPU3,3)-P(IPU4,3)  
        P(IPU5,4)=SHPR-P(IPU3,4)-P(IPU4,4)  
        K(I,1)=21   
        K(I,2)=KFRES    
        DO 410 J=1,5    
  410   P(I,J)=P(IPU5,J)    
        N=IPU5  
        MINT(23)=KFRES  
    
      ELSEIF(IDOC.EQ.12) THEN   
C...Z0 and W+/- scattering: store bosons and outgoing partons.  
        PHI(1)=PARU(2)*RLU(0)   
        PHI(2)=PHI(1)-PHIR  
        DO 420 JT=1,2   
        I=MINT(84)+2+JT 
        K(I,1)=1    
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3   
        K(I,2)=MINT(20+JT)  
        K(I,3)=MINT(83)+IDOC+JT-2   
        P(I,5)=ULMASS(K(I,2))   
        IF(0.5*SHPR*Z(JT).LE.P(I,5)) P(I,5)=0.  
        PABS=SQRT(MAX(0.,(0.5*SHPR*Z(JT))**2-P(I,5)**2))    
        PTABS=PABS*SQRT(MAX(0.,1.-CTHE(JT)**2)) 
        P(I,1)=PTABS*COS(PHI(JT))   
        P(I,2)=PTABS*SIN(PHI(JT))   
        P(I,3)=PABS*CTHE(JT)*(-1)**(JT+1)   
        P(I,4)=0.5*SHPR*Z(JT)   
        IZW=MINT(83)+6+JT   
        K(IZW,1)=21 
        IF(MINT(14+JT).EQ.MINT(20+JT)) THEN 
          K(IZW,2)=23   
        ELSE    
          K(IZW,2)=ISIGN(24,LUCHGE(MINT(14+JT))-LUCHGE(MINT(20+JT)))    
        ENDIF   
        K(IZW,3)=IZW-2  
        P(IZW,1)=-P(I,1)    
        P(IZW,2)=-P(I,2)    
        P(IZW,3)=(0.5*SHPR-PABS*CTHE(JT))*(-1)**(JT+1)  
        P(IZW,4)=0.5*SHPR*(1.-Z(JT))    
        P(IZW,5)=-SQRT(MAX(0.,P(IZW,3)**2+PTABS**2-P(IZW,4)**2))    
        IPU=MINT(84)+4+JT   
        K(IPU,1)=3  
        K(IPU,2)=KFPR(ISUB,JT)  
        K(IPU,3)=MINT(83)+8+JT  
        IF(IABS(K(IPU,2)).LE.10.OR.K(IPU,2).EQ.21) THEN 
          P(IPU,5)=ULMASS(K(IPU,2)) 
        ELSE    
          P(IPU,5)=SQRT(VINT(63+MOD(JS+JT,2)))  
        ENDIF   
        MINT(22+JT)=K(IZW,2)    
  420   CONTINUE    
        IF(ISUB.EQ.72) K(MINT(84)+4+INT(1.5+RLU(0)),2)=-24  
C...Find rotation and boost for hard scattering subsystem.  
        I1=MINT(83)+7   
        I2=MINT(83)+8   
        BEXCM=(P(I1,1)+P(I2,1))/(P(I1,4)+P(I2,4))   
        BEYCM=(P(I1,2)+P(I2,2))/(P(I1,4)+P(I2,4))   
        BEZCM=(P(I1,3)+P(I2,3))/(P(I1,4)+P(I2,4))   
        GAMCM=(P(I1,4)+P(I2,4))/SHR 
        BEPCM=BEXCM*P(I1,1)+BEYCM*P(I1,2)+BEZCM*P(I1,3) 
        PX=P(I1,1)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEXCM 
        PY=P(I1,2)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEYCM 
        PZ=P(I1,3)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEZCM 
        THECM=ULANGL(PZ,SQRT(PX**2+PY**2))  
        PHICM=ULANGL(PX,PY) 
C...Store hard scattering subsystem. Rotate and boost it.   
        SQLAM=(SH-P(IPU5,5)**2-P(IPU6,5)**2)**2-4.*P(IPU5,5)**2*    
     &  P(IPU6,5)**2    
        PABS=SQRT(MAX(0.,SQLAM/(4.*SH)))    
        CTHWZ=VINT(23)  
        STHWZ=SQRT(MAX(0.,1.-CTHWZ**2)) 
        PHIWZ=VINT(24)-PHICM    
        P(IPU5,1)=PABS*STHWZ*COS(PHIWZ) 
        P(IPU5,2)=PABS*STHWZ*SIN(PHIWZ) 
        P(IPU5,3)=PABS*CTHWZ    
        P(IPU5,4)=SQRT(PABS**2+P(IPU5,5)**2)    
        P(IPU6,1)=-P(IPU5,1)    
        P(IPU6,2)=-P(IPU5,2)    
        P(IPU6,3)=-P(IPU5,3)    
        P(IPU6,4)=SQRT(PABS**2+P(IPU6,5)**2)    
        CALL LUDBRB(IPU5,IPU6,THECM,PHICM,DBLE(BEXCM),DBLE(BEYCM),  
     &  DBLE(BEZCM))    
        DO 430 JT=1,2   
        I1=MINT(83)+8+JT    
        I2=MINT(84)+4+JT    
        K(I1,1)=21  
        K(I1,2)=K(I2,2) 
        DO 430 J=1,5    
  430   P(I1,J)=P(I2,J) 
        N=IPU6  
        MINT(7)=MINT(83)+9  
        MINT(8)=MINT(83)+10 
      ENDIF 
    
      IF(IDOC.GE.8) THEN    
C...Store colour connection indices.    
        DO 440 J=1,2    
        JC=J    
        IF(KCS.EQ.-1) JC=3-J    
        IF(ICOL(KCC,1,JC).NE.0.AND.K(IPU1,1).EQ.14) K(IPU1,J+3)=    
     &  K(IPU1,J+3)+MINT(84)+ICOL(KCC,1,JC) 
        IF(ICOL(KCC,2,JC).NE.0.AND.K(IPU2,1).EQ.14) K(IPU2,J+3)=    
     &  K(IPU2,J+3)+MINT(84)+ICOL(KCC,2,JC) 
        IF(ICOL(KCC,3,JC).NE.0.AND.K(IPU3,1).EQ.3) K(IPU3,J+3)= 
     &  MSTU(5)*(MINT(84)+ICOL(KCC,3,JC))   
  440   IF(ICOL(KCC,4,JC).NE.0.AND.K(IPU4,1).EQ.3) K(IPU4,J+3)= 
     &  MSTU(5)*(MINT(84)+ICOL(KCC,4,JC))   
    
C...Copy outgoing partons to documentation lines.   
        DO 450 I=1,2    
        I1=MINT(83)+IDOC-2+I    
        I2=MINT(84)+2+I 
        K(I1,1)=21  
        K(I1,2)=K(I2,2) 
        IF(IDOC.LE.9) K(I1,3)=0 
        IF(IDOC.GE.11) K(I1,3)=MINT(83)+2+I 
        DO 450 J=1,5    
  450   P(I1,J)=P(I2,J) 
      ENDIF 
      MINT(52)=N    
    
C...Low-pT events: remove gluons used for string drawing purposes.  
      IF(ISUB.EQ.95) THEN   
        K(IPU3,1)=K(IPU3,1)+10  
        K(IPU4,1)=K(IPU4,1)+10  
        DO 460 J=41,66  
  460   VINT(J)=0.  
        DO 470 I=MINT(83)+5,MINT(83)+8  
        DO 470 J=1,5    
  470   P(I,J)=0.   
      ENDIF 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYSSPA(IPU1,IPU2)  
    
C...Generates spacelike parton showers. 
      IMPLICIT DOUBLE PRECISION(D)  
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      DIMENSION KFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVS(2),ROBO(5),   
     &XFS(2,-6:6),XFA(-6:6),XFB(-6:6),XFN(-6:6),WTAP(-6:6),WTSF(-6:6),  
     &THE2(2),ALAM(2),DQ2(3),DPC(3),DPD(4),DPB(4)   
    
C...Calculate maximum virtuality and check that evolution possible. 
      IPUS1=IPU1    
      IPUS2=IPU2    
      ISUB=MINT(1)  
      Q2E=VINT(52)  
      IF(ISET(ISUB).EQ.1) THEN  
        Q2E=Q2E/PARP(67)    
      ELSEIF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN   
        Q2E=PMAS(23,1)**2   
        IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77) Q2E=PMAS(24,1)**2 
      ENDIF 
      TMAX=LOG(PARP(67)*PARP(63)*Q2E/PARP(61)**2)   
      IF(PARP(67)*Q2E.LT.MAX(PARP(62)**2,2.*PARP(61)**2).OR.    
     &TMAX.LT.0.2) RETURN   
    
C...Common constants and initial values. Save normal Lambda value.  
      XE0=2.*PARP(65)/VINT(1)   
      ALAMS=PARU(111)   
      PARU(111)=PARP(61)    
      NS=N  
  100 N=NS  
      DO 110 JT=1,2 
      KFLS(JT)=MINT(14+JT)  
      KFLS(JT+2)=KFLS(JT)   
      XS(JT)=VINT(40+JT)    
      ZS(JT)=1. 
      Q2S(JT)=PARP(67)*Q2E  
      TEVS(JT)=TMAX 
      ALAM(JT)=PARP(61) 
      THE2(JT)=100. 
      DO 110 KFL=-6,6   
  110 XFS(JT,KFL)=XSFX(JT,KFL)  
      DSH=VINT(44)  
      IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) DSH=VINT(26)*VINT(2)   
    
C...Pick up leg with highest virtuality.    
  120 N=N+1 
      JT=1  
      IF(N.GT.NS+1.AND.Q2S(2).GT.Q2S(1)) JT=2   
      KFLB=KFLS(JT) 
      XB=XS(JT) 
      DO 130 KFL=-6,6   
  130 XFB(KFL)=XFS(JT,KFL)  
      DSHR=2D0*SQRT(DSH)    
      DSHZ=DSH/DBLE(ZS(JT)) 
      XE=MAX(XE0,XB*(1./(1.-PARP(66))-1.))  
      IF(XB+XE.GE.0.999) THEN   
        Q2B=0.  
        GOTO 220    
      ENDIF 
    
C...Maximum Q2 without or with Q2 ordering. Effective Lambda and n_f.   
      IF(MSTP(62).LE.1) THEN    
        Q2B=0.5*(1./ZS(JT)+1.)*Q2S(JT)+0.5*(1./ZS(JT)-1.)*(Q2S(3-JT)-   
     &  SNGL(DSH)+SQRT((SNGL(DSH)+Q2S(1)+Q2S(2))**2+8.*Q2S(1)*Q2S(2)*   
     &  ZS(JT)/(1.-ZS(JT))))    
        TEVB=LOG(PARP(63)*Q2B/ALAM(JT)**2)  
      ELSE  
        Q2B=Q2S(JT) 
        TEVB=TEVS(JT)   
      ENDIF 
      ALSDUM=ULALPS(PARP(63)*Q2B)   
      TEVB=TEVB+2.*LOG(ALAM(JT)/PARU(117))  
      TEVBSV=TEVB   
      ALAM(JT)=PARU(117)    
      B0=(33.-2.*MSTU(118))/6.  
    
C...Calculate Altarelli-Parisi and structure function weights.  
      DO 140 KFL=-6,6   
      WTAP(KFL)=0.  
  140 WTSF(KFL)=0.  
      IF(KFLB.EQ.21) THEN   
        WTAPQ=16.*(1.-SQRT(XB+XE))/(3.*SQRT(XB))    
        DO 150 KFL=-MSTP(54),MSTP(54)   
        IF(KFL.EQ.0) WTAP(KFL)=6.*LOG((1.-XB)/XE)   
  150   IF(KFL.NE.0) WTAP(KFL)=WTAPQ    
      ELSE  
        WTAP(0)=0.5*XB*(1./(XB+XE)-1.)  
        WTAP(KFLB)=8.*LOG((1.-XB)*(XB+XE)/XE)/3.    
      ENDIF 
  160 WTSUM=0.  
      IF(KFLB.NE.21) XFBO=XFB(KFLB) 
      IF(KFLB.EQ.21) XFBO=XFB(0)
C***************************************************************
C**********ERROR HAS OCCURED HERE
      IF(XFBO.EQ.0.0) THEN
                WRITE(MSTU(11),1000)
                WRITE(MSTU(11),1001) KFLB,XFB(FKLB)
		XFBO=0.00001
      ENDIF
C****************************************************************    
      DO 170 KFL=-MSTP(54),MSTP(54) 
      WTSF(KFL)=XFB(KFL)/XFBO   
  170 WTSUM=WTSUM+WTAP(KFL)*WTSF(KFL)   
      WTSUM=MAX(0.0001,WTSUM)   
    
C...Choose new t: fix alpha_s, alpha_s(Q2), alpha_s(k_T2).  
  180 IF(MSTP(64).LE.0) THEN    
        TEVB=TEVB+LOG(RLU(0))*PARU(2)/(PARU(111)*WTSUM) 
      ELSEIF(MSTP(64).EQ.1) THEN    
        TEVB=TEVB*EXP(MAX(-100.,LOG(RLU(0))*B0/WTSUM))  
      ELSE  
        TEVB=TEVB*EXP(MAX(-100.,LOG(RLU(0))*B0/(5.*WTSUM))) 
      ENDIF 
  190 Q2REF=ALAM(JT)**2*EXP(TEVB)   
      Q2B=Q2REF/PARP(63)    
    
C...Evolution ended or select flavour for branching parton. 
      IF(Q2B.LT.PARP(62)**2) THEN   
        Q2B=0.  
      ELSE  
        WTRAN=RLU(0)*WTSUM  
        KFLA=-MSTP(54)-1    
  200   KFLA=KFLA+1 
        WTRAN=WTRAN-WTAP(KFLA)*WTSF(KFLA)   
        IF(KFLA.LT.MSTP(54).AND.WTRAN.GT.0.) GOTO 200   
        IF(KFLA.EQ.0) KFLA=21   
    
C...Choose z value and corrective weight.   
        IF(KFLB.EQ.21.AND.KFLA.EQ.21) THEN  
          Z=1./(1.+((1.-XB)/XB)*(XE/(1.-XB))**RLU(0))   
          WTZ=(1.-Z*(1.-Z))**2  
        ELSEIF(KFLB.EQ.21) THEN 
          Z=XB/(1.-RLU(0)*(1.-SQRT(XB+XE)))**2  
          WTZ=0.5*(1.+(1.-Z)**2)*SQRT(Z)    
        ELSEIF(KFLA.EQ.21) THEN 
          Z=XB*(1.+RLU(0)*(1./(XB+XE)-1.))  
          WTZ=1.-2.*Z*(1.-Z)    
        ELSE    
          Z=1.-(1.-XB)*(XE/((XB+XE)*(1.-XB)))**RLU(0)   
          WTZ=0.5*(1.+Z**2) 
        ENDIF   
    
C...Option with resummation of soft gluon emission as effective z shift.    
        IF(MSTP(65).GE.1) THEN  
          RSOFT=6.  
          IF(KFLB.NE.21) RSOFT=8./3.    
          Z=Z*(TEVB/TEVS(JT))**(RSOFT*XE/((XB+XE)*B0))  
          IF(Z.LE.XB) GOTO 180  
        ENDIF   
    
C...Option with alpha_s(k_T2)Q2): demand k_T2 > cutoff, reweight.   
        IF(MSTP(64).GE.2) THEN  
          IF((1.-Z)*Q2B.LT.PARP(62)**2) GOTO 180    
          ALPRAT=TEVB/(TEVB+LOG(1.-Z))  
          IF(ALPRAT.LT.5.*RLU(0)) GOTO 180  
          IF(ALPRAT.GT.5.) WTZ=WTZ*ALPRAT/5.    
        ENDIF   
    
C...Option with angular ordering requirement.   
        IF(MSTP(62).GE.3) THEN  
          THE2T=(4.*Z**2*Q2B)/(VINT(2)*(1.-Z)*XB**2)    
          IF(THE2T.GT.THE2(JT)) GOTO 180    
        ENDIF   
    
C...Weighting with new structure functions. 
        CALL PYSTFU(MINT(10+JT),XB,Q2REF,XFN,JT)   
        IF(KFLB.NE.21) XFBN=XFN(KFLB)   
        IF(KFLB.EQ.21) XFBN=XFN(0)  
        IF(XFBN.LT.1E-20) THEN  
          IF(KFLA.EQ.KFLB) THEN 
            TEVB=TEVBSV 
            WTAP(KFLB)=0.   
            GOTO 160    
          ELSEIF(TEVBSV-TEVB.GT.0.2) THEN   
            TEVB=0.5*(TEVBSV+TEVB)  
            GOTO 190    
          ELSE  
            XFBN=1E-10  
          ENDIF 
        ENDIF   
        DO 210 KFL=-MSTP(54),MSTP(54)   
  210   XFB(KFL)=XFN(KFL)   
        XA=XB/Z 
        CALL PYSTFU(MINT(10+JT),XA,Q2REF,XFA,JT)   
        IF(KFLA.NE.21) XFAN=XFA(KFLA)   
        IF(KFLA.EQ.21) XFAN=XFA(0)  
        IF(XFAN.LT.1E-20) GOTO 160  
        IF(KFLA.NE.21) WTSFA=WTSF(KFLA) 
        IF(KFLA.EQ.21) WTSFA=WTSF(0)    
        IF(WTZ*XFAN/XFBN.LT.RLU(0)*WTSFA) GOTO 160  
      ENDIF 
    
C...Define two hard scatterers in their CM-frame.   
  220 IF(N.EQ.NS+2) THEN    
        DQ2(JT)=Q2B 
        DPLCM=SQRT((DSH+DQ2(1)+DQ2(2))**2-4D0*DQ2(1)*DQ2(2))/DSHR   
        DO 240 JR=1,2   
        I=NS+JR 
        IF(JR.EQ.1) IPO=IPUS1   
        IF(JR.EQ.2) IPO=IPUS2   
        DO 230 J=1,5    
        K(I,J)=0    
        P(I,J)=0.   
  230   V(I,J)=0.   
        K(I,1)=14   
        K(I,2)=KFLS(JR+2)   
        K(I,4)=IPO  
        K(I,5)=IPO  
        P(I,3)=DPLCM*(-1)**(JR+1)   
        P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR 
        P(I,5)=-SQRT(SNGL(DQ2(JR))) 
        K(IPO,1)=14 
        K(IPO,3)=I  
        K(IPO,4)=MOD(K(IPO,4),MSTU(5))+MSTU(5)*I    
  240   K(IPO,5)=MOD(K(IPO,5),MSTU(5))+MSTU(5)*I    
    
C...Find maximum allowed mass of timelike parton.   
      ELSEIF(N.GT.NS+2) THEN    
        JR=3-JT 
        DQ2(3)=Q2B  
        DPC(1)=P(IS(1),4)   
        DPC(2)=P(IS(2),4)   
        DPC(3)=0.5*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))    
        DPD(1)=DSH+DQ2(JR)+DQ2(JT)  
        DPD(2)=DSHZ+DQ2(JR)+DQ2(3)  
        DPD(3)=SQRT(DPD(1)**2-4D0*DQ2(JR)*DQ2(JT))  
        DPD(4)=SQRT(DPD(2)**2-4D0*DQ2(JR)*DQ2(3))   
        IKIN=0  
        IF(Q2S(JR).GE.(0.5*PARP(62))**2.AND.DPD(1)-DPD(3).GE.   
     &  1D-10*DPD(1)) IKIN=1    
        IF(IKIN.EQ.0) DMSMA=(DQ2(JT)/DBLE(ZS(JT))-DQ2(3))*(DSH/ 
     &  (DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))    
        IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/(2.*  
     &  DQ2(JR))-DQ2(JT)-DQ2(3) 
    
C...Generate timelike parton shower (if required).  
        IT=N    
        DO 250 J=1,5    
        K(IT,J)=0   
        P(IT,J)=0.  
  250   V(IT,J)=0.  
        K(IT,1)=3   
        K(IT,2)=21  
        IF(KFLB.EQ.21.AND.KFLS(JT+2).NE.21) K(IT,2)=-KFLS(JT+2) 
        IF(KFLB.NE.21.AND.KFLS(JT+2).EQ.21) K(IT,2)=KFLB    
        P(IT,5)=ULMASS(K(IT,2)) 
        IF(SNGL(DMSMA).LE.P(IT,5)**2) GOTO 100  
        IF(MSTP(63).GE.1) THEN  
          P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR    
          P(IT,3)=SQRT(P(IT,4)**2-P(IT,5)**2)   
          IF(MSTP(63).EQ.1) THEN    
            Q2TIM=DMSMA 
          ELSEIF(MSTP(63).EQ.2) THEN    
            Q2TIM=MIN(SNGL(DMSMA),PARP(71)*Q2S(JT)) 
          ELSE  
C'''Here remains to introduce angular ordering in first branching.  
            Q2TIM=DMSMA 
          ENDIF 
          CALL LUSHOW(IT,0,SQRT(Q2TIM)) 
          IF(N.GE.IT+1) P(IT,5)=P(IT+1,5)   
        ENDIF   
    
C...Reconstruct kinematics of branching: timelike parton shower.    
        DMS=P(IT,5)**2  
        IF(IKIN.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/(DSH+DQ2(JT))  
        IF(IKIN.EQ.1) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*DPD(2)+0.5*DPD(3)*   
     &  DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/(4.*DSH*DPC(3)**2) 
        IF(DPT2.LT.0.) GOTO 100 
        DPB(1)=(0.5*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/  
     &  DSHR)/DPC(3)-DPC(3) 
        P(IT,1)=SQRT(SNGL(DPT2))    
        P(IT,3)=DPB(1)*(-1)**(JT+1) 
        P(IT,4)=(DSHZ-DSH-DMS)/DSHR 
        IF(N.GE.IT+1) THEN  
          DPB(1)=SQRT(DPB(1)**2+DPT2)   
          DPB(2)=SQRT(DPB(1)**2+DMS)    
          DPB(3)=P(IT+1,3)  
          DPB(4)=SQRT(DPB(3)**2+DMS)    
          DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)* 
     &    DPB(1))   
          CALL LUDBRB(IT+1,N,0.,0.,0D0,0D0,DBEZ)    
          THE=ULANGL(P(IT,3),P(IT,1))   
          CALL LUDBRB(IT+1,N,THE,0.,0D0,0D0,0D0)    
        ENDIF   
    
C...Reconstruct kinematics of branching: spacelike parton.  
        DO 260 J=1,5    
        K(N+1,J)=0  
        P(N+1,J)=0. 
  260   V(N+1,J)=0. 
        K(N+1,1)=14 
        K(N+1,2)=KFLB   
        P(N+1,1)=P(IT,1)    
        P(N+1,3)=P(IT,3)+P(IS(JT),3)    
        P(N+1,4)=P(IT,4)+P(IS(JT),4)    
        P(N+1,5)=-SQRT(SNGL(DQ2(3)))    
    
C...Define colour flow of branching.    
        K(IS(JT),3)=N+1 
        K(IT,3)=N+1 
        ID1=IT  
        IF((K(N+1,2).GT.0.AND.K(N+1,2).NE.21.AND.K(ID1,2).GT.0.AND. 
     &  K(ID1,2).NE.21).OR.(K(N+1,2).LT.0.AND.K(ID1,2).EQ.21).OR.   
     &  (K(N+1,2).EQ.21.AND.K(ID1,2).EQ.21.AND.RLU(0).GT.0.5).OR.   
     &  (K(N+1,2).EQ.21.AND.K(ID1,2).LT.0)) ID1=IS(JT)  
        ID2=IT+IS(JT)-ID1   
        K(N+1,4)=K(N+1,4)+ID1   
        K(N+1,5)=K(N+1,5)+ID2   
        K(ID1,4)=K(ID1,4)+MSTU(5)*(N+1) 
        K(ID1,5)=K(ID1,5)+MSTU(5)*ID2   
        K(ID2,4)=K(ID2,4)+MSTU(5)*ID1   
        K(ID2,5)=K(ID2,5)+MSTU(5)*(N+1) 
        N=N+1   
    
C...Boost to new CM-frame.  
        CALL LUDBRB(NS+1,N,0.,0.,-DBLE((P(N,1)+P(IS(JR),1))/(P(N,4)+    
     &  P(IS(JR),4))),0D0,-DBLE((P(N,3)+P(IS(JR),3))/(P(N,4)+   
     &  P(IS(JR),4))))  
        IR=N+(JT-1)*(IS(1)-N)   
        CALL LUDBRB(NS+1,N,-ULANGL(P(IR,3),P(IR,1)),PARU(2)*RLU(0), 
     &  0D0,0D0,0D0)    
      ENDIF 
    
C...Save quantities, loop back. 
      IS(JT)=N  
      Q2S(JT)=Q2B   
      DQ2(JT)=Q2B   
      IF(MSTP(62).GE.3) THE2(JT)=THE2T  
      DSH=DSHZ  
      IF(Q2B.GE.(0.5*PARP(62))**2) THEN 
        KFLS(JT+2)=KFLS(JT) 
        KFLS(JT)=KFLA   
        XS(JT)=XA   
        ZS(JT)=Z    
        DO 270 KFL=-6,6 
  270   XFS(JT,KFL)=XFA(KFL)    
        TEVS(JT)=TEVB   
      ELSE  
        IF(JT.EQ.1) IPU1=N  
        IF(JT.EQ.2) IPU2=N  
      ENDIF 
      IF(N.GT.MSTU(4)-MSTU(32)-10) THEN 
        CALL LUERRM(11,'(PYSSPA:) no more memory left in LUJETS')   
        IF(MSTU(21).GE.1) N=NS  
        IF(MSTU(21).GE.1) RETURN    
      ENDIF 
      IF(MAX(Q2S(1),Q2S(2)).GE.(0.5*PARP(62))**2.OR.N.LE.NS+1) GOTO 120 
    
C...Boost hard scattering partons to frame of shower initiators.    
      DO 280 J=1,3  
  280 ROBO(J+2)=(P(NS+1,J)+P(NS+2,J))/(P(NS+1,4)+P(NS+2,4)) 
      DO 290 J=1,5  
  290 P(N+2,J)=P(NS+1,J)    
      ROBOT=ROBO(3)**2+ROBO(4)**2+ROBO(5)**2    
      IF(ROBOT.GE.0.999999) THEN    
        ROBOT=1.00001*SQRT(ROBOT)   
        ROBO(3)=ROBO(3)/ROBOT   
        ROBO(4)=ROBO(4)/ROBOT   
        ROBO(5)=ROBO(5)/ROBOT   
      ENDIF 
      CALL LUDBRB(N+2,N+2,0.,0.,-DBLE(ROBO(3)),-DBLE(ROBO(4)),  
     &-DBLE(ROBO(5)))   
      ROBO(2)=ULANGL(P(N+2,1),P(N+2,2)) 
      ROBO(1)=ULANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))    
      CALL LUDBRB(MINT(83)+5,NS,ROBO(1),ROBO(2),DBLE(ROBO(3)),  
     &DBLE(ROBO(4)),DBLE(ROBO(5)))  
    
C...Store user information. Reset Lambda value. 
      K(IPU1,3)=MINT(83)+3  
      K(IPU2,3)=MINT(83)+4  
      DO 300 JT=1,2 
      MINT(12+JT)=KFLS(JT)  
  300 VINT(140+JT)=XS(JT)   
      PARU(111)=ALAMS   
 1000 FORMAT(5X,'structure function has a zero point here')
 1001 FORMAT(5X,'xf(x,i=',I5,')=',F10.5)

      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYMULT(MMUL)   
    
C...Initializes treatment of multiple interactions, selects kinematics  
C...of hardest interaction if low-pT physics included in run, and   
C...generates all non-hardest interactions. 
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      DIMENSION NMUL(20),SIGM(20),KSTR(500,2)   
      SAVE XT2,XT2FAC,XC2,XTS,IRBIN,RBIN,NMUL,SIGM  
    
C...Initialization of multiple interaction treatment.   
      IF(MMUL.EQ.1) THEN    
        IF(MSTP(122).GE.1) WRITE(MSTU(11),1000) MSTP(82)    
        ISUB=96 
        MINT(1)=96  
        VINT(63)=0. 
        VINT(64)=0. 
        VINT(143)=1.    
        VINT(144)=1.    
    
C...Loop over phase space points: xT2 choice in 20 bins.    
  100   SIGSUM=0.   
        DO 120 IXT2=1,20    
        NMUL(IXT2)=MSTP(83) 
        SIGM(IXT2)=0.   
        DO 110 ITRY=1,MSTP(83)  
        RSCA=0.05*((21-IXT2)-RLU(0))    
        XT2=VINT(149)*(1.+VINT(149))/(VINT(149)+RSCA)-VINT(149) 
        XT2=MAX(0.01*VINT(149),XT2) 
        VINT(25)=XT2    
    
C...Choose tau and y*. Calculate cos(theta-hat).    
        IF(RLU(0).LE.COEF(ISUB,1)) THEN 
          TAUP=(2.*(1.+SQRT(1.-XT2))/XT2-1.)**RLU(0)    
          TAU=XT2*(1.+TAUP)**2/(4.*TAUP)    
        ELSE    
          TAU=XT2*(1.+TAN(RLU(0)*ATAN(SQRT(1./XT2-1.)))**2) 
        ENDIF   
        VINT(21)=TAU    
        CALL PYKLIM(2)  
        RYST=RLU(0) 
        MYST=1  
        IF(RYST.GT.COEF(ISUB,7)) MYST=2 
        IF(RYST.GT.COEF(ISUB,7)+COEF(ISUB,8)) MYST=3    
        CALL PYKMAP(2,MYST,RLU(0))  
        VINT(23)=SQRT(MAX(0.,1.-XT2/TAU))*(-1)**INT(1.5+RLU(0)) 
    
C...Calculate differential cross-section.   
        VINT(71)=0.5*VINT(1)*SQRT(XT2)  
        CALL PYSIGH(NCHN,SIGS)  
  110   SIGM(IXT2)=SIGM(IXT2)+SIGS  
  120   SIGSUM=SIGSUM+SIGM(IXT2)    
        SIGSUM=SIGSUM/(20.*MSTP(83))    
    
C...Reject result if sigma(parton-parton) is smaller than hadronic one. 
        IF(SIGSUM.LT.1.1*VINT(106)) THEN    
          IF(MSTP(122).GE.1) WRITE(MSTU(11),1100) PARP(82),SIGSUM   
          PARP(82)=0.9*PARP(82) 
          VINT(149)=4.*PARP(82)**2/VINT(2)  
          GOTO 100  
        ENDIF   
        IF(MSTP(122).GE.1) WRITE(MSTU(11),1200) PARP(82), SIGSUM    
    
C...Start iteration to find k factor.   
        YKE=SIGSUM/VINT(106)    
        SO=0.5  
        XI=0.   
        YI=0.   
        XK=0.5  
        IIT=0   
  130   IF(IIT.EQ.0) THEN   
          XK=2.*XK  
        ELSEIF(IIT.EQ.1) THEN   
          XK=0.5*XK 
        ELSE    
          XK=XI+(YKE-YI)*(XF-XI)/(YF-YI)    
        ENDIF   
    
C...Evaluate overlap integrals. 
        IF(MSTP(82).EQ.2) THEN  
          SP=0.5*PARU(1)*(1.-EXP(-XK))  
          SOP=SP/PARU(1)    
        ELSE    
          IF(MSTP(82).EQ.3) DELTAB=0.02 
          IF(MSTP(82).EQ.4) DELTAB=MIN(0.01,0.05*PARP(84))  
          SP=0. 
          SOP=0.    
          B=-0.5*DELTAB 
  140     B=B+DELTAB    
          IF(MSTP(82).EQ.3) THEN    
            OV=EXP(-B**2)/PARU(2)   
          ELSE  
            CQ2=PARP(84)**2 
            OV=((1.-PARP(83))**2*EXP(-MIN(100.,B**2))+2.*PARP(83)*  
     &      (1.-PARP(83))*2./(1.+CQ2)*EXP(-MIN(100.,B**2*2./(1.+CQ2)))+ 
     &      PARP(83)**2/CQ2*EXP(-MIN(100.,B**2/CQ2)))/PARU(2)   
          ENDIF 
          PACC=1.-EXP(-MIN(100.,PARU(1)*XK*OV)) 
          SP=SP+PARU(2)*B*DELTAB*PACC   
          SOP=SOP+PARU(2)*B*DELTAB*OV*PACC  
          IF(B.LT.1..OR.B*PACC.GT.1E-6) GOTO 140    
        ENDIF   
        YK=PARU(1)*XK*SO/SP 
    
C...Continue iteration until convergence.   
        IF(YK.LT.YKE) THEN  
          XI=XK 
          YI=YK 
          IF(IIT.EQ.1) IIT=2    
        ELSE    
          XF=XK 
          YF=YK 
          IF(IIT.EQ.0) IIT=1    
        ENDIF   
        IF(ABS(YK-YKE).GE.1E-5*YKE) GOTO 130    
    
C...Store some results for subsequent use.  
        VINT(145)=SIGSUM    
        VINT(146)=SOP/SO    
        VINT(147)=SOP/SP    
    
C...Initialize iteration in xT2 for hardest interaction.    
      ELSEIF(MMUL.EQ.2) THEN    
        IF(MSTP(82).LE.0) THEN  
        ELSEIF(MSTP(82).EQ.1) THEN  
          XT2=1.    
          XT2FAC=XSEC(96,1)/VINT(106)*VINT(149)/(1.-VINT(149))  
        ELSEIF(MSTP(82).EQ.2) THEN  
          XT2=1.    
          XT2FAC=VINT(146)*XSEC(96,1)/VINT(106)*VINT(149)*(1.+VINT(149))    
        ELSE    
          XC2=4.*CKIN(3)**2/VINT(2) 
          IF(CKIN(3).LE.CKIN(5).OR.MINT(82).GE.2) XC2=0.    
        ENDIF   
    
      ELSEIF(MMUL.EQ.3) THEN    
C...Low-pT or multiple interactions (first semihard interaction):   
C...choose xT2 according to dpT2/pT2**2*exp(-(sigma above pT2)/norm)    
C...or (MSTP(82)>=2) dpT2/(pT2+pT0**2)**2*exp(-....).   
        ISUB=MINT(1)    
        IF(MSTP(82).LE.0) THEN  
          XT2=0.    
        ELSEIF(MSTP(82).EQ.1) THEN  
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(RLU(0)))   
        ELSEIF(MSTP(82).EQ.2) THEN  
          IF(XT2.LT.1..AND.EXP(-XT2FAC*XT2/(VINT(149)*(XT2+ 
     &    VINT(149)))).GT.RLU(0)) XT2=1.    
          IF(XT2.GE.1.) THEN    
            XT2=(1.+VINT(149))*XT2FAC/(XT2FAC-(1.+VINT(149))*LOG(1.-    
     &      RLU(0)*(1.-EXP(-XT2FAC/(VINT(149)*(1.+VINT(149)))))))-  
     &      VINT(149)   
          ELSE  
            XT2=-XT2FAC/LOG(EXP(-XT2FAC/(XT2+VINT(149)))+RLU(0)*    
     &      (EXP(-XT2FAC/VINT(149))-EXP(-XT2FAC/(XT2+VINT(149)))))- 
     &      VINT(149)   
          ENDIF 
          XT2=MAX(0.01*VINT(149),XT2)   
        ELSE    
          XT2=(XC2+VINT(149))*(1.+VINT(149))/(1.+VINT(149)- 
     &    RLU(0)*(1.-XC2))-VINT(149)    
          XT2=MAX(0.01*VINT(149),XT2)   
        ENDIF   
        VINT(25)=XT2    
    
C...Low-pT: choose xT2, tau, y* and cos(theta-hat) fixed.   
        IF(MSTP(82).LE.1.AND.XT2.LT.VINT(149)) THEN 
          IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)-1   
          IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)-1 
          ISUB=95   
          MINT(1)=ISUB  
          VINT(21)=0.01*VINT(149)   
          VINT(22)=0.   
          VINT(23)=0.   
          VINT(25)=0.01*VINT(149)   
    
        ELSE    
C...Multiple interactions (first semihard interaction). 
C...Choose tau and y*. Calculate cos(theta-hat).    
          IF(RLU(0).LE.COEF(ISUB,1)) THEN   
            TAUP=(2.*(1.+SQRT(1.-XT2))/XT2-1.)**RLU(0)  
            TAU=XT2*(1.+TAUP)**2/(4.*TAUP)  
          ELSE  
            TAU=XT2*(1.+TAN(RLU(0)*ATAN(SQRT(1./XT2-1.)))**2)   
          ENDIF 
          VINT(21)=TAU  
          CALL PYKLIM(2)    
          RYST=RLU(0)   
          MYST=1    
          IF(RYST.GT.COEF(ISUB,7)) MYST=2   
          IF(RYST.GT.COEF(ISUB,7)+COEF(ISUB,8)) MYST=3  
          CALL PYKMAP(2,MYST,RLU(0))    
          VINT(23)=SQRT(MAX(0.,1.-XT2/TAU))*(-1)**INT(1.5+RLU(0))   
        ENDIF   
        VINT(71)=0.5*VINT(1)*SQRT(VINT(25)) 
    
C...Store results of cross-section calculation. 
      ELSEIF(MMUL.EQ.4) THEN    
        ISUB=MINT(1)    
        XTS=VINT(25)    
        IF(ISET(ISUB).EQ.1) XTS=VINT(21)    
        IF(ISET(ISUB).EQ.2) XTS=(4.*VINT(48)+2.*VINT(63)+2.*VINT(64))/  
     &  VINT(2) 
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) XTS=VINT(26) 
        RBIN=MAX(0.000001,MIN(0.999999,XTS*(1.+VINT(149))/  
     &  (XTS+VINT(149))))   
        IRBIN=INT(1.+20.*RBIN)  
        IF(ISUB.EQ.96) NMUL(IRBIN)=NMUL(IRBIN)+1    
        IF(ISUB.EQ.96) SIGM(IRBIN)=SIGM(IRBIN)+VINT(153)    
    
C...Choose impact parameter.    
      ELSEIF(MMUL.EQ.5) THEN    
        IF(MSTP(82).EQ.3) THEN  
          VINT(148)=RLU(0)/(PARU(2)*VINT(147))  
        ELSE    
          RTYPE=RLU(0)  
          CQ2=PARP(84)**2   
          IF(RTYPE.LT.(1.-PARP(83))**2) THEN    
            B2=-LOG(RLU(0)) 
          ELSEIF(RTYPE.LT.1.-PARP(83)**2) THEN  
            B2=-0.5*(1.+CQ2)*LOG(RLU(0))    
          ELSE  
            B2=-CQ2*LOG(RLU(0)) 
          ENDIF 
          VINT(148)=((1.-PARP(83))**2*EXP(-MIN(100.,B2))+2.*PARP(83)*   
     &    (1.-PARP(83))*2./(1.+CQ2)*EXP(-MIN(100.,B2*2./(1.+CQ2)))+ 
     &    PARP(83)**2/CQ2*EXP(-MIN(100.,B2/CQ2)))/(PARU(2)*VINT(147))   
        ENDIF   
    
C...Multiple interactions (variable impact parameter) : reject with 
C...probability exp(-overlap*cross-section above pT/normalization). 
        RNCOR=(IRBIN-20.*RBIN)*NMUL(IRBIN)  
        SIGCOR=(IRBIN-20.*RBIN)*SIGM(IRBIN) 
        DO 150 IBIN=IRBIN+1,20  
        RNCOR=RNCOR+NMUL(IBIN)  
  150   SIGCOR=SIGCOR+SIGM(IBIN)    
        SIGABV=(SIGCOR/RNCOR)*VINT(149)*(1.-XTS)/(XTS+VINT(149))    
        VINT(150)=EXP(-MIN(100.,VINT(146)*VINT(148)*SIGABV/VINT(106)))  
    
C...Generate additional multiple semihard interactions. 
      ELSEIF(MMUL.EQ.6) THEN    
    
C...Reconstruct strings in hard scattering. 
        ISUB=MINT(1)    
        NMAX=MINT(84)+4 
        IF(ISET(ISUB).EQ.1) NMAX=MINT(84)+2 
        NSTR=0  
        DO 170 I=MINT(84)+1,NMAX    
        KCS=KCHG(LUCOMP(K(I,2)),2)*ISIGN(1,K(I,2))  
        IF(KCS.EQ.0) GOTO 170   
        DO 160 J=1,4    
        IF(KCS.EQ.1.AND.(J.EQ.2.OR.J.EQ.4)) GOTO 160    
        IF(KCS.EQ.-1.AND.(J.EQ.1.OR.J.EQ.3)) GOTO 160   
        IF(J.LE.2) THEN 
          IST=MOD(K(I,J+3)/MSTU(5),MSTU(5)) 
        ELSE    
          IST=MOD(K(I,J+1),MSTU(5)) 
        ENDIF   
        IF(IST.LT.MINT(84).OR.IST.GT.I) GOTO 160    
        IF(KCHG(LUCOMP(K(IST,2)),2).EQ.0) GOTO 160  
        NSTR=NSTR+1 
        IF(J.EQ.1.OR.J.EQ.4) THEN   
          KSTR(NSTR,1)=I    
          KSTR(NSTR,2)=IST  
        ELSE    
          KSTR(NSTR,1)=IST  
          KSTR(NSTR,2)=I    
        ENDIF   
  160   CONTINUE    
  170   CONTINUE    
    
C...Set up starting values for iteration in xT2.    
        XT2=VINT(25)    
        IF(ISET(ISUB).EQ.1) XT2=VINT(21)    
        IF(ISET(ISUB).EQ.2) XT2=(4.*VINT(48)+2.*VINT(63)+2.*VINT(64))/  
     &  VINT(2) 
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) XT2=VINT(26) 
        ISUB=96 
        MINT(1)=96  
        IF(MSTP(82).LE.1) THEN  
          XT2FAC=XSEC(ISUB,1)*VINT(149)/((1.-VINT(149))*VINT(106))  
        ELSE    
          XT2FAC=VINT(146)*VINT(148)*XSEC(ISUB,1)/VINT(106)*    
     &    VINT(149)*(1.+VINT(149))  
        ENDIF   
        VINT(63)=0. 
        VINT(64)=0. 
        VINT(151)=0.    
        VINT(152)=0.    
        VINT(143)=1.-VINT(141)  
        VINT(144)=1.-VINT(142)  
    
C...Iterate downwards in xT2.   
  180   IF(MSTP(82).LE.1) THEN  
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(RLU(0)))   
          IF(XT2.LT.VINT(149)) GOTO 220 
        ELSE    
          IF(XT2.LE.0.01*VINT(149)) GOTO 220    
          XT2=XT2FAC*(XT2+VINT(149))/(XT2FAC-(XT2+VINT(149))*   
     &    LOG(RLU(0)))-VINT(149)    
          IF(XT2.LE.0.) GOTO 220    
          XT2=MAX(0.01*VINT(149),XT2)   
        ENDIF   
        VINT(25)=XT2    
    
C...Choose tau and y*. Calculate cos(theta-hat).    
        IF(RLU(0).LE.COEF(ISUB,1)) THEN 
          TAUP=(2.*(1.+SQRT(1.-XT2))/XT2-1.)**RLU(0)    
          TAU=XT2*(1.+TAUP)**2/(4.*TAUP)    
        ELSE    
          TAU=XT2*(1.+TAN(RLU(0)*ATAN(SQRT(1./XT2-1.)))**2) 
        ENDIF   
        VINT(21)=TAU    
        CALL PYKLIM(2)  
        RYST=RLU(0) 
        MYST=1  
        IF(RYST.GT.COEF(ISUB,7)) MYST=2 
        IF(RYST.GT.COEF(ISUB,7)+COEF(ISUB,8)) MYST=3    
        CALL PYKMAP(2,MYST,RLU(0))  
        VINT(23)=SQRT(MAX(0.,1.-XT2/TAU))*(-1)**INT(1.5+RLU(0)) 
    
C...Check that x not used up. Accept or reject kinematical variables.   
        X1M=SQRT(TAU)*EXP(VINT(22)) 
        X2M=SQRT(TAU)*EXP(-VINT(22))    
        IF(VINT(143)-X1M.LT.0.01.OR.VINT(144)-X2M.LT.0.01) GOTO 180 
        VINT(71)=0.5*VINT(1)*SQRT(XT2)  
        CALL PYSIGH(NCHN,SIGS)  
        IF(SIGS.LT.XSEC(ISUB,1)*RLU(0)) GOTO 180    
    
C...Reset K, P and V vectors. Select some variables.    
        DO 190 I=N+1,N+2    
        DO 190 J=1,5    
        K(I,J)=0    
        P(I,J)=0.   
  190   V(I,J)=0.   
        RFLAV=RLU(0)    
        PT=0.5*VINT(1)*SQRT(XT2)    
        PHI=PARU(2)*RLU(0)  
        CTH=VINT(23)    
    
C...Add first parton to event record.   
        K(N+1,1)=3  
        K(N+1,2)=21 
        IF(RFLAV.GE.MAX(PARP(85),PARP(86))) K(N+1,2)=   
     &  1+INT((2.+PARJ(2))*RLU(0))  
        P(N+1,1)=PT*COS(PHI)    
        P(N+1,2)=PT*SIN(PHI)    
        P(N+1,3)=0.25*VINT(1)*(VINT(41)*(1.+CTH)-VINT(42)*(1.-CTH)) 
        P(N+1,4)=0.25*VINT(1)*(VINT(41)*(1.+CTH)+VINT(42)*(1.-CTH)) 
        P(N+1,5)=0. 
    
C...Add second parton to event record.  
        K(N+2,1)=3  
        K(N+2,2)=21 
        IF(K(N+1,2).NE.21) K(N+2,2)=-K(N+1,2)   
        P(N+2,1)=-P(N+1,1)  
        P(N+2,2)=-P(N+1,2)  
        P(N+2,3)=0.25*VINT(1)*(VINT(41)*(1.-CTH)-VINT(42)*(1.+CTH)) 
        P(N+2,4)=0.25*VINT(1)*(VINT(41)*(1.-CTH)+VINT(42)*(1.+CTH)) 
        P(N+2,5)=0. 
    
        IF(RFLAV.LT.PARP(85).AND.NSTR.GE.1) THEN    
C....Choose relevant string pieces to place gluons on.  
          DO 210 I=N+1,N+2  
          DMIN=1E8  
          DO 200 ISTR=1,NSTR    
          I1=KSTR(ISTR,1)   
          I2=KSTR(ISTR,2)   
          DIST=(P(I,4)*P(I1,4)-P(I,1)*P(I1,1)-P(I,2)*P(I1,2)-   
     &    P(I,3)*P(I1,3))*(P(I,4)*P(I2,4)-P(I,1)*P(I2,1)-   
     &    P(I,2)*P(I2,2)-P(I,3)*P(I2,3))/MAX(1.,P(I1,4)*P(I2,4)-    
     &    P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)-P(I1,3)*P(I2,3))  
          IF(ISTR.EQ.1.OR.DIST.LT.DMIN) THEN    
            DMIN=DIST   
            IST1=I1 
            IST2=I2 
            ISTM=ISTR   
          ENDIF 
  200     CONTINUE  
    
C....Colour flow adjustments, new string pieces.    
          IF(K(IST1,4)/MSTU(5).EQ.IST2) K(IST1,4)=MSTU(5)*I+    
     &    MOD(K(IST1,4),MSTU(5))    
          IF(MOD(K(IST1,5),MSTU(5)).EQ.IST2) K(IST1,5)= 
     &    MSTU(5)*(K(IST1,5)/MSTU(5))+I 
          K(I,5)=MSTU(5)*IST1   
          K(I,4)=MSTU(5)*IST2   
          IF(K(IST2,5)/MSTU(5).EQ.IST1) K(IST2,5)=MSTU(5)*I+    
     &    MOD(K(IST2,5),MSTU(5))    
          IF(MOD(K(IST2,4),MSTU(5)).EQ.IST1) K(IST2,4)= 
     &    MSTU(5)*(K(IST2,4)/MSTU(5))+I 
          KSTR(ISTM,2)=I    
          KSTR(NSTR+1,1)=I  
          KSTR(NSTR+1,2)=IST2   
  210     NSTR=NSTR+1   
    
C...String drawing and colour flow for gluon loop.  
        ELSEIF(K(N+1,2).EQ.21) THEN 
          K(N+1,4)=MSTU(5)*(N+2)    
          K(N+1,5)=MSTU(5)*(N+2)    
          K(N+2,4)=MSTU(5)*(N+1)    
          K(N+2,5)=MSTU(5)*(N+1)    
          KSTR(NSTR+1,1)=N+1    
          KSTR(NSTR+1,2)=N+2    
          KSTR(NSTR+2,1)=N+2    
          KSTR(NSTR+2,2)=N+1    
          NSTR=NSTR+2   
    
C...String drawing and colour flow for q-qbar pair. 
        ELSE    
          K(N+1,4)=MSTU(5)*(N+2)    
          K(N+2,5)=MSTU(5)*(N+1)    
          KSTR(NSTR+1,1)=N+1    
          KSTR(NSTR+1,2)=N+2    
          NSTR=NSTR+1   
        ENDIF   
    
C...Update remaining energy; iterate.   
        N=N+2   
        IF(N.GT.MSTU(4)-MSTU(32)-10) THEN   
          CALL LUERRM(11,'(PYMULT:) no more memory left in LUJETS') 
          IF(MSTU(21).GE.1) RETURN  
        ENDIF   
        MINT(31)=MINT(31)+1 
        VINT(151)=VINT(151)+VINT(41)    
        VINT(152)=VINT(152)+VINT(42)    
        VINT(143)=VINT(143)-VINT(41)    
        VINT(144)=VINT(144)-VINT(42)    
        IF(MINT(31).LT.240) GOTO 180    
  220   CONTINUE    
      ENDIF 
    
C...Format statements for printout. 
 1000 FORMAT(/1X,'****** PYMULT: initialization of multiple inter', 
     &'actions for MSTP(82) =',I2,' ******')    
 1100 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,    
     &E9.2,' mb: rejected') 
 1200 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,    
     &E9.2,' mb: accepted') 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYREMN(IPU1,IPU2)  
    
C...Adds on target remnants (one or two from each side) and 
C...includes primordial kT. 
      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/
      COMMON/HISTRNG/NFP(300,15),PPHI(300,15),NFT(300,15),PTHI(300,15)
      SAVE  /HISTRNG/
C...COMMON BLOCK FROM HIJING
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(6),IS(2),ROBO(5)   
    
C...Special case for lepton-lepton interaction. 
      IF(MINT(43).EQ.1) THEN    
        DO 100 JT=1,2   
        I=MINT(83)+JT+2 
        K(I,1)=21   
        K(I,2)=K(I-2,2) 
        K(I,3)=I-2  
        DO 100 J=1,5    
  100   P(I,J)=P(I-2,J) 
      ENDIF 
    
C...Find event type, set pointers.  
      IF(IPU1.EQ.0.AND.IPU2.EQ.0) RETURN    
      ISUB=MINT(1)  
      ILEP=0    
      IF(IPU1.EQ.0) ILEP=1  
      IF(IPU2.EQ.0) ILEP=2  
      IF(ISUB.EQ.95) ILEP=-1    
      IF(ILEP.EQ.1) IQ=MINT(84)+1   
      IF(ILEP.EQ.2) IQ=MINT(84)+2   
      IP=MAX(IPU1,IPU2) 
      ILEPR=MINT(83)+5-ILEP 
      NS=N  
    
C...Define initial partons, including primordial kT.    
  110 DO 130 JT=1,2 
      I=MINT(83)+JT+2   
      IF(JT.EQ.1) IPU=IPU1  
      IF(JT.EQ.2) IPU=IPU2  
      K(I,1)=21 
      K(I,3)=I-2    
      IF(ISUB.EQ.95) THEN   
        K(I,2)=21   
        SHS=0.  
      ELSEIF(MINT(40+JT).EQ.1.AND.IPU.NE.0) THEN    
        K(I,2)=K(IPU,2) 
        P(I,5)=P(IPU,5) 
        P(I,1)=0.   
        P(I,2)=0.   
        PMS(JT)=P(I,5)**2   
      ELSEIF(IPU.NE.0) THEN 
        K(I,2)=K(IPU,2) 
        P(I,5)=P(IPU,5) 
C...No primordial kT or chosen according to truncated Gaussian or   
C...exponential.    
        IF(IHPR2(5).LE.0) THEN
120	     IF(MSTP(91).LE.0) THEN
               PT=0. 
             ELSEIF(MSTP(91).EQ.1) THEN
               PT=PARP(91)*SQRT(-LOG(RLU(0)))
             ELSE    
               RPT1=RLU(0)   
               RPT2=RLU(0)   
               PT=-PARP(92)*LOG(RPT1*RPT2)   
             ENDIF   
             IF(PT.GT.PARP(93)) GOTO 120 
	     PHI=PARU(2)*RLU(0)  
	     RPT1=PT*COS(PHI)  
	     RPT2=PT*SIN(PHI)
	ELSE IF(IHPR2(5).EQ.1) THEN
	     IF(JT.EQ.1) JPT=NFP(IHNT2(11),11)
	     IF(JT.EQ.2) JPT=NFT(IHNT2(12),11)
1205	     PTGS=PARP(91)*SQRT(-LOG(RLU(0)))
	     IF(PTGS.GT.PARP(93)) GO TO 1205
	     PHI=2.0*HIPR1(40)*RLU(0)
	     RPT1=PTGS*COS(PHI)
	     RPT2=PTGS*SIN(PHI)
	     DO 1210 I_INT=1,JPT-1
		PKCMX=16.0
		PKCSQ=4.0*HIRND2(3,0.0,PKCMX)
		PHI=2.0*HIPR1(40)*RLU(0)
		PTABS=SQRT(ABS(PKCSQ))
		RPT1=RPT1+PTABS*COS(PHI)
		RPT2=RPT2+PTABS*SIN(PHI)
1210	     CONTINUE
	ENDIF
C			********When initial interaction among soft partons is
C				assumed the primordial pt comes from the sum of
C				pt of JPT-1 number of initial interaction, JPT
C				is the number of interaction including present
C				one that nucleon hassuffered 
	P(I,1)=RPT1
	P(I,2)=RPT2  
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2   
      ELSE  
        K(I,2)=K(IQ,2)  
        Q2=VINT(52) 
        P(I,5)=-SQRT(Q2)    
        PMS(JT)=-Q2 
        SHS=(1.-VINT(43-JT))*Q2/VINT(43-JT)+VINT(5-JT)**2   
      ENDIF 
  130 CONTINUE  
    
C...Kinematics construction for initial partons.    
      I1=MINT(83)+3 
      I2=MINT(83)+4 
      IF(ILEP.EQ.0) SHS=VINT(141)*VINT(142)*VINT(2)+    
     &(P(I1,1)+P(I2,1))**2+(P(I1,2)+P(I2,2))**2 
      SHR=SQRT(MAX(0.,SHS)) 
      IF(ILEP.EQ.0) THEN    
        IF((SHS-PMS(1)-PMS(2))**2-4.*PMS(1)*PMS(2).LE.0.) GOTO 110  
        P(I1,4)=0.5*(SHR+(PMS(1)-PMS(2))/SHR)   
        P(I1,3)=SQRT(MAX(0.,P(I1,4)**2-PMS(1))) 
        P(I2,4)=SHR-P(I1,4) 
        P(I2,3)=-P(I1,3)    
      ELSEIF(ILEP.EQ.1) THEN    
        P(I1,4)=P(IQ,4) 
        P(I1,3)=P(IQ,3) 
        P(I2,4)=P(IP,4) 
        P(I2,3)=P(IP,3) 
      ELSEIF(ILEP.EQ.2) THEN    
        P(I1,4)=P(IP,4) 
        P(I1,3)=P(IP,3) 
        P(I2,4)=P(IQ,4) 
        P(I2,3)=P(IQ,3) 
      ENDIF 
      IF(MINT(43).EQ.1) RETURN  
    
C...Transform partons to overall CM-frame (not for leptoproduction).    
      IF(ILEP.EQ.0) THEN    
        ROBO(3)=(P(I1,1)+P(I2,1))/SHR   
        ROBO(4)=(P(I1,2)+P(I2,2))/SHR   
        CALL LUDBRB(I1,I2,0.,0.,-DBLE(ROBO(3)),-DBLE(ROBO(4)),0D0)  
        ROBO(2)=ULANGL(P(I1,1),P(I1,2)) 
        CALL LUDBRB(I1,I2,0.,-ROBO(2),0D0,0D0,0D0)  
        ROBO(1)=ULANGL(P(I1,3),P(I1,1)) 
        CALL LUDBRB(I1,I2,-ROBO(1),0.,0D0,0D0,0D0)  
        NMAX=MAX(MINT(52),IPU1,IPU2)    
        CALL LUDBRB(I1,NMAX,ROBO(1),ROBO(2),DBLE(ROBO(3)),DBLE(ROBO(4)),    
     &  0D0)    
        ROBO(5)=MAX(-0.999999,MIN(0.999999,(VINT(141)-VINT(142))/   
     &  (VINT(141)+VINT(142)))) 
        CALL LUDBRB(I1,NMAX,0.,0.,0D0,0D0,DBLE(ROBO(5)))    
      ENDIF 
    
C...Check invariant mass of remnant system: 
C...hadronic events or leptoproduction. 
      IF(ILEP.LE.0) THEN    
        IF(MSTP(81).LE.0.OR.MSTP(82).LE.0.OR.ISUB.EQ.95) THEN   
          VINT(151)=0.  
          VINT(152)=0.  
        ENDIF   
        PEH=P(I1,4)+P(I2,4)+0.5*VINT(1)*(VINT(151)+VINT(152))   
        PZH=P(I1,3)+P(I2,3)+0.5*VINT(1)*(VINT(151)-VINT(152))   
        SHH=(VINT(1)-PEH)**2-(P(I1,1)+P(I2,1))**2-(P(I1,2)+P(I2,2))**2- 
     &  PZH**2  
        PMMIN=P(MINT(83)+1,5)+P(MINT(83)+2,5)+ULMASS(K(I1,2))+  
     &  ULMASS(K(I2,2)) 
        IF(SHR.GE.VINT(1).OR.SHH.LE.(PMMIN+PARP(111))**2) THEN  
          MINT(51)=1    
          RETURN    
        ENDIF   
        SHR=SQRT(SHH+(P(I1,1)+P(I2,1))**2+(P(I1,2)+P(I2,2))**2) 
      ELSE  
        PEI=P(IQ,4)+P(IP,4) 
        PZI=P(IQ,3)+P(IP,3) 
        PMS(ILEP)=MAX(0.,PEI**2-PZI**2) 
        PMMIN=P(ILEPR-2,5)+ULMASS(K(ILEPR,2))+SQRT(PMS(ILEP))   
        IF(SHR.LE.PMMIN+PARP(111)) THEN 
          MINT(51)=1    
          RETURN    
        ENDIF   
      ENDIF 
    
C...Subdivide remnant if necessary, store first parton. 
  140 I=NS  
      DO 190 JT=1,2 
      IF(JT.EQ.ILEP) GOTO 190   
      IF(JT.EQ.1) IPU=IPU1  
      IF(JT.EQ.2) IPU=IPU2  
      CALL PYSPLI(MINT(10+JT),MINT(12+JT),KFLCH(JT),KFLSP(JT))  
      I=I+1 
      IS(JT)=I  
      DO 150 J=1,5  
      K(I,J)=0  
      P(I,J)=0. 
  150 V(I,J)=0. 
      K(I,1)=3  
      K(I,2)=KFLSP(JT)  
      K(I,3)=MINT(83)+JT    
      P(I,5)=ULMASS(K(I,2)) 
    
C...First parton colour connections and transverse mass.    
      KFLS=(3-KCHG(LUCOMP(KFLSP(JT)),2)*ISIGN(1,KFLSP(JT)))/2   
      K(I,KFLS+3)=IPU   
      K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I    
      IF(KFLCH(JT).EQ.0) THEN   
        P(I,1)=-P(MINT(83)+JT+2,1)  
        P(I,2)=-P(MINT(83)+JT+2,2)  
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2   
    
C...When extra remnant parton or hadron: find relative pT, store.   
      ELSE  
        CALL LUPTDI(1,P(I,1),P(I,2))    
        PMS(JT+2)=P(I,5)**2+P(I,1)**2+P(I,2)**2 
        I=I+1   
        DO 160 J=1,5    
        K(I,J)=0    
        P(I,J)=0.   
  160   V(I,J)=0.   
        K(I,1)=1    
        K(I,2)=KFLCH(JT)    
        K(I,3)=MINT(83)+JT  
        P(I,5)=ULMASS(K(I,2))   
        P(I,1)=-P(MINT(83)+JT+2,1)-P(I-1,1) 
        P(I,2)=-P(MINT(83)+JT+2,2)-P(I-1,2) 
        PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2 
C...Relative distribution of energy for particle into two jets. 
        IMB=1   
        IF(MOD(MINT(10+JT)/1000,10).NE.0) IMB=2 
        IF(IABS(KFLCH(JT)).LE.10.OR.KFLCH(JT).EQ.21) THEN   
          CHIK=PARP(92+2*IMB)   
          IF(MSTP(92).LE.1) THEN    
            IF(IMB.EQ.1) CHI(JT)=RLU(0) 
            IF(IMB.EQ.2) CHI(JT)=1.-SQRT(RLU(0))    
          ELSEIF(MSTP(92).EQ.2) THEN    
            CHI(JT)=1.-RLU(0)**(1./(1.+CHIK))   
          ELSEIF(MSTP(92).EQ.3) THEN    
            CUT=2.*0.3/VINT(1)  
  170       CHI(JT)=RLU(0)**2   
            IF((CHI(JT)**2/(CHI(JT)**2+CUT**2))**0.25*(1.-CHI(JT))**CHIK    
     &      .LT.RLU(0)) GOTO 170    
          ELSE  
            CUT=2.*0.3/VINT(1)  
            CUTR=(1.+SQRT(1.+CUT**2))/CUT   
  180       CHIR=CUT*CUTR**RLU(0)   
            CHI(JT)=(CHIR**2-CUT**2)/(2.*CHIR)  
            IF((1.-CHI(JT))**CHIK.LT.RLU(0)) GOTO 180   
          ENDIF 
C...Relative distribution of energy for particle into jet plus particle.    
        ELSE    
          IF(MSTP(92).LE.1) THEN    
            IF(IMB.EQ.1) CHI(JT)=RLU(0) 
            IF(IMB.EQ.2) CHI(JT)=1.-SQRT(RLU(0))    
          ELSE  
            CHI(JT)=1.-RLU(0)**(1./(1.+PARP(93+2*IMB))) 
          ENDIF 
          IF(MOD(KFLCH(JT)/1000,10).NE.0) CHI(JT)=1.-CHI(JT)    
        ENDIF   
        PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/(1.-CHI(JT))    
        KFLS=KCHG(LUCOMP(KFLCH(JT)),2)*ISIGN(1,KFLCH(JT))   
        IF(KFLS.NE.0) THEN  
          K(I,1)=3  
          KFLS=(3-KFLS)/2   
          K(I,KFLS+3)=IPU   
          K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I    
        ENDIF   
      ENDIF 
  190 CONTINUE  
      IF(SHR.LE.SQRT(PMS(1))+SQRT(PMS(2))) GOTO 140 
      N=I   
    
C...Reconstruct kinematics of remnants. 
      DO 200 JT=1,2 
      IF(JT.EQ.ILEP) GOTO 200   
      PE=0.5*(SHR+(PMS(JT)-PMS(3-JT))/SHR)  
      PZ=SQRT(PE**2-PMS(JT))    
      IF(KFLCH(JT).EQ.0) THEN   
        P(IS(JT),4)=PE  
        P(IS(JT),3)=PZ*(-1)**(JT-1) 
      ELSE  
        PW1=CHI(JT)*(PE+PZ) 
        P(IS(JT)+1,4)=0.5*(PW1+PMS(JT+4)/PW1)   
        P(IS(JT)+1,3)=0.5*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)  
        P(IS(JT),4)=PE-P(IS(JT)+1,4)    
        P(IS(JT),3)=PZ*(-1)**(JT-1)-P(IS(JT)+1,3)   
      ENDIF 
  200 CONTINUE  
    
C...Hadronic events: boost remnants to correct longitudinal frame.  
      IF(ILEP.LE.0) THEN    
        CALL LUDBRB(NS+1,N,0.,0.,0D0,0D0,-DBLE(PZH/(VINT(1)-PEH)))  
C...Leptoproduction events: boost colliding subsystem.  
      ELSE  
        NMAX=MAX(IP,MINT(52))   
        PEF=SHR-PE  
        PZF=PZ*(-1)**(ILEP-1)   
        PT2=P(ILEPR,1)**2+P(ILEPR,2)**2 
        PHIPT=ULANGL(P(ILEPR,1),P(ILEPR,2)) 
        CALL LUDBRB(MINT(84)+1,NMAX,0.,-PHIPT,0D0,0D0,0D0)  
        RQP=P(IQ,3)*(PT2+PEI**2)-P(IQ,4)*PEI*PZI    
        SINTH=P(IQ,4)*SQRT(PT2*(PT2+PEI**2)/(RQP**2+PT2*    
     &  P(IQ,4)**2*PZI**2))*SIGN(1.,-RQP)   
        CALL LUDBRB(MINT(84)+1,NMAX,ASIN(SINTH),0.,0D0,0D0,0D0) 
        BETAX=(-PEI*PZI*SINTH+SQRT(PT2*(PT2+PEI**2-(PZI*SINTH)**2)))/   
     &  (PT2+PEI**2)    
        CALL LUDBRB(MINT(84)+1,NMAX,0.,0.,DBLE(BETAX),0D0,0D0)  
        CALL LUDBRB(MINT(84)+1,NMAX,0.,PHIPT,0D0,0D0,0D0)   
        PEM=P(IQ,4)+P(IP,4) 
        PZM=P(IQ,3)+P(IP,3) 
        BETAZ=(-PEM*PZM+PZF*SQRT(PZF**2+PEM**2-PZM**2))/(PZF**2+PEM**2) 
        CALL LUDBRB(MINT(84)+1,NMAX,0.,0.,0D0,0D0,DBLE(BETAZ))  
        CALL LUDBRB(I1,I2,ASIN(SINTH),0.,DBLE(BETAX),0D0,0D0)   
        CALL LUDBRB(I1,I2,0.,PHIPT,0D0,0D0,DBLE(BETAZ)) 
      ENDIF 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYRESD 
    
C...Allows resonances to decay (including parton showers for hadronic   
C...channels).  
      IMPLICIT DOUBLE PRECISION(D)  
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      DIMENSION IREF(10,6),KDCY(2),KFL1(2),KFL2(2),NSD(2),ILIN(6),  
     &COUP(6,4),PK(6,4),PKK(6,6),CTHE(2),PHI(2),WDTP(0:40), 
     &WDTE(0:40,0:5)    
      COMPLEX FGK,HA(6,6),HC(6,6)   
    
C...The F, Xi and Xj functions of Gunion and Kunszt 
C...(Phys. Rev. D33, 665, plus errata from the authors).    
      FGK(I1,I2,I3,I4,I5,I6)=4.*HA(I1,I3)*HC(I2,I6)*(HA(I1,I5)* 
     &HC(I1,I4)+HA(I3,I5)*HC(I3,I4))    
      DIGK(DT,DU)=-4.*D34*D56+DT*(3.*DT+4.*DU)+DT**2*(DT*DU/(D34*D56)-  
     &2.*(1./D34+1./D56)*(DT+DU)+2.*(D34/D56+D56/D34))  
      DJGK(DT,DU)=8.*(D34+D56)**2-8.*(D34+D56)*(DT+DU)-6.*DT*DU-    
     &2.*DT*DU*(DT*DU/(D34*D56)-2.*(1./D34+1./D56)*(DT+DU)+ 
     &2.*(D34/D56+D56/D34)) 
    
C...Define initial two objects, initialize loop.    
      ISUB=MINT(1)  
      SH=VINT(44)   
      IREF(1,5)=0   
      IREF(1,6)=0   
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        IREF(1,1)=MINT(84)+2+ISET(ISUB) 
        IREF(1,2)=0 
        IREF(1,3)=MINT(83)+6+ISET(ISUB) 
        IREF(1,4)=0 
      ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        IREF(1,1)=MINT(84)+1+ISET(ISUB) 
        IREF(1,2)=MINT(84)+2+ISET(ISUB) 
        IREF(1,3)=MINT(83)+5+ISET(ISUB) 
        IREF(1,4)=MINT(83)+6+ISET(ISUB) 
      ENDIF 
      NP=1  
      IP=0  
  100 IP=IP+1   
      NINH=0    
    
C...Loop over one/two resonances; reset decay rates.    
      JTMAX=2   
      IF(IP.EQ.1.AND.(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3)) JTMAX=1  
      DO 140 JT=1,JTMAX 
      KDCY(JT)=0    
      KFL1(JT)=0    
      KFL2(JT)=0    
      NSD(JT)=IREF(IP,JT)   
      ID=IREF(IP,JT)    
      IF(ID.EQ.0) GOTO 140  
      KFA=IABS(K(ID,2)) 
      IF(KFA.LT.23.OR.KFA.GT.40) GOTO 140   
      IF(MDCY(KFA,1).NE.0) THEN 
        IF(ISUB.EQ.1.OR.ISUB.EQ.141) MINT(61)=1 
        CALL PYWIDT(KFA,P(ID,5),WDTP,WDTE)  
        IF(KCHG(KFA,3).EQ.0) THEN   
          IPM=2 
        ELSE    
          IPM=(5+ISIGN(1,K(ID,2)))/2    
        ENDIF   
        IF(JTMAX.EQ.1.OR.IABS(K(IREF(IP,1),2)).NE.IABS(K(IREF(IP,2),2)))    
     &  THEN    
          I12=4 
        ELSE    
          IF(JT.EQ.1) I12=INT(4.5+RLU(0))   
          I12=9-I12 
        ENDIF   
        RKFL=(WDTE(0,1)+WDTE(0,IPM)+WDTE(0,I12))*RLU(0) 
        DO 120 I=1,MDCY(KFA,3)  
        IDC=I+MDCY(KFA,2)-1 
        KFL1(JT)=KFDP(IDC,1)*ISIGN(1,K(ID,2))   
        KFL2(JT)=KFDP(IDC,2)*ISIGN(1,K(ID,2))   
        RKFL=RKFL-(WDTE(I,1)+WDTE(I,IPM)+WDTE(I,I12))   
        IF(RKFL.LE.0.) GOTO 130 
  120   CONTINUE    
  130   CONTINUE    
      ENDIF 
    
C...Summarize result on decay channel chosen.   
      IF((KFA.EQ.23.OR.KFA.EQ.24).AND.KFL1(JT).EQ.0) NINH=NINH+1    
      IF(KFL1(JT).EQ.0) GOTO 140    
      KDCY(JT)=2    
      IF(IABS(KFL1(JT)).LE.10.OR.KFL1(JT).EQ.21) KDCY(JT)=1 
      IF((IABS(KFL1(JT)).GE.23.AND.IABS(KFL1(JT)).LE.25).OR.    
     &(IABS(KFL1(JT)).EQ.37)) KDCY(JT)=3    
      NSD(JT)=N 
    
C...Fill decay products, prepared for parton showers for quarks.    
      IF(KDCY(JT).EQ.1) THEN    
        CALL LU2ENT(-(N+1),KFL1(JT),KFL2(JT),P(ID,5))   
      ELSE  
        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))  
      ENDIF 
      IF(JTMAX.EQ.1) THEN   
        CTHE(JT)=VINT(13)+(VINT(33)-VINT(13)+VINT(34)-VINT(14))*RLU(0)  
        IF(CTHE(JT).GT.VINT(33)) CTHE(JT)=CTHE(JT)+VINT(14)-VINT(33)    
        PHI(JT)=VINT(24)    
      ELSE  
        CTHE(JT)=2.*RLU(0)-1.   
        PHI(JT)=PARU(2)*RLU(0)  
      ENDIF 
  140 CONTINUE  
      IF(MINT(3).EQ.1.AND.IP.EQ.1) THEN 
        MINT(25)=KFL1(1)    
        MINT(26)=KFL2(1)    
      ENDIF 
      IF(JTMAX.EQ.1.AND.KDCY(1).EQ.0) GOTO 530  
      IF(JTMAX.EQ.2.AND.KDCY(1).EQ.0.AND.KDCY(2).EQ.0) GOTO 530 
      IF(MSTP(45).LE.0.OR.IREF(IP,2).EQ.0.OR.NINH.GE.1) GOTO 500    
      IF(K(IREF(1,1),2).EQ.25.AND.IP.EQ.1) GOTO 500 
      IF(K(IREF(1,1),2).EQ.25.AND.KDCY(1)*KDCY(2).EQ.0) GOTO 500    
    
C...Order incoming partons and outgoing resonances. 
      ILIN(1)=MINT(84)+1    
      IF(K(MINT(84)+1,2).GT.0) ILIN(1)=MINT(84)+2   
      IF(K(ILIN(1),2).EQ.21) ILIN(1)=2*MINT(84)+3-ILIN(1)   
      ILIN(2)=2*MINT(84)+3-ILIN(1)  
      IMIN=1    
      IF(IREF(IP,5).EQ.25) IMIN=3   
      IMAX=2    
      IORD=1    
      IF(K(IREF(IP,1),2).EQ.23) IORD=2  
      IF(K(IREF(IP,1),2).EQ.24.AND.K(IREF(IP,2),2).EQ.-24) IORD=2   
      IF(IABS(K(IREF(IP,IORD),2)).EQ.25) IORD=3-IORD    
      IF(KDCY(IORD).EQ.0) IORD=3-IORD   
    
C...Order decay products of resonances. 
      DO 390 JT=IORD,3-IORD,3-2*IORD    
      IF(KDCY(JT).EQ.0) THEN    
        ILIN(IMAX+1)=NSD(JT)    
        IMAX=IMAX+1 
      ELSEIF(K(NSD(JT)+1,2).GT.0) THEN  
        ILIN(IMAX+1)=N+2*JT-1   
        ILIN(IMAX+2)=N+2*JT 
        IMAX=IMAX+2 
        K(N+2*JT-1,2)=K(NSD(JT)+1,2)    
        K(N+2*JT,2)=K(NSD(JT)+2,2)  
      ELSE  
        ILIN(IMAX+1)=N+2*JT 
        ILIN(IMAX+2)=N+2*JT-1   
        IMAX=IMAX+2 
        K(N+2*JT-1,2)=K(NSD(JT)+1,2)    
        K(N+2*JT,2)=K(NSD(JT)+2,2)  
      ENDIF 
  390 CONTINUE  
    
C...Find charge, isospin, left- and righthanded couplings.  
      XW=PARU(102)  
      DO 410 I=IMIN,IMAX    
      DO 400 J=1,4  
  400 COUP(I,J)=0.  
      KFA=IABS(K(ILIN(I),2))    
      IF(KFA.GT.20) GOTO 410    
      COUP(I,1)=LUCHGE(KFA)/3.  
      COUP(I,2)=(-1)**MOD(KFA,2)    
      COUP(I,4)=-2.*COUP(I,1)*XW    
      COUP(I,3)=COUP(I,2)+COUP(I,4) 
  410 CONTINUE  
      SQMZ=PMAS(23,1)**2    
      GZMZ=PMAS(23,1)*PMAS(23,2)    
      SQMW=PMAS(24,1)**2    
      GZMW=PMAS(24,1)*PMAS(24,2)    
      SQMZP=PMAS(32,1)**2   
      GZMZP=PMAS(32,1)*PMAS(32,2)   
    
C...Select random angles; construct massless four-vectors.  
  420 DO 430 I=N+1,N+4  
      K(I,1)=1  
      DO 430 J=1,5  
  430 P(I,J)=0. 
      DO 440 JT=1,JTMAX 
      IF(KDCY(JT).EQ.0) GOTO 440    
      ID=IREF(IP,JT)    
      P(N+2*JT-1,3)=0.5*P(ID,5) 
      P(N+2*JT-1,4)=0.5*P(ID,5) 
      P(N+2*JT,3)=-0.5*P(ID,5)  
      P(N+2*JT,4)=0.5*P(ID,5)   
      CTHE(JT)=2.*RLU(0)-1. 
      PHI(JT)=PARU(2)*RLU(0)    
      CALL LUDBRB(N+2*JT-1,N+2*JT,ACOS(CTHE(JT)),PHI(JT),   
     &DBLE(P(ID,1)/P(ID,4)),DBLE(P(ID,2)/P(ID,4)),DBLE(P(ID,3)/P(ID,4)))    
  440 CONTINUE  
    
C...Store incoming and outgoing momenta, with random rotation to    
C...avoid accidental zeroes in HA expressions.  
      DO 450 I=1,IMAX   
      K(N+4+I,1)=1  
      P(N+4+I,4)=SQRT(P(ILIN(I),1)**2+P(ILIN(I),2)**2+P(ILIN(I),3)**2+  
     &P(ILIN(I),5)**2)  
      P(N+4+I,5)=P(ILIN(I),5)   
      DO 450 J=1,3  
  450 P(N+4+I,J)=P(ILIN(I),J)   
      THERR=ACOS(2.*RLU(0)-1.)  
      PHIRR=PARU(2)*RLU(0)  
      CALL LUDBRB(N+5,N+4+IMAX,THERR,PHIRR,0D0,0D0,0D0) 
      DO 460 I=1,IMAX   
      DO 460 J=1,4  
  460 PK(I,J)=P(N+4+I,J)    
    
C...Calculate internal products.    
      IF(ISUB.EQ.22.OR.ISUB.EQ.23.OR.ISUB.EQ.25) THEN   
        DO 470 I1=IMIN,IMAX-1   
        DO 470 I2=I1+1,IMAX 
        HA(I1,I2)=SQRT((PK(I1,4)-PK(I1,3))*(PK(I2,4)+PK(I2,3))/ 
     &  (1E-20+PK(I1,1)**2+PK(I1,2)**2))*CMPLX(PK(I1,1),PK(I1,2))-  
     &  SQRT((PK(I1,4)+PK(I1,3))*(PK(I2,4)-PK(I2,3))/   
     &  (1E-20+PK(I2,1)**2+PK(I2,2)**2))*CMPLX(PK(I2,1),PK(I2,2))   
        HC(I1,I2)=CONJG(HA(I1,I2))  
        IF(I1.LE.2) HA(I1,I2)=CMPLX(0.,1.)*HA(I1,I2)    
        IF(I1.LE.2) HC(I1,I2)=CMPLX(0.,1.)*HC(I1,I2)    
        HA(I2,I1)=-HA(I1,I2)    
  470   HC(I2,I1)=-HC(I1,I2)    
      ENDIF 
      DO 480 I=1,2  
      DO 480 J=1,4  
  480 PK(I,J)=-PK(I,J)  
      DO 490 I1=IMIN,IMAX-1 
      DO 490 I2=I1+1,IMAX   
      PKK(I1,I2)=2.*(PK(I1,4)*PK(I2,4)-PK(I1,1)*PK(I2,1)-   
     &PK(I1,2)*PK(I2,2)-PK(I1,3)*PK(I2,3))  
  490 PKK(I2,I1)=PKK(I1,I2) 
    
      IF(IREF(IP,5).EQ.25) THEN 
C...Angular weight for H0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons 
        WT=16.*PKK(3,5)*PKK(4,6)    
        IF(IP.EQ.1) WTMAX=SH**2 
        IF(IP.GE.2) WTMAX=P(IREF(IP,6),5)**4    
    
      ELSEIF(ISUB.EQ.1) THEN    
        IF(KFA.NE.37) THEN  
C...Angular weight for gamma*/Z0 -> 2 quarks/leptons    
          EI=KCHG(IABS(MINT(15)),1)/3.  
          AI=SIGN(1.,EI+0.1)    
          VI=AI-4.*EI*XW    
          EF=KCHG(KFA,1)/3. 
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          GG=1. 
          GZ=1./(8.*XW*(1.-XW))*SH*(SH-SQMZ)/((SH-SQMZ)**2+GZMZ**2) 
          ZZ=1./(16.*XW*(1.-XW))**2*SH**2/((SH-SQMZ)**2+GZMZ**2)    
          IF(MSTP(43).EQ.1) THEN    
C...Only gamma* production included 
            GZ=0.   
            ZZ=0.   
          ELSEIF(MSTP(43).EQ.2) THEN    
C...Only Z0 production included 
            GG=0.   
            GZ=0.   
          ENDIF 
          ASYM=2.*(EI*AI*GZ*EF*AF+4.*VI*AI*ZZ*VF*AF)/(EI**2*GG*EF**2+   
     &    EI*VI*GZ*EF*VF+(VI**2+AI**2)*ZZ*(VF**2+AF**2))    
          WT=1.+ASYM*CTHE(JT)+CTHE(JT)**2   
          WTMAX=2.+ABS(ASYM)    
        ELSE    
C...Angular weight for gamma*/Z0 -> H+ + H- 
          WT=1.-CTHE(JT)**2 
          WTMAX=1.  
        ENDIF   
    
      ELSEIF(ISUB.EQ.2) THEN    
C...Angular weight for W+/- -> 2 quarks/leptons 
        WT=(1.+CTHE(JT))**2 
        WTMAX=4.    
    
      ELSEIF(ISUB.EQ.15.OR.ISUB.EQ.19) THEN 
C...Angular weight for f + fb -> gluon/gamma + Z0 ->    
C...-> gluon/gamma + 2 quarks/leptons   
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)* 
     &  (PKK(1,3)**2+PKK(2,4)**2)+((COUP(1,3)*COUP(3,4))**2+    
     &  (COUP(1,4)*COUP(3,3))**2)*(PKK(1,4)**2+PKK(2,3)**2) 
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*  
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2) 
    
      ELSEIF(ISUB.EQ.16.OR.ISUB.EQ.20) THEN 
C...Angular weight for f + fb' -> gluon/gamma + W+/- -> 
C...-> gluon/gamma + 2 quarks/leptons   
        WT=PKK(1,3)**2+PKK(2,4)**2  
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2 
    
      ELSEIF(ISUB.EQ.22) THEN   
C...Angular weight for f + fb -> Z0 + Z0 -> 4 quarks/leptons    
        S34=P(IREF(IP,IORD),5)**2   
        S56=P(IREF(IP,3-IORD),5)**2 
        TI=PKK(1,3)+PKK(1,4)+S34    
        UI=PKK(1,5)+PKK(1,6)+S56    
        WT=COUP(1,3)**4*((COUP(3,3)*COUP(5,3)*ABS(FGK(1,2,3,4,5,6)/ 
     &  TI+FGK(1,2,5,6,3,4)/UI))**2+(COUP(3,4)*COUP(5,3)*ABS(   
     &  FGK(1,2,4,3,5,6)/TI+FGK(1,2,5,6,4,3)/UI))**2+(COUP(3,3)*    
     &  COUP(5,4)*ABS(FGK(1,2,3,4,6,5)/TI+FGK(1,2,6,5,3,4)/UI))**2+ 
     &  (COUP(3,4)*COUP(5,4)*ABS(FGK(1,2,4,3,6,5)/TI+FGK(1,2,6,5,4,3)/  
     &  UI))**2)+COUP(1,4)**4*((COUP(3,3)*COUP(5,3)*ABS(    
     &  FGK(2,1,5,6,3,4)/TI+FGK(2,1,3,4,5,6)/UI))**2+(COUP(3,4)*    
     &  COUP(5,3)*ABS(FGK(2,1,6,5,3,4)/TI+FGK(2,1,3,4,6,5)/UI))**2+ 
     &  (COUP(3,3)*COUP(5,4)*ABS(FGK(2,1,5,6,4,3)/TI+FGK(2,1,4,3,5,6)/  
     &  UI))**2+(COUP(3,4)*COUP(5,4)*ABS(FGK(2,1,6,5,4,3)/TI+   
     &  FGK(2,1,4,3,6,5)/UI))**2)   
        WTMAX=4.*S34*S56*(COUP(1,3)**4+COUP(1,4)**4)*(COUP(3,3)**2+ 
     &  COUP(3,4)**2)*(COUP(5,3)**2+COUP(5,4)**2)*4.*(TI/UI+UI/TI+  
     &  2.*SH*(S34+S56)/(TI*UI)-S34*S56*(1./TI**2+1./UI**2))    
    
      ELSEIF(ISUB.EQ.23) THEN   
C...Angular weight for f + fb' -> Z0 + W +/- -> 4 quarks/leptons    
        D34=P(IREF(IP,IORD),5)**2   
        D56=P(IREF(IP,3-IORD),5)**2 
        DT=PKK(1,3)+PKK(1,4)+D34    
        DU=PKK(1,5)+PKK(1,6)+D56    
        CAWZ=COUP(2,3)/SNGL(DT)-2.*(1.-XW)*COUP(1,2)/(SH-SQMW)  
        CBWZ=COUP(1,3)/SNGL(DU)+2.*(1.-XW)*COUP(1,2)/(SH-SQMW)  
        WT=COUP(5,3)**2*ABS(CAWZ*FGK(1,2,3,4,5,6)+CBWZ* 
     &  FGK(1,2,5,6,3,4))**2+COUP(5,4)**2*ABS(CAWZ* 
     &  FGK(1,2,3,4,6,5)+CBWZ*FGK(1,2,6,5,3,4))**2  
        WTMAX=4.*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*(CAWZ**2*  
     &  DIGK(DT,DU)+CBWZ**2*DIGK(DU,DT)+CAWZ*CBWZ*DJGK(DT,DU))  
    
      ELSEIF(ISUB.EQ.24) THEN   
C...Angular weight for f + fb -> Z0 + H0 -> 2 quarks/leptons + H0   
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)* 
     &  PKK(1,3)*PKK(2,4)+((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)* 
     &  COUP(3,3))**2)*PKK(1,4)*PKK(2,3)    
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*  
     &  (PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4)) 
    
      ELSEIF(ISUB.EQ.25) THEN   
C...Angular weight for f + fb -> W+ + W- -> 4 quarks/leptons    
        D34=P(IREF(IP,IORD),5)**2   
        D56=P(IREF(IP,3-IORD),5)**2 
        DT=PKK(1,3)+PKK(1,4)+D34    
        DU=PKK(1,5)+PKK(1,6)+D56    
        CDWW=(COUP(1,3)*SQMZ/(SH-SQMZ)+COUP(1,2))/SH    
        CAWW=CDWW+0.5*(COUP(1,2)+1.)/SNGL(DT)   
        CBWW=CDWW+0.5*(COUP(1,2)-1.)/SNGL(DU)   
        CCWW=COUP(1,4)*SQMZ/(SH-SQMZ)/SH    
        WT=ABS(CAWW*FGK(1,2,3,4,5,6)-CBWW*FGK(1,2,5,6,3,4))**2+ 
     &  CCWW**2*ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))**2   
        WTMAX=4.*D34*D56*(CAWW**2*DIGK(DT,DU)+CBWW**2*DIGK(DU,DT)-CAWW* 
     &  CBWW*DJGK(DT,DU)+CCWW**2*(DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))) 
    
      ELSEIF(ISUB.EQ.26) THEN   
C...Angular weight for f + fb' -> W+/- + H0 -> 2 quarks/leptons + H0    
        WT=PKK(1,3)*PKK(2,4)    
        WTMAX=(PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))   
    
      ELSEIF(ISUB.EQ.30) THEN   
C...Angular weight for f + g -> f + Z0 -> f + 2 quarks/leptons  
        IF(K(ILIN(1),2).GT.0) WT=((COUP(1,3)*COUP(3,3))**2+ 
     &  (COUP(1,4)*COUP(3,4))**2)*(PKK(1,4)**2+PKK(3,5)**2)+    
     &  ((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*COUP(3,3))**2)*    
     &  (PKK(1,3)**2+PKK(4,5)**2)   
        IF(K(ILIN(1),2).LT.0) WT=((COUP(1,3)*COUP(3,3))**2+ 
     &  (COUP(1,4)*COUP(3,4))**2)*(PKK(1,3)**2+PKK(4,5)**2)+    
     &  ((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*COUP(3,3))**2)*    
     &  (PKK(1,4)**2+PKK(3,5)**2)   
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*  
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2) 
    
      ELSEIF(ISUB.EQ.31) THEN   
C...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons  
        IF(K(ILIN(1),2).GT.0) WT=PKK(1,4)**2+PKK(3,5)**2    
        IF(K(ILIN(1),2).LT.0) WT=PKK(1,3)**2+PKK(4,5)**2    
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2 
    
      ELSEIF(ISUB.EQ.141) THEN  
C...Angular weight for gamma*/Z0/Z'0 -> 2 quarks/leptons    
        EI=KCHG(IABS(MINT(15)),1)/3.    
        AI=SIGN(1.,EI+0.1)  
        VI=AI-4.*EI*XW  
        API=SIGN(1.,EI+0.1) 
        VPI=API-4.*EI*XW    
        EF=KCHG(KFA,1)/3.   
        AF=SIGN(1.,EF+0.1)  
        VF=AF-4.*EF*XW  
        APF=SIGN(1.,EF+0.1) 
        VPF=APF-4.*EF*XW    
        GG=1.   
        GZ=1./(8.*XW*(1.-XW))*SH*(SH-SQMZ)/((SH-SQMZ)**2+GZMZ**2)   
        GZP=1./(8.*XW*(1.-XW))*SH*(SH-SQMZP)/((SH-SQMZP)**2+GZMZP**2)   
        ZZ=1./(16.*XW*(1.-XW))**2*SH**2/((SH-SQMZ)**2+GZMZ**2)  
        ZZP=2./(16.*XW*(1.-XW))**2* 
     &  SH**2*((SH-SQMZ)*(SH-SQMZP)+GZMZ*GZMZP)/    
     &  (((SH-SQMZ)**2+GZMZ**2)*((SH-SQMZP)**2+GZMZP**2))   
        ZPZP=1./(16.*XW*(1.-XW))**2*SH**2/((SH-SQMZP)**2+GZMZP**2)  
        IF(MSTP(44).EQ.1) THEN  
C...Only gamma* production included 
          GZ=0. 
          GZP=0.    
          ZZ=0. 
          ZZP=0.    
          ZPZP=0.   
        ELSEIF(MSTP(44).EQ.2) THEN  
C...Only Z0 production included 
          GG=0. 
          GZ=0. 
          GZP=0.    
          ZZP=0.    
          ZPZP=0.   
        ELSEIF(MSTP(44).EQ.3) THEN  
C...Only Z'0 production included    
          GG=0. 
          GZ=0. 
          GZP=0.    
          ZZ=0. 
          ZZP=0.    
        ELSEIF(MSTP(44).EQ.4) THEN  
C...Only gamma*/Z0 production included  
          GZP=0.    
          ZZP=0.    
          ZPZP=0.   
        ELSEIF(MSTP(44).EQ.5) THEN  
C...Only gamma*/Z'0 production included 
          GZ=0. 
          ZZ=0. 
          ZZP=0.    
        ELSEIF(MSTP(44).EQ.6) THEN  
C...Only Z0/Z'0 production included 
          GG=0. 
          GZ=0. 
          GZP=0.    
        ENDIF   
        ASYM=2.*(EI*AI*GZ*EF*AF+EI*API*GZP*EF*APF+4.*VI*AI*ZZ*VF*AF+    
     &  (VI*API+VPI*AI)*ZZP*(VF*APF+VPF*AF)+4.*VPI*API*ZPZP*VPF*APF)/   
     &  (EI**2*GG*EF**2+EI*VI*GZ*EF*VF+EI*VPI*GZP*EF*VPF+   
     &  (VI**2+AI**2)*ZZ*(VF**2+AF**2)+(VI*VPI+AI*API)*ZZP* 
     &  (VF*VPF+AF*APF)+(VPI**2+API**2)*ZPZP*(VPF**2+APF**2))   
        WT=1.+ASYM*CTHE(JT)+CTHE(JT)**2 
        WTMAX=2.+ABS(ASYM)  
    
      ELSE  
        WT=1.   
        WTMAX=1.    
      ENDIF 
C...Obtain correct angular distribution by rejection techniques.    
      IF(WT.LT.RLU(0)*WTMAX) GOTO 420   
    
C...Construct massive four-vectors using angles chosen. Mark decayed    
C...resonances, add documentation lines. Shower evolution.  
  500 DO 520 JT=1,JTMAX 
      IF(KDCY(JT).EQ.0) GOTO 520    
      ID=IREF(IP,JT)    
      CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,ACOS(CTHE(JT)),PHI(JT),   
     &DBLE(P(ID,1)/P(ID,4)),DBLE(P(ID,2)/P(ID,4)),DBLE(P(ID,3)/P(ID,4)))    
      K(ID,1)=K(ID,1)+10    
      K(ID,4)=NSD(JT)+1 
      K(ID,5)=NSD(JT)+2 
      IDOC=MINT(83)+MINT(4) 
      DO 510 I=NSD(JT)+1,NSD(JT)+2  
      MINT(4)=MINT(4)+1 
      I1=MINT(83)+MINT(4)   
      K(I,3)=I1 
      K(I1,1)=21    
      K(I1,2)=K(I,2)    
      K(I1,3)=IREF(IP,JT+2) 
      DO 510 J=1,5  
  510 P(I1,J)=P(I,J)    
      IF(JTMAX.EQ.1) THEN   
        MINT(7)=MINT(83)+6+2*ISET(ISUB) 
        MINT(8)=MINT(83)+7+2*ISET(ISUB) 
      ENDIF 
      IF(MSTP(71).GE.1.AND.KDCY(JT).EQ.1) CALL LUSHOW(NSD(JT)+1,    
     &NSD(JT)+2,P(ID,5))    
    
C...Check if new resonances were produced, loop back if needed. 
      IF(KDCY(JT).NE.3) GOTO 520    
      NP=NP+1   
      IREF(NP,1)=NSD(JT)+1  
      IREF(NP,2)=NSD(JT)+2  
      IREF(NP,3)=IDOC+1 
      IREF(NP,4)=IDOC+2 
      IREF(NP,5)=K(IREF(IP,JT),2)   
      IREF(NP,6)=IREF(IP,JT)    
  520 CONTINUE  
  530 IF(IP.LT.NP) GOTO 100 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYDIFF 
    
C...Handles diffractive and elastic scattering. 
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
    
C...Reset K, P and V vectors. Store incoming particles. 
      DO 100 JT=1,MSTP(126)+10  
      I=MINT(83)+JT 
      DO 100 J=1,5  
      K(I,J)=0  
      P(I,J)=0. 
  100 V(I,J)=0. 
      N=MINT(84)    
      MINT(3)=0 
      MINT(21)=0    
      MINT(22)=0    
      MINT(23)=0    
      MINT(24)=0    
      MINT(4)=4 
      DO 110 JT=1,2 
      I=MINT(83)+JT 
      K(I,1)=21 
      K(I,2)=MINT(10+JT)    
      P(I,5)=VINT(2+JT) 
      P(I,3)=VINT(5)*(-1)**(JT+1)   
  110 P(I,4)=SQRT(P(I,3)**2+P(I,5)**2)  
      MINT(6)=2 
    
C...Subprocess; kinematics. 
      ISUB=MINT(1)  
      SQLAM=(VINT(2)-VINT(63)-VINT(64))**2-4.*VINT(63)*VINT(64) 
      PZ=SQRT(SQLAM)/(2.*VINT(1))   
      DO 150 JT=1,2 
      I=MINT(83)+JT 
      PE=(VINT(2)+VINT(62+JT)-VINT(65-JT))/(2.*VINT(1)) 
    
C...Elastically scattered particle. 
      IF(MINT(16+JT).LE.0) THEN 
        N=N+1   
        K(N,1)=1    
        K(N,2)=K(I,2)   
        K(N,3)=I+2  
        P(N,3)=PZ*(-1)**(JT+1)  
        P(N,4)=PE   
        P(N,5)=P(I,5)   
    
C...Diffracted particle: valence quark kicked out.  
      ELSEIF(MSTP(101).EQ.1) THEN   
        N=N+2   
        K(N-1,1)=2  
        K(N,1)=1    
        K(N-1,3)=I+2    
        K(N,3)=I+2  
        CALL PYSPLI(K(I,2),21,K(N,2),K(N-1,2))  
        P(N-1,5)=ULMASS(K(N-1,2))   
        P(N,5)=ULMASS(K(N,2))   
        SQLAM=(VINT(62+JT)-P(N-1,5)**2-P(N,5)**2)**2-   
     &  4.*P(N-1,5)**2*P(N,5)**2    
        P(N-1,3)=(PE*SQRT(SQLAM)+PZ*(VINT(62+JT)+P(N-1,5)**2-   
     &  P(N,5)**2))/(2.*VINT(62+JT))*(-1)**(JT+1)   
        P(N-1,4)=SQRT(P(N-1,3)**2+P(N-1,5)**2)  
        P(N,3)=PZ*(-1)**(JT+1)-P(N-1,3) 
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)    
    
C...Diffracted particle: gluon kicked out.  
      ELSE  
        N=N+3   
        K(N-2,1)=2  
        K(N-1,1)=2  
        K(N,1)=1    
        K(N-2,3)=I+2    
        K(N-1,3)=I+2    
        K(N,3)=I+2  
        CALL PYSPLI(K(I,2),21,K(N,2),K(N-2,2))  
        K(N-1,2)=21 
        P(N-2,5)=ULMASS(K(N-2,2))   
        P(N-1,5)=0. 
        P(N,5)=ULMASS(K(N,2))   
C...Energy distribution for particle into two jets. 
  120   IMB=1   
        IF(MOD(K(I,2)/1000,10).NE.0) IMB=2  
        CHIK=PARP(92+2*IMB) 
        IF(MSTP(92).LE.1) THEN  
          IF(IMB.EQ.1) CHI=RLU(0)   
          IF(IMB.EQ.2) CHI=1.-SQRT(RLU(0))  
        ELSEIF(MSTP(92).EQ.2) THEN  
          CHI=1.-RLU(0)**(1./(1.+CHIK)) 
        ELSEIF(MSTP(92).EQ.3) THEN  
          CUT=2.*0.3/VINT(1)    
  130     CHI=RLU(0)**2 
          IF((CHI**2/(CHI**2+CUT**2))**0.25*(1.-CHI)**CHIK.LT.  
     &    RLU(0)) GOTO 130  
        ELSE    
          CUT=2.*0.3/VINT(1)    
          CUTR=(1.+SQRT(1.+CUT**2))/CUT 
  140     CHIR=CUT*CUTR**RLU(0) 
          CHI=(CHIR**2-CUT**2)/(2.*CHIR)    
          IF((1.-CHI)**CHIK.LT.RLU(0)) GOTO 140 
        ENDIF   
        IF(CHI.LT.P(N,5)**2/VINT(62+JT).OR.CHI.GT.1.-P(N-2,5)**2/   
     &  VINT(62+JT)) GOTO 120   
        SQM=P(N-2,5)**2/(1.-CHI)+P(N,5)**2/CHI  
        IF((SQRT(SQM)+PARJ(32))**2.GE.VINT(62+JT)) GOTO 120 
        PZI=(PE*(VINT(62+JT)-SQM)+PZ*(VINT(62+JT)+SQM))/    
     &  (2.*VINT(62+JT))    
        PEI=SQRT(PZI**2+SQM)    
        PQQP=(1.-CHI)*(PEI+PZI) 
        P(N-2,3)=0.5*(PQQP-P(N-2,5)**2/PQQP)*(-1)**(JT+1)   
        P(N-2,4)=SQRT(P(N-2,3)**2+P(N-2,5)**2)  
        P(N-1,3)=(PZ-PZI)*(-1)**(JT+1)  
        P(N-1,4)=ABS(P(N-1,3))  
        P(N,3)=PZI*(-1)**(JT+1)-P(N-2,3)    
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)    
      ENDIF 
    
C...Documentation lines.    
      K(I+2,1)=21   
      IF(MINT(16+JT).EQ.0) K(I+2,2)=MINT(10+JT) 
      IF(MINT(16+JT).NE.0) K(I+2,2)=10*(MINT(10+JT)/10) 
      K(I+2,3)=I    
      P(I+2,3)=PZ*(-1)**(JT+1)  
      P(I+2,4)=PE   
      P(I+2,5)=SQRT(VINT(62+JT))    
  150 CONTINUE  
    
C...Rotate outgoing partons/particles using cos(theta). 
      CALL LUDBRB(MINT(83)+3,N,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0) 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYFRAM(IFRAME) 
    
C...Performs transformations between different coordinate frames.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
    
      IF(IFRAME.LT.1.OR.IFRAME.GT.2) THEN   
        WRITE(MSTU(11),1000) IFRAME,MINT(6) 
        RETURN  
      ENDIF 
      IF(IFRAME.EQ.MINT(6)) RETURN  
    
      IF(MINT(6).EQ.1) THEN 
C...Transform from fixed target or user specified frame to  
C...CM-frame of incoming particles. 
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))  
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)   
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)   
        MINT(6)=2   
    
      ELSE  
C...Transform from particle CM-frame to fixed target or user specified  
C...frame.  
        CALL LUROBO(VINT(6),VINT(7),VINT(8),VINT(9),VINT(10))   
        MINT(6)=1   
      ENDIF 
      MSTI(6)=MINT(6)   
    
 1000 FORMAT(1X,'Error: illegal values in subroutine PYFRAM.',1X,   
     &'No transformation performed.'/1X,'IFRAME =',1X,I5,'; MINT(6) =', 
     &1X,I5)    
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYWIDT(KFLR,RMAS,WDTP,WDTE)    
    
C...Calculates full and partial widths of resonances.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      DIMENSION WDTP(0:40),WDTE(0:40,0:5)   
    
C...Some common constants.  
      KFLA=IABS(KFLR)   
      SQM=RMAS**2   
      AS=ULALPS(SQM)    
      AEM=PARU(101) 
      XW=PARU(102)  
      RADC=1.+AS/PARU(1)    
    
C...Reset width information.    
      DO 100 I=0,40 
      WDTP(I)=0.    
      DO 100 J=0,5  
  100 WDTE(I,J)=0.  
    
      IF(KFLA.EQ.21) THEN   
C...QCD:    
        DO 110 I=1,MDCY(21,3)   
        IDC=I+MDCY(21,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 110  
        IF(I.LE.8) THEN 
C...QCD -> q + qb   
          WDTP(I)=(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  110   CONTINUE    
    
      ELSEIF(KFLA.EQ.23) THEN   
C...Z0: 
        IF(MINT(61).EQ.1) THEN  
          EI=KCHG(IABS(MINT(15)),1)/3.  
          AI=SIGN(1.,EI)    
          VI=AI-4.*EI*XW    
          SQMZ=PMAS(23,1)**2    
          GZMZ=PMAS(23,2)*PMAS(23,1)    
          GGI=EI**2 
          GZI=EI*VI/(8.*XW*(1.-XW))*SQM*(SQM-SQMZ)/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          ZZI=(VI**2+AI**2)/(16.*XW*(1.-XW))**2*SQM**2/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          IF(MSTP(43).EQ.1) THEN    
C...Only gamma* production included 
            GZI=0.  
            ZZI=0.  
          ELSEIF(MSTP(43).EQ.2) THEN    
C...Only Z0 production included 
            GGI=0.  
            GZI=0.  
          ENDIF 
        ELSEIF(MINT(61).EQ.2) THEN  
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(114)=0.  
        ENDIF   
        DO 120 I=1,MDCY(23,3)   
        IDC=I+MDCY(23,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 120  
        IF(I.LE.8) THEN 
C...Z0 -> q + qb    
          EF=KCHG(I,1)/3.   
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=3.*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=3.*((GGI*EF**2+GZI*EF*VF+ZZI*VF**2)*    
     &      (1.+2.*RM1)+ZZI*AF**2*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=3.*EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            GZF=3.*EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            ZZF=3.*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ENDIF 
          WID2=1.   
        ELSEIF(I.LE.16) THEN    
C...Z0 -> l+ + l-, nu + nub 
          EF=KCHG(I+2,1)/3. 
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          WDTP(I)=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*    
     &    SQRT(MAX(0.,1.-4.*RM1))   
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=((GGI*EF**2+GZI*EF*VF+ZZI*VF**2)*   
     &      (1.+2.*RM1)+ZZI*AF**2*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            GZF=EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            ZZF=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ENDIF 
          WID2=1.   
        ELSE    
C...Z0 -> H+ + H-   
          CF=2.*(1.-2.*XW)  
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=0.25*CF**2*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))  
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=0.25*(GGI+GZI*CF+ZZI*CF**2)*(1.-4.*RM1)*    
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=0.25*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))    
            GZF=0.25*CF*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1)) 
            ZZF=0.25*CF**2*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))  
          ENDIF 
          WID2=WIDS(37,1)   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
          VINT(111)=VINT(111)+GGF*WID2  
          VINT(112)=VINT(112)+GZF*WID2  
          VINT(114)=VINT(114)+ZZF*WID2  
        ENDIF   
  120   CONTINUE    
        IF(MSTP(43).EQ.1) THEN  
C...Only gamma* production included 
          VINT(112)=0.  
          VINT(114)=0.  
        ELSEIF(MSTP(43).EQ.2) THEN  
C...Only Z0 production included 
          VINT(111)=0.  
          VINT(112)=0.  
        ENDIF   
    
      ELSEIF(KFLA.EQ.24) THEN   
C...W+/-:   
        DO 130 I=1,MDCY(24,3)   
        IDC=I+MDCY(24,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 130  
        IF(I.LE.16) THEN    
C...W+/- -> q + qb' 
          WDTP(I)=3.*(2.-RM1-RM2-(RM1-RM2)**2)* 
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))* 
     &    VCKM((I-1)/4+1,MOD(I-1,4)+1)*RADC 
          WID2=1.   
        ELSE    
C...W+/- -> l+/- + nu   
          WDTP(I)=(2.-RM1-RM2-(RM1-RM2)**2)*    
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))  
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  130   CONTINUE    
    
      ELSEIF(KFLA.EQ.25) THEN   
C...H0: 
        DO 170 I=1,MDCY(25,3)   
        IDC=I+MDCY(25,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 170  
        IF(I.LE.8) THEN 
C...H0 -> q + qb    
          WDTP(I)=3.*RM1*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
          WID2=1.   
        ELSEIF(I.LE.12) THEN    
C...H0 -> l+ + l-   
          WDTP(I)=RM1*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
          WID2=1.   
        ELSEIF(I.EQ.13) THEN    
C...H0 -> g + g; quark loop contribution only   
          ETARE=0.  
          ETAIM=0.  
          DO 140 J=1,2*MSTP(1)  
          EPS=(2.*PMAS(J,1)/RMAS)**2    
          IF(EPS.LE.1.) THEN    
            IF(EPS.GT.1.E-4) THEN   
              ROOT=SQRT(1.-EPS) 
              RLN=LOG((1.+ROOT)/(1.-ROOT))  
            ELSE    
              RLN=LOG(4./EPS-2.)    
            ENDIF   
            PHIRE=0.25*(RLN**2-PARU(1)**2)  
            PHIIM=0.5*PARU(1)*RLN   
          ELSE  
            PHIRE=-(ASIN(1./SQRT(EPS)))**2  
            PHIIM=0.    
          ENDIF 
          ETARE=ETARE+0.5*EPS*(1.+(EPS-1.)*PHIRE)   
          ETAIM=ETAIM+0.5*EPS*(EPS-1.)*PHIIM    
  140     CONTINUE  
          ETA2=ETARE**2+ETAIM**2    
          WDTP(I)=(AS/PARU(1))**2*ETA2  
          WID2=1.   
        ELSEIF(I.EQ.14) THEN    
C...H0 -> gamma + gamma; quark, charged lepton and W loop contributions 
          ETARE=0.  
          ETAIM=0.  
          DO 150 J=1,3*MSTP(1)+1    
          IF(J.LE.2*MSTP(1)) THEN   
            EJ=KCHG(J,1)/3. 
            EPS=(2.*PMAS(J,1)/RMAS)**2  
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            JL=2*(J-2*MSTP(1))-1    
            EJ=KCHG(10+JL,1)/3. 
            EPS=(2.*PMAS(10+JL,1)/RMAS)**2  
          ELSE  
            EPS=(2.*PMAS(24,1)/RMAS)**2 
          ENDIF 
          IF(EPS.LE.1.) THEN    
            IF(EPS.GT.1.E-4) THEN   
              ROOT=SQRT(1.-EPS) 
              RLN=LOG((1.+ROOT)/(1.-ROOT))  
            ELSE    
              RLN=LOG(4./EPS-2.)    
            ENDIF   
            PHIRE=0.25*(RLN**2-PARU(1)**2)  
            PHIIM=0.5*PARU(1)*RLN   
          ELSE  
            PHIRE=-(ASIN(1./SQRT(EPS)))**2  
            PHIIM=0.    
          ENDIF 
          IF(J.LE.2*MSTP(1)) THEN   
            ETARE=ETARE+0.5*3.*EJ**2*EPS*(1.+(EPS-1.)*PHIRE)    
            ETAIM=ETAIM+0.5*3.*EJ**2*EPS*(EPS-1.)*PHIIM 
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            ETARE=ETARE+0.5*EJ**2*EPS*(1.+(EPS-1.)*PHIRE)   
            ETAIM=ETAIM+0.5*EJ**2*EPS*(EPS-1.)*PHIIM    
          ELSE  
            ETARE=ETARE-0.5-0.75*EPS*(1.+(EPS-2.)*PHIRE)    
            ETAIM=ETAIM+0.75*EPS*(EPS-2.)*PHIIM 
          ENDIF 
  150     CONTINUE  
          ETA2=ETARE**2+ETAIM**2    
          WDTP(I)=(AEM/PARU(1))**2*0.5*ETA2 
          WID2=1.   
        ELSEIF(I.EQ.15) THEN    
C...H0 -> gamma + Z0; quark, charged lepton and W loop contributions    
          ETARE=0.  
          ETAIM=0.  
          DO 160 J=1,3*MSTP(1)+1    
          IF(J.LE.2*MSTP(1)) THEN   
            EJ=KCHG(J,1)/3. 
            AJ=SIGN(1.,EJ+0.1)  
            VJ=AJ-4.*EJ*XW  
            EPS=(2.*PMAS(J,1)/RMAS)**2  
            EPSP=(2.*PMAS(J,1)/PMAS(23,1))**2   
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            JL=2*(J-2*MSTP(1))-1    
            EJ=KCHG(10+JL,1)/3. 
            AJ=SIGN(1.,EJ+0.1)  
            VJ=AI-4.*EJ*XW  
            EPS=(2.*PMAS(10+JL,1)/RMAS)**2  
            EPSP=(2.*PMAS(10+JL,1)/PMAS(23,1))**2   
          ELSE  
            EPS=(2.*PMAS(24,1)/RMAS)**2 
            EPSP=(2.*PMAS(24,1)/PMAS(23,1))**2  
          ENDIF 
          IF(EPS.LE.1.) THEN    
            ROOT=SQRT(1.-EPS)   
            IF(EPS.GT.1.E-4) THEN   
              RLN=LOG((1.+ROOT)/(1.-ROOT))  
            ELSE    
              RLN=LOG(4./EPS-2.)    
            ENDIF   
            PHIRE=0.25*(RLN**2-PARU(1)**2)  
            PHIIM=0.5*PARU(1)*RLN   
            PSIRE=-(1.+0.5*ROOT*RLN)    
            PSIIM=0.5*PARU(1)*ROOT  
          ELSE  
            PHIRE=-(ASIN(1./SQRT(EPS)))**2  
            PHIIM=0.    
            PSIRE=-(1.+SQRT(EPS-1.)*ASIN(1./SQRT(EPS))) 
            PSIIM=0.    
          ENDIF 
          IF(EPSP.LE.1.) THEN   
            ROOT=SQRT(1.-EPSP)  
            IF(EPSP.GT.1.E-4) THEN  
              RLN=LOG((1.+ROOT)/(1.-ROOT))  
            ELSE    
              RLN=LOG(4./EPSP-2.)   
            ENDIF   
            PHIREP=0.25*(RLN**2-PARU(1)**2) 
            PHIIMP=0.5*PARU(1)*RLN  
            PSIREP=-(1.+0.5*ROOT*RLN)   
            PSIIMP=0.5*PARU(1)*ROOT 
          ELSE  
            PHIREP=-(ASIN(1./SQRT(EPSP)))**2    
            PHIIMP=0.   
            PSIREP=-(1.+SQRT(EPSP-1.)*ASIN(1./SQRT(EPSP)))  
            PSIIMP=0.   
          ENDIF 
          FXYRE=EPS*EPSP/(8.*(EPS-EPSP))*(1.-EPS*EPSP/(EPS-EPSP)*(PHIRE-    
     &    PHIREP)+2.*EPS/(EPS-EPSP)*(PSIRE-PSIREP)) 
          FXYIM=EPS*EPSP/(8.*(EPS-EPSP))*(-EPS*EPSP/(EPS-EPSP)*(PHIIM-  
     &    PHIIMP)+2.*EPS/(EPS-EPSP)*(PSIIM-PSIIMP)) 
          F1RE=EPS*EPSP/(2.*(EPS-EPSP))*(PHIRE-PHIREP)  
          F1IM=EPS*EPSP/(2.*(EPS-EPSP))*(PHIIM-PHIIMP)  
          IF(J.LE.2*MSTP(1)) THEN   
            ETARE=ETARE-3.*EJ*VJ*(FXYRE-0.25*F1RE)  
            ETAIM=ETAIM-3.*EJ*VJ*(FXYIM-0.25*F1IM)  
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            ETARE=ETARE-EJ*VJ*(FXYRE-0.25*F1RE) 
            ETAIM=ETAIM-EJ*VJ*(FXYIM-0.25*F1IM) 
          ELSE  
            ETARE=ETARE-SQRT(1.-XW)*(((1.+2./EPS)*XW/SQRT(1.-XW)-   
     &      (5.+2./EPS))*FXYRE+(3.-XW/SQRT(1.-XW))*F1RE)    
            ETAIM=ETAIM-SQRT(1.-XW)*(((1.+2./EPS)*XW/SQRT(1.-XW)-   
     &      (5.+2./EPS))*FXYIM+(3.-XW/SQRT(1.-XW))*F1IM)    
          ENDIF 
  160     CONTINUE  
          ETA2=ETARE**2+ETAIM**2    
          WDTP(I)=(AEM/PARU(1))**2*(1.-(PMAS(23,1)/RMAS)**2)**3/XW*ETA2 
          WID2=WIDS(23,2)   
        ELSE    
C...H0 -> Z0 + Z0, W+ + W-  
          WDTP(I)=(1.-4.*RM1+12.*RM1**2)*SQRT(MAX(0.,1.-4.*RM1))/   
     &    (2.*(18-I))   
          WID2=WIDS(7+I,1)  
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  170   CONTINUE    
    
      ELSEIF(KFLA.EQ.32) THEN   
C...Z'0:    
        IF(MINT(61).EQ.1) THEN  
          EI=KCHG(IABS(MINT(15)),1)/3.  
          AI=SIGN(1.,EI)    
          VI=AI-4.*EI*XW    
          SQMZ=PMAS(23,1)**2    
          GZMZ=PMAS(23,2)*PMAS(23,1)    
          API=SIGN(1.,EI)   
          VPI=API-4.*EI*XW  
          SQMZP=PMAS(32,1)**2   
          GZPMZP=PMAS(32,2)*PMAS(32,1)  
          GGI=EI**2 
          GZI=EI*VI/(8.*XW*(1.-XW))*SQM*(SQM-SQMZ)/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          GZPI=EI*VPI/(8.*XW*(1.-XW))*SQM*(SQM-SQMZP)/  
     &    ((SQM-SQMZP)**2+GZPMZP**2)    
          ZZI=(VI**2+AI**2)/(16.*XW*(1.-XW))**2*SQM**2/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          ZZPI=2.*(VI*VPI+AI*API)/(16.*XW*(1.-XW))**2*  
     &    SQM**2*((SQM-SQMZ)*(SQM-SQMZP)+GZMZ*GZPMZP)/  
     &    (((SQM-SQMZ)**2+GZMZ**2)*((SQM-SQMZP)**2+GZPMZP**2))  
          ZPZPI=(VPI**2+API**2)/(16.*XW*(1.-XW))**2*SQM**2/ 
     &    ((SQM-SQMZP)**2+GZPMZP**2)    
          IF(MSTP(44).EQ.1) THEN    
C...Only gamma* production included 
            GZI=0.  
            GZPI=0. 
            ZZI=0.  
            ZZPI=0. 
            ZPZPI=0.    
          ELSEIF(MSTP(44).EQ.2) THEN    
C...Only Z0 production included 
            GGI=0.  
            GZI=0.  
            GZPI=0. 
            ZZPI=0. 
            ZPZPI=0.    
          ELSEIF(MSTP(44).EQ.3) THEN    
C...Only Z'0 production included    
            GGI=0.  
            GZI=0.  
            GZPI=0. 
            ZZI=0.  
            ZZPI=0. 
          ELSEIF(MSTP(44).EQ.4) THEN    
C...Only gamma*/Z0 production included  
            GZPI=0. 
            ZZPI=0. 
            ZPZPI=0.    
          ELSEIF(MSTP(44).EQ.5) THEN    
C...Only gamma*/Z'0 production included 
            GZI=0.  
            ZZI=0.  
            ZZPI=0. 
          ELSEIF(MSTP(44).EQ.6) THEN    
C...Only Z0/Z'0 production included 
            GGI=0.  
            GZI=0.  
            GZPI=0. 
          ENDIF 
        ELSEIF(MINT(61).EQ.2) THEN  
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ENDIF   
        DO 180 I=1,MDCY(32,3)   
        IDC=I+MDCY(32,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 180  
        IF(I.LE.8) THEN 
C...Z'0 -> q + qb   
          EF=KCHG(I,1)/3.   
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          APF=SIGN(1.,EF+0.1)   
          VPF=APF-4.*EF*XW  
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=3.*(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=3.*((GGI*EF**2+GZI*EF*VF+GZPI*EF*VPF+ZZI*VF**2+ 
     &      ZZPI*VF*VPF+ZPZPI*VPF**2)*(1.+2.*RM1)+(ZZI*AF**2+   
     &      ZZPI*AF*APF+ZPZPI*APF**2)*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=3.*EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            GZF=3.*EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            GZPF=3.*EF*VPF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC 
            ZZF=3.*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
            ZZPF=3.*(VF*VPF*(1.+2.*RM1)+AF*APF*(1.-4.*RM1))*    
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
            ZPZPF=3.*(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ENDIF 
          WID2=1.   
        ELSE    
C...Z'0 -> l+ + l-, nu + nub    
          EF=KCHG(I+2,1)/3. 
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          APF=SIGN(1.,EF+0.1)   
          VPF=API-4.*EF*XW  
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*    
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=((GGI*EF**2+GZI*EF*VF+GZPI*EF*VPF+ZZI*VF**2+    
     &      ZZPI*VF*VPF+ZPZPI*VPF**2)*(1.+2.*RM1)+(ZZI*AF**2+   
     &      ZZPI*AF*APF+ZPZPI*APF**2)*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            GZF=EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            GZPF=EF*VPF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1)) 
            ZZF=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
            ZZPF=(VF*VPF*(1.+2.*RM1)+AF*APF*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1)) 
            ZPZPF=(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ENDIF 
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
          VINT(111)=VINT(111)+GGF   
          VINT(112)=VINT(112)+GZF   
          VINT(113)=VINT(113)+GZPF  
          VINT(114)=VINT(114)+ZZF   
          VINT(115)=VINT(115)+ZZPF  
          VINT(116)=VINT(116)+ZPZPF 
        ENDIF   
  180   CONTINUE    
        IF(MSTP(44).EQ.1) THEN  
C...Only gamma* production included 
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ELSEIF(MSTP(44).EQ.2) THEN  
C...Only Z0 production included 
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ELSEIF(MSTP(44).EQ.3) THEN  
C...Only Z'0 production included    
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
        ELSEIF(MSTP(44).EQ.4) THEN  
C...Only gamma*/Z0 production included  
          VINT(113)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ELSEIF(MSTP(44).EQ.5) THEN  
C...Only gamma*/Z'0 production included 
          VINT(112)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
        ELSEIF(MSTP(44).EQ.6) THEN  
C...Only Z0/Z'0 production included 
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
        ENDIF   
    
      ELSEIF(KFLA.EQ.37) THEN   
C...H+/-:   
        DO 190 I=1,MDCY(37,3)   
        IDC=I+MDCY(37,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 190  
        IF(I.LE.4) THEN 
C...H+/- -> q + qb' 
          WDTP(I)=3.*((RM1*PARU(121)+RM2/PARU(121))*    
     &    (1.-RM1-RM2)-4.*RM1*RM2)* 
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))*RADC 
          WID2=1.   
        ELSE    
C...H+/- -> l+/- + nu   
          WDTP(I)=((RM1*PARU(121)+RM2/PARU(121))*   
     &    (1.-RM1-RM2)-4.*RM1*RM2)* 
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))  
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  190   CONTINUE    
    
      ELSEIF(KFLA.EQ.40) THEN   
C...R:  
        DO 200 I=1,MDCY(40,3)   
        IDC=I+MDCY(40,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 200  
        IF(I.LE.4) THEN 
C...R -> q + qb'    
          WDTP(I)=3.*RADC   
          WID2=1.   
        ELSE    
C...R -> l+ + l'-   
          WDTP(I)=1.    
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  200   CONTINUE    
    
      ENDIF 
      MINT(61)=0    
    
      RETURN    
      END   
    
C***********************************************************************    
    
      SUBROUTINE PYKLIM(ILIM)   
    
C...Checks generated variables against pre-set kinematical limits;  
C...also calculates limits on variables used in generation. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
    
C...Common kinematical expressions. 
      ISUB=MINT(1)  
      IF(ISUB.EQ.96) GOTO 110   
      SQM3=VINT(63) 
      SQM4=VINT(64) 
      IF(ILIM.NE.1) THEN    
        TAU=VINT(21)    
        RM3=SQM3/(TAU*VINT(2))  
        RM4=SQM4/(TAU*VINT(2))  
        BE34=SQRT((1.-RM3-RM4)**2-4.*RM3*RM4)   
      ENDIF 
      PTHMIN=CKIN(3)    
      IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) PTHMIN=MAX(CKIN(3),CKIN(5))  
    
      IF(ILIM.EQ.0) THEN    
C...Check generated values of tau, y*, cos(theta-hat), and tau' against 
C...pre-set kinematical limits. 
        YST=VINT(22)    
        CTH=VINT(23)    
        TAUP=VINT(26)   
        IF(ISET(ISUB).LE.2) THEN    
          X1=SQRT(TAU)*EXP(YST) 
          X2=SQRT(TAU)*EXP(-YST)    
        ELSE    
          X1=SQRT(TAUP)*EXP(YST)    
          X2=SQRT(TAUP)*EXP(-YST)   
        ENDIF   
        XF=X1-X2    
        IF(TAU*VINT(2).LT.CKIN(1)**2) MINT(51)=1    
        IF(CKIN(2).GE.0..AND.TAU*VINT(2).GT.CKIN(2)**2) MINT(51)=1  
        IF(X1.LT.CKIN(21).OR.X1.GT.CKIN(22)) MINT(51)=1 
        IF(X2.LT.CKIN(23).OR.X2.GT.CKIN(24)) MINT(51)=1 
        IF(XF.LT.CKIN(25).OR.XF.GT.CKIN(26)) MINT(51)=1 
        IF(YST.LT.CKIN(7).OR.YST.GT.CKIN(8)) MINT(51)=1 
        IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN 
          PTH=0.5*BE34*SQRT(TAU*VINT(2)*(1.-CTH**2))    
          Y3=YST+0.5*LOG((1.+RM3-RM4+BE34*CTH)/(1.+RM3-RM4-BE34*CTH))   
          Y4=YST+0.5*LOG((1.+RM4-RM3-BE34*CTH)/(1.+RM4-RM3+BE34*CTH))   
          YLARGE=MAX(Y3,Y4) 
          YSMALL=MIN(Y3,Y4) 
          ETALAR=10.    
          ETASMA=-10.   
          STH=SQRT(1.-CTH**2)   
          IF(STH.LT.1.E-6) GOTO 100 
          EXPET3=((1.+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH+    
     &    SQRT(((1.+RM3-RM4)*COSH(YST)+BE34*SINH(YST)*CTH)**2-4.*RM3))/ 
     &    (BE34*STH)    
          EXPET4=((1.-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH+    
     &    SQRT(((1.-RM3+RM4)*COSH(YST)-BE34*SINH(YST)*CTH)**2-4.*RM4))/ 
     &    (BE34*STH)    
          ETA3=LOG(MIN(1.E10,MAX(1.E-10,EXPET3)))   
          ETA4=LOG(MIN(1.E10,MAX(1.E-10,EXPET4)))   
          ETALAR=MAX(ETA3,ETA4) 
          ETASMA=MIN(ETA3,ETA4) 
  100     CTS3=((1.+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH)/ 
     &    SQRT(((1.+RM3-RM4)*COSH(YST)+BE34*SINH(YST)*CTH)**2-4.*RM3)   
          CTS4=((1.-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH)/ 
     &    SQRT(((1.-RM3+RM4)*COSH(YST)-BE34*SINH(YST)*CTH)**2-4.*RM4)   
          CTSLAR=MAX(CTS3,CTS4) 
          CTSSMA=MIN(CTS3,CTS4) 
          IF(PTH.LT.PTHMIN) MINT(51)=1  
          IF(CKIN(4).GE.0..AND.PTH.GT.CKIN(4)) MINT(51)=1   
          IF(YLARGE.LT.CKIN(9).OR.YLARGE.GT.CKIN(10)) MINT(51)=1    
          IF(YSMALL.LT.CKIN(11).OR.YSMALL.GT.CKIN(12)) MINT(51)=1   
          IF(ETALAR.LT.CKIN(13).OR.ETALAR.GT.CKIN(14)) MINT(51)=1   
          IF(ETASMA.LT.CKIN(15).OR.ETASMA.GT.CKIN(16)) MINT(51)=1   
          IF(CTSLAR.LT.CKIN(17).OR.CTSLAR.GT.CKIN(18)) MINT(51)=1   
          IF(CTSSMA.LT.CKIN(19).OR.CTSSMA.GT.CKIN(20)) MINT(51)=1   
          IF(CTH.LT.CKIN(27).OR.CTH.GT.CKIN(28)) MINT(51)=1 
        ENDIF   
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          IF(TAUP*VINT(2).LT.CKIN(31)**2) MINT(51)=1    
          IF(CKIN(32).GE.0..AND.TAUP*VINT(2).GT.CKIN(32)**2) MINT(51)=1 
        ENDIF   
    
      ELSEIF(ILIM.EQ.1) THEN    
C...Calculate limits on tau 
C...0) due to definition    
        TAUMN0=0.   
        TAUMX0=1.   
C...1) due to limits on subsystem mass  
        TAUMN1=CKIN(1)**2/VINT(2)   
        TAUMX1=1.   
        IF(CKIN(2).GE.0.) TAUMX1=CKIN(2)**2/VINT(2) 
C...2) due to limits on pT-hat (and non-overlapping rapidity intervals) 
        TM3=SQRT(SQM3+PTHMIN**2)    
        TM4=SQRT(SQM4+PTHMIN**2)    
        YDCOSH=1.   
        IF(CKIN(9).GT.CKIN(12)) YDCOSH=COSH(CKIN(9)-CKIN(12))   
        TAUMN2=(TM3**2+2.*TM3*TM4*YDCOSH+TM4**2)/VINT(2)    
        TAUMX2=1.   
C...3) due to limits on pT-hat and cos(theta-hat)   
        CTH2MN=MIN(CKIN(27)**2,CKIN(28)**2) 
        CTH2MX=MAX(CKIN(27)**2,CKIN(28)**2) 
        TAUMN3=0.   
        IF(CKIN(27)*CKIN(28).GT.0.) TAUMN3= 
     &  (SQRT(SQM3+PTHMIN**2/(1.-CTH2MN))+  
     &  SQRT(SQM4+PTHMIN**2/(1.-CTH2MN)))**2/VINT(2)    
        TAUMX3=1.   
        IF(CKIN(4).GE.0..AND.CTH2MX.LT.1.) TAUMX3=  
     &  (SQRT(SQM3+CKIN(4)**2/(1.-CTH2MX))+ 
     &  SQRT(SQM4+CKIN(4)**2/(1.-CTH2MX)))**2/VINT(2)   
C...4) due to limits on x1 and x2   
        TAUMN4=CKIN(21)*CKIN(23)    
        TAUMX4=CKIN(22)*CKIN(24)    
C...5) due to limits on xF  
        TAUMN5=0.   
        TAUMX5=MAX(1.-CKIN(25),1.+CKIN(26)) 
        VINT(11)=MAX(TAUMN0,TAUMN1,TAUMN2,TAUMN3,TAUMN4,TAUMN5) 
        VINT(31)=MIN(TAUMX0,TAUMX1,TAUMX2,TAUMX3,TAUMX4,TAUMX5) 
        IF(MINT(43).EQ.1.AND.(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.2)) THEN 
          VINT(11)=0.99999  
          VINT(31)=1.00001  
        ENDIF   
        IF(VINT(31).LE.VINT(11)) MINT(51)=1 
    
      ELSEIF(ILIM.EQ.2) THEN    
C...Calculate limits on y*  
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) TAU=VINT(26) 
        TAURT=SQRT(TAU) 
C...0) due to kinematics    
        YSTMN0=LOG(TAURT)   
        YSTMX0=-YSTMN0  
C...1) due to explicit limits   
        YSTMN1=CKIN(7)  
        YSTMX1=CKIN(8)  
C...2) due to limits on x1  
        YSTMN2=LOG(MAX(TAU,CKIN(21))/TAURT) 
        YSTMX2=LOG(MAX(TAU,CKIN(22))/TAURT) 
C...3) due to limits on x2  
        YSTMN3=-LOG(MAX(TAU,CKIN(24))/TAURT)    
        YSTMX3=-LOG(MAX(TAU,CKIN(23))/TAURT)    
C...4) due to limits on xF  
        YEPMN4=0.5*ABS(CKIN(25))/TAURT  
        YSTMN4=SIGN(LOG(SQRT(1.+YEPMN4**2)+YEPMN4),CKIN(25))    
        YEPMX4=0.5*ABS(CKIN(26))/TAURT  
        YSTMX4=SIGN(LOG(SQRT(1.+YEPMX4**2)+YEPMX4),CKIN(26))    
C...5) due to simultaneous limits on y-large and y-small    
        YEPSMN=(RM3-RM4)*SINH(CKIN(9)-CKIN(11)) 
        YEPSMX=(RM3-RM4)*SINH(CKIN(10)-CKIN(12))    
        YDIFMN=ABS(LOG(SQRT(1.+YEPSMN**2)-YEPSMN))  
        YDIFMX=ABS(LOG(SQRT(1.+YEPSMX**2)-YEPSMX))  
        YSTMN5=0.5*(CKIN(9)+CKIN(11)-YDIFMN)    
        YSTMX5=0.5*(CKIN(10)+CKIN(12)+YDIFMX)   
C...6) due to simultaneous limits on cos(theta-hat) and y-large or  
C...   y-small  
        CTHLIM=SQRT(1.-4.*PTHMIN**2/(BE34*TAU*VINT(2))) 
        RZMN=BE34*MAX(CKIN(27),-CTHLIM) 
        RZMX=BE34*MIN(CKIN(28),CTHLIM)  
        YEX3MX=(1.+RM3-RM4+RZMX)/MAX(1E-10,1.+RM3-RM4-RZMX) 
        YEX4MX=(1.+RM4-RM3-RZMN)/MAX(1E-10,1.+RM4-RM3+RZMN) 
        YEX3MN=MAX(1E-10,1.+RM3-RM4+RZMN)/(1.+RM3-RM4-RZMN) 
        YEX4MN=MAX(1E-10,1.+RM4-RM3-RZMX)/(1.+RM4-RM3+RZMX) 
        YSTMN6=CKIN(9)-0.5*LOG(MAX(YEX3MX,YEX4MX))  
        YSTMX6=CKIN(12)-0.5*LOG(MIN(YEX3MN,YEX4MN)) 
        VINT(12)=MAX(YSTMN0,YSTMN1,YSTMN2,YSTMN3,YSTMN4,YSTMN5,YSTMN6)  
        VINT(32)=MIN(YSTMX0,YSTMX1,YSTMX2,YSTMX3,YSTMX4,YSTMX5,YSTMX6)  
        IF(MINT(43).EQ.1) THEN  
          VINT(12)=-0.00001 
          VINT(32)=0.00001  
        ELSEIF(MINT(43).EQ.2) THEN  
          VINT(12)=0.99999*YSTMX0   
          VINT(32)=1.00001*YSTMX0   
        ELSEIF(MINT(43).EQ.3) THEN  
          VINT(12)=-1.00001*YSTMX0  
          VINT(32)=-0.99999*YSTMX0  
        ENDIF   
        IF(VINT(32).LE.VINT(12)) MINT(51)=1 
    
      ELSEIF(ILIM.EQ.3) THEN    
C...Calculate limits on cos(theta-hat)  
        YST=VINT(22)    
C...0) due to definition    
        CTNMN0=-1.  
        CTNMX0=0.   
        CTPMN0=0.   
        CTPMX0=1.   
C...1) due to explicit limits   
        CTNMN1=MIN(0.,CKIN(27)) 
        CTNMX1=MIN(0.,CKIN(28)) 
        CTPMN1=MAX(0.,CKIN(27)) 
        CTPMX1=MAX(0.,CKIN(28)) 
C...2) due to limits on pT-hat  
        CTNMN2=-SQRT(1.-4.*PTHMIN**2/(BE34**2*TAU*VINT(2))) 
        CTPMX2=-CTNMN2  
        CTNMX2=0.   
        CTPMN2=0.   
        IF(CKIN(4).GE.0.) THEN  
          CTNMX2=-SQRT(MAX(0.,1.-4.*CKIN(4)**2/(BE34**2*TAU*VINT(2))))  
          CTPMN2=-CTNMX2    
        ENDIF   
C...3) due to limits on y-large and y-small 
        CTNMN3=MIN(0.,MAX((1.+RM3-RM4)/BE34*TANH(CKIN(11)-YST), 
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(10)-YST))) 
        CTNMX3=MIN(0.,(1.+RM3-RM4)/BE34*TANH(CKIN(12)-YST), 
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(9)-YST))   
        CTPMN3=MAX(0.,(1.+RM3-RM4)/BE34*TANH(CKIN(9)-YST),  
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(12)-YST))  
        CTPMX3=MAX(0.,MIN((1.+RM3-RM4)/BE34*TANH(CKIN(10)-YST), 
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(11)-YST))) 
        VINT(13)=MAX(CTNMN0,CTNMN1,CTNMN2,CTNMN3)   
        VINT(33)=MIN(CTNMX0,CTNMX1,CTNMX2,CTNMX3)   
        VINT(14)=MAX(CTPMN0,CTPMN1,CTPMN2,CTPMN3)   
        VINT(34)=MIN(CTPMX0,CTPMX1,CTPMX2,CTPMX3)   
        IF(VINT(33).LE.VINT(13).AND.VINT(34).LE.VINT(14)) MINT(51)=1    
    
      ELSEIF(ILIM.EQ.4) THEN    
C...Calculate limits on tau'    
C...0) due to kinematics    
        TAPMN0=TAU  
        TAPMX0=1.   
C...1) due to explicit limits   
        TAPMN1=CKIN(31)**2/VINT(2)  
        TAPMX1=1.   
        IF(CKIN(32).GE.0.) TAPMX1=CKIN(32)**2/VINT(2)   
        VINT(16)=MAX(TAPMN0,TAPMN1) 
        VINT(36)=MIN(TAPMX0,TAPMX1) 
        IF(MINT(43).EQ.1) THEN  
          VINT(16)=0.99999  
          VINT(36)=1.00001  
        ENDIF   
        IF(VINT(36).LE.VINT(16)) MINT(51)=1 
    
      ENDIF 
      RETURN    
    
C...Special case for low-pT and multiple interactions:  
C...effective kinematical limits for tau, y*, cos(theta-hat).   
  110 IF(ILIM.EQ.0) THEN    
      ELSEIF(ILIM.EQ.1) THEN    
        IF(MSTP(82).LE.1) VINT(11)=4.*PARP(81)**2/VINT(2)   
        IF(MSTP(82).GE.2) VINT(11)=PARP(82)**2/VINT(2)  
        VINT(31)=1. 
      ELSEIF(ILIM.EQ.2) THEN    
        VINT(12)=0.5*LOG(VINT(21))  
        VINT(32)=-VINT(12)  
      ELSEIF(ILIM.EQ.3) THEN    
        IF(MSTP(82).LE.1) ST2EFF=4.*PARP(81)**2/(VINT(21)*VINT(2))  
        IF(MSTP(82).GE.2) ST2EFF=0.01*PARP(82)**2/(VINT(21)*VINT(2))    
        VINT(13)=-SQRT(MAX(0.,1.-ST2EFF))   
        VINT(33)=0. 
        VINT(14)=0. 
        VINT(34)=-VINT(13)  
      ENDIF 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYKMAP(IVAR,MVAR,VVAR) 
    
C...Maps a uniform distribution into a distribution of a kinematical    
C...variable according to one of the possibilities allowed. It is   
C...assumed that kinematical limits have been set by a PYKLIM call. 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
    
C...Convert VVAR to tau variable.   
      ISUB=MINT(1)  
      IF(IVAR.EQ.1) THEN    
        TAUMIN=VINT(11) 
        TAUMAX=VINT(31) 
        IF(MVAR.EQ.3.OR.MVAR.EQ.4) THEN 
          TAURE=VINT(73)    
          GAMRE=VINT(74)    
        ELSEIF(MVAR.EQ.5.OR.MVAR.EQ.6) THEN 
          TAURE=VINT(75)    
          GAMRE=VINT(76)    
        ENDIF   
        IF(MINT(43).EQ.1.AND.(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.2)) THEN 
          TAU=1.    
        ELSEIF(MVAR.EQ.1) THEN  
          TAU=TAUMIN*(TAUMAX/TAUMIN)**VVAR  
        ELSEIF(MVAR.EQ.2) THEN  
          TAU=TAUMAX*TAUMIN/(TAUMIN+(TAUMAX-TAUMIN)*VVAR)   
        ELSEIF(MVAR.EQ.3.OR.MVAR.EQ.5) THEN 
          RATGEN=(TAURE+TAUMAX)/(TAURE+TAUMIN)*TAUMIN/TAUMAX    
          TAU=TAURE*TAUMIN/((TAURE+TAUMIN)*RATGEN**VVAR-TAUMIN) 
        ELSE    
          AUPP=ATAN((TAUMAX-TAURE)/GAMRE)   
          ALOW=ATAN((TAUMIN-TAURE)/GAMRE)   
          TAU=TAURE+GAMRE*TAN(ALOW+(AUPP-ALOW)*VVAR)    
        ENDIF   
        VINT(21)=MIN(TAUMAX,MAX(TAUMIN,TAU))    
    
C...Convert VVAR to y* variable.    
      ELSEIF(IVAR.EQ.2) THEN    
        YSTMIN=VINT(12) 
        YSTMAX=VINT(32) 
        IF(MINT(43).EQ.1) THEN  
          YST=0.    
        ELSEIF(MINT(43).EQ.2) THEN  
          IF(ISET(ISUB).LE.2) YST=-0.5*LOG(VINT(21))    
          IF(ISET(ISUB).GE.3) YST=-0.5*LOG(VINT(26))    
        ELSEIF(MINT(43).EQ.3) THEN  
          IF(ISET(ISUB).LE.2) YST=0.5*LOG(VINT(21)) 
          IF(ISET(ISUB).GE.3) YST=0.5*LOG(VINT(26)) 
        ELSEIF(MVAR.EQ.1) THEN  
          YST=YSTMIN+(YSTMAX-YSTMIN)*SQRT(VVAR) 
        ELSEIF(MVAR.EQ.2) THEN  
          YST=YSTMAX-(YSTMAX-YSTMIN)*SQRT(1.-VVAR)  
        ELSE    
          AUPP=ATAN(EXP(YSTMAX))    
          ALOW=ATAN(EXP(YSTMIN))    
          YST=LOG(TAN(ALOW+(AUPP-ALOW)*VVAR))   
        ENDIF   
        VINT(22)=MIN(YSTMAX,MAX(YSTMIN,YST))    
    
C...Convert VVAR to cos(theta-hat) variable.    
      ELSEIF(IVAR.EQ.3) THEN    
        RM34=2.*VINT(63)*VINT(64)/(VINT(21)*VINT(2))**2 
        RSQM=1.+RM34    
        IF(2.*VINT(71)**2/(VINT(21)*VINT(2)).LT.0.0001) RM34=MAX(RM34,  
     &  2.*VINT(71)**2/(VINT(21)*VINT(2)))  
        CTNMIN=VINT(13) 
        CTNMAX=VINT(33) 
        CTPMIN=VINT(14) 
        CTPMAX=VINT(34) 
        IF(MVAR.EQ.1) THEN  
          ANEG=CTNMAX-CTNMIN    
          APOS=CTPMAX-CTPMIN    
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=CTNMIN+(CTNMAX-CTNMIN)*VCTN 
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=CTPMIN+(CTPMAX-CTPMIN)*VCTP 
          ENDIF 
        ELSEIF(MVAR.EQ.2) THEN  
          RMNMIN=MAX(RM34,RSQM-CTNMIN)  
          RMNMAX=MAX(RM34,RSQM-CTNMAX)  
          RMPMIN=MAX(RM34,RSQM-CTPMIN)  
          RMPMAX=MAX(RM34,RSQM-CTPMAX)  
          ANEG=LOG(RMNMIN/RMNMAX)   
          APOS=LOG(RMPMIN/RMPMAX)   
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=RSQM-RMNMIN*(RMNMAX/RMNMIN)**VCTN   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=RSQM-RMPMIN*(RMPMAX/RMPMIN)**VCTP   
          ENDIF 
        ELSEIF(MVAR.EQ.3) THEN  
          RMNMIN=MAX(RM34,RSQM+CTNMIN)  
          RMNMAX=MAX(RM34,RSQM+CTNMAX)  
          RMPMIN=MAX(RM34,RSQM+CTPMIN)  
          RMPMAX=MAX(RM34,RSQM+CTPMAX)  
          ANEG=LOG(RMNMAX/RMNMIN)   
          APOS=LOG(RMPMAX/RMPMIN)   
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=RMNMIN*(RMNMAX/RMNMIN)**VCTN-RSQM   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=RMPMIN*(RMPMAX/RMPMIN)**VCTP-RSQM   
          ENDIF 
        ELSEIF(MVAR.EQ.4) THEN  
          RMNMIN=MAX(RM34,RSQM-CTNMIN)  
          RMNMAX=MAX(RM34,RSQM-CTNMAX)  
          RMPMIN=MAX(RM34,RSQM-CTPMIN)  
          RMPMAX=MAX(RM34,RSQM-CTPMAX)  
          ANEG=1./RMNMAX-1./RMNMIN  
          APOS=1./RMPMAX-1./RMPMIN  
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=RSQM-1./(1./RMNMIN+ANEG*VCTN)   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=RSQM-1./(1./RMPMIN+APOS*VCTP)   
          ENDIF 
        ELSEIF(MVAR.EQ.5) THEN  
          RMNMIN=MAX(RM34,RSQM+CTNMIN)  
          RMNMAX=MAX(RM34,RSQM+CTNMAX)  
          RMPMIN=MAX(RM34,RSQM+CTPMIN)  
          RMPMAX=MAX(RM34,RSQM+CTPMAX)  
          ANEG=1./RMNMIN-1./RMNMAX  
          APOS=1./RMPMIN-1./RMPMAX  
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=1./(1./RMNMIN-ANEG*VCTN)-RSQM   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=1./(1./RMPMIN-APOS*VCTP)-RSQM   
          ENDIF 
        ENDIF   
        IF(CTH.LT.0.) CTH=MIN(CTNMAX,MAX(CTNMIN,CTH))   
        IF(CTH.GT.0.) CTH=MIN(CTPMAX,MAX(CTPMIN,CTH))   
        VINT(23)=CTH    
    
C...Convert VVAR to tau' variable.  
      ELSEIF(IVAR.EQ.4) THEN    
        TAU=VINT(11)    
        TAUPMN=VINT(16) 
        TAUPMX=VINT(36) 
        IF(MINT(43).EQ.1) THEN  
          TAUP=1.   
        ELSEIF(MVAR.EQ.1) THEN  
          TAUP=TAUPMN*(TAUPMX/TAUPMN)**VVAR 
        ELSE    
          AUPP=(1.-TAU/TAUPMX)**4   
          ALOW=(1.-TAU/TAUPMN)**4   
          TAUP=TAU/(1.-(ALOW+(AUPP-ALOW)*VVAR)**0.25)   
        ENDIF   
        VINT(26)=MIN(TAUPMX,MAX(TAUPMN,TAUP))   
      ENDIF 
    
      RETURN    
      END   
    
C***********************************************************************    
    
      SUBROUTINE PYSIGH(NCHN,SIGS)  
    
C...Differential matrix elements for all included subprocesses. 
C...Note that what is coded is (disregarding the COMFAC factor) 
C...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,  
C...when d(sigma-hat) is given in the zero-width limit, the delta   
C...function in tau is replaced by a Breit-Wigner:  
C...1/pi*(s*m_res*Gamma_res)/((s*tau-m_res^2)^2+(m_res*Gamma_res)^2);   
C...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);   
C...i.e., dimensionless quantities. COMFAC contains the factor  
C...pi/s and the conversion factor from GeV^-2 to mb.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      DIMENSION X(2),XPQ(-6:6),KFAC(2,-40:40),WDTP(0:40),WDTE(0:40,0:5) 
    
C...Reset number of channels and cross-section. 
      NCHN=0    
      SIGS=0.   
    
C...Read kinematical variables and limits.  
      ISUB=MINT(1)  
      TAUMIN=VINT(11)   
      YSTMIN=VINT(12)   
      CTNMIN=VINT(13)   
      CTPMIN=VINT(14)   
      XT2MIN=VINT(15)   
      TAUPMN=VINT(16)   
      TAU=VINT(21)  
      YST=VINT(22)  
      CTH=VINT(23)  
      XT2=VINT(25)  
      TAUP=VINT(26) 
      TAUMAX=VINT(31)   
      YSTMAX=VINT(32)   
      CTNMAX=VINT(33)   
      CTPMAX=VINT(34)   
      XT2MAX=VINT(35)   
      TAUPMX=VINT(36)   
    
C...Derive kinematical quantities.  
      IF(ISET(ISUB).LE.2.OR.ISET(ISUB).EQ.5) THEN   
        X(1)=SQRT(TAU)*EXP(YST) 
        X(2)=SQRT(TAU)*EXP(-YST)    
      ELSE  
        X(1)=SQRT(TAUP)*EXP(YST)    
        X(2)=SQRT(TAUP)*EXP(-YST)   
      ENDIF 
      IF(MINT(43).EQ.4.AND.ISET(ISUB).GE.1.AND. 
     &(X(1).GT.0.999.OR.X(2).GT.0.999)) RETURN  
      SH=TAU*VINT(2)    
      SQM3=VINT(63) 
      SQM4=VINT(64) 
      RM3=SQM3/SH   
      RM4=SQM4/SH   
      BE34=SQRT((1.-RM3-RM4)**2-4.*RM3*RM4) 
      RPTS=4.*VINT(71)**2/SH    
      BE34L=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4-RPTS))   
      RM34=2.*RM3*RM4   
      RSQM=1.+RM34  
      RTHM=(4.*RM3*RM4+RPTS)/(1.-RM3-RM4+BE34L) 
      TH=-0.5*SH*MAX(RTHM,1.-RM3-RM4-BE34*CTH)  
      UH=-0.5*SH*MAX(RTHM,1.-RM3-RM4+BE34*CTH)  
      SQPTH=0.25*SH*BE34**2*(1.-CTH**2) 
      SH2=SH**2 
      TH2=TH**2 
      UH2=UH**2 
    
C...Choice of Q2 scale. 
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        Q2=SH   
      ELSEIF(MOD(ISET(ISUB),2).EQ.0.OR.ISET(ISUB).EQ.5) THEN    
        IF(MSTP(32).EQ.1) THEN  
          Q2=2.*SH*TH*UH/(SH**2+TH**2+UH**2)    
        ELSEIF(MSTP(32).EQ.2) THEN  
          Q2=SQPTH+0.5*(SQM3+SQM4)  
        ELSEIF(MSTP(32).EQ.3) THEN  
          Q2=MIN(-TH,-UH)   
        ELSEIF(MSTP(32).EQ.4) THEN  
          Q2=SH 
        ENDIF   
        IF(ISET(ISUB).EQ.5.AND.MSTP(82).GE.2) Q2=Q2+PARP(82)**2 
      ENDIF 
    
C...Store derived kinematical quantities.   
      VINT(41)=X(1) 
      VINT(42)=X(2) 
      VINT(44)=SH   
      VINT(43)=SQRT(SH) 
      VINT(45)=TH   
      VINT(46)=UH   
      VINT(48)=SQPTH    
      VINT(47)=SQRT(SQPTH)  
      VINT(50)=TAUP*VINT(2) 
      VINT(49)=SQRT(MAX(0.,VINT(50)))   
      VINT(52)=Q2   
      VINT(51)=SQRT(Q2) 
    
C...Calculate parton structure functions.   
      IF(ISET(ISUB).LE.0) GOTO 145  
      IF(MINT(43).GE.2) THEN    
        Q2SF=Q2 
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          Q2SF=PMAS(23,1)**2    
          IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77) Q2SF=PMAS(24,1)**2  
        ENDIF   
        DO 100 I=3-MINT(41),MINT(42)    
        XSF=X(I)    
        IF(ISET(ISUB).EQ.5) XSF=X(I)/VINT(142+I)    
        CALL PYSTFU(MINT(10+I),XSF,Q2SF,XPQ,I)    
        DO 100 KFL=-6,6 
  100   XSFX(I,KFL)=XPQ(KFL)
      ENDIF 
    
C...Calculate alpha_strong and K-factor.    
      IF(MSTP(33).NE.3) AS=ULALPS(Q2)   
      FACK=1.   
      FACA=1.   
      IF(MSTP(33).EQ.1) THEN    
        FACK=PARP(31)   
      ELSEIF(MSTP(33).EQ.2) THEN    
        FACK=PARP(31)   
        FACA=PARP(32)/PARP(31)  
      ELSEIF(MSTP(33).EQ.3) THEN    
        Q2AS=PARP(33)*Q2    
        IF(ISET(ISUB).EQ.5.AND.MSTP(82).GE.2) Q2AS=Q2AS+    
     &  PARU(112)*PARP(82)  
        AS=ULALPS(Q2AS) 
      ENDIF 
      RADC=1.+AS/PARU(1)    
    
C...Set flags for allowed reacting partons/leptons. 
      DO 130 I=1,2  
      DO 110 J=-40,40   
  110 KFAC(I,J)=0   
      IF(MINT(40+I).EQ.1) THEN  
        KFAC(I,MINT(10+I))=1    
      ELSE  
        DO 120 J=-40,40 
        KFAC(I,J)=KFIN(I,J) 
        IF(ABS(J).GT.MSTP(54).AND.J.NE.21) KFAC(I,J)=0  
        IF(ABS(J).LE.6) THEN    
          IF(XSFX(I,J).LT.1.E-10) KFAC(I,J)=0   
        ELSEIF(J.EQ.21) THEN    
          IF(XSFX(I,0).LT.1.E-10) KFAC(I,21)=0  
        ENDIF   
  120   CONTINUE    
      ENDIF 
  130 CONTINUE  
    
C...Lower and upper limit for flavour loops.    
      MIN1=0    
      MAX1=0    
      MIN2=0    
      MAX2=0    
      DO 140 J=-20,20   
      IF(KFAC(1,-J).EQ.1) MIN1=-J   
      IF(KFAC(1,J).EQ.1) MAX1=J 
      IF(KFAC(2,-J).EQ.1) MIN2=-J   
      IF(KFAC(2,J).EQ.1) MAX2=J 
  140 CONTINUE  
      MINA=MIN(MIN1,MIN2)   
      MAXA=MAX(MAX1,MAX2)   
    
C...Common conversion factors (including Jacobian) for subprocesses.    
      SQMZ=PMAS(23,1)**2    
      GMMZ=PMAS(23,1)*PMAS(23,2)    
      SQMW=PMAS(24,1)**2    
      GMMW=PMAS(24,1)*PMAS(24,2)    
      SQMH=PMAS(25,1)**2    
      GMMH=PMAS(25,1)*PMAS(25,2)    
      SQMZP=PMAS(32,1)**2   
      GMMZP=PMAS(32,1)*PMAS(32,2)   
      SQMHC=PMAS(37,1)**2   
      GMMHC=PMAS(37,1)*PMAS(37,2)   
      SQMR=PMAS(40,1)**2    
      GMMR=PMAS(40,1)*PMAS(40,2)    
      AEM=PARU(101) 
      XW=PARU(102)  
    
C...Phase space integral in tau and y*. 
      COMFAC=PARU(1)*PARU(5)/VINT(2)    
      IF(MINT(43).EQ.4) COMFAC=COMFAC*FACK  
      IF((MINT(43).GE.2.OR.ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4).AND. 
     &ISET(ISUB).NE.5) THEN 
        ATAU0=LOG(TAUMAX/TAUMIN)    
        ATAU1=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)   
        H1=COEF(ISUB,1)+(ATAU0/ATAU1)*COEF(ISUB,2)/TAU  
        IF(MINT(72).GE.1) THEN  
          TAUR1=VINT(73)    
          GAMR1=VINT(74)    
          ATAU2=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1  
          ATAU3=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/    
     &    GAMR1 
          H1=H1+(ATAU0/ATAU2)*COEF(ISUB,3)/(TAU+TAUR1)+ 
     &    (ATAU0/ATAU3)*COEF(ISUB,4)*TAU/((TAU-TAUR1)**2+GAMR1**2)  
        ENDIF   
        IF(MINT(72).EQ.2) THEN  
          TAUR2=VINT(75)    
          GAMR2=VINT(76)    
          ATAU4=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2  
          ATAU5=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/    
     &    GAMR2 
          H1=H1+(ATAU0/ATAU4)*COEF(ISUB,5)/(TAU+TAUR2)+ 
     &    (ATAU0/ATAU5)*COEF(ISUB,6)*TAU/((TAU-TAUR2)**2+GAMR2**2)  
        ENDIF   
        COMFAC=COMFAC*ATAU0/(TAU*H1)    
      ENDIF 
      IF(MINT(43).EQ.4.AND.ISET(ISUB).NE.5) THEN    
        AYST0=YSTMAX-YSTMIN 
        AYST1=0.5*(YSTMAX-YSTMIN)**2    
        AYST2=AYST1 
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))  
        H2=(AYST0/AYST1)*COEF(ISUB,7)*(YST-YSTMIN)+(AYST0/AYST2)*   
     &  COEF(ISUB,8)*(YSTMAX-YST)+(AYST0/AYST3)*COEF(ISUB,9)/COSH(YST)  
        COMFAC=COMFAC*AYST0/H2  
      ENDIF 
    
C...2 -> 1 processes: reduction in angular part of phase space integral 
C...for case of decaying resonance. 
      ACTH0=CTNMAX-CTNMIN+CTPMAX-CTPMIN 
      IF((ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3).AND.  
     &MDCY(KFPR(ISUB,1),1).EQ.1) THEN   
        IF(KFPR(ISUB,1).EQ.25.OR.KFPR(ISUB,1).EQ.37) THEN   
          COMFAC=COMFAC*0.5*ACTH0   
        ELSE    
          COMFAC=COMFAC*0.125*(3.*ACTH0+CTNMAX**3-CTNMIN**3+    
     &    CTPMAX**3-CTPMIN**3)  
        ENDIF   
    
C...2 -> 2 processes: angular part of phase space integral. 
      ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        ACTH1=LOG((MAX(RM34,RSQM-CTNMIN)*MAX(RM34,RSQM-CTPMIN))/    
     &  (MAX(RM34,RSQM-CTNMAX)*MAX(RM34,RSQM-CTPMAX)))  
        ACTH2=LOG((MAX(RM34,RSQM+CTNMAX)*MAX(RM34,RSQM+CTPMAX))/    
     &  (MAX(RM34,RSQM+CTNMIN)*MAX(RM34,RSQM+CTPMIN)))  
        ACTH3=1./MAX(RM34,RSQM-CTNMAX)-1./MAX(RM34,RSQM-CTNMIN)+    
     &  1./MAX(RM34,RSQM-CTPMAX)-1./MAX(RM34,RSQM-CTPMIN)   
        ACTH4=1./MAX(RM34,RSQM+CTNMIN)-1./MAX(RM34,RSQM+CTNMAX)+    
     &  1./MAX(RM34,RSQM+CTPMIN)-1./MAX(RM34,RSQM+CTPMAX)   
        H3=COEF(ISUB,10)+   
     &  (ACTH0/ACTH1)*COEF(ISUB,11)/MAX(RM34,RSQM-CTH)+ 
     &  (ACTH0/ACTH2)*COEF(ISUB,12)/MAX(RM34,RSQM+CTH)+ 
     &  (ACTH0/ACTH3)*COEF(ISUB,13)/MAX(RM34,RSQM-CTH)**2+  
     &  (ACTH0/ACTH4)*COEF(ISUB,14)/MAX(RM34,RSQM+CTH)**2   
        COMFAC=COMFAC*ACTH0*0.5*BE34/H3 
      ENDIF 
    
C...2 -> 3, 4 processes: phace space integral in tau'.  
      IF(MINT(43).GE.2.AND.(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4)) THEN   
        ATAUP0=LOG(TAUPMX/TAUPMN)   
        ATAUP1=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU) 
        H4=COEF(ISUB,15)+   
     &  ATAUP0/ATAUP1*COEF(ISUB,16)/TAUP*(1.-TAU/TAUP)**3   
        IF(1.-TAU/TAUP.GT.1.E-4) THEN   
          FZW=(1.+TAU/TAUP)*LOG(TAUP/TAU)-2.*(1.-TAU/TAUP)  
        ELSE    
          FZW=1./6.*(1.-TAU/TAUP)**3*TAU/TAUP   
        ENDIF   
        COMFAC=COMFAC*ATAUP0*FZW/H4 
      ENDIF 
    
C...Phase space integral for low-pT and multiple interactions.  
      IF(ISET(ISUB).EQ.5) THEN  
        COMFAC=PARU(1)*PARU(5)*FACK*0.5*VINT(2)/SH2 
        ATAU0=LOG(2.*(1.+SQRT(1.-XT2))/XT2-1.)  
        ATAU1=2.*ATAN(1./XT2-1.)/SQRT(XT2)  
        H1=COEF(ISUB,1)+(ATAU0/ATAU1)*COEF(ISUB,2)/SQRT(TAU)    
        COMFAC=COMFAC*ATAU0/H1  
        AYST0=YSTMAX-YSTMIN 
        AYST1=0.5*(YSTMAX-YSTMIN)**2    
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))  
        H2=(AYST0/AYST1)*COEF(ISUB,7)*(YST-YSTMIN)+(AYST0/AYST1)*   
     &  COEF(ISUB,8)*(YSTMAX-YST)+(AYST0/AYST3)*COEF(ISUB,9)/COSH(YST)  
        COMFAC=COMFAC*AYST0/H2  
        IF(MSTP(82).LE.1) COMFAC=COMFAC*XT2**2*(1./VINT(149)-1.)    
C...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is 
C...introduced to make cross-section finite for xT2 -> 0.   
        IF(MSTP(82).GE.2) COMFAC=COMFAC*XT2**2/(VINT(149)*  
     &  (1.+VINT(149))) 
      ENDIF 
    
C...A: 2 -> 1, tree diagrams.   
    
  145 IF(ISUB.LE.10) THEN   
      IF(ISUB.EQ.1) THEN    
C...f + fb -> gamma*/Z0.    
        MINT(61)=2  
        CALL PYWIDT(23,SQRT(SH),WDTP,WDTE)  
        FACZ=COMFAC*AEM**2*4./3.    
        DO 150 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 150    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        FACF=1. 
        IF(IABS(I).LE.10) FACF=FACA/3.  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACF*FACZ*(EI**2*VINT(111)+EI*VI/(8.*XW*(1.-XW))*    
     &  SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)*VINT(112)+(VI**2+AI**2)/    
     &  (16.*XW*(1.-XW))**2*SH2/((SH-SQMZ)**2+GMMZ**2)*VINT(114))   
  150   CONTINUE    
    
      ELSEIF(ISUB.EQ.2) THEN    
C...f + fb' -> W+/-.    
        CALL PYWIDT(24,SQRT(SH),WDTP,WDTE)  
        FACW=COMFAC*(AEM/XW)**2*1./24*SH2/((SH-SQMW)**2+GMMW**2)    
        DO 170 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 170   
        IA=IABS(I)  
        DO 160 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 160   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 160  
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 160 
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FACF=1. 
        IF(IA.LE.10) FACF=VCKM((IA+1)/2,(JA+1)/2)*FACA/3.   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACF*FACW*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))   
  160   CONTINUE    
  170   CONTINUE    
    
      ELSEIF(ISUB.EQ.3) THEN    
C...f + fb -> H0.   
        CALL PYWIDT(25,SQRT(SH),WDTP,WDTE)  
        FACH=COMFAC*(AEM/XW)**2*1./48.*(SH/SQMW)**2*    
     &  SH2/((SH-SQMH)**2+GMMH**2)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))  
        DO 180 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 180    
        RMQ=PMAS(IABS(I),1)**2/SH   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*RMQ*SQRT(MAX(0.,1.-4.*RMQ)) 
  180   CONTINUE    
    
      ELSEIF(ISUB.EQ.4) THEN    
C...gamma + W+/- -> W+/-.   
    
      ELSEIF(ISUB.EQ.5) THEN    
C...Z0 + Z0 -> H0.  
        CALL PYWIDT(25,SQRT(SH),WDTP,WDTE)  
        FACH=COMFAC*1./(128.*PARU(1)**2*16.*(1.-XW)**3)*(AEM/XW)**4*    
     &  (SH/SQMW)**2*SH2/((SH-SQMH)**2+GMMH**2)*    
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,4)) 
        DO 200 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 200   
        DO 190 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 190   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*(VI**2+AI**2)*(VJ**2+AJ**2) 
  190   CONTINUE    
  200   CONTINUE    
    
      ELSEIF(ISUB.EQ.6) THEN    
C...Z0 + W+/- -> W+/-.  
    
      ELSEIF(ISUB.EQ.7) THEN    
C...W+ + W- -> Z0.  
    
      ELSEIF(ISUB.EQ.8) THEN    
C...W+ + W- -> H0.  
        CALL PYWIDT(25,SQRT(SH),WDTP,WDTE)  
        FACH=COMFAC*1./(128*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2*   
     &  SH2/((SH-SQMH)**2+GMMH**2)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))  
        DO 220 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 220   
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)    
        DO 210 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 210   
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)    
        IF(EI*EJ.GT.0.) GOTO 210    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*VINT(180+I)*VINT(180+J) 
  210   CONTINUE    
  220   CONTINUE    
      ENDIF 
    
C...B: 2 -> 2, tree diagrams.   
    
      ELSEIF(ISUB.LE.20) THEN   
      IF(ISUB.EQ.11) THEN   
C...f + f' -> f + f'.   
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2 
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-  
     &  MSTP(34)*2./3.*UH2/(SH*TH)) 
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-   
     &  MSTP(34)*2./3.*SH2/(TH*UH)) 
        DO 240 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 240   
        DO 230 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 230   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQ1   
        IF(I.EQ.-J) SIGH(NCHN)=FACQQB   
        IF(I.EQ.J) THEN 
          SIGH(NCHN)=0.5*SIGH(NCHN) 
          NCHN=NCHN+1   
          ISIG(NCHN,1)=I    
          ISIG(NCHN,2)=J    
          ISIG(NCHN,3)=2    
          SIGH(NCHN)=0.5*FACQQ2 
        ENDIF   
  230   CONTINUE    
  240   CONTINUE    
    
      ELSEIF(ISUB.EQ.12) THEN   
C...f + fb -> f' + fb' (q + qb -> q' + qb' only).   
        CALL PYWIDT(21,SQRT(SH),WDTP,WDTE)  
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+   
     &  WDTE(0,3)+WDTE(0,4))    
        DO 250 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 250    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQB   
  250   CONTINUE    
    
      ELSEIF(ISUB.EQ.13) THEN   
C...f + fb -> g + g (q + qb -> g + g only). 
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2) 
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2) 
        DO 260 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 260    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=0.5*FACGG2   
  260   CONTINUE    
    
      ELSEIF(ISUB.EQ.14) THEN   
C...f + fb -> g + gamma (q + qb -> g + gamma only). 
        FACGG=COMFAC*AS*AEM*8./9.*(TH2+UH2)/(TH*UH) 
        DO 270 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 270    
        EI=KCHG(IABS(I),1)/3.   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGG*EI**2  
  270   CONTINUE    
    
      ELSEIF(ISUB.EQ.15) THEN   
C...f + fb -> g + Z0 (q + qb -> g + Z0 only).   
        FACZG=COMFAC*AS*AEM/(XW*(1.-XW))*1./18.*    
     &  (TH2+UH2+2.*SQM4*SH)/(TH*UH)    
        FACZG=FACZG*WIDS(23,2)  
        DO 280 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 280    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZG*(VI**2+AI**2)  
  280   CONTINUE    
    
      ELSEIF(ISUB.EQ.16) THEN   
C...f + fb' -> g + W+/- (q + qb' -> g + W+/- only). 
        FACWG=COMFAC*AS*AEM/XW*2./9.*(TH2+UH2+2.*SQM4*SH)/(TH*UH)   
        DO 300 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 300   
        IA=IABS(I)  
        DO 290 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 290   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 290  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACWG*FCKM*WIDS(24,(5-KCHW)/2)   
  290   CONTINUE    
  300   CONTINUE    
    
      ELSEIF(ISUB.EQ.17) THEN   
C...f + fb -> g + H0 (q + qb -> g + H0 only).   
    
      ELSEIF(ISUB.EQ.18) THEN   
C...f + fb -> gamma + gamma.    
        FACGG=COMFAC*FACA*AEM**2*1./3.*(TH2+UH2)/(TH*UH)    
        DO 310 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 310    
        EI=KCHG(IABS(I),1)/3.   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGG*EI**4  
  310   CONTINUE    
    
      ELSEIF(ISUB.EQ.19) THEN   
C...f + fb -> gamma + Z0.   
        FACGZ=COMFAC*FACA*AEM**2/(XW*(1.-XW))*1./24.*   
     &  (TH2+UH2+2.*SQM4*SH)/(TH*UH)    
        FACGZ=FACGZ*WIDS(23,2)  
        DO 320 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 320    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGZ*EI**2*(VI**2+AI**2)    
  320   CONTINUE    
    
      ELSEIF(ISUB.EQ.20) THEN   
C...f + fb' -> gamma + W+/-.    
        FACGW=COMFAC*FACA*AEM**2/XW*1./6.*  
     &  ((2.*UH-TH)/(3.*(SH-SQM4)))**2*(TH2+UH2+2.*SQM4*SH)/(TH*UH) 
        DO 340 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 340   
        IA=IABS(I)  
        DO 330 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 330   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 330  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGW*FCKM*WIDS(24,(5-KCHW)/2)   
  330   CONTINUE    
  340   CONTINUE    
      ENDIF 
    
      ELSEIF(ISUB.LE.30) THEN   
      IF(ISUB.EQ.21) THEN   
C...f + fb -> gamma + H0.   
    
      ELSEIF(ISUB.EQ.22) THEN   
C...f + fb -> Z0 + Z0.  
        FACZZ=COMFAC*FACA*(AEM/(XW*(1.-XW)))**2*1./768.*    
     &  (UH/TH+TH/UH+2.*(SQM3+SQM4)*SH/(TH*UH)- 
     &  SQM3*SQM4*(1./TH2+1./UH2))  
        FACZZ=FACZZ*WIDS(23,1)  
        DO 350 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 350    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZZ*(VI**4+6.*VI**2*AI**2+AI**4)   
  350   CONTINUE    
    
      ELSEIF(ISUB.EQ.23) THEN   
C...f + fb' -> Z0 + W+/-.   
        FACZW=COMFAC*FACA*(AEM/XW)**2*1./6. 
        FACZW=FACZW*WIDS(23,2)  
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        DO 370 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 370   
        IA=IABS(I)  
        DO 360 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 360   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 360  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        EI=KCHG(IA,1)/3.    
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        EJ=KCHG(JA,1)/3.    
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        IF(VI+AI.GT.0) THEN 
          VISAV=VI  
          AISAV=AI  
          VI=VJ 
          AI=AJ 
          VJ=VISAV  
          AJ=AISAV  
        ENDIF   
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZW*FCKM*(1./(SH-SQMW)**2* 
     &  ((9.-8.*XW)/4.*THUH+(8.*XW-6.)/4.*SH*(SQM3+SQM4))+  
     &  (THUH-SH*(SQM3+SQM4))/(2.*(SH-SQMW))*((VJ+AJ)/TH-(VI+AI)/UH)+   
     &  THUH/(16.*(1.-XW))*((VJ+AJ)**2/TH2+(VI+AI)**2/UH2)+ 
     &  SH*(SQM3+SQM4)/(8.*(1.-XW))*(VI+AI)*(VJ+AJ)/(TH*UH))*   
     &  WIDS(24,(5-KCHW)/2) 
  360   CONTINUE    
  370   CONTINUE    
    
      ELSEIF(ISUB.EQ.24) THEN   
C...f + fb -> Z0 + H0.  
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        FACHZ=COMFAC*FACA*(AEM/(XW*(1.-XW)))**2*1./96.* 
     &  (THUH+2.*SH*SQMZ)/(SH-SQMZ)**2  
        FACHZ=FACHZ*WIDS(23,2)*WIDS(25,2)   
        DO 380 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 380    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHZ*(VI**2+AI**2)  
  380   CONTINUE    
    
      ELSEIF(ISUB.EQ.25) THEN   
C...f + fb -> W+ + W-.  
        FACWW=COMFAC*FACA*(AEM/XW)**2*1./12.    
        FACWW=FACWW*WIDS(24,1)  
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        DO 390 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 390    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        DSIGWW=THUH/SH2*(3.-(SH-3.*(SQM3+SQM4))/(SH-SQMZ)*  
     &  (VI+AI)/(2.*AI*(1.-XW))+(SH/(SH-SQMZ))**2*  
     &  (1.-2.*(SQM3+SQM4)/SH+12.*SQM3*SQM4/SH2)*(VI**2+AI**2)/ 
     &  (8.*(1.-XW)**2))-2.*SQMZ/(SH-SQMZ)*(VI+AI)/AI+  
     &  SQMZ*SH/(SH-SQMZ)**2*(1.-2.*(SQM3+SQM4)/SH)*(VI**2+AI**2)/  
     &  (2.*(1.-XW))    
        IF(KCHG(IABS(I),1).LT.0) THEN   
          DSIGWW=DSIGWW+2.*(1.+SQMZ/(SH-SQMZ)*(VI+AI)/(2.*AI))* 
     &    (THUH/(SH*TH)-(SQM3+SQM4)/TH)+THUH/TH2    
        ELSE    
          DSIGWW=DSIGWW+2.*(1.+SQMZ/(SH-SQMZ)*(VI+AI)/(2.*AI))* 
     &    (THUH/(SH*UH)-(SQM3+SQM4)/UH)+THUH/UH2    
        ENDIF   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACWW*DSIGWW 
  390   CONTINUE    
    
      ELSEIF(ISUB.EQ.26) THEN   
C...f + fb' -> W+/- + H0.   
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        FACHW=COMFAC*FACA*(AEM/XW)**2*1./24.*(THUH+2.*SH*SQMW)/ 
     &  (SH-SQMW)**2    
        FACHW=FACHW*WIDS(25,2)  
        DO 410 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 410   
        IA=IABS(I)  
        DO 400 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(1,J).EQ.0) GOTO 400   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 400  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHW*FCKM*WIDS(24,(5-KCHW)/2)   
  400   CONTINUE    
  410   CONTINUE    
    
      ELSEIF(ISUB.EQ.27) THEN   
C...f + fb -> H0 + H0.  
    
      ELSEIF(ISUB.EQ.28) THEN   
C...f + g -> f + g (q + g -> q + g only).   
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*  
     &  FACA    
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)   
        DO 430 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 430 
        DO 420 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 420    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 420    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQG1   
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=FACQG2   
  420   CONTINUE    
  430   CONTINUE    
    
      ELSEIF(ISUB.EQ.29) THEN   
C...f + g -> f + gamma (q + g -> q + gamma only).   
        FGQ=COMFAC*FACA*AS*AEM*1./3.*(SH2+UH2)/(-SH*UH) 
        DO 450 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 450 
        EI=KCHG(IABS(I),1)/3.   
        FACGQ=FGQ*EI**2 
        DO 440 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 440    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 440    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGQ    
  440   CONTINUE    
  450   CONTINUE    
    
      ELSEIF(ISUB.EQ.30) THEN   
C...f + g -> f + Z0 (q + g -> q + Z0 only). 
        FZQ=COMFAC*FACA*AS*AEM/(XW*(1.-XW))*1./48.* 
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)   
        FZQ=FZQ*WIDS(23,2)  
        DO 470 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 470 
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        FACZQ=FZQ*(VI**2+AI**2) 
        DO 460 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 460    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 460    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZQ    
  460   CONTINUE    
  470   CONTINUE    
      ENDIF 
    
      ELSEIF(ISUB.LE.40) THEN   
      IF(ISUB.EQ.31) THEN   
C...f + g -> f' + W+/- (q + g -> q' + W+/- only).   
        FACWQ=COMFAC*FACA*AS*AEM/XW*1./12.* 
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)   
        DO 490 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 490 
        IA=IABS(I)  
        KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I)) 
        DO 480 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 480    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 480    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACWQ*VINT(180+I)*WIDS(24,(5-KCHW)/2)    
  480   CONTINUE    
  490   CONTINUE    
    
      ELSEIF(ISUB.EQ.32) THEN   
C...f + g -> f + H0 (q + g -> q + H0 only). 
    
      ELSEIF(ISUB.EQ.33) THEN   
C...f + gamma -> f + g (q + gamma -> q + g only).   
    
      ELSEIF(ISUB.EQ.34) THEN   
C...f + gamma -> f + gamma. 
    
      ELSEIF(ISUB.EQ.35) THEN   
C...f + gamma -> f + Z0.    
    
      ELSEIF(ISUB.EQ.36) THEN   
C...f + gamma -> f' + W+/-. 
    
      ELSEIF(ISUB.EQ.37) THEN   
C...f + gamma -> f + H0.    
    
      ELSEIF(ISUB.EQ.38) THEN   
C...f + Z0 -> f + g (q + Z0 -> q + g only). 
    
      ELSEIF(ISUB.EQ.39) THEN   
C...f + Z0 -> f + gamma.    
    
      ELSEIF(ISUB.EQ.40) THEN   
C...f + Z0 -> f + Z0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.50) THEN   
      IF(ISUB.EQ.41) THEN   
C...f + Z0 -> f' + W+/-.    
    
      ELSEIF(ISUB.EQ.42) THEN   
C...f + Z0 -> f + H0.   
    
      ELSEIF(ISUB.EQ.43) THEN   
C...f + W+/- -> f' + g (q + W+/- -> q' + g only).   
    
      ELSEIF(ISUB.EQ.44) THEN   
C...f + W+/- -> f' + gamma. 
    
      ELSEIF(ISUB.EQ.45) THEN   
C...f + W+/- -> f' + Z0.    
    
      ELSEIF(ISUB.EQ.46) THEN   
C...f + W+/- -> f' + W+/-.  
    
      ELSEIF(ISUB.EQ.47) THEN   
C...f + W+/- -> f' + H0.    
    
      ELSEIF(ISUB.EQ.48) THEN   
C...f + H0 -> f + g (q + H0 -> q + g only). 
    
      ELSEIF(ISUB.EQ.49) THEN   
C...f + H0 -> f + gamma.    
    
      ELSEIF(ISUB.EQ.50) THEN   
C...f + H0 -> f + Z0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.60) THEN   
      IF(ISUB.EQ.51) THEN   
C...f + H0 -> f' + W+/-.    
    
      ELSEIF(ISUB.EQ.52) THEN   
C...f + H0 -> f + H0.   
    
      ELSEIF(ISUB.EQ.53) THEN   
C...g + g -> f + fb (g + g -> q + qb only). 
        CALL PYWIDT(21,SQRT(SH),WDTP,WDTE)  
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 500 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQ1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=FACQQ2   
  500   CONTINUE    
    
      ELSEIF(ISUB.EQ.54) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.55) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.56) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.57) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.58) THEN   
C...gamma + gamma -> f + fb.    
    
      ELSEIF(ISUB.EQ.59) THEN   
C...gamma + Z0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.60) THEN   
C...gamma + W+/- -> f + fb'.    
      ENDIF 
    
      ELSEIF(ISUB.LE.70) THEN   
      IF(ISUB.EQ.61) THEN   
C...gamma + H0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.62) THEN   
C...Z0 + Z0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.63) THEN   
C...Z0 + W+/- -> f + fb'.   
    
      ELSEIF(ISUB.EQ.64) THEN   
C...Z0 + H0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.65) THEN   
C...W+ + W- -> f + fb.  
    
      ELSEIF(ISUB.EQ.66) THEN   
C...W+/- + H0 -> f + fb'.   
    
      ELSEIF(ISUB.EQ.67) THEN   
C...H0 + H0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.68) THEN   
C...g + g -> g + g. 
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+    
     &  TH2/SH2)*FACA   
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+    
     &  SH2/UH2)*FACA   
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2) 
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 510 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=0.5*FACGG2   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=3  
        SIGH(NCHN)=0.5*FACGG3   
  510   CONTINUE    
    
      ELSEIF(ISUB.EQ.69) THEN   
C...gamma + gamma -> W+ + W-.   
    
      ELSEIF(ISUB.EQ.70) THEN   
C...gamma + W+/- -> gamma + W+/-.   
      ENDIF 
    
      ELSEIF(ISUB.LE.80) THEN   
      IF(ISUB.EQ.71) THEN   
C...Z0 + Z0 -> Z0 + Z0. 
        BE2=1.-4.*SQMZ/SH   
        TH=-0.5*SH*BE2*(1.-CTH) 
        UH=-0.5*SH*BE2*(1.+CTH) 
        SHANG=1./(1.-XW)*SQMW/SQMZ*(1.+BE2)**2  
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        THANG=1./(1.-XW)*SQMW/SQMZ*(BE2-CTH)**2 
        ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG    
        ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG    
        UHANG=1./(1.-XW)*SQMW/SQMZ*(BE2+CTH)**2 
        AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG    
        AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG    
        FACH=0.5*COMFAC*1./(4096.*PARU(1)**2*16.*(1.-XW)**2)*   
     &  (AEM/XW)**4*(SH/SQMW)**2*((ASHRE+ATHRE+AUHRE)**2+   
     &  (ASHIM+ATHIM+AUHIM)**2)*SQMZ/SQMW   
        DO 530 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 530   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        AVI=AI**2+VI**2 
        DO 520 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 520   
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        AVJ=AJ**2+VJ**2 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*AVI*AVJ 
  520   CONTINUE    
  530   CONTINUE    
    
      ELSEIF(ISUB.EQ.72) THEN   
C...Z0 + Z0 -> W+ + W-. 
        BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))   
        CTH2=CTH**2 
        TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)   
        UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)   
        SHANG=4.*SQRT(SQMW/(SQMZ*(1.-XW)))*(1.-2.*SQMW/SH)* 
     &  (1.-2.*SQMZ/SH) 
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        ATWRE=(1.-XW)/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        ATWIM=0.    
        AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        AUWIM=0.    
        A4RE=2.*(1.-XW)/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)    
        A4IM=0. 
        FACH=COMFAC*1./(4096.*PARU(1)**2*16.*(1.-XW)**2)*(AEM/XW)**4*   
     &  (SH/SQMW)**2*((ASHRE+ATWRE+AUWRE+A4RE)**2+  
     &  (ASHIM+ATWIM+AUWIM+A4IM)**2)*SQMZ/SQMW  
        DO 550 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 550   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        AVI=AI**2+VI**2 
        DO 540 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 540   
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        AVJ=AJ**2+VJ**2 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*AVI*AVJ 
  540   CONTINUE    
  550   CONTINUE    
    
      ELSEIF(ISUB.EQ.73) THEN   
C...Z0 + W+/- -> Z0 + W+/-. 
        BE2=1.-2.*(SQMZ+SQMW)/SH+((SQMZ-SQMW)/SH)**2    
        EP1=1.+(SQMZ-SQMW)/SH   
        EP2=1.-(SQMZ-SQMW)/SH   
        TH=-0.5*SH*BE2*(1.-CTH) 
        UH=(SQMZ-SQMW)**2/SH-0.5*SH*BE2*(1.+CTH)    
        THANG=SQRT(SQMW/(SQMZ*(1.-XW)))*(BE2-EP1*CTH)*(BE2-EP2*CTH) 
        ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG    
        ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG    
        ASWRE=(1.-XW)/SQMZ*SH/(SH-SQMW)*(-BE2*(EP1+EP2)**4*CTH+ 
     &  1./4.*(BE2+EP1*EP2)**2*((EP1-EP2)**2-4.*BE2*CTH)+   
     &  2.*BE2*(BE2+EP1*EP2)*(EP1+EP2)**2*CTH-  
     &  1./16.*SH/SQMW*(EP1**2-EP2**2)**2*(BE2+EP1*EP2)**2) 
        ASWIM=0.    
        AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*(-BE2*(EP2+EP1*CTH)*    
     &  (EP1+EP2*CTH)*(BE2+EP1*EP2)+BE2*(EP2+EP1*CTH)*  
     &  (BE2+EP1*EP2*CTH)*(2.*EP2-EP2*CTH+EP1)-BE2*(EP2+EP1*CTH)**2*    
     &  (BE2-EP2**2*CTH)-1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+  
     &  2.*BE2*(1.-CTH))+1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*   
     &  (EP1**2-EP2**2)**2-BE2*(EP1+EP2*CTH)*(EP2+EP1*CTH)* 
     &  (BE2+EP1*EP2)+BE2*(EP1+EP2*CTH)*(BE2+EP1*EP2*CTH)*  
     &  (2.*EP1-EP1*CTH+EP2)-BE2*(EP1+EP2*CTH)**2*(BE2-EP1**2*CTH)- 
     &  1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+2.*BE2*(1.-CTH))+  
     &  1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*(EP1**2-EP2**2)**2) 
        AUWIM=0.    
        A4RE=(1.-XW)/SQMZ*(EP1**2*EP2**2*(CTH**2-1.)-   
     &  2.*BE2*(EP1**2+EP2**2+EP1*EP2)*CTH-2.*BE2*EP1*EP2)  
        A4IM=0. 
        FACH=COMFAC*1./(4096.*PARU(1)**2*4.*(1.-XW))*(AEM/XW)**4*   
     &  (SH/SQMW)**2*((ATHRE+ASWRE+AUWRE+A4RE)**2+  
     &  (ATHIM+ASWIM+AUWIM+A4IM)**2)*SQRT(SQMZ/SQMW)    
        DO 570 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 570   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        AVI=AI**2+VI**2 
        DO 560 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 560   
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AI-4.*EJ*XW  
        AVJ=AJ**2+VJ**2 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*(AVI*VINT(180+J)+VINT(180+I)*AVJ)   
  560   CONTINUE    
  570   CONTINUE    
    
      ELSEIF(ISUB.EQ.75) THEN   
C...W+ + W- -> gamma + gamma.   
    
      ELSEIF(ISUB.EQ.76) THEN   
C...W+ + W- -> Z0 + Z0. 
        BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))   
        CTH2=CTH**2 
        TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)   
        UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)   
        SHANG=4.*SQRT(SQMW/(SQMZ*(1.-XW)))*(1.-2.*SQMW/SH)* 
     &  (1.-2.*SQMZ/SH) 
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        ATWRE=(1.-XW)/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        ATWIM=0.    
        AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        AUWIM=0.    
        A4RE=2.*(1.-XW)/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)    
        A4IM=0. 
        FACH=0.5*COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2* 
     &  ((ASHRE+ATWRE+AUWRE+A4RE)**2+(ASHIM+ATWIM+AUWIM+A4IM)**2)   
        DO 590 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 590   
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)    
        DO 580 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 580   
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)    
        IF(EI*EJ.GT.0.) GOTO 580    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*VINT(180+I)*VINT(180+J) 
  580   CONTINUE    
  590   CONTINUE    
    
      ELSEIF(ISUB.EQ.77) THEN   
C...W+/- + W+/- -> W+/- + W+/-. 
        BE2=1.-4.*SQMW/SH   
        BE4=BE2**2  
        CTH2=CTH**2 
        CTH3=CTH**3 
        TH=-0.5*SH*BE2*(1.-CTH) 
        UH=-0.5*SH*BE2*(1.+CTH) 
        SHANG=(1.+BE2)**2   
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        THANG=(BE2-CTH)**2  
        ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG    
        ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG    
        SGZANG=1./SQMW*BE2*(3.-BE2)**2*CTH  
        ASGRE=XW*SGZANG 
        ASGIM=0.    
        ASZRE=(1.-XW)*SH/(SH-SQMZ)*SGZANG   
        ASZIM=0.    
        TGZANG=1./SQMW*(BE2*(4.-2.*BE2+BE4)+BE2*(4.-10.*BE2+BE4)*CTH+   
     &  (2.-11.*BE2+10.*BE4)*CTH2+BE2*CTH3) 
        ATGRE=0.5*XW*SH/TH*TGZANG   
        ATGIM=0.    
        ATZRE=0.5*(1.-XW)*SH/(TH-SQMZ)*TGZANG   
        ATZIM=0.    
        A4RE=1./SQMW*(1.+2.*BE2-6.*BE2*CTH-CTH2)    
        A4IM=0. 
        FACH=COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2* 
     &  ((ASHRE+ATHRE+ASGRE+ASZRE+ATGRE+ATZRE+A4RE)**2+ 
     &  (ASHIM+ATHIM+ASGIM+ASZIM+ATGIM+ATZIM+A4IM)**2)  
        DO 610 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 610   
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)    
        DO 600 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 600   
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)    
        IF(EI*EJ.GT.0.) GOTO 600    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*VINT(180+I)*VINT(180+J) 
  600   CONTINUE    
  610   CONTINUE    
    
      ELSEIF(ISUB.EQ.78) THEN   
C...W+/- + H0 -> W+/- + H0. 
    
      ELSEIF(ISUB.EQ.79) THEN   
C...H0 + H0 -> H0 + H0. 
    
      ENDIF 
    
C...C: 2 -> 2, tree diagrams with masses.   
    
      ELSEIF(ISUB.LE.90) THEN   
      IF(ISUB.EQ.81) THEN   
C...q + qb -> Q + QB.   
        FACQQB=COMFAC*AS**2*4./9.*(((TH-SQM3)**2+   
     &  (UH-SQM3)**2)/SH2+2.*SQM3/SH)   
        IF(MSTP(35).GE.1) THEN  
          IF(MSTP(35).EQ.1) THEN    
            ALSSG=PARP(35)  
          ELSE  
            MST115=MSTU(115)    
            MSTU(115)=MSTP(36)  
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))   
            ALSSG=ULALPS(Q2BN)  
            MSTU(115)=MST115    
          ENDIF 
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))   
          FREPU=XREPU/(EXP(MIN(100.,XREPU))-1.) 
          PARI(81)=FREPU    
          FACQQB=FACQQB*FREPU   
        ENDIF   
        DO 620 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 620    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQB   
  620   CONTINUE    
    
      ELSEIF(ISUB.EQ.82) THEN   
C...g + g -> Q + QB.    
        FACQQ1=COMFAC*FACA*AS**2*1./6.*((UH-SQM3)/(TH-SQM3)-    
     &  2.*(UH-SQM3)**2/SH2+4.*SQM3/SH*(TH*UH-SQM3**2)/(TH-SQM3)**2)    
        FACQQ2=COMFAC*FACA*AS**2*1./6.*((TH-SQM3)/(UH-SQM3)-    
     &  2.*(TH-SQM3)**2/SH2+4.*SQM3/SH*(TH*UH-SQM3**2)/(UH-SQM3)**2)    
        IF(MSTP(35).GE.1) THEN  
          IF(MSTP(35).EQ.1) THEN    
            ALSSG=PARP(35)  
          ELSE  
            MST115=MSTU(115)    
            MSTU(115)=MSTP(36)  
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))   
            ALSSG=ULALPS(Q2BN)  
            MSTU(115)=MST115    
          ENDIF 
          XATTR=4.*PARU(1)*ALSSG/(3.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))    
          FATTR=XATTR/(1.-EXP(-MIN(100.,XATTR)))    
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))   
          FREPU=XREPU/(EXP(MIN(100.,XREPU))-1.) 
          FATRE=(2.*FATTR+5.*FREPU)/7.  
          PARI(81)=FATRE    
          FACQQ1=FACQQ1*FATRE   
          FACQQ2=FACQQ2*FATRE   
        ENDIF   
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 630 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQ1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=FACQQ2   
  630   CONTINUE    
    
      ENDIF 
    
C...D: Mimimum bias processes.  
    
      ELSEIF(ISUB.LE.100) THEN  
      IF(ISUB.EQ.91) THEN   
C...Elastic scattering. 
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.92) THEN   
C...Single diffractive scattering.  
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.93) THEN   
C...Double diffractive scattering.  
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.94) THEN   
C...Central diffractive scattering. 
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.95) THEN   
C...Low-pT scattering.  
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.96) THEN   
C...Multiple interactions: sum of QCD processes.    
        CALL PYWIDT(21,SQRT(SH),WDTP,WDTE)  
    
C...q + q' -> q + q'.   
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2 
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-  
     &  MSTP(34)*2./3.*UH2/(SH*TH)) 
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-   
     &  MSTP(34)*2./3.*SH2/(TH*UH)) 
        DO 650 I=-3,3   
        IF(I.EQ.0) GOTO 650 
        DO 640 J=-3,3   
        IF(J.EQ.0) GOTO 640 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=111    
        SIGH(NCHN)=FACQQ1   
        IF(I.EQ.-J) SIGH(NCHN)=FACQQB   
        IF(I.EQ.J) THEN 
          SIGH(NCHN)=0.5*SIGH(NCHN) 
          NCHN=NCHN+1   
          ISIG(NCHN,1)=I    
          ISIG(NCHN,2)=J    
          ISIG(NCHN,3)=112  
          SIGH(NCHN)=0.5*FACQQ2 
        ENDIF   
  640   CONTINUE    
  650   CONTINUE    
    
C...q + qb -> q' + qb' or g + g.    
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+   
     &  WDTE(0,3)+WDTE(0,4))    
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2) 
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2) 
        DO 660 I=-3,3   
        IF(I.EQ.0) GOTO 660 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=121    
        SIGH(NCHN)=FACQQB   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=131    
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=132    
        SIGH(NCHN)=0.5*FACGG2   
  660   CONTINUE    
    
C...q + g -> q + g. 
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*  
     &  FACA    
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)   
        DO 680 I=-3,3   
        IF(I.EQ.0) GOTO 680 
        DO 670 ISDE=1,2 
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=281    
        SIGH(NCHN)=FACQG1   
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=282    
        SIGH(NCHN)=FACQG2   
  670   CONTINUE    
  680   CONTINUE    
    
C...g + g -> q + qb or g + g.   
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+    
     &  TH2/SH2)*FACA   
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+    
     &  SH2/UH2)*FACA   
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2) 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=531    
        SIGH(NCHN)=FACQQ1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=532    
        SIGH(NCHN)=FACQQ2   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=681    
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=682    
        SIGH(NCHN)=0.5*FACGG2   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=683    
        SIGH(NCHN)=0.5*FACGG3   
      ENDIF 
    
C...E: 2 -> 1, loop diagrams.   
    
      ELSEIF(ISUB.LE.110) THEN  
      IF(ISUB.EQ.101) THEN  
C...g + g -> gamma*/Z0. 
    
      ELSEIF(ISUB.EQ.102) THEN  
C...g + g -> H0.    
        CALL PYWIDT(25,SQRT(SH),WDTP,WDTE)  
        ETARE=0.    
        ETAIM=0.    
        DO 690 I=1,2*MSTP(1)    
        EPS=4.*PMAS(I,1)**2/SH  
        IF(EPS.LE.1.) THEN  
          IF(EPS.GT.1.E-4) THEN 
            ROOT=SQRT(1.-EPS)   
            RLN=LOG((1.+ROOT)/(1.-ROOT))    
          ELSE  
            RLN=LOG(4./EPS-2.)  
          ENDIF 
          PHIRE=0.25*(RLN**2-PARU(1)**2)    
          PHIIM=0.5*PARU(1)*RLN 
        ELSE    
          PHIRE=-(ASIN(1./SQRT(EPS)))**2    
          PHIIM=0.  
        ENDIF   
        ETARE=ETARE+0.5*EPS*(1.+(EPS-1.)*PHIRE) 
        ETAIM=ETAIM+0.5*EPS*(EPS-1.)*PHIIM  
  690   CONTINUE    
        ETA2=ETARE**2+ETAIM**2  
        FACH=COMFAC*FACA*(AS/PARU(1)*AEM/XW)**2*1./512.*    
     &  (SH/SQMW)**2*ETA2*SH2/((SH-SQMH)**2+GMMH**2)*   
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,4)) 
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 700 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH 
  700   CONTINUE    
    
      ENDIF 
    
C...F: 2 -> 2, box diagrams.    
    
      ELSEIF(ISUB.LE.120) THEN  
      IF(ISUB.EQ.111) THEN  
C...f + fb -> g + H0 (q + qb -> g + H0 only).   
        A5STUR=0.   
        A5STUI=0.   
        DO 710 I=1,2*MSTP(1)    
        SQMQ=PMAS(I,1)**2   
        EPSS=4.*SQMQ/SH 
        EPSH=4.*SQMQ/SQMH   
        A5STUR=A5STUR+SQMQ/SQMH*(4.+4.*SH/(TH+UH)*(PYW1AU(EPSS,1)-  
     &  PYW1AU(EPSH,1))+(1.-4.*SQMQ/(TH+UH))*(PYW2AU(EPSS,1)-   
     &  PYW2AU(EPSH,1)))    
        A5STUI=A5STUI+SQMQ/SQMH*(4.*SH/(TH+UH)*(PYW1AU(EPSS,2)- 
     &  PYW1AU(EPSH,2))+(1.-4.*SQMQ/(TH+UH))*(PYW2AU(EPSS,2)-   
     &  PYW2AU(EPSH,2)))    
  710   CONTINUE    
        FACGH=COMFAC*FACA/(144.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW* 
     &  SQMH/SH*(UH**2+TH**2)/(UH+TH)**2*(A5STUR**2+A5STUI**2)  
        FACGH=FACGH*WIDS(25,2)  
        DO 720 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 720    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGH    
  720   CONTINUE    
    
      ELSEIF(ISUB.EQ.112) THEN  
C...f + g -> f + H0 (q + g -> q + H0 only). 
        A5TSUR=0.   
        A5TSUI=0.   
        DO 730 I=1,2*MSTP(1)    
        SQMQ=PMAS(I,1)**2   
        EPST=4.*SQMQ/TH 
        EPSH=4.*SQMQ/SQMH   
        A5TSUR=A5TSUR+SQMQ/SQMH*(4.+4.*TH/(SH+UH)*(PYW1AU(EPST,1)-  
     &  PYW1AU(EPSH,1))+(1.-4.*SQMQ/(SH+UH))*(PYW2AU(EPST,1)-   
     &  PYW2AU(EPSH,1)))    
        A5TSUI=A5TSUI+SQMQ/SQMH*(4.*TH/(SH+UH)*(PYW1AU(EPST,2)- 
     &  PYW1AU(EPSH,2))+(1.-4.*SQMQ/(SH+UH))*(PYW2AU(EPST,2)-   
     &  PYW2AU(EPSH,2)))    
  730   CONTINUE    
        FACQH=COMFAC*FACA/(384.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW* 
     &  SQMH/(-TH)*(UH**2+SH**2)/(UH+SH)**2*(A5TSUR**2+A5TSUI**2)   
        FACQH=FACQH*WIDS(25,2)  
        DO 750 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 750 
        DO 740 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 740    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 740    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQH    
  740   CONTINUE    
  750   CONTINUE    
    
      ELSEIF(ISUB.EQ.113) THEN  
C...g + g -> g + H0.    
        A2STUR=0.   
        A2STUI=0.   
        A2USTR=0.   
        A2USTI=0.   
        A2TUSR=0.   
        A2TUSI=0.   
        A4STUR=0.   
        A4STUI=0.   
        DO 760 I=6,2*MSTP(1)    
C'''Only t-quarks yet included  
        SQMQ=PMAS(I,1)**2   
        EPSS=4.*SQMQ/SH 
        EPST=4.*SQMQ/TH 
        EPSU=4.*SQMQ/UH 
        EPSH=4.*SQMQ/SQMH   
        IF(EPSH.LT.1.E-6) GOTO 760  
        BESTU=0.5*(1.+SQRT(1.+EPSS*TH/UH))  
        BEUST=0.5*(1.+SQRT(1.+EPSU*SH/TH))  
        BETUS=0.5*(1.+SQRT(1.+EPST*UH/SH))  
        BEUTS=BESTU 
        BETSU=BEUST 
        BESUT=BETUS 
        W3STUR=PYI3AU(BESTU,EPSH,1)-PYI3AU(BESTU,EPSS,1)-   
     &  PYI3AU(BESTU,EPSU,1)    
        W3STUI=PYI3AU(BESTU,EPSH,2)-PYI3AU(BESTU,EPSS,2)-   
     &  PYI3AU(BESTU,EPSU,2)    
        W3SUTR=PYI3AU(BESUT,EPSH,1)-PYI3AU(BESUT,EPSS,1)-   
     &  PYI3AU(BESUT,EPST,1)    
        W3SUTI=PYI3AU(BESUT,EPSH,2)-PYI3AU(BESUT,EPSS,2)-   
     &  PYI3AU(BESUT,EPST,2)    
        W3TSUR=PYI3AU(BETSU,EPSH,1)-PYI3AU(BETSU,EPST,1)-   
     &  PYI3AU(BETSU,EPSU,1)    
        W3TSUI=PYI3AU(BETSU,EPSH,2)-PYI3AU(BETSU,EPST,2)-   
     &  PYI3AU(BETSU,EPSU,2)    
        W3TUSR=PYI3AU(BETUS,EPSH,1)-PYI3AU(BETUS,EPST,1)-   
     &  PYI3AU(BETUS,EPSS,1)    
        W3TUSI=PYI3AU(BETUS,EPSH,2)-PYI3AU(BETUS,EPST,2)-   
     &  PYI3AU(BETUS,EPSS,2)    
        W3USTR=PYI3AU(BEUST,EPSH,1)-PYI3AU(BEUST,EPSU,1)-   
     &  PYI3AU(BEUST,EPST,1)    
        W3USTI=PYI3AU(BEUST,EPSH,2)-PYI3AU(BEUST,EPSU,2)-   
     &  PYI3AU(BEUST,EPST,2)    
        W3UTSR=PYI3AU(BEUTS,EPSH,1)-PYI3AU(BEUTS,EPSU,1)-   
     &  PYI3AU(BEUTS,EPSS,1)    
        W3UTSI=PYI3AU(BEUTS,EPSH,2)-PYI3AU(BEUTS,EPSU,2)-   
     &  PYI3AU(BEUTS,EPSS,2)    
        B2STUR=SQMQ/SQMH**2*(SH*(UH-SH)/(SH+UH)+2.*TH*UH*(UH+2.*SH)/    
     &  (SH+UH)**2*(PYW1AU(EPST,1)-PYW1AU(EPSH,1))+(SQMQ-SH/4.)*    
     &  (0.5*PYW2AU(EPSS,1)+0.5*PYW2AU(EPSH,1)-PYW2AU(EPST,1)+W3STUR)+  
     &  SH**2*(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(PYW2AU(EPST,1)- 
     &  PYW2AU(EPSH,1))+0.5*TH*UH/SH*(PYW2AU(EPSH,1)-2.*PYW2AU(EPST,1))+    
     &  0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*W3TSUR) 
        B2STUI=SQMQ/SQMH**2*(2.*TH*UH*(UH+2.*SH)/(SH+UH)**2*    
     &  (PYW1AU(EPST,2)-PYW1AU(EPSH,2))+(SQMQ-SH/4.)*   
     &  (0.5*PYW2AU(EPSS,2)+0.5*PYW2AU(EPSH,2)-PYW2AU(EPST,2)+W3STUI)+  
     &  SH**2*(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(PYW2AU(EPST,2)- 
     &  PYW2AU(EPSH,2))+0.5*TH*UH/SH*(PYW2AU(EPSH,2)-2.*PYW2AU(EPST,2))+    
     &  0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*W3TSUI) 
        B2SUTR=SQMQ/SQMH**2*(SH*(TH-SH)/(SH+TH)+2.*UH*TH*(TH+2.*SH)/    
     &  (SH+TH)**2*(PYW1AU(EPSU,1)-PYW1AU(EPSH,1))+(SQMQ-SH/4.)*    
     &  (0.5*PYW2AU(EPSS,1)+0.5*PYW2AU(EPSH,1)-PYW2AU(EPSU,1)+W3SUTR)+  
     &  SH**2*(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(PYW2AU(EPSU,1)- 
     &  PYW2AU(EPSH,1))+0.5*UH*TH/SH*(PYW2AU(EPSH,1)-2.*PYW2AU(EPSU,1))+    
     &  0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*W3USTR) 
        B2SUTI=SQMQ/SQMH**2*(2.*UH*TH*(TH+2.*SH)/(SH+TH)**2*    
     &  (PYW1AU(EPSU,2)-PYW1AU(EPSH,2))+(SQMQ-SH/4.)*   
     &  (0.5*PYW2AU(EPSS,2)+0.5*PYW2AU(EPSH,2)-PYW2AU(EPSU,2)+W3SUTI)+  
     &  SH**2*(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(PYW2AU(EPSU,2)- 
     &  PYW2AU(EPSH,2))+0.5*UH*TH/SH*(PYW2AU(EPSH,2)-2.*PYW2AU(EPSU,2))+    
     &  0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*W3USTI) 
        B2TSUR=SQMQ/SQMH**2*(TH*(UH-TH)/(TH+UH)+2.*SH*UH*(UH+2.*TH)/    
     &  (TH+UH)**2*(PYW1AU(EPSS,1)-PYW1AU(EPSH,1))+(SQMQ-TH/4.)*    
     &  (0.5*PYW2AU(EPST,1)+0.5*PYW2AU(EPSH,1)-PYW2AU(EPSS,1)+W3TSUR)+  
     &  TH**2*(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(PYW2AU(EPSS,1)- 
     &  PYW2AU(EPSH,1))+0.5*SH*UH/TH*(PYW2AU(EPSH,1)-2.*PYW2AU(EPSS,1))+    
     &  0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*W3STUR) 
        B2TSUI=SQMQ/SQMH**2*(2.*SH*UH*(UH+2.*TH)/(TH+UH)**2*    
     &  (PYW1AU(EPSS,2)-PYW1AU(EPSH,2))+(SQMQ-TH/4.)*   
     &  (0.5*PYW2AU(EPST,2)+0.5*PYW2AU(EPSH,2)-PYW2AU(EPSS,2)+W3TSUI)+  
     &  TH**2*(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(PYW2AU(EPSS,2)- 
     &  PYW2AU(EPSH,2))+0.5*SH*UH/TH*(PYW2AU(EPSH,2)-2.*PYW2AU(EPSS,2))+    
     &  0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*W3STUI) 
        B2TUSR=SQMQ/SQMH**2*(TH*(SH-TH)/(TH+SH)+2.*UH*SH*(SH+2.*TH)/    
     &  (TH+SH)**2*(PYW1AU(EPSU,1)-PYW1AU(EPSH,1))+(SQMQ-TH/4.)*    
     &  (0.5*PYW2AU(EPST,1)+0.5*PYW2AU(EPSH,1)-PYW2AU(EPSU,1)+W3TUSR)+  
     &  TH**2*(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(PYW2AU(EPSU,1)- 
     &  PYW2AU(EPSH,1))+0.5*UH*SH/TH*(PYW2AU(EPSH,1)-2.*PYW2AU(EPSU,1))+    
     &  0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*W3UTSR) 
        B2TUSI=SQMQ/SQMH**2*(2.*UH*SH*(SH+2.*TH)/(TH+SH)**2*    
     &  (PYW1AU(EPSU,2)-PYW1AU(EPSH,2))+(SQMQ-TH/4.)*   
     &  (0.5*PYW2AU(EPST,2)+0.5*PYW2AU(EPSH,2)-PYW2AU(EPSU,2)+W3TUSI)+  
     &  TH**2*(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(PYW2AU(EPSU,2)- 
     &  PYW2AU(EPSH,2))+0.5*UH*SH/TH*(PYW2AU(EPSH,2)-2.*PYW2AU(EPSU,2))+    
     &  0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*W3UTSI) 
        B2USTR=SQMQ/SQMH**2*(UH*(TH-UH)/(UH+TH)+2.*SH*TH*(TH+2.*UH)/    
     &  (UH+TH)**2*(PYW1AU(EPSS,1)-PYW1AU(EPSH,1))+(SQMQ-UH/4.)*    
     &  (0.5*PYW2AU(EPSU,1)+0.5*PYW2AU(EPSH,1)-PYW2AU(EPSS,1)+W3USTR)+  
     &  UH**2*(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(PYW2AU(EPSS,1)- 
     &  PYW2AU(EPSH,1))+0.5*SH*TH/UH*(PYW2AU(EPSH,1)-2.*PYW2AU(EPSS,1))+    
     &  0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*W3SUTR) 
        B2USTI=SQMQ/SQMH**2*(2.*SH*TH*(TH+2.*UH)/(UH+TH)**2*    
     &  (PYW1AU(EPSS,2)-PYW1AU(EPSH,2))+(SQMQ-UH/4.)*   
     &  (0.5*PYW2AU(EPSU,2)+0.5*PYW2AU(EPSH,2)-PYW2AU(EPSS,2)+W3USTI)+  
     &  UH**2*(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(PYW2AU(EPSS,2)- 
     &  PYW2AU(EPSH,2))+0.5*SH*TH/UH*(PYW2AU(EPSH,2)-2.*PYW2AU(EPSS,2))+    
     &  0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*W3SUTI) 
        B2UTSR=SQMQ/SQMH**2*(UH*(SH-UH)/(UH+SH)+2.*TH*SH*(SH+2.*UH)/    
     &  (UH+SH)**2*(PYW1AU(EPST,1)-PYW1AU(EPSH,1))+(SQMQ-UH/4.)*    
     &  (0.5*PYW2AU(EPSU,1)+0.5*PYW2AU(EPSH,1)-PYW2AU(EPST,1)+W3UTSR)+  
     &  UH**2*(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(PYW2AU(EPST,1)- 
     &  PYW2AU(EPSH,1))+0.5*TH*SH/UH*(PYW2AU(EPSH,1)-2.*PYW2AU(EPST,1))+    
     &  0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*W3TUSR) 
        B2UTSI=SQMQ/SQMH**2*(2.*TH*SH*(SH+2.*UH)/(UH+SH)**2*    
     &  (PYW1AU(EPST,2)-PYW1AU(EPSH,2))+(SQMQ-UH/4.)*   
     &  (0.5*PYW2AU(EPSU,2)+0.5*PYW2AU(EPSH,2)-PYW2AU(EPST,2)+W3UTSI)+  
     &  UH**2*(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(PYW2AU(EPST,2)- 
     &  PYW2AU(EPSH,2))+0.5*TH*SH/UH*(PYW2AU(EPSH,2)-2.*PYW2AU(EPST,2))+    
     &  0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*W3TUSI) 
        B4STUR=SQMQ/SQMH*(-2./3.+(SQMQ/SQMH-1./4.)*(PYW2AU(EPSS,1)- 
     &  PYW2AU(EPSH,1)+W3STUR)) 
        B4STUI=SQMQ/SQMH*(SQMQ/SQMH-1./4.)*(PYW2AU(EPSS,2)- 
     &  PYW2AU(EPSH,2)+W3STUI)  
        B4TUSR=SQMQ/SQMH*(-2./3.+(SQMQ/SQMH-1./4.)*(PYW2AU(EPST,1)- 
     &  PYW2AU(EPSH,1)+W3TUSR)) 
        B4TUSI=SQMQ/SQMH*(SQMQ/SQMH-1./4.)*(PYW2AU(EPST,2)- 
     &  PYW2AU(EPSH,2)+W3TUSI)  
        B4USTR=SQMQ/SQMH*(-2./3.+(SQMQ/SQMH-1./4.)*(PYW2AU(EPSU,1)- 
     &  PYW2AU(EPSH,1)+W3USTR)) 
        B4USTI=SQMQ/SQMH*(SQMQ/SQMH-1./4.)*(PYW2AU(EPSU,2)- 
     &  PYW2AU(EPSH,2)+W3USTI)  
        A2STUR=A2STUR+B2STUR+B2SUTR 
        A2STUI=A2STUI+B2STUI+B2SUTI 
        A2USTR=A2USTR+B2USTR+B2UTSR 
        A2USTI=A2USTI+B2USTI+B2UTSI 
        A2TUSR=A2TUSR+B2TUSR+B2TSUR 
        A2TUSI=A2TUSI+B2TUSI+B2TSUI 
        A4STUR=A4STUR+B4STUR+B4USTR+B4TUSR  
        A4STUI=A4STUI+B4STUI+B4USTI+B4TUSI  
  760   CONTINUE    
        FACGH=COMFAC*FACA*3./(128.*PARU(1)**2)*AEM/XW*AS**3*    
     &  SQMH/SQMW*SQMH**3/(SH*TH*UH)*(A2STUR**2+A2STUI**2+A2USTR**2+    
     &  A2USTI**2+A2TUSR**2+A2TUSI**2+A4STUR**2+A4STUI**2)  
        FACGH=FACGH*WIDS(25,2)  
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 770 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGH    
  770   CONTINUE    
    
      ELSEIF(ISUB.EQ.114) THEN  
C...g + g -> gamma + gamma. 
        ASRE=0. 
        ASIM=0. 
        DO 780 I=1,2*MSTP(1)    
        EI=KCHG(IABS(I),1)/3.   
        SQMQ=PMAS(I,1)**2   
        EPSS=4.*SQMQ/SH 
        EPST=4.*SQMQ/TH 
        EPSU=4.*SQMQ/UH 
        IF(EPSS+ABS(EPST)+ABS(EPSU).LT.3.E-6) THEN  
          A0STUR=1.+(TH-UH)/SH*LOG(TH/UH)+0.5*(TH2+UH2)/SH2*    
     &    (LOG(TH/UH)**2+PARU(1)**2)    
          A0STUI=0. 
          A0TSUR=1.+(SH-UH)/TH*LOG(-SH/UH)+0.5*(SH2+UH2)/TH2*   
     &    LOG(-SH/UH)**2    
          A0TSUI=-PARU(1)*((SH-UH)/TH+(SH2+UH2)/TH2*LOG(-SH/UH))    
          A0UTSR=1.+(TH-SH)/UH*LOG(-TH/SH)+0.5*(TH2+SH2)/UH2*   
     &    LOG(-TH/SH)**2    
          A0UTSI=PARU(1)*((TH-SH)/UH+(TH2+SH2)/UH2*LOG(-TH/SH)) 
          A1STUR=-1.    
          A1STUI=0. 
          A2STUR=-1.    
          A2STUI=0. 
        ELSE    
          BESTU=0.5*(1.+SQRT(1.+EPSS*TH/UH))    
          BEUST=0.5*(1.+SQRT(1.+EPSU*SH/TH))    
          BETUS=0.5*(1.+SQRT(1.+EPST*UH/SH))    
          BEUTS=BESTU   
          BETSU=BEUST   
          BESUT=BETUS   
          A0STUR=1.+(1.+2.*TH/SH)*PYW1AU(EPST,1)+(1.+2.*UH/SH)* 
     &    PYW1AU(EPSU,1)+0.5*((TH2+UH2)/SH2-EPSS)*(PYW2AU(EPST,1)+  
     &    PYW2AU(EPSU,1))-0.25*EPST*(1.-0.5*EPSS)*(PYI3AU(BESUT,EPSS,1)+    
     &    PYI3AU(BESUT,EPST,1))-0.25*EPSU*(1.-0.5*EPSS)*    
     &    (PYI3AU(BESTU,EPSS,1)+PYI3AU(BESTU,EPSU,1))+  
     &    0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU)* 
     &    (PYI3AU(BETSU,EPST,1)+PYI3AU(BETSU,EPSU,1))   
          A0STUI=(1.+2.*TH/SH)*PYW1AU(EPST,2)+(1.+2.*UH/SH)*    
     &    PYW1AU(EPSU,2)+0.5*((TH2+UH2)/SH2-EPSS)*(PYW2AU(EPST,2)+  
     &    PYW2AU(EPSU,2))-0.25*EPST*(1.-0.5*EPSS)*(PYI3AU(BESUT,EPSS,2)+    
     &    PYI3AU(BESUT,EPST,2))-0.25*EPSU*(1.-0.5*EPSS)*    
     &    (PYI3AU(BESTU,EPSS,2)+PYI3AU(BESTU,EPSU,2))+  
     &    0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU)* 
     &    (PYI3AU(BETSU,EPST,2)+PYI3AU(BETSU,EPSU,2))   
          A0TSUR=1.+(1.+2.*SH/TH)*PYW1AU(EPSS,1)+(1.+2.*UH/TH)* 
     &    PYW1AU(EPSU,1)+0.5*((SH2+UH2)/TH2-EPST)*(PYW2AU(EPSS,1)+  
     &    PYW2AU(EPSU,1))-0.25*EPSS*(1.-0.5*EPST)*(PYI3AU(BETUS,EPST,1)+    
     &    PYI3AU(BETUS,EPSS,1))-0.25*EPSU*(1.-0.5*EPST)*    
     &    (PYI3AU(BETSU,EPST,1)+PYI3AU(BETSU,EPSU,1))+  
     &    0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU)* 
     &    (PYI3AU(BESTU,EPSS,1)+PYI3AU(BESTU,EPSU,1))   
          A0TSUI=(1.+2.*SH/TH)*PYW1AU(EPSS,2)+(1.+2.*UH/TH)*    
     &    PYW1AU(EPSU,2)+0.5*((SH2+UH2)/TH2-EPST)*(PYW2AU(EPSS,2)+  
     &    PYW2AU(EPSU,2))-0.25*EPSS*(1.-0.5*EPST)*(PYI3AU(BETUS,EPST,2)+    
     &    PYI3AU(BETUS,EPSS,2))-0.25*EPSU*(1.-0.5*EPST)*    
     &    (PYI3AU(BETSU,EPST,2)+PYI3AU(BETSU,EPSU,2))+  
     &    0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU)* 
     &    (PYI3AU(BESTU,EPSS,2)+PYI3AU(BESTU,EPSU,2))   
          A0UTSR=1.+(1.+2.*TH/UH)*PYW1AU(EPST,1)+(1.+2.*SH/UH)* 
     &    PYW1AU(EPSS,1)+0.5*((TH2+SH2)/UH2-EPSU)*(PYW2AU(EPST,1)+  
     &    PYW2AU(EPSS,1))-0.25*EPST*(1.-0.5*EPSU)*(PYI3AU(BEUST,EPSU,1)+    
     &    PYI3AU(BEUST,EPST,1))-0.25*EPSS*(1.-0.5*EPSU)*    
     &    (PYI3AU(BEUTS,EPSU,1)+PYI3AU(BEUTS,EPSS,1))+  
     &    0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS)* 
     &    (PYI3AU(BETUS,EPST,1)+PYI3AU(BETUS,EPSS,1))   
          A0UTSI=(1.+2.*TH/UH)*PYW1AU(EPST,2)+(1.+2.*SH/UH)*    
     &    PYW1AU(EPSS,2)+0.5*((TH2+SH2)/UH2-EPSU)*(PYW2AU(EPST,2)+  
     &    PYW2AU(EPSS,2))-0.25*EPST*(1.-0.5*EPSU)*(PYI3AU(BEUST,EPSU,2)+    
     &    PYI3AU(BEUST,EPST,2))-0.25*EPSS*(1.-0.5*EPSU)*    
     &    (PYI3AU(BEUTS,EPSU,2)+PYI3AU(BEUTS,EPSS,2))+  
     &    0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS)* 
     &    (PYI3AU(BETUS,EPST,2)+PYI3AU(BETUS,EPSS,2))   
          A1STUR=-1.-0.25*(EPSS+EPST+EPSU)*(PYW2AU(EPSS,1)+ 
     &    PYW2AU(EPST,1)+PYW2AU(EPSU,1))+0.25*(EPSU+0.5*EPSS*EPST)* 
     &    (PYI3AU(BESUT,EPSS,1)+PYI3AU(BESUT,EPST,1))+  
     &    0.25*(EPST+0.5*EPSS*EPSU)*(PYI3AU(BESTU,EPSS,1)+  
     &    PYI3AU(BESTU,EPSU,1))+0.25*(EPSS+0.5*EPST*EPSU)*  
     &    (PYI3AU(BETSU,EPST,1)+PYI3AU(BETSU,EPSU,1))   
          A1STUI=-0.25*(EPSS+EPST+EPSU)*(PYW2AU(EPSS,2)+PYW2AU(EPST,2)+ 
     &    PYW2AU(EPSU,2))+0.25*(EPSU+0.5*EPSS*EPST)*    
     &    (PYI3AU(BESUT,EPSS,2)+PYI3AU(BESUT,EPST,2))+  
     &    0.25*(EPST+0.5*EPSS*EPSU)*(PYI3AU(BESTU,EPSS,2)+  
     &    PYI3AU(BESTU,EPSU,2))+0.25*(EPSS+0.5*EPST*EPSU)*  
     &    (PYI3AU(BETSU,EPST,2)+PYI3AU(BETSU,EPSU,2))   
          A2STUR=-1.+0.125*EPSS*EPST*(PYI3AU(BESUT,EPSS,1)+ 
     &    PYI3AU(BESUT,EPST,1))+0.125*EPSS*EPSU*(PYI3AU(BESTU,EPSS,1)+  
     &    PYI3AU(BESTU,EPSU,1))+0.125*EPST*EPSU*(PYI3AU(BETSU,EPST,1)+  
     &    PYI3AU(BETSU,EPSU,1)) 
          A2STUI=0.125*EPSS*EPST*(PYI3AU(BESUT,EPSS,2)+ 
     &    PYI3AU(BESUT,EPST,2))+0.125*EPSS*EPSU*(PYI3AU(BESTU,EPSS,2)+  
     &    PYI3AU(BESTU,EPSU,2))+0.125*EPST*EPSU*(PYI3AU(BETSU,EPST,2)+  
     &    PYI3AU(BETSU,EPSU,2)) 
        ENDIF   
        ASRE=ASRE+EI**2*(A0STUR+A0TSUR+A0UTSR+4.*A1STUR+A2STUR) 
        ASIM=ASIM+EI**2*(A0STUI+A0TSUI+A0UTSI+4.*A1STUI+A2STUI) 
  780   CONTINUE    
        FACGG=COMFAC*FACA/(8.*PARU(1)**2)*AS**2*AEM**2*(ASRE**2+ASIM**2)    
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 790 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGG    
  790   CONTINUE    
    
      ELSEIF(ISUB.EQ.115) THEN  
C...g + g -> gamma + Z0.    
    
      ELSEIF(ISUB.EQ.116) THEN  
C...g + g -> Z0 + Z0.   
    
      ELSEIF(ISUB.EQ.117) THEN  
C...g + g -> W+ + W-.   
    
      ENDIF 
    
C...G: 2 -> 3, tree diagrams.   
    
      ELSEIF(ISUB.LE.140) THEN  
      IF(ISUB.EQ.121) THEN  
C...g + g -> f + fb + H0.   
    
      ENDIF 
    
C...H: 2 -> 1, tree diagrams, non-standard model processes. 
    
      ELSEIF(ISUB.LE.160) THEN  
      IF(ISUB.EQ.141) THEN  
C...f + fb -> gamma*/Z0/Z'0.    
        MINT(61)=2  
        CALL PYWIDT(32,SQRT(SH),WDTP,WDTE)  
        FACZP=COMFAC*AEM**2*4./9.   
        DO 800 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 800    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        API=SIGN(1.,EI) 
        VPI=API-4.*EI*XW    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZP*(EI**2*VINT(111)+EI*VI/(8.*XW*(1.-XW))*    
     &  SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)*VINT(112)+EI*VPI/(8.*XW*    
     &  (1.-XW))*SH*(SH-SQMZP)/((SH-SQMZP)**2+GMMZP**2)*VINT(113)+  
     &  (VI**2+AI**2)/(16.*XW*(1.-XW))**2*SH2/((SH-SQMZ)**2+GMMZ**2)*   
     &  VINT(114)+2.*(VI*VPI+AI*API)/(16.*XW*(1.-XW))**2*SH2*   
     &  ((SH-SQMZ)*(SH-SQMZP)+GMMZ*GMMZP)/(((SH-SQMZ)**2+GMMZ**2)*  
     &  ((SH-SQMZP)**2+GMMZP**2))*VINT(115)+(VPI**2+API**2)/    
     &  (16.*XW*(1.-XW))**2*SH2/((SH-SQMZP)**2+GMMZP**2)*VINT(116)) 
  800   CONTINUE    
    
      ELSEIF(ISUB.EQ.142) THEN  
C...f + fb' -> H+/-.    
        CALL PYWIDT(37,SQRT(SH),WDTP,WDTE)  
        FHC=COMFAC*(AEM/XW)**2*1./48.*(SH/SQMW)**2*SH2/ 
     &  ((SH-SQMHC)**2+GMMHC**2)    
C'''No construction yet for leptons 
        DO 840 I=1,MSTP(54)/2   
        IL=2*I-1    
        IU=2*I  
        RMQL=PMAS(IL,1)**2/SH   
        RMQU=PMAS(IU,1)**2/SH   
        FACHC=FHC*((RMQL*PARU(121)+RMQU/PARU(121))*(1.-RMQL-RMQU)-  
     &  4.*RMQL*RMQU)/SQRT(MAX(0.,(1.-RMQL-RMQU)**2-4.*RMQL*RMQU))  
        IF(KFAC(1,IL)*KFAC(2,-IU).EQ.0) GOTO 810    
        KCHHC=(KCHG(IL,1)-KCHG(IU,1))/3 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=IL 
        ISIG(NCHN,2)=-IU    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  810   IF(KFAC(1,-IL)*KFAC(2,IU).EQ.0) GOTO 820    
        KCHHC=(-KCHG(IL,1)+KCHG(IU,1))/3    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=-IL    
        ISIG(NCHN,2)=IU 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  820   IF(KFAC(1,IU)*KFAC(2,-IL).EQ.0) GOTO 830    
        KCHHC=(KCHG(IU,1)-KCHG(IL,1))/3 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=IU 
        ISIG(NCHN,2)=-IL    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  830   IF(KFAC(1,-IU)*KFAC(2,IL).EQ.0) GOTO 840    
        KCHHC=(-KCHG(IU,1)+KCHG(IL,1))/3    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=-IU    
        ISIG(NCHN,2)=IL 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  840   CONTINUE    
    
      ELSEIF(ISUB.EQ.143) THEN  
C...f + fb -> R.    
        CALL PYWIDT(40,SQRT(SH),WDTP,WDTE)  
        FACR=COMFAC*(AEM/XW)**2*1./9.*SH2/((SH-SQMR)**2+GMMR**2)    
        DO 860 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 860   
        IA=IABS(I)  
        DO 850 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 850   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.IABS(IA-JA).NE.2) GOTO 850   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACR*(WDTE(0,1)+WDTE(0,(10-(I+J))/4)+WDTE(0,4))  
  850   CONTINUE    
  860   CONTINUE    
    
      ENDIF 
    
C...I: 2 -> 2, tree diagrams, non-standard model processes. 
    
      ELSE  
      IF(ISUB.EQ.161) THEN  
C...f + g -> f' + H+/- (q + g -> q' + H+/- only).   
        FHCQ=COMFAC*FACA*AS*AEM/XW*1./24    
        DO 900 I=1,MSTP(54) 
        IU=I+MOD(I,2)   
        SQMQ=PMAS(IU,1)**2  
        FACHCQ=FHCQ/PARU(121)*SQMQ/SQMW*(SH/(SQMQ-UH)+  
     &  2.*SQMQ*(SQMHC-UH)/(SQMQ-UH)**2+(SQMQ-UH)/SH+   
     &  2.*SQMQ/(SQMQ-UH)+2.*(SQMHC-UH)/(SQMQ-UH)*(SQMHC-SQMQ-SH)/SH)   
        IF(KFAC(1,-I)*KFAC(2,21).EQ.0) GOTO 870 
        KCHHC=ISIGN(1,-KCHG(I,1))   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=-I 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  870   IF(KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 880  
        KCHHC=ISIGN(1,KCHG(I,1))    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  880   IF(KFAC(1,21)*KFAC(2,-I).EQ.0) GOTO 890 
        KCHHC=ISIGN(1,-KCHG(I,1))   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  890   IF(KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 900  
        KCHHC=ISIGN(1,KCHG(I,1))    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=I  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  900   CONTINUE    
    
      ENDIF 
      ENDIF 
    
C...Multiply with structure functions.  
      IF(ISUB.LE.90.OR.ISUB.GE.96) THEN 
        DO 910 ICHN=1,NCHN  
        IF(MINT(41).EQ.2) THEN  
          KFL1=ISIG(ICHN,1) 
          IF(KFL1.EQ.21) KFL1=0 
          SIGH(ICHN)=SIGH(ICHN)*XSFX(1,KFL1)    
        ENDIF   
        IF(MINT(42).EQ.2) THEN  
          KFL2=ISIG(ICHN,2) 
          IF(KFL2.EQ.21) KFL2=0 
          SIGH(ICHN)=SIGH(ICHN)*XSFX(2,KFL2)    
        ENDIF   
  910   SIGS=SIGS+SIGH(ICHN)    
      ENDIF 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYSTFU(KF,X,Q2,XPQ,JBT)    
    
C			*******JBT specifies beam or target of the particle
C...Gives proton and pi+ parton structure functions according to a few  
C...different parametrizations. Note that what is coded is x times the  
C...probability distribution, i.e. xq(x,Q2) etc.    
      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/
      COMMON/HIJCRDN/YP(3,300),YT(3,300)
      SAVE  /HIJCRDN/
C			********COMMON BLOCK FROM HIJING
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      DIMENSION XPQ(-6:6),XQ(6),TX(6),TT(6),TS(6),NEHLQ(8,2),   
     &CEHLQ(6,6,2,8,2),CDO(3,6,5,2),COW(3,5,4,2)    
    
C...The following data lines are coefficients needed in the 
C...Eichten, Hinchliffe, Lane, Quigg proton structure function  
C...parametrizations, see below.    
C...Powers of 1-x in different cases.   
      DATA NEHLQ/3,4,7,5,7,7,7,7,3,4,7,6,7,7,7,7/   
C...Expansion coefficients for up valence quark distribution.   
      DATA (((CEHLQ(IX,IT,NX,1,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 7.677E-01,-2.087E-01,-3.303E-01,-2.517E-02,-1.570E-02,-1.000E-04,    
     2-5.326E-01,-2.661E-01, 3.201E-01, 1.192E-01, 2.434E-02, 7.620E-03,    
     3 2.162E-01, 1.881E-01,-8.375E-02,-6.515E-02,-1.743E-02,-5.040E-03,    
     4-9.211E-02,-9.952E-02, 1.373E-02, 2.506E-02, 8.770E-03, 2.550E-03,    
     5 3.670E-02, 4.409E-02, 9.600E-04,-7.960E-03,-3.420E-03,-1.050E-03,    
     6-1.549E-02,-2.026E-02,-3.060E-03, 2.220E-03, 1.240E-03, 4.100E-04,    
     1 2.395E-01, 2.905E-01, 9.778E-02, 2.149E-02, 3.440E-03, 5.000E-04,    
     2 1.751E-02,-6.090E-03,-2.687E-02,-1.916E-02,-7.970E-03,-2.750E-03,    
     3-5.760E-03,-5.040E-03, 1.080E-03, 2.490E-03, 1.530E-03, 7.500E-04,    
     4 1.740E-03, 1.960E-03, 3.000E-04,-3.400E-04,-2.900E-04,-1.800E-04,    
     5-5.300E-04,-6.400E-04,-1.700E-04, 4.000E-05, 6.000E-05, 4.000E-05,    
     6 1.700E-04, 2.200E-04, 8.000E-05, 1.000E-05,-1.000E-05,-1.000E-05/    
      DATA (((CEHLQ(IX,IT,NX,1,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 7.237E-01,-2.189E-01,-2.995E-01,-1.909E-02,-1.477E-02, 2.500E-04,    
     2-5.314E-01,-2.425E-01, 3.283E-01, 1.119E-01, 2.223E-02, 7.070E-03,    
     3 2.289E-01, 1.890E-01,-9.859E-02,-6.900E-02,-1.747E-02,-5.080E-03,    
     4-1.041E-01,-1.084E-01, 2.108E-02, 2.975E-02, 9.830E-03, 2.830E-03,    
     5 4.394E-02, 5.116E-02,-1.410E-03,-1.055E-02,-4.230E-03,-1.270E-03,    
     6-1.991E-02,-2.539E-02,-2.780E-03, 3.430E-03, 1.720E-03, 5.500E-04,    
     1 2.410E-01, 2.884E-01, 9.369E-02, 1.900E-02, 2.530E-03, 2.400E-04,    
     2 1.765E-02,-9.220E-03,-3.037E-02,-2.085E-02,-8.440E-03,-2.810E-03,    
     3-6.450E-03,-5.260E-03, 1.720E-03, 3.110E-03, 1.830E-03, 8.700E-04,    
     4 2.120E-03, 2.320E-03, 2.600E-04,-4.900E-04,-3.900E-04,-2.300E-04,    
     5-6.900E-04,-8.200E-04,-2.000E-04, 7.000E-05, 9.000E-05, 6.000E-05,    
     6 2.400E-04, 3.100E-04, 1.100E-04, 0.000E+00,-2.000E-05,-2.000E-05/    
C...Expansion coefficients for down valence quark distribution. 
      DATA (((CEHLQ(IX,IT,NX,2,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 3.813E-01,-8.090E-02,-1.634E-01,-2.185E-02,-8.430E-03,-6.200E-04,    
     2-2.948E-01,-1.435E-01, 1.665E-01, 6.638E-02, 1.473E-02, 4.080E-03,    
     3 1.252E-01, 1.042E-01,-4.722E-02,-3.683E-02,-1.038E-02,-2.860E-03,    
     4-5.478E-02,-5.678E-02, 8.900E-03, 1.484E-02, 5.340E-03, 1.520E-03,    
     5 2.220E-02, 2.567E-02,-3.000E-05,-4.970E-03,-2.160E-03,-6.500E-04,    
     6-9.530E-03,-1.204E-02,-1.510E-03, 1.510E-03, 8.300E-04, 2.700E-04,    
     1 1.261E-01, 1.354E-01, 3.958E-02, 8.240E-03, 1.660E-03, 4.500E-04,    
     2 3.890E-03,-1.159E-02,-1.625E-02,-9.610E-03,-3.710E-03,-1.260E-03,    
     3-1.910E-03,-5.600E-04, 1.590E-03, 1.590E-03, 8.400E-04, 3.900E-04,    
     4 6.400E-04, 4.900E-04,-1.500E-04,-2.900E-04,-1.800E-04,-1.000E-04,    
     5-2.000E-04,-1.900E-04, 0.000E+00, 6.000E-05, 4.000E-05, 3.000E-05,    
     6 7.000E-05, 8.000E-05, 2.000E-05,-1.000E-05,-1.000E-05,-1.000E-05/    
      DATA (((CEHLQ(IX,IT,NX,2,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 3.578E-01,-8.622E-02,-1.480E-01,-1.840E-02,-7.820E-03,-4.500E-04,    
     2-2.925E-01,-1.304E-01, 1.696E-01, 6.243E-02, 1.353E-02, 3.750E-03,    
     3 1.318E-01, 1.041E-01,-5.486E-02,-3.872E-02,-1.038E-02,-2.850E-03,    
     4-6.162E-02,-6.143E-02, 1.303E-02, 1.740E-02, 5.940E-03, 1.670E-03,    
     5 2.643E-02, 2.957E-02,-1.490E-03,-6.450E-03,-2.630E-03,-7.700E-04,    
     6-1.218E-02,-1.497E-02,-1.260E-03, 2.240E-03, 1.120E-03, 3.500E-04,    
     1 1.263E-01, 1.334E-01, 3.732E-02, 7.070E-03, 1.260E-03, 3.400E-04,    
     2 3.660E-03,-1.357E-02,-1.795E-02,-1.031E-02,-3.880E-03,-1.280E-03,    
     3-2.100E-03,-3.600E-04, 2.050E-03, 1.920E-03, 9.800E-04, 4.400E-04,    
     4 7.700E-04, 5.400E-04,-2.400E-04,-3.900E-04,-2.400E-04,-1.300E-04,    
     5-2.600E-04,-2.300E-04, 2.000E-05, 9.000E-05, 6.000E-05, 4.000E-05,    
     6 9.000E-05, 1.000E-04, 2.000E-05,-2.000E-05,-2.000E-05,-1.000E-05/    
C...Expansion coefficients for up and down sea quark distributions. 
      DATA (((CEHLQ(IX,IT,NX,3,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 6.870E-02,-6.861E-02, 2.973E-02,-5.400E-03, 3.780E-03,-9.700E-04,    
     2-1.802E-02, 1.400E-04, 6.490E-03,-8.540E-03, 1.220E-03,-1.750E-03,    
     3-4.650E-03, 1.480E-03,-5.930E-03, 6.000E-04,-1.030E-03,-8.000E-05,    
     4 6.440E-03, 2.570E-03, 2.830E-03, 1.150E-03, 7.100E-04, 3.300E-04,    
     5-3.930E-03,-2.540E-03,-1.160E-03,-7.700E-04,-3.600E-04,-1.900E-04,    
     6 2.340E-03, 1.930E-03, 5.300E-04, 3.700E-04, 1.600E-04, 9.000E-05,    
     1 1.014E+00,-1.106E+00, 3.374E-01,-7.444E-02, 8.850E-03,-8.700E-04,    
     2 9.233E-01,-1.285E+00, 4.475E-01,-9.786E-02, 1.419E-02,-1.120E-03,    
     3 4.888E-02,-1.271E-01, 8.606E-02,-2.608E-02, 4.780E-03,-6.000E-04,    
     4-2.691E-02, 4.887E-02,-1.771E-02, 1.620E-03, 2.500E-04,-6.000E-05,    
     5 7.040E-03,-1.113E-02, 1.590E-03, 7.000E-04,-2.000E-04, 0.000E+00,    
     6-1.710E-03, 2.290E-03, 3.800E-04,-3.500E-04, 4.000E-05, 1.000E-05/    
      DATA (((CEHLQ(IX,IT,NX,3,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 1.008E-01,-7.100E-02, 1.973E-02,-5.710E-03, 2.930E-03,-9.900E-04,    
     2-5.271E-02,-1.823E-02, 1.792E-02,-6.580E-03, 1.750E-03,-1.550E-03,    
     3 1.220E-02, 1.763E-02,-8.690E-03,-8.800E-04,-1.160E-03,-2.100E-04,    
     4-1.190E-03,-7.180E-03, 2.360E-03, 1.890E-03, 7.700E-04, 4.100E-04,    
     5-9.100E-04, 2.040E-03,-3.100E-04,-1.050E-03,-4.000E-04,-2.400E-04,    
     6 1.190E-03,-1.700E-04,-2.000E-04, 4.200E-04, 1.700E-04, 1.000E-04,    
     1 1.081E+00,-1.189E+00, 3.868E-01,-8.617E-02, 1.115E-02,-1.180E-03,    
     2 9.917E-01,-1.396E+00, 4.998E-01,-1.159E-01, 1.674E-02,-1.720E-03,    
     3 5.099E-02,-1.338E-01, 9.173E-02,-2.885E-02, 5.890E-03,-6.500E-04,    
     4-3.178E-02, 5.703E-02,-2.070E-02, 2.440E-03, 1.100E-04,-9.000E-05,    
     5 8.970E-03,-1.392E-02, 2.050E-03, 6.500E-04,-2.300E-04, 2.000E-05,    
     6-2.340E-03, 3.010E-03, 5.000E-04,-3.900E-04, 6.000E-05, 1.000E-05/    
C...Expansion coefficients for gluon distribution.  
      DATA (((CEHLQ(IX,IT,NX,4,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 9.482E-01,-9.578E-01, 1.009E-01,-1.051E-01, 3.456E-02,-3.054E-02,    
     2-9.627E-01, 5.379E-01, 3.368E-01,-9.525E-02, 1.488E-02,-2.051E-02,    
     3 4.300E-01,-8.306E-02,-3.372E-01, 4.902E-02,-9.160E-03, 1.041E-02,    
     4-1.925E-01,-1.790E-02, 2.183E-01, 7.490E-03, 4.140E-03,-1.860E-03,    
     5 8.183E-02, 1.926E-02,-1.072E-01,-1.944E-02,-2.770E-03,-5.200E-04,    
     6-3.884E-02,-1.234E-02, 5.410E-02, 1.879E-02, 3.350E-03, 1.040E-03,    
     1 2.948E+01,-3.902E+01, 1.464E+01,-3.335E+00, 5.054E-01,-5.915E-02,    
     2 2.559E+01,-3.955E+01, 1.661E+01,-4.299E+00, 6.904E-01,-8.243E-02,    
     3-1.663E+00, 1.176E+00, 1.118E+00,-7.099E-01, 1.948E-01,-2.404E-02,    
     4-2.168E-01, 8.170E-01,-7.169E-01, 1.851E-01,-1.924E-02,-3.250E-03,    
     5 2.088E-01,-4.355E-01, 2.239E-01,-2.446E-02,-3.620E-03, 1.910E-03,    
     6-9.097E-02, 1.601E-01,-5.681E-02,-2.500E-03, 2.580E-03,-4.700E-04/    
      DATA (((CEHLQ(IX,IT,NX,4,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 2.367E+00, 4.453E-01, 3.660E-01, 9.467E-02, 1.341E-01, 1.661E-02,    
     2-3.170E+00,-1.795E+00, 3.313E-02,-2.874E-01,-9.827E-02,-7.119E-02,    
     3 1.823E+00, 1.457E+00,-2.465E-01, 3.739E-02, 6.090E-03, 1.814E-02,    
     4-1.033E+00,-9.827E-01, 2.136E-01, 1.169E-01, 5.001E-02, 1.684E-02,    
     5 5.133E-01, 5.259E-01,-1.173E-01,-1.139E-01,-4.988E-02,-2.021E-02,    
     6-2.881E-01,-3.145E-01, 5.667E-02, 9.161E-02, 4.568E-02, 1.951E-02,    
     1 3.036E+01,-4.062E+01, 1.578E+01,-3.699E+00, 6.020E-01,-7.031E-02,    
     2 2.700E+01,-4.167E+01, 1.770E+01,-4.804E+00, 7.862E-01,-1.060E-01,    
     3-1.909E+00, 1.357E+00, 1.127E+00,-7.181E-01, 2.232E-01,-2.481E-02,    
     4-2.488E-01, 9.781E-01,-8.127E-01, 2.094E-01,-2.997E-02,-4.710E-03,    
     5 2.506E-01,-5.427E-01, 2.672E-01,-3.103E-02,-1.800E-03, 2.870E-03,    
     6-1.128E-01, 2.087E-01,-6.972E-02,-2.480E-03, 2.630E-03,-8.400E-04/    
C...Expansion coefficients for strange sea quark distribution.  
      DATA (((CEHLQ(IX,IT,NX,5,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 4.968E-02,-4.173E-02, 2.102E-02,-3.270E-03, 3.240E-03,-6.700E-04,    
     2-6.150E-03,-1.294E-02, 6.740E-03,-6.890E-03, 9.000E-04,-1.510E-03,    
     3-8.580E-03, 5.050E-03,-4.900E-03,-1.600E-04,-9.400E-04,-1.500E-04,    
     4 7.840E-03, 1.510E-03, 2.220E-03, 1.400E-03, 7.000E-04, 3.500E-04,    
     5-4.410E-03,-2.220E-03,-8.900E-04,-8.500E-04,-3.600E-04,-2.000E-04,    
     6 2.520E-03, 1.840E-03, 4.100E-04, 3.900E-04, 1.600E-04, 9.000E-05,    
     1 9.235E-01,-1.085E+00, 3.464E-01,-7.210E-02, 9.140E-03,-9.100E-04,    
     2 9.315E-01,-1.274E+00, 4.512E-01,-9.775E-02, 1.380E-02,-1.310E-03,    
     3 4.739E-02,-1.296E-01, 8.482E-02,-2.642E-02, 4.760E-03,-5.700E-04,    
     4-2.653E-02, 4.953E-02,-1.735E-02, 1.750E-03, 2.800E-04,-6.000E-05,    
     5 6.940E-03,-1.132E-02, 1.480E-03, 6.500E-04,-2.100E-04, 0.000E+00,    
     6-1.680E-03, 2.340E-03, 4.200E-04,-3.400E-04, 5.000E-05, 1.000E-05/    
      DATA (((CEHLQ(IX,IT,NX,5,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 6.478E-02,-4.537E-02, 1.643E-02,-3.490E-03, 2.710E-03,-6.700E-04,    
     2-2.223E-02,-2.126E-02, 1.247E-02,-6.290E-03, 1.120E-03,-1.440E-03,    
     3-1.340E-03, 1.362E-02,-6.130E-03,-7.900E-04,-9.000E-04,-2.000E-04,    
     4 5.080E-03,-3.610E-03, 1.700E-03, 1.830E-03, 6.800E-04, 4.000E-04,    
     5-3.580E-03, 6.000E-05,-2.600E-04,-1.050E-03,-3.800E-04,-2.300E-04,    
     6 2.420E-03, 9.300E-04,-1.000E-04, 4.500E-04, 1.700E-04, 1.100E-04,    
     1 9.868E-01,-1.171E+00, 3.940E-01,-8.459E-02, 1.124E-02,-1.250E-03,    
     2 1.001E+00,-1.383E+00, 5.044E-01,-1.152E-01, 1.658E-02,-1.830E-03,    
     3 4.928E-02,-1.368E-01, 9.021E-02,-2.935E-02, 5.800E-03,-6.600E-04,    
     4-3.133E-02, 5.785E-02,-2.023E-02, 2.630E-03, 1.600E-04,-8.000E-05,    
     5 8.840E-03,-1.416E-02, 1.900E-03, 5.800E-04,-2.500E-04, 1.000E-05,    
     6-2.300E-03, 3.080E-03, 5.500E-04,-3.700E-04, 7.000E-05, 1.000E-05/    
C...Expansion coefficients for charm sea quark distribution.    
      DATA (((CEHLQ(IX,IT,NX,6,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 9.270E-03,-1.817E-02, 9.590E-03,-6.390E-03, 1.690E-03,-1.540E-03,    
     2 5.710E-03,-1.188E-02, 6.090E-03,-4.650E-03, 1.240E-03,-1.310E-03,    
     3-3.960E-03, 7.100E-03,-3.590E-03, 1.840E-03,-3.900E-04, 3.400E-04,    
     4 1.120E-03,-1.960E-03, 1.120E-03,-4.800E-04, 1.000E-04,-4.000E-05,    
     5 4.000E-05,-3.000E-05,-1.800E-04, 9.000E-05,-5.000E-05,-2.000E-05,    
     6-4.200E-04, 7.300E-04,-1.600E-04, 5.000E-05, 5.000E-05, 5.000E-05,    
     1 8.098E-01,-1.042E+00, 3.398E-01,-6.824E-02, 8.760E-03,-9.000E-04,    
     2 8.961E-01,-1.217E+00, 4.339E-01,-9.287E-02, 1.304E-02,-1.290E-03,    
     3 3.058E-02,-1.040E-01, 7.604E-02,-2.415E-02, 4.600E-03,-5.000E-04,    
     4-2.451E-02, 4.432E-02,-1.651E-02, 1.430E-03, 1.200E-04,-1.000E-04,    
     5 1.122E-02,-1.457E-02, 2.680E-03, 5.800E-04,-1.200E-04, 3.000E-05,    
     6-7.730E-03, 7.330E-03,-7.600E-04,-2.400E-04, 1.000E-05, 0.000E+00/    
      DATA (((CEHLQ(IX,IT,NX,6,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 9.980E-03,-1.945E-02, 1.055E-02,-6.870E-03, 1.860E-03,-1.560E-03,    
     2 5.700E-03,-1.203E-02, 6.250E-03,-4.860E-03, 1.310E-03,-1.370E-03,    
     3-4.490E-03, 7.990E-03,-4.170E-03, 2.050E-03,-4.400E-04, 3.300E-04,    
     4 1.470E-03,-2.480E-03, 1.460E-03,-5.700E-04, 1.200E-04,-1.000E-05,    
     5-9.000E-05, 1.500E-04,-3.200E-04, 1.200E-04,-6.000E-05,-4.000E-05,    
     6-4.200E-04, 7.600E-04,-1.400E-04, 4.000E-05, 7.000E-05, 5.000E-05,    
     1 8.698E-01,-1.131E+00, 3.836E-01,-8.111E-02, 1.048E-02,-1.300E-03,    
     2 9.626E-01,-1.321E+00, 4.854E-01,-1.091E-01, 1.583E-02,-1.700E-03,    
     3 3.057E-02,-1.088E-01, 8.022E-02,-2.676E-02, 5.590E-03,-5.600E-04,    
     4-2.845E-02, 5.164E-02,-1.918E-02, 2.210E-03,-4.000E-05,-1.500E-04,    
     5 1.311E-02,-1.751E-02, 3.310E-03, 5.100E-04,-1.200E-04, 5.000E-05,    
     6-8.590E-03, 8.380E-03,-9.200E-04,-2.600E-04, 1.000E-05,-1.000E-05/    
C...Expansion coefficients for bottom sea quark distribution.   
      DATA (((CEHLQ(IX,IT,NX,7,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 9.010E-03,-1.401E-02, 7.150E-03,-4.130E-03, 1.260E-03,-1.040E-03,    
     2 6.280E-03,-9.320E-03, 4.780E-03,-2.890E-03, 9.100E-04,-8.200E-04,    
     3-2.930E-03, 4.090E-03,-1.890E-03, 7.600E-04,-2.300E-04, 1.400E-04,    
     4 3.900E-04,-1.200E-03, 4.400E-04,-2.500E-04, 2.000E-05,-2.000E-05,    
     5 2.600E-04, 1.400E-04,-8.000E-05, 1.000E-04, 1.000E-05, 1.000E-05,    
     6-2.600E-04, 3.200E-04, 1.000E-05,-1.000E-05, 1.000E-05,-1.000E-05,    
     1 8.029E-01,-1.075E+00, 3.792E-01,-7.843E-02, 1.007E-02,-1.090E-03,    
     2 7.903E-01,-1.099E+00, 4.153E-01,-9.301E-02, 1.317E-02,-1.410E-03,    
     3-1.704E-02,-1.130E-02, 2.882E-02,-1.341E-02, 3.040E-03,-3.600E-04,    
     4-7.200E-04, 7.230E-03,-5.160E-03, 1.080E-03,-5.000E-05,-4.000E-05,    
     5 3.050E-03,-4.610E-03, 1.660E-03,-1.300E-04,-1.000E-05, 1.000E-05,    
     6-4.360E-03, 5.230E-03,-1.610E-03, 2.000E-04,-2.000E-05, 0.000E+00/    
      DATA (((CEHLQ(IX,IT,NX,7,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 8.980E-03,-1.459E-02, 7.510E-03,-4.410E-03, 1.310E-03,-1.070E-03,    
     2 5.970E-03,-9.440E-03, 4.800E-03,-3.020E-03, 9.100E-04,-8.500E-04,    
     3-3.050E-03, 4.440E-03,-2.100E-03, 8.500E-04,-2.400E-04, 1.400E-04,    
     4 5.300E-04,-1.300E-03, 5.600E-04,-2.700E-04, 3.000E-05,-2.000E-05,    
     5 2.000E-04, 1.400E-04,-1.100E-04, 1.000E-04, 0.000E+00, 0.000E+00,    
     6-2.600E-04, 3.200E-04, 0.000E+00,-3.000E-05, 1.000E-05,-1.000E-05,    
     1 8.672E-01,-1.174E+00, 4.265E-01,-9.252E-02, 1.244E-02,-1.460E-03,    
     2 8.500E-01,-1.194E+00, 4.630E-01,-1.083E-01, 1.614E-02,-1.830E-03,    
     3-2.241E-02,-5.630E-03, 2.815E-02,-1.425E-02, 3.520E-03,-4.300E-04,    
     4-7.300E-04, 8.030E-03,-5.780E-03, 1.380E-03,-1.300E-04,-4.000E-05,    
     5 3.460E-03,-5.380E-03, 1.960E-03,-2.100E-04, 1.000E-05, 1.000E-05,    
     6-4.850E-03, 5.950E-03,-1.890E-03, 2.600E-04,-3.000E-05, 0.000E+00/    
C...Expansion coefficients for top sea quark distribution.  
      DATA (((CEHLQ(IX,IT,NX,8,1),IX=1,6),IT=1,6),NX=1,2)/  
     1 4.410E-03,-7.480E-03, 3.770E-03,-2.580E-03, 7.300E-04,-7.100E-04,    
     2 3.840E-03,-6.050E-03, 3.030E-03,-2.030E-03, 5.800E-04,-5.900E-04,    
     3-8.800E-04, 1.660E-03,-7.500E-04, 4.700E-04,-1.000E-04, 1.000E-04,    
     4-8.000E-05,-1.500E-04, 1.200E-04,-9.000E-05, 3.000E-05, 0.000E+00,    
     5 1.300E-04,-2.200E-04,-2.000E-05,-2.000E-05,-2.000E-05,-2.000E-05,    
     6-7.000E-05, 1.900E-04,-4.000E-05, 2.000E-05, 0.000E+00, 0.000E+00,    
     1 6.623E-01,-9.248E-01, 3.519E-01,-7.930E-02, 1.110E-02,-1.180E-03,    
     2 6.380E-01,-9.062E-01, 3.582E-01,-8.479E-02, 1.265E-02,-1.390E-03,    
     3-2.581E-02, 2.125E-02, 4.190E-03,-4.980E-03, 1.490E-03,-2.100E-04,    
     4 7.100E-04, 5.300E-04,-1.270E-03, 3.900E-04,-5.000E-05,-1.000E-05,    
     5 3.850E-03,-5.060E-03, 1.860E-03,-3.500E-04, 4.000E-05, 0.000E+00,    
     6-3.530E-03, 4.460E-03,-1.500E-03, 2.700E-04,-3.000E-05, 0.000E+00/    
      DATA (((CEHLQ(IX,IT,NX,8,2),IX=1,6),IT=1,6),NX=1,2)/  
     1 4.260E-03,-7.530E-03, 3.830E-03,-2.680E-03, 7.600E-04,-7.300E-04,    
     2 3.640E-03,-6.050E-03, 3.030E-03,-2.090E-03, 5.900E-04,-6.000E-04,    
     3-9.200E-04, 1.710E-03,-8.200E-04, 5.000E-04,-1.200E-04, 1.000E-04,    
     4-5.000E-05,-1.600E-04, 1.300E-04,-9.000E-05, 3.000E-05, 0.000E+00,    
     5 1.300E-04,-2.100E-04,-1.000E-05,-2.000E-05,-2.000E-05,-1.000E-05,    
     6-8.000E-05, 1.800E-04,-5.000E-05, 2.000E-05, 0.000E+00, 0.000E+00,    
     1 7.146E-01,-1.007E+00, 3.932E-01,-9.246E-02, 1.366E-02,-1.540E-03,    
     2 6.856E-01,-9.828E-01, 3.977E-01,-9.795E-02, 1.540E-02,-1.790E-03,    
     3-3.053E-02, 2.758E-02, 2.150E-03,-4.880E-03, 1.640E-03,-2.500E-04,    
     4 9.200E-04, 4.200E-04,-1.340E-03, 4.600E-04,-8.000E-05,-1.000E-05,    
     5 4.230E-03,-5.660E-03, 2.140E-03,-4.300E-04, 6.000E-05, 0.000E+00,    
     6-3.890E-03, 5.000E-03,-1.740E-03, 3.300E-04,-4.000E-05, 0.000E+00/    
    
C...The following data lines are coefficients needed in the 
C...Duke, Owens proton structure function parametrizations, see below.  
C...Expansion coefficients for (up+down) valence quark distribution.    
      DATA ((CDO(IP,IS,1,1),IS=1,6),IP=1,3)/    
     1 4.190E-01, 3.460E+00, 4.400E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2 4.000E-03, 7.240E-01,-4.860E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     3-7.000E-03,-6.600E-02, 1.330E+00, 0.000E+00, 0.000E+00, 0.000E+00/    
      DATA ((CDO(IP,IS,1,2),IS=1,6),IP=1,3)/    
     1 3.740E-01, 3.330E+00, 6.030E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2 1.400E-02, 7.530E-01,-6.220E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     3 0.000E+00,-7.600E-02, 1.560E+00, 0.000E+00, 0.000E+00, 0.000E+00/    
C...Expansion coefficients for down valence quark distribution. 
      DATA ((CDO(IP,IS,2,1),IS=1,6),IP=1,3)/    
     1 7.630E-01, 4.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2-2.370E-01, 6.270E-01,-4.210E-01, 0.000E+00, 0.000E+00, 0.000E+00,    
     3 2.600E-02,-1.900E-02, 3.300E-02, 0.000E+00, 0.000E+00, 0.000E+00/    
      DATA ((CDO(IP,IS,2,2),IS=1,6),IP=1,3)/    
     1 7.610E-01, 3.830E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2-2.320E-01, 6.270E-01,-4.180E-01, 0.000E+00, 0.000E+00, 0.000E+00,    
     3 2.300E-02,-1.900E-02, 3.600E-02, 0.000E+00, 0.000E+00, 0.000E+00/    
C...Expansion coefficients for (up+down+strange) sea quark distribution.    
      DATA ((CDO(IP,IS,3,1),IS=1,6),IP=1,3)/    
     1 1.265E+00, 0.000E+00, 8.050E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2-1.132E+00,-3.720E-01, 1.590E+00, 6.310E+00,-1.050E+01, 1.470E+01,    
     3 2.930E-01,-2.900E-02,-1.530E-01,-2.730E-01,-3.170E+00, 9.800E+00/    
      DATA ((CDO(IP,IS,3,2),IS=1,6),IP=1,3)/    
     1 1.670E+00, 0.000E+00, 9.150E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2-1.920E+00,-2.730E-01, 5.300E-01, 1.570E+01,-1.010E+02, 2.230E+02,    
     3 5.820E-01,-1.640E-01,-7.630E-01,-2.830E+00, 4.470E+01,-1.170E+02/    
C...Expansion coefficients for charm sea quark distribution.    
      DATA ((CDO(IP,IS,4,1),IS=1,6),IP=1,3)/    
     1 0.000E+00,-3.600E-02, 6.350E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2 1.350E-01,-2.220E-01, 3.260E+00,-3.030E+00, 1.740E+01,-1.790E+01,    
     3-7.500E-02,-5.800E-02,-9.090E-01, 1.500E+00,-1.130E+01, 1.560E+01/    
       DATA ((CDO(IP,IS,4,2),IS=1,6),IP=1,3)/   
     1 0.000E+00,-1.200E-01, 3.510E+00, 0.000E+00, 0.000E+00, 0.000E+00,    
     2 6.700E-02,-2.330E-01, 3.660E+00,-4.740E-01, 9.500E+00,-1.660E+01,    
     3-3.100E-02,-2.300E-02,-4.530E-01, 3.580E-01,-5.430E+00, 1.550E+01/    
C...Expansion coefficients for gluon distribution.  
      DATA ((CDO(IP,IS,5,1),IS=1,6),IP=1,3)/    
     1 1.560E+00, 0.000E+00, 6.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,    
     2-1.710E+00,-9.490E-01, 1.440E+00,-7.190E+00,-1.650E+01, 1.530E+01,    
     3 6.380E-01, 3.250E-01,-1.050E+00, 2.550E-01, 1.090E+01,-1.010E+01/    
      DATA ((CDO(IP,IS,5,2),IS=1,6),IP=1,3)/    
     1 8.790E-01, 0.000E+00, 4.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,    
     2-9.710E-01,-1.160E+00, 1.230E+00,-5.640E+00,-7.540E+00,-5.960E-01,    
     3 4.340E-01, 4.760E-01,-2.540E-01,-8.170E-01, 5.500E+00, 1.260E-01/    
    
C...The following data lines are coefficients needed in the 
C...Owens pion structure function parametrizations, see below.  
C...Expansion coefficients for up and down valence quark distributions. 
      DATA ((COW(IP,IS,1,1),IS=1,5),IP=1,3)/    
     1  4.0000E-01,  7.0000E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00, 
     2 -6.2120E-02,  6.4780E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00, 
     3 -7.1090E-03,  1.3350E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/ 
      DATA ((COW(IP,IS,1,2),IS=1,5),IP=1,3)/    
     1  4.0000E-01,  6.2800E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00, 
     2 -5.9090E-02,  6.4360E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00, 
     3 -6.5240E-03,  1.4510E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/ 
C...Expansion coefficients for gluon distribution.  
      DATA ((COW(IP,IS,2,1),IS=1,5),IP=1,3)/    
     1  8.8800E-01,  0.0000E+00,  3.1100E+00,  6.0000E+00,  0.0000E+00, 
     2 -1.8020E+00, -1.5760E+00, -1.3170E-01,  2.8010E+00, -1.7280E+01, 
     3  1.8120E+00,  1.2000E+00,  5.0680E-01, -1.2160E+01,  2.0490E+01/ 
      DATA ((COW(IP,IS,2,2),IS=1,5),IP=1,3)/    
     1  7.9400E-01,  0.0000E+00,  2.8900E+00,  6.0000E+00,  0.0000E+00, 
     2 -9.1440E-01, -1.2370E+00,  5.9660E-01, -3.6710E+00, -8.1910E+00, 
     3  5.9660E-01,  6.5820E-01, -2.5500E-01, -2.3040E+00,  7.7580E+00/ 
C...Expansion coefficients for (up+down+strange) quark sea distribution.    
      DATA ((COW(IP,IS,3,1),IS=1,5),IP=1,3)/    
     1  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00, 
     2 -2.4280E-01, -2.1200E-01,  8.6730E-01,  1.2660E+00,  2.3820E+00, 
     3  1.3860E-01,  3.6710E-03,  4.7470E-02, -2.2150E+00,  3.4820E-01/ 
      DATA ((COW(IP,IS,3,2),IS=1,5),IP=1,3)/    
     1  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00, 
     2 -1.4170E-01, -1.6970E-01, -2.4740E+00, -2.5340E+00,  5.6210E-01, 
     3 -1.7400E-01, -9.6230E-02,  1.5750E+00,  1.3780E+00, -2.7010E-01/ 
C...Expansion coefficients for charm quark sea distribution.    
      DATA ((COW(IP,IS,4,1),IS=1,5),IP=1,3)/    
     1  0.0000E+00, -2.2120E-02,  2.8940E+00,  0.0000E+00,  0.0000E+00, 
     2  7.9280E-02, -3.7850E-01,  9.4330E+00,  5.2480E+00,  8.3880E+00, 
     3 -6.1340E-02, -1.0880E-01, -1.0852E+01, -7.1870E+00, -1.1610E+01/ 
      DATA ((COW(IP,IS,4,2),IS=1,5),IP=1,3)/    
     1  0.0000E+00, -8.8200E-02,  1.9240E+00,  0.0000E+00,  0.0000E+00, 
     2  6.2290E-02, -2.8920E-01,  2.4240E-01, -4.4630E+00, -8.3670E-01, 
     3 -4.0990E-02, -1.0820E-01,  2.0360E+00,  5.2090E+00, -4.8400E-02/ 
    
C...Euler's beta function, requires ordinary Gamma function 
      EULBET(X,Y)=PYGAMM(X)*PYGAMM(Y)/PYGAMM(X+Y)   
    
C...Reset structure functions, check x and hadron flavour.  
      ALAM=0.   
      DO 100 KFL=-6,6   
  100 XPQ(KFL)=0.   
      IF(X.LT.0..OR.X.GT.1.) THEN   
        WRITE(MSTU(11),1000) X  
        RETURN  
      ENDIF 
      KFA=IABS(KF)  
      IF(KFA.NE.211.AND.KFA.NE.2212.AND.KFA.NE.2112) THEN   
        WRITE(MSTU(11),1100) KF 
        RETURN  
      ENDIF 
    
C...Call user-supplied structure function. Select proton/neutron/pion.  
      IF(MSTP(51).EQ.0.OR.MSTP(52).GE.2) THEN   
        KFE=KFA 
        IF(KFA.EQ.2112) KFE=2212    
        CALL PYSTFE(KFE,X,Q2,XPQ)   
        GOTO 230    
      ENDIF 
      IF(KFA.EQ.211) GOTO 200   
    
      IF(MSTP(51).EQ.1.OR.MSTP(51).EQ.2) THEN   
C...Proton structure functions from Eichten, Hinchliffe, Lane, Quigg.   
C...Allowed variable range: 5 GeV2 < Q2 < 1E8 GeV2; 1E-4 < x < 1    
    
C...Determine set, Lamdba and x and t expansion variables.  
        NSET=MSTP(51)   
        IF(NSET.EQ.1) ALAM=0.2  
        IF(NSET.EQ.2) ALAM=0.29 
        TMIN=LOG(5./ALAM**2)    
        TMAX=LOG(1E8/ALAM**2)   
        IF(MSTP(52).EQ.0) THEN  
          T=TMIN    
        ELSE    
          T=LOG(Q2/ALAM**2) 
        ENDIF   
        VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))    
        NX=1    
        IF(X.LE.0.1) NX=2   
        IF(NX.EQ.1) VX=(2.*X-1.1)/0.9   
        IF(NX.EQ.2) VX=MAX(-1.,(2.*LOG(X)+11.51293)/6.90776)    
        CXS=1.  
        IF(X.LT.1E-4.AND.ABS(PARP(51)-1.).GT.0.01) CXS= 
     &  (1E-4/X)**(PARP(51)-1.) 
    
C...Chebyshev polynomials for x and t expansion.    
        TX(1)=1.    
        TX(2)=VX    
        TX(3)=2.*VX**2-1.   
        TX(4)=4.*VX**3-3.*VX    
        TX(5)=8.*VX**4-8.*VX**2+1.  
        TX(6)=16.*VX**5-20.*VX**3+5.*VX 
        TT(1)=1.    
        TT(2)=VT    
        TT(3)=2.*VT**2-1.   
        TT(4)=4.*VT**3-3.*VT    
        TT(5)=8.*VT**4-8.*VT**2+1.  
        TT(6)=16.*VT**5-20.*VT**3+5.*VT 
    
C...Calculate structure functions.  
        DO 120 KFL=1,6  
        XQSUM=0.    
        DO 110 IT=1,6   
        DO 110 IX=1,6   
  110   XQSUM=XQSUM+CEHLQ(IX,IT,NX,KFL,NSET)*TX(IX)*TT(IT)  
  120   XQ(KFL)=XQSUM*(1.-X)**NEHLQ(KFL,NSET)*CXS   
    
C...Put into output array.  
        XPQ(0)=XQ(4)    
        XPQ(1)=XQ(2)+XQ(3)  
        XPQ(2)=XQ(1)+XQ(3)  
        XPQ(3)=XQ(5)    
        XPQ(4)=XQ(6)    
        XPQ(-1)=XQ(3)   
        XPQ(-2)=XQ(3)   
        XPQ(-3)=XQ(5)   
        XPQ(-4)=XQ(6)   
    
C...Special expansion for bottom (threshold effects).   
        IF(MSTP(54).GE.5) THEN  
          IF(NSET.EQ.1) TMIN=8.1905 
          IF(NSET.EQ.2) TMIN=7.4474 
          IF(T.LE.TMIN) GOTO 140    
          VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))  
          TT(1)=1.  
          TT(2)=VT  
          TT(3)=2.*VT**2-1. 
          TT(4)=4.*VT**3-3.*VT  
          TT(5)=8.*VT**4-8.*VT**2+1.    
          TT(6)=16.*VT**5-20.*VT**3+5.*VT   
          XQSUM=0.  
          DO 130 IT=1,6 
          DO 130 IX=1,6 
  130     XQSUM=XQSUM+CEHLQ(IX,IT,NX,7,NSET)*TX(IX)*TT(IT)  
          XPQ(5)=XQSUM*(1.-X)**NEHLQ(7,NSET)    
          XPQ(-5)=XPQ(5)    
  140     CONTINUE  
        ENDIF   
    
C...Special expansion for top (threshold effects).  
        IF(MSTP(54).GE.6) THEN  
          IF(NSET.EQ.1) TMIN=11.5528    
          IF(NSET.EQ.2) TMIN=10.8097    
          TMIN=TMIN+2.*LOG(PMAS(6,1)/30.)   
          TMAX=TMAX+2.*LOG(PMAS(6,1)/30.)   
          IF(T.LE.TMIN) GOTO 160    
          VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))  
          TT(1)=1.  
          TT(2)=VT  
          TT(3)=2.*VT**2-1. 
          TT(4)=4.*VT**3-3.*VT  
          TT(5)=8.*VT**4-8.*VT**2+1.    
          TT(6)=16.*VT**5-20.*VT**3+5.*VT   
          XQSUM=0.  
          DO 150 IT=1,6 
          DO 150 IX=1,6 
  150     XQSUM=XQSUM+CEHLQ(IX,IT,NX,8,NSET)*TX(IX)*TT(IT)  
          XPQ(6)=XQSUM*(1.-X)**NEHLQ(8,NSET)    
          XPQ(-6)=XPQ(6)    
  160     CONTINUE  
        ENDIF   
    
      ELSEIF(MSTP(51).EQ.3.OR.MSTP(51).EQ.4) THEN   
C...Proton structure functions from Duke, Owens.    
C...Allowed variable range: 4 GeV2 < Q2 < approx 1E6 GeV2.  
    
C...Determine set, Lambda and s expansion parameter.    
        NSET=MSTP(51)-2 
        IF(NSET.EQ.1) ALAM=0.2  
        IF(NSET.EQ.2) ALAM=0.4  
        IF(MSTP(52).LE.0) THEN  
          SD=0. 
        ELSE    
          SD=LOG(LOG(MAX(Q2,4.)/ALAM**2)/LOG(4./ALAM**2))   
        ENDIF   
    
C...Calculate structure functions.  
        DO 180 KFL=1,5  
        DO 170 IS=1,6   
  170   TS(IS)=CDO(1,IS,KFL,NSET)+CDO(2,IS,KFL,NSET)*SD+    
     &  CDO(3,IS,KFL,NSET)*SD**2    
        IF(KFL.LE.2) THEN   
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)*(1.+TS(3)*X)/(EULBET(TS(1),    
     &    TS(2)+1.)*(1.+TS(3)*TS(1)/(TS(1)+TS(2)+1.)))  
        ELSE    
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2+  
     &    TS(6)*X**3)   
        ENDIF   
  180   CONTINUE    
    
C...Put into output arrays. 
        XPQ(0)=XQ(5)    
        XPQ(1)=XQ(2)+XQ(3)/6.   
        XPQ(2)=3.*XQ(1)-XQ(2)+XQ(3)/6.  
        XPQ(3)=XQ(3)/6. 
        XPQ(4)=XQ(4)    
        XPQ(-1)=XQ(3)/6.    
        XPQ(-2)=XQ(3)/6.    
        XPQ(-3)=XQ(3)/6.    
        XPQ(-4)=XQ(4)   
    
C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli. 
C...These are accessed via PYSTFE since the files needed may not always 
C...available.  
      ELSEIF(MSTP(51).GE.11.AND.MSTP(51).LE.13) THEN    
        CALL PYSTFE(2212,X,Q2,XPQ)  
    
C...Unknown proton parametrization. 
      ELSE  
        WRITE(MSTU(11),1200) MSTP(51)   
      ENDIF 
      GOTO 230  
    
  200 IF((MSTP(51).GE.1.AND.MSTP(51).LE.4).OR.  
     &(MSTP(51).GE.11.AND.MSTP(51).LE.13)) THEN 
C...Pion structure functions from Owens.    
C...Allowed variable range: 4 GeV2 < Q2 < approx 2000 GeV2. 
    
C...Determine set, Lambda and s expansion variable. 
        NSET=1  
        IF(MSTP(51).EQ.2.OR.MSTP(51).EQ.4.OR.MSTP(51).EQ.13) NSET=2 
        IF(NSET.EQ.1) ALAM=0.2  
        IF(NSET.EQ.2) ALAM=0.4  
        IF(MSTP(52).LE.0) THEN  
          SD=0. 
        ELSE    
          SD=LOG(LOG(MAX(Q2,4.)/ALAM**2)/LOG(4./ALAM**2))   
        ENDIF   
    
C...Calculate structure functions.  
        DO 220 KFL=1,4  
        DO 210 IS=1,5   
  210   TS(IS)=COW(1,IS,KFL,NSET)+COW(2,IS,KFL,NSET)*SD+    
     &  COW(3,IS,KFL,NSET)*SD**2    
        IF(KFL.EQ.1) THEN   
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/EULBET(TS(1),TS(2)+1.) 
        ELSE    
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)  
        ENDIF   
  220   CONTINUE    
    
C...Put into output arrays. 
        XPQ(0)=XQ(2)    
        XPQ(1)=XQ(3)/6. 
        XPQ(2)=XQ(1)+XQ(3)/6.   
        XPQ(3)=XQ(3)/6. 
        XPQ(4)=XQ(4)    
        XPQ(-1)=XQ(1)+XQ(3)/6.  
        XPQ(-2)=XQ(3)/6.    
        XPQ(-3)=XQ(3)/6.    
        XPQ(-4)=XQ(4)   
    
C...Unknown pion parametrization.   
      ELSE  
        WRITE(MSTU(11),1200) MSTP(51)   
      ENDIF 
    
C...Isospin conjugation for neutron, charge conjugation for antipart.   
  230 IF(KFA.EQ.2112) THEN  
        XPS=XPQ(1)  
        XPQ(1)=XPQ(2)   
        XPQ(2)=XPS  
        XPS=XPQ(-1) 
        XPQ(-1)=XPQ(-2) 
        XPQ(-2)=XPS 
      ENDIF 
      IF(KF.LT.0) THEN  
        DO 240 KFL=1,4  
        XPS=XPQ(KFL)    
        XPQ(KFL)=XPQ(-KFL)  
  240   XPQ(-KFL)=XPS   
      ENDIF 
    
C...Check positivity and reset above maximum allowed flavour.   
      DO 250 KFL=-6,6   
      XPQ(KFL)=MAX(0.,XPQ(KFL)) 
  250 IF(IABS(KFL).GT.MSTP(54)) XPQ(KFL)=0. 

C...consider nuclear effect on the structure function
      	IF((JBT.NE.1.AND.JBT.NE.2).OR.IHPR2(6).EQ.0) GO TO 400
      	ATNM=IHNT2(2*JBT-1)
      	IF(ATNM.LE.1.0) GO TO 400
      	IF(JBT.EQ.1) THEN
	  BBR2=(YP(1,IHNT2(11))**2+YP(2,IHNT2(11))**2)/1.44/ATNM**0.66666
      	ELSEIF(JBT.EQ.2) THEN
	  BBR2=(YT(1,IHNT2(12))**2+YT(2,IHNT2(12))**2)/1.44/ATNM**0.66666
      	ENDIF
      	BBR2=MIN(1.0,BBR2)
	ABX=(ATNM**0.33333333-1.0)
      	APX=HIPR1(6)*4.0/3.0*ABX*SQRT(1.0-BBR2)
      	AAX=1.192*ALOG(ATNM)**0.1666666
      	RRX=AAX*(X**3-1.2*X**2+0.21*X)+1.0
     &		-(APX-1.079*ABX*SQRT(X)/ALOG(ATNM+1.0))*EXP(-X**2.0/0.01)
      	DO 300 KFL=-6,6
		XPQ(KFL)=XPQ(KFL)*RRX
300   	CONTINUE

C			********consider the nuclear effect on the structure
C				fucntion which also depends on the impact
C				parameter of the nuclear reaction

400   	CONTINUE    
C...Formats for error printouts.    
 1000 FORMAT(' Error: x value outside physical range, x =',1P,E12.3)    
 1100 FORMAT(' Error: illegal particle code for structure function,',   
     &' KF =',I5)   
 1200 FORMAT(' Error: bad value of parameter MSTP(51) in PYSTFU,',  
     &' MSTP(51) =',I5) 
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYSPLI(KF,KFLIN,KFLCH,KFLSP)   
    
C...In case of a hadron remnant which is more complicated than just a   
C...quark or a diquark, split it into two (partons or hadron + parton). 
      DIMENSION KFL(3)  
    
C...Preliminaries. Parton composition.  
      KFA=IABS(KF)  
      KFS=ISIGN(1,KF)   
      KFL(1)=MOD(KFA/1000,10)   
      KFL(2)=MOD(KFA/100,10)    
      KFL(3)=MOD(KFA/10,10) 
      KFLR=KFLIN*KFS    
      KFLCH=0   
    
C...Subdivide meson.    
      IF(KFL(1).EQ.0) THEN  
        KFL(2)=KFL(2)*(-1)**KFL(2)  
        KFL(3)=-KFL(3)*(-1)**IABS(KFL(2))   
        IF(KFLR.EQ.KFL(2)) THEN 
          KFLSP=KFL(3)  
        ELSEIF(KFLR.EQ.KFL(3)) THEN 
          KFLSP=KFL(2)  
        ELSEIF(IABS(KFLR).EQ.21.AND.RLU(0).GT.0.5) THEN 
          KFLSP=KFL(2)  
          KFLCH=KFL(3)  
        ELSEIF(IABS(KFLR).EQ.21) THEN   
          KFLSP=KFL(3)  
          KFLCH=KFL(2)  
        ELSEIF(KFLR*KFL(2).GT.0) THEN   
          CALL LUKFDI(-KFLR,KFL(2),KFDUMP,KFLCH)    
          KFLSP=KFL(3)  
        ELSE    
          CALL LUKFDI(-KFLR,KFL(3),KFDUMP,KFLCH)    
          KFLSP=KFL(2)  
        ENDIF   
    
C...Subdivide baryon.   
      ELSE  
        NAGR=0  
        DO 100 J=1,3    
  100   IF(KFLR.EQ.KFL(J)) NAGR=NAGR+1  
        IF(NAGR.GE.1) THEN  
          RAGR=0.00001+(NAGR-0.00002)*RLU(0)    
          IAGR=0    
          DO 110 J=1,3  
          IF(KFLR.EQ.KFL(J)) RAGR=RAGR-1.   
  110     IF(IAGR.EQ.0.AND.RAGR.LE.0.) IAGR=J   
        ELSE    
          IAGR=1.00001+2.99998*RLU(0)   
        ENDIF   
        ID1=1   
        IF(IAGR.EQ.1) ID1=2 
        IF(IAGR.EQ.1.AND.KFL(3).GT.KFL(2)) ID1=3    
        ID2=6-IAGR-ID1  
        KSP=3   
        IF(MOD(KFA,10).EQ.2.AND.KFL(1).EQ.KFL(2)) THEN  
          IF(IAGR.NE.3.AND.RLU(0).GT.0.25) KSP=1    
        ELSEIF(MOD(KFA,10).EQ.2.AND.KFL(2).GE.KFL(3)) THEN  
          IF(IAGR.NE.1.AND.RLU(0).GT.0.25) KSP=1    
        ELSEIF(MOD(KFA,10).EQ.2) THEN   
          IF(IAGR.EQ.1) KSP=1   
          IF(IAGR.NE.1.AND.RLU(0).GT.0.75) KSP=1    
        ENDIF   
        KFLSP=1000*KFL(ID1)+100*KFL(ID2)+KSP    
        IF(KFLIN.EQ.21) THEN    
          KFLCH=KFL(IAGR)   
        ELSEIF(NAGR.EQ.0.AND.KFLR.GT.0) THEN    
          CALL LUKFDI(-KFLR,KFL(IAGR),KFDUMP,KFLCH) 
        ELSEIF(NAGR.EQ.0) THEN  
          CALL LUKFDI(10000+KFLSP,-KFLR,KFDUMP,KFLCH)   
          KFLSP=KFL(IAGR)   
        ENDIF   
      ENDIF 
    
C...Add on correct sign for result. 
      KFLCH=KFLCH*KFS   
      KFLSP=KFLSP*KFS   
    
      RETURN    
      END   
    
C*********************************************************************  
    
      FUNCTION PYGAMM(X)    
    
C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;    
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions 
C...(Dover, 1965) 6.1.36.   
      DIMENSION B(8)    
      DATA B/-0.577191652,0.988205891,-0.897056937,0.918206857, 
     &-0.756704078,0.482199394,-0.193527818,0.035868343/    
    
      NX=INT(X) 
      DX=X-NX   
    
      PYGAMM=1. 
      DO 100 I=1,8  
  100 PYGAMM=PYGAMM+B(I)*DX**I  
      IF(X.LT.1.) THEN  
        PYGAMM=PYGAMM/X 
      ELSE  
        DO 110 IX=1,NX-1    
  110   PYGAMM=(X-IX)*PYGAMM    
      ENDIF 
    
      RETURN    
      END   
    
C***********************************************************************    
    
      FUNCTION PYW1AU(EPS,IREIM)    
    
C...Calculates real and imaginary parts of the auxiliary function W1;   
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,   
C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987    
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
    
      ASINH(X)=LOG(X+SQRT(X**2+1.)) 
      ACOSH(X)=LOG(X+SQRT(X**2-1.)) 
    
      IF(EPS.LT.0.) THEN    
        W1RE=2.*SQRT(1.-EPS)*ASINH(SQRT(-1./EPS))   
        W1IM=0. 
      ELSEIF(EPS.LT.1.) THEN    
        W1RE=2.*SQRT(1.-EPS)*ACOSH(SQRT(1./EPS))    
        W1IM=-PARU(1)*SQRT(1.-EPS)  
      ELSE  
        W1RE=2.*SQRT(EPS-1.)*ASIN(SQRT(1./EPS)) 
        W1IM=0. 
      ENDIF 
    
      IF(IREIM.EQ.1) PYW1AU=W1RE    
      IF(IREIM.EQ.2) PYW1AU=W1IM    
    
      RETURN    
      END   
    
C***********************************************************************    
    
      FUNCTION PYW2AU(EPS,IREIM)    
    
C...Calculates real and imaginary parts of the auxiliary function W2;   
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,   
C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987    
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
    
      ASINH(X)=LOG(X+SQRT(X**2+1.)) 
      ACOSH(X)=LOG(X+SQRT(X**2-1.)) 
    
      IF(EPS.LT.0.) THEN    
        W2RE=4.*(ASINH(SQRT(-1./EPS)))**2   
        W2IM=0. 
      ELSEIF(EPS.LT.1.) THEN    
        W2RE=4.*(ACOSH(SQRT(1./EPS)))**2-PARU(1)**2 
        W2IM=-4.*PARU(1)*ACOSH(SQRT(1./EPS))    
      ELSE  
        W2RE=-4.*(ASIN(SQRT(1./EPS)))**2    
        W2IM=0. 
      ENDIF 
    
      IF(IREIM.EQ.1) PYW2AU=W2RE    
      IF(IREIM.EQ.2) PYW2AU=W2IM    
    
      RETURN    
      END   
    
C***********************************************************************    
    
      FUNCTION PYI3AU(BE,EPS,IREIM) 
    
C...Calculates real and imaginary parts of the auxiliary function I3;   
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,   
C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987    
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
    
      IF(EPS.LT.1.) GA=0.5*(1.+SQRT(1.-EPS))    
    
      IF(EPS.LT.0.) THEN    
        F3RE=PYSPEN((GA-1.)/(GA+BE-1.),0.,1)-PYSPEN(GA/(GA+BE-1.),0.,1)+    
     &  PYSPEN((BE-GA)/BE,0.,1)-PYSPEN((BE-GA)/(BE-1.),0.,1)+   
     &  (LOG(BE)**2-LOG(BE-1.)**2)/2.+LOG(GA)*LOG((GA+BE-1.)/BE)+   
     &  LOG(GA-1.)*LOG((BE-1.)/(GA+BE-1.))  
        F3IM=0. 
      ELSEIF(EPS.LT.1.) THEN    
        F3RE=PYSPEN((GA-1.)/(GA+BE-1.),0.,1)-PYSPEN(GA/(GA+BE-1.),0.,1)+    
     &  PYSPEN(GA/(GA-BE),0.,1)-PYSPEN((GA-1.)/(GA-BE),0.,1)+   
     &  LOG(GA/(1.-GA))*LOG((GA+BE-1.)/(BE-GA)) 
        F3IM=-PARU(1)*LOG((GA+BE-1.)/(BE-GA))   
      ELSE  
        RSQ=EPS/(EPS-1.+(2.*BE-1.)**2)  
        RCTHE=RSQ*(1.-2.*BE/EPS)    
        RSTHE=SQRT(RSQ-RCTHE**2)    
        RCPHI=RSQ*(1.+2.*(BE-1.)/EPS)   
        RSPHI=SQRT(RSQ-RCPHI**2)    
        R=SQRT(RSQ) 
        THE=ACOS(RCTHE/R)   
        PHI=ACOS(RCPHI/R)   
        F3RE=PYSPEN(RCTHE,RSTHE,1)+PYSPEN(RCTHE,-RSTHE,1)-  
     &  PYSPEN(RCPHI,RSPHI,1)-PYSPEN(RCPHI,-RSPHI,1)+   
     &  (PHI-THE)*(PHI+THE-PARU(1)) 
        F3IM=PYSPEN(RCTHE,RSTHE,2)+PYSPEN(RCTHE,-RSTHE,2)-  
     &  PYSPEN(RCPHI,RSPHI,2)-PYSPEN(RCPHI,-RSPHI,2)    
      ENDIF 
    
      IF(IREIM.EQ.1) PYI3AU=2./(2.*BE-1.)*F3RE  
      IF(IREIM.EQ.2) PYI3AU=2./(2.*BE-1.)*F3IM  
    
      RETURN    
      END   
    
C***********************************************************************    
    
      FUNCTION PYSPEN(XREIN,XIMIN,IREIM)    
    
C...Calculates real and imaginary part of Spence function; see  
C...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.    
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      DIMENSION B(0:14) 
    
      DATA B/   
     & 1.000000E+00,        -5.000000E-01,         1.666667E-01,    
     & 0.000000E+00,        -3.333333E-02,         0.000000E+00,    
     & 2.380952E-02,         0.000000E+00,        -3.333333E-02,    
     & 0.000000E+00,         7.575757E-02,         0.000000E+00,    
     &-2.531135E-01,         0.000000E+00,         1.166667E+00/    
    
      XRE=XREIN 
      XIM=XIMIN 
      IF(ABS(1.-XRE).LT.1.E-6.AND.ABS(XIM).LT.1.E-6) THEN   
        IF(IREIM.EQ.1) PYSPEN=PARU(1)**2/6. 
        IF(IREIM.EQ.2) PYSPEN=0.    
        RETURN  
      ENDIF 
    
      XMOD=SQRT(XRE**2+XIM**2)  
      IF(XMOD.LT.1.E-6) THEN    
        IF(IREIM.EQ.1) PYSPEN=0.    
        IF(IREIM.EQ.2) PYSPEN=0.    
        RETURN  
      ENDIF 
    
      XARG=SIGN(ACOS(XRE/XMOD),XIM) 
      SP0RE=0.  
      SP0IM=0.  
      SGN=1.    
      IF(XMOD.GT.1.) THEN   
        ALGXRE=LOG(XMOD)    
        ALGXIM=XARG-SIGN(PARU(1),XARG)  
        SP0RE=-PARU(1)**2/6.-(ALGXRE**2-ALGXIM**2)/2.   
        SP0IM=-ALGXRE*ALGXIM    
        SGN=-1. 
        XMOD=1./XMOD    
        XARG=-XARG  
        XRE=XMOD*COS(XARG)  
        XIM=XMOD*SIN(XARG)  
      ENDIF 
      IF(XRE.GT.0.5) THEN   
        ALGXRE=LOG(XMOD)    
        ALGXIM=XARG 
        XRE=1.-XRE  
        XIM=-XIM    
        XMOD=SQRT(XRE**2+XIM**2)    
        XARG=SIGN(ACOS(XRE/XMOD),XIM)   
        ALGYRE=LOG(XMOD)    
        ALGYIM=XARG 
        SP0RE=SP0RE+SGN*(PARU(1)**2/6.-(ALGXRE*ALGYRE-ALGXIM*ALGYIM))   
        SP0IM=SP0IM-SGN*(ALGXRE*ALGYIM+ALGXIM*ALGYRE)   
        SGN=-SGN    
      ENDIF 
    
      XRE=1.-XRE    
      XIM=-XIM  
      XMOD=SQRT(XRE**2+XIM**2)  
      XARG=SIGN(ACOS(XRE/XMOD),XIM) 
      ZRE=-LOG(XMOD)    
      ZIM=-XARG 
    
      SPRE=0.   
      SPIM=0.   
      SAVERE=1. 
      SAVEIM=0. 
      DO 100 I=0,14 
      TERMRE=(SAVERE*ZRE-SAVEIM*ZIM)/FLOAT(I+1) 
      TERMIM=(SAVERE*ZIM+SAVEIM*ZRE)/FLOAT(I+1) 
      SAVERE=TERMRE 
      SAVEIM=TERMIM 
      SPRE=SPRE+B(I)*TERMRE 
  100 SPIM=SPIM+B(I)*TERMIM 
    
      IF(IREIM.EQ.1) PYSPEN=SP0RE+SGN*SPRE  
      IF(IREIM.EQ.2) PYSPEN=SP0IM+SGN*SPIM  
    
      RETURN    
      END   
    
*********************************************************************** 
    
      SUBROUTINE PYTEST(MTEST)  
    
C...Purpose: to provide a simple program (disguised as a subroutine) to 
C...run at installation as a check that the program works as intended.  
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
    
C...Common initial values. Loop over initiating conditions. 
      MSTP(122)=1   
      IF(MTEST.LE.0) MSTP(122)=0    
      MDCY(LUCOMP(111),1)=0 
      NERR=0    
      DO 130 IPROC=1,7  
    
C...Reset process type, kinematics cuts, and the flags used.    
      MSEL=0    
      DO 100 ISUB=1,200 
  100 MSUB(ISUB)=0  
      CKIN(1)=2.    
      CKIN(3)=0.    
      MSTP(2)=1 
      MSTP(33)=0    
      MSTP(81)=1    
      MSTP(82)=1    
      MSTP(111)=1   
      MSTP(131)=0   
      MSTP(133)=0   
      PARP(131)=0.01    
    
C...Prompt photon production at fixed target.   
      IF(IPROC.EQ.1) THEN   
        PZSUM=300.  
        PESUM=SQRT(PZSUM**2+ULMASS(211)**2)+ULMASS(2212)    
        PQSUM=2.    
        MSEL=10 
        CKIN(3)=5.  
        CALL PYINIT('FIXT','pi+','p',PZSUM) 
    
C...QCD processes at ISR energies.  
      ELSEIF(IPROC.EQ.2) THEN   
        PESUM=63.   
        PZSUM=0.    
        PQSUM=2.    
        MSEL=1  
        CKIN(3)=5.  
        CALL PYINIT('CMS','p','p',PESUM)    
    
C...W production + multiple interactions at CERN Collider.  
      ELSEIF(IPROC.EQ.3) THEN   
        PESUM=630.  
        PZSUM=0.    
        PQSUM=0.    
        MSEL=12 
        CKIN(1)=20. 
        MSTP(82)=4  
        MSTP(2)=2   
        MSTP(33)=3  
        CALL PYINIT('CMS','p','pbar',PESUM) 
    
C...W/Z gauge boson pairs + overlayed events at the Tevatron.   
      ELSEIF(IPROC.EQ.4) THEN   
        PESUM=1800. 
        PZSUM=0.    
        PQSUM=0.    
        MSUB(22)=1  
        MSUB(23)=1  
        MSUB(25)=1  
        CKIN(1)=200.    
        MSTP(111)=0 
        MSTP(131)=1 
        MSTP(133)=2 
        PARP(131)=0.04  
        CALL PYINIT('CMS','p','pbar',PESUM) 
    
C...Higgs production at LHC.    
      ELSEIF(IPROC.EQ.5) THEN   
        PESUM=17000.    
        PZSUM=0.    
        PQSUM=0.    
        MSEL=16 
        PMAS(25,1)=300. 
        CKIN(1)=200.    
        MSTP(81)=0  
        MSTP(111)=0 
        CALL PYINIT('CMS','p','pbar',PESUM) 
    
C...Z' production at SSC.   
      ELSEIF(IPROC.EQ.6) THEN   
        PESUM=40000.    
        PZSUM=0.    
        PQSUM=0.    
        MSEL=21 
        PMAS(32,1)=600. 
        CKIN(1)=400.    
        MSTP(81)=0  
        MSTP(111)=0 
        CALL PYINIT('CMS','p','pbar',PESUM) 
    
C...W pair production at 1 TeV e+e- collider.   
      ELSEIF(IPROC.EQ.7) THEN   
        PESUM=1000. 
        PZSUM=0.    
        PQSUM=0.    
        MSUB(25)=1  
        CALL PYINIT('CMS','e+','e-',PESUM)  
      ENDIF 
    
C...Generate 20 events of each required type.   
      DO 120 IEV=1,20   
      CALL PYTHIA   
      PESUMM=PESUM  
      IF(IPROC.EQ.4) PESUMM=MSTI(41)*PESUM  
    
C...Check conservation of energy/momentum/flavour.  
      MERR=0    
      DEVE=ABS(PLU(0,4)-PESUMM)+ABS(PLU(0,3)-PZSUM) 
      DEVT=ABS(PLU(0,1))+ABS(PLU(0,2))  
      DEVQ=ABS(PLU(0,6)-PQSUM)  
      IF(DEVE.GT.1E-3*PESUM.OR.DEVT.GT.MAX(0.01,1E-5*PESUM).OR. 
     &DEVQ.GT.0.1) MERR=1   
      IF(MERR.NE.0) WRITE(MSTU(11),1000) IPROC,IEV  
    
C...Check that all KF codes are known ones, and that partons/particles  
C...satisfy energy-momentum-mass relation.  
      DO 110 I=1,N  
      IF(K(I,1).GT.20) GOTO 110 
      IF(LUCOMP(K(I,2)).EQ.0) THEN  
        WRITE(MSTU(11),1100) I  
        MERR=MERR+1 
      ENDIF 
      PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2* 
     &SIGN(1.,P(I,5))   
      IF(ABS(PD).GT.MAX(0.1,0.002*P(I,4)**2,0.002*P(I,5)**2).OR.    
     &(P(I,5).GE.0..AND.P(I,4).LT.0.)) THEN 
        WRITE(MSTU(11),1200) I  
        MERR=MERR+1 
      ENDIF 
  110 CONTINUE  
    
C...Listing of erronoeus events, and first event of each type.  
      IF(MERR.GE.1) NERR=NERR+1 
      IF(NERR.GE.10) THEN   
        WRITE(MSTU(11),1300)    
        CALL LULIST(1)  
        STOP    
      ENDIF 
      IF(MTEST.GE.1.AND.(MERR.GE.1.OR.IEV.EQ.1)) THEN   
        IF(MERR.GE.1) WRITE(MSTU(11),1400)  
        CALL LULIST(1)  
      ENDIF 
  120 CONTINUE  
    
C...List statistics for each process type.  
      IF(MTEST.GE.1) CALL PYSTAT(1) 
  130 CONTINUE  
    
C...Summarize result of run.    
      IF(NERR.EQ.0) WRITE(MSTU(11),1500)    
      IF(NERR.GT.0) WRITE(MSTU(11),1600) NERR   
      RETURN    
    
C...Formats for information.    
 1000 FORMAT(/5X,'Energy/momentum/flavour nonconservation for process', 
     &I2,', event',I4)  
 1100 FORMAT(/5X,'Entry no.',I4,' in following event not known code')   
 1200 FORMAT(/5X,'Entry no.',I4,' in following event has faulty ',  
     &'kinematics') 
 1300 FORMAT(/5X,'This is the tenth error experienced! Something is ',  
     &'wrong.'/5X,'Execution will be stopped after listing of event.')  
 1400 FORMAT(5X,'Faulty event follows:')    
 1500 FORMAT(//5X,'End result of run: no errors detected.') 
 1600 FORMAT(//5X,'End result of run:',I2,' errors detected.'/  
     &5X,'This should not have happened!')  
      END   
    
C*********************************************************************  
    
      BLOCK DATA PYDATA 
    
C...Give sensible default values to all status codes and parameters.    
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      COMMON/PYINT6/PROC(0:200) 
      CHARACTER PROC*28 
      SAVE /PYINT6/ 
    
C...Default values for allowed processes and kinematics constraints.    
      DATA MSEL/1/  
      DATA MSUB/200*0/  
      DATA ((KFIN(I,J),J=-40,40),I=1,2)/40*1,0,80*1,0,40*1/ 
      DATA CKIN/    
     &   2.0, -1.0,  0.0, -1.0,  1.0,  1.0, -10.,  10., -10.,  10., 
     1  -10.,  10., -10.,  10., -10.,  10., -1.0,  1.0, -1.0,  1.0, 
     2   0.0,  1.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0,   0.,   0., 
     3   2.0, -1.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     4   160*0./    
    
C...Default values for main switches and parameters. Reset information. 
      DATA (MSTP(I),I=1,100)/   
     &     3,    1,    2,    0,    0,    0,    0,    0,    0,    0, 
     1     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     2     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     3     1,    2,    0,    0,    0,    2,    0,    0,    0,    0, 
     4     1,    0,    3,    7,    1,    0,    0,    0,    0,    0, 
     5     1,    1,   20,    6,    0,    0,    0,    0,    0,    0, 
     6     1,    2,    2,    2,    1,    0,    0,    0,    0,    0, 
     7     1,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     8     1,    1,  100,    0,    0,    0,    0,    0,    0,    0, 
     9     1,    4,    0,    0,    0,    0,    0,    0,    0,    0/ 
      DATA (MSTP(I),I=101,200)/ 
     &     1,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     1     1,    1,    1,    0,    0,    0,    0,    0,    0,    0, 
     2     0,    1,    2,    1,    1,   20,    0,    0,    0,    0, 
     3     0,    4,    0,    1,    0,    0,    0,    0,    0,    0, 
     4     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     5     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     6     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     7     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     8     5,    3, 1989,   11,   24,    0,    0,    0,    0,    0, 
     9     0,    0,    0,    0,    0,    0,    0,    0,    0,    0/ 
      DATA (PARP(I),I=1,100)/   
     &  0.25,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     1    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     2    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     3   1.5,  2.0, 0.075,  0.,  0.2,   0.,   0.,   0.,   0.,   0., 
     4    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     5   1.0, 2.26, 1.E4, 1.E-4,  0.,   0.,   0.,   0.,   0.,   0., 
     6  0.25,  1.0, 0.25,  1.0,  2.0, 1.E-3, 4.0,   0.,   0.,   0., 
     7   4.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     8   1.6, 1.85,  0.5,  0.2, 0.33, 0.66,  0.7,  0.5,   0.,   0., 
     9  0.44, 0.44,  2.0,  1.0,   0.,  3.0,  1.0, 0.75,   0.,   0./ 
      DATA (PARP(I),I=101,200)/ 
     & -0.02,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     1   2.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     2   0.4,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     3  0.01,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     4    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     5    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     6    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     7    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     8    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     9    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0./ 
      DATA MSTI/200*0/  
      DATA PARI/200*0./ 
      DATA MINT/400*0/  
      DATA VINT/400*0./ 
    
C...Constants for the generation of the various processes.  
      DATA (ISET(I),I=1,100)/   
     &    1,    1,    1,   -1,    3,   -1,   -1,    3,   -2,   -2,  
     1    2,    2,    2,    2,    2,    2,   -1,    2,    2,    2,  
     2   -1,    2,    2,    2,    2,    2,   -1,    2,    2,    2,  
     3    2,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  
     4   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  
     5   -1,   -1,    2,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  
     6   -1,   -1,   -1,   -1,   -1,   -1,   -1,    2,   -1,   -1,  
     7    4,    4,    4,   -1,   -1,    4,    4,   -1,   -1,   -2,  
     8    2,    2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     9    0,    0,    0,   -1,    0,    5,   -2,   -2,   -2,   -2/  
      DATA (ISET(I),I=101,200)/ 
     &   -1,    1,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     1    2,    2,    2,    2,   -1,   -1,   -1,   -2,   -2,   -2,  
     2   -1,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     3   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     4    1,    1,    1,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     5   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     6    2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     7   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     8   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,  
     9   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2/  
      DATA ((KFPR(I,J),J=1,2),I=1,50)/  
     &   23,    0,   24,    0,   25,    0,   24,    0,   25,    0,  
     &   24,    0,   23,    0,   25,    0,    0,    0,    0,    0,  
     1    0,    0,    0,    0,   21,   21,   21,   22,   21,   23,  
     1   21,   24,   21,   25,   22,   22,   22,   23,   22,   24,  
     2   22,   25,   23,   23,   23,   24,   23,   25,   24,   24,  
     2   24,   25,   25,   25,    0,   21,    0,   22,    0,   23,  
     3    0,   24,    0,   25,    0,   21,    0,   22,    0,   23,  
     3    0,   24,    0,   25,    0,   21,    0,   22,    0,   23,  
     4    0,   24,    0,   25,    0,   21,    0,   22,    0,   23,  
     4    0,   24,    0,   25,    0,   21,    0,   22,    0,   23/  
      DATA ((KFPR(I,J),J=1,2),I=51,100)/    
     5    0,   24,    0,   25,    0,    0,    0,    0,    0,    0,  
     5    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     6    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     6    0,    0,    0,    0,   21,   21,   24,   24,   22,   24,  
     7   23,   23,   24,   24,   23,   24,   23,   25,   22,   22,  
     7   23,   23,   24,   24,   24,   25,   25,   25,    0,    0,  
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/  
      DATA ((KFPR(I,J),J=1,2),I=101,150)/   
     &   23,    0,   25,    0,    0,    0,    0,    0,    0,    0,  
     &    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     1   21,   25,    0,   25,   21,   25,   22,   22,   22,   23,  
     1   23,   23,   24,   24,    0,    0,    0,    0,    0,    0,  
     2    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     2    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     3    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     3    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     4   32,    0,   37,    0,   40,    0,    0,    0,    0,    0,  
     4    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/  
      DATA ((KFPR(I,J),J=1,2),I=151,200)/   
     5    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     5    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     6    0,   37,    0,    0,    0,    0,    0,    0,    0,    0,  
     6    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     7    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     7    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/  
      DATA COEF/4000*0./    
      DATA (((ICOL(I,J,K),K=1,2),J=1,4),I=1,40)/    
     1 4,0,3,0,2,0,1,0,3,0,4,0,1,0,2,0,2,0,0,1,4,0,0,3,3,0,0,4,1,0,0,2, 
     2 3,0,0,4,1,4,3,2,4,0,0,3,4,2,1,3,2,0,4,1,4,0,2,3,4,0,3,4,2,0,1,2, 
     3 3,2,1,0,1,4,3,0,4,3,3,0,2,1,1,0,3,2,1,4,1,0,0,2,2,4,3,1,2,0,0,1, 
     4 3,2,1,4,1,4,3,2,4,2,1,3,4,2,1,3,3,4,4,3,1,2,2,1,2,0,3,1,2,0,0,0, 
     5 4,2,1,0,0,0,1,0,3,0,0,3,1,2,0,0,4,0,0,4,0,0,1,2,2,0,0,1,4,4,3,3, 
     6 2,2,1,1,4,4,3,3,3,3,4,4,1,1,2,2,3,2,1,3,1,2,0,0,4,2,1,4,0,0,1,2, 
     7 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
     8 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
     9 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/ 
    
C...Character constants: name of processes. 
      DATA PROC(0)/                    'All included subprocesses   '/  
      DATA (PROC(I),I=1,20)/    
     1'f + fb -> gamma*/Z0         ',  'f + fb'' -> W+/-             ', 
     2'f + fb -> H0                ',  'gamma + W+/- -> W+/-        ',  
     3'Z0 + Z0 -> H0               ',  'Z0 + W+/- -> W+/-           ',  
     4'                            ',  'W+ + W- -> H0               ',  
     5'                            ',  '                            ',  
     6'f + f'' -> f + f''            ','f + fb -> f'' + fb''          ',    
     7'f + fb -> g + g             ',  'f + fb -> g + gamma         ',  
     8'f + fb -> g + Z0            ',  'f + fb'' -> g + W+/-         ', 
     9'f + fb -> g + H0            ',  'f + fb -> gamma + gamma     ',  
     &'f + fb -> gamma + Z0        ',  'f + fb'' -> gamma + W+/-     '/ 
      DATA (PROC(I),I=21,40)/   
     1'f + fb -> gamma + H0        ',  'f + fb -> Z0 + Z0           ',  
     2'f + fb'' -> Z0 + W+/-        ', 'f + fb -> Z0 + H0           ',  
     3'f + fb -> W+ + W-           ',  'f + fb'' -> W+/- + H0        ', 
     4'f + fb -> H0 + H0           ',  'f + g -> f + g              ',  
     5'f + g -> f + gamma          ',  'f + g -> f + Z0             ',  
     6'f + g -> f'' + W+/-          ', 'f + g -> f + H0             ',  
     7'f + gamma -> f + g          ',  'f + gamma -> f + gamma      ',  
     8'f + gamma -> f + Z0         ',  'f + gamma -> f'' + W+/-      ', 
     9'f + gamma -> f + H0         ',  'f + Z0 -> f + g             ',  
     &'f + Z0 -> f + gamma         ',  'f + Z0 -> f + Z0            '/  
      DATA (PROC(I),I=41,60)/   
     1'f + Z0 -> f'' + W+/-         ', 'f + Z0 -> f + H0            ',  
     2'f + W+/- -> f'' + g          ', 'f + W+/- -> f'' + gamma      ', 
     3'f + W+/- -> f'' + Z0         ', 'f + W+/- -> f'' + W+/-       ', 
     4'f + W+/- -> f'' + H0         ', 'f + H0 -> f + g             ',  
     5'f + H0 -> f + gamma         ',  'f + H0 -> f + Z0            ',  
     6'f + H0 -> f'' + W+/-         ', 'f + H0 -> f + H0            ',  
     7'g + g -> f + fb             ',  'g + gamma -> f + fb         ',  
     8'g + Z0 -> f + fb            ',  'g + W+/- -> f + fb''         ', 
     9'g + H0 -> f + fb            ',  'gamma + gamma -> f + fb     ',  
     &'gamma + Z0 -> f + fb        ',  'gamma + W+/- -> f + fb''     '/ 
      DATA (PROC(I),I=61,80)/   
     1'gamma + H0 -> f + fb        ',  'Z0 + Z0 -> f + fb           ',  
     2'Z0 + W+/- -> f + fb''        ', 'Z0 + H0 -> f + fb           ',  
     3'W+ + W- -> f + fb           ',  'W+/- + H0 -> f + fb''        ', 
     4'H0 + H0 -> f + fb           ',  'g + g -> g + g              ',  
     5'gamma + gamma -> W+ + W-    ',  'gamma + W+/- -> gamma + W+/-',  
     6'Z0 + Z0 -> Z0 + Z0          ',  'Z0 + Z0 -> W+ + W-          ',  
     7'Z0 + W+/- -> Z0 + W+/-      ',  'Z0 + Z0 -> Z0 + H0          ',  
     8'W+ + W- -> gamma + gamma    ',  'W+ + W- -> Z0 + Z0          ',  
     9'W+/- + W+/- -> W+/- + W+/-  ',  'W+/- + H0 -> W+/- + H0      ',  
     &'H0 + H0 -> H0 + H0          ',  '                            '/  
      DATA (PROC(I),I=81,100)/  
     1'q + qb -> Q + QB, massive   ',  'g + g -> Q + QB, massive    ',  
     2'                            ',  '                            ',  
     3'                            ',  '                            ',  
     4'                            ',  '                            ',  
     5'                            ',  '                            ',  
     6'Elastic scattering          ',  'Single diffractive          ',  
     7'Double diffractive          ',  'Central diffractive         ',  
     8'Low-pT scattering           ',  'Semihard QCD 2 -> 2         ',  
     9'                            ',  '                            ',  
     &'                            ',  '                            '/  
      DATA (PROC(I),I=101,120)/ 
     1'g + g -> gamma*/Z0          ',  'g + g -> H0                 ',  
     2'                            ',  '                            ',  
     3'                            ',  '                            ',  
     4'                            ',  '                            ',  
     5'                            ',  '                            ',  
     6'f + fb -> g + H0            ',  'q + g -> q + H0             ',  
     7'g + g -> g + H0             ',  'g + g -> gamma + gamma      ',  
     8'g + g -> gamma + Z0         ',  'g + g -> Z0 + Z0            ',  
     9'g + g -> W+ + W-            ',  '                            ',  
     &'                            ',  '                            '/  
      DATA (PROC(I),I=121,140)/ 
     1'g + g -> f + fb + H0        ',  '                            ',  
     2'                            ',  '                            ',  
     3'                            ',  '                            ',  
     4'                            ',  '                            ',  
     5'                            ',  '                            ',  
     6'                            ',  '                            ',  
     7'                            ',  '                            ',  
     8'                            ',  '                            ',  
     9'                            ',  '                            ',  
     &'                            ',  '                            '/  
      DATA (PROC(I),I=141,160)/ 
     1'f + fb -> gamma*/Z0/Z''0     ', 'f + fb'' -> H+/-             ', 
     2'f + fb -> R                 ',  '                            ',  
     3'                            ',  '                            ',  
     4'                            ',  '                            ',  
     5'                            ',  '                            ',  
     6'                            ',  '                            ',  
     7'                            ',  '                            ',  
     8'                            ',  '                            ',  
     9'                            ',  '                            ',  
     &'                            ',  '                            '/  
      DATA (PROC(I),I=161,180)/ 
     1'f + g -> f'' + H+/-          ', '                            ',  
     2'                            ',  '                            ',  
     3'                            ',  '                            ',  
     4'                            ',  '                            ',  
     5'                            ',  '                            ',  
     6'                            ',  '                            ',  
     7'                            ',  '                            ',  
     8'                            ',  '                            ',  
     9'                            ',  '                            ',  
     &'                            ',  '                            '/  
      DATA (PROC(I),I=181,200)/     20*'                            '/  
    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYKCUT(MCUT)   
    
C...Dummy routine, which the user can replace in order to make cuts on  
C...the kinematics on the parton level before the matrix elements are   
C...evaluated and the event is generated. The cross-section estimates   
C...will automatically take these cuts into account, so the given   
C...values are for the allowed phase space region only. MCUT=0 means    
C...that the event has passed the cuts, MCUT=1 that it has failed.  
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
    
      MCUT=0    
    
      RETURN    
      END   
    
C*********************************************************************  
    
      SUBROUTINE PYSTFE(KF,X,Q2,XPQ)    
    
C...This is a dummy routine, where the user can introduce an interface  
C...to his own external structure function parametrization. 
C...Arguments in:   
C...KF : 2212 for p, 211 for pi+; isospin conjugation for n and charge  
C...    conjugation for pbar, nbar or pi- is performed by PYSTFU.   
C...X : x value.    
C...Q2 : Q^2 value. 
C...Arguments out:  
C...XPQ(-6:6) : x * f(x,Q2), with index according to KF code,   
C...    except that gluon is placed in 0. Thus XPQ(0) = xg, 
C...    XPQ(1) = xd, XPQ(-1) = xdbar, XPQ(2) = xu, XPQ(-2) = xubar, 
C...    XPQ(3) = xs, XPQ(-3) = xsbar, XPQ(4) = xc, XPQ(-4) = xcbar, 
C...    XPQ(5) = xb, XPQ(-5) = xbbar, XPQ(6) = xt, XPQ(-6) = xtbar. 
C...    
C...One such interface, to the Diemos, Ferroni, Longo, Martinelli   
C...proton structure functions, already comes with the package. What    
C...the user needs here is external files with the three routines   
C...FXG160, FXG260 and FXG360 of the authors above, plus the    
C...interpolation routine FINT, which is part of the CERN library   
C...KERNLIB package. To avoid problems with unresolved external 
C...references, the external calls are commented in the current 
C...version. To enable this option, remove the C* at the beginning  
C...of the relevant lines.  
C...    
C...Alternatively, the routine can be used as an interface to the   
C...structure function evolution program of Tung. This can be achieved  
C...by removing C* at the beginning of some of the lines below. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      DIMENSION XPQ(-6:6),XFDFLM(9) 
      CHARACTER CHDFLM(9)*5,HEADER*40   
      DATA CHDFLM/'UPVAL','DOVAL','GLUON','QBAR ','UBAR ','SBAR ',  
     &'CBAR ','BBAR ','TBAR '/  
      DATA HEADER/'Tung evolution package has been invoked'/    
      DATA INIT/0/  
    
C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli. 
C...Allowed variable range 10 GeV2 < Q2 < 1E8 GeV2, 5E-5 < x < .95. 
      IF(MSTP(51).GE.11.AND.MSTP(51).LE.13.AND.MSTP(52).LE.1) THEN  
        XDFLM=MAX(0.51E-4,X)    
        Q2DFLM=MAX(10.,MIN(1E8,Q2)) 
        IF(MSTP(52).EQ.0) Q2DFLM=10.    
        DO 100 J=1,9    
        IF(MSTP(52).EQ.1.AND.J.EQ.9) THEN   
          Q2DFLM=Q2DFLM*(40./PMAS(6,1))**2  
          Q2DFLM=MAX(10.,MIN(1E8,Q2))   
        ENDIF   
        XFDFLM(J)=0.    
C...Remove C* on following three lines to enable the DFLM options.  
C*      IF(MSTP(51).EQ.11) CALL FXG160(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))    
C*      IF(MSTP(51).EQ.12) CALL FXG260(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))    
C*      IF(MSTP(51).EQ.13) CALL FXG360(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))    
  100   CONTINUE    
        IF(X.LT.0.51E-4.AND.ABS(PARP(51)-1.).GT.0.01) THEN  
          CXS=(0.51E-4/X)**(PARP(51)-1.)    
          DO 110 J=1,7  
  110     XFDFLM(J)=XFDFLM(J)*CXS   
        ENDIF   
        XPQ(0)=XFDFLM(3)    
        XPQ(1)=XFDFLM(2)+XFDFLM(5)  
        XPQ(2)=XFDFLM(1)+XFDFLM(5)  
        XPQ(3)=XFDFLM(6)    
        XPQ(4)=XFDFLM(7)    
        XPQ(5)=XFDFLM(8)    
        XPQ(6)=XFDFLM(9)    
        XPQ(-1)=XFDFLM(5)   
        XPQ(-2)=XFDFLM(5)   
        XPQ(-3)=XFDFLM(6)   
        XPQ(-4)=XFDFLM(7)   
        XPQ(-5)=XFDFLM(8)   
        XPQ(-6)=XFDFLM(9)   
    
C...Proton structure function evolution from Wu-Ki Tung: parton 
C...distribution functions incorporating heavy quark mass effects.  
C...Allowed variable range: PARP(52) < Q < PARP(53); PARP(54) < x < 1.  
      ELSE  
        IF(INIT.EQ.0) THEN  
          I1=0  
          IF(MSTP(52).EQ.4) I1=1    
          IHDRN=1   
          NU=MSTP(53)   
          I2=MSTP(51)   
          IF(MSTP(51).GE.11) I2=MSTP(51)-3  
          I3=0  
          IF(MSTP(52).EQ.3) I3=1    
    
C...Convert to Lambda in CWZ scheme (approximately linear relation).    
          ALAM=0.75*PARP(1) 
          TPMS=PMAS(6,1)    
          QINI=PARP(52) 
          QMAX=PARP(53) 
          XMIN=PARP(54) 
    
C...Initialize evolution (perform calculation or read results from  
C...file).  
C...Remove C* on following two lines to enable Tung initialization. 
C*        CALL PDFSET(I1,IHDRN,ALAM,TPMS,QINI,QMAX,XMIN,NU,HEADER,  
C*   &    I2,I3,IRET,IRR)   
          INIT=1    
        ENDIF   
    
C...Put into output array.  
        Q=SQRT(Q2)  
        DO 200 I=-6,6   
        FIXQ=0. 
C...Remove C* on following line to enable structure function call.  
C*      FIXQ=MAX(0.,PDF(10,1,I,X,Q,IR)) 
  200   XPQ(I)=X*FIXQ   
    
C...Change order of u and d quarks from Tung to PYTHIA convention.  
        XPS=XPQ(1)  
        XPQ(1)=XPQ(2)   
        XPQ(2)=XPS  
        XPS=XPQ(-1) 
        XPQ(-1)=XPQ(-2) 
        XPQ(-2)=XPS 
      ENDIF 
    
      RETURN    
      END   
