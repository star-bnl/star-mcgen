C: definitions from /star/u2c/nevski/bin/geant3.def
*******************************************************************************
      SUBROUTINE HEPTUP                                                   
      CALL HERMES(0)                                                     
      CALL HEPHELP                                                        
      END                                                                 
*******************************************************************************
      SUBROUTINE HEPEXAMPLE                                               
      INTEGER MM(2)/0,0/,DD(2)/0,0/,IW(2)/90,91/,PIPE,P/0/                
      REAL PP(3),VV(3)                                                    
*   call HEPfat
*   call HEPdense
      NP=12000                                                            
C *                                                                       
      DO 5011 J=1,3                                                       
         CALL HEPEVENT ('hijing',12,NP, 3.,1.5,100.,0.1, 197.,97.,197.,   
     *   97.)                                                             
C    *                                                                    
         DO 5021 I=1,NP                                                   
            PP(1)=1                                                       
            PP(2)=2                                                       
            PP(3)=3*RNDM()                                                
            VV(1)=0                                                       
            VV(2)=0                                                       
            VV(3)=.01*RNDM()                                              
            CALL HEPPART (I,1,421,MM,DD,PP,10.,1.,VV,0.)                  
5021     CONTINUE                                                         
5022     CONTINUE                                                         
5011  CONTINUE                                                            
5012  CONTINUE                                                            
      CALL HEPEND('z')                                                    
      END                                                                 
*******************************************************************************
      SUBROUTINE HEPHELP                                                  
      PRINT *,'*********************************************************  
     ******************'                                                  
      PRINT *,'* A utility set to write a standard HEPEVNT n-tuple 999 i  
     *n evgen.run.nt  *'                                                  
      PRINT *,'*********************************************************  
     ******************'                                                  
      PRINT *,'*          mandatory Calles:                               
     *                *'                                                  
      PRINT *,'* HEPEvent (generator, run, Npart, B,F,Et,At, A1,Z1,A2,Z2  
     *) - new event   *'                                                  
      PRINT *,'* HEPPart  (ipa,ist,pdg, moth,idau,pp, Ep,Am,vv,vt) - wri  
     *te new particle *'                                                  
      PRINT *,'* HEPEnd   (option) - close ntuple and compress it on "z"  
     * option         *'                                                  
      PRINT *,'*          optional Calls:                                 
     *                *'                                                  
      PRINT *,'* HEPdens  - dense packing: no mother-daughter relations,  
     * no vertex info *'                                                  
      PRINT *,'* HEPfat   - fat packing: precise vertex info              
     *                *'                                                  
      PRINT *,'* HEPnormal- return to default packing: vertex limited wi  
     *thin 1 mk       *'                                                  
      PRINT *,'*          experts Call:                                   
     *                *'                                                  
      PRINT *,'* HEPmax (IPdg, IRef, NPart, Vxyzt, Nbit) - set limits on  
     * HEP variables  *'                                                  
      PRINT *,'*********************************************************  
     ******************'                                                  
      END                                                                
*************************************************************************
      SUBROUTINE HEPEVENT (GENERATOR, RUN, NPART, B,F,ET,AT, A1,Z1,A2,    
     *Z2)                                                                 
      IMPLICIT NONE                                                       
      INTEGER IVER/11/,IPMX/1000000/,MXRF/1/,MXPA/32000/,NV/16/           
      INTEGER MAXIP,MAXRF,MAXPA,MAXNV,K,IC/0/,ID/999/                     
      REAL VXMAX,VXMX/0.001/,VXMM,VXRM                                    
* Input parameters:
      CHARACTER GENERATOR*(*)                                             
      INTEGER RUN,NPART,IPA,IST,PDG,MOTH(2),IDAU(2)                       
      REAL B,F,ET,AT,A1,Z1,A2,Z2,PP(3),EP,AM,VV(3),VT                     
* Cernlib related:
      INTEGER NWPAW,IPAW,LENOCC,SYSTEMF                                   
      PARAMETER (NWPAW=1000000)                                           
      COMMON /PAWC/ IPAW(NWPAW)                                           
* Hepevnt related:
      INTEGER NTRACK,ITYPE,IS,L,MREF                                      
      CHARACTER CR*8, OPTION*1, GENER*20, FILE*20                         
      INTEGER NP,IDRUN,IEVT,IDAT,ITIM,IGEN                                
      COMMON /HEP_HEAD/ NP,IDRUN,IEVT,IDAT,ITIM,IGEN                      
      INTEGER IP,ISTAT,IPDG,MOT1,MOT2,IDA1,IDA2                           
      REAL PXYZ,ENER,MASS,VXYZ,VTIME                                      
      COMMON /HEP_PART/ IP,ISTAT,IPDG,MOT1,MOT2,IDA1,IDA2,PXYZ(3),ENER,   
     *MASS, VXYZ(3),VTIME                                                 
* Local:
      PARAMETER (K=7)                                                     
      INTEGER I,CC(K)                                                     
      CHARACTER*20 GG(K),FF(K)                                            
      DATA (GG(I),FF(I),CC(I),I=1,K) / 'nexus' , 'nexus.inp' , 7
     *, 'venus' , 'optns.dat' , 6, 'hijing' , 'hijev.inp' , 5, 'mevsim' 
     *, 'mult_gen.in', 4, 'rqmd' , 'rqmd.inp' , 3, 'pythia' 
     *, 'pythia.data', 1, 'user' , 'user.input' , 0/  
      integer iquest
      COMMON/QUEST/IQUEST(100)
      integer IQUEST10_SAVED
      LOGICAL FIRST/.TRUE./                                               
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Check FIRST                                                             
      IF (FIRST) THEN                                                     
      FIRST=.FALSE.                                                       
      MREF = MXRF*MXPA                                                    
      CALL VZERO(IP,16)                                                   
      GENER=GENERATOR                                                     
      CALL CUTOL(GENER)                                                   
*   Is HBOOK and memory initialised ?
C *                                                                       
C    Check IPAW(1)==0                                                     
         IF (IPAW(1).EQ.0) THEN                                           
         PRINT *,' HBOOK initialised for HEP'                             
         CALL HLIMIT(NWPAW)                                               
      END IF                                                              
      CALL HEPNUMBER (RUN,CR)                                             
      IDRUN=RUN                                                           
      VXMM=0                                                              
C Check VXMX>0                                                            
      IF (VXMX.GT.0) VXMM=-VXMX                                           
      FILE='evgen.'//CR(1:LENOCC(CR))//'.nt'                              
      IQUEST10_SAVED = IQUEST(10) 
      IQUEST(10) = 65000
      CALL HROPEN(99,'HEPEVNT',FILE,'QN7',4096,IS)
      IQUEST(10) = IQUEST10_SAVED 

*yf   CALL HROPEN (99,'HEPEVNT',FILE,'N',1024,IS)                         
      CALL RZCDIR ('//HEPEVNT', ' ')                                      
      CALL HBSET ('BSIZE',4096, IS)                                       
      CALL HBNT (ID,'HEPEVNT',' ')                                        
      CALL HEPBNAME (ID,IP, 'itrac' , 0, -1, MXPA)                        
      CALL HEPBNAME (ID,ISTAT,'istat' , 0, -1, IVER)                      
      CALL HEPBNAME (ID,IPDG, 'ipdg' , 0, -IPMX,IPMX)                     
      CALL HEPBNAME (ID,MOT1, 'moth1' , 0, -1, MREF)                      
      CALL HEPBNAME (ID,MOT2, 'moth2' , 0, -MREF, 1)                      
      CALL HEPBNAME (ID,IDA1, 'idau1' , 0, -1, MREF)                      
      CALL HEPBNAME (ID,IDA2, 'idau2' , 0, -1, MREF)                      
      CALL HEPBNAME (ID,PXYZ, 'Pxyz(3)' , 0, 0, 0)                        
      CALL HEPBNAME (ID,ENER, 'ener' , 0, 0, 0)                           
      CALL HEPBNAME (ID,MASS, 'mass:R:' ,16, -1, 10)                      
C  mm                                                                     
      CALL HEPBNAME (ID,VXYZ, 'Vxyz(3):R:' ,NV, VXMM, VXMX)               
C mm/c                                                                    
      CALL HEPBNAME (ID,VTIME,'Vtime:R:' ,NV, 0, VXMX)                    
*   1 mm/c=0.33 ns;   ct=3.e11: tmax=5000 -> 17 ns
C *                                                                       
C    Loop here                                                            
         DO 5011 IGEN=1,K-1                                               
         L=MIN(LENOCC(GENER),LENOCC(GG(IGEN)))                            
         IF (GENER(1:L).EQ.GG(IGEN)(1:L))GO TO 5012                       
5011  CONTINUE                                                            
5012  CONTINUE                                                            
      CALL HEPINPUT(FF(IGEN))                                             
      ENDIF                                                               
*
      IEVT=IEVT+1                                                         
      IP=NPART                                                            
      ISTAT=IVER                                                          
      IPDG=IPMX                                                          
      CALL DATIME (IDAT,ITIM)                                             
      IPDG=IPDG-1                                                         
      PXYZ(1)=RUN                                                         
      PXYZ(2)=IEVT                                                        
      PXYZ(3)=IDAT                                                        
      ENER=ITIM                                                           
      MASS=CC(IGEN)                                                       
      CALL HFNT(ID)                                                       
      IPDG=IPDG-1                                                         
      PXYZ(1)=B                                                           
      PXYZ(2)=F                                                           
      PXYZ(3)=ET                                                          
      ENER=AT                                                             
      MASS=1                                                              
      CALL HFNT(ID)                                                       
      IPDG=IPDG-1                                                         
      PXYZ(1)=A1                                                          
      PXYZ(2)=Z1                                                          
      PXYZ(3)=A2                                                          
      ENER=Z2                                                             
      MASS=2                                                              
      CALL HFNT(ID)                                                       
      NP=NPART                                                            
      RETURN                                                              
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENTRY HEPPART (IPA,IST,PDG,MOTH,IDAU,PP,EP,AM,VV,VT)                
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL RZCDIR('//HEPEVNT',' ')                                        
      IP = IPA                                                            
      ISTAT = IST                                                         
      IPDG = PDG                                                          
      MOT1 = MOTH(1)                                                      
      MOT2 = MOTH(2)                                                      
      IDA1 = IDAU(1)                                                      
      IDA2 = IDAU(2)                                                      
      CALL UCOPY(PP,PXYZ,3)                                               
      CALL UCOPY(VV,VXYZ,3)                                               
      VTIME = VT                                                          
      MASS = AM                                                           
      ENER = EP                                                           
*  if (ipa==Np) Ip=-1
      CALL HFNT(ID)                                                       
      RETURN                                                              
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENTRY HEPEND(OPTION)                                                
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL HROUT(0,IC,'NT')                                               
      CALL HREND('HEPEVNT')                                               
C Check OPTION=='z' | OPTION=='Z'                                         
      IF (OPTION.EQ.'z' .OR. OPTION.EQ.'Z') I=SYSTEMF('gzip -f '//        
     *FILE(1:LENOCC(FILE)))                                               
      RETURN                                                              
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENTRY HEPNORMAL                                                     
      MXRF=1                                                              
      NV=16                                                               
      RETURN                                                              
      ENTRY HEPDENSE                                                      
      MXRF=0                                                              
      NV= 1                                                               
      RETURN                                                              
      ENTRY HEPFAT                                                        
      MXRF=1                                                              
      VXMX=0                                                              
      RETURN                                                              
      ENTRY HEPMAX (MAXIP, MAXRF, MAXPA, VXMAX, MAXNV)                    
      IPMX=MAXIP                                                          
      MREF=MAXRF                                                          
      MXPA=MAXPA                                                          
      VXMX=VXMAX                                                          
      NV=MAXNV                                                            
      RETURN                                                              
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END                                                                 
*************************************************************************
      SUBROUTINE HEPINPUT (INPUT)                                         
      CHARACTER INPUT*(*),LINE*128                                        
      INTEGER LENOCC,LI/98/,ID/998/                                       
      CLOSE (LI)                                                          
      CALL HBNT (ID,'HEPinput',' ')                                       
      CALL HBNAMC (ID,'HEPinput',LINE, 'line(4):C*32:')                   
      OPEN (LI,FILE=INPUT(1:LENOCC(INPUT)),STATUS='OLD',ERR=  5010)       
C *                                                                       
C    Loop here                                                            
5021     CONTINUE                                                         
         READ (LI,'(a)',ERR=5010,END=5010) LINE                           
         CALL HFNT(ID)                                                    
      GO TO 5021                                                          
5022  CONTINUE                                                            
5010  CLOSE (LI)                                                          
      END                                                                 
*************************************************************************
      SUBROUTINE HEPBNAME(ID,VAR,FORM,NB,IA,IB)                           
      IMPLICIT NONE                                                       
      INTEGER LENOCC,INDEX,NB,ID,IA,IB,VAR,L                              
      CHARACTER C*8,CC*80,FORM*(*)                                        
      CC=FORM                                                             
      L=INDEX(FORM,':')                                                   
C Check L>0                                                               
      IF (L.GT.0) CC=FORM(1:L-1)                                          
C *                                                                       
C    Check IA!=0 | IB!=0                                                  
         IF (IA.NE.0 .OR. IB.NE.0) THEN                                   
         CC = FORM                                                        
         CALL HEPNUMBER(NB,C)                                             
C    Check NB>0                                                           
         IF (NB.GT.0) CC=CC(1:LENOCC(CC))//C(1:LENOCC(C))                 
C    Check INDEX(CC,':')>0                                                
         IF (INDEX(CC,':').GT.0) CC=CC(1:LENOCC(CC))//':'                 
         CALL HEPNUMBER(IA,C)                                             
         CC=CC(1:LENOCC(CC))//'['//C(1:LENOCC(C))                         
         CALL HEPNUMBER(IB,C)                                             
         CC=CC(1:LENOCC(CC))//','//C(1:LENOCC(C))//']'                    
      END IF                                                              
      CALL HBNAME(ID,'particle',VAR,CC(1:LENOCC(CC)))                     
      END                                                                 
*************************************************************************
      SUBROUTINE HEPNUMBER(NUM,CNUM)                                      
      IMPLICIT NONE                                                       
      CHARACTER CNUM*(*),S*10                                             
      INTEGER NUM,L,I,I1,I2                                               
C *                                                                       
C    Check ABS(NUM)<=1000000                                              
         IF (ABS(NUM).LE.1000000) THEN                                    
         WRITE (S, * ) NUM                                                
      ELSE                                                                
         WRITE (S,'(f10.6)') NUM                                          
      END IF                                                              
      I1=10                                                               
      I2=1                                                                
C *                                                                       
C    Loop here                                                            
         DO 5011 I=1,10                                                   
C    Skip Unless S(I:I)!=' '                                              
         IF (S(I:I).EQ.' ')GO TO 5011                                     
         I1=MIN(I1,I)                                                     
         I2=MAX(I2,I)                                                     
5011  CONTINUE                                                            
5012  CONTINUE                                                            
      CNUM=S(I1:I2)                                                       
      L=I2-I1+1                                                           
1     CONTINUE                                                            
      END                                                                 
