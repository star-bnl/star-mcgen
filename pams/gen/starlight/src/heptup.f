C: definitions from /star/u2c/nevski/bin/geant3.def
*******************************************************************************
      SUBROUTINE HEPTUP                                                   2
      CALL HERMES (0)                                                     3
      CALL HEPHELP                                                        4
      END                                                                 5
*******************************************************************************
      SUBROUTINE HEPEXAMPLE                                               8
      INTEGER MM(2)/0,0/,DD(2)/0,0/,IW(2)/90,91/,PIPE,P/0/                9
      REAL PP(3),VV(3)                                                    10
*   call HEPfat
*   call HEPdense
      INTEGER N/1/,K                                                      13
C *                                                                       15
      DO 5011 K=1,N                                                       15
         NP=12000                                                         16
C    *                                                                    17
         DO 5021 J=1,10                                                   17
            CALL HEPEVENT ('hijing',0,NP, 3.,1.5,100.,0.1, 197.,97.,      18
     *      197.,97.)                                                     18
C       *                                                                 19
            DO 5031 I=1,NP                                                19
               PP(1)=1                                                    19
               PP(2)=2                                                    19
               PP(3)=3*RNDM()                                             19
               VV(1)=0                                                    19
               VV(2)=0                                                    20
               VV(3)=.01*RNDM()                                           20
               CALL HEPPART (I,1,421,MM,DD,PP,10.,1.,VV,0.)               21
5031        CONTINUE                                                      22
5032        CONTINUE                                                      22
5021     CONTINUE                                                         23
5022     CONTINUE                                                         23
5011  CONTINUE                                                            24
5012  CONTINUE                                                            24
      N=2*N                                                               25
*   call hepend('z')
      END                                                                 27
*******************************************************************************
      SUBROUTINE HEPLIGHT                                                 30
      IMPLICIT NONE                                                       31
      REAL PP(3),VV(3)/0,0,0/,E1,MDEC,VT/0/                               32
      INTEGER MM(2)/0,0/,DD(2)/0,0/                                       33
      CALL HEPEVENT ('starlight',1,2,99.,0.5,100.,0.1,197.,97.,197.,      35
     *97.)                                                                35
      E1=0.126805                                                         35
      MDEC=0.105658                                                       35
      PP(1)=-5.99571E-02                                                  36
      PP(2)=-2.52263E-02                                                  36
      PP(3)=2.61670E-02                                                   36
      CALL HEPPART (1,1,13,MM,DD,PP,E1,MDEC,VV,VT)                        36
      PP(1)=6.16919E-02                                                   38
      PP(2)=2.06726E-02                                                   38
      PP(3)=0.585180                                                      38
      CALL HEPPART (2,1,-13,MM,DD,PP,E1,MDEC,VV,VT)                       38
      END                                                                 42
*******************************************************************************
      SUBROUTINE HEPHELP                                                  45
      PRINT *,'*********************************************************  47
     ******************'                                                  47
      PRINT *,'* A utility set to write a standard HEPEVNT n-tuple 999 i  48
     *n evgen.run.nt  *'                                                  48
      PRINT *,'*********************************************************  49
     ******************'                                                  49
      PRINT *,'*          mandatory Calles:                               50
     *                *'                                                  50
      PRINT *,'* HEPEvent (generator, run, Npart, B,F,Et,At, A1,Z1,A2,Z2  51
     *) - new event   *'                                                  51
      PRINT *,'* HEPPart  (ipa,ist,pdg, moth,idau,pp, Ep,Am,vv,vt) - wri  52
     *te new particle *'                                                  52
      PRINT *,'* HEPEnd   (option) - close ntuple and compress it on "z"  53
     * option         *'                                                  53
      PRINT *,'*          optional Calls:                                 54
     *                *'                                                  54
      PRINT *,'* HEPdens  - dense packing: no mother-daughter relations,  55
     * no vertex info *'                                                  55
      PRINT *,'* HEPfat   - fat packing: precise vertex info              56
     *                *'                                                  56
      PRINT *,'* HEPnormal- return to default packing: vertex limited wi  57
     *thin 1 mk       *'                                                  57
      PRINT *,'*          experts Call:                                   58
     *                *'                                                  58
      PRINT *,'* HEPmax (IPdg, IRef, NPart, Vxyzt, Nbit) - set limits on  59
     * HEP variables  *'                                                  59
      PRINT *,'*********************************************************  60
     ******************'                                                  60
      END                                                                 61
*******************************************************************************
      SUBROUTINE HEPEVENT (GENERATOR, RUN, NPART, B,F,ET,AT, A1,Z1,A2,    65
     *Z2)                                                                 65
      IMPLICIT NONE                                                       66
      INTEGER IVER/11/,IPMX/1000000/,MXRF/1/,MXPA/32000/,NV/16/           67
      INTEGER MAXIP,MAXRF,MAXPA,MAXNV,K,IC/0/,ID/999/                     68
      REAL VXMAX,VXMX/0.001/,VXMM,VXRM                                    69
* Input parameters:
      CHARACTER GENERATOR*(*)                                             71
      INTEGER RUN,NPART,IPA,IST,PDG,MOTH(2),IDAU(2)                       72
      REAL B,F,ET,AT,A1,Z1,A2,Z2,PP(3),EP,AM,VV(3),VT                     73
* Cernlib related:
      INTEGER NWPAW,IPAW, LENOCC, SYSTEMF, LREC,BSIZE                     76
      PARAMETER (NWPAW=2 000 000, LREC=4096, BSIZE=LREC)                  77
      COMMON /PAWC/ IPAW(NWPAW)                                           78
* Hepevnt related:
      INTEGER NTRACK,ITYPE,IS,L,MREF,LUH/99/,LUR/98/,LUN/97/              81
      CHARACTER CR*8, OPTION*1, GENER*20, FILE*20                         82
      INTEGER NP,IDRUN,IEVT,IDAT,ITIM,IGEN                                83
      COMMON /HEP_HEAD/ NP,IDRUN,IEVT,IDAT,ITIM,IGEN                      84
      INTEGER IP,ISTAT,IPDG,MOT1,MOT2,IDA1,IDA2                           85
      REAL PXYZ,ENER,MASS,VXYZ,VTIME                                      86
      COMMON /HEP_PART/ IP,ISTAT,IPDG,MOT1,MOT2,IDA1,IDA2,PXYZ(3),ENER,   88
     *MASS, VXYZ(3),VTIME                                                 88
* Local:
      PARAMETER (K=7)                                                     90
      INTEGER I,CC(K)                                                     91
      CHARACTER*20 GG(K),FF(K)                                            92
      DATA (GG(I),FF(I),CC(I),I=1,K) / 'starlight', 'starlight.in', 7, '  100
     *venus' , 'optns.dat' , 6, 'hijing' , 'hijev.inp' , 5, 'mevsim' , '  100
     *mult_gen.in' , 4, 'rqmd' , 'rqmd.inp' , 3, 'pythia' , 'pythia.data  100
     *' , 1, 'user' , 'user.input' , 0/                                   100
      LOGICAL FIRST/.TRUE./                                               101
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Check FIRST                                                             103
      IF (FIRST) THEN                                                     103
      FIRST=.FALSE.                                                       104
      MREF = MXRF*MXPA                                                    105
      CALL VZERO(IP,16)                                                   106
      GENER=GENERATOR                                                     107
      CALL CUTOL(GENER)                                                   108
*   Is HBOOK and memory initialised ?
C *                                                                       110
C    Check IPAW(1)==0                                                     110
         IF (IPAW(1).EQ.0) THEN                                           110
         PRINT *,' HBOOK initialised for HEP'                             110
         CALL HLIMIT(NWPAW)                                               111
      END IF                                                              111
*   print *,' hbook initialised with len = ',Ipaw(1)
      IDRUN=RUN                                                           114
C Check RUN==0                                                            115
      IF (RUN.EQ.0) THEN                                                  115
5010  CALL HEPNUMBER (RUN,CR)                                             116
      OPEN (LUN,FILE='next',STATUS='NEW',ERR=5010)                        117
      OPEN (LUR,FILE='run', STATUS='OLD',ERR=  5020)                      118
      READ (LUR,*) IDRUN                                                  119
      CLOSE (LUR)                                                         120
5020  IDRUN = IDRUN+1                                                     121
      WRITE (LUN,*) IDRUN                                                 122
      CLOSE (LUN)                                                         123
      IS = SYSTEMF ('mv -f next run')                                     124
      ENDIF                                                               125
      CALL HEPNUMBER (IDRUN,CR)                                           127
      VXMM=0                                                              127
C Check VXMX>0                                                            128
      IF (VXMX.GT.0) VXMM=-VXMX                                           128
      FILE='evgen.'//CR(1:LENOCC(CR))//'.nt'                              129
      CALL HROPEN (LUH,'HEPEVNT',FILE,'N7',LREC,IS)                       130
C Check IS!=0                                                             131
      IF (IS.NE.0) STOP ' HEPTUPLE: Can not open output file '            131
      CALL RZCDIR ('//HEPEVNT', ' ')                                      132
      CALL HBSET ('BSIZE',BSIZE,IS)                                       133
C Check IS!=0                                                             134
      IF (IS.NE.0) STOP ' HEPTUPLE: Can not set buffer size '             134
      CALL HBNT (ID,'HEPEVNT',' ')                                        136
      CALL HEPBNAME (ID,IP, 'itrac' , 0, -1, MXPA)                        137
      CALL HEPBNAME (ID,ISTAT,'istat' , 0, -1, IVER)                      138
      CALL HEPBNAME (ID,IPDG, 'ipdg' , 0, -IPMX,IPMX)                     139
      CALL HEPBNAME (ID,MOT1, 'moth1' , 0, -1, MREF)                      140
      CALL HEPBNAME (ID,MOT2, 'moth2' , 0, -MREF, 1)                      141
      CALL HEPBNAME (ID,IDA1, 'idau1' , 0, -1, MREF)                      142
      CALL HEPBNAME (ID,IDA2, 'idau2' , 0, -1, MREF)                      143
      CALL HEPBNAME (ID,PXYZ, 'Pxyz(3)' , 0, 0, 0)                        144
      CALL HEPBNAME (ID,ENER, 'ener' , 0, 0, 0)                           145
      CALL HEPBNAME (ID,MASS, 'mass:R:' ,16, -1, 10)                      146
C  mm                                                                     147
      CALL HEPBNAME (ID,VXYZ, 'Vxyz(3):R:' ,NV, VXMM, VXMX)               147
C mm/c                                                                    148
      CALL HEPBNAME (ID,VTIME,'Vtime:R:' ,NV, 0, VXMX)                    148
*   1 mm/c=0.33 ns;   ct=3.e11: tmax=5000 -> 17 ns
C *                                                                       152
C    Loop here                                                            152
         DO 5031 IGEN=1,K-1                                               152
         L=MIN(LENOCC(GENER),LENOCC(GG(IGEN)))                            152
         IF (GENER(1:L).EQ.GG(IGEN)(1:L))GO TO 5032                       153
5031  CONTINUE                                                            154
5032  CONTINUE                                                            154
      CALL HEPINPUT(FF(IGEN))                                             154
      ENDIF                                                               155
*
      IEVT=IEVT+1                                                         156
      IP=NPART                                                            156
      ISTAT=IVER                                                          156
      IPDG=IPMX                                                           156
      CALL DATIME (IDAT,ITIM)                                             156
      IPDG=IPDG-1                                                         157
      PXYZ(1)=IDRUN                                                       157
      PXYZ(2)=IEVT                                                        157
      PXYZ(3)=IDAT                                                        157
      PXYZ(4)=ITIM                                                        157
      PXYZ(5)=CC(IGEN)                                                    157
      CALL HFNT(ID)                                                       157
      IPDG=IPDG-1                                                         158
      PXYZ(1)=B                                                           158
      PXYZ(2)=F                                                           158
      PXYZ(3)=ET                                                          158
      PXYZ(4)=AT                                                          158
      PXYZ(5)=1                                                           158
      CALL HFNT(ID)                                                       158
      IPDG=IPDG-1                                                         159
      PXYZ(1)=A1                                                          159
      PXYZ(2)=Z1                                                          159
      PXYZ(3)=A2                                                          159
      PXYZ(4)=Z2                                                          159
      PXYZ(5)=2                                                           159
      CALL HFNT(ID)                                                       159
      NP=NPART                                                            160
      RETURN                                                              162
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENTRY HEPPART (IPA,IST,PDG,MOTH,IDAU,PP,EP,AM,VV,VT)                165
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL RZCDIR('//HEPEVNT',' ')                                        168
      IP = IPA                                                            169
      ISTAT = IST                                                         170
      IPDG = PDG                                                          171
      MOT1 = MOTH(1)                                                      172
      MOT2 = MOTH(2)                                                      173
      IDA1 = IDAU(1)                                                      174
      IDA2 = IDAU(2)                                                      175
      CALL UCOPY(PP,PXYZ,3)                                               176
      CALL UCOPY(VV,VXYZ,3)                                               177
      VTIME = VT                                                          178
      MASS = AM                                                           179
      ENER = EP                                                           180
*  if (ipa==Np) Ip=-1
      CALL HFNT(ID)                                                       182
      RETURN                                                              183
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENTRY HEPEND(OPTION)                                                186
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL HROUT(0,IC,'NT')                                               188
      CALL HREND('HEPEVNT')                                               190
C Check OPTION=='z' | OPTION=='Z'                                         191
      IF (OPTION.EQ.'z' .OR. OPTION.EQ.'Z') I=SYSTEMF('gzip -f '//        191
     *FILE(1:LENOCC(FILE)))                                               191
      RETURN                                                              192
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENTRY HEPNORMAL                                                     194
      MXRF=1                                                              194
      NV=16                                                               194
      RETURN                                                              195
      ENTRY HEPDENSE                                                      195
      MXRF=0                                                              195
      NV= 1                                                               195
      RETURN                                                              196
      ENTRY HEPFAT                                                        196
      MXRF=1                                                              196
      VXMX=0                                                              196
      RETURN                                                              197
      ENTRY HEPMAX (MAXIP, MAXRF, MAXPA, VXMAX, MAXNV)                    198
      IPMX=MAXIP                                                          198
      MREF=MAXRF                                                          198
      MXPA=MAXPA                                                          198
      VXMX=VXMAX                                                          198
      NV=MAXNV                                                            198
      RETURN                                                              200
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END                                                                 202
*************************************************************************
      SUBROUTINE HEPINPUT (INPUT)                                         206
      CHARACTER INPUT*(*),LINE*128                                        207
      INTEGER LENOCC,LI/98/,ID/998/                                       212
      CLOSE (LI)                                                          212
      CALL HBNT (ID,'HEPinput',' ')                                       212
      CALL HBNAMC (ID,'HEPinput',LINE, 'line(4):C*32:')                   212
      OPEN (LI,FILE=INPUT(1:LENOCC(INPUT)),STATUS='OLD',ERR=  5010)       212
C *                                                                       212
C    Loop here                                                            212
5021     CONTINUE                                                         212
         READ (LI,'(a)',ERR=5010,END=5010) LINE                           212
         CALL HFNT(ID)                                                    213
      GO TO 5021                                                          213
5022  CONTINUE                                                            213
5010  CLOSE (LI)                                                          214
      END                                                                 215
*************************************************************************
      SUBROUTINE HEPBNAME(ID,VAR,FORM,NB,IA,IB)                           218
      IMPLICIT NONE                                                       219
      INTEGER LENOCC,INDEX,NB,ID,IA,IB,VAR,L                              220
      CHARACTER C*8,CC*80,FORM*(*)                                        221
      CC=FORM                                                             221
      L=INDEX(FORM,':')                                                   221
C Check L>0                                                               222
      IF (L.GT.0) CC=FORM(1:L-1)                                          222
C *                                                                       224
C    Check IA!=0 | IB!=0                                                  224
         IF (IA.NE.0 .OR. IB.NE.0) THEN                                   224
         CC = FORM                                                        225
         CALL HEPNUMBER(NB,C)                                             226
C    Check NB>0                                                           227
         IF (NB.GT.0) CC=CC(1:LENOCC(CC))//C(1:LENOCC(C))                 227
C    Check INDEX(CC,':')>0                                                228
         IF (INDEX(CC,':').GT.0) CC=CC(1:LENOCC(CC))//':'                 228
         CALL HEPNUMBER(IA,C)                                             228
         CC=CC(1:LENOCC(CC))//'['//C(1:LENOCC(C))                         229
         CALL HEPNUMBER(IB,C)                                             229
         CC=CC(1:LENOCC(CC))//','//C(1:LENOCC(C))//']'                    230
      END IF                                                              231
      CALL HBNAME(ID,'particle',VAR,CC(1:LENOCC(CC)))                     232
      END                                                                 233
*************************************************************************
      SUBROUTINE HEPNUMBER(NUM,CNUM)                                      237
      IMPLICIT NONE                                                       238
      CHARACTER CNUM*(*),S*14                                             239
      INTEGER ANUM,NUM,L,I,I1,I2                                          240
      REAL RNUM                                                           241
      EQUIVALENCE (RNUM,ANUM)                                             242
      ANUM=NUM                                                            244
C *                                                                       244
C    Check ABS(NUM)<=1000000                                              244
         IF (ABS(NUM).LE.1000000) THEN                                    244
         WRITE (S, * ) ANUM                                               244
      ELSE                                                                245
         WRITE (S,'(f14.6)') RNUM                                         245
      END IF                                                              246
      I1=14                                                               246
      I2=1                                                                246
C *                                                                       247
C    Loop here                                                            247
         DO 5011 I=1,14                                                   247
C    Skip Unless S(I:I)!=' '                                              247
         IF (S(I:I).EQ.' ')GO TO 5011                                     247
         I1=MIN(I1,I)                                                     247
         I2=MAX(I2,I)                                                     247
5011  CONTINUE                                                            248
5012  CONTINUE                                                            248
      CNUM=S(I1:I2)                                                       248
      L=I2-I1+1                                                           249
1     CONTINUE                                                            250
      END                                                                 250
*CMZ :  2.00/02 19/01/2000  23.56.38  by  Pavel Nevski
*-- Author :
      SUBROUTINE RZMAKE(LUNIN,CHDIR,NWKEY,CHFORM,CHTAG,NRECPIN,CHOPT)     254
*
************************************************************************
*
*           Routine to create a new RZ file
*           To use an already existing file CALL RZFILE
* Input:
*   LUNIN   Logical unit number associated with  the RZ file.   A FORTRAN
*           OPEN statement must precede the call to RZFILE.
*           Starting address of the memory area which will contain the RZ
*           information ('M' option)
*   CHDIR   Character variable specifying  the name of the  top directory
*           to be associated with unit LUN.
*   NWKEY   Number of words associated to a key (maximum 5)
*   CHFORM  Character variable describing each element  of the key vector
*           'B' Bit string but not zero
*           'H' Hollerith (4 characters)
*           'I' Integer (nonzero)
*           Ex: CHFORM='IIH' for NWKEY=3 and the 2 first keys are integer
*               and the third one is Hollerith
*   CHTAG   Character array defined as CHARACTER*8 CHTAG(NWKEY).
*           Each  element of  the  array allows  the  description of  the
*           corresponding element in the key vector with a tag of up to 8
*           characters.
*   NRECP   Number of physical records for primary allocation
*   CHOPT   Character variable specifying the selected options.
*           medium
*             default
*                   Disk
*             'M'   Memory
*                   In this  case the user  must have allocated  at least
*                   NRecp*LRecp words of memory starting at address LUNIN.
*                   LRecp is passes as the first word in this area
*           mode
*             default
*                   Native mode
*             'X'   Exchange mode
*           other
*             'F'   Format NRECP records (unless 'M')
*             'C'   C I/O (unless 'M')
*                   LRECL (words) taken from IQUEST(10)
*             'O'   OLD format for Cycle information (default is NEW)
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 14.09.93 No longer force exchange mode for LINUX
*          : 09.03.94 S.Banerjee (Change in cycle structure)
*          : 30.01.95 J.Shiers. Permit nrecp>65000 for new format
*          : 10.12.97 P.Nevski  Default is NEW
*          : 01.03.99 Perevozchikov  default NRECP for NEW is 1M
************************************************************************
*
*KEEP,ZUNIT.
      COMMON /ZUNIT/ IQREAD,IQPRNT,IQPR2,IQLOG,IQPNCH,IQTTIN,IQTYPE       309
      COMMON /ZUNITZ/IQDLUN,IQFLUN,IQHLUN, NQUSED                         310
*KEEP,ZSTATE.
      COMMON /ZSTATE/QVERSN,NQPHAS,IQDBUG,NQDCUT,NQWCUT,NQERR , NQLOGD,   313
     *NQLOGM,NQLOCK,NQDEVZ,NQOPTS(6)                                      313
*KEEP,RZCL.
      PARAMETER (IQDROP=25)                                               315
      PARAMETER (IQMARK=26)                                               315
      PARAMETER (IQCRIT=27)                                               315
      PARAMETER (IQSYSX=28)                                               315
      INTEGER IQUEST                                                      316
      COMMON/QUEST/IQUEST(100)                                            317
      COMMON /ZVFAUT/IQVID(2),IQVSTA,IQVLOG,IQVTHR(2),IQVREM(2,6)         318
      COMMON /ZEBQ/ IQFENC(4), LQ(100)                                    319
      DIMENSION IQ(92), Q(92)                                             320
      EQUIVALENCE (IQ(1),LQ(9)), (Q(1),IQ(1))                             321
      COMMON /MZCA/ NQSTOR,NQOFFT(16),NQOFFS(16),NQALLO(16), NQIAM ,      325
     *LQATAB,LQASTO,LQBTIS, LQWKTB,NQWKTB,LQWKFZ , MQKEYS(3),NQINIT,      325
     *NQTSYS,NQM99,NQPERM,NQFATA,NQCASE , NQTRAC,MQTRAC(48)               325
      EQUIVALENCE (KQSP,NQOFFS(1))                                        326
      COMMON /MZCB/ JQSTOR,KQT,KQS, JQDIVI,JQDIVR , JQKIND,JQMODE,        330
     *JQDIVN,JQSHAR,JQSHR1,JQSHR2,NQRESV , LQSTOR,NQFEND,NQSTRU,NQREF,    330
     *NQLINK,NQMINR,LQ2END , JQDVLL,JQDVSY,NQLOGL,NQSNAM(6)               330
      DIMENSION IQCUR(16)                                                 331
      EQUIVALENCE (IQCUR(1),LQSTOR)                                       332
      COMMON /MZCC/ LQPSTO,NQPFEN,NQPSTR,NQPREF,NQPLK,NQPMIN,LQP2E ,      340
     *JQPDVL,JQPDVS,NQPLOG,NQPNAM(6) , LQSYSS(10), LQSYSR(10),            340
     *IQTDUM(22) , LQSTA(21), LQEND(20), NQDMAX(20),IQMODE(20) ,          340
     *IQKIND(20),IQRCU(20), IQRTO(20), IQRNO(20) , NQDINI(20),            340
     *NQDWIP(20),NQDGAU(20),NQDGAF(20) , NQDPSH(20),NQDRED(20),           340
     *NQDSIZ(20) , IQDN1(20), IQDN2(20), KQFT, LQFSTA(21)                 340
      DIMENSION IQTABV(16)                                                341
      EQUIVALENCE (IQTABV(1),LQPSTO)                                      342
C
      COMMON /RZCL/ LTOP,LRZ0,LCDIR,LRIN,LROUT,LFREE,LUSED,LPURG ,        345
     *LTEMP,LCORD,LFROM                                                   345
      EQUIVALENCE (LQRS,LQSYSS(7))                                        346
C
*KEEP,RZDIR.
      PARAMETER (NLPATM=100)                                              349
      COMMON /RZDIRN/NLCDIR,NLNDIR,NLPAT                                  350
      COMMON /RZDIRC/CHCDIR(NLPATM),CHNDIR(NLPATM),CHPAT(NLPATM)          351
      CHARACTER*16 CHNDIR, CHCDIR, CHPAT                                  352
C
*KEEP,RZCLUN.
      COMMON /RZCLUN/LUN,LREC,ISAVE,IMODEX,IRELAT,NHPWD,IHPWD(2) ,        356
     *IZRECL,IMODEC,IMODEH                                                356
C
*KEEP,RZK.
      PARAMETER (KUP=5)                                                   363
      PARAMETER (KPW1=7)                                                  363
      PARAMETER (KNCH=9)                                                  363
      PARAMETER (KDATEC=10)                                               363
      PARAMETER (KDATEM=11)                                               363
      PARAMETER (KQUOTA=12)                                               363
      PARAMETER (KRUSED=13)                                               363
      PARAMETER (KWUSED=14)                                               363
      PARAMETER (KMEGA=15)                                                363
      PARAMETER (KRZVER=16)                                               363
      PARAMETER (KIRIN=17)                                                363
      PARAMETER (KIROUT=18)                                               363
      PARAMETER (KRLOUT=19)                                               363
      PARAMETER (KIP1=20)                                                 363
      PARAMETER (KNFREE=22)                                               363
      PARAMETER (KNSD=23)                                                 363
      PARAMETER (KLD=24)                                                  363
      PARAMETER (KLB=25)                                                  363
      PARAMETER (KLS=26)                                                  363
      PARAMETER (KLK=27)                                                  363
      PARAMETER (KLF=28)                                                  363
      PARAMETER (KLC=29)                                                  363
      PARAMETER (KLE=30)                                                  363
      PARAMETER (KNKEYS=31)                                               363
      PARAMETER (KNWKEY=32)                                               363
      PARAMETER (KKDES=33)                                                363
      PARAMETER (KNSIZE=253)                                              363
      PARAMETER (KEX=6)                                                   363
      PARAMETER (KNMAX=100)                                               363
C
*KEEP,RZCYCLE.
*
*     Pointers to cycle content
*
*     KLCYCL : length of cycle block (4,7)
*     KPPCYC : pointer to previous cycle
*     KFRCYC : first record number
*     KSRCYC : secord record number
*     KFLCYC : creation date/time and other stuff
*     KORCYC : offset in first record to data
*     KCNCYC : cycle number
*     KNWCYC : number of words in d/s
*     KKYCYC : key number to which this cycle belongs (only for version 1)
*     KVSCYC : version of RZ cycles structure (0, 1)
*
      INTEGER KLCYCL, KPPCYC, KFRCYC, KSRCYC, KFLCYC, KORCYC, KCNCYC,     381
     *KNWCYC, KKYCYC, KVSCYC                                              381
      COMMON/RZCYCLE/KLCYCL, KPPCYC, KFRCYC, KSRCYC, KFLCYC, KORCYC,      383
     *KCNCYC, KNWCYC, KKYCYC, KVSCYC                                      383
*KEND.
* Nrecp: If non zero, overwrites default NHrecp in AGI version of RZMAKE.
*        Otherwise it is 32k in standard cernlib and 1M in AGI.
*        In principal it may be passed in IQUEST(10) by calling HROPEN
*        with Copt='Q', but HRFILE anyway resets it to 100<NQUOTA<65000.
      INTEGER NHRECP                                                      389
      COMMON /RZNHRECP/ NHRECP                                            390
      CHARACTER CHOPT*(*),CHDIR*(*),CHFORM*(*)                            391
      CHARACTER*16 CHTOP                                                  392
      CHARACTER*(*) CHTAG(*)                                              393
      DIMENSION IOPTV(6),IHDIR(2)                                         394
      EQUIVALENCE (IOPTM,IOPTV(1)), (IOPTX,IOPTV(2)) , (IOPTF,IOPTV(3)),  397
     * (IOPTC,IOPTV(4)) , (IOPTN,IOPTV(5)), (IOPTO,IOPTV(6))              397
      INTEGER HMEM/0/                                                     398
*
*-----------------------------------------------------------------------
*
*KEEP,Q$JBIT.
      JBIT(IZW,IZP) = AND(ISHFTR(IZW,IZP-1),1)                            403
      JBYT(IZW,IZP,NZB) = ISHFTR(LSHIFT(IZW,33-IZP-NZB),32-NZB)           404
*KEND.
      IQUEST(1)=0                                                         406
      LOGLV = MIN(NQLOGD,4)                                               407
      LOGLV = MAX(LOGLV,-3)                                               408
*
      CALL UOPTC(CHOPT,'MXFCNO',IOPTV)                                    410
      IOPTN = 1-IOPTO                                                     413
      IMODEX = IOPTX                                                      414
      IMODEC = IOPTC                                                      415
      IRELAT = 0                                                          416
      LUNP = LUNIN                                                        417
      IF (IOPTC.NE.0) LUNP = IQUEST(11)                                   418
*
**              VP correction on NRECP default is 1M for new format
*
      NRECP = NRECPIN                                                     423
C Check IOPTN.NE.0                                                        424
C NEW RZ FORMAT                                                           424
      IF (IOPTN.NE.0) THEN                                                424
      M = 1 000 000                                                       425
C Check NHRECP.GT.0                                                       426
      IF (NHRECP.GT.0) M=NHRECP*1024                                      426
      NRECP = MAX(NRECPIN,M)                                              427
      ENDIF                                                               428
*
*                Check NWKEY and NRECP
*
      IF (NWKEY.LE.0.OR.NWKEY.GT.KNMAX) THEN                              432
      IF (LOGLV.GE.-2) WRITE(IQLOG,9010)                                  433
9010  FORMAT(' RZMAKE. NWKEY input value is invalid')                     434
      IQUEST(1) =1                                                        435
      IQUEST(11)=NWKEY                                                    436
      GO TO 99                                                            437
      ENDIF                                                               438
      IF (NRECP.LT.2.OR.(NRECP.GT.65000.AND.IOPTN.EQ.0)) THEN             440
      IF (LOGLV.GE.-2) WRITE(IQLOG,9011)                                  441
9011  FORMAT(' RZMAKE. NRECP input value is invalid')                     442
      IQUEST(1) =1                                                        443
      IQUEST(11)=NRECP                                                    444
      GO TO 99                                                            445
      ENDIF                                                               446
*
*     Find record length (as specified in the OPEN statement)
*
      IF (IOPTM.NE.0) THEN                                                450
*
*          A, Memory option. LUNIN itself IS the buffer and contains the
*                            value of LUNP as the block length (what a mess)
      LUNP = HMEM+1                                                       454
      LRECP = LUNIN                                                       455
      IF (LRECP.LT.100.OR.LRECP.GT.10000) LRECP=1024                      456
      ELSE                                                                458
*
*          B, Standard option DISK. Use information as specified
*             in the Fortran OPEN statement
*
C Check IOPTC.EQ.0                                                        463
      IF (IOPTC.EQ.0) THEN                                                463
      INQUIRE(UNIT=LUNP,RECL=LRECB)                                       464
*
      LRECP=LRECB/4                                                       466
      ELSE                                                                467
*
*     Take LRECL from IQUEST(10) in case of C I/O option
*
      LRECP = IQUEST(10)                                                  471
      ENDIF                                                               472
      ENDIF                                                               473
*
      LUN = LUNP                                                          475
      IZRECL = LRECP                                                      476
      IF (LUN.LE.0.AND.IOPTM.EQ.0) THEN                                   477
      IF (LOGLV.GE.-2) WRITE(IQLOG,9012)                                  478
9012  FORMAT(' RZMAKE. LUN input value is invalid')                       479
      IQUEST(1) =1                                                        480
      IQUEST(11)=LUN                                                      481
      GO TO 99                                                            482
      ENDIF                                                               483
      IF (LRECP.LT.50) THEN                                               484
      IF (LOGLV.GE.-2) WRITE(IQLOG,9013)                                  485
9013  FORMAT(' RZMAKE. LRECP input value less than 50')                   486
      IQUEST(1) =1                                                        487
      IQUEST(11)=LRECP                                                    488
      GO TO 99                                                            489
      ENDIF                                                               490
*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*          Save existing material (if any)
*
      CALL RZSAVE                                                         495
*
      IF (LOGLV.GE.0) WRITE(IQLOG,9014) LUNP,NRECP,LRECP,CHOPT            497
9014  FORMAT(' RZMAKE. Unit ',I6,' Initializing Nrec=',I6, ' with LREC=   499
     *',I6,', OPT= ',A)                                                   499
      CALL MZSDIV (0,-7)                                                  500
*
*           Check if LUN not already defined
*
      LRZ=LQRS                                                            504
10    IF (LRZ.NE.0) THEN                                                  505
      IF (IQ(KQSP+LRZ-5).EQ.LUN) THEN                                     506
      IF (LOGLV.GE.-2) WRITE(IQLOG,9015)                                  507
9015  FORMAT(' RZMAKE. Logical unit number already in use')               508
      IQUEST(1) =1                                                        509
      IQUEST(11)=LUN                                                      510
      GO TO 99                                                            511
      ELSE                                                                512
      LRZ=LQ(KQSP+LRZ)                                                    513
      GO TO 10                                                            514
      ENDIF                                                               515
      ENDIF                                                               516
*
*            First call to RZMAKE, create link area
*
      IF (LQRS.EQ.0) THEN                                                 520
      CALL MZLINK(JQPDVS,'RZCL',LTOP,LTOP,LFROM)                          521
      CALL MZBOOK(JQPDVS,LRZ0,LQRS,1,'RZ0 ',2,2,36,2,0)                   522
      IQ(KQSP+LRZ0-5)=0                                                   523
      ISAVE = 1                                                           524
      NHPWD = 0                                                           525
      CALL VBLANK(IHPWD,2)                                                526
      ENDIF                                                               527
      NCHD = LEN(CHDIR)                                                   528
      IF (NCHD.GT.16) NCHD=16                                             529
      CHTOP = CHDIR(1:NCHD)                                               530
*
*            Create control bank
*
      IDTIME=0                                                            534
      CALL RZDATE(IDTIME,IDATE,ITIME,2)                                   535
      KTAGS = KKDES+(NWKEY-1)/10+1                                        536
      NREC = NRECP                                                        537
      LREC = LRECP                                                        538
      NWREC = (NREC-1)/32 +1                                              539
      NW = 50+NWREC                                                       540
      NRD = (NW-1)/LREC +1                                                541
      NWL = NRD*LREC                                                      542
      LD = KTAGS+2*NWKEY                                                  543
      LB = LD+NRD+1                                                       544
      LS = LB+3+NWREC                                                     545
      LK = LS                                                             546
      LF = LS                                                             547
*
      CALL MZBOOK (JQPDVS,LTOP,LQRS,1,'RZ  ',10,9,NWL,2,0)                549
*
*            Disk or memory
*
      IF (IOPTM.EQ.0) THEN                                                553
      IQ(KQSP+LTOP-5) = LUN                                               554
*
*            C I/O?
      IF (IOPTC.NE.0) CALL SBIT1(IQ(KQSP+LTOP),5)                         557
      ELSE                                                                558
      NMEM=IQ(KQSP+LRZ0)+1                                                559
      IQ(KQSP+LRZ0)=NMEM                                                  560
      IQ(KQSP+LTOP-5)=-NMEM                                               561
      IF (2*NMEM.GT.IQ(KQSP+LRZ0-1)) THEN                                 562
      CALL MZPUSH(JQPDVS,LRZ0,0,10,'I')                                   563
      ENDIF                                                               564
      IQ(KQSP+LRZ0+2*NMEM-1)=LOCF(LUNIN)-LOCF(IQ(1))+1                    565
      IQ(KQSP+LRZ0+2*NMEM )=LRECP                                         566
      LUN=-NMEM                                                           567
      ENDIF                                                               568
*
*            Pre-format file
*
      IF ((IOPTF.NE.0).AND.(IOPTM.EQ.0)) THEN                             572
      DO 100 I=2,NRECP                                                    573
100   CALL RZIODO(LUN,LREC,I,IQ(KQSP+LTOP+1),2)                           574
      IF (IQUEST(1).NE.0) THEN                                            575
      IF (LOGLV.GE.-1) WRITE(IQLOG,1000) I-1                              576
1000  FORMAT(' RZMAKE. Could only pre-format',I6,' records')              577
      IQUEST(1)=0                                                         578
      ENDIF                                                               579
      ENDIF                                                               580
*
*            Write empty record for locks
*
      CALL RZIODO(LUN,LREC,1,IQ(KQSP+LTOP+1),2)                           584
      IF (IQUEST(1).NE.0) GO TO 99                                        585
*
*            Build top-directory parameters
*
      CALL SBIT1(IQ(KQSP+LTOP),2)                                         589
      CALL VBLANK(IQ(KQSP+LTOP+1),4)                                      590
      CALL UCTOH(CHDIR,IQ(KQSP+LTOP+1),4,NCHD)                            591
      CALL ZHTOI(IQ(KQSP+LTOP+1),IQ(KQSP+LTOP+1),4)                       592
*
      NHPWD = 0                                                           594
      CALL VBLANK(IHPWD,2)                                                595
      CALL UCOPY(IHPWD,IQ(KQSP+LTOP+KPW1),2)                              596
      IQ(KQSP+LTOP+KPW1+2) = NCHD                                         597
      IF (IMODEX.GT.0) THEN                                               598
      CALL SBIT1(IQ(KQSP+LTOP+KPW1+2),12)                                 599
      ENDIF                                                               600
      IQ(KQSP+LTOP+KDATEC) = IDTIME                                       601
      IQ(KQSP+LTOP+KDATEM) = IDTIME                                       602
      IQ(KQSP+LTOP+KQUOTA) = NREC                                         603
      IQ(KQSP+LTOP+KRUSED) = NRD                                          604
      IQ(KQSP+LTOP+KWUSED) = NWL                                          605
C Check IOPTN.NE.0                                                        606
      IF (IOPTN.NE.0) THEN                                                606
      WRITE(IQLOG,7001) CHDIR                                             607
7001  FORMAT(' RZMAKE. new RZ format selected for ',A)                    608
*     +        ' This file will not be readable with versions',
*     +            ' of RZ prior to release 94B')
      IQ(KQSP+LTOP+KRZVER) = 1                                            611
      ELSE                                                                612
      WRITE(IQLOG,7007) CHDIR                                             613
7007  FORMAT(' *****  RZMAKE. OLD RZ format selected for ',A)             614
*     +        ' This file will have the limit on the number of',
*     +           ' blocks < 64 K')
      IQ(KQSP+LTOP+KRZVER) = 0                                            617
      ENDIF                                                               618
      IQ(KQSP+LTOP+KIP1) = 2                                              619
      IQ(KQSP+LTOP+KNFREE) = NWL-LF                                       620
      IQ(KQSP+LTOP+KLD) = LD                                              621
      IQ(KQSP+LTOP+KLB) = LB                                              622
      IQ(KQSP+LTOP+KLS) = LS                                              623
      IQ(KQSP+LTOP+KLK) = LK                                              624
      IQ(KQSP+LTOP+KLF) = LF                                              625
      IQ(KQSP+LTOP+KLC) = NWL+1                                           626
      IQ(KQSP+LTOP+KLE) = NWL                                             627
      IQ(KQSP+LTOP+KNWKEY) = NWKEY                                        628
      IQ(KQSP+LTOP+LD) = NRD                                              629
      IQ(KQSP+LTOP+LB) = NWREC                                            630
      IQ(KQSP+LTOP+LB+1) = LREC                                           631
      IQ(KQSP+LTOP+LB+2) = IDTIME                                         632
*
      NCHF=LEN(CHFORM)                                                    634
      NCH =LEN(CHTAG(1))                                                  635
      IF (NCH.GT.8) NCH=8                                                 636
      DO 20 I=1,NWKEY                                                     637
      IF (NCH.LT.8) CALL VBLANK(IHDIR,2)                                  638
      CALL UCTOH(CHTAG(I),IHDIR,4,NCH)                                    639
      CALL UCOPY(IHDIR,IQ(KQSP+LTOP+KTAGS+2*(I-1)),2)                     640
      IFORM=2                                                             641
      IF (I.LE.NCHF) THEN                                                 642
      IF (CHFORM(I:I).EQ.'B') IFORM=1                                     643
      IF (CHFORM(I:I).EQ.'H') IFORM=3                                     644
      IF (CHFORM(I:I).EQ.'A') IFORM=4                                     645
      ENDIF                                                               646
      IKDES=(I-1)/10                                                      647
      IKBIT1=3*I-30*IKDES-2                                               648
      CALL SBYT(IFORM,IQ(KQSP+LTOP+KKDES+IKDES),IKBIT1,3)                 649
20    CONTINUE                                                            650
      CALL ZHTOI(IQ(KQSP+LTOP+KTAGS),IQ(KQSP+LTOP+KTAGS),2*NWKEY)         651
      DO 30 I=1,NRD                                                       652
      IQ(KQSP+LTOP+LD+I)=I+1                                              653
      CALL SBIT1(IQ(KQSP+LTOP+LB+3),I+1)                                  654
30    CONTINUE                                                            655
*
*            Store default LOG level
*
      LOGL = LOGLV + 3                                                    659
      CALL SBYT(LOGL,IQ(KQSP+LTOP),15,3)                                  660
      CALL RZVCYC(LTOP)                                                   661
*
*            Allocate free records
*
      CALL MZBOOK(JQPDVS,LFREE,LTOP,-2,'RZFR',0,0,3,2,0)                  665
      IQ(KQSP+LFREE-5)=LUN                                                666
      IQ(KQSP+LFREE+1)=1                                                  667
      IQ(KQSP+LFREE+2)=NRD+2                                              668
      IQ(KQSP+LFREE+3)=NREC                                               669
*
*            Allocate space for used records
*
      CALL MZBOOK(JQPDVS,LUSED,LTOP,-3,'RZUS',0,0,21,2,0)                 673
*
      IQ(KQSP+LUSED-5)=LUN                                                675
      LRIN = 0                                                            676
      LPURG = 0                                                           677
      LROUT = 0                                                           678
      LCDIR = LTOP                                                        679
      NLCDIR= 1                                                           680
      NLNDIR= 1                                                           681
      NLPAT = 1                                                           682
      CHCDIR(1)=CHTOP                                                     683
      CHNDIR(1)=CHTOP                                                     684
      IQUEST(1)=0                                                         685
*
99    RETURN                                                              687
C     prevent "never used" warning
99999 I=JBIT(1,2)+JBYT(1,2,3)                                             689
      END                                                                 690
