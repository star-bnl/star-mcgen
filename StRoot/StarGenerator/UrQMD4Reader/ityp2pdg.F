cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function pdgid (ityp, iso3)
c
c     Revision : 1.0
c
coutput pdgid  : Particle-ID according to Particle Data Group  
c 
c     converts UrQMD-Id to PDG-Id 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.inc'

      integer pdgid
      integer ityp
      integer iso3

      integer tab_size
      parameter (TAB_SIZE = 181)

      logical anti
      integer abs_ityp
      integer norm_iso3
      integer idtab(3,TAB_SIZE)
      integer first
      integer last
      integer next
      integer isoit

      data idtab/
c Neutron
     .       1, -1,  2112,  
c Proton
     .       1,  1,  2212,
c N*
     .       2, -1, 12112,        2,  1, 12212,
     .       3, -1,  1214,        3,  1,  2124, 
     .       4, -1, 22112,        4,  1, 22212,
     .       5, -1, 32112,        5,  1, 32212,
     .       6, -1,  2116,        6,  1,  2216,
     .       7, -1, 12116,        7,  1, 12216,
     .       8, -1, 21214,        8,  1, 22124,
     .       9, -1, 42112,        9,  1, 42212, 
     .      10, -1, 31214,       10,  1, 32124, 
     .      14, -1,  1218,       14,  1,  2128, 
     .      23, -1,  1218,       23,  1,  2128,
c Delta
     .      24, -3,  1114,  24, -1,  2114,  24, 1,  2214,  24, 3,  2224,
     .      25, -3, 31114,  25, -1, 32114,  25, 1, 32214,  25, 3, 32224,
     .      26, -3,  1112,  26, -1,  1212,  26, 1,  2122,  26, 3,  2222,
     .      27, -3, 11114,  27, -1, 12114,  27, 1, 12214,  27, 3, 12224,
     .      28, -3, 11112,  28, -1, 11212,  28, 1, 12122,  28, 3, 12222,
     .      29, -3,  1116,  29, -1,  1216,  29, 1,  2126,  29, 3,  2226,
     .      30, -3, 21112,  30, -1, 21212,  30, 1, 22122,  30, 3, 22222,
     .      31, -3, 21114,  31, -1, 22114,  31, 1, 22214,  31, 3, 22224,
     .      32, -3, 11116,  32, -1, 11216,  32, 1, 12126,  32, 3, 12226,
     .      33, -3,  1118,  33, -1,  2118,  33, 1,  2218,  33, 3,  2228,
     .      40, -3,  1118,  40, -1,  2118,  40, 1,  2218,  40, 3,  2228,
c Lambda
     .      41,  0,  3122,
     .      42,  0, 13122,   
     .      43,  0,  3124,   
     .      44,  0, 23122,   
     .      45,  0, 33122,
     .      46,  0, 13124,
     .      47,  0, 43122,   
     .      48,  0, 53122,   
     .      49,  0,  3126,   
     .      50,  0, 13126,   
     .      51,  0, 23124,   
     .      52,  0,  3128,   
     .      53,  0, 23126,   
c Sigma
     .      54, -2,  3112,    54,  0,  3212,    54,  2,  3222,
     .      55, -2,  3114,    55,  0,  3214,    55,  2,  3224,
     .      56, -2, 13112,    56,  0, 13212,    56,  2, 13222,
     .      57, -2, 13114,    57,  0, 13214,    57,  2, 13224,
     .      58, -2, 23112,    58,  0, 23212,    58,  2, 23222,
     .      59, -2,  3116,    59,  0,  3216,    59,  2,  3226,
     .      60, -2, 13116,    60,  0, 13216,    60,  2, 13226,
     .      61, -2, 23114,    61,  0, 23214,    61,  2, 23224,
     .      62, -2,  3118,    62,  0,  3218,    62,  2,  3228,
c Xi
     .      63, -1,  3312,     63,  1,  3322,
     .      64, -1,  3314,     64,  1,  3324,
     .      66, -1, 13314,     66,  1, 13324,
c Omega
     .      69,  0,  3334,
c  Lambda_C
     .      70,  0,  4122,
c gamma
     .     100,  0,    22, 
c pion
     .     101, -2,  -211,    101,  0,   111,    101,  2,   211,
c eta
     .     102,  0,   221,
c omega
     .     103,  0,   223,
c rho
     .     104, -2,  -213,    104,  0,   113,    104,  2,   213,
c f0(980)
     .     105,  0, 10221,
c kaon
     .     106, -1,   311,    106,  1,   321,
c eta'
     .     107,  0,   331,
c k*(892)
     .     108, -1,   313,    108,  1,   323,
c phi
     .     109,  0,   333,
c k0*(1430)
     .     110, -1, 10311,    110,  1, 10321,
c a0(980)
     .     111, -2,-10211,    111,  0, 10111,    111,  2, 10211,
c f0(1370)
     .     112,  0, 20221,
c k1(1270)
     .     113, -1, 10313,    113,  1, 10323,
c a1(1260)
     .     114, -2,-20213,    114,  0, 20113,    114,  2, 20213,
c f1(1285)
     .     115,  0, 20223,
c f1'(1510)
     .     116,  0, 40223,
c k2*(1430)
     .     117, -1,   315,    117,  1,   325,
c a2(1329)
     .     118, -2,  -215,    118,  0,   115,    118,  2,   215,
c f2(1270)
     .     119,  0,   225,
c f2'(1525)
     .     120,  0,   335,
c k1(1400)
     .     121, -1, 20313,    121,  1, 20323,
c b1
     .     122, -2,-10213,    122,  0, 10113,    122,  2, 10213,
c h1
     .     123,  0, 10223,
c K* (1410)
     .     125, -1, 30313,    125,  1, 30323,
c rho (1450)
     .     126, -2,-40213,    126,  0, 40113,    126,  2, 40213,
c omega (1420)
     .     127,  0, 50223,
c phi(1680)
     .     128,  0, 10333,
c k*(1680)
     .     129, -1, 40313,    129,  1, 40323,
c rho(1700)
     .     130, -2,-30213,    130,  0, 30113,    130,  2, 30213,
c omega(1600)
     .     131,  0, 60223,
c phi(1850)     
     .     132,  0,   337,
c D
     .     133,  -1,   421,    133,   1,   411,  
c D*
     .     134,  -1, 10421,    134,   1, 10411,  
c J/Psi
     .     135,  0, 443,
c Chi_c
     .     136,  0, 10441,
c Psi'
     .     137,  0, 100443,
c Ds
     .     138,   0,   431,  
c Ds*
     .     139,   0,   433 /


c PYTHIA pdgid's (only need to subtract the offset!)
      if(abs(ityp).gt.1000) then
         if(ityp.gt.0) then
            pdgid=ityp-1000
         else
            pdgid=-1*(abs(ityp)-1000)
         endif
         return
      endif

cb check for antiparticles
      if (ityp.lt.0) then
cl its an antiparticle         
         anti = .true.
         abs_ityp = abs(ityp)
         norm_iso3 = -iso3
cl only mesons with odd isospin can have a negative ITYPE
         if ((abs_ityp.gt.minmes).and.
     .        (mod(isoit(abs_ityp),2).eq.0)) then
            call error ('pdgid','Negative ITYP not allowed',
     .           dble(ityp),3)
            pdgid = 0
            return
         endif
      else
         anti = .false.
         abs_ityp = ityp
         norm_iso3 = iso3
      endif

cb search for the ITYP in IDTAB

      first = 1
      last = TAB_SIZE
      if (idtab(1,first).eq.abs_ityp) then 
         next = first 
         goto 200
      endif
      if (idtab(1,last).eq.abs_ityp) then
         next = last 
         goto 200
      endif

 100  continue

cl ITYP not found in IDTAB
      if (last.le.(first+1)) then
         pdgid = 0
         return
      endif

      next = (first+last)/2

      if (idtab(1,next).lt.abs_ityp) then 
         first = next
         goto 100
      elseif (idtab(1,next).gt.abs_ityp) then
         last = next
         goto 100
      endif

 200  continue

cl calculate the entry with the wanted ISO3
      next = next - (idtab(2,next)-norm_iso3)/2
      
cl check if we found the correct values in IDTAB
      if ((idtab(1,next).eq.abs_ityp).and.
     .    (idtab(2,next).eq.norm_iso3)) then
         if (anti) then
            pdgid = -idtab(3,next)
         else
            pdgid = idtab(3,next)
         endif
      else
         call error ('pdgid','Error in tablelookup',dble(next),3)
         pdgid = 0
      endif
      
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function partname (ityp)
c
c     Revision : 1.0
c
coutput partname : Name of the Particle as a character string 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.inc'

      character*15 partname
      integer ityp

      integer abs_ityp
      character*7 baryon_names(maxbar-minbar+1) 
      character*11 meson_names(maxmes-minmes+1) 
      character*1 prefix

      data baryon_names/
     .     'Nukleon',
     .     'N(1440)',
     .     'N(1520)',
     .     'N(1535)',
     .     'N(1650)',
     .     'N(1675)',
     .     'N(1680)',
     .     'N(1700)',
     .     'N(1710)',
     .     'N(1720)',
     .     'N(1900)',
     .     'N(1990)',
     .     'N(2080)',
     .     'N(2190)',
     .     'N(2220)',
     .     'N(2250)',
     .     'N(2570)',
     .     'N(1875)',
     .     'N(2700)',
     .     'N(3100)',
     .     'N(3500)',
     .     'N(3800)',
     .     'N(4100)',
     .     'D(1232)',
     .     'D(1600)',
     .     'D(1620)',
     .     'D(1700)',
     .     'D(1900)',
     .     'D(1905)',
     .     'D(1910)',
     .     'D(1920)',
     .     'D(1930)',
     .     'D(1950)',
     .     'D(2420)',
     .     'D(2750)',
     .     'D(2950)',
     .     'D(3300)',
     .     'D(3500)',
     .     'D(3700)',
     .     'D(4100)',
     .     'Lambda',
     .     'L(1405)',
     .     'L(1520)',
     .     'L(1600)',
     .     'L(1670)',
     .     'L(1690)',
     .     'L(1800)',
     .     'L(1810)',
     .     'L(1820)',
     .     'L(1830)',
     .     'L(1890)',
     .     'L(2100)',
     .     'L(2110)',
     .     'Sigma',
     .     'S(1385)',
     .     'S(1660)',
     .     'S(1670)',
     .     'S(1750)',
     .     'S(1775)',
     .     'S(1915)',
     .     'S(1940)',
     .     'S(2030)',
     .     'Xi',
     .     'X(1530)',
     .     'X(1690)',
     .     'X(1820)',
     .     'X(1950)',
     .     'X(2030)',
     .     'Omega',
     .     'LambdaC'/
     

      data meson_names/
     .     'gamma',
     .     'pion',
     .     'eta',
     .     'omega',
     .     'rho',
     .     'f0(980)',
     .     'kaon',
     .     'eta''',
     .     'k*(892)',
     .     'phi',
     .     'k0(1430)',
     .     'a0(980)',
     .     'f0(1370)',
     .     'k1(1270)',
     .     'a1(1260)',
     .     'f1(1285)',
     .     'f1(1510)',
     .     'k2*(1430)',
     .     'a2(1329)',
     .     'f2(1270)',
     .     'f2''(1525)',
     .     'k1(1400)',
     .     'b1',
     .     'h1',
     .     'h1''',
     .     'k*(1410)',
     .     'rho(1450)',
     .     'omega(1420)',
     .     'phi(1680)',
     .     'k*(1680)',
     .     'rho(1700)',
     .     'omega(1600)',
     .     'phi3(1850)', 
     .     'D',
     .     'D*',
     .     'J/Psi',
     .     'Chi_c',
     .     'Psi Prime',
     .     'D_s',
     .     'D_s*'/

      abs_ityp = abs(ityp)

c    set the prefix for anti-particles      
      if (ityp.lt.0) then
         prefix = '*'
      else 
         prefix = ' '
      endif

      if ((abs_ityp.ge.minbar).and.(abs_ityp.le.maxbar)) then
         partname = prefix//baryon_names (abs_ityp)
      elseif ((abs_ityp.ge.minmes).and.(abs_ityp.le.maxmes)) then
         partname = prefix//meson_names (abs_ityp-99)
      else
         call error ('partname','ITYP out of range',dble(ityp),3)
         partname = '---'
      endif

      return
      end







