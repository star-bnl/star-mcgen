c-----------------------------------------------------------------------
c                  venus 4.12       06 / 09-apr-1993
c-----------------------------------------------------------------------
c
c
c  author:
c  -------
c
c  klaus werner
c
c  univ. heidelberg, inst. f. theoretische physik
c  philosophenweg 19, 6900 heidelberg, germany
c
c  werner@dhdmpi5.bitnet, werner@hobbit.mpi-hd.mpg.de, 28877::werner
c
c
c  venus:        
c  ------
c
c  venus is a monte carlo procedure to simulate hadronic interactions at
c  ultrarelativistic energies (hadron-hadron, hadron-nucleus, nucleus-
c  nucleus scattering), and also interactions involving leptons (e+e- 
c  astore_egz.Fastore_egz.Fastore_egz.Fannihilation, lepton-nucleon, lepton-nucleus scattering).
c  venus is based on gribov-regge theory (of multiple pomeron exchange)
c  and classical relativistic string dynamics. a detailed description can 
c  be found in the univ. heidelberg preprint hd-tvp-93-1 (270 pages), to
c  be published in physics reports.
c
c
c  distribution:
c  -------------
c
c  This code should not be distributed without notifying the author, in-
c  dicating what the code is going to be used for. depending on the type
c  of reaction and the kind of observables analysed,  the systematic un-
c  certainties of the venus simulations vary strongly,  and this should 
c  be clarified before using venus.
c
c
c  important features:
c  -------------------
c
c  covariant treatment of secondary interactions. each produced particle
c  is allowed to reinteract with other produced particles or with 
c  spectators. important for hadron-nucleus, nucleus-nucleus and 
c  lepton-nucleus scattering.
c  no final state interaction, if 'radiac' is set zero.
c
c  participation of antiquarks (in addition to quarks) in the
c  colour exchange mechanism to form strings.
c
c  possibility of diquark breakup, leading to multi-strings, like
c  a forward quark linked via two (!) strings to two backward
c  quarks (double-string). such strings fragment differently than
c  quark-diquark strings. in case of the double-string, the forward
c  quark will fragment via two breaks into a leading baryon.
c
c  sophisticated fragmentation procedure. since space-time
c  evolution is an important issue concerning final state inter-
c  actions, it is not enough to have a fragmentation model which
c  works, it should have the right space-time description! there-
c  for the field-feynman model of earlier versions (<3.00) has been
c  abandoned and replaced by a very powerful and appealing procedure
c  suggested by artru/mennessier.
c
c  very large resonance table, including for example all nucleon
c  resonances up to 2 gev. the decay table includes also kshort
c  and lambda decays. decay suppression is governed by the parameters
c  'ndecay' and 'ndecax'. Setting 'ndecay' to 1 suppresses all decays.
c  for further details search for 'decay suppression' in subr. ainitl
c
c  for version >= 4.01, we use gribovs cut-pomeron probabilities to
c  determine the weights for multi-colour-exchange.
c
c
c  file definitions
c  ----------------
c
c  the program unit 'subr. afiles' contains all file definitions.
c  the files are:
c
c     optns: input file (may be the terminal).
c     check: contains information and error messages.
c     histo: used for analysis (not described here).
c     data : collects all produced particles.
c
c
c  input to run venus
c  ------------------
c
c  to simulate a nucleon-nucleon, nucleon-nucleus, nucleus-nucleus
c  collision, the following input is required:
c
c     'had'
c     nevent ndisplay energy zproj aproj ztarg atarg i
c     'param1' value1
c     'param2' value2
c         ...
c     'parami' valuei
c     'stp'
c
c  where 'had' and 'stp' has to be taken litterally. the variables are:
c
c     nevent : number of events
c     ndisplay: a message occurs after ndisplay events being done
c     energy: energy per incident nucleon
c     zproj: number of projectile protons
c     aproj: number of projectile nucleons
c     ztarg: number of target protons
c     atarg: number of target nucleons
c     i: number of parameters or options to be changed
c
c  if i>0 one has to provide i lines each containing
c
c          'param' value
c
c  param is the name of the parameter or option, value is the new value.
c  a list of parameters/options and there default values can be found
c  in 'subr. ainitl'.
c
c  example1:
c
c  o + pb at 200 agev, 5000 evts, message after 100,
c  default parameter values:
c
c     'had'
c     5000 100    200.0    8 16 82 207    0
c     'stp'
c
c  example2:
c
c  p + pb at 200 gev, 5000 evts, no message ,
c  fragmentation mean transverse momentum changed to .450 gev:
c
c     'had'
c     5000 10000    200.0    1 1 82 207    1
c     'ptf   ' .450
c     'stp'
c
c
c  input for hadron-nucleus for arbitrary hadrons
c  ----------------------------------------------
c
c  a hadron (arbitrary) nucleus collision requires:
c
c     'had'
c     nevent ndisplay energy -1 1 ztarg atarg i
c     idptl
c     'param1' value1
c     'param2' value2
c         ...
c     'parami' valuei
c     'stp'
c
c  where idptl is the particle id of the hadron (id-description later).
c  for meson-nucleon scattering the diffractive probabilities 'wproj'
c  and 'wtarg' have to be changed.
c  warning: meson-nucleon and (even more) meson-nucleus
c  scattering has been tested much less than purely nucleonic reactions,
c  so the parameters are not neccessarily optimized, and the author
c  would appreciate (in any case) comments and suggestions!
c
c  example3:
c
c  pion + p at 200 gev, 5000 evts, no message, 'wproj','wtarg'
c  changed:
c
c     'had'
c     5000 10000    200.0   -1 1 1 1       2
c     120
c     'wproj' 0.04
c     'wtarg' 0.04
c     'stp'
c
c  example4:
c
c  pion + au at 200 gev, 5000 evts, no message, 'wproj','wtarg'
c  changed:
c
c     'had'
c     5000 10000    200.0   -1 1 79 197      3
c     120
c     'wproj' 0.04
c     'wtarg' 0.04
c     'stp'
c
c  (120 is the id-code for pion+).
c
c
c  input for interactions involving leptons
c  ----------------------------------------
c
c  contact the author (klaus werner).
c
c
c  output
c  ------
c
c  the results are stored in the file 'ifdt', defined in 'subr.
c  afiles'. the storage is performed in 'subr. astore', which is
c  written in a very tranparent way to make changes by users
c  very easy. in order to learn what kind of variables are stored in
c  which way, the user is referred to 'subr. astore' (subroutines
c  are written in alphabetic order).
c  
c
c  notes
c  -----
c
c
c     *** the following dimensions are used (def with parameter statemt)
c            max # of strings per event:   mxstr=20000
c            max # of particles per event: mxptl=150000
c
c
c  updates
c  -------
c
c  venus 2.00  string flip involving antiquarks.
c  venus 3.00  fragmentation according to classical string theory.
c              extended resonance table.
c              rescattering of produced particles among themselves.
c              rescattering of produced particles with spectators.
c  venus 3.01  removing a slight impact-parameter trigger.
c  venus 3.02  adding analysis routines 'wimi', 'wimino'.
c  venus 3.04  adding phase space factor for qq-qqbar production.
c              removing an impact parameter trigger at ags energies.
c              no rescattering for iorptl=-1.
c              removing 'tmefrz' parameter, adding 'ndecay' instead.
c              adding parameter 'jorptl'. changing dimension mxptl.
c  venus 3.05  also rescattering if iorptl=0.
c              removing 'tmeiac' parameter, adding 'tmefor' instead.
c              separating stmt numbers and if stmts.
c  venus 3.06  lorentz contraction of rescattering distance.
c              radiac+radias as interaction radius.
c              meaning of amsiac modified (amsiac->m+amsiac).
c              multiple rescattering per time step included.
c              cluster decay acc. to phase space.
c              option 'bminim' added.
c  venus 3.07  option 'labsys' added.
c              analysis routines 'wspa' and 'wspafi' added.
c              nucleon hard core properly treated.
c              string origin as average collision location.
c  venus 3.10  rho0 -> rhophi*rho0 + (1-rhophi)*omega .
c              new sr idquac, wclu, wclufi.
c              changes to allow lhc energies (26m gev),
c                 new sr utlob2, utrot2, utchm,
c                 new sr utclea to clean /cptl/ .
c              covariant rescattering:
c                 taus (hyperbola and tangent) replaces t.
c                 new sr/fu uttain, uttaus.
c              covariant treatment in sr jamfra. new sr uttaix.
c              iaaptl=0 for parents of interact. ptls.
c              wspa: iacptl=0 for parents and childr of counted ptls.
c              new sr wspa, wspafi: defs acc kinetic theory.
c              no mass limit'n for cluster i.a. any more.
c              new parameter pcomp.
c              sr jclude considerably modified.
c              nrptla and nstptl removed from /cptl/ .
c              ibptl considered in elastic scattering.
c              structure functions (st...) reorganized.
c              new features regarding muon-nucleus scattering.
c                 new sr/fu: lep..., stx... (from g. jancso).
c              call idtr4 only within idquac.
c              idquac: special treatment for !id!<100.
c              hmucex changed and renamed racpro.
c              hard component for quark-pt (fctn sptq).
c  venus 3.12  ish redefined.
c              no string mass requirement any more at string formatn.
c              no en adjustment any more in sr hasi, only pt is chgd.
c              const sea probability possible (parameter prosea).
c              ish=16 prints calculated sea prob.
c              problems involving >10 quarks of the same flavour fixed.
c              update of subr. pvjpsf.
c              sr jamfra: split string: split not poss for ns=1 -> bl
c              sr wspa(fi): new option iwcent.
c              seedc determined in sr astore (so far in sr avenus).
c              qmust and ptq log energy dependent as default.
c              exact proton and neutron number (not average any more)
c              new option ioctet for octet-octet fragmentation.
c                 ioctet=2 default.
c              change in idtr4 to make ic=1,1 and 2,2 equally probabl
c  venus 3.14  new cluster decay acc to phase space (new sr jclude).
c  venus 3.15  leading strings treated as multistrings.
c                 changes in sr's hahabs, hastfs, hastfl, hastpr.
c                 new sr jsplit.
c                 new fn utamnu (formerly part of fn sjcgam).
c  venus 4.01  sr racpro: gribovs cut-pomeron probabilities.
c  venus 4.02  sr wjpsi, wjpsif for jpsi analysis added.
c              sr wspa: some storage-saving fixes (also in jfrade).
c              sr ainitl: energy dependent phard (as default)
c  venus 4.10  sr jamfra: split into jamfra and jambr1
c              alternative to jambr1: jambr2 = SAMBA string breaking proc
c              new general analysis sr's wzzzz, wzzzzf
c              low mass q-qbar str: fusion, pt transfer, only sea
c              new value w=.32
c              diffractive val + sea
c              phard=.30+12(lg s - lg 30^2)
c              new value pth=1.0
c              hard scatt: e > pt
c              new value ptq1=.26
c              new value iopt=2
c              new value cutmss=0.001
c              new parameter cutmsq=2.0
c              new valence quark count
c              all legs absorbed: random excitation
c              x>1-.85/s: redo h-h
c              x>1-10/s and only qqbar and no q-qbar: cluster
c              rescale strings from multistring
c              absorption of small mass legs
c              covariant string formation
c              new value rstras=0
c              sr idmimj removed
c              random formation time
c              arbitrary large entropy possible
c              high energy rescattering: colour exchange
c                 new parameter delrap=1.5
c              new value iopent=5
c              new parameter thaged=.23
c              new parameter kentro=100000 (repl kzentr)
c              new parameter uentro=4
c              new parameter ebamin=1.5
c              split (covariant) of leading particle
c              epart=elead/nsccx
c              ssplit(x)=z^3-z^5, z=(x+1)/2
c              remnant (after cutting legs): 33% excited, 33% s=3/2 acc
c              small mass str dec
c              new value pdiqua=.12
c              pdiq=0
c              new value raphdi=0
c              sr idresi: rho0=omega, rho0=rho+, no a,f, rho -> .950
c              new value parea=.60
c              pudx=pud*0.945
c              sr hastfr and sr jsplit more accurate (double precision)
c              sr jamfra: changes reg small mass str and backw +jets
c              sr veanly: selection of nondiffractive events possible
c              sr jclude: changes to avoid redo (infinite loops possible)
c  venus 4.11  sr jclude: modifications concerning redo cluster decay     
c              sr jfrade: split into jfrade and jintac   
c              sr jintac: split into jintac, jintcc, jintfs, jintce,
c                  jintch, jintpa ,jintel, jintfu
c              sr jintce: no iaaptl()=2
c  venus 4.12  fctn uttaus -> sr uttaus
c              rename sr jintac  to  sr jinta1
c              new sr jinta2 as alternative to jinta1 (option iojint)
c                default: iojint=2 --> call jinta2
c                iojint=1 --> call jinta1 --> equiv 4.10+ 
c                      (4.10+ = 4.10 + redo improvements) 
c              new sr jintfp
c              new sr jintcl
c
c  list of subroutines and functions (incl 'program aamain')
c  ---------------------------------------------------------
c
c     aamain  afiles  ainitl  astore  atermi  avenus
c     hahabs  hapapa  hasi    hastfc  hastfl  hastfr  hastfs
c     hastpr  hdecmp  hrescl
c     idchrg  idcomi  idcomj  idcomp  idcon   iddeco
c     idenco  idflav  idlabl  idmass  
c     idmix   idquac  idres   idresi  idsgl   idspin
c     idtau   idtr4   idtra   idtrb
c     jambr1  jambr2  jamfra  jcentd  jcentp  jcentr  jclude  jdeca
c     jdecay  jdecin  jestpr  jetgen  jfrade  jinta1  jinta2  
c     jintcc  jintce  jintch  jintcl  jintel  jintfp  jintfs  
c     jintfu  jintpa  jrescl
c     jsplit  jstrij  jstrip
c     lepexp  lepstr  leptar
c     nucini  nuclco  nucoge  nucoll  nucstr
c     pvjpsf
c     racpro  randfl  rndm    ranfgt  rangen  ranhvy
c     ransig  ranstc  ranxq
c     sbet    sdensi  sgam    sgau    sjcen0  sjcen4  sjcent  sjcgam
c     sjctem  smass   smassi  smassp  smasst  smassu
c     spoc    sptf    spth    sptj    sptq    ssplit  ssplix
c     staa    stgdv   stgus   stguv   stodpr  stoglu
c     stospi  stospr  stovpi  stovpr  stxd    stxs
c     stxu    stxus   stxzne  stxzpr  stzzse  stzzvl
c     utacos  utamnu  utamnx  utamny  utamst  utamsy  
c     utar    utaxis  utchm
c     utclea  utfost  uthist  uthsea  utinvt  utkprj  utksix
c     utkstr  utlobo  utloc   utlow   utlow3  utlow6
c     utmom   utovel  utpage  utpart  utpcm   utquad
c     utquaf  utremb  utrepl  utrota  utsqsi  utsqst
c     utsstr  utstop  uttain  uttaix  uttaus  uttima
c     uttucl
c     veanif  veanly
c     wclu    wclufi  wimi    wimino  wjpsi   wjpsif
c     wspa    wspafi  wzzzz   wzzzzf
c
c
c  particle id (same as isajet - see isajet.doc, f.e.paige)
c  --------------------------------------------------------
c
c  quarks and leptons are numbered in order of mass:
c        up=1 dn=2 st=3 ch=4 bt=5 tp=6
c        nue=11 e-=12 num=13 mu-=14 nut=15 tau-=16
c  with a negative sign for antiparticles.
c  arbitrary conventions are:
c        gl=9 gm=10 ks=20 kl=-20 w+=80 z0=90 .
c  the code for a meson is an integer +-jkl, where j.le.k are the
c  quarks and l is the spin. the sign is for the j quark.
c  flavor singlet mesons are ordered by mass:
c        pi0=110 eta=220 etap=330 etac=440 .
c  similarly, the code for a baryon is a compound integer +-ijkl
c  formed from the three quarks i,j,k and a spin label l=0,1.
c  the code for a diquark is +-ij00.
c
c  list of ident codes:
c
c     ident     label           mass    charge
c         1     UP        .30000E+00       .67
c        -1     UB        .30000E+00      -.67
c         2     DN        .30000E+00      -.33
c        -2     DB        .30000E+00       .33
c         3     ST        .50000E+00      -.33
c        -3     SB        .50000E+00       .33
c         4     CH        .16000E+01       .67
c        -4     CB        .16000E+01      -.67
c         5     BT        .49000E+01      -.33
c        -5     BB        .49000E+01       .33
c         6     TP        .30000E+02       .67
c        -6     TB        .30000E+02      -.67
c
c         9     GL       0.               0.00
c
c        10     GM       0.               0.00
c
c        11     NUE      0.               0.00
c       -11     ANUE     0.               0.00
c        12     E-        .51100E-03     -1.00
c       -12     E+        .51100E-03      1.00
c        13     NUM      0.               0.00
c       -13     ANUM     0.               0.00
c        14     MU-       .10566E+00     -1.00
c       -14     MU+       .10566E+00      1.00
c        15     NUT      0.               0.00
c       -15     ANUT     0.               0.00
c        16     TAU-      .18070E+01     -1.00
c       -16     TAU+      .18070E+01      1.00
c
c        20     KS        .49767E+00      0.00
c       -20     KL        .49767E+00      0.00
c
c        80     W+        SIN2W=.23       1.00
c       -80     W-        SIN2W=.23      -1.00
c        90     Z0        SIN2W=.23       0.00
c
c       110     PI0       .13496E+00      0.00
c       120     PI+       .13957E+00      1.00
c      -120     PI-       .13957E+00     -1.00
c       220     ETA       .54880E+00      0.00
c       130     K+        .49367E+00      1.00
c      -130     K-        .49367E+00     -1.00
c       230     K0        .49767E+00      0.00
c      -230     AK0       .49767E+00      0.00
c       330     ETAP      .95760E+00      0.00
c       140     AD0       .18633E+01      0.00
c      -140     D0        .18633E+01      0.00
c       240     D-        .18683E+01     -1.00
c      -240     D+        .18683E+01      1.00
c       340     F-        .20300E+01     -1.00
c      -340     F+        .20300E+01      1.00
c       440     ETAC      .29760E+01      0.00
c
c       111     RHO0      .77000E+00      0.00
c       121     RHO+      .77000E+00      1.00
c      -121     RHO-      .77000E+00     -1.00
c       221     OMEG      .78260E+00      0.00
c       131     K*+       .88810E+00      1.00
c      -131     K*-       .88810E+00     -1.00
c       231     K*0       .89220E+00      0.00
c      -231     AK*0      .89220E+00      0.00
c       331     PHI       .10196E+01      0.00
c       141     AD*0      .20060E+01      0.00
c      -141     D*0       .20060E+01      0.00
c       241     D*-       .20086E+01     -1.00
c      -241     D*+       .20086E+01      1.00
c       341     F*-       .21400E+01     -1.00
c      -341     F*+       .21400E+01      1.00
c       441     JPSI      .30970E+01      0.00
c
c      1120     P         .93828E+00      1.00
c     -1120     AP        .93828E+00     -1.00
c      1220     N         .93957E+00      0.00
c     -1220     AN        .93957E+00      0.00
c      1130     S+        .11894E+01      1.00
c     -1130     AS-       .11894E+01     -1.00
c      1230     S0        .11925E+01      0.00
c     -1230     AS0       .11925E+01      0.00
c      2130     L         .11156E+01      0.00
c     -2130     AL        .11156E+01      0.00
c      2230     S-        .11974E+01     -1.00
c     -2230     AS+       .11974E+01      1.00
c      1330     XI0       .13149E+01      0.00
c     -1330     AXI0      .13149E+01      0.00
c      2330     XI-       .13213E+01     -1.00
c     -2330     AXI+      .13213E+01      1.00
c      1140     SC++      .24300E+01      2.00
c     -1140     ASC--     .24300E+01     -2.00
c      1240     SC+       .24300E+01      1.00
c     -1240     ASC-      .24300E+01     -1.00
c      2140     LC+       .22600E+01      1.00
c     -2140     ALC-      .22600E+01     -1.00
c      2240     SC0       .24300E+01      0.00
c     -2240     ASC0      .24300E+01      0.00
c      1340     USC.      .25000E+01      1.00
c     -1340     AUSC.     .25000E+01     -1.00
c      3140     SUC.      .24000E+01      1.00
c     -3140     ASUC.     .24000E+01     -1.00
c      2340     DSC.      .25000E+01      0.00
c     -2340     ADSC.     .25000E+01      0.00
c      3240     SDC.      .24000E+01      0.00
c     -3240     ASDC.     .24000E+01      0.00
c      3340     SSC.      .26000E+01      0.00
c     -3340     ASSC.     .26000E+01      0.00
c      1440     UCC.      .35500E+01      2.00
c     -1440     AUCC.     .35500E+01     -2.00
c      2440     DCC.      .35500E+01      1.00
c     -2440     ADCC.     .35500E+01     -1.00
c      3440     SCC.      .37000E+01      1.00
c     -3440     ASCC.     .37000E+01     -1.00
c
c      1111     DL++      .12320E+01      2.00
c     -1111     ADL--     .12320E+01     -2.00
c      1121     DL+       .12320E+01      1.00
c     -1121     ADL-      .12320E+01     -1.00
c      1221     DL0       .12320E+01      0.00
c     -1221     ADL0      .12320E+01      0.00
c      2221     DL-       .12320E+01     -1.00
c     -2221     ADL+      .12320E+01      1.00
c      1131     S*+       .13823E+01      1.00
c     -1131     AS*-      .13823E+01     -1.00
c      1231     S*0       .13820E+01      0.00
c     -1231     AS*0      .13820E+01      0.00
c      2231     S*-       .13875E+01     -1.00
c     -2231     AS*+      .13875E+01      1.00
c      1331     XI*0      .15318E+01      0.00
c     -1331     AXI*0     .15318E+01      0.00
c      2331     XI*-      .15350E+01     -1.00
c     -2331     AXI*+     .15350E+01      1.00
c      3331     OM-       .16722E+01     -1.00
c     -3331     AOM+      .16722E+01      1.00
c      1141     UUC*      .26300E+01      2.00
c     -1141     AUUC*     .26300E+01     -2.00
c      1241     UDC*      .26300E+01      1.00
c     -1241     AUDC*     .26300E+01     -1.00
c      2241     DDC*      .26300E+01      0.00
c     -2241     ADDC*     .26300E+01      0.00
c      1341     USC*      .27000E+01      1.00
c     -1341     AUSC*     .27000E+01     -1.00
c      2341     DSC*      .27000E+01      0.00
c     -2341     ADSC*     .27000E+01      0.00
c      3341     SSC*      .28000E+01      0.00
c     -3341     ASSC*     .28000E+01      0.00
c      1441     UCC*      .37500E+01      2.00
c     -1441     AUCC*     .37500E+01     -2.00
c      2441     DCC*      .37500E+01      1.00
c     -2441     ADCC*     .37500E+01     -1.00
c      3441     SCC*      .39000E+01      1.00
c     -3441     ASCC*     .39000E+01     -1.00
c      4441     CCC*      .48000E+01      2.00
c     -4441     ACCC*     .48000E+01     -2.00
 
c-----------------------------------------------------------------------
      subroutine afiles
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      character c*8
      ifop=5
      ifmt=6
      ifch=31
      ifhi=35
      ifdt=51
      open(unit=ifop,file='optns.ven',status='old')
      return
      end
 
c-----------------------------------------------------------------------
      subroutine ainitl(jread)
c     performs all initializations
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      common/paro4/grigam,grirsq,gridel,grislo,gricel
      common/paro5/delvol,deleps
      common/nevnt/nevnt
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter(nprbms=20)
      common /cprbms/prbms(nprbms)
      parameter (nptj=129)
      common /cptj/xptj(nptj),qptj(nptj),wptj(nptj)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      parameter (nptf=129)
      common /cptf/xptf(nptf),wptf(nptf)
     *,qptfu(nptf),qptfs(nptf),qptfuu(nptf),qptfus(nptf),qptfss(nptf)
      parameter (ngau=129)
      common /cgau/xgau(ngau),qgau(ngau),wgau(ngau)
      parameter (nsplit=129)
      common /csplit/xsplit(nsplit),qsplit(nsplit),qsplix(nsplit)
     *,wsplit(nsplit)
      parameter (nstru=2049)
      common/stru/xstru(nstru),wstru(nstru)
     *,qvap(nstru),qvat(nstru),qsep(nstru),qset(nstru)
      character choic*3
      parameter(mxpar=99)
      character*6 parcha(mxpar)
      double precision parval(mxpar)
      common/paramc/npar,idmmy,parval,parcha
      common /cvsn/iversn /cden/massnr,rmx
      external sdensi,stzzvl,stzzse
     *,sptj,spth,sptq,sptf,sgau,ssplit,ssplix
      common /cptlu/nptlu
      common /cquama/quama /cipio/ipio /clep/icinpu,idscat
      parameter (mxtau=4,mxvol=10,mxeps=20)
      common/cjintc/clust(mxtau,mxvol,mxeps)
      common/outfgdt/ifgdt
c  ifgdt=1 text output for input to nuclear fragmentis program
      COMMON/HEADPSS/IFIRST,PSSHEP(5),VSSHEP(4)
      INTEGER IFIRST
      DOUBLE PRECISION PSSHEP,VSSHEP
      SAVE /HEADPSS/
c     version number
c     --------------
         iversn=412
 
c     fragmentation parameters/options
c     --------------------------------
c     prob. for u or d quark production ( =(1-p_strange)/2 ):
         pud=0.455 
c     qq-qqbar probability
         pdiqua=.12
c     spin probabilities (for light and heavy flavours):
         pspinl=.50
         pspinh=.75
c     isospin probability:
         pispn=.50
c     option for p_t distribution (1=exponential,2=Gaussian):
         ioptf=1
c     average p_transverse
         ptf=.40
c     string tension:
         tensn=1.0
c     string decay parameter
         parea=.60
c     threshold resonance -> string
         delrem=1.
c     cutoff for kmaxor beyond which pdiq=0 in sr jsplit
         kutdiq=4
c     option for breaking procedure (1=amor,2=samba)
         iopbrk=1
 
c     proton-proton parameters/options
c     --------------------------------
c     option for quark p_t distribution (1=exponential,2=Gaussian,3=powe
         ioptq=2
c     mean transverse momentum of quarks
c      (q1+q2*ln(e)+q3*ln(e)**2, e=sqrt(s)):
         ptq1=0.260
         ptq2=0.
         ptq3=0.
c     probability for semihard interaction (not used if negative):
         phard=-1
c     cutoff parameter for p_t distr. for semihard interactions:
         pth=1.0
c     effective ratio of strange sea over u sea:
         rstras=0.
c     diffractive probability:
         wproj=.32
         wtarg=.32
c     effective cutoff mass in structure functions:
         cutmsq=2.0    ! for sea-valence ratio
         cutmss=.001
c     valence quark fraction in case of diffractive interaction
         pvalen=0.30
c     phase space parameters:
         delmss=.300
 
c     gribov-regge-theory parameters
c     ------------------------------
c     gamma (in fm**2):
         grigam=3.64*0.04
c     r**2(in fm**2):
         grirsq=3.56*0.04
c     delta=intercept of regge trajectory-1:
         gridel=0.07
c     slope of regge trajectory (in fm**2):
         grislo=0.25*0.04
c     c (determines relative weight of elastic and diffr cross sctn):
         gricel=1.5
 
c     nucleus-nucleus parameters
c     --------------------------
c     hard core distance:
         core=0.8
c     jpsi nucleon cross section [fm**2]:
         sigj=0.2
 
c     rescattering parameters
c     -----------------------
c     reaction time:
         taurea=1.5
c     overlap parameter (not used if negative)
         overlp=-1.
c     baryon radius:
         radiac=0.65
c     meson radius:
         radias=0.35
c     interaction mass:
         amsiac=.8
c     option to call jinta1 (1) or jinta2 (2)
         iojint=1
c     print options
         amprif=0
         delvol=1
         deleps=1

c     cluster decay parameters/options
c     --------------------------------
c     option for entropy calculation:
c      iopent=0: zero entropy
c      iopent=1: oscillator model (0 for k.le.uentro)
c      iopent=2: fermi gas with const volume (0 for k.le.uentro)
c      iopent=3: fermi gas with const density (0 for k.le.uentro)
c      iopent=4: fermi gas with const volume - new (0 for k.le.uentro)
c      iopent=5: resonance gas (hagedorn) (0 for u.le.uentro)
         iopent=5
         uentro=4.
         kentro=100000
c     decay time (comoving frame):
         taunll=1.0
c     oscillator quantum
         omega=.500
 
c     presently not used
c     ------------------
 
c     cluster decay initializations
c     -----------------------------
c     average hadron masses, two lowest multiplets (if possible):
c     n/delta,lambda/sigma,xi,omega,pi/rho,kaon,delta:
         asuhax(1)=1.134
         asuhax(2)=1.301
         asuhax(3)=1.461
         asuhax(4)=1.673
         asuhax(5)=0.6125
         asuhax(6)=0.7915
         asuhax(7)=1.2320
c     lowest masses:
         asuhay(1)=0.940
         asuhay(2)=1.200
         asuhay(3)=1.322
         asuhay(4)=1.673
         asuhay(5)=0.1400
         asuhay(6)=0.4977
         asuhay(7)=1.2320
 
c     technical parameters
c     --------------------
c     min tau for space-time evolution:
         taumin=0
c     max tau for space-time evolution
         taumax=10
c     tau steps for space-time evotution
         numtau=51
c     range for pt distribution
         ptmx=6
c     range for gauss distribution
         gaumx=8
c     parameter determining range for density distribution
         fctrmx=10
c     try-again parameter
         ntrymx=10
c     max time for jpsi evolution
         taumx=20
c     time steps for jpsi evolution
         nsttau=100
 
c     options
c     -------
c     option for minimum energy in sjcgam:
c       iopenu = 1 : sum of hadron masses
c       iopenu = 2 : bag model curve with minimum at nonzero strangen.
         iopenu=1
c     parameter theta in berger/jaffe mass formula
         themas=.51225
c     sea probability (if .lt. 0. then calculated from structure fncts)
         prosea=-1.
c     inelastic pp cross section [fm**2]
c      (if negative: calculated from gribov-regge-theory):
         sigppi=-1
c     multistring parameter (q1+q2*ln(e)+q3*ln(e)**2, e=sqrt(s)):
c      (not used if racpro is called with 'gri'-option (default))
         qmust1=0.50
         qmust2=0.
         qmust3=0.
c     entro() calculated (1) or from data (2)
         ientro=2
c     dual parton model (1) or not (else)
         idpm=0
c     antiquark color exchange (1) or not (0):
         iaqu=1
c     minimum number of valence quarks:
         neqmn=-5
c     maximum number of valence quarks:
         neqmx=5
c     upper limit for rapidity interval for intermittency analysis
         ymximi=2
c     clean /cptl/ if nclean > 0 (every nclean_th time step)
         nclean=0
c     trafo from pp-cm into lab-system (1) or not (.ne.1)
         labsys=1
c     maximum number of collisions:
         ncolmx=100000
c     maximum resonance spin (spin in a genaral sense: mod(/id/,10))
         maxres=99999
c     momentum rescaling (1=yes):
         irescl=1
c     nue energy
         elepti=43.00
c     mue  energy
         elepto=26.24
c     mue angle
         angmue=3.9645/360.*2*3.14159
c     jpsi to be produced (1) or not (0):
         jpsi=0
c     jpsi final state interaction (1) or not (0):
         jpsifi=0
c     collision trigger (only coll between ko1 and ko2 are used):
         ko1ko2=00009999
c     print option:
c     ish=14: call uttima
c     ish=15: prints ptls read from data file in sr veanly
c     ish=16: prints sea prob.
c     ish=17: prints random numbers
c     ish=18: sr jclude, no-phase-space clusters
c     ish=19: sr ainitl, call smassp
c     ish=20: sr veanly, prints evt nr if evt is accepted
c     ish=21: creates histogram for sea distribution
c     ish=22: sr jfrade, msg after call utclea
c     ish=23: call jintfp
c     ish=24: call jintcl
c     ish=90,91,92,93,94,95: more and more detailed messages.
         ish=0
c     print option:
c     ishsub=ijmn, ij specifies location where ish=mn.
c     ij=01: sr jclude
c     ij=02: sr jetgen
c     ij=03: sr jfrade, starting before fragmentation
c     ij=04: sr jdecay
c     ij=05: sr jdeca
c     ij=06: sr nucoll
c     ij=07: sr nucoge+-
c     ij=08: sr astore
c     ij=09: sr jfrade, starting after fragmentation
c     ij=10: sr jfrade, starting before decay
c     ij=11: sr jfrade, starting after interactions
c     ij=12: sr jcentr, entro() in data format
c     ij=13: sr jcentp
c     ij=14: sr jdeca if cluster decay
c     ij=15: sr jsplit
c     ij=16: sr jfrade
c     ij=17: sr racpro
c     ij=18: sr utclea
c     ij=19: sr jinta1, jinta2, after call utclea
         ishsub=0
c     print option:
c     if ishevt.ne.0: for evt#.ne.ishevt ish is set to 0
         ishevt=0
c     print marks between whom ish is set to ish(init):
         ipagi=0
c     verify option for input reading:
         ivi=1
c     maximum impact parameter (bmaxim=0=>central):
         bmaxim=10000.
c     minimum impact parameter:
         bminim=0.
c     store only stable ptl (0) or also parents (1):
         istmax=0
c     suppression (1) or not of messages
         isup=0
c     suppression of calling jfrade (0). jfrade=fragm+decay+rescattering
         ifrade=1
c..   decay suppression. ndecay specifies which resonances are not decay
c..   0000001 : all resonances
c..   0000010 : k_short/long (+-20)
c..   0000100 : lambda (+-2130)
c..   0001000 : sigma (+-1130,+-2230)
c..   0010000 : cascade (+-2330,+-1330)
c..   0100000 : omega (+-3331)
c..   1000000 : pi0 (110)
         ndecay=0
c..   decay suppression. ndecax specifies which resonances are not decay
c..   0000001 : jpsi 	
c..   0000010 : k_zero (+-230)
c..   0000100 : delta (+-1111,+-1121,+-1221,+-2221)
c..   0001000 : rho,omega,phi (111,+-121,221,331)
c..   0010000 : eta (220)
c..   0100000 : etaprime (330)
c..   1000000 : a0 (112), a+- (+-122)
         ndecax=0
c..   decay suppression. ndecaw specifies which resonances are not decay
c..   0000001 : f0 (332)
c..   0000010 : k* (+-131,+-231)
         ndecaw=0
c     fill zzzz histograms (1) or not (0)
         iwzzzz=0
c     fill intermittency histograms (1) or not (0)
         imihis=0
c     fill space-time histograms (1) or not (0)
         isphis=0
c     fill cluster histograms (1) or not (0)
         iclhis=0
c     fill jpsi histograms (1) or not (0)
         ijphis=0
c     rho/rho+phi ratio
         rhophi=.5
c     wspa: all ptls (1) or only interacting ptls (else)
         ispall=1
c     tmin in wspa
         wtmini=-3
c     t-step in wspa
         wtstep=1
c     only central point (1) or longitudinal distribution (else) in wspa
         iwcent=0
c     quark masses
         smas=0.
         uumas=0.
         usmas=0.
         ssmas=0.
 
c     constants
c     ---------
      pi=3.1415927
      prom=.94
      piom=.14
      ainfin=1e30
 
c     initializations
c     ---------------
      laproj=0
      maproj=0
      latarg=0
      maproj=0
      idproj=1120
      idtarg=1120
      do 6 i=1,99
      prob(i)=0
      do 6 j=1,2
      icbac(i,j)=0
6     icfor(i,j)=0
      pnll=0
      do 2 i=1,mxpar
      parval(i)=0
2     parcha(i)='      '
      imsg=0
      jerr=0
      ntevt=0
      nrevt=0
      naevt=0
      nrstr=0
      nrptl=0
      inoiac=0
      ilamas=0
      nptlu=0
      do itau=1,mxtau
      do ivol=1,mxvol
      do ieps=1,mxeps
      clust(itau,ivol,ieps)=0
      enddo
      enddo
      enddo
 
c     read from ifop
c     ---------------
      if(jread.ne.1)return
      if(isup.ne.1)write(ifmt,105)float(iversn)/100
105   format(
     * ' !-----------------------------------------------------!'
     */' !            venus',f5.2,6x,'-     k. werner          !'
     */' !-----------------------------------------------------!')
1     continue
      if(isup.ne.1)write(ifmt,*)'choice? [a3]      type hlp for help'
      read(ifop,*,err=1)choic
      call utlow3(choic)
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)choic
      if(choic.eq.'hlp')write(ifmt,100)
      if(choic.eq.'hlp')goto1
      if(choic.eq.'sup')then
      isup=1
      goto1
      endif
100   format(
     * ' !----------------------------------------------------!'
     */' !  choice [format:character*3]                       !'
     */' !   str ... string fragmentation                     !'
     */' !   had ... hadronic collisions (hh,ha,aa)           !'
     */' !   col ... only collision analysis                  !'
     */' !   ans ... only analysis (string frg.)              !'
     */' !   anh ... only analysis (hadr. coll.)              !'
     */' !----------------------------------------------------!')
      if(choic.eq.'str')ichoic=1
      if(choic.eq.'had')ichoic=2
      if(choic.eq.'col')ichoic=3
      if(choic.eq.'ans')ichoic=4
      if(choic.eq.'anh')ichoic=5
 
      if(ichoic.eq.1.or.ichoic.eq.4)then
      if(isup.ne.1)write(ifmt,*)'#events?[i]   modsho?[i]   '
     *,'engy?[f]   #str?[i])   tg?[2i]   #par?[i]'
      read(ifop,*)nevent,modsho,engy,ns,latarg,matarg,npar
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,103)
     *nevent,modsho,engy,ns,latarg,matarg,npar
103   format(1x,i6,i5,f9.3,3x,i6,3x,2i6,3x,i6)
      if(ns.eq.0)stop'nr of strings must not be 0.'
      if(ns.gt.99)stop'nr of strings has to be < 100.'
      if(ns.gt.0)then
      icinpu=0
      do 4 n=1,ns
      if(isup.ne.1)write(ifmt,*)n,'. string?[fiiii]'
      read(ifop,*)prob(n)
     *,icbac(n,1),icbac(n,2),icfor(n,1),icfor(n,2)
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,106)prob(n)
     *,icbac(n,1),icbac(n,2),icfor(n,1),icfor(n,2)
4     continue
106   format(9x,f5.2,4x,2i8,4x,2i8)
      else
      icinpu=-ns
      endif
 
      elseif(ichoic.eq.2.or.ichoic.eq.3.or.ichoic.eq.5)then
      if(isup.ne.1)write(ifmt,*)'#events?[i]   modsho?[i]   '
     *,'p(lab)?[f]   pj,tg?[4i]   #par?[i]'
C      read(ifop,*)nevent,modsho,pnll
C     *,laproj,maproj,latarg,matarg,npar
      read(ifop,*)pnll,laproj,maproj,latarg,matarg
       nevent=1
       modsho=1
       npar=8
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,104)nevent,modsho,pnll
     *,laproj,maproj,latarg,matarg,npar
104   format(1x,i6,i5,f10.0,3x,4i4,i6)
      if(laproj.lt.0.or.latarg.lt.0)
     *call utstop('ainitl: code only good for p-p          ')
      if(laproj.lt.0.and.latarg.lt.0)then
      if(isup.ne.1)write(ifmt,*)'id_proj? [i]  id_targ? [i]'
      read(ifop,*)idproj,idtarg
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)idproj,idtarg
      elseif(laproj.lt.0.and.latarg.ge.0)then
      if(isup.ne.1)write(ifmt,*)'id_proj? [i]'
      read(ifop,*)idproj
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)idproj
      elseif(laproj.ge.0.and.latarg.lt.0)then
      if(isup.ne.1)write(ifmt,*)'id_targ? [i]'
      read(ifop,*)idtarg
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)idtarg
      endif
      if(iabs(idproj).lt.1000.or.iabs(idtarg).lt.1000)sigppi=2
 
      else
      stop'wrong choice'
      endif
      PSSHEP(4)=abs(pnll)
      aproj=maproj 
      deast=10.
      if(aproj.gt.9) deast=100.
      if(aproj.gt.99) deast=1000.
      aproj=aproj/deast
      atarg=matarg
      PSSHEP(5) = atarg + aproj
      if(npar.gt.0.and.npar.le.mxpar)then
      do 3 n=1,npar
      if(isup.ne.1)write(ifmt,101)n
101   format(7x,i1,'. parameter? [a6,f]')
      read(ifop,*)parcha(n),parval(n)
      call utlow6(parcha(n))
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)parcha(n),parval(n)
      if(parcha(n).eq.'iversn')iversn=parval(n)
      if(parcha(n).eq.'pud   ')pud   =parval(n)
      if(parcha(n).eq.'pspinl')pspinl=parval(n)
      if(parcha(n).eq.'pspinh')pspinh=parval(n)
      if(parcha(n).eq.'pispn ')pispn =parval(n)
      if(parcha(n).eq.'taunll')taunll=parval(n)
      if(parcha(n).eq.'overlp')overlp=parval(n)
      if(parcha(n).eq.'iopenu')iopenu=parval(n)
      if(parcha(n).eq.'maxres')maxres=parval(n)
      if(parcha(n).eq.'ptf   ')ptf   =parval(n)
      if(parcha(n).eq.'ioptq ')ioptq =parval(n)
      if(parcha(n).eq.'pvalen')pvalen=parval(n)
      if(parcha(n).eq.'ioptf ')ioptf =parval(n)
      if(parcha(n).eq.'irescl')irescl=parval(n)
      if(parcha(n).eq.'tensn ')tensn =parval(n)
      if(parcha(n).eq.'ifrade')ifrade=nint(parval(n))
      if(parcha(n).eq.'ndecay')ndecay=nint(parval(n))
      if(parcha(n).eq.'ndecax')ndecax=nint(parval(n))
      if(parcha(n).eq.'ientro')ientro=nint(parval(n))
      if(parcha(n).eq.'pdiqua')pdiqua=parval(n)
      if(parcha(n).eq.'smas  ')smas  =parval(n)
      if(parcha(n).eq.'uumas ')uumas =parval(n)
      if(parcha(n).eq.'usmas ')usmas =parval(n)
      if(parcha(n).eq.'ssmas ')ssmas =parval(n)
      if(parcha(n).eq.'parea ')parea =parval(n)
      if(parcha(n).eq.'delrem')delrem=parval(n)
      if(parcha(n).eq.'ptq1  ')ptq1  =parval(n)
      if(parcha(n).eq.'ptq2  ')ptq2  =parval(n)
      if(parcha(n).eq.'ptq3  ')ptq3  =parval(n)
      if(parcha(n).eq.'pth   ')pth   =parval(n)
      if(parcha(n).eq.'phard ')phard =parval(n)
      if(parcha(n).eq.'neqmn ')neqmn =parval(n)
      if(parcha(n).eq.'neqmx ')neqmx =parval(n)
      if(parcha(n).eq.'iaqu  ')iaqu  =parval(n)
      if(parcha(n).eq.'wproj ')wproj =parval(n)
      if(parcha(n).eq.'wtarg ')wtarg =parval(n)
      if(parcha(n).eq.'qmust ')stop'*** qmust not used anymore! ***'
      if(parcha(n).eq.'qmust1')qmust1=parval(n)
      if(parcha(n).eq.'qmust2')qmust2=parval(n)
      if(parcha(n).eq.'qmust3')qmust3=parval(n)
      if(parcha(n).eq.'sigppi')sigppi=parval(n)
      if(parcha(n).eq.'grigam')grigam=parval(n)
      if(parcha(n).eq.'grirsq')grirsq=parval(n)
      if(parcha(n).eq.'gridel')gridel=parval(n)
      if(parcha(n).eq.'grislo')grislo=parval(n)
      if(parcha(n).eq.'gricel')gricel=parval(n)
      if(parcha(n).eq.'core  ')core  =parval(n)
      if(parcha(n).eq.'ntrymx')ntrymx=parval(n)
      if(parcha(n).eq.'taumx ')taumx =parval(n)
      if(parcha(n).eq.'nsttau')nsttau=parval(n)
      if(parcha(n).eq.'idpm  ')idpm  =parval(n)
      if(parcha(n).eq.'ymximi')ymximi=parval(n)
      if(parcha(n).eq.'nclean')nclean=parval(n)
      if(parcha(n).eq.'labsys')labsys=parval(n)
      if(parcha(n).eq.'ncolmx')ncolmx=parval(n)
      if(parcha(n).eq.'sigj  ')sigj  =parval(n)
      if(parcha(n).eq.'ptmx  ')ptmx  =parval(n)
      if(parcha(n).eq.'gaumx ')gaumx =parval(n)
      if(parcha(n).eq.'cutmss')cutmss=parval(n)
      if(parcha(n).eq.'prosea')prosea=parval(n)
      if(parcha(n).eq.'rstras')rstras=parval(n)
      if(parcha(n).eq.'fctrmx')fctrmx=parval(n)
      if(parcha(n).eq.'delmss')delmss=parval(n)
      if(parcha(n).eq.'iopbrk')iopbrk=parval(n)
      if(parcha(n).eq.'ndecaw')ndecaw=parval(n)
      if(parcha(n).eq.'jpsi  ')jpsi  =parval(n)
      if(parcha(n).eq.'jpsifi')jpsifi=parval(n)
      if(parcha(n).eq.'ko1ko2')ko1ko2=parval(n)
      if(parcha(n).eq.'ish   ')ish   =parval(n)
      if(parcha(n).eq.'ishsub')ishsub=parval(n)
      if(parcha(n).eq.'ishevt')ishevt=parval(n)
      if(parcha(n).eq.'ipagi ')ipagi =parval(n)
      if(parcha(n).eq.'ivi   ')ivi   =parval(n)
      if(parcha(n).eq.'bmaxim')bmaxim=parval(n)
      if(parcha(n).eq.'bminim')bminim=parval(n)
      if(parcha(n).eq.'istmax')istmax=parval(n)
      if(parcha(n).eq.'radiac')radiac=parval(n)
      if(parcha(n).eq.'radias')radias=parval(n)
      if(parcha(n).eq.'themas')themas=parval(n)
      if(parcha(n).eq.'cutmsq')cutmsq=parval(n)
      if(parcha(n).eq.'taurea')taurea=parval(n)
      if(parcha(n).eq.'kutdiq')kutdiq=parval(n)
      if(parcha(n).eq.'amsiac')amsiac=parval(n)
      if(parcha(n).eq.'iopent')iopent=parval(n)
      if(parcha(n).eq.'iojint')iojint=parval(n)
      if(parcha(n).eq.'amprif')amprif=parval(n)
      if(parcha(n).eq.'delvol')delvol=parval(n)
      if(parcha(n).eq.'deleps')deleps=parval(n)
      if(parcha(n).eq.'uentro')uentro=parval(n)
      if(parcha(n).eq.'taumin')taumin=parval(n)
      if(parcha(n).eq.'taumax')taumax=parval(n)
      if(parcha(n).eq.'numtau')numtau=parval(n)
      if(parcha(n).eq.'elepti')elepti=parval(n)
      if(parcha(n).eq.'elepto')elepto=parval(n)
      if(parcha(n).eq.'angmue')angmue=parval(n)
      if(parcha(n).eq.'iwzzzz')iwzzzz=parval(n)
      if(parcha(n).eq.'imihis')imihis=parval(n)
      if(parcha(n).eq.'isphis')isphis=parval(n)
      if(parcha(n).eq.'iclhis')iclhis=parval(n)
      if(parcha(n).eq.'ijphis')ijphis=parval(n)
      if(parcha(n).eq.'rhophi')rhophi=parval(n)
      if(parcha(n).eq.'ispall')ispall=parval(n)
      if(parcha(n).eq.'wtmini')wtmini=parval(n)
      if(parcha(n).eq.'wtstep')wtstep=parval(n)
      if(parcha(n).eq.'iwcent')iwcent=nint(parval(n))
      if(parcha(n).eq.'kentro')kentro=parval(n)
      if(parcha(n).eq.'ifgdt ')ifgdt=nint(parval(n))
      if(parcha(n).eq.'gen#  ')PSSHEP(3)=parval(n)
      if(parcha(n).eq.'run#  ')VSSHEP(1)=parval(n)
      if(parcha(n).eq.'start#')VSSHEP(2)=parval(n)
3     continue
      endif
 
      if(isup.ne.1)write(ifmt,*)'initializations ...'
      nevnt=nevent
      if(isphis.eq.1)labsys=0
      if(icinpu.ge.1)elepti=engy
      if(idpm.eq.1)iaqu=0
      if(idpm.eq.1)neqmn=2
      if(idpm.eq.1)neqmx=3
      if(iopenu.eq.2)call smassi(themas)
      if(iopenu.eq.2.and.ish.eq.19)call smassp
      if(iopenu.eq.2.and.ish.eq.19)stop
      call sjcen0
 
      if(ichoic.eq.1.or.ichoic.eq.4)then
      yhaha=0
      ypjtl=0
      endif
 
           if(ichoic.eq.2.or.ichoic.eq.3.or.ichoic.eq.5)then

      call idmass(idproj,amproj)
      call idmass(idtarg,amtarg)
         if(pnll.lt.0)then
      engy=-pnll
      pnll=sqrt( ((engy**2-amproj**2-amtarg**2)/2/amtarg)**2-amproj**2 )
         else
      engy=sqrt( 2*sqrt(pnll**2+amproj**2)*amtarg+amtarg**2+amproj**2 )
         endif
      if(engy.lt.12.)
     *call utstop('ainitl: incident energy too small       ')
      engyi=engy
      pnlli=pnll
         if(pnll/amproj.lt.1e2)then
      trm=sqrt(pnll**2+amproj**2)
      engy=sqrt((trm+amtarg-pnll)*(trm+amtarg+pnll))
         else
      trm=amproj**2/2/pnll+amtarg
      engy=sqrt(trm*(2*pnll+trm))
         endif
      d1=abs(pnlli-pnll)/pnll
      d2=abs(engyi-engy)/engy
         if(d1.gt.1e-3.or.d2.gt.1e-3)then
      if(ish.ge.0)then
      call utmsg('ainitl')
      write(ifch,*)'*****  pnll,pnlli:',pnll,pnlli
      write(ifch,*)'*****  engy,engyi:',engy,engyi
      call utmsgf
      endif
         endif
      s=engy**2
      pnllx=utpcm(engy,amproj,amtarg)
      yhaha=alog((sqrt(pnll**2+s)+pnll)/sqrt(s))
      ypjtl=alog((sqrt(pnll**2+amproj**2)+pnll)/amproj)

           endif

      if(ichoic.eq.4.or.ichoic.eq.5)goto1000

      if(ichoic.eq.1)then
      s=12.**2
      endif

      qmust=qmust1+qmust2*alog(engy)+qmust3*alog(engy)**2
      ptq=ptq1+ptq2*alog(engy)+ptq3*alog(engy)**2
      if(phard.lt.0.)then
c-c   phard=1-exp(-0.023329-0.000351*sqrt(s))
c-c   phard=0
c-c   if(sqrt(s).gt.100.)
c-c  *phard=1-exp(-0.000275*(sqrt(s)-100.))
      phard=0.030+0.12*(alog10(s)-alog10(30.**2))
      phard=amin1(1.,phard)
      phard=amax1(0.030,phard)
      endif
 
           if(ichoic.eq.1.or.ichoic.eq.2)then

      call jdecin(.false.)
      c=ptmx
      quama=0.
      call utquaf(sptf,nptf,xptf,qptfu,wptf,0.,.33*c,.66*c,c)
      quama=smas
      call utquaf(sptf,nptf,xptf,qptfs,wptf,0.,.33*c,.66*c,c)
      quama=uumas
      call utquaf(sptf,nptf,xptf,qptfuu,wptf,0.,.33*c,.66*c,c)
      quama=usmas
      call utquaf(sptf,nptf,xptf,qptfus,wptf,0.,.33*c,.66*c,c)
      quama=ssmas
      call utquaf(sptf,nptf,xptf,qptfss,wptf,0.,.33*c,.66*c,c)
      c=gaumx
      call utquaf(sgau,ngau,xgau,qgau,wgau,0.,.33*c,.66*c,c)
      call utquaf(ssplit,nsplit,xsplit,qsplit,wsplit,0.,.33,.66,1.)
      call utquaf(ssplix,nsplit,xsplit,qsplix,wsplit,0.,.33,.66,1.)
 
      xcut=cutmsq/sqrt(s)
      b=amin1(.05,xcut*500)
      a=amin1(.2*b,xcut*100)
      ipio=0
      if(iabs(idproj).lt.1000)ipio=1
      call utquaf(stzzvl,nstru,xstru,qvap,wstru,0.,a,b,1.)
      qvapc=qvap(nstru)
      call utquaf(stzzse,nstru,xstru,qsep,wstru,0.,a,b,1.)
      qsepc=qsep(nstru)
      ipio=0
      if(iabs(idtarg).lt.1000)ipio=1
      call utquaf(stzzvl,nstru,xstru,qvat,wstru,0.,a,b,1.)
      qvatc=qvat(nstru)
      call utquaf(stzzse,nstru,xstru,qset,wstru,0.,a,b,1.)
      qsetc=qset(nstru)
      if(ish.eq.16)write(ifch,*)' '
      if(ish.eq.16)write(ifch,*)'sr anitl: qvapc,qsepc,qvatc,qsetc:'
      if(ish.eq.16)write(ifch,*)qvapc,qsepc,qvatc,qsetc
      if(ish.eq.16)write(ifch,*)' '
 
      if(prosea.ge.0.)then
      qvapc=1.
      qsepc=prosea
      qvatc=1.
      qsetc=prosea
      endif
 
      xcut=cutmss/sqrt(s)
      b=amin1(.05,xcut*500)
      a=amin1(.2*b,xcut*100)
      ipio=0
      if(iabs(idproj).lt.1000)ipio=1
      call utquaf(stzzvl,nstru,xstru,qvap,wstru,0.,a,b,1.)
      call utquaf(stzzse,nstru,xstru,qsep,wstru,0.,a,b,1.)
      ipio=0
      if(iabs(idtarg).lt.1000)ipio=1
      call utquaf(stzzvl,nstru,xstru,qvat,wstru,0.,a,b,1.)
      call utquaf(stzzse,nstru,xstru,qset,wstru,0.,a,b,1.)
 
      if(ish.eq.21)call uthsea
      if(ish.eq.21)stop
 
      c=ptmx
      call utquaf(sptj,nptj,xptj,qptj,wptj,0.,.33*c,.66*c,c)
      c=ptmx
      call utquaf(spth,nptq,xptq,qpth,wptq,0.,.33*c,.66*c,c)
      call utquaf(sptq,nptq,xptq,qptq,wptq,0.,.33*c,.66*c,c)
      call racpro('gri',qmust,nprbms,prbms)
 
           endif
 
      massnr=maproj
      rmproj=0
      if(ichoic.le.3.and.massnr.gt.1)then
      c=1.19*massnr**(1./3.)-1.61*massnr**(-1./3.)+fctrmx*.54
      rmproj=c
      call utquaf(sdensi,ndep,xdep,qdep,wdep,0.,.33*c,.66*c,c)
      endif
 
      massnr=matarg
      rmtarg=0
      if(ichoic.le.3.and.massnr.gt.1)then
      c=1.19*massnr**(1./3.)-1.61*massnr**(-1./3.)+fctrmx*.54
      rmtarg=c
      call utquaf(sdensi,ndet,xdet,qdet,wdet,0.,.33*c,.66*c,c)
      endif
 
      if(idpm.eq.1)qsepc=0
      if(idpm.eq.1)qsetc=0
 
      bmax=rmproj+rmtarg

      if(iwzzzz.eq.1)call wzzzzi
 
      rewind ifdt
      if(isup.ne.1)write(ifmt,*)'generate events ...'
 
1000  continue
      if(ish.eq.17)write(ifch,*)'ainitl (exit)'
      return
 
      entry shopar
      if(isup.ne.1.and.ish.ge.90)write(ifmt,102)('*',i=1,69)
     * ,taunll,maxres,ptf,ptq,ioptq,ptmx,neqmn,iaqu,wtarg
     * ,wproj,qmust,sigppi,core,fctrmx,ncolmx,labsys,irescl,overlp
     * ,ntrymx,delmss,gaumx,bmaxim,pud,pspinl,pspinh,pispn
     * ,nclean,jpsi,cutmss,rstras,istmax,tensn,neqmx,ipagi,ndecay
     * ,pdiqua,parea,delrem,taumx,nsttau,sigj,jpsifi,bminim
     * ,radiac,taumin,taumax,numtau,themas,amsiac,elepti
102   format(' ',69a1
     */' *  taunll=',f7.3,'  maxres=',i7,'  ptf   =',f7.3
     *,'  ptq   =',f7.3,'   *'/' *  ioptq =',i7,'  ptmx  =',f7.3
     *,'  neqmn =',i7,'  iaqu  =',i7,'   *'/' *  wtarg =',f7.3
     *,'  wproj =',f7.3,'  qmust =',f7.3,'  sigppi=',f7.3,'   *'
     */' *  core  =',f7.3,'  fctrmx=',f7.3,'  ncolmx=',i7
     *,'  labsys=',  i7,'   *'/' *  irescl=',i7,'  overlp=',f7.3
     *,'  ntrymx=',i7,'  delmss=',f7.3,'   *'/
     *,'  gaumx =',f7.3,'  bmaxim=',f7.1,'   *'/' *  pud   =',f7.3
     *,'  pspinl=',f7.3,'  pspinh=',f7.3,'  pispn =',f7.3,'   *'
     */' *  nclean=',i7,'  jpsi  =',i7,'  cutmss=',f7.3
     *,'  rstras=',f7.3,'   *'/' *  istmax=',i7,'  tensn =',f7.3
     *,'  neqmx =',i7,'  ipagi =',i7,'   *'/' *  ndecay=',i7
     *,'  pdiqua=',f7.3,'  parea =',f7.3,'  delrem=',f7.3,'   *'
     */' *  taumx =',f7.3,'  nsttau=',i7,'  sigj  =',f7.3
     *,'  jpsifi=',i7,'   *'/' *  bminim=',f7.1,'  radiac=',f7.3
     *,'  taumin=',f7.3,'  taumax=',f7.3,'   *'/' *  numtau=',i7
     *,'  themas=',f7.3,'  amsiac=',f7.3,'  elepti=',f7.3,'   *')
      if(isup.ne.1.and.ish.ge.90)write(ifmt,108)
     *  elepto,angmue,smas,uumas,usmas,ssmas,iopbrk,ndecaw
     * ,imihis,kentro,radias,isphis,rhophi,iclhis,ishsub,iopenu
     * ,iopent,cutmsq,kutdiq,idpm,taurea,ispall,ymximi,pth,phard
     * ,ioptf,prosea,ndecax,wtmini,wtstep,iwcent,ishevt,pvalen
     * ,ifrade,ientro,grigam,grirsq,gridel,grislo,gricel,ijphis
     * ,uentro,iwzzzz,iojint,amprif,delvol,deleps
     * ,('*',i=1,69)
108   format(
     * ' *  elepto=',f7.3,'  angmue=',f7.3,'  smas  =',f7.3
     *,'  uumas =',f7.3,'   *'/' *  usmas =',f7.3,'  ssmas =',f7.3
     *,'  iopbrk=',i7,'  ndecaw=',i7  ,'   *'/' *  imihis=',i7
     *,'  kentro=',i7,'  radias=',f7.3,'  isphis=',i7  ,'   *'
     */' *  rhophi=',f7.3,'  iclhis=',i7  ,'  ishsub=',i7
     *,'  iopenu=',i7  ,'   *'/' *  iopent=',i7,'  cutmsq=',f7.3
     *,'  kutdiq=',i7  ,'  idpm  =',i7  ,'   *'/' *  taurea=',f7.3
     *,'  ispall=',i7  ,'  ymximi=',f7.3,'  pth   =',f7.3,'   *'
     */' *  phard =',f7.3,'  ioptf =',i7  ,'  prosea=',f7.3
     *,'  ndecax=',i7  ,'   *'/' *  wtmini=',f7.3,'  wtstep=',f7.3
     *,'  iwcent=',i7  ,'  ishevt=',i7  ,'   *'/' *  pvalen=',f7.3
     *,'  ifrade=',i7  ,'  ientro=',i7  ,'  grigam=',f7.3,'   *'
     */' *  grirsq=',f7.3,'  gridel=',f7.3,'  grislo=',f7.3
     *,'  gricel=',f7.3,'   *'/' *  ijphis=',i7  ,'  uentro=',f7.3
     *,'  iwzzzz=',i7  ,'  iojint=',i7  ,'   *'/' *  amprif=',f7.3
     *,'  delvol=',f7.3,'  deleps=',f7.3,'         ',7x  ,'   *'
     */' ',69a1)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine astore
c     writes the results of a simulation on tape.
c     contains a description of the stored variables.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      dimension nrptla(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      common /cptlu/nptlu

c -- dhwi
c put in call to zebra output

      call astore_egz
      return

c -- dhwi

      ish0=ish
      ishs0=ishsub
      if(ishevt.ne.0.and.nrevt.ne.ishevt)ish=0
      if(ishevt.ne.0.and.nrevt.ne.ishevt)ishsub=0
      if(ishsub/100.eq.8)ish=mod(ishsub,100)
      if(ish.eq.14)call uttima('*** astore *** ')
      if(ish.eq.17)write(ifch,*)' '
      if(ish.eq.17)write(ifch,*)'astore (entry)'
      if(ichoic.eq.4.or.ichoic.eq.5)goto1000
      if(nevt.ne.1)goto1000
      nptlu=max0(nptl,nptlu)
      if(iwzzzz.eq.1)goto9999
      if(imihis.eq.1)goto9999
      if(isphis.eq.1)goto9999
      if(iclhis.eq.1)goto9999
      if(ijphis.eq.1)goto9999
      levt=1
      lptl=3
 
      nptls=0
      do 1 i=1,nptl
      nrptla(i)=-999
      if(istptl(i).gt.istmax)goto1
      nptls=nptls+1
      nrptla(i)=nptls
1     continue
 
c     store event variables:
 
      write(ifdt)levt,nrevt,nptls,bimevt,kolevt,colevt,pmxevt
     *,egyevt,npjevt,ntgevt
 
c     levt................... record label (levt=1)
c     nrevt.................. event number
c     nptls ................. number of (stored!) particles per event
c     bimevt ................ impact parameter
c     kolevt,colevt ......... real/effective # of collisions
c     pmxevt ................ reference momentum
c     egyevt ................ pp cm energy (had) or string energy (str)
c     npjevt,ntgevt ......... # of proj/targ participants
 
      if(ichoic.eq.3)goto1000
 
      do 5 i=1,nptl
      if(nrptla(i).le.0)goto5
      nrptl=nrptl+1
      if(iorptl(i).gt.0)then
      ioptl=nrptla(iorptl(i))
      else
      ioptl=iorptl(i)
      endif
      if(jorptl(i).gt.0)then
      joptl=nrptla(jorptl(i))
      else
      joptl=jorptl(i)
      endif
 
c     store particle variables:
 
      write(ifdt)lptl,nrevt,nrptla(i),i,idptl(i),(pptl(l,i),l=1,5),ioptl
     *,joptl,istptl(i),(xorptl(l,i),l=1,4),(tivptl(l,i),l=1,2)
     *,nqjptl(i)
 
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
 
5     continue
4     continue
      goto1000
 
c     fill histograms
c     ---------------
9999  continue
      if(iwzzzz.eq.1)call wzzzz
      if(imihis.eq.1)call wimi
      if(isphis.eq.1)call wspa
      if(iclhis.eq.1)call wclu
      if(ijphis.eq.1)call wjpsi
 
1000  continue
      call uttimx(timaa)
      if(isup.ne.1.and.mod(nrevt,modsho).eq.0)
     *write(ifmt,*)nrevt,timaa
      if(ish.eq.14)call uttima('    astore f   ')
      if(ish.eq.17)write(ifch,*)'astore (exit)'
      ish=ish0
      ishsub=ishs0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine atermi
c     terminates run
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cifipr/ifipr /cptlu/nptlu
      character dcn*3
c MAB
      integer istat, egzend
c MAB

c>>>>>dhwi + MAB
c>>>>>ensure clean ending of zebra file
      istat=egzend(0)
c>>>>>dhwi + MAB

      if(ish.eq.14)call uttima('*** atermi *** ')
      call veanly
      if(ichoic.eq.4.or.ichoic.eq.5)stop
      if(ifipr.eq.1)then
      call shopar
      if(inoiac.gt.0.or.ilamas.gt.0)then
      if(ish.ge.90)then
      call utmsg('atermi')
      if(inoiac.gt.0)write(ifch,*)'*****  inoiac = ',inoiac
      if(ilamas.gt.0)write(ifch,*)'*****  ilamas = ',ilamas
      call utmsgf
      endif
      endif
      if(isup.ne.1.and.imsg.gt.0)write(ifmt,106)
     *('*',i=1,69),imsg,ifch,('*',i=1,69)
106   format(' ',69a1/' *',67x,'*'
     */' *',12x,i5,' message(s) - written to file ',i2,18x,'*'
     */' *',67x,'*'/' ',69a1)
      endif
      if(isup.ne.1)write(ifmt,*)
     *'stop execution (stp) or run again (agn)?'
      read(ifop,*)dcn
      call utlow3(dcn)
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)dcn
      if(isup.ne.1)write(ifmt,*)'nptl_max:',nptlu
      if(dcn.eq.'stp')stop
      if(dcn.eq.'agn')goto1000
      write(ifmt,*)'dont understand. stp assumed.'
      stop
 
1000  continue
      if(ish.eq.14)call uttima('    atermi f   ')
      return
      end
 
c-----------------------------------------------------------------------
      subroutine avenus
c     generates one venus event
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      ish0=ish
      ishs0=ishsub
      if(ishevt.ne.0.and.nrevt+1.ne.ishevt)ish=0
      if(ishevt.ne.0.and.nrevt+1.ne.ishevt)ishsub=0
      if(nrevt.eq.0.and.(ish.eq.13.or.ish.eq.14))call uttimt
      if(ish.eq.14)call uttima('*** avenus *** ')
      if(ish.eq.17)write(ifch,*)' '
      if(ish.eq.17)write(ifch,*)'avenus (entry)'
      ipage=0
      call utpage
      if(ish.ge.91)
     *write(ifch,113)('-',l=1,79),ipage,nrevt+1,('-',l=1,79)
113   format(/1x,79a1/1x,i5,'.page            '
     *,5x,'evt:',i4,5x,/1x,79a1/)
      if(ichoic.eq.4.or.ichoic.eq.5)goto1000
      bimevt=-1
      ntevt0=ntevt
3     ntevt=ntevt0
2     ntevt=ntevt+1
      if(ichoic.eq.1)then
      call jetgen(ier)
      if(ier.eq.1)goto3
      elseif(ichoic.eq.2.or.ichoic.eq.3)then
      call nucoll
      endif
      if(ichoic.eq.1.or.ichoic.eq.2)then
      call jfrade(ier)
      if(ier.eq.1)goto3
      endif
      if((ichoic.eq.2.or.ichoic.eq.3).and.nevt.eq.0)then
      bimevt=-1
      goto2
      endif
      if(jerr.gt.0)
     *call utstop('avenus: jerr > 0                        ')
      nrevt=nrevt+1
 
1000  continue
      if(ish.eq.14)call uttima('    avenus f   ')
      if(ish.eq.17)write(ifch,*)'avenus (exit)'
      ish=ish0
      ishsub=ishs0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hahabs(proj,targ,iap,iat,iskip,irethh)
c     performs a basic (one color exchange) hadron-hadron collision
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      common/ckol/kol
      parameter (nsi=6,nsix=100,nvmax=10)
      parameter (nha=nsix+2)
      real proj(nsi,nha),targ(nsi,nha)
     *,strap(nsi,nsix+1),strat(nsi,nsix+1),str0(nsi,nsix+1)
     *,strlp(nsi,nsix+1),strlt(nsi,nsix+1)
     *,strpwo(nsi,2),strtwo(nsi,2)
      parameter (mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      real siap(nsi),sibp(nsi),silp(nsi),sixp(nsi,nsix),
     *siat(nsi),sibt(nsi),silt(nsi),sixt(nsi,nsix)
     *,siapwo(nsi),sibpwo(nsi),siatwo(nsi),sibtwo(nsi)
      common /cnew/kotri,newico,newcol
     */clead/lead,cooav3,cooav4
     */cnce/nces,ncolex
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      character xflap*3,xflbp*3,xflat*3,xflbt*3
      parameter (nflav=6)
      integer icvp(2),icvt(2),jcvp(nflav,2),jcvt(nflav,2)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cstsh/nstsh
      common/cpzstr/pzstrl,estrl,istrl,isea
      real psum(5),psumx(5)
      r=rangen()
      ntry=0
      irethh=0
      call utremb(proj,targ,2)
      goto9
 
9994  iret=0
      call utrest(proj,targ,2)
9     do 8 n=1,nsi
      do 8 nx=1,nsix
      sixp(n,nx)=0
8     sixt(n,nx)=0
      call utkprj(proj,kmaxp)
      call utkprj(targ,kmaxt)
      if(ish.ge.91)then
      write(ifch,100)((proj(i,j),i=1,nsi),j=2,kmaxp+2)
      write(ifch,102)((targ(i,j),i=1,nsi),j=2,kmaxt+2)
      endif
      if(iskip.eq.1)goto9999
 
9997  call hdecmp(proj,silp,sixp)
      call hdecmp(targ,silt,sixt)
1     icvp(1)=nint(proj(5,1))
      icvp(2)=nint(proj(6,1))
      icvt(1)=nint(targ(5,1))
      icvt(2)=nint(targ(6,1))
      call hapapa(silp,iflap,xflap,ptap,ptahp,iflbp,xflbp,ptbp,ptbhp
     *,nqaqp,icvp)
      call hapapa(silt,iflat,xflat,ptat,ptaht,iflbt,xflbt,ptbt,ptbht
     *,nqaqt,icvt)
      if(nqaqp*nqaqt.eq.-1)goto1
 
9998  ntry=ntry+1
      if(ntry.gt.ntrymx+1)then
      iskip=1
      goto9999
      endif
 
      call hdecmp(proj,silp,sixp)
      call hdecmp(targ,silt,sixt)
      call hastfs(silp,sixp,iflap,xflap,ptap,ptahp
     *,silt,sixt,iflbt,xflbt,ptbt,ptbht,strap,iret,1)
      if(iret.eq.1)then
      iskip=1
      goto9999
      endif
      if(iret.eq.2)goto9997
      if(iret.eq.3)then
      if(xflbp(1:2).eq.'va')xflbp(1:2)='se'
      if(xflat(1:2).eq.'va')xflat(1:2)='se'
      iflat=iflbt
      iflbp=iflap
      goto9998
      endif
      if(iret.eq.5)then
      if(xflap(1:2).eq.'se'.and.xflbp(1:2).eq.'se')then
      iflbp=-iabs(iflbp)
      iflap=iflbp
      endif
      if(xflat(1:2).eq.'se'.and.xflbt(1:2).eq.'se')then
      iflat=-iabs(iflat)
      iflbt=iflat
      endif
      goto9998
      endif
      call hastfs(silt,sixt,iflat,xflat,ptat,ptaht
     *,silp,sixp,iflbp,xflbp,ptbp,ptbhp,strat,iret,2)
      if(iret.eq.1)then
      iskip=1
      goto9999
      endif
      if(iret.eq.2)goto9997
      if(iret.eq.3)then
      if(xflap(1:2).eq.'va')xflap(1:2)='se'
      if(xflbt(1:2).eq.'va')xflbt(1:2)='se'
      iflbt=iflat
      iflap=iflbp
      goto9998
      endif
      if(iret.eq.5)then
      if(xflap(1:2).eq.'se'.and.xflbp(1:2).eq.'se')then
      iflap=-iabs(iflap)
      iflbp=iflap
      endif
      if(xflat(1:2).eq.'se'.and.xflbt(1:2).eq.'se')then
      iflbt=-iabs(iflbt)
      iflat=iflbt
      endif
      goto9998
      endif
 
      call hastfc(silp,sixp,iret)
      if(iret.eq.0)
     *call hastfc(silt,sixt,iret)
           if(iret.eq.1)then
      if(ish.ge.91)write(ifch,*)'light string --> redo hastfs'
      goto9998
           elseif(iret.eq.2)then
      if(ish.ge.91)write(ifch,*)'jc>9 --> redo hapapa'
      goto9997
           endif
 
9999  continue
 
      if(iskip.eq.1)then
      if(ish.ge.91)write(ifch,*)'skip'
      colevt=colevt-1./ncolex
      call hdecmp(proj,silp,sixp)
      call hdecmp(targ,silt,sixt)
      call hastfw(silp,sixp)
      call hastfw(silt,sixt)
      endif
 
      call hastfl(silp,sixp,strlp)
      call hastfl(silt,sixt,strlt)
 
      if(iskip.eq.1)goto9995
 
      call iddeco(icvp,jcvp)
      call iddeco(icvt,jcvt)
      if(xflap(1:2).eq.'va')jcvp(iabs(iflap),1)=jcvp(iabs(iflap),1)-1
      if(xflbp(1:2).eq.'va')jcvp(iabs(iflbp),2)=jcvp(iabs(iflbp),2)-1
      if(xflat(1:2).eq.'va')jcvt(iabs(iflat),1)=jcvt(iabs(iflat),1)-1
      if(xflbt(1:2).eq.'va')jcvt(iabs(iflbt),2)=jcvt(iabs(iflbt),2)-1
      call idenco(jcvp,icvp,ireten)
      if(ireten.eq.1)
     *call utstop('hapapa: idenco ret code = 1             ')
      call idenco(jcvt,icvt,ireten)
      if(ireten.eq.1)
     *call utstop('hapapa: idenco ret code = 1             ')
 
      proj(5,1)=icvp(1)
      proj(6,1)=icvp(2)
      targ(5,1)=icvt(1)
      targ(6,1)=icvt(2)
      call utksix(sixp,kmaxpn)
      call utksix(sixt,kmaxtn)
      if(kmaxpn+1.gt.nsix)
     *call utstop('hahabs: dimension nsix too small        ')
      if(kmaxtn+1.gt.nsix)
     *call utstop('hahabs: dimension nsix too small        ')
      kmxp=kmaxpn+1
      kmxt=kmaxtn+1
      do 5 n=1,nsi
      proj(n,2)=strlp(n,1)
      targ(n,2)=strlt(n,1)
      do 6 j=1,kmxp
6     proj(n,2+j)=sixp(n,j)
      do 7 j=1,kmxt
7     targ(n,2+j)=sixt(n,j)
      kproja(2,iproj)=kol
      if(kmaxpn.gt.kmaxp)then
      do 2 k=kmaxp+1,kmaxpn
2     kproja(2+k,iproj)=kol
      endif
      ktarga(2,itarg)=kol
      if(kmaxtn.gt.kmaxt)then
      do 3 k=kmaxt+1,kmaxtn
3     ktarga(2+k,itarg)=kol
      endif
5     continue
 
      nstsh=0
      lead=0
      isplt=0
11    call hastpr(strap,isplt)
      if(isplt.ne.0)
     *call utstop('hahabs: isplt /= 0                      ')
      isplt=0
12    call hastpr(strat,isplt)
      if(isplt.ne.0)
     *call utstop('hahabs: isplt /= 0                      ')
 
      if(kmaxpn.gt.kmaxp)then
      proj(3,1)=proj(3,1)+coord(3,kol)
      proj(4,1)=proj(4,1)+coord(4,kol)
      endif
      if(kmaxtn.gt.kmaxt)then
      targ(3,1)=targ(3,1)+coord(3,kol)
      targ(4,1)=targ(4,1)+coord(4,kol)
      endif
      kmaxp=kmaxpn
      kmaxt=kmaxtn
 
9995  lead=1
 
c     write leading string (proj)
c     ---------------------------
 
           if(iap.eq.1)then
 
      if(kmaxp.gt.0)then
      cooav3=proj(3,1)/kmaxp
      cooav4=proj(4,1)/kmaxp
      else
      cooav3=coord(3,kol)
      cooav4=coord(4,kol)
      endif
      call utkstr(strlp,kmaxor)
      pzstrl=0
      estrl=0
      isea=1
      do i=1,4
      psum(i)=0
      enddo
           do k=1,kmaxor
      do i=1,4
      psum(i)=psum(i)+strlp(i,k)
      enddo
      ampr=0
c-c   if(k.eq.1)ampr=.94
      pzstrl=pzstrl+strlp(3,k)
      estrl=estrl+sqrt(strlp(4,k)**2+ampr**2)
      if(k.ge.2.and.strlp(5,k).gt.0..and.strlp(6,k).gt.0.)isea=0
           enddo
      psum(5)=psum(4)**2-psum(3)**2-psum(2)**2-psum(1)**2
      if(psum(5).gt.0.)psum(5)=sqrt(psum(5))
      if(pzstrl/pnllx.gt.1.-0.850/engy**2.and.iskip.ne.1)goto1002
      istrl=0
      if(nstsh.eq.0.and.isea.eq.1
     *.and.pzstrl/pnllx.gt.1.-10.000/engy**2)istrl=1
      do i=1,nsi
      do k=1,kmaxor
      str0(i,k)=strlp(i,k)
      enddo
      enddo
      nstr0=nstr
17    nstr=nstr0
      isplt=0
15    continue
      do i=1,nsi
      do k=1,kmaxor
      strlp(i,k)=str0(i,k)
      enddo
      enddo
13    call hastpr(strlp,isplt)
      if(isplt.eq.-1.and.iskip.eq.1)goto9996
      if(isplt.eq.-1)goto9994
      if(isplt.eq.-3)goto15
      if(isplt.eq.-4)goto1001
      if(isplt.eq.-5)goto17
      if(isplt.gt.0)goto13
           if(nstr.gt.nstr0+1)then
      do i=1,4
      psumx(i)=0
      do j=nstr0+1,nstr
      psumx(i)=psumx(i)+pstr(i,j)
      enddo
      enddo
      psumx(5)=psumx(4)**2-psumx(3)**2-psumx(2)**2-psumx(1)**2
      if(psumx(5).gt.0.)psumx(5)=sqrt(psumx(5))
c-c   write(6,*)' '
c-c   do j=nstr0+1,nstr
c-c   write(6,109)j,(icstr(k,j)/100,k=1,4)
c-c  *,sqrt(pstr(1,j)**2+pstr(2,j)**2),pstr(3,j)/pnllx,pstr(5,j)
c-c  *,irlstr(j)
c-c   enddo
c-c   write(6,*)psum
c-c   write(6,*)psumx
      do i=1,4
      do j=nstr0+1,nstr
      pstr(i,j+nstr-nstr0)=pstr(i,j)
      enddo
      enddo
      call hrescl(nstr0+1,nstr,psum,ifail)
      if(ifail.gt.0)then
c-c   write(6,*)'rescale failed'
      do i=1,4
      do j=nstr0+1,nstr
      pstr(i,j)=pstr(i,j+nstr-nstr0)
      enddo
      enddo
      endif
      do i=1,4
      psumx(i)=0
      do j=nstr0+1,nstr
      psumx(i)=psumx(i)+pstr(i,j)
      enddo
      enddo
      psumx(5)=psumx(4)**2-psumx(3)**2-psumx(2)**2-psumx(1)**2
      if(psumx(5).gt.0.)psumx(5)=sqrt(psumx(5))
c-c   do j=nstr0+1,nstr
c-c   write(6,109)j,(icstr(k,j)/100,k=1,4)
c-c  *,sqrt(pstr(1,j)**2+pstr(2,j)**2),pstr(3,j)/pnllx,pstr(5,j)
c-c  *,irlstr(j)
c-c   enddo
c-c   write(6,*)psum
c-c   write(6,*)psumx
           endif
109   format(' /cstr/',i4,3x,4i5,2x,3(1x,e9.3),i4)
 
           endif
 
c     write leading string (targ)
c     ---------------------------
 
           if(iat.eq.1)then
 
      if(kmaxt.gt.0)then
      cooav3=targ(3,1)/kmaxt
      cooav4=targ(4,1)/kmaxt
      else
      cooav3=coord(3,kol)
      cooav4=coord(4,kol)
      endif
      call utkstr(strlt,kmaxor)
      pzstrl=0
      estrl=0
      isea=1
      do i=1,4
      psum(i)=0
      enddo
           do k=1,kmaxor
      do i=1,4
      psum(i)=psum(i)+strlt(i,k)
      enddo
      ampr=0
c-c   if(k.eq.1)ampr=.94
      pzstrl=pzstrl+strlt(3,k)
      estrl=estrl+sqrt(strlt(4,k)**2+ampr**2)
      if(k.ge.2.and.strlt(5,k).gt.0..and.strlt(6,k).gt.0.)isea=0
           enddo
      psum(5)=psum(4)**2-psum(3)**2-psum(2)**2-psum(1)**2
      if(psum(5).gt.0.)psum(5)=sqrt(psum(5))
      if(-pzstrl/pnllx.gt.1.-0.850/engy**2.and.iskip.ne.1)goto1002
      istrl=0
      if(nstsh.eq.0.and.isea.eq.1
     *.and.-pzstrl/pnllx.gt.1.-10.000/engy**2)istrl=1
      do i=1,nsi
      do k=1,kmaxor
      str0(i,k)=strlt(i,k)
      enddo
      enddo
      nstr0=nstr
18    nstr=nstr0
      isplt=0
16    continue
      do i=1,nsi
      do k=1,kmaxor
      strlt(i,k)=str0(i,k)
      enddo
      enddo
14    call hastpr(strlt,isplt)
      if(isplt.eq.-1.and.iskip.eq.1)goto9996
      if(isplt.eq.-1)goto9994
      if(isplt.eq.-3)goto16
      if(isplt.eq.-4)goto1001
      if(isplt.eq.-5)goto18
      if(isplt.gt.0)goto14
           if(nstr.gt.nstr0+1)then
      do i=1,4
      psumx(i)=0
      do j=nstr0+1,nstr
      psumx(i)=psumx(i)+pstr(i,j)
      enddo
      enddo
      psumx(5)=psumx(4)**2-psumx(3)**2-psumx(2)**2-psumx(1)**2
      if(psumx(5).gt.0.)psumx(5)=sqrt(psumx(5))
c-c   write(6,*)' '
c-c   do j=nstr0+1,nstr
c-c   write(6,109)j,(icstr(k,j)/100,k=1,4)
c-c  *,sqrt(pstr(1,j)**2+pstr(2,j)**2),pstr(3,j)/pnllx,pstr(5,j)
c-c  *,irlstr(j)
c-c   enddo
c-c   write(6,*)psum
c-c   write(6,*)psumx
      do i=1,4
      do j=nstr0+1,nstr
      pstr(i,j+nstr-nstr0)=pstr(i,j)
      enddo
      enddo
      call hrescl(nstr0+1,nstr,psum,ifail)
      if(ifail.gt.0)then
c-c   write(6,*)'rescale failed'
      do i=1,4
      do j=nstr0+1,nstr
      pstr(i,j)=pstr(i,j+nstr-nstr0)
      enddo
      enddo
      endif
      do i=1,4
      psumx(i)=0
      do j=nstr0+1,nstr
      psumx(i)=psumx(i)+pstr(i,j)
      enddo
      enddo
      psumx(5)=psumx(4)**2-psumx(3)**2-psumx(2)**2-psumx(1)**2
      if(psumx(5).gt.0.)psumx(5)=sqrt(psumx(5))
c-c   do j=nstr0+1,nstr
c-c   write(6,109)j,(icstr(k,j)/100,k=1,4)
c-c  *,sqrt(pstr(1,j)**2+pstr(2,j)**2),pstr(3,j)/pnllx,pstr(5,j)
c-c  *,irlstr(j)
c-c   enddo
c-c   write(6,*)psum
c-c   write(6,*)psumx
           endif
 
           endif
 
c     exit
c     ----
 
      if(.not.(ish.ge.91))return
 
      write(ifch,100)((proj(i,j),i=1,nsi),j=2,kmaxp+2)
      write(ifch,102)((targ(i,j),i=1,nsi),j=2,kmaxt+2)
100   format(' proj:',4f11.5,2f8.0/
     *   ,20('      ',4f11.5,2f8.0/))
102   format(' targ:',4f11.5,2f8.0/
     *   ,50('      ',4f11.5,2f8.0/))
101   format('  ',i3,'. trial')
      return
 
9996  continue
      if(ish.ge.91)then
      write(ifch,*)' '
      write(ifch,*)'skip impossible, strl cant be stored'
      endif
      iskip=2
      return
 
1001  continue
      if(ish.ge.91)then
      write(ifch,*)' '
      write(ifch,*)
     *'sr hahabs: not justified complete absorption -> irethh=1'
      endif
      irethh=1
      if(iskip.eq.1)goto9996
      return
 
1002  continue
      if(ish.ge.91)then
      write(ifch,*)' '
      write(ifch,*)'sr hahabs: too fast multi-string -> irethh=1'
      endif
      irethh=1
      if(iskip.eq.1)
     *call utstop('hahabs: iskip=1 should not happen       ')
      return
 
      end
 
c-----------------------------------------------------------------------
      subroutine hapapa(sil,ifla,xfla,pta,ptah,iflb,xflb,ptb,ptbh,nqaq
     *,icv)
c     determines momenta and flavor of participating partons in a hadron
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (nsi=6,nsix=100,nflav=6)
      parameter (nha=nsix+2)
      real sia(nsi),sib(nsi),sil(nsi),probab(nflav),probsu(nflav+1)
      integer jv0(nflav),jc(nflav,2),jcx(nflav),id(2)
     *,icv(2),jcv(nflav,2)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      parameter (nstru=2049)
      common/stru/xstru(nstru),wstru(nstru)
     *,qvap(nstru),qvat(nstru),qsep(nstru),qset(nstru)
      common /cnew/kotri,newico,newcol
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      character xfla*3,xflb*3,chpt*1
 
      sgn=sign(1.,sil(3))
      if(sgn.gt.0.)qvalc=qvapc
      if(sgn.gt.0.)qseac=qsepc
      if(sgn.lt.0.)qvalc=qvatc
      if(sgn.lt.0.)qseac=qsetc
      if(sgn.gt.0.)chpt='p'
      if(sgn.lt.0.)chpt='t'
      if(sgn.gt.0.)w=wproj
      if(sgn.lt.0.)w=wtarg
      call iddeco(icv,jcv)
      nvq=0
      nva=0
      do 12 i=1,nflav
           if(ish.ge.0.and. (jcv(i,1).lt.0.or.jcv(i,2).lt.0) )then
      call utmsg('hahabs')
      write(ifch,*)'*****  negative jcv'
      write(ifch,*)'jcv:'
      write(ifch,*)jcv
      call utmsgf
           endif
      nvq=nvq+jcv(i,1)
12    nva=nva+jcv(i,2)
      do 13 i=1,2
13    id(i)=nint(sil(4+i))
      call iddeco(id,jc)
      nqu=0
      naq=0
      do 11 i=1,nflav
      nqu=nqu+jc(i,1)
11    naq=naq+jc(i,2)
      neq=nqu-naq
 
c     determine nqaq,ival1,ival2
c     --------------------------
      nqaq=0
      ival=0
            if(rangen().gt.w)then
      pq=1
      pa=1
      if(neq.le.neqmn)pq=0
      if(neq.ge.neqmx)pa=0
      sum=(nvq*pq+nva*pa*iaqu)*qvalc+qseac*(.5*pq+.5*pa*iaqu)
      if(sum.eq.0.)goto14
      pvq=nvq*pq*qvalc/sum
      pva=nva*pa*iaqu*qvalc/sum
      psq=.5*pq*qseac/sum
      psa=.5*pa*iaqu*qseac/sum
      r=rangen()
      if(r.lt.pvq+pva)ival=1
      nqaq=-1
      if(ival.eq.1.and.r.lt.pvq)nqaq=1
      if(ival.eq.0.and.r.lt.pvq+pva+psq)nqaq=1
           if(nqu-naq-nqaq.lt.neqmn)then
      if(ish.ge.90)then
      call utmsg('hapapa')
      write(ifch,*)'*****  neq < neqmn'
      write(ifch,*)'pvq/a psq/a:',pvq,pva,psq,psa
      write(ifch,*)'nqaq:',nqaq
      call utmsgf
      endif
      nqaq=0
           endif
           if(nqu-naq-nqaq.gt.neqmx)then
      if(ish.ge.91)write(ifch,*)' '
      if(ish.ge.90)then
      call utmsg('hapapa')
      write(ifch,*)'*****  neq > neqmx'
      write(ifch,*)'pvq/a psq/a:',pvq,pva,psq,psa
      write(ifch,*)'nqaq:',nqaq
      endif
      if(ish.ge.91)write(ifch,*)' '
      nqaq=0
           endif
      if(nqu-naq-nqaq.lt.neqmn)
     *call utstop('hapapa: neq.lt.neqmn                    ')
      if(nqu-naq-nqaq.gt.neqmx)
     *call utstop('hapapa: neq.gt.neqmx                    ')
           endif
14    continue
      ival1=ival
      ival2=ival
           if(nqaq.eq.0)then
      ival1=0
      sum=nvq*qvalc+qseac*.5
      if(sum.eq.0.)goto15
      pvq=nvq*qvalc/sum
      psq=.5*qseac/sum
      r=rangen()
      if(r.lt.pvq)ival1=1
   15 continue
c-c   if(nvq.gt.0.and.rangen().lt.pvalen)ival1=1
      ival2=0
      sum=nva*iaqu*qvalc+qseac*.5*iaqu
      if(sum.eq.0.)goto16
      pva=nva*iaqu*qvalc/sum
      psa=.5*iaqu*qseac/sum
      r=rangen()
      if(r.lt.pva)ival2=1
   16 continue
c-c   if(nva.gt.0.and.rangen().lt.pvalen)ival2=1
      if(ival1.eq.1.and.ival2.eq.1)then
      r=rangen()
      if(r.lt..5)then
      ival1=0
      else
      ival2=0
      endif
      endif
           endif
 
c     quark
c     -----
      xfla='---'
           if(nqaq.ge.0)then
      su=0
      do 20 i=1,4
      probab(i)=0
      if(ival1.eq.1)probab(i)=jcv(i,1)
      if(ival1.eq.0.and.i.le.2)probab(i)=1
      if(ival1.eq.0.and.i.eq.3)probab(i)=rstras
20    su=su+probab(i)
      probsu(1)=0
      do 30 i=1,4
30    probsu(i+1)=probsu(i)+probab(i)/su
      r=rangen()
      do 50 i=1,4
      if(probsu(i).le.r.and.r.le.probsu(i+1))ifla=i
50    continue
      if(ifla.gt.4)ifla=4
      if(ival1.eq.1)xfla='va'//chpt
      if(ival1.eq.0)xfla='se'//chpt
      pta=utinvt(nptq,xptq,qptq,rangen()*qptq(nptq))
      ptah=utinvt(nptq,xptq,qpth,rangen()*qpth(nptq))
           endif
 
c     antiquark
c     ---------
      xflb='---'
           if(nqaq.le.0)then
      su=0
      do 60 i=1,4
      probab(i)=0
      if(ival2.eq.1)probab(i)=jcv(i,2)
      if(ival2.eq.0.and.i.le.2)probab(i)=1
      if(ival2.eq.0.and.i.eq.3)probab(i)=rstras
60    su=su+probab(i)
      probsu(1)=0
      do 70 i=1,4
70    probsu(i+1)=probsu(i)+probab(i)/su
      r=rangen()
      do 80 i=1,4
      if(probsu(i).le.r.and.r.le.probsu(i+1))iflb=i
80    continue
      if(iflb.gt.4)iflb=4
      if(nqaq.eq.0.and.ival1.eq.0.and.ival2.eq.0)iflb=ifla
c-c   if(nqaq.eq.0.and.ival2.eq.0)iflb=ifla
c-c   if(nqaq.eq.0.and.ival1.eq.0)ifla=iflb
      if(ival2.eq.1)xflb='va'//chpt
      if(ival2.eq.0)xflb='se'//chpt
      ptb=utinvt(nptq,xptq,qptq,rangen()*qptq(nptq))
      ptbh=utinvt(nptq,xptq,qpth,rangen()*qpth(nptq))
           endif
 
      if(nqaq.eq.0.and.ival1.eq.0.and.ival2.eq.0)ifla=-ifla
      if(nqaq.eq.0.and.ival1.eq.0.and.ival2.eq.0)iflb=-iflb
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hasi(qaq,sil,ifl,xfl,pt0,pt0h,si,iret,jord,ixflab
     *               ,ptdiff)
c     determines string ingredients (=jets)
c iret=0: ok
c iret=1: remnant changes direction
c iret=2: jc(,)=10
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cnsta/pi,prom,piom,ainfin
      character xfl*3,qaq*5
      parameter (nsi=6,nflav=6,nsix=100)
      real sil(nsi),si(nsi),six(nsi,nsix)
      integer jc(nflav,2),id(2),ic(2)
      real ptdiff(2)
      iret=0
      iretx=0
      ntry=0
      ptdiff(1)=0.
      ptdiff(2)=0.
      ifla=iabs(ifl)
           if(xfl.eq.'000')then
      do 11 j=3,nsi
11    si(j)=0
           else
      do 10 j=1,nsi
10    si(j)=0
           endif
      if(xfl.eq.'---')goto1000
      if(jord.eq.1.and.qaq.eq.'quark')ams=amproj
      if(jord.eq.2.and.qaq.eq.'antiq')ams=amproj
      if(jord.eq.2.and.qaq.eq.'quark')ams=amtarg
      if(jord.eq.1.and.qaq.eq.'antiq')ams=amtarg
c-c   enll=sqrt(ams**2+pnllx**2)
      enll=sil(4)
      sgn=sign(1.,sil(3))
      ihard=0
      if(rangen().le.phard)ihard=1
      do 1 i=1,2
1     id(i)=nint(sil(4+i))
      call iddeco(id,jc)
9999  ntry=ntry+1
      if(ntry.gt.20)then
      iret=1
      goto1000
      endif
           if(ihard.eq.0)then
      pt=pt0
      if(pt.gt.enll)pt=rangen()*enll
      enmin=0
           elseif(ihard.eq.1)then
      pt=pt0h
      if(pt.gt.enll)pt=rangen()*enll
      enmin=pt
           endif
      ptfull=pt
      lo=0
51    lo=lo+1
      if(xfl.eq.'000')then
      en=0.
      pt=0.
      else
      en=ranstc(xfl,enmin/enll)*enll
      endif
      if(pt.gt.en.and.rangen().lt.0.5.and.lo.le.10)goto51
      if(pt.gt.en)pt=rangen()*en
      phi=2.*pi*rangen()
      if(xfl.ne.'000')then
      si(1)=pt*cos(phi)
      si(2)=pt*sin(phi)
      si(3)=sgn*sqrt(en**2-pt**2)
      if(si(1).eq.0..and.si(2).eq.0..and.si(3).eq.0.)then
      iretx=1
      iret=1
      goto1000
      endif
c-c   ptdiff(1)=(ptfull-pt)*cos(phi)
c-c   ptdiff(2)=(ptfull-pt)*sin(phi)
      endif
      if((sil(3)-si(3))*sil(3).lt.0.)then
      if(ish.ge.94)write(ifch,*)'sil,si,enll,en',sil,si,enll,en
      goto9999
      endif
      si(4)=sqrt(si(1)**2+si(2)**2+si(3)**2)
      if(qaq.eq.'quark')si(5)=10**(nflav-ifla)
      if(qaq.eq.'antiq')si(6)=10**(nflav-ifla)
 
      if(ixflab.eq.0)then
      sil(1)=sil(1)-si(1)
      sil(2)=sil(2)-si(2)
      sil(3)=sil(3)-si(3)
      sil(4)=sqrt(sil(1)**2+sil(2)**2+sil(3)**2)
      endif
 
      l1=0
      l2=0
      if(jc(ifla,1).gt.0)l1=1
      if(jc(ifla,2).gt.0)l2=1
      if(qaq.eq.'quark'.and.l1.eq.0)jc(ifla,2)=jc(ifla,2)+1
      if(qaq.eq.'quark'.and.l1.eq.1)jc(ifla,1)=jc(ifla,1)-1
      if(qaq.eq.'antiq'.and.l2.eq.0)jc(ifla,1)=jc(ifla,1)+1
      if(qaq.eq.'antiq'.and.l2.eq.1)jc(ifla,2)=jc(ifla,2)-1
      if(jc(ifla,1).eq.10.or.jc(ifla,2).eq.10)goto9998
      call idenco(jc,ic,ireten)
      if(ireten.eq.1)
     *call utstop('hasi  : idenco ret code = 1             ')
      sil(5)=ic(1)
      sil(6)=ic(2)
      goto1000
 
9998  iret=2
      if(ish.ge.90)then
      call utmsg('hasi  ')
      write(ifch,*)'*****  jc(,)=10'
      write(ifch,*)jc
      call utmsgf
      endif
 
1000  continue
           if(iretx.eq.1)then
      call utmsg('hasi')
      write(ifch,*)'*****  si(1/2/3)=0'
      write(ifch,*)si
      call utmsgf
           endif
      if(ish.ge.93)then
      if(iret.ne.0.or.ntry.gt.1)
     *write(ifch,*)'iret=',iret,'  ntry=',ntry
      if(ntry.gt.0)
     *write(ifch,100)xfl,pt0,pt,en,en/enll
100   format(1x,'hasi: xfl=',a3
     *,' pt0=',e8.3,' pt=',e8.3,' en=',e8.3,' x=',e8.3)
      if(ntry.eq.0)
     *write(ifch,101)xfl
101   format(1x,'hasi: xfl=',a3)
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hastfc(sil,six,iret)
c     checks long strings
c     iret=0: ok
c     iret=1: mass .lt. minimal mass
c     iret=2: jc(,) .gt. 9
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common /cnce/nces,ncolex
      parameter (nsi=6,nflav=6,nsix=100)
      real sil(nsi),six(nsi,nsix),si(nsi)
      integer ic(2),jc(nflav,2),jcp(nflav,2)
 
      iret=0
 
      call utksix(six,kmax)
      if(kmax.lt.1)goto9000
 
      do 110 n=1,4
      si(n)=sil(n)
      do 110 k=1,kmax
110   si(n)=si(n)+six(n,k)
      s=si(4)**2-si(1)**2-si(2)**2-si(3)**2
 
      ic(1)=nint(sil(5))
      ic(2)=nint(sil(6))
      call iddeco(ic,jc)
      do 130 k=1,kmax
      ic(1)=nint(abs(six(5,k)))
      ic(2)=nint(abs(six(6,k)))
      call iddeco(ic,jcp)
      do 120 i=1,nflav
      do 120 j=1,2
120   jc(i,j)=jc(i,j)+jcp(i,j)
130   continue
 
      do 125 i=1,nflav
      if(jc(i,1)-jc(i,2).ge.10)iret=2
      if(jc(i,2)-jc(i,1).ge.10)iret=2
125   continue
 
      if(iret.eq.2)goto9000
 
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      amin=utamnu(keu,ked,kes,kec,4)/2
 
c-c   if(s.lt.amin**2)iret=1
 
9000  continue
 
           if(ish.ge.94)then
      write(ifch,*)' '
      if(iret.eq.0)
     *write(ifch,*)('-',l=1,79)
      if(iret.ne.0)
     *write(ifch,*)('#',l=1,79)
      write(ifch,*)'iret= ',iret ,'  kmax= ',kmax,'  nrevt= ',nrevt
     *,'  nces= ',nces
      write(ifch,8004)sil
8004  format(' sil: ',4f11.5,2f8.0)
           if(kmax.gt.0)then
      write(ifch,8007)(six(i,1),i=1,nsi)
8007  format(' six: ',4f11.5,2f8.0)
      if(kmax.gt.1)then
      do 1 j=2,kmax
1     write(ifch,8008)(six(i,j),i=1,nsi)
8008  format('      ',4f11.5,2f8.0)
      endif
      write(ifch,8005)(si(i),i=1,4)
8005  format('  si: ',4f11.5)
      write(ifch,8006)(si(i)**2,i=1,4)
8006  format(' si^2:',4f11.5)
      write(ifch,*)'jc:'
      write(ifch,*)jc
      write(ifch,*)'keu,ked,kes,kec: ',keu,ked,kes,kec
      write(ifch,*)'s= ',s,'    amin**2= ',amin**2,'   amin= ',amin
           endif
      if(iret.ne.0)
     *write(ifch,*)('#',l=1,79)
      if(iret.eq.0)
     *write(ifch,*)('-',l=1,79)
      write(ifch,*)' '
           endif
 
      if(iret.ne.0)return
 
      entry hastfw(sil,six)
      if(.not.(ish.ge.91))return
      call utksix(six,kmax)
      write(ifch,103)sil,((six(i,j),i=1,nsi),j=1,kmax)
103   format(/' sil: ',4f11.5,2f8.0/,' six: ',4f11.5,2f8.0/
     *,50('      ',4f11.5,2f8.0/))
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hastfl(sil,six,strl)
c     forms a leading string
c-----------------------------------------------------------------------
      parameter (nsi=6,nsix=100)
      real sil(nsi),six(nsi,nsix),strl(nsi,nsix+1)
      do 100 n=1,nsi
      strl(n,1)=sil(n)
100   strl(n,2)=0
      call utksix(six,kmax)
      if(kmax.eq.0)return
      do 110 k=1,kmax
      do 110 n=1,nsi
110   strl(n,k+1)=six(n,k)
      do 120 n=1,nsi
120   strl(n,kmax+2)=0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hastfr(sila,silb,sia,sib,iret)
c-----------------------------------------------------------------------
c  subtracts sia-sib from sila-silb (only momentum components)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax

      parameter (nsi=6)
      real sila(nsi),silb(nsi),sia(nsi),sib(nsi)
      double precision a,d,b,put,pux,puy,pot,pox,poy,pep,pem,pyp,pym
     *,pip,pim,pi3,pi4,pap,pam,pa3,pa4,sla,slb
      double precision ssila(nsi),ssilb(nsi),ssia(nsi),ssib(nsi)
     *,plab(4),pab(4)
 
      iret=0

      do i=1,4
      ssila(i)=sila(i)
      ssilb(i)=silb(i)
      ssia(i)=sia(i)
      ssib(i)=sib(i)
      plab(i)=ssila(i)+ssilb(i)
      pab(i)=ssia(i)+ssib(i)
      enddo
      sla=ssila(3)/dabs(ssila(3))
      slb=-sla
 
      if(ssia(3)*ssila(3).lt.0.)then
      call utmsg('hastfr')
      write(ifch,*)'*****  ssia(3)*ssila(3)<0'
      write(ifch,*)ssia(3),ssila(3)
      call utmsgf
      endif
 
      pep=0
      pem=0
      pot=0
      put=0
      a=0
      d=0
      pyp=0
      pym=0
      pip=0
      pam=0
 
      pox=ssila(1)-ssia(1)
      poy=ssila(2)-ssia(2)
      pot=dsqrt(pox**2+poy**2)
      pux=ssilb(1)-ssib(1)
      puy=ssilb(2)-ssib(2)
      put=dsqrt(pux**2+puy**2)
      pep= plab(4)+sla*plab(3) - (pab(4)+sla*pab(3))
      pem= plab(4)-sla*plab(3) - (pab(4)-sla*pab(3))
      if(pep.lt.0.)goto9001
      if(pem.lt.0.)goto9001
      a=(pem*pep-put**2-pot**2)/2.
      d=put*pot
      if(a**2-d**2.lt.0.)goto9001
      b=dsqrt(a**2-d**2)
      pyp=(a+put**2-b)/pem
      pym=(a+pot**2-b)/pep
      if(pyp.lt.0.)goto9001
      if(pym.lt.0.)goto9001
      pip=pep-pyp
      if(pip.lt.0.)goto9001
      pim=pym
      pap=pyp
      pam=pem-pym
      if(pam.lt.0.)goto9001
      pi3=(pip-pim)/2*sla
      pi4=(pip+pim)/2
      pa3=(pap-pam)/2*sla
      pa4=(pap+pam)/2
      if(pi3*ssila(3).lt.0.)goto9001
      if(pa3*ssilb(3).lt.0.)goto9001
      ssila(1)=pox
      ssila(2)=poy
      ssila(3)=pi3
      ssila(4)=pi4
      ssilb(1)=pux
      ssilb(2)=puy
      ssilb(3)=pa3
      ssilb(4)=pa4
 
      do i=1,4
      sila(i)=ssila(i)
      silb(i)=ssilb(i)
      sia(i)=ssia(i)
      sib(i)=ssib(i)
      enddo

      if(ish.lt.90)goto9000
 
      if(dabs(pip*pim-pot**2).gt.1e-4)then
      call utmsg('hastfr')
      write(ifch,*)'*****  pip*pim /= pot**2'
      write(ifch,*)'pip*pim=',pip*pim
      write(ifch,*)'pot**2=',pot**2
      write(ifch,*)'pip=',pip
      write(ifch,*)'pim=',pim
      write(ifch,*)'pot=',pot
      call utmsgf
      endif
      if(dabs(pap*pam-put**2).gt.1e-4)then
      call utmsg('hastfr')
      write(ifch,*)'*****  pap*pam /= put**2'
      write(ifch,*)'pap*pam=',pap*pam
      write(ifch,*)'put**2=',put**2
      write(ifch,*)'pap=',pap
      write(ifch,*)'pam=',pam
      write(ifch,*)'put=',put
      call utmsgf
      endif
      if(dabs(ssila(4)**2
     *-ssila(1)**2-ssila(2)**2-ssila(3)**2).gt.1e-4)then
      call utmsg('hastfr')
      write(ifch,*)'*****  mass**2 of ssila nonzero'
      write(ifch,*)'mass**2=',ssila(4)**2
     *-ssila(1)**2-ssila(2)**2-ssila(3)**2
      call utmsgf
      endif
      if(dabs(ssilb(4)**2
     *-ssilb(1)**2-ssilb(2)**2-ssilb(3)**2).gt.1e-4)then
      call utmsg('hastfr')
      write(ifch,*)'*****  mass**2 of ssilb nonzero'
      write(ifch,*)'mass**2=',ssilb(4)**2
     *-ssilb(1)**2-ssilb(2)**2-ssilb(3)**2
      call utmsgf
      endif
      do 14 n=1,4
      if(dabs(plab(n)
     *-ssila(n)-ssilb(n)-ssia(n)-ssib(n)).gt.1e-4)goto15
14    continue
      goto16
15    continue
      call utmsg('hastfr')
      write(ifch,*)'***** ssila+ssilb /= ssila_new+ssilb_new+ssia+ssib'
      write(ifch,*)'ssila_old+ssilb_old:'
      write(ifch,*)plab
      write(ifch,*)'ssila_new+ssilb_new+ssia+ssib:'
      write(ifch,*)((ssila(k)+ssilb(k)+ssia(k)+ssib(k)),k=1,4)
      write(ifch,*)'ssila_new:'
      write(ifch,*)(ssila(n),n=1,4)
      write(ifch,*)'ssilb_new:'
      write(ifch,*)(ssilb(n),n=1,4)
      write(ifch,*)'ssia:'
      write(ifch,*)(ssia(n),n=1,4)
      write(ifch,*)'ssib:'
      write(ifch,*)(ssib(n),n=1,4)
      call utmsgf
16    continue
 
9000  return
 
9001  iret=1
      if(ish.lt.90)return
      call utmsg('hastfr')
      write(ifch,*)'string subtraction not possible'
      write(ifch,*)'sngl(ssila):'
      write(ifch,*)(sngl(ssila(n)),n=1,4)
      write(ifch,*)'sngl(ssilb):'
      write(ifch,*)(sngl(ssilb(n)),n=1,4)
      write(ifch,*)'sngl(ssia):'
      write(ifch,*)(sngl(ssia(n)),n=1,4)
      write(ifch,*)'sngl(ssib):'
      write(ifch,*)(sngl(ssib(n)),n=1,4)
      write(ifch,*)'pep,pem:',pep,pem
      write(ifch,*)'pot,put:',pot,put
      write(ifch,*)'a,d:',a,d
      write(ifch,*)'a**2-d**2:',a**2-d**2
      write(ifch,*)'pyp,pym:',pyp,pym
      write(ifch,*)'pip,pam:',pip,pam
      write(ifch,*)'pi3,sila(3):',pi3,sila(3)
      write(ifch,*)'pa3,silb(3):',pa3,silb(3)
      call utmsgf
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hastfs(sila,sixa,ifla,xfla,pta,ptah
     *                 ,silb,sixb,iflb,xflb,ptb,ptbh  ,str,iret  ,jord)
c     forms a short (=q-qbar) string
c iret=0: ok
c iret=1: in sr hasi: remnant changes direction or zero si().
c iret=2: in sr hasi: jc(,)=10 .
c iret=3: string mass too small. not 4 , 5 .
c iret=4: string mass too small. equal flavour,zero momentum,xfl='000'.
c iret=5: string mass too small. valence quarks involved.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (nsi=6,nsix=100)
      real sila(nsi),sixa(nsi,nsix),silb(nsi),sixb(nsi,nsix)
     *,sia(nsi),sib(nsi),str(nsi,nsix+1),sila0(nsi),silb0(nsi),si(nsi)
     *,six(nsi),sts(nsi,2)
      character xfla*3,xflb*3
      integer ic(2)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      real ptdiff(2)
 
      do 50 n=1,nsi
      str(n,1)=0
      str(n,2)=0
      str(n,3)=0
50    continue
 
      call utksix(sixa,kmaxa)
      call utksix(sixb,kmaxb)
 
      pt=utinvt(nptq,xptq,qptq,rangen()*qptq(nptq))
      phi=2.*pi*rangen()
      sia(1)=pt*cos(phi)
      sia(2)=pt*sin(phi)
      sib(1)=-sia(1)
      sib(2)=-sia(2)
 
      ixflab=0
      if(  xfla.ne.'---'.and.xflb.ne.'---'
     *.and.(xfla.ne.'000'.or.xflb.ne.'000'))ixflab=1
      call hasi('quark',sila,ifla,xfla,pta,ptah,sia,iret,jord,ixflab
     *,ptdiff)
      if(iret.eq.1)goto9001
      if(iret.eq.2)goto9002
      sila(1)=sila(1)-ptdiff(1)
      sila(2)=sila(2)-ptdiff(2)
      silb(1)=silb(1)+ptdiff(1)
      silb(2)=silb(2)+ptdiff(2)
      call hasi('antiq',silb,iflb,xflb,ptb,ptbh,sib,iret,jord,ixflab
     *,ptdiff)
      if(iret.eq.1)goto9001
      if(iret.eq.2)goto9002
      silb(1)=silb(1)-ptdiff(1)
      silb(2)=silb(2)-ptdiff(2)
      sila(1)=sila(1)+ptdiff(1)
      sila(2)=sila(2)+ptdiff(2)
      if(sila(5).eq.0..and.sila(6).eq.0.)then
      sila(5)=100000
      sila(6)=100000
      endif
      if(silb(5).eq.0..and.silb(6).eq.0.)then
      silb(5)=100000
      silb(6)=100000
      endif
           if(xfla.eq.'000'.and.xflb.eq.'000')then
      goto9003
           elseif(xfla.ne.'---'.and.xflb.ne.'---')then
      call hastfr(sila,silb,sia,sib,iret)
      if(iret.eq.1)goto9003
      call utfost(sia,sib,sts)
      call utamst(sts,am,amin,iret)
      if(iret.ne.0)goto9003
      do 70 n=1,nsi
      str(n,1)=sia(n)
70    str(n,2)=sib(n)
           elseif(xfla.ne.'---'.and.xflb.eq.'---')then
      if(kmaxb+1.gt.nsix)
     *call utstop('hastfs: nsix too small                  ')
      do 60 n=1,nsi
      sixb(n,kmaxb+1)=sia(n)
      if(kmaxb+2.le.nsix)sixb(n,kmaxb+2)=0
60    continue
      if(ifla.lt.0)sixb(5,kmaxb+1)=-sixb(5,kmaxb+1)
      if(ifla.lt.0)sixb(6,kmaxb+1)=-sixb(6,kmaxb+1)
      kmaxb=kmaxb+1
           elseif(xfla.eq.'---'.and.xflb.ne.'---')then
      if(kmaxa+1.gt.nsix)
     *call utstop('hastfs: nsix too small                  ')
      do 80 n=1,nsi
      sixa(n,kmaxa+1)=sib(n)
      if(kmaxa+2.le.nsix)sixa(n,kmaxa+2)=0
80    continue
      if(iflb.lt.0)sixa(5,kmaxa+1)=-sixa(5,kmaxa+1)
      if(iflb.lt.0)sixa(6,kmaxa+1)=-sixa(6,kmaxa+1)
      kmaxa=kmaxa+1
           elseif(xfla.eq.'---'.and.xflb.eq.'---')then
c     no action
           else
      call utstop('hastfs: if/else error                   ')
           endif
 
9000  iret=0
      if(ish.ge.91)write(ifch,102)sia,xfla,sib,xflb
102   format(' sia: ',4f11.5,2f8.0,2x,a3
     *      /' sib: ',4f11.5,2f8.0,2x,a3)
      return
 
9001  iret=1
      if(ish.ge.91)write(ifch,102)sia,xfla,sib,xflb
      if(ish.ge.91)write(ifch,*)'remnant changes direction'
      if(ish.ge.91)write(ifch,*)' or zero si(1/2/3)'
      if(ish.ge.91)write(ifch,*)' '
      return
 
9002  iret=2
c     jc(,)=10
      if(ish.ge.91)write(ifch,102)sia,xfla,sib,xflb
      return
 
9003  iret=3
      if(xfla(1:2).eq.'va'.or.xflb(1:2).eq.'va')goto9005
      wt=abs(sia(1)+sib(1))+abs(sia(2)+sib(2))
     *+abs(sia(3)+sib(3))
      if(iabs(ifla).eq.iabs(iflb).and.wt.lt.1e-3
     *.and.xfla.eq.'000'.and.xflb.eq.'000')goto9004
      iflb=iabs(ifla)*iflb/iabs(iflb)
      if(ish.ge.91)write(ifch,102)sia,xfla,sib,xflb
      if(ish.ge.91)write(ifch,*)'string mass too small (3)'
      if(ish.ge.91)write(ifch,*)' '
      xfla='000'
      xflb='000'
      return
 
9004  iret=4
      if(ish.ge.91)write(ifch,102)sia,xfla,sib,xflb
      if(ish.ge.91)write(ifch,*)'string mass too small (4)'
      if(ish.ge.91)write(ifch,*)' '
      return
 
9005  iret=5
      if(ish.ge.91)write(ifch,102)sia,xfla,sib,xflb
      if(ish.ge.91)write(ifch,*)'string mass too small (5)'
      if(ish.ge.91)write(ifch,*)' '
      xfla(1:2)='se'
      xflb(1:2)='se'
      return
 
      end
 
c-----------------------------------------------------------------------
      subroutine hastpr(stro,isplt)
c-----------------------------------------------------------------------
c processes a string
c out: isplt=0 : processing of stro finished
c           >0 : hastpr to be redone to process reduced string
c           -1 : error
c           -3 : hastpr to be redone with new epart,apart
c           -4 : haha collision to be redone bec of too fast ld ptl
c                             (only for ichoic=2)
c           -5 : hastpr to be redone because of produced s=3/2 particle
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (nsi=6,nsix=100)
      real str(nsi,2),strp(nsi,2),stro(nsi,nsix+1)
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter (nflav=6)
      integer jc1(nflav,2),jc2(nflav,2),jc(nflav,2),ic1(2),ic2(2),ic(2)
     *,icx(2),icy(2),ic1x(2)
      parameter (ngau=129)
      common /cgau/xgau(ngau),qgau(ngau),wgau(ngau)
     */clead/lead,cooav3,cooav4 /ckol/kol
      real p1(4),p2(4)
      double precision ppx,ppy,ppz,ppe,ppm,pp1(4),pp2(4),ppp
      parameter (nha=nsix+2,mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      common/cjspli/jpart,nscc,nsccx,apart,epart,alead,elead,sgnsil
      common/cnsuc/nsuc
      common/cstsh/nstsh
      common/cpzstr/pzstrl,estrl,istrl,isea
      zero=0.
      if(isplt.eq.-3)then
      jpart=1
      isplt=0
      elseif(isplt.eq.0)then
      jpart=0
      endif
      isplt0=isplt
 
      call utkstr(stro,kmax)
 
c     zero string
c     -----------
      if(kmax.eq.0)then
      if(ish.ge.91)write(ifch,*)'zero string'
      if(ish.ge.91)write(ifch,*)' '
      isplt=0
      return
      endif
 
c     print
c     -----
      if(ish.ge.91)then
      write(ifch,105)(stro(i,1),i=1,4),(nint(stro(i,1)),i=5,6)
      if(kmax.gt.1)then
      do 8 k=2,kmax
8     write(ifch,104)(stro(i,k),i=1,4),(nint(stro(i,k)),i=5,6)
      endif
105   format(' str: ',4f11.5,2i8)
104   format('      ',4f11.5,2i8)
      write(ifch,*)' '
      endif
 
c     check leading object
c     --------------------
      do 3 i=1,2
3     ic1(i)=nint(stro(4+i,1))
      call idcomi(ic1,ic1x)
      do 4 i=1,2
4     stro(4+i,1)=ic1x(i)
      call iddeco(ic1x,jc1)
      npa=0
      do 7 n=1,nflav
      do 7 i=1,2
7     npa=npa+jc1(n,i)
 
c     split string
c     ------------
           if(kmax.gt.2.or.kmax.eq.2.and.isplt.gt.0
     *                 .or.kmax.eq.2.and.npa.gt.3
     *                 .or.kmax.eq.2.and.lead.eq.1)then
           if(isplt.eq.0)then
      nsuc=0
      kmaxor=kmax
      sgnsil=1
      if(stro(3,1).lt.0.)sgnsil=-1
      do 17 n=1,nflav
      do 17 i=1,2
17    jc2(n,i)=0
      do 15 k=1,kmax
      do 14 i=1,2
14    ic1(i)=nint(abs(stro(4+i,k)))
      call iddeco(ic1,jc1)
      do 16 n=1,nflav
      do 16 i=1,2
      jc2(n,i)=jc2(n,i)+jc1(n,i)
      if(n.gt.4.and.jc2(n,i).ne.0)
     *call utstop('hastpr: flavour > 4                     ')
16    continue
15    continue
      keu=jc2(1,1)-jc2(1,2)
      ked=jc2(2,1)-jc2(2,2)
      kes=jc2(3,1)-jc2(3,2)
      kec=jc2(4,1)-jc2(4,2)
      alead=utamnu(keu,ked,kes,kec,0)
c-c   alead=0
      elead=stro(4,1)
      if(jpart.eq.0)then
      nsccx=kmaxor-1
      nstr0=nstr
      else
      nsccx=max(1,nscc)
      nstr=nstr0
      endif
      apart=alead/nsccx
      epart=elead/nsccx
      if(alead.gt.elead)then
      if(ish.ge.90)then
      call utmsg('hastpr')
      write(ifch,*)'*****  alead > elead  ==> ret code = -1'
      write(ifch,*)'alead=',alead,'   elead=',elead
      call utmsgf
      endif
      isplt=-1
      return
      endif
           endif
      call jsplit(stro,str,kolsp,ier,kmaxor)
      if(ier.eq.1)then
      isplt=-1
      return
      endif
      if(ier.eq.3)then
      isplt=-3
      return
      endif
      isplt=isplt+1
      if(ier.eq.2)then
      if(ish.ge.91)write(ifch,*)'multistring: absorption of one leg'
      if(ish.ge.91)write(ifch,*)' '
      return
      endif
      kolz=kolsp
      if(kolz.lt.1.or.kolz.gt.kollmx)write(ifch,*)'kolz=',kolz
      if(kolz.lt.1.or.kolz.gt.kollmx)
     *call utstop('hastpr: kolz out of range (1)           ')
      kin=1
           else
      kin=0
      if(isplt.gt.0)kin=1
      do 5 n=1,nsi
      do 5 i=1,2
      str(n,i)=stro(n,i)
      if(n.ge.5)str(n,i)=abs(str(n,i))
5     continue
      if(isplt.gt.0)then
      kolz=kproja(2,iproj)
      if(sgnsil.lt.0.)kolz=ktarga(2,itarg)
      elseif(lead.eq.1.and.kmax.eq.2)then
      kolz=kproja(2,iproj)
      if(str(3,1).lt.0.)kolz=ktarga(2,itarg)
      else
      kolz=kol
      endif
      if(kolz.lt.1.or.kolz.gt.kollmx)then
      call utmsg('hastpr')
      write(ifch,*)'*****  kolz out of range;   kolz:',kolz
     *,'   set to:',kol
      write(ifch,105)(str(i,1),i=1,4),(nint(str(i,1)),i=5,6)
      write(ifch,104)(str(i,2),i=1,4),(nint(str(i,2)),i=5,6)
      call utmsgf
      kolz=kol
      endif
      isplt=0
           endif
 
c     add two jets
c     ------------
      do 1 i=1,2
      ic1(i)=nint(str(4+i,1))
1     ic2(i)=nint(str(4+i,2))
      ia1=ic1(1)+ic1(2)
      ia2=ic2(1)+ic2(2)
           if(ia1.eq.0.and.ia2.eq.0)then
      write(ifch,*)'stro:'
      write(ifch,105)(stro(i,1),i=1,4),(nint(stro(i,1)),i=5,6)
      if(kmax.gt.1)then
      do 18 k=2,kmax
18    write(ifch,104)(stro(i,k),i=1,4),(nint(stro(i,k)),i=5,6)
      endif
      write(ifch,*)'str:'
      write(ifch,105)(str(i,1),i=1,4),(nint(str(i,1)),i=5,6)
      write(ifch,104)(str(i,2),i=1,4),(nint(str(i,2)),i=5,6)
      call utstop('hastpr: str=0                           ')
           endif
      call iddeco(ic1,jc1)
      call iddeco(ic2,jc2)
      n1=0
      n2=0
      do 9 n=1,nflav
      do 9 i=1,2
      n1=n1+jc1(n,i)
      n2=n2+jc2(n,i)
      jc(n,i)=jc1(n,i)+jc2(n,i)
      if(n.gt.4.and.jc(n,i).ne.0)
     *call utstop('hastpr: flavour > 4                     ')
9     continue
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      amstr=utamnu(keu,ked,kes,kec,6)
      call idenco(jc,ic,ireten)
      if(ireten.eq.1)
     *call utstop('hastpr: idenco ret code = 1             ')
      idstr=idtra(ic,0,0,3)
 
c     determine rot
c     -------------
      do 10 i=1,4
      pp1(i)=str(i,1)
10    pp2(i)=str(i,2)
      ppx=pp1(1)+pp2(1)
      ppy=pp1(2)+pp2(2)
      ppz=pp1(3)+pp2(3)
      ppp=dsqrt(ppx**2+ppy**2+ppz**2)
           if(ia1.ne.0.and.ia2.ne.0)then
      pp1(4)=dsqrt(pp1(1)**2+pp1(2)**2+pp1(3)**2)
      pp2(4)=dsqrt(pp2(1)**2+pp2(2)**2+pp2(3)**2)
      ppe=pp1(4)+pp2(4)
      ppm=dsqrt((ppe-ppp)*(ppe+ppp))
      call utlob2(1,ppx,ppy,ppz,ppe,ppm,pp1(1),pp1(2),pp1(3),pp1(4))
      call utlob2(1,ppx,ppy,ppz,ppe,ppm,pp2(1),pp2(2),pp2(3),pp2(4))
           else
      ppe=pp1(4)+pp2(4)
      if((ppe-ppp)*(ppe+ppp).gt.0d0)then
      ppm=dsqrt((ppe-ppp)*(ppe+ppp))
      else
      ppm=0
      ppe=ppp
      endif
      do 12 i=1,4
      pp1(i)=0.
12    pp2(i)=0.
           endif
      px=ppx
      py=ppy
      pz=ppz
      e=ppe
      am=ppm
      p=ppp
      do 6 i=1,4
      p1(i)=pp1(i)
6     p2(i)=pp2(i)
      if(ish.ge.92)then
      write(ifch,107)(p1(i),i=1,4),ic1
      write(ifch,107)(p2(i),i=1,4),ic2
107   format(' jet: ',4f11.5,2i8)
      write(ifch,*)' '
      endif
      if(p1(3).ne.0.)then
      dircn=p1(3)/abs(p1(3))
      else
      dircn=1
      endif
 
c     mass adjustment
c     ---------------
      am0=am
      call idres(idstr,am,idstrr,iadj)
 
c diffractive particle:
c-c   if((ia1.eq.0.or.ia2.eq.0).and.isplt0.eq.0)then
c-c   if(mod(idstrr,10).eq.0.and.rangen().lt.0.1)then
c-c   am=0
c-c   call idres(idstr,am,idstrr,iadj)
c-c   am=am+.300
c-c   call idres(idstr,am,idstrr,iadj)
c-c   endif
c-c   endif
 
c remnant:
      if((ia1.eq.0.or.ia2.eq.0).and.isplt0.gt.0.and.nsuc.gt.0)then
      am=0
      call idres(idstr,am,idstrr,iadj)
      if(n1+n2.eq.2)then
      if(rangen().lt.0.33)am=am+0.400
      elseif(mod(idstrr,10).eq.0)then
      if(rangen().lt.0.33)am=am+0.300
      else
      if(rangen().gt.0.33)then
      isplt=-5
      return
      endif
      endif
      call idres(idstr,am,idstrr,iadj)
      endif
 
c all legs absorbed
      if((ia1.eq.0.or.ia2.eq.0).and.isplt0.gt.0.and.nsuc.eq.0)then
      if(isea.eq.0)then
      isplt=-4
      return
      endif
      am=0
      if(istrl.eq.1)am=sqrt( engy**2*(1-abs(pzstrl)/pnllx)+.94**2 )
      call idres(idstr,am,idstrr,iadj)
      if(istrl.ne.1)then
      if(mod(idstrr,10).eq.0)am=am+rangen()
      if(mod(idstrr,10).ne.0)am=am+rangen()-.300
c-c   if(mod(idstrr,10).eq.0.and.rangen().lt.0.33)am=am+0.300
      call idres(idstr,am,idstrr,iadj)
      endif
      endif
 
c leading string = hadron:
      if(lead.eq.1.and.ia1.ne.0.and.ia2.ne.0
     *.and.isplt0.eq.0.and.isplt.eq.0.and.idstrr.ne.0)then
      am=0
      call idres(idstr,am,idstrr,iadj)
      if(mod(idstrr,10).eq.0.and.rangen().lt.0.33)am=am+0.300
      call idres(idstr,am,idstrr,iadj)
      endif
 
      idstr=idstrr
      if(idstr.eq.0.and.am.lt.amstr)am=amstr
      ppm=am
      e=dsqrt(ppp**2+ppm**2)
 
c     write /cstr/
c     ------------
      nstr=nstr+1
      if(nstr.gt.mxstr)
     *call utstop('hastpr: nstr>mxstr                      ')
      if(lead.eq.0)nstsh=1
      nsuc=nsuc+1
      irlstr(nstr)=nsuc
      if(isplt0.eq.0.and.isplt.eq.0)irlstr(nstr)=0
      icstr(1,nstr)=ic1(1)
      icstr(2,nstr)=ic1(2)
      icstr(3,nstr)=ic2(1)
      icstr(4,nstr)=ic2(2)
      pstr(1,nstr)=px
      pstr(2,nstr)=py
      pstr(3,nstr)=pz
      pstr(4,nstr)=e
      pstr(5,nstr)=am
      rotstr(1,nstr)=p1(1)
      rotstr(2,nstr)=p1(2)
      rotstr(3,nstr)=p1(3)
      if(p1(1).eq.0..and.p1(2).eq.0..and.p1(3).eq.0.)
     *rotstr(3,nstr)=1.
      do 2 i=1,4
2     xorstr(i,nstr)=coord(i,kolz)
      iorstr(nstr)=-kolz
      amsac=amsac+am
      if(ish.ge.91)then
      if(isplt.gt.0)then
      write(ifch,*)'split off string:'
      elseif(isplt0.gt.0)then
      write(ifch,*)'remainder:'
      elseif(lead.eq.1)then
      write(ifch,*)'ordinary baryonic string:'
      else
      write(ifch,*)'ordinary mesonic string:'
      endif
      write(ifch,106)nstr,(icstr(k,nstr)/100,k=1,4)
     *,(pstr(i,nstr),i=3,5)
106   format(/' /cstr/',i4,3x,4i5,2x,3(1x,e9.3)/)
      endif
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hdecmp(bar,sil,six)
c     decomposes bar into sil,six
c-----------------------------------------------------------------------
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      real bar(nsi,nha),sil(nsi),six(nsi,nsix)
      do 1 n=1,nsi
1     sil(n)=bar(n,2)
      do 3 m=1,nsix
      sixsqr=0
      do 2 n=1,nsi
      six(n,m)=bar(n,2+m)
2     sixsqr=sixsqr+six(n,m)**2
      if(sixsqr.le.1e-5)return
3     continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine hrescl(j1,j2,psum,ifail)
c     rescales string momenta of strings j1-j2 to have total mom psum.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      dimension psum(5),beta(3)
      double precision ppsum(5),padd(5),pp(5)
      data errlim/.001/

      ifail=1

      do k=1,5
      ppsum(k)=psum(k)
      enddo

      if(j1.ge.j2)
     *call utstop('hrescl: j1 .ge. j2                      ')

      do 100 k=1,5
100   padd(k)=0.
      do 110 j=j1,j2
      do 110 k=1,5
      padd(k)=padd(k)+pstr(k,j)
110   continue
           if(padd(5).ge.ppsum(5))then
      if(ish.ge.90)then
      call utmsg('hrescl')
      write(ifch,*)'*****  sum of str masses .ge. psum(5)'
      do 1 j=j1,j2
1     write(ifch,109)j,(icstr(k,j)/100,k=1,4)
     *,(pstr(l,j),l=3,5)
109   format(' /cstr/',i4,3x,4i5,2x,3(1x,e9.3))
      write(ifch,*)'ppsum(345):',(sngl(ppsum(k)),k=3,5)
      call utmsgf
      endif
      return
           endif
      padd(5)=padd(4)**2-padd(1)**2-padd(2)**2-padd(3)**2
           if(padd(5).le.0)then
      do 2 j=j1,j2
2     write(ifch,108)j,(pstr(l,j),l=1,5)
108   format(' /cstr/',i4,2x,5(1x,e9.3))
      call utstop('hrescl: mass**2 of string-sum negative  ')
           endif
      padd(5)=sqrt(padd(5))
 
c     boost strings to rest
c     ---------------------
           do j=j1,j2
      do k=1,4
      pp(k)=pstr(k,j)
      enddo 
      call utlob2(1,padd(1),padd(2),padd(3),padd(4),padd(5)    
     *,pp(1),pp(2),pp(3),pp(4))
      do k=1,4
      pstr(k,j)=pp(k)
      enddo 
           enddo  
 
c     rescale momenta in rest frame
c     -----------------------------
      scal=1.
      do 301 ipass=1,200
      sum=0.
      do 310 j=j1,j2
      do 320 k=1,3
320   pstr(k,j)=scal*pstr(k,j)
      pstr(4,j)=sqrt(pstr(1,j)**2+pstr(2,j)**2+pstr(3,j)**2
     *+pstr(5,j)**2)
      sum=sum+pstr(4,j)
310   continue
      scal=psum(5)/sum
      if(abs(scal-1.).le.errlim) goto300
301   continue
      if(ish.ge.90)then
      call utmsg('hrescl')
      write(ifch,*)'*****  scal=',scal
      call utmsgf
      endif
300   continue
 
c     boost back with ppsum
c     ---------------------
           do j=j1,j2
      do k=1,4
      pp(k)=pstr(k,j)
      enddo 
      call utlob2(-1,ppsum(1),ppsum(2),ppsum(3),ppsum(4),ppsum(5)    
     *,pp(1),pp(2),pp(3),pp(4))
      do k=1,4
      pstr(k,j)=pp(k)
      enddo 
           enddo  

      ifail=0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idchrg(id,chrg)
c     computes charge of particle with ident code id
c     ichrg must be dimensioned nqlep+12
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      dimension ichrg(53),ifl(3)
      data ichrg/0,2,-1,-1,2,-1,2,-1,2,0,0,0,-3,0,-3,0,-3,0,-3,3*0
     *,2,-1,-1,2,-1,2,-1,2,0,0,0,-3,0,-3,0,-3,0,-3,3,0
     *,3,0,0,0,3,3,3,6,6,6,0/
      idabs=iabs(id)
      call idflav(id,ifl(1),ifl(2),ifl(3),jspin,index)
      if(idabs.lt.100) goto200
      isum=0
      do 100 i=1,3
      isum=isum+ichrg(iabs(ifl(i))+1)*isign(1,ifl(i))
  100 continue
      chrg=isum/3.
      return
200   chrg=ichrg(index+1)*isign(1,id)
      chrg=chrg/3.
      return
      end

c-----------------------------------------------------------------------
      subroutine idcomk(ic)
c     compactifies ic
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      call idcomp(ic,icx,jc,1)
      ic(1)=icx(1)
      ic(2)=icx(2)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idcomi(ic,icx)
c     compactifies ic
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      call idcomp(ic,icx,jc,1)
      return
      end
  
c-----------------------------------------------------------------------
      subroutine idcomj(jc)
c     compactifies jc
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      call idcomp(ic,icx,jc,2)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idcomp(ic,icx,jc,im)
c-----------------------------------------------------------------------
c     compactifies ic,jc
c     input: im (1 or 2)
c            ic (if im=1)
c            jc (if im=2) 
c     output: icx (if im=1)
c             jc  
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      if(im.eq.1)call iddeco(ic,jc)
      icx(1)=0
      icx(2)=0
           do n=1,nflav
           do j=1,2
      if(jc(n,j).ne.0)goto1
           enddo
           enddo
      return
1     l=0
           do n=1,nflav
      k=min0(jc(n,1),jc(n,2))
      jc(n,1)=jc(n,1)-k
      jc(n,2)=jc(n,2)-k
      if(jc(n,1).lt.0.or.jc(n,2).lt.0)
     *call utstop('idcomp: jc negative                     ')
      l=l+jc(n,1)+jc(n,2)
           enddo
           if(l.eq.0)then
      jc(1,1)=1              
      jc(1,2)=1
           endif
           if(im.eq.1)then
      call idenco(jc,icx,ireten)
      if(ireten.eq.1)
     *call utstop('idcomp: idenco ret code = 1             ')
           endif
      return
      end
 
c-----------------------------------------------------------------------
      integer function idcon(ifl1,ifl2,ifail)
c     constructs particle identification code for ifl1+ifl2
c-----------------------------------------------------------------------
      common/didib/ndidib
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      real pmix1(3,2),pmix2(3,2),pspin1(8)
      data pmix1/.25,.25,.5,0.,.5,1./,pmix2/.5,.5,1.,0.,0.,1./
      do 1 i=1,8
      if(i.le.3)pspin1(i)=pspinl
      if(i.gt.3)pspin1(i)=pspinh
1     continue
      ifail=0
 
c               meson:
                if(iabs(ifl1).le.6.and.iabs(ifl2).le.6)then
      if(ifl1*ifl2.gt.0)
     *call utstop('idcon : ifl1*ifl2.gt.0                  ')
      ihigh=max0(iabs(ifl1),iabs(ifl2))
      jspin=int(rangen()+pspin1(ihigh))
      id1=ifl1
      id2=ifl2
           if(id1+id2.eq.0)then
      rnd=rangen()
      id1=iabs(id1)
      id1=int(pmix1(id1,jspin+1)+rnd)+int(pmix2(id1,jspin+1)+rnd)+1
      id2=-id1
           endif
           if(iabs(id1).gt.iabs(id2))then
      isave=id1
      id1=id2
      id2=isave
           endif
      idcon=isign(100*iabs(id1)+10*iabs(id2)+jspin,id1)
 
c               baryon:
                elseif(mod(ifl1,100).eq.0.and.iabs(ifl2).le.6
     *             .or.mod(ifl2,100).eq.0.and.iabs(ifl1).le.6)then
           if(mod(ifl1,100).eq.0)then
      id3=mod(ifl1/100,10)
      id2=ifl1/1000
      id1=ifl2
           elseif(mod(ifl2,100).eq.0)then
      id3=mod(ifl2/100,10)
      id2=ifl2/1000
      id1=ifl1
           endif
           if(iabs(id1).gt.iabs(id2))then
      iswap=id1
      id1=id2
      id2=iswap
           endif
           if(iabs(id2).gt.iabs(id3))then
      iswap=id2
      id2=id3
      id3=iswap
           endif
           if(iabs(id1).gt.iabs(id2))then
      iswap=id1
      id1=id2
      id2=iswap
           endif
      jspin=1
           if(.not.(id1.eq.id2.and.id2.eq.id3))then
      ihigh=iabs(id3)
      jspin=int(rangen()+pspin1(ihigh))
           endif
           if(jspin.ne.1.and.id1.ne.id2.and.id2.ne.id3
     *      .and.rangen().le.pispn)then
      iswap=id1
      id1=id2
      id2=iswap
           endif
      idcon=1000*iabs(id1)+100*iabs(id2)+10*iabs(id3)+jspin
      idcon=isign(idcon,ifl1)
 
               else
      ifail=1
      r=rangen()
      if(r.le..5)idcon=121
      if(r.gt..5)idcon=-121
                endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine iddeco(ic,jc)
c     decode particle id
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)
      ici=ic(1)
      jc(6,1)=mod(ici,10)
      jc(5,1)=mod(ici/10,10)
      jc(4,1)=mod(ici/100,10)
      jc(3,1)=mod(ici/1000,10)
      jc(2,1)=mod(ici/10000,10)
      jc(1,1)=mod(ici/100000,10)
      ici=ic(2)
      jc(6,2)=mod(ici,10)
      jc(5,2)=mod(ici/10,10)
      jc(4,2)=mod(ici/100,10)
      jc(3,2)=mod(ici/1000,10)
      jc(2,2)=mod(ici/10000,10)
      jc(1,2)=mod(ici/100000,10)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idenco(jc,ic,ireten)
c     encode particle id
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)
      ireten=0
      ic(1)=0
      do 20 i=1,nflav
      if(jc(i,1).ge.10)goto22
20    ic(1)=ic(1)+jc(i,1)*10**(nflav-i)
      ic(2)=0
      do 21 i=1,nflav
      if(jc(i,2).ge.10)goto22
21    ic(2)=ic(2)+jc(i,2)*10**(nflav-i)
      return
22    ireten=1
      ic(1)=0
      ic(2)=0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idflav(id,ifl1,ifl2,ifl3,jspin,index)
c     unpacks the ident code id=+/-ijkl
c
c          mesons--
c          i=0, j<=k, +/- is sign for j
c          id=110 for pi0, id=220 for eta, etc.
c
c          baryons--
c          i<=j<=k in general
c          j<i<k for second state antisymmetric in (i,j), eg. l = 2130
c
c          other--
c          id=1,...,6 for quarks
c          id=9 for gluon
c          id=10 for photon
c          id=11,...,16 for leptons
c          id=20 for ks, id=-20 for kl
c
c          i=21...26 for scalar quarks
c          i=29 for gluino
c          i=30 for photino
c          i=31...36 for scalar leptons
c          i=39 for wino
c          i=40 for zino
c
c          id=80 for w+
c          id=81,...,89 for higgs mesons
c          id=90 for z0
c
c          diquarks--
c          id=+/-ij00, i<j for diquark composed of i,j.
c
c          index is a sequence number used internally
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      parameter ( nqlep=41,nmes=2,nbary=2)
      idabs=iabs(id)
      i=idabs/1000
      j=mod(idabs/100,10)
      k=mod(idabs/10,10)
      jspin=mod(idabs,10)
      if(id.ne.0.and.mod(id,100).eq.0) goto300
      if(j.eq.0) goto200
      if(i.eq.0) goto100
c          baryons
c          only x,y baryons are qqx, qqy, q=u,d,s.
      ifl1=isign(i,id)
      ifl2=isign(j,id)
      ifl3=isign(k,id)
      if(k.le.6) then
        index=max0(i-1,j-1)**2+i+max0(i-j,0)+(k-1)*k*(2*k-1)/6
     1  +109*jspin+36*nmes+nqlep+11
      else
        index=max0(i-1,j-1)**2+i+max0(i-j,0)+9*(k-7)+91
     1  +109*jspin+36*nmes+nqlep+11
      endif
      return
c          mesons
100   continue
      ifl1=0
      ifl2=isign(j,id)
      ifl3=isign(k,-id)
      index=j+k*(k-1)/2+36*jspin+nqlep
      index=index+11
      return
c          quarks, leptons, etc
200   continue
      ifl1=0
      ifl2=0
      ifl3=0
      jspin=0
      index=idabs
      if(idabs.lt.20) return
c          define index=20 for ks, index=21 for kl
      index=idabs+1
      if(id.eq.20) index=20
c          index=nqlep+1,...,nqlep+11 for w+, higgs, z0
      if(idabs.lt.80) return
      index=nqlep+idabs-79
      return
300   ifl1=isign(i,id)
      ifl2=isign(j,id)
      ifl3=0
      jspin=0
      index=0
      return
      end
 
c-----------------------------------------------------------------------
      function idlabl(id)
c     returns the character*8 label for the particle id
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      parameter ( nqlep=41,nmes=2,nbary=2)
c
      character*8 idlabl
      character*8 llep,lmes0,lmes1,lbar0,labar0,lbar1,labar1
      character*8 lqq,laqq
      dimension llep(104)
      dimension lmes0(64),lmes1(64)
      dimension lbar0(109),labar0(109),lbar1(109),labar1(109)
      dimension lqq(21),laqq(21)
c          diquark labels
      data lqq/
     1'uu0. ','ud0. ','dd0. ','us0. ','ds0. ','ss0. ','uc0. ','dc0. ',
     2'sc0. ','cc0. ','ub0. ','db0. ','sb0. ','cb0. ','bb0. ','ut0. ',
     3'dt0. ','st0. ','ct0. ','bt0. ','tt0. '/
      data laqq/
     1'auu0.','aud0.','add0.','aus0.','ads0.','ass0.','auc0.','adc0.',
     2'asc0.','acc0.','aub0.','adb0.','asb0.','acb0.','abb0.','aut0.',
     3'adt0.','ast0.','act0.','abt0.','att0.'/
c          quark and lepton labels
      data llep/
     *'     ','up   ','ub   ','dn   ','db   ','st   ','sb   ','ch   ',
     *'cb   ','bt   ','bb   ','tp   ','tb   ','y    ','yb   ','x    ',
     *'xb   ','gl   ','err  ','gm   ','err  ','nue  ','anue ','e-   ',
     *'e+   ','num  ','anum ','mu-  ','mu+  ','nut  ','anut ','tau- ',
     *'tau+ ','err  ','err  ','err  ','err  ','err  ','err  ','ks   ',
     *'err  ','err  ','kl   ',
     *'upss ','ubss ','dnss ','dbss ','stss ','sbss ','chss ','cbss ',
     *'btss ','bbss ','tpss ','tbss ','err  ','err  ','err  ','err  ',
     *'glss ','err  ','gmss ','err  ','ness ','aness','e-ss ','e+ss ',
     *'nmss ','anmss','mu-ss','mu+ss','ntss ','antss','t-ss ','t+ss ',
     *'err  ','err  ','err  ','err  ','w+ss ','w-ss ','z0ss ','err  ',
     *'w+   ','w-   ','h10  ','ah10 ','h20  ','ah20 ','h30  ','ah30 ',
     *'h4+  ','h4-  ','h5+  ','h5-  ','h6+  ','h6-  ','h7++ ','h7-- ',
     *'h8++ ','h8-- ','h9++ ','h9-- ','z0   '/
c          0- meson labels
      data lmes0/
     1'pi0  ','pi+  ','eta  ','pi-  ','k+   ','k0   ','etap ','ak0  ',
     2'k-   ','ad0  ','d-   ','f-   ','etac ','f+   ','d+   ','d0   ',
     2'ub.  ','db.  ','sb.  ','cb.  ','bb.  ','bc.  ','bs.  ','bd.  ',
     3'bu.  ','ut.  ','dt.  ','st.  ','ct.  ','bt.  ','tt.  ','tb.  ',
     4'tc.  ','ts.  ','td.  ','tu.  ','uy.  ','dy.  ','sy.  ','cy.  ',
     5'by.  ','ty.  ','yy.  ','yt.  ','yb.  ','yc.  ','ys.  ','yd.  ',
     6'yu.  ','ux.  ','dx.  ','sx.  ','cx.  ','bx.  ','tx.  ','yx.  ',
     7'xx.  ','xy.  ','xt.  ','xb.  ','xc.  ','xs.  ','xd.  ','xu.  '/
c          1- meson labels
      data lmes1/
     1'rho0 ','rho+ ','omeg ','rho- ','k*+  ','k*0  ','phi  ','ak*0 ',
     2'k*-  ','ad*0 ','d*-  ','f*-  ','jpsi ','f*+  ','d*+  ','d*0  ',
     3'ub*  ','db*  ','sb*  ','cb*  ','upsl ','bc*  ','bs*  ','bd*  ',
     4'bu*  ','ut*  ','dt*  ','st*  ','ct*  ','bt*  ','tt*  ','tb*  ',
     5'tc*  ','ts*  ','td*  ','tu*  ','uy*  ','dy*  ','sy*  ','cy*  ',
     6'by*  ','ty*  ','yy*  ','yt*  ','yb*  ','yc*  ','ys*  ','yd*  ',
     7'yu*  ','ux*  ','dx*  ','sx*  ','cx*  ','bx*  ','tx*  ','yx*  ',
     8'xx*  ','xy*  ','xt*  ','xb*  ','xc*  ','xs*  ','xd*  ','xu*  '/
c          1/2+ baryon labels
      data lbar0/
     1'err  ','p    ','n    ','err  ','err  ','s+   ','s0   ','s-   ',
     2'l    ','xi0  ','xi-  ','err  ','err  ','err  ','sc++ ','sc+  ',
     3'sc0  ','lc+  ','usc. ','dsc. ','ssc. ','sdc. ','suc. ','ucc. ',
     4'dcc. ','scc. ','err  ','err  ','err  ','err  ','uub. ','udb. ',
     5'ddb. ','dub. ','usb. ','dsb. ','ssb. ','sdb. ','sub. ','ucb. ',
     6'dcb. ','scb. ','ccb. ','csb. ','cdb. ','cub. ','ubb. ','dbb. ',
     7'sbb. ','cbb. ','err  ','err  ','err  ','err  ','err  ','utt. ',
     8'udt. ','ddt. ','dut. ','ust. ','dst. ','sst. ','sdt. ','sut. ',
     9'uct. ','dct. ','sct. ','cct. ','cst. ','cdt. ','cut. ','ubt. ',
     1'dbt. ','sbt. ','cbt. ','bbt. ','bct. ','bst. ','bdt. ','but. ',
     2'utt. ','dtt. ','stt. ','ctt. ','btt. ','err  ','err  ','err  ',
     3'err  ','err  ','err  ','uuy. ','udy. ','ddy. ','duy. ','usy. ',
     4'dsy. ','ssy. ','sdy. ','suy. ','uux. ','udx. ','ddx. ','dux. ',
     5'usx. ','dsx. ','ssx. ','sdx. ','sux. '/
      data labar0/
     1'err  ','ap   ','an   ','err  ','err  ','as-  ','as0  ','as+  ',
     2'al   ','axi0 ','axi+ ','err  ','err  ','err  ','asc--','asc- ',
     3'asc0 ','alc- ','ausc.','adsc.','assc.','asdc.','asuc.','aucc.',
     4'adcc.','ascc.','err  ','err  ','err  ','err  ','auub.','audb.',
     5'addb.','adub.','ausb.','adsb.','assb.','asdb.','asub.','aucb.',
     6'adcb.','ascb.','accb.','acsb.','acdb.','acub.','aubb.','adbb.',
     7'asbb.','acbb.','err  ','err  ','err  ','err  ','err  ','autt.',
     8'audt.','addt.','adut.','aust.','adst.','asst.','asdt.','asut.',
     9'auct.','adct.','asct.','acct.','acst.','acdt.','acut.','aubt.',
     1'adbt.','asbt.','acbt.','abbt.','abct.','abst.','abdt.','abut.',
     2'autt.','adtt.','astt.','actt.','abtt.','err  ','err  ','err  ',
     3'err  ','err  ','err  ','auuy.','audy.','addy.','aduy.','ausy.',
     4'adsy.','assy.','asdy.','asuy.','auux.','audx.','addx.','adux.',
     5'ausx.','adsx.','assx.','asdx.','asux.'/
c          3/2+ baryon labels
      data lbar1/
     1'dl++ ','dl+  ','dl0  ','dl-  ','err  ','s*+  ','s*0  ','s*-  ',
     2'err  ','xi*0 ','xi*- ','om-  ','err  ','err  ','uuc* ','udc* ',
     3'ddc* ','err  ','usc* ','dsc* ','ssc* ','err  ','err  ','ucc* ',
     4'dcc* ','scc* ','ccc* ','err  ','err  ','err  ','uub* ','udb* ',
     5'ddb* ','err  ','usb* ','dsb* ','ssb* ','err  ','err  ','ucb* ',
     6'dcb* ','scb* ','ccb* ','err  ','err  ','err  ','ubb* ','dbb* ',
     7'sbb* ','cbb* ','bbb* ','err  ','err  ','err  ','err  ','utt* ',
     8'udt* ','ddt* ','err  ','ust* ','dst* ','sst* ','err  ','err  ',
     9'uct* ','dct* ','sct* ','cct* ','err  ','err  ','err  ','ubt* ',
     1'dbt* ','sbt* ','cbt* ','bbt* ','err  ','err  ','err  ','err  ',
     2'utt* ','dtt* ','stt* ','ctt* ','btt* ','ttt* ','err  ','err  ',
     3'err  ','err  ','err  ','uuy* ','udy* ','ddy* ','err  ','usy* ',
     4'dsy* ','ssy* ','err  ','err  ','uux* ','udx* ','ddx* ','err  ',
     5'usx* ','dsx* ','ssx* ','err  ','err  '/
      data labar1/
     1'adl--','adl- ','adl0 ','adl+ ','err  ','as*- ','as*0 ','as*+ ',
     2'err  ','axi*0','axi*+','aom+ ','err  ','err  ','auuc*','audc*',
     3'addc*','err  ','ausc*','adsc*','assc*','err  ','err  ','aucc*',
     4'adcc*','ascc*','accc*','err  ','err  ','err  ','auub*','audb*',
     5'addb*','err  ','ausb*','adsb*','assb*','err  ','err  ','aucb*',
     6'adcb*','ascb*','accb*','err  ','err  ','err  ','aubb*','adbb*',
     7'asbb*','acbb*','abbb*','err  ','err  ','err  ','err  ','autt*',
     8'audt*','addt*','err  ','aust*','adst*','asst*','err  ','err  ',
     9'auct*','adct*','asct*','acct*','err  ','err  ','err  ','aubt*',
     1'adbt*','asbt*','acbt*','abbt*','err  ','err  ','err  ','err  ',
     2'autt*','adtt*','astt*','actt*','abtt*','attt*','err  ','err  ',
     3'err  ','err  ','err  ','auuy*','audy*','addy*','err  ','ausy*',
     4'adsy*','assy*','err  ','err  ','auux*','audx*','addx*','err  ',
     5'ausx*','adsx*','assx*','err  ','err  '/
c          entry
      call idflav(id,ifl1,ifl2,ifl3,jspin,index)
      if(iabs(id).lt.100) goto200
      if(iabs(id).lt.1000) goto100
      if(id.ne.0.and.mod(id,100).eq.0) goto300
c          baryons
      index=index-109*jspin-36*nmes-nqlep
      index=index-11
      if(jspin.eq.0.and.id.gt.0) idlabl=lbar0(index)
      if(jspin.eq.0.and.id.lt.0) idlabl=labar0(index)
      if(jspin.eq.1.and.id.gt.0) idlabl=lbar1(index)
      if(jspin.eq.1.and.id.lt.0) idlabl=labar1(index)
      return
c          mesons
100   continue
      i=max0(ifl2,ifl3)
      j=-min0(ifl2,ifl3)
      index=max0(i-1,j-1)**2+i+max0(i-j,0)
      if(jspin.eq.0) idlabl=lmes0(index)
      if(jspin.eq.1) idlabl=lmes1(index)
      return
c          quarks, leptons, etc.
200   continue
      index=2*index
      if(id.le.0) index=index+1
      idlabl=llep(index)
      return
300   i=iabs(ifl1)
      j=iabs(ifl2)
      index=i+j*(j-1)/2
      if(id.gt.0) idlabl=lqq(index)
      if(id.lt.0) idlabl=laqq(index)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idmass(id,amass)
c     returns the mass of the particle with ident code id.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      dimension ammes0(10),ammes1(10),ambar0(30),ambar1(30)
      dimension amlep(52)
      parameter ( nqlep=41,nmes=2,nbary=2)
      data amlep/.3,.3,.5,1.6,4.9,30.,-1.,-1.,0.,0.,
     *0.,.511003e-3,0.,.105661,0.,1.807,3*-1.,.49767,.49767,
     *100.3,100.3,100.5,101.6,104.9,130.,2*-1.,100.,0.,
     *100.,100.005,100.,100.1,100.,101.8,2*-1.,100.,100.,
     *11*0./
c          0- meson mass table
      data ammes0/.13496,.13957,.5488,.49367,.49767,.9576,1.8633
     1,1.8683,2.030,2.976/
c          1- meson mass table
      data ammes1/.770,.770,.7826,.8881,.8922,1.0196,2.006,2.0086
     1,2.140,3.097/
c          1/2+ baryon mass table
      data ambar0/-1.,.93828,.93957,2*-1.,1.1894,1.1925,1.1974
     1,1.1156,1.3149,1.3213,3*-1.,2.43,2.43,2.43,2.26
     2,2.50,2.50,2.60,2.40,2.40,3.55,3.55,3.70,4*-1./
c          3/2+ baryon mass table
      data ambar1/1.232,1.232,1.232,1.232,-1.,1.3823,1.3820
     1,1.3875,-1.,1.5318,1.5350,1.6722,2*-1.
     2,2.63,2.63,2.63,-1.,2.70,2.70,2.80,2*-1.,3.75,3.75
     3,3.90,4.80,3*-1./
c          entry
      call idflav(id,ifl1,ifl2,ifl3,jspin,index)
      if(id.ne.0.and.mod(id,100).eq.0) goto400
      if(iabs(ifl1).gt.4.or.iabs(ifl2).gt.4.or.iabs(ifl3).gt.4)
     1goto300
      if(ifl2.eq.0) goto200
      if(ifl1.eq.0) goto100
c          baryons
      index=index-109*jspin-36*nmes-nqlep
      index=index-11
      amass=(1-jspin)*ambar0(index)+jspin*ambar1(index)
      return
c          mesons
100   continue
      index=index-36*jspin-nqlep
      index=index-11
      amass=(1-jspin)*ammes0(index)+jspin*ammes1(index)
      return
c          quarks and leptons
200   continue
      amass=amlep(index)
      return
c          b and t particles
300   continue
      amass=amlep(iabs(ifl2))+amlep(iabs(ifl3))-.03+.04*jspin
      if(ifl1.ne.0) amass=amass+amlep(iabs(ifl1))
      return
c          diquarks
400   amass=amlep(iabs(ifl1))+amlep(iabs(ifl2))
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idmix(ic,jspin,icm,idm)
c     accounts for flavour mixing
c-----------------------------------------------------------------------
      parameter (nflav=6)
      real pmix1(3,2),pmix2(3,2)
      integer ic(2),icm(2)
      data pmix1/.25,.25,.5,0.,.5,1./,pmix2/.5,.5,1.,0.,0.,1./
      icm(1)=0
      icm(2)=0
      idm=0
      i=ic(1)
      if(i.ne.ic(2))return
      id=0
      if(i.eq.100000)id=1
      if(i.eq. 10000)id=2
      if(i.eq.  1000)id=3
      if(id.eq.0)return
      rnd=rangen()
      idm=int(pmix1(id,jspin+1)+rnd)+int(pmix2(id,jspin+1)+rnd)+1
      icm(1)=10**(nflav-idm)
      icm(2)=ic(1)
      idm=idm*100+idm*10+jspin
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idquac(i,nq,ns,na,jc)
c     returns quark content of ptl i from /cptl/ .
c        nq = # quarks - # antiquarks
c        ns = # strange quarks - # strange antiquarks
c        na = # quarks + # antiquarks
c        jc(nflav,2) = jc-type particle identification code.
c-----------------------------------------------------------------------
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)
 
      if(iabs(idptl(i)).eq.20)then
      idptl(i)=230
      if(rangen().lt..5)idptl(i)=-230
      goto9999
      endif
 
      if(iabs(idptl(i)).lt.100)then
      nq=0
      ns=0
      do 1 n=1,nflav
      jc(n,1)=0
1     jc(n,2)=0
      return
      endif
 
9999  if(idptl(i)/10**8.ne.7)then
      call idtr4(idptl(i),ic)
      call iddeco(ic,jc)
      else
      call idtrb(ibptl(1,i),ibptl(2,i),ibptl(3,i),ibptl(4,i),jc)
      endif
      na=0
      nq=0
      do 53 n=1,nflav
      na=na+jc(n,1)+jc(n,2)
53    nq=nq+jc(n,1)-jc(n,2)
      ns=   jc(3,1)-jc(3,2)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idres(id,am,idr,iadj)
c     returns resonance id idr corresponding to mass am.
c     performs mass adjustment, if necessary (if so iadj=1, 0 else)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
 
      if(am.eq.0.)am=1e-5
 
      idi=id
      ami=am
      idr=0
      iadj=0
 
      if(id.eq.0)goto9999
 
      do 3 k=1,mxmx
      do 3 m=2,mxma
      if(iabs(id).eq.idmx(m,k))id=idmx(1,k)*10*id/iabs(id)
3     continue
 
      ix=iabs(id)/10
      if(ix.lt.1.or.ix.gt.mxindx)
     *call utstop('idres: ix out of range.                 ')
      i=indx(ix)
      if(i.lt.1.or.i.gt.mxre)
     *call utstop('idres: particle not in table            ')
      do 1 j=1,mxma-1
      if(am.ge.rema(i,j).and.am.le.rema(i,j+1))then
      if(j-1.gt.9)
     *call utstop('idres: spin > 9                         ')
      idr=id/10*10+(j-1)*id/iabs(id)
      goto2
      endif
1     continue
      goto9999
2     continue
 
      do 4 k=1,mxmx
      if(ix.eq.idmx(1,k))then
      if(j.lt.1.or.j.gt.mxma-1)
     *call utstop('idres: index j out of range             ')
      if(idmx(j+1,k).ne.0)idr=idmx(j+1,k)*id/iabs(id)
      endif
4     continue
 
      iy=mod(iabs(idr),10)
      if(iy.gt.maxres)then
      iadj=0
      idr=0
      goto9999
      endif
 
      if(iy.ne.0.and.iy.ne.1)goto9999
 
      call idmass(idr,am)
      if(am.lt.0.)then
      write(ifch,*)' '
      write(ifch,*)'*****  error in idres: '
     *,'neg mass returned from idmass'
      write(ifch,*)'id,am(input):',idi,ami
      write(ifch,*)'idr,am:',idr,am
      call utstop('idres: neg mass returned from idmass    ')
      endif
      iadj=1
 
9999  id=idi
      if(.not.(ish.ge.93))return
      write(ifch,*)' '
      write(ifch,*)'return from idres. id,ami,am,idr,iadj:'
      write(ifch,*)id,ami,am,idr,iadj
      write(ifch,*)' '
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idresi
c     initializes /crema/
c-----------------------------------------------------------------------
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
      parameter (n=29)
      dimension remai(n,mxma),rewii(n,mxma),idmxi(mxma,mxmx)
     *,icrei(n,2*mxma)
      do 3 i=1,mxindx
3     indx(i)=0
      do 4 k=1,mxre
      do 4 m=1,mxma
4     rema(k,m)=0
 
      data (idmxi(j,1),j=1,mxma)/ 11, 110, 111,   0,   0,   0,   0, 4*0/
      data (idmxi(j,2),j=1,mxma)/ 22, 220, 330, 331,   0,   0,   0, 4*0/
      data (idmxi(j,3),j=1,mxma)/123,2130,1230,1231,   0,   0,   0, 4*0/
      data (idmxi(j,4),j=1,mxma)/124,2140,1240,1241,   0,   0,   0, 4*0/
      data (idmxi(j,5),j=1,mxma)/134,3140,1340,1341,   0,   0,   0, 4*0/
      data (idmxi(j,6),j=1,mxma)/234,3240,2340,2341,   0,   0,   0, 4*0/
 
      data ((icrei(k,m),m=1,2*mxma),k=1,10)/
     *111,000000, 9*300000,    11*0,
     *222,000000, 9*030000,    11*0,
     *112,       10*210000,    11*0,
     *122,       10*120000,    11*0,
     *113,       10*201000,    11*0,
     *223,       10*021000,    11*0,
     *123,       10*111000,    11*0,
     *133,       10*102000,    11*0,
     *233,       10*012000,    11*0,
     *333,000000, 9*003000,    11*0/
      data ((icrei(k,m),m=1,2*mxma),k=11,20)/
     *114,       10*200100,    11*0,
     *124,       10*110100,    11*0,
     *224,       10*020100,    11*0,
     *134,       10*101100,    11*0,
     *234,       10*011100,    11*0,
     *334,       10*002100,    11*0,
     *144,       10*100200,    11*0,
     *244,       10*010200,    11*0,
     *344,       10*001200,    11*0,
     *444,000000, 9*000300,    11*0/
      data ((icrei(k,m),m=1,2*mxma),k=21,29)/
     * 11,  10*100000,    0,   10*100000,
     * 22,  10*001000,    0,   10*001000,
     * 12,  10*100000,    0,   10*010000,
     * 13,  10*100000,    0,   10*001000,
     * 23,  10*010000,    0,   10*001000,
     * 14,  10*100000,    0,   10*000100,
     * 24,  10*010000,    0,   10*000100,
     * 34,  10*001000,    0,   10*000100,
     * 44,  10*000100,    0,   10*000100/
 
      data ((remai(k,m),m=1,mxma),k=1,10)/
     *111.,0.000,1.425,1.660,1.825,2.000,0.000,0.000,0.000,0.000,0.000,
     *222.,0.000,1.425,1.660,1.825,2.000,0.000,0.000,0.000,0.000,0.000,
     *112.,1.080,1.315,1.485,1.575,1.645,1.685,1.705,1.825,2.000,0.000,
     *122.,1.080,1.315,1.485,1.575,1.645,1.685,1.705,1.825,2.000,0.000,
     *113.,1.300,1.500,1.700,1.850,2.000,0.000,0.000,0.000,0.000,0.000,
     *223.,1.300,1.500,1.700,1.850,2.000,0.000,0.000,0.000,0.000,0.000,
     *123.,1.117,1.300,1.395,1.465,1.540,1.655,1.710,1.800,1.885,2.000,
     *133.,1.423,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *233.,1.428,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *333.,0.000,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/
      data ((remai(k,m),m=1,mxma),k=11,20)/
     *114.,2.530,2.730,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *124.,2.345,2.530,2.730,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *224.,2.530,2.730,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *134.,2.450,2.600,2.800,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *234.,2.450,2.600,2.800,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *334.,2.700,2.900,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *144.,3.650,3.850,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *244.,3.650,3.850,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *344.,3.800,4.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *444.,0.000,5.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/
      data ((remai(k,m),m=1,mxma),k=21,29)/
     * 11.,0.450,0.950,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 22.,0.750,0.965,1.080,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 12.,0.450,0.950,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 13.,0.700,1.050,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 23.,0.700,1.050,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 14.,1.935,2.077,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 24.,1.938,2.079,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 34.,2.085,2.195,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 44.,3.037,3.158,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/
 
      data ((rewii(k,m),m=1,mxma),k=1,5)/
     *111.,0.000e+00,0.115e+00,0.140e+00,0.250e+00,0.250e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *222.,0.000e+00,0.115e+00,0.140e+00,0.250e+00,0.250e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *112.,0.000e+00,0.115e+00,0.200e+00,0.140e+00,0.140e+00,
     *     0.145e+00,0.250e+00,0.140e+00,0.250e+00,0.000e+00,
     *122.,0.000e+00,0.115e+00,0.200e+00,0.140e+00,0.140e+00,
     *     0.145e+00,0.250e+00,0.140e+00,0.250e+00,0.000e+00,
     *113.,0.824e-14,0.036e+00,0.080e+00,0.100e+00,0.170e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00/
      data ((rewii(k,m),m=1,mxma),k=6,10)/
     *223.,0.445e-14,0.039e+00,0.080e+00,0.100e+00,0.170e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *123.,0.250e-14,0.890e-05,0.036e+00,0.040e+00,0.016e+00,
     *     0.090e+00,0.080e+00,0.100e+00,0.145e+00,0.170e+00,
     *133.,0.227e-14,0.009e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *233.,0.400e-14,0.010e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *333.,0.000e+00,0.800e-14,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00/
      data ((rewii(k,m),m=1,mxma),k=11,15)/
     *114.,0.400e-11,0.010e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *124.,0.400e-11,0.400e-11,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *224.,0.400e-11,0.010e+00,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *134.,0.150e-11,0.400e-11,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *234.,0.150e-11,0.400e-11,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00/
      data ((rewii(k,m),m=1,mxma),k=16,20)/
     *334.,0.400e-11,0.010e+00,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *144.,0.400e-11,0.010e+00,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *244.,0.400e-11,0.010e+00,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *344.,0.400e-11,0.010e+00,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *444.,0.400e-11,0.010e+00,0.010e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00/
      data ((rewii(k,m),m=1,mxma),k=21,25)/
     * 11.,0.757e-08,0.153e+00,0.057e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 22.,0.105e-05,0.210e-03,0.034e+00,0.004e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 12.,0.000e+00,0.153e+00,0.057e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 13.,0.000e+00,0.051e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 23.,0.197e-02,0.051e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00/
      data ((rewii(k,m),m=1,mxma),k=26,29)/
     * 14.,0.154e-11,0.002e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 24.,0.615e-12,0.002e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 34.,0.150e-11,0.020e+00,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     * 44.,0.010e+00,0.068e-03,0.000e+00,0.000e+00,0.000e+00,
     *     0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00/
 
      do 2 j=1,mxma
      do 2 i=1,mxmx
2     idmx(j,i)=idmxi(j,i)
 
      if(n.gt.mxre)
     *call utstop('idresi: dimension mxre too small        ')
      do 1 k=1,n
      ix=nint(remai(k,1))
      ix2=nint(rewii(k,1))
      ix3=icrei(k,1)
      if(ix.ne.ix2)
     *call utstop('idresi: ix /= ix2                       ')
      if(ix.ne.ix3)
     *call utstop('idresi: ix /= ix3                       ')
      if(ix.lt.1.or.ix.gt.mxindx)
     *call utstop('idresi: ix out of range.                ')
      indx(ix)=k
      rema(k,1)=0.
      rewi(k,1)=0.
      icre1(k,1)=0
      icre2(k,1)=0
      do 1 m=2,mxma
      rema(k,m)=remai(k,m)
      rewi(k,m)=rewii(k,m)
      icre1(k,m)=icrei(k,m)
1     icre2(k,m)=icrei(k,mxma+m)
 
      indx(33) =indx(22)
      indx(213)=indx(123)
      indx(214)=indx(124)
      indx(314)=indx(134)
      indx(324)=indx(234)
 
      return
      end
 
c-----------------------------------------------------------------------
      integer function idsgl(ic,gen,cmp)
c     returns 1 for singlets (qqq or qqbar) 0 else.
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),jcx(nflav,2),icx(2)
      character gen*6,cmp*6
      idsgl=0
      if(cmp.eq.'cmp-ys')then
      call idcomi(ic,icx)
      else
      icx(1)=ic(1)
      icx(2)=ic(2)
      endif
      call iddeco(icx,jcx)
      nq=0
      na=0
      do 1 i=1,nflav
      nq=nq+jcx(i,1)
1     na=na+jcx(i,2)
      if(nq.eq.0.and.na.eq.0)return
      if(gen.eq.'gen-no')then
      if(nq.eq.3.and.na.eq.0.or.nq.eq.1.and.na.eq.1
     *.or.nq.eq.0.and.na.eq.3)idsgl=1
      elseif(gen.eq.'gen-ys')then
      if(mod(nq-na,3).eq.0)idsgl=1
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idspin(ii,ic,id,jspin)
c     determines particle spin
c-----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      parameter (nflav=6)
      integer ic(2),jc(nflav,2)
      real pspin1(8)
      do 3 i=1,8
      if(i.le.3)pspin1(i)=pspinl
      if(i.gt.3)pspin1(i)=pspinh
3     continue
      call iddeco(ic,jc)
      ihigh=0
      if(ii.eq.1)then
      do 4 i=1,nflav
      if(jc(i,1).ne.0)ihigh=i
      if(jc(i,2).ne.0)ihigh=i
4     continue
      else
      call idflav(id,i1,i2,i3,idu1,idu2)
      do 5 i=1,nflav
      if(i.eq.iabs(i1).or.i.eq.iabs(i2).or.i.eq.iabs(i3))ihigh=i
5     continue
      endif
      jspin=int(rangen()+pspin1(ihigh))
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idtau(id,p4,p5,taugm)
c     returns lifetime*gamma for id with energy p4, mass p5
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
           if(iabs(id).lt.100.and.id.ne.20)then
      wi=0
           elseif(id.eq.20)then
      wi=.197/2.675e13
           elseif(iabs(id).lt.1e8)then
      ix=iabs(id)/10
      if(ix.lt.1.or.ix.gt.mxindx)
     *call utstop('idtau: ix out of range.                 ')
      ii=indx(ix)
      jj=mod(iabs(id),10)+2
      do 75 imx=1,mxmx
      do 75 ima=2,mxma
      if(iabs(id).eq.idmx(ima,imx))jj=ima
75    continue
      if(ii.lt.1.or.ii.gt.mxre.or.jj.lt.1.or.jj.gt.mxma)then
      write(ifch,*)' '
      write(ifch,*)'id,ii,jj:',id,'   ',ii,jj
      call utstop('idtau: ii or jj out of range            ')
      endif
      wi=rewi(ii,jj)
           else
      tauz=taunll
c-c   tauz=amin1(9./p5**2,tauz)
c-c   tauz=amax1(.2,tauz)
      wi=.197/tauz
           endif
      if(wi.eq.0.)then
      tau=ainfin
      else
      tau=.197/wi
      endif
      if(p5.ne.0.)then
      gm=p4/p5
      else
      gm=ainfin
      endif
      if(tau.ge.ainfin.or.gm.ge.ainfin)then
      taugm=ainfin
      else
      taugm=tau*gm
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idtr4(id,ic)
c     transforms generalized paige_id -> werner_id  (for < 4 flv)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
      integer ic(2)
           if(iabs(id).lt.1e8)then
      ix=iabs(id)/10
      if(ix.lt.1.or.ix.gt.mxindx)goto9999
      ii=indx(ix)
      if(ii.eq.0)goto9998
      jj=mod(iabs(id),10)+2
      do 27 imx=1,mxmx
      do 27 ima=2,mxma
      if(iabs(id).eq.idmx(ima,imx))jj=ima
27    continue
      if(id.gt.0)then
      ic(1)=icre1(ii,jj)
      ic(2)=icre2(ii,jj)
      else
      ic(2)=icre1(ii,jj)
      ic(1)=icre2(ii,jj)
      endif
      if(ic(1).eq.100000.and.ic(2).eq.100000.and.rangen().lt.0.5)then
      ic(1)=010000
      ic(2)=010000
      endif
           elseif(mod(id/10**8,10).eq.8)then
      ic(1)=mod(id,10**8)/10000*100
      ic(2)=mod(id,10**4)*100
           else
      write(ifch,*)' '
      write(ifch,*)'***** id: ',id
      call utstop('idtr4: unrecognized id                  ')
           endif
      return
 
9998  write(ifch,*)' '
      write(ifch,*)'id: ',id
      call utstop('idtr4: indx=0.                          ')
 
9999  write(ifch,*)' '
      write(ifch,*)'id: ',id
      call utstop('idtr4: ix out of range.                 ')
      end
 
c-----------------------------------------------------------------------
      integer function idtra(ic,ier,ires,imix)
c     tranforms from werner-id to paige-id
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (nidt=44,nflav=6)
      integer idt(3,nidt),ic(2),jc(nflav,2),icm(2)
      data idt/
     * 100000,000000,   1   ,010000,000000,   2   ,001000,000000,   3
     *,000100,000000,   4
     *,200000,000000,1100   ,110000,000000,1200   ,020000,000000,2200
     *,101000,000000,1300   ,011000,000000,2300   ,002000,000000,3300
     *,100100,000000,1400   ,010100,000000,2400   ,001100,000000,3400
     *,000200,000000,4400
     *,100000,100000, 110   ,100000,010000, 120   ,010000,010000, 220
     *,100000,001000, 130   ,010000,001000, 230   ,001000,001000, 330
     *,100000,000100, 140   ,010000,000100, 240   ,001000,000100, 340
     *,000100,000100, 440
     *,300000,000000,1111   ,210000,000000,1120   ,120000,000000,1220
     *,030000,000000,2221   ,201000,000000,1130   ,111000,000000,1230
     *,021000,000000,2230   ,102000,000000,1330   ,012000,000000,2330
     *,003000,000000,3331   ,200100,000000,1140   ,110100,000000,1240
     *,020100,000000,2240   ,101100,000000,1340   ,011100,000000,2340
     *,002100,000000,3340   ,100200,000000,1440   ,010200,000000,2440
     *,001200,000000,3440   ,000300,000000,4441/
 
      idtra=0
      if(ic(1).eq.0.and.ic(2).eq.0)return
      do 1 i=1,nidt
      if(ic(2).eq.idt(1,i).and.ic(1).eq.idt(2,i))idtra=-idt(3,i)
      if(ic(1).eq.idt(1,i).and.ic(2).eq.idt(2,i))idtra=idt(3,i)
1     continue
      isi=1
      if(idtra.ne.0)isi=idtra/iabs(idtra)
 
      jspin=0
      if(idtra.ne.0.and.mod(idtra,10).eq.0.and.ires.eq.1)
     *call idspin(1,ic,idu,jspin)
 
      if(imix.eq.1)then
      call idmix(ic,jspin,icm,idtram)
      if(idtram.ne.0)idtra=idtram
      if(jspin.eq.0)then
      if(iabs(idtra).eq.1230.and.rangen().lt.pispn)idtra=isi*2130
      if(iabs(idtra).eq.1240.and.rangen().lt.pispn)idtra=isi*2140
      if(iabs(idtra).eq.1340.and.rangen().lt.pispn)idtra=isi*3140
      if(iabs(idtra).eq.2340.and.rangen().lt.pispn)idtra=isi*3240
      endif
      elseif(imix.eq.2)then
      if(idtra.eq.220)idtra=110
      if(idtra.eq.330)idtra=110
      elseif(imix.eq.3)then
      if(idtra.eq.220)idtra=110
      if(idtra.eq.330)idtra=220
      endif
 
      if(idtra.ne.0)idtra=idtra+jspin*isi
 
      if(idtra.ne.0)return
      if(ier.ne.1)return
      jerr=jerr+1
      write(ifch,*)' '
      write(ifch,*)'***** error in idtra: unknown code'
      write(ifch,*)'ic = ',ic
      write(ifch,*)' '
      return
 
      entry idtrai(num,id,ier)
      idtrai=0
      j=0
      do 2 i=1,nidt
      if(iabs(id).eq.idt(3,i))j=i
2     continue
      if(j.ne.0)then
      if(id.lt.0)then
      idtrai=idt(3-num,j)
      else
      idtrai=idt(num,j)
      endif
      return
      endif
      if(ier.ne.1)return
      jerr=jerr+1
      write(ifch,*)' '
      write(ifch,*)'***** error in idtrai: unknown code'
      write(ifch,*)'id = ',id
      write(ifch,*)' '
      return
      end
 
c-----------------------------------------------------------------------
      subroutine idtrb(ib1,ib2,ib3,ib4,jc)
c     id transformation ib -> jc
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      parameter (nflav=6)
      integer ib(4),jc(nflav,2)
      jc(1,1)=ib1/10**4
      jc(2,1)=ib2/10**4
      jc(3,1)=ib3/10**4
      jc(4,1)=ib4/10**4
      jc(5,1)=0
      jc(6,1)=0
      jc(1,2)=mod(ib1,10**4)
      jc(2,2)=mod(ib2,10**4)
      jc(3,2)=mod(ib3,10**4)
      jc(4,2)=mod(ib4,10**4)
      jc(5,2)=0
      jc(6,2)=0
      return
      entry idtrbi(jc,ib1,ib2,ib3,ib4)
      ib1=jc(1,1)*10**4+jc(1,2)
      ib2=jc(2,1)*10**4+jc(2,2)
      ib3=jc(3,1)*10**4+jc(3,2)
      ib4=jc(4,1)*10**4+jc(4,2)
      ib5=jc(5,1)*10**4+jc(5,2)
      ib6=jc(6,1)*10**4+jc(6,2)
      if(ib5.ne.0.or.ib6.ne.0)then
      write(ifch,*)' '
      write(ifch,*)'***** error in idtrbi: bottom or top quarks'
      write(ifch,*)'jc:'
      write(ifch,*)jc
      call utstop('idtrbi: bottom or top quarks            ')
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jambr1(ns,nsg,iret)
c-----------------------------------------------------------------------
c "AMOR" (koch+werner, 89)
c breaks string ns according to a-m model.
c ns: points to the current fragmenting string.
c nsg: points to the last produced substring (succ increased)
c     input:
c psg(,ns): momentum of string ns in pp-cms
c pjt(,2*ns-1),pjt(,2*ns): momentum of end of string ns in pp-cms
c     output:
c xorsg(,nsg): origin of substring nsg in pp-cms
c pjt(,2*nsg-1),pjt(,2*nsg): momentum of end of substring nsg in pp-cms
c icjt(,2*nsg-1),icjt(,2*nsg): ic-code  of end of substrg nsg in pp-cms
c psg(,nsg): momentum of substring nsg in pp-cms
c xbkptl(,nptlc-nptl2): breakpoint of ns in pp-cms
c isptl(,nptlc-nptl2): 1 if string ns breaks
c-----------------------------------------------------------------------
      parameter (mxpc=500,mxsg=500)
      parameter (mxjt=2*mxsg)
      double precision xbkptl
      common/cpc/xbkptl(2,mxpc),isptl(mxpc)
      double precision psg,rotsg,pjt,xorsg
      common/csg/psg(5,mxsg),rotsg(3,mxsg),pjt(5,mxjt),xorsg(4,mxsg)
     *,isg(mxsg),icjt(2,mxjt),iorsg(mxsg)
      common/cjambr/nptl2,nptlc
 
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,taufos
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/cdelre/delrer
      parameter (nptf=129)
      common /cptf/xptf(nptf),wptf(nptf)
     *,qptfu(nptf),qptfs(nptf),qptfuu(nptf),qptfus(nptf),qptfss(nptf)
 
      parameter (nflav=6)
      integer jcp(nflav,2),jcm(nflav,2),jcpm(nflav,2),jcmp(nflav,2)
     *,jcps(nflav,2),jcms(nflav,2),icp(2),icm(2),icpm(2),icmp(2)
     *,icps(2),icms(2),ic(2),ics(2),jc(nflav,2),icx(2),jcx(nflav,2)
     *,icy(2),jcy(nflav,2),icj1(2),icj2(2),jcj1(nflav,2),jcj2(nflav,2)
     *,jcz(nflav,2),icz(2),icmz(2),icpz(2),jcpz(nflav,2),jcmz(nflav,2)
      double precision xor1(4),xor2(4),xbk(4),tor,zor
      double precision arp(4),arm(4),p0p,p0m,pp,pm,pt,px,py,phi,eeta
 
      ncnt3=0
      nsg0=nsg
9993  nsg=nsg0
 
      iret=0
      do i=1,2
      icp(i)=icjt(i,2*ns-1)
      icm(i)=icjt(i,2*ns)
      enddo
      call iddeco(icp,jcp)
      call iddeco(icm,jcm)
      np=0
      do 7 nf=1,nflav
      jc(nf,1)=jcp(nf,1)+jcm(nf,1)
      jc(nf,2)=jcp(nf,2)+jcm(nf,2)
7     np=np+jcp(nf,1)-jcp(nf,2)
      call idenco(jc,ic,ireten)
      id=idtra(ic,0,0,3)
      amms=utamnx(jcp,jcm)
      am=psg(5,ns)
 
c     split string
c     ------------
      p0p=0
      p0m=0
      do 5 j=2*ns-1,2*ns
      p0p=p0p+pjt(4,j)+pjt(3,j)
5     p0m=p0m+pjt(4,j)-pjt(3,j)
9994  ncnt3=ncnt3+1
           if(ncnt3.gt.100)then
           if(ish.ge.90)then
           call utmsg('jambr1')
           write(ifch,*)'*****  split kinem not possible.'
           write(ifch,112)
     *     (icjt(j,2*ns-1),j=1,2),(icjt(j,2*ns),j=1,2),psg(5,ns)
112        format(1x,2i8,4x,2i8,4x,f7.2)
           call utmsgf
           endif
           if(ns.eq.1)iret=9999
           if(ns.eq.1)return
           iret=9996
           return
           endif
      do 17 nf=1,nflav
      do 17 i=1,2
      jcpm(nf,i)=0
17    jcmp(nf,i)=0
 
c     determine flavour
c     -----------------
      nqu=0
      if(ish.ge.93)write(ifch,*)'ordinary string fragmentation'
      if(ish.ge.93)write(ifch,*)' '
      nqu=1
      if(rangen().le.pdiqua)nqu=2
      ii=0
      if(mod(np+nqu,3).eq.0)ii=1
      if(mod(np-nqu,3).eq.0)ii=2
           if(ii.eq.0)
     *     call utstop('jambr1: no singlet construction possible')
      ifltt=0
      do 8 n=1,nqu
      ifl=int(rangen()/pud)+1
      ifltt=ifltt*10+(ifl+1)/2
      jcpm(ifl,ii)=jcpm(ifl,ii)+1
8     jcmp(ifl,3-ii)=jcmp(ifl,3-ii)+1
      call idenco(jcpm,icpm,ireten)
           if(ireten.eq.1)
     *     call utstop('jambr1: idenco(jcpm... ret.code=1       ')
      call idenco(jcmp,icmp,ireten)
           if(ireten.eq.1)
     *     call utstop('jambr1: idenco(jcmp... ret.code=1       ')
      do 25 nf=1,nflav
      do 25 i=1,2
      jcps(nf,i)=jcp(nf,i)+jcpm(nf,i)
25    jcms(nf,i)=jcm(nf,i)+jcmp(nf,i)
      call idenco(jcps,icps,ireten)
      call idenco(jcms,icms,ireten)
      call idcomk(icps)
      call idcomk(icms)

c     calculate p+,p-,pt of string breaking
c     -------------------------------------
           if(ish.ge.93)write(ifch,109)
     *     icm(1),icmp(1),icpm(1),icp(1)
     *     ,icm(2),icmp(2),icpm(2),icp(2)
109        format(1x,'flavors:',2(i11,i7)/9x,2(i11,i7)/)
           if(ish.ge.93)write(ifch,*)'ifltt:',ifltt
           if(ish.ge.93)write(ifch,*)' '
      idp=idtra(icps,0,0,3)
      idm=idtra(icms,0,0,3)
      ammp=utamny(jcp,jcpm)
      ammm=utamny(jcm,jcmp)
      if(ifltt.eq.1)then
      pt=utinvt(nptf,xptf,qptfu ,rangen()*qptfu(nptf))
      elseif(ifltt.eq.2)then
      pt=utinvt(nptf,xptf,qptfs ,rangen()*qptfs(nptf))
      elseif(ifltt.eq.11)then
      pt=utinvt(nptf,xptf,qptfuu,rangen()*qptfuu(nptf))
      elseif(ifltt.eq.12.or.ifltt.eq.21)then
      pt=utinvt(nptf,xptf,qptfus,rangen()*qptfus(nptf))
      elseif(ifltt.eq.22)then
      pt=utinvt(nptf,xptf,qptfss,rangen()*qptfss(nptf))
      endif
      lo=1
      goto48
47    lo=lo+1
      pt=rangen()*pt
48    continue
      phi=2.*pi*rangen()
      px=pt*dcos(phi)
      py=pt*dsin(phi)
      tmmp=sqrt(pt**2+ammp**2)
      tmmm=sqrt(pt**2+ammm**2)
      area=-alog(rangen())/parea
      etam=dlog((tmmm**2+area)/(sqrt(area)*p0m))
      etap=dlog((sqrt(area)*p0p)/(tmmp**2+area))
      if(etam.gt.etap.and.lo.lt.5)goto47
      if(etam.gt.etap)goto9994
      if(nqu.eq.2)then
      endif
      eta=etam+rangen()*(etap-etam)
      amp=sqrt(p0p*sqrt(area)*exp(-eta)-area-pt**2)
      amm=sqrt(p0m*sqrt(area)*exp( eta)-area-pt**2)
      call idres(idp,amp,idpr,iadjp)
      call idres(idm,amm,idmr,iadjm)
      r=rangen()
           if(idpr.eq.110.and.r.lt.0.5)then
      idp=220
      amp=.549
      if(r.lt.0.6667)amp=.958
      call idres(idp,amp,idpr,iadjp)
      iadjp=1
           endif
      r=rangen()
           if(idmr.eq.110.and.r.lt.0.5)then
      idm=220
      amm=.549
      if(r.lt.0.6667)amm=.958
      call idres(idm,amm,idmr,iadjm)
      iadjm=1
           endif
      tmp=sqrt(pt**2+amp**2)
      tmm=sqrt(pt**2+amm**2)
           if(iadjp.eq.1.and.iadjm.ne.1)then
      eta=dlog((sqrt(area)*p0p)/(tmp**2+area))
      if(eta.lt.etam)goto9994
      amm=sqrt(p0m*sqrt(area)*exp(eta)-area-pt**2)
      call idres(idm,amm,idmr,iadjm)
      tmm=sqrt(pt**2+amm**2)
           endif
           if(iadjm.eq.1.and.iadjp.ne.1)then
      eta=dlog((tmm**2+area)/(sqrt(area)*p0m))
      if(eta.gt.etap)goto9994
      amp=sqrt(p0p*sqrt(area)*exp(-eta)-area-pt**2)
      call idres(idp,amp,idpr,iadjp)
      tmp=sqrt(pt**2+amp**2)
           endif
           if(iadjp.eq.1.and.iadjm.eq.1)then
      tm=(p0p*p0m-tmm**2-tmp**2)/2
      if(tm.lt.0.)goto9994
      if((tm-tmp*tmm)*(tm+tmp*tmm).lt.0.)goto9994
      area=tm-sqrt((tm-tmp*tmm)*(tm+tmp*tmm))
      eeta=p0p*sqrt(area)/(tmp**2+area)
      if(eeta.le.0.)goto9994
      eta=dlog(eeta)
           endif
      pp=sqrt(area)*exp(eta)
      pm=sqrt(area)*exp(-eta)
           if(p0p-pp-pt**2/pm.lt.0..or.p0m-pm-pt**2/pp.lt.0.)then
           if(ish.ge.90)then
           call utmsg('jambr1')
           write(ifch,*)'*****  negative jet energy. split again.'
           if(p0p-pp-pt**2/pm.lt.0.)
     *     write(ifch,*)'p0p: ',p0p,'   pp+pt**2/pm: ',pp+pt**2/pm
           if(p0m-pm-pt**2/pp.lt.0.)
     *     write(ifch,*)'p0m: ',p0m,'   pm+pt**2/pp: ',pm+pt**2/pp
           call utmsgf
           endif
           goto9994
           endif
 
c     determine break pnt and new origins
c     -----------------------------------
      xor1(1)=0
      xor1(2)=0
      xor1(3)=pp/2/tensn
      xor1(4)=pp/2/tensn
           if(ish.ge.93)write(ifch,118)xor1
118        format(' new origin +: ',13x,4f7.2)
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,xor1(1),xor1(2),xor1(3))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,xor1(1),xor1(2),xor1(3),xor1(4))
      do 44 j=1,4
44    xor1(j)=xor1(j)+xorsg(j,ns)
           if(ish.ge.93)write(ifch,117)xor1
      xbk(1)=0
      xbk(2)=0
      xbk(3)=(pp-pm)/2/tensn
      xbk(4)=(pp+pm)/2/tensn
           if(ish.ge.93)write(ifch,116)xbk
116        format(' breaking pnt: ',13x,4f7.2)
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,xbk(1),xbk(2),xbk(3))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,xbk(1),xbk(2),xbk(3),xbk(4))
      do 45 j=1,4
45    xbk(j)=xbk(j)+xorsg(j,ns)
      if(ish.ge.93)write(ifch,117)xbk
117   format(28x,4f7.2/)
      xor2(1)=0
      xor2(2)=0
      xor2(3)=-pm/2/tensn
      xor2(4)= pm/2/tensn
      if(ish.ge.93)write(ifch,119)xor2
119   format(' new origin -: ',13x,4f7.2)
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,xor2(1),xor2(2),xor2(3))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,xor2(1),xor2(2),xor2(3),xor2(4))
      do 46 j=1,4
46    xor2(j)=xor2(j)+xorsg(j,ns)
           if(ish.ge.93)write(ifch,117)xor2
 
c     store the two substrings
c     ------------------------
      nsgb=nsg
      do 9998 m=1,2
      nsg=nsg+1
           if(nsg.gt.mxsg)
     *     call utstop('jambr1: nsg>mxsg                        ')
      isg(nsg)=nptlc
      if(m.eq.1)then !-------------------
      do 41 j=1,4
41    xorsg(j,nsg)=xor1(j)
      pjt(1,2*nsg-1)=0
      pjt(2,2*nsg-1)=0
      pjt(3,2*nsg-1)=(p0p-pp-pt**2/pm)/2
      pjt(4,2*nsg-1)=(p0p-pp-pt**2/pm)/2
      pjt(5,2*nsg-1)=0
      do 10 j=1,2
10    icjt(j,2*nsg-1)=icjt(j,2*ns-1)
      pjt(1,2*nsg)=px
      pjt(2,2*nsg)=py
      pjt(3,2*nsg)=(pt**2/pm-pm)/2
      pjt(4,2*nsg)=(pt**2/pm+pm)/2
      pjt(5,2*nsg)=0
      do 11 j=1,2
      icjt(j,2*nsg)=icpm(j)
      do 11 nf=1,nflav
      jc(nf,j)=jcps(nf,j)
11    continue
      else !-------------------------
      do 42 j=1,4
42    xorsg(j,nsg)=xor2(j)
      pjt(1,2*nsg-1)=-px
      pjt(2,2*nsg-1)=-py
      pjt(3,2*nsg-1)=(pp-pt**2/pp)/2
      pjt(4,2*nsg-1)=(pp+pt**2/pp)/2
      pjt(5,2*nsg-1)=0
      do 12 j=1,2
      icjt(j,2*nsg-1)=icmp(j)
      do 12 nf=1,nflav
      jc(nf,j)=jcms(nf,j)
12    continue
      pjt(1,2*nsg)=0
      pjt(2,2*nsg)=0
      pjt(3,2*nsg)=-(p0m-pm-pt**2/pp)/2
      pjt(4,2*nsg)= (p0m-pm-pt**2/pp)/2
      pjt(5,2*nsg)=0
      do 13 j=1,2
13    icjt(j,2*nsg)=icjt(j,2*ns)
      endif !---------------------
      do 18 i=1,2
      icx(i)=icjt(i,2*nsg-1)
18    icy(i)=icjt(i,2*nsg)
      call iddeco(icx,jcx)
      call iddeco(icy,jcy)
      do 28 n=1,nflav
      do 28 i=1,2
      k=jcx(n,i)+jcy(n,i)-jc(n,i)
      if(k.le.0)goto28
      do 19 l=1,k
      jx=jcx(n,i)
      jy=jcy(n,i)
      if(jx.gt.0.and.jy.eq.0)then
      jcx(n,i)=jcx(n,i)-1
      elseif(jx.eq.0.and.jy.gt.0)then
      jcy(n,i)=jcy(n,i)-1
      elseif(jx.gt.0.and.jy.gt.0)then
      if(rangen().lt..5)then
      jcx(n,i)=jcx(n,i)-1
      else
      jcy(n,i)=jcy(n,i)-1
      endif
      else
           write(ifch,*)('*',lp=1,71)
           write(ifch,*)'*****  ic;  ',ic
           write(ifch,*)'*****  icx: ',icx
           write(ifch,*)'*****  icy: ',icy
           write(ifch,*)'*****  n,i,k,l: ',n,i,k,l
           write(ifch,*)'*****  jx,jy: ',jx,jy
           write(ifch,*)('*',lp=1,71)
           call utstop('jambr1: error during jet compactificatn ')
      endif
19    continue
28    continue
      call idenco(jcx,icx,ireten)
           if(ireten.eq.1)
     *     call utstop('jambr1: idenco(jcx... ret.code=1        ')
      call idenco(jcy,icy,ireten)
           if(ireten.eq.1)
     *     call utstop('jambr1: idenco(jcy... ret.code=1        ')
      do 20 i=1,2
      icjt(i,2*nsg-1)=icx(i)
20    icjt(i,2*nsg)=icy(i)
           if(ish.ge.93)write(ifch,108)nsg,ns
     *     ,(icjt(j,2*nsg-1),j=1,2),(sngl(pjt(j,2*nsg-1)),j=1,5)
108        format(2x,i3,i4,2x,2i7,3x,5(e8.2,1x))
           if(ish.ge.93)write(ifch,108)nsg,ns
     *     ,(icjt(j,2*nsg  ),j=1,2),(sngl(pjt(j,2*nsg  )),j=1,5)
      do 14 j=1,4
14    psg(j,nsg)=pjt(j,2*nsg-1)+pjt(j,2*nsg)
      psg(5,nsg)=
     *sqrt(psg(4,nsg)**2-psg(3,nsg)**2-psg(2,nsg)**2-psg(1,nsg)**2)
           if(ish.ge.93)write(ifch,110)nsg,ns,(psg(j,nsg),j=1,5)
110        format(2x,i3,i4,19x,5(e8.2,1x))
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,psg(1,nsg),psg(2,nsg),psg(3,nsg))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,psg(1,nsg),psg(2,nsg),psg(3,nsg),psg(4,nsg))
      iorsg(nsg)=ns
           if(ish.ge.93)
     *     write(ifch,102)nsg,iorsg(nsg),(psg(j,nsg),j=1,5)
102        format(2x,i3,i4,19x,5(e8.2,1x))
           if(ish.ge.93)write(ifch,*)' '
c-c   do 56 i=1,2
c-c   icpz(i)=icjt(i,2*nsg-1)
c-c56 icmz(i)=icjt(i,2*nsg)
c-c   call iddeco(icpz,jcpz)
c-c   call iddeco(icmz,jcmz)
c-c   do 57 nf=1,nflav
c-c   jcz(nf,1)=jcpz(nf,1)+jcmz(nf,1)
c-c57 jcz(nf,2)=jcpz(nf,2)+jcmz(nf,2)
c-c   call idenco(jcz,icz,ireten)
c-c   idz=idtra(icz,0,0,3)
c-c   amz=psg(5,nsg)
c-c   call idres(idz,amz,idrz,iadj)
c-c   if(idrz.eq.110.and.rangen().lt.0.5)goto9993
9998  continue
 
      xbkptl(1,nptlc-nptl2)=xbk(3)
      xbkptl(2,nptlc-nptl2)=xbk(4)
      isptl(nptlc-nptl2)=1
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jambr2(ns,nsg,iret)
c-----------------------------------------------------------------------
c "SAMBA" (scholten+werner, mar 92)
c breaks string ns according to a-m model.
c ns: points to the current fragmenting string.
c nsg: points to the last produced substring (succ increased)
c     input:
c psg(,ns): momentum of string ns in pp-cms
c pjt(,2*ns-1),pjt(,2*ns): momentum of end of string ns in pp-cms
c     output:
c xorsg(,nsg): origin of substring nsg in pp-cms
c pjt(,2*nsg-1),pjt(,2*nsg): momentum of end of substring nsg in pp-cms
c icjt(,2*nsg-1),icjt(,2*nsg): ic-code  of end of substrg nsg in pp-cms
c psg(,nsg): momentum of substring nsg in pp-cms
c xbkptl(,nptlc-nptl2): breakpoint of ns in pp-cms
c isptl(,nptlc-nptl2): 1 if string ns breaks
c-----------------------------------------------------------------------
      parameter (mxpc=500,mxsg=500)
      parameter (mxjt=2*mxsg)
      double precision xbkptl
      common/cpc/xbkptl(2,mxpc),isptl(mxpc)
      double precision psg,rotsg,pjt,xorsg
      common/csg/psg(5,mxsg),rotsg(3,mxsg),pjt(5,mxjt),xorsg(4,mxsg)
     *,isg(mxsg),icjt(2,mxjt),iorsg(mxsg)
      common/cjambr/nptl2,nptlc
 
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,taufos
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/cdelre/delrer
      parameter (nptf=129)
      common /cptf/xptf(nptf),wptf(nptf)
     *,qptfu(nptf),qptfs(nptf),qptfuu(nptf),qptfus(nptf),qptfss(nptf)
 
      parameter (nflav=6)
      integer jcp(nflav,2),jcm(nflav,2),jcpm(nflav,2),jcmp(nflav,2)
     *,jcps(nflav,2),jcms(nflav,2),icp(2),icm(2),icpm(2),icmp(2)
     *,icps(2),icms(2),ic(2),ics(2),jc(nflav,2),icx(2),jcx(nflav,2)
     *,icy(2),jcy(nflav,2),icj1(2),icj2(2),jcj1(nflav,2),jcj2(nflav,2)
     *,jcz(nflav,2),icz(2),icmz(2),icpz(2),jcpz(nflav,2),jcmz(nflav,2)
      double precision xor1(4),xor2(4),xbk(4),tor,zor,pw(4)
      double precision arp(4),arm(4),p0p,p0m,pp,pm,pt,px,py,phi,eeta
 
      DATA DONE,DTWO,DFOUR,fyst/1.D0,2.D0,4.D0,0.0001/
 
      DIMENSION XBREAK(100),YBREAK(100),PXBREAK(100),PYBREAK(100),
     & ITRD(100),ICPMBR(2,100),ICMPBR(2,100)
      PARAMETER (NBRM=99)
      REAL*8 DY,DR,DALFA,DYST,DONE,AMAX,RA,ANEW,DY1,DY2,DYT,AMAM,
     & RMX,RMY,A,XP,YP,XBR,YBR,YTD,DTWO,DFOUR,DA,YNEW,
     & DATWID,X2
      LOGICAL LAST
 
      ncnt3=0
      nsg0=nsg
9993  nsg=nsg0
 
c initialization
c --------------
 
      iret=0
      do 6 i=1,2
      icp(i)=icjt(i,2*ns-1)
      icm(i)=icjt(i,2*ns)
      ICPMBR(i,1)=ICM(I)	!flavor of left side,=break#1
    6 CONTINUE
 
      call iddeco(icp,jcp)
      call iddeco(icm,jcm)
      np=0
      do 7 nf=1,nflav
      np=np+jcp(nf,1)-jcp(nf,2)
7     continue
      AMMS=UTAMNX(JCP,JCM)
      am=psg(5,ns)
C
      p0p=0
      p0m=0
      do 5 j=2*ns-1,2*ns
      p0p=p0p+pjt(4,j)+pjt(3,j)
5     p0m=p0m+pjt(4,j)-pjt(3,j)
      XT=P0M
      YT=P0P
C
      DALFA=DBLE(PAREA)
C	initialize
      XBREAK(1)=XT
      YBREAK(1)=0.
      PXBREAK(1)=0
      PYBREAK(1)=0
 
c  redo
c  ----
 
      ncnt3=0
9994  ncnt3=ncnt3+1
      IBR=1
      XP=XT
      YP=0.
      YTD=YT
           if(ncnt3.gt.100)then
           if(ish.ge.90)then
           call utmsg('jambr2')
           write(ifch,*)'*****  split kinem not possible.'
           write(ifch,112)
     *     (icjt(j,2*ns-1),j=1,2),(icjt(j,2*ns),j=1,2),psg(5,ns)
112        format(1x,2i8,4x,2i8,4x,f7.2)
           call utmsgf
           endif
           if(ns.eq.1)iret=9999
           if(ns.eq.1)return
           iret=9996
           return
           endif
 
c  search for breakpoints
c  ----------------------
 
    9 CONTINUE
      do 17 nf=1,nflav
      do 17 i=1,2
      jcpm(nf,i)=0
17    jcmp(nf,i)=0
 
c     ordinary strings
c     ----------------
      nqu=0
      if(ish.ge.92)write(ifch,*)'ordinary string fragmentation'
      if(ish.ge.92)write(ifch,*)' '
      nqu=1
      if(rangen().le.pdiqua)nqu=2
      ii=0
      if(mod(np+nqu,3).eq.0)ii=1
      if(mod(np-nqu,3).eq.0)ii=2
           if(ii.eq.0)
     *     call utstop('jambr2: no singlet construction possible')
      ifltt=0
      do 18 n=1,nqu
      ifl=int(rangen()/pud)+1
      ifltt=ifltt*10+(ifl+1)/2
      jcpm(ifl,ii)=jcpm(ifl,ii)+1
18    jcmp(ifl,3-ii)=jcmp(ifl,3-ii)+1
      call idenco(jcpm,icpm,ireten)
           if(ireten.eq.1)
     *     call utstop('jambr2: idenco(jcpm... ret.code=1       ')
      call idenco(jcmp,icmp,ireten)
           if(ireten.eq.1)
     *     call utstop('jambr2: idenco(jcmp... ret.code=1       ')
 
c  calculate pt of string breaking
c  -------------------------------------
           if(ish.ge.92)write(ifch,109)
     *     icm(1),icmp(1),icpm(1),icp(1)
     *     ,icm(2),icmp(2),icpm(2),icp(2)
109        format(1x,'flavors:',2(i11,i7)/9x,2(i11,i7)/)
           if(ish.ge.92)write(ifch,*)'ifltt:',ifltt
           if(ish.ge.92)write(ifch,*)' '
      ammp=UTAMNY(JCP,JCPM)
      ammm=UTAMNY(JCMP,JCM)
      if(ifltt.eq.1)then
      pt=utinvt(nptf,xptf,qptfu ,rangen()*qptfu(nptf))
      elseif(ifltt.eq.2)then
      pt=utinvt(nptf,xptf,qptfs ,rangen()*qptfs(nptf))
      elseif(ifltt.eq.11)then
      pt=utinvt(nptf,xptf,qptfuu,rangen()*qptfuu(nptf))
      elseif(ifltt.eq.12.or.ifltt.eq.21)then
      pt=utinvt(nptf,xptf,qptfus,rangen()*qptfus(nptf))
      elseif(ifltt.eq.22)then
      pt=utinvt(nptf,xptf,qptfss,rangen()*qptfss(nptf))
      endif
 
      RMX=0.001
      IF(IBR.EQ.1) RMX=AMMM**2+PT**2
      RMY=AMMP**2+PT**2
      IF(RMY.GT.(XP*(YTD-YP))) GOTO 8
 
C     Single Artru-Mennessier Break (SAMB)
c     ------------------------------------
C     in: XP,YP,YTD,RMX,RMY
C     in: FYST,DALFA
C     out: XBR,YBR
 
      DYT=YTD-YP
C	Y-stepping range
      AMAM=DONE+(RMX-RMY)/(XP*DYT)
      A=(DONE-DFOUR*RMX/(AMAM*AMAM*XP*DYT))
      IF(A .LT. 0.) GOTO 8	! string cannot break any further
      A=SQRT(A)
      DY1=DYT*AMAM*(DONE-A)/DTWO
      DY2=DYT*AMAM*(DONE+A)/DTWO
      DYST=DBLE(FYST)*(DYT-DY2)
    3 CONTINUE
      DR=DBLE(1.-RANgen())
      AMAX=XP*DYT*A*AMAM+RMX*LOG(DY1/DY2)+RMY*LOG((DYT-DY2)/(DYT-DY1))
      RA=-DLOG(DR)/DALFA
      IF(RA.GT.AMAX .AND. IBR.GT.1) GOTO 8
      RA=MOD(RA,AMAX)
      DY=DYT*DSQRT(RMX)/(DSQRT(RMX)+DSQRT(RMY))	! second derivative area
      DATWID=XP*(DY-DY1)+RMX*LOG(DY1/DY)+RMY*LOG((DYT-DY)/(DYT-DY1))
      ANEW=RA-DATWID	! DATWID= area for which second derivative=0
      ICOUNT=0
      IF(ANEW .GT. 0.) GOTO 1
    2 CONTINUE
      DA=XP - RMX/DY - RMY/(DYT-DY)	! first derivative of area
      YNEW=DY+ANEW/DA	! always approach from above
      DY=INT(YNEW/DYST)*DYST
      IF(DY.LT.DY1) THEN
      	WRITE(IFCH,*) 'DY,DY1',DY,DY1
      	GOTO 4
      ENDIF
      ANEW=RA-XP*(DY-DY1)-RMX*LOG(DY1/DY)-RMY*LOG((DYT-DY)/(DYT-DY1))	!
      IF(ANEW .LT. 0.) GOTO 2
      GOTO 4
    1 CONTINUE
      DA=XP - RMX/DY - RMY/(DYT-DY)
      YNEW=DY+ANEW/DA	! always approach from BELOW
      DY=(INT(YNEW/DYST)+1)*DYST
      ANEW=RA-XP*(DY-DY1)-RMX*LOG(DY1/DY)-RMY*LOG((DYT-DY)/(DYT-DY1))
      IF(ANEW .GT. 0.) GOTO 1
      DY=DY-DYST
      ANEW=RA-XP*(DY-DY1)-RMX*LOG(DY1/DY)-RMY*LOG((DYT-DY)/(DYT-DY1))
    4 CONTINUE
      YBR=DY+RANgen()*DYST
      IF(YBR .GT. DY2) YBR=DY2
      XBR=ANEW/DYST+RMY/(DYT-YBR)
      X2=(XP-RMX/YBR)
      IF(XBR .GT. X2) XBR=X2
C	BETTER:	solve for YBR from ANEW=0
C		find XBR from homogeneous (not exp) distr, x1<x<x2
C      X1=RMY/(DYT-YBR)
      YBR=YBR+YP
C
C	end SAMB
C
      IF(IBR.GE.NBRM)
     *     call utstop('jambr2: IBR>NBRM                        ')
      IBR=IBR+1
      XBREAK(IBR)=XBR
      YBREAK(IBR)=YBR
      phi=2.*pi*rangen()
      pxBREAK(IBR)=pt*dcos(phi)
      pyBREAK(IBR)=pt*dsin(phi)
      DO I=1,2
      ICPMBR(i,IBR)=ICPM(I)
      ICMPBR(i,IBR)=ICMP(I)
      ENDDO
      XP=XBR	! better, but longer: =XT
      YP=YBR
      amleft=sqrt(ybr*(xt-xbr))
	amright=sqrt(xbr*(yt-ybr))
      GOTO 9	! search for next break point
    8 CONTINUE  ! ibr loop finished
 
      IF(IBR.EQ.1) GOTO 9994	! no breakpoint found
 
C  initialize tail end
C  -------------------
      XBREAK(IBR+1)=0.
      YBREAK(IBR+1)=YT
      PXBREAK(IBR+1)=0
      PYBREAK(IBR+1)=0
      DO I=1,2
      ICmpBR(i,IBR+1)=ICP(I)	!flavor of right side,=break#last
      ENDDO
 
c          print
c          -----
           if(ish.ge.92)then
           write(ifch,*)' '
           write(ifch,*)'icmpbr(1/2,)  x/ybreak  px/ybreak:'
           do 52 ib=1,ibr+1
           write(ifch,100)icmpbr(1,ib),icmpbr(2,ib)
     *     ,xbreak(ib),ybreak(ib),pxbreak(ib),pybreak(ib)
52         write(ifch,100)icpmbr(1,ib),icpmbr(2,ib)
100        format(1x,i10,i10,2e11.3,7x,2e11.3)
           write(ifch,*)' '
           endif
 
C  Time order breakpoints
C  ----------------------
 
      T1=0.
      DO J=2,IBR
      T2=2.
      DO 20 I=2,IBR
C      T=XBREAK(I)/XT+YBREAK(I)/YT	! time
      T=XBREAK(I)*YBREAK(I)/(XT*YT)	! area past
      IF(T.LE.T1 .OR. T.GT.T2) GOTO 20
      T2=T
      NT=I
   20 CONTINUE
      T1=T2
      ITRD(J)=NT
      ENDDO	! LOOP J
 
c          print
c          -----
           if(ish.ge.92)then
           write(ifch,*)' '
           write(ifch,*)'itrd:'
           write(ifch,*)(itrd(ib),ib=2,ibr)
           write(ifch,*)' '
           endif
 
C  Check acceptance criterium
c  --------------------------
 
      nbr=0
      DO 12 J=2,IBR	! check points according to T ordering
      I=ITRD(J)
      XBRI=XBREAK(I)
      YBRI=YBREAK(I)
C	Find neighboring earlier break points
      IR=IBR+1
      IL=1
      DO JN=2,J
      IN=ITRD(JN)
      IF(IN.GT.I .AND. IN.LT.IR) IR=IN
      IF(IN.LT.I .AND. IN.GT.IL) IL=IN
      ENDDO	! loop JN
      PML=(XBREAK(IL)-XBRI)
      PPL=(YBRI-YBREAK(IL))
      PMR=(XBRI-XBREAK(IR))
      PPR=(YBREAK(IR)-YBRI)
      PX=PXBREAK(I)
      PY=PYBREAK(I)
      PXL=PX-PXBREAK(IL)
      PYL=PY-PYBREAK(IL)
      ammm=PML*PPL-(PXL*PXL+PYL*PYL)
      PXR=PXBREAK(IR)-PX
      PYR=PYBREAK(IR)-PY
      ammp=PMR*PPR-(PXR*PXR+PYR*PYR)
      if(ammp.le.0.) goto 13
      if(ammm.le.0.) goto 13
      ammp=sqrt(ammp)
      ammm=sqrt(ammm)
      amp=ammp
      amm=ammm
      call utresm(icmpbr(1,ir),icmpbr(2,ir)
     *,icpmbr(1,i),icpmbr(2,i),amp,idpr,iadjp,ireten)
      call utresm(icpmbr(1,il),icpmbr(2,il)
     *,icmpbr(1,i),icmpbr(2,i),amm,idmr,iadjm,ireten)
      IF(amp.gt.ammp) GOTO 13
      IF(amm.gt.ammm) GOTO 13
      pt2=px*px+py*py
      ptl2=pxbreak(il)**2+pybreak(il)**2
      ptr2=pxbreak(ir)**2+pybreak(ir)**2
      d=(pml*ppl-ptl2-pt2)/2.
      d=d*d-pt2*ptl2
      if(d.le.0.)goto13
      d=(pmr*ppr-ptr2-pt2)/2.
      d=d*d-pt2*ptr2
      if(d.le.0.)goto13
      nbr=nbr+1
      GOTO 12
   13 CONTINUE
      ITRD(J)=-1
      XBREAK(I)=-1
      YBREAK(I)=-1
   12 CONTINUE	! loop J
      if(nbr.eq.0) goto 9994
 
c          print
c          -----
           if(ish.ge.92)then
           write(ifch,*)' '
           write(ifch,*)'icmpbr(1/2,)  x/ybreak  px/ybreak:'
           do 53 ib=1,ibr+1
           write(ifch,100)icmpbr(1,ib),icmpbr(2,ib)
     *     ,xbreak(ib),ybreak(ib),pxbreak(ib),pybreak(ib)
53         write(ifch,100)icpmbr(1,ib),icpmbr(2,ib)
           write(ifch,*)' '
           endif
 
c  Build new strings
c  -----------------
 
      NBREAK=0
      TAUAVE=0.
      IL=1
      IN=IL
      IR=1
      LAST=.FALSE.
   11 CONTINUE
      IR=IR+1
      IF(IR.EQ.(IBR+1)) LAST=.TRUE.
      IF(XBREAK(IR).LT.0.) GOTO 11	! breakpoint was not acceptable
      IF(IN.EQ.1) GOTO 14
   15 CONTINUE
           if(ish.ge.92)then
           write(ifch,*)' '
           write(ifch,*)'il in ir:  ',il,in,ir
           write(ifch,*)' '
           endif
      PML=XBREAK(IL)-XBREAK(IN)
      PPL=YBREAK(IN)-YBREAK(IL)
      PMR=(XBREAK(IN)-XBREAK(IR))
      PPR=(YBREAK(IR)-YBREAK(IN))
      PX=PXBREAK(IN)
      PY=PYBREAK(IN)
      PXL=PX-PXBREAK(IL)
      PYL=PY-PYBREAK(IL)
      PTL2=(PXL*PXL+PYL*PYL)
      AMM=SQRT(PML*PPL-PTL2)
      PXR=PXBREAK(IR)-PX
      PYR=PYBREAK(IR)-PY
      PTR2=(PXR*PXR+PYR*PYR)
      AMP=SQRT(PMR*PPR-PTR2)
      call utresm(icmpbr(1,ir),icmpbr(2,ir)
     *,icpmbr(1,in),icpmbr(2,in),amp,idpr,iadjp,ireten)
      call utresm(icpmbr(1,il),icpmbr(2,il)
     *,icmpbr(1,in),icmpbr(2,in),amm,idmr,iadjm,ireten)
      AREA=PPL*PMR
      P0P=PPL+PPR
      P0M=PML+PMR
      tmm=sqrt(ptL2+amm**2)
      tmP=sqrt(ptR2+amP**2)
           if(tmp*tmp.gt.ppr*pmr+1e-4)then
           call utmsg('jambr2')
           write(ifch,*)'*****  tmp*tmp.gt.ppr*pmr'
           write(ifch,*)'tmp*tmp ppr*pmr: ',tmp*tmp,ppr*pmr
           call utmsgf
           endif
	   if(tmm*tmm.gt.ppl*pml+1e-4)then
           call utmsg('jambr2')
           write(ifch,*)'*****  tmm*tmm.gt.ppl*pml'
           write(ifch,*)'tmm*tmm ppl*pml: ',tmm*tmm,ppl*pml
           call utmsgf
           endif
      ETA=0.5*ALOG(PPL/PMR)
c ------ adjust p not m
      if(iadjp.eq.1.and.iadjm.ne.1)then
           if(ish.ge.92)write(ifch,*)'adjust p not m'
      etaNEW=dlog((sqrt(area)*p0p)/(tmp**2+area))
      if(etaNEW.lt.eta-1e-4) THEN	!should never happen
           call utmsg('jambr2')
           write(ifch,*)'*****  etanew.lt.eta'
           write(ifch,*)'etanew eta: ',etanew,eta
           call utmsgf
      ETANEW=ETA
      ENDIF
      xbreak(IN)=SQRT(AREA)*EXP(-ETANEW)+XBREAK(IR)
      YBREAK(IN)=SQRT(AREA)*EXP(ETANEW)+YBREAK(IL)
c ------ adjust m not p
      ELSEif(iadjM.eq.1.and.iadjP.ne.1)then
           if(ish.ge.92)write(ifch,*)'adjust m not p'
      etaNEW=-dlog((sqrt(area)*p0M)/(tmM**2+area))
      if(etaNEW.Gt.eta+1e-4) THEN	!should never happen
           call utmsg('jambr2')
           write(ifch,*)'*****  etanew.gt.eta'
           write(ifch,*)'etanew eta: ',etanew,eta
           call utmsgf
      ETANEW=ETA
      ENDIF
      xbreak(IN)=SQRT(AREA)*EXP(-ETANEW)+XBREAK(IR)
      YBREAK(IN)=SQRT(AREA)*EXP(ETANEW)+YBREAK(IL)
c ------ adjust both
      ELSEif(iadjp.eq.1.and.iadjm.eq.1)then
           if(ish.ge.92)write(ifch,*)'adjust both'
      tmp2=tmp**2
      tmm2=tmm**2
      d=(p0p*p0m-tmp2-tmm2)**2-4*tmp2*tmm2
           if(d.lt.0.)then
           call utmsg('jambr2')
           write(ifch,*)'***** negative d'
           write(ifch,*)'d: ',d
           call utmsgf
           d=0.
           endif
      d=sqrt(d)
      t=p0p*p0m-tmp2+tmm2
      p1=(t+d)/(2*p0m)
      p2=(t-d)/(2*p0m)
      if((p1+p2)/2.gt.ppl)then
      p=p2
      else
      p=p1
      endif
      xbreak(in)=p0m-tmm2/p+xbreak(ir)
      ybreak(in)=p+ybreak(il)
      endif
c ------
           if(ish.ge.92)then
           write(ifch,*)' '
           write(ifch,*)'x/ybreak:',xbreak(in),ybreak(in)
           write(ifch,*)' '
           endif
 
c     write sg
c     --------
      NBREAK=NBREAK+1
      TAUAVE=TAUAVE+XBREAK(IN)*YBREAK(IN)
      MM=1
      IF(LAST) MM=2
      do 9998 m=1,MM
      nsg=nsg+1
           if(nsg.gt.mxsg)
     *     call utstop('jambr2: nsg>mxsg                        ')
           if(ish.ge.92)write(ifch,*)'nsg last:  ',nsg,last
      isg(nsg)=nptlc
      IF(M.EQ.2) THEN
      IL=IN
      IN=IR
      ENDIF
      xor2(1)=0
      xor2(2)=0
      xor2(3)=(YBREAK(IL)-XBREAK(IN))/2/tensn
      xor2(4)=(YBREAK(IL)+XBREAK(IN))/2/tensn
      if(ish.ge.92)write(ifch,119)xor2
119   format(/' origin: ',13x,4f7.2/)
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,xor2(1),xor2(2),xor2(3))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,xor2(1),xor2(2),xor2(3),xor2(4))
      do 46 j=1,4
46    xor2(j)=xor2(j)+xorsg(j,ns)
c          if(ish.ge.92)write(ifch,117)xor2
      do 42 j=1,4
42    xorsg(j,nsg)=xor2(j)
      PP=YBREAK(IN)-YBREAK(IL)
      PM=XBREAK(IL)-XBREAK(IN)
      psg(1,nsg)=pxbreak(in)-pxbreak(il)
      psg(2,nsg)=pybreak(in)-pybreak(il)
      psg(3,nsg)=(pp-pm)/2
      psg(4,nsg)=(pp+pm)/2
      ss=
     *psg(4,nsg)**2-psg(3,nsg)**2-psg(2,nsg)**2-psg(1,nsg)**2
      psg(5,nsg)=sqrt(ss)
           if(ish.ge.92)write(ifch,110)nsg,ns,(psg(j,nsg),j=1,5)
110        format(2x,i3,i4,19x,5(e8.2,1x))
      do 55 k=1,4
55    pw(k)=psg(k,nsg)
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,psg(1,nsg),psg(2,nsg),psg(3,nsg))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,psg(1,nsg),psg(2,nsg),psg(3,nsg),psg(4,nsg))
      iorsg(nsg)=ns
           if(ish.ge.92)
     *     write(ifch,102)nsg,iorsg(nsg),(psg(j,nsg),j=1,5)
102        format(2x,i3,i4,19x,5(e8.2,1x))
           if(ish.ge.92)write(ifch,*)' '
 
c     write jt
c     --------
      do 41 j=1,2
      icjt(j,2*nsg-1)=icmpBR(j,IN)
      do 41 nf=1,nflav
      jc(nf,j)=jcms(nf,j)
41    continue
      do 43 j=1,2
43    icjt(j,2*nsg)=icPMBR(j,IL)
      idx=iadjm
      if(m.eq.2)idx=iadjp
           if(idx.eq.0)then
      pjt(1,2*nsg-1)=pxBREAK(IN)            ! right jet
      pjt(2,2*nsg-1)=pyBREAK(IN)
      PTJR2=PXBREAK(IN)**2+PYBREAK(IN)**2
      PTJL2=PXBREAK(IL)**2+PYBREAK(IL)**2
      AA=PM*PP-PTJL2-PTJR2
      PRM=(AA/2+PTJR2-SQRT(AA*AA/4-PTJL2*PTJR2))/PP
      PLP=(AA/2+PTJL2-SQRT(AA*AA/4-PTJL2*PTJR2))/PM
      pjt(3,2*nsg-1)=(PP-PLP-PRM)/2
      pjt(4,2*nsg-1)=(pp-PLP+PRM)/2
      pjt(5,2*nsg-1)=0
      pjt(1,2*nsg)=-PXBREAK(IL)             ! left jet
      pjt(2,2*nsg)=-PYBREAK(IL)
      pjt(3,2*nsg)=(PLP-PM+PRM)/2
      pjt(4,2*nsg)= (PLP+PM-PRM)/2
      pjt(5,2*nsg)=0
           if(ish.ge.92)then
      write(ifch,108)nsg,ns
     *,(icjt(j,2*nsg-1),j=1,2),(sngl(pjt(j,2*nsg-1)),j=1,5)
108   format(2x,i3,i4,2x,2i7,3x,5(e8.2,1x))
      write(ifch,108)nsg,ns
     *,(icjt(j,2*nsg  ),j=1,2),(sngl(pjt(j,2*nsg  )),j=1,5)
           endif
      err=0
      do 44 j=1,4
44    err=err+(pw(j)-pjt(j,2*nsg-1)-pjt(j,2*nsg))**2
           endif
      do 56 i=1,2
      icpz(i)=icjt(i,2*nsg-1)
56    icmz(i)=icjt(i,2*nsg)
      call iddeco(icpz,jcpz)
      call iddeco(icmz,jcmz)
      do 57 nf=1,nflav
      jcz(nf,1)=jcpz(nf,1)+jcmz(nf,1)
57    jcz(nf,2)=jcpz(nf,2)+jcmz(nf,2)
      call idenco(jcz,icz,ireten)
      idz=idtra(icz,0,0,3)
      amz=psg(5,nsg)
      call idres(idz,amz,idrz,iadj)
      if(idrz.eq.110.and.rangen().lt.0.5)goto9993
9998  CONTINUE
 
   14 CONTINUE
      IF(LAST) GOTO 10
      IL=IN
      IN=IR
      GOTO 11
   10 CONTINUE
c
c     determine break pnt
c     -----------------------------------
      TAUAVE=TAUAVE/NBREAK
      xbk(1)=0
      xbk(2)=0
      xbk(3)=(SQRT(TAUAVE*YT/XT)-SQRT(TAUAVE*XT/YT))/2/tensn
      xbk(4)=(SQRT(TAUAVE*YT/XT)+SQRT(TAUAVE*XT/YT))/2/tensn
           if(ish.ge.92)write(ifch,116)xbk
116        format(' breaking pnt: ',13x,4f7.2)
      call utrot2(-1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,xbk(1),xbk(2),xbk(3))
      call utlob2(-1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,xbk(1),xbk(2),xbk(3),xbk(4))
      do 45 j=1,4
45    xbk(j)=xbk(j)+xorsg(j,ns)
      if(ish.ge.92)write(ifch,117)xbk
117   format(28x,4f7.2/)
C
      xbkptl(1,nptlc-nptl2)=xbk(3)
      xbkptl(2,nptlc-nptl2)=xbk(4)
      isptl(nptlc-nptl2)=1
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jamfra(js,newevt)
c-----------------------------------------------------------------------
c     fragments string js according to a-m model.
c     version mar-92 (calls jambr1 or jambr2)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,taufos
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cdelre/delrer
 
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
 
      parameter (mxpc=500,mxsg=500)
      parameter (mxjt=2*mxsg)
      double precision xbkptl
      common/cpc/xbkptl(2,mxpc),isptl(mxpc)
      double precision psg,rotsg,pjt,xorsg
      common/csg/psg(5,mxsg),rotsg(3,mxsg),pjt(5,mxjt),xorsg(4,mxsg)
     *,isg(mxsg),icjt(2,mxjt),iorsg(mxsg)
      common/cjambr/nptl2,nptlc
 
      parameter (nptf=129)
      common /cptf/xptf(nptf),wptf(nptf)
     *,qptfu(nptf),qptfs(nptf),qptfuu(nptf),qptfus(nptf),qptfss(nptf)
      parameter (ngau=129)
      common /cgau/xgau(ngau),qgau(ngau),wgau(ngau)
 
      parameter (nflav=6)
      integer jcp(nflav,2),jcm(nflav,2),jcpm(nflav,2),jcmp(nflav,2)
     *,jcps(nflav,2),jcms(nflav,2)
     *,icp(2),icm(2),icpm(2),icmp(2),icup(2),icum(2)
     *,icps(2),icms(2),ic(2),ics(2),jc(nflav,2),icx(2),jcx(nflav,2)
     *,icy(2),jcy(nflav,2)
     *,icj1(2),icj2(2),jcj1(nflav,2),jcj2(nflav,2)
      double precision xor1(4),xor2(4),xbk(4),tor,zor
      double precision arp(4),arm(4),p0p,p0m,pp,pm,pt,px,py,phi,eeta
 
      delrec=0.600
      delrex=0.050
      zero=0
      newevt=0
 
      psg(1,1)=pstr(1,js)
      psg(2,1)=pstr(2,js)
      psg(3,1)=pstr(3,js)
      psg(4,1)=pstr(4,js)
      psg(5,1)=pstr(5,js)
      rotsg(1,1)=rotstr(1,js)
      rotsg(2,1)=rotstr(2,js)
      rotsg(3,1)=rotstr(3,js)
      xorsg(1,1)=xorstr(1,js)
      xorsg(2,1)=xorstr(2,js)
      xorsg(3,1)=xorstr(3,js)
      xorsg(4,1)=xorstr(4,js)
      s=pstr(5,js)**2
      isg(1)=iorstr(js)
      iorsg(1)=0
      do 1 i=1,2
      pjt(1,i)=0
      pjt(2,i)=0
      pjt(3,i)=sqrt(s)/2*(3-2*i)
      pjt(4,i)=sqrt(s)/2
      pjt(5,i)=0
      j0=2*i-2
      if(rotstr(3,js).lt.0.)j0=4-2*i
      do 3 j=1,2
3     icjt(j,i)=icstr(j0+j,js)
1     continue
      icj1(1)=icjt(1,1)
      icj1(2)=icjt(2,1)
      icj2(1)=icjt(1,2)
      icj2(2)=icjt(2,2)
      call iddeco(icj1,jcj1)
      call iddeco(icj2,jcj2)
      nqj1=0
      nqj2=0
      do 4 nf=1,nflav
      nqj1=nqj1+jcj1(nf,1)-jcj1(nf,2)
4     nqj2=nqj2+jcj2(nf,1)-jcj2(nf,2)
      if(nqj1.ge.0.and.nqj2.ge.0)isi=0
      if(nqj1.ge.0.and.nqj2.lt.0)isi=1
      if(nqj1.lt.0.and.nqj2.ge.0)isi=2
      if(nqj1.lt.0.and.nqj2.lt.0)isi=3
      nqjstr=isi*10**6+iabs(nqj1)*10**3+iabs(nqj2)
 
c     entry string fragmentation
c     --------------------------
      nptl2=nptl
      ncnt2=0
9996  ncnt2=ncnt2+1
      nptl=nptl2
      nsg=1
      if(ncnt2.gt.1000)goto1001
 
c     entry substring processing
c     --------------------------
      ns=0
9999  ns=ns+1
      if(ns.gt.nsg)goto9997
      if(ish.ge.92)write(ifch,*)'entry substring processing'
      if(ish.ge.92)write(ifch,*)' '
      if(ish.ge.92)write(ifch,107)ns
     *,(icjt(j,2*ns-1 ),j=1,2),(sngl(pjt(j,2*ns-1 )),j=1,5)
107   format(2x,i3,3x,3x,2i7,3x,5(e8.2,1x))
      if(ish.ge.92)write(ifch,107)ns
     *,(icjt(j,2*ns   ),j=1,2),(sngl(pjt(j,2*ns   )),j=1,5)
      if(ish.ge.92)write(ifch,101)ns,(psg(j,ns),j=1,5)
101   format(2x,i3,23x,5(e8.2,1x)/)
114   format(' left breaking pnt: ',8x,4f7.2/)
 
c     determine id and min.mass
c     -------------------------
      delrer=delrex+rangen()*(delrem-delrex)
      do 6 i=1,2
      icp(i)=icjt(i,2*ns-1)
6     icm(i)=icjt(i,2*ns)
      call idcomk(icp)
      call idcomk(icm)
      do 29 i=1,2
      icjt(i,2*ns-1)=icp(i)
29    icjt(i,2*ns)=icm(i)
      call iddeco(icp,jcp)
      call iddeco(icm,jcm)
      amms=utamnx(jcm,jcp)
      nubar=0
      do 7 nf=1,nflav
      jc(nf,1)=jcp(nf,1)+jcm(nf,1)
      jc(nf,2)=jcp(nf,2)+jcm(nf,2)
      if(nf.gt.4.and.(jc(nf,1).ne.0.or.jc(nf,2).ne.0))
     *call utstop('jamfra: flavour > 4                     ')
7     nubar=nubar+jc(nf,1)-jc(nf,2)
      if(ns.eq.1)numbar=nubar
      call idenco(jc,ic,ireten)
      call idcomk(ic)
      id=idtra(ic,0,0,3)
      idk=id
      am=psg(5,ns)
      iremn=0
      if(icp(1)+icp(2).eq.0.or.icm(1)+icm(2).eq.0)iremn=1
      if(am.gt.amms+delrer.and.iremn.eq.0)id=0
 
c     mass adjustment => fragment again
c     ---------------------------------
      amc=psg(5,ns)
      call idres(id,amc,idr,iadj)
      if(abs(amc-psg(5,ns)).gt.1e-3)then
      if(ish.ge.90)then
      call utmsg('jamfra')
      write(ifch,*)'*****  mass changed. fragment again.'
      write(ifch,*)'mass before, after: ',psg(5,ns),amc
      write(ifch,*)'ic,id,idr: ',ic,id,idr
      call utmsgf
      endif
      goto9996
      endif
 
c     if mass < min.mass ==> again
c     ----------------------------
c-c        if(idr.eq.0.and.am.lt.amms-1e-3)then
c-c   if(ish.ge.90)then
c-c   call utmsg('jamfra')
c-c   write(ifch,*)'*****  mass < min mass. ',am,amms
c-c   write(ifch,*)'idk,id,idr,psg(5,ns):'
c-c   write(ifch,*)idk,id,idr,psg(5,ns)
c-c   write(ifch,*)'p_jet:'
c-c   write(ifch,107)ns
c-c  *,(icjt(j,2*ns-1 ),j=1,2),(sngl(pjt(j,2*ns-1 )),j=1,5)
c-c   write(ifch,107)ns
c-c  *,(icjt(j,2*ns   ),j=1,2),(sngl(pjt(j,2*ns   )),j=1,5)
c-c   write(ifch,*)'p_str:'
c-c   write(ifch,101)ns,(psg(j,ns),j=1,5)
c-c   call utmsgf
c-c   endif
c-c   if(ns.eq.1)goto1001
c-c   goto9996
c-c        endif
 
c     write on /cptl/ (1)
c     -------------------
      iori=isg(ns)
      nptlb=nptl
9995  nptl=nptl+1
      nptlc=nptl
      if(nptl.gt.mxptl)
     *call utstop('jamfra: nptl>mxptl                      ')
      if(nptl-nptl2.gt.mxpc)
     *call utstop('jamfra: nptl-nptl2>mxpc                 ')
      do 40 j=1,5
40    pptl(j,nptl)=psg(j,ns)
      do 38 j=1,4
38    xorptl(j,nptl)=xorsg(j,ns)
      xbkptl(1,nptl-nptl2)=ainfin
      xbkptl(2,nptl-nptl2)=ainfin
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=ainfin
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      iclptl(nptl)=0
      isptl(nptl-nptl2)=0
      iorptl(nptl)=iori
      jorptl(nptl)=0
      if(nptl.eq.nptl2+1)nqjptl(nptl)=nqjstr
      if(nptl.gt.nptl2+1)nqjptl(nptl)=-nqjstr
      idptl(nptl)=idr
           if(idr.eq.0)then
      if(ic(1).eq.0.and.ic(2).eq.0)then
      idptl(nptl)=7*10**8
      call idtrbi(jc,ibptl(1,nptl),ibptl(2,nptl)
     *,ibptl(3,nptl),ibptl(4,nptl))
      else
      call idtrbi(jc,idum1,idum2,idum3,idum4)
      idptl(nptl)=8*10**8+ic(1)*100+ic(2)/100
      endif
           endif
           if(ish.ge.91)then
      write(ifch,115)nptl,iorptl(nptl),idptl(nptl)
     *,(pptl(j,nptl),j=1,5)
115   format(1x,'/cptl/',i4,i6,i10,3x,5(e8.2,1x))
      if(idptl(nptl).eq.7e8)write(ifch,*)(ibptl(i,nptl),i=1,4)
      if(ish.ge.92)write(ifch,*)' '
           endif
      if(idr.ne.0)goto9999
 
      nqp=0
      nqm=0
      nap=0
      nam=0
      np=0
      nm=0
      jp=0
      jm=0
      do 23 nf=1,nflav
      nqp=nqp+jcp(nf,1)
      nqm=nqm+jcm(nf,1)
      nap=nap+jcp(nf,2)
      nam=nam+jcm(nf,2)
      np=np+jcp(nf,1)-jcp(nf,2)
      nm=nm+jcm(nf,1)-jcm(nf,2)
      jp=jp+jcp(nf,1)+jcp(nf,2)
23    jm=jm+jcm(nf,1)+jcm(nf,2)
 
c     quark-cluster
c     -------------
      if(icp(1)+icp(2).eq.0.or.icm(1)+icm(2).eq.0)goto78
      if(am.gt.amms+delrer)goto77
      if(idk.eq.0.and.am.gt.amms+delrec)goto77
      if(idk.ne.0.and.am.gt.amms+delrex)goto77
78    continue
      if(ish.ge.92)write(ifch,*)'quark-cluster'
      if(ish.ge.92)write(ifch,*)' '
      if(icp(1)+icp(2).eq.0.and.icm(1)+icm(2).eq.0)
     *call utstop('jamfra: zero string.                    ')
      goto9999
77    continue
 
c     jet trafos
c     ----------
      if(ish.ge.92)write(ifch,*)'jet trafos'
      if(ish.ge.92)write(ifch,*)' '
      do 75 k=1,4
      arp(k)=pjt(k,2*ns-1)
75    arm(k)=pjt(k,2*ns)
      if(ish.ge.90)call utchm(arp,arm,1)
      nso=iorsg(ns)
      if(nso.gt.0)then
      call utrot2(-1,rotsg(1,nso),rotsg(2,nso),rotsg(3,nso)
     *,arp(1),arp(2),arp(3))
      call utlob2(-1,psg(1,nso),psg(2,nso),psg(3,nso),psg(4,nso)
     *,psg(5,nso),arp(1),arp(2),arp(3),arp(4))
      call utrot2(-1,rotsg(1,nso),rotsg(2,nso),rotsg(3,nso)
     *,arm(1),arm(2),arm(3))
      call utlob2(-1,psg(1,nso),psg(2,nso),psg(3,nso),psg(4,nso)
     *,psg(5,nso),arm(1),arm(2),arm(3),arm(4))
      if(ish.ge.90)call utchm(arp,arm,2)
      call utlob2(1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,arp(1),arp(2),arp(3),arp(4))
      call utlob2(1,psg(1,ns),psg(2,ns),psg(3,ns),psg(4,ns),psg(5,ns)
     *,arm(1),arm(2),arm(3),arm(4))
      if(ish.ge.90)call utchm(arp,arm,3)
      rotsg(1,ns)=(arp(1)-arm(1))/2
      rotsg(2,ns)=(arp(2)-arm(2))/2
      rotsg(3,ns)=(arp(3)-arm(3))/2
      call utrot2(1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,arp(1),arp(2),arp(3))
      call utrot2(1,rotsg(1,ns),rotsg(2,ns),rotsg(3,ns)
     *,arm(1),arm(2),arm(3))
      endif
      do i=1,2
      icup(i)=icjt(i,2*ns-1)
      icum(i)=icjt(i,2*ns)
      enddo
           if(arp(3).lt.0)then
      do i=1,2
      icjt(i,2*ns-1)=icum(i)
      icjt(i,2*ns)=icup(i)
      enddo
      do k=1,4
      pjt(k,2*ns-1)=arm(k)
      pjt(k,2*ns)=arp(k)
      enddo
           else
      do i=1,2
      icjt(i,2*ns-1)=icup(i)
      icjt(i,2*ns)=icum(i)
      enddo
      do k=1,4
      pjt(k,2*ns-1)=arp(k)
      pjt(k,2*ns)=arm(k)
      enddo
           endif
      s=psg(5,ns)**2
           if(abs(sngl(pjt(3,2*ns-1)-pjt(3,2*ns))-sqrt(s)).gt.1e-2
     *.and.abs(sngl(pjt(3,2*ns-1)-pjt(3,2*ns))-sqrt(s)).gt.
     *1e-2*sqrt(s))then
      if(ish.ge.90)then
      call utmsg('jamfra')
      write(ifch,*)'*****  pz1-pz2-sqrt(s) nonzero'
      write(ifch,*)'value:   ',sngl(pjt(3,2*ns-1)-pjt(3,2*ns))-sqrt(s)
      write(ifch,*)'sqrt(s): ',sqrt(s)
      write(ifch,*)(sngl(pjt(i,2*ns-1)),i=1,4)
      write(ifch,*)(sngl(pjt(i,2*ns)),i=1,4)
      call utmsgf
      endif
           endif
      if(ish.ge.92)write(ifch,103)ns,nso
     *,(icjt(j,2*ns-1 ),j=1,2),(sngl(pjt(j,2*ns-1 )),j=1,5)
103   format(2x,i3,i3,3x,2i7,3x,5(e8.2,1x))
      if(ish.ge.92)write(ifch,104)
     * (icjt(j,2*ns   ),j=1,2),(sngl(pjt(j,2*ns   )),j=1,5)
104   format(2x,9x,2i7,3x,5(e8.2,1x)/)
 
      if(iopbrk.eq.1)call jambr1(ns,nsg,iret)   ! AMOR
      if(iopbrk.eq.2)call jambr2(ns,nsg,iret)   ! SAMBA
      if(iret.eq.9996)goto9996
      goto9999
 
c     write on /cptl/ (2)
c     -------------------
9997  continue
      if(nptl.le.nptl2)
     *call utstop('jamfra: nptl<=nptl2                     ')
      if(ish.ge.91)write(ifch,*)' '
      do 67 i=nptl2+1,nptl
      istptl(i)=isptl(i-nptl2)
      io=iorptl(i)
      if(io.gt.0)then
      if(ifrptl(1,io).eq.0)ifrptl(1,io)=i
      ifrptl(2,io)=i
      else
      iorptl(i)=0
      endif
      zor=xorsg(3,1)
      tor=xorsg(4,1)
      r=rangen()
      tauran=-taurea*alog(r)
      call uttaix(i,tauran,zor,tor,zis,tis)
      tivptl(1,i)=amax1(tis,xorptl(4,i))
           if(istptl(i).ne.0)then
      taubr=
     *dsqrt((xbkptl(2,i-nptl2)-tor)**2-(xbkptl(1,i-nptl2)-zor)**2)
      call uttaix(i,taubr,zor,tor,zis,tis)
      tivptl(2,i)=tis
           elseif(io.gt.0)then
      taubr=
     *dsqrt((xbkptl(2,io-nptl2)-tor)**2-(xbkptl(1,io-nptl2)-zor)**2)
      call uttaix(i,taubr,xorsg(3,1),xorsg(4,1),zisbr,tisbr)
      call idtau(idptl(i),pptl(4,i),pptl(5,i),taugm)
      tivptl(2,i)=tisbr+taugm
           else
      call idtau(idptl(i),pptl(4,i),pptl(5,i),taugm)
      tivptl(2,i)=xorptl(4,i)+taugm
           endif
      if(ish.ge.92)write(ifch,120)i,iorptl(i),idptl(i)
     *,(sngl(xbkptl(j,i-nptl2)),j=1,2)
     *,(xorptl(j,i),j=3,4),(tivptl(j,i),j=1,2)
120   format(1x,'/cptl/',i4,i5,i10
     *,2x,e8.2,1x,e8.2,1x,e8.2,1x,e8.2,1x,e8.2,1x,e8.2)
67    continue
 
1000  return
 
1001  newevt=1
      goto1000
 
      end
 
c----------------------------------------------------------------------
      subroutine jcentd
c----------------------------------------------------------------------
      parameter (kparx=15,nquax=12)
      common/centro/entro(1+kparx,1+nquax)
      real entrox(1+kparx,1+nquax)
 
      do 1 k=1,1+kparx
      do 1 n=1,1+nquax
1     entro(k,n)=entrox(k,n)
 
      if(kparx.ne.15)
     *call utstop('icentd: insufficient initialization;   k')
      if(nquax.ne.12)
     *call utstop('icentd: insufficient initialization;   n')
      data (entrox( 1,1+n),n=0,12)/
     *.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,
     *.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,
     *.00000E+00/
      data (entrox( 4,1+n),n=0,12)/
     *.40254E+01,.59349E+01,.74855E+01,.87464E+01,.98267E+01,.10770E+02,
     *.11611E+02,.12368E+02,.13058E+02,.13691E+02,.14277E+02,.14822E+02,
     *.15332E+02/
      data (entrox( 7,1+n),n=0,12)/
     *.61944E+01,.89306E+01,.11114E+02,.12969E+02,.14595E+02,.16054E+02,
     *.17380E+02,.18599E+02,.19728E+02,.20780E+02,.21767E+02,.22696E+02,
     *.23574E+02/
      data (entrox(10,1+n),n=0,12)/
     *.68876E+01,.10183E+02,.12855E+02,.15159E+02,.17201E+02,.19046E+02,
     *.20734E+02,.22296E+02,.23751E+02,.25115E+02,.26400E+02,.27616E+02,
     *.28771E+02/
      data (entrox(13,1+n),n=0,12)/
     *.61944E+01,.99602E+01,.13068E+02,.15784E+02,.18208E+02,.20406E+02,
     *.22425E+02,.24295E+02,.26041E+02,.27681E+02,.29228E+02,.30696E+02,
     *.32092E+02/
      data (entrox(16,1+n),n=0,12)/
     *.40254E+01,.82375E+01,.11781E+02,.14923E+02,.17745E+02,.20311E+02,
     *.22669E+02,.24853E+02,.26890E+02,.28803E+02,.30609E+02,.32320E+02,
     *.33948E+02/
 
      return
      end
 
c----------------------------------------------------------------------
      subroutine jcentp
c     plots entro(,) and ftn sjcent
c----------------------------------------------------------------------
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (kparx=15,nquax=12)
      common/centro/entro(1+kparx,1+nquax)
      real xplot(101),yplot(101)
      character text*50
      common/cjcent/nsymx,igx
      integer iplot(5)
      data iplot/0,0,0,1,1/
      ish0=ish
      if(ishsub/100.eq.13)ish=mod(ishsub,100)
      if(.not.(ish.ge.95))goto1000
 
      if(iplot(1).eq.1)then
      text='title ig=     nsym=     k= 3,6,9,12,15$       '
      write(text(10:11),122)igx
      write(text(20:21),122)nsymx
122   format(i2)
      do 54 kx=1,5
      k=3*kx
      do 55 n=0,12
      xplot(1+n)=n
55    yplot(1+n)=entro(1+k,1+n)
      call uthist(0.,12.,0.,40.,1+12,xplot,yplot,'lin','linlin'
     *,'xaxis energy / omega            $                 '
     *,'yaxis entropy                   $                 ',text)
54    continue
      endif
 
      if(iplot(2).eq.1)then
      text='title ig=     nsym=     k= 3,6,9,12,15$       '
      write(text(10:11),122)igx
      write(text(20:21),122)nsymx
      do 64 kx=1,5
      k=3*kx
      do 65 n=0,12
      xplot(1+n)=n*1./k
65    yplot(1+n)=entro(1+k,1+n)/k
      call uthist(0.,5.,0.,6.,1+12,xplot,yplot,'lin','linlin'
     *,'xaxis energy / omega / particle $                 '
     *,'yaxis entropy / particle        $                 ',text)
64    continue
      endif
 
      if(iplot(3).eq.1)then
      do 62 kx=1,10
      k=3*kx
      do 63 n=0,100
      x=n/20.
      xplot(1+n)=x
63    yplot(1+n)=sjcent(k,k,k*x)/k
      call uthist(0.,5.,0.,6.,101,xplot,yplot,'lin','linlin'
     *,'xaxis energy / omega / particle $                 '
     *,'yaxis entropy / particle        $                 '
     *,'title inter(extra)polated       $                 ')
62    continue
      endif
 
      if(iplot(4).eq.1)then
      do 58 kx=1,10
      k=3*kx
      do 59 n=0,100
      x=n/2.
      xplot(1+n)=x
59    yplot(1+n)=sjcent(k,k,x)
      call uthist(0.,50.,0.,80.,101,xplot,yplot,'lin','linlin'
     *,'xaxis energy / omega            $                 '
     *,'yaxis entropy                   $                 '
     *,'title inter(extra)polated  ku=k $                 ')
58    continue
      endif
 
      if(iplot(5).eq.1)then
      do 48 kx=1,10
      k=3*kx
      do 49 n=0,100
      x=n/2.
      xplot(1+n)=x
49    yplot(1+n)=sjcent(k,0,x)
      call uthist(0.,50.,0.,80.,101,xplot,yplot,'lin','linlin'
     *,'xaxis energy / omega            $                 '
     *,'yaxis entropy                   $                 '
     *,'title inter(extra)polated  ku=0 $                 ')
48    continue
      endif
 
1000  continue
      ish=ish0
      return
      end
 
c----------------------------------------------------------------------
      subroutine jcentr(nsym,ig,idi,init)
c     fills array degen(1+k,1+n) containing
c       the number of mixed symmetric (box of width nsym)
c       k-particle states with energy n (units of omega), based on
c       oscillator (dim: idi) wavefunctions (addit. degeneracy: ig) .
c     fills array entro(1+k,1+n) = alog( degen(1+k,1+n) )  .
c     entro(,) written in form of data for jcentd if ish=95
c       (ishsub=12.. selects this).
c     init must be set 1 for the first call of jcentr
c       (init=1 calls uttucl, utpart, utovel)
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (kparx=15,nquax=12)
      parameter (nymax=1000)
      parameter (jpamax=nquax*nquax*nquax,kkpmax=nquax*nquax*2)
      common/cdegen/degen(1+kparx,1+nquax)
      common/centro/entro(1+kparx,1+nquax)
      character*1 deli(1+kparx,1+nquax)
      real xplot(1+nquax),yplot(1+nquax)
      common/cparta/parta(nquax),ipart(nquax,jpamax)
      integer iyo(nymax,kparx),kkk(kparx),iyox(kparx)
     *,jyo(nymax)
      character*1 cyox(1+nquax,kparx)
      double precision tucl
      parameter (ktumax=100,ntumax=100)
      common/ctucl/tucl(1+ktumax,1+ntumax)
      parameter (itamax=1000,nsymax=20)
      integer iyoz(itamax,nsymax,2),iyol(itamax),iyor(nsymax)
     *,iyoo(itamax),iyom(itamax)
      real yofa(itamax)
      parameter (iovmax=100,jovmax=100)
      common/covel/ovel(1+iovmax,1+jovmax)
      character ciga*7
      common/cjcent/nsymx,igx
 
      if(init.eq.1)then
      call uttucl
      call utpart
      call utovel
      endif
 
      write(ifmt,*)'execute sr jcentr ...'
      if(nsym.ge.nsymax)
     *call utstop('jcentr: dimension nsymax too small      ')
      if(mod(kparx,nsym).ne.0)
     *call utstop('jcentr: kparx should be multiple of nsym')
      mmax=1+(kparx-1)/nsym
      izero=0
      nsymx=nsym
      igx=ig
 
c     loop over energy
c     ----------------
      do 6001 n=0,nquax
      write(ifmt,*)'sr jcentr: n=',n
      if(n.eq.0)then
      kkp=1
      else
      kkp=-ipart(n,1)
      if(kkp.gt.kkpmax)
     *call utstop('jcentr: dimension kkpmax too small      ')
      if(-ipart(n,1).ne.nint(parta(n)))
     *call utstop('jcentr: # of partitions disagree        ')
      endif
      if(ish.ge.93)then
      write(ifch,103)('-',ic=1,79),n,kkp,('-',ic=1,79)
103   format(/79a1/7x,'n = ',i2,'   --->   ',i4,' partition(s)'/79a1/)
      jj=2
      do 8 kk=1,kkp
      if(n.gt.0)ll=-ipart(n,jj)
      if(n.eq.0)write(ifch,102)kk,izero
      if(n.gt.0)
     *write(ifch,102)kk,(ipart(n,i),i=jj+1,jj+min0(20,ll))
102   format(2x,i2,'. partition:',3x,20i3)
8     jj=jj+ll+1
      endif
 
c     zero particles
c     --------------
      degen(1,1+n)=1
      entro(1,1+n)=0
      deli(1,1+n)=','
      if(n.eq.nquax)deli(1,1+n)='/'
 
c     loop over particle number
c     -------------------------
      do 6001 kx=1,kparx/3
      k=3*kx
      if(k.gt.kparx)
     *call utstop('jcentr: dimension kparx too small       ')
      degen(1+k,1+n)=0
      deli(1+k,1+n)=','
      if(n.eq.nquax)deli(1+k,1+n)='/'
      kkk(k)=kkp
      mmaxk=1+(k-1)/nsym
      if(k.lt.kparx)then
      do 5 ny=1,nymax
      do 5 l=k+1,kparx
5     iyo(ny,l)=-1
      endif
      if(ish.ge.93)write(ifch,100)('-',ic=1,11),k,('-',ic=1,11)
100   format(/11a1/3x,'k = ',i2/11a1/)
 
c     loop over partitions
c     --------------------
      jj=2
      do 6002 kk=1,kkp
      if(n.gt.0)ll=-ipart(n,jj)
      if(n.eq.0)ll=1
      if(ll.gt.k)then
      kkk(k)=kk-1
      goto6003
      endif
      if(ish.ge.93.and.n.eq.0)
     *write(ifch,102)kk,izero
      if(ish.ge.93.and.n.gt.0)
     *write(ifch,102)kk,(ipart(n,i),i=jj+1,jj+min0(20,ll))
      if(ish.ge.93)write(ifch,*)' '
 
c     construct young tableaus
c     ------------------------
      ny1=1
      ny2=1
      jyo(ny1)=1
      do 46 l=1,k
46    iyo(ny1,l)=0
           if(n.gt.0)then
      iyo(ny1,k)=ipart(n,jj+1)
      ny1=ny1-1
45    ny1=ny1+1
      if(jyo(ny1).eq.ll)goto51
      do 43 lx=1,nsym
      do 44 mx=1,mmaxk
      l=nsym+1-lx
      m=mmaxk+1-mx
           if(iyo(ny1,(m-1)*nsym+l).eq.0)then
      if(l.eq.nsym.or.(l.ne.nsym.and.iyo(ny1,(m-1)*nsym+l+1).ne.0))then
      ny2=ny2+1
      if(ny2.gt.nymax)
     *call utstop('jcentr: dimension nymax too small       ')
      jy=jyo(ny1)+1
      jyo(ny2)=jy
      do 50 lp=1,k
50    iyo(ny2,lp)=iyo(ny1,lp)
      iyo(ny2,(m-1)*nsym+l)=ipart(n,jj+jy)
      if(ny2.gt.1)then
      do 47 ny3=1,ny2-1
      if(jyo(ny3).ne.jyo(ny2))goto47
      do 48 lp=1,k
48    if(iyo(ny3,lp).ne.iyo(ny2,lp))goto47
      ny2=ny2-1
      goto49
47    continue
49    continue
      endif
      endif
      goto43
           endif
44    continue
43    continue
      goto45
51    continue
          endif
 
c     loop over young tableaus
c     ------------------------
      do 6005 ny=ny1,ny2
 
      if(ish.ge.93.and.nsym.eq.3)then
      write(ifch,117)ny-ny1+1,((iyo(ny,(m-1)*nsym+i),i=1,nsym),m=1,1)
117   format(2x,i2,'. tableau:',5x,3i2)
      if(mmaxk.gt.1)
     *write(ifch,110)((iyo(ny,(m-1)*nsym+i),i=1,nsym),m=2,mmaxk)
110   format(19x,3i2)
      write(ifch,*)' '
      endif
 
c     select nn-subtableaus
c     ---------------------
      ditab=1
      do 6004 nn=0,n
      iga=ig*nint(dexp(tucl(1+idi,1+nn)))
 
      do 15 l=1,kparx
15    iyox(l)=-1
      do 11 m=1,mmaxk
      memp=m-1
      ms=(m-1)*nsym
      do 11 i=1,nsym
      if(iyo(ny,ms+i).eq.nn)goto12
11    continue
12    continue
      mmaxks=mmaxk-memp
      do 13 l=ms+1,k
      if(iyo(ny,l).eq.nn)iyox(l-ms)=iyo(ny,l)
13    continue
      if(nsym.gt.1)then
      do 17 i=1,nsym-1
      do 16 m=1,mmaxk
      ms=(m-1)*nsym
      if(iyox(ms+1).ne.-1)goto18
16    continue
      do 19 m=1,mmaxk
      ms=(m-1)*nsym
      do 20 l=1,nsym-1
20    iyox(ms+l)=iyox(ms+l+1)
      iyox(ms+nsym)=-1
19    continue
17    continue
18    continue
      endif
      do 38 l=1,kparx
      cyox(1+nn,l)=' '
      if(iyox(l).ge.0)cyox(1+nn,l)='X'
38    continue
 
c     skip for empty tableaus
c     -----------------------
      do 21 l=1,k
      if(iyox(l).ne.-1)goto22
21    continue
      goto6004
22    continue
 
c     print
c     -----
      if(ish.ge.93.and.nsym.eq.3)then
      if(iga.ge.  1.and.iga.lt.  10)write(ciga,107)iga
      if(iga.gt. 10.and.iga.lt. 100)write(ciga,108)iga
      if(iga.gt.100.and.iga.lt.1000)write(ciga,109)iga
107   format('GL(',i1,')  ')
108   format('GL(',i2,') ')
109   format('GL(',i3,')')
      write(ifch,116)nn,((cyox(1+nn,(m-1)*nsym+i),i=1,nsym),m=1,1)
     *,ciga
116   format(2x,i2,'-subtableau:',3x,3(1x,a1),3x,a7)
      if(mmaxks.gt.1)then
      do 39 m=2,mmaxks
39    write(ifch,106)(cyox(1+nn,(m-1)*nsym+i),i=1,nsym)
106   format(19x,3(1x,a1))
      endif
      if(ish.ge.94)write(ifch,*)' '
      endif
 
c     reductions of tableaus 1+iyox()
c     -------------------------------
      disuta=0
      n1=1
      n3=1
      n300=1
      do 24 i=1,nsym
      ma=0
      me=0
      do 25 m=1,mmaxk
      if(1+iyox((m-1)*nsym+i).le.0)goto25
      if(ma.eq.0)ma=m
      me=m
25    continue
      iyoz(n1,i,1)=ma
      iyoz(n1,i,2)=me
24    continue
      iyol(n1)=0
      yofa(n1)=1
      iyoo(n1)=0
      iyom(n1)=1
      n1=n1-1
5000  n1=n1+1
      n30=n3
      if(n1.gt.1.and.iyol(n1).gt.iyol(n1-1))n300=n3
      n2=n1
      do 26 i=1,nsym
26    iyor(i)=max0(0,iyoz(n1,i,2)-1)
      lev=iyol(n1)
      fac=yofa(n1)
      if(iga-1.le.iovmax.and.lev.le.jovmax.and.lev.le.iga-1)then
      facx=exp(ovel(iga,1+lev))
      if(abs(facx-fac).gt.1e-5*fac)then
      write(ifch,*)' '
      write(ifch,*)'n=',n,'    k=',k,'    kk=',kk,'    ny=',ny
      write(ifch,*)'facx=',facx,'     fac=',fac
      call utstop('jcentr: binomials differ                ')
      endif
      endif
      if(lev.gt.iga-1)goto5003
      iheim=0
      do 32 i=1,nsym
      ihei=iyoz(n1,i,2)-iyoz(n1,i,1)
      if(ihei.gt.iheim)iheim=ihei
32    continue
      if(iheim.eq.0)disuta=disuta+fac*iyom(n1)
           if(ish.ge.94.and.nsym.eq.3)then
      write(ifch,112)n1,(iyoz(n1,i,1),i=1,nsym),iyoo(n1)
     *,lev,fac,iyom(n1),iheim,disuta
     *,(iyoz(n1,i,2),i=1,nsym)
112   format(3x,i3,2x,3i2,2x,'ori:',i3,3x,'lev:',i2,3x,'fac:',f8.1
     *,3x,'mul:',i2,3x,'hei:',i2,3x,'sum:',f8.1/8x,3i2/)
           endif
      n2=n2-1
5001  n2=n2+1
      do 27 i=1,nsym
      if(iyoz(n2,i,2)-1.lt.iyor(i))goto27
      if(i.lt.nsym.and.iyoz(n2,i,2)-1.lt.iyoz(n2,i+1,2))goto27
      n3=n3+1
      if(n3.gt.itamax)
     *call utstop('jcentr: dimension itamax too small      ')
      iyol(n3)=iyol(n1)+1
      yofa(n3)=( yofa(n1)*(iga-iyol(n3)) )/iyol(n3)
      iyoo(n3)=n1
      iyom(n3)=iyom(n1)
      do 28 j=1,nsym
      iyoz(n3,j,1)=iyoz(n2,j,1)
28    iyoz(n3,j,2)=iyoz(n2,j,2)
      iyoz(n3,i,2)=iyoz(n3,i,2)-1
      if(iyoz(n3,i,2).lt.iyoz(n3,i,1))then
      iyoz(n3,i,1)=0
      iyoz(n3,i,2)=0
      endif
      if(n30+1.lt.n3)then
      do 30 nch=n30+1,n3-1
      do 31 ij=1,nsym
      if(iyoz(nch,ij,1).ne.iyoz(n3,ij,1)
     *.or.iyoz(nch,ij,2).ne.iyoz(n3,ij,2))goto30
31    continue
      n3=n3-1
      goto27
30    continue
      endif
           if(ish.ge.95.and.nsym.eq.3)then
      write(ifch,113)n1,(iyoz(n1,ij,1),ij=1,nsym)
     *,n2,(iyoz(n2,ij,1),ij=1,nsym),n3,(iyoz(n3,ij,1),ij=1,nsym)
     *,(iyoz(n1,ij,2),ij=1,nsym)
     *,(iyoz(n2,ij,2),ij=1,nsym),(iyoz(n3,ij,2),ij=1,nsym)
113   format(3x,'n1: ',i2,3x,3i2,4x,'n2: ',i2,3x,3i2,4x,'n3: ',i2,3x,3i2
     */12x,3i2,13x,3i2,13x,3i2/)
           endif
27    continue
      if(n2.eq.n1)n2=n30
      if(n2.lt.n3)goto5001
           if(n30.lt.n3.and.n300.lt.n30)then
      n3s=n3
      n3m=n3
      n3=n30
      do 33 nch3=n30+1,n3s
      n3=n3+1
      do 34 nch=n300+1,n30
      do 35 ij=1,nsym
      if(iyoz(nch,ij,1).ne.iyoz(n3,ij,1)
     *.or.iyoz(nch,ij,2).ne.iyoz(n3,ij,2))goto34
35    continue
      iyom(nch)=iyom(nch)+iyom(n3)
      n3=n3-1
      n3m=n3m-1
      if(n3+1.le.n3m)then
      do 36 nm=n3+1,n3m
      iyol(nm)=iyol(nm+1)
      yofa(nm)=yofa(nm+1)
      iyoo(nm)=iyoo(nm+1)
      iyom(nm)=iyom(nm+1)
      do 37 j=1,nsym
      iyoz(nm,j,1)=iyoz(nm+1,j,1)
37    iyoz(nm,j,2)=iyoz(nm+1,j,2)
36    continue
      endif
      goto33
34    continue
33    continue
           endif
      if(n1.lt.n3)goto5000
5003  continue
 
      dialt=-99999.
      if(iyoz(1,1,1).eq.1)then
      dialt=1
      do 40 i=1,nsym
      if(iyoz(1,i,1).gt.1)
     *call utstop('jcentr: iyoz(,,1).gt.1                  ')
      if(iyoz(1,i,1).eq.0)goto40
      do 41 m=1,iyoz(1,i,2)
      hak=iyoz(1,i,2)-m+1
      if(i.lt.nsym)then
      do 42 j=i+1,nsym
      if(iyoz(1,j,2).ge.m)hak=hak+1
42    continue
      endif
      dialt=dialt*(iga+i-m)/hak
41    continue
40    continue
      endif
 
      if(dialt.ne.-99999.
     *.and.abs(dialt-disuta).gt.1e-5*dialt)then
      write(ifch,*)' '
      write(ifch,*)'n=',n,'    k=',k,'    kk=',kk,'    ny=',ny
      write(ifch,*)'disuta=',disuta,'     dialt=',dialt
      call utstop('jcentr: dimensions differ               ')
      endif
 
      if(ish.ge.93)write(ifch,118)disuta
118   format(28x,'subtab.-dimension:',f8.1/)
 
      ditab=ditab*disuta
 
6004  continue
 
      if(ish.ge.93)write(ifch,119)ditab
119   format(35x,'tableau-dimension:',f11.1/)
 
      degen(1+k,1+n)=degen(1+k,1+n)+ditab
 
6005  continue
 
6002  jj=jj+ll+1
6003  continue
      if(kkk(k).eq.0)
     *call utstop('jcentr: no allowed partition            ')
 
      entro(1+k,1+n)=-9999999.
      if(degen(1+k,1+n).gt.0.)entro(1+k,1+n)=alog(degen(1+k,1+n))
 
      if(ish.ge.93)write(ifch,120)degen(1+k,1+n)
120   format(49x,'degeneracy:',f11.1/)
 
6001  continue
 
      if(ish.ge.92)then
121   format(79a1)
      write(ifch,*)' '
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   degeneracy(k,n)         ig=',ig,'          '
     *,'nsym=',nsym
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   k:','   3','   6','   9','   12'
      write(ifch,121)('-',ic=1,79)
      do 52 n=0,nquax
52    write(ifch,*)n,(degen(1+3*k,1+n),k=1,4)
      write(ifch,*)' '
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   entropy(k,n)            ig=',ig,'          '
     *,'nsym=',nsym
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   k:','   3','   6','   9','   12'
      write(ifch,121)('-',ic=1,79)
      do 53 n=0,nquax
53    write(ifch,*)n,(entro(1+3*k,1+n),k=1,4)
      endif
 
      ish0=ish
      if(ishsub/100.eq.12)ish=mod(ishsub,100)
      if(ish.ge.95)then
      write(ifch,*)' '
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   entropy(k,n)            ig=',ig
     *,'          nsym=',nsym
      write(ifch,121)('-',ic=1,79)
      write(ifch,125)kparx
125   format(6x,'if(kparx.ne.',i2,')'/5x,
     *56h*call utstop('jcentd: insufficient initialization;   k'))
      write(ifch,126)nquax
126   format(6x,'if(nquax.ne.',i2,')'/5x,
     *56h*call utstop('jcentd: insufficient initialization;   n'))
      do 58 kx=0,kparx/3
      k=3*kx
      write(ifch,123)k+1,nquax
123   format(6x,'data (entro(',i2,',1+n),n=0,',i2,')/')
      write(ifch,124)(entro(1+k,1+n),deli(1+k,1+n), n=0,nquax)
124   format(10(5x,'*',6(e10.5,a1)/))
58    continue
      endif
      ish=ish0
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jclude(ip,iret)
c-----------------------------------------------------------------------
c  decays cluster ip from /cptl/ .
c  reqires jcentr or jcentd to be called before.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter (nflav=6)
      integer ica(2),ice(2),ico(2),ica0(2)
     *,jca(nflav,2),jce(nflav,2),jco(nflav,2),jca0(nflav,2)
      common/csjcga/amnull,nsuha(7),asuha(7),entrpy,amegam,nopha
      common/centex/entexp
      common/cscal/scal
      common/ctimel/ntc
      common/cutinv/lutinv
  
      parameter (moxmax=30,iomaxm=25,iomax=54)
      parameter (momax=moxmax*iomax)
      integer idol(iomax),idolis(iomax),ifolis(iomax)
      real pol(iomax)
      real xo(1+momax),qo(1+momax),yo(1+momax),oo(1+momax),go(1+momax)
      real pa(5),pa0(5),pe(5),u(3)
      integer ko(1+momax)

      data idolis/
     * 110,120,-120,220,130,-130,230,-230,330
     *,1120,-1120,1220,-1220,1130,-1130,1230,-1230,2130,-2130
     *,2230,-2230,1330,-1330,2330,-2330
     *,1111,-1111,2221,-2221,3331,-3331
     *,111,121,-121,221,131,-131,231,-231,331
     *,1121,-1121,1221,-1221,1131,-1131,1231,-1231
     *,2231,-2231,1331,-1331,2331,-2331/
      data ifolis/9*1,16*2,6*4,9*3,14*4/

      ish0=ish
      if(ishsub/100.eq.1)ish=mod(ishsub,100)
      if(ish.ge.92)then
      write(ifch,*)' '
      write(ifch,*)('-',l=1,79)
      write(ifch,*)'   cluster decay of',ip,idptl(ip),pptl(5,ip)
      write(ifch,*)('-',l=1,79)
      write(ifch,*)' '
      endif
      delpox=.01
      nptlb=nptl
      loop=0
      iret=0
      ipoxra=0
      istpfl=0
      ebamin=1.5 
      loopmx=20
      nophax=20

c     original cluster --> pa,jca
c     ---------------------------

6010  continue
      do 1 j=1,5
1     pa(j)=pptl(j,ip)
      call idquac(ip,ndu,ndu,ndu,jca)
      if(ish.ge.92)then
      write(ifch,*)'initial jca:'
      write(ifch,*)jca
      endif
      call idcomj(jca)
      if(ish.ge.92)then
      write(ifch,*)'initial jca after compactification:'
      write(ifch,*)jca
      endif
 
      do 18 j=1,5
18    pa0(j)=pa(j)
      nqn=0
      do 2 nf=1,nflav
      nqn=nqn+jca(nf,1)-jca(nf,2)
      do 2 j=1,2
2     jca0(nf,j)=jca(nf,j)
      nba=iabs(nqn)/3
      eba=ainfin
      if(nba.gt.0)eba=pa(5)/nba
      ish00=ish

c     initialization for redo
c     -----------------------

6001  loop=loop+1
           if(ish00.eq.90)then
      if(loop.eq.loopmx)then
      ish=92
      write(ifch,117)('-',l=1,75),loop,ish
117   format(//1x,75a1/1x,i2,'. attempt to decay this cluster!'
     */1x,'ish set to: ',i2)
      endif 
           endif
      if(loop.gt.loopmx)then
      ish=ish00
      goto1001
      endif
      nopha=0
      ismall=0
      nptl=nptlb
      do 20 j=1,5
20    pa(j)=pa0(j)
      do 3 nf=1,nflav
      do 3 j=1,2
3     jca(nf,j)=jca0(nf,j)
 
c     loop over sequential decays
c     ---------------------------

6002  nptl=nptl+1
      if(nptl.gt.mxptl)
     *call utstop('jclude: nptl>mxptl                      ')
      isjca=0
      do nf=1,nflav
      do j=1,2
      isjca=isjca+iabs(jca(nf,j))
      enddo
      enddo
      if(isjca.eq.0)then
      jca(1,1)=1
      jca(1,2)=1
      endif
      call idcomj(jca)
      call idenco(jca,ica,ireten)
           if(ish.ge.92)then
      write(ifch,*)' '
      write(ifch,100)(pa(k),k=1,5)
100   format(1x,'pa:',20x,5(e10.3,1x))
      write(ifch,116)jca
116   format(1x,'jca: ',6i5/6x,6i5)
           endif

      kau=jca(1,1)-jca(1,2)
      kad=jca(2,1)-jca(2,2)
      kas=jca(3,1)-jca(3,2)
      kac=jca(4,1)-jca(4,2)
      amicl=utamnu(kau,kad,kas,kac,5)

           if(nopha.eq.nophax.or.nopha.gt.1.and.pa(5).gt.2*amicl)then
      if(ish.ge.90)then
      call utmsg('jclude')
      if(ismall.eq.0)write(ifch,*)'*****  no phase space --> redo decay'
      if(ismall.eq.1)write(ifch,*)'*****  small phase sp --> redo decay'
      if(ismall.ne.0.and.ismall.ne.1)
     *write(ifch,*)'***** small/no phase space --> redo decay'
      if(ismall.ne.0.and.ismall.ne.1)write(ifch,*)'ismall=',ismall
      write(ifch,*)'nt=',ntc,'   loop=',loop,'   nopha=',nopha
      if(ismall.eq.1)
     *write(ifch,*)'yo_hit=',yohit,'   yo_max=',yomax
      write(ifch,*)'m_min=',amicl,'   m=',pa(5)
      write(ifch,100)(pa(k),k=1,5)
      write(ifch,116)jca
      call utmsgf
      endif
      goto6001
           endif

           if(nopha.gt.1)then
      if(ish.ge.90)then
      call utmsg('jclude')
      if(ismall.eq.0)write(ifch,*)'*****  no phase space --> incrs mass'
      if(ismall.eq.1)write(ifch,*)'*****  small phase sp --> incrs mass'
      if(ismall.ne.0.and.ismall.ne.1)
     *write(ifch,*)'***** small/no phase space --> incrs mass'
      if(ismall.ne.0.and.ismall.ne.1)write(ifch,*)'ismall=',ismall
      write(ifch,*)'nt=',ntc,'   loop=',loop,'   nopha=',nopha
      if(ismall.eq.1)
     *write(ifch,*)'yo_hit=',yohit,'   yo_max=',yomax
      write(ifch,*)'m_min=',amicl,'   m=',pa(5)
      write(ifch,100)(pa(k),k=1,5)
      write(ifch,116)jca
      endif
      pa(5)=pa(5)*1.05
      pa(4)=sqrt(pa(1)**2+pa(2)**2+pa(3)**2+pa(5)**2)
      if(ish.ge.90)then
      write(ifch,100)(pa(k),k=1,5)
      call utmsgf
      endif
           endif

      idar=0
           if(nopha.eq.0)then
      ich=0
      ida=0
      if(ireten.eq.0)ida=idtra(ica,0,0,3)
      ama=pa(5)
      call idres(ida,ama,idar,iadj)
      if(idar.ne.0.and.eba.lt.ebamin)then
      ama=0
      call idres(ida,ama,idar,iadj)
      if(mod(idar,10).ne.0.and.loop.lt.loopmx/2)goto6001 
      endif
      if(idar.ne.idptl(ip))ich=1
           endif

           if(idar.ne.0)then
      if(nptl.gt.nptlb+1.or.ich.eq.1)then
      idptl(nptl)=idar
      pa(5)=ama
      pa(4)=sqrt(ama**2+pa(1)**2+pa(2)**2+pa(3)**2)
      do 15 j=1,5
15    pptl(j,nptl)=pa(j)
      if(ish.ge.92)
     *write(ifch,101)nptl,idptl(nptl),(pptl(k,nptl),k=1,5)
101   format(1x,'/cptl/',i6,i11,2x,5(e10.3,1x))
      else
      nptl=nptl-1
      endif
      goto7001
           endif

      ama=pa(5)
 
c     loop over hadrons
c     -----------------

      mo=1
      po=0
      xo(1)=0
      qo(1)=0
      iom=0
      qol=0
      if(ish.ge.92)write(ifch,106)
106   format(' partial decay widths:')

      do 6003 io=1,iomax
      if(eba.lt.ebamin.and.io.gt.iomaxm)goto6003

      ido=idolis(io)
      iom=iom+1
      pol(iom)=po
      idol(iom)=ido
      call idmass(ido,amo)
      idptl(nptl)=ido
      call idquac(nptl,ndu,ndu,ndu,jco)
      do 12 nf=1,nflav
      do 12 j=1,2
12    jce(nf,j)=jca(nf,j)-jco(nf,j)
      do 22 nf=1,nflav
      do 22 j=1,2
      if(jce(nf,j).lt.0)then
      jce(nf,3-j)=jce(nf,3-j)-jce(nf,j)
      jce(nf,j)=0
      endif
22    continue
      ke=0
      do 13 i=5,nflav
      if(jce(i,1).ne.0.or.jce(i,2).ne.0.or.istpfl.gt.0)then
      if(istpfl.eq.2)
     *call utstop('jclude: more than 4 flavours            ')
      istpfl=istpfl+1
      ish=93
      if(istpfl.eq.2)goto6010
      goto6001
      endif
13    continue
      keu=jce(1,1)-jce(1,2)
      ked=jce(2,1)-jce(2,2)
      kes=jce(3,1)-jce(3,2)
      kec=jce(4,1)-jce(4,2)
 
c     loop over hadron momenta
c     ------------------------

      if(mo+moxmax.gt.1+momax)
     *call utstop('jclude: dimension momax too small       ')
      pox=0.
      mom=mo
      pom=po
      qom=qo(mo)
      gamxm=0.
      poxm=0.
      gamy=0.
      gamx=0.
      delpo=delpox
      yoint=0

      do 6004 mox=1,moxmax

      delpo=delpo*1.2
      mo=mo+1
      po=po+delpo
      ko(mo)=ido
      if(mox.gt.1)ko(mo)=0
      xo(mo)=po
      pox=pox+delpo
      gamx=sjcgam(keu,ked,kes,kec,ama,amo,pox,mox)
      go(mo)=gamx
      yo(mo)=ifolis(io)*delpo*gamx
      oo(mo)=entexp 
      yoint=yoint+ifolis(io)*delpo*(gamx+gamy)/2.         
      i6005=0
      if(gamx.eq.0.)i6005=1
      if(mox.gt.1.and.gamxm.gt.0..and.gamx.lt.1e-2*gamxm)i6005=1
           if((i6005.eq.0.or.mox.gt.1).and.ish.ge.93)then
      if(mox.eq.1)
     *write(ifch,109)keu,ked,kes,nsuha,ama,asuha,amnull,ido,amo
      write(ifch,105)mo,xo(mo),pox,sqrt(pox**2+amo**2)
     *,ama-sqrt(pox**2+amo**2),entrpy,gamx
     *,yo(mo),oo(mo)
109   format(/' u_d_s:',3i3,'  n:',7i5,11x,'  a:',f10.2
     */18x,'m:',7f5.2,'  m0:',f6.2,'  o:',i5,f5.2
     */'   mo     xo   pox      eo ama-eo     entrpy     '
     *,'  gamx         yo         oo ')
105   format(i5,1x,f6.2,f6.2,2x,f6.2,f6.2
     *,3x,e9.3,2x,e9.3,2x,e9.3,2x,e9.3)
           endif
      if(i6005.eq.1)goto6005
      if(gamx.gt.gamxm)then
      gamxm=gamx
      poxm=pox
      endif
      gamy=gamx

6004  continue

      if(ish.ge.90.and.ipoxra.eq.0)then
      ipoxra=1
      call utmsg('jclude')
      write(ifch,*)'*****  pox-range too small'
      write(ifch,100)(pa(k),k=1,5)
      write(ifch,*)'jca:'
      write(ifch,*)jca
      write(ifch,*)'pox= ',pox, '   gamx= ',gamx
      write(ifch,*)'poxm=',poxm,'   gamxm=',gamxm
      call utmsgf
      endif

6005  continue
 
      if(yoint.eq.0.)then
      if(ish.ge.94)write(ifch,*)' '
      if(ish.ge.94)write(ifch,*)'no phase space for',ido
      mo=mom
      po=pom
      iom=iom-1
      else
      if(ish.ge.92)write(ifch,107)ido,poxm,yoint,oo(mo)
107   format(' ido,poxm,yoint,oo:',i6,3x,f5.2,3x,e9.3,3x,e9.3)
      endif
 
6003  continue
      if(ish.ge.92)write(ifch,*)' '

c     no phase space
c     --------------
      if(iom.eq.0)then
      if(ish.ge.92)write(ifch,*)'no phase space'
      nptl=nptl-1
      nopha=nopha+1
      goto6002
      endif

c     determine qo(m)
c     ---------------
      oomax=0
      do m=1,mo
      if(oo(m).gt.oomax)oomax=oo(m)
      enddo
      oox=oomax-10
      do m=1,mo
      if(oo(m).lt.oox)then
      yo(m)=0
      else
      yo(m)=yo(m)*exp(oo(m)-oox)
      endif
      enddo       
      qo(1)=0        
      do m=2,mo
      if(ko(m).ne.0)yom=0
      qo(m)=qo(m-1)+(yom+yo(m))/2.
      yom=yo(m)
      enddo
      yomax=0
      do m=1,mo
      if(yo(m).gt.yomax)yomax=yo(m)
      enddo

      if(ish.ge.93)then
      ido=0
      xox=0
      write(ifch,*)' '
      write(ifch,*)'accumulated decay probabilities:'
      write(ifch,*)' '
      do m=1,mo
      if(ko(m).ne.0)write(ifch,*)' '
      if(ko(m).ne.0)ido=ko(m)
      write(ifch,110)m,ido,xo(m)-xox,qo(m)
110   format(1x,'m,ido,xo,qo: ',i6,i6,3x,f6.2,3x,e9.3)  
      if(m.lt.mo.and.ko(m+1).ne.0)xox=xo(m)
      enddo
      endif

c     select randomly hadron+momentum
c     -------------------------------

      lo=0
24    lo=lo+1
           if(lo.gt.10)then
      if(ish.ge.92)write(ifch,*)'small phase space:',yohit,yomax
      nptl=nptl-1
      nopha=nopha+1
      ismall=1
      goto6002
           endif

      xos=utinvt(mo,xo,qo,rangen()*qo(mo))
      mhit=lutinv
      do 28 i=2,iom
      ix=i-1
      if(xos.lt.pol(i))goto29
28    continue
      ix=iom
29    continue
      pos=xos-pol(ix)
      ids=idol(ix)
      call idmass(ids,ams)
      idptl(nptl)=ids

      call idquac(nptl,ndu,ndu,ndu,jco)
      do 25 nf=1,nflav
      do 25 j=1,2
25    jce(nf,j)=jca(nf,j)-jco(nf,j)
      do 26 nf=1,nflav
      do 26 j=1,2
      if(jce(nf,j).lt.0)then
      jce(nf,3-j)=jce(nf,3-j)-jce(nf,j)
      jce(nf,j)=0
      endif
26    continue
      ke=0
      do 27 i=5,nflav
      if(jce(i,1).ne.0.or.jce(i,2).ne.0)
     *call utstop('jclude: flavour > 4                     ')
27    continue

      keu=jce(1,1)-jce(1,2)
      ked=jce(2,1)-jce(2,2)
      kes=jce(3,1)-jce(3,2)
      kec=jce(4,1)-jce(4,2)
      gamxht=sjcgam(keu,ked,kes,kec,ama,ams,pos,1)
      yohit=gamxht*(xo(mhit+1)-xo(mhit))*4
      if(yohit.lt.1e-5.and.yohit.lt.1e-5*yomax)goto24
 
      pptl(5,nptl)=ams
      u(3)=2.*rangen()-1.
      phi=2.*pi*rangen()
      u(1)=sqrt(1.-u(3)**2)*cos(phi)
      u(2)=sqrt(1.-u(3)**2)*sin(phi)
      do 5 j=1,3
5     pptl(j,nptl)=pos*u(j)
      pptl(4,nptl)=sqrt(pptl(5,nptl)**2+pos**2)
      if(ish.ge.93)
     *write(ifch,101)nptl,idptl(nptl),(pptl(k,nptl),k=1,5)
      call utlobo(-1,pa(1),pa(2),pa(3),pa(4),pa(5)
     *,pptl(1,nptl),pptl(2,nptl),pptl(3,nptl),pptl(4,nptl))
      if(ish.ge.92)
     *write(ifch,101)nptl,idptl(nptl),(pptl(k,nptl),k=1,5)
 
      pe(5)=amegam
      do 16 j=1,3
16    pe(j)=-pos*u(j)
      pe(4)=sqrt(pe(5)**2+pos**2)
      call utlobo(-1,pa(1),pa(2),pa(3),pa(4),pa(5)
     *,pe(1),pe(2),pe(3),pe(4))
 
      do 17 j=1,5
17    pa(j)=pe(j)
      do 23 nf=1,nflav
      do 23 j=1,2
23    jca(nf,j)=jce(nf,j)
      nopha=0
 
      goto6002
 
7001  continue
 
c     check energy conservation
c     -------------------------
      if(nptl.le.nptlb)goto1000

      ifail=1
           if(nptl.gt.nptlb+1)then                          
      ishrsc=ish
      ish=0
      call jrescl(nptlb+1,nptl,pa0,ifail)
      ish=ishrsc
c-c   if(ifail.ne.0.and.ish.ge.90)then
c-c   call utmsg('jclude')
c-c   write(ifch,*)'*****  ifail_jrescl=',ifail
c-c   call utmsgf
c-c   endif
           endif

c-c   do n=nptlb+1,nptl
c-c   p=pptl(3,n)
c-c   e=pptl(4,n)
c-c   y=100
c-c   dy=3
c-c   if(e-p.ne.0..and.e+p.ne.0.)y=.5*alog((e+p)/(e-p))+dy
c-c   ida=iabs(idptl(n))
c-c   if(ida.gt.1000.and.mod(ida,10).ne.0.and.y.lt.10.)then
c-c   write(6,*)('-',k=1,69)
c-c   p=pptl(3,ip)
c-c   e=pptl(4,ip)
c-c   y=100
c-c   dy=3
c-c   if(e-p.ne.0..and.e+p.ne.0.)y=.5*alog((e+p)/(e-p))+dy
c-c   write(6,115)ip,idptl(ip)
c-c  *,(pptl(k,ip),k=3,5),y
c-c   write(6,*)'-------> ' 
c-c   do m=nptlb+1,nptl
c-c   p=pptl(3,m)
c-c   e=pptl(4,m)
c-c   y=100
c-c   dy=3
c-c   if(e-p.ne.0..and.e+p.ne.0.)y=.5*alog((e+p)/(e-p))+dy
c-c   write(6,115)m,idptl(m)
c-c  *,(pptl(k,m),k=3,5),y
c-c   enddo
c-c   endif
c-c   enddo
 115  format(1x,'/cptl/',i4,i10
     *,3x,3(e8.2,1x),2x,2(e8.2,1x))

      if(ish.ge. 0)then
      px=0
      py=0
      pz=0
      e=0
      do 10 n=nptlb+1,nptl
      px=px+pptl(1,n)
      py=py+pptl(2,n)
      pz=pz+pptl(3,n)
10    e=e+  pptl(4,n)
      endif
 
      if(ish.ge.92)then
      write(ifch,*)' '
      write(ifch,*)'check energy conservation'
      do 11 n=nptlb+1,nptl
11    write(ifch,101)n,idptl(n),(pptl(k,n),k=1,5)
      am=sqrt(e**2-px**2-py**2-pz**2)
      write(ifch,103)px,py,pz,e,am
      write(ifch,104)(pa0(k),k=1,5)
103   format(1x,'p_sum:  ',15x,5(e10.3,1x))
104   format(1x,'p_clu:  ',15x,5(e10.3,1x))
      endif
 
      if(ish.ge.90.and.(
     *ifail.eq.0.and. 
     *   (abs(px-pa0(1)).gt.1e-2*abs(px).and.abs(px-pa0(1)).gt.1e-2
     *.or.abs(py-pa0(2)).gt.1e-2*abs(py).and.abs(py-pa0(2)).gt.1e-2
     *.or.abs(pz-pa0(3)).gt.1e-2*abs(pz).and.abs(pz-pa0(3)).gt.1e-2
     *.or.abs(e -pa0(4)).gt.1e-2*abs(e ).and.abs(e -pa0(4)).gt.1e-2)
     *.or.
     *ifail.ne.0.and.
     *   (abs(px-pa0(1)).gt.1e-2*abs(px).and.abs(px-pa0(1)).gt.1e-2
     *.or.abs(py-pa0(2)).gt.1e-2*abs(py).and.abs(py-pa0(2)).gt.1e-2
     *.or.abs(pz-pa0(3)).gt.1e-2*abs(pz).and.abs(pz-pa0(3)).gt.1e-2
     *.or.abs(e -pa0(4)).gt.35e-1*abs(e ).and.abs(e -pa0(4)).gt.35e-1)
     *   ))then
      call utmsg('jclude')
      write(ifch,*)'*****  p_sum /= p_clu'
      write(ifch,*)'ifail_jrescl:',ifail,'   scal:',scal
      do 30 n=nptlb+1,nptl
30    write(ifch,101)n,idptl(n),(pptl(k,n),k=1,5)
      write(ifch,103)px,py,pz,e
      write(ifch,104)(pa0(k),k=1,4)
      call utmsgf
      endif
  
1000  continue
      if(ish.ge.92)write(ifch,*)' '
      if(ish.ge.92)write(ifch,*)('-',l=1,25)
      if(ish.ge.92)write(ifch,*)'   return from jclude   '
      if(ish.ge.92)write(ifch,*)('-',l=1,25)
      if(ish.ge.92)write(ifch,*)' '
      ish=ish0
      return
 
1001  iret=1
      if(ish.ge.90)then
      call utmsg('jclude')
      write(ifch,*)'*****  cluster decay not possible  --> iret=1'
      write(ifch,100)(pa0(k),k=1,5)
      write(ifch,*)'jca:'
      write(ifch,*)jca0
      call utmsgf
      endif
      return
 
      end
 
c-----------------------------------------------------------------------
      subroutine jdeca(i,iret)
c     decays i (calls jdecay)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat,zor,tor
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat 
      common/cclude/kclude
      parameter (nflav=6)
      integer jcdu(nflav,2)
      iret=0
      if(istptl(i).ne.0)return
      ish0=ish
      if(ishsub/100.eq.5)ish=mod(ishsub,100)
      if(ish.ge.93)then
      write(ifch,*)('-',l=1,79)
      write(ifch,*)'entry jdeca. decay of',i,idptl(i),pptl(5,i)
      endif
      call idmass(111,amrho0)
      call idmass(221,amomeg)
      ioi=iorptl(i)
      if(.not.(ioi.gt.0.and.iabs(idptl(ioi)).lt.10000
     *.and.jorptl(i).eq.0))then
      if(idptl(i).eq.111)idptl(i)=221
      if(idptl(i).eq.221.and.rangen().gt.0.5)idptl(i)=111
      endif
      ida=iabs(idptl(i))
      if(mod(ndecay        ,10).eq.1)goto1000
      if(mod(ndecay/10     ,10).eq.1.and.ida.eq.  20)goto1000
      if(mod(ndecay/100    ,10).eq.1.and.ida.eq.2130)goto1000
      if(mod(ndecay/1000   ,10).eq.1.and.ida.eq.1130)goto1000
      if(mod(ndecay/1000   ,10).eq.1.and.ida.eq.2230)goto1000
      if(mod(ndecay/10000  ,10).eq.1.and.ida.eq.2330)goto1000
      if(mod(ndecay/10000  ,10).eq.1.and.ida.eq.1330)goto1000
      if(mod(ndecay/100000 ,10).eq.1.and.ida.eq.3331)goto1000
      if(mod(ndecay/1000000,10).eq.1.and.ida.eq. 110)goto1000
      if(mod(ndecax        ,10).eq.1.and.ida.eq. 441)goto1000
      if(mod(ndecax/10     ,10).eq.1.and.ida.eq. 230)goto1000
      if(mod(ndecax/100    ,10).eq.1.and.ida.eq.1111)goto1000
      if(mod(ndecax/100    ,10).eq.1.and.ida.eq.1121)goto1000
      if(mod(ndecax/100    ,10).eq.1.and.ida.eq.1221)goto1000
      if(mod(ndecax/100    ,10).eq.1.and.ida.eq.2221)goto1000
      if(mod(ndecax/1000   ,10).eq.1.and.ida.eq. 111)goto1000
      if(mod(ndecax/1000   ,10).eq.1.and.ida.eq. 121)goto1000
      if(mod(ndecax/1000   ,10).eq.1.and.ida.eq. 221)goto1000
      if(mod(ndecax/1000   ,10).eq.1.and.ida.eq. 331)goto1000
      if(mod(ndecax/10000  ,10).eq.1.and.ida.eq. 220)goto1000
      if(mod(ndecax/100000 ,10).eq.1.and.ida.eq. 330)goto1000
      if(mod(ndecax/1000000,10).eq.1.and.ida.eq. 112)goto1000
      if(mod(ndecax/1000000,10).eq.1.and.ida.eq. 122)goto1000
      if(mod(ndecaw        ,10).eq.1.and.ida.eq. 332)goto1000
      if(mod(ndecaw/10     ,10).eq.1.and.ida.eq. 131)goto1000
      if(mod(ndecaw/10     ,10).eq.1.and.ida.eq.-131)goto1000
      if(mod(ndecaw/10     ,10).eq.1.and.ida.eq. 231)goto1000
      if(mod(ndecaw/10     ,10).eq.1.and.ida.eq.-231)goto1000
      t=tivptl(2,i)
      nptlb=nptl
      if(nptl.gt.mxptl-10)
     *call utstop('jdeca: mxptl too small                  ')
      ish=ish0
      call jdecay(i,iret)
      if(ishsub/100.eq.5)ish=mod(ishsub,100)
      if(iret.eq.1)goto1000
      if(nptl.le.nptlb)goto1000
      ish00=ish
      if(ishsub/100.eq.14.and.kclude.eq.1)ish=mod(ishsub,100)
      if(ish.ge.91)write(ifch,101)sngl(ttaus)
101   format(1x,'decay  at tau =',e10.3)
      if(ish.ge.91)write(ifch,115)i,idptl(i)
     *,(pptl(k,i),k=3,5),(tivptl(k,i),k=1,2)
115   format(1x,'/cptl/',i4,i10
     *,3x,3(e8.2,1x),2x,2(e8.2,1x))
      istptl(i)=1
      ifrptl(1,i)=nptlb+1
      ifrptl(2,i)=nptl
      x=xorptl(1,i)+(t-xorptl(4,i))*pptl(1,i)/pptl(4,i)
      y=xorptl(2,i)+(t-xorptl(4,i))*pptl(2,i)/pptl(4,i)
      z=xorptl(3,i)+(t-xorptl(4,i))*pptl(3,i)/pptl(4,i)
      if(ish.ge.93)write(ifch,*)
     *'loop over decay products ',nptlb+1,' - ',nptl,' :'
      do 20 n=nptlb+1,nptl
      if(ish.ge.93)write(ifch,*)'particle: ',n,idptl(n)
      iorptl(n)=i
      jorptl(n)=0
      istptl(n)=0
      ifrptl(1,n)=0
      ifrptl(2,n)=0
      xorptl(1,n)=x
      xorptl(2,n)=y
      xorptl(3,n)=z
      xorptl(4,n)=t
      nqjptl(n)=nqjptl(i)
      io=n
1     io=iorptl(io)
      if(ish.ge.93)write(ifch,*)'io = ',io
      if(iorptl(io).gt.0)goto1
      if(ish.ge.93)write(ifch,*)'origin: ',io,idptl(io)
      zor=xorptl(3,io)
      tor=xorptl(4,io)
      call idquac(io,nq,ndu,ndu,jcdu)
      r=rangen()
      tauran=-taurea*alog(r)
      call uttaix(n,tauran,zor,tor,zis,tis)
      tivptl(1,n)=amax1(t,tis)
      call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
      tivptl(2,n)=t+taugm
      iclptl(n)=1
      if(ish.ge.91)write(ifch,115)n,idptl(n)
     *,(pptl(k,n),k=3,5),(tivptl(k,n),k=1,2)
20    continue
      ish=ish00
 
1000  continue
      if(ish.ge.93)write(ifch,*)('-',l=1,79)
      ish=ish0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jdecay(ip,iret)
c     decays particle ip from /cptl/
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/wco/wmass2,wgam2
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      common/cnsta/pi,prom,piom,ainfin
      dimension pgen(5,5),rnd(5),u(3),beta(3),idqk(3),rot(3,3),psave(3)
     1,reduce(5)
      dimension pold(4),pnew(4)
      dimension prest(4,5)
      common/cclude/kclude
      data reduce/1.,1.,2.,5.,15./
      data twome/1.022006e-3/
 
c          fctn definitions
      dot(i1,i2)=prest(4,i1)*prest(4,i2)-prest(1,i1)*prest(1,i2)
     *-prest(2,i1)*prest(2,i2)-prest(3,i1)*prest(3,i2)
c          charged w propagator.
      wprop(z)=(z-wmass2**2)**2+(wmass2*wgam2)**2
 
      ish0=ish
      if(ishsub/100.eq.4)ish=mod(ishsub,100)
      if(ish.ge.93)then
      write(ifch,*)('-',l=1,79)
      write(ifch,*)'decay of',ip,idptl(ip),pptl(5,ip)
      endif
      ish=ish0
 
      iret=0
      kclude=0
 
c     no k_long decay
c     ---------------
      if(idptl(ip).eq.-20)goto1000
 
c     cluster decay
c     -------------
      if(iabs(idptl(ip)).gt.1e8)then
      kclude=1
      call jclude(ip,iret)
      goto1000
      endif
 
      if(ishsub/100.eq.4)ish=mod(ishsub,100)
      if(ish.ge.93)write(ifch,*)'ordinary decay'
 
c     select decay mode
c     -----------------
      ntry=0
2     ntry=ntry+1
           if(ntry.gt.100)then
      if(ish.ge.90)then
      call utmsg('jdecay')
      write(ifch,*)'*****  decay not possible. iret = 1.'
      write(ifch,*)'id,mass: ',idptl(ip),pptl(5,ip)
      call utmsgf
      endif
      iret=1
      goto1000
           endif
      idlv1=idptl(ip)
      amss=pptl(5,ip)
1     continue
      ipoint=look(iabs(idlv1))-1
      if(ipoint.lt.0) goto1000
      try=rangen()
100   ipoint=ipoint+1
      if(ish.ge.93)write(ifch,*)'ipoint,cbr,try',ipoint,cbr(ipoint),try
      if(try.gt.cbr(ipoint)) goto100
      nadd=0
      sum=0.
      nstart=nptl+1
      do 110 i=1,5
      if(mode(i,ipoint).eq.0) goto110
      if(nptl+nadd+1.gt.mxptl) goto9999
      nadd=nadd+1
      new=nptl+nadd
      idptl(new)=mode(i,ipoint)
      idlv1=idptl(new)
      call idmass(idlv1,pptl(5,new))
      sum=sum+pptl(5,new)
110   continue
      if(nadd.ne.1.and.sum+1e-2.ge.amss)goto2
      nadd1=nadd-1
      do 120 j=1,5
      pgen(j,1)=pptl(j,ip)
120   continue
      pgen(5,nadd)=pptl(5,nptl+nadd)
      if(nadd.eq.1) goto700
      if(nadd.eq.2) goto400
 
c     use kroll-wada distribution for dalitz decays.
c     ----------------------------------------------
      if(ish.ge.93)write(ifch,*)'>= 3 body decay'
      if(.not.((idptl(ip).eq.110.or.idptl(ip).eq.220).and.
     1iabs(idptl(nptl+2)).eq.12)) goto130
      ntry=0
125   ntry=ntry+1
           if(ntry.gt.10)then
      if(ish.ge. 0)then
      call utmsg('jdecay')
      write(ifch,*)'*****  ntry > 10. iret = 1.'
      write(ifch,*)'amee,ree,wtee',amee,ree,wtee                        
      call utmsgf
      endif
      iret=1
      goto1000
           endif
      amee=twome*(pptl(5,ip)/twome)**rangen()
      ree=(twome/amee)**2
      wtee=(1.-(amee/pptl(5,ip))**2)**3*sqrt(1.-ree)*(1.+.5*ree)
      if(wtee.lt.rangen()) goto125
      pgen(5,2)=amee
      goto400
130   continue
 
c     calculate maximum phase-space weight
c     ------------------------------------
      wtmax=1./reduce(nadd)
      sum1=pgen(5,1)
      sum2=sum-pptl(5,nptl+1)
      do 200 i=1,nadd1
      wtmax=wtmax*utpcm(sum1,sum2,pptl(5,nptl+i))
      sum1=sum1-pptl(5,nptl+i)
      sum2=sum2-pptl(5,nptl+i+1)
200   continue
 
c     generate uniform nadd-body phase space
c     --------------------------------------
      ntry=0
300   ntry=ntry+1
           if(ntry.gt.10000)then
      if(ish.ge. 0)then
      call utmsg('jdecay')
      write(ifch,*)'*****  infinite loop (2). iret = 1.'
      write(ifch,*)'ip,idptl(ip),pptl(5,ip):',ip,idptl(ip),pptl(5,ip)
      write(ifch,*)'wt,wtmax:',wt,wtmax
      write(ifch,*)'i,pgen(5,i),pptl(5,nptl+i),idptl(nptl+i):'
      do i=1,nadd
      write(ifch,*)i,pgen(5,i),pptl(5,nptl+i),idptl(nptl+i)
      enddo
      call utmsgf
      endif
      iret=1
      goto1000
           endif
      rnd(1)=1.
      do 310 i=2,nadd1
      rnew=rangen()
      i1=i-1
      do 320 jj1=1,i1
      j=i-jj1
      jsave=j+1
      if(rnew.le.rnd(j)) goto310
      rnd(jsave)=rnd(j)
320   continue
310   rnd(jsave)=rnew
      rnd(nadd)=0.
      wt=1.
      sum1=sum
      do 330 i=2,nadd
      sum1=sum1-pptl(5,nptl+i-1)
      pgen(5,i)=sum1+rnd(i)*(pgen(5,1)-sum)
      a=pgen(5,i-1)
      b=pgen(5,i)
      c=pptl(5,nptl+i-1)
      wt=wt*utpcm(a,b,c)
330   continue
      if(wt.lt.rangen()*wtmax) goto300
 
c     carry out two-body decays in pgen frames
c     ----------------------------------------
400   continue
      if(ish.ge.93)write(ifch,*)'2 body decay'
      do 410 i=1,nadd1
      qcm=utpcm(pgen(5,i),pgen(5,i+1),pptl(5,nptl+i))
      u(3)=2.*rangen()-1.
      phi=2.*pi*rangen()
      u(1)=sqrt(1.-u(3)**2)*cos(phi)
      u(2)=sqrt(1.-u(3)**2)*sin(phi)
      do 420 j=1,3
      pptl(j,nptl+i)=qcm*u(j)
      pgen(j,i+1)=-pptl(j,nptl+i)
420   continue
      pptl(4,nptl+i)=sqrt(qcm**2+pptl(5,nptl+i)**2)
      pgen(4,i+1)=sqrt(qcm**2+pgen(5,i+1)**2)
410   continue
      do 430 j=1,4
      pptl(j,nptl+nadd)=pgen(j,nadd)
430   continue
 
c     boost pgen frames to lab frame
c          also save momenta in rest frame (last frame)
c     -------------------------------------------------
      do 500 ii=1,nadd1
      i=nadd-ii
      do 510 j=1,3
      beta(j)=pgen(j,i)/pgen(4,i)
510   continue
      gamma=pgen(4,i)/pgen(5,i)
      do 520 k=i,nadd
      k1=nptl+k
      bp=beta(1)*pptl(1,k1)+beta(2)*pptl(2,k1)+beta(3)*pptl(3,k1)
      do 530 j=1,3
      prest(j,k)=pptl(j,k1)
      pptl(j,k1)=pptl(j,k1)+gamma*beta(j)*(pptl(4,k1)
     1+bp*gamma/(gamma+1.))
530   continue
      prest(4,k)=pptl(4,k1)
      pptl(4,k1)=gamma*(pptl(4,k1)+bp)
520   continue
500   continue
 
c     matrix elements
c     ---------------
      if(nadd.eq.3.and.(idptl(ip).eq.221.or.idptl(ip).eq.331)) goto610
      if(nadd.eq.3.and.iabs(idptl(nptl+1)).lt.20.and.
     1idptl(nptl+1).ne.10) goto620
      goto800
 
c     omeg and phi decay
c          use vectors in rest frame
c     ------------------------------
610   wt=(pptl(5,nptl+1)*pptl(5,nptl+2)*pptl(5,nptl+3))**2
     1-(pptl(5,nptl+1)*dot(2,3))**2
     2-(pptl(5,nptl+2)*dot(1,3))**2
     3-(pptl(5,nptl+3)*dot(1,2))**2
     4+2.*dot(1,2)*dot(2,3)*dot(1,3)
      if(wt.lt.rangen()*pptl(5,ip)**6/108.) goto300
      goto800
 
c     semileptonic and quark decays
c          use vectors in rest frame, where ip has (m,0,0,0)
c          include w propagator
c     ------------------------------------------------------
620   wt=(pptl(5,ip)*prest(4,2))*dot(1,3)
      s12=pptl(5,nptl+1)**2+pptl(5,nptl+2)**2+2.*dot(1,2)
      s12max=pptl(5,ip)**2
      wt=wt*wprop(s12max)/wprop(s12)
      if(wt.lt.rangen()*pptl(5,ip)**4/16.) goto300
      goto800
 
c     one-particle decays
c     -------------------
700   continue
      do 710 j=1,5
      pptl(j,nptl+1)=pptl(j,ip)
710   continue
 
c     swap particles and antiparticles if idptl(ip)<0
c     -----------------------------------------------
800   continue
      if(idptl(ip).ge.0.or.iabs(idptl(ip)).eq.20) goto900
      do 810 i=1,nadd
      idabs=iabs(idptl(nptl+i))
      ifl1=idabs/1000
      ifl2=mod(idabs/100,10)
      ifl3=mod(idabs/10,10)
      if(ifl1.eq.0.and.ifl2.ne.0.and.ifl2.eq.ifl3) goto810
      if(idabs.eq.9.or.idabs.eq.10.or.idabs.eq.20) goto810
      if(idabs.eq.29.or.idabs.eq.30.or.idabs.eq.40) goto810
      idptl(nptl+i)=-idptl(nptl+i)
810   continue
 
900   continue
      nptl=nptl+nadd
      if(nptl.gt.mxptl)
     *call utstop('jdecay: nptl>mxptl                      ')
      nqk=0
      if(iabs(idptl(nptl)).lt.10.or.mod(idptl(nptl),100).eq.0)
     *call utstop('jdecay: decay ptcl is parton            ')
 
1000  continue
      if(ish.ge.93)write(ifch,*)('-',l=1,79)
      ish=ish0
      return
 
9999  call utstop('jdecay: mxptl too small                 ')
      end
 
c-----------------------------------------------------------------------
      subroutine jdecin(lprint)
c     sets up /dkytab/
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/wco/wmass2,wgam2
      dimension imode(6)
      character*8 idlabl,lmode(6),lres
      character*8 iquit,iblank,lread(10)
      logical lprint
      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      common/nodcay/nodcay,noeta,nopi0,nonunu,noevol,nohadr
      logical nodcay,noeta,nopi0,nonunu,noevol,nohadr
      parameter (ndectb=1171)
      real dectab(7,ndectb)
 
      if(ientro.eq.1)call jcentr(3,6,3,1)
      if(ientro.eq.2)call jcentd
      call jcentp
 
      call idresi
 
c     determine wmass2,wgam2
c     ----------------------
      alfa=1./137.036
      gf=1.16570e-5
      sin2w=.215
      sinw=sqrt(sin2w)
      cosw=sqrt(1.-sin2w)
      amw=sqrt(pi*alfa/(.9304*sqrt(2.)*gf))/sinw
      wmass2=amw
      call idmass(5,amlep5)
      call idmass(6,amlep6)
      ngam=12
      if(amlep5+amlep6.gt.amw) ngam=9
      wgam2=gf*amw**3/(6.*pi*sqrt(2.))*ngam
 
      data iquit/'////'/,iblank/' '/
      ird=0
      do 1 i=1,mxlook
1     look(i)=0
      do 2 i=1,mxdky
      do 3 j=1,5
3     mode(j,i)=0
2     cbr(i)=0.
      nodcay=.false.
      noeta=.false.
      nopi0=.false.
      nonunu=.false.
      noevol=.false.
      nohadr=.false.
      if(lprint) write(ifch,10)
10    format('1',30('*')/' *',28x,'*'/
     1' *',5x,'isajet decay table',5x,'*'/
     2' *',28x,'*'/' ',30('*')//
     36x,'part',18x,'decay mode',19x,'cum br',15x,'ident',17x,
     4'decay ident'/)
      loop=0
      iold=0
      if(nodcay) return
 
200   loop=loop+1
      if(loop.gt.mxdky) goto9999
220   do 210 i=1,5
      imode(i)=0
      lmode(i)=iblank
210   continue
      ird=ird+1
      if(ird.gt.ndectb)return
      ires=nint(dectab(1,ird))
      br=dectab(2,ird)
      do 215 i=1,5
215   imode(i)=nint(dectab(2+i,ird))
      if(nopi0.and.ires.eq.110) goto220
      if(noeta.and.ires.eq.220) goto220
      if(ires.eq.iold) goto230
      if(ires.lt.0.or.ires.gt.mxlook)
     *call utstop('jdecin: ires out of range               ')
      look(ires)=loop
230   iold=ires
      cbr(loop)=br
      do 240 i=1,5
      mode(i,loop)=imode(i)
      if(imode(i).ne.0) lmode(i)=idlabl(imode(i))
240   continue
      lres=idlabl(ires)
      if(lprint) write(ifch,20) lres,(lmode(k),k=1,5),
     1br,ires,(imode(k),k=1,5)
20    format(6x,a5,6x,5(a5,2x),3x,f8.5,15x,i5,4x,5(i5,2x))
      goto200
 
9999  write(ifch,*)'loop=', loop
      call utstop('jdecin: loop > mxdky                    ')
 
      data ((dectab(i,j),i=1,7),j=  1, 18)/
     *  110., .98850,  10.,  10.,   0.,   0.,   0.
     *, 110.,1.00000,  10.,  12., -12.,   0.,   0.
     *, 220., .38000,  10.,  10.,   0.,   0.,   0.
     *, 220., .71000, 110., 110., 110.,   0.,   0.
     *, 220., .94600, 120.,-120., 110.,   0.,   0.
     *, 220., .99500, 120.,-120.,  10.,   0.,   0.
     *, 220.,1.00000,  10.,  12., -12.,   0.,   0.
     *, 330., .44100, 220., 120.,-120.,   0.,   0.
     *, 330., .66100, 220., 110., 110.,   0.,   0.
     *, 330., .95900, 111.,  10.,   0.,   0.,   0.
     *, 330., .98000, 221.,  10.,   0.,   0.,   0.
     *, 330.,1.00000,  10.,  10.,   0.,   0.,   0.
     *, 121.,1.00000, 120., 110.,   0.,   0.,   0.
     *, 111., .99989, 120.,-120.,   0.,   0.,   0.
     *, 111., .99993,  12., -12.,   0.,   0.,   0.
     *, 111.,1.00000,  14., -14.,   0.,   0.,   0.
     *, 221., .89900, 120.,-120., 110.,   0.,   0.
     *, 221., .91200, 120.,-120.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j= 19, 36)/
     *  221., .99992, 110.,  10.,   0.,   0.,   0.
     *, 221.,1.00000,  12., -12.,   0.,   0.,   0.
     *, 331., .48600, 130.,-130.,   0.,   0.,   0.
     *, 331., .83700, 230.,-230.,   0.,   0.,   0.
     *, 331., .98400, 120.,-120., 110.,   0.,   0.
     *, 331., .99944, 220.,  10.,   0.,   0.,   0.
     *, 331., .99975,  12., -12.,   0.,   0.,   0.
     *, 331.,1.00000,  14., -14.,   0.,   0.,   0.
     *, 230., .50000,  20.,   0.,   0.,   0.,   0.
     *, 230.,1.00000, -20.,   0.,   0.,   0.,   0.
     *, 131., .66670, 230., 120.,   0.,   0.,   0.
     *, 131.,1.00000, 130., 110.,   0.,   0.,   0.
     *, 231., .66670, 130.,-120.,   0.,   0.,   0.
     *, 231.,1.00000, 230., 110.,   0.,   0.,   0.
     *, 240., .11000,  12., -11., 230.,   0.,   0.
     *, 240., .17000,  12., -11., 231.,   0.,   0.
     *, 240., .28000,  14., -13., 230.,   0.,   0.
     *, 240., .34000,  14., -13., 231.,   0.,   0./
      data ((dectab(i,j),i=1,7),j= 37, 54)/
     *  240., .37800, 230.,-120.,   0.,   0.,   0.
     *, 240., .56300, 230.,-121.,   0.,   0.,   0.
     *, 240., .60800, 231.,-120.,   0.,   0.,   0.
     *, 240., .62100, 230.,-120., 110.,   0.,   0.
     *, 240., .71000, 130.,-120.,-120.,   0.,   0.
     *, 240., .80100, 230.,-120.,-120., 120.,   0.
     *, 240., .87900, 130.,-120.,-120., 110.,   0.
     *, 240., .95400, 230.,-120., 110., 110.,   0.
     *, 240., .96600, 230.,-130.,   0.,   0.,   0.
     *, 240., .97600, 331.,-120.,   0.,   0.,   0.
     *, 240., .98800,-130., 231.,   0.,   0.,   0.
     *, 240.,1.00000,-131., 230.,   0.,   0.,   0.
     *, 140., .04500,  12., -11., 130.,   0.,   0.
     *, 140., .07500,  12., -11., 131.,   0.,   0.
     *, 140., .12000,  14., -13., 130.,   0.,   0.
     *, 140., .15000,  14., -13., 131.,   0.,   0.
     *, 140., .20300, 130.,-120.,   0.,   0.,   0.
     *, 140., .22700, 230., 110.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j= 55, 72)/
     *  140., .24700, 230., 220.,   0.,   0.,   0.
     *, 140., .28900, 230., 221.,   0.,   0.,   0.
     *, 140., .45100, 130.,-121.,   0.,   0.,   0.
     *, 140., .53600, 131.,-120.,   0.,   0.,   0.
     *, 140., .56200, 231., 110.,   0.,   0.,   0.
     *, 140., .57600, 230., 111.,   0.,   0.,   0.
     *, 140., .58700, 130.,-120., 110.,   0.,   0.
     *, 140., .60300, 230.,-120., 120.,   0.,   0.
     *, 140., .72700, 130.,-120.,-120., 120.,   0.
     *, 140., .87600, 230.,-120., 120., 110.,   0.
     *, 140., .96900, 130.,-120., 110., 110.,   0.
     *, 140.,1.00000, 230., 110., 110., 110.,   0.
     *, 340., .03250,  12., -11., 220.,   0.,   0.
     *, 340., .06500,  12., -11., 331.,   0.,   0.
     *, 340., .09750,  14., -13., 220.,   0.,   0.
     *, 340., .13000,  14., -13., 331.,   0.,   0.
     *, 340., .17900,-130., 230.,   0.,   0.,   0.
     *, 340., .22800,-120., 220.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j= 73, 90)/
     *  340., .33800,-131., 230.,   0.,   0.,   0.
     *, 340., .44800,-130., 231.,   0.,   0.,   0.
     *, 340., .55800,-120., 331.,   0.,   0.,   0.
     *, 340., .57500,-130., 230., 110.,   0.,   0.
     *, 340., .59200,-230., 230.,-120.,   0.,   0.
     *, 340., .69400,-130., 230.,-120., 120.,   0.
     *, 340., .79600,-130., 230., 110., 110.,   0.
     *, 340., .89800,-130., 130.,-120., 110.,   0.
     *, 340.,1.00000,-230., 230.,-120., 110.,   0.
     *, 241., .64000, 140.,-120.,   0.,   0.,   0.
     *, 241., .92000, 240., 110.,   0.,   0.,   0.
     *, 241.,1.00000, 240.,  10.,   0.,   0.,   0.
     *, 141., .55000, 140., 110.,   0.,   0.,   0.
     *, 141.,1.00000, 140.,  10.,   0.,   0.,   0.
     *, 341.,1.00000, 340.,  10.,   0.,   0.,   0.
     *, 441., .07400,  12., -12.,   0.,   0.,   0.
     *, 441., .14800,  14., -14.,   0.,   0.,   0.
     *, 441., .15210,-121., 120.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j= 91,108)/
     *  441., .15620, 111., 110.,   0.,   0.,   0.
     *, 441., .16020, 121.,-120.,   0.,   0.,   0.
     *, 441., .16300,-121., 111., 120.,   0.,   0.
     *, 441., .16580, 121.,-121., 110.,   0.,   0.
     *, 441., .16860, 121., 111.,-120.,   0.,   0.
     *, 441., .28740, 120.,-120., 130.,-130.,   0.
     *, 441., .40620, 110., 110., 130.,-130.,   0.
     *, 441., .52500, 120.,-120., 120.,-120.,   0.
     *, 441., .64380, 120.,-120., 110., 110.,   0.
     *, 441., .76260, 110., 110., 110., 110.,   0.
     *, 441., .88130, 120.,-120., 230.,-230.,   0.
     *, 441.,1.00000, 110., 110., 230., 230.,   0.
     *, 150., .06000, -12.,  11., 140.,   0.,   0.
     *, 150., .12000, -12.,  11., 141.,   0.,   0.
     *, 150., .18000, -14.,  13., 140.,   0.,   0.
     *, 150., .24000, -14.,  13., 141.,   0.,   0.
     *, 150., .25500, -16.,  15., 140.,   0.,   0.
     *, 150., .27000, -16.,  15., 141.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=109,126)/
     *  150., .28050, 140., 120.,   0.,   0.,   0.
     *, 150., .29100, 140., 121.,   0.,   0.,   0.
     *, 150., .30150, 141., 120.,   0.,   0.,   0.
     *, 150., .31200, 141., 121.,   0.,   0.,   0.
     *, 150., .32650, 140.,-340.,   0.,   0.,   0.
     *, 150., .34100, 140.,-341.,   0.,   0.,   0.
     *, 150., .35550, 141.,-340.,   0.,   0.,   0.
     *, 150., .37000, 141.,-341.,   0.,   0.,   0.
     *, 150., .39800, 140., 120., 110.,   0.,   0.
     *, 150., .42600, 140., 120., 220.,   0.,   0.
     *, 150., .45400, 140., 120., 111.,   0.,   0.
     *, 150., .48200, 140., 120., 221.,   0.,   0.
     *, 150., .51000, 140., 121., 110.,   0.,   0.
     *, 150., .53800, 140., 121., 220.,   0.,   0.
     *, 150., .56600, 140., 121., 111.,   0.,   0.
     *, 150., .59400, 140., 121., 221.,   0.,   0.
     *, 150., .62200, 141., 120., 110.,   0.,   0.
     *, 150., .65000, 141., 120., 220.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=127,144)/
     *  150., .67800, 141., 120., 111.,   0.,   0.
     *, 150., .70600, 141., 120., 221.,   0.,   0.
     *, 150., .73400, 141., 121., 110.,   0.,   0.
     *, 150., .76200, 141., 121., 220.,   0.,   0.
     *, 150., .79000, 141., 121., 111.,   0.,   0.
     *, 150., .81800, 141., 121., 221.,   0.,   0.
     *, 150., .83200, 140., 130.,-230.,   0.,   0.
     *, 150., .84600, 140., 130.,-231.,   0.,   0.
     *, 150., .86000, 140., 131.,-230.,   0.,   0.
     *, 150., .87400, 140., 131.,-231.,   0.,   0.
     *, 150., .88800, 141., 130.,-230.,   0.,   0.
     *, 150., .90200, 141., 130.,-231.,   0.,   0.
     *, 150., .91600, 141., 131.,-230.,   0.,   0.
     *, 150., .93000, 141., 131.,-231.,   0.,   0.
     *, 150., .93300, 140.,-140., 130.,   0.,   0.
     *, 150., .93600, 140.,-140., 131.,   0.,   0.
     *, 150., .93900, 140.,-141., 130.,   0.,   0.
     *, 150., .94200, 140.,-141., 131.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=145,162)/
     *  150., .94500, 141.,-140., 130.,   0.,   0.
     *, 150., .94800, 141.,-140., 131.,   0.,   0.
     *, 150., .95100, 141.,-141., 130.,   0.,   0.
     *, 150., .95400, 141.,-141., 131.,   0.,   0.
     *, 150., .95700, 140.,-240., 230.,   0.,   0.
     *, 150., .96000, 140.,-240., 231.,   0.,   0.
     *, 150., .96300, 140.,-241., 230.,   0.,   0.
     *, 150., .96600, 140.,-241., 231.,   0.,   0.
     *, 150., .96900, 141.,-240., 230.,   0.,   0.
     *, 150., .97200, 141.,-240., 231.,   0.,   0.
     *, 150., .97500, 141.,-241., 230.,   0.,   0.
     *, 150., .97800, 141.,-241., 231.,   0.,   0.
     *, 150., .97950, 140.,-340., 330.,   0.,   0.
     *, 150., .98100, 140.,-340., 331.,   0.,   0.
     *, 150., .98250, 140.,-341., 331.,   0.,   0.
     *, 150., .98400, 140.,-341., 331.,   0.,   0.
     *, 150., .98550, 141.,-340., 330.,   0.,   0.
     *, 150., .98700, 141.,-340., 331.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=163,180)/
     *  150., .98850, 141.,-341., 331.,   0.,   0.
     *, 150., .99000, 141.,-341., 331.,   0.,   0.
     *, 150., .99200, 441., 130., 110.,   0.,   0.
     *, 150., .99400, 441., 131., 110.,   0.,   0.
     *, 150., .99600, 441., 230., 120.,   0.,   0.
     *, 150., .99800, 441., 231., 120.,   0.,   0.
     *, 150., .99900, 441., 330., 130.,   0.,   0.
     *, 150.,1.00000, 441., 331., 130.,   0.,   0.
     *, 250., .06000, -12.,  11., 240.,   0.,   0.
     *, 250., .12000, -12.,  11., 241.,   0.,   0.
     *, 250., .18000, -14.,  13., 240.,   0.,   0.
     *, 250., .24000, -14.,  13., 241.,   0.,   0.
     *, 250., .25500, -16.,  15., 240.,   0.,   0.
     *, 250., .27000, -16.,  15., 241.,   0.,   0.
     *, 250., .28050, 240., 120.,   0.,   0.,   0.
     *, 250., .29100, 240., 121.,   0.,   0.,   0.
     *, 250., .30150, 241., 120.,   0.,   0.,   0.
     *, 250., .31200, 241., 121.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=181,198)/
     *  250., .32650, 240.,-340.,   0.,   0.,   0.
     *, 250., .34100, 240.,-341.,   0.,   0.,   0.
     *, 250., .35550, 241.,-340.,   0.,   0.,   0.
     *, 250., .37000, 241.,-341.,   0.,   0.,   0.
     *, 250., .39800, 240., 120., 110.,   0.,   0.
     *, 250., .42600, 240., 120., 220.,   0.,   0.
     *, 250., .45400, 240., 120., 111.,   0.,   0.
     *, 250., .48200, 240., 120., 221.,   0.,   0.
     *, 250., .51000, 240., 121., 110.,   0.,   0.
     *, 250., .53800, 240., 121., 220.,   0.,   0.
     *, 250., .56600, 240., 121., 111.,   0.,   0.
     *, 250., .59400, 240., 121., 221.,   0.,   0.
     *, 250., .62200, 241., 120., 110.,   0.,   0.
     *, 250., .65000, 241., 120., 220.,   0.,   0.
     *, 250., .67800, 241., 120., 111.,   0.,   0.
     *, 250., .70600, 241., 120., 221.,   0.,   0.
     *, 250., .73400, 241., 121., 110.,   0.,   0.
     *, 250., .76200, 241., 121., 220.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=199,216)/
     *  250., .79000, 241., 121., 111.,   0.,   0.
     *, 250., .81800, 241., 121., 221.,   0.,   0.
     *, 250., .83200, 240., 130.,-230.,   0.,   0.
     *, 250., .84600, 240., 130.,-231.,   0.,   0.
     *, 250., .86000, 240., 131.,-230.,   0.,   0.
     *, 250., .87400, 240., 131.,-231.,   0.,   0.
     *, 250., .88800, 241., 130.,-230.,   0.,   0.
     *, 250., .90200, 241., 130.,-231.,   0.,   0.
     *, 250., .91600, 241., 131.,-230.,   0.,   0.
     *, 250., .93000, 241., 131.,-231.,   0.,   0.
     *, 250., .93300, 240.,-140., 130.,   0.,   0.
     *, 250., .93600, 240.,-140., 131.,   0.,   0.
     *, 250., .93900, 240.,-141., 130.,   0.,   0.
     *, 250., .94200, 240.,-141., 131.,   0.,   0.
     *, 250., .94500, 241.,-140., 130.,   0.,   0.
     *, 250., .94800, 241.,-140., 131.,   0.,   0.
     *, 250., .95100, 241.,-141., 130.,   0.,   0.
     *, 250., .95400, 241.,-141., 131.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=217,234)/
     *  250., .95700, 240.,-240., 230.,   0.,   0.
     *, 250., .96000, 240.,-240., 231.,   0.,   0.
     *, 250., .96300, 240.,-241., 230.,   0.,   0.
     *, 250., .96600, 240.,-241., 231.,   0.,   0.
     *, 250., .96900, 241.,-240., 230.,   0.,   0.
     *, 250., .97200, 241.,-240., 231.,   0.,   0.
     *, 250., .97500, 241.,-241., 230.,   0.,   0.
     *, 250., .97800, 241.,-241., 231.,   0.,   0.
     *, 250., .97950, 240.,-340., 330.,   0.,   0.
     *, 250., .98100, 240.,-340., 331.,   0.,   0.
     *, 250., .98250, 240.,-341., 331.,   0.,   0.
     *, 250., .98400, 240.,-341., 331.,   0.,   0.
     *, 250., .98550, 241.,-340., 330.,   0.,   0.
     *, 250., .98700, 241.,-340., 331.,   0.,   0.
     *, 250., .98850, 241.,-341., 331.,   0.,   0.
     *, 250., .99000, 241.,-341., 331.,   0.,   0.
     *, 250., .99200, 441., 130.,-120.,   0.,   0.
     *, 250., .99400, 441., 131.,-120.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=235,252)/
     *  250., .99600, 441., 230., 220.,   0.,   0.
     *, 250., .99800, 441., 231., 221.,   0.,   0.
     *, 250., .99900, 441., 330., 230.,   0.,   0.
     *, 250.,1.00000, 441., 331., 230.,   0.,   0.
     *, 350., .06000, -12.,  11., 340.,   0.,   0.
     *, 350., .12000, -12.,  11., 341.,   0.,   0.
     *, 350., .18000, -14.,  13., 340.,   0.,   0.
     *, 350., .24000, -14.,  13., 341.,   0.,   0.
     *, 350., .25500, -16.,  15., 340.,   0.,   0.
     *, 350., .27000, -16.,  15., 341.,   0.,   0.
     *, 350., .28050, 340., 120.,   0.,   0.,   0.
     *, 350., .29100, 340., 121.,   0.,   0.,   0.
     *, 350., .30150, 341., 120.,   0.,   0.,   0.
     *, 350., .31200, 341., 121.,   0.,   0.,   0.
     *, 350., .32650, 340.,-340.,   0.,   0.,   0.
     *, 350., .34100, 340.,-341.,   0.,   0.,   0.
     *, 350., .35550, 341.,-340.,   0.,   0.,   0.
     *, 350., .37000, 341.,-341.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=253,270)/
     *  350., .39800, 340., 120., 110.,   0.,   0.
     *, 350., .42600, 340., 120., 220.,   0.,   0.
     *, 350., .45400, 340., 120., 111.,   0.,   0.
     *, 350., .48200, 340., 120., 221.,   0.,   0.
     *, 350., .51000, 340., 121., 110.,   0.,   0.
     *, 350., .53800, 340., 121., 220.,   0.,   0.
     *, 350., .56600, 340., 121., 111.,   0.,   0.
     *, 350., .59400, 340., 121., 221.,   0.,   0.
     *, 350., .62200, 341., 120., 110.,   0.,   0.
     *, 350., .65000, 341., 120., 220.,   0.,   0.
     *, 350., .67800, 341., 120., 111.,   0.,   0.
     *, 350., .70600, 341., 120., 221.,   0.,   0.
     *, 350., .73400, 341., 121., 110.,   0.,   0.
     *, 350., .76200, 341., 121., 220.,   0.,   0.
     *, 350., .79000, 341., 121., 111.,   0.,   0.
     *, 350., .81800, 341., 121., 221.,   0.,   0.
     *, 350., .83200, 340., 130.,-230.,   0.,   0.
     *, 350., .84600, 340., 130.,-231.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=271,288)/
     *  350., .86000, 340., 131.,-230.,   0.,   0.
     *, 350., .87400, 340., 131.,-231.,   0.,   0.
     *, 350., .88800, 341., 130.,-230.,   0.,   0.
     *, 350., .90200, 341., 130.,-231.,   0.,   0.
     *, 350., .91600, 341., 131.,-230.,   0.,   0.
     *, 350., .93000, 341., 131.,-231.,   0.,   0.
     *, 350., .93300, 340.,-140., 130.,   0.,   0.
     *, 350., .93600, 340.,-140., 131.,   0.,   0.
     *, 350., .93900, 340.,-141., 130.,   0.,   0.
     *, 350., .94200, 340.,-141., 131.,   0.,   0.
     *, 350., .94500, 341.,-140., 130.,   0.,   0.
     *, 350., .94800, 341.,-140., 131.,   0.,   0.
     *, 350., .95100, 341.,-141., 130.,   0.,   0.
     *, 350., .95400, 341.,-141., 131.,   0.,   0.
     *, 350., .95700, 340.,-240., 230.,   0.,   0.
     *, 350., .96000, 340.,-240., 231.,   0.,   0.
     *, 350., .96300, 340.,-241., 230.,   0.,   0.
     *, 350., .96600, 340.,-241., 231.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=289,306)/
     *  350., .96900, 341.,-240., 230.,   0.,   0.
     *, 350., .97200, 341.,-240., 231.,   0.,   0.
     *, 350., .97500, 341.,-241., 230.,   0.,   0.
     *, 350., .97800, 341.,-241., 231.,   0.,   0.
     *, 350., .97950, 340.,-340., 330.,   0.,   0.
     *, 350., .98100, 340.,-340., 331.,   0.,   0.
     *, 350., .98250, 340.,-341., 331.,   0.,   0.
     *, 350., .98400, 340.,-341., 331.,   0.,   0.
     *, 350., .98550, 341.,-340., 330.,   0.,   0.
     *, 350., .98700, 341.,-340., 331.,   0.,   0.
     *, 350., .98850, 341.,-341., 331.,   0.,   0.
     *, 350., .99000, 341.,-341., 331.,   0.,   0.
     *, 350., .99200, 441., 130.,-130.,   0.,   0.
     *, 350., .99400, 441., 131.,-130.,   0.,   0.
     *, 350., .99600, 441., 230.,-230.,   0.,   0.
     *, 350., .99800, 441., 231.,-230.,   0.,   0.
     *, 350., .99900, 441., 330., 330.,   0.,   0.
     *, 350.,1.00000, 441., 331., 331.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=307,324)/
     *  160., .33330,  -1.,   2.,  -5.,   0.,   0.
     *, 160., .66660,  -4.,   3.,  -5.,   0.,   0.
     *, 160., .77770, -11.,  12.,  -5.,   0.,   0.
     *, 160., .88880, -13.,  14.,  -5.,   0.,   0.
     *, 160.,1.00000, -15.,  16.,  -5.,   0.,   0.
     *, 260., .33330,  -1.,   2.,  -5.,   0.,   0.
     *, 260., .66660,  -4.,   3.,  -5.,   0.,   0.
     *, 260., .77770, -11.,  12.,  -5.,   0.,   0.
     *, 260., .88880, -13.,  14.,  -5.,   0.,   0.
     *, 260.,1.00000, -15.,  16.,  -5.,   0.,   0.
     *, 360., .33330,  -1.,   2.,  -5.,   0.,   0.
     *, 360., .66660,  -4.,   3.,  -5.,   0.,   0.
     *, 360., .77770, -11.,  12.,  -5.,   0.,   0.
     *, 360., .88880, -13.,  14.,  -5.,   0.,   0.
     *, 360.,1.00000, -15.,  16.,  -5.,   0.,   0.
     *, 151.,1.00000, 150.,  10.,   0.,   0.,   0.
     *, 251.,1.00000, 250.,  10.,   0.,   0.,   0.
     *, 351.,1.00000, 350.,  10.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=325,342)/
     *  161.,1.00000, 160.,  10.,   0.,   0.,   0.
     *, 261.,1.00000, 260.,  10.,   0.,   0.,   0.
     *, 361.,1.00000, 360.,  10.,   0.,   0.,   0.
     *,1230.,1.00000,2130.,  10.,   0.,   0.,   0.
     *,1111.,1.00000,1120., 120.,   0.,   0.,   0.
     *,1121., .66670,1120., 110.,   0.,   0.,   0.
     *,1121.,1.00000,1220., 120.,   0.,   0.,   0.
     *,1221., .66670,1220., 110.,   0.,   0.,   0.
     *,1221.,1.00000,1120.,-120.,   0.,   0.,   0.
     *,2221.,1.00000,1220.,-120.,   0.,   0.,   0.
     *,1131., .88000,2130., 120.,   0.,   0.,   0.
     *,1131., .94000,1130., 110.,   0.,   0.,   0.
     *,1131.,1.00000,1230., 120.,   0.,   0.,   0.
     *,1231., .88000,2130., 110.,   0.,   0.,   0.
     *,1231., .94000,1130.,-120.,   0.,   0.,   0.
     *,1231.,1.00000,2230., 120.,   0.,   0.,   0.
     *,2231., .88000,2130.,-120.,   0.,   0.,   0.
     *,2231., .94000,1230.,-120.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=343,360)/
     * 2231.,1.00000,2230., 110.,   0.,   0.,   0.
     *,1331., .66670,2330., 120.,   0.,   0.,   0.
     *,1331.,1.00000,1330., 110.,   0.,   0.,   0.
     *,2331., .66670,1330.,-120.,   0.,   0.,   0.
     *,2331.,1.00000,2330., 110.,   0.,   0.,   0.
     *,  16., .18000,  12., -11.,  15.,   0.,   0.
     *,  16., .36000,  14., -13.,  15.,   0.,   0.
     *,  16., .45100,-120.,  15.,   0.,   0.,   0.
     *,  16., .66000,-121.,  15.,   0.,   0.,   0.
     *,  16., .78000, 110., 110.,-120.,  15.,   0.
     *,  16., .83600, 120.,-120.,-120.,  15.,   0.
     *,  16.,1.00000, 120., 110.,-120.,-120.,  15.
     *,2140., .03750, -12.,  11.,2130.,   0.,   0.
     *,2140., .07500, -12.,  11.,1231.,   0.,   0.
     *,2140., .11250, -14.,  13.,2130.,   0.,   0.
     *,2140., .15000, -14.,  13.,1231.,   0.,   0.
     *,2140., .18200,2130., 120.,   0.,   0.,   0.
     *,2140., .21300,1230., 110.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=361,378)/
     * 2140., .24400,1120.,-230.,   0.,   0.,   0.
     *,2140., .29500,1131., 110.,   0.,   0.,   0.
     *,2140., .34600,1231., 120.,   0.,   0.,   0.
     *,2140., .39700,1121.,-230.,   0.,   0.,   0.
     *,2140., .44800,1111.,-130.,   0.,   0.,   0.
     *,2140., .49900,1130., 111.,   0.,   0.,   0.
     *,2140., .55000,1230., 121.,   0.,   0.,   0.
     *,2140., .60100,1120.,-231.,   0.,   0.,   0.
     *,2140., .65800,1120.,-230., 120.,-120.,   0.
     *,2140., .71500,1120.,-230., 110., 110.,   0.
     *,2140., .77200,1120.,-130., 120., 110.,   0.
     *,2140., .82900,1220.,-230., 120., 110.,   0.
     *,2140., .88600,1220.,-130., 120., 120.,   0.
     *,2140., .94300,2130., 120., 120.,-120.,   0.
     *,2140.,1.00000,2130., 120., 110., 110.,   0.
     *,1140.,1.00000,2140., 120.,   0.,   0.,   0.
     *,1240.,1.00000,2140., 110.,   0.,   0.,   0.
     *,2240.,1.00000,2140.,-120.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=379,396)/
     * 1340., .03750, -12.,  11.,1330.,   0.,   0.
     *,1340., .07500, -12.,  11.,1331.,   0.,   0.
     *,1340., .11250, -14.,  13.,1330.,   0.,   0.
     *,1340., .15000, -14.,  13.,1331.,   0.,   0.
     *,1340., .19900,1330., 120.,   0.,   0.,   0.
     *,1340., .24800,1231.,-230.,   0.,   0.,   0.
     *,1340., .28800,1330., 120.,   0.,   0.,   0.
     *,1340., .32800,1131.,-230.,   0.,   0.,   0.
     *,1340., .36800,1330., 121.,   0.,   0.,   0.
     *,1340., .40800,1130.,-230.,   0.,   0.,   0.
     *,1340., .44800,1330., 120., 110.,   0.,   0.
     *,1340., .48800,2330., 120., 120.,   0.,   0.
     *,1340., .52800,1130.,-130., 120.,   0.,   0.
     *,1340., .56800,1130.,-230., 110.,   0.,   0.
     *,1340., .60800,1230.,-230., 120.,   0.,   0.
     *,1340., .66400,2130.,-230., 120., 110.,   0.
     *,1340., .72000,2130.,-130., 120., 120.,   0.
     *,1340., .77600,1130.,-230., 120., 120.,   0./
      data ((dectab(i,j),i=1,7),j=397,414)/
     * 1340., .83200,1130.,-230., 110., 110.,   0.
     *,1340., .88800,1330., 120., 120.,-120.,   0.
     *,1340., .94400,1330., 120., 110., 110.,   0.
     *,1340.,1.00000,2330., 120., 120., 110.,   0.
     *,3140., .03750, -12.,  11.,1330.,   0.,   0.
     *,3140., .07500, -12.,  11.,1331.,   0.,   0.
     *,3140., .11250, -14.,  13.,1330.,   0.,   0.
     *,3140., .15000, -14.,  13.,1331.,   0.,   0.
     *,3140., .19900,1330., 120.,   0.,   0.,   0.
     *,3140., .24800,1231.,-230.,   0.,   0.,   0.
     *,3140., .28800,1330., 120.,   0.,   0.,   0.
     *,3140., .32800,1131.,-230.,   0.,   0.,   0.
     *,3140., .36800,1330., 121.,   0.,   0.,   0.
     *,3140., .40800,1130.,-230.,   0.,   0.,   0.
     *,3140., .44800,1330., 120., 110.,   0.,   0.
     *,3140., .48800,2330., 120., 120.,   0.,   0.
     *,3140., .52800,1130.,-130., 120.,   0.,   0.
     *,3140., .56800,1130.,-230., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=415,432)/
     * 3140., .60800,1230.,-230., 120.,   0.,   0.
     *,3140., .66400,2130.,-230., 120., 110.,   0.
     *,3140., .72000,2130.,-130., 120., 120.,   0.
     *,3140., .77600,1130.,-230., 120., 120.,   0.
     *,3140., .83200,1130.,-230., 110., 110.,   0.
     *,3140., .88800,1330., 120., 120.,-120.,   0.
     *,3140., .94400,1330., 120., 110., 110.,   0.
     *,3140.,1.00000,2330., 120., 120., 110.,   0.
     *,2340., .03750, -12.,  11.,2330.,   0.,   0.
     *,2340., .07500, -12.,  11.,2331.,   0.,   0.
     *,2340., .11250, -14.,  13.,2330.,   0.,   0.
     *,2340., .15000, -14.,  13.,2331.,   0.,   0.
     *,2340., .17500,2330., 120.,   0.,   0.,   0.
     *,2340., .20000,1330., 110.,   0.,   0.,   0.
     *,2340., .22500,1130.,-130.,   0.,   0.,   0.
     *,2340., .25000,1230.,-230.,   0.,   0.,   0.
     *,2340., .29500,2331., 120.,   0.,   0.,   0.
     *,2340., .34000,1331., 110.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=433,450)/
     * 2340., .38500,1131.,-130.,   0.,   0.,   0.
     *,2340., .43000,1231.,-230.,   0.,   0.,   0.
     *,2340., .47500,2330., 121.,   0.,   0.,   0.
     *,2340., .52000,1330., 111.,   0.,   0.,   0.
     *,2340., .56500,1130.,-131.,   0.,   0.,   0.
     *,2340., .61000,1230.,-231.,   0.,   0.,   0.
     *,2340., .64900,2130.,-230., 120.,-120.,   0.
     *,2340., .68800,2130.,-230., 110., 110.,   0.
     *,2340., .72700,2130.,-130., 120., 110.,   0.
     *,2340., .76600,1130.,-230.,-120., 110.,   0.
     *,2340., .80500,1130.,-130., 120.,-120.,   0.
     *,2340., .84400,1130.,-130., 110., 110.,   0.
     *,2340., .88300,1330., 120.,-120., 110.,   0.
     *,2340., .92200,1330., 110., 110., 110.,   0.
     *,2340., .96100,2330., 120., 120.,-120.,   0.
     *,2340.,1.00000,2330., 120., 110., 110.,   0.
     *,3240., .03750, -12.,  11.,2330.,   0.,   0.
     *,3240., .07500, -12.,  11.,2331.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=451,468)/
     * 3240., .11250, -14.,  13.,2330.,   0.,   0.
     *,3240., .15000, -14.,  13.,2331.,   0.,   0.
     *,3240., .17500,2330., 120.,   0.,   0.,   0.
     *,3240., .20000,1330., 110.,   0.,   0.,   0.
     *,3240., .22500,1130.,-130.,   0.,   0.,   0.
     *,3240., .25000,1230.,-230.,   0.,   0.,   0.
     *,3240., .29500,2331., 120.,   0.,   0.,   0.
     *,3240., .34000,1331., 110.,   0.,   0.,   0.
     *,3240., .38500,1131.,-130.,   0.,   0.,   0.
     *,3240., .43000,1231.,-230.,   0.,   0.,   0.
     *,3240., .47500,2330., 121.,   0.,   0.,   0.
     *,3240., .52000,1330., 111.,   0.,   0.,   0.
     *,3240., .56500,1130.,-131.,   0.,   0.,   0.
     *,3240., .61000,1230.,-231.,   0.,   0.,   0.
     *,3240., .64900,2130.,-230., 120.,-120.,   0.
     *,3240., .68800,2130.,-230., 110., 110.,   0.
     *,3240., .72700,2130.,-130., 120., 110.,   0.
     *,3240., .76600,1130.,-230.,-120., 110.,   0./
      data ((dectab(i,j),i=1,7),j=469,486)/
     * 3240., .80500,1130.,-130., 120.,-120.,   0.
     *,3240., .84400,1130.,-130., 110., 110.,   0.
     *,3240., .88300,1330., 120.,-120., 110.,   0.
     *,3240., .92200,1330., 110., 110., 110.,   0.
     *,3240., .96100,2330., 120., 120.,-120.,   0.
     *,3240.,1.00000,2330., 120., 110., 110.,   0.
     *,3340., .07500, -12.,  11.,3331.,   0.,   0.
     *,3340., .15000, -12.,  11.,3331.,   0.,   0.
     *,3340., .25000,1330.,-230.,   0.,   0.,   0.
     *,3340., .31000,3331., 120.,   0.,   0.,   0.
     *,3340., .37000,1331.,-230.,   0.,   0.,   0.
     *,3340., .43000,1330.,-231.,   0.,   0.,   0.
     *,3340., .49000,2330.,-230., 120.,   0.,   0.
     *,3340., .55000,1330.,-230., 110.,   0.,   0.
     *,3340., .61000,1330.,-130., 120.,   0.,   0.
     *,3340., .67500,3331., 120., 120.,-120.,   0.
     *,3340., .74000,3331., 120., 110., 110.,   0.
     *,3340., .80500,1330.,-230., 120.,-120.,   0./
      data ((dectab(i,j),i=1,7),j=487,504)/
     * 3340., .87000,1330.,-230., 110., 110.,   0.
     *,3340., .93500,2330.,-230., 120., 110.,   0.
     *,3340.,1.00000,2330.,-130., 120., 120.,   0.
     *,1141.,1.00000,2140., 120.,   0.,   0.,   0.
     *,1241.,1.00000,2140., 110.,   0.,   0.,   0.
     *,2241.,1.00000,2140.,-120.,   0.,   0.,   0.
     *,1341., .66670,2340., 120.,   0.,   0.,   0.
     *,1341.,1.00000,1340., 110.,   0.,   0.,   0.
     *,2341., .66670,1340.,-120.,   0.,   0.,   0.
     *,2341.,1.00000,2340., 110.,   0.,   0.,   0.
     *,3341.,1.00000,3340., 110.,   0.,   0.,   0.
     *,1150., .06000,  12., -11.,1140.,   0.,   0.
     *,1150., .12000,  12., -11.,1141.,   0.,   0.
     *,1150., .18000,  14., -13.,1140.,   0.,   0.
     *,1150., .24000,  14., -13.,1141.,   0.,   0.
     *,1150., .25500,  16., -15.,1140.,   0.,   0.
     *,1150., .27000,  16., -15.,1141.,   0.,   0.
     *,1150., .28925,1140.,-120.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=505,522)/
     * 1150., .30850,1140.,-121.,   0.,   0.,   0.
     *,1150., .32775,1141.,-120.,   0.,   0.,   0.
     *,1150., .34700,1141.,-121.,   0.,   0.,   0.
     *,1150., .35775,1140., 340.,   0.,   0.,   0.
     *,1150., .36850,1140., 341.,   0.,   0.,   0.
     *,1150., .37925,1141., 340.,   0.,   0.,   0.
     *,1150., .39000,1141., 341.,   0.,   0.,   0.
     *,1150., .42050,1140.,-120., 110.,   0.,   0.
     *,1150., .45100,1140.,-120., 220.,   0.,   0.
     *,1150., .48150,1140.,-120., 111.,   0.,   0.
     *,1150., .51200,1140.,-120., 221.,   0.,   0.
     *,1150., .54250,1140.,-121., 110.,   0.,   0.
     *,1150., .57300,1140.,-121., 220.,   0.,   0.
     *,1150., .60350,1140.,-121., 111.,   0.,   0.
     *,1150., .63400,1140.,-121., 221.,   0.,   0.
     *,1150., .66450,1141.,-120., 110.,   0.,   0.
     *,1150., .69500,1141.,-120., 220.,   0.,   0.
     *,1150., .72550,1141.,-120., 111.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=523,540)/
     * 1150., .75600,1141.,-120., 221.,   0.,   0.
     *,1150., .78650,1141.,-121., 110.,   0.,   0.
     *,1150., .81700,1141.,-121., 220.,   0.,   0.
     *,1150., .84750,1141.,-121., 111.,   0.,   0.
     *,1150., .87800,1141.,-121., 221.,   0.,   0.
     *,1150., .89325,1140.,-130., 230.,   0.,   0.
     *,1150., .90850,1140.,-130., 231.,   0.,   0.
     *,1150., .92375,1140.,-131., 230.,   0.,   0.
     *,1150., .93900,1140.,-131., 231.,   0.,   0.
     *,1150., .95425,1141.,-130., 230.,   0.,   0.
     *,1150., .96950,1141.,-130., 231.,   0.,   0.
     *,1150., .98475,1141.,-131., 230.,   0.,   0.
     *,1150.,1.00000,1141.,-131., 231.,   0.,   0.
     *,1250., .06000,  12., -11.,1240.,   0.,   0.
     *,1250., .12000,  12., -11.,1241.,   0.,   0.
     *,1250., .18000,  14., -13.,1240.,   0.,   0.
     *,1250., .24000,  14., -13.,1241.,   0.,   0.
     *,1250., .25500,  16., -15.,1240.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=541,558)/
     * 1250., .27000,  16., -15.,1241.,   0.,   0.
     *,1250., .28925,1240.,-120.,   0.,   0.,   0.
     *,1250., .30850,1240.,-121.,   0.,   0.,   0.
     *,1250., .32775,1241.,-120.,   0.,   0.,   0.
     *,1250., .34700,1241.,-121.,   0.,   0.,   0.
     *,1250., .35775,1240., 340.,   0.,   0.,   0.
     *,1250., .36850,1240., 341.,   0.,   0.,   0.
     *,1250., .37925,1241., 340.,   0.,   0.,   0.
     *,1250., .39000,1241., 341.,   0.,   0.,   0.
     *,1250., .42050,1240.,-120., 110.,   0.,   0.
     *,1250., .45100,1240.,-120., 220.,   0.,   0.
     *,1250., .48150,1240.,-120., 111.,   0.,   0.
     *,1250., .51200,1240.,-120., 221.,   0.,   0.
     *,1250., .54250,1240.,-121., 110.,   0.,   0.
     *,1250., .57300,1240.,-121., 220.,   0.,   0.
     *,1250., .60350,1240.,-121., 111.,   0.,   0.
     *,1250., .63400,1240.,-121., 221.,   0.,   0.
     *,1250., .66450,1241.,-120., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=559,576)/
     * 1250., .69500,1241.,-120., 220.,   0.,   0.
     *,1250., .72550,1241.,-120., 111.,   0.,   0.
     *,1250., .75600,1241.,-120., 221.,   0.,   0.
     *,1250., .78650,1241.,-121., 110.,   0.,   0.
     *,1250., .81700,1241.,-121., 220.,   0.,   0.
     *,1250., .84750,1241.,-121., 111.,   0.,   0.
     *,1250., .87800,1241.,-121., 221.,   0.,   0.
     *,1250., .89325,1240.,-130., 230.,   0.,   0.
     *,1250., .90850,1240.,-130., 231.,   0.,   0.
     *,1250., .92375,1240.,-131., 230.,   0.,   0.
     *,1250., .93900,1240.,-131., 231.,   0.,   0.
     *,1250., .95425,1241.,-130., 230.,   0.,   0.
     *,1250., .96950,1241.,-130., 231.,   0.,   0.
     *,1250., .98475,1241.,-131., 230.,   0.,   0.
     *,1250.,1.00000,1241.,-131., 231.,   0.,   0.
     *,1350., .06000,  12., -11.,1340.,   0.,   0.
     *,1350., .12000,  12., -11.,1341.,   0.,   0.
     *,1350., .18000,  14., -13.,1340.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=577,594)/
     * 1350., .24000,  14., -13.,1341.,   0.,   0.
     *,1350., .25500,  16., -15.,1340.,   0.,   0.
     *,1350., .27000,  16., -15.,1341.,   0.,   0.
     *,1350., .28925,1340.,-120.,   0.,   0.,   0.
     *,1350., .30850,1340.,-121.,   0.,   0.,   0.
     *,1350., .32775,1341.,-120.,   0.,   0.,   0.
     *,1350., .34700,1341.,-121.,   0.,   0.,   0.
     *,1350., .35775,1340., 340.,   0.,   0.,   0.
     *,1350., .36850,1340., 341.,   0.,   0.,   0.
     *,1350., .37925,1341., 340.,   0.,   0.,   0.
     *,1350., .39000,1341., 341.,   0.,   0.,   0.
     *,1350., .42050,1340.,-120., 110.,   0.,   0.
     *,1350., .45100,1340.,-120., 220.,   0.,   0.
     *,1350., .48150,1340.,-120., 111.,   0.,   0.
     *,1350., .51200,1340.,-120., 221.,   0.,   0.
     *,1350., .54250,1340.,-121., 110.,   0.,   0.
     *,1350., .57300,1340.,-121., 220.,   0.,   0.
     *,1350., .60350,1340.,-121., 111.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=595,612)/
     * 1350., .63400,1340.,-121., 221.,   0.,   0.
     *,1350., .66450,1341.,-120., 110.,   0.,   0.
     *,1350., .69500,1341.,-120., 220.,   0.,   0.
     *,1350., .72550,1341.,-120., 111.,   0.,   0.
     *,1350., .75600,1341.,-120., 221.,   0.,   0.
     *,1350., .78650,1341.,-121., 110.,   0.,   0.
     *,1350., .81700,1341.,-121., 220.,   0.,   0.
     *,1350., .84750,1341.,-121., 111.,   0.,   0.
     *,1350., .87800,1341.,-121., 221.,   0.,   0.
     *,1350., .89325,1340.,-130., 230.,   0.,   0.
     *,1350., .90850,1340.,-130., 231.,   0.,   0.
     *,1350., .92375,1340.,-131., 230.,   0.,   0.
     *,1350., .93900,1340.,-131., 231.,   0.,   0.
     *,1350., .95425,1341.,-130., 230.,   0.,   0.
     *,1350., .96950,1341.,-130., 231.,   0.,   0.
     *,1350., .98475,1341.,-131., 230.,   0.,   0.
     *,1350.,1.00000,1341.,-131., 231.,   0.,   0.
     *,2150., .06000,  12., -11.,2140.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=613,630)/
     * 2150., .12000,  12., -11.,1241.,   0.,   0.
     *,2150., .18000,  14., -13.,2140.,   0.,   0.
     *,2150., .24000,  14., -13.,1241.,   0.,   0.
     *,2150., .25500,  16., -15.,2140.,   0.,   0.
     *,2150., .27000,  16., -15.,1241.,   0.,   0.
     *,2150., .28925,2140.,-120.,   0.,   0.,   0.
     *,2150., .30850,2140.,-121.,   0.,   0.,   0.
     *,2150., .32775,1241.,-120.,   0.,   0.,   0.
     *,2150., .34700,1241.,-121.,   0.,   0.,   0.
     *,2150., .35775,2140., 340.,   0.,   0.,   0.
     *,2150., .36850,2140., 341.,   0.,   0.,   0.
     *,2150., .37925,1241., 340.,   0.,   0.,   0.
     *,2150., .39000,1241., 341.,   0.,   0.,   0.
     *,2150., .42050,2140.,-120., 110.,   0.,   0.
     *,2150., .45100,2140.,-120., 220.,   0.,   0.
     *,2150., .48150,2140.,-120., 111.,   0.,   0.
     *,2150., .51200,2140.,-120., 221.,   0.,   0.
     *,2150., .54250,2140.,-121., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=631,648)/
     * 2150., .57300,2140.,-121., 220.,   0.,   0.
     *,2150., .60350,2140.,-121., 111.,   0.,   0.
     *,2150., .63400,2140.,-121., 221.,   0.,   0.
     *,2150., .66450,1241.,-120., 110.,   0.,   0.
     *,2150., .69500,1241.,-120., 220.,   0.,   0.
     *,2150., .72550,1241.,-120., 111.,   0.,   0.
     *,2150., .75600,1241.,-120., 221.,   0.,   0.
     *,2150., .78650,1241.,-121., 110.,   0.,   0.
     *,2150., .81700,1241.,-121., 220.,   0.,   0.
     *,2150., .84750,1241.,-121., 111.,   0.,   0.
     *,2150., .87800,1241.,-121., 221.,   0.,   0.
     *,2150., .89325,2140.,-130., 230.,   0.,   0.
     *,2150., .90850,2140.,-130., 231.,   0.,   0.
     *,2150., .92375,2140.,-131., 230.,   0.,   0.
     *,2150., .93900,2140.,-131., 231.,   0.,   0.
     *,2150., .95425,1241.,-130., 230.,   0.,   0.
     *,2150., .96950,1241.,-130., 231.,   0.,   0.
     *,2150., .98475,1241.,-131., 230.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=649,666)/
     * 2150.,1.00000,1241.,-131., 231.,   0.,   0.
     *,2250., .06000,  12., -11.,2240.,   0.,   0.
     *,2250., .12000,  12., -11.,2241.,   0.,   0.
     *,2250., .18000,  14., -13.,2240.,   0.,   0.
     *,2250., .24000,  14., -13.,2241.,   0.,   0.
     *,2250., .25500,  16., -15.,2240.,   0.,   0.
     *,2250., .27000,  16., -15.,2241.,   0.,   0.
     *,2250., .28925,2240.,-120.,   0.,   0.,   0.
     *,2250., .30850,2240.,-121.,   0.,   0.,   0.
     *,2250., .32775,2241.,-120.,   0.,   0.,   0.
     *,2250., .34700,2241.,-121.,   0.,   0.,   0.
     *,2250., .35775,2240., 340.,   0.,   0.,   0.
     *,2250., .36850,2240., 341.,   0.,   0.,   0.
     *,2250., .37925,2241., 340.,   0.,   0.,   0.
     *,2250., .39000,2241., 341.,   0.,   0.,   0.
     *,2250., .42050,2240.,-120., 110.,   0.,   0.
     *,2250., .45100,2240.,-120., 220.,   0.,   0.
     *,2250., .48150,2240.,-120., 111.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=667,684)/
     * 2250., .51200,2240.,-120., 221.,   0.,   0.
     *,2250., .54250,2240.,-121., 110.,   0.,   0.
     *,2250., .57300,2240.,-121., 220.,   0.,   0.
     *,2250., .60350,2240.,-121., 111.,   0.,   0.
     *,2250., .63400,2240.,-121., 221.,   0.,   0.
     *,2250., .66450,2241.,-120., 110.,   0.,   0.
     *,2250., .69500,2241.,-120., 220.,   0.,   0.
     *,2250., .72550,2241.,-120., 111.,   0.,   0.
     *,2250., .75600,2241.,-120., 221.,   0.,   0.
     *,2250., .78650,2241.,-121., 110.,   0.,   0.
     *,2250., .81700,2241.,-121., 220.,   0.,   0.
     *,2250., .84750,2241.,-121., 111.,   0.,   0.
     *,2250., .87800,2241.,-121., 221.,   0.,   0.
     *,2250., .89325,2240.,-130., 230.,   0.,   0.
     *,2250., .90850,2240.,-130., 231.,   0.,   0.
     *,2250., .92375,2240.,-131., 230.,   0.,   0.
     *,2250., .93900,2240.,-131., 231.,   0.,   0.
     *,2250., .95425,2241.,-130., 230.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=685,702)/
     * 2250., .96950,2241.,-130., 231.,   0.,   0.
     *,2250., .98475,2241.,-131., 230.,   0.,   0.
     *,2250.,1.00000,2241.,-131., 231.,   0.,   0.
     *,2350., .06000,  12., -11.,2340.,   0.,   0.
     *,2350., .12000,  12., -11.,2341.,   0.,   0.
     *,2350., .18000,  14., -13.,2340.,   0.,   0.
     *,2350., .24000,  14., -13.,2341.,   0.,   0.
     *,2350., .25500,  16., -15.,2340.,   0.,   0.
     *,2350., .27000,  16., -15.,2341.,   0.,   0.
     *,2350., .28925,2340.,-120.,   0.,   0.,   0.
     *,2350., .30850,2340.,-121.,   0.,   0.,   0.
     *,2350., .32775,2341.,-120.,   0.,   0.,   0.
     *,2350., .34700,2341.,-121.,   0.,   0.,   0.
     *,2350., .35775,2340., 340.,   0.,   0.,   0.
     *,2350., .36850,2340., 341.,   0.,   0.,   0.
     *,2350., .37925,2341., 340.,   0.,   0.,   0.
     *,2350., .39000,2341., 341.,   0.,   0.,   0.
     *,2350., .42050,2340.,-120., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=703,720)/
     * 2350., .45100,2340.,-120., 220.,   0.,   0.
     *,2350., .48150,2340.,-120., 111.,   0.,   0.
     *,2350., .51200,2340.,-120., 221.,   0.,   0.
     *,2350., .54250,2340.,-121., 110.,   0.,   0.
     *,2350., .57300,2340.,-121., 220.,   0.,   0.
     *,2350., .60350,2340.,-121., 111.,   0.,   0.
     *,2350., .63400,2340.,-121., 221.,   0.,   0.
     *,2350., .66450,2341.,-120., 110.,   0.,   0.
     *,2350., .69500,2341.,-120., 220.,   0.,   0.
     *,2350., .72550,2341.,-120., 111.,   0.,   0.
     *,2350., .75600,2341.,-120., 221.,   0.,   0.
     *,2350., .78650,2341.,-121., 110.,   0.,   0.
     *,2350., .81700,2341.,-121., 220.,   0.,   0.
     *,2350., .84750,2341.,-121., 111.,   0.,   0.
     *,2350., .87800,2341.,-121., 221.,   0.,   0.
     *,2350., .89325,2340.,-130., 230.,   0.,   0.
     *,2350., .90850,2340.,-130., 231.,   0.,   0.
     *,2350., .92375,2340.,-131., 230.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=721,738)/
     * 2350., .93900,2340.,-131., 231.,   0.,   0.
     *,2350., .95425,2341.,-130., 230.,   0.,   0.
     *,2350., .96950,2341.,-130., 231.,   0.,   0.
     *,2350., .98475,2341.,-131., 230.,   0.,   0.
     *,2350.,1.00000,2341.,-131., 231.,   0.,   0.
     *,3150., .06000,  12., -11.,3140.,   0.,   0.
     *,3150., .12000,  12., -11.,1341.,   0.,   0.
     *,3150., .18000,  14., -13.,3140.,   0.,   0.
     *,3150., .24000,  14., -13.,1341.,   0.,   0.
     *,3150., .25500,  16., -15.,3140.,   0.,   0.
     *,3150., .27000,  16., -15.,1341.,   0.,   0.
     *,3150., .28925,3140.,-120.,   0.,   0.,   0.
     *,3150., .30850,3140.,-121.,   0.,   0.,   0.
     *,3150., .32775,1341.,-120.,   0.,   0.,   0.
     *,3150., .34700,1341.,-121.,   0.,   0.,   0.
     *,3150., .35775,3140., 340.,   0.,   0.,   0.
     *,3150., .36850,3140., 341.,   0.,   0.,   0.
     *,3150., .37925,1341., 340.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=739,756)/
     * 3150., .39000,1341., 341.,   0.,   0.,   0.
     *,3150., .42050,3140.,-120., 110.,   0.,   0.
     *,3150., .45100,3140.,-120., 220.,   0.,   0.
     *,3150., .48150,3140.,-120., 111.,   0.,   0.
     *,3150., .51200,3140.,-120., 221.,   0.,   0.
     *,3150., .54250,3140.,-121., 110.,   0.,   0.
     *,3150., .57300,3140.,-121., 220.,   0.,   0.
     *,3150., .60350,3140.,-121., 111.,   0.,   0.
     *,3150., .63400,3140.,-121., 221.,   0.,   0.
     *,3150., .66450,1341.,-120., 110.,   0.,   0.
     *,3150., .69500,1341.,-120., 220.,   0.,   0.
     *,3150., .72550,1341.,-120., 111.,   0.,   0.
     *,3150., .75600,1341.,-120., 221.,   0.,   0.
     *,3150., .78650,1341.,-121., 110.,   0.,   0.
     *,3150., .81700,1341.,-121., 220.,   0.,   0.
     *,3150., .84750,1341.,-121., 111.,   0.,   0.
     *,3150., .87800,1341.,-121., 221.,   0.,   0.
     *,3150., .89325,3140.,-130., 230.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=757,774)/
     * 3150., .90850,3140.,-130., 231.,   0.,   0.
     *,3150., .92375,3140.,-131., 230.,   0.,   0.
     *,3150., .93900,3140.,-131., 231.,   0.,   0.
     *,3150., .95425,1341.,-130., 230.,   0.,   0.
     *,3150., .96950,1341.,-130., 231.,   0.,   0.
     *,3150., .98475,1341.,-131., 230.,   0.,   0.
     *,3150.,1.00000,1341.,-131., 231.,   0.,   0.
     *,3250., .06000,  12., -11.,3240.,   0.,   0.
     *,3250., .12000,  12., -11.,2341.,   0.,   0.
     *,3250., .18000,  14., -13.,3240.,   0.,   0.
     *,3250., .24000,  14., -13.,2341.,   0.,   0.
     *,3250., .25500,  16., -15.,3240.,   0.,   0.
     *,3250., .27000,  16., -15.,2341.,   0.,   0.
     *,3250., .28925,3240.,-120.,   0.,   0.,   0.
     *,3250., .30850,3240.,-121.,   0.,   0.,   0.
     *,3250., .32775,2341.,-120.,   0.,   0.,   0.
     *,3250., .34700,2341.,-121.,   0.,   0.,   0.
     *,3250., .35775,3240., 340.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=775,792)/
     * 3250., .36850,3240., 341.,   0.,   0.,   0.
     *,3250., .37925,2341., 340.,   0.,   0.,   0.
     *,3250., .39000,2341., 341.,   0.,   0.,   0.
     *,3250., .42050,3240.,-120., 110.,   0.,   0.
     *,3250., .45100,3240.,-120., 220.,   0.,   0.
     *,3250., .48150,3240.,-120., 111.,   0.,   0.
     *,3250., .51200,3240.,-120., 221.,   0.,   0.
     *,3250., .54250,3240.,-121., 110.,   0.,   0.
     *,3250., .57300,3240.,-121., 220.,   0.,   0.
     *,3250., .60350,3240.,-121., 111.,   0.,   0.
     *,3250., .63400,3240.,-121., 221.,   0.,   0.
     *,3250., .66450,2341.,-120., 110.,   0.,   0.
     *,3250., .69500,2341.,-120., 220.,   0.,   0.
     *,3250., .72550,2341.,-120., 111.,   0.,   0.
     *,3250., .75600,2341.,-120., 221.,   0.,   0.
     *,3250., .78650,2341.,-121., 110.,   0.,   0.
     *,3250., .81700,2341.,-121., 220.,   0.,   0.
     *,3250., .84750,2341.,-121., 111.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=793,810)/
     * 3250., .87800,2341.,-121., 221.,   0.,   0.
     *,3250., .89325,3240.,-130., 230.,   0.,   0.
     *,3250., .90850,3240.,-130., 231.,   0.,   0.
     *,3250., .92375,3240.,-131., 230.,   0.,   0.
     *,3250., .93900,3240.,-131., 231.,   0.,   0.
     *,3250., .95425,2341.,-130., 230.,   0.,   0.
     *,3250., .96950,2341.,-130., 231.,   0.,   0.
     *,3250., .98475,2341.,-131., 230.,   0.,   0.
     *,3250.,1.00000,2341.,-131., 231.,   0.,   0.
     *,3350., .06000,  12., -11.,3340.,   0.,   0.
     *,3350., .12000,  12., -11.,3341.,   0.,   0.
     *,3350., .18000,  14., -13.,3340.,   0.,   0.
     *,3350., .24000,  14., -13.,3341.,   0.,   0.
     *,3350., .25500,  16., -15.,3340.,   0.,   0.
     *,3350., .27000,  16., -15.,3341.,   0.,   0.
     *,3350., .28925,3340.,-120.,   0.,   0.,   0.
     *,3350., .30850,3340.,-121.,   0.,   0.,   0.
     *,3350., .32775,3341.,-120.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=811,828)/
     * 3350., .34700,3341.,-121.,   0.,   0.,   0.
     *,3350., .35775,3340., 340.,   0.,   0.,   0.
     *,3350., .36850,3340., 341.,   0.,   0.,   0.
     *,3350., .37925,3341., 340.,   0.,   0.,   0.
     *,3350., .39000,3341., 341.,   0.,   0.,   0.
     *,3350., .42050,3340.,-120., 110.,   0.,   0.
     *,3350., .45100,3340.,-120., 220.,   0.,   0.
     *,3350., .48150,3340.,-120., 111.,   0.,   0.
     *,3350., .51200,3340.,-120., 221.,   0.,   0.
     *,3350., .54250,3340.,-121., 110.,   0.,   0.
     *,3350., .57300,3340.,-121., 220.,   0.,   0.
     *,3350., .60350,3340.,-121., 111.,   0.,   0.
     *,3350., .63400,3340.,-121., 221.,   0.,   0.
     *,3350., .66450,3341.,-120., 110.,   0.,   0.
     *,3350., .69500,3341.,-120., 220.,   0.,   0.
     *,3350., .72550,3341.,-120., 111.,   0.,   0.
     *,3350., .75600,3341.,-120., 221.,   0.,   0.
     *,3350., .78650,3341.,-121., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=829,846)/
     * 3350., .81700,3341.,-121., 220.,   0.,   0.
     *,3350., .84750,3341.,-121., 111.,   0.,   0.
     *,3350., .87800,3341.,-121., 221.,   0.,   0.
     *,3350., .89325,3340.,-130., 230.,   0.,   0.
     *,3350., .90850,3340.,-130., 231.,   0.,   0.
     *,3350., .92375,3340.,-131., 230.,   0.,   0.
     *,3350., .93900,3340.,-131., 231.,   0.,   0.
     *,3350., .95425,3341.,-130., 230.,   0.,   0.
     *,3350., .96950,3341.,-130., 231.,   0.,   0.
     *,3350., .98475,3341.,-131., 230.,   0.,   0.
     *,3350.,1.00000,3341.,-131., 231.,   0.,   0.
     *,1160., .33300,   1.,  -2.,1500.,   0.,   0.
     *,1160., .66700,   4.,  -3.,1500.,   0.,   0.
     *,1160., .77800, -12.,  11.,1500.,   0.,   0.
     *,1160., .88900, -14.,  13.,1500.,   0.,   0.
     *,1160.,1.00000, -16.,  15.,1500.,   0.,   0.
     *,1260., .33300,   1.,  -2.,2500.,   0.,   0.
     *,1260., .66700,   4.,  -3.,2500.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=847,864)/
     * 1260., .77800, -12.,  11.,2500.,   0.,   0.
     *,1260., .88900, -14.,  13.,2500.,   0.,   0.
     *,1260.,1.00000, -16.,  15.,2500.,   0.,   0.
     *,2260., .33300,   1.,  -2.,2500.,   0.,   0.
     *,2260., .66700,   4.,  -3.,2500.,   0.,   0.
     *,2260., .77800, -12.,  11.,2500.,   0.,   0.
     *,2260., .88900, -14.,  13.,2500.,   0.,   0.
     *,2260.,1.00000, -16.,  15.,2500.,   0.,   0.
     *,2160., .33300,   1.,  -2.,1500.,   0.,   0.
     *,2160., .66700,   4.,  -3.,1500.,   0.,   0.
     *,2160., .77800, -12.,  11.,1500.,   0.,   0.
     *,2160., .88900, -14.,  13.,1500.,   0.,   0.
     *,2160.,1.00000, -16.,  15.,1500.,   0.,   0.
     *,1360., .33300,   1.,  -2.,3500.,   0.,   0.
     *,1360., .66700,   4.,  -3.,3500.,   0.,   0.
     *,1360., .77800, -12.,  11.,3500.,   0.,   0.
     *,1360., .88900, -14.,  13.,3500.,   0.,   0.
     *,1360.,1.00000, -16.,  15.,3500.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=865,882)/
     * 2360., .33300,   1.,  -2.,3500.,   0.,   0.
     *,2360., .66700,   4.,  -3.,3500.,   0.,   0.
     *,2360., .77800, -12.,  11.,3500.,   0.,   0.
     *,2360., .88900, -14.,  13.,3500.,   0.,   0.
     *,2360.,1.00000, -16.,  15.,3500.,   0.,   0.
     *,3360., .33300,   1.,  -2.,3500.,   0.,   0.
     *,3360., .66700,   4.,  -3.,3500.,   0.,   0.
     *,3360., .77800, -12.,  11.,3500.,   0.,   0.
     *,3360., .88900, -14.,  13.,3500.,   0.,   0.
     *,3360.,1.00000, -16.,  15.,3500.,   0.,   0.
     *,1151.,1.00000,1150.,  10.,   0.,   0.,   0.
     *,1251.,1.00000,1250.,  10.,   0.,   0.,   0.
     *,2251.,1.00000,2250.,  10.,   0.,   0.,   0.
     *,1351.,1.00000,1350.,  10.,   0.,   0.,   0.
     *,2351.,1.00000,2350.,  10.,   0.,   0.,   0.
     *,3351.,1.00000,3350.,  10.,   0.,   0.,   0.
     *,1161.,1.00000,1160.,  10.,   0.,   0.,   0.
     *,1261.,1.00000,1260.,  10.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=883,886)/
     * 2261.,1.00000,2260.,  10.,   0.,   0.,   0.
     *,1361.,1.00000,1360.,  10.,   0.,   0.,   0.
     *,2361.,1.00000,2360.,  10.,   0.,   0.,   0.
     *,3361.,1.00000,3360.,  10.,   0.,   0.,   0./
c    *---------------------------------------------
c    *    delta++ resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=887,900)/
c    *--dl++(1620)---------------------------------
     * 1112., .30000,1120., 120.,   0.,   0.,   0.
     *,1112., .66000,1111., 110.,   0.,   0.,   0.
     *,1112., .90000,1121., 120.,   0.,   0.,   0.
     *,1112.,1.00000,1120., 120., 110.,   0.,   0.
c    *--dl++(1700)---------------------------------
     *,1113., .15000,1120., 120.,   0.,   0.,   0.
     *,1113., .51000,1111., 110.,   0.,   0.,   0.
     *,1113., .75000,1121., 120.,   0.,   0.,   0.
     *,1113.,1.00000,1120., 120., 110.,   0.,   0.
c    *--dl++(1925)---------------------------------
     *,1114., .28000,1120., 120.,   0.,   0.,   0.
     *,1114., .40600,1111., 110.,   0.,   0.,   0.
     *,1114., .49000,1121., 120.,   0.,   0.,   0.
     *,1114., .69000,1120., 121.,   0.,   0.,   0.
     *,1114., .70000,1130., 130.,   0.,   0.,   0.
     *,1114.,1.00000,1122., 120.,   0.,   0.,   0./
c    *---------------------------------------------
c    *    delta- resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=901,914)/
c    *--dl-(1620)----------------------------------
     * 2222., .30000,1220.,-120.,   0.,   0.,   0.
     *,2222., .66000,2221., 110.,   0.,   0.,   0.
     *,2222., .90000,1221.,-120.,   0.,   0.,   0.
     *,2222.,1.00000,1220., 110.,-120.,   0.,   0.
c    *--dl-(1700)----------------------------------
     *,2223., .15000,1220.,-120.,   0.,   0.,   0.
     *,2223., .51000,2221., 110.,   0.,   0.,   0.
     *,2223., .75000,1221.,-120.,   0.,   0.,   0.
     *,2223.,1.00000,1220., 110.,-120.,   0.,   0.
c    *--dl-(1925)----------------------------------
     *,2224., .28000,1220.,-120.,   0.,   0.,   0.
     *,2224., .40600,2221., 110.,   0.,   0.,   0.
     *,2224., .49000,1221.,-120.,   0.,   0.,   0.
     *,2224., .69000,1220.,-121.,   0.,   0.,   0.
     *,2224., .70000,2230., 230.,   0.,   0.,   0.
     *,2224.,1.00000,1222.,-120.,   0.,   0.,   0./
c    *---------------------------------------------
c    *    n*+ resonances + delta+ resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=915,931)/
c    *--n*+(1440)----------------------------------
     * 1122., .20000,1120., 110.,   0.,   0.,   0.
     *,1122., .60000,1220., 120.,   0.,   0.,   0.
     *,1122., .68000,1111.,-120.,   0.,   0.,   0.
     *,1122., .73000,1121., 110.,   0.,   0.,   0.
     *,1122., .76000,1221., 120.,   0.,   0.,   0.
     *,1122., .84000,1120., 120.,-120.,   0.,   0.
     *,1122., .87000,1120., 110., 110.,   0.,   0.
     *,1122.,1.00000,1220., 120., 110.,   0.,   0.
c    *--n*+(1530)----------------------------------
     *,1123., .17000,1120., 110.,   0.,   0.,   0.
     *,1123., .51000,1220., 120.,   0.,   0.,   0.
     *,1123., .57000,1111.,-120.,   0.,   0.,   0.
     *,1123., .61000,1121., 110.,   0.,   0.,   0.
     *,1123., .63000,1221., 120.,   0.,   0.,   0.
     *,1123., .67000,1120., 120.,-120.,   0.,   0.
     *,1123., .68000,1120., 110., 110.,   0.,   0.
     *,1123., .75000,1220., 120., 110.,   0.,   0.
     *,1123.,1.00000,1120., 220.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=932,948)/
c    *--dl+(1620)----------------------------------
     * 1124., .20000,1120., 110.,   0.,   0.,   0.
     *,1124., .30000,1220., 120.,   0.,   0.,   0.
     *,1124., .54000,1111.,-120.,   0.,   0.,   0.
     *,1124., .58000,1121., 110.,   0.,   0.,   0.
     *,1124., .90000,1221., 120.,   0.,   0.,   0.
     *,1124., .96000,1120., 120.,-120.,   0.,   0.
     *,1124.,1.00000,1220., 120., 110.,   0.,   0.
c    *--n*+(1665)----------------------------------
     *,1125., .16700,1120., 110.,   0.,   0.,   0.
     *,1125., .49970,1220., 120.,   0.,   0.,   0.
     *,1125., .62470,1111.,-120.,   0.,   0.,   0.
     *,1125., .70800,1121., 110.,   0.,   0.,   0.
     *,1125., .74970,1221., 120.,   0.,   0.,   0.
     *,1125., .82080,1120., 120.,-120.,   0.,   0.
     *,1125., .85190,1120., 110., 110.,   0.,   0.
     *,1125., .96300,1220., 120., 110.,   0.,   0.
     *,1125., .97300,1120., 220.,   0.,   0.,   0.
     *,1125.,1.00000,2130., 130.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=949,955)/
c    *--dl+(1700)----------------------------------
     * 1126., .10000,1120., 110.,   0.,   0.,   0.
     *,1126., .15000,1220., 120.,   0.,   0.,   0.
     *,1126., .39000,1111.,-120.,   0.,   0.,   0.
     *,1126., .43000,1121., 110.,   0.,   0.,   0.
     *,1126., .75000,1221., 120.,   0.,   0.,   0.
     *,1126., .91500,1120., 120.,-120.,   0.,   0.
     *,1126.,1.00000,1220., 120., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=956,969)/
c    *--n*+(1710)----------------------------------
     * 1127., .04430,1120., 110.,   0.,   0.,   0.
     *,1127., .13290,1220., 120.,   0.,   0.,   0.
     *,1127., .23790,1111.,-120.,   0.,   0.,   0.
     *,1127., .30790,1121., 110.,   0.,   0.,   0.
     *,1127., .34290,1221., 120.,   0.,   0.,   0.
     *,1127., .41190,1120., 120.,-120.,   0.,   0.
     *,1127., .48090,1120., 110., 110.,   0.,   0.
     *,1127., .54990,1220., 120., 110.,   0.,   0.
     *,1127., .66070,1120., 220.,   0.,   0.,   0.
     *,1127., .72800,2130., 130.,   0.,   0.,   0.
     *,1127., .74930,1230., 130.,   0.,   0.,   0.
     *,1127., .76000,1130., 230.,   0.,   0.,   0.
     *,1127., .84000,1120., 111.,   0.,   0.,   0.
     *,1127.,1.00000,1220., 121.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=970,980)/
c    *--dl+(1925)----------------------------------
     * 1128., .18700,1120., 110.,   0.,   0.,   0.
     *,1128., .28000,1220., 120.,   0.,   0.,   0.
     *,1128., .36400,1111.,-120.,   0.,   0.,   0.
     *,1128., .37800,1121., 110.,   0.,   0.,   0.
     *,1128., .49000,1221., 120.,   0.,   0.,   0.
     *,1128., .62300,1120., 111.,   0.,   0.,   0.
     *,1128., .69000,1220., 121.,   0.,   0.,   0.
     *,1128., .69350,1130., 230.,   0.,   0.,   0.
     *,1128., .69900,1230., 130.,   0.,   0.,   0.
     *,1128., .89900,1122., 110.,   0.,   0.,   0.
     *,1128.,1.00000,1222., 120.,   0.,   0.,   0./
c    *---------------------------------------------
c    *    n*0  resonances + delta0 resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=981,997)/
c    *----------n*0(1440)--------------------------
     * 1222., .20000,1220., 110.,   0.,   0.,   0.
     *,1222., .60000,1120.,-120.,   0.,   0.,   0.
     *,1222., .68000,2221., 120.,   0.,   0.,   0.
     *,1222., .73000,1221., 110.,   0.,   0.,   0.
     *,1222., .76000,1121.,-120.,   0.,   0.,   0.
     *,1222., .84000,1220., 120.,-120.,   0.,   0.
     *,1222., .87000,1220., 110., 110.,   0.,   0.
     *,1222.,1.00000,1120.,-120., 110.,   0.,   0.
c    *----------n*0(1530)--------------------------
     *,1223., .17000,1220., 110.,   0.,   0.,   0.
     *,1223., .51000,1120.,-120.,   0.,   0.,   0.
     *,1223., .57000,2221., 120.,   0.,   0.,   0.
     *,1223., .61000,1221., 110.,   0.,   0.,   0.
     *,1223., .63000,1121.,-120.,   0.,   0.,   0.
     *,1223., .67000,1220., 120.,-120.,   0.,   0.
     *,1223., .68000,1220., 110., 110.,   0.,   0.
     *,1223., .75000,1120.,-120., 110.,   0.,   0.
     *,1223.,1.00000,1220., 220.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=998,1014)/
c    *----------dl0(1620)--------------------------
     * 1224., .20000,1220., 110.,   0.,   0.,   0.
     *,1224., .30000,1120.,-120.,   0.,   0.,   0.
     *,1224., .54000,2221., 120.,   0.,   0.,   0.
     *,1224., .58000,1221., 110.,   0.,   0.,   0.
     *,1224., .90000,1121.,-120.,   0.,   0.,   0.
     *,1224., .96500,1220., 120.,-120.,   0.,   0.
     *,1224.,1.00000,1120.,-120., 110.,   0.,   0.
c    *----------n*0(1665)--------------------------
     *,1225., .16700,1220., 110.,   0.,   0.,   0.
     *,1225., .49970,1120.,-120.,   0.,   0.,   0.
     *,1225., .62470,2221., 120.,   0.,   0.,   0.
     *,1225., .70800,1221., 110.,   0.,   0.,   0.
     *,1225., .74970,1121.,-120.,   0.,   0.,   0.
     *,1225., .82080,1220., 120.,-120.,   0.,   0.
     *,1225., .85190,1220., 110., 110.,   0.,   0.
     *,1225., .96300,1120.,-120., 110.,   0.,   0.
     *,1225., .97300,1220., 220.,   0.,   0.,   0.
     *,1225.,1.00000,2130., 230.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1015,1021)/
c    *----------dl0(1700)--------------------------
     * 1226., .10000,1220., 110.,   0.,   0.,   0.
     *,1226., .15000,1120.,-120.,   0.,   0.,   0.
     *,1226., .39000,2221., 120.,   0.,   0.,   0.
     *,1226., .43000,1221., 110.,   0.,   0.,   0.
     *,1226., .75000,1121.,-120.,   0.,   0.,   0.
     *,1226., .91500,1220., 120.,-120.,   0.,   0.
     *,1226.,1.00000,1120.,-120., 110.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1022,1035)/
c    *----------n*0(1710)--------------------------
     * 1227., .04430,1220., 110.,   0.,   0.,   0.
     *,1227., .13290,1120.,-120.,   0.,   0.,   0.
     *,1227., .23790,2221., 120.,   0.,   0.,   0.
     *,1227., .30790,1221., 110.,   0.,   0.,   0.
     *,1227., .34290,1121.,-120.,   0.,   0.,   0.
     *,1227., .41190,1220., 120.,-120.,   0.,   0.
     *,1227., .48090,1220., 110., 110.,   0.,   0.
     *,1227., .54990,1120.,-120., 110.,   0.,   0.
     *,1227., .66070,1220., 220.,   0.,   0.,   0.
     *,1227., .72800,2130., 230.,   0.,   0.,   0.
     *,1227., .73870,1230., 230.,   0.,   0.,   0.
     *,1227., .76000,2230., 130.,   0.,   0.,   0.
     *,1227., .92000,1120.,-121.,   0.,   0.,   0.
     *,1227.,1.00000,1220., 111.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1036,1046)/
c    *----------dl0(1925)--------------------------
     * 1228., .18700,1220., 110.,   0.,   0.,   0.
     *,1228., .28000,1120.,-120.,   0.,   0.,   0.
     *,1228., .36400,2221., 120.,   0.,   0.,   0.
     *,1228., .37800,1221., 110.,   0.,   0.,   0.
     *,1228., .49000,1121.,-120.,   0.,   0.,   0.
     *,1228., .55700,1220., 111.,   0.,   0.,   0.
     *,1228., .69000,1120.,-121.,   0.,   0.,   0.
     *,1228., .69350,2230., 130.,   0.,   0.,   0.
     *,1228., .70000,1230., 230.,   0.,   0.,   0.
     *,1228., .80000,1122.,-120.,   0.,   0.,   0.
     *,1228.,1.00000,1222., 110.,   0.,   0.,   0./
c    *---------------------------------------------
c    *   lambda resonances + sigma0 resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=1047,1059)/
c    *----------lambda(1405)-----------------------
     * 1233., .33000,1230., 110.,   0.,   0.,   0.
     *,1233., .66000,2230., 120.,   0.,   0.,   0.
     *,1233.,1.00000,1130.,-120.,   0.,   0.,   0.
c    *----------lambda(1520)-----------------------
     *,1234., .22500,1120.,-130.,   0.,   0.,   0.
     *,1234., .48000,1220.,-230.,   0.,   0.,   0.
     *,1234., .62000,1230., 110.,   0.,   0.,   0.
     *,1234., .76000,2230., 120.,   0.,   0.,   0.
     *,1234., .90000,1130.,-120.,   0.,   0.,   0.
     *,1234., .96000,2130., 120.,-120.,   0.,   0.
     *,1234., .99000,2130., 110., 110.,   0.,   0.
     *,1234., .99330,1130.,-120., 110.,   0.,   0.
     *,1234., .99660,2230., 120., 110.,   0.,   0.
     *,1234.,1.00000,1230., 120.,-120.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1060,1075)/
c    *----------lambda(1645)-----------------------
     * 1235., .10000,1120.,-130.,   0.,   0.,   0.
     *,1235., .20000,1220.,-230.,   0.,   0.,   0.
     *,1235., .35000,1230., 110.,   0.,   0.,   0.
     *,1235., .50000,2230., 120.,   0.,   0.,   0.
     *,1235., .65000,1130.,-120.,   0.,   0.,   0.
     *,1235., .75000,2130., 120.,-120.,   0.,   0.
     *,1235., .80000,2130., 110., 110.,   0.,   0.
     *,1235., .84500,1130.,-120., 110.,   0.,   0.
     *,1235., .89000,2230., 120., 110.,   0.,   0.
     *,1235., .93500,1230., 120.,-120.,   0.,   0.
     *,1235.,1.00000,2130., 220.,   0.,   0.,   0.
c    *----------sigma0(1665)-----------------------
     *,1236., .10000,1120.,-130.,   0.,   0.,   0.
     *,1236., .20000,1220.,-230.,   0.,   0.,   0.
     *,1236., .40000,2230., 120.,   0.,   0.,   0.
     *,1236., .60000,1130.,-120.,   0.,   0.,   0.
     *,1236.,1.00000,2130., 110.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1076,1084)/
c    *----------sigma0(1776)-----------------------
     * 1237., .17500,1120.,-130.,   0.,   0.,   0.
     *,1237., .35000,1220.,-230.,   0.,   0.,   0.
     *,1237., .38750,2230., 120.,   0.,   0.,   0.
     *,1237., .42500,1130.,-120.,   0.,   0.,   0.
     *,1237., .57500,2130., 110.,   0.,   0.,   0.
     *,1237., .60000,2231., 120.,   0.,   0.,   0.
     *,1237., .62500,1131.,-120.,   0.,   0.,   0.
     *,1237., .75000,1234., 110.,   0.,   0.,   0.
     *,1237.,1.00000,1230., 220.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1085,1094)/
c    *----------lambda(1845)-----------------------
     * 1238., .17000,1120.,-130.,   0.,   0.,   0.
     *,1238., .34000,1220.,-230.,   0.,   0.,   0.
     *,1238., .44000,1230., 110.,   0.,   0.,   0.
     *,1238., .54000,2230., 120.,   0.,   0.,   0.
     *,1238., .64000,1130.,-120.,   0.,   0.,   0.
     *,1238., .70000,1231., 110.,   0.,   0.,   0.
     *,1238., .76000,2231., 120.,   0.,   0.,   0.
     *,1238., .82000,1131.,-120.,   0.,   0.,   0.
     *,1238., .91000,1120.,-131.,   0.,   0.,   0.
     *,1238.,1.00000,1220.,-231.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1095,1106)/
c    *----------sigma0(1930)-----------------------
     * 1239., .07500,1120.,-130.,   0.,   0.,   0.
     *,1239., .15000,1220.,-230.,   0.,   0.,   0.
     *,1239., .20000,1121.,-130.,   0.,   0.,   0.
     *,1239., .25000,1221.,-230.,   0.,   0.,   0.
     *,1239., .32500,1120.,-131.,   0.,   0.,   0.
     *,1239., .40000,1220.,-231.,   0.,   0.,   0.
     *,1239., .47500,2230., 120.,   0.,   0.,   0.
     *,1239., .55000,1130.,-120.,   0.,   0.,   0.
     *,1239., .70000,2130., 110.,   0.,   0.,   0.
     *,1239., .77500,2231., 120.,   0.,   0.,   0.
     *,1239., .85000,1131.,-120.,   0.,   0.,   0.
     *,1239.,1.00000,1234., 110.,   0.,   0.,   0./
c    *---------------------------------------------
c    *            sigma+ resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=1107,1118)/
c    *----------sigma+(1665)-----------------------
     * 1132., .20000,1120.,-230.,   0.,   0.,   0.
     *,1132., .40000,1130., 110.,   0.,   0.,   0.
     *,1132., .60000,1230., 120.,   0.,   0.,   0.
     *,1132.,1.00000,2130., 120.,   0.,   0.,   0.
c    *----------sigma+(1776)-----------------------
     *,1133., .35000,1120.,-230.,   0.,   0.,   0.
     *,1133., .38750,1130., 110.,   0.,   0.,   0.
     *,1133., .42500,1230., 120.,   0.,   0.,   0.
     *,1133., .57500,2130., 120.,   0.,   0.,   0.
     *,1133., .60000,1131., 110.,   0.,   0.,   0.
     *,1133., .62500,1231., 120.,   0.,   0.,   0.
     *,1133., .75000,1234., 120.,   0.,   0.,   0.
     *,1133.,1.00000,1130., 220.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1119,1128)/
c    *----------sigma+(1930)-----------------------
     * 1134., .15000,1120.,-230.,   0.,   0.,   0.
     *,1134., .22500,1111.,-130.,   0.,   0.,   0.
     *,1134., .25000,1121.,-230.,   0.,   0.,   0.
     *,1134., .40000,1120.,-231.,   0.,   0.,   0.
     *,1134., .47500,1130., 110.,   0.,   0.,   0.
     *,1134., .55000,1230., 120.,   0.,   0.,   0.
     *,1134., .70000,2130., 120.,   0.,   0.,   0.
     *,1134., .77500,1131., 110.,   0.,   0.,   0.
     *,1134., .85000,1231., 120.,   0.,   0.,   0.
     *,1134.,1.00000,1234., 120.,   0.,   0.,   0./
c    *---------------------------------------------
c    *            sigma- resonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=1129,1140)/
c    *----------sigma-(1665)-----------------------
     * 2232., .20000,1220.,-130.,   0.,   0.,   0.
     *,2232., .40000,2230., 110.,   0.,   0.,   0.
     *,2232., .60000,1230.,-120.,   0.,   0.,   0.
     *,2232.,1.00000,2130.,-120.,   0.,   0.,   0.
c    *----------sigma-(1776)-----------------------
     *,2233., .35000,1220.,-130.,   0.,   0.,   0.
     *,2233., .38750,2230., 110.,   0.,   0.,   0.
     *,2233., .42500,1230.,-120.,   0.,   0.,   0.
     *,2233., .57500,2130.,-120.,   0.,   0.,   0.
     *,2233., .60000,2231., 110.,   0.,   0.,   0.
     *,2233., .62500,1231.,-120.,   0.,   0.,   0.
     *,2233., .75000,1234.,-120.,   0.,   0.,   0.
     *,2233.,1.00000,2230., 220.,   0.,   0.,   0./
      data ((dectab(i,j),i=1,7),j=1141,1150)/
c    *----------sigma-(1930)-----------------------
     * 2234., .15000,1220.,-130.,   0.,   0.,   0.
     *,2234., .17500,1221.,-130.,   0.,   0.,   0.
     *,2234., .25000,2221.,-230.,   0.,   0.,   0.
     *,2234., .40000,1220.,-131.,   0.,   0.,   0.
     *,2234., .47500,2230., 110.,   0.,   0.,   0.
     *,2234., .55000,1230.,-120.,   0.,   0.,   0.
     *,2234., .70000,2130.,-120.,   0.,   0.,   0.
     *,2234., .77500,2231., 110.,   0.,   0.,   0.
     *,2234., .85000,1231.,-120.,   0.,   0.,   0.
     *,2234.,1.00000,1234.,-120.,   0.,   0.,   0./
c    *---------------------------------------------
c    *      additional mesonresonances
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=1151,1159)/
c    *-----------f0(975)---------------------------
     *  332., .50000, 120.,-120.,   0.,   0.,   0.
     *, 332., .75000, 110., 110.,   0.,   0.,   0.
     *, 332., .87500, 130.,-130.,   0.,   0.,   0.
     *, 332.,1.00000, 230.,-230.,   0.,   0.,   0.
c    *-----------a0(980)---------------------------
     *, 112., .56000, 110., 220.,   0.,   0.,   0.
     *, 112., .78000, 130.,-130.,   0.,   0.,   0.
     *, 112.,1.00000, 230.,-230.,   0.,   0.,   0.
c    *-----------a+(980)---------------------------
     *, 122., .60000, 120., 220.,   0.,   0.,   0.
     *, 122.,1.00000, 130.,-230.,   0.,   0.,   0./
c    *---------------------------------------------
c    *      weak baryon decays
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=1160,1169)/
c    *-----------lambda(1116)----------------------
     * 2130.,0.64200,1120.,-120.,   0.,   0.,   0.
     *,2130.,1.00000,1220., 110.,   0.,   0.,   0.
c    *-----------sigma+(1180)----------------------
     *,1130.,0.51580,1120., 110.,   0.,   0.,   0.
     *,1130.,1.00000,1220., 120.,   0.,   0.,   0.
c    *-----------sigma-(1180)----------------------
     *,2230.,1.00000,1220.,-120.,   0.,   0.,   0.
c    *---------kaskade-(1360)----------------------
     *,2330.,1.00000,2130.,-120.,   0.,   0.,   0.
c    *---------kaskade0(1360)----------------------
     *,1330.,1.00000,2130., 110.,   0.,   0.,   0.
c    *---------omega-(1680)------------------------
     *,3331.,0.68000,2130.,-130.,   0.,   0.,   0.
     *,3331.,0.82000,1330.,-120.,   0.,   0.,   0.
     *,3331.,1.00000,2330., 110.,   0.,   0.,   0./
c    *---------------------------------------------
c    *      weak meson decays
c    *---------------------------------------------
      data ((dectab(i,j),i=1,7),j=1170,1171)/
c    *-----------k0s(975)--------------------------
     *   20., .68610, 120.,-120.,   0.,   0.,   0.
     *,  20.,1.00000, 110., 110.,   0.,   0.,   0./
c    *---------------------------------------------
      end
 
c-----------------------------------------------------------------------
      subroutine jestpr(ic1,ic2,am,ier)
c     processes strings
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      logical revers
      parameter (nflav=6)
      integer ic1(2),ic2(2),jc(nflav,2),jc1(nflav,2),jc2(nflav,2)
      parameter (nsi=6,nsix=100)
      real stro(nsi,nsix+1)
      common/ckol/kol /clead/lead,cooav3,cooav4
      parameter (nha=nsix+2,mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      common/cpzstr/pzstrl,estrl,istrl,isea
 
      ier=0
 
c     pmax
c     ----
      call iddeco(ic1,jc)
      nq=0
      do 7 nf=1,nflav
7     nq=nq+jc(nf,1)-jc(nf,2)
      am1min=piom
      if(iabs(nq).ge.2)am1min=prom
      call iddeco(ic2,jc)
      nq=0
      do 8 nf=1,nflav
8     nq=nq+jc(nf,1)-jc(nf,2)
      am2min=piom
      if(iabs(nq).ge.2)am2min=prom
      if(am.le.am1min+am2min)then
      pmax=am/2
      else
      pmax=utpcm(am,am1min,am2min)
      endif
 
c     hastpr
c     ------
      iproj=1          !this leads to iorstr()=-1
      itarg=1          !
      kproja(2,1)=1    !
      ktarga(2,1)=1    !
      kproja(3,1)=1    !
      ktarga(3,1)=1    !
      kol=1            !
      do 1 i=1,4         !this leads to xorstr()=0
1     coord(i,1)=0       !
      call iddeco(ic1,jc1)
      npa1=0
      do 2 n=1,nflav
      do 2 i=1,2
2     npa1=npa1+jc1(n,i)
      call iddeco(ic2,jc2)
      npa2=0
      do 3 n=1,nflav
      do 3 i=1,2
3     npa2=npa2+jc2(n,i)
      ii=1
      if(npa2.gt.npa1)ii=2
      lead=0
      if(npa1.gt.1.or.npa2.gt.1)lead=1
 
      nstr0=nstr
17    nstr=nstr0
 
      do 15 n=1,nsi
      do 15 i=1,3
15    stro(n,i)=0
      stro(3,ii)=am/2
      stro(4,ii)=am/2
      stro(5,ii)=ic1(1)
      stro(6,ii)=ic1(2)
      stro(3,3-ii)=-am/2
      stro(4,3-ii)=am/2
      stro(5,3-ii)=ic2(1)
      stro(6,3-ii)=ic2(2)
 
      pzstrl=stro(3,1)
      estrl=stro(4,1)
      pzstrl=pzstrl+stro(3,2)
      estrl=estrl+stro(4,2)
      istrl=0
 
      call utpage
      if(ish.ge.91)then
      write(ifch,110)('-',l=1,79),ipage,('-',l=1,79)
110   format(1x,79a1/1x,i5,'.page            '
     *,'string generation'/1x,79a1/)
      write(ifch,105)(stro(i,1),i=1,4),(nint(stro(i,1)),i=5,6)
      write(ifch,104)(stro(i,2),i=1,4),(nint(stro(i,2)),i=5,6)
105   format(' str: ',4f11.5,2i8)
104   format('      ',4f11.5,2i8/)
      endif
 
      isea=1
      isplt=0
14    call hastpr(stro,isplt)
      if(isplt.eq.-1)goto9001
      if(isplt.eq.-3)goto9001
      if(isplt.eq.-4)
     *call utstop('jestpr: isplt=-4                        ')
      if(isplt.eq.-5)goto17
      if(isplt.gt.0)goto14
 
      pmxevt=pmax
      egyevt=am
 
      return
 
9001  ier=1
      return
 
      end
 
c-----------------------------------------------------------------------
      subroutine jetgen(ier)
c     generates strings
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      common/partnr/summas,iptnr,ic4,px4,py4,js4,pex,pey,pezet,pe0,nps
     */cntevm/ntevm
      common /clep/icinpu,idscat
      integer ic1(2),ic2(2),ic4(2)
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      parameter (mamx=250)
      real targa(nsi,nha,mamx)
      integer icstri(4,18)
 
      data ((icstri(i,j),i=1,4),j=1,8)/
     *100000,000000,110000,000000,
     *010000,000000,200000,000000,
     *100000,000000,210000,100000,
     *010000,000000,210000,010000,
     *001000,000000,210000,001000,
     *000000,100000,310000,000000,
     *000000,010000,220000,000000,
     *000000,001000,211000,000000/
      data ((icstri(i,j),i=1,4),j=11,18)/
     *100000,000000,020000,000000,
     *010000,000000,110000,000000,
     *100000,000000,120000,100000,
     *010000,000000,120000,010000,
     *001000,000000,120000,001000,
     *000000,100000,220000,000000,
     *000000,010000,130000,000000,
     *000000,001000,121000,000000/
 
      ish0=ish
      if(ishsub/100.eq.2)ish=mod(ishsub,100)
      if(ish.ge.93)then
      write(ifch,*)('-',l=1,79)
      write(ifch,*)'string generation. sr jetgen.'
      endif
 
      ier=0
 
      nevt=1
      nstr=0
      nptl=0
 
           if(icinpu.ge.1)then
3     call lepexp(xbj,qsq)
      rnu=qsq/(2.* prom * xbj)
      wsq=prom**2 + 2.* prom * rnu - qsq
      amst=sqrt(wsq)
      elepto=elepti-rnu
      cosang=1-qsq/2/elepti/elepto
      if(ish.ge.93)write(ifch,*)'elepti,elepto,rnu: ',elepti,elepto,rnu
      if(ish.ge.93)write(ifch,*)'xbj,qsq,cosang=1-qsq/2/elepti/elepto: '
     *,xbj,qsq,cosang
      if(rnu.gt.elepti)then
      if(ish.ge.93)write(ifch,*)'*****  q0 too large: ',rnu
      goto3
      endif
      if(cosang.gt.1.and.cosang.lt.1+1e-3)then
      cosang=1.
      elseif(cosang.lt.-1..or.cosang.gt.1.)then
      cosano=cosang
      if(cosang.ge.0.)cosang=1.
      if(cosang.lt.0.)cosang=-1.
      call utmsg('jetgen')
      write(ifch,*)'*****  cosang out of range'
      write(ifch,*)'elepti,elepto,rnu: ',elepti,elepto,rnu
      write(ifch,*)'xbj,qsq,cosang=1-qsq/2/elepti/elepto: '
     *,xbj,qsq,cosano
      write(ifch,*)'cosang_new: ',cosang
      call utmsgf
      endif
      angmue=acos(cosang)
      call leptar(xbj,qsq,matarg,latarg,idscat)
           endif
 
      if(matarg.gt.0)then
      call nucoge
      nptl=0
      call nucini('str',targa,latarg,matarg,-1)
      call nucstr(ier)
      if(ier.eq.1)goto99999
      endif
 
           if(icinpu.eq.0)then
      r=rangen()
      ps=0
      do 1 k=1,99
      ps=ps+prob(k)
      if(r.le.ps)goto2
1     continue
      call utstop('jetgen: no k found                      ')
2     ic1(1)=icfor(k,1)
      ic1(2)=icfor(k,2)
      ic2(1)=icbac(k,1)
      ic2(2)=icbac(k,2)
      amst=engy
           elseif(icinpu.gt.0)then
      call lepstr(idscat,xbj,qsq,ids)
      ic1(1)=icstri(1,ids)
      ic1(2)=icstri(2,ids)
      ic2(1)=icstri(3,ids)
      ic2(2)=icstri(4,ids)
           endif
      call jestpr(ic1,ic2,amst,ier)
 
99999 ish=ish0
      return
      end

c-----------------------------------------------------------------------
      subroutine jfrade(ier)
c-----------------------------------------------------------------------
c  performs string fragmentation/decay and fin. state interactions
c-----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer iacptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      real xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt

      common/partnr/summas,iptnr,ic4,px4,py4,js4,pex,pey,pezet,pe0,nps
     */didib/ndidib /cnfr/nrfra
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc
      parameter(nprbms=20)
      common /cprbms/prbms(nprbms) /cnce/nces,ncolex /ckol/kol 
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      parameter (mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      double precision ppfus(5)
      
      real psum(5),u(5),pei(5),pej(5)
      parameter (nflav=6)
      integer jci(nflav,2),jcj(nflav,2),jc(nflav,2)
     *,ici(2),icj(2),ic(2),ic4(2),jcx(nflav,2),jcdu(nflav,2)
      character dash*1,chadec*5

c     initialization
c     --------------
      if(ish.eq.14)call uttima('*** jfrade *** ')
      if(ish.eq.13)call uttima('*** jfrade *** ')
      ier=0
      iret=0
      ish0=ish
      if(ishsub/100.eq.16)ish=mod(ishsub,100)
      if(nevt.ne.1.or.ifrade.eq.0)goto1000
      nptlpt=iabs(maproj)+iabs(matarg)
      if(ichoic.eq.1.or.ichoic.eq.4)then
      ycmmax=alog(egyevt/.4)
      etapro=ycmmax*.6667
      etatar=-etapro
      else
      etapro=(ypjtl-yhaha)*.6667
      etatar=-yhaha*.6667
      endif
      detap=etapro
      detat=etatar
      tpro=dcosh(detap)
      zpro=dsinh(detap)
      ttar=dcosh(detat)
      ztar=dsinh(detat)
 
c     write
c     -----
      call utpage
      if(ish.ge.91)write(ifch,110)('-',l=1,79),ipage,('-',l=1,79)
110   format(1x,79a1/1x,i5,'.page            '
     *,'strings before rescaling'/1x,79a1/)
      do 9 j=1,nstr
      if(ish.ge.91)write(ifch,109)j,(icstr(k,j)/100,k=1,4)
     *,sqrt(pstr(1,j)**2+pstr(2,j)**2),pstr(3,j),pstr(5,j)
     *,irlstr(j)
9     continue
109   format(' /cstr/',i4,3x,4i5,2x,3(1x,e9.3),i4)
 
c     string rescaling
c     ----------------
      if(ichoic.le.2.and.nstr.gt.1.and.irescl.eq.1)then
      psum(1)=0
      psum(2)=0
      if(ichoic.eq.1)then
      psum(3)=0
      psum(4)=egyevt
      else
      psum(3)=(npjevt-ntgevt)*pnllx
      psum(4)=npjevt*sqrt(amproj**2+pnllx**2)
     *+ntgevt*sqrt(amtarg**2+pnllx**2)
      endif
      psum(5)=sqrt(psum(4)**2-psum(3)**2)
      call hrescl(1,nstr,psum,ifail)
      if(ifail.ne.0)goto1001
      endif
 
c     write
c     -----
      call utpage
      if(ish.ge.91)write(ifch,111)('-',l=1,79),ipage,('-',l=1,79)
111   format(/1x,79a1/1x,i5,'.page            '
     *,'strings after rescaling'/1x,79a1/)
      do 10 j=1,nstr
      if(ish.ge.91)write(ifch,109)j,(icstr(k,j)/100,k=1,4)
     *,sqrt(pstr(1,j)**2+pstr(2,j)**2),pstr(3,j),pstr(5,j)
     *,irlstr(j)
10    continue
 
c     fragmentation
c     -------------
      if(ishsub/100.eq.3)ish=mod(ishsub,100)
      if(ish.eq.13)call uttima('               ')
      do 3 j=1,nstr
      call utpage
      if(ish.ge.91)write(ifch,102)('-',l=1,79),ipage,j
     *,(icstr(k,j),k=1,4),sqrt(pstr(1,j)**2+pstr(2,j)**2
     *+pstr(3,j)**2),pstr(4,j),pstr(5,j),('-',l=1,79)
102   format(/1x,79a1
     */1x,i5,'.page  str:',i3,4i7,1x,3(e8.2,1x)/1x,79a1/)
      call jamfra(j,newevt)
      if(newevt.eq.1)goto10011
3     continue
      if(ish.eq.13)call uttima('fragmentation  ')
 
c     print /cptl/
c     ------------
      if(ishsub/100.eq.9)ish=mod(ishsub,100)
      call utpage
      if(ish.ge.91)then
      write(ifch,104)('-',l=1,79),ipage,('-',l=1,79)
104   format(/1x,79a1/1x,i5,'.page            '
     *,'ptls after fragmentation'/1x,79a1/)
      do 18 n=1,nptl
      if(n.gt.nptlpt.or.istptl(n).ne.0)
     *write(ifch,116)iorptl(n),jorptl(n),n,ifrptl(1,n),ifrptl(2,n)
     *,idptl(n),pptl(3,n),pptl(4,n),pptl(5,n),istptl(n)
18    continue
      endif
 
c     initial decay
c     -------------
      if(ishsub/100.eq.10)ish=mod(ishsub,100)
      call utpage
      if(ish.eq.13)call uttima('               ')
      iacn=1
      if(radiac.gt.0.)then
      ttaus=taumin
      else
      ttaus=ainfin
      iacn=0
      endif
      if(ish.ge.91)then
      write(ifch,119)('-',l=1,79),ipage,sngl(ttaus),('-',l=1,79)
119   format(/1x,79a1/1x,i5,'.page            '
     *,'decay before tau =',e10.3/1x,79a1/)
      endif
      np1=1
21    np2=nptl
      do 5 i=np1,np2
      zfi=xorptl(3,i)+(tivptl(2,i)-xorptl(4,i))*pptl(3,i)/pptl(4,i)
      call uttaus(zfi,tzfi,szfi)
      tme=tzfi
      if(tivptl(2,i).le.tme.or.iacn.eq.0)call jdeca(i,iret)
      if(iret.eq.1)goto1001
5     continue
      np1=np2+1
      if(np1.le.nptl) goto21
      if(ish.eq.13)call uttima('decay ini      ')

c     interaction and decay 
c     ---------------------
      if(iacn.eq.0)goto5000
      call utpage
           if(ish.ge.91)then
      write(ifch,118)('-',l=1,79),ipage,('-',l=1,79)
118   format(/1x,79a1/1x,i5,'.page            '
     *,'interactions and decay'/1x,79a1/)
           endif
      do 23 nt=1,numtau
      ntc=nt
      dtaus=1./(numtau-1.)*(taumax-taumin)
      ttaus=taumin+(nt-1.)*dtaus
      if(iojint.eq.1)call jinta1
      if(iojint.eq.2)call jinta2
      ttaus=ttaus+dtaus
      np1=1
36    np2=nptl
      do 37 ip=np1,np2
      zfi=xorptl(3,ip)+(tivptl(2,ip)-xorptl(4,ip))*pptl(3,ip)/pptl(4,ip)
      call uttaus(zfi,tzfi,szfi)
      tme=tzfi
      if(tivptl(2,ip).le.tme)call jdeca(ip,iret)
      if(iret.eq.1)goto1001
37    continue
      np1=np2+1
      if(np1.le.nptl)then
      do 4 ip=np1,nptl
4     iaaptl(ip)=1        
      goto36
      endif
      if(ish.eq.13)call uttima('decay          ')
23    continue

c     final decay
c     -----------
      call utpage
      if(ish.ge.91)then
      write(ifch,120)('-',l=1,79),ipage,('-',l=1,79)
120   format(/1x,79a1/1x,i5,'.page            '
     *,'final decay'/1x,79a1/)
      endif
      nbef=nptl
      n0bef=nptl0
      naft=nptl
      n0aft=nptl0
           if(nclean.gt.0)then
      call utclea(nptl0)
      naft=nptl
      n0aft=nptl0
           endif
      ishnpt=ish
      if(ishsub/100.eq.19)ish=mod(ishsub,100)
           if(ish.eq.22)then
      write(ifmt,131)nbef,naft
131   format(1x,'bef fin de: nbef=',i8,4x,'naft=',i8)
           endif
      ish=ishnpt
      np1=1
41    np2=nptl
      do 42 ip=np1,np2
      call jdeca(ip,iret)
      if(iret.eq.1)goto1001
42    continue
      np1=np2+1
      if(np1.le.nptl)then
      do 6 ip=np1,nptl
6     iaaptl(ip)=1
      goto41
      endif
      ishnpt=ish
      if(ishsub/100.eq.19)ish=mod(ishsub,100)
           if(ish.eq.22)then
      write(ifmt,132)nptl      
132   format(1x,'aft fin de: nptl=',i8)                  
           endif
      ish=ishnpt
      if(ish.eq.13)call uttima('decay fin      ')
 
c     print /cptl/
c     ------------
      call utpage
      if(ishsub/100.eq.11)ish=mod(ishsub,100)
      if(ish.ge.91)then
      write(ifch,117)('-',l=1,79),ipage,('-',l=1,79)
117   format(/1x,79a1/1x,i5,'.page            '
     *,'ptls after performing interactions'/1x,79a1/
     *1x,' ior',' jor',4x,'   n',4x,'ifr1ifr2',10x,'id',3x,
     *5x,'ptr',7x,'pz',4x,'mass ','ist'/)
      do 34 n=1,nptl
      if(n.gt.nptlpt.or.istptl(n).ne.0)
     *write(ifch,116)iorptl(n),jorptl(n),n,ifrptl(1,n),ifrptl(2,n)
     *,idptl(n),sqrt(pptl(1,n)**2+pptl(2,n)**2),pptl(3,n),pptl(5,n)
     *,istptl(n),iaaptl(n)
34    continue
116   format(1x,i4,i4,4x,i4,4x,i4,i4,i12,3x,3(e8.2,1x),i2,i2)
      endif

c     trafo -> lab cm
c     ---------------
5000  continue
      if(labsys.eq.1)then
      do 7 i=1,nptl
      amt=sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2)
      pz=pptl(3,i)
      e=pptl(4,i)
      yi=sign(1.,pz)*alog((e+abs(pz))/amt)
      y=yi+yhaha
      pptl(3,i)=amt*sinh(y)
7     pptl(4,i)=amt*cosh(y)
      endif
 
c     finish
c     ------

      call utpage
      if(ish.ge.91)write(ifch,113)('-',l=1,79),ipage,('-',l=1,79)
113   format(/1x,79a1/1x,i5,'.page            '
     *,'particle production finished'/1x,79a1/)
      if(ish.eq.13)call uttima('               ')
      goto1000
 
10011 continue
      if(ish.ge.90)then
      call utmsg('jfrade')
      write(ifch,*)'*****  fragmentation not possible'
      write(ifch,112)j,(icstr(k,j)/100,k=1,4)
     *,(pstr(l,j),l=1,5)
112   format( ' str:',3x,i3,3x,3x,4i5,5f7.2,f7.2)
      call utmsgf
      endif
      ier=1
      goto1000
 
1001  continue
      ier=1
      goto1000
 
1000  continue
           if(ier.eq.1)then
      call utmsg('jfrade')
      write(ifmt,*)'*****  redo event',nrevt+1
      write(ifch,*)'*****  redo event',nrevt+1
      call utmsgf
           endif
      if(ish.eq.14)call uttima('    jfrade f   ')
      ish=ish0

      return
      end
 
c----------------------------------------------------------------------
      subroutine jinta1
c----------------------------------------------------------------------
c  secondary interactions
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer iacptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      parameter (mxifr=mxptl)
      common/cifrij/ifrij(mxifr)
      real xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt

      common/partnr/summas,iptnr,ic4,px4,py4,js4,pex,pey,pezet,pe0,nps
     */didib/ndidib /cnfr/nrfra
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc
      parameter(nprbms=20)
      common /cprbms/prbms(nprbms) /cnce/nces,ncolex /ckol/kol 
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      parameter (mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      common/cjint/nptl0,rnup,rnut,velp,velt,bx,by,xaver(4)

      double precision ppfus(5)
      real psum(5),u(5),pei(5),pej(5)
      parameter (nflav=6)
      integer jci(nflav,2),jcj(nflav,2),jc(nflav,2)
     *,ici(2),icj(2),ic(2),ic4(2),jcx(nflav,2),jcdu(nflav,2)
      integer ifr(mxifr)
      real proj(nsi,nha),targ(nsi,nha)
      double precision arp(5),arm(5),arq(5),rot(3),boo(5)
      character chadec*5

c     initialization for ntc=1
c     ------------------------
           if(ntc.eq.1)then
         if(overlp.ge.0..and.maproj.ne.0.and.matarg.ne.0)then
      rnup=1.19*maproj**(1./3.)-1.61*maproj**(-1./3.)
      rnut=1.19*matarg**(1./3.)-1.61*matarg**(-1./3.)
      if(maproj.eq.1)rnup=0
      if(matarg.eq.1)rnut=0
      rnup=rnup+overlp
      rnut=rnut+overlp
      velp=pptl(3,1)/pptl(4,1)
      velt=pptl(3,iabs(maproj)+1)/pptl(4,iabs(maproj)+1)
      bx=cos(phievt)*bimevt
      by=sin(phievt)*bimevt
         endif
      nptl0=nptl
      do 51 i=1,nptl0
      iaaptl(i)=1
      call idquac(i,nqi,ndummy,ndummy,jci)
      if(nqi.ne.0)radptl(i)=radiac
      if(nqi.eq.0)radptl(i)=radias
      amiptl(i)=prom+amsiac
      if(nqi.eq.0)amiptl(i)=piom+amsiac
51    continue
           endif

c     initialization for each ntc
c     ---------------------------
      nt=ntc
      call utpage
      if(ish.eq.13)call uttima('               ')
      nbef=nptl
      n0bef=nptl0
      naft=nptl
      n0aft=nptl0
           if(nclean.gt.0.and.mod(nt-1,nclean).eq.0)then
      call utclea(nptl0)
      naft=nptl
      n0aft=nptl0
           endif
      taus=ttaus
      nptlpt=iabs(maproj)+iabs(matarg)
           if(nptl.gt.nptl0)then
      do 52 i=nptl0+1,nptl
      call idquac(i,nqi,ndummy,ndummy,jci)
      if(nqi.ne.0)radptl(i)=radiac
      if(nqi.eq.0)radptl(i)=radias
      amiptl(i)=prom+amsiac
      if(nqi.eq.0)amiptl(i)=piom+amsiac
52    continue
           endif
      do 44 i=1,nptl
      iacptl(i)=1
      if(iclptl(i).eq.0.and.pptl(5,i).gt.amiptl(i))iacptl(i)=0
      if(iorptl(i).eq.-1)iacptl(i)=0
      call uttain(i,x,y,z,t,n,0)
      if(n.eq.1.or.n.eq.2.or.n.eq.9)iacptl(i)=0
      if(nclean.gt.0
     *.and.mod(nt,nclean).eq.0.and.n.eq.2)istptl(i)=2
      if(istptl(i).eq.2)iacptl(i)=0
      if(iaaptl(i).eq.0)iacptl(i)=0
      if(iabs(idptl(i)).lt.100.and.iabs(idptl(i)).ne.20)iacptl(i)=0
      if(idptl(i).eq.441.and.jpsifi.eq.0)iacptl(i)=0
      xptl(i)=x
      yptl(i)=y
      zptl(i)=z
      tptl(i)=t
      if(overlp.ge.0..and.maproj.ne.0.and.matarg.ne.0)then
      if((x-bx/2)**2+(y-by/2)**2+(z-velp*taus)**2.lt.rnup**2.and.
     *(x+bx/2)**2+(y+by/2)**2+(z-velt*taus)**2.lt.rnut**2)iacptl(i)=0
      endif
44    continue
      nptl0=nptl
      i=0
 
c     i loop --> 24
c     -------------
9999  i=i+1
      if(iacptl(i).ne.1)goto24
      j0=nptlpt+1
      if(i.gt.nptlpt)j0=i+1
      if(i.gt.nptl0)j0=1
 
c     j loop --> 25
c     -------------
      do 25 j=j0,nptl
      if(iacptl(j).ne.1)goto25
      radsqr=(radptl(i)+radptl(j))**2
      if((zptl(i)-zptl(j))**2.gt.radsqr)goto25
      if((yptl(i)-yptl(j))**2.gt.radsqr)goto25
      if((xptl(i)-xptl(j))**2.gt.radsqr)goto25
      if(i.eq.j)goto25
      if(iorptl(i).gt.0.and.iorptl(j).eq.iorptl(i))goto25
      if(iorptl(i).eq.j)goto25
      if(iorptl(j).eq.i)goto25
      pde=(pptl(3,i)+pptl(3,j))/(pptl(4,i)+pptl(4,j))
      gam2i=1-pde**2
      if(gam2i.eq.0.)goto25
      if((zptl(i)-zptl(j)-(tptl(i)-tptl(j))*pde)**2
     *.gt.radsqr*gam2i)goto25
      if((xptl(i)-xptl(j))**2+(yptl(i)-yptl(j))**2+
     *1/gam2i*(zptl(i)-zptl(j)-(tptl(i)-tptl(j))*pde)**2
     *.gt.radsqr)goto25
      call jintcc(i,j,iret)  !check children
      if(iret.eq.1)goto25                   

      iactn=0
      nptl00=nptl
      nstr00=nstr
      xaver(1)=(xptl(i)+xptl(j))/2    
      xaver(2)=(yptl(i)+yptl(j))/2    
      xaver(3)=(zptl(i)+zptl(j))/2    
      xaver(4)=(tptl(i)+tptl(j))/2    

      call jintfs(i,j,nqifus,jc,amim,iret)    !fusion (1)        
      if(iret.eq.1)goto25                 

      call jintce(i,j,amim,iactn,iret)        !colour exchange (1)
      if(iret.eq.25)goto25             

      tivptl(2,i)=tptl(i)
      tivptl(2,j)=tptl(j)
      istptl(i)=1
      istptl(j)=1
      iacptl(i)=0
      iacptl(j)=0

      call jintch(i,j,kmax)  
      do k=1,kmax
      n=ifrij(k)            ! n: i,j child
      istptl(n)=2
      iacptl(n)=0
      enddo
      call jintpa(i,j,kmax) ! also: iaaptl(n)=0 
      do k=1,kmax
      n=ifrij(k)            ! n: i,j parent
      iacptl(n)=0
      enddo

      if(iactn.eq.1)then                     !colour exchange (2)
      do n=nptl00+1,nptl
      iaaptl(n)=1
      call idquac(n,nqi,ndummy,ndummy,jcdu)
      if(nqi.eq.0)then
      radptl(n)=radias
      amiptl(n)=piom+amsiac
      else
      radptl(n)=radiac
      amiptl(n)=prom+amsiac
      endif
      iacptl(n)=1
      if(pptl(5,n).gt.amiptl(n))iacptl(n)=0
      call uttain(n,x,y,z,t,k,0)
      if(k.eq.1.or.k.eq.2.or.k.eq.9)iacptl(n)=0
      if(iabs(idptl(n)).lt.100.and.iabs(idptl(n)).ne.20)iacptl(n)=0
      xptl(n)=x
      yptl(n)=y
      zptl(n)=z
      tptl(n)=t
      if(overlp.ge.0..and.maproj.ne.0.and.matarg.ne.0)then
      if((x-bx/2)**2+(y-by/2)**2+(z-velp*taus)**2.lt.rnup**2.and.
     *(x+bx/2)**2+(y+by/2)**2+(z-velt*taus)**2.lt.rnut**2)iacptl(n)=0
      endif
      enddo
      goto24                   
      endif

      call jintel(i,j,amim,iactn)            !elastic scattering
      if(iactn.eq.2)then
      do n=nptl00+1,nptl
      if(n.eq.nptl00+1)ij=i
      if(n.eq.nptl00+2)ij=j
      radptl(n)=radptl(ij)
      amiptl(n)=amiptl(ij)
      xptl(n)=xaver(1)
      yptl(n)=xaver(2)
      zptl(n)=xaver(3)                    
      tptl(n)=xaver(4)
      iacptl(n)=1
      iaaptl(n)=1
      enddo
      goto24
      endif

      call jintfu(i,j,jc,iactn)               !fusion (2)                     
      do n=nptl00+1,nptl
      if(nqifus.eq.0)then
      radptl(n)=radias
      amiptl(n)=piom+amsiac
      else
      radptl(n)=radiac
      amiptl(n)=prom+amsiac
      endif
      xptl(n)=xaver(1)
      yptl(n)=xaver(2)
      zptl(n)=xaver(3)                    
      tptl(n)=xaver(4)
      iacptl(n)=1
      iaaptl(n)=1
      enddo
      goto24

25    continue
 
24    continue
      if(i.lt.nptl-1)goto9999

      ishnpt=ish
      if(ishsub/100.eq.19)ish=mod(ishsub,100)
      if(ish.eq.22)then
      write(ifmt,131)nt,nbef,naft,nptl
131   format(1x,'nt=',i5,4x,'nbef=',i8,4x,'naft=',i8,4x,'nptl=',i8)
      endif
      ish=ishnpt
      if(ish.eq.13)call uttima('interactions   ')
      return
      end

c----------------------------------------------------------------------
      subroutine jinta2
c----------------------------------------------------------------------
c  secondary interactions: percolation--cluster model
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer iacptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      parameter (mxifr=mxptl)
      common/cifrij/ifrij(mxifr)
      real xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
     *,optl(mxptl),uptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt

      common/partnr/summas,iptnr,ic4,px4,py4,js4,pex,pey,pezet,pe0,nps
     */didib/ndidib /cnfr/nrfra
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc
      parameter(nprbms=20)
      common /cprbms/prbms(nprbms) /cnce/nces,ncolex /ckol/kol 
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      parameter (mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      common/cjint/nptl0,rnup,rnut,velp,velt,bx,by,xaver(4)

      double precision ppfus(5)
      real psum(5),u(5),pei(5),pej(5)
      parameter (nflav=6)
      integer jci(nflav,2),jcj(nflav,2),jc(nflav,2)
     *,ici(2),icj(2),ic(2),ic4(2),jcx(nflav,2),jcdu(nflav,2)
      integer ifr(mxifr)
      real proj(nsi,nha),targ(nsi,nha)
      double precision arp(5),arm(5),arq(5),rot(3),boo(5)
      character chadec*5

c     initialization for ntc=1
c     ------------------------
           if(ntc.eq.1)then
         if(overlp.ge.0..and.maproj.ne.0.and.matarg.ne.0)then
      rnup=1.19*maproj**(1./3.)-1.61*maproj**(-1./3.)
      rnut=1.19*matarg**(1./3.)-1.61*matarg**(-1./3.)
      if(maproj.eq.1)rnup=0
      if(matarg.eq.1)rnut=0
      rnup=rnup+overlp
      rnut=rnut+overlp
      velp=pptl(3,1)/pptl(4,1)
      velt=pptl(3,iabs(maproj)+1)/pptl(4,iabs(maproj)+1)
      bx=cos(phievt)*bimevt
      by=sin(phievt)*bimevt
         endif
      nptl0=nptl
      do 51 i=1,nptl0
      iaaptl(i)=1
      call idquac(i,nqi,ndummy,ndummy,jci)
      if(nqi.eq.0)then
      radptl(i)=radias
      desptl(i)=radias
      dezptl(i)=0          
      amiptl(i)=piom+amsiac
      else
      radptl(i)=radiac
      desptl(i)=radiac
      dezptl(i)=0          
      amiptl(i)=prom+amsiac
      endif
51    continue
           endif

c     initialization for each ntc
c     ---------------------------
      nt=ntc
      call utpage
      if(ish.eq.13)call uttima('               ')
      nbef=nptl
      n0bef=nptl0
      naft=nptl
      n0aft=nptl0
           if(nclean.gt.0.and.mod(nt-1,nclean).eq.0)then
      call utclea(nptl0)
      naft=nptl
      n0aft=nptl0
           endif
      taus=ttaus
      nptlpt=iabs(maproj)+iabs(matarg)
           if(nptl.gt.nptl0)then
      do 52 i=nptl0+1,nptl
      call idquac(i,nqi,ndummy,ndummy,jci)
      if(nqi.eq.0)then
      radptl(i)=radias
      desptl(i)=radias
      dezptl(i)=0          
      amiptl(i)=piom+amsiac
      else
      radptl(i)=radiac
      desptl(i)=radiac
      dezptl(i)=0          
      amiptl(i)=prom+amsiac
      endif
52    continue
           endif
      do 44 i=1,nptl
      iacptl(i)=1
      if(iclptl(i).eq.0.and.pptl(5,i).gt.amiptl(i))iacptl(i)=0
      if(iorptl(i).eq.-1)iacptl(i)=0
      call uttain(i,x,y,z,t,n,0)
      if(n.eq.1.or.n.eq.2.or.n.eq.9)iacptl(i)=0
      if(nclean.gt.0
     *.and.mod(nt,nclean).eq.0.and.n.eq.2)istptl(i)=2
      if(istptl(i).eq.2)iacptl(i)=0
      if(iaaptl(i).eq.0)iacptl(i)=0
      if(iabs(idptl(i)).lt.100.and.iabs(idptl(i)).ne.20)iacptl(i)=0
      if(idptl(i).eq.441.and.jpsifi.eq.0)iacptl(i)=0
      xptl(i)=x
      yptl(i)=y
      zptl(i)=z
      tptl(i)=t
      call uttaus(z,tz,sz)
      optl(i)=sz+desptl(i)
      uptl(i)=sz-desptl(i)
      if(overlp.ge.0..and.maproj.ne.0.and.matarg.ne.0)then
      if((x-bx/2)**2+(y-by/2)**2+(z-velp*taus)**2.lt.rnup**2.and.
     *(x+bx/2)**2+(y+by/2)**2+(z-velt*taus)**2.lt.rnut**2)iacptl(i)=0
      endif
44    continue
      nptl0=nptl
      i=0
 
c     i loop --> 24
c     -------------
9999  i=i+1
      if(iacptl(i).ne.1)goto24
      j0=nptlpt+1
      if(i.gt.nptlpt)j0=i+1
      if(i.gt.nptl0)j0=1
 
c     j loop --> 25
c     -------------
      do 25 j=j0,nptl
      if(iacptl(j).ne.1)goto25
      if(optl(i).lt.uptl(j))goto25
      if(optl(j).lt.uptl(i))goto25
      radsqr=(radptl(i)+radptl(j))**2
      if((yptl(i)-yptl(j))**2.gt.radsqr)goto25
      if((xptl(i)-xptl(j))**2.gt.radsqr)goto25
      if(i.eq.j)goto25
      if(iorptl(i).gt.0.and.iorptl(j).eq.iorptl(i))goto25
      if(iorptl(i).eq.j)goto25
      if(iorptl(j).eq.i)goto25
      if((xptl(i)-xptl(j))**2+(yptl(i)-yptl(j))**2.gt.radsqr)goto25
      call jintcc(i,j,iret)  !check children
      if(iret.eq.1)goto25                   

      iactn=0
      nptl00=nptl
      nstr00=nstr
      xaver(1)=(xptl(i)+xptl(j))/2    
      xaver(2)=(yptl(i)+yptl(j))/2    
      xaver(3)=(zptl(i)+zptl(j))/2    
      xaver(4)=(tptl(i)+tptl(j))/2    

      call jintfs(i,j,nqifus,jc,amim,iret)    !fusion (1)        
      if(iret.eq.1)goto25                 

      call jintce(i,j,amim,iactn,iret)        !colour exchange (1)
      if(iret.eq.25)goto25             

      tivptl(2,i)=tptl(i)
      tivptl(2,j)=tptl(j)
      istptl(i)=1
      istptl(j)=1
      iacptl(i)=0
      iacptl(j)=0

      call jintch(i,j,kmax)  
      do k=1,kmax
      n=ifrij(k)            ! n: i,j child
      istptl(n)=2
      iacptl(n)=0
      enddo
      call jintpa(i,j,kmax) ! also: iaaptl(n)=0 
      do k=1,kmax
      n=ifrij(k)            ! n: i,j parent
      iacptl(n)=0
      enddo

      if(iactn.eq.1)then                     !colour exchange (2)
      do n=nptl00+1,nptl
      iaaptl(n)=1
      call idquac(n,nqi,ndummy,ndummy,jcdu)
      if(nqi.eq.0)then
      radptl(n)=radias
      desptl(n)=radias
      dezptl(n)=0          
      amiptl(n)=piom+amsiac
      else
      radptl(n)=radiac
      desptl(n)=radiac
      dezptl(n)=0          
      amiptl(n)=prom+amsiac
      endif
      iacptl(n)=1
      if(pptl(5,n).gt.amiptl(n))iacptl(n)=0
      call uttain(n,x,y,z,t,k,0)
      if(k.eq.1.or.k.eq.2.or.k.eq.9)iacptl(n)=0
      if(iabs(idptl(n)).lt.100.and.iabs(idptl(n)).ne.20)iacptl(n)=0
      xptl(n)=x
      yptl(n)=y
      zptl(n)=z
      tptl(n)=t
      call uttaus(z,tz,sz)
      optl(n)=sz+desptl(n)
      uptl(n)=sz-desptl(n)
      if(overlp.ge.0..and.maproj.ne.0.and.matarg.ne.0)then
      if((x-bx/2)**2+(y-by/2)**2+(z-velp*taus)**2.lt.rnup**2.and.
     *(x+bx/2)**2+(y+by/2)**2+(z-velt*taus)**2.lt.rnut**2)iacptl(n)=0
      endif
      enddo
      goto24                   
      endif

      call jintel(i,j,amim,iactn)            !elastic scattering
      if(iactn.eq.2)then
      do n=nptl00+1,nptl
      if(n.eq.nptl00+1)ij=i
      if(n.eq.nptl00+2)ij=j
      radptl(n)=radptl(ij)
      desptl(n)=desptl(ij)
      dezptl(n)=dezptl(ij) 
      amiptl(n)=amiptl(ij)
      xptl(n)=xaver(1)
      yptl(n)=xaver(2)
      zptl(n)=xaver(3)                    
      tptl(n)=xaver(4)
      z=zptl(n)
      call uttaus(z,tz,sz)
      optl(n)=sz+desptl(n)
      uptl(n)=sz-desptl(n)
      iacptl(n)=1
      iaaptl(n)=1
      enddo
      goto24
      endif

      call jintfu(i,j,jc,iactn)               !fusion (2)                     
      do n=nptl00+1,nptl
      if(nqifus.eq.0)then
      radptl(n)=radias
      desptl(n)=radias
      dezptl(n)=0          
      amiptl(n)=piom+amsiac
      else
      radptl(n)=radiac
      desptl(n)=radiac
      dezptl(n)=0          
      amiptl(n)=prom+amsiac
      endif
      xptl(n)=xaver(1)
      yptl(n)=xaver(2)
      zptl(n)=xaver(3)                    
      tptl(n)=xaver(4)
      z=zptl(n)
      call uttaus(z,tz,sz)
      optl(n)=sz+desptl(n)
      uptl(n)=sz-desptl(n)
      iacptl(n)=1
      iaaptl(n)=1
      enddo
           if(ish.eq.23.and.pptl(5,nptl).ge.amprif)then
      n=nptl
      call jintfp(i,j,n,
     *xptl(i),yptl(i),radptl(i),optl(i),uptl(i),
     *xptl(j),yptl(j),radptl(j),optl(j),uptl(j),
     *xptl(n),yptl(n),radptl(n),optl(n),uptl(n))
           endif
      goto24

25    continue
 
24    continue
      if(i.lt.nptl-1)goto9999

      if(ish.eq.24)then
      if(taus.eq.1..or.taus.eq.2..or.taus.eq.4..or.taus.eq.8.)then
      do 45 i=1,nptl
      if(iacptl(i).eq.0.and.i.ne.nptl)goto45
      call jintcl(i,xptl(i),yptl(i),radptl(i),optl(i),uptl(i),iacptl(i))
45    continue
      endif
      endif

      ishnpt=ish
      if(ishsub/100.eq.19)ish=mod(ishsub,100)
      if(ish.eq.22)then
      write(ifmt,131)nt,nbef,naft,nptl
131   format(1x,'nt=',i5,4x,'nbef=',i8,4x,'naft=',i8,4x,'nptl=',i8)
      endif
      ish=ishnpt
      if(ish.eq.13)call uttima('interactions   ')
      return
      end

c---------------------------------------------------------------------
      subroutine jintcc(i,j,iret)
c---------------------------------------------------------------------
c  iret=1 if i = child of j
c  iret=1 if j = child of i
c  iret=0 else
c---------------------------------------------------------------------
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer iacptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      real xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
      parameter (mxifr=mxptl)
      integer ifr(mxifr),ifrx(mxifr)

      iret=0    
           do kk=1,2
      k1=0
      k2=0
      if(kk.eq.1)n=i
      if(kk.eq.2)n=j
      if(kk.eq.1)m=j
      if(kk.eq.2)m=i
      loop=0
55    loop=loop+1
           if(ifrptl(1,n).gt.0)then
      do 56 k=ifrptl(1,n),ifrptl(2,n)
      if(k.eq.m)goto1001
      k2=k2+1
      if(k2.gt.mxifr)
     *call utstop('jintcc: k2 > mxifr                      ')
56    ifr(k2)=k
           endif
      k1=k1+1
      if(k1.le.k2)then
      n=ifr(k1)
      goto55
      endif
           enddo
      goto1000

1001  iret=1

1000  return
      end

c----------------------------------------------------------------------
      subroutine jintce(i,j,amim,iactn,iret)
c----------------------------------------------------------------------
c  colour exchange interaction of ptls i,j
c  input:
c  i,j: ptl numbers; amf: mass, amim: min mass, of fused obj 
c  iactn=1: ce done
c  iret=0: ce done  or  criteria for ce not fulfilled
c  iret=25: skip
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer iacptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      real xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt

      common/partnr/summas,iptnr,ic4,px4,py4,js4,pex,pey,pezet,pe0,nps
     */didib/ndidib /cnfr/nrfra
      parameter (mxindx=1000,mxre=100,mxma=11,mxmx=6)
      common/crema/indx(mxindx),rema(mxre,mxma),rewi(mxre,mxma)
     *,idmx(mxma,mxmx),icre1(mxre,mxma),icre2(mxre,mxma)
      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc
      parameter(nprbms=20)
      common /cprbms/prbms(nprbms) /cnce/nces,ncolex /ckol/kol 
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      parameter (mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      common/cjint/nptl0,rnup,rnut,velp,velt,bx,by,xaver(4)

      double precision ppfus(5)
      real psum(5),u(5),pei(5),pej(5)
      parameter (nflav=6)
      integer jci(nflav,2),jcj(nflav,2),jc(nflav,2)
     *,ici(2),icj(2),ic(2),ic4(2),jcx(nflav,2),jcdu(nflav,2)
      parameter (mxifr=mxptl)
      integer ifr(mxifr),ifrx(mxifr)
      real proj(nsi,nha),targ(nsi,nha)
      double precision arp(5),arm(5),arq(5),rot(3),boo(5)
      character chadec*5

c     initialization
c     --------------      
      iret=0
      delrap=1.5
      delamf=1.0
      ish00=ish
      pnllx0=pnllx
      taus=ttaus
      amf=pptl(5,nptl+1)

c     check whether ce criteria fulfilled
c     -----------------------------------
      vei=pptl(3,i)/pptl(4,i)
      vej=pptl(3,j)/pptl(4,j)
      rai=0.
      raj=0.
      if(abs(vei).lt.1..and.abs(vej).lt.1.)then
      rai=0.5*alog((1+vei)/(1-vei))
      raj=0.5*alog((1+vej)/(1-vej))
      endif
      if(.not.(amf.gt.amim+delamf.and.abs(rai-raj).gt.delrap
     *.and.iabs(idptl(i)).lt.10000
     *.and.iabs(idptl(j)).lt.10000))goto1000
      iactn=1

c     print
c     -----
      if(ish.ge.92)write(ifch,*)' '
      if(ish.ge.91)write(ifch,101)ntc,taus
      if(ish.ge.92)write(ifch,*)' '
      if(ish.ge.91)write(ifch,115)i,idptl(i)
     *,(pptl(k,i),k=3,5),(xorptl(k,i),k=3,4)
      if(ish.ge.91)write(ifch,115)j,idptl(j)
     *,(pptl(k,j),k=3,5),(xorptl(k,j),k=3,4)
      if(ish.ge.92)write(ifch,*)' '
      if(ish.ge.92)write(ifch,*)'amf,amim: ',amf,amim
      if(ish.ge.92)write(ifch,*)'rai,raj:  ',rai,raj
      if(ish.ge.92)write(ifch,*)' '
101   format(1x,'co.ex. interaction --- nt=',i3,' --- taus=',f6.2)
115   format(1x,'/cptl/',i6,i10
     *,3x,2(e8.2,1x),e8.3,3x,2(e8.2,1x),2x,i4)

c     initialization      
c     --------------
      nptl00=nptl
      nstr00=nstr

      iproj=1          !this leads to iorstr()=-1
      itarg=1          !
      do l=1,nha       !
      kproja(l,1)=1    !
      ktarga(l,1)=1    !
      enddo            !
      kol=1            !
      coord(1,1)=xaver(1)    ! to determine xorstr()
      coord(2,1)=xaver(2)    !
      coord(3,1)=xaver(3)    !
      coord(4,1)=xaver(4)    !

      idp=idptl(i)
      idm=idptl(j)
      do l=1,5
      arp(l)=pptl(l,i)
      arm(l)=pptl(l,j)
      boo(l)=pptl(l,nptl+1)
      enddo
      s=pptl(5,nptl+1)**2

c     boosts into i-j cm
c     ------------------
      call utlob2(1,boo(1),boo(2),boo(3),boo(4),boo(5)    
     *,arp(1),arp(2),arp(3),arp(4))
      call utlob2(1,boo(1),boo(2),boo(3),boo(4),boo(5)    
     *,arm(1),arm(2),arm(3),arm(4))
      if(arp(3).lt.0.)then
      idq=idm
      idm=idp
      idp=idq
      do l=1,5
      arq(l)=arm(l)
      arm(l)=arp(l)
      arp(l)=arq(l)
      enddo
      endif
      pnllx=arp(3)
      rot(1)=(arp(1)-arm(1))/2
      rot(2)=(arp(2)-arm(2))/2
      rot(3)=(arp(3)-arm(3))/2
      call utrot2(1,rot(1),rot(2),rot(3)
     *,arp(1),arp(2),arp(3))
      call utrot2(1,rot(1),rot(2),rot(3)
     *,arm(1),arm(2),arm(3))

c     checks
c     ------
           if(arp(3).lt.0)then
      if(ish.ge. 0)then
      call utmsg('jintce')
      write(ifch,*)'*****  z-comp of +jet < 0.'
      write(ifch,*)(sngl(arp(l)),l=1,4)
      write(ifch,*)(sngl(arm(l)),l=1,4)
      call utmsgf
      endif
      goto10025
           endif
           if(abs(sngl(arp(4)+arm(4))-sqrt(s)).gt.
     *2.e-2*sqrt(s))then
      if(ish.ge. 0)then
      call utmsg('jintce')
      write(ifch,*)'*****  arp(4)+arm(4)-sqrt(s) nonzero'
      write(ifch,*)'value:   ',sngl(arp(4)+arm(4))-sqrt(s)
      write(ifch,*)'sqrt(s): ',sqrt(s)
      write(ifch,*)(sngl(arp(i)),i=1,4)
      write(ifch,*)(sngl(arm(i)),i=1,4)
      call utmsgf
      endif
           endif

c     fill proj, targ
c     ---------------
      call idtr4(idp,ic)
      do m=1,nha
      do n=1,nsi
      proj(n,m)=0
      enddo
      enddo
      proj(5,1)=ic(1)
      proj(6,1)=ic(2)
      proj(5,2)=ic(1)
      proj(6,2)=ic(2)
      do l=1,4
      proj(l,2)=arp(l)
      enddo
      call idtr4(idm,ic)
      do m=1,nha
      do n=1,nsi
      targ(n,m)=0
      enddo
      enddo
      targ(5,1)=ic(1)
      targ(6,1)=ic(2)
      targ(5,2)=ic(1)
      targ(6,2)=ic(2)
      do l=1,4
      targ(l,2)=arm(l)
      enddo

c     redo
c     ----
      call utremb(proj,targ,1)
      loo=1
      goto58
57    continue
      loo=loo+1
      if(loo.gt.20)goto10025
      if(ish.ge.91)write(ifch,*)'redo hh collision'
      call utrest(proj,targ,1)
58    continue
 
c     number of ce's
c     --------------
      ncolex=1
c-c   lo=0
c-c16 lo=lo+1
c-c        if(lo.eq.3)then
c-c   if(ish.ge.90)then
c-c   call utmsg('jintce')
c-c   write(ifch,*)'*****  lo=3'
c-c   call utmsgf
c-c   endif
c-c        endif
c-c   r=rangen()
c-c   ncolex=0
c-c15 ncolex=ncolex+1
c-c   if(ncolex.gt.nprbms)goto16
c-c   if(r.gt.prbms(ncolex))goto15
 
c     do ncolex ce's
c     --------------
      iskip=0
      do nce=1,ncolex
      if(ish.ge.92.and.nce.gt.1)write(ifch,*)' '
      if(ish.ge.92.and.nce.gt.1)write(ifch,*)nce,'. colour exchange'
      if(ish.ge.92.and.nce.gt.1)write(ifch,*)' '
      nces=nce
      ish=ish-2
      call hahabs(proj,targ,nce/ncolex,nce/ncolex,iskip,irethh)
      ish=ish+2
      if(iskip.ne.0)goto10025
      if(irethh.eq.1)goto57
      enddo

c     fragmentation
c     -------------
      if(nstr.le.nstr00)goto10050

      if(ish.ge.92)write(ifch,123)
123   format(/1x,'strings before rescaling'/)
      do l=nstr00+1,nstr
      if(ish.ge.92)write(ifch,109)l,(icstr(k,l)/100,k=1,4)
     *,sqrt(pstr(1,l)**2+pstr(2,l)**2),pstr(3,l),pstr(5,l)
     *,irlstr(l)
      enddo
109   format(' /cstr/',i4,3x,4i5,2x,3(1x,e9.3),i4)
 
      if(irescl.eq.1)then
      psum(1)=0
      psum(2)=0
      psum(3)=0
      psum(4)=sqrt(s)
      psum(5)=sqrt(s)
      ish=ish-2
      call hrescl(nstr00+1,nstr,psum,ifail)
      ish=ish+2
      if(ifail.ne.0)goto57
      endif
 
      if(ish.ge.92)write(ifch,124)
124   format(/1x,'strings after rescaling'/)
      do l=nstr00+1,nstr
      if(ish.ge.92)write(ifch,109)l,(icstr(k,l)/100,k=1,4)
     *,sqrt(pstr(1,l)**2+pstr(2,l)**2),pstr(3,l),pstr(5,l)
     *,irlstr(l)
      enddo
  
      do l=nstr00+1,nstr
      if(ish.ge.92)write(ifch,122)l
     *,(icstr(k,l),k=1,4)
     *,sqrt(pstr(1,l)**2+pstr(2,l)**2),pstr(3,l),pstr(5,l)
122   format(/'str:',i7,4i7,1x,3(e8.2,1x)/)
      ish=ish-2
      call jamfra(l,newevt)
      ish=ish+2
      if(newevt.eq.1)goto10025
      enddo

c     boost particles from string fragmentation
c     -----------------------------------------
      if(nptl.le.nptl00)goto10051

      iorptl(nptl00+1)=i
      jorptl(nptl00+1)=j
      do l=1,5
      psum(l)=0
      enddo

           do n=nptl00+1,nptl

      do l=1,5
      arp(l)=pptl(l,n)
      enddo
      call utrot2(-1,rot(1),rot(2),rot(3)
     *,arp(1),arp(2),arp(3))
      call utlob2(-1,boo(1),boo(2),boo(3),boo(4),boo(5)    
     *,arp(1),arp(2),arp(3),arp(4))
      do l=1,5
      pptl(l,n)=arp(l)
      enddo
      nqjptl(n)=0

      do l=1,4
      if(istptl(n).eq.0)psum(l)=psum(l)+pptl(l,n)
      enddo

      if(ish.ge.91)then
      write(ifch,115)n,idptl(n)
     *,(pptl(k,n),k=3,5),(xorptl(k,n),k=3,4),istptl(n)
      endif

           enddo

      if(ish.ge.92)then
      write(ifch,*)'p_fin:',(psum(k),k=1,4)
      write(ifch,*)'p_ini:',(boo(k),k=1,4)
      endif
               
      if(ish.ge.90.and.
     *   (abs(boo(1)-psum(1)).gt.5.e-3*abs(boo(1))
     *.or.abs(boo(2)-psum(2)).gt.5.e-3*abs(boo(2))
     *.or.abs(boo(3)-psum(3)).gt.5.e-3*abs(boo(3))
     *.or.abs(boo(4)-psum(4)).gt.1.e-1*abs(boo(4)))
     *.and.
     *   (abs(boo(1)-psum(1)).gt.5.e-3
     *.or.abs(boo(2)-psum(2)).gt.5.e-3
     *.or.abs(boo(3)-psum(3)).gt.5.e-3
     *.or.abs(boo(4)-psum(4)).gt.1.e-1))then
      call utmsg('jintce')
      write(ifch,*)'*****  p_ini /= p_fin'
      write(ifch,*)'final particles:'
      do n=nptl00+1,nptl
      write(ifch,125)n,idptl(n),(pptl(k,n),k=1,5),istptl(n)                                   
125   format(1x,'/cptl/',i6,i10,3x,5(e9.3,1x),i2)
      enddo 
      write(ifch,*)'p_fin:',(psum(k),k=1,4)
      write(ifch,*)'p_ini:',(sngl(boo(k)),k=1,4)
      call utmsgf
      endif
      goto1000

c     finish
c     ------

10025 continue
      iret=25
      if(ish.ge.91)write(ifch,*)'skip'
      nptl=nptl00
      nstr=nstr00
      goto1000

10050 continue      
      call utmsg('jintce')
      write(ifch,*)'incident particles:'
      write(ifch,115)i,idptl(i),(pptl(k,i),k=3,5)      
      write(ifch,115)j,idptl(j),(pptl(k,j),k=3,5)
      write(ifch,*)'nstr=',nstr,'    nstr00=',nstr00
      call utmsgf
      call utstop('jintce: nstr .le. nstr00                ')

10051 continue      
      call utmsg('jintce')
      write(ifch,*)'incident particles:'
      write(ifch,115)i,idptl(i),(pptl(k,i),k=3,5)      
      write(ifch,115)j,idptl(j),(pptl(k,j),k=3,5)
      write(ifch,*)'nstr=',nstr,'    nstr00=',nstr00
      write(ifch,*)'nptl=',nptl,'    nptl00=',nptl00
      call utmsgf
      call utstop('jintce: nptl .le. nptl00                ')
  
1000  continue
      pnllx=pnllx0
      ish=ish00
      return
      end

c----------------------------------------------------------------------
      subroutine jintch(i,j,kmax)
c----------------------------------------------------------------------
c  writes children of i,j to ifrij(1-kmax)
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter (mxifr=mxptl)
      common/cifrij/ifrij(mxifr)

      kcur=0
      kmax=0
      loop=0
31    loop=loop+1
      if(loop.eq.1)n=i
      if(loop.eq.2)n=j
           if(ifrptl(1,n).gt.0)then
      if(ish.ge.92)write(ifch,*)n,' ---> ',ifrptl(1,n),ifrptl(2,n)
      do 30 k=ifrptl(1,n),ifrptl(2,n)
      kmax=kmax+1
      if(kmax.gt.mxifr)
     *call utstop('jintch: kmax > mxifr                    ')
30    ifrij(kmax)=k
           endif
      if(loop.eq.1)goto31
      kcur=kcur+1
      if(kcur.le.kmax)then
      n=ifrij(kcur)
      goto31
      endif

      return
      end

c---------------------------------------------------------------------
      subroutine jintcl(i,x,y,rad,o,u,iac)
c---------------------------------------------------------------------
c  fills histogram concerning cluster characteristics
c---------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/paro5/delvol,deleps
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/nevnt/nevnt /cvsn/iversn

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)

      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat 
c MAB      parameter (mxtau=4,mxvol=10,mxeps=10)
      parameter (mxtau=4,mxvol=10,mxeps=20)
      common/cjintc/clust(mxtau,mxvol,mxeps)

      character ax*2

      tau=sngl(ttaus)
      itau=1+alog(tau)/alog(2.)

      if(iac.eq.0)goto1

      vol=(o-u)*rad**2*pi
      ivol=1+vol/delvol
      eps=pptl(5,i)/vol
      ieps=1+eps/deleps
      if(itau.ge.1.and.itau.le.mxtau
     *.and.ivol.ge.1.and.ivol.le.mxvol
     *.and.ieps.ge.1.and.ieps.le.mxeps)
     *clust(itau,ivol,ieps)=clust(itau,ivol,ieps)+1

1     continue
      if(i.lt.nptl.or.nrevt+1.lt.nevnt)goto1000

      if(itau.eq.1)write(ifch,105)maproj,matarg,engy,iversn/100.
105   format(/1x,'proj=',i3,3x,'targ=',i3,3x,'engy=',f7.2
     *,3x,'VENUS ',f4.2)
      write(ifch,100)tau,itau,nevnt,delvol,deleps
100   format(/1x,'tau=',f5.2,3x,'itau=',i1
     *,3x,'nevnt=',i4,3x,'delvol=',f5.1,3x,'deleps=',f5.1/)
      write(ifch,101)
101   format(9x,'ivol=1 ivol=2 ivol=3 ivol=4 ivol=5 '
     *,'ivol=6 ivol=7 ivol=8 ivol=9 ivol=10 ')
      do ie=1,mxeps
      if(ie.le.9)write(ax,102)ie
102   format('0',i1)
      if(ie.gt.9)write(ax,103)ie
103   format(i2)
      write(ifch,104)ax,(nint(clust(itau,iv,ie)),iv=1,mxvol)
104   format(' ieps=',a2,10i7)
      enddo

1000  continue
      return
      end

c----------------------------------------------------------------------
      subroutine jintel(i,j,amim,iactn)
c----------------------------------------------------------------------
c  elastic scattering of ptls i,j
c  el scatt requir not fulf: return without action, iactn unchanged
c    else: elastic scattering done, iactn=2
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common/cjint/nptl0,rnup,rnut,velp,velt,bx,by,xaver(4)

      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc

      parameter (nflav=6)
      integer jci(nflav,2),jcj(nflav,2),jc(nflav,2),jcx(nflav,2)
      real p(5),u(5),pei(5),pej(5)
      double precision ppfus(5)

c     check
c     -----
      amf=pptl(5,nptl+1)
      if(amf.ge.amim)goto1000

c     print
c     -----
      if(ish.ge.91)then
      write(ifch,101)ntc,sngl(ttaus)
101   format(1x,'elastic scattering --- nt=',i3,' --- taus=',f6.2)
      write(ifch,115)i,idptl(i)
     *,(pptl(k,i),k=3,5),(xorptl(k,i),k=3,4)
      write(ifch,115)j,idptl(j)
     *,(pptl(k,j),k=3,5),(xorptl(k,i),k=3,4)
115   format(1x,'/cptl/',i6,i10
     *,3x,2(e8.2,1x),e8.3,3x,2(e8.2,1x),2x,i4)
      endif

c     initialization
c     --------------
      iactn=2
      do k=1,5
      p(k)=pptl(k,nptl+1)
      enddo

c     determine momenta of outgoing particles (pei,pej)
c     -------------------------------------------------
           if(p(5).le.(pptl(5,i)+pptl(5,j))*.99)then
      if(ish.ge.90)then
      call utmsg('jintel')
      write(ifch,132)p(5),pptl(5,i)+pptl(5,j)
132   format(1x,'*****  m_fus < m_i+m_j ---> qcm set zero    ( '
     *,2f10.3,' )')
      write(ifch,133)'p_i:  ',(pptl(k,i),k=1,5)
      write(ifch,133)'p_j:  ',(pptl(k,j),k=1,5)
      write(ifch,133)'p_fus:',(p(k),k=1,5)
133   format(1x,a6,3x,5f10.3) 
      call utmsgf
      endif
      qcm=0.
           elseif(p(5).le.pptl(5,i)+pptl(5,j))then
      qcm=0.
           else
      qcm=utpcm(p(5),pptl(5,i),pptl(5,j))
           endif

c isotropic
      u(3)=2.*rangen()-1.
      phi=2.*pi*rangen()
      u(1)=sqrt(1.-u(3)**2)*cos(phi)
      u(2)=sqrt(1.-u(3)**2)*sin(phi)
      do 47 k=1,3
      pei(k)= qcm*u(k)
47    pej(k)=-qcm*u(k)

c nonisotropic
c-c   pt=utinvt(nptq,xptq,qptq,rangen()*qptq(nptq))
c-c   if(pt.ge.qcm)pt=rangen()*qcm
c-c   qpl=sqrt(qcm**2-pt**2)
c-c   u(3)=qpl
c-c   phi=2.*pi*rangen()
c-c   u(1)=pt*cos(phi)
c-c   u(2)=pt*sin(phi)
c-c   call utaxis(i,j,a1,a2,a3)
c-c   ivt=1
c-c   if(a3.lt.0.)then
c-c   a1=-a1
c-c   a2=-a2
c-c   a3=-a3
c-c   ivt=-1
c-c   endif 
c-c   call utrota(-1,a1,a2,a3,u(1),u(2),u(3))
c-c   do 47 k=1,3
c-c   pei(k)= u(k)*ivt
c-c47 pej(k)=-u(k)*ivt

      pei(4)=sqrt(qcm**2+pptl(5,i)**2)
      pej(4)=sqrt(qcm**2+pptl(5,j)**2)
      pei(5)=pptl(5,i)
      pej(5)=pptl(5,j)
      call utlobo(-1,p(1),p(2),p(3),p(4),p(5)
     *,pei(1),pei(2),pei(3),pei(4))
      call utlobo(-1,p(1),p(2),p(3),p(4),p(5)
     *,pej(1),pej(2),pej(3),pej(4))

c     fill /cptl/
c     -----------
      do 49 lo=1,2
      nptl=nptl+1
      if(lo.eq.1)ij=i
      if(lo.eq.2)ij=j
      do 48 k=1,5
      if(lo.eq.1)pptl(k,nptl)=pei(k)
      if(lo.eq.2)pptl(k,nptl)=pej(k)
48    continue
      istptl(nptl)=0
      idptl(nptl)=idptl(ij)
      ibptl(1,nptl)=ibptl(1,ij)
      ibptl(2,nptl)=ibptl(2,ij)
      ibptl(3,nptl)=ibptl(3,ij)
      ibptl(4,nptl)=ibptl(4,ij)
      xorptl(1,nptl)=xaver(1)           
      xorptl(2,nptl)=xaver(2)
      xorptl(3,nptl)=xaver(3)
      xorptl(4,nptl)=xaver(4)
      iorptl(nptl)=i
      jorptl(nptl)=j
      tivptl(1,nptl)=xaver(4)                
      call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
      tivptl(2,nptl)=tivptl(1,nptl)+taugm
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      iclptl(nptl)=iclptl(ij)
      nqjptl(nptl)=nqjptl(ij)
      if(ish.ge.91)write(ifch,115)nptl,idptl(nptl)
     *,(pptl(k,nptl),k=3,5),(xorptl(k,nptl),k=3,4)
49    continue

1000  return
      end

c----------------------------------------------------------------------
      subroutine jintfp(i,j,n,
     *xpl1,ypl1,radpl1,opl1,upl1,
     *xpl2,ypl2,radpl2,opl2,upl2,
     *xpl3,ypl3,radpl3,opl3,upl3)
c----------------------------------------------------------------------
c  printout
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/cnsta/pi,prom,piom,ainfin

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)

      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat 

      real xpl(3),ypl(3),radpl(3),opl(3),upl(3)
      character ax*10,lin*59,mark*1,txt*8

      xpl(1)  =xpl1
      ypl(1)  =ypl1
      radpl(1)=radpl1
      opl(1)  =opl1
      upl(1)  =upl1
      xpl(2)  =xpl2
      ypl(2)  =ypl2
      radpl(2)=radpl2
      opl(2)  =opl2
      upl(2)  =upl2
      xpl(3)  =xpl3
      ypl(3)  =ypl3
      radpl(3)=radpl3
      opl(3)  =opl3
      upl(3)  =upl3

      v1=(opl1-upl1)*radpl1**2*pi
      v2=(opl2-upl2)*radpl2**2*pi
      v3=(opl3-upl3)*radpl3**2*pi

      write(ifch,103)sngl(ttaus)
103   format(' monitoring fusion interaction at ttaus=',f5.2)
      tx=ttaus
      tx=amax1(tx,1.) 
      write(txt(3:7),104)tx
104   format(f5.2)
      write(ifch,105)'incoming:',i,idptl(i),pptl(5,i),v1,pptl(5,i)/v1
      write(ifch,105)'incoming:',j,idptl(j),pptl(5,j),v2,pptl(5,j)/v2
      write(ifch,105)'fused:   ',n,idptl(n),pptl(5,n),v3,pptl(5,n)/v3
105   format(1x,a9,2x,'nr=',i6,2x,'id=',i9,2x,'mass=',f5.1,2x
     *,'vol=',f5.1,2x,'eps=',f4.1)

      do mm=1,3
      if(mm.eq.1)txt(1:2)='s/'
      if(mm.eq.1)txt(8:8)=':'
      if(mm.eq.2)txt='x:      '
      if(mm.eq.3)txt='y:      '
      ax='---------!'
      write(ifch,100)txt,(ax,l=1,6)
100   format(2x,a8,3x,'-3',8x,'-2',8x,'-1',9x,'0',9x,'1',9x,'2',9x,'3'
     *      /13x,' !',6a10)
      do k=1,3
      do l=1,59
      wo=-3+l*0.1 + 0.05
      wu=-3+l*0.1 - 0.05
           if(mm.eq.1)then
      o=opl(k)/tx
      u=upl(k)/tx
           elseif(mm.eq.2)then
      o=xpl(k)+radpl(k)
      u=xpl(k)-radpl(k)
           elseif(mm.eq.3)then
      o=ypl(k)+radpl(k)
      u=ypl(k)-radpl(k)
           endif
      mark='O'
      if(k.eq.3)mark='X'
      lin(l:l)=' '
      if(wu.le.o.and.wo.ge.u)lin(l:l)=mark
      enddo
      write(ifch,101)u,o,lin
      enddo
101   format(1x,2f6.2,' !',a59,'!')
      write(ifch,102)('-',l=1,59)
102   format(13x,' !',59a1,'!')
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine jintfs(i,j,nqi,jc,amim,iret) 
c----------------------------------------------------------------------
c  input: ptl numbers i,j   
c  output: ppfus(5): momentum , nqi: net quark number,                  
c              jc: jc-code, amim: minimum mass,    of fused ptl
c          iret=0 if ok, 1 else
c          ppfus() written to pptl(,nptl+1)
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common/cjint/nptl0,rnup,rnut,velp,velt,bx,by,xaver(4)
      parameter (nflav=6)
      integer jci(nflav,2),jcj(nflav,2),jc(nflav,2),jcx(nflav,2)
      real p(5)
      double precision ppfus(5),pp52

      iret=0

      do 35 k=1,4
      p(k)=pptl(k,i)+pptl(k,j)
35    ppfus(k)=p(k)
      pp52=ppfus(4)**2-ppfus(3)**2-ppfus(2)**2-ppfus(1)**2
           if(pp52.le.0.)then
      if(ish.ge.90)then
      call utmsg('jintfs')
      write(ifch,*)'*****  mfus**2 < 0    (',pp52,' )'
      write(ifch,*)(ppfus(m),m=1,4)
      call utmsgf
      endif
      goto1001
           endif
      ppfus(5)=dsqrt(pp52)
      if(nptl+1.gt.mxptl)
     *call utstop('jintfs: nptl>mxptl                      ')
      do 36 k=1,5
36    pptl(k,nptl+1)=ppfus(k)
      amf=ppfus(5)

      istptl(nptl+1)=0

      call idquac(i,ndummy,ndummy,ndummy,jci)
      call idquac(j,ndummy,ndummy,ndummy,jcj)
      do 29 n=1,nflav
      do 29 k=1,2
29    jc(n,k)=jci(n,k)+jcj(n,k)
      nqi=0
      do 54 n=1,nflav
54    nqi=nqi+jc(n,1)-jc(n,2)
                       
      call idcomj(jc)
      amim=utamnz(jc,5)+.200
      if(amf.lt.amim.and.i.gt.nptl0.and.j.gt.nptl0)goto1001
      goto1000

1001  iret=1

1000  return
      end

c----------------------------------------------------------------------
      subroutine jintfu(i,j,jc,iactn)
c----------------------------------------------------------------------
c  fusion of ptls i,j: determine cluster
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)

      parameter (mxlook=10000,mxdky=2000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      common/cjint/nptl0,rnup,rnut,velp,velt,bx,by,xaver(4)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc

      parameter (nflav=6)
      integer jc(nflav,2),ic(2),ib(4)

c     print
c     -----
      if(ish.ge.91)then
      write(ifch,101)ntc,sngl(ttaus)
101   format(1x,'fusion interaction --- nt=',i3,' --- taus=',f6.2)
      write(ifch,115)i,idptl(i)
     *,(pptl(k,i),k=3,5),(xorptl(k,i),k=3,4)
      write(ifch,115)j,idptl(j)
     *,(pptl(k,j),k=3,5),(xorptl(k,i),k=3,4)
115   format(1x,'/cptl/',i6,i10
     *,3x,2(e8.2,1x),e8.3,3x,2(e8.2,1x),2x,i4)
      endif

c     initialization
c     --------------
      iactn=3
      amf=pptl(5,nptl+1)
      nptl=nptl+1

c     determine idr, ib(1-4)
c     ----------------------
      idr=0
      do 40 nf=1,nflav
      do 40 ij=1,2
      if(jc(nf,ij).ge.10)idr=7*10**8
40    continue
           if(idr/10**8.ne.7)then
      call idenco(jc,ic,ireten)
      if(ireten.eq.1)
     *call utstop('jintfu: idenco ret code = 1             ')
      id=idtra(ic,0,0,3)
43    amc=amf
      call idres(id,amc,idr,iadj)
           if(idr.ne.0)then
      lid=look(iabs(idr))
      if((lid.le.0.or.lid.gt.0.and.mode(2,lid).eq.0)
     *.and.pptl(5,nptl).gt.amc+1e-3)then
      amf=amf+0.010
      goto43
      endif
      if((lid.le.0.or.lid.gt.0.and.mode(2,lid).eq.0)
     *.and.abs(amc-pptl(5,nptl)).gt.1e-3)then
      if(ish.ge.90)then
      call utmsg('jintfu')
      write(ifch,*)'*****  not on mass shell after fusion: '
     *,pptl(5,nptl),amc
      call utmsgf
      endif
      endif
           endif
           if(idr.eq.0)then
      if(mod(ic(1),100).ne.0.or.mod(ic(2),100).ne.0)then
      idr=9*10**8
      else
      idr=8*10**8+ic(1)*100+ic(2)/100
      endif
           endif
           else
      call idtrbi(jc,ib(1),ib(2),ib(3),ib(4))
      idr=idr
     *+mod(jc(1,1)+jc(2,1)+jc(3,1)+jc(4,1),10**4)*10**4
     *+mod(jc(1,2)+jc(2,2)+jc(3,2)+jc(4,2),10**4)
      if(ish.ge.93)write(ifch,*)'ib:',(ib(kk),kk=1,4)
           endif

c     fill /cptl/
c     -----------
      idptl(nptl)=idr
      do k=1,4
      ibptl(k,nptl)=ib(k)
      enddo
      xorptl(1,nptl)=xaver(1)
      xorptl(2,nptl)=xaver(2)
      xorptl(3,nptl)=xaver(3)
      xorptl(4,nptl)=xaver(4)
      iorptl(nptl)=i
      jorptl(nptl)=j
      tivptl(1,nptl)=xaver(4) 
      call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
      tivptl(2,nptl)=tivptl(1,nptl)+taugm
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      iclptl(nptl)=1
      nqjptl(nptl)=0

c     print + return
c     --------------
      if(ish.ge.91)then
      n=nptl
      write(ifch,115)n,idptl(n)
     *,(pptl(k,n),k=3,5),(xorptl(k,n),k=3,4)
      endif
      return
      end

c----------------------------------------------------------------------
      subroutine jintpa(i,j,kmax)
c----------------------------------------------------------------------
c  writes parents of i,j to ifrij(1-kmax)
c  sets iaaptl()=0 for parents
c----------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax

      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      parameter (mxifr=mxptl)
      common/cifrij/ifrij(mxifr)

      kcur=0
      kmax=0
      loop=0
12    loop=loop+1
      if(loop.eq.1)n=i
      if(loop.eq.2)n=j
           if(iorptl(n).gt.0)then
      if(ish.ge.92)write(ifch,*)n,' <--- ',iorptl(n),jorptl(n)
      if(kmax+2.gt.mxifr)
     *call utstop('jintpa: kmax+2 > mxifr                  ')
      if(iaaptl(iorptl(n)).ne.0)then
      kmax=kmax+1
      ifrij(kmax)=iorptl(n)
      iaaptl(iorptl(n))=0
      endif
      if(jorptl(n).gt.0)then
      if(iaaptl(jorptl(n)).ne.0)then
      kmax=kmax+1
      ifrij(kmax)=jorptl(n)
      iaaptl(jorptl(n))=0
      endif
      endif
           endif
      if(loop.eq.1)goto12
8     kcur=kcur+1
      if(kcur.le.kmax)then
      n=ifrij(kcur)
      if(iaaptl(n).ne.0)
     *call utstop('jintpa: should not happen               ')
      goto12
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine jrescl(j1,j2,psum,ifail)
c     rescales ptl momenta of ptls j1-j2 to have total mom psum.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      dimension psum(5),beta(3)
      double precision ppsum(5),padd(5),pp(5)
      data errlim/.001/
      common/cscal/scal

      ifail=1

      do k=1,5
      ppsum(k)=psum(k)
      enddo

      if(j1.ge.j2)
     *call utstop('jrescl: j1 .ge. j2                      ')

      do 100 k=1,5
100   padd(k)=0.
      do 110 j=j1,j2
      do 110 k=1,5
      padd(k)=padd(k)+pptl(k,j)
110   continue
           if(padd(5).ge.ppsum(5))then
      if(ish.ge.90)then
      call utmsg('jrescl')
      write(ifch,*)'*****  sum of ptl masses .ge. ppsum(5)'
      do 1 j=j1,j2
1     write(ifch,109)j,idptl(j)
     *,(pptl(l,j),l=3,5)
109   format(' /cptl/',i6,i10,3x,3(1x,e9.3))
      write(ifch,*)'ppsum(345):',(sngl(ppsum(k)),k=3,5)
      call utmsgf
      endif
      return
           endif
      padd(5)=padd(4)**2-padd(1)**2-padd(2)**2-padd(3)**2
           if(padd(5).le.0)then
      do 2 j=j1,j2
2     write(ifch,108)j,(pptl(l,j),l=1,5)
108   format(' /cptl/',i4,2x,5(1x,e9.3))
      call utstop('jrescl: mass**2 of string-sum negative  ')
           endif
      padd(5)=sqrt(padd(5))

c     boost ptls to rest
c     ------------------
           do j=j1,j2
      do k=1,4
      pp(k)=pptl(k,j)
      enddo 
      call utlob2(1,padd(1),padd(2),padd(3),padd(4),padd(5)    
     *,pp(1),pp(2),pp(3),pp(4))
      do k=1,4
      pptl(k,j)=pp(k)
      enddo 
           enddo  
 
c     rescale momenta in rest frame
c     -----------------------------
      scal=1.
      do 301 ipass=1,200
      sum=0.
      do 310 j=j1,j2
      do 320 k=1,3
320   pptl(k,j)=scal*pptl(k,j)
      pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2
     *+pptl(5,j)**2)
      sum=sum+pptl(4,j)
310   continue
      scal=psum(5)/sum
      if(abs(scal-1.).le.errlim) goto300
301   continue
      if(ish.ge.90)then
      call utmsg('jrescl')
      write(ifch,*)'*****  scal=',scal
      call utmsgf
      endif
300   continue
 
c     boost back with ppsum
c     ---------------------
           do j=j1,j2
      do k=1,4
      pp(k)=pptl(k,j)
      enddo 
      call utlob2(-1,ppsum(1),ppsum(2),ppsum(3),ppsum(4),ppsum(5)    
     *,pp(1),pp(2),pp(3),pp(4))
      do k=1,4
      pptl(k,j)=pp(k)
      enddo 
           enddo  

      ifail=0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jsplit(stro,str,kolsp,ier,kmaxor)
c-----------------------------------------------------------------------
c  splits string stro into q-qbar string str and remainder (->stro)
c  dimensions: stro(nsi,nsix+1),str(nsi,2)
c  ier=0: ok ; ier=1: error ; ier=2: absorption ;
c  ier=3: again with new apart, epart
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (nptf=129)
      common /cptf/xptf(nptf),wptf(nptf)
     *,qptfu(nptf),qptfs(nptf),qptfuu(nptf),qptfus(nptf),qptfss(nptf)
      parameter (nptq=129)
      common /cptq/xptq(nptq),qptq(nptq),qpth(nptq),wptq(nptq)
      parameter (nsi=6,nsix=100)
      common/cnfusn/nfusn(nsix+1)
      parameter (nflav=6)
      integer jc(nflav,2),ic(2),jcx(nflav,2),icx(2)
      parameter (nha=nsix+2,mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      common/cjspli/jpart,nscc,nsccx,apart,epart,alead,elead,sgnsil
      parameter (nsplit=129)
      common /csplit/xsplit(nsplit),qsplit(nsplit),qsplix(nsplit)
     *,wsplit(nsplit)
      common/cpzstr/pzstrl,estrl,istrl,isea

      real stro(nsi,nsix+1),stro0(nsi,nsix+1),str(nsi,2)
      double precision sstr(nsi,2),po(5),po123,a,d,pnl3,pnl4
     *,pox,poy,pot,pop,pom,pux,puy,put,pup,pum,pap,pam,pep,pem,pyp,pym
     *,pip,pim,pp,pm

      pudx=pud*.945
      ish0=ish
      if(ishsub/100.eq.15)ish=mod(ishsub,100)
 
      call utkstr(stro,kmax)
           if(kmax.eq.kmaxor)then
      nscc=0
      if(jpart.eq.0)then
      do i=1,kmaxor
      nfusn(i)=0
      enddo
      endif
           endif
      kmax0=kmax
      do 10 i=1,nsi
      do 10 k=1,kmax
10    stro0(i,k)=stro(i,k)
c-c   pdiq=pdiqua/2
c-c   if(kmaxor.gt.kutdiq)pdiq=0.
c-c   pdiq=pdiqua*0.5*nsccx/(kmaxor-1.)
      pdiq=0
      loop=0
      ncorn=0
      xunter=0.
      xober=1.
      messco=0
      kolsp=kproja(kmax+1,iproj)
      if(sgnsil.lt.0.)kolsp=ktarga(kmax+1,itarg)
      pap=2*epart
      pam=0
      pux=stro(1,kmax)
      puy=stro(2,kmax)
      put=dsqrt(pux**2+puy**2)
      pup=stro(4,kmax)-abs(stro(3,kmax))
      pum=stro(4,kmax)+abs(stro(3,kmax))
 
c     print
c     -----
      if(ish.ge.92)then
      if(ish.ge.93)write(ifch,*)('-',l=1,79)
      write(ifch,*)'multi-string detected. split off q-qbar string:'
      if(ish.ge.93)write(ifch,*)('-',l=1,79)
      write(ifch,*)' '
      if(ish.ge.93)then
      write(ifch,*)'input string stro:'
      write(ifch,*)' '
      write(ifch,105)(stro(i,1),i=1,4),(nint(stro(i,1)),i=5,6)
      do 8 k=2,kmax
8     write(ifch,104)(stro(i,k),i=1,4),(nint(stro(i,k)),i=5,6)
      write(ifch,*)' '
      endif
      endif
 
c     absorption
c     ----------
      if(istrl.eq.1)goto1002
 
c     reset
c     -----
5001  loop=0
5000  loop=loop+1
      if(loop.ge.5)goto1002
      kmax=kmax0
      do 11 i=1,nsi
      do 11 k=1,kmax
11    stro(i,k)=stro0(i,k)
 
c     split off hadron + determine remainder string (momenta)
c     -------------------------------------------------------
      nqu=1
      if(rangen().lt.pdiq)nqu=2
      ifltt=0
      iflto=0
      do 24 i=1,nqu
      ifl=int(rangen()/pudx)+1
      iflto=iflto*10+ifl
24    ifltt=ifltt*10+(ifl+1)/2
      if(ifltt.eq.1)then
      pt=utinvt(nptf,xptf,qptfu ,rangen()*qptfu(nptf))
      elseif(ifltt.eq.2)then
      pt=utinvt(nptf,xptf,qptfs ,rangen()*qptfs(nptf))
      elseif(ifltt.eq.11)then
      pt=utinvt(nptf,xptf,qptfuu,rangen()*qptfuu(nptf))
      elseif(ifltt.eq.12.or.ifltt.eq.21)then
      pt=utinvt(nptf,xptf,qptfus,rangen()*qptfus(nptf))
      elseif(ifltt.eq.22)then
      pt=utinvt(nptf,xptf,qptfss,rangen()*qptfss(nptf))
      endif
      at=sqrt(apart**2+pt**2)
      r=rangen()
c-c   x=sqrt(1-sqrt(r))                !f(x)=x-x**3
      x=utinvt(nsplit,xsplit,qsplit,r*qsplit(nsplit))
      if(nsccx.ge.2)
     *x=utinvt(nsplit,xsplit,qsplix,r*qsplix(nsplit))
      ystrl=.5*alog((estrl+pzstrl)/(estrl-pzstrl))
      xunter=abs(0.94*sinh(ystrl))/epart
      if(xunter.gt.1)xunter=1.
      xunter=xunter-(xober-xunter)
      x=xunter+x*(xober-xunter)
c-c   en=apart+x*(epart-apart)
      pz=x*epart
      sgnpo=1
      if(pz.lt.0)sgnpo=-1
      en=sqrt(pz**2+at**2)
c-c   if(at.gt.en)at=apart+rangen()*(en-apart)
c-c   pt=sqrt((at-apart)*(at+apart))
      phi=2.*pi*rangen()
      po(1)=pt*cos(phi)
      po(2)=pt*sin(phi)
c-c   po(3)=sgnsil*sqrt(en**2-at**2)
      po(3)=sgnsil*pz
      po(4)=en
      po123=sqrt(po(3)**2+po(2)**2+po(1)**2)
           if(po(4)-po123.gt.0.)then
      po(5)=sqrt((po(4)-po123)*(po(4)+po123))
           else
      po(5)=0
      if(po(4)-po123.lt.-1e-4*po(4))then
      call utmsg('jsplit')
      write(ifch,*)'*****  !po123! > po4'
      write(ifch,*)'po123,po4:',po123,po(4)
      write(ifch,*)'po1:',po(1)
      write(ifch,*)'po2:',po(2)
      write(ifch,*)'po3:',po(3)
      write(ifch,*)'po5:',po(5)
      call utmsgf
      endif
           endif
      pox=po(1)
      poy=po(2)
      pot=dsqrt(pox**2+poy**2)
      pop=po(4)+dabs(po(3))*sgnpo
      pom=po(4)-dabs(po(3))*sgnpo
      stro(3,1)=stro(3,1)-epart*sgnsil
      stro(4,1)=stro(4,1)-epart
      do 2 i=1,4
      stro(i,1)=stro(i,1)+po(i)
2     stro(i,kmax)=0.
      yld=100
      if(po(4)-dabs(po(3)).ne.0.)
     *yld=0.5*dlog((po(4)+dabs(po(3)))/(po(4)-dabs(po(3))))
 
c     determine split string (momenta)
c     --------------------------------
      pp=pap+pup
      pm=pam+pum
      pep=pp-pop
      pem=pm-pom
      if(pep.le.0.)goto5000
      if(pem.le.0.)goto5000
      a=(pem*pep-put**2-pot**2)/2.
      d=put*pot
      if(a**2-d**2.lt.0.)goto5000
      pyp=(a+put**2-dsqrt(a**2-d**2))
      if(pyp.lt.0..and.pyp.gt.-1e-6)pyp=0.
      pyp=pyp/pem
      pym=(a+pot**2-dsqrt(a**2-d**2))
      if(pym.lt.0..and.pym.gt.-1e-6)pym=0.
      pym=pym/pep
      if(pyp.lt.0.)goto5000
      if(pym.lt.0.)goto5000
      pip=pep-pyp
      if(pip.lt.0.)goto5000
      pim=pym
      pap=pyp
      pam=pem-pym
      if(pam.lt.0.)goto5000
      sstr(1,1)=-pox
      sstr(2,1)=-poy
      sstr(3,1)=sgnsil*(pip-pim)/2
      sstr(4,1)=(pip+pim)/2
      sstr(1,2)=pux
      sstr(2,2)=puy
      sstr(3,2)=sgnsil*(pap-pam)/2
      sstr(4,2)=(pap+pam)/2
      pnl3=sstr(3,1)+sstr(3,2)
      pnl4=sstr(4,1)+sstr(4,2)
      ynl=100
      if(pnl4-dabs(pnl3).ne.0.)
     *ynl=0.5*dlog((pnl4+dabs(pnl3))/(pnl4-dabs(pnl3)))
      do n=1,4
      do i=1,2
      str(n,i)=sstr(n,i)
      enddo 
      enddo 

c     checks
c     ------
      if(ish.ge.90)then
      if(dabs(pip*pim-pot**2).gt.1e-4)then
      call utmsg('jsplit')
      write(ifch,*)'*****  pip*pim /= pot**2'
      write(ifch,*)'pip*pim=',pip*pim
      write(ifch,*)'pot**2=',pot**2
      write(ifch,*)'pip=',pip
      write(ifch,*)'pim=',pim
      write(ifch,*)'pot=',pot
      call utmsgf
      endif
      if(dabs(pap*pam-put**2).gt.1e-4)then
      call utmsg('jsplit')
      write(ifch,*)'*****  pap*pam /= put**2'
      write(ifch,*)'pap*pam=',pap*pam
      write(ifch,*)'put**2=',put**2
      write(ifch,*)'pap=',pap
      write(ifch,*)'pam=',pam
      write(ifch,*)'put=',put
      call utmsgf
      endif
      if(dabs(sstr(4,1)**2
     *-sstr(1,1)**2-sstr(2,1)**2-sstr(3,1)**2).gt.1e-4)then
      call utmsg('jsplit')
      write(ifch,*)'*****  mass**2 of string end 1 nonzero'
      write(ifch,*)'mass**2=',sstr(4,1)**2
     *-sstr(1,1)**2-sstr(2,1)**2-sstr(3,1)**2
      call utmsgf
      endif
      if(dabs(sstr(4,2)**2
     *-sstr(1,2)**2-sstr(2,2)**2-sstr(3,2)**2).gt.1e-4)then
      call utmsg('jsplit')
      write(ifch,*)'*****  mass**2 of string end 2 nonzero'
      write(ifch,*)'mass**2=',sstr(4,2)**2
     *-sstr(1,2)**2-sstr(2,2)**2-sstr(3,2)**2
      call utmsgf
      endif
      do 14 n=1,4
      if(abs(str(n,1)+stro(n,1)-stro0(n,1)
     *+str(n,2)+stro(n,kmax)-stro0(n,kmax)).gt.1e-4
     *.and.abs(str(n,1)+stro(n,1)-stro0(n,1)
     *+str(n,2)+stro(n,kmax)-stro0(n,kmax))
     *.gt.1e-4*abs(stro0(n,1)+stro0(n,kmax)))goto15
14    continue
      goto16
15    continue
      call utmsg('jsplit')
      write(ifch,*)'*****  p_str + p_stro /= p_stro0'
      write(ifch,*)'p_str + p_stro:'
      write(ifch,104)((str(n,1)+str(n,2)+stro(n,1)+stro(n,kmax)),n=1,4)
      write(ifch,*)'p_stro0:'
      write(ifch,104)((stro0(n,1)+stro0(n,kmax)),n=1,4)
      write(ifch,*)'str:'
      write(ifch,104)(str(n,1),n=1,4)
      write(ifch,104)(str(n,2),n=1,4)
      write(ifch,*)'stro:'
      write(ifch,104)(stro(n,1),n=1,4)
      write(ifch,104)(stro(n,kmax),n=1,4)
      write(ifch,*)'stro0:'
      write(ifch,104)(stro0(n,1),n=1,4)
      write(ifch,104)(stro0(n,kmax),n=1,4)
      call utmsgf
16    continue
      if(ish.ge.93)then
      write(ifch,100)apart,epart,en,pt
100   format(3x,'apart,epart,en,pt:',4f11.5)
      write(ifch,101)(sngl(po(k)),k=1,5)
101   format(3x,'po:',5f11.5/)
      if(nqu.eq.2)then
      write(ifch,*)'diquark-antidiquark break'
      write(ifch,*)' '
      endif
      endif
      endif
 
c     flavour
c     -------
      do 1 i=1,2
1     ic(i)=nint(abs(stro(4+i,kmax)))
      do 3 i=1,2
3     str(4+i,2)=ic(i)
      do 4 n=1,nflav
      do 4 i=1,2
4     jcx(n,i)=0
      m=0
      if(ic(1).gt.0.and.nqu.eq.1)m=2
      if(ic(2).gt.0.and.nqu.eq.1)m=1
      if(ic(1).gt.0.and.nqu.eq.2)m=1
      if(ic(2).gt.0.and.nqu.eq.2)m=2
      if(m.eq.0)
     *call utstop('jsplit: m = 0                           ')
      do 25 n=1,nqu
      if(n.eq.1)nfl=mod(iflto,10)
      if(n.eq.2)nfl=iflto/10
25    jcx(nfl,m)=jcx(nfl,m)+1
      call idenco(jcx,icx,ireten)
      if(ireten.eq.1)
     *call utstop('jsplit: idenco ret code = 1             ')
      do 5 i=1,2
5     str(4+i,1)=icx(i)
           call utamst(str,am,amin,iret)
           if(iret.ne.0)goto1002
           if(nfusn(kmax).eq.1)goto1002
      do 6 i=1,2
6     ic(i)=nint(stro(4+i,1))
      call iddeco(ic,jc)
      do 26 n=1,nqu
      if(n.eq.1)nfl=mod(iflto,10)
      if(n.eq.2)nfl=iflto/10
      if(jc(nfl,m).gt.0)then
      jc(nfl,m)=jc(nfl,m)-1
      else
      jc(nfl,3-m)=jc(nfl,3-m)+1
      endif
26    continue
      nn=0
      do 27 n=1,nflav
      do 27 i=1,2
27    nn=nn+jc(n,i)
      if(nn.eq.0)then
      nfl=int(rangen()/pudx)+1
      jc(nfl,1)=1
      jc(nfl,2)=1
      endif
      call idenco(jc,ic,ireten)
           if(ireten.eq.1)then
      if(ish.ge.90.and.messco.eq.0)then
      call utmsg('jsplit')
      write(ifch,*)'*****  idenco ret code = 1.   redo jsplit'
      write(ifch,*)'jc:'
      write(ifch,*)jc
      call utmsgf
      messco=1
      endif
      goto5000
           endif
      do 7 i=1,2
      stro(4+i,1)=ic(i)
7     stro(4+i,kmax)=0.
 
c     ok
c     --
1000  ier=0
      nscc=nscc+1
      goto10002
 
c     error
c     -----
1001  ier=1
      if(ish.ge.90)then
      call utmsg('jsplit')
      write(ifch,*)'*****  split not possible'
      call utmsgf
      endif
      goto10001
 
c     absorption
c     ----------
1002  ier=2
      nfusn(kmax)=1
      if(ish.ge.92)then
      write(ifch,*)'absorption of stro(,kmax)'
      write(ifch,*)' '
      endif
      do 18 i=1,nsi
      str(i,1)=0
18    str(i,2)=0
      do 17 i=1,nsi
      do 17 k=1,kmax
17    stro(i,k)=stro0(i,k)
      do 19 i=1,4
      stro(i,1)=stro(i,1)+stro(i,kmax)
19    stro(i,kmax)=0.
      do 20 i=1,2
20    ic(i)=nint(stro(4+i,1))
      call iddeco(ic,jc)
      do 21 i=1,2
21    icx(i)=nint(abs(stro(4+i,kmax)))
      call iddeco(icx,jcx)
      do 22 nf=1,nflav
      do 22 i=1,2
22    jc(nf,i)=jc(nf,i)+jcx(nf,i)
      call idenco(jc,ic,ireten)
           if(ireten.eq.1)then
      if(ish.ge.90)then
      call utmsg('jsplit')
      write(ifch,*)'*****  idenco ret code = 1'
     *,'   (after absorption)'
      write(ifch,*)'jc:'
      write(ifch,*)jc
      call utmsgf
      endif
      goto1001
           endif
      do 23 i=1,2
      stro(4+i,1)=ic(i)
23    stro(4+i,kmax)=0.
      goto10002
 
c     check nscc
c     ----------
10002 continue
      if(kmaxor.gt.2.and.kmax.eq.2.and.max(1,nscc).ne.nsccx)then
      ier=3
      if(ish.ge.91)then
      write(ifch,*)'redo string processing with new apart, epart'
      write(ifch,*)' '
      endif
      goto10001
      endif
 
c     print
c     -----
      if(ish.ge.92)then
      if(ish.ge.93)then
      write(ifch,*)'remainder string:'
      write(ifch,*)' '
      write(ifch,105)(stro(i,1),i=1,4),(nint(stro(i,1)),i=5,6)
      do 9 k=2,kmax
9     write(ifch,104)(stro(i,k),i=1,4),(nint(stro(i,k)),i=5,6)
      write(ifch,*)' '
      write(ifch,*)'split string:'
      write(ifch,*)' '
      endif
      write(ifch,105)(str(i,1),i=1,4),(nint(str(i,1)),i=5,6)
      write(ifch,104)(str(i,2),i=1,4),(nint(str(i,2)),i=5,6)
105   format(' str: ',4f11.5,2i8)
104   format('      ',4f11.5,2i8)
      write(ifch,*)' '
      endif
10001 continue
      if(ish.ge.93)then
      write(ifch,*)('-',l=1,25)
      write(ifch,*)'   return from jsplit   '
      write(ifch,*)('-',l=1,25)
      write(ifch,*)' '
      endif
      ish=ish0
      return
 
      end
 
c-----------------------------------------------------------------------
      subroutine jstrij(ipa1,ipa2,jc,idhad,jcx,ier,isame,ierrc)
c     strips hadron idhad off jet jc. isame=1: zero flavor required
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (nflav=6)
      integer jc(nflav,2),ic(2),jcx(nflav,2),icx(2),icp(2),jcp(nflav,2)
     *,ic1(2),jc1(nflav,2),ic2(2),jc2(nflav,2)
      ier=0
      call randfl(jc,ipa1,ii,ic1,isame)
      if(ii.eq.0)goto9999
      call iddeco(ic1,jc1)
      do 8 n=1,nflav
      do 8 i=1,2
      jcp(n,i)=jc1(n,i)
8     jcx(n,i)=jc(n,i)-jc1(n,i)
      if(isame.ne.1)then
      call randfl(jcx,ipa2,ii,ic2,0)
      if(ii.eq.0)goto9999
      call iddeco(ic2,jc2)
      do 4 n=1,nflav
      do 4 i=1,2
      jcp(n,i)=jcp(n,i)+jc2(n,i)
4     jcx(n,i)=jcx(n,i)-jc2(n,i)
      endif
      call idenco(jcp,icp,ireten)
      if(ireten.eq.1)
     *call utstop('jstrij: idenco ret code = 1             ')
      idhad=idtra(icp,1,0,3)
      if(idhad.eq.0)goto9998
      return
9998  jerr=jerr+1
      if(ish.ge.91)write(ifch,*)' '
      write(ifch,*)'***** error in jstrij: idhad = 0'
      write(ifch,*)'code bef strip:'
      write(ifch,*)jc
      write(ifch,*)'jc1:'
      write(ifch,*)jc1
      write(ifch,*)'jc2:'
      write(ifch,*)jc2
      write(ifch,*)'ptl code:'
      write(ifch,*)jcp
      if(ish.ge.91)write(ifch,*)' '
      ier=1
      return
9999  ier=1
      if(ierrc.eq.0)return
      jerr=jerr+1
      if(ish.ge.91)write(ifch,*)' '
      write(ifch,*)'***** error in jstrij: no parton found'
      write(ifch,*)'code bef strip:'
      write(ifch,*)jc
      if(ish.ge.91)write(ifch,*)' '
      return
      end
 
c-----------------------------------------------------------------------
      subroutine jstrip(ipa1,ipa2,ic,idhad,icx,ier,isame,ierrc)
c     strips hadron idhad off jet ic. isame=1: zero flavor required
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2),jcx(nflav,2),icx(2)
      call iddeco(ic,jc)
      call jstrij(ipa1,ipa2,jc,idhad,jcx,ier,isame,ierrc)
      call idenco(jcx,icx,ireten)
      if(ireten.eq.1)
     *call utstop('jstrip: idenco ret code = 1             ')
      return
      end
 
c-----------------------------------------------------------------------
      subroutine lepexp(rxbj,rqsq)
c     generates x_bjorken and q**2 according to an experimental
c     distribution ( given in array xq(nxbj,nqsq) ).
c-----------------------------------------------------------------------
      parameter (nxbj=10,nqsq=10)
      parameter (xbjmin=0.,qsqmin=4.)
      parameter (xbjwid=0.025, qsqwid=4.)
      dimension xq(nxbj,nqsq),vxq(nxbj*nqsq)
      equivalence (xq(1,1),vxq(1))
 
      data vxq/
     &         1304.02,   366.40,    19.84,    10.79,     6.42,
     &            4.54,     4.15,     3.38,     2.03,     1.56,
     &          241.63,  1637.26,   427.36,   164.51,    73.72,
     &           43.07,    20.73,    12.78,     9.34,     5.83,
     &            0.01,   724.66,   563.79,   275.08,   176.13,
     &          106.44,    85.82,    54.52,    37.12,    28.65,
     &            0.01,   202.40,   491.10,   245.13,   157.07,
     &          104.43,    61.05,    49.42,    37.84,    26.79,
     &            0.01,     3.77,   316.38,   226.92,   133.45,
     &           90.30,    63.67,    48.42,    35.73,    28.04,
     &            0.01,     0.01,   153.74,   213.09,   114.14,
     &           76.26,    60.02,    43.15,    43.47,    25.60,
     &            0.01,     0.01,    39.31,   185.74,   108.56,
     &           88.40,    47.29,    39.35,    31.80,    22.91,
     &            0.01,     0.01,     0.01,   104.61,   107.01,
     &           66.24,    45.34,    37.45,    33.44,    23.78,
     &            0.01,     0.01,     0.01,    56.58,    99.39,
     &           67.78,    43.28,    35.98,    34.63,    18.31,
     &            0.01,     0.01,     0.01,    13.56,    76.25,
     &           64.30,    42.80,    28.56,    21.19,    20.75 /
 
      data init/0/
      init=init+1
      if(init.eq.1) then
      n=nxbj*nqsq
      sum=0.
      do 1 i=1,n
      sum=sum+vxq(i)
1     continue
      do 2 i=2,n
2     vxq(i)=vxq(i)+vxq(i-1)
      do 3 i=1,n
3     vxq(i)=vxq(i)/sum
      endif
 
      n=nxbj*nqsq
      r=rangen()
      call utloc(vxq,n,r,iloc)
      if(iloc.ge.n) iloc=iloc-1
      i=mod(iloc,nxbj)+1
      if(i.eq.0) i=nxbj
      j=iloc/nxbj + 1
      dxint=vxq(1)
      if(iloc.gt.0) dxint=vxq(iloc+1)-vxq(iloc)
      dxbj=xbjwid*abs(r-vxq(iloc+1))/dxint
      dy=qsqwid*rangen()
      rxbj=xbjmin+xbjwid*float(i-1)+dxbj
      rqsq=qsqmin+qsqwid*float(j-1)+dy
      return
      end
 
c-----------------------------------------------------------------------
      subroutine lepstr(nuclon,xbj,qsq,nstrng)
c     returns string code nstrng for given nuclon, xbj, qsq.
c     nuclon  : the struck nucleon (1120/1220 = proton/neutron)
c     xbj,qsq : x-bjorken and q**2
c     nstrng  : string code:
c        proton strings:          neutron strings:
c        1 : u  ---  ud           11 : u  ---  dd
c        2 : d  ---  uu           12 : d  ---  ud
c        3 : u  ---  ub(uud)      13 : u  ---  ub(udd)
c        4 : d  ---  db(uud)      14 : d  ---  db(udd)
c        5 : s  ---  sb(uud)      15 : s  ---  sb(udd)
c        6 : ub ---  u(uud)       16 : ub ---  u(udd)
c        7 : db ---  d(uud)       17 : db ---  d(udd)
c        8 : sb ---  s(uud)       18 : sb ---  s(udd)
c-----------------------------------------------------------------------
      dimension quarks(9)
      nstrng=0
      if(nuclon.ne.1120.and.nuclon.ne.1220) return
 
c     proton-string (valence part)
      if(nuclon.eq.1120) then
      quarks(1) =  4.*  stxu(xbj,qsq)
      quarks(2) =       stxd(xbj,qsq)
      endif
 
c     neutron-string (valence part)
      if(nuclon.eq.1220) then
      quarks(1) =  4.*  stxd(xbj,qsq)
      quarks(2) =       stxu(xbj,qsq)
      endif
 
c     the sea contributions (proton/neutron)
      udsea     =       stxus(xbj,qsq)
      ssea      =       stxs(xbj,qsq)
      quarks(3) =  4.*  udsea
      quarks(4) =       udsea
      quarks(5) =       ssea
      quarks(6) =  4.*  udsea
      quarks(7) =       udsea
      quarks(8) =       ssea
 
      quarks(9) = 0.
      do 11 i=1,8
11    quarks(9) = quarks(9)+quarks(i)
 
      r=rangen() * quarks(9)
      sumq=0.
      do 12 i=1,8
      nstrng=i
      sumq=sumq+quarks(i)
      if(r.le.sumq) go to 13
12    continue
13    continue
 
      if(nuclon.eq.1220) nstrng=nstrng+10
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine leptar(xbj,qsq,matarg,latarg,nuclon)
c     returns nuclon = id of hit target nucleon.
c     xbj   : x bjorken
c     qsq   : q squared
c     matarg : a  of target
c     latarg : z  of target
c     nuclon: id of target nucleon (1120/1220 for proton/neutron)
c-----------------------------------------------------------------------
c     proton part of the structure function:
      f2p = latarg * stxzpr(xbj,qsq)
c     neutron part of the structure function:
      f2n = (matarg-latarg) * stxzne(xbj,qsq)
c     structure function of nucleus:
      f2a = f2p + f2n
c     select  the target-nucleon  ( proton or neutron ):
      rn = rangen()
      f2rndm = f2a * rn
      nuclon = 1120
      if(f2rndm .lt. f2n) nuclon = 1220
      return
      end
 
c-----------------------------------------------------------------------
      subroutine nucini(opt,nuc,la,ma,isi)
c     initializes nucleon-momenta.
c     writes nucleons on /cptl/ (p,ifr,icl).
c-----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter (nsi=6,nsix=100,mamx=250)
      parameter (nha=nsix+2)
      real nuc(nsi,nha,mamx)
      parameter (mamx2=mamx*2)
      common/cnnn/nnnptl(mamx2)
      character opt*3
      nhax=3
      if(nha.lt.nhax)
     *call utstop('nucini: nha too small                   ')
      las=0
      mas=0
 
      do 1 l=1,ma
 
      do 2 m=1,nha
      do 2 n=1,nsi
2     nuc(n,m,l)=0
 
      if(opt.eq.'nuc')then
      if(isi.gt.0.and.laproj.lt.0)then
      id=idproj
      elseif(isi.lt.0.and.latarg.lt.0)then
      id=idtarg
      else
      id=1220
      if(rangen().le.(la-las)*1./(ma-mas))id=1120
      if(id.eq.1120)las=las+1
      mas=mas+1
      endif
      call idmass(id,ams)
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      p1=0
      p2=0
      p3=isi*pnllx
      p0=sqrt(pnllx**2+ams**2)
      p0x=pnllx
      p5=ams
 
      elseif(opt.eq.'str')then
      id=1220
      if(rangen().le.(la-las)*1./(ma-mas))id=1120
      if(id.eq.1120)las=las+1
      mas=mas+1
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      p1=0
      p2=0
      p3=0
      p0=prom
      p5=prom
      endif
 
      nuc(5,1,l)=ic1
      nuc(6,1,l)=ic2
      nuc(3,2,l)=p3
      nuc(4,2,l)=p0x
      nuc(5,2,l)=ic1
      nuc(6,2,l)=ic2
      nptl=nptl+1
      if(nptl.ne.nnnptl(nptl))
     *call utstop('nucini: nptl and nnnptl dont match      ')
      idptl(nptl)=id
      pptl(1,nptl)=p1
      pptl(2,nptl)=p2
      pptl(3,nptl)=p3
      pptl(4,nptl)=p0
      pptl(5,nptl)=p5
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      iclptl(nptl)=1
      nqjptl(nptl)=0
1     continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine nuclco(massnr,n,x,y,z,ynuc)
c     calculates coordinates of the nucleons in a nucleus.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common /col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      real x(n),y(n),z(n)
      if(massnr.eq.0)return
      if(massnr.gt.n)
     *call utstop('nuclco: massnr.gt.n                     ')
      if(massnr.eq.1)then
      x(1)=0
      y(1)=0
      z(1)=0
      return
      endif
      do 3 i=1,massnr
10    continue
      if(massnr.eq.maproj)r=utinvt(ndep,xdep,qdep,rangen()*qdep(ndep))
      if(massnr.eq.matarg)r=utinvt(ndet,xdet,qdet,rangen()*qdet(ndet))
      if(massnr.ne.maproj.and.massnr.ne.matarg)
     *call utstop('nuclco: nucleus neither proj nor targ   ')
      theta=acos(1-2*rangen())
      phi=-pi+2*pi*rangen()
      x(i)=r*sin(theta)*cos(phi)
      y(i)=r*sin(theta)*sin(phi)
      z(i)=r*cos(theta)
      if(i.eq.1)goto3
      if(core.eq.0.)goto3
      do 2 j=1,i-1
      if((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2.lt.core**2)goto10
2     continue
3     continue
      if(ish.ge.93)write(ifch,*)' '
      if(ish.ge.93)write(ifch,*)'nucleon coordinates:'
      do 4 i=1,massnr
      z(i)=z(i)/cosh(ynuc)
      if(ish.ge.93)write(ifch,*)'i x y z: ',i,x(i),y(i),z(i)
4     continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine nucoge
c     writes nucleons on /cptl/ (xor,tiv,ist,ior,jor,nst).
c     determines sequence of collisions according to nuclear geometry
c       if maproj>0.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (ndep=129,ndet=129,kollmx=5000)
      common /col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      parameter (mamx=250)
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      parameter (mamx2=mamx*2)
      common/cnnn/nnnptl(mamx2)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      if(ish.eq.17)write(ifch,*)' '
      if(ish.eq.17)write(ifch,*)'nucoge (entry)'
      vel=tanh(ypjtl-yhaha)+tanh(yhaha)
      nproj=1
      ntarg=1
      do 15 n=1,kollmx
      nord(n)=n
      do 16 k=1,4
16    coord(k,n)=0
15    continue
 
      if(matarg.lt.0)then
      koll=-matarg
      ntarg=koll
      bimp=0
      nptl=nptl+1
      do 3 i=1,4
3     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
      nnnptl(nptl)=nptl
      do 1 k=1,koll
      nrproj(k)=1
      nrtarg(k)=k
      nptl=nptl+1
      do 4 i=1,4
4     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
      nnnptl(nptl)=nptl
1     continue
      goto9999
 
      elseif(maproj.eq.1.and.matarg.eq.1)then
      koll=1
      bimp=0
      nrproj(1)=1
      nrtarg(1)=1
      do 5 ii=1,2
      nptl=nptl+1
      do 2 i=1,4
2     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
5     nnnptl(nptl)=nptl
      goto9999
      endif
 
      call nuclco(maproj,mamx,xproj,yproj,zproj,ypjtl-yhaha)
      call nuclco(matarg,mamx,xtarg,ytarg,ztarg,yhaha)
      bx=0
      by=0
      if(maproj.gt.0)then
      if(bimevt.lt.0)then
      b1=bminim
      b2=amin1(rmproj+rmtarg,bmaxim)
      if(b1.gt.b2)
     *call utstop('nucoge: bmin > bmax                     ')
      bimp=sqrt(b1**2+(b2**2-b1**2)*rangen())
      phi=rangen()*2*pi
      else
      phi=phievt
      bimp=bimevt
      endif
      bx=cos(phi)*bimp
      by=sin(phi)*bimp
      do 6 i=1,maproj
      nptl=nptl+1
      xorptl(1,nptl)=xproj(i)+bx/2
      xorptl(2,nptl)=yproj(i)+by/2
      xorptl(3,nptl)=zproj(i)
      xorptl(4,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)= ainfin
      istptl(nptl)=0
      iorptl(nptl)=0
      jorptl(nptl)=0
6     nnnptl(nptl)=nptl
      endif
      do 7 i=1,matarg
      nptl=nptl+1
      xorptl(1,nptl)=xtarg(i)-bx/2
      xorptl(2,nptl)=ytarg(i)-by/2
      xorptl(3,nptl)=ztarg(i)
      xorptl(4,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)= ainfin
      istptl(nptl)=0
      iorptl(nptl)=0
      jorptl(nptl)=0
7     nnnptl(nptl)=nptl
      if(maproj.eq.0)goto1000
      koll=0
      do 12 i=1,maproj
      do 12 j=1,matarg
      if((xproj(i)+bx-xtarg(j))**2
     *+(yproj(i)+by-ytarg(j))**2.gt.sigppi/pi)goto12
           if(koll.gt.0)then
      kp=0
      kt=0
      do kl=1,koll
      if(nrproj(kl).eq.i)kp=1
      if(nrtarg(kl).eq.j)kt=1
      enddo
      ftrsig=1.0
      ftr=1.0
      if(kp.eq.1)ftr=ftr*ftrsig
      if(kt.eq.1)ftr=ftr*ftrsig
      if(kp.eq.1.or.kt.eq.1)then
      if((xproj(i)+bx-xtarg(j))**2
     *+(yproj(i)+by-ytarg(j))**2.gt.ftr*sigppi/pi)goto12
      endif
           endif
      koll=koll+1
      if(koll.gt.kollmx)
     *call utstop('nucoge: kollmx too small                ')
      nrproj(koll)=i
      nrtarg(koll)=j
      distce(koll)=ztarg(j)-zproj(i)
      coord(1,koll)=(xproj(i)+xtarg(j))/2
      coord(2,koll)=(yproj(i)+ytarg(j))/2
      coord(3,koll)=(zproj(i)+ztarg(j))/2
      coord(4,koll)=distce(koll)/vel
      istptl(i)=1
      iorptl(i)=-1
      tivptl(2,i)=coord(4,koll)
      istptl(maproj+j)=1
      iorptl(maproj+j)=-1
      tivptl(2,maproj+j)=coord(4,koll)
12    continue
      if(koll.le.1)goto9999
 
      do 21 n=2,koll
      do 22 m=1,n-1
      if(nrproj(m).eq.nrproj(n))goto21
22    continue
      nproj=nproj+1
21    continue
      do 23 n=2,koll
      do 24 m=1,n-1
      if(nrtarg(m).eq.nrtarg(n))goto23
24    continue
      ntarg=ntarg+1
23    continue
 
      do 20 n=1,koll-1
      do 20 m=n+1,koll
      if(distce(nord(m)).lt.distce(nord(n)))then
      nordm=nord(m)
      nord(m)=nord(n)
      nord(n)=nordm
      endif
20    continue
 
9999  continue
      if(koll.le.0)goto1000
      if(koll.lt.ko1ko2/10000.or.koll.gt.mod(ko1ko2,10000))goto1000
      nevt=1
      bimevt=bimp
      phievt=phi
      kolevt=koll
      npjevt=nproj
      ntgevt=ntarg
      pmxevt=pnll
      egyevt=engy
 
1000  continue
      if(ish.eq.17)write(ifch,*)'nucoge (exit)'
      return
      end
 
c-----------------------------------------------------------------------
      subroutine nucoll
c     performs a  nucleus-nucleus collision (incl. nucleon-nucleon)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      real kolbad(kollmx)
      common/civ/ivp0,ivp,ivt0,ivt /cntevm/ntevm
     */cnew/kotri,newico,newcol
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      real proj(nsi,nha),targ(nsi,nha)
     *,proj2(nsi,nha),targ2(nsi,nha)
      parameter(nprbms=20)
      common /cprbms/prbms(nprbms)
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      common/partnr/summas,iptnr,ic4,px4,py4,js4,pex,pey,pezet,pe0,nps
     */cnce/nces,ncolex /cnfr/nrfra /ckol/kol
      parameter (mamx=250)
      real proja(nsi,nha,mamx),targa(nsi,nha,mamx)
      integer ncopj(mamx),ncotg(mamx),ic4(2)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      character dash*1
      dash='-'
      ish0=ish
      if(ishsub/100.eq.6)ish=mod(ishsub,100)
      if(ish.eq.14)call uttima('*** nucoll *** ')
      if(ish.eq.17)write(ifch,*)' '
      if(ish.eq.17)write(ifch,*)'nucoll (entry)'
 
      if(ishsub/100.eq.7)ish=mod(ishsub,100)
      nevt=0
      nptl=0
      call nucoge
      if(ichoic.eq.3.or.koll.eq.0)goto1000
      naevt=naevt+1
      if(koll.lt.ko1ko2/10000.or.koll.gt.mod(ko1ko2,10000))goto1000
      colevt=koll
      amsev=sqrt((npjevt*sqrt(amproj**2+pnllx**2)
     *+ntgevt*sqrt(amtarg**2+pnllx**2))**2
     *-((npjevt-ntgevt)*pnllx)**2)
      if(ishsub/100.eq.7)ish=ish0
 
      itry=0
38    continue
      itry=itry+1
      nstr=0
      nptl=0
      sumpt2=0
      amsac=0
      call nucini('nuc',proja,laproj,maproj,1)
      call nucini('nuc',targa,latarg,iabs(matarg),-1)
 
      kolran=rangen()*koll+1
      kolran=min0(kolran,koll)
      kolran=max0(kolran,1)
 
c--------------------------
      do 2 kols=1,koll
c--------------------------
 
      kol=kols
      iskip=0
      i=nrproj(nord(kol))
      j=nrtarg(nord(kol))
      iproj=i
      itarg=j
      nrecop=0
      nrecot=0
      if(kol.lt.koll)then
      do 33 k=kol+1,koll
      if(nrproj(nord(k)).eq.i)nrecop=nrecop+1
      if(nrtarg(nord(k)).eq.j)nrecot=nrecot+1
33    continue
      endif
      ncop=0
      ncot=0
      if(kol.gt.1)then
      do 34 k=1,kol-1
      if(nrproj(nord(k)).eq.i)ncop=ncop+1
      if(nrtarg(nord(k)).eq.j)ncot=ncot+1
34    continue
      endif
      if(ncop.ge.ncolmx.or.ncot.ge.ncolmx)iskip=1
 
      do 56 m=1,nha
      sm=0
      do 57 n=1,nsi
      sm=sm+proja(n,m,i)**2
57    proj(n,m)=proja(n,m,i)
      if(m.ge.3.and.sm.lt.1e-5)goto58
56    continue
58    continue
      do 53 m=1,nha
      sm=0
      do 54 n=1,nsi
      sm=sm+targa(n,m,j)**2
54    targ(n,m)=targa(n,m,j)
      if(m.ge.3.and.sm.lt.1e-5)goto55
53    continue
55    continue
 
      if(kol.eq.kolran.and.jpsi.eq.1)then
      call pvjpsf(proj,targ)
      iskip=1
      endif
 
c     nr of colour exchanges
c     ----------------------
 
      amsac0=amsac
      iskip0=iskip
      smpt2=sumpt2
      call utremb(proj,targ,1)
      goto4
3     continue
      if(ish.ge.91)write(ifch,*)'redo hh collision'
      call utrest(proj,targ,1)
      sumpt2=smpt2
      iskip=iskip0
      amsac=amsac0
4     continue
 
      ncolex=1
      if(iskip.ne.1)then
      lo=0
16    lo=lo+1
      if(lo.eq.3)then
      if(ish.ge.90)then
      call utmsg('nucoll')
      write(ifch,*)'*****  lo=3'
      call utmsgf
      endif
      endif
      r=rangen()
      ncolex=0
15    ncolex=ncolex+1
      if(ncolex.gt.nprbms)goto16
      if(r.gt.prbms(ncolex))goto15
      endif
 
c     hadron-hadron collision
c     -----------------------
 
      do 31 nce=1,ncolex
      nces=nce
      call utpage
      if(ish.ge.91)then
      write(ifch,101)(dash,l=1,79),ipage,kol,nce
     *,nrproj(nord(kol)),nrtarg(nord(kol)),(dash,l=1,79)
101   format(/1x,79a1/1x,i5,'.page             col: ',i2
     *,'   cex: ',i2,'   prj: ',i3,'   trg: ',i3/1x,79a1/)
      endif
      istorp=0
      istort=0
      if(nrecop.eq.0)istorp=1
      if(nrecot.eq.0)istort=1
      call hahabs(proj,targ
     *,istorp*(nce/ncolex),istort*(nce/ncolex),iskip,irethh)
      if(iskip.ge.2)goto9997
      if(irethh.eq.1)goto3
      if(irescl.eq.1.and.amsac.gt.amsev)goto9998
      sumpt2=sumpt2+
     *proj(1,2)**2+proj(2,2)**2+targ(1,2)**2+targ(2,2)**2
      if(ish.eq.11)write(ifch,*)'sumpt2:',nrevt,kol,nce,sumpt2
      iskip=0
31    continue
 
      do 40 m=1,nha
      sm=0
      do 41 n=1,nsi
      sm=sm+proj(n,m)**2
41    proja(n,m,i)=proj(n,m)
      if(m.ge.3.and.sm.lt.1e-5)goto42
40    continue
42    continue
      do 43 m=1,nha
      sm=0
      do 44 n=1,nsi
      sm=sm+targ(n,m)**2
44    targa(n,m,j)=targ(n,m)
      if(m.ge.3.and.sm.lt.1e-5)goto45
43    continue
45    continue
 
c------------------
  2   continue
c------------------
 
      if(sumpt2.lt.1e-5)goto9999
 
1000  continue
      if(ish.eq.14)call uttima('    nucoll f   ')
      if(ish.eq.17)write(ifch,*)'nucoll (exit)'
      ish=ish0
      return
 
9999  inoiac=inoiac+1
      if(ish.ge.91
     *.or.ish.eq.11)then
      call utmsg('nucoll')
      write(ifch,*)'*****  no interaction. redo nucoll'
      call utmsgf
      endif
      goto38
 
9998  ilamas=ilamas+1
      if(ish.ge.91
     *.or.ish.eq.12)then
      call utmsg('nucoll')
      write(ifch,*)'*****  amsac>amsev: ',amsac,amsev
     *,' . redo nucoll'
      do 1 ns=1,nstr
1     write(ifch,109)ns,(icstr(l,ns)/100,l=1,4)
     *,(pstr(l,ns),l=3,5)
109   format(1x,i3,3x,4i5,2x,3(1x,e9.3))
      call utmsgf
      endif
      goto38
 
9997  continue
      if(ish.ge.90)then
      call utmsg('nucoll')
      write(ifch,*)'*****  iskip>=2. redo nucoll'
      call utmsgf
      endif
      goto38
 
      end
 
c-----------------------------------------------------------------------
      subroutine nucstr(ier)
c     performes x and p trafos for nucleons for string decay in nucleus
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common /clep/icinpu,idscat
      if(matarg.ne.nptl)
     *call utstop('nucstr: matarg /= nptl                  ')
      ncnt=0
3     ncnt=ncnt+1
      if(ncnt.ge.10)goto1001
      n0=min0(1+int(rangen()*nptl),nptl)
      if(icinpu.gt.0.and.idptl(n0).ne.idscat)goto3
      istptl(n0)=1
      iorptl(n0)=-1
      tivptl(2,n0)=0
      do 2 n=1,nptl
      do 1 j=1,3
1     xorptl(j,n)=xorptl(j,n)-xorptl(j,n0)
      phi=rangen()*2*pi
      p1=-elepto*sin(angmue)*sin(phi)
      p2=-elepto*sin(angmue)*cos(phi)
      p3=elepti-elepto*cos(angmue)
      call utrota(1,p1,p2,p3,xorptl(1,n),xorptl(2,n),xorptl(3,n))
      p3=sqrt(p1**2+p2**2+p3**2)
      p1=0
      p2=0
      p4=prom+elepti-elepto
      call utlobo(1,p1,p2,p3,p4,sqrt(p4**2-p3**2-p2**2-p1**2)
     *,xorptl(1,n),xorptl(2,n),xorptl(3,n),xorptl(4,n))
2     call utlobo(1,p1,p2,p3,p4,sqrt(p4**2-p3**2-p2**2-p1**2)
     *,pptl(1,n),pptl(2,n),pptl(3,n),pptl(4,n))
 
      ier=0
      return
 
1001  ier=1
      if(ish.ge.90)then
      call utmsg('nucstr')
      write(ifch,*)'*****  idscat not possible ==> redo event.'
      call utmsgf
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine pvjpsf(proj,targ)
c     forms a jpsi
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/cnsta/pi,prom,piom,ainfin
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      parameter (ndep=129,ndet=129,kollmx=5000)
      common/col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      common/ckol/kol
      parameter (nptj=129)
      common /cptj/xptj(nptj),qptj(nptj),wptj(nptj)
      parameter (ngau=129)
      common /cgau/xgau(ngau),qgau(ngau),wgau(ngau)
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      real proj(nsi,nha),targ(nsi,nha),psum(5),pr(4),tg(4)
      parameter (mamx=250)
      common/cproja/iproj,itarg,kproja(nha,mamx),ktarga(nha,mamx)
      pax=0.19
      prox=proj(3,2)
      tarx=targ(3,2)
      loop=0
      do 3 n=1,4
      pr(n)=proj(n,2)
3     tg(n)=targ(n,2)
 
 
      call utpage
      if(ish.ge.91)write(ifch,110)('-',l=1,79),ipage,('-',l=1,79)
110   format(1x,79a1/1x,i5,'.page            '
     *,'jpsi formation'/1x,79a1/)
 
5000  loop=loop+1
      if(loop.gt.100)then
      call utmsg('pvjpsf')
      write(ifch,*)'*****  jpsi formation not possible'
      call utmsgf
      goto1000
      endif
      do 4 n=1,4
      proj(n,2)=pr(n)
4     targ(n,2)=tg(n)
 
c     jpsi momenta
c     ------------
      id=441
      call idmass(id,am)
      s=am**2
      pt=utinvt(nptj,xptj,qptj,rangen()*qptj(nptj))
      phi=2.*pi*rangen()
      px=pt*cos(phi)
      py=pt*sin(phi)
      lo=0
1     lo=lo+1
      if(lo.gt.10)
     *call utstop('pvjpsf: lo > 10                         ')
      z=pax*utinvt(ngau,xgau,qgau,rangen()*qgau(ngau))
      if(z.gt.1.)goto1
      pz=z*pnllx*ransig()
      e=sqrt(s+px**2+py**2+pz**2)
      pp=e+pz
      pm=e-pz
 
c     proj and targ momenta
c     ---------------------
      r=rangen()
      pox=proj(1,2)-r*px
      poy=proj(2,2)-r*py
      pot=sqrt(pox**2+poy**2)
      r=1-r
      pux=targ(1,2)-r*px
      puy=targ(2,2)-r*py
      put=sqrt(pux**2+puy**2)
      pop=proj(4,2)+proj(3,2)
      pom=proj(4,2)-proj(3,2)
      pup=targ(4,2)+targ(3,2)
      pum=targ(4,2)-targ(3,2)
      pep=pop+pup-pp
      pem=pom+pum-pm
      a=(pem*pep-put**2-pot**2)/2.
      d=put*pot
      if(a**2-d**2.lt.0.)goto5000
      pyp=(a+put**2-sqrt(a**2-d**2))/pem
      pym=(a+pot**2-sqrt(a**2-d**2))/pep
      if(pyp.lt.0.)goto5000
      if(pym.lt.0.)goto5000
      pip=pep-pyp
      if(pip.lt.0.)goto5000
      pim=pym
      pap=pyp
      pam=pem-pym
      if(pam.lt.0.)goto5000
      proj(1,2)=pox
      proj(2,2)=poy
      proj(3,2)=(pip-pim)/2
      proj(4,2)=(pip+pim)/2
      targ(1,2)=pux
      targ(2,2)=puy
      targ(3,2)=(pap-pam)/2
      targ(4,2)=(pap+pam)/2
      kproja(2,iproj)=kol
      ktarga(2,itarg)=kol
      if(proj(3,2)*prox.lt.0.)goto5000
      if(targ(3,2)*tarx.lt.0.)goto5000
 
c     jpsi string
c     -----------
      nstr=nstr+1
      if(nstr.gt.mxstr)
     *call utstop('pvjpsf: nstr>mxstr                      ')
      iorstr(nstr)=-kol
      icstr(1,nstr)=000100
      icstr(2,nstr)=0
      icstr(3,nstr)=0
      icstr(4,nstr)=000100
      pstr(1,nstr)=px
      pstr(2,nstr)=py
      pstr(3,nstr)=pz
      pstr(4,nstr)=e
      pstr(5,nstr)=am
      rotstr(1,nstr)=0
      rotstr(2,nstr)=0
      rotstr(3,nstr)=1
      do 2 i=1,4
2     xorstr(i,nstr)=coord(i,kol)
           if(ish.ge.91)then
      j=nstr
      write(ifch,100)j,(icstr(k,j)/100,k=1,4)
     *,pstr(3,j),pstr(4,j),pstr(5,j)
100   format(' /cstr/',i4,3x,4i5,2x,3(1x,e9.3))
           endif
 
c     checks
c     ------
      if(ish.ge.90)then
      if(abs(pip*pim-pot**2).gt.1e-4)then
      call utmsg('pvjpsf')
      write(ifch,*)'*****  pip*pim /= pot**2'
      write(ifch,*)'pip*pim=',pip*pim,'pot**2=',pot**2
      write(ifch,*)'pip=',pip,'   pim=',pim,'   pot=',pot
      call utmsgf
      endif
      if(abs(pap*pam-put**2).gt.1e-4)then
      call utmsg('pvjpsf')
      write(ifch,*)'*****  pap*pam /= put**2'
      write(ifch,*)'pap*pam=',pap*pam,'put**2=',put**2
      write(ifch,*)'pap=',pap,'   pam=',pam,'   put=',put
      call utmsgf
      endif
      if(abs(proj(4,2)**2
     *-proj(1,2)**2-proj(2,2)**2-proj(3,2)**2).gt.1e-4)then
      call utmsg('pvjpsf')
      write(ifch,*)'*****  mass**2 of proj nonzero'
      write(ifch,*)'mass**2=',proj(4,2)**2
     *-proj(1,2)**2-proj(2,2)**2-proj(3,2)**2
      call utmsgf
      endif
      if(abs(targ(4,2)**2
     *-targ(1,2)**2-targ(2,2)**2-targ(3,2)**2).gt.1e-4)then
      call utmsg('pvjpsf')
      write(ifch,*)'*****  mass**2 of targ nonzero'
      write(ifch,*)'mass**2=',targ(4,2)**2
     *-targ(1,2)**2-targ(2,2)**2-targ(3,2)**2
      call utmsgf
      endif
      do 14 n=1,4
      if(abs(pr(n)+tg(n)
     *-proj(n,2)-targ(n,2)-pstr(n,nstr)).gt.1e-4)goto15
14    continue
      goto16
15    continue
      call utmsg('pvjpsf')
      write(ifch,*)'*****  proj + targ /= proj_new + targ_new +jpsi'
      write(ifch,*)'proj,targ:'
      write(ifch,*)pr
      write(ifch,*)tg
      write(ifch,*)'proj_new,targ_new,jpsi:'
      write(ifch,*)(proj(n,2),n=1,4)
      write(ifch,*)(targ(n,2),n=1,4)
      write(ifch,*)(pstr(n,nstr),n=1,4)
      call utmsgf
16    continue
      endif
 
      amsac=amsac+am
 
1000  return
      end
 
c-----------------------------------------------------------------------
      subroutine racpro(typ,qmu,n,acprob)
c     returns the array acprob containing accumulated prob for:
c       exponential or poisson distribution (for typ = exp or poi)
c       multi pomeron cuts acc to gribov (for typ = gri).
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      common/paro4/grigam,grirsq,gridel,grislo,gricel
      common/cnsta/pi,prom,piom,ainfin
      real acprob(n)
      character typ*3
      ish0=ish
      if(ishsub/100.eq.17)ish=mod(ishsub,100)
 
c     gribov-regge
c     ------------
           if(typ.eq.'gri')then
      xi=alog(engy**2)
      sig=8*pi*grigam*exp(xi*gridel)
      z=2*gricel*grigam/(grirsq+grislo*xi)*exp(xi*gridel)
 
c     0-cut-pomeron probability
c     -------------------------
      qz=1/(-z)
      qzh=1/(-z/2)
      fz=0
      fzh=0
      do 10 i=1,20
      qz=qz*(-z)/i
      qzh=qzh*(-z/2)/i
      fz=fz+qz/i
10    fzh=fzh+qzh/i
      reler=abs(qz/i/fz)
      relerh=abs(qzh/i/fzh)
      sigppz=sig*(fzh-fz)
      sigppe=sigppz/gricel
      sigppd=sigppz/gricel*(gricel-1)
 
c     n-cut-pomeron probabilities, n>0
c     --------------------------------
      r=0
      do 21 i=1,n
      if(i.eq.1)q=1
      if(i.gt.1)q=q*z/(i-1)
      r=r+q
      si=sig/i/z*(1-exp(-z)*r)
      if(i.eq.1)acprob(i)=si
      if(i.gt.1)acprob(i)=acprob(i-1)+si
21    continue
      if(sigppi.lt.0.)sigppi=acprob(n)
      do 22 i=1,n
22    acprob(i)=acprob(i)/acprob(n)
 
c     warnings
c     --------
           if(reler.gt.1e-3)then
      call utmsg('racpro')
      write(ifch,*)'*****  reler too large'
      write(ifch,*)'reler,qz/i,fz:',reler,qz/i,fz
      call utmsgf
           endif
           if(relerh.gt.1e-3)then
      call utmsg('racpro')
      write(ifch,*)'*****  relerh too large'
      write(ifch,*)'relerh,qzh/i,fzh:',relerh,qzh/i,fzh
      call utmsgf
           endif
           if(sigppz.lt.0.)then
      call utmsg('racpro')
      write(ifch,*)'*****  negative sigppz'
      write(ifch,*)'sigppz,sig,fzh,fz:',sigppz,sig,fzh,fz
      call utmsgf
           endif
      do 24 i=1,n
      if(acprob(i).lt.0)then
      call utmsg('racpro')
      write(ifch,*)'*****  negative acprob(i)'
      write(ifch,*)'i,acprob(i):',i,acprob(i)
      call utmsgf
      endif
24    continue
      relern=acprob(n)-acprob(n-1)
           if(relern.gt.1e-3)then
      call utmsg('racpro')
      write(ifch,*)'*****  relern too large'
      write(ifch,*)'relern:',relern
      call utmsgf
           endif
 
c     print
c     -----
      if(ish.ge.92)then
100   format(79a1)
      write(ifch,*)' '
      write(ifch,100)('-',ic=1,79)
      write(ifch,*)'   cross sections and',
     *' cut-pomeron weights acc to gribov-regge-theory'
      write(ifch,100)('-',ic=1,79)
      write(ifch,*)'   cms-energy [gev]:',engy
      write(ifch,100)('-',ic=1,79)
      write(ifch,*)' '
      write(ifch,*)'   parameters: '
      write(ifch,*)' '
      write(ifch,*)'     gamma: ',grigam
      write(ifch,*)'     r**2:  ',grirsq
      write(ifch,*)'     delta: ',gridel
      write(ifch,*)'     slope: ',grislo
      write(ifch,*)'     c:     ',gricel
      write(ifch,*)' '
      write(ifch,*)'   cross sections:'
      write(ifch,*)' '
      write(ifch,*)'     elastic:     ',sigppe
      write(ifch,*)'     diffractive: ',sigppd
      write(ifch,*)'     inelastic:   ',sigppi
      write(ifch,*)' '
      write(ifch,*)'   weights w(n) of n cut pomerons:'
      write(ifch,*)'      ( n - w(n) - w(n)_exp )'
      write(ifch,*)' '
      a=qmu/(qmu+1)
      i=1
      write(ifch,*)i,acprob(i),(1-a)*a**(i-1)
      do 25 i=2,n
25    write(ifch,*)i,acprob(i)-acprob(i-1),(1-a)*a**(i-1)
      write(ifch,*)' '
      endif
 
c     poisson
c     -------
           elseif(typ.eq.'poi')then
      z=qmu
      if(z.ge.n-1.)
     *call utstop('racpro: z >= n-1.                       ')
      k=amax1(0.,z-(n-1-z))+1
      if(k.ge.n)
     *call utstop('racpro: k >= n.                         ')
           if(k.eq.1)then
      prbab=exp(-z)
      acprob(1)=prbab
      do 1 i=2,n
      prbab=prbab*z/(i-1)
1     acprob(i)=acprob(i-1)+prbab
           elseif(k.gt.1)then
      x=z*exp(-z/(k-1.))
      prbab=1
      do 3 i=1,k-1
      prbab=prbab*x/i
3     acprob(i)=0
      acprob(k)=prbab
      do 4 i=k+1,n
      prbab=prbab*z/(i-1)
4     acprob(i)=acprob(i-1)+prbab
           else
      call utstop('racpro: k <= 0.                         ')
           endif
 
c     exponential
c     -----------
           elseif(typ.eq.'exp')then
      a=qmu/(qmu+1)
      prbab=1-a
      acprob(1)=prbab
      do 2 i=2,n
      prbab=prbab*a
2     acprob(i)=acprob(i-1)+prbab
           endif
 
      ish=ish0
      return
      end
 
c-----------------------------------------------------------------------
      subroutine randfl(jc,iqa0,iflav,ic,isame)
c     returns random flavour ic(2) (iqa0=1:quark,2:antiquark,11:diquark)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (nflav=6)
      real probab(nflav),probsu(nflav+1)
      integer jc(nflav,2),jc0(nflav,2),ic(2)
      character dash*1
      dash='-'
      if(ish.ge.95)then
      write(ifch,*)(dash,l=1,79)
      write(ifch,*)'subroutine randfl'
      write(ifch,*)'iqa0:',iqa0
      write(ifch,*)'jc:'
      write(ifch,*)jc
      endif
      iflav=0
      ic(1)=0
      ic(2)=0
      do 10 n=1,nflav
      do 10 i=1,2
10    jc0(n,i)=0
      iqa1=iqa0*10
9999  iqa1=iqa1/10
      if(iqa1.eq.0)goto9998
      iqa=mod(iqa1,10)
      su=0
      do 20 i=1,nflav
      probab(i)=jc(i,iqa)-jc0(i,iqa)
      if(isame.eq.1)probab(i)=probab(i)*(jc(i,3-iqa)-jc0(i,3-iqa))
20    su=su+probab(i)
      if(su.lt..5)then
      iflav=0
      ic(1)=0
      ic(2)=0
      goto9998
      endif
      probsu(1)=0.
      do 30 i=1,nflav
      probsu(i+1)=probsu(i)+probab(i)/su
      if(probsu(i+1)-probsu(i).lt.1e-5)probsu(i+1)=probsu(i)
30    continue
      r=rangel(probsu(nflav+1))
      do 50 i=1,nflav
      if(probsu(i).le.r.and.r.lt.probsu(i+1))iflav=i
50    continue
      jc0(iflav,iqa)=jc0(iflav,iqa)+1
      if(isame.eq.1)jc0(iflav,3-iqa)=jc0(iflav,3-iqa)+1
      call idenco(jc0,ic,ireten)
      if(ireten.eq.1)
     *call utstop('randfl: idenco ret code = 1             ')
      if(ish.ge.95)then
      write(ifch,*)'probab:'
      write(ifch,*)probab
      write(ifch,*)'probsu:'
      write(ifch,*)probsu
      write(ifch,*)'ran#:',r,'   flav:',iflav
      endif
      goto9999
9998  continue
      if(ish.ge.95)write(ifch,*)(dash,l=1,79)
      return
      end
 
c-----------------------------------------------------------------------
      function rangen()
c     generates a random number
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
1     rangen=rndm(-1)
      if(rangen.le.0.)goto1
      if(rangen.ge.1.)goto1
      if(ish.eq.17)write(ifch,*)'   rangen()= ',rangen
     *,'    es: ',nrevt+1
      return
      entry rangel(xmax)
2     rangel=rndm(-1)*xmax
      if(rangel.le.0.)goto2
      if(rangel.ge.xmax)goto2
      return
      end
 
c-----------------------------------------------------------------------
      subroutine ranhvy(x,eps)
c     generates x for heavy particle fragmentation according to
c     the peterson form
c          d(x)=1/(x*(1-1/x-eps/(1-x))**2)
c              =d0(x)*d1(x)*d2(x)
c          d0(x)=(1-x)**2/((1-x)**2+eps)**2
c          d1(x)=x
c          d2(x)=(((1-x)**2+eps)/((1-x)**2+eps*x))**2
c     using x=1-y**pow
c     generates flat in x if eps>1.
c-----------------------------------------------------------------------
      data aln4/1.3863/
      if(eps.lt.1.) then
        pow=alog((3.+eps)/eps)/aln4
        ymx=(eps*(3.*pow-1.)/(pow+1.))**(.5/pow)
        zmx=1-ymx**pow
        d0mx=(1-zmx)**2/((1.-zmx)**2+eps)**2*pow*ymx**(pow-1.)
        d2mx=2./(2.-sqrt(eps))
      else
        pow=1.
        zmx=0.
        d0mx=(1.-zmx)**2/((1.-zmx)**2+eps)**2
        d2mx=1.+eps
      endif
c
c          generate z according to (1-z)**2/((1-z)**2+eps*z)**2
1     continue
      y=rangen()
      z=1.-y**pow
c
      d0z=(1.-z)**2/((1.-z)**2+eps)**2*pow*y**(pow-1.)
      if(d0z.lt.rangen()*d0mx) goto1
c
c          check remaining factors
      d1=z
      d2=(((1.-z)**2+eps)/((1.-z)**2+eps*z))**2
      if(d1*d2.lt.rangen()*d2mx) goto1
c
c          good x
      x=z
      return
      end
 
c-----------------------------------------------------------------------
      function ransig()
c     returns randomly +1 or -1
c-----------------------------------------------------------------------
      ransig=1
      if(rangen().gt.0.5)ransig=-1
      return
      end
 
c-----------------------------------------------------------------------
      function ranstc(xfl,xmin)
c     returns random number according to a quark structure fctn
c       with x>=xmin
c-----------------------------------------------------------------------
      parameter (nstru=2049)
      common/stru/xstru(nstru),wstru(nstru)
     *,qvap(nstru),qvat(nstru),qsep(nstru),qset(nstru)
      character xfl*3
      if(xfl.eq.'vap')ranstc=ranxq(nstru,xstru,qvap,xmin)
      if(xfl.eq.'sep')ranstc=ranxq(nstru,xstru,qsep,xmin)
      if(xfl.eq.'vat')ranstc=ranxq(nstru,xstru,qvat,xmin)
      if(xfl.eq.'set')ranstc=ranxq(nstru,xstru,qset,xmin)
      return
      end
 
c-----------------------------------------------------------------------
      function ranxq(n,x,q,xmin)
c     returns random number according to x(i) q(i) with x>=xmin
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      real x(n),q(n)
      imin=1
      if(xmin.eq.0.)goto3
      i1=1
      i2=n
1     i=i1+(i2-i1)/2
      if(x(i).lt.xmin)then
      i1=i
      elseif(x(i).gt.xmin)then
      i2=i
      else
      imin=i
      goto3
      endif
      if(i2-i1.gt.1)goto1
      imin=i2
3     continue
      if(q(imin).gt.q(n)*.9999)then
      ranxq=xmin
      goto4
      endif
      qran=q(imin)+rangen()*(q(n)-q(imin))
      ranxq=utinvt(n,x,q,qran)
4     continue
 
      if(ranxq.lt.xmin)then
      call utmsg('ranxq ')
      write(ifch,*)'*****  ranxq=',ranxq,' <       xmin=',xmin
      write(ifch,*)'q(imin) q q(n):',q(imin),qran,q(n)
      write(ifch,*)'x(imin) x x(n):',x(imin),ranxq,x(n)
      call utmsgf
      ranxq=xmin
      endif
 
      return
      end
 
c-----------------------------------------------------------------------
      function sbet(z,w)
c-----------------------------------------------------------------------
      sbet=sgam(z)*sgam(w)/sgam(z+w)
      return
      end
 
c-----------------------------------------------------------------------
      function sdensi(r)
c     nuclear density
c-----------------------------------------------------------------------
      common /cden/massnr,rmx
      a=0.54
      r0=1.19*massnr**(1./3.)-1.61*massnr**(-1./3.)
      sdensi=r**2/(1+exp((r-r0)/a))
      return
      end
 
c-----------------------------------------------------------------------
      function sgam(x)
c     gamma fctn
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      dimension c(13)
      data c
     1/ 0.00053 96989 58808, 0.00261 93072 82746, 0.02044 96308 23590,
     2  0.07309 48364 14370, 0.27964 36915 78538, 0.55338 76923 85769,
     3  0.99999 99999 99998,-0.00083 27247 08684, 0.00469 86580 79622,
     4  0.02252 38347 47260,-0.17044 79328 74746,-0.05681 03350 86194,
     5  1.13060 33572 86556/
      z=x
      if(x .gt. 0.0) goto1
      if(x .eq. aint(x)) goto5
      z=1.0-z
    1 f=1.0/z
      if(z .le. 1.0) goto4
      f=1.0
    2 continue
      if(z .lt. 2.0) goto3
      z=z-1.0
      f=f*z
      goto2
    3 z=z-1.0
    4 sgam=
     1 f*((((((c(1)*z+c(2))*z+c(3))*z+c(4))*z+c(5))*z+c(6))*z+c(7))/
     2   ((((((c(8)*z+c(9))*z+c(10))*z+c(11))*z+c(12))*z+c(13))*z+1.0)
      if(x .gt. 0.0) return
      sgam=3.141592653589793/(sin(3.141592653589793*x)*sgam)
      return
    5 write(ifch,10)x
   10 format(1x,'argument of gamma function = ',e20.5)
      call utstop('sgam : negative integer argument        ')
      end
 
c-----------------------------------------------------------------------
      function sgau(x)
c     returns gaussian distribution (not normalized)
c-----------------------------------------------------------------------
      sgau=exp(-x**2/2)
      return
      end
 
c----------------------------------------------------------------------
      subroutine sjcen0
c     initializes /csjcen/
c----------------------------------------------------------------------
      common/csjcen/ent(16000)	
      real en(16000)	
      data (en(i),i=     1,   171)/
     *  0.00,  8.90, 13.15, 16.70, 19.94, 22.97, 25.83, 28.56, 31.18,
     * 33.77, 36.25, 38.70, 41.05, 43.39, 45.71, 48.00, 50.17, 52.31,
     * 54.50, 56.64, 58.74, 60.80, 62.82, 64.86, 66.84, 68.85, 70.84,
     * 72.72, 74.72, 76.54, 78.52, 80.37, 82.19, 84.03, 85.90, 87.73,
     * 89.48, 91.34, 93.04, 94.81, 96.61,  0.00,  9.28, 13.30, 16.72,
     * 19.89, 22.86, 25.71, 28.41, 31.06, 33.62, 36.12, 38.53, 40.96,
     * 43.27, 45.57, 47.80, 50.02, 52.18, 54.39, 56.51, 58.59, 60.65,
     * 62.79, 64.75, 66.71, 68.78, 70.72, 72.65, 74.57, 76.46, 78.37,
     * 80.29, 82.15, 84.06, 85.80, 87.60, 89.39, 91.25, 93.01, 94.72,
     * 96.49,  0.00,  9.39, 13.39, 16.84, 20.05, 23.06, 25.90, 28.66,
     * 31.28, 33.88, 36.44, 38.88, 41.25, 43.59, 45.87, 48.21, 50.36,
     * 52.64, 54.79, 56.89, 58.97, 61.07, 63.15, 65.13, 67.18, 69.09,
     * 70.90, 73.07, 74.93, 76.90, 78.75, 80.67, 82.46, 84.31, 86.18,
     * 88.01, 89.77, 91.60, 93.30, 95.13, 96.86,  0.00,  9.14, 13.57,
     * 17.25, 20.58, 23.68, 26.57, 29.40, 32.04, 34.69, 37.20, 39.63,
     * 42.04, 44.39, 46.76, 49.01, 51.25, 53.40, 55.53, 57.68, 59.80,
     * 61.84, 63.89, 65.89, 67.97, 69.90, 71.84, 73.77, 75.79, 77.62,
     * 79.53, 81.38, 83.25, 85.10, 86.98, 88.80, 90.56, 92.38, 94.14,
     * 95.98, 97.65,  0.00, 10.73, 15.45, 19.30, 22.72, 25.89, 28.81/
      data (en(i),i=   172,   342)/
     * 31.64, 34.31, 36.91, 39.42, 41.83, 44.26, 46.58, 48.87, 51.10,
     * 53.30, 55.47, 57.60, 59.68, 61.76, 63.79, 65.87, 67.88, 69.86,
     * 71.79, 73.78, 75.66, 77.54, 79.49, 81.34, 83.21, 85.03, 86.87,
     * 88.71, 90.48, 92.38, 94.09, 95.86, 97.64, 99.38,  0.00, 11.36,
     * 15.90, 19.62, 22.92, 25.99, 28.88, 31.60, 34.23, 36.77, 39.28,
     * 41.68, 44.03, 46.39, 48.56, 50.88, 53.07, 55.17, 57.36, 59.43,
     * 61.54, 63.53, 65.62, 67.55, 69.52, 71.51, 73.45, 75.49, 77.25,
     * 79.18, 81.03, 82.96, 84.79, 86.58, 88.45, 90.17, 92.03, 93.85,
     * 95.50, 97.26, 99.12,  0.00, 11.76, 16.34, 20.00, 23.24, 26.25,
     * 29.09, 31.78, 34.41, 36.93, 39.40, 41.78, 44.11, 46.44, 48.72,
     * 50.90, 53.15, 55.27, 57.40, 59.48, 61.57, 63.60, 65.71, 67.67,
     * 69.58, 71.64, 73.54, 75.43, 77.36, 79.27, 81.11, 82.97, 84.85,
     * 86.59, 88.38, 90.21, 92.06, 93.85, 95.71, 97.39, 99.14,  0.00,
     * 11.95, 16.59, 20.22, 23.43, 26.44, 29.27, 31.98, 34.57, 37.11,
     * 39.58, 41.98, 44.28, 46.59, 48.89, 51.12, 53.31, 55.45, 57.61,
     * 59.63, 61.77, 63.82, 65.84, 67.87, 69.82, 71.80, 73.75, 75.71,
     * 77.53, 79.49, 81.33, 83.19, 85.05, 86.96, 88.64, 90.45, 92.32,
     * 94.06, 95.78, 97.62, 99.37,  0.00, 12.00, 16.65, 20.31, 23.55,
     * 26.55, 29.41, 32.14, 34.75, 37.33, 39.84, 42.21, 44.58, 46.93/
      data (en(i),i=   343,   513)/
     * 49.22, 51.47, 53.65, 55.82, 57.93, 60.02, 62.14, 64.22, 66.21,
     * 68.27, 70.16, 72.11, 74.05, 76.04, 77.95, 79.79, 81.77, 83.53,
     * 85.39, 87.18, 89.11, 90.94, 92.61, 94.44, 96.27, 98.04, 99.73,
     *  0.00, 11.85, 16.47, 20.18, 23.53, 26.65, 29.56, 32.36, 35.06,
     * 37.66, 40.22, 42.66, 45.01, 47.35, 49.71, 51.94, 54.15, 56.31,
     * 58.46, 60.60, 62.66, 64.72, 66.73, 68.77, 70.79, 72.77, 74.71,
     * 76.65, 78.57, 80.46, 82.34, 84.18, 86.10, 87.89, 89.73, 91.49,
     * 93.34, 95.07, 96.81, 98.62,100.33,  0.00, 11.22, 16.20, 20.26,
     * 23.82, 27.06, 30.11, 32.95, 35.68, 38.37, 40.91, 43.41, 45.80,
     * 48.16, 50.49, 52.73, 54.96, 57.16, 59.34, 61.41, 63.56, 65.65,
     * 67.62, 69.72, 71.62, 73.67, 75.62, 77.53, 79.47, 81.37, 83.22,
     * 85.08, 86.91, 88.82, 90.56, 92.35, 94.13, 95.96, 97.72, 99.51,
     *101.20,  0.00, 12.10, 17.25, 21.36, 24.99, 28.30, 31.38, 34.28,
     * 37.04, 39.73, 42.26, 44.77, 47.23, 49.59, 51.88, 54.11, 56.33,
     * 58.54, 60.68, 62.78, 64.91, 66.98, 68.97, 70.98, 72.98, 74.91,
     * 76.87, 78.79, 80.68, 82.64, 84.52, 86.39, 88.14, 89.98, 91.87,
     * 93.67, 95.42, 97.18, 98.98,100.71,102.48,  0.00, 12.82, 17.88,
     * 21.83, 25.35, 28.56, 31.56, 34.39, 37.09, 39.71, 42.25, 44.71,
     * 47.10, 49.49, 51.74, 54.01, 56.17, 58.35, 60.53, 62.61, 64.62/
      data (en(i),i=   514,   684)/
     * 66.73, 68.73, 70.66, 72.71, 74.58, 76.59, 78.56, 80.36, 82.34,
     * 84.14, 86.04, 87.91, 89.66, 91.52, 93.25, 95.04, 96.88, 98.61,
     *100.33,102.11,  0.00, 13.26, 18.39, 22.39, 25.86, 28.99, 31.97,
     * 34.73, 37.44, 40.00, 42.48, 44.93, 47.24, 49.58, 51.89, 54.10,
     * 56.30, 58.40, 60.56, 62.56, 64.66, 66.64, 68.70, 70.75, 72.69,
     * 74.66, 76.62, 78.48, 80.40, 82.19, 84.09, 85.96, 87.77, 89.62,
     * 91.48, 93.20, 94.99, 96.80, 98.50,100.28,102.05,  0.00, 13.54,
     * 18.77, 22.80, 26.23, 29.41, 32.31, 35.08, 37.69, 40.25, 42.68,
     * 45.12, 47.45, 49.74, 52.01, 54.25, 56.41, 58.54, 60.71, 62.78,
     * 64.77, 66.84, 68.88, 70.86, 72.77, 74.73, 76.68, 78.62, 80.52,
     * 82.35, 84.24, 86.02, 87.90, 89.69, 91.46, 93.28, 95.12, 96.81,
     * 98.57,100.30,102.04,  0.00, 13.73, 19.04, 23.09, 26.55, 29.67,
     * 32.57, 35.30, 37.95, 40.49, 42.96, 45.38, 47.70, 49.99, 52.24,
     * 54.37, 56.62, 58.72, 60.85, 62.96, 65.02, 67.02, 69.11, 71.02,
     * 72.99, 74.95, 76.83, 78.73, 80.68, 82.48, 84.35, 86.27, 88.08,
     * 89.94, 91.64, 93.46, 95.15, 97.05, 98.73,100.49,102.23,  0.00,
     * 13.86, 19.17, 23.24, 26.69, 29.84, 32.70, 35.50, 38.09, 40.66,
     * 43.15, 45.52, 47.85, 50.21, 52.44, 54.65, 56.81, 58.98, 61.05,
     * 63.14, 65.21, 67.22, 69.30, 71.23, 73.25, 75.15, 77.02, 79.03/
      data (en(i),i=   685,   855)/
     * 80.96, 82.83, 84.70, 86.45, 88.35, 90.21, 92.00, 93.74, 95.47,
     * 97.21, 99.02,100.83,102.45,  0.00, 13.89, 19.20, 23.29, 26.74,
     * 29.89, 32.79, 35.60, 38.25, 40.82, 43.30, 45.69, 48.09, 50.37,
     * 52.65, 54.86, 57.08, 59.24, 61.34, 63.45, 65.53, 67.52, 69.55,
     * 71.53, 73.61, 75.53, 77.43, 79.39, 81.19, 83.05, 85.05, 86.88,
     * 88.71, 90.55, 92.30, 94.18, 95.96, 97.69, 99.38,101.16,102.87,
     *  0.00, 13.81, 19.11, 23.17, 26.66, 29.81, 32.78, 35.62, 38.33,
     * 40.92, 43.43, 45.91, 48.33, 50.64, 52.91, 55.22, 57.37, 59.61,
     * 61.71, 63.86, 65.93, 67.97, 70.07, 72.02, 73.98, 75.93, 77.86,
     * 79.77, 81.72, 83.57, 85.48, 87.27, 89.19, 90.91, 92.85, 94.58,
     * 96.39, 98.07, 99.90,101.66,103.42,  0.00, 13.59, 18.81, 22.87,
     * 26.43, 29.73, 32.79, 35.68, 38.47, 41.08, 43.73, 46.21, 48.65,
     * 51.03, 53.33, 55.59, 57.89, 60.00, 62.23, 64.35, 66.43, 68.47,
     * 70.56, 72.55, 74.54, 76.53, 78.52, 80.35, 82.34, 84.19, 86.00,
     * 87.91, 89.81, 91.59, 93.41, 95.21, 97.07, 98.86,100.56,102.37,
     *104.00,  0.00, 12.85, 18.26, 22.64, 26.46, 29.92, 33.12, 36.10,
     * 38.97, 41.69, 44.32, 46.88, 49.33, 51.76, 54.08, 56.40, 58.68,
     * 60.84, 63.03, 65.18, 67.26, 69.35, 71.42, 73.41, 75.39, 77.42,
     * 79.33, 81.31, 83.26, 85.09, 86.97, 88.89, 90.68, 92.52, 94.38/
      data (en(i),i=   856,  1026)/
     * 96.13, 98.00, 99.71,101.46,103.28,104.95,  0.00, 13.22, 18.77,
     * 23.14, 26.94, 30.35, 33.55, 36.58, 39.45, 42.19, 44.82, 47.38,
     * 49.90, 52.29, 54.67, 56.95, 59.18, 61.41, 63.63, 65.77, 67.87,
     * 69.93, 72.02, 74.04, 76.03, 78.03, 79.97, 81.89, 83.78, 85.71,
     * 87.56, 89.43, 91.31, 93.11, 94.91, 96.78, 98.48,100.36,102.07,
     *103.80,105.60,  0.00, 13.97, 19.47, 23.74, 27.43, 30.77, 33.91,
     * 36.83, 39.63, 42.35, 44.95, 47.48, 49.88, 52.25, 54.57, 56.86,
     * 59.14, 61.34, 63.45, 65.56, 67.66, 69.78, 71.83, 73.76, 75.71,
     * 77.77, 79.72, 81.60, 83.50, 85.38, 87.25, 89.19, 91.00, 92.74,
     * 94.64, 96.37, 98.14,100.08,101.72,103.48,105.15,  0.00, 14.42,
     * 20.05, 24.35, 28.00, 31.32, 34.40, 37.31, 40.08, 42.73, 45.23,
     * 47.74, 50.13, 52.51, 54.79, 57.08, 59.29, 61.43, 63.62, 65.73,
     * 67.75, 69.79, 71.83, 73.81, 75.79, 77.75, 79.66, 81.63, 83.44,
     * 85.32, 87.18, 89.12, 90.99, 92.77, 94.51, 96.31, 98.14, 99.81,
     *101.65,103.32,105.07,  0.00, 14.74, 20.47, 24.80, 28.49, 31.79,
     * 34.88, 37.76, 40.45, 43.09, 45.62, 48.08, 50.46, 52.83, 55.07,
     * 57.29, 59.50, 61.65, 63.76, 65.89, 67.92, 69.94, 71.98, 73.94,
     * 75.88, 77.81, 79.87, 81.72, 83.62, 85.47, 87.34, 89.18, 90.95,
     * 92.82, 94.56, 96.35, 98.21, 99.93,101.72,103.43,105.08,  0.00/
      data (en(i),i=  1027,  1197)/
     * 14.95, 20.86, 25.24, 28.96, 32.26, 35.29, 38.15, 40.81, 43.45,
     * 45.98, 48.37, 50.78, 53.06, 55.36, 57.55, 59.75, 61.92, 64.00,
     * 66.13, 68.10, 70.14, 72.22, 74.15, 76.16, 78.09, 80.00, 81.83,
     * 83.69, 85.63, 87.45, 89.34, 91.08, 92.95, 94.71, 96.50, 98.30,
     *100.05,101.77,103.49,105.24,  0.00, 15.13, 21.02, 25.44, 29.20,
     * 32.54, 35.55, 38.43, 41.12, 43.71, 46.23, 48.64, 51.01, 53.29,
     * 55.55, 57.76, 59.94, 62.08, 64.20, 66.31, 68.30, 70.34, 72.30,
     * 74.39, 76.34, 78.25, 80.20, 82.01, 83.90, 85.79, 87.63, 89.47,
     * 91.33, 93.16, 94.85, 96.77, 98.43,100.30,101.98,103.67,105.40,
     *  0.00, 15.31, 21.19, 25.64, 29.44, 32.73, 35.79, 38.65, 41.33,
     * 43.92, 46.42, 48.85, 51.19, 53.51, 55.78, 57.98, 60.15, 62.31,
     * 64.40, 66.50, 68.55, 70.58, 72.58, 74.51, 76.48, 78.47, 80.38,
     * 82.26, 84.16, 86.06, 87.85, 89.64, 91.53, 93.30, 95.15, 96.88,
     * 98.65,100.44,102.21,103.94,105.65,  0.00, 15.35, 21.30, 25.75,
     * 29.54, 32.88, 35.94, 38.75, 41.49, 44.10, 46.56, 49.00, 51.40,
     * 53.70, 55.99, 58.21, 60.39, 62.49, 64.62, 66.74, 68.72, 70.83,
     * 72.78, 74.79, 76.70, 78.69, 80.59, 82.52, 84.35, 86.30, 88.06,
     * 89.98, 91.72, 93.51, 95.39, 97.14, 98.90,100.65,102.41,104.20,
     *105.83,  0.00, 15.38, 21.28, 25.80, 29.56, 32.91, 35.98, 38.83/
      data (en(i),i=  1198,  1368)/
     * 41.53, 44.17, 46.69, 49.12, 51.53, 53.87, 56.12, 58.33, 60.54,
     * 62.71, 64.85, 66.90, 69.00, 71.00, 72.99, 75.06, 77.03, 78.88,
     * 80.85, 82.80, 84.71, 86.58, 88.41, 90.27, 92.13, 93.91, 95.68,
     * 97.42, 99.28,101.19,102.73,104.47,106.07,  0.00, 15.34, 21.24,
     * 25.70, 29.49, 32.84, 35.89, 38.81, 41.53, 44.17, 46.71, 49.23,
     * 51.63, 53.93, 56.25, 58.49, 60.72, 62.88, 65.07, 67.16, 69.19,
     * 71.30, 73.34, 75.28, 77.27, 79.28, 81.64, 83.09, 85.07, 86.93,
     * 88.81, 90.69, 92.53, 94.22, 96.05, 97.90, 99.65,101.41,103.22,
     *104.93,106.70,  0.00, 15.25, 21.07, 25.51, 29.30, 32.64, 35.76,
     * 38.70, 41.48, 44.17, 46.79, 49.26, 51.75, 54.11, 56.46, 58.72,
     * 61.03, 63.16, 65.34, 67.44, 69.55, 71.68, 73.70, 75.74, 77.66,
     * 79.62, 81.59, 83.61, 85.42, 87.29, 89.24, 91.05, 92.95, 94.80,
     * 96.56, 98.40,100.08,101.95,103.64,105.44,107.21,  0.00, 14.99,
     * 20.73, 25.13, 28.92, 32.37, 35.56, 38.61, 41.46, 44.25, 46.90,
     * 49.44, 52.00, 54.44, 56.79, 59.14, 61.34, 63.57, 65.76, 67.94,
     * 70.05, 72.14, 74.23, 76.27, 78.19, 80.20, 82.22, 84.20, 86.11,
     * 87.98, 89.76, 91.64, 93.56, 95.30, 97.20, 99.01,100.80,102.51,
     *104.32,105.99,107.80,  0.00, 14.19, 20.02, 24.67, 28.70, 32.35,
     * 35.70, 38.87, 41.85, 44.66, 47.40, 50.01, 52.53, 55.05, 57.43/
      data (en(i),i=  1369,  1539)/
     * 59.75, 62.03, 64.28, 66.48, 68.68, 70.78, 72.96, 74.98, 77.01,
     * 79.05, 81.02, 83.00, 84.97, 86.92, 88.80, 90.71, 92.59, 94.39,
     * 96.27, 98.08, 99.88,101.75,103.52,105.29,106.98,108.72,  0.00,
     * 14.19, 20.05, 24.68, 28.68, 32.24, 35.54, 38.64, 41.60, 44.42,
     * 47.16, 49.77, 52.31, 54.78, 57.20, 59.54, 61.88, 64.10, 66.36,
     * 68.52, 70.65, 72.77, 74.79, 76.87, 78.91, 80.89, 82.85, 84.80,
     * 86.73, 88.68, 90.60, 92.43, 94.28, 96.12, 97.94, 99.78,101.59,
     *103.35,105.17,106.88,108.67,  0.00, 14.94, 20.82, 25.37, 29.23,
     * 32.71, 35.97, 38.96, 41.89, 44.66, 47.38, 49.91, 52.46, 54.87,
     * 57.19, 59.61, 61.87, 64.03, 66.26, 68.40, 70.52, 72.65, 74.66,
     * 76.72, 78.77, 80.70, 82.66, 84.59, 86.57, 88.40, 90.21, 92.09,
     * 94.00, 95.81, 97.62, 99.47,101.26,102.99,104.83,106.50,108.34,
     *  0.00, 15.41, 21.41, 26.01, 29.88, 33.16, 36.55, 39.57, 42.42,
     * 45.17, 47.75, 50.31, 52.80, 55.19, 57.52, 59.80, 62.06, 64.30,
     * 66.42, 68.57, 70.69, 72.79, 74.80, 76.85, 78.84, 80.79, 82.76,
     * 84.64, 86.57, 88.45, 90.14, 92.16, 94.02, 95.84, 97.68, 99.45,
     *101.20,103.00,104.79,106.48,108.24,  0.00, 15.73, 21.89, 26.55,
     * 30.44, 33.92, 37.10, 40.09, 42.95, 45.65, 48.22, 50.72, 53.19,
     * 55.56, 57.87, 60.12, 62.40, 64.55, 66.69, 68.82, 70.95, 73.01/
      data (en(i),i=  1540,  1710)/
     * 75.02, 76.98, 78.98, 80.98, 82.81, 84.83, 86.67, 88.55, 90.40,
     * 92.27, 94.13, 95.98, 97.74, 99.54,101.08,103.02,104.80,106.54,
     *108.21,  0.00, 16.00, 22.20, 27.11, 30.86, 34.36, 37.56, 40.55,
     * 43.33, 45.88, 48.66, 51.14, 53.53, 55.92, 58.26, 60.47, 62.66,
     * 64.84, 67.03, 69.05, 71.17, 73.21, 75.21, 77.20, 79.17, 81.19,
     * 83.09, 84.98, 86.87, 88.68, 90.61, 92.42, 94.22, 96.03, 97.85,
     * 99.64,101.38,103.12,104.91,106.70,108.36,  0.00, 16.19, 22.51,
     * 27.26, 31.22, 34.77, 37.95, 40.95, 43.73, 46.42, 48.96, 51.49,
     * 53.95, 56.30, 58.55, 60.79, 62.96, 65.08, 67.27, 69.37, 71.45,
     * 73.46, 75.49, 77.45, 79.50, 81.36, 83.26, 85.24, 87.06, 88.86,
     * 90.80, 92.58, 94.40, 96.14, 98.03, 99.81,101.58,103.29,104.99,
     *106.73,108.49,  0.00, 16.35, 22.69, 27.53, 31.50, 35.04, 38.24,
     * 41.23, 44.06, 46.73, 49.31, 51.79, 54.18, 56.53, 58.80, 61.06,
     * 63.25, 65.43, 67.53, 69.64, 71.65, 73.69, 75.65, 77.64, 79.62,
     * 81.50, 83.45, 85.36, 87.25, 89.08, 90.97, 92.80, 94.62, 96.48,
     * 98.23, 99.97,101.68,103.53,105.29,106.94,108.67,  0.00, 16.44,
     * 22.88, 27.71, 31.71, 35.28, 38.50, 41.51, 44.31, 47.00, 49.54,
     * 52.05, 54.42, 56.76, 59.06, 61.31, 63.48, 65.60, 67.72, 69.84,
     * 71.89, 73.86, 75.90, 77.88, 79.84, 81.74, 83.67, 85.51, 87.50/
      data (en(i),i=  1711,  1881)/
     * 89.38, 91.13, 92.97, 94.83, 96.63, 98.42,100.23,101.88,103.66,
     *105.44,107.18,108.93,  0.00, 16.58, 23.01, 27.83, 31.85, 35.45,
     * 38.67, 41.66, 44.48, 47.17, 49.75, 52.24, 54.61, 56.97, 59.23,
     * 61.46, 63.67, 65.86, 67.91, 70.02, 72.09, 74.12, 76.09, 78.08,
     * 80.02, 82.01, 83.89, 85.76, 87.66, 89.52, 91.38, 93.19, 95.09,
     * 96.83, 98.68,100.41,102.15,103.93,105.69,107.40,109.12,  0.00,
     * 16.62, 23.04, 27.91, 31.97, 35.54, 38.80, 41.80, 44.63, 47.30,
     * 49.85, 52.34, 54.75, 57.09, 59.41, 61.65, 63.83, 65.98, 68.10,
     * 70.21, 72.29, 74.31, 76.31, 78.24, 80.26, 82.21, 84.06, 85.97,
     * 87.84, 89.71, 91.66, 93.44, 95.22, 97.04, 98.83,100.72,102.37,
     *104.19,105.89,107.62,109.31,  0.00, 16.66, 23.05, 27.93, 32.01,
     * 35.60, 38.84, 41.81, 44.67, 47.33, 49.91, 52.40, 54.85, 57.20,
     * 59.48, 61.73, 63.95, 66.11, 68.25, 70.31, 72.47, 74.47, 76.53,
     * 78.53, 80.49, 82.44, 84.34, 86.17, 88.18, 89.97, 91.86, 93.68,
     * 95.56, 97.29, 99.09,100.91,102.70,104.48,106.27,107.91,109.59,
     *  0.00, 16.62, 23.02, 27.89, 31.93, 35.51, 38.78, 41.80, 44.61,
     * 47.36, 49.90, 52.45, 54.87, 57.24, 59.52, 61.84, 64.05, 66.22,
     * 68.45, 70.51, 72.62, 74.60, 76.67, 78.72, 80.71, 82.61, 84.63,
     * 86.52, 88.27, 90.30, 92.18, 94.01, 95.78, 97.61, 99.46,101.16/
      data (en(i),i=  1882,  2052)/
     *102.93,104.81,106.46,108.21,110.01,  0.00, 16.57, 22.92, 27.73,
     * 31.82, 35.37, 38.61, 41.64, 44.55, 47.25, 49.88, 52.41, 54.85,
     * 57.26, 59.66, 61.95, 64.21, 66.39, 68.54, 70.67, 72.80, 74.86,
     * 76.93, 78.91, 80.95, 82.93, 84.82, 86.85, 88.68, 90.64, 92.46,
     * 94.31, 96.22, 98.03, 99.81,101.60,103.28,105.22,106.93,108.68,
     *110.46,  0.00, 16.47, 22.71, 27.51, 31.54, 35.11, 38.38, 41.43,
     * 44.37, 47.09, 49.82, 52.40, 54.97, 57.35, 59.76, 62.09, 64.36,
     * 66.60, 68.77, 70.98, 73.10, 75.19, 77.21, 79.29, 81.34, 83.21,
     * 85.21, 87.15, 89.04, 90.97, 92.89, 94.73, 96.56, 98.45,100.33,
     *102.04,103.90,105.61,107.36,109.22,110.86,  0.00, 16.19, 22.33,
     * 27.09, 31.11, 34.68, 38.05, 41.22, 44.65, 47.06, 49.82, 52.47,
     * 55.04, 57.50, 59.97, 62.35, 64.64, 66.92, 69.18, 71.36, 73.52,
     * 75.56, 77.66, 79.76, 81.78, 83.77, 85.74, 87.75, 89.66, 91.52,
     * 93.45, 95.36, 97.23, 99.09,100.81,102.71,104.51,106.33,108.03,
     *109.83,111.59,  0.00, 15.39, 21.54, 26.45, 30.70, 34.51, 38.01,
     * 41.31, 44.39, 47.36, 50.17, 52.85, 55.51, 58.02, 60.50, 62.90,
     * 65.27, 67.56, 69.77, 72.01, 74.15, 76.32, 78.46, 80.47, 82.52,
     * 84.57, 86.58, 88.55, 90.47, 92.40, 94.26, 96.17, 98.06, 99.87,
     *101.78,103.53,105.32,107.15,108.94,110.67,112.51,  0.00, 15.07/
      data (en(i),i=  2053,  2223)/
     * 21.20, 26.05, 30.23, 33.94, 37.37, 40.57, 43.57, 46.49, 49.28,
     * 51.94, 54.57, 57.09, 59.57, 61.98, 64.30, 66.62, 68.88, 71.08,
     * 73.28, 75.38, 77.50, 79.59, 81.66, 83.65, 85.66, 87.58, 89.58,
     * 91.53, 93.41, 95.32, 97.16, 99.04,100.91,102.70,104.50,106.29,
     *108.06,109.87,111.56,  0.00, 15.86, 22.00, 26.78, 30.90, 34.36,
     * 37.96, 40.96, 44.02, 46.79, 49.54, 52.20, 54.71, 57.27, 59.70,
     * 62.03, 64.39, 66.64, 68.85, 71.07, 73.20, 75.38, 77.46, 79.52,
     * 81.48, 83.50, 85.58, 87.46, 89.35, 91.32, 93.24, 95.06, 96.98,
     * 98.84,100.60,102.53,104.25,106.13,107.84,109.59,111.31,  0.00,
     * 16.29, 22.60, 27.45, 31.53, 35.18, 38.48, 41.56, 44.51, 47.35,
     * 50.04, 52.70, 55.20, 57.66, 60.03, 62.41, 64.65, 66.92, 69.15,
     * 71.32, 73.46, 75.49, 77.52, 79.62, 81.63, 83.62, 85.64, 87.62,
     * 89.48, 91.36, 93.31, 95.10, 96.96, 98.77,100.67,102.46,104.20,
     *106.03,107.73,109.52,111.27,  0.00, 16.63, 23.07, 27.97, 32.11,
     * 35.77, 39.05, 42.15, 45.08, 47.88, 50.58, 53.16, 55.67, 58.10,
     * 60.40, 62.78, 65.05, 67.26, 69.45, 71.57, 73.69, 75.82, 77.90,
     * 79.82, 81.90, 83.87, 85.78, 87.72, 89.61, 91.50, 93.44, 95.22,
     * 97.10, 98.91,100.77,102.53,104.35,106.07,107.88,109.60,111.40,
     *  0.00, 16.93, 23.44, 28.45, 32.59, 36.27, 39.55, 42.67, 45.62/
      data (en(i),i=  2224,  2394)/
     * 48.37, 51.06, 53.63, 56.19, 58.55, 60.87, 63.21, 65.44, 67.62,
     * 69.77, 71.92, 74.05, 76.14, 78.11, 80.15, 82.10, 84.08, 86.08,
     * 88.00, 89.85, 91.82, 93.57, 95.40, 97.29, 99.11,100.91,102.71,
     *104.49,106.21,108.03,109.75,111.47,  0.00, 17.11, 23.75, 28.78,
     * 33.01, 36.67, 40.03, 43.11, 46.05, 48.83, 51.50, 54.06, 56.56,
     * 58.96, 61.30, 63.53, 65.80, 68.02, 70.11, 72.30, 74.34, 76.40,
     * 78.40, 80.47, 82.40, 84.33, 86.32, 88.22, 90.17, 91.95, 93.90,
     * 95.76, 97.47, 99.28,101.12,102.89,104.63,106.42,108.19,110.00,
     *111.61,  0.00, 17.27, 24.00, 29.05, 33.29, 37.04, 40.38, 43.51,
     * 46.43, 49.23, 51.88, 54.46, 56.87, 59.31, 61.62, 63.93, 66.09,
     * 68.28, 70.43, 72.59, 74.61, 76.65, 78.77, 80.69, 82.72, 84.60,
     * 86.55, 88.43, 90.36, 92.20, 94.07, 95.88, 97.74, 99.46,101.33,
     *103.14,104.90,106.69,108.36,110.13,111.87,  0.00, 17.73, 24.36,
     * 29.30, 33.56, 37.32, 40.68, 43.94, 46.76, 49.53, 52.19, 54.75,
     * 57.22, 59.68, 62.02, 64.26, 66.45, 68.61, 70.78, 72.94, 74.88,
     * 76.93, 78.95, 81.00, 82.91, 84.87, 86.88, 88.78, 90.63, 92.48,
     * 94.29, 96.20, 98.05, 99.81,101.59,103.31,105.13,106.95,108.68,
     *110.24,112.00,  0.00, 17.56, 24.32, 29.48, 33.81, 37.53, 40.96,
     * 44.08, 47.02, 49.83, 52.47, 55.02, 57.48, 59.85, 62.20, 64.46/
      data (en(i),i=  2395,  2565)/
     * 66.70, 68.88, 70.98, 73.08, 75.18, 77.20, 79.17, 81.17, 83.20,
     * 85.09, 86.98, 88.93, 90.81, 92.63, 94.54, 96.32, 98.13, 99.97,
     *101.72,103.55,105.24,107.01,108.72,110.55,112.24,  0.00, 17.67,
     * 24.45, 29.65, 33.95, 37.72, 41.15, 44.28, 47.25, 50.00, 52.69,
     * 55.22, 57.67, 60.01, 62.38, 64.66, 66.90, 69.06, 71.18, 73.32,
     * 75.36, 77.38, 79.43, 81.37, 83.40, 85.34, 87.28, 89.08, 90.99,
     * 92.88, 94.72, 96.57, 98.39,100.19,101.91,103.72,105.44,107.22,
     *108.92,110.68,112.41,  0.00, 17.73, 24.53, 29.73, 34.09, 37.86,
     * 41.28, 44.46, 47.40, 50.17, 52.85, 55.39, 57.87, 60.27, 62.53,
     * 64.84, 67.05, 69.25, 71.40, 73.45, 75.56, 77.60, 79.58, 81.59,
     * 83.56, 85.48, 87.45, 89.33, 91.19, 93.04, 94.95, 96.79, 98.54,
     *100.42,102.17,103.96,105.65,107.48,109.28,110.94,112.70,  0.00,
     * 17.74, 24.58, 29.81, 34.11, 37.94, 41.36, 44.54, 47.49, 50.27,
     * 52.92, 55.49, 57.94, 60.34, 62.69, 64.96, 67.22, 69.42, 71.53,
     * 73.62, 75.74, 77.76, 79.77, 81.78, 83.73, 85.71, 87.60, 89.54,
     * 91.39, 93.23, 94.95, 96.92, 98.82,100.56,102.41,104.11,105.94,
     *107.63,109.43,111.15,112.85,  0.00, 17.79, 24.62, 29.80, 34.15,
     * 37.96, 41.40, 44.56, 47.53, 50.32, 53.00, 55.56, 58.04, 60.41,
     * 62.75, 65.04, 67.28, 69.50, 71.62, 73.73, 75.86, 77.88, 79.94/
      data (en(i),i=  2566,  2736)/
     * 81.89, 83.90, 85.84, 87.77, 89.69, 91.59, 93.44, 95.31, 97.16,
     * 98.97,100.84,102.61,104.35,106.17,107.95,109.68,111.34,113.15,
     *  0.00, 17.83, 24.58, 29.75, 34.11, 37.90, 41.38, 44.52, 47.49,
     * 50.30, 52.93, 55.55, 58.01, 60.48, 62.77, 65.05, 67.34, 69.54,
     * 71.70, 73.88, 75.92, 78.04, 80.08, 82.14, 84.03, 86.05, 87.99,
     * 89.90, 91.82, 93.65, 95.57, 97.43, 99.22,101.05,102.82,104.66,
     *106.41,108.17,109.99,111.67,113.42,  0.00, 17.74, 24.58, 29.67,
     * 34.00, 37.81, 41.25, 44.41, 47.36, 50.19, 52.90, 55.48, 57.97,
     * 60.40, 62.80, 65.11, 67.36, 69.60, 71.79, 73.91, 76.06, 78.14,
     * 80.16, 82.25, 84.28, 86.24, 88.11, 90.06, 92.01, 93.93, 95.81,
     * 97.68, 99.53,101.36,103.16,104.95,106.65,108.51,110.25,112.00,
     *113.67,  0.00, 17.64, 24.38, 29.51, 33.81, 37.59, 41.04, 44.21,
     * 47.18, 50.02, 52.75, 55.34, 57.88, 60.38, 62.80, 65.12, 67.47,
     * 69.64, 71.87, 74.06, 76.19, 78.30, 80.39, 82.41, 84.43, 86.43,
     * 88.38, 90.32, 92.27, 94.19, 96.07, 97.92, 99.82,101.66,103.52,
     *105.28,107.01,108.85,110.63,112.37,114.07,  0.00, 17.49, 24.15,
     * 29.25, 33.52, 37.30, 40.73, 43.92, 46.93, 49.81, 52.53, 55.25,
     * 57.80, 60.34, 62.77, 65.18, 67.48, 69.75, 71.97, 74.19, 76.42,
     * 78.53, 80.58, 82.70, 84.75, 86.69, 88.72, 90.68, 92.65, 94.53/
      data (en(i),i=  2737,  2907)/
     * 96.50, 98.31,100.19,102.08,103.91,105.73,107.55,109.28,111.06,
     *112.86,114.62,  0.00, 17.24, 23.78, 28.81, 33.02, 36.84, 40.28,
     * 43.56, 46.67, 49.60, 52.43, 55.18, 57.83, 60.35, 62.88, 65.32,
     * 67.70, 70.02, 72.30, 74.50, 76.68, 78.83, 80.95, 83.06, 85.13,
     * 87.14, 89.20, 91.22, 93.14, 95.01, 97.01, 98.84,100.69,102.56,
     *104.38,106.27,108.08,109.88,111.65,113.40,115.26,  0.00, 16.45,
     * 22.92, 28.05, 32.48, 36.44, 40.09, 43.50, 46.71, 49.76, 52.69,
     * 55.47, 58.17, 60.77, 63.35, 65.79, 68.18, 70.58, 72.84, 75.10,
     * 77.36, 79.53, 81.63, 83.79, 85.84, 87.86, 89.91, 91.90, 93.84,
     * 95.84, 97.78, 99.66,101.57,103.36,105.23,107.07,108.90,110.74,
     *112.48,114.23,116.01,  0.00, 15.89, 22.26, 27.29, 31.62, 35.46,
     * 39.03, 42.35, 45.46, 48.41, 51.22, 53.99, 56.68, 59.25, 61.75,
     * 64.22, 66.60, 68.92, 71.24, 73.46, 75.70, 77.84, 80.00, 82.08,
     * 84.17, 86.24, 88.26, 90.29, 92.19, 94.21, 96.10, 98.02, 99.90,
     *101.85,103.66,105.55,107.35,109.14,110.97,112.77,114.48,  0.00,
     * 16.48, 23.07, 28.06, 32.34, 36.13, 39.59, 42.81, 45.88, 48.77,
     * 52.68, 54.54, 56.93, 59.43, 61.94, 64.32, 66.69, 69.07, 71.27,
     * 73.48, 75.69, 77.84, 79.95, 82.01, 84.08, 86.18, 88.14, 90.18,
     * 92.13, 94.03, 95.93, 97.87, 99.79,101.64,103.40,105.31,107.06/
      data (en(i),i=  2908,  3078)/
     *108.93,110.63,112.44,114.11,  0.00, 17.09, 23.71, 28.79, 33.05,
     * 36.83, 40.26, 43.49, 46.47, 49.37, 52.16, 54.80, 57.46, 59.93,
     * 62.40, 64.74, 67.04, 69.41, 71.62, 73.82, 76.04, 78.14, 80.24,
     * 82.26, 84.29, 86.36, 88.32, 90.27, 92.29, 94.20, 96.11, 97.94,
     * 99.85,101.71,103.56,105.29,107.12,108.99,110.66,112.49,114.21,
     *  0.00, 17.44, 24.16, 29.32, 33.63, 37.47, 40.92, 44.07, 47.04,
     * 49.95, 52.72, 55.38, 57.91, 60.40, 62.87, 65.20, 67.58, 69.83,
     * 72.00, 74.16, 76.32, 78.52, 80.50, 82.58, 84.65, 86.61, 88.58,
     * 90.54, 92.47, 94.37, 96.25, 98.14, 99.99,101.83,103.68,105.52,
     *107.25,109.05,110.89,112.54,114.29,  0.00, 17.70, 24.51, 29.74,
     * 34.10, 37.99, 41.43, 44.63, 47.58, 50.49, 53.25, 55.87, 58.43,
     * 60.93, 63.32, 65.67, 67.94, 70.23, 72.42, 74.61, 76.69, 78.87,
     * 80.83, 82.90, 84.95, 86.93, 88.88, 90.81, 92.72, 94.66, 96.49,
     * 98.44,100.28,102.08,103.91,105.74,107.49,109.21,111.03,112.72,
     *114.53,  0.00, 17.94, 24.83, 30.14, 34.56, 38.43, 41.90, 45.11,
     * 48.12, 51.00, 53.69, 56.36, 58.91, 61.36, 63.74, 66.12, 68.42,
     * 70.65, 72.84, 74.98, 77.11, 79.18, 81.25, 83.24, 85.24, 87.20,
     * 89.19, 91.09, 93.09, 94.97, 96.82, 98.68,100.49,102.36,104.17,
     *105.90,107.69,109.47,111.25,112.99,114.67,  0.00, 18.11, 25.10/
      data (en(i),i=  3079,  3249)/
     * 30.44, 34.92, 38.78, 42.29, 45.52, 48.55, 51.46, 54.17, 56.82,
     * 59.39, 61.76, 64.17, 66.51, 68.81, 71.02, 73.17, 75.31, 77.46,
     * 79.52, 81.56, 83.56, 85.63, 87.56, 89.46, 91.43, 93.28, 95.22,
     * 97.10, 98.97,100.80,102.58,104.34,106.13,107.91,109.71,111.45,
     *113.19,114.87,  0.00, 18.26, 25.30, 30.69, 35.20, 39.12, 42.62,
     * 45.87, 48.92, 51.80, 54.56, 57.18, 59.72, 62.17, 64.55, 66.88,
     * 69.14, 71.39, 73.53, 75.67, 77.81, 79.91, 81.91, 83.91, 85.91,
     * 87.83, 89.79, 91.66, 93.61, 95.47, 97.35, 99.16,100.99,102.77,
     *104.66,106.44,108.18,109.93,111.68,113.36,115.11,  0.00, 18.43,
     * 25.49, 30.92, 35.44, 39.41, 42.94, 46.21, 49.27, 52.14, 54.90,
     * 57.52, 60.06, 62.53, 64.87, 67.20, 69.49, 71.69, 73.87, 76.03,
     * 78.08, 80.11, 82.19, 84.20, 86.18, 88.09, 90.07, 91.96, 93.81,
     * 95.69, 97.57, 99.40,101.23,103.04,104.88,106.58,108.39,110.13,
     *111.90,113.69,115.34,  0.00, 18.53, 25.65, 31.08, 35.65, 39.64,
     * 43.17, 46.44, 49.54, 52.44, 55.18, 57.78, 60.33, 62.80, 65.15,
     * 67.48, 69.77, 71.96, 74.16, 76.23, 78.37, 80.39, 82.46, 84.47,
     * 86.44, 88.41, 90.27, 92.17, 94.13, 95.97, 97.82, 99.65,101.41,
     *103.34,105.14,106.89,108.59,110.43,112.12,113.82,115.57,  0.00,
     * 18.62, 25.75, 31.25, 35.81, 39.79, 43.37, 46.66, 49.75, 52.69/
      data (en(i),i=  3250,  3420)/
     * 55.43, 58.04, 60.57, 63.06, 65.42, 67.73, 70.00, 72.24, 74.40,
     * 76.48, 78.55, 80.63, 82.68, 84.65, 86.58, 88.58, 90.52, 92.43,
     * 94.34, 96.21, 98.05, 99.83,101.71,103.54,105.29,107.10,108.86,
     *110.62,112.31,114.01,115.71,  0.00, 18.71, 25.88, 31.36, 35.95,
     * 39.92, 43.52, 46.85, 49.92, 52.82, 55.59, 58.23, 60.77, 63.26,
     * 65.62, 67.92, 70.14, 72.38, 74.54, 76.69, 78.77, 80.84, 82.84,
     * 84.85, 86.84, 88.79, 90.70, 92.55, 94.59, 96.37, 98.26,100.00,
     *101.86,103.66,105.53,107.31,109.00,110.74,112.54,114.21,115.96,
     *  0.00, 18.73, 25.96, 31.45, 36.03, 40.01, 43.62, 46.96, 50.04,
     * 52.96, 55.75, 58.35, 60.90, 63.36, 65.74, 68.07, 70.31, 72.57,
     * 74.73, 76.81, 78.95, 80.98, 83.04, 85.04, 87.00, 88.99, 90.96,
     * 92.75, 94.70, 96.57, 98.42,100.23,102.08,103.86,105.67,107.40,
     *109.16,110.89,112.73,114.44,116.17,  0.00, 18.81, 25.97, 31.47,
     * 36.05, 40.07, 43.70, 47.04, 50.15, 53.04, 55.82, 58.45, 60.98,
     * 63.46, 65.86, 68.15, 70.46, 72.67, 74.90, 76.93, 79.04, 81.17,
     * 83.13, 85.15, 87.20, 89.12, 91.07, 92.93, 94.83, 96.72, 98.63,
     *100.41,102.25,104.09,105.91,107.66,109.40,111.19,112.93,114.62,
     *116.33,  0.00, 18.82, 25.99, 31.73, 36.07, 40.08, 43.73, 47.08,
     * 50.14, 53.07, 55.86, 58.48, 61.04, 63.49, 65.89, 68.21, 70.49/
      data (en(i),i=  3421,  3591)/
     * 72.74, 74.91, 77.07, 79.16, 81.20, 83.24, 85.31, 87.26, 89.22,
     * 91.15, 93.10, 95.00, 96.86, 98.78,100.56,102.46,104.25,106.02,
     *107.83,109.65,111.39,113.14,114.82,116.59,  0.00, 18.80, 25.97,
     * 31.48, 36.06, 40.08, 43.70, 47.02, 50.13, 53.04, 55.81, 58.49,
     * 61.02, 63.51, 65.90, 68.26, 70.52, 72.75, 74.95, 77.11, 79.23,
     * 81.33, 83.34, 85.41, 87.41, 89.39, 91.35, 93.24, 95.18, 97.03,
     * 98.94,100.75,102.69,104.41,106.27,108.01,109.84,111.57,113.39,
     *115.08,116.88,  0.00, 18.78, 25.91, 31.38, 35.95, 40.00, 43.63,
     * 46.95, 50.04, 52.93, 55.73, 58.39, 60.96, 63.46, 65.86, 68.23,
     * 70.51, 72.80, 74.95, 77.11, 79.27, 81.36, 83.48, 85.48, 87.45,
     * 89.47, 91.43, 93.44, 95.34, 97.18, 99.06,101.04,102.82,104.60,
     *106.48,108.23,110.09,111.84,113.56,115.38,117.12,  0.00, 18.74,
     * 25.81, 31.26, 35.83, 39.85, 43.46, 46.75, 49.86, 52.79, 55.59,
     * 58.25, 60.83, 63.36, 65.77, 68.16, 70.50, 72.74, 74.97, 77.20,
     * 79.33, 81.45, 83.55, 85.56, 87.65, 89.57, 91.58, 93.52, 95.53,
     * 97.39, 99.32,101.15,103.02,104.91,106.72,108.55,110.33,112.08,
     *113.79,115.58,117.41,  0.00, 18.65, 25.68, 31.04, 35.61, 39.62,
     * 43.21, 46.53, 49.62, 52.54, 55.39, 58.03, 60.69, 63.20, 65.68,
     * 68.09, 70.43, 72.74, 75.01, 77.21, 79.37, 81.51, 83.59, 85.64/
      data (en(i),i=  3592,  3762)/
     * 87.79, 89.80, 91.78, 93.74, 95.70, 97.65, 99.53,101.43,103.29,
     *105.19,106.99,108.76,110.64,112.51,114.11,115.98,117.80,  0.00,
     * 18.48, 25.44, 30.80, 35.28, 39.25, 42.84, 46.16, 49.28, 52.27,
     * 55.10, 57.87, 60.51, 63.09, 65.61, 68.04, 70.44, 72.74, 75.06,
     * 77.32, 79.50, 81.65, 83.78, 85.89, 87.91, 90.00, 92.05, 94.04,
     * 96.00, 97.90, 99.87,101.76,103.60,105.54,107.37,109.16,110.97,
     *112.79,114.66,116.46,118.16,  0.00, 18.20, 25.05, 30.34, 34.77,
     * 38.72, 42.33, 45.72, 48.88, 51.95, 54.89, 57.70, 60.42, 63.03,
     * 65.57, 68.11, 70.53, 72.93, 75.20, 77.47, 79.70, 81.94, 84.07,
     * 86.22, 88.30, 90.36, 92.44, 94.44, 96.37, 98.40,100.28,102.19,
     *104.11,105.94,107.87,109.73,111.57,113.32,115.19,116.96,118.66,
     *  0.00, 17.44, 24.17, 29.50, 34.12, 38.24, 42.06, 45.55, 48.87,
     * 52.00, 54.99, 57.89, 60.69, 63.34, 65.96, 68.52, 70.96, 73.36,
     * 75.72, 78.02, 80.31, 82.54, 84.68, 86.89, 88.93, 91.00, 93.10,
     * 95.09, 97.11, 99.11,101.02,102.95,104.90,106.81,108.66,110.49,
     *112.39,114.19,115.92,117.80,119.61,  0.00, 16.63, 23.23, 28.45,
     * 32.93, 36.93, 40.58, 43.98, 47.19, 50.21, 53.09, 55.88, 58.58,
     * 61.25, 63.82, 66.28, 68.75, 71.14, 73.46, 75.76, 78.01, 80.18,
     * 82.37, 84.54, 86.63, 88.68, 90.73, 92.80, 94.73, 96.74, 98.72/
      data (en(i),i=  3763,  3933)/
     *100.69,102.52,104.42,106.34,108.16,110.05,111.83,113.72,115.45,
     *117.24,  0.00, 17.35, 24.11, 29.26, 33.64, 37.61, 41.19, 44.52,
     * 47.62, 50.62, 53.44, 56.23, 58.95, 61.52, 63.36, 66.50, 68.93,
     * 71.21, 73.55, 75.86, 78.04, 80.19, 82.42, 84.52, 86.57, 88.67,
     * 90.64, 92.75, 94.67, 96.34, 98.58,100.53,102.43,104.25,106.19,
     *108.00,109.87,111.70,113.42,115.31,116.99,  0.00, 17.95, 24.68,
     * 29.93, 34.41, 38.33, 41.91, 45.21, 48.28, 51.26, 54.03, 56.80,
     * 59.43, 62.08, 64.54, 66.99, 69.35, 71.69, 73.98, 76.20, 78.40,
     * 80.60, 82.67, 84.78, 86.86, 88.94, 90.97, 92.93, 94.83, 96.81,
     * 98.72,100.63,102.58,104.47,106.24,108.15,109.90,111.71,113.49,
     *115.31,117.05,  0.00, 18.31, 25.16, 30.50, 35.06, 38.98, 42.57,
     * 45.85, 49.03, 51.90, 54.74, 57.43, 60.06, 62.66, 65.11, 67.46,
     * 69.86, 72.16, 74.45, 76.65, 78.84, 80.98, 83.11, 85.26, 87.22,
     * 89.33, 91.29, 93.24, 95.18, 97.06, 99.03,100.95,102.79,104.60,
     *106.44,108.29,110.06,111.95,113.69,115.49,117.34,  0.00, 18.53,
     * 25.54, 30.95, 35.52, 39.52, 43.11, 46.45, 49.52, 52.45, 55.27,
     * 57.99, 60.62, 63.12, 65.62, 68.04, 70.36, 72.61, 74.93, 77.12,
     * 79.25, 81.37, 83.43, 85.59, 87.57, 89.60, 91.58, 93.55, 95.46,
     * 97.39, 99.34,101.22,103.05,104.92,106.70,108.51,110.34,112.12/
      data (en(i),i=  3934,  4104)/
     *113.98,115.66,117.47,  0.00, 18.68, 26.29, 31.36, 36.24, 39.97,
     * 43.62, 46.95, 50.27, 52.98, 55.79, 58.53, 61.15, 63.65, 66.12,
     * 68.65, 70.84, 73.04, 75.30, 77.48, 79.64, 81.77, 83.90, 86.05,
     * 87.98, 90.13, 91.89, 93.90, 95.78, 97.76, 99.58,101.49,103.52,
     *105.22,107.03,108.84,110.62,112.43,114.08,115.92,117.67,  0.00,
     * 18.91, 26.20, 31.75, 36.34, 40.39, 44.03, 47.42, 50.54, 53.48,
     * 56.30, 59.02, 61.68, 64.09, 66.56, 68.95, 71.29, 73.51, 75.73,
     * 77.99, 80.15, 82.25, 84.24, 86.31, 88.30, 90.36, 92.27, 94.22,
     * 96.21, 98.08, 99.94,101.81,103.67,105.51,107.34,109.13,110.83,
     *112.65,114.45,116.12,117.99,  0.00, 19.00, 26.32, 31.92, 36.64,
     * 40.72, 44.38, 47.73, 50.91, 53.85, 56.72, 59.44, 62.04, 64.54,
     * 66.93, 69.34, 71.69, 73.93, 76.17, 78.32, 80.49, 82.56, 84.61,
     * 86.65, 88.66, 90.64, 92.63, 94.57, 96.48, 98.43,100.27,102.15,
     *103.92,105.76,107.52,109.36,111.11,112.87,114.74,116.37,118.11,
     *  0.00, 19.15, 26.53, 32.20, 36.87, 41.03, 44.72, 48.13, 51.24,
     * 54.26, 57.05, 59.80, 62.42, 64.88, 67.39, 69.76, 72.06, 74.30,
     * 76.52, 78.71, 80.79, 82.90, 84.94, 87.00, 89.02, 91.00, 92.95,
     * 94.88, 96.74, 98.66,100.55,102.32,104.20,106.02,107.90,109.68,
     *111.38,113.16,114.97,116.70,118.36,  0.00, 19.29, 26.69, 32.39/
      data (en(i),i=  4105,  4275)/
     * 37.12, 41.24, 45.03, 48.40, 51.56, 54.56, 57.40, 60.13, 62.77,
     * 65.27, 67.71, 70.08, 72.40, 74.62, 76.82, 79.02, 81.09, 83.21,
     * 85.27, 87.28, 89.36, 91.30, 93.23, 95.15, 97.08, 98.96,100.80,
     *102.63,104.54,106.27,108.15,109.86,111.65,113.45,115.25,116.96,
     *118.68,  0.00, 19.45, 26.82, 32.54, 37.31, 41.47, 45.25, 48.65,
     * 51.77, 54.84, 57.66, 60.43, 63.05, 65.59, 67.98, 70.39, 72.70,
     * 74.92, 77.12, 79.29, 81.40, 83.55, 85.60, 87.63, 89.54, 91.62,
     * 93.49, 95.41, 97.42, 99.19,101.05,102.99,104.73,106.54,108.33,
     *110.18,111.87,113.65,115.43,117.17,118.86,  0.00, 19.52, 26.98,
     * 32.72, 37.50, 41.70, 45.43, 48.86, 52.05, 55.11, 57.96, 60.67,
     * 63.28, 65.83, 68.27, 70.62, 72.96, 75.23, 77.40, 79.55, 81.68,
     * 83.80, 85.80, 87.80, 89.81, 91.80, 93.81, 95.68, 97.51, 99.40,
     *101.28,103.18,105.00,106.77,108.62,110.31,112.11,113.86,115.61,
     *117.36,119.15,  0.00, 19.63, 27.09, 32.82, 37.65, 41.82, 45.59,
     * 49.01, 52.22, 55.30, 58.13, 60.86, 63.50, 66.01, 68.48, 70.82,
     * 73.16, 75.43, 77.57, 79.76, 81.91, 83.95, 86.04, 88.07, 90.06,
     * 92.03, 93.96, 95.86, 97.78, 99.61,101.50,103.42,105.20,106.98,
     *108.85,110.53,112.32,114.15,115.77,117.61,119.25,  0.00, 19.65,
     * 27.17, 32.89, 37.74, 41.95, 45.73, 49.16, 52.38, 55.45, 58.33/
      data (en(i),i=  4276,  4446)/
     * 61.03, 63.67, 66.20, 68.65, 71.01, 73.33, 75.56, 77.80, 79.97,
     * 82.07, 84.18, 86.23, 88.24, 90.26, 92.21, 94.10, 96.10, 97.99,
     * 99.89,101.73,103.57,105.36,107.18,108.97,110.83,112.64,114.25,
     *116.04,117.81,119.48,  0.00, 19.73, 27.24, 32.95, 37.81, 42.02,
     * 45.76, 49.27, 52.48, 55.56, 58.42, 61.14, 63.79, 66.33, 68.79,
     * 71.15, 73.46, 75.73, 77.94, 80.11, 82.23, 84.31, 86.34, 88.37,
     * 90.35, 92.33, 94.28, 96.26, 98.16,100.02,101.90,103.71,105.58,
     *107.35,109.11,110.98,112.75,114.49,116.20,117.93,119.69,  0.00,
     * 19.74, 27.24, 33.02, 37.85, 42.06, 45.81, 49.31, 52.59, 55.62,
     * 58.50, 61.25, 63.83, 66.42, 68.87, 71.20, 73.55, 75.80, 78.04,
     * 80.17, 82.34, 84.39, 86.45, 88.47, 90.47, 92.50, 94.47, 96.42,
     * 98.35,100.13,102.05,103.85,105.74,107.61,109.31,111.08,112.93,
     *114.60,116.41,118.14,119.88,  0.00, 19.79, 27.24, 32.99, 37.86,
     * 42.05, 45.86, 49.33, 52.58, 55.62, 58.52, 61.24, 63.88, 66.43,
     * 68.86, 71.25, 73.59, 75.87, 78.10, 80.24, 82.40, 84.50, 86.58,
     * 88.63, 90.58, 92.57, 94.57, 96.44, 98.42,100.35,102.15,103.97,
     *105.88,107.68,109.47,111.24,113.09,114.79,116.51,118.36,120.08,
     *  0.00, 19.75, 27.27, 33.02, 37.81, 42.01, 45.82, 49.31, 52.54,
     * 55.59, 58.46, 61.21, 63.88, 66.38, 68.87, 71.25, 73.58, 75.88/
      data (en(i),i=  4447,  4617)/
     * 78.09, 80.26, 82.45, 84.54, 86.61, 88.69, 90.68, 92.64, 94.65,
     * 96.57, 98.49,100.37,102.25,104.13,105.97,107.86,109.60,111.48,
     *113.25,115.00,116.78,118.47,120.23,  0.00, 19.75, 27.21, 32.94,
     * 37.73, 41.96, 45.78, 49.26, 52.47, 55.54, 58.41, 61.14, 63.81,
     * 66.36, 68.84, 71.21, 73.57, 75.83, 78.12, 80.27, 82.42, 84.60,
     * 86.67, 88.71, 90.73, 92.76, 94.69, 96.67, 98.64,100.47,102.44,
     *104.29,106.16,107.98,109.81,111.60,113.41,115.16,116.94,118.73,
     *120.45,  0.00, 19.72, 27.12, 32.87, 37.65, 41.83, 45.63, 49.10,
     * 52.36, 55.38, 58.27, 61.01, 63.64, 66.25, 68.69, 71.11, 73.47,
     * 75.76, 78.04, 80.23, 82.43, 84.57, 86.65, 88.67, 90.78, 92.78,
     * 94.75, 96.78, 98.68,100.66,102.53,104.44,106.28,108.15,109.93,
     *111.71,113.57,115.38,117.12,118.86,120.72,  0.00, 19.62, 27.02,
     * 32.69, 37.47, 41.69, 45.49, 48.92, 52.16, 55.15, 58.16, 60.84,
     * 63.52, 66.08, 68.58, 71.05, 73.37, 75.71, 78.00, 80.20, 82.44,
     * 84.59, 86.71, 88.78, 90.85, 92.82, 94.91, 96.85, 98.82,100.77,
     *102.65,104.64,106.41,108.34,110.13,111.96,113.80,115.63,117.44,
     *119.17,120.95,  0.00, 19.54, 26.89, 32.51, 37.22, 41.41, 45.15,
     * 48.66, 51.86, 54.92, 57.80, 60.56, 63.25, 65.86, 68.37, 70.82,
     * 73.24, 75.57, 77.94, 80.14, 82.40, 84.59, 86.73, 88.80, 90.91/
      data (en(i),i=  4618,  4788)/
     * 92.98, 94.98, 96.99, 99.01,100.87,102.89,104.75,106.67,108.57,
     *110.38,112.24,114.02,115.88,117.65,119.42,121.20,  0.00, 19.38,
     * 26.62, 32.21, 36.90, 41.03, 44.77, 48.22, 51.47, 54.54, 57.47,
     * 60.29, 63.00, 65.65, 68.22, 70.71, 73.14, 75.53, 77.89, 80.20,
     * 82.41, 84.63, 86.82, 88.95, 91.08, 93.09, 95.15, 97.22, 99.24,
     *101.20,103.17,105.01,107.00,108.82,110.72,112.52,114.35,116.17,
     *118.03,119.87,121.65,  0.00, 19.08, 26.21, 31.74, 36.39, 40.48,
     * 44.24, 47.69, 51.00, 54.14, 57.11, 60.02, 62.85, 65.55, 68.13,
     * 70.69, 73.21, 75.60, 78.03, 80.35, 82.61, 84.81, 87.00, 89.17,
     * 91.32, 93.43, 95.47, 97.53, 99.54,101.57,103.57,105.44,107.36,
     *109.31,111.17,112.99,114.91,116.68,118.58,120.32,122.12,  0.00,
     * 18.33, 25.33, 30.88, 35.63, 39.89, 43.83, 47.46, 50.88, 54.09,
     * 57.21, 60.17, 63.04, 65.76, 68.43, 71.04, 73.55, 76.01, 78.42,
     * 80.81, 83.08, 85.32, 87.56, 89.78, 91.89, 94.02, 96.12, 98.18,
     *100.24,102.19,104.25,106.15,108.06,110.03,111.88,113.76,115.65,
     *117.45,119.31,121.12,122.86,  0.00, 17.32, 24.15, 29.50, 34.13,
     * 38.25, 41.96, 45.52, 48.78, 51.88, 54.85, 57.77, 60.49, 63.15,
     * 65.71, 68.32, 70.76, 73.17, 75.60, 77.89, 80.19, 82.40, 84.64,
     * 86.78, 88.95, 91.05, 93.13, 95.18, 97.16, 99.18,101.18,103.08/
      data (en(i),i=  4789,  4959)/
     *105.05,106.99,108.90,110.78,112.60,114.45,116.34,118.10,119.96,
     *  0.00, 18.06, 24.87, 30.34, 34.89, 38.93, 42.62, 46.08, 49.29,
     * 52.36, 55.26, 58.06, 60.80, 64.11, 66.03, 68.51, 70.74, 73.36,
     * 75.76, 78.09, 80.33, 82.52, 84.67, 86.79, 88.93, 91.11, 93.14,
     * 95.18, 97.21, 99.37,101.11,103.06,105.00,106.88,108.80,110.65,
     *112.45,114.26,116.11,117.91,119.74,  0.00, 18.62, 25.65, 31.11,
     * 35.73, 39.80, 43.41, 46.84, 50.03, 53.06, 55.89, 58.76, 61.40,
     * 64.05, 66.54, 68.95, 71.46, 73.75, 76.14, 78.49, 80.73, 82.95,
     * 85.06, 87.20, 89.34, 91.37, 93.39, 95.36, 97.39, 99.42,101.38,
     *103.24,105.17,106.99,108.95,110.84,112.62,114.51,116.27,118.05,
     *119.81,  0.00, 18.84, 26.01, 31.59, 36.26, 40.34, 44.05, 47.52,
     * 50.67, 53.67, 56.58, 59.30, 62.01, 64.63, 67.15, 69.61, 72.02,
     * 74.40, 76.68, 78.90, 81.14, 83.30, 85.45, 87.58, 89.68, 91.68,
     * 93.74, 95.69, 97.71, 99.68,101.57,103.54,105.44,107.28,109.13,
     *110.99,112.88,114.72,116.45,118.28,119.99,  0.00, 19.16, 26.45,
     * 32.02, 36.76, 40.91, 44.65, 48.10, 51.30, 54.33, 57.15, 59.92,
     * 62.59, 65.17, 67.72, 70.20, 72.55, 74.91, 77.21, 79.41, 81.57,
     * 83.75, 85.89, 88.05, 90.11, 92.11, 94.16, 96.09, 98.15,100.05,
     *101.98,103.88,105.80,107.59,109.48,111.28,113.08,114.94,116.68/
      data (en(i),i=  4960,  5130)/
     *118.53,120.21,  0.00, 19.38, 26.76, 32.43, 37.25, 41.40, 45.15,
     * 48.68, 51.83, 54.87, 57.75, 60.46, 63.15, 65.75, 68.25, 70.67,
     * 73.09, 75.40, 77.66, 79.92, 82.14, 84.31, 86.33, 88.48, 90.51,
     * 92.52, 94.52, 96.57, 98.49,100.38,102.33,104.17,106.05,107.90,
     *109.81,111.64,113.41,115.29,117.03,118.81,120.50,  0.00, 19.56,
     * 27.00, 32.77, 37.58, 41.83, 45.61, 49.09, 52.30, 55.34, 58.23,
     * 60.97, 63.69, 66.27, 68.73, 71.17, 73.59, 75.87, 78.17, 80.40,
     * 82.55, 84.73, 86.78, 88.89, 90.94, 92.93, 94.92, 96.89, 98.86,
     *100.75,102.66,104.53,106.44,108.30,110.09,111.98,113.71,115.57,
     *117.27,119.04,120.75,  0.00, 19.82, 27.24, 33.07, 37.91, 42.16,
     * 45.99, 49.51, 52.76, 55.76, 58.66, 61.44, 64.16, 66.74, 69.25,
     * 71.70, 74.07, 76.36, 78.58, 80.86, 82.97, 85.15, 87.17, 89.28,
     * 91.24, 93.37, 95.31, 97.27, 99.21,101.11,103.02,104.92,106.84,
     *108.56,110.49,112.24,114.05,115.80,117.57,119.30,121.10,  0.00,
     * 19.89, 27.50, 33.30, 38.20, 42.52, 46.36, 49.85, 53.10, 56.18,
     * 35.51, 40.10, 44.19, 67.16, 51.27, 72.06, 57.42, 60.28, 63.00,
     * 81.21, 83.38, 70.71, 87.62, 89.72, 77.75, 79.95, 82.15, 97.67,
     * 86.44,101.50,103.34, 92.55, 94.53,108.92,110.78,100.40,114.33,
     *116.11,106.05,119.62,121.38,  0.00, 20.04, 27.64, 33.49, 38.46/
      data (en(i),i=  5131,  5301)/
     * 42.79, 46.62, 50.17, 53.43, 56.50, 59.46, 62.23, 64.94, 67.53,
     * 70.05, 72.48, 74.79, 77.09, 79.38, 81.59, 83.77, 85.90, 87.93,
     * 90.01, 92.06, 94.04, 96.03, 97.95, 99.89,101.81,103.71,103.77,
     *107.45,109.25,111.10,112.78,114.60,116.44,118.13,119.86,121.66,
     *  0.00, 20.17, 27.83, 33.70, 38.66, 43.00, 46.90, 50.44, 53.75,
     * 56.79, 59.78, 61.30, 65.26, 67.85, 70.36, 72.82, 75.16, 77.49,
     * 79.70, 81.91, 84.11, 86.19, 88.31, 90.34, 92.36, 94.38, 96.34,
     * 98.27,100.22,102.13,103.95,105.84,107.66,109.53,111.36,113.16,
     *114.91,116.65,118.41,120.19,121.96,  0.00, 20.25, 27.95, 33.89,
     * 38.89, 43.22, 47.13, 50.68, 53.99, 57.09, 60.04, 62.86, 65.59,
     * 68.18, 70.68, 73.10, 75.46, 77.83, 79.99, 82.20, 84.41, 86.50,
     * 88.59, 90.61, 92.69, 94.68, 96.63, 98.58,100.48,102.43,104.27,
     *106.11,107.96,109.76,111.57,113.40,115.11,116.89,118.70,120.37,
     *122.14,  0.00, 20.37, 28.07, 34.02, 39.01, 43.42, 47.33, 50.93,
     * 54.19, 57.32, 60.31, 63.10, 65.81, 68.44, 70.95, 73.37, 75.76,
     * 78.03, 80.33, 82.49, 84.65, 86.79, 88.84, 90.88, 92.88, 94.93,
     * 96.90, 98.83,100.70,102.59,104.56,106.41,108.22,109.98,111.81,
     *113.56,115.39,117.17,118.97,120.67,122.32,  0.00, 20.65, 28.18,
     * 34.14, 39.18, 43.55, 47.46, 51.08, 54.39, 57.56, 60.51, 63.33/
      data (en(i),i=  5302,  5472)/
     * 66.07, 68.64, 71.20, 73.61, 75.96, 78.25, 80.50, 82.70, 84.88,
     * 87.00, 89.07, 91.21, 93.16, 95.15, 97.12, 99.07,101.02,102.89,
     *104.78,106.59,108.43,110.22,112.04,113.90,115.61,117.34,119.08,
     *120.89,122.55,  0.00, 20.54, 28.25, 34.27, 39.27, 43.69, 47.61,
     * 51.19, 54.56, 57.71, 60.70, 63.55, 66.23, 68.83, 71.37, 73.77,
     * 76.15, 78.45, 80.72, 82.89, 85.08, 87.22, 89.29, 91.38, 93.34,
     * 95.37, 97.28, 99.28,101.14,103.05,104.92,106.77,108.67,110.48,
     *112.22,113.98,115.84,117.60,119.37,121.03,122.81,  0.00, 20.57,
     * 28.33, 34.32, 39.38, 43.74, 47.72, 51.33, 54.69, 57.84, 60.82,
     * 63.65, 66.36, 69.00, 71.49, 73.94, 76.29, 78.63, 80.89, 83.10,
     * 85.24, 87.39, 89.43, 91.53, 93.55, 95.50, 97.46, 99.41,101.35,
     *103.25,105.09,106.96,108.85,110.66,112.46,114.28,116.03,117.82,
     *119.48,121.25,122.99,  0.00, 20.59, 28.41, 34.41, 39.41, 43.83,
     * 47.78, 51.41, 54.76, 57.94, 60.91, 63.76, 66.51, 69.10, 71.61,
     * 74.06, 76.41, 78.73, 81.01, 83.20, 85.36, 87.42, 89.53, 91.65,
     * 93.64, 95.66, 97.64, 99.53,101.52,103.43,105.29,107.15,109.01,
     *110.73,112.63,114.36,116.20,118.00,119.66,121.39,123.15,  0.00,
     * 20.67, 28.37, 34.43, 39.51, 43.83, 47.84, 51.45, 54.80, 57.97,
     * 60.99, 63.83, 66.52, 69.14, 71.64, 74.13, 76.53, 78.82, 81.10/
      data (en(i),i=  5473,  5643)/
     * 83.27, 85.44, 87.57, 89.64, 91.74, 93.76, 95.80, 97.74, 99.66,
     *101.61,103.53,105.42,107.27,109.13,110.95,112.78,114.56,116.30,
     *118.13,119.88,121.57,123.33,  0.00, 20.69, 28.45, 34.44, 39.45,
     * 43.87, 47.83, 51.46, 54.84, 58.03, 61.01, 63.80, 66.57, 69.16,
     * 71.69, 74.17, 76.52, 78.90, 81.13, 83.35, 85.54, 87.65, 89.73,
     * 91.77, 93.81, 95.80, 97.80, 99.80,101.70,103.63,105.49,107.39,
     *109.24,111.06,112.91,114.67,116.48,118.23,119.99,121.73,123.50,
     *  0.00, 20.68, 28.39, 34.43, 39.47, 43.84, 47.80, 51.44, 54.83,
     * 57.96, 60.99, 63.85, 66.53, 69.15, 71.68, 74.15, 76.55, 78.83,
     * 81.11, 83.36, 85.56, 87.70, 89.77, 91.85, 93.88, 95.90, 97.88,
     * 99.88,101.80,103.67,105.58,107.48,109.36,111.20,112.96,114.82,
     *116.56,118.42,120.10,121.94,123.59,  0.00, 20.67, 28.38, 34.40,
     * 39.39, 43.80, 47.71, 51.38, 54.76, 57.95, 60.89, 63.75, 66.49,
     * 69.07, 71.61, 74.08, 76.50, 78.84, 81.08, 83.30, 85.55, 87.64,
     * 89.76, 91.83, 93.89, 95.93, 97.89, 99.89,101.83,103.76,105.69,
     *107.57,109.45,111.29,113.07,114.89,116.76,118.47,120.26,122.04,
     *123.73,  0.00, 20.60, 28.32, 34.29, 39.31, 43.68, 47.62, 51.28,
     * 54.64, 57.79, 60.77, 63.63, 66.38, 68.99, 71.51, 73.99, 76.38,
     * 78.76, 81.05, 83.28, 85.45, 87.64, 89.74, 91.74, 93.90, 95.91/
      data (en(i),i=  5644,  5814)/
     * 97.95, 99.95,101.89,103.81,105.72,107.68,109.55,111.39,113.24,
     *115.02,116.83,118.64,120.39,122.13,123.96,  0.00, 20.57, 28.25,
     * 34.20, 39.15, 43.52, 47.47, 51.11, 54.48, 57.62, 60.61, 63.49,
     * 66.17, 68.82, 71.36, 73.84, 76.28, 78.64, 80.96, 83.15, 85.40,
     * 87.56, 89.69, 91.81, 93.91, 95.92, 97.97,100.00,101.93,103.93,
     *105.82,107.75,109.59,111.47,113.29,115.21,116.93,118.77,120.54,
     *122.34,124.09,  0.00, 20.50, 28.16, 34.06, 38.97, 43.34, 47.28,
     * 50.88, 54.23, 57.37, 60.37, 63.24, 65.95, 68.59, 71.17, 73.71,
     * 76.08, 78.47, 80.81, 83.11, 85.30, 87.54, 89.65, 91.81, 93.88,
     * 95.91, 97.95,100.02,102.00,103.94,105.93,107.86,109.74,111.56,
     *113.45,115.34,117.14,118.93,120.75,122.53,124.32,  0.00, 20.37,
     * 27.99, 33.80, 38.72, 43.09, 46.99, 50.58, 53.92, 57.07, 60.05,
     * 62.92, 65.67, 68.38, 70.97, 73.48, 75.95, 78.32, 80.64, 82.99,
     * 85.23, 87.43, 89.67, 91.76, 93.90, 96.00, 98.00,100.04,102.08,
     *104.11,106.02,107.94,109.92,111.79,113.69,115.50,117.37,119.23,
     *121.01,122.85,124.55,  0.00, 20.20, 27.70, 33.50, 38.41, 42.66,
     * 46.55, 50.17, 53.51, 56.64, 59.68, 62.55, 65.37, 68.07, 70.69,
     * 73.23, 75.75, 78.18, 80.58, 82.88, 85.22, 87.46, 89.66, 91.86,
     * 93.96, 96.12, 98.14,100.27,102.27,104.27,106.20,108.15,110.10/
      data (en(i),i=  5815,  5985)/
     *112.04,113.89,115.83,117.64,119.55,121.29,123.17,124.94,  0.00,
     * 19.96, 27.08, 33.02, 37.90, 42.14, 45.99, 49.59, 52.98, 56.19,
     * 59.23, 62.20, 65.09, 67.87, 70.50, 73.16, 75.67, 78.17, 80.61,
     * 82.95, 85.30, 87.64, 89.81, 92.04, 94.23, 96.32, 98.44,100.52,
     *102.54,104.55,106.57,108.54,110.50,112.43,114.35,116.22,118.03,
     *119.99,121.75,123.64,125.38,  0.00, 19.20, 26.43, 32.17, 37.06,
     * 41.45, 45.46, 49.24, 52.74, 56.07, 59.23, 62.26, 65.18, 68.01,
     * 70.77, 73.41, 75.99, 78.51, 80.98, 83.40, 85.74, 88.05, 90.35,
     * 92.52, 94.69, 96.91, 98.99,101.06,103.15,105.19,107.20,109.21,
     *111.16,113.11,115.06,116.88,118.79,120.64,122.55,124.31,126.13,
     *  0.00, 17.94, 24.97, 30.52, 35.23, 39.48, 43.35, 46.93, 50.32,
     * 53.51, 56.54, 59.43, 62.27, 64.97, 67.61, 70.16, 72.70, 75.13,
     * 77.59, 79.94, 82.25, 84.53, 86.75, 88.98, 91.16, 93.27, 95.32,
     * 97.48, 99.54,101.58,103.52,105.54,107.43,109.44,111.32,113.25,
     *115.10,116.96,118.87,120.64,122.47,  0.00, 18.70, 25.83, 31.32,
     * 36.04, 40.24, 44.05, 47.56, 50.90, 54.00, 57.04, 59.89, 62.65,
     * 65.30, 67.93, 70.30, 72.95, 75.28, 77.76, 80.11, 82.43, 84.62,
     * 86.89, 89.12, 91.17, 93.34, 95.42, 97.45, 99.55,101.51,103.50,
     *105.46,107.48,109.26,111.25,113.13,115.06,116.90,118.69,120.49/
      data (en(i),i=  5986,  6156)/
     *122.31,  0.00, 19.13, 26.43, 31.03, 36.80, 40.25, 44.09, 47.65,
     * 51.59, 54.73, 57.68, 60.54, 62.73, 65.87, 68.61, 70.36, 73.44,
     * 75.38, 77.64, 80.58, 82.87, 85.10, 86.81, 89.06, 91.56, 93.65,
     * 95.77, 97.42, 99.91,101.37,103.74,105.69,107.67,109.18,111.44,
     *113.27,115.21,117.06,118.85,120.26,122.14,  0.00, 19.48, 26.91,
     * 32.61, 37.44, 41.66, 45.51, 49.02, 52.31, 55.39, 58.38, 61.19,
     * 63.92, 66.51, 69.14, 71.59, 74.15, 76.48, 78.83, 81.15, 83.36,
     * 85.53, 87.72, 89.95, 91.99, 94.08, 96.11, 98.12,100.24,102.19,
     *104.06,106.09,107.99,109.83,111.71,113.65,115.44,117.23,119.07,
     *120.90,122.78,  0.00, 19.78, 27.28, 33.04, 37.94, 42.21, 46.03,
     * 49.64, 52.94, 56.06, 58.98, 61.79, 64.49, 67.12, 69.72, 72.20,
     * 74.68, 77.02, 79.36, 81.62, 83.90, 86.06, 88.25, 90.33, 92.50,
     * 94.47, 96.59, 98.64,100.65,102.54,104.47,106.38,108.36,110.23,
     *112.11,113.94,115.74,117.58,119.37,121.17,122.98,  0.00, 20.04,
     * 27.62, 33.46, 38.40, 42.68, 46.57, 50.19, 53.51, 56.60, 59.52,
     * 62.40, 65.04, 67.71, 70.29, 72.78, 75.23, 77.60, 79.89, 82.20,
     * 84.42, 86.61, 88.75, 90.87, 92.98, 94.92, 97.05, 99.00,101.05,
     *102.95,104.85,106.85,108.72,110.53,112.47,114.32,116.14,117.94,
     *119.75,121.46,123.27,  0.00, 20.22, 27.86, 33.86, 38.78, 43.12/
      data (en(i),i=  6157,  6327)/
     * 47.09, 50.63, 54.07, 57.16, 60.06, 62.90, 65.64, 68.31, 70.82,
     * 73.35, 75.75, 78.15, 80.38, 82.63, 84.92, 87.12, 89.26, 91.30,
     * 93.45, 95.47, 97.52, 99.44,101.43,103.34,105.33,107.21,109.10,
     *111.03,112.76,114.69,116.53,118.23,120.03,121.81,123.56,  0.00,
     * 20.39, 28.20, 34.10, 39.10, 43.51, 47.51, 51.11, 54.49, 57.61,
     * 60.53, 63.37, 66.10, 68.80, 71.36, 73.82, 76.25, 78.58, 80.89,
     * 83.22, 85.37, 87.56, 89.71, 91.83, 93.87, 95.83, 97.86, 99.81,
     *101.84,103.77,105.69,107.57,109.50,111.30,113.19,115.02,116.75,
     *118.61,120.38,122.17,123.90,  0.00, 20.57, 28.37, 34.37, 39.41,
     * 43.83, 47.82, 51.49, 54.84, 58.00, 60.96, 63.85, 66.55, 69.20,
     * 71.77, 74.27, 76.70, 79.06, 81.35, 83.59, 85.80, 87.98, 90.08,
     * 92.16, 94.26, 96.29, 98.30,100.29,102.19,104.16,106.06,107.97,
     *109.87,111.70,113.53,115.32,117.14,118.90,120.67,122.45,124.13,
     *  0.00, 20.70, 28.55, 34.57, 39.68, 44.12, 48.11, 51.80, 55.21,
     * 58.34, 61.34, 64.21, 66.98, 69.61, 72.18, 74.70, 77.08, 79.43,
     * 81.74, 83.99, 86.23, 88.35, 90.47, 92.58, 94.66, 96.67, 98.66,
     *100.65,102.58,104.53,106.44,108.27,110.21,112.08,113.86,115.72,
     *117.47,119.22,120.99,122.82,124.53,  0.00, 20.95, 28.54, 34.85,
     * 39.94, 44.37, 48.43, 52.13, 55.50, 58.67, 61.71, 64.60, 67.31/
      data (en(i),i=  6328,  6498)/
     * 70.03, 72.58, 75.05, 77.49, 79.79, 82.11, 84.36, 86.56, 88.77,
     * 90.85, 93.00, 94.99, 97.04, 99.03,100.98,102.96,104.90,106.72,
     *108.67,110.50,112.40,114.13,115.97,117.76,119.63,121.36,123.11,
     *124.85,  0.00, 20.93, 28.81, 34.98, 40.10, 44.61, 48.65, 52.37,
     * 55.75, 59.00, 61.99, 64.90, 67.66, 70.30, 72.88, 75.42, 77.82,
     * 80.18, 82.51, 84.73, 86.92, 89.05, 91.20, 93.32, 95.37, 97.38,
     * 99.31,101.33,103.25,105.18,107.10,108.98,110.80,112.67,114.50,
     *116.35,118.09,119.87,121.69,123.41,125.07,  0.00, 21.07, 29.01,
     * 35.17, 40.28, 44.80, 48.88, 52.58, 56.03, 59.27, 62.28, 65.15,
     * 67.95, 70.64, 73.23, 75.71, 78.15, 80.53, 82.80, 85.03, 87.30,
     * 89.39, 91.53, 93.64, 95.63, 97.65, 99.65,101.65,103.54,105.47,
     *107.37,109.24,111.11,112.93,114.74,116.56,118.32,120.14,121.94,
     *123.68,125.38,  0.00, 21.16, 29.12, 35.25, 40.44, 44.98, 49.05,
     * 52.84, 56.23, 59.49, 62.47, 65.42, 68.21, 70.93, 73.50, 76.00,
     * 78.40, 80.76, 83.07, 85.33, 87.52, 89.69, 91.77, 93.90, 95.96,
     * 97.92, 99.90,101.89,103.81,105.77,107.65,109.49,111.39,113.17,
     *114.98,116.81,118.64,120.41,122.11,123.83,125.64,  0.00, 21.29,
     * 29.22, 35.41, 40.59, 45.14, 49.24, 52.99, 56.43, 59.65, 62.69,
     * 65.63, 68.45, 71.15, 73.74, 76.21, 78.64, 81.00, 83.30, 85.55/
      data (en(i),i=  6499,  6669)/
     * 87.75, 89.92, 92.07, 94.15, 96.16, 98.18,100.21,102.13,104.07,
     *105.97,107.86,109.74,111.58,113.47,115.26,117.02,118.88,120.58,
     *122.37,124.15,125.93,  0.00, 21.32, 29.34, 35.48, 40.72, 45.27,
     * 49.40, 53.14, 56.56, 59.81, 62.89, 65.81, 68.62, 71.35, 73.92,
     * 76.41, 78.88, 81.21, 83.51, 85.77, 87.99, 90.15, 92.24, 94.36,
     * 96.39, 98.40,100.43,102.37,104.27,106.21,108.07,109.98,111.82,
     *113.70,115.51,117.25,119.06,120.85,122.58,124.39,126.07,  0.00,
     * 21.38, 29.36, 35.59, 40.79, 45.40, 49.47, 53.24, 56.70, 59.98,
     * 63.04, 65.98, 68.81, 71.52, 74.07, 76.62, 79.01, 81.42, 83.69,
     * 85.94, 88.16, 90.32, 92.46, 94.51, 96.58, 98.59,100.61,102.56,
     *104.44,106.37,108.32,110.20,112.02,113.81,115.63,117.50,119.25,
     *121.04,122.80,124.52,126.22,  0.00, 21.45, 29.46, 35.68, 40.88,
     * 45.45, 49.58, 53.32, 56.82, 60.08, 63.15, 66.11, 68.91, 71.62,
     * 74.20, 76.73, 79.15, 81.52, 83.81, 86.11, 88.28, 90.48, 92.58,
     * 94.66, 96.79, 98.70,100.71,102.66,104.66,106.58,108.42,110.32,
     *112.16,113.96,115.81,117.64,119.32,121.21,122.96,124.76,126.42,
     *  0.00, 21.47, 29.50, 35.75, 40.91, 45.53, 49.59, 53.37, 56.90,
     * 60.15, 63.28, 66.18, 69.00, 71.71, 74.31, 76.85, 79.25, 81.63,
     * 83.96, 86.18, 88.38, 90.60, 92.72, 94.77, 96.82, 98.85,100.82/
      data (en(i),i=  6670,  6840)/
     *102.81,104.79,106.64,108.55,110.45,112.33,114.16,115.97,117.81,
     *119.55,121.30,123.13,124.82,126.57,  0.00, 21.48, 29.49, 35.71,
     * 40.98, 45.53, 49.64, 53.39, 56.89, 60.17, 63.31, 66.23, 69.09,
     * 71.75, 74.34, 76.87, 79.34, 81.70, 84.02, 86.29, 88.46, 90.67,
     * 92.76, 94.88, 96.89, 98.92,100.96,102.92,104.78,106.78,108.75,
     *110.59,112.44,114.26,116.06,117.86,119.68,121.49,123.26,125.02,
     *126.73,  0.00, 21.52, 29.53, 35.76, 40.97, 45.56, 49.66, 53.42,
     * 56.94, 60.23, 63.28, 66.26, 69.08, 71.78, 74.40, 76.90, 79.36,
     * 81.73, 84.02, 86.32, 88.49, 90.69, 92.86, 94.92, 97.02, 99.00,
     *101.06,102.98,104.96,106.91,108.73,110.70,112.49,114.41,116.24,
     *117.95,119.84,121.64,123.39,125.14,126.83,  0.00, 21.48, 29.52,
     * 35.76, 40.96, 45.51, 49.65, 53.39, 56.90, 60.17, 63.31, 66.26,
     * 69.08, 71.75, 74.36, 76.91, 79.36, 81.72, 84.03, 86.32, 88.50,
     * 90.69, 92.86, 94.93, 96.97, 99.07,101.06,103.05,105.02,106.96,
     *108.86,110.74,112.58,114.43,116.29,118.13,119.95,121.67,123.49,
     *125.27,126.93,  0.00, 21.49, 29.47, 35.71, 40.93, 45.49, 49.57,
     * 53.34, 56.89, 60.17, 63.24, 66.17, 69.04, 71.71, 74.32, 76.88,
     * 79.31, 81.65, 83.96, 86.25, 88.49, 90.70, 92.81, 95.00, 97.05,
     * 99.04,101.06,103.08,105.00,106.99,108.92,110.82,112.69,114.52/
      data (en(i),i=  6841,  7011)/
     *116.41,118.16,119.97,121.81,123.61,125.39,127.15,  0.00, 21.50,
     * 29.49, 35.62, 40.84, 45.42, 49.49, 53.25, 56.76, 60.02, 63.16,
     * 66.10, 68.90, 71.60, 74.25, 76.75, 79.20, 81.62, 83.90, 86.23,
     * 88.41, 90.65, 92.79, 94.87, 97.01, 98.92,100.96,103.09,105.04,
     *106.97,108.90,110.84,112.67,114.61,116.47,118.22,120.07,121.87,
     *123.75,125.46,127.21,  0.00, 21.42, 29.41, 35.57, 40.71, 45.28,
     * 49.34, 53.12, 56.62, 59.91, 62.96, 65.95, 68.74, 71.47, 74.07,
     * 76.61, 79.07, 81.47, 83.82, 86.09, 88.32, 90.50, 92.74, 94.80,
     * 96.95, 98.99,101.02,103.03,105.06,106.98,108.89,110.82,112.78,
     *114.64,116.51,118.35,120.18,121.97,123.84,125.53,127.34,  0.00,
     * 21.37, 29.33, 35.45, 40.57, 45.10, 49.19, 52.99, 56.42, 59.71,
     * 62.45, 65.72, 68.56, 71.28, 73.89, 76.44, 78.88, 81.29, 83.69,
     * 85.99, 88.18, 90.44, 92.62, 94.76, 96.82, 98.94,100.99,103.05,
     *105.05,107.04,108.94,110.87,112.83,114.75,116.58,118.46,120.29,
     *122.06,123.92,125.70,127.51,  0.00, 21.30, 29.14, 35.27, 40.41,
     * 44.87, 48.98, 52.72, 56.15, 59.46, 62.53, 65.45, 68.29, 71.01,
     * 73.65, 76.22, 78.67, 81.05, 83.49, 85.81, 88.10, 90.33, 92.51,
     * 94.69, 96.77, 98.92,100.97,103.02,104.99,106.99,109.02,110.97,
     *112.83,114.78,116.71,118.54,120.43,122.22,124.06,125.87,127.64/
      data (en(i),i=  7012,  7182)/
     *  0.00, 21.16, 28.98, 35.04, 40.14, 44.58, 48.65, 52.36, 55.84,
     * 59.07, 62.17, 65.13, 67.98, 70.67, 73.35, 75.92, 78.43, 80.90,
     * 83.28, 85.61, 87.91, 90.25, 92.40, 94.62, 96.76, 98.91,100.95,
     *103.06,105.09,107.05,109.07,111.04,112.94,114.91,116.84,118.72,
     *120.56,122.42,124.23,126.03,127.85,  0.00, 20.98, 28.75, 34.76,
     * 39.76, 44.19, 48.24, 51.96, 55.42, 58.61, 61.76, 64.74, 67.58,
     * 70.36, 73.00, 75.66, 78.17, 80.70, 83.09, 85.50, 87.83, 90.10,
     * 92.40, 94.57, 96.80, 98.90,101.03,103.10,105.15,107.26,109.25,
     *111.20,113.22,115.10,117.00,118.89,120.78,122.64,124.54,126.37,
     *128.14,  0.00, 20.71, 28.32, 34.28, 39.19, 43.69, 47.64, 51.33,
     * 54.84, 58.11, 61.30, 64.31, 67.26, 70.05, 72.75, 75.48, 78.07,
     * 80.61, 83.12, 85.53, 87.89, 90.18, 92.49, 94.71, 96.88, 99.13,
     *101.24,103.37,105.42,107.46,109.50,111.54,113.50,115.44,117.44,
     *119.27,121.21,123.09,124.87,126.75,128.62,  0.00, 19.96, 27.47,
     * 33.32, 38.39, 42.90, 47.03, 50.89, 54.47, 57.90, 61.17, 64.24,
     * 67.24, 70.13, 72.95, 75.65, 78.31, 80.86, 83.37, 85.85, 88.25,
     * 90.59, 92.94, 95.16, 97.42, 99.60,101.79,103.84,105.96,108.09,
     *110.11,112.10,114.09,116.04,118.00,119.91,121.83,123.73,125.62,
     *127.48,129.33,  0.00, 18.59, 25.81, 31.47, 36.33, 40.63, 44.62/
      data (en(i),i=  7183,  7353)/
     * 48.28, 51.74, 54.97, 58.09, 61.11, 63.98, 66.75, 69.40, 72.02,
     * 74.58, 77.02, 79.49, 81.85, 84.21, 86.54, 88.77, 91.06, 93.23,
     * 95.38, 97.51, 99.69,101.76,103.77,105.75,107.79,109.79,111.77,
     *113.73,115.57,117.52,119.40,121.23,123.05,124.96,  0.00, 19.31,
     * 26.65, 32.32, 37.15, 41.48, 45.33, 48.96, 52.33, 55.55, 58.66,
     * 61.57, 64.36, 67.10, 69.74, 72.38, 74.85, 77.31, 79.76, 82.11,
     * 84.43, 86.77, 89.00, 91.16, 93.36, 95.47, 97.67, 99.66,101.75,
     *103.80,105.85,107.81,109.77,111.72,113.65,115.58,117.42,119.33,
     *121.18,122.97,124.78,  0.00, 19.76, 27.20, 33.02, 37.89, 42.29,
     * 46.16, 49.72, 53.09, 56.28, 59.35, 62.24, 65.02, 67.76, 70.42,
     * 72.91, 75.44, 77.89, 80.27, 82.65, 84.88, 87.24, 89.36, 91.61,
     * 93.75, 95.91, 97.99,100.03,102.13,104.10,106.16,108.04,110.07,
     *112.03,113.91,115.76,117.66,119.52,121.33,123.20,124.98,  0.00,
     * 20.07, 27.72, 33.53, 38.49, 42.90, 46.80, 50.44, 53.61, 56.97,
     * 60.02, 62.90, 63.43, 68.41, 70.99, 73.49, 75.97, 78.47, 80.88,
     * 83.19, 85.48, 87.77, 89.90, 92.07, 94.27, 96.35, 98.40,100.45,
     *102.46,104.50,106.51,108.49,110.37,112.33,114.15,116.05,117.93,
     *119.84,121.65,123.44,125.20,  0.00, 20.37, 28.09, 34.06, 39.05,
     * 43.41, 47.40, 51.11, 54.51, 57.66, 60.67, 63.64, 66.44, 69.03/
      data (en(i),i=  7354,  7524)/
     * 71.60, 74.17, 76.61, 79.06, 81.45, 83.78, 86.04, 88.27, 90.47,
     * 92.56, 94.72, 96.79, 98.88,101.00,102.96,104.95,106.96,108.93,
     *110.83,112.72,114.59,116.49,118.33,120.16,121.99,123.75,125.57,
     *  0.00, 20.59, 28.38, 34.43, 39.45, 43.93, 47.92, 51.64, 55.04,
     * 58.24, 61.26, 64.21, 66.92, 69.63, 72.18, 74.76, 77.24, 79.64,
     * 82.02, 84.31, 86.55, 88.80, 90.95, 93.14, 95.23, 97.35, 99.35,
     *101.43,103.39,105.41,107.35,109.27,111.19,113.16,115.03,116.85,
     *118.67,120.55,122.27,124.19,125.90,  0.00, 20.75, 28.69, 34.72,
     * 39.85, 44.38, 48.36, 52.13, 55.56, 58.77, 61.81, 64.68, 67.52,
     * 70.15, 72.75, 75.30, 77.76, 80.25, 82.54, 84.83, 87.11, 89.34,
     * 91.49, 93.62, 95.74, 97.79, 99.85,101.92,103.85,105.85,107.82,
     *109.74,111.65,113.49,115.42,117.29,119.08,120.92,122.69,124.57,
     *126.33,  0.00, 21.05, 28.94, 35.08, 40.20, 44.73, 48.78, 52.52,
     * 55.98, 59.23, 62.31, 65.23, 67.96, 70.68, 73.26, 75.79, 78.31,
     * 80.68, 83.04, 85.37, 87.62, 89.85, 91.99, 94.10, 96.25, 98.27,
     *100.27,102.33,104.33,106.27,108.20,110.17,112.00,113.94,115.85,
     *117.60,119.51,121.28,123.11,124.90,126.61,  0.00, 21.19, 29.27,
     * 35.42, 40.45, 45.08, 49.13, 52.92, 56.44, 59.63, 62.73, 65.61,
     * 68.45, 71.16, 73.78, 76.34, 78.77, 81.20, 83.48, 85.82, 88.10/
      data (en(i),i=  7525,  7695)/
     * 90.27, 92.45, 94.55, 96.70, 98.67,100.81,102.72,104.77,106.72,
     *108.68,110.64,112.49,114.32,116.19,118.03,119.77,121.63,123.40,
     *125.22,126.97,  0.00, 21.40, 29.35, 35.57, 40.74, 45.35, 49.50,
     * 53.29, 56.77, 60.09, 63.15, 66.06, 68.89, 71.55, 74.20, 76.75,
     * 79.25, 81.64, 83.95, 86.28, 88.49, 90.75, 92.87, 94.99, 97.11,
     * 99.17,101.15,103.21,105.18,107.11,109.07,110.93,112.84,114.70,
     *116.58,118.38,120.22,121.96,123.76,125.58,127.37,  0.00, 21.53,
     * 29.54, 35.79, 41.00, 45.67, 49.79, 53.63, 57.14, 60.42, 63.50,
     * 66.44, 69.27, 71.95, 74.61, 77.19, 79.68, 82.06, 84.38, 86.67,
     * 88.91, 91.11, 93.30, 95.40, 97.51, 99.53,101.61,103.54,105.54,
     *107.47,109.38,111.37,113.16,115.10,116.88,118.79,120.54,122.37,
     *124.13,125.83,127.70,  0.00, 21.63, 29.72, 35.96, 41.25, 45.90,
     * 50.05, 53.89, 57.43, 60.73, 63.79, 66.78, 69.60, 72.32, 74.98,
     * 77.54, 79.98, 82.42, 84.77, 87.08, 89.30, 91.54, 93.69, 95.78,
     * 97.83, 99.90,101.98,103.96,105.91,107.82,109.73,111.72,113.51,
     *115.38,117.22,119.07,120.84,122.64,124.40,126.21,127.99,  0.00,
     * 21.76, 29.90, 36.12, 41.49, 46.07, 50.32, 54.15, 57.69, 61.02,
     * 64.12, 67.08, 69.89, 72.66, 75.34, 77.89, 80.40, 82.75, 85.15,
     * 87.41, 89.68, 91.87, 93.99, 96.16, 98.23,100.23,102.21,104.30/
      data (en(i),i=  7696,  7866)/
     *106.24,108.19,110.12,112.00,113.86,115.75,117.56,119.40,121.25,
     *122.96,124.74,126.58,128.21,  0.00, 21.87, 29.97, 36.32, 41.65,
     * 46.33, 50.52, 54.36, 57.96, 61.28, 64.40, 67.34, 70.21, 72.99,
     * 75.59, 78.17, 80.65, 83.08, 85.44, 87.70, 90.00, 92.19, 94.33,
     * 96.42, 98.54,100.55,102.56,104.59,106.57,108.52,110.44,112.30,
     *114.22,116.05,117.88,119.60,121.46,123.30,125.03,126.78,128.49,
     *  0.00, 21.95, 30.11, 36.46, 41.80, 46.49, 50.70, 54.57, 58.14,
     * 61.47, 64.62, 67.60, 70.43, 73.21, 75.92, 78.47, 80.92, 83.38,
     * 85.70, 88.01, 90.24, 92.45, 94.63, 96.76, 98.82,100.89,102.87,
     *104.82,106.86,108.75,110.69,112.57,114.47,116.33,118.16,119.94,
     *121.73,123.57,125.32,127.08,128.86,  0.00, 22.01, 30.19, 36.59,
     * 41.91, 46.59, 50.88, 54.77, 58.31, 61.65, 64.82, 67.83, 70.70,
     * 73.48, 76.16, 78.69, 81.18, 83.58, 85.97, 88.31, 90.50, 92.74,
     * 94.88, 97.02, 99.10,101.12,103.17,105.16,107.07,109.01,110.93,
     *112.87,114.65,116.58,118.36,120.23,122.04,123.74,125.51,127.28,
     *129.00,  0.00, 23.22, 31.13, 36.68, 42.09, 46.77, 50.99, 54.89,
     * 58.90, 61.82, 64.99, 68.00, 70.87, 73.66, 76.36, 79.21, 81.38,
     * 83.85, 86.21, 88.75, 91.01, 93.18, 95.14, 97.49, 99.31,101.60,
     *103.59,105.37,107.50,109.27,111.15,113.03,115.14,116.74,118.83/
      data (en(i),i=  7867,  8037)/
     *120.62,122.23,123.99,126.00,127.57,129.24,  0.00, 22.26, 30.49,
     * 36.84, 42.20, 46.92, 51.21, 55.09, 58.62, 61.95, 65.14, 68.15,
     * 71.06, 73.85, 76.52, 79.11, 81.62, 84.01, 86.43, 88.73, 90.90,
     * 93.15, 95.30, 97.45, 99.49,101.55,103.63,105.56,107.56,109.48,
     *111.36,113.29,115.15,117.05,118.82,120.64,122.43,124.19,126.00,
     *127.73,129.44,  0.00, 22.15, 30.48, 36.82, 42.23, 46.94, 51.18,
     * 55.13, 58.73, 62.10, 65.24, 68.32, 71.16, 73.98, 76.66, 79.24,
     * 81.74, 84.18, 86.56, 88.82, 91.07, 93.25, 95.45, 97.55, 99.64,
     *101.74,103.73,105.71,107.68,109.65,111.56,113.45,115.30,117.10,
     *118.96,120.78,122.58,124.42,126.13,127.89,129.63,  0.00, 22.20,
     * 30.52, 36.89, 42.32, 47.07, 51.29, 55.22, 58.79, 62.17, 65.35,
     * 68.44, 71.29, 74.10, 76.76, 79.35, 81.88, 84.30, 86.62, 88.94,
     * 91.21, 93.46, 95.63, 97.71, 99.74,101.84,103.89,105.90,107.85,
     *109.78,111.69,113.59,115.50,117.33,119.14,120.95,122.77,124.50,
     *126.31,128.10,129.86,  0.00, 22.20, 30.50, 36.95, 42.34, 47.10,
     * 51.34, 55.25, 58.85, 62.20, 65.45, 68.46, 71.40, 74.15, 76.82,
     * 79.46, 81.90, 84.38, 86.76, 89.06, 91.30, 93.55, 95.73, 97.83,
     * 99.91,101.94,104.00,105.97,107.96,109.88,111.79,113.69,115.58,
     *117.47,119.31,121.11,122.93,124.65,126.48,128.16,130.00,  0.00/
      data (en(i),i=  8038,  8208)/
     * 22.28, 30.42, 36.88, 42.35, 47.03, 51.36, 55.28, 58.90, 62.25,
     * 65.50, 68.53, 71.46, 74.25, 76.88, 79.49, 81.99, 84.43, 86.79,
     * 89.11, 91.40, 93.61, 95.74, 97.91, 99.95,102.03,104.07,106.03,
     *108.15,109.98,111.86,113.78,115.66,117.51,119.37,121.17,123.00,
     *124.80,126.53,128.32,130.07,  0.00, 22.31, 30.57, 37.02, 42.39,
     * 47.14, 51.38, 55.30, 58.88, 62.27, 65.51, 68.54, 71.48, 74.29,
     * 76.94, 79.49, 82.02, 84.47, 86.85, 89.13, 91.41, 93.65, 95.79,
     * 97.95, 99.99,102.11,104.11,106.10,108.12,110.07,111.98,113.86,
     *115.80,117.65,119.48,121.29,123.08,124.92,126.68,128.46,130.18,
     *  0.00, 22.28, 30.56, 36.99, 42.38, 47.10, 51.37, 55.27, 58.84,
     * 62.25, 65.45, 68.51, 71.42, 74.24, 76.94, 79.49, 82.03, 84.45,
     * 86.84, 89.18, 91.39, 93.62, 95.79, 97.91,100.01,102.11,104.17,
     *106.15,108.13,110.12,112.00,113.90,115.81,117.67,119.52,121.39,
     *123.26,124.94,126.81,128.54,130.30,  0.00, 22.31, 30.57, 36.92,
     * 42.34, 47.04, 51.30, 55.21, 58.80, 62.25, 65.46, 68.48, 71.36,
     * 74.19, 76.89, 79.42, 81.97, 84.44, 86.79, 89.14, 91.37, 93.61,
     * 95.83, 97.90,100.06,102.10,104.15,106.15,108.17,110.10,112.02,
     *113.94,115.84,117.68,119.56,121.48,123.22,125.09,126.89,128.64,
     *130.38,  0.00, 22.27, 30.50, 36.93, 42.25, 47.03, 51.19, 55.14/
      data (en(i),i=  8209,  8379)/
     * 58.76, 62.16, 65.35, 68.41, 71.26, 74.07, 76.76, 79.40, 81.90,
     * 84.31, 86.73, 89.05, 91.32, 93.55, 95.69, 97.89, 99.99,102.03,
     *104.07,106.17,108.04,110.10,112.07,114.00,115.91,117.79,119.62,
     *121.46,123.34,125.06,126.95,128.71,130.46,  0.00, 22.27, 30.48,
     * 36.62, 42.16, 46.90, 51.13, 54.89, 58.65, 62.00, 65.10, 68.15,
     * 71.17, 73.92, 76.61, 79.24, 81.78, 84.13, 86.57, 88.90, 91.20,
     * 93.38, 95.63, 97.76, 99.93,102.02,104.06,106.11,108.06,109.98,
     *111.99,113.95,115.83,117.76,119.63,121.51,123.35,125.10,126.84,
     *128.76,130.57,  0.00, 22.18, 30.42, 36.75, 42.10, 46.72, 50.96,
     * 54.87, 58.48, 61.89, 65.02, 68.11, 71.01, 73.74, 76.47, 79.07,
     * 81.56, 84.01, 86.44, 88.73, 91.09, 93.35, 95.54, 97.68, 99.77,
     *101.91,104.02,105.99,108.08,110.05,111.99,113.92,115.90,117.81,
     *119.68,121.55,123.39,125.19,126.97,128.80,130.61,  0.00, 22.12,
     * 30.26, 36.60, 41.90, 46.56, 50.80, 54.66, 58.27, 61.66, 64.85,
     * 67.84, 70.73, 73.57, 76.22, 78.84, 81.37, 83.81, 86.24, 88.56,
     * 90.88, 93.13, 95.40, 97.53, 99.70,101.80,103.90,105.95,107.97,
     *109.99,112.00,113.92,115.90,117.77,119.65,121.60,123.45,125.27,
     *127.08,128.95,130.73,  0.00, 22.31, 30.11, 36.42, 41.70, 46.34,
     * 50.54, 54.38, 57.98, 61.45, 64.55, 67.55, 70.47, 73.24, 75.97/
      data (en(i),i=  8380,  8550)/
     * 78.62, 81.13, 83.67, 86.00, 88.39, 90.69, 93.04, 95.30, 97.50,
     * 99.65,101.75,103.80,105.95,107.91,109.97,112.00,113.99,115.88,
     *117.85,119.74,121.63,123.50,125.41,127.24,129.03,130.96,  0.00,
     * 21.94, 29.93, 36.17, 41.41, 46.07, 50.24, 54.08, 57.61, 61.00,
     * 64.14, 67.16, 70.12, 72.88, 75.59, 78.24, 80.84, 83.32, 85.70,
     * 88.17, 90.49, 92.76, 95.05, 97.26, 99.49,101.67,103.76,105.84,
     *107.90,109.99,112.02,113.97,115.98,117.92,119.81,121.71,123.60,
     *125.45,127.32,129.16,131.08,  0.00, 21.79, 29.72, 35.88, 41.04,
     * 45.71, 49.82, 53.64, 57.19, 60.51, 63.72, 66.75, 69.70, 72.53,
     * 75.23, 77.88, 80.54, 83.04, 85.55, 87.96, 90.35, 92.72, 94.96,
     * 97.28, 99.47,101.67,103.78,105.93,108.02,110.07,112.08,114.09,
     *116.11,118.01,120.06,121.96,123.79,125.69,127.56,129.43,131.33,
     *  0.00, 21.52, 29.33, 35.40, 40.53, 45.09, 49.08, 53.02, 56.58,
     * 59.94, 63.18, 66.29, 69.23, 72.18, 74.92, 77.69, 80.32, 82.92,
     * 85.42, 87.90, 90.34, 92.67, 94.99, 97.31, 99.56,101.73,103.96,
     *106.06,108.19,110.29,112.33,114.36,116.38,118.34,120.27,122.25,
     *124.23,126.04,127.99,129.80,131.66,  0.00, 20.76, 28.40, 34.48,
     * 39.67, 44.30, 48.49, 52.45, 56.16, 59.67, 62.96, 66.13, 69.25,
     * 72.17, 75.04, 77.82, 80.51, 83.09, 85.65, 88.16, 90.66, 93.02/
      data (en(i),i=  8551,  8721)/
     * 95.37, 97.73, 99.97,102.23,104.37,106.57,108.68,110.77,112.88,
     *114.89,116.90,118.89,120.91,122.84,124.79,126.68,128.57,130.43,
     *132.33,  0.00, 19.17, 26.54, 32.33, 37.33, 41.71, 45.90, 49.55,
     * 53.08, 56.45, 59.58, 62.66, 65.57, 68.43, 71.15, 73.80, 76.42,
     * 78.88, 81.32, 83.74, 86.09, 88.45, 90.73, 93.01, 95.30, 97.45,
     * 99.64,101.73,103.91,105.91,107.96,110.04,112.06,113.98,115.96,
     *117.93,119.77,121.66,123.60,125.45,127.34,  0.00, 19.90, 27.39,
     * 33.23, 38.16, 42.60, 46.53, 50.24, 53.71, 57.01, 60.17, 63.12,
     * 66.03, 68.81, 71.49, 74.11, 76.63, 79.20, 81.57, 84.01, 86.34,
     * 88.68, 90.95, 93.23, 95.32, 97.60, 99.71,101.86,103.92,105.97,
     *108.06,110.05,112.06,113.92,115.89,117.91,119.83,121.69,123.48,
     *125.44,127.20,  0.00, 20.31, 27.97, 33.89, 38.95, 43.34, 47.34,
     * 51.04, 54.47, 57.75, 60.84, 63.91, 66.76, 69.48, 72.18, 74.76,
     * 77.28, 79.72, 82.19, 84.55, 86.84, 89.21, 91.47, 93.70, 95.81,
     * 98.01,100.15,102.27,104.23,106.33,108.40,110.35,112.38,114.25,
     *116.25,118.09,120.04,121.90,123.80,125.56,127.42,  0.00, 20.66,
     * 28.43, 34.45, 39.52, 44.03, 48.04, 51.78, 55.22, 58.50, 61.64,
     * 64.55, 67.45, 70.21, 72.83, 75.35, 77.88, 80.34, 82.80, 85.15,
     * 87.48, 89.70, 91.98, 94.20, 96.37, 98.45,100.66,102.70,104.71/
      data (en(i),i=  8722,  8892)/
     *106.81,108.82,110.77,112.68,114.72,116.60,118.47,120.42,122.22,
     *124.12,125.90,127.77,  0.00, 20.97, 28.84, 34.93, 40.07, 44.61,
     * 48.61, 52.43, 55.88, 59.16, 62.31, 65.24, 68.11, 70.77, 73.43,
     * 76.00, 78.51, 80.99, 83.41, 85.78, 88.11, 90.27, 92.57, 94.73,
     * 96.84, 99.05,101.12,103.16,105.27,107.26,109.26,111.18,113.16,
     *115.08,117.03,118.86,120.78,122.60,124.45,126.25,128.08,  0.00,
     * 21.23, 29.13, 35.38, 40.50, 45.05, 49.20, 52.99, 56.45, 59.75,
     * 62.91, 65.89, 68.73, 71.45, 74.08, 76.64, 79.12, 81.55, 84.02,
     * 86.33, 88.66, 90.94, 93.07, 95.31, 97.48, 99.60,101.61,103.68,
     *105.71,107.76,109.70,111.71,113.56,115.57,117.45,119.31,121.21,
     *123.09,124.82,126.68,128.54,  0.00, 21.40, 29.47, 35.67, 40.87,
     * 45.52, 49.63, 53.45, 57.00, 60.34, 63.46, 66.41, 69.27, 71.98,
     * 74.61, 77.20, 79.68, 82.17, 84.57, 86.91, 89.15, 91.46, 93.71,
     * 95.86, 98.00,100.07,102.21,104.23,106.24,108.22,110.20,112.16,
     *114.03,115.97,117.85,119.77,121.60,123.43,125.26,127.08,128.84,
     *  0.00, 21.65, 29.70, 35.96, 41.26, 45.91, 50.06, 53.89, 57.45,
     * 60.79, 63.94, 66.93, 69.81, 72.56, 75.20, 77.75, 80.27, 82.69,
     * 85.13, 87.46, 89.73, 92.00, 94.21, 96.33, 98.51,100.61,102.66,
     *104.72,106.71,108.72,110.70,112.60,114.55,116.46,118.28,120.14/
      data (en(i),i=  8893,  9063)/
     *122.05,123.88,125.65,127.40,129.21,  0.00, 21.78, 29.91, 36.25,
     * 41.58, 46.26, 50.46, 54.31, 57.89, 61.25, 64.42, 67.38, 70.27,
     * 73.01, 75.67, 78.22, 80.74, 83.24, 85.57, 87.94, 90.22, 92.49,
     * 94.70, 96.83, 98.98,101.08,103.15,105.15,107.15,109.13,111.14,
     *113.05,114.93,116.90,118.77,120.57,122.41,124.31,126.10,127.88,
     *129.68,  0.00, 21.94, 30.09, 36.51, 41.82, 46.60, 51.03, 54.88,
     * 58.29, 61.65, 64.86, 67.82, 70.72, 73.61, 76.08, 78.71, 81.24,
     * 83.67, 86.09, 88.37, 90.70, 93.05, 95.13, 97.32, 99.40,101.51,
     *103.61,105.62,107.71,109.54,111.56,113.43,115.39,117.38,119.13,
     *120.97,122.90,124.56,126.45,128.20,130.04,  0.00, 22.13, 30.34,
     * 36.71, 42.12, 46.81, 51.13, 55.04, 58.64, 62.02, 65.22, 68.24,
     * 71.12, 73.91, 76.54, 79.10, 81.65, 84.15, 86.51, 88.89, 91.11,
     * 93.40, 95.55, 97.71, 99.83,101.92,104.02,106.05,107.99,110.02,
     *111.94,113.85,115.74,117.69,119.48,121.34,123.18,125.05,126.85,
     *128.67,130.36,  0.00, 22.25, 30.50, 36.92, 42.32, 47.08, 51.38,
     * 55.30, 58.93, 62.38, 65.59, 68.61, 71.46, 74.26, 76.96, 79.51,
     * 82.13, 84.52, 86.92, 89.27, 91.60, 93.83, 96.00, 98.16,100.24,
     *102.39,104.37,106.47,108.46,110.40,112.41,114.35,116.26,118.09,
     *119.94,121.73,123.55,125.43,127.17,129.00,130.71,  0.00, 22.38/
      data (en(i),i=  9064,  9234)/
     * 30.63, 37.12, 42.51, 47.34, 51.65, 55.53, 59.23, 62.62, 65.88,
     * 68.95, 71.81, 74.56, 77.25, 79.91, 82.46, 84.91, 87.28, 89.69,
     * 91.88, 94.17, 96.32, 98.54,100.64,102.75,104.76,106.82,108.81,
     *110.77,112.69,114.61,116.57,118.43,120.23,122.06,123.96,125.75,
     *127.49,129.29,131.09,  0.00, 22.50, 30.80, 37.25, 42.65, 47.60,
     * 51.87, 55.84, 59.48, 62.91, 66.15, 69.21, 72.11, 74.89, 77.63,
     * 80.25, 82.76, 85.23, 87.64, 90.02, 92.30, 94.54, 96.73, 98.88,
     *101.02,103.07,105.11,107.15,109.11,111.06,113.05,114.94,116.81,
     *118.70,120.57,122.41,124.22,126.03,127.81,129.62,131.41,  0.00,
     * 22.61, 30.89, 37.42, 42.90, 47.71, 52.07, 56.03, 59.73, 63.14,
     * 66.37, 69.44, 72.36, 75.19, 77.88, 80.54, 83.07, 85.56, 87.98,
     * 90.29, 92.62, 94.87, 97.07, 99.21,101.32,103.38,105.46,107.48,
     *109.50,111.41,113.39,115.27,117.14,119.03,120.89,122.71,124.54,
     *126.30,128.18,129.92,131.64,  0.00, 22.72, 31.07, 37.61, 43.06,
     * 47.94, 52.21, 56.23, 59.96, 63.44, 66.61, 69.71, 72.62, 75.47,
     * 78.18, 80.77, 83.38, 85.88, 88.24, 90.59, 92.90, 95.18, 97.32,
     * 99.46,101.62,103.71,105.72,107.75,109.73,111.69,113.63,115.59,
     *117.44,119.32,121.18,123.01,124.86,126.68,128.40,130.16,132.01,
     *  0.00, 22.74, 31.16, 37.69, 43.13, 48.03, 52.43, 56.38, 60.08/
      data (en(i),i=  9235,  9405)/
     * 63.56, 66.80, 69.89, 72.81, 75.66, 78.33, 81.03, 83.64, 86.08,
     * 88.52, 90.85, 93.13, 95.37, 97.61, 99.75,101.89,103.95,105.97,
     *108.03,110.02,111.94,113.88,115.86,117.70,119.53,121.46,123.30,
     *125.12,126.88,128.71,130.46,132.22,  0.00, 22.82, 31.19, 37.75,
     * 43.27, 48.11, 52.54, 56.54, 60.24, 63.73, 66.94, 70.07, 73.03,
     * 75.87, 78.62, 81.31, 83.79, 86.32, 88.75, 91.05, 93.35, 95.64,
     * 97.83,100.00,102.10,104.18,106.25,108.29,110.22,112.24,114.17,
     *116.03,117.97,119.81,121.68,123.50,125.30,127.10,128.93,130.68,
     *132.37,  0.00, 22.85, 31.33, 37.91, 43.39, 48.27, 52.66, 56.69,
     * 60.40, 63.87, 67.15, 70.20, 73.21, 76.07, 78.82, 81.48, 84.06,
     * 86.48, 88.93, 91.31, 93.58, 95.83, 98.05,100.20,102.30,104.36,
     *106.44,108.45,110.48,112.45,114.38,116.31,118.19,120.04,121.86,
     *123.73,125.57,127.36,129.13,130.86,132.59,  0.00, 23.11, 31.33,
     * 37.97, 43.53, 48.35, 52.77, 56.80, 60.53, 63.99, 67.24, 70.41,
     * 73.35, 76.20, 78.96, 81.65, 84.20, 86.70, 89.09, 91.46, 93.76,
     * 96.03, 98.20,100.38,102.48,104.58,106.67,108.62,110.64,112.60,
     *114.53,116.43,118.36,120.25,122.09,123.93,125.77,127.51,129.30,
     *131.10,132.80,  0.00, 22.99, 31.45, 38.03, 43.58, 48.48, 52.88,
     * 56.91, 60.59, 64.09, 67.38, 70.48, 73.43, 76.32, 79.08, 81.74/
      data (en(i),i=  9406,  9576)/
     * 84.32, 86.82, 89.27, 91.60, 93.89, 96.13, 98.40,100.54,102.63,
     *104.71,106.79,108.83,110.77,112.79,114.68,116.63,118.53,120.38,
     *122.21,124.11,125.92,127.74,129.47,131.32,133.04,  0.00, 23.02,
     * 31.51, 38.05, 43.64, 48.50, 52.91, 56.96, 60.66, 64.16, 67.43,
     * 70.56, 73.58, 76.45, 79.20, 81.88, 84.44, 86.94, 89.33, 91.78,
     * 94.04, 96.26, 98.49,100.63,102.72,104.83,106.93,108.98,110.96,
     *112.85,114.83,116.77,118.66,120.54,122.40,124.19,126.06,127.85,
     *129.64,131.40,133.13,  0.00, 23.07, 31.51, 38.10, 43.65, 48.57,
     * 52.99, 56.97, 60.70, 64.19, 67.50, 70.66, 73.62, 76.51, 79.22,
     * 81.91, 84.49, 87.03, 89.46, 91.79, 94.12, 96.38, 98.55,100.73,
     *102.82,104.93,107.05,109.04,111.05,113.00,114.93,116.90,118.80,
     *120.64,122.53,124.34,126.13,128.00,129.72,131.49,133.33,  0.00,
     * 23.06, 31.57, 38.10, 43.72, 48.61, 53.05, 57.06, 60.81, 64.28,
     * 67.53, 70.71, 73.68, 76.61, 79.36, 82.00, 84.55, 87.05, 89.49,
     * 91.87, 94.21, 96.40, 98.62,100.80,102.93,104.99,107.09,109.13,
     *111.14,113.12,114.99,116.95,118.88,120.80,122.59,124.49,126.30,
     *128.09,129.83,131.67,133.37,  0.00, 23.09, 31.51, 38.15, 43.71,
     * 48.61, 53.01, 57.03, 60.74, 64.21, 67.57, 70.71, 73.71, 76.59,
     * 79.33, 82.03, 84.60, 87.07, 89.49, 91.84, 94.19, 96.47, 98.67/
      data (en(i),i=  9577,  9747)/
     *100.84,102.96,105.10,107.15,109.19,111.21,113.13,115.11,117.04,
     *118.89,120.81,122.66,124.51,126.34,128.13,129.97,131.70,133.50,
     *  0.00, 23.06, 31.55, 38.10, 43.67, 48.59, 52.97, 57.02, 60.75,
     * 64.19, 67.50, 70.66, 73.66, 76.57, 79.30, 82.00, 84.54, 87.10,
     * 89.49, 91.85, 94.19, 96.45, 98.66,100.83,103.02,105.10,107.13,
     *109.21,111.17,113.16,115.12,117.08,119.01,120.84,122.67,124.54,
     *126.42,128.25,130.00,131.76,133.51,  0.00, 23.06, 31.55, 38.14,
     * 43.65, 48.53, 52.95, 56.97, 60.64, 64.19, 67.51, 70.62, 73.61,
     * 76.53, 79.29, 81.93, 84.50, 87.00, 89.47, 91.84, 94.13, 96.40,
     * 98.61,100.80,102.98,105.05,107.16,109.16,111.17,113.14,115.11,
     *117.07,118.95,120.89,122.76,124.62,126.41,128.23,130.02,131.82,
     *133.64,  0.00, 23.04, 31.48, 38.06, 43.59, 48.51, 52.88, 56.87,
     * 60.61, 64.13, 67.40, 70.58, 73.55, 76.40, 79.22, 81.87, 84.47,
     * 86.94, 89.36, 91.76, 94.09, 96.35, 98.59,100.75,102.92,105.00,
     *107.06,109.12,111.18,113.12,115.13,117.06,119.01,120.83,122.79,
     *124.65,126.42,128.33,130.12,131.91,133.67,  0.00, 23.01, 31.46,
     * 37.95, 43.55, 48.43, 52.80, 56.75, 60.44, 64.01, 67.25, 70.44,
     * 73.46, 76.30, 79.07, 81.76, 84.34, 86.86, 89.25, 91.63, 93.94,
     * 96.21, 98.48,100.65,102.86,104.91,107.02,109.04,110.29,113.09/
      data (en(i),i=  9748,  9918)/
     *115.04,117.01,118.94,120.87,122.73,124.64,126.45,128.32,130.10,
     *131.94,133.69,  0.00, 22.99, 31.40, 37.98, 43.48, 48.31, 52.66,
     * 56.64, 60.37, 63.87, 67.21, 70.28, 73.31, 76.16, 78.90, 81.56,
     * 84.18, 86.66, 89.10, 91.46, 93.81, 96.13, 98.29,100.55,102.67,
     *104.79,106.91,108.96,111.04,113.06,115.03,116.99,118.94,120.85,
     *122.73,124.57,126.49,128.29,130.14,131.98,133.77,  0.00, 22.91,
     * 31.29, 37.81, 43.30, 48.12, 52.50, 56.49, 60.19, 63.68, 67.02,
     * 70.10, 73.07, 75.94, 78.70, 81.32, 83.99, 86.44, 88.92, 91.34,
     * 93.63, 95.92, 98.15,100.39,102.60,104.70,106.79,108.90,110.95,
     *112.96,114.94,116.90,118.84,120.77,122.67,124.57,126.49,128.36,
     *130.12,132.01,133.83,  0.00, 22.85, 31.24, 37.70, 43.14, 47.97,
     * 52.25, 56.29, 59.97, 63.45, 66.73, 69.88, 72.83, 75.72, 78.48,
     * 81.13, 83.71, 86.23, 88.68, 91.11, 93.44, 95.73, 98.06,100.17,
     *102.37,104.58,106.70,108.76,110.79,112.83,114.89,116.79,118.84,
     *120.81,122.67,124.57,126.51,128.32,130.13,132.01,133.91,  0.00,
     * 22.76, 31.04, 37.50, 42.95, 47.73, 52.02, 56.05, 59.68, 63.15,
     * 66.37, 69.54, 72.48, 75.39, 78.14, 80.83, 83.38, 85.99, 88.39,
     * 90.81, 93.19, 95.49, 97.74,100.00,102.31,104.40,106.47,108.64,
     *110.77,112.81,114.81,116.77,118.75,120.79,122.72,124.58,126.48/
      data (en(i),i=  9919, 10089)/
     *128.40,130.21,132.13,133.93,  0.00, 22.59, 30.89, 37.23, 42.63,
     * 47.43, 51.73, 55.68, 59.33, 62.76, 66.05, 69.13, 72.16, 74.99,
     * 77.79, 80.44, 83.06, 85.60, 88.10, 90.56, 92.97, 95.26, 97.54,
     * 99.88,102.05,104.22,106.45,108.60,110.67,112.74,114.80,116.83,
     *118.77,120.74,122.73,124.70,126.52,128.48,130.36,132.24,134.06,
     *  0.00, 22.47, 30.61, 36.95, 42.28, 47.05, 51.28, 55.17, 58.83,
     * 62.28, 65.55, 68.69, 71.64, 74.57, 77.32, 80.07, 82.74, 85.30,
     * 87.86, 90.30, 92.69, 95.13, 97.46, 99.73,101.96,104.19,106.39,
     *108.55,110.69,112.76,114.85,116.90,118.88,120.85,122.87,124.78,
     *126.78,128.69,130.55,132.38,134.22,  0.00, 22.27, 30.26, 36.46,
     * 41.72, 46.39, 50.64, 54.59, 58.24, 61.70, 64.93, 68.14, 71.25,
     * 74.12, 77.00, 79.75, 82.42, 85.09, 87.65, 90.14, 92.64, 95.04,
     * 97.40, 99.76,102.00,104.31,106.46,108.72,110.79,112.91,114.99,
     *117.06,119.14,121.14,123.07,125.11,126.98,128.98,130.86,132.83,
     *134.67,  0.00, 21.48, 29.34, 35.51, 40.87, 45.58, 49.93, 53.97,
     * 57.75, 61.29, 64.69, 67.92, 71.09, 74.08, 77.03, 79.85, 82.62,
     * 85.25, 87.86, 90.44, 92.95, 95.37, 97.76,100.09,102.39,104.65,
     *106.93,109.12,111.26,113.42,115.50,117.59,119.63,121.65,123.69,
     *125.61,127.59,129.54,131.43,133.37,135.23,  0.00, 19.71, 27.23/
      data (en(i),i= 10090, 10260)/
     * 33.17, 38.31, 42.79, 46.94, 50.77, 54.37, 57.78, 61.00, 64.14,
     * 67.15, 70.00, 72.77, 75.50, 78.11, 80.68, 83.13, 85.65, 88.00,
     * 90.34, 92.64, 94.92, 97.21, 99.44,101.59,103.76,105.90,107.96,
     *110.06,112.11,114.10,116.14,118.16,120.10,122.04,123.99,125.92,
     *127.75,129.64,  0.00, 20.46, 28.12, 34.09, 39.14, 43.63, 47.73,
     * 51.51, 55.04, 58.40, 61.58, 64.66, 67.56, 70.42, 73.20, 75.83,
     * 78.49, 81.02, 83.41, 85.84, 88.20, 90.58, 92.88, 95.20, 97.41,
     * 99.55,101.72,103.89,106.04,108.05,110.18,112.18,114.20,116.18,
     *118.19,120.10,122.00,123.88,125.75,127.68,129.59,  0.00, 20.90,
     * 20.52, 34.74, 39.91, 39.72, 48.52, 52.34, 55.85, 59.21, 62.35,
     * 65.40, 68.31, 71.11, 73.86, 76.50, 78.01, 81.61, 83.99, 85.69,
     * 87.75, 88.72, 93.35, 93.31, 95.58, 97.76, 99.99,102.08,106.37,
     *108.48,110.50,112.54,114.58,116.55,116.47,120.37,122.32,122.27,
     *126.14,128.01,129.86,  0.00, 21.20, 29.12, 35.29, 40.51, 45.08,
     * 49.27, 53.04, 56.59, 59.93, 63.16, 66.12, 69.08, 71.84, 74.54,
     * 77.19, 79.69, 82.20, 84.65, 86.98, 89.34, 91.71, 93.99, 96.22,
     * 98.36,100.52,102.77,104.79,106.87,108.97,111.08,112.96,114.94,
     *116.98,118.90,120.80,122.67,124.61,126.44,128.39,130.18,  0.00,
     * 21.75, 29.52, 35.76, 41.16, 45.67, 49.80, 53.72, 57.40, 60.60/
      data (en(i),i= 10261, 10431)/
     * 63.78, 66.84, 69.81, 72.49, 75.34, 77.91, 80.42, 82.92, 85.28,
     * 87.69, 90.11, 92.35, 94.59, 96.90, 98.95,101.19,103.19,105.45,
     *107.38,109.47,111.50,113.44,115.42,117.41,119.35,121.26,123.15,
     *124.99,126.85,128.72,130.57,  0.00, 21.81, 29.91, 36.16, 41.44,
     * 46.18, 50.38, 54.24, 57.85, 61.22, 64.44, 67.49, 70.42, 73.15,
     * 75.49, 78.47, 80.96, 83.45, 85.91, 88.28, 90.65, 92.91, 95.15,
     * 97.40, 99.53,101.73,103.82,105.86,107.96,109.94,112.04,113.93,
     *115.96,117.89,119.81,121.68,123.60,125.46,127.31,129.17,130.93,
     *  0.00, 22.00, 30.13, 36.51, 41.87, 46.56, 50.87, 54.73, 58.38,
     * 61.76, 64.98, 68.02, 70.99, 73.75, 76.42, 79.03, 81.56, 84.01,
     * 86.50, 88.86, 91.22, 93.53, 95.68, 97.94,100.14,102.22,104.36,
     *106.39,108.50,110.53,112.56,114.47,116.42,118.36,120.25,122.16,
     *124.06,125.88,127.68,129.57,131.31,  0.00, 22.19, 30.43, 36.78,
     * 42.22, 46.99, 51.36, 55.22, 58.85, 62.27, 65.50, 68.54, 71.46,
     * 74.29, 77.03, 79.59, 82.16, 84.63, 87.02, 89.43, 91.78, 94.03,
     * 96.29, 98.46,100.65,102.78,104.88,106.97,109.00,111.00,113.00,
     *114.94,116.89,118.81,120.68,122.56,124.49,126.30,128.24,129.94,
     *131.76,  0.00, 22.37, 30.65, 37.16, 42.56, 47.38, 51.71, 55.62,
     * 59.29, 62.71, 65.99, 69.12, 71.96, 74.80, 77.54, 80.14, 82.65/
      data (en(i),i= 10432, 10602)/
     * 85.10, 87.61, 89.92, 92.29, 94.59, 96.81, 99.00,101.17,103.31,
     *105.40,107.45,109.45,111.51,113.47,115.47,117.33,119.28,121.26,
     *123.11,124.97,126.78,128.59,130.44,132.20,  0.00, 22.53, 30.89,
     * 37.39, 42.86, 47.67, 52.09, 55.99, 59.69, 63.15, 66.39, 69.50,
     * 72.41, 75.25, 77.94, 80.61, 83.11, 85.60, 88.04, 90.45, 92.75,
     * 95.06, 97.29, 99.50,101.63,103.73,105.82,107.88,109.96,111.96,
     *113.91,115.83,117.79,119.72,121.63,123.43,125.37,127.22,129.02,
     *130.77,132.59,  0.00, 22.71, 31.03, 37.59, 43.12, 47.97, 52.30,
     * 56.33, 60.04, 63.51, 66.73, 69.90, 72.79, 75.67, 78.39, 80.98,
     * 83.61, 86.05, 88.53, 90.89, 93.28, 95.49, 97.78, 99.93,102.08,
     *104.21,106.33,108.37,110.37,112.40,114.38,116.31,118.25,120.15,
     *121.98,123.94,125.73,127.57,129.46,131.24,133.04,  0.00, 22.79,
     * 31.17, 37.75, 43.36, 48.17, 52.57, 56.64, 60.31, 63.85, 67.10,
     * 70.27, 73.23, 76.05, 78.79, 81.41, 84.02, 86.48, 88.92, 91.33,
     * 93.66, 95.91, 98.19,100.41,102.56,104.63,106.76,108.83,110.85,
     *112.79,114.74,116.75,118.68,120.55,122.40,124.30,126.13,127.99,
     *129.75,131.59,133.39,  0.00, 22.92, 31.51, 37.98, 43.57, 48.44,
     * 52.85, 56.89, 60.62, 64.16, 67.46, 70.56, 73.55, 76.45, 79.16,
     * 81.77, 84.40, 86.92, 89.35, 91.74, 94.07, 96.35, 98.62,100.76/
      data (en(i),i= 10603, 10773)/
     *102.97,105.05,107.17,109.26,111.20,113.21,115.18,117.05,119.00,
     *120.97,122.86,124.68,126.56,128.39,130.11,131.97,133.77,  0.00,
     * 23.03, 31.53, 38.15, 43.77, 48.67, 53.12, 57.19, 60.92, 64.43,
     * 67.77, 70.88, 73.90, 76.70, 79.52, 82.15, 84.76, 87.26, 89.73,
     * 92.08, 94.48, 96.74, 98.93,101.19,103.29,105.42,107.51,109.55,
     *111.58,113.65,115.54,117.49,119.40,121.31,123.15,125.03,126.91,
     *128.76,130.54,132.34,134.09,  0.00, 23.19, 31.70, 38.32, 43.92,
     * 49.17, 53.33, 57.36, 61.21, 64.69, 68.04, 71.20, 74.21, 77.04,
     * 79.79, 82.43, 85.11, 87.61, 90.04, 92.42, 94.78, 97.07, 99.35,
     *101.54,103.70,105.77,107.89,109.92,111.96,113.95,115.92,117.79,
     *119.80,121.65,123.58,125.37,127.25,129.08,130.84,132.66,134.42,
     *  0.00, 23.26, 31.82, 38.45, 44.12, 49.05, 53.56, 57.58, 61.39,
     * 64.96, 68.28, 71.42, 74.47, 77.31, 80.09, 82.80, 85.40, 87.94,
     * 90.40, 92.73, 95.19, 97.40, 99.66,101.88,104.02,106.13,108.23,
     *110.29,112.28,114.28,116.25,118.19,120.10,121.96,123.79,125.77,
     *127.53,129.29,131.17,132.96,134.73,  0.00, 23.35, 31.91, 38.60,
     * 44.26, 49.24, 53.67, 57.82, 61.61, 65.13, 68.49, 71.66, 74.67,
     * 77.57, 80.30, 83.05, 85.64, 88.23, 90.64, 93.07, 95.42, 97.73,
     * 99.92,102.13,104.29,106.40,108.54,110.59,112.55,114.56,116.56/
      data (en(i),i= 10774, 10944)/
     *118.46,120.38,122.28,124.11,126.03,127.81,129.63,131.39,133.20,
     *134.99,  0.00, 23.43, 32.01, 38.73, 44.37, 49.39, 53.87, 57.95,
     * 61.76, 65.34, 68.71, 71.89, 74.89, 77.79, 80.55, 83.35, 85.90,
     * 88.43, 90.95, 93.36, 95.69, 97.95,100.18,102.43,104.60,106.70,
     *108.79,110.84,112.85,114.85,116.77,118.77,120.63,122.58,124.46,
     *126.31,128.09,129.90,131.75,133.47,135.28,  0.00, 23.52, 32.14,
     * 38.78, 44.50, 49.48, 53.99, 58.13, 61.92, 65.50, 68.87, 72.08,
     * 75.08, 78.00, 80.81, 83.52, 86.14, 88.70, 91.13, 93.57, 95.94,
     * 98.20,100.42,102.65,104.84,106.93,109.02,111.13,113.06,115.09,
     *117.06,118.95,120.89,122.75,124.73,126.48,128.36,130.17,131.97,
     *133.79,135.53,  0.00, 23.53, 32.21, 38.93, 44.59, 49.63, 54.11,
     * 58.28, 62.09, 65.69, 69.00, 72.19, 75.25, 78.17, 80.99, 83.68,
     * 86.34, 88.89, 91.35, 93.75, 96.15, 98.43,100.69,102.90,105.04,
     *107.17,109.22,111.31,113.33,115.28,117.23,119.19,121.14,123.02,
     *124.87,126.78,128.57,130.36,132.16,133.96,135.75,  0.00, 23.58,
     * 32.21, 39.01, 44.67, 49.69, 54.25, 58.36, 62.20, 65.79, 69.15,
     * 72.34, 75.41, 78.31, 81.15, 83.88, 86.53, 89.07, 91.53, 93.94,
     * 96.33, 98.57,100.89,103.04,105.23,107.33,109.46,111.52,113.49,
     *115.49,117.43,119.40,121.33,123.24,125.10,126.94,128.75,130.62/
      data (en(i),i= 10945, 11115)/
     *132.40,134.17,135.94,  0.00, 23.67, 32.36, 39.12, 44.72, 49.78,
     * 54.30, 58.45, 62.32, 65.92, 69.30, 72.45, 75.51, 78.48, 81.24,
     * 84.00, 86.68, 89.22, 91.68, 94.11, 96.45, 98.79,100.99,103.25,
     *105.41,107.51,109.53,111.64,113.72,115.58,117.64,119.55,121.54,
     *123.42,125.27,127.12,128.98,130.80,132.53,134.34,136.13,  0.00,
     * 23.72, 32.35, 39.16, 44.83, 49.84, 54.40, 58.55, 62.40, 65.97,
     * 69.35, 72.52, 75.61, 78.58, 81.41, 84.15, 86.77, 89.31, 91.85,
     * 94.27, 96.60, 98.89,101.17,103.39,105.54,107.66,109.77,111.83,
     *113.79,115.80,117.77,119.74,121.63,123.54,125.39,127.30,129.15,
     *130.95,132.77,134.49,136.26,  0.00, 23.75, 32.41, 39.21, 44.89,
     * 49.89, 54.49, 58.63, 62.45, 66.01, 69.44, 72.61, 75.69, 78.67,
     * 81.48, 84.19, 86.87, 89.47, 91.87, 94.36, 96.70, 99.00,101.24,
     *103.49,105.63,107.76,109.82,111.95,113.93,115.97,117.90,119.82,
     *121.81,123.67,125.55,127.35,129.26,131.08,132.88,134.62,136.41,
     *  0.00, 23.80, 32.43, 39.20, 44.94, 49.92, 54.49, 58.63, 62.52,
     * 66.08, 69.48, 72.65, 75.74, 78.73, 81.58, 84.30, 86.93, 89.47,
     * 92.01, 94.40, 96.79, 99.10,101.34,103.57,105.74,107.87,109.95,
     *112.01,114.03,116.00,118.03,119.93,121.90,123.79,125.68,127.50,
     *129.39,131.23,132.96,134.77,136.55,  0.00, 23.78, 32.44, 39.24/
      data (en(i),i= 11116, 11286)/
     * 44.93, 49.98, 54.52, 58.67, 62.50, 66.10, 69.49, 72.70, 75.78,
     * 78.75, 81.61, 84.36, 87.00, 89.57, 92.05, 94.49, 96.85, 99.14,
     *101.42,103.60,105.77,107.94,109.88,112.04,114.06,116.08,118.12,
     *120.05,122.02,123.84,125.77,127.60,129.45,131.30,133.04,134.87,
     *136.67,  0.00, 23.80, 32.47, 39.29, 44.98, 49.99, 54.51, 58.68,
     * 62.50, 66.14, 69.53, 72.73, 75.81, 78.79, 81.62, 84.33, 86.99,
     * 89.56, 92.05, 94.52, 96.84, 99.16,101.42,103.66,105.80,107.95,
     *110.03,112.12,114.14,116.16,118.14,120.06,121.97,123.95,125.80,
     *127.68,129.53,131.35,133.12,134.97,136.74,  0.00, 23.84, 32.46,
     * 39.24, 44.96, 49.97, 54.53, 58.64, 62.47, 66.11, 69.49, 72.74,
     * 75.82, 78.76, 81.57, 84.32, 87.00, 89.62, 92.09, 94.49, 96.87,
     * 99.15,101.45,103.62,105.83,107.92,110.06,112.13,114.17,116.14,
     *118.15,120.13,122.05,123.97,125.82,127.73,129.58,131.42,133.23,
     *134.98,136.80,  0.00, 23.85, 32.45, 39.24, 44.93, 49.95, 54.47,
     * 58.64, 62.48, 66.02, 69.41, 72.71, 75.77, 78.71, 81.56, 84.28,
     * 86.96, 89.50, 92.03, 94.45, 96.84, 99.11,101.41,103.64,105.82,
     *107.93,110.04,112.09,114.12,116.14,118.12,120.14,122.04,123.99,
     *125.88,127.73,129.55,131.43,133.23,135.07,136.83,  0.00, 23.78,
     * 32.43, 39.16, 44.87, 49.88, 54.39, 58.59, 62.35, 65.96, 69.33/
      data (en(i),i= 11287, 11457)/
     * 72.61, 75.69, 78.63, 81.51, 84.22, 86.91, 89.43, 91.96, 94.39,
     * 96.76, 99.07,101.34,103.55,105.72,107.92,109.96,112.07,114.13,
     *116.15,118.12,120.05,122.02,123.88,125.84,127.74,129.54,131.43,
     *133.21,135.02,136.82,  0.00, 23.80, 32.38, 39.17, 44.85, 49.80,
     * 54.33, 58.49, 62.30, 65.89, 69.34, 72.54, 75.59, 78.55, 81.39,
     * 84.11, 86.72, 89.36, 91.82, 94.31, 96.62, 98.94,101.26,103.52,
     *105.67,107.82,109.93,112.00,114.00,116.09,118.03,119.99,121.97,
     *123.87,125.82,127.66,129.61,131.41,133.23,135.07,136.89,  0.00,
     * 23.73, 32.37, 39.06, 44.77, 49.76, 54.22, 58.34, 62.19, 65.78,
     * 69.17, 72.43, 75.45, 78.40, 81.26, 83.99, 86.67, 89.20, 91.71,
     * 94.13, 96.53, 98.84,101.13,103.32,105.54,107.71,109.76,111.92,
     *113.96,116.01,117.96,119.94,121.92,123.86,125.77,127.68,129.54,
     *131.41,133.25,135.03,136.85,  0.00, 23.72, 32.30, 38.94, 44.63,
     * 49.64, 54.05, 58.19, 62.03, 65.62, 68.99, 72.23, 75.29, 78.20,
     * 81.08, 83.78, 86.45, 89.00, 91.51, 93.95, 96.31, 98.66,100.92,
     *103.19,105.35,107.52,109.67,111.72,113.85,115.88,117.89,119.84,
     *121.81,123.81,125.72,127.56,129.48,131.32,133.18,135.02,136.89,
     *  0.00, 23.66, 32.22, 38.87, 44.52, 49.47, 53.90, 58.03, 61.84,
     * 65.42, 68.83, 71.99, 75.08, 77.99, 80.87, 83.53, 86.19, 88.81/
      data (en(i),i= 11458, 11628)/
     * 91.29, 93.76, 96.06, 98.43,100.75,103.01,105.20,107.34,109.47,
     *111.53,113.68,115.72,117.75,119.76,121.69,123.72,125.62,127.54,
     *129.47,131.34,133.16,135.08,136.84,  0.00, 23.54, 32.06, 38.74,
     * 44.33, 49.28, 53.73, 57.79, 61.58, 65.19, 68.56, 71.73, 74.75,
     * 77.75, 80.54, 83.32, 85.97, 88.50, 91.04, 93.48, 95.85, 98.23,
     *100.53,102.75,105.02,107.17,109.35,111.44,113.53,115.63,117.67,
     *119.67,121.67,123.60,125.62,127.50,129.37,131.26,133.18,135.01,
     *136.83,  0.00, 23.44, 31.90, 38.50, 44.08, 49.01, 53.42, 57.48,
     * 61.31, 64.88, 68.20, 71.44, 74.45, 77.41, 80.24, 82.95, 85.63,
     * 88.22, 90.70, 93.17, 95.57, 97.95,100.25,102.51,104.75,106.98,
     *109.18,111.29,113.39,115.52,117.52,119.54,121.63,123.54,125.54,
     *127.44,129.38,131.32,133.21,135.07,136.89,  0.00, 23.36, 31.70,
     * 38.31, 43.78, 48.69, 53.08, 57.17, 60.94, 64.49, 67.85, 70.96,
     * 74.04, 77.01, 79.81, 82.56, 85.26, 87.87, 90.39, 92.85, 95.27,
     * 97.70,100.03,102.32,104.59,106.78,108.98,111.13,113.27,115.38,
     *117.50,119.52,121.57,123.53,125.55,127.52,129.38,131.37,133.22,
     *135.07,137.00,  0.00, 23.16, 31.52, 37.94, 43.48, 48.29, 52.70,
     * 56.68, 60.46, 63.96, 67.29, 70.50, 73.59, 76.51, 79.36, 82.16,
     * 84.82, 87.47, 90.05, 92.55, 95.03, 97.44, 99.80,102.12,104.46/
      data (en(i),i= 11629, 11799)/
     *106.67,108.89,111.10,113.23,115.37,117.45,119.52,121.64,123.58,
     *125.59,127.56,129.51,131.48,133.45,135.30,137.16,  0.00, 22.91,
     * 31.07, 37.34, 42.79, 47.69, 51.96, 56.07, 59.73, 63.26, 66.69,
     * 69.85, 73.00, 75.97, 78.90, 81.77, 84.48, 87.19, 89.86, 92.35,
     * 94.85, 97.35, 99.73,102.05,104.41,106.73,109.02,111.17,113.38,
     *115.43,117.57,119.65,121.71,123.78,125.85,127.83,129.74,131.77,
     *133.68,135.56,137.35,  0.00, 22.15, 30.18, 36.56, 41.97, 46.86,
     * 51.31, 55.35, 59.22, 62.88, 66.37, 69.67, 72.85, 75.95, 78.92,
     * 81.77, 84.57, 87.30, 90.02, 92.56, 95.08, 97.56,100.01,102.42,
     *104.72,107.06,109.29,111.59,113.79,115.95,118.10,120.19,122.28,
     *124.37,126.32,128.40,130.31,132.31,134.25,136.18,138.03,  0.00,
     * 20.24, 27.94, 33.99, 39.20, 43.80, 48.01, 51.94, 55.61, 59.07,
     * 62.39, 65.53, 68.56, 71.52, 74.32, 77.06, 79.75, 82.37, 84.96,
     * 87.37, 89.83, 92.22, 94.52, 96.79, 99.09,101.31,103.53,105.75,
     *107.85,110.02,112.07,114.20,116.25,118.22,120.22,122.26,124.14,
     *126.11,128.08,129.93,131.88,  0.00, 20.99, 28.73, 34.86, 40.05,
     * 44.69, 48.75, 52.82, 56.29, 59.75, 63.92, 66.13, 69.11, 72.04,
     * 74.76, 77.51, 80.14, 82.77, 85.24, 87.56, 90.06, 92.43, 94.77,
     * 97.02, 99.33,101.47,103.69,105.89,107.96,110.14,112.24,114.29/
      data (en(i),i= 11800, 11970)/
     *116.33,118.32,120.26,122.32,124.23,126.15,128.02,129.93,131.87,
     *  0.00, 21.37, 29.33, 35.57, 40.78, 45.00, 49.67, 53.53, 57.08,
     * 60.53, 63.76, 66.89, 69.87, 72.72, 75.53, 78.22, 80.82, 83.31,
     * 85.84, 88.28, 90.70, 92.97, 95.24, 97.59, 99.78,102.00,104.15,
     *106.36,108.41,110.55,112.59,114.62,116.66,118.67,120.66,122.63,
     *124.60,126.47,128.33,130.22,132.15,  0.00, 21.74, 29.84, 36.09,
     * 41.44, 46.14, 50.40, 54.28, 57.90, 61.33, 64.57, 67.67, 70.56,
     * 73.48, 76.23, 78.81, 81.47, 84.06, 86.50, 88.91, 91.26, 93.63,
     * 95.86, 98.11,100.35,102.53,104.75,106.82,108.97,111.03,113.15,
     *115.15,117.18,119.13,121.13,123.04,124.92,126.89,128.83,130.64,
     *132.48,  0.00, 22.06, 30.23, 36.59, 41.97, 46.72, 50.98, 54.93,
     * 58.60, 62.04, 65.22, 68.32, 71.29, 74.20, 76.91, 79.60, 82.12,
     * 84.67, 87.12, 89.51, 91.86, 94.20, 96.55, 98.74,100.93,103.17,
     *105.38,107.48,109.55,111.62,113.60,115.60,117.68,119.60,121.58,
     *123.45,125.47,127.22,129.21,131.03,132.88,  0.00, 22.30, 30.60,
     * 36.99, 42.40, 47.18, 51.51, 55.45, 59.16, 62.60, 65.85, 68.99,
     * 71.98, 74.83, 77.53, 80.21, 82.79, 85.28, 87.73, 90.11, 92.49,
     * 94.83, 97.06, 99.41,101.56,103.78,105.84,108.02,110.06,112.16,
     *114.12,116.20,118.14,120.06,122.09,123.92,125.90,127.76,129.60/
      data (en(i),i= 11971, 12141)/
     *131.48,133.35,  0.00, 22.57, 30.83, 37.49, 42.80, 47.65, 52.07,
     * 56.08, 59.74, 63.25, 66.49, 69.55, 72.54, 75.46, 78.10, 80.86,
     * 83.41, 85.96, 88.34, 90.78, 93.12, 95.42, 97.69, 99.99,102.18,
     *104.32,106.50,108.61,110.64,112.69,114.67,116.67,118.69,120.66,
     *122.62,124.51,126.38,128.25,130.17,132.03,133.75,  0.00, 22.74,
     * 31.11, 37.65, 43.12, 48.06, 52.38, 56.43, 60.16, 63.67, 66.97,
     * 70.11, 73.05, 75.97, 78.67, 81.41, 83.98, 86.49, 88.95, 91.33,
     * 93.65, 96.03, 98.27,100.54,102.76,104.86,107.00,109.12,111.20,
     *113.16,115.25,117.22,119.20,121.09,123.02,124.99,126.82,128.75,
     *130.62,132.37,134.24,  0.00, 22.92, 31.33, 37.90, 43.49, 48.36,
     * 52.78, 56.82, 60.57, 64.13, 67.42, 70.57, 73.57, 76.47, 79.19,
     * 81.88, 84.49, 87.00, 89.45, 91.83, 94.22, 96.53, 98.84,101.09,
     *103.27,105.37,107.55,109.60,111.67,113.70,115.70,117.73,119.65,
     *121.65,123.52,125.44,127.35,129.20,131.09,132.90,134.61,  0.00,
     * 23.11, 31.58, 38.20, 43.76, 48.71, 53.15, 57.23, 60.99, 64.52,
     * 67.86, 71.06, 74.06, 76.92, 79.71, 82.42, 85.00, 87.49, 89.94,
     * 92.38, 94.79, 97.08, 99.38,101.55,103.76,105.91,108.05,110.11,
     *112.13,114.20,116.23,118.23,120.14,122.07,124.02,125.85,127.77,
     *129.68,131.49,133.30,135.09,  0.00, 23.24, 31.84, 38.49, 44.13/
      data (en(i),i= 12142, 12312)/
     * 48.99, 53.46, 57.57, 61.36, 64.93, 68.31, 71.45, 74.47, 77.41,
     * 80.17, 82.87, 85.43, 87.94, 90.45, 92.86, 95.24, 97.61, 99.86,
     *102.13,104.23,106.41,108.50,110.56,112.70,114.72,116.65,118.65,
     *120.60,122.54,124.45,126.35,128.22,130.05,131.91,133.77,135.53,
     *  0.00, 23.36, 31.96, 38.61, 44.29, 49.28, 53.75, 57.87, 61.72,
     * 65.25, 68.64, 71.80, 74.89, 77.76, 80.57, 83.27, 85.87, 88.39,
     * 90.90, 93.31, 95.67, 98.01,100.30,102.49,104.68,106.82,109.00,
     *111.09,113.08,115.17,117.16,119.18,121.08,122.97,124.90,126.76,
     *128.64,130.45,132.36,134.16,135.89,  0.00, 23.49, 32.17, 38.83,
     * 44.45, 49.51, 54.00, 58.18, 61.98, 65.57, 68.99, 72.15, 75.23,
     * 78.14, 80.93, 83.67, 86.31, 88.80, 91.31, 93.76, 96.09, 98.47,
     *100.76,102.93,105.12,107.29,109.43,111.46,113.51,115.58,117.59,
     *119.49,121.45,123.40,125.29,127.22,129.07,130.94,132.70,134.57,
     *136.36,  0.00, 23.64, 32.12, 39.10, 44.74, 49.73, 54.24, 58.45,
     * 62.30, 65.88, 69.30, 72.50, 75.55, 78.52, 81.34, 84.11, 86.69,
     * 89.17, 91.74, 94.12, 96.53, 98.86,101.13,103.32,105.55,107.73,
     *109.83,111.91,113.96,116.02,117.90,119.93,121.85,123.74,125.76,
     *127.57,129.45,131.30,133.14,134.87,136.72,  0.00, 23.79, 32.40,
     * 39.15, 44.92, 49.90, 54.51, 58.71, 62.58, 66.17, 69.58, 72.78/
      data (en(i),i= 12313, 12483)/
     * 75.90, 78.86, 81.66, 84.34, 86.98, 89.55, 92.08, 94.52, 96.89,
     * 99.20,101.45,103.70,105.93,108.06,110.17,112.23,114.33,116.31,
     *118.30,120.26,122.20,124.11,126.05,127.93,129.80,131.61,133.42,
     *135.26,137.09,  0.00, 23.86, 32.57, 39.37, 45.06, 50.13, 54.72,
     * 58.88, 62.80, 66.43, 69.81, 73.08, 76.17, 79.14, 81.96, 84.65,
     * 87.23, 89.87, 92.35, 94.84, 97.24, 99.54,101.84,104.05,106.29,
     *108.41,110.54,112.60,114.64,116.74,118.68,120.63,122.57,124.54,
     *126.41,128.33,130.20,132.02,133.83,135.64,137.37,  0.00, 23.91,
     * 32.63, 39.46, 45.21, 50.29, 54.99, 59.10, 62.96, 66.62, 70.10,
     * 73.35, 76.42, 79.36, 82.23, 84.94, 87.60, 90.17, 92.66, 95.10,
     * 97.57, 99.90,102.18,104.38,106.53,108.72,110.88,112.91,114.97,
     *117.06,118.97,120.92,122.94,124.84,126.72,128.58,130.45,132.31,
     *134.13,135.91,137.67,  0.00, 24.05, 33.08, 39.61, 45.41, 50.45,
     * 55.06, 59.26, 63.22, 66.86, 70.26, 73.54, 76.63, 79.55, 82.43,
     * 85.21, 87.87, 90.44, 92.98, 95.40, 97.84,100.15,102.44,104.63,
     *106.86,109.02,111.15,113.28,115.29,117.35,119.27,121.29,123.16,
     *125.15,127.07,128.88,130.75,132.63,134.41,136.20,138.03,  0.00,
     * 24.08, 32.91, 39.74, 45.53, 50.56, 55.21, 59.40, 63.35, 67.00,
     * 70.49, 73.74, 76.84, 79.80, 82.64, 85.43, 88.14, 90.70, 93.23/
      data (en(i),i= 12484, 12654)/
     * 95.72, 98.13,100.41,102.70,104.97,107.15,109.36,111.46,113.52,
     *115.60,117.63,119.57,121.52,123.50,125.40,127.27,129.20,131.03,
     *132.87,134.70,136.54,138.30,  0.00, 24.21, 32.98, 39.86, 45.67,
     * 50.74, 55.39, 59.58, 63.47, 67.18, 70.65, 73.91, 77.06, 80.02,
     * 82.90, 85.61, 88.36, 90.94, 93.49, 95.94, 98.32,100.70,102.98,
     *105.24,107.41,109.55,111.69,113.74,115.84,117.84,119.87,121.78,
     *123.74,125.70,127.61,129.41,131.27,133.17,134.99,136.72,138.53,
     *  0.00, 24.25, 33.00, 39.90, 45.74, 50.80, 55.49, 59.71, 63.66,
     * 67.31, 70.81, 74.07, 77.20, 80.20, 83.07, 85.83, 88.56, 91.17,
     * 93.68, 96.16, 98.59,100.91,103.19,105.40,107.64,109.78,111.95,
     *114.01,116.05,118.04,120.04,122.02,124.01,125.92,127.81,129.71,
     *131.49,133.39,135.17,136.95,138.76,  0.00, 24.42, 33.10, 40.00,
     * 45.86, 51.01, 55.63, 59.87, 63.83, 67.48, 71.00, 74.30, 77.41,
     * 80.32, 83.23, 86.05, 88.73, 91.36, 93.87, 96.31, 98.74,101.10,
     *103.42,105.65,107.85,110.09,112.15,114.23,116.26,118.36,120.30,
     *122.29,124.28,126.12,128.00,129.88,131.77,133.57,135.42,137.26,
     *138.98,  0.00, 24.26, 33.18, 40.05, 45.96, 51.01, 55.70, 59.95,
     * 63.91, 67.58, 71.07, 74.35, 77.48, 80.50, 83.37, 86.14, 88.88,
     * 91.52, 94.04, 96.49, 98.93,101.26,103.57,105.85,108.05,110.21/
      data (en(i),i= 12655, 12825)/
     *112.30,114.41,116.42,118.47,120.65,122.41,124.42,126.29,128.19,
     *130.10,131.94,133.81,135.55,137.45,139.17,  0.00, 24.46, 33.22,
     * 40.15, 45.98, 51.12, 55.81, 60.08, 64.01, 67.67, 71.18, 74.49,
     * 77.62, 80.60, 83.49, 86.30, 89.04, 91.61, 94.17, 96.66, 99.08,
     *101.39,103.74,106.00,108.19,110.36,112.49,114.56,116.59,118.66,
     *120.66,122.64,124.54,126.49,128.37,130.27,132.12,133.97,135.74,
     *137.54,139.39,  0.00, 24.50, 33.32, 40.22, 46.05, 51.25, 55.87,
     * 60.16, 64.04, 67.75, 71.26, 74.54, 77.68, 80.69, 83.64, 86.35,
     * 89.11, 91.73, 94.27, 96.77, 99.19,101.52,103.85,106.05,108.31,
     *110.48,112.60,114.68,116.73,118.80,120.75,122.77,124.73,126.63,
     *128.53,130.41,132.22,134.05,135.94,137.71,139.51,  0.00, 24.51,
     * 33.28, 40.26, 46.06, 51.22, 55.92, 60.14, 64.10, 67.81, 71.27,
     * 74.59, 77.77, 80.77, 83.64, 86.51, 89.18, 91.81, 94.40, 96.84,
     * 99.27,101.66,103.96,106.21,108.43,110.55,112.69,114.80,116.91,
     *118.88,120.91,122.86,124.81,126.74,128.63,130.53,132.39,134.18,
     *135.99,137.80,139.58,  0.00, 24.55, 33.43, 40.30, 46.16, 51.30,
     * 55.97, 60.25, 64.13, 67.87, 71.35, 74.66, 77.78, 80.82, 83.77,
     * 86.55, 89.25, 91.88, 94.44, 96.93, 99.32,101.73,103.98,106.29,
     *108.52,110.63,112.80,114.89,116.94,118.98,120.95,122.93,124.93/
      data (en(i),i= 12826, 12996)/
     *126.82,128.70,130.65,132.50,134.33,136.11,137.92,139.78,  0.00,
     * 24.51, 33.39, 40.33, 46.14, 51.31, 55.96, 60.20, 64.19, 67.90,
     * 71.42, 74.68, 77.82, 80.85, 83.78, 86.59, 89.30, 91.91, 94.46,
     * 96.93, 99.38,101.75,104.10,106.33,108.50,110.72,112.87,114.90,
     *117.00,119.04,121.01,123.03,124.99,126.92,128.81,130.63,132.58,
     *134.39,136.24,137.99,139.77,  0.00, 24.54, 33.42, 40.31, 46.20,
     * 51.32, 55.99, 60.20, 64.21, 67.87, 71.37, 74.67, 77.85, 80.90,
     * 83.77, 86.57, 89.36, 91.92, 94.54, 97.01, 99.43,101.75,104.07,
     *106.35,108.58,110.70,112.90,114.98,117.02,119.08,121.03,123.03,
     *125.01,126.93,128.86,130.76,132.64,134.45,136.33,138.07,139.87,
     *  0.00, 24.50, 33.36, 40.33, 46.19, 51.22, 55.97, 60.25, 64.21,
     * 67.87, 71.37, 74.65, 77.81, 80.85, 83.81, 86.57, 89.28, 91.92,
     * 94.47, 97.00, 99.43,101.80,104.04,106.36,108.59,110.75,112.85,
     *114.92,117.06,119.23,121.03,123.05,125.00,126.93,128.87,130.78,
     *132.58,134.51,136.27,138.11,139.96,  0.00, 24.55, 33.37, 40.30,
     * 46.11, 51.30, 55.94, 60.18, 64.12, 67.84, 71.30, 74.63, 77.80,
     * 80.86, 83.76, 86.52, 89.28, 91.90, 94.47, 96.95, 99.36,101.70,
     *104.03,106.31,108.51,110.67,112.82,114.95,117.02,119.04,121.04,
     *123.08,125.05,126.94,128.85,130.74,132.62,134.48,136.36,138.13/
      data (en(i),i= 12997, 13167)/
     *140.00,  0.00, 24.41, 33.35, 40.22, 46.07, 51.28, 55.86, 60.18,
     * 64.05, 67.79, 71.24, 74.57, 77.71, 80.76, 83.71, 86.48, 89.22,
     * 91.85, 94.37, 96.88, 99.27,101.63,103.97,106.28,108.47,110.64,
     *112.82,114.90,117.01,118.98,121.04,123.06,125.02,126.97,128.84,
     *130.72,132.63,134.53,136.31,138.11,139.95,  0.00, 24.68, 33.33,
     * 40.27, 46.04, 51.20, 55.83, 60.07, 64.01, 67.69, 71.16, 74.51,
     * 77.63, 80.68, 83.58, 86.42, 89.12, 91.74, 94.33, 96.80, 99.20,
     *101.57,103.87,106.13,108.37,110.59,112.67,114.85,116.91,118.98,
     *120.95,122.96,124.95,126.86,128.81,130.70,132.57,134.48,136.29,
     *138.19,139.92,  0.00, 24.45, 33.29, 40.17, 45.99, 51.07, 55.75,
     * 59.92, 63.87, 67.56, 71.01, 74.39, 77.53, 80.58, 83.45, 86.26,
     * 88.99, 91.63, 94.14, 96.67, 99.07,101.38,103.74,106.02,108.25,
     *110.44,112.59,114.74,116.79,118.88,120.91,122.92,124.88,126.80,
     *128.73,130.65,132.57,134.36,136.28,138.09,139.97,  0.00, 24.44,
     * 33.19, 40.08, 45.86, 50.98, 55.59, 59.84, 63.71, 67.41, 70.87,
     * 74.20, 77.36, 80.41, 83.29, 86.08, 88.78, 91.44, 94.00, 96.48,
     * 98.92,101.29,103.61,105.89,108.10,110.29,112.43,114.58,116.66,
     *118.71,120.77,122.74,124.75,126.69,128.68,130.60,132.46,134.39,
     *136.22,138.07,139.86,  0.00, 24.40, 33.22, 39.97, 45.83, 50.85/
      data (en(i),i= 13168, 13338)/
     * 55.47, 59.65, 63.62, 67.29, 70.74, 74.01, 77.21, 80.24, 83.10,
     * 85.90, 88.58, 91.24, 93.81, 96.28, 98.75,101.08,103.39,105.67,
     *107.96,110.16,112.30,114.39,116.54,118.60,120.62,122.68,124.67,
     *126.65,128.58,130.48,132.44,134.26,136.19,138.07,139.87,  0.00,
     * 24.32, 33.03, 39.84, 45.58, 50.63, 55.29, 59.49, 63.38, 67.04,
     * 70.50, 73.83, 76.97, 79.93, 82.83, 85.66, 88.37, 90.95, 93.53,
     * 96.06, 98.45,100.84,103.21,105.47,107.68,109.96,112.14,114.24,
     *116.32,118.45,120.48,122.52,124.49,126.48,128.49,130.42,132.36,
     *134.22,136.12,137.91,139.77,  0.00, 24.23, 32.91, 39.72, 45.42,
     * 50.46, 55.01, 59.22, 63.12, 66.78, 70.26, 73.53, 76.64, 79.69,
     * 82.54, 85.37, 88.03, 90.67, 93.21, 95.75, 98.16,100.55,102.90,
     *105.21,107.49,109.71,111.85,114.04,116.12,118.25,120.33,122.35,
     *124.36,126.41,128.37,130.27,132.26,134.12,136.10,137.89,139.74,
     *  0.00, 24.10, 32.73, 39.51, 45.21, 50.26, 54.77, 58.99, 62.86,
     * 66.50, 69.92, 73.17, 76.30, 79.33, 82.17, 84.98, 87.75, 90.35,
     * 92.94, 95.42, 97.91,100.26,102.66,104.95,107.22,109.51,111.67,
     *113.81,116.00,118.10,120.16,122.21,124.24,126.25,128.28,130.28,
     *132.22,134.10,136.04,137.97,139.76,  0.00, 23.96, 32.56, 39.27,
     * 44.92, 49.90, 54.43, 58.56, 62.49, 66.06, 69.50, 72.79, 75.87/
      data (en(i),i= 13339, 13509)/
     * 78.85, 81.79, 84.57, 87.30, 90.00, 92.53, 95.10, 97.58, 99.98,
     *102.34,104.67,106.97,109.24,111.48,113.69,115.80,117.95,120.09,
     *122.19,124.23,126.20,128.25,130.17,132.16,134.06,136.08,138.00,
     *139.88,  0.00, 23.81, 32.36, 38.90, 44.54, 49.52, 53.98, 58.11,
     * 61.93, 65.56, 68.98, 72.28, 75.41, 78.39, 81.28, 84.11, 86.87,
     * 89.53, 92.13, 94.70, 97.21, 99.75,102.06,104.46,106.78,109.03,
     *111.38,113.58,115.73,117.91,119.99,122.08,124.20,126.28,128.27,
     *130.31,132.27,134.18,136.15,138.08,140.01,  0.00, 23.54, 31.93,
     * 38.48, 44.06, 48.96, 53.36, 57.49, 61.33, 64.90, 68.33, 71.62,
     * 74.81, 77.85, 80.83, 83.67, 86.48, 89.21, 91.91, 94.46, 97.01,
     * 99.54,101.95,104.36,106.72,109.06,111.37,113.55,115.82,117.98,
     *120.09,122.21,124.36,126.42,128.44,130.42,132.48,134.43,136.36,
     *138.35,140.22,  0.00, 22.84, 31.04, 37.54, 43.05, 48.06, 52.55,
     * 56.75, 60.69, 64.42, 67.97, 71.33, 74.56, 77.72, 80.76, 83.69,
     * 86.51, 89.27, 91.97, 94.62, 97.19, 99.75,102.23,104.65,107.02,
     *109.33,111.69,113.93,116.19,118.38,120.56,122.67,124.78,126.84,
     *128.89,130.94,132.91,134.93,136.97,138.84,140.84,  0.00, 20.76,
     * 28.61, 34.79, 40.09, 44.76, 49.11, 53.09, 56.88, 60.35, 63.66,
     * 66.90, 69.95, 72.97, 75.88, 78.60, 81.39, 84.03, 86.62, 89.11/
      data (en(i),i= 13510, 13680)/
     * 91.56, 93.94, 96.30, 98.66,100.91,103.20,105.41,107.57,109.79,
     *111.91,114.04,116.15,118.24,120.31,122.33,124.30,126.28,128.23,
     *130.13,132.12,133.97,  0.00, 21.51, 29.42, 35.68, 40.95, 45.65,
     * 49.87, 53.85, 57.50, 61.02, 64.34, 67.49, 70.51, 73.52, 76.30,
     * 79.06, 81.67, 84.36, 86.97, 88.36, 91.57, 94.24, 96.64, 98.86,
     *101.10,103.38,105.56,107.81,109.93,112.14,114.18,116.23,118.36,
     *120.39,122.39,124.38,126.32,128.29,130.25,132.15,133.99,  0.00,
     * 13.22, 30.05, 36.37, 36.45, 41.79, 50.75, 54.65, 58.31, 58.40,
     * 65.16, 68.29, 71.31, 74.24, 74.29, 77.07, 82.41, 82.46, 87.59,
     * 90.00, 92.41, 94.82, 96.64, 99.38,101.58,101.68,103.98,108.32,
     *110.44,112.48,114.71,116.75,116.75,120.75,122.71,124.68,126.70,
     *128.57,130.51,130.56,134.35,  0.00, 22.19, 21.19, 29.68, 36.25,
     * 47.07, 46.63, 55.40, 55.00, 62.56, 65.89, 69.04, 71.99, 75.05,
     * 77.80, 80.57, 83.14, 82.90, 85.48, 90.74, 93.11, 92.82, 97.76,
     *100.01, 99.80,102.05,106.65,106.48,108.61,113.14,115.17,117.16,
     *119.24,121.24,123.26,125.17,124.98,126.92,130.90,130.70,132.66,
     *  0.00, 22.71, 30.91, 37.42, 42.91, 47.68, 52.01, 56.05, 59.78,
     * 63.27, 66.72, 69.73, 72.79, 75.67, 78.55, 81.19, 83.85, 86.42,
     * 88.89, 91.31, 93.69, 96.04, 98.35,100.70,102.91,105.07,107.31/
      data (en(i),i= 13681, 13851)/
     *109.43,111.59,113.66,115.67,117.71,119.78,121.78,123.75,125.71,
     *127.65,129.52,131.44,133.25,135.16,  0.00, 22.77, 31.18, 37.73,
     * 43.30, 48.04, 52.53, 56.61, 60.38, 63.87, 67.22, 70.40, 73.45,
     * 76.32, 79.11, 81.85, 84.56, 87.04, 89.63, 92.01, 94.35, 96.68,
     * 98.92,101.28,103.48,105.72,107.86,109.99,112.12,114.25,116.30,
     *118.29,120.28,122.23,124.22,126.22,128.15,129.98,131.89,133.71,
     *135.62,  0.00, 23.00, 31.48, 38.11, 43.63, 48.56, 53.04, 57.12,
     * 60.86, 64.47, 67.78, 71.04, 74.07, 76.97, 79.76, 82.47, 85.16,
     * 87.66, 90.21, 92.61, 94.91, 97.30, 99.56,101.90,104.14,106.32,
     *108.50,110.63,112.67,114.79,116.78,118.87,120.80,122.83,124.79,
     *126.70,128.62,130.59,132.40,134.25,136.11,  0.00, 23.19, 31.74,
     * 38.42, 44.03, 49.02, 53.47, 57.55, 61.37, 64.93, 68.33, 71.56,
     * 74.64, 77.52, 80.35, 83.08, 85.73, 88.31, 90.78, 93.17, 95.54,
     * 97.86,100.23,102.45,104.66,106.90,108.99,111.18,113.28,115.30,
     *117.39,119.34,121.33,123.31,125.31,127.22,129.12,131.05,132.88,
     *134.77,136.64,  0.00, 23.38, 32.03, 38.70, 44.49, 49.33, 53.83,
     * 58.00, 61.89, 65.48, 68.83, 72.08, 75.12, 78.06, 80.88, 83.64,
     * 86.26, 88.85, 91.28, 93.80, 96.17, 98.50,100.83,102.99,105.27,
     *107.42,109.58,111.73,113.86,115.87,117.92,119.94,121.94,123.83/
      data (en(i),i= 13852, 14022)/
     *125.89,127.67,129.59,131.49,133.35,135.19,137.01,  0.00, 23.58,
     * 32.24, 38.96, 44.69, 49.66, 54.25, 58.40, 62.22, 65.86, 69.23,
     * 72.48, 75.61, 78.52, 81.35, 84.13, 86.74, 89.32, 91.85, 94.21,
     * 96.62, 98.98,101.32,103.56,105.81,107.99,110.13,112.22,114.37,
     *116.41,118.44,120.46,122.38,124.37,126.32,128.21,130.17,132.04,
     *133.88,135.65,137.55,  0.00, 23.70, 32.43, 39.24, 44.95, 49.95,
     * 54.52, 58.72, 62.61, 66.25, 69.63, 72.94, 76.03, 78.99, 81.85,
     * 84.57, 87.27, 89.80, 92.29, 94.66, 97.19, 99.50,101.80,104.11,
     *106.27,108.50,110.64,112.70,114.77,116.88,118.87,120.90,122.86,
     *124.82,126.80,128.71,130.56,132.44,134.30,136.10,137.93,  0.00,
     * 23.86, 32.66, 39.46, 45.17, 50.27, 54.83, 59.07, 62.96, 66.62,
     * 69.99, 73.28, 76.42, 79.43, 82.27, 85.02, 87.65, 90.26, 92.80,
     * 95.17, 97.62, 99.98,102.28,104.53,106.76,108.95,111.08,113.22,
     *115.28,117.34,119.42,121.33,123.33,125.30,127.25,129.18,131.03,
     *132.88,134.72,136.56,138.43,  0.00, 24.03, 32.81, 39.65, 45.38,
     * 50.51, 55.11, 59.34, 63.22, 66.97, 70.42, 73.66, 76.76, 79.77,
     * 82.62, 85.37, 88.08, 90.67, 93.21, 95.60, 98.06,100.42,102.73,
     *104.98,107.21,109.39,111.55,113.64,115.78,117.79,119.78,121.81,
     *123.79,125.71,127.68,129.52,131.46,133.34,135.22,137.02,138.85/
      data (en(i),i= 14023, 14193)/
     *  0.00, 24.15, 32.96, 39.82, 45.64, 50.77, 55.37, 59.67, 63.53,
     * 67.25, 70.75, 74.04, 77.15, 80.17, 83.04, 85.82, 88.49, 91.01,
     * 93.60, 96.07, 98.47,100.85,103.16,105.42,107.63,109.86,111.96,
     *114.08,116.15,118.24,120.23,122.27,124.24,126.18,128.06,129.97,
     *131.90,133.77,135.61,137.37,139.24,  0.00, 24.37, 33.09, 40.02,
     * 45.87, 50.99, 55.63, 59.93, 63.85, 67.58, 70.99, 74.38, 77.47,
     * 80.47, 83.34, 86.17, 88.81, 91.42, 93.93, 96.40, 98.89,101.22,
     *103.55,105.83,108.06,110.22,112.36,114.51,116.60,118.61,120.66,
     *122.66,124.64,126.56,128.46,130.38,132.24,134.14,135.99,137.76,
     *139.60,  0.00, 24.40, 33.25, 40.17, 45.99, 51.14, 55.85, 60.07,
     * 64.04, 67.82, 71.28, 74.60, 77.78, 80.77, 83.67, 86.48, 89.18,
     * 91.75, 94.30, 96.75, 99.24,101.56,103.93,106.18,108.44,110.61,
     *112.79,114.88,117.01,118.99,121.02,123.02,124.98,126.92,128.86,
     *130.72,132.61,134.50,136.32,138.12,139.96,  0.00, 24.33, 33.38,
     * 40.29, 46.08, 51.33, 56.00, 60.28, 64.27, 68.05, 71.56, 74.86,
     * 78.03, 80.99, 84.00, 86.75, 89.47, 92.09, 94.59, 97.08, 99.53,
     *101.92,104.22,106.49,108.73,110.90,113.14,115.23,117.24,119.37,
     *121.39,123.37,125.34,127.32,129.15,131.12,132.99,134.80,136.73,
     *138.51,140.34,  0.00, 24.68, 33.48, 40.46, 46.31, 51.56, 56.22/
      data (en(i),i= 14194, 14364)/
     * 60.55, 64.46, 68.24, 71.81, 75.11, 78.31, 81.38, 84.27, 87.07,
     * 89.72, 92.35, 94.93, 97.45, 99.85,102.23,104.59,106.85,109.06,
     *111.32,113.47,115.57,117.62,119.68,121.68,123.70,125.65,127.63,
     *129.57,131.46,133.29,135.17,137.04,138.88,140.68,  0.00, 24.68,
     * 33.56, 40.63, 46.50, 51.69, 56.38, 60.68, 64.73, 68.46, 71.95,
     * 75.33, 78.57, 81.54, 84.51, 87.32, 89.97, 92.66, 95.19, 97.70,
     *100.14,102.53,104.88,107.15,109.41,111.62,113.74,115.89,117.97,
     *119.99,122.05,124.01,126.00,127.96,129.85,131.73,133.59,135.51,
     *137.34,139.15,140.96,  0.00, 22.28, 31.99, 39.25, 46.61, 51.84,
     * 56.51, 60.90, 64.86, 68.66, 72.18, 75.51, 78.07, 81.78, 84.70,
     * 87.53, 90.24, 92.87, 95.48, 97.94, 99.89,102.29,105.20,106.96,
     *109.18,111.90,114.03,116.14,118.22,120.25,121.82,124.32,126.25,
     *127.76,130.16,132.06,133.93,135.35,137.16,139.07,141.21,  0.00,
     * 24.81, 33.76, 40.76, 46.74, 51.98, 56.69, 61.04, 65.04, 68.81,
     * 72.32, 75.73, 78.90, 81.97, 84.89, 87.69, 90.42, 93.13, 95.72,
     * 98.21,100.70,103.05,105.39,107.70,109.94,112.11,114.27,116.44,
     *118.51,120.49,122.57,124.60,126.53,128.50,130.38,132.27,134.18,
     *136.03,137.83,139.71,141.48,  0.00, 24.87, 33.83, 40.87, 46.85,
     * 52.08, 56.85, 61.14, 65.19, 68.94, 72.49, 75.84, 79.12, 82.13/
      data (en(i),i= 14365, 14535)/
     * 85.09, 87.94, 90.65, 93.32, 95.94, 98.41,100.91,103.28,105.63,
     *107.93,110.18,112.39,114.50,116.61,118.74,120.78,122.81,124.77,
     *126.80,128.73,130.62,132.56,134.40,136.30,138.12,139.94,141.77,
     *  0.00, 24.94, 34.24, 40.99, 46.94, 52.18, 56.99, 61.26, 65.34,
     * 69.08, 72.67, 76.03, 79.23, 82.27, 85.25, 88.04, 90.79, 93.47,
     * 96.11, 98.60,101.08,103.46,105.79,108.09,110.35,112.57,114.77,
     *116.84,118.94,121.01,123.05,125.01,127.06,128.96,130.85,132.74,
     *134.59,136.52,138.33,140.16,141.93,  0.00, 25.05, 34.00, 41.07,
     * 47.03, 52.26, 57.03, 61.35, 65.48, 69.18, 72.73, 76.12, 79.38,
     * 82.40, 85.38, 88.22, 90.99, 93.67, 96.26, 98.77,101.24,103.67,
     *105.98,108.29,110.54,112.77,114.90,117.06,119.13,121.18,123.21,
     *125.17,127.18,129.12,131.09,132.99,134.83,136.68,138.57,140.33,
     *142.19,  0.00, 25.05, 34.05, 41.08, 47.17, 52.36, 57.10, 61.52,
     * 65.54, 69.29, 72.86, 76.23, 79.48, 82.60, 85.51, 88.32, 91.10,
     * 93.82, 96.43, 98.92,101.42,103.81,106.15,108.50,110.75,112.95,
     *115.09,117.23,119.30,121.36,123.42,125.41,127.39,129.35,131.19,
     *133.16,134.97,136.85,138.74,140.52,142.36,  0.00, 25.06, 34.18,
     * 41.18, 47.15, 52.40, 57.13, 61.54, 65.59, 69.42, 73.01, 76.37,
     * 79.57, 82.65, 85.61, 88.43, 91.26, 93.93, 96.53, 99.09,101.52/
      data (en(i),i= 14536, 14706)/
     *103.91,106.34,108.62,110.84,113.08,115.26,117.35,119.50,121.50,
     *123.61,125.53,127.52,129.44,131.40,133.33,135.10,137.04,138.94,
     *140.71,142.46,  0.00, 25.16, 34.16, 41.24, 47.20, 52.46, 57.22,
     * 61.64, 65.64, 69.43, 73.04, 76.41, 79.63, 82.74, 85.84, 88.56,
     * 91.33, 94.02, 96.65, 99.17,101.65,104.04,106.42,108.73,110.97,
     *113.19,115.36,117.49,119.60,121.66,123.69,125.78,127.63,129.60,
     *131.54,133.41,135.35,137.14,139.02,140.86,142.59,  0.00, 25.15,
     * 34.15, 41.23, 47.23, 52.56, 57.26, 61.65, 65.68, 69.55, 73.09,
     * 76.52, 79.71, 82.81, 85.72, 88.66, 91.44, 94.31, 96.68, 99.25,
     *101.77,104.12,106.50,108.82,111.07,113.25,115.45,117.57,119.67,
     *121.73,123.74,125.74,127.78,129.72,131.58,133.55,135.40,137.26,
     *139.15,140.97,142.76,  0.00, 25.13, 34.18, 41.27, 47.25, 52.56,
     * 57.28, 61.72, 65.77, 69.59, 73.14, 76.53, 79.73, 82.81, 85.79,
     * 88.69, 91.49, 94.14, 96.77, 99.33,101.80,104.24,106.56,108.89,
     *111.16,113.36,115.50,117.64,119.78,121.83,123.84,125.84,127.85,
     *129.81,131.72,133.66,135.48,137.36,139.19,141.03,142.87,  0.00,
     * 25.20, 34.20, 41.28, 47.26, 52.61, 57.33, 61.75, 65.80, 69.56,
     * 73.18, 76.53, 79.76, 82.89, 85.82, 88.69, 91.54, 94.22, 96.86,
     * 99.40,101.84,104.28,106.63,108.90,111.17,113.35,115.59,117.74/
      data (en(i),i= 14707, 14877)/
     *119.85,121.87,123.92,125.93,127.87,129.83,131.76,133.71,135.61,
     *137.48,139.28,141.16,142.91,  0.00, 25.17, 34.22, 41.28, 47.29,
     * 52.56, 57.36, 61.71, 65.77, 69.58, 73.17, 76.56, 79.73, 82.90,
     * 85.87, 88.74, 91.54, 94.20, 96.80, 99.39,101.88,104.28,106.65,
     *108.96,111.23,113.44,115.62,117.74,119.85,121.90,123.93,125.98,
     *127.95,129.84,131.86,133.74,135.64,137.48,139.34,141.13,142.97,
     *  0.00, 25.19, 34.30, 41.31, 47.32, 52.62, 57.36, 61.74, 65.75,
     * 69.55, 73.11, 76.55, 79.76, 82.84, 85.88, 88.74, 91.54, 94.19,
     * 96.84, 99.41,101.86,104.29,106.60,108.93,111.28,113.43,115.61,
     *117.70,119.82,121.95,123.93,125.99,127.98,129.94,131.86,133.73,
     *135.67,137.53,139.38,141.17,142.96,  0.00, 25.25, 34.20, 41.33,
     * 47.28, 52.57, 57.29, 61.68, 65.76, 69.51, 73.07, 76.46, 79.73,
     * 82.86, 85.80, 88.72, 91.53, 94.16, 96.81, 99.35,101.84,104.23,
     *106.60,108.93,111.19,113.37,115.54,117.70,119.83,121.89,123.91,
     *125.96,127.90,129.93,131.73,133.71,135.64,137.54,139.38,141.16,
     *143.00,  0.00, 25.20, 34.24, 41.29, 47.20, 52.51, 57.30, 61.62,
     * 65.69, 69.51, 73.01, 76.40, 79.71, 82.76, 85.80, 88.63, 91.44,
     * 94.12, 96.73, 99.23,101.77,104.21,106.52,108.87,111.07,113.34,
     *115.49,117.64,119.76,121.83,123.87,125.91,127.92,129.83,131.77/
      data (en(i),i= 14878, 15048)/
     *133.70,135.64,137.48,139.35,141.20,143.04,  0.00, 25.05, 34.16,
     * 41.20, 47.29, 52.47, 57.20, 61.51, 65.64, 69.40, 72.91, 76.32,
     * 79.57, 82.70, 85.69, 88.54, 91.39, 94.06, 96.63, 99.19,101.68,
     *104.06,106.42,108.72,111.07,113.23,115.46,117.54,119.68,121.77,
     *123.85,125.85,127.82,129.77,131.73,133.61,135.58,137.46,139.34,
     *141.13,143.00,  0.00, 25.17, 34.18, 41.21, 47.15, 52.42, 57.11,
     * 61.50, 65.53, 69.27, 72.79, 76.23, 79.52, 82.61, 85.59, 88.44,
     * 91.23, 93.95, 96.51, 99.11,101.55,103.99,106.35,108.64,110.88,
     *113.13,115.35,117.47,119.61,121.64,123.75,125.74,127.74,129.69,
     *131.63,133.59,135.51,137.37,139.24,141.16,142.91,  0.00, 25.10,
     * 34.06, 41.01, 47.02, 52.27, 56.94, 61.34, 65.41, 69.13, 72.73,
     * 76.15, 79.31, 82.47, 85.44, 88.34, 91.05, 93.74, 96.40, 98.94,
     *101.40,103.82,106.14,108.49,110.74,113.02,115.20,117.35,119.47,
     *121.54,123.62,125.62,127.63,129.65,131.58,133.54,135.45,137.31,
     *139.21,141.08,142.93,  0.00, 25.09, 34.06, 41.04, 46.92, 52.20,
     * 56.94, 61.19, 65.28, 69.00, 72.57, 75.95, 79.20, 82.31, 85.26,
     * 88.10, 90.93, 93.61, 96.19, 98.75,101.23,103.61,105.98,108.31,
     *110.60,112.82,115.03,117.16,119.27,121.43,123.40,125.47,127.46,
     *129.48,131.43,133.40,135.31,137.25,139.11,141.00,142.87,  0.00/
      data (en(i),i= 15049, 15219)/
     * 25.01, 33.92, 40.93, 46.80, 52.07, 56.76, 61.05, 65.06, 68.83,
     * 72.41, 75.73, 78.97, 82.08, 85.07, 87.87, 90.68, 93.36, 95.94,
     * 98.51,101.01,103.43,105.76,108.07,110.38,112.58,114.79,117.00,
     *119.10,121.20,123.30,125.33,127.38,129.39,131.34,133.33,135.18,
     *137.11,139.08,140.87,142.78,  0.00, 24.96, 33.84, 40.80, 46.68,
     * 51.87, 56.59, 60.83, 64.81, 68.57, 72.17, 75.52, 78.70, 81.80,
     * 84.79, 87.63, 90.38, 93.10, 95.71, 98.23,100.71,103.14,105.51,
     *107.82,110.15,112.39,114.55,116.76,118.94,121.05,123.15,125.18,
     *127.19,129.23,131.18,133.22,135.09,136.94,138.97,140.80,142.67,
     *  0.00, 24.92, 33.73, 40.67, 46.47, 51.63, 56.35, 60.59, 64.55,
     * 68.35, 71.84, 75.22, 78.46, 81.52, 84.50, 87.32, 90.05, 92.79,
     * 95.37, 97.97,100.41,102.84,105.28,107.61,109.88,112.09,114.35,
     *116.56,118.72,120.83,122.89,124.95,127.01,129.06,131.02,133.05,
     *135.00,136.94,138.83,140.76,142.56,  0.00, 24.74, 33.56, 40.48,
     * 46.25, 51.40, 56.03, 60.29, 64.27, 67.99, 71.57, 74.86, 78.08,
     * 81.17, 84.12, 86.95, 89.69, 92.43, 95.03, 97.58,100.07,102.52,
     *104.93,107.31,109.62,111.86,114.13,116.24,118.50,120.60,122.77,
     *124.79,126.85,128.91,130.91,132.90,134.87,136.88,138.76,140.69,
     *142.66,  0.00, 24.62, 33.36, 40.16, 45.98, 51.05, 55.68, 59.97/
      data (en(i),i= 15220, 15390)/
     * 63.87, 67.64, 71.10, 74.45, 77.67, 80.69, 83.69, 86.53, 89.30,
     * 92.02, 94.64, 97.20, 99.74,102.21,104.56,106.94,109.30,111.56,
     *113.83,116.03,118.27,120.40,122.61,124.71,126.72,128.82,130.88,
     *132.88,134.88,136.79,138.78,140.74,142.56,  0.00, 24.39, 33.13,
     * 39.89, 45.65, 50.71, 55.29, 59.50, 63.38, 67.13, 70.70, 73.93,
     * 77.13, 80.18, 83.10, 86.03, 88.86, 91.53, 94.24, 96.79, 99.33,
     *101.85,104.31,106.64,109.03,111.36,113.69,115.88,118.14,120.34,
     *122.45,124.64,126.74,128.76,130.87,132.90,134.88,136.85,138.79,
     *140.76,142.73,  0.00, 24.16, 32.73, 39.42, 45.08, 50.10, 54.63,
     * 58.86, 62.72, 66.40, 69.90, 73.24, 76.53, 79.63, 82.58, 85.52,
     * 88.39, 91.15, 93.87, 96.48, 99.06,101.58,104.14,106.57,108.91,
     *111.29,113.63,115.90,118.12,120.36,122.50,124.66,126.81,128.96,
     *131.00,133.02,135.01,137.06,139.04,141.01,142.95,  0.00, 23.49,
     * 31.81, 38.47, 44.14, 49.16, 53.81, 58.06, 62.07, 65.89, 69.45,
     * 72.94, 76.24, 79.38, 82.50, 85.48, 88.40, 91.16, 93.93, 96.59,
     * 99.25,101.79,104.33,106.81,109.21,111.59,113.92,116.16,118.46,
     *120.66,122.89,125.05,127.17,129.34,131.43,133.49,135.46,137.54,
     *139.47,141.50,143.46,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/
      do 1 i=1,16000
1     ent(i)=en(i)
      return
      end
 
c----------------------------------------------------------------------
      function sjcen4(k,ku,u)
c     returns total entropy.
c     constant volume 15 fm**3.
c     input: total quark number k; up and down quarks ku;
c     total excitation energy u[gev].
c     initial call of sjcen0 required!!!
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/csjcen/ent(16000)	
 
      sjcen4 = 0
 
      if(k.lt.ku)then
           if(ish.ge.90)then
      call utmsg('sjcen4')
      write(ifch,*)'*****  k < ku'
      write(ifch,*)'*****  k: ',k,'   ku: ',ku
      call utmsgf
           endif
      return
      endif
 
      if(k.gt.45)then
      ko=k
      k=45
           if(ish.ge.90)then
      call utmsg('sjcen4')
      write(ifch,*)'*****  k > 45'
      write(ifch,*)'*****  k: ',ko,'   k_new: ',k
      call utmsgf
           endif
      endif
 
      if(mod(k,3).ne.0)
     *call utstop('sjcen4: noninteger baryon number        ')
 
      ik=0
      do 1 i=3,k-3,3
1     ik=ik+i+1       		
      if(u.lt.10)then
      du=mod(u,.25)
      iu=int(u/.25)+1
      ientry=(ik+ku)*41+iu
      sjcen4=ent(ientry)+(ent(ientry+1)-ent(ientry))*du/.25
      else
      iu=40
      ientry=(ik+ku)*41+iu
      sjcen4=ent(ientry)+(ent(ientry+1)-ent(ientry))/.25*(u-9.75)
      endif
 
      return
      end
 
c----------------------------------------------------------------------
      function sjcent(k,ku,u)
c     returns entropy.
c     input: quark number k; energy u [gev].
c     iopent=1: oscillator model.
c       interpolates and extrapolates entro(1+k,1+n)
c       from sr jcentr (jcentd).
c       sr jcentr (jcentd) has to be called before!!
c     iopent=2,3: fermi gas model; joerg  aichelin.
c       iopent=2: const volume, iopent=3: const density
c     iopent=4: fermi gas (new)
c     iopent=5: hagedorn
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (kparx=15,nquax=12)
      common/centro/entro(1+kparx,1+nquax)
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      iop=mod(iopent,10)
 
      sjcent=0
 
           if(iop.eq.0)then
 
      return
 
           elseif(iop.eq.1)then
 
      if(k.le.kentro)return
      if(mod(k,3).ne.0)
     *call utstop('sjcent: k must be multiple of 3         ')
      e=u/omega
      if(k.gt.kparx)goto5001
      n=int(e)
      if(n.lt.0)
     *call utstop('sjcent: negative energy                 ')
      if(n.ge.nquax)n=nquax-1
      sjcent=entro(1+k,1+n)+(e-n)*(entro(1+k,1+n+1)-entro(1+k,1+n))
      if(k.le.3)return
5002  e3=e*3./k
      n3=int(e3)
      if(n3.ge.nquax)n3=nquax-1
      sentr3=entro(1+3,1+n3)+(e3-n3)*(entro(1+3,1+n3+1)-entro(1+3,1+n3))
      sjcent=amin1(sjcent,k/3.*sentr3)
      return
5001  continue
      l=kparx
      el=(e*l)/k
      nl=int(el)
      if(nl.ge.nquax)nl=nquax-1
      sentrl=entro(1+l,1+nl)+(el-nl)*(entro(1+l,1+nl+1)-entro(1+l,1+nl))
      sjcent=(k*1.)/l*sentrl
      goto5002
 
           elseif(iop.eq.2)then
 
c     constant volume 5 fm**3
      if(k.le.kentro)return
      sjcent = 0
      x1=12.96*k**(-.315)
      if(k.lt.19)then
      x2=.785 + .005*k
      else
      x2=.88
      endif
      sjcent=k*x1*(u/k)**x2
      return
 
           elseif(iop.eq.3)then
 
c     constant density
      if(k.le.kentro)return
      sjcent = 0
      x1=9.785
      x2=.7926
      sjcent=k*x1*(u/k)**x2
      return
 
           elseif(iop.eq.4)then
 
c     constant volume 15 fm**3
      if(k.le.kentro)return
      sjcent=sjcen4(k,ku,u)
      return

           elseif(iop.eq.5)then
 
c     hagedorn
      if(u.le.uentro.and.k.le.kentro)return
      thaged=.250
      sjcent = u/thaged   
      return
 
           else
 
      call utstop('sjcent: invalid option iop              ')
 
           endif
 
      end
 
c----------------------------------------------------------------------
      function sjcgam(keux,kedx,kesx,kecx,ama,amo,po,mox)
c     returns partial decay width dgamma = phase space * density
c       for decay of cluster into cluster and hadron.
c     ke*x: net quark number
c     ama: cluster mass;  amo: hadron mass;  po: hadron momentum
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/cnsta/pi,prom,piom,ainfin
      common/csjcga/amnull,nsuha(7),asuha(7),entrpy,amegam,nopha
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      common/centex/entexp
      sjcgam=0.
      entrpy=-99.9
      amegam=0.
      if(mox.eq.1)entexp=0.
      entrpy=entrpy-entexp
      ke=iabs(keux+kedx+kesx+kecx)
      ku=iabs(keux+kedx)
      if(mox.eq.1)amnull=utamnu(keux,kedx,kesx,kecx,iopenu)
      eo=sqrt(amo**2+po**2)
      if(eo.gt.ama)return
      ame2=(ama-eo)**2-po**2
      if(ame2.lt.0.)return
      amegam=sqrt(ame2)
      e=amegam-amnull
      if(e.lt.0.)return
      entrpy=sjcent(ke,ku,e)
      if(mox.eq.1)entexp=entrpy
      entrpy=entrpy-entexp
      sjcgam=1./8./pi**2/ama*
     *exp(entrpy)*
     *po**2/sqrt(amo**2+po**2)
      return
      end
 
c-----------------------------------------------------------------------
      function sjctem(k,ku,u)
c     returns temperature. input: quark number k; u,d qurks ku; energy u
c-----------------------------------------------------------------------
      common/cnsta/pi,prom,piom,ainfin
      du=0.001
      beta=(sjcent(k,ku,u+du)-sjcent(k,ku,u))/du
      if(beta.lt.1/ainfin)then
      sjctem=ainfin
      else
      sjctem=1/beta
      endif
      return
      end
 
c-----------------------------------------------------------------------
      function smass(a,y,z)
c     returns cluster mass (in gev) (per cluster, not (!) per nucleon)
c     according to berger/jaffe mass formula, prc35(1987)213 eq.2.31,
c     see also c. dover, BNL-46322, intersections-meeting, tucson, 91.
c     a: massnr, y: hypercharge, z: charge,
c-----------------------------------------------------------------------
      common/cmass/thet,epsi,as,ac,dy,dz,ym,cz,zm,sigma,rzero
      ymin=ym*a
      zmin=cz/(dz/a+zm/a**(1./3.))
      smass=epsi*a+as*a**(2./3.)+(ac/a**(1./3.)+dz/a/2.)*(z-zmin)**2
     *+dy/a/2.*(y-ymin)**2
      return
      end
 
c-----------------------------------------------------------------------
      subroutine smassi(theta)
c     initialization for smass.
c     calculates parameters for berger/jaffe mass formula
c     (prc35(1987)213 eq.2.31, see also c. dover, BNL-46322).
c     theta: parameter that determines all parameters in mass formula.
c-----------------------------------------------------------------------
      common/cmass/thet,epsi,as,ac,dy,dz,ym,cz,zm,sigma,rzero
      thet=theta
 
      astr=.150
      pi=3.14159
      alp=1./137.
 
      co=cos(theta)
      si=sin(theta)
      bet=(1+co**3)/2.
      rzero=si/astr/(  2./3./pi*(1+co**3)  )**(1./3.)
      cs=astr/si
      cz=-astr/si*(  (  .5*(1+co**3)  )**(1./3.)-1  )
      sigma=6./8./pi*(astr/si)**3*(co**2/6.-si**2*(1-si)/3.-
     *1./3./pi*(pi/2.-theta-sin(2*theta)+si**3*alog((1+co)/si)))
 
      epsi=astr*((.5*(1+co**3))**(1./3.)+2)/si
      as=4*pi*sigma*rzero**2
      ac=3./5.*alp/rzero
      dz=astr/si*bet**(1./3.)*co**2*
     *(co**4*(1+bet**(2./3.))+(1+bet)**2)/
     *(  (2*co**2+bet**(1./3.))*(co**4*(1+bet**(2./3.))+(1+bet)**2)-
     *(co**4+bet**(1./3.)*(1+bet))*((2*bet**(2./3.)-1)*co**2+1+bet)  )
      dy=astr/6.*(1+co**3)**3/si*
     *(  1+(1+co)/(4*(1+co**3))**(2./3.)  )/
     *(co**6+co+co*(.5*(1+co**3))**(4./3.))
      zm=6*alp/(5*rzero)
      ym=(1-co**3)/(1+co**3)
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine smassp
c     prints smass.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/cmass/thet,epsi,as,ac,dy,dz,ym,cz,zm,sigma,rzero
      real eng(14),ymi(14),zmi(14)
      pi=3.14159
      write(ifch,*)'parameters of mass formula:'
      write(ifch,*)'---------------------------'
      write(ifch,*)' '
      write(ifch,*)'theta=',thet,'   epsi=',epsi
      write(ifch,*)'as=',as,'   ac=',ac
      write(ifch,*)'dy=',dy,'   dz=',dz
      write(ifch,*)'ym=',ym
      write(ifch,*)'cz dz zm=',cz,dz,zm
      write(ifch,*)'sigma**1/3=',sigma**(1./3.),'   rzero=',rzero
      write(ifch,*)' '
      write(ifch,*)'mass:'
      write(ifch,*)'-----'
      write(ifch,5000)(j,j=1,14)
5000  format(/5x,'a:',14i5/)
      do 4 j=1,14
      a=j
      ymi(j)=ym*a
4     zmi(j)=cz/(dz/a+zm/a**(1./3.))
      write(ifch,5002)(ymi(j),j=1,14)
5002  format(1x,'ymin: ',14f5.2/)
      write(ifch,5003)(zmi(j),j=1,14)
5003  format(1x,'zmin: ',14f5.2/)
      do 2 i=1,15
      ns=11-i
      do 3 j=1,14
      a=j
      y=a-ns
      z=0.
3     eng(j)=smass(a,y,z)/a
      write(ifch,5001)ns,(eng(j),j=1,14)
5001  format(1x,'s=',i2,2x,14f5.2)
2     continue
      write(ifch,*)' '
      write(ifch,*)'mass-mass(free):'
      write(ifch,*)'----------------'
      write(ifch,5000)(j,j=1,14)
      do 5 i=1,15
      ns=11-i
      do 6 j=1,14
      a=j
      y=a-ns
      z=0.
      call smassu(a,y,z,ku,kd,ks,kc)
6     eng(j)=(smass(a,y,z)-utamnu(ku,kd,ks,kc,3))/a
      write(ifch,5001)ns,(eng(j),j=1,14)
5     continue
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine smasst(kux,kdx,ksx,kcx,a,y,z)
c     input: kux,kdx,ksx,kcx = net quark numbers (for u,d,s,c quarks).
c     output: massnr a, hypercharge y and charge z.
c-----------------------------------------------------------------------
      sg=1
      if(kux+kdx+ksx+kcx.lt.0.)sg=-1
      ku=sg*kux
      kd=sg*kdx
      ks=sg*ksx
      kc=sg*kcx
      k=ku+kd+ks+kc
      if(mod(k,3).ne.0)stop'noninteger baryon number'
      a=k/3
      y=a-ks
      nz=2*ku-kd-ks+2*kc
      if(mod(nz,3).ne.0)stop'noninteger charge'
      z=nz/3
      return
      end
 
c-----------------------------------------------------------------------
      subroutine smassu(ax,yx,zx,ku,kd,ks,kc)
c     input: massnr ax, hypercharge yx and charge zx.
c     output: ku,kd,ks,kc = net quark numbers (for u,d,s,c quarks).
c-----------------------------------------------------------------------
      sg=1
      if(ax.lt.0.)sg=-1
      a=sg*ax
      y=sg*yx
      z=sg*zx
      ku=nint(a+z)
      kd=nint(a-z+y)
      ks=nint(a-y)
      kc=0
      return
      end
 
c-----------------------------------------------------------------------
      function spoc(a,b,c,d,x)
c     power fctn with cutoff
c-----------------------------------------------------------------------
      spoc=0
      if(a.eq.0..and.b.eq.0.)return
      spoc =a+b*x**c
      spoc0=a+b*d**c
      spoc=amin1(spoc,spoc0)
      spoc=amax1(0.,spoc)
      return
      end
 
c-----------------------------------------------------------------------
      function sptf(x)
c     returns pt-distribution for fragmentation
c-----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/cnsta/pi,prom,piom,ainfin
      common /cquama/quama
      if(ioptf.eq.1)sptf=x*exp(-2/ptf*sqrt(x**2+quama**2))
      if(ioptf.eq.2)sptf=x*exp(-pi/4*(x**2+quama**2)/ptf**2)
      return
      end
 
c-----------------------------------------------------------------------
      function spth(x)
c     returns pt-distribution for hard scattering
c-----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      spth=pth**2*2.*x/(x**2+pth**2)**3
      return
      end
 
c-----------------------------------------------------------------------
      function sptj(x)
c     jpsi pt-distribution in 200 gev pp
c-----------------------------------------------------------------------
      a=0.95
      c=1/0.363
      z=x/a
      sptj=1/a*c**c/sgam(c)*z**(c-1)*exp(-c*z)
      return
      end
 
c-----------------------------------------------------------------------
      function sptq(x)
c     returns pt-distribution of quarks in nucleons
c-----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/cnsta/pi,prom,piom,ainfin
      if(ioptq.eq.1)sptq=4./ptq**2*x*exp(-2/ptq*x)
      if(ioptq.eq.2)sptq=pi/2./ptq**2*x*exp(-pi/4*x**2/ptq**2)
      if(ioptq.eq.3)sptq=ptq**2*2.*x/(x**2+ptq**2)**2
      return
      end
 
c-----------------------------------------------------------------------
      function ssplit(x)
c     returns splitting function
c-----------------------------------------------------------------------
c-c   ssplit= ( 0.3 + 10.0*x**7 - 10.3*x**15 ) /sqrt(x**2+.2**2)
c-c   ssplit=   0.2 +  2.4*x**2 -  2.6*x**5
c-c   ssplit=   0.4 +  2.0*x    -  2.4*x**3
c-c   ssplit= ( 0.8 +  2.0*x    -  2.8*x**3  ) /sqrt(x**2+.2**2)
      ssplit= x**3 - x**5
      return
      end
 
c-----------------------------------------------------------------------
      function ssplix(x)
c     returns splitting function
c-----------------------------------------------------------------------
c-c   ssplix= (1-x) - (1-x)**7
c-c   ssplix= 1
c-c   z=2*x-.5
c-c   ssplix= z**3 - z**5
c-c   if(z.lt.0.)ssplix=0.
c-c   if(z.gt.1.)ssplix=0.
      ssplix= x**3 - x**5
      return
      end
 
c-----------------------------------------------------------------------
c     structure functions.
c-----------------------------------------------------------------------
      subroutine staa(x,q2i,z,s)
c     initialization
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      z=sqrt(x**2+xcut**2)
      q2=amax1(4.00001,q2i)
      s=alog(alog(q2/.2**2)/alog(4/.2**2))
      return
      end
 
      function stgdv(x,q2)
c     d-valence structure fctn
      stgdv=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      stgdv=.541*x**.364*(1-x**2)**5.09/z
      return
      end
 
      function stgus(x,q2)
      stgus=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      stgus=0.25*(1-x)**7/z
      return
      end
 
      function stguv(x,q2)
      stguv=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      stguv=1.23*x**.421*(1-x**2)**3.37/z
      return
      end
 
      function stodpr(x,q2)
      stodpr=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      eta3=.763-.237*s+.026*abs(s)**2
      eta4=4.0+.627*s-.019*abs(s)**2
      gd=-.421*s+.033*abs(s)**2
      b=sbet(eta3,eta4+1)*(1+gd*eta3/(eta3+eta4+1))
      stodpr=1/b*x**eta3*(1-x)**eta4*(1+gd*x)/z
      return
      end
 
      function stoglu(x,q2)
c     glue structure fctn
      stoglu=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      stoglu=(1.56-1.71*s+.638*abs(s)**2)*
     *x**(-.949*s+.325*abs(s)**2)*
     *(1-x)**(6.0+1.44*s-1.05*abs(s)**2)*
     *(1+(9.0-7.19*s+.255*abs(s)**2)*x
     *+(-16.5*s+10.9*abs(s)**2)*x**2
     *+(15.3*s-10.1*abs(s)**2)*x**3)
     */z
      return
      end
 
      function stospi(x,q2)
      stospi=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      a=0.9
      b=5
      stospi=a*(1-x)**b/z
      return
      end
 
      function stospr(x,q2)
      stospr=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      stospr=(1.265-1.132*s+.293*abs(s)**2)*
     *x**(-.372*s-.029*abs(s)**2)*
     *(1-x)**(8.05+1.59*s-.153*abs(s)**2)*
     *(1+(6.31*s-.273*abs(s)**2)*x
     *+(-10.5*s-3.17*abs(s)**2)*x**2
     *+(14.7*s+9.80*abs(s)**2)*x**3)
     */z
      return
      end
 
      function stovpi(x,q2)
      stovpi=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      a=0.4
      b=0.7
      stovpi=x**a*(1-x)**b/sbet(a,b+1)/z
      return
      end
 
      function stovpr(x,q2)
      stovpr=0
      if(x.eq.0.)return
      call staa(x,q2,z,s)
      eta1=.419+.004*s-.007*abs(s)**2
      eta2=3.46+.724*s-.066*abs(s)**2
      gud=4.4-4.86*s+1.33*abs(s)**2
      a=sbet(eta1,eta2+1)*(1+gud*eta1/(eta1+eta2+1))
      stovpr=3/a*x**eta1*(1-x)**eta2*(1+gud*x)/z
      return
      end
 
      function stxd(x,q2)
c     valence d-quark distribution.
c     from glueck, hoffmann, reya, z. phys. c13 (1982) 119.
      call staa(x,q2,z,s)
      a=.364-.0368*s
      c=2.-.5414*s**.8
      d=5.09+.3463*s
      stxd=c*x**a*(1.-x**c)**d/sbet(d+1.,a/c)
      return
      end
 
      function stxs(x,q2)
      call staa(x,q2,z,s)
c     s-sea distribution.
c     from glueck, hoffmann, reya, z. phys. c13 (1982) 119.
      a=.0625+.1132*s**1.3
      b=12.64*s-51.70*s**1.8+38.02*s**2
      c=4.448*s
      d=7.+1.562*s
      e=.3081*s**2.5
      f=47.24+67.91*s
      stxs=a*(1.+b*x+c*x**2)*(1.-x)**d + e*exp(-f*x)
      return
      end
 
      function stxu(x,q2)
c     valence u-quark distribution.
c     from glueck, hoffmann, reya, z. phys. c13 (1982) 119.
      call staa(x,q2,z,s)
      a=.421-.0412*s
      c=2.-.6223*s**.8
      d=3.37+.4319*s
      stxu=2.*c*x**a*(1.-x**c)**d/sbet(d+1.,a/c)
      return
      end
 
      function stxus(x,q2)
c     u-sea distribution.
c     from glueck, hoffmann, reya, z. phys. c13 (1982) 119.
      call staa(x,q2,z,s)
      a=.25+.088*s**1.3
      b=.8128*s-2.003*s**1.8+.0831*s**2
      c=3.97*s
      d=7.+1.666*s
      e=.2487*s**2.5
      f=27.8+59.68*s
      stxus=a*(1.+b*x+c*x**2)*(1.-x)**d + e*exp(-f*x)
      return
      end
 
      function stxzne(x,q2)
c     structure function of neutron
      stxzne=(4.*stxd(x,q2) + stxu(x,q2) +
     *10.*stxus(x,q2) + 2.*stxs(x,q2))/9.
      return
      end
 
      function stxzpr(x,q2)
c     structure function of proton
      stxzpr=(4.*stxu(x,q2) + stxd(x,q2) +
     *10.*stxus(x,q2) + 2.*stxs(x,q2))/9.
      return
      end
 
      function stzzse(x)
c     sea structure fctn
      common /cipio/ipio
      if(ipio.eq.0)then
      stzzse=stospr(x,0.)
      elseif(ipio.eq.1)then
      stzzse=stospi(x,0.)
      else
      call utstop('stzzse: ipio must be 0 or 1!            ')
      endif
      return
      end
 
      function stzzvl(x)
c     valence structure fctn
      common /cipio/ipio
      if(ipio.eq.0)then
      stzzvl=stovpr(x,0.)/3
      elseif(ipio.eq.1)then
      stzzvl=stovpi(x,0.)/3
      else
      call utstop('stzzvl: ipio must be 0 or 1!            ')
      endif
      return
      end
 
c-----------------------------------------------------------------------
      function utacos(x)
c     returns acos(x) for -1 <= x <= 1 , acos(+-1) else
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      argum=x
      if(x.lt.-1.)then
      if(ish.ge.90)then
      call utmsg('utacos')
      write(ifch,*)'*****  argum = ',argum,' set -1'
      call utmsgf
      endif
      argum=-1.
      elseif(x.gt.1.)then
      if(ish.ge.90)then
      call utmsg('utacos')
      write(ifch,*)'*****  argum = ',argum,' set 1'
      call utmsgf
      endif
      argum=1.
      endif
      utacos=acos(argum)
      return
      end
 
c----------------------------------------------------------------------
      function utamnu(keux,kedx,kesx,kecx,modus)
c     returns min mass of cluster with given u,d,s,c content
c     keux: net u quark number
c     kedx: net d quark number
c     kesx: net s quark number
c     kecx: net c quark number
c     modus:0,1,2,3,4,5,6
c----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/cnsta/pi,prom,piom,ainfin
      common/csjcga/amnull,nsuha(7),asuha(7),entrpy,amegam,nopha
      common/paro3/iwzzzz,uentro,sigppe,sigppd,asuhax(7),asuhay(7),omega
      amnull=0.
 
      do 1 i=1,7
      if(modus.eq.0)asuha(i)=asuhay(i)
      if(modus.eq.1.and.nopha.ge.1)asuha(i)=asuhay(i)
      if(modus.eq.1.and.nopha.eq.0)asuha(i)=asuhax(i)
      if(modus.eq.2)asuha(i)=asuhay(i)
      if(modus.eq.3)asuha(i)=asuhay(i)
      if(modus.eq.4)asuha(i)=asuhax(i)    !two lowest multiplets
      if(modus.eq.5)asuha(i)=asuhay(i)    !lowest multiplet
      if(modus.eq.6)asuha(i)=sqrt(asuhay(i)**2+delmss**2)
1     nsuha(i)=0
      ke=iabs(keux+kedx+kesx+kecx)
 
      if(keux+kedx+kesx+kecx.ge.0)then
      keu=keux
      ked=kedx
      kes=kesx
      kec=kecx
      else
      keu=-keux
      ked=-kedx
      kes=-kesx
      kec=-kecx
      endif
 
      if(kec.ne.0)then
10    continue
      if(kec.lt.0)then
      kec=kec+1
      if(keu.gt.ked)then
      keu=keu-1
      else
      ked=ked-1
      endif
      amnull=amnull+1.87
      goto10
      endif
11    continue
      if(kec.gt.0)then
      kec=kec-1
      if(keu.lt.ked)then
      keu=keu+1
      else
      ked=ked+1
      endif
      amnull=amnull+1.87
      goto11
      endif
      endif
 
5     continue
      if(kes.lt.0)then
      amnull=amnull+asuha(6)
      if(keu.ge.ked)then
      keu=keu-1
      else
      ked=ked-1
      endif
      kes=kes+1
      goto5
      endif
 
6     continue
      if(ked.lt.0)then
      if(keu.ge.kes)then
      amnull=amnull+asuha(5)
      keu=keu-1
      else
      amnull=amnull+asuha(6)
      kes=kes-1
      endif
      ked=ked+1
      goto6
      endif
 
7     continue
      if(keu.lt.0)then
      if(ked.ge.kes)then
      amnull=amnull+asuha(5)
      ked=ked-1
      else
      amnull=amnull+asuha(6)
      kes=kes-1
      endif
      keu=keu+1
      goto7
      endif
 
      if(keu+ked+kes+kec.ne.ke)
     *call utstop('utamnu: sum_kei /= ke                   ')
      keq=keu+ked
      keqx=keq
 
      if(modus.eq.2.and.ke.gt.3)then
      call smasst(keu,ked,kes,kec,a,y,z)
      amnuz=smass(a,y,z)
      endif
 
      amnux=0
 
      i=4
2     i=i-1
3     continue
      if((4-i)*kes.gt.(i-1)*keq)then
      amnux=amnux+asuha(1+i)
      kes=kes-i
      keq=keq-3+i
      if(kes.lt.0)
     *call utstop('utamnu: negative kes                     ')
      if(keq.lt.0)
     *call utstop('utamnu: negative keq                     ')
      goto3
      endif
      if(i.gt.1)goto2
 
      if(keqx.gt.keq)then
      do 8 k=1,keqx-keq
      if(keu.ge.ked)then
      keu=keu-1
      else
      ked=ked-1
      endif
8     continue
      endif
      if(keu+ked.ne.keq)
     *call utstop('utamnu: keu+ked /= keq                   ')
 
9     continue
      if(keu.gt.2*ked)then
      amnux=amnux+asuha(7)
      keu=keu-3
      if(keu.lt.0)
     *call utstop('utamnu: negative keu                     ')
      goto9
      endif
      if(ked.gt.2*keu)then
      amnux=amnux+asuha(7)
      ked=ked-3
      if(ked.lt.0)
     *call utstop('utamnu: negative ked                     ')
      goto9
      endif
 
      keq=keu+ked
      if(mod(keq,3).ne.0)
     *call utstop('utamnu: mod(keq,3) /= 0                  ')
      amnux=amnux+asuha(1)*keq/3
 
      amnu=amnux
      if(modus.eq.2.and.ke.gt.3)amnu=amin1(amnux,amnuz)
      amnull=amnull+amnu
 
      if(amnull.eq.0..and.modus.gt.0)amnull=asuha(5)
 
1000  utamnu=amnull
      return
      end
 
c-----------------------------------------------------------------------
      function utamnx(jcp,jcm)
c-----------------------------------------------------------------------
c returns minimum mass for the decay of jcp---jcm (by calling utamnu).
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jcp(nflav,2),jcm(nflav,2)
 
      do i=1,nflav
      do j=1,2
      if(jcp(i,j).ne.0)goto1
      enddo
      enddo
      keu=jcm(1,1)-jcm(1,2)
      ked=jcm(2,1)-jcm(2,2)
      kes=jcm(3,1)-jcm(3,2)
      kec=jcm(4,1)-jcm(4,2)
      utamnx=utamnu(keu,ked,kes,kec,5)
      return
1     continue
 
      do i=1,nflav
      do j=1,2
      if(jcm(i,j).ne.0)goto2
      enddo
      enddo
      keu=jcp(1,1)-jcp(1,2)
      ked=jcp(2,1)-jcp(2,2)
      kes=jcp(3,1)-jcp(3,2)
      kec=jcp(4,1)-jcp(4,2)
      utamnx=utamnu(keu,ked,kes,kec,5)
      return
2     continue
 
      keu=jcp(1,1)-jcp(1,2)
      ked=jcp(2,1)-jcp(2,2)
      kes=jcp(3,1)-jcp(3,2)
      kec=jcp(4,1)-jcp(4,2)
      ke=keu+ked+kes+kec
      if(mod(ke+1,3).eq.0)then
      keu=keu+1
      elseif(mod(ke-1,3).eq.0)then
      keu=keu-1
      else
      call utstop('utamnx: no singlet possible (1)         ')
      endif
      amms=utamnu(keu,ked,kes,kec,5)
      keu=jcm(1,1)-jcm(1,2)
      ked=jcm(2,1)-jcm(2,2)
      kes=jcm(3,1)-jcm(3,2)
      kec=jcm(4,1)-jcm(4,2)
      ke=keu+ked+kes+kec
      if(mod(ke+1,3).eq.0)then
      keu=keu+1
      elseif(mod(ke-1,3).eq.0)then
      keu=keu-1
      else
      call utstop('utamnx: no singlet possible (2)         ')
      endif
      utamnx=amms+utamnu(keu,ked,kes,kec,5)
      return
      end
 
c-----------------------------------------------------------------------
      function utamny(jcp,jcm)
c-----------------------------------------------------------------------
c returns minimum mass of jcp+jcm (by calling utamnu).
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jcp(nflav,2),jcm(nflav,2),jc(nflav,2)
      do 7 nf=1,nflav
      jc(nf,1)=jcp(nf,1)+jcm(nf,1)
7     jc(nf,2)=jcp(nf,2)+jcm(nf,2)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      utamny=utamnu(keu,ked,kes,kec,5)
      return
      end

c-----------------------------------------------------------------------
      function utamnz(jc,modus)
c-----------------------------------------------------------------------
c returns minimum mass of jc (by calling utamnu).
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      utamnz=utamnu(keu,ked,kes,kec,modus)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utamst(sts,am,amin,iret)
c     input:  sts   = string  (single)
c     output: am    = mass
c             amin  = minimum mass
c             iret  = return code (=3 if am**2.lt.amin**2, 0 else)
c-----------------------------------------------------------------------
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      parameter (nsi=6,nflav=6)
      real sts(nsi,2)
      integer ic1(2),ic2(2),jc1(nflav,2),jc2(nflav,2),jc(nflav,2)
      iret=0
      am2= (sts(4,1)+sts(4,2))**2 - (sts(3,1)+sts(3,2))**2 -
     * (sts(2,1)+sts(2,2))**2 - (sts(1,1)+sts(1,2))**2
      am=-am2
      do 1 i=1,2
      ic1(i)=nint(sts(4+i,1))
1     ic2(i)=nint(sts(4+i,2))
      call iddeco(ic1,jc1)
      call iddeco(ic2,jc2)
      do 2 nf=1,nflav
      jc(nf,1)=jc1(nf,1)+jc2(nf,1)
      jc(nf,2)=jc1(nf,2)+jc2(nf,2)
      if(nf.gt.4.and.(jc(nf,1).ne.0.or.jc(nf,2).ne.0))
     *call utstop('utamst: flavour > 4                     ')
2     continue
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      amin=utamnu(keu,ked,kes,kec,4)
      if(am2.lt.amin**2)then
      iret=3
      return
      else
      am=sqrt(am2)
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utar(i1,i2,i3,x0,x1,x2,x3,xx)
c     returns the array xx with xx(1)=x0 <= xx(i) <= xx(i3)=x3
c-----------------------------------------------------------------------
      real xx(i3)
      do 1 i=1,i1-1
  1   xx(i)=x0+(i-1.)/(i1-1.)*(x1-x0)
      do 2 i=i1,i2-1
  2   xx(i)=x1+(i-i1*1.)/(i2-i1*1.)*(x2-x1)
      do 3 i=i2,i3
  3   xx(i)=x2+(i-i2*1.)/(i3-i2*1.)*(x3-x2)
      return
      end

c---------------------------------------------------------------------
      subroutine utaxis(i,j,a1,a2,a3)
c     calculates the axis defined by the ptls i,j in the i,j cm system
c---------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt                                     
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      a1=0
      a2=0
      a3=1
      pi1=pptl(1,i) 
      pi2=pptl(2,i) 
      pi3=pptl(3,i) 
      pi4=pptl(4,i) 
      pj1=pptl(1,j) 
      pj2=pptl(2,j) 
      pj3=pptl(3,j) 
      pj4=pptl(4,j) 
      p1=pi1+pj1
      p2=pi2+pj2
      p3=pi3+pj3
      p4=pi4+pj4
      p5=sqrt(p4**2-p3**2-p2**2-p1**2)
      call utlobo(1,p1,p2,p3,p4,p5,pi1,pi2,pi3,pi4)                                
      call utlobo(1,p1,p2,p3,p4,p5,pj1,pj2,pj3,pj4)                                
           err=(pi1+pj1)**2+(pi2+pj2)**2+(pi3+pj3)**2 
           if(err.gt.1e-3)then
      call utmsg('utaxis')                                                      
      write(ifch,*)'*****  err=',err
      write(ifch,*)'pi:',pi1,pi2,pi3,pi4
      write(ifch,*)'pj:',pj1,pj2,pj3,pj4
      call utmsgf                                                               
           endif
      a=sqrt( (pj1-pi1)**2 + (pj2-pi2)**2 + (pj3-pi3)**2 )
      if(a.eq.0.)return
      a1=(pi1-pj1)/a
      a2=(pi2-pj2)/a
      a3=(pi3-pj3)/a
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utchm(arp,arm,ii)
c     checks whether arp**2=0 and arm**2=0.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      double precision arp(4),arm(4),difp,difm
      difp=arp(4)**2-arp(1)**2-arp(2)**2-arp(3)**2
      difm=arm(4)**2-arm(1)**2-arm(2)**2-arm(3)**2
      if(dabs(difp).gt.1e-3*arp(4)**2
     *.or.dabs(difm).gt.1e-3*arm(4)**2)then
      call utmsg('utchm ')
      write(ifch,*)'*****  mass non zero  -  ',ii
      write(ifch,*)'jet-mass**2`s:    ',difp,difm
      write(ifch,*)'energy**2`s:      ',arp(4)**2,arm(4)**2
      write(ifch,*)(sngl(arp(i)),i=1,4)
      write(ifch,*)(sngl(arm(i)),i=1,4)
      call utmsgf
      endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utclea(nptl0)
c     overwrites istptl=2 particles in /cptl/, reduces so nptl.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer newptl(mxptl),oldptl(mxptl),ii(mxptl)
 
      ish0=ish
      if(ishsub/100.eq.18)ish=mod(ishsub,100)
 
      if(ish.ge.92)then
      write(ifch,*)('-',l=1,68)
      write(ifch,*)'sr utclea. initial.'
      write(ifch,*)('-',l=1,68)
      do 34 n=1,nptl
      write(ifch,116)iorptl(n),jorptl(n),n,ifrptl(1,n),ifrptl(2,n)
     *,idptl(n),sqrt(pptl(1,n)**2+pptl(2,n)**2),pptl(3,n),pptl(5,n)
     *,istptl(n)
34    continue
116   format(1x,i4,i4,4x,i4,4x,i4,i4,i12,3x,3(e8.2,1x),i2,i2)
      endif
 
      i=0
1     i=i+1
      if(i.gt.nptl)goto1000
      if(istptl(i).eq.2)goto2
      newptl(i)=i
      oldptl(i)=i
      goto1
 
2     i=i-1
      j=i
3     i=i+1
4     j=j+1
      if(j.gt.nptl)goto5
      newptl(j)=0
      if(istptl(j).eq.2)goto4
      newptl(j)=i
      oldptl(i)=j
      call utrepl(i,j)
      goto3
 
5     nptl=i-1
      if(nptl.eq.0)nptl0=0  
      if(nptl.eq.0)goto1000

20    n0=newptl(nptl0)
      if(n0.gt.0)then
      nptl0=n0
      else
      nptl0=nptl0-1
      if(nptl0.gt.0)goto20
      endif

      do 11 k=1,nptl
      io=iorptl(k)
      if(io.le.0)ii(k)=io
      if(io.gt.0)ii(k)=newptl(io)
11    continue
      do 12 k=1,nptl
12    iorptl(k)=ii(k)
 
      do 13 k=1,nptl
      jo=jorptl(k)
      if(jo.le.0)ii(k)=jo
      if(jo.gt.0)ii(k)=newptl(jo)
13    continue
      do 14 k=1,nptl
14    jorptl(k)=ii(k)
 
      do 15 k=1,nptl
      if1=ifrptl(1,k)
      if(if1.le.0)ii(k)=if1
      if(if1.gt.0)ii(k)=newptl(if1)
15    continue
      do 16 k=1,nptl
16    ifrptl(1,k)=ii(k)
 
      do 17 k=1,nptl
      if2=ifrptl(2,k)
      if(if2.le.0)ii(k)=if2
      if(if2.gt.0)ii(k)=newptl(if2)
17    continue
      do 18 k=1,nptl
18    ifrptl(2,k)=ii(k)
 
      do 19 k=1,nptl
      if(ifrptl(1,k).eq.0.and.ifrptl(2,k).gt.0)ifrptl(1,k)=ifrptl(2,k)
      if(ifrptl(2,k).eq.0.and.ifrptl(1,k).gt.0)ifrptl(2,k)=ifrptl(1,k)
19    continue
 
1000  continue
 
      if(ish.ge.92)then
      write(ifch,*)('-',l=1,68)
      write(ifch,*)'sr utclea. final.'
      write(ifch,*)('-',l=1,68)
      do 35 n=1,nptl
      write(ifch,116)iorptl(n),jorptl(n),n,ifrptl(1,n),ifrptl(2,n)
     *,idptl(n),sqrt(pptl(1,n)**2+pptl(2,n)**2),pptl(3,n),pptl(5,n)
     *,istptl(n)
35    continue
      write(ifch,*)('-',l=1,79)
      endif
 
      ish=ish0
      return
      end
 
c----------------------------------------------------------------------
      subroutine utfost(sia,sib,sts)
c     input :  sia = string end
c              sib = string end
c     output : sts = string (single)
c----------------------------------------------------------------------
      parameter (nsi=6)
      real sia(nsi),sib(nsi),sts(nsi,2)
      do 1 n=1,nsi
      sts(n,1)=sia(n)
1     sts(n,2)=sib(n)
      return
      end
 
c----------------------------------------------------------------------
      subroutine uthist(x1,x2,y1,y2,n,x,y,lity,lilo,text1,text2,text3)
c     writes arrays x,y in histo-format
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      real x(n),y(n)
      character*50 text1,text2,text3
      character lity*3,lilo*6
      write(ifch,100)('-',i=1,69)
      write(ifch,100)('-',i=1,69)
100   format(1x,1h',69a1,1h')
      write(ifch,101)x1,x2,0.,y1,y2,0.,'_____',0
     *,'___',lity,'___','_______',lilo,0.,3,n,0.
101   format(1x,6f8.2,2x,1h',a5,1h',i11
     */4x,1h',a3,1h',3x,1h',a3,a3,1h',3x,1h',a7,a6,1h',f8.2,2i6,f10.4)
      write(ifch,100)('-',i=1,69)
      write(ifch,102)0.,0.,text1
      write(ifch,102)0.,0.,text2
      write(ifch,102)0.,0.,text3
102   format(4x,f4.2,2x,f4.2,2x,1h',a50,1h')
      write(ifch,100)('-',i=1,69)
      do 36 k=1,n
36    write(ifch,103)k,x(k),0.,y(k),0.,0.
103   format(4x,i5,2x,5e12.4)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine uthsea
c  creates histogram of sea structure function
c-----------------------------------------------------------------------
      parameter (nstru=2049)
      common/stru/xstru(nstru),wstru(nstru)
     *,qvap(nstru),qvat(nstru),qsep(nstru),qset(nstru)
      real xar(1000),yar(1000)
      dx=1
      x1=-1
      do 3 k=1,2
      dx=dx/1000.
      x1=x1-3
      do 2 n=1,1000
      yar(n)=0
2     xar(n)=-dx/2+n*dx
      do 1 i=1,100000
      x=utinvt(nstru,xstru,qsep,rangen()*qsep(nstru))
      n=1+x/dx
      if(n.le.1000)yar(n)=yar(n)+1
1     continue
      call uthist(x1,x1+4,0.,5.,1000,xar,yar,'poc','loglog'
     *,'xaxis x                         $                 '
     *,'yaxis counts                    $                 '
     *,'title sea quark structure function                ')
3     continue
      return
      end
 
c-----------------------------------------------------------------------
      function utinvt(n,x,q,y)
c     returns x with y=q(x)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cutinv/lutinv
      real x(n),q(n)
      if(q(n).eq.0.)
     *call utstop('utinvt: q(n)=0                          ')
           if(y.lt.0.)then
      if(ish.ge.90)then
      call utmsg('utinvt')
      write(ifch,*)'*****  y=',y,' < 0'
      call utmsgf
      endif
      y=0.
           elseif(y.gt.q(n))then
      if(ish.ge.90)then
      call utmsg('utinvt')
      write(ifch,*)'*****  y=',y,' > ',q(n)
      call utmsgf
      endif
      y=q(n)
           endif
      lu=1
      lo=n
1     lz=(lo+lu)/2
      if((q(lu).le.y).and.(y.le.q(lz)))then
      lo=lz
      elseif((q(lz).lt.y).and.(y.le.q(lo)))then
      lu=lz
      else
      write(ifch,*)' '
      write(ifch,*)'q(1),y,q(n):',q(1),y,q(n)
      write(ifch,*)'lu,lz,lo:',lu,lz,lo
      write(ifch,*)'q(lu),q(lz),q(lo):',q(lu),q(lz),q(lo)
      call utstop('utinvt: no interval found               ')
      endif
      if((lo-lu).ge.2) goto1
      if(lo.le.lu)
     *call utstop('utinvt: lo.le.lu                        ')
      utinvt=x(lu)+(y-q(lu))*(x(lo)-x(lu))/(q(lo)-q(lu))
      lutinv=lu
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utkprj(proj,kmax)
c     returns kmax for proj
c-----------------------------------------------------------------------
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      real proj(nsi,nha),sil(nsi),six(nsi,nsix)
      call hdecmp(proj,sil,six)
      call utksix(six,kmax)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utksix(six,kmax)
c     returns kmax for six
c-----------------------------------------------------------------------
      parameter (nsi=6,nsix=100)
      real six(nsi,nsix)
      do 10 k=0,nsix-1
      kmax=k
      sixsqr=0
      do 20 i=1,nsi
20    sixsqr=sixsqr+six(i,k+1)**2
      if(sixsqr.le.1e-5)return
10    continue
      call utstop('utksix: dimension nsix too small        ')
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utkstr(str,kmax)
c     returns kmax for str
c-----------------------------------------------------------------------
      parameter (nsi=6,nsix=100)
      real str(nsi,nsix+1)
      do 10 k=0,nsix
      kmax=k
      strsqr=0
      do 20 i=1,nsi
20    strsqr=strsqr+str(i,k+1)**2
      if(strsqr.le.1e-5)return
10    continue
      call utstop('utkstr: dimension nsix too small        ')
      end
 
c-----------------------------------------------------------------------
      subroutine utlob2(isig,p1,p2,p3,p4,p5,x1,x2,x3,x4)
c     performs a lorentz boost, double prec.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      double precision beta(4),z(4),p1,p2,p3,p4,p5,pp,bp,x1,x2,x3,x4
     *,xx0,xx1,x10,x20,x30,x40,x4x
           if(ish.ge.90)then
      if(ish.ge.93)then
      write(ifch,*)' '
      write(ifch,101)sngl(x1),sngl(x2),sngl(x3),sngl(x4)
     *,sngl(x4**2-x3**2-x2**2-x1**2)
101   format(' utlob2:',f9.5,4f11.5)
      endif
      pp=p4**2-p3**2-p2**2-p1**2
      if(dabs(pp-p5**2).gt.1e-3*p4**2.and.dabs(pp-p5**2).gt.1e-3)then
      call utmsg('utlob2')
      write(ifch,*)'*****  p**2 .ne. p5**2'
      write(ifch,*)'p**2,p5**2: ',pp,p5**2
      write(ifch,*)'p: ',p1,p2,p3,p4,p5
      call utmsgf
      endif
      x10=x1
      x20=x2
      x30=x3
      x40=x4
           endif
      xx0=x4**2-x3**2-x2**2-x1**2
      if(p5.le.0.)then
      call utmsg('utlob2')
      write(ifch,*)'*****  p5 negative.'
      write(ifch,*)'p(5): ',p1,p2,p3,p4,p5
      call utmsgf
      call utstop('utlob2: p5 negative.                    ')
      endif
      z(1)=x1
      z(2)=x2
      z(3)=x3
      z(4)=x4
      beta(1)=-p1/p5
      beta(2)=-p2/p5
      beta(3)=-p3/p5
      beta(4)= p4/p5
      bp=0.
      do 220 k=1,3
220   bp=bp+z(k)*isig*beta(k)
      do 230 k=1,3
230   z(k)=z(k)+isig*beta(k)*z(4)
     *+isig*beta(k)*bp/(beta(4)+1.)
      z(4)=beta(4)*z(4)+bp
      x1=z(1)
      x2=z(2)
      x3=z(3)
      x4=z(4)
      if(ish.ge.93)
     *write(ifch,101)sngl(x1),sngl(x2),sngl(x3),sngl(x4)
     *,sngl(x4**2-x3**2-x2**2-x1**2)
      x4x=x4
      x4=dsqrt(xx0+x1**2+x2**2+x3**2)
      if(ish.ge.93)then
      write(ifch,101)sngl(x1),sngl(x2),sngl(x3),sngl(x4)
     *,sngl(x4**2-x3**2-x2**2-x1**2)
      write(ifch,*)' '
      endif
           if(ish.ge.90)then
      if(dabs(x4-x4x).gt.1e-2*dabs(x4).and.dabs(x4-x4x).gt.1e-2)then
      call utmsg('utlob2')
      write(ifch,*)'*****  x**2_ini .ne. x**2_fin.'
      write(ifch,*)'x1 x2 x3 x4 x**2 (initial/final/corrected):'
      write(ifch,101)sngl(x10),sngl(x20),sngl(x30),sngl(x40)
     *,sngl(x40**2-x30**2-x20**2-x10**2)
      write(ifch,101)sngl(x1),sngl(x2),sngl(x3),sngl(x4x)
     *,sngl(x4x**2-x3**2-x2**2-x1**2)
      write(ifch,101)sngl(x1),sngl(x2),sngl(x3),sngl(x4)
     *,sngl(x4**2-x3**2-x2**2-x1**2)
      call utmsgf
      endif
           endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utlobo(isig,p1,p2,p3,p4,p5,x1,x2,x3,x4)
c     performs a lorentz boost
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      real beta(4),z(4)
      if(p5.le.0.)then
      call utmsg('utlobo')
      write(ifch,*)'*****  mass <= 0.'
      write(ifch,*)'p(5): ',p1,p2,p3,p4,p5
      call utmsgf
      call utstop('utlobo: mass <= 0.                      ')
      endif
      z(1)=x1
      z(2)=x2
      z(3)=x3
      z(4)=x4
      beta(1)=-p1/p5
      beta(2)=-p2/p5
      beta(3)=-p3/p5
      beta(4)= p4/p5
      bp=0.
      do 220 k=1,3
220   bp=bp+z(k)*isig*beta(k)
      do 230 k=1,3
230   z(k)=z(k)+isig*beta(k)*z(4)
     *+isig*beta(k)*bp/(beta(4)+1.)
      z(4)=beta(4)*z(4)+bp
      x1=z(1)
      x2=z(2)
      x3=z(3)
      x4=z(4)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utloc(ar,n,a,l)
c-----------------------------------------------------------------------
      real ar(n)
      do 1 i=1,n
      l=i-1
      if(a.lt.ar(i))return
1     continue
      l=n
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utlow(cone)
c-----------------------------------------------------------------------
      character*1 cone
      if(cone.eq.'A')cone='a'
      if(cone.eq.'B')cone='b'
      if(cone.eq.'C')cone='c'
      if(cone.eq.'D')cone='d'
      if(cone.eq.'E')cone='e'
      if(cone.eq.'F')cone='f'
      if(cone.eq.'G')cone='g'
      if(cone.eq.'H')cone='h'
      if(cone.eq.'I')cone='i'
      if(cone.eq.'J')cone='j'
      if(cone.eq.'K')cone='k'
      if(cone.eq.'L')cone='l'
      if(cone.eq.'M')cone='m'
      if(cone.eq.'N')cone='n'
      if(cone.eq.'O')cone='o'
      if(cone.eq.'P')cone='p'
      if(cone.eq.'Q')cone='q'
      if(cone.eq.'R')cone='r'
      if(cone.eq.'S')cone='s'
      if(cone.eq.'T')cone='t'
      if(cone.eq.'U')cone='u'
      if(cone.eq.'V')cone='v'
      if(cone.eq.'W')cone='w'
      if(cone.eq.'X')cone='x'
      if(cone.eq.'Y')cone='y'
      if(cone.eq.'Z')cone='z'
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utlow3(cthree)
c-----------------------------------------------------------------------
      character cthree*3
      do 1 i=1,3
1     call utlow(cthree(i:i))
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utlow6(csix)
c-----------------------------------------------------------------------
      character csix*6
      do 1 i=1,6
1     call utlow(csix(i:i))
      return
      end
 
c-----------------------------------------------------------------------
      function utmom(k,n,x,q)
c     calculates kth moment for f(x) with q(i)=int[0,x(i)]f(z)dz
c-----------------------------------------------------------------------
      real x(n),q(n)
      if(n.lt.2)
     *call utstop('utmom : dimension too small             ')
      utmom=0
      do 1 i=2,n
1     utmom=utmom+((x(i)+x(i-1))/2)**k*(q(i)-q(i-1))
      utmom=utmom/q(n)
      return
      end
 
c----------------------------------------------------------------------
      subroutine utovel
c     fills array ovel(1+i,1+j) containing the logarithm (ln) of
c       i_over_j:
c     ovel(1+i,1+j)=alog(i!/j!/(j-i)!)      i>=0 j>=0
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (iovmax=100,jovmax=100)
      common/covel/ovel(1+iovmax,1+jovmax)
 
      do 1 i=0,iovmax
      ovel(1+i,1)=0.
1     ovel(1+i,1+i)=0.
      do 2 j=1,jovmax-1
      do 3 i=j+1,iovmax
3     ovel(1+i,1+j)=ovel(1+i,j)+alog((i-j+1.)/j)
2     continue
 
      if(ish.ge.90)then
      do 5 j=1,49
      do 5 i=j+1,50
      relerr=abs((exp(ovel(1+i,1+j))
     *-exp(ovel(i,j))-exp(ovel(i,1+j)))/exp(ovel(1+i,1+j)))
      if(relerr.gt.1e-4)then
      call utmsg('utovel')
      write(ifch,*)'*****  ovel(,) violates recurrence relation'
      write(ifch,*)exp(ovel(1+i,1+j)),exp(ovel(i,j))+exp(ovel(i,1+j))
      call utmsgf
      endif
5     continue
      endif
 
      if(ish.ge.95)then
      write(ifch,*)' '
      write(ifch,*)'   exp( ovel(1+i,1+j) )'
      write(ifch,*)' '
      do 6 i=0,15
6     write(ifch,*)(exp( ovel(1+i,1+j) ),j=0,min0(4,i))
      endif
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utpage
c     increases ipage by 1, changes ish
c-----------------------------------------------------------------------
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cishi/ishi
      ipage=ipage+1
      if(ipagi.le.0)return
      if(ipage.eq.1)ishi=ish
      ish=0
      if(ipage.ge.ipagi/10000.and.ipage.le.mod(ipagi,10000))ish=ishi
      return
      end
 
c----------------------------------------------------------------------
      subroutine utpart
c     fills array partx(k,n) containing the number
c       of partitions of n into at most k integers (k>=1 n>=1)  .
c     fills array part(k,n) containing  the number
c       of partitions of n into k integers (k>=1 n>=1)  .
c     fills array ipart(n,j) containing the partitions of n.
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (nquax=12)
      parameter (kpamax=nquax)
      parameter (jpamax=nquax*nquax*nquax)
      common/cpart/part(kpamax,nquax)
      common/cparta/parta(nquax),ipart(nquax,jpamax)
      common/cpartx/partx(kpamax,nquax)
      integer iswap(jpamax)
      write(ifmt,*)'execute sr utpart ...'
      do 10 n=1,nquax
      do 10 j=1,jpamax
10    ipart(n,j)=0
 
      do 6 n=1,nquax
      partx(1,n)=1
6     part(1,n)=1
      do 1 k=2,kpamax
      do 2 n=1,nquax
      u=0
      do 3 i=0,n/k
      if(n-i*k.eq.0)u=u+1
3     if(n-i*k.ne.0)u=u+partx(k-1,n-i*k)
      partx(k,n)=u
2     continue
1     continue
      do 7 n=1,nquax
      parta(n)=part(1,n)
      do 7 k=2,kpamax
      part(k,n)=partx(k,n)-partx(k-1,n)
7     parta(n)=parta(n)+part(k,n)
 
      if(ish.ge.93)then
121   format(79a1)
      write(ifch,*)' '
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   partitions of n into at most k integers'
      write(ifch,121)('-',ic=1,79)
      write(ifch,101)((k),k=1,10),kpamax
      write(ifch,121)('-',ic=1,79)
      do 8 n=1,nquax
8     write(ifch,100)n,(nint(partx(k,n)),k=1,10),nint(parta(n))
100   format(3x,'n = ',i2,3x,10i5,3x,i5)
101   format(9x,'k =',10i5,3x,i5)
      endif
 
c     n=1
c     ---
      ipart(1,1)=-1
      ipart(1,2)=-1
      ipart(1,3)=1
 
c     n>1
c     ---
      do 11 n=2,nquax
      write(ifmt,*)'sr utpart: n=',n
      if(ish.ge.95)write(ifch,*)' '
      if(ish.ge.95)write(ifch,*)'n=',n
      ipart(n,1)=-1
      ipart(n,2)=-1
      ipart(n,3)=n
      ii=3
 
c     n=n1+n2, multiply partitions of n1 and n2
c     -----------------------------------------
      do 12 n1=1,n/2
      n2=n-n1
      if(ish.ge.95)write(ifch,*)' '
      if(ish.ge.95)write(ifch,*)'n1,n2=',n1,n2
      k1max=-ipart(n1,1)
      k2max=-ipart(n2,1)
      if(k1max.lt.0.or.k2max.lt.0)
     *call utstop('utpart: kimax negative                  ')
      j1=2
      do 13 k1=1,k1max
      l1=-ipart(n1,j1)
      if(l1.lt.0)
     *call utstop('utpart: l1 negative (1)                 ')
      j2=2
      do 15 k2=1,k2max
      l2=-ipart(n2,j2)
      if(l2.lt.0)
     *call utstop('utpart: l2 negative (1)                 ')
 
      ipart(n,1)=ipart(n,1)-1
      ii=ii+1
      if(ii.gt.jpamax)goto1000
      ipart(n,ii)=-l1-l2
      ii0=ii+1
      do 17 i1=1,l1
      ii=ii+1
      if(ii.gt.jpamax)goto1000
17    ipart(n,ii)=ipart(n1,j1+i1)
      do 18 i2=1,l2
      ii=ii+1
      if(ii.gt.jpamax)goto1000
18    ipart(n,ii)=ipart(n2,j2+i2)
 
      if(ish.ge.95)then
      write(ifch,*)' '
      write(ifch,*)'k1,k2=',k1,k2
      write(ifch,103)n,-ipart(n,1),nint(parta(n))
103   format(/3x,'n = ',i2,'   p(n) = ',i4,'   p0(n) = ',i4/)
      jy=2
      do 29 ky=1,-ipart(n,1)
      ly=-ipart(n,jy)
      write(ifch,102)ky,ly,(ipart(n,i),i=jy+1,jy+min0(20,ly))
102   format(3x,i2,3x,i2,3x,20i3)
29    jy=jy+ly+1
      endif
 
      ipri=0
      do 19 i=ii0,ii-1
      do 19 j=i+1,ii
      if(ipart(n,i).lt.ipart(n,j))then
      ipri=1
      isti=ipart(n,i)
      ipart(n,i)=ipart(n,j)
      ipart(n,j)=isti
      endif
19    continue
      if(ipri.eq.1.and.ish.ge.95)then
      write(ifch,*)' '
      write(ifch,*)'sequence changed'
      endif
 
      jx=2
      do 28 kx=1,-ipart(n,1)-1
      lx=-ipart(n,jx)
      if(lx.lt.0)
     *call utstop('utpart: lx negative                     ')
      if(lx.eq.l1+l2)then
      do 22 l=1,lx
      if(ipart(n,jx+l).ne.ipart(n,ii0-1+l))goto23
22    continue
      if(ish.ge.95)then
      write(ifch,*)' '
      write(ifch,*)'exists already'
      endif
      ipri=1
      do 33 i=ii0-1,ii
33    ipart(n,i)=0
      ii=ii0-2
      ipart(n,1)=ipart(n,1)+1
      goto30
23    continue
      endif
28    jx=jx+lx+1
30    continue
 
      if(ipri.eq.1.and.ish.ge.95)then
      write(ifch,103)n,-ipart(n,1),nint(parta(n))
      jy=2
      do 32 ky=1,-ipart(n,1)
      ly=-ipart(n,jy)
      write(ifch,102)ky,ly,(ipart(n,i),i=jy+1,jy+min0(20,ly))
32    jy=jy+ly+1
      endif
 
15    j2=j2+l2+1
13    j1=j1+l1+1
12    continue
 
c     ordering
c     --------
      km=-ipart(n,1)
      if(km.lt.0)
     *call utstop('utpart: km negative                     ')
      j1=2
      do 20 k1=1,km-1
      l1=-ipart(n,j1)
      if(l1.lt.0)
     *call utstop('utpart: l1 negative (2)                 ')
      j2=2
      do 21 k2=1,km
      l2=-ipart(n,j2)
      if(l2.lt.0)
     *call utstop('utpart: l2 negative (2)                 ')
      if(k2.le.k1)goto21
      if(l1.gt.l2)then
      do 24 l=1,l2+1
24    iswap(l)=ipart(n,j2-1+l)
      do 25 l=1,j2-j1
      i=j2-l
25    ipart(n,i+l2+1)=ipart(n,i)
      do 26 l=1,l2+1
26    ipart(n,j1-1+l)=iswap(l)
      l1=-ipart(n,j1)
      if(l1.lt.0)
     *call utstop('utpart: l1 negative (2)                 ')
      if(ish.ge.95)then
      write(ifch,*)' '
      write(ifch,*)'order changed.      k1,k2=',k1,k2
      write(ifch,103)n,-ipart(n,1),nint(parta(n))
      jy=2
      do 31 ky=1,-ipart(n,1)
      ly=-ipart(n,jy)
      write(ifch,102)ky,ly,(ipart(n,i),i=jy+1,jy+min0(20,ly))
31    jy=jy+ly+1
      endif
      endif
21    j2=j2+l2+1
20    j1=j1+l1+1
 
      if(ish.ge.93.and.n.le.8)then
      write(ifch,113)('-',ic=1,79),n,-ipart(n,1),('-',ic=1,79)
113   format(/79a1/7x,'n = ',i2,'   --->   ',i4,' partitions'/79a1)
      jy=2
      do 27 ky=1,-ipart(n,1)
      ly=-ipart(n,jy)
      write(ifch,112)ky,(ipart(n,i),i=jy+1,jy+min0(20,ly))
112   format(2x,i2,'. partition:',3x,20i3)
27    jy=jy+ly+1
      endif
      if(-ipart(n,1).ne.nint(parta(n)))
     *call utstop('utpart: # of partitions wrong           ')
11    continue
 
      return
 
1000  write(ifch,*)' '
      write(ifch,*)('*',j=1,79)
      write(ifch,*)'***** n=',n
      write(ifch,*)'***** jpamax=',jpamax
      write(ifch,*)('*',j=1,79)
      call utstop('utpart: dimension jpamax too small.     ')
 
      end
 
c-----------------------------------------------------------------------
      function utpcm(a,b,c)
c     calculates cm momentum for a-->b+c
c-----------------------------------------------------------------------
      val=(a**2-b**2-c**2)**2-(2.*b*c)**2
      if(val.lt.0..and.val.gt.-1e-4)then
      utpcm=0
      return
      endif
      utpcm=sqrt(val)/(2.*a)
      return
      end
 
c-----------------------------------------------------------------------
      function utquad(m,x,f,k)
c     performs an integration according to simpson
c-----------------------------------------------------------------------
      real x(m),f(m)
      utquad=0
      do 1 i=1,k-1
  1   utquad=utquad+(f(i)+f(i+1))/2*(x(i+1)-x(i))
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utquaf(fu,n,x,q,wo,x0,x1,x2,x3)
c     returns q(i) = integral [x(1)->x(i)] fu(x) dx
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
c MAB
      external fu
      real x(n),q(n),wo(n)
      parameter (m=10)
      real xa(m),fa(m)
      if(x1.lt.x0.or.x2.lt.x1.or.x3.lt.x2)then
      if(ish.ge.90)then
      call utmsg('utquaf')
      write(ifch,*)'   xi=',x0,x1,x2,x3
      call utmsgf
      endif
      endif
      call utar(n/3,n*2/3,n,x0,x1,x2,x3,x)
      q(1)=0
      do 2 i=2,n
      do 3 k=1,m
      z=x(i-1)+(k-1.)/(m-1.)*(x(i)-x(i-1))
      xa(k)=z
3     fa(k)=fu(z)
      q(i)=q(i-1)+utquad(m,xa,fa,m)
2     continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utremb(proj,targ,ii)
c     remembers nevt,nstr,nptl,proj,targ
c-----------------------------------------------------------------------
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      parameter(mxstr=20000)
      common/cstr/nstr,pstr(5,mxstr),icstr(4,mxstr),iorstr(mxstr)
     *,irlstr(mxstr),rotstr(3,mxstr),xorstr(4,mxstr)
      common/cevt/phievt,nevt,bimevt,kolevt,colevt,pmxevt,egyevt
     *,npjevt,ntgevt
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      real proj(nsi,nha),targ(nsi,nha)
      common/cremb/nevtre(2),nstrre(2),nptlre(2)
     *,projre(2,nsi,nha),targre(2,nsi,nha)
 
      nevtre(ii)=nevt
      nstrre(ii)=nstr
      nptlre(ii)=nptl
      do 56 m=1,nha
      sm=0
      do 57 n=1,nsi
      sm=sm+proj(n,m)**2
57    projre(ii,n,m)=proj(n,m)
      if(m.ge.3.and.sm.lt.1e-5)goto58
56    continue
58    continue
      do 53 m=1,nha
      sm=0
      do 54 n=1,nsi
      sm=sm+targ(n,m)**2
54    targre(ii,n,m)=targ(n,m)
      if(m.ge.3.and.sm.lt.1e-5)goto55
53    continue
55    continue
      return
 
      entry utrest(proj,targ,ii)
      nevt=nevtre(ii)
      nstr=nstrre(ii)
      nptl=nptlre(ii)
      do 66 m=1,nha
      sm=0
      do 67 n=1,nsi
      sm=sm+projre(ii,n,m)**2
67    proj(n,m)=projre(ii,n,m)
      if(m.ge.3.and.sm.lt.1e-5)goto68
66    continue
68    continue
      do 63 m=1,nha
      sm=0
      do 64 n=1,nsi
      sm=sm+targre(ii,n,m)**2
64    targ(n,m)=targre(ii,n,m)
      if(m.ge.3.and.sm.lt.1e-5)goto65
63    continue
65    continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utrepl(i,j)
c     i is replaced by j in /cptl/
c-----------------------------------------------------------------------
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common/c2ptl/iaaptl(mxptl),radptl(mxptl),amiptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      do 1 k=1,5
1     pptl(k,i)  =pptl(k,j)
      iorptl(i)  =iorptl(j)
      idptl(i)   =idptl(j)
      istptl(i)  =istptl(j)
      do 2 k=1,2
2     tivptl(k,i)=tivptl(k,j)
      do 3 k=1,2
3     ifrptl(k,i)=ifrptl(k,j)
      jorptl(i)  =jorptl(j)
      do 4 k=1,4
4     xorptl(k,i)=xorptl(k,j)
      do 5 k=1,4
5     ibptl(k,i) =ibptl(k,j)
      iclptl(i)  =iclptl(j)
      nqjptl(i)  =nqjptl(j)
      iaaptl(i)  =iaaptl(j)
      radptl(i)  =radptl(j)
      amiptl(i)  =amiptl(j)

      return
      end
 
c-----------------------------------------------------------------------
      subroutine utresm(icp1,icp2,icm1,icm2,amp,idpr,iadj,ireten)
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer icm(2),icp(2),jcm(nflav,2),jcp(nflav,2)
      icm(1)=icm1
      icm(2)=icm2
      icp(1)=icp1
      icp(2)=icp2
      CALL IDDECO(ICM,JCM)
      CALL IDDECO(ICP,JCP)
      do 37 nf=1,nflav
      do 37 k=1,2
37    jcP(nf,k)=jcp(nf,k)+jcm(nf,k)
      CALL IDENCO(JCP,ICP,IRETEN)
      IDP=IDTRA(ICP,0,0,3)
      call idres(idp,amp,idpr,iadj)
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utrot2(isig,ax,ay,az,x,y,z)
c     performs a rotation, double prec.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      double precision ax,ay,az,x,y,z,rx,ry,rz
     *,alp,bet,cosa,sina,cosb,sinb,xs,ys,zs
         if(ax**2.eq.0.and.ay**2.eq.0.and.az**2.eq.0.)then
      write(ifch,*)' '
      write(ifch,*)'ax**2,ay**2,az**2:',ax**2,ay**2,az**2
      write(ifch,*)'ax,ay,az:',ax,ay,az
      call utstop('utrot2: zero vector.                    ')
         endif
         if(az.ge.0.)then
      rx=ax
      ry=ay
      rz=az
         else
      rx=-ax
      ry=-ay
      rz=-az
         endif
      if(rz**2+ry**2.ne.0.)then
      alp=dabs(dacos(rz/dsqrt(rz**2+ry**2)))*sign(1.,sngl(ry))
      bet=
     *dabs(dacos(dsqrt(rz**2+ry**2)/dsqrt(rz**2+ry**2+rx**2)))*
     *sign(1.,sngl(rx))
      else
      alp=3.1415927/2
      bet=3.1415927/2
      endif
      cosa=dcos(alp)
      sina=dsin(alp)
      cosb=dcos(bet)
      sinb=dsin(bet)
           if(isig.gt.0)then
      xs=x*cosb-y*sina*sinb-z*cosa*sinb
      ys=       y*cosa     -z*sina
      zs=x*sinb+y*sina*cosb+z*cosa*cosb
           elseif(isig.lt.0)then
      xs= x*cosb            +z*sinb
      ys=-x*sinb*sina+y*cosa+z*cosb*sina
      zs=-x*sinb*cosa-y*sina+z*cosb*cosa
           endif
      x=xs
      y=ys
      z=zs
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utrota(isig,ax,ay,az,x,y,z)
c     performs a rotation
c-----------------------------------------------------------------------
         if(az.ge.0.)then
      rx=ax
      ry=ay
      rz=az
         else
      rx=-ax
      ry=-ay
      rz=-az
         endif
      alp=abs(utacos(rz/sqrt(rz**2+ry**2)))*sign(1.,ry)
      bet=
     *abs(utacos(sqrt(rz**2+ry**2)/sqrt(rz**2+ry**2+rx**2)))*sign(1.,rx)
      cosa=cos(alp)
      sina=sin(alp)
      cosb=cos(bet)
      sinb=sin(bet)
           if(isig.gt.0)then
      xs=x*cosb-y*sina*sinb-z*cosa*sinb
      ys=       y*cosa     -z*sina
      zs=x*sinb+y*sina*cosb+z*cosa*cosb
           elseif(isig.lt.0)then
      xs= x*cosb            +z*sinb
      ys=-x*sinb*sina+y*cosa+z*cosb*sina
      zs=-x*sinb*cosa-y*sina+z*cosb*cosa
           endif
      x=xs
      y=ys
      z=zs
      return
      end
 
c-----------------------------------------------------------------------
      function utsqsi(sil)
c-----------------------------------------------------------------------
      parameter (nsi=6)
      real sil(nsi)
      utsqsi=0
      do 1 n=1,nsi
1     utsqsi=utsqsi+sil(n)**2
      return
      end
 
c-----------------------------------------------------------------------
      function utsqst(str)
c-----------------------------------------------------------------------
      parameter (nsi=6,nsix=100)
      parameter (nha=nsix+2)
      real str(nsi,2)
      utsqst=0
      do 1 n=1,nsi
1     utsqst=utsqst+str(n,1)**2+str(n,2)**2
      return
      end
 
c-----------------------------------------------------------------------
      function utsstr(str)
c     calculates the mass of string str
c-----------------------------------------------------------------------
      parameter (nsi=6)
      real str(nsi,2)
      utsstr=(str(4,1)+str(4,2))**2-(str(1,1)+str(1,2))**2
     *-(str(2,1)+str(2,2))**2-(str(3,1)+str(3,2))**2
      return
      end
 
c-----------------------------------------------------------------------
      subroutine utstop(text)
c     returns error message and stops execution
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      character text*40,star*1,txt*6
      star='*'
      do 1 i=1,2
      if(i.eq.1)ifi=ifch
      if(i.eq.2)ifi=ifmt
      write(ifi,*)' '
      write(ifi,*)(star,j=1,77)
      write(ifi,*)'***** stop in ',text
      write(ifi,*)'***** current event number: ',nrevt+1
      write(ifi,*)'***** current page number: ',ipage
      write(ifi,*)(star,j=1,77)
1     write(ifi,*)' '
      stop
      entry utmsg(txt)
      imsg=imsg+1
      if(ish.ne.90.or.ishsub.ne.0)write(ifch,*)' '
      write(ifch,*)('-',j=1,77)
      write(ifch,100)txt,nrevt+1,ipage
100   format(1x,'***** msg from ',a6,'.   eps:',i7,i5)
      return
      entry utmsgf
      if(ish.eq.90.and.ishsub.eq.0)return
      write(ifch,*)('-',j=1,77)
      write(ifch,*)' '
      end
 
c-----------------------------------------------------------------------
      subroutine uttain(i,x,y,z,t,n,iopt)
c     returns intersection of ptl-i-trajectory with taus-line.
c     n=0 if ok, n=1 if ptl lives later, n=2 if earlier,
c       n=9 if tiv1>tiv2, n=3,4,5 else.
c     iopt=0: formation time considered, iopt=1 if not
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      double precision vv,zza,zz,tt,xo3,xo4,ti1,ti2,derr,dd
      double precision ttp,zzp,ttt,zzt,vvt,vvp
 
      derr=1e-2
      ttp=tpro*ttaus
      zzp=zpro*ttaus
      ttt=ttar*ttaus
      zzt=ztar*ttaus
      vv=pptl(3,i)/pptl(4,i)
      xo3=xorptl(3,i)
      xo4=xorptl(4,i)
      zza=xo3-xo4*vv
      if(iopt.eq.0)ti1=tivptl(1,i)
      if(iopt.eq.1)ti1=xo4
      ti2=tivptl(2,i)
 
      if(ti1.gt.ti2)goto1009
 
      zfi=xo3+(ti2-xo4)*vv
      call uttaus(zfi,tzfi,szfi)
      tfi=tzfi
      if(tfi.ge.sngl(ti2))goto1002
 
      zin=xo3+(ti1-xo4)*vv
      call uttaus(zin,tzin,szin)
      tin=tzin
      if(tin.le.sngl(ti1))goto1001
 
           if(sngl(ttaus).le.0.)then
      tt=ttaus
      zz=xo3+(tt-xo4)*vv
      if(tt.lt.ti1.or.tt.ge.ti2)goto10031
      goto1000
           else
      vvt=zzt/ttt
      vvp=zzp/ttp
      tt=(ttt+(zza-zzt)*vvt)/(1-vv*vvt)
      zz=xo3+(tt-xo4)*vv
      if(zz.le.zzt)then
      if(tt.lt.ti1.or.tt.ge.ti2)goto10032
      goto1000
      endif
      tt=(ttp+(zza-zzp)*vvp)/(1-vv*vvp)
      zz=xo3+(tt-xo4)*vv
      if(zz.ge.zzp)then
      if(tt.lt.ti1.or.tt.ge.ti2)goto10033
      goto1000
      endif
      dd=1-vv**2
      if(sngl(dd).eq.0..and.vv.gt.0.)then
      tt=-(ttaus**2+zza**2)/2./zza
      elseif(sngl(dd).eq.0..and.vv.lt.0.)then
      tt=(ttaus**2+zza**2)/2./zza
      else
      tt=(zza*vv+dsqrt(zza**2+ttaus**2*dd))/dd
      endif
      zz=xo3+(tt-xo4)*vv
      if(tt.lt.ti1.or.tt.ge.ti2)goto10034
      if(tt.lt.0.)goto10035
      if(zz.le.zzt.or.zz.ge.zzp)goto1004
      if(dabs(ttaus**2-(tt+zz)*(tt-zz)).gt.derr*ttaus**2)goto1005
      goto1000
           endif
 
1000  n=0
10001 t=tt
      z=zz
      x=xorptl(1,i)+(t-xorptl(4,i))*pptl(1,i)/pptl(4,i)
      y=xorptl(2,i)+(t-xorptl(4,i))*pptl(2,i)/pptl(4,i)
      return
1001  n=1
      return
1002  n=2
      return
10031 n=31
      goto1003
10032 n=32
      goto1003
10033 n=33
      goto1003
10034 n=34
1003  if(dabs(tt-ti1).le.derr*dabs(tt))goto1000
      if(dabs(tt-ti2).le.derr*dabs(tt))goto1000
      if(ish.ge.90)then
      call utmsg('uttain')
      write(ifch,*)'*****  ti1 < tt < ti2   not fulfilled - ',n
      write(ifch,*)sngl(ti1),sngl(tt),sngl(ti2)
      call utmsgf
      endif
      goto10001
10035 continue
      if(ish.ge.90)then
      call utmsg('uttain')
      write(ifch,*)'*****  tt < 0     ( ',tt,' )'
      write(ifch,*)'vv,dd:',vv,dd
      write(ifch,*)'zza,ttaus:',zza,ttaus
      call utmsgf
      endif
      goto10001
1004  n=4
      if(dabs(zz-zzt).le.derr*dabs(zz))goto1000
      if(dabs(zz-zzp).le.derr*dabs(zz))goto1000
      if(ish.ge.90)then
      call utmsg('uttain')
      write(ifch,*)'*****  zzt < zz < zzp   not fulfilled'
      write(ifch,*)sngl(zzt),sngl(zz),sngl(zzp)
      call utmsgf
      endif
      goto10001
1005  n=5
      if(dabs(ttaus**2-(tt+zz)*(tt-zz)).le.derr)goto1000
      if(ish.ge.90)then
      call utmsg('uttain')
      write(ifch,*)'*****  ttaus**2 .ne. (tt+zz)*(tt-zz)'
      write(ifch,*)sngl(ttaus**2),sngl((tt+zz)*(tt-zz))
      call utmsgf
      endif
      goto10001
1009  n=9
      return
      end
 
c-----------------------------------------------------------------------
      subroutine uttaix(i,tau,zor,tor,z,t)
c     returns intersection z,t of ptl-i-trajectory with hyperbola h.
c        h: (t-tor)**2-(z-zor)**2=tau**2 .
c        zor, tor double precision.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      double precision tor,zor,tors,zors,vv,cc,dd,ttau,derr,tt,zz
      common/cnsta/pi,prom,piom,ainfin
      derr=1e-3
      ttau=tau
      zors=xorptl(3,i)-zor
      tors=xorptl(4,i)-tor
      vv=pptl(3,i)/pptl(4,i)
           if(dabs(vv).gt.1d0)then
      if(ish.ge.90.and.dabs(vv).gt.1.001d0)then
      call utmsg('uttaix')
      write(ifch,*)'*****  !v! > 1'
      write(ifch,*)'v: ',vv
      write(ifch,*)'p,e: ',pptl(3,i),pptl(4,i)
      call utmsgf
      endif
      vv=vv/dabs(vv)
           endif
      cc=zors-tors*vv
      dd=1d0-vv**2
           if(dd.eq.0d0.and.cc.eq.0d0)then
      if(tau.eq.0.)tt=0.
      if(tau.ne.0.)tt=ainfin
      zz=tt
      goto1000
           elseif(dd.eq.0d0)then
      tt=-(ttau**2+cc**2)/2./cc/vv
           elseif(dd.lt.1e-8)then
      tt=-(ttau**2+cc**2)/2./cc/vv
      call utmsg('uttaix')
      write(ifch,*)'*****  dd = ',dd,'    treated as zero'
      call utmsgf
           else
      tt=(cc*vv+dsqrt(cc**2+ttau**2*dd))
      tt=tt/dd
           endif
      zz=cc+tt*vv
      if(dabs(ttau**2-(tt+zz)*(tt-zz)).gt.derr*ttau**2
     *.and.dabs(ttau**2-(tt+zz)*(tt-zz)).gt.derr
     *.and.tors**2-zors**2.lt.1e6)then
      if(ish.ge.90)then
      call utmsg('uttaix')
      write(ifch,*)'*****  ttau**2 .ne. (tt+zz)*(tt-zz)'
      write(ifch,*)sngl(ttau**2),sngl((tt+zz)*(tt-zz))
      write(ifch,*)'tau,t,z:'
      write(ifch,*)tau,tt,zz
      write(ifch,*)'#,id(ptl):',i,idptl(i)
      write(ifch,*)'zor,tor(str):',zor,tor
      write(ifch,*)'zors,tors,p,e(ptl):'
      write(ifch,*)sngl(zors),sngl(tors),pptl(3,i),pptl(4,i)
      call utmsgf
      endif
      endif
1000  z=zz+zor
      t=tt+tor
      return
      end
 
c-----------------------------------------------------------------------
      subroutine uttaus(z,tz,sz)
c     returns time tz and inv z-coord sz corresponding to ttaus and z.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat 
      double precision zz
      double precision ttp,zzp,ttt,zzt
      common/cnsta/pi,prom,piom,ainfin
      ttp=tpro*ttaus
      zzp=zpro*ttaus
      ttt=ttar*ttaus
      zzt=ztar*ttaus
      zz=z
      tz=ttaus
      sz=z
      if(ttaus.le.0.)return
           if(zz.le.zzt)then
      tz=ttt+(zz-zzt)*zzt/ttt
      sz=ttaus*detat+(zz-zzt)/ttar
           elseif(zz.ge.zzp)then
      tz=ttp+(zz-zzp)*zzp/ttp
      sz=ttaus*detap+(zz-zzp)/tpro
           else
      if(sngl(ttaus).ge.ainfin)then
      tz=ttaus
      sz=0
      if(ish.ge.90)then
      call utmsg('uttaus')
      write(ifch,*)'*****  large ttaus; set tz=ttaus, sz=0'
      write(ifch,*)'ttaus=',ttaus,'zz=',zz
      call utmsgf
      endif
      else
      tz=dsqrt(ttaus**2+zz**2)
      sz=ttaus*0.5*dlog((tz+zz)/(tz-zz))
      endif
           endif
      return
      end
 
c-----------------------------------------------------------------------
      subroutine uttima(text)
c     returns time.
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      character*15 text
      timaa=0
      tima=0
c-c   call timax(timaa)
c-c   call timad(tima)
      if(text.eq.'               ')return
      if(ish.ge.91)write(ifch,*)' '
      write(ifch,100)text,tima/5.,timaa/5.
      if(ish.ge.91)write(ifch,*)' '
100   format(1x,a15,5x,f12.5,5x,f12.5)
      return
 
      entry uttimt
c-c   call timast(1e10)
      return
 
      entry uttimx(timaax)
      timaax=0
c-c   call timax(timaax)
      return
      end
 
c----------------------------------------------------------------------
      subroutine uttucl
c     fills array tucl(1+k,1+n) containing the logarithm (ln) of
c       the number of k-tupels n_i with sum_i n_i = n:
c       tucl(1+k,1+n)=dlog((n+1)*(n+2)...*(n+k-1)/(k-1)!)  k>=0 n>=0 .
c     double precision tucl.
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      double precision tucl,dd
      parameter (ktumax=100,ntumax=100)
      common/ctucl/tucl(1+ktumax,1+ntumax)
 
      do 6 n=0,ntumax
      tucl(1,1+n)=-100
6     tucl(2,1+n)=0
      do 7 k=2,ktumax
7     tucl(1+k,1)=0
      do 2 n=1,ntumax
      do 1 k=2,ktumax
      dd=(n+k-1.)/(k-1.)
      tucl(1+k,1+n)=tucl(k,1+n)+dlog(dd)
1     continue
2     continue
 
      if(ish.ge.93)then
      do 5 k=2,50
      do 5 n=1,50
      relerr=dabs((dexp(tucl(1+k,1+n))
     *-dexp(tucl(1+k,n))-dexp(tucl(k,1+n)))/dexp(tucl(1+k,1+n)))
      if(relerr.gt.1e-4)then
      call utmsg('uttucl')
      write(ifch,*)'*****  tucl(,) violates recurrence relation'
      write(ifch,*)sngl(dexp(tucl(1+k,1+n)))
     *,sngl(dexp(tucl(1+k,n))+dexp(tucl(k,1+n)))
      call utmsgf
      endif
5     continue
      endif
 
      if(ish.ge.93)then
121   format(79a1)
      write(ifch,*)' '
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   degeneracy of energy levels of k-dimensional'
     *,' oscillator'
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   =='
     *,'     number of k-tupels of length n (sum_i n_i = n)'
      write(ifch,121)('-',ic=1,79)
      write(ifch,*)'   k:','   1','   2','   3','   4'
      write(ifch,121)('-',ic=1,79)
      do 8 n=0,50
8     write(ifch,*)n,(sngl(dexp( tucl(1+k,1+n) )),k=1,4)
      endif
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine veanif(in,io,im,ic,ih,id,si,ls,te,en,pl
     *,px,yp,is,ii,nt,na,bm,pp,pr,ai,lp,mp
     *,lt,mt,pb,jb,jf,iv,iu,np,pv,pc,ne)
c     transfers variables to veanly
c-----------------------------------------------------------------------
      common/cvsn/iversn
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (ndep=129,ndet=129,kollmx=5000)
      common /col/rmproj,rmtarg,bmax,bimp,koll,nproj,ntarg
     *,xdep(ndep),qdep(ndep),wdep(ndep),xdet(ndet),qdet(ndet),wdet(ndet)
     *,nrproj(kollmx),nrtarg(kollmx),distce(kollmx),nord(kollmx)
     *,coord(4,kollmx)
      common/cnsta/pi,prom,piom,ainfin
      parameter(mxpar=99)
      character*6 parcha(mxpar),pc(mxpar)
      double precision parval(mxpar),pv(mxpar)
      common/paramc/npar,idmmy,parval,parcha
      dimension pb(99),jb(99,2),jf(99,2)
 
      in=iversn
      io=ifop
      im=ifmt
      ic=ifch
      ih=ifhi
      id=ifdt
      si=sigppi
      ls=labsys
      te=tensn
      en=engy
      pl=pnll
      px=pnllx
      yp=yhaha
      is=ish
      ii=ichoic
      nt=ntevt
      na=naevt
      bm=bmax
      pp=pi
      pr=prom
      ai=ainfin
      lp=laproj
      mp=maproj
      lt=latarg
      mt=matarg
      do 1 i=1,99
      pb(i)=prob(i)
      do 1 j=1,2
      jb(i,j)=icbac(i,j)
1     jf(i,j)=icfor(i,j)
      iv=ivi
      iu=isup
      np=npar
      do 2 n=1,mxpar
      pv(n)=parval(n)
2     pc(n)=parcha (n)
      ne=nrevt
 
      if(ichoic.eq.4.or.ichoic.eq.5)bm=0
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine veanly
c     analyses venus results
c-----------------------------------------------------------------------
      common/cifipr/ifipr
      parameter(mxpar=99,nflav=6,kmaxmx=10000,mxscat=5000)
      parameter (mxtxt=mxpar,mxtrig=mxpar,mxidp=mxpar)
      character das*1,expmt*5,xobsrv*6,yobsrv*6,ili*6,lilo*3,mde*3
     *,cone*1,fch*1
      character*6 parcha(mxpar)
      double precision parval(mxpar)
      real parfi(mxpar)
      common/canar/ar(kmaxmx,5)
      real yacevt(kmaxmx),yactot(kmaxmx)
      real trig1(mxtrig),trig2(mxtrig)
      character*6 trig(mxtrig),trg
      integer idptl(mxidp),nacevt(kmaxmx)
      real zactot(kmaxmx)
      real arscat(mxscat,2)
      integer jc1(nflav,2),jc2(nflav,2)
      character trigh*22,idptlh*6,parh*14,modana*3,muopt*13
      character*50 txt(mxtxt)
      real xt(mxtxt),yt(mxtxt)
      real facmo(5)
      dimension prob(99),icbac(99,2),icfor(99,2)
      call veanif(iversn,ifop,ifmt,ifch,ifhi,ifdt,sigppi,labsys,tensn
     *,engy,pnll,pnllx,yhaha,ish,ichoic,ntevt,naevt,bmax,pi,prom,ainfin
     *,laproj,maproj,latarg,matarg,prob,icbac,icfor,ivi,isup,npar,parval
     *,parcha,nrevt)
      ish0=ish
      ifhi0=ifhi
      das='-'
      fch='-'
      zero=0
      izero=0
      levt=1
      lptl=3
 
99    continue
      ish=ish0
      ifhi=ifhi0
      ioptn=0
      crofac=1
      pltfac=1
      if(isup.ne.1)write(ifmt,*)
     *'analysis: sup or nos or ret or ren or stp?'
      read(ifop,*)modana
      call utlow3(modana)
      if(isup.ne.1.and.ivi.eq.1)write(ifmt,*)modana
      if(modana.eq.'stp')stop
      ifipr=1
      if(modana.eq.'ret')return
      ifipr=0
      if(modana.eq.'ren')return
      do 10 k=1,kmaxmx
      zactot(k)=0
10    yactot(k)=0
      do 1 k=1,kmaxmx
      do 1 l=1,5
1     ar(k,l)=0
      do 6 i=1,mxidp
6     idptl(i)=0
      do 38 i=1,mxtrig
      trig1(i)=0
      trig2(i)=0
38    trig(i)='      '
 
c     options
c     -------
3     continue
      if(modana.ne.'sup')write(ifmt,*)'x-observ?[a6;hlp]  '
     *,'y-observ?[a6;hlp]  #trig?[i]  #had?[i]'
      read(ifop,fmt=*,err=3)xobsrv,yobsrv,ntrig,nidp
      do 90 i=1,6
      call utlow(xobsrv(i:i))
90    call utlow(yobsrv(i:i))
      if(ivi.eq.1.and.modana.ne.'sup')
     *write(ifmt,1030)xobsrv,yobsrv,ntrig,nidp
1030  format(1x,a6,1x,a6,i2,i2)
      if(xobsrv.eq.'stp   '.or.xobsrv.eq.'stop  ')return
      if(xobsrv.eq.'hlp   '.or.xobsrv.eq.'help  ')write(ifmt,104)
      if(xobsrv.eq.'hlp   '.or.xobsrv.eq.'help  ')goto3
      if(xobsrv.eq.'hlpptl')write(ifmt,108)
      if(xobsrv.eq.'hlpptl')goto3
      if(xobsrv.eq.'hlpevt')write(ifmt,109)
      if(xobsrv.eq.'hlpevt')goto3
      if(xobsrv.eq.'hlppar')write(ifmt,120)
      if(xobsrv.eq.'hlppar')goto3
      if(yobsrv.eq.'stp   '.or.yobsrv.eq.'stop  ')return
      if(yobsrv.eq.'hlp   '.or.yobsrv.eq.'help  ')write(ifmt,104)
      if(yobsrv.eq.'hlp   '.or.yobsrv.eq.'help  ')goto3
      if(yobsrv.eq.'hlpptl')write(ifmt,108)
      if(yobsrv.eq.'hlpptl')write(ifmt,107)
      if(yobsrv.eq.'hlpptl')goto3
      if(yobsrv.eq.'hlpevt')write(ifmt,109)
      if(yobsrv.eq.'hlpevt')goto3
      if(yobsrv.eq.'hlppar')write(ifmt,120)
      if(yobsrv.eq.'hlppar')goto3
      if(xobsrv(4:6).ne.yobsrv(4:6))write(ifmt,*)
     *'x-observable and y-observable have to be of the same type'
      if(xobsrv(4:6).ne.yobsrv(4:6))goto3
104   format(
     * ' !----------------------------------------------------------!'
     */' ! observ/trig [format:character(a6)]                       !'
     */' !----------------------------------------------------------!'
     */' !   stp/stop ..... stop  (for observ)                      !'
     */' !   xxxptl ....... particle observ. (type hlpptl for help) !'
     */' !   xxxevt ....... event observabl. (type hlpevt for help) !'
     */' !   parxxx ....... parameters       (type hlppar for help) !'
     */' !   trjptl ....... ptl trajectories (xobsrv, yobsrv=xxxptl)!'
     */' !----------------------------------------------------------!')
108   format(
     * ' !----------------------------------------------------------!'
     */' !   numptl ... number of ptls  (=1)                        !'
     */' !   nufptl ... number of ptls in forw hemisph              !'
     */' !   nubptl ... number of ptls in backw hemisph             !'
     */' !   nldptl ... number of lambdas                           !'
     */' !   nalptl ... number of alambdas                          !'
     */' !   nksptl ... number of k_shorts                          !'
     */' !   nngptl ... number of negatives                         !'
     */' !   npmptl ... number of ptls * 0+1-1 (zero/pos/neg chrg)  !'
     */' !   npaptl ... number of ptls * +1-1 (ptl/antiptl)         !'
     */' !   r.cptl ... hadronization point in pp cm (x,y,z,0)      !'
     */' !   yinptl ... space-time rapidity at hadr. point          !'
     */' !   tinptl ... proper time at hadr. point                  !'
     */' !   tmeptl ... time                                        !'
     */' !   tauptl ... eigentime                                   !'
     */' !   r.tptl ... space time (x,y,z,0)                        !'
     */' !   rttptl ... transverse space time                       !'
     */' !   xcmptl ... longitud. momentum fraction of ptl in pp cm !'
     */' !   xenptl ... energy fraction of ptl in pp cm             !'
     */' !   xelptl ... energy fraction of ptl in lab sys           !'
     */' !   xpcptl ... momentum fraction in cm                     !'
     */' !   xpaptl ... momentum fraction absolute in cm            !'
     */' !   xpeptl ... z-momentum (cm) / engy/2                    !'
     */' !   sobptl ... s / beta of ptl in pp cm                    !')
107   format(
     * ' !   xlbptl ... longitud. momentum fraction of ptl in lab   !'
     */' !   ycmptl ... rapdity of ptl in pp cm                     !'
     */' !   ycnptl ... -rapdity of ptl in pp cm                    !'
     */' !   yacptl ... /rapdity/ of ptl in pp cm                   !'
     */' !   yrcptl ... space time rapdity of ptl in pp cm          !'
     */' !   ylbptl ... rapdity of ptl in lab                       !'
     */' !   gcmptl ... pseudo rapdity of ptl in pp cm              !'
     */' !   glbptl ... pseudo rapdity of ptl in lab                !'
     */' !   ptrptl ... transverse momentum of ptl                  !'
     */' !   plbptl ... momentum of ptl (lab)                       !'
     */' !   etlptl ... transverse energy of ptl                    !'
     */' !   pt2ptl ... transverse momentum squared of ptl          !'
     */' !   amtptl ... transverse mass of ptl                      !'
     */' !   amsptl ... mass of ptl                                 !'
     */' !   amlptl ... log10(mass) of ptl                          !'
     */' !   phlptl ... azimuthal angle (lab)                       !'
     */' !   alsptl ... log s of evt                                !'
     */' !   nqjptl ... quark numbers of jets                       !'
     */' !   iorptl ... origin of ptl                               !'
     */' !----------------------------------------------------------!')
109   format(
     * ' !----------------------------------------------------------!'
     */' !   numevt ... number of evts  (=1)                        !'
     */' !   xlbevt ... momentum fraction of evt in pp cm           !'
     */' !   enlevt ... total energy of evt                         !'
     */' !   ev4evt ... zero degree energy of evt, eta>4            !'
     */' !   ev6evt ... zero degree energy of evt, eta>6            !'
     */' !   eplevt ... total eff energy of evt + proj.spect.energy !'
     */' !   etoevt ... total energy of evt + spect.energy          !'
     */' !   alsevt ... log s                                       !'
     */' !   etlevt ... transverve energy of evt                    !'
     */' !   ezlevt ... neutral energy of evt                       !'
     */' !   mulevt ... multiplicity of evt                         !'
     */' !   mufevt ... forw multiplicity of evt (rin1<y<rin2)      !'
     */' !   mubevt ... backw multiplicity of evt (-rin2<y<-rin1)   !'
     */' !   mbfevt ... mubevt*mufevt                               !'
     */' !   mldevt ... lambda multiplicity of evt                  !'
     */' !   malevt ... alambda multiplicity of evt                 !'
     */' !   mksevt ... k_short multiplicity of evt                 !'
     */' !   mngevt ... neg ptl multiplicity of evt (full phas sp)  !'
     */' !   rlnevt ... lambda/neg multiplicity of evt              !'
     */' !   ranevt ... alambda/neg multiplicity of evt             !'
     */' !   rknevt ... k_short/neg multiplicity of evt             !'
     */' !   bimevt ... impact parameter                            !'
     */' !   kolevt ... number of collisions                        !'
     */' !   colevt ... effective number of collisions              !'
     */' !   egyevt ... string energy (for choic=str)               !'
     */' !   npjevt ... number of involved projectile nucleons      !'
     */' !   ntgevt ... number of involved target nucleons          !'
     */' !----------------------------------------------------------!')
120   format(
     * ' !----------------------------------------------------------!'
     */' !   parp12 ... parameter par1,par2                         !'
     */' !   parp34 ... parameter par3,par4                         !'
     */' !   parpbi ... individual (pp) probability                 !'
     */' !   parrin ... rapidity interval                           !'
     */' !   parish ... ish                                         !'
     */' !----------------------------------------------------------!')
      par1=0
      par2=0
      par3=0
      par4=0
      ipbi=0
      parpbi=0
      rin1=0
      rin2=1e+10
      if(ntrig.gt.0)then
      if(ntrig.gt.mxtrig)then
      write(ifmt,*)'too many triggers'
      goto3
      endif
      do 4 nt=1,ntrig
5     continue
      if(modana.ne.'sup')write(ifmt,111)nt
111   format(7x,'trigger',i2,' ? [a3;hlp]    interval? [2f]')
      read(ifop,fmt=*,err=5)trig(nt),trig1(nt),trig2(nt)
      do 92 i=1,6
92    call utlow(trig(nt)(i:i))
      if(ivi.eq.1.and.modana.ne.'sup')
     *write(ifmt,1031)trig(nt),trig1(nt),trig2(nt)
1031  format(7x,a6,f7.2,f7.2)
      if(trig(nt).eq.'hlp')write(ifmt,101)
      if(trig(nt).eq.'hlp')goto5
      trg=trig(nt)
      if(trg(1:3).eq.'par')then
      if(trg.eq.'parp12')par1=trig1(nt)
      if(trg.eq.'parp12')par2=trig2(nt)
      if(trg.eq.'parp34')par3=trig1(nt)
      if(trg.eq.'parp34')par4=trig2(nt)
      if(trg.eq.'parpbi')ipbi=1
      if(trg.eq.'parpbi')parpbi=trig1(nt)
      if(trg.eq.'parrin')rin1=trig1(nt)
      if(trg.eq.'parrin')rin2=trig2(nt)
      if(trg.eq.'parish')ish=nint(trig1(nt))
      endif
4     continue
      endif
      if(nidp.gt.0)then
      if(nidp.gt.mxidp)then
      write(ifmt,*)'too many ptcl ids'
      goto3
      endif
      do 7 ni=1,nidp
2     continue
      if(modana.ne.'sup')write(ifmt,112)ni
112   format(7x,'particle id',i2,' ? [i]')
      read(ifop,fmt=*,err=2)idptl(ni)
      if(ivi.eq.1.and.modana.ne.'sup')write(ifmt,1032)idptl(ni)
1032  format(7x,i5)
7     continue
      endif
      if(ish.eq.20)ifhi=ifch
 
8     continue
      if(modana.ne.'sup')
     *write(ifmt,*)'interv?[2f]   #steps?[i]   mode?[a3]   axes?[6f]'
      read(ifop,fmt=*,err=8)x1,x2,kmax,mde,z1,z2,zi,y1,y2,yi
      do 91 i=1,3
91    call utlow(mde(i:i))
      if(ivi.eq.1.and.modana.ne.'sup')
     *write(ifmt,1033)x1,x2,kmax,mde,z1,z2,zi,y1,y2,yi
1033  format(1x,2f7.2,i4,3x,a3,6f7.2)
      if(kmax.lt.0)then
      write(ifmt,*)'kmax negative'
      goto8
      elseif(kmax.gt.kmaxmx)then
      write(ifmt,*)'kmax too large'
      goto8
      elseif(kmax.gt.0)then
      xinc=(x2-x1)/kmax
      endif
 
      if(mde.eq.'sho')then
      write(ifmt,*)'which event do you want to see? [i]'
      read(ifop,*)nevtsh
      endif
 
      if(mde(1:1).eq.'p')then
9     continue
      if(modana.ne.'sup')write(ifmt,*)'opt?[i,<0]  fct?[f]  '
     *,'line?[a6]  opt/dsp?[a13]  shft?[f]  #txt?[i]  expmt?[a5;i]'
      read(ifop,fmt=*,err=9)
     *ioptn,crofac,ili,muopt,shift,ntxt,expmt,idpexp
      do 93 i=1,6
93    call utlow(ili(i:i))
      do 94 i=1,13
94    call utlow(muopt(i:i))
      do 89 i=1,5
89    call utlow(expmt(i:i))
      if(ivi.eq.1.and.modana.ne.'sup')write(ifmt,1034)
     *ioptn,crofac,ili,muopt,shift,ntxt,expmt,idpexp
1034  format(1x,i3,f6.2,2x,a6,2x,a13,f5.1,i3,2x,a5,i9)
      if(ioptn.lt.0)write(ifmt,105)
      if(ioptn.lt.0)write(ifmt,115)
      if(ioptn.lt.0)goto9
105   format(
     * ' !----------------------------------------------------------!'
     */' ! option [format:integer]                                  !'
     */' !----------------------------------------------------------!'
     */' !   __0   *1                                               !'
     */' !   __1   / # of events                                    !'
     */' !   __2   / # of valid events                              !'
     */' !   __3   / # of counts (normalized to one)                !'
     */' !   __4   / # of counts per bin                            !'
     */' !   __5   / # of all events (incl trigger exclusions)      !'
     */' !   __6   / integral(y) (normalized to one)                !'
     */' !----------------------------------------------------------!')
115   format(
     * ' !   _0_   *1                                               !'
     */' !   _1_   / bin size                                       !'
     */' !   _2_   / bin size * sigppi(prt-ncl)                     !'
     */' !----------------------------------------------------------!'
     */' !   0__   *1                                               !'
     */' !   1__   * x                                              !'
     */' !   2__   / (2pi*x)                                        !'
     */' !   3__   * <x> ,  x/<x>  (kno-representation)             !'
     */' !   4__   scatter plot                                     !'
     */' !----------------------------------------------------------!')
      ioptn1=mod(ioptn,10)
      ioptn2=mod(ioptn/10,10)
      ioptn3=mod(ioptn/100,10)
      if(ipbi.eq.1)then
      if(ioptn1.eq.1.or.ioptn1.eq.5)stop'ipbi=1 => ioptn1 .ne. 1,5 !'
      endif
      if(crofac.lt.0.)then
      pltfac=-crofac
      crofac=1.
      endif
      if(expmt(1:2).eq.'fi')then
      if(idpexp.gt.mxpar)then
      write(ifmt,*)'too many parameters'
      goto9
      endif
      read(ifop,*)(parfi(i),i=1,idpexp)
      endif
      if(ntxt.gt.0)then
      if(ntxt.gt.mxtxt)then
      write(ifmt,*)'too many text inputs'
      goto9
      endif
      do 12 n=1,ntxt
      read(ifop,*)xt(n),yt(n),txt(n)
      if(ivi.eq.1.and.modana.ne.'sup')
     *write(ifmt,1007)xt(n),yt(n),txt(n)
12    continue
      endif
      endif
 
c     read from ifdt and analyse
c     --------------------------

      if(xobsrv.eq.'noanly'.or.yobsrv.eq.'noanly')goto95
      if(xobsrv(4:6).eq.'his'.or.yobsrv(4:6).eq.'his')goto95
      if(ish.ge.91
     *.or.ish.eq.15.or.ish.eq.20)write(ifch,106)
106   format(/' veanly: read from ifdt'/)
      rewind ifdt
      nevt=0
      nptl=0
      zevt0=0
      nscat=0
      do 44 i=1,5
44    facmo(i)=0
96    nevt=nevt+1
      do 29 k=1,kmax
      nacevt(k)=0
29    yacevt(k)=0
      read(ifdt,end=98)levts,nevts,nptevt,bimevt,kolevt,colevt,pmxevt
     *,egyevt,npjevt,ntgevt
      if(ish.ge.91)write(ifch,103)nevts,nptevt,bimevt
     *,kolevt,npjevt,ntgevt
103   format(' evt:',i2,i4,f9.2,3i5)
      if(levts.ne.levt)stop'***** label mismatch during read evt'
      if(nevts.ne.nevt)stop'***** evt# mismatch during read evt'

      pxlevt=0
      pylevt=0
      pzlevt=0
      enlevt=0
      ev4evt=0
      ev6evt=0
      eplevt=0
      etlevt=0
      ezlevt=0
      mulevt=0
      mufevt=0
      mubevt=0
      mldevt=0
      malevt=0
      mksevt=0
      mngevt=0
      rlnevt=0
      ranevt=0
      rknevt=0

      nptl=0
      do 97 npts=1,nptevt
      nptl=nptl+1
      read(ifdt)lptls,nevts,nptls,nptlo,idp
     *,pxlptl,pylptl,pzlptl,p0lptl,amsptl,iorptl,jorptl,istptl
     *,xo1ptl,xo2ptl,xo3ptl,xo4ptl,t1ptl,t2ptl,nqjptl
      if(lptls.ne.lptl)stop'***** label mismatch during read ptl'
      if(nevts.ne.nevt)stop'***** evt# mismatch during read ptl'
      if(nptls.ne.nptl)stop'***** ptl# mismatch during read ptl'
      numptl=1
      nldptl=0
      if(idp.eq.2130)nldptl=1
      nalptl=0
      if(idp.eq.-2130)nalptl=1
      nksptl=0
      if(idp.eq.20)nksptl=1
      nngptl=0
      if(idp.eq.-120.or.idp.eq.-130.or.idp.eq.-1120.or.idp.eq.-1130
     *.or.idp.eq.2230.or.idp.eq.2330.or.idp.eq.3331)nngptl=1
      chrg=0
      if(iabs(idp).le.9999.and.mod(iabs(idp),10).le.1)
     *call idchrg(idp,chrg)
      if(chrg.eq.0.)then
      npmptl=0
      else
      npmptl=sign(1.,chrg)
      endif
      npaptl=isign(1,idp)
      if(pxlptl.eq.0.)phlptl=pi/2
      if(pxlptl.ne.0.)phlptl=abs(atan(pylptl/pxlptl))
      if(pxlptl.lt.0..and.pylptl.ge.0.)phlptl=pi-phlptl
      if(pxlptl.lt.0..and.pylptl.lt.0.)phlptl=pi+phlptl
      if(pxlptl.ge.0..and.pylptl.lt.0.)phlptl=2*pi-phlptl
      phlptl=phlptl/pi*180
      plbptl=sqrt(pxlptl**2+pylptl**2+pzlptl**2)
      ptrptl=sqrt(pxlptl**2+pylptl**2)
      pt2ptl=ptrptl**2
      if(ptrptl.le.1e-6)then
      glbptl=sign(1.,pzlptl)*1e+6
      else
      glbptl=sign(1.,pzlptl)*
     *alog((sqrt(pzlptl**2+ptrptl**2)+abs(pzlptl))/ptrptl)
      endif
      amlptl=-1
      if(amsptl.gt.0.)amlptl=alog10(amsptl)
      amtptl=sqrt(amsptl**2+ptrptl**2)
      ylbptl=sign(1.,pzlptl)*alog((p0lptl+abs(pzlptl))/amtptl)
      ycmptl=ylbptl-yhaha*labsys
      ycnptl=-ycmptl
      yacptl=abs(ycmptl)
      pc1ptl=pxlptl
      pc2ptl=pylptl
      pc3ptl=amtptl*sinh(ycmptl)
      pc4ptl=amtptl*cosh(ycmptl)
      if(ptrptl.le.1e-6)then
      gcmptl=sign(1.,pc3ptl)*1e+6
      else
      gcmptl=sign(1.,pc3ptl)*
     *alog((sqrt(pc3ptl**2+ptrptl**2)+abs(pc3ptl))/ptrptl)
      endif
      xlbptl=pzlptl/pmxevt
      xcmptl=xlbptl
      if(ichoic.ne.1.and.ichoic.ne.4)xcmptl=0
      if(ichoic.ne.1.and.ichoic.ne.4.and.pnllx.ne.0.)xcmptl=pc3ptl/pnllx
      xenptl=2*pc4ptl/engy
      xelptl=xenptl
      if(ichoic.ne.1.and.ichoic.ne.4)xelptl=0.
      if(ichoic.ne.1.and.ichoic.ne.4.and.pnll.ne.0.)xelptl=p0lptl/pnll
      xpcptl=2*sqrt(pc3ptl**2+ptrptl**2)/engy
      xpaptl=abs(xpcptl)
      xpeptl=2*pc3ptl/engy
      sobptl=ainfin
      if(sqrt(pc3ptl**2+pc2ptl**2+pc1ptl**2).ne.0)
     *sobptl=engy**2/sqrt(pc3ptl**2+pc2ptl**2+pc1ptl**2)*pc4ptl
      eefptl=p0lptl
      if(idp.ge.1000)eefptl=eefptl-prom
      if(idp.le.-1000)eefptl=eefptl+prom
      if(glbptl.gt.4.)ev4evt=ev4evt+eefptl
      if(glbptl.gt.6.)ev6evt=ev6evt+eefptl
      etlptl=0
      if(ptrptl.ne.0)
     *etlptl=eefptl*ptrptl/sqrt(pzlptl**2+ptrptl**2)
      vel=sqrt(pzlptl**2+ptrptl**2)/p0lptl
      if(idp.eq.10)ezlptl=p0lptl
      if(abs(chrg).gt.0.99.and.vel.gt.0.80)ezlptl=.500
      alsptl=alog10(engy**2)
      nufptl=0
      nubptl=0
      if(idp.eq.120 .or.idp.eq.-120.or.idp.eq.130 .or.idp.eq.-130
     *.or.idp.eq.1120.or.idp.eq.-1120.or.idp.eq.1130.or.idp.eq.-1130
     *.or.idp.eq.2230.or.idp.eq.-2230.or.idp.eq.2330.or.idp.eq.-2330
     *.or.idp.eq.3331.or.idp.eq.-3331)then
      if(ycmptl.ge.rin1.and.ycmptl.le.rin2)nufptl=1
      if(ycmptl.ge.-rin2.and.ycmptl.le.-rin1)nubptl=1
      endif
      mufevt=mufevt+nufptl
      mubevt=mubevt+nubptl
      mbfevt=mubevt*mufevt

           if(nidp.gt.0)then
      if(idptl(1).eq.9990)then
      elseif(idptl(1).eq.9980)then
      if(idp.ne.20 .and.idp.ne.-20.and.iabs(idp).lt.100)goto20
      elseif(idptl(1).eq.9970)then
      if(idp.ne.120 .and.idp.ne.-120.and.idp.ne.130 .and.idp.ne.-130
     *.and.idp.ne.1120.and.idp.ne.-1120.and.idp.ne.1130.and.idp.ne.-1130
     *.and.idp.ne.2230.and.idp.ne.-2230.and.idp.ne.2330.and.idp.ne.-2330
     *.and.idp.ne.3331.and.idp.ne.-3331)goto20
      elseif(idptl(1).eq.9960)then
      if(idp.ne.120 .and.idp.ne.130 .and.idp.ne.1120.and.idp.ne.1130
     *.and.idp.ne.-2230.and.idp.ne.-2330.and.idp.ne.-3331)goto20
      elseif(idptl(1).eq.-9960)then
      if(idp.ne.-120.and.idp.ne.-130.and.idp.ne.-1120.and.idp.ne.-1130
     *.and.idp.ne.2230.and.idp.ne.2330.and.idp.ne.3331)goto20
      elseif(idptl(1).eq.9950)then
      if(istptl.ne.0)goto20
      if(iabs(idp).lt.1000)goto20
      elseif(idptl(1).eq.9940)then
      if(istptl.ne.0)goto20
      if(idp.lt.1000)goto20
      elseif(idptl(1).eq.-9940)then
      if(istptl.ne.0)goto20
      if(idp.gt.-1000)goto20
      elseif(idptl(1).eq.9930)then
      if(istptl.ne.0)goto20
      if(idp.ne.20 .and.idp.ne.-20.and.iabs(idp).lt.100)goto20
      if(iabs(idp).gt.1000)goto20
      else
      do 39 n=1,nidp
      if(idp.eq.idptl(n))goto40
39    continue
      goto20
40    continue
      endif
           endif

           if(ntrig.gt.0)then
      do 22 nt=1,ntrig
      trg=trig(nt)
      if(trg(4:6).ne.'ptl')goto22
      tr=0.5*(trig1(nt)+trig2(nt))
      if(trig(nt).eq.'numptl')tr=numptl
      if(trig(nt).eq.'nufptl')tr=nufptl
      if(trig(nt).eq.'nubptl')tr=nubptl
      if(trig(nt).eq.'npmptl')tr=npmptl
      if(trig(nt).eq.'npaptl')tr=npaptl
      if(trig(nt).eq.'xcmptl')tr=xcmptl
      if(trig(nt).eq.'xenptl')tr=xenptl
      if(trig(nt).eq.'xelptl')tr=xelptl
      if(trig(nt).eq.'xpcptl')tr=xpcptl
      if(trig(nt).eq.'xpaptl')tr=xpaptl
      if(trig(nt).eq.'xpeptl')tr=xpeptl
      if(trig(nt).eq.'sobptl')tr=sobptl
      if(trig(nt).eq.'xlbptl')tr=xlbptl
      if(trig(nt).eq.'ycmptl')tr=ycmptl
      if(trig(nt).eq.'ycnptl')tr=ycnptl
      if(trig(nt).eq.'yacptl')tr=yacptl
      if(trig(nt).eq.'yrcptl')tr=yrcptl
      if(trig(nt).eq.'ylbptl')tr=ylbptl
      if(trig(nt).eq.'gcmptl')tr=gcmptl
      if(trig(nt).eq.'glbptl')tr=glbptl
      if(trig(nt).eq.'ptrptl')tr=ptrptl
      if(trig(nt).eq.'plbptl')tr=plbptl
      if(trig(nt).eq.'etlptl')tr=etlptl
      if(trig(nt).eq.'pt2ptl')tr=pt2ptl
      if(trig(nt).eq.'amtptl')tr=amtptl
      if(trig(nt).eq.'amsptl')tr=amsptl
      if(trig(nt).eq.'amlptl')tr=amlptl
      if(trig(nt).eq.'phlptl')tr=phlptl
      if(trig(nt).eq.'alsptl')tr=alsptl
      if(trig(nt).eq.'nqjptl')tr=nqjptl
      if(trig(nt).eq.'iorptl')tr=iorptl
      if(tr.lt.trig1(nt).or.tr.gt.trig2(nt))goto20
22    continue
           endif

      if(xobsrv.eq.'numptl')x=numptl
      if(xobsrv.eq.'nufptl')x=nufptl
      if(xobsrv.eq.'nubptl')x=nubptl
      if(xobsrv.eq.'npmptl')x=npmptl
      if(xobsrv.eq.'npaptl')x=npaptl
      if(xobsrv.eq.'xcmptl')x=xcmptl
      if(xobsrv.eq.'xenptl')x=xenptl
      if(xobsrv.eq.'xelptl')x=xelptl
      if(xobsrv.eq.'xpcptl')x=xpcptl
      if(xobsrv.eq.'xpaptl')x=xpaptl
      if(xobsrv.eq.'xpeptl')x=xpeptl
      if(xobsrv.eq.'sobptl')x=sobptl
      if(xobsrv.eq.'xlbptl')x=xlbptl
      if(xobsrv.eq.'ycmptl')x=ycmptl
      if(xobsrv.eq.'ycnptl')x=ycnptl
      if(xobsrv.eq.'yacptl')x=yacptl
      if(xobsrv.eq.'yrcptl')x=yrcptl
      if(xobsrv.eq.'ylbptl')x=ylbptl
      if(xobsrv.eq.'gcmptl')x=gcmptl
      if(xobsrv.eq.'glbptl')x=glbptl
      if(xobsrv.eq.'ptrptl')x=ptrptl
      if(xobsrv.eq.'plbptl')x=plbptl
      if(xobsrv.eq.'etlptl')x=etlptl
      if(xobsrv.eq.'pt2ptl')x=pt2ptl
      if(xobsrv.eq.'amtptl')x=amtptl
      if(xobsrv.eq.'amsptl')x=amsptl
      if(xobsrv.eq.'amlptl')x=amlptl
      if(xobsrv.eq.'phlptl')x=phlptl
      if(xobsrv.eq.'alsptl')x=alsptl
      if(xobsrv.eq.'nqjptl')x=nqjptl
      if(xobsrv.eq.'iorptl')x=iorptl
      if(yobsrv.eq.'numptl')y=numptl
      if(yobsrv.eq.'nufptl')y=nufptl
      if(yobsrv.eq.'nubptl')y=nubptl
      if(yobsrv.eq.'npmptl')y=npmptl
      if(yobsrv.eq.'npaptl')y=npaptl
      if(yobsrv.eq.'xcmptl')y=xcmptl
      if(yobsrv.eq.'xenptl')y=xenptl
      if(yobsrv.eq.'xelptl')y=xelptl
      if(yobsrv.eq.'xpcptl')y=xpcptl
      if(yobsrv.eq.'xpaptl')y=xpaptl
      if(yobsrv.eq.'xpeptl')y=xpeptl
      if(yobsrv.eq.'sobptl')y=sobptl
      if(yobsrv.eq.'xlbptl')y=xlbptl
      if(yobsrv.eq.'ycmptl')y=ycmptl
      if(yobsrv.eq.'ycnptl')y=ycnptl
      if(yobsrv.eq.'yacptl')y=yacptl
      if(yobsrv.eq.'yrcptl')y=yrcptl
      if(yobsrv.eq.'ylbptl')y=ylbptl
      if(yobsrv.eq.'gcmptl')y=gcmptl
      if(yobsrv.eq.'glbptl')y=glbptl
      if(yobsrv.eq.'ptrptl')y=ptrptl
      if(yobsrv.eq.'plbptl')y=plbptl
      if(yobsrv.eq.'etlptl')y=etlptl
      if(yobsrv.eq.'pt2ptl')y=pt2ptl
      if(yobsrv.eq.'amtptl')y=amtptl
      if(yobsrv.eq.'amsptl')y=amsptl
      if(yobsrv.eq.'amlptl')y=amlptl
      if(yobsrv.eq.'phlptl')y=phlptl
      if(yobsrv.eq.'alsptl')y=alsptl
      if(yobsrv.eq.'nqjptl')y=nqjptl
      if(yobsrv.eq.'iorptl')y=iorptl
      pxlevt=pxlevt+pxlptl
      pylevt=pylevt+pylptl
      pzlevt=pzlevt+pzlptl
      enlevt=enlevt+p0lptl
      eplevt=eplevt+eefptl
      etlevt=etlevt+etlptl
      ezlevt=ezlevt+ezlptl
      mulevt=mulevt+numptl
      mldevt=mldevt+nldptl
      malevt=malevt+nalptl
      mksevt=mksevt+nksptl
      itrptl=1
      if(ish.ge.91
     *.or.ish.eq.15)write(ifch,101)
     *nevt,nptlo,idp,ycmptl,amsptl,istptl,iorptl
101   format(1x,i5,i5,i11,f7.1,2x,e8.2,i4,i10)
      mscat=min(mxscat,nint(mxscat*crofac))
      if(xobsrv(4:6).eq.'ptl'.and.ioptn3.eq.4.and.nscat.lt.mscat)then
      if(xobsrv.eq.'trjptl')then
      nscat=nscat+1
      t1=amin1(t1ptl,y2)
      t1=amax1(t1,y1)
      arscat(nscat,1)=xo3ptl+pc3ptl/pc4ptl*(t1-xo4ptl)
      arscat(nscat,2)=t1
      nscat=nscat+1
      t2=amin1(t2ptl,y2)
      t2=amax1(t2,y1)
      arscat(nscat,1)=xo3ptl+pc3ptl/pc4ptl*(t2-xo4ptl)
      arscat(nscat,2)=t2
      else
      nscat=nscat+1
      arscat(nscat,1)=x
      arscat(nscat,2)=y
      endif
      endif
      goto21
20    itrptl=0
21    continue
      mngevt=mngevt+nngptl
 
c     last ptl in event
c     -----------------
      lptevt=0
      if(npts.eq.nptevt)then
      lptevt=1
      numevt=1
      if(mngevt.gt.0.)rlnevt=mldevt/real(mngevt)
      if(mngevt.gt.0.)ranevt=malevt/real(mngevt)
      if(mngevt.gt.0.)rknevt=mksevt/real(mngevt)
      enlevt=enlevt/1000.
      eplevt=eplevt/1000.
      ev4evt=ev4evt/1000.
      ev6evt=ev6evt/1000.
      alsevt=alog10(engy**2)
      if(ntrig.gt.0)then
      do 25 nt=1,ntrig
      trg=trig(nt)
      if(trg(4:6).ne.'evt')goto25
      tr=0.5*(trig1(nt)+trig2(nt))
      if(trig(nt).eq.'numevt')tr=numevt
      if(trig(nt).eq.'xlbevt')tr=pzlevt/pmxevt
      if(trig(nt).eq.'enlevt')tr=enlevt
      if(trig(nt).eq.'eplevt')tr=eplevt
      if(trig(nt).eq.'ev4evt')tr=ev4evt
      if(trig(nt).eq.'ev6evt')tr=ev6evt
      if(trig(nt).eq.'etoevt')tr=etoevt
      if(trig(nt).eq.'alsevt')tr=alsevt
      if(trig(nt).eq.'etlevt')tr=etlevt
      if(trig(nt).eq.'ezlevt')tr=ezlevt
      if(trig(nt).eq.'mulevt')tr=mulevt
      if(trig(nt).eq.'mufevt')tr=mufevt
      if(trig(nt).eq.'mubevt')tr=mubevt
      if(trig(nt).eq.'mbfevt')tr=mbfevt
      if(trig(nt).eq.'mldevt')tr=mldevt
      if(trig(nt).eq.'malevt')tr=malevt
      if(trig(nt).eq.'mksevt')tr=mksevt
      if(trig(nt).eq.'mngevt')tr=mngevt
      if(trig(nt).eq.'rlnevt')tr=rlnevt
      if(trig(nt).eq.'ranevt')tr=ranevt
      if(trig(nt).eq.'rknevt')tr=rknevt
      if(trig(nt).eq.'bimevt')tr=bimevt
      if(trig(nt).eq.'kolevt')tr=kolevt
      if(trig(nt).eq.'colevt')tr=colevt
      if(trig(nt).eq.'egyevt')tr=egyevt
      if(trig(nt).eq.'npjevt')tr=npjevt
      if(trig(nt).eq.'ntgevt')tr=ntgevt
      if(tr.lt.trig1(nt).or.tr.gt.trig2(nt))goto26
25    continue
      endif
      if(xobsrv.eq.'numevt')x=numevt
      if(xobsrv.eq.'xlbevt')x=pzlevt/pmxevt
      if(xobsrv.eq.'enlevt')x=enlevt
      if(xobsrv.eq.'eplevt')x=eplevt
      if(xobsrv.eq.'ev4evt')x=ev4evt
      if(xobsrv.eq.'ev6evt')x=ev6evt
      if(xobsrv.eq.'etoevt')x=etoevt
      if(xobsrv.eq.'alsevt')x=alsevt
      if(xobsrv.eq.'etlevt')x=etlevt
      if(xobsrv.eq.'ezlevt')x=ezlevt
      if(xobsrv.eq.'mulevt')x=mulevt
      if(xobsrv.eq.'mufevt')x=mufevt
      if(xobsrv.eq.'mubevt')x=mubevt
      if(xobsrv.eq.'mbfevt')x=mbfevt
      if(xobsrv.eq.'mldevt')x=mldevt
      if(xobsrv.eq.'malevt')x=malevt
      if(xobsrv.eq.'mksevt')x=mksevt
      if(xobsrv.eq.'mngevt')x=mngevt
      if(xobsrv.eq.'rlnevt')x=rlnevt
      if(xobsrv.eq.'ranevt')x=ranevt
      if(xobsrv.eq.'rknevt')x=rknevt
      if(xobsrv.eq.'bimevt')x=bimevt
      if(xobsrv.eq.'kolevt')x=kolevt
      if(xobsrv.eq.'colevt')x=colevt
      if(xobsrv.eq.'egyevt')x=egyevt
      if(xobsrv.eq.'npjevt')x=npjevt
      if(xobsrv.eq.'ntgevt')x=ntgevt
      if(yobsrv.eq.'numevt')y=numevt
      if(yobsrv.eq.'xlbevt')y=pzlevt/pmxevt
      if(yobsrv.eq.'enlevt')y=enlevt
      if(yobsrv.eq.'eplevt')y=eplevt
      if(yobsrv.eq.'ev4evt')y=ev4evt
      if(yobsrv.eq.'ev6evt')y=ev6evt
      if(yobsrv.eq.'etoevt')y=etoevt
      if(yobsrv.eq.'alsevt')y=alsevt
      if(yobsrv.eq.'etlevt')y=etlevt
      if(yobsrv.eq.'ezlevt')y=ezlevt
      if(yobsrv.eq.'mulevt')y=mulevt
      if(yobsrv.eq.'mufevt')y=mufevt
      if(yobsrv.eq.'mubevt')y=mubevt
      if(yobsrv.eq.'mbfevt')y=mbfevt
      if(yobsrv.eq.'mldevt')y=mldevt
      if(yobsrv.eq.'malevt')y=malevt
      if(yobsrv.eq.'mksevt')y=mksevt
      if(yobsrv.eq.'mngevt')y=mngevt
      if(yobsrv.eq.'rlnevt')y=rlnevt
      if(yobsrv.eq.'ranevt')y=ranevt
      if(yobsrv.eq.'rknevt')y=rknevt
      if(yobsrv.eq.'bimevt')y=bimevt
      if(yobsrv.eq.'kolevt')y=kolevt
      if(yobsrv.eq.'colevt')y=colevt
      if(yobsrv.eq.'egyevt')y=egyevt
      if(yobsrv.eq.'npjevt')y=npjevt
      if(yobsrv.eq.'ntgevt')y=ntgevt
      if(xobsrv(4:6).eq.'evt'.and.ioptn3.eq.4.and.nscat.lt.mxscat)then
      nscat=nscat+1
      arscat(nscat,1)=x
      arscat(nscat,2)=y
      endif
      if(ish.eq.20)then
      write(ifch,*)'evt:',nevt
      endif
      itrevt=1
      goto27
26    itrevt=0
27    continue
      endif
 
      if((xobsrv(4:6).eq.'ptl'.and.itrptl.eq.1
     *.or.xobsrv(4:6).eq.'evt'.and.itrevt.eq.1.and.lptevt.eq.1)
     *.and.ioptn3.ne.4)then
      if(x.ge.x1.and.x.le.x2)then
      if(ioptn3.eq.1)y=y*x
      if(ioptn3.eq.2)y=y/(2*pi*x)
      endif
      k=1+(x-x1)/xinc
      if(k.ge.1.and.k.le.kmax)nacevt(k)=nacevt(k)+1
      if(k.ge.1.and.k.le.kmax)yacevt(k)=yacevt(k)+y
      endif
      if(lptevt.eq.1.and.itrevt.eq.1)then
      if(ipbi.ne.1)wei=1
      if(ipbi.eq.1)wei=parpbi*kolevt
      zevt0=zevt0+wei
      do 28 k=1,kmax
      zactot(k)=zactot(k)+nacevt(k)*wei
28    yactot(k)=yactot(k)+yacevt(k)*wei
      endif
 
97    continue
      goto96
98    continue
      nevt=nevt-1
 
      zcount=0
      if(ioptn3.eq.4)then
      first=0
      varnce=0
      goto95
      endif
      yacsum=0
      do 30 k=1,kmax
      ar(k,1)=x1+(k-.5)*xinc
      if(ioptn1.ne.4)ar(k,3)=yactot(k)
      if(ioptn1.eq.4.and.yactot(k).ne.0.)ar(k,3)=yactot(k)/zactot(k)
      ar(k,5)=zactot(k)
      yacsum=yacsum+yactot(k)
30    zcount=zcount+zactot(k)
      if(bmax.eq.0.)sgm=sigppi*10
      if(bmax.ne.0.)sgm=pi*bmax**2*nevt/ntevt *10
      if(ioptn1.eq.0)cnorm=1.
      if(ioptn1.eq.1)cnorm=1./nevt
      if(ioptn1.eq.2.and.zevt0.ne.0.)cnorm=1/zevt0
      if(ioptn1.eq.3.and.zcount.ne.0.)cnorm=1/zcount
      if(ioptn1.eq.4)cnorm=1.
      if(ioptn1.eq.5.and.naevt.ne.0)cnorm=1./naevt
      if(ioptn1.eq.5.and.naevt.eq.0)cnorm=0
      if(ioptn1.eq.6.and.yacsum.ne.0.)cnorm=1./yacsum
      if(ioptn2.eq.1)cnorm=cnorm/xinc
      if(ioptn2.eq.2)cnorm=cnorm/xinc*sgm
      first=0
      secnd=0
      do 31 k=1,kmax
      if(zcount.ne.0.)first=first+ar(k,1)*ar(k,5)/zcount
      if(zcount.ne.0.)secnd=secnd+ar(k,1)**2*ar(k,5)/zcount
31    continue
      varnce=secnd-first**2
      do 32 k=1,kmax
      ar(k,3)=ar(k,3)*cnorm*crofac
      if(ioptn3.eq.3)ar(k,3)=ar(k,3)*first
      if(ioptn3.eq.3)ar(k,1)=ar(k,1)/first
      if(ar(k,5).ne.0.)ar(k,4)=ar(k,3)/sqrt(ar(k,5))
32    continue
 
c     write head
c     ----------
95    write(ifhi,1013)fch,
     *(das,i=1,68),float(iversn)/100,engy,laproj,maproj,latarg,matarg
1013  format(1x,1h',69a1,1h'/1x,2h'h,4x,f4.2,4x,f9.3,4x,4i4,27x,1h')
      if(maproj.eq.0)then
      do 48 n=1,99
      if(prob(n).eq.0.)goto47
48    write(ifhi,1014)prob(n)*100
     *,icbac(n,1),icbac(n,2),icfor(n,1),icfor(n,2)
1014  format(1x,1h',5x,f6.2,4x,2i8,4x,2i8,18x,1h')
      endif
47    write(ifhi,1012)xobsrv,yobsrv,ntrig,nidp,npar
1012  format(1x,1h',5x,a6,4x,a6,4x,i2,4x,i2,4x,i2,30x,1h')
      nmx=max0(ntrig,nidp)
      nmx=max0(nmx,npar)
      if(nmx.gt.0)then
      do 34 n=1,nmx
      parv=parval(n)
      if(parval(n).ge.1000.)parv=parval(n)/1000
      if(n.le.ntrig.and.n.le.nidp.and.n.le.npar)then
      write(ifhi,1021)trig(n),trig1(n),trig2(n),idptl(n),parcha(n)
     *,parv
1021  format(1x,1h',5x,a6,1x,e10.4,1x,e10.4,4x,i5,4x,a6,f8.3,9x,1h')
      elseif(n.le.ntrig.and.n.le.nidp)then
      write(ifhi,1022)trig(n),trig1(n),trig2(n),idptl(n)
1022  format(1x,1h',5x,a6,1x,e10.4,1x,e10.4,4x,i5,4x,14x    ,9x,1h')
      elseif(n.le.ntrig.and.n.le.npar)then
      write(ifhi,1023)trig(n),trig1(n),trig2(n),parcha(n),parv
1023  format(1x,1h',5x,a6,1x,e10.4,1x,e10.4,4x,5x,4x,a6,f8.3,9x,1h')
      elseif(n.le.nidp.and.n.le.npar)then
      write(ifhi,1024)idptl(n),parcha(n),parv
1024  format(1x,1h',5x,28x         ,4x,i5,4x,a6,f8.3,9x,1h')
      elseif(n.le.ntrig)then
      write(ifhi,1025)trig(n),trig1(n),trig2(n)
1025  format(1x,1h',5x,a6,1x,e10.4,1x,e10.4,4x,5x,4x,14x    ,9x,1h')
      elseif(n.le.nidp)then
      write(ifhi,1026)idptl(n)
1026  format(1x,1h',5x,28x         ,4x,i5,4x,14x    ,9x,1h')
      elseif(n.le.npar)then
      write(ifhi,1027)parcha(n),parv
1027  format(1x,1h',5x,28x         ,4x,5x,4x,a6,f8.3,9x,1h')
      endif
34    continue
      endif
      write(ifhi,1003)1.*nrevt,zevt0,zcount,ioptn,crofac,first,varnce
1003  format(1x,2h'm,4x,e9.4,1x,e9.4,1x,e9.4,2x,i3,f8.2,2f9.3,4x,1h')
      if(modana.ne.'sup')write(ifmt,1002)1.*nevt,zevt0,zcount
     *,first,varnce
1002  format(1x,'=========>   ',e9.4,1x,e9.4,1x,e9.4,1x,f9.3,1x,f9.3
     *,'   <=========')
      fch='+'
      if(mde(1:1).ne.'p')goto99
      fch='-'
 
c     fill histogram
c     --------------
      kmx=kmax
      if(ioptn3.eq.4)kmx=nscat
      if(xobsrv.eq.'iwzzzz'
     *.or.yobsrv.eq.'iwzzzz')call wzzzzf(ioptn,kmx)
      if(xobsrv.eq.'imihis'
     *.or.yobsrv.eq.'imihis')call wimino(ioptn,kmx)
      if(xobsrv.eq.'isphis'
     *.or.yobsrv.eq.'isphis')call wspafi(ioptn,kmx)
      if(xobsrv.eq.'iclhis'
     *.or.yobsrv.eq.'iclhis')call wclufi(ioptn,kmx)
      if(xobsrv.eq.'ijphis'
     *.or.yobsrv.eq.'ijphis')call wjpsif(ioptn,kmx)
      write(ifhi,1004)(das,i=1,69)
     *,z1,z2,zi,y1,y2,yi,expmt,idpexp
     *,mde,ili,muopt,shift,ntxt,kmx,pltfac
1004  format(1x,1h',69a1,1h'/1x,6f8.2,2x,1h',a5,1h',i11
     */4x,1h',a3,1h',3x,1h',a6,1h',3x,1h',a13,1h',f8.2,2i6,f10.4)
      if(expmt(1:2).eq.'fi')then
      write(ifhi,1006)(das,i=1,69)
1006  format(1x,1h',69a1,1h')
      write(ifhi,1008)(parfi(i),i=1,idpexp)
1008  format(4x,5e13.5)
      endif
      if(ntxt.gt.0)then
      write(ifhi,1006)(das,i=1,69)
      do 11 n=1,ntxt
11    write(ifhi,1007)xt(n),yt(n),txt(n)
1007  format(4x,f4.2,2x,f4.2,2x,1h',a50,1h')
      endif
      if(kmx.gt.0)then
      write(ifhi,1006)(das,i=1,69)
      zer=0
      do 36 k=1,kmx
      if(ioptn3.eq.4
     *.and.xobsrv(4:6).ne.'his'.and.yobsrv(4:6).ne.'his')then
      write(ifhi,1005)k,arscat(k,1),zer,arscat(k,2),zer,zer
      else
      write(ifhi,1005)k,(ar(k,j),j=1,5)
      endif
36    continue
1005  format(4x,i5,2x,5e12.4)
      endif
      goto99
      end
 
c-----------------------------------------------------------------------
      subroutine wclu
c     analysis. cluster (resonance) formation and decay.
c     fills arclu(n,k); k:mass bin.
c     n=1,12: cluster formation (m=meson, b=baryon, ab=antibaryon
c       s=strange, as=antistrange, ns=nonstrange).
c       arclu contains the number of clusters for the foll. reactions:
c          n=1:   m+m
c          n=2:   m+b
c          n=3:   b+b
c          n=4:   m+ab
c          n=5:   b+ab
c          n=6:   p+p
c          n=7:   p+lda
c          n=8:   m/s+b/ns
c          n=9:   m/ns+b/s
c          n=10:  m/ns+ab/as
c          n=11:  m/as+ab/ns
c          n=12:  m/ns+ab/ns
c     n=13,60: cluster decay.
c          n=13,16:  # of all, k+, k-, lda   from m
c          n=17,20:  # of all, k+, k-, lda   from b
c          n=21,24:  # of all, k+, k-, lda   from bb (6 quarks)
c          n=25,36:  like 13-24, however nonstrange clusters required
c          n=37,60:  # of clusters corr to n=13,36 (for norm.)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer ifuptl(mxptl)
      parameter(nclu=60,kclu=25)
      common/cclu/arclu(nclu,kclu),aclu(kclu)
      real tar(nclu)
      parameter (ncomx=100)
      integer iar(ncomx)
      parameter (nflav=6)
      integer jcdu(nflav,2)
      amin=0
      amax=5
 
      if(nrevt.eq.1)then
      do 44 k=1,kclu
      aclu(k)=amin+(k-.5)*(amax-amin)/kclu
      do 44 n=1,nclu
44    arclu(n,k)=0
      endif
 
      do 1 i=1,nptl
      ifuptl(i)=0
      if(istptl(i).ne.1.or.iorptl(i).eq.-1)goto1
         if(ifrptl(1,i).ne.0.and.istptl(ifrptl(1,i)).ne.2)then
      if(ifrptl(2,i).lt.ifrptl(1,i))
     *call utstop('wclu: ifrptl2<ifrptl1                   ')
      en=pptl(4,i)
      do 2 j=ifrptl(1,i),ifrptl(2,i)
2     en=en-pptl(4,j)
      if(abs(en/pptl(4,i)).gt.1e-3)then
      if(ish.ge.90)then
      call utmsg('jclu  ')
      write(ifch,*)'*****  en_cluster.ne.en_fragments.'
      write(ifch,*)'en_cl, en_cl-en_fr:',pptl(4,i),en
      call utmsgf
      endif
      endif
         else
      ifuptl(i)=2
         endif
1     continue
 
      if(ish.eq.10)then
      write(ifch,101)('-',i=1,16)
101   format(/1x,'ptl list in wclu'/1x,16a1/)
      do 10 i=1,nptl
10    write(ifch,100)i,idptl(i),ifuptl(i)
     *,iorptl(i),jorptl(i),ifrptl(1,i),ifrptl(2,i),istptl(i)
100   format(1x,i4,i10,i6,7x,2i4,7x,2i4,i8)
      write(ifch,*)' '
      endif
 
      do 3 i=1,nptl
      if(istptl(i).eq.2)goto3
      io=iorptl(i)
      jo=jorptl(i)
      if(jo.eq.0)goto3
      ifuptl(io)=ifuptl(io)-2
      ifuptl(jo)=ifuptl(jo)-2
      if(abs(pptl(4,i)-pptl(4,io)-pptl(4,jo)).gt.1e-3)then
      ifuptl(io)=ifuptl(io)+1
      ifuptl(jo)=ifuptl(jo)+1
      goto3
      endif
      idi=idptl(io)
      idj=idptl(jo)
      call idquac(io,nqi,nsi,ndummy,jcdu)
      call idquac(jo,nqj,nsj,ndummy,jcdu)
      k=1+(pptl(5,i)-amin)/(amax-amin)*kclu
         if(k.ge.1.and.k.le.kclu)then
      if(nqi.eq.0.and.nqj.eq.0)       arclu( 1,k)=arclu( 1,k)+1
      if(nqi.eq.0.and.nqj.eq.3
     *.or.nqj.eq.0.and.nqi.eq.3)      arclu( 2,k)=arclu( 2,k)+1
      if(nqi.eq.3.and.nqj.eq.3)       arclu( 3,k)=arclu( 3,k)+1
      if(nqi.eq.0.and.nqj.eq.-3
     *.or.nqj.eq.0.and.nqi.eq.-3)     arclu( 4,k)=arclu( 4,k)+1
      if(nqi.eq.3.and.nqj.eq.-3
     *.or.nqj.eq.3.and.nqi.eq.-3)     arclu( 5,k)=arclu( 5,k)+1
      if(idi.eq.1120.and.idj.eq.1120) arclu( 6,k)=arclu( 6,k)+1
      if(idi.eq.1120.and.idj.eq.2130
     *.or.idj.eq.1120.and.idi.eq.2130)arclu( 7,k)=arclu( 7,k)+1
      if(nqi.eq.0.and.nsi.eq.1
     *.and.nqj.eq.3.and.nsj.eq.0
     *.or.nqj.eq.0.and.nsj.eq.1
     *.and.nqi.eq.3.and.nsi.eq.0)     arclu( 8,k)=arclu( 8,k)+1
      if(nqi.eq.0.and.nsi.eq.0
     *.and.nqj.eq.3.and.nsj.eq.1
     *.or.nqj.eq.0.and.nsj.eq.0
     *.and.nqi.eq.3.and.nsi.eq.1)     arclu( 9,k)=arclu( 9,k)+1
      if(nqi.eq.0.and.nsi.eq.0
     *.and.nqj.eq.-3.and.nsj.eq.-1
     *.or.nqj.eq.0.and.nsj.eq.0
     *.and.nqi.eq.-3.and.nsi.eq.-1)   arclu(10,k)=arclu(10,k)+1
      if(nqi.eq.0.and.nsi.eq.-1
     *.and.nqj.eq.-3.and.nsj.eq.0
     *.or.nqj.eq.0.and.nsj.eq.-1
     *.and.nqi.eq.-3.and.nsi.eq.0)    arclu(11,k)=arclu(11,k)+1
      if(nqi.eq.0.and.nsi.eq.0
     *.and.nqj.eq.-3.and.nsj.eq.0
     *.or.nqj.eq.0.and.nsj.eq.0
     *.and.nqi.eq.-3.and.nsi.eq.0)    arclu(12,k)=arclu(12,k)+1
         endif
3     continue
 
      do 4 i=1,nptl
      if(ifuptl(i).eq.0)goto4
      if(ish.ge.90
     *.or.ish.eq.10)then
      call utmsg('jclu  ')
      write(ifch,*)'*****  ifuptl .ne. 0.'
      write(ifch,*)'i,ifuptl(i): ',i,ifuptl(i)
      call utmsgf
      endif
4     continue
 
      do 5 i=1,nptl
      if(istptl(i).ne.1.or.iorptl(i).eq.-1)goto5
      if(.not.(ifrptl(1,i).ne.0.and.istptl(ifrptl(1,i)).ne.2))goto5
      k=1+(pptl(5,i)-amin)/(amax-amin)*kclu
      if(k.lt.1.or.k.gt.kclu)goto5
      call idquac(i,nq,ns,ndummy,jcdu)
      do 45 n=1,nclu
45    tar(n)=0
      iar(1)=i
      nco=1
      npo=0
7     npo=npo+1
      ipo=iar(npo)
         if(istptl(ipo).eq.2)goto5
         if(istptl(ipo).eq.0)then
      id=idptl(ipo)
      if(nq.eq.0)               tar(13)=tar(13)+1
      if(nq.eq.0.and.id.eq.130) tar(14)=tar(14)+1
      if(nq.eq.0.and.id.eq.-130)tar(15)=tar(15)+1
      if(nq.eq.0.and.id.eq.2130)tar(16)=tar(16)+1
      if(nq.eq.3)               tar(17)=tar(17)+1
      if(nq.eq.3.and.id.eq.130) tar(18)=tar(18)+1
      if(nq.eq.3.and.id.eq.-130)tar(19)=tar(19)+1
      if(nq.eq.3.and.id.eq.2130)tar(20)=tar(20)+1
      if(nq.eq.6)               tar(21)=tar(21)+1
      if(nq.eq.6.and.id.eq.130) tar(22)=tar(22)+1
      if(nq.eq.6.and.id.eq.-130)tar(23)=tar(23)+1
      if(nq.eq.6.and.id.eq.2130)tar(24)=tar(24)+1
      if(ns.eq.0.and.nq.eq.0)               tar(25)=tar(25)+1
      if(ns.eq.0.and.nq.eq.0.and.id.eq.130) tar(26)=tar(26)+1
      if(ns.eq.0.and.nq.eq.0.and.id.eq.-130)tar(27)=tar(27)+1
      if(ns.eq.0.and.nq.eq.0.and.id.eq.2130)tar(28)=tar(28)+1
      if(ns.eq.0.and.nq.eq.3)               tar(29)=tar(29)+1
      if(ns.eq.0.and.nq.eq.3.and.id.eq.130) tar(30)=tar(30)+1
      if(ns.eq.0.and.nq.eq.3.and.id.eq.-130)tar(31)=tar(31)+1
      if(ns.eq.0.and.nq.eq.3.and.id.eq.2130)tar(32)=tar(32)+1
      if(ns.eq.0.and.nq.eq.6)               tar(33)=tar(33)+1
      if(ns.eq.0.and.nq.eq.6.and.id.eq.130) tar(34)=tar(34)+1
      if(ns.eq.0.and.nq.eq.6.and.id.eq.-130)tar(35)=tar(35)+1
      if(ns.eq.0.and.nq.eq.6.and.id.eq.2130)tar(36)=tar(36)+1
      if(npo.lt.nco)goto7
         elseif(ifrptl(1,ipo).ne.0.and.istptl(ifrptl(1,ipo)).ne.2)then
      do 6 j=ifrptl(1,ipo),ifrptl(2,ipo)
      nco=nco+1
      if(nco.gt.ncomx)
     *call utstop('wclu: ncomx too small. increase dim.    ')
6     iar(nco)=j
      goto7
         else
      goto5
         endif
      if(nq.eq.0)            tar(37)=tar(37)+1
      if(nq.eq.3)            tar(41)=tar(41)+1
      if(nq.eq.6)            tar(45)=tar(45)+1
      if(ns.eq.0.and.nq.eq.0)tar(49)=tar(49)+1
      if(ns.eq.0.and.nq.eq.3)tar(53)=tar(53)+1
      if(ns.eq.0.and.nq.eq.6)tar(57)=tar(57)+1
      do 8 n=13,36
8     arclu(n,k)=arclu(n,k)+tar(n)
      do 9 n=37,60
9     arclu(n,k)=arclu(n,k)+tar(1+(n-1)/4*4)
5     continue
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wclufi(n,kmx)
c     analysis. cluster (resonance) formation and decay. normalization.
c     0 <= n <= 36, meaning of n see initial comment in subr wclu.
c-----------------------------------------------------------------------
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter(nclu=60,kclu=25)
      common/cclu/arclu(nclu,kclu),aclu(kclu)
      parameter(kmaxmx=10000)
      common/canar/ar(kmaxmx,5)
      kmx=kclu
      da=aclu(2)-aclu(1)
      do 1 k=1,kclu
      ar(k,1)=aclu(k)
      if(n.ge.1.and.n.le.12.and.arclu(n,k).gt.0.)then
      ar(k,3)=arclu(n,k)/nrevt/da
      ar(k,4)=arclu(n,k)/nrevt/da  /   sqrt(arclu(n,k))
      elseif(n.ge.13.and.n.le.36.and.arclu(n,k).gt.0.)then
      ar(k,3)=arclu(n,k)/arclu(n+24,k)
      ar(k,4)=arclu(n,k)/arclu(n+24,k)  /  sqrt(arclu(n,k))
      elseif(arclu(n,k).gt.0.)then
      call utstop('wclufi: only optn 1-36 possible.        ')
      endif
1     continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wimi
c     analysis. calculates factorial moments.
c-----------------------------------------------------------------------
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      integer nhis(32)
      double precision afacmo,dfacmo
      common /cfacmo/afacmo(5,32)
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      common/cwimiy/ymin,ymax
      integer ifacmo(5)
      ymin=-ymximi
      ymax=ymximi
 
      if(nrevt.eq.1)then
      do 44 m=1,5
      do 44 kmax=1,32
44    afacmo(m,kmax)=0
      endif
 
      do 45 kmax=1,32
      do 40 m=1,5
40    ifacmo(m)=0
      do 47 k=1,32
47    nhis(k)=0
      do 46 i=1,nptl
      if(istptl(i).ne.0)goto46
      amtptl=sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2)
      ylbptl=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
      k=1+(ylbptl-ymin)/(ymax-ymin)*kmax
      if(k.ge.1.and.k.le.kmax)nhis(k)=nhis(k)+1
46    continue
      do 41 k=1,kmax
      n=nhis(k)
      ifacmo(1)=ifacmo(1)+n
      ifacmo(2)=ifacmo(2)+n*(n-1)
      ifacmo(3)=ifacmo(3)+n*(n-1)*(n-2)
      ifacmo(4)=ifacmo(4)+n*(n-1)*(n-2)*(n-3)
      ifacmo(5)=ifacmo(5)+n*(n-1)*(n-2)*(n-3)*(n-4)
41    continue
      do 42 m=1,5
42    afacmo(m,kmax)=afacmo(m,kmax)+1d0*ifacmo(m)/kmax
45    continue
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wimino(mom,kmx)
c     analysis. normalizes factorial moments.
c-----------------------------------------------------------------------
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cnsta/pi,prom,piom,ainfin
      double precision afacmo,f1,fm
      common /cfacmo/afacmo(5,32)
      parameter(kmaxmx=10000)
      common/canar/ar(kmaxmx,5)
      common/cwimiy/ymin,ymax
      kmx=32
      do 45 kmax=1,32
      f1=afacmo(  1,kmax)/nrevt
      fm=afacmo(mom,kmax)/nrevt/f1**mom
      ar(kmax,1)=-alog((ymax-ymin)/kmax)
      if(fm.le.0d0)ar(kmax,3)=-ainfin
      if(fm.gt.0d0)ar(kmax,3)=dlog(fm)
45    continue
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wjpsi
c     jpsi analysis.
c-----------------------------------------------------------------------
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter (njpsi=50)
      common /cjpsi/xjpsi(njpsi),yjpsi(njpsi),zjpsi(njpsi)
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain
      detap=(ypjtl-yhaha)*.6667
      detat=-yhaha*.6667
      tpro=dcosh(detap)
      zpro=dsinh(detap)
      ttar=dcosh(detat)
      ztar=dsinh(detat)
 
      radsqr=1
      taumin=0
      taumax=20
 
      if(nrevt.eq.1)then
      do 4 nt=1,njpsi
      taus=taumin+(nt-1.)/(njpsi-1.)*(taumax-taumin)
      xjpsi(nt)=taus
      yjpsi(nt)=0
4     zjpsi(nt)=0
      endif
 
      do 1 i=1,nptl
      if(idptl(i).eq.441)j=i
1     continue
 
      do 2 nt=1,njpsi
      ttaus=xjpsi(nt)
      call uttain(j,xj,yj,zj,tj,n,1)
      if(n.eq.1.or.n.eq.2.or.n.eq.9)goto2
      zjpsi(nt)=zjpsi(nt)+1
      do 3 i=1,nptl
      if(i.eq.j)goto3
      call uttain(i,x,y,z,t,n,1)
      if(n.eq.1.or.n.eq.2.or.n.eq.9)goto3
      pde=(pptl(3,i)+pptl(3,j))/(pptl(4,i)+pptl(4,j))
      gam2i=1-pde**2
      if(gam2i.eq.0.)goto3
      if((x-xj)**2+(y-yj)**2+1/gam2i*(z-zj-(t-tj)*pde)**2
     *.gt.radsqr)goto3
      yjpsi(nt)=yjpsi(nt)+1
3     continue
2     continue
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wjpsif(ise,kmx)
c     jpsi analysis. normalized histogram -> /canar/
c-----------------------------------------------------------------------
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter(kmaxmx=10000)
      common/canar/ar(kmaxmx,5)
      parameter (njpsi=50)
      common /cjpsi/xjpsi(njpsi),yjpsi(njpsi),zjpsi(njpsi)
 
      do 1 i=1,njpsi
      ar(i,1)=xjpsi(i)
      ar(i,2)=0
      if(zjpsi(i).eq.0.)ar(i,3)=0.
      if(zjpsi(i).ne.0.)ar(i,3)=yjpsi(i)/zjpsi(i)
      ar(i,4)=0
      ar(i,5)=zjpsi(i)/nrevt
1     continue
 
      kmx=njpsi
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wspa
c     analysis. space-time evolution.
c     fills ar(i,j,k) with
c        ar: dspa:sum{1}, sspa:sum{strangn}, cspa:sum{p3/p4},
c            espa:sum{p4}, pspa:sum{p3}, tspa:sum{p3**2/p4}
c            aspa:sum{p5}, px2spa:sum{p1**2}, py2spa:sum{p2**2}
c        i: z-bin number  (for iwcent=1:baryon#-bin)
c        j: t-bin number
c        k: meson (1) or baryon (2)  (for iwcent=1: only k=1 filled)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      common/cnsta/pi,prom,piom,ainfin
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/paro2/nevent,modsho,engy,pnll,pnllx,yhaha,ish,ichoic,prosea
     *,laproj,maproj,latarg,matarg,prob(99),icbac(99,2),icfor(99,2),ivi
     *,istmax,ipagi,jpsi,jpsifi,isup,ndecay,elepti,elepto,angmue,ishsub
     *,idproj,idtarg,amproj,amtarg,ypjtl,imihis,isphis,rhophi,iclhis
     *,ijphis,ientro,kutdiq,ishevt,idpm,taurea,ispall,ymximi,ndecax
      parameter (nspaz=51,nspat=9)
      common /cspati/dspa(nspaz,nspat,2),cspa(nspaz,nspat,2)
     *,espa(nspaz,nspat,2),pspa(nspaz,nspat,2),tspa(nspaz,nspat,2)
     *,aspa(nspaz,nspat,2),px2spa(nspaz,nspat,2),py2spa(nspaz,nspat,2)
     *,zspa(nspaz),sspa(nspaz,nspat,2),atra,dznull
      parameter (mxptl=150000)
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),iclptl(mxptl),nqjptl(mxptl)
      integer ibaptl(mxptl),iacptl(mxptl),ntbptl(mxptl),ntsptl(mxptl)
      real amiptl(mxptl)
      parameter (nflav=6)
      integer jci(nflav,2),ici(2)
      integer jcdu(nflav,2)
      parameter (mxifr=mxptl)
      integer ifr(mxifr),ifrx(mxifr)
      mami=min0(maproj,matarg)
      r0=1.19*mami**(1./3.)-1.61*mami**(-1./3.)
      atra=(r0/2)**2*pi
      tmin=wtmini
      tstep=wtstep
      zmin=-6
      zmax=6
      iz0=(nspaz-1)/2+1
 
      if(nrevt.eq.1)then
      do 44 i=1,nspaz
      zspa(i)=zmin+(i-.5)*(zmax-zmin)/nspaz
      if(i.eq.iz0+4)zplus=zspa(i)
      if(i.eq.iz0-5)zminus=zspa(i)
      if(iwcent.eq.1)zspa(i)=i-25
      do 44 j=1,nspat
      do 44 k=1,2
      dspa(i,j,k)=0.
      sspa(i,j,k)=0.
      cspa(i,j,k)=0.
      espa(i,j,k)=0.
      pspa(i,j,k)=0.
      tspa(i,j,k)=0.
      aspa(i,j,k)=0.
      px2spa(i,j,k)=0.
44    py2spa(i,j,k)=0.
      dznull=zplus-zminus
      endif
 
      do 52 i=1,nptl
      call idquac(i,nq,ns,ndummy,jcdu)
      ibaptl(i)=1
      if(nq.ne.0)ibaptl(i)=2
      ntbptl(i)=nq/3
      ntsptl(i)=ns
      amiptl(i)=prom+amsiac
      if(nq.eq.0)amiptl(i)=piom+amsiac
52    continue
 
      nr=0
      do 45 it=1,nspat
      t=tmin+(it-1)*tstep
      do 43 i=1,nptl
43    iacptl(i)=1
      do 46 i=1,nptl
      if(istptl(i).gt.1)goto46
      if(iacptl(i).eq.0)goto46
      if(iabs(idptl(i)).lt.100.and.iabs(idptl(i)).ne.20)goto46
      if(ispall.eq.1)then
      if(t.lt.xorptl(4,i).or.t.ge.tivptl(2,i))goto46
      else
      if(t.lt.tivptl(1,i).or.t.ge.tivptl(2,i))goto46
      if(iclptl(i).eq.0.and.pptl(5,i).gt.amiptl(i))goto46
      endif
      x=xorptl(1,i)+(t-xorptl(4,i))*pptl(1,i)/pptl(4,i)
      y=xorptl(2,i)+(t-xorptl(4,i))*pptl(2,i)/pptl(4,i)
      z=xorptl(3,i)+(t-xorptl(4,i))*pptl(3,i)/pptl(4,i)
      if(x**2+y**2.gt.atra/pi)goto46
      iz=1+(z-zmin)/(zmax-zmin)*nspaz
      if(iz.lt.1.or.iz.gt.nspaz)goto46
      ib=ibaptl(i)
      if(iwcent.eq.1)then
      if(iz.lt.iz0-4.or.iz.gt.iz0+4)goto46
      iz=25+ntbptl(i)
      if(iz.lt.1.or.iz.gt.nspaz)goto46
      ib=1
      endif
      nr=nr+1
      dspa(iz,it,ib)=dspa(iz,it,ib)+1
      sspa(iz,it,ib)=sspa(iz,it,ib)+ntsptl(i)
      cspa(iz,it,ib)=cspa(iz,it,ib)+pptl(3,i)/pptl(4,i)
      espa(iz,it,ib)=espa(iz,it,ib)+pptl(4,i)
      pspa(iz,it,ib)=pspa(iz,it,ib)+pptl(3,i)
      tspa(iz,it,ib)=tspa(iz,it,ib)+pptl(3,i)**2/pptl(4,i)
      aspa(iz,it,ib)=aspa(iz,it,ib)+pptl(5,i)
      px2spa(iz,it,ib)=px2spa(iz,it,ib)+pptl(1,i)**2
      py2spa(iz,it,ib)=py2spa(iz,it,ib)+pptl(2,i)**2
      write(ifdt)nr,iz,it,ib,pptl(3,i),pptl(4,i)
 
c     children: iacptl=0
c     ------------------
      k1=0
      k2=0
      loop=0
31    loop=loop+1
      if(loop.eq.1)n=i
      if(ifrptl(1,n).gt.0)then
      if(ish.ge.92)write(ifch,*)n,' ---> ',ifrptl(1,n),ifrptl(2,n)
      do 30 k=ifrptl(1,n),ifrptl(2,n)
      if(iacptl(k).eq.0)goto30
      iacptl(k)=0
      k2=k2+1
           if(k2.gt.mxifr/10)then
      if(ish.ge.90.and.k2-1.eq.mxifr/10)then
      call utmsg('wspa  ')
      write(ifch,*)'*****  k2= ',k2,' > mxifr/10 = ',mxifr/10
      do 55 kk=1,k2
55    write(ifch,*)kk,'        ',ifrx(kk),'  --->  ',ifr(kk)
      call utmsgf
      endif
      if(k2.gt.mxifr)then
      call utmsg('wspa  ')
      write(ifch,*)'*****  k2= ',k2,' > mxifr = ',mxifr
      call utstop('wspa  : k2 > mxifr                      ')
      endif
           endif
      ifrx(k2)=0
      if(k.eq.ifrptl(1,n))ifrx(k2)=n
      ifr(k2)=k
30    continue
      endif
      k1=k1+1
      if(k1.le.k2)then
      n=ifr(k1)
      if(iacptl(n).ne.0)
     *call utstop('wspa  : should not happen               ')
      goto31
      endif
 
c     parents: iacptl=0
c     -----------------
      k1=0
      k2=0
      loop=0
12    loop=loop+1
      if(loop.eq.1)n=i
      if(iorptl(n).gt.0)then
      if(ish.ge.92)write(ifch,*)n,' <--- ',iorptl(n),jorptl(n)
           if(k2+2.gt.mxifr/10)then
      if(ish.ge.90.and.(k2+1.eq.mxifr/10.or.k2.eq.mxifr/10))then
      call utmsg('wspa  ')
      write(ifch,*)'*****  k2+2= ',k2+2,' > mxifr/10 = ',mxifr/10
      do 56 kk=1,k2
56    write(ifch,*)kk,'        ',ifrx(kk),'  <---  ',ifr(kk)
      call utmsgf
      endif
      if(k2+2.gt.mxifr)then
      call utmsg('wspa  ')
      write(ifch,*)'*****  k2+2= ',k2+2,' > mxifr = ',mxifr
      call utstop('wspa  : k2+2 > mxifr                    ')
      endif
           endif
      if(iacptl(iorptl(n)).ne.0)then
      k2=k2+1
      ifr(k2)=iorptl(n)
      ifrx(k2)=n
      iacptl(iorptl(n))=0
      endif
      if(jorptl(n).gt.0)then
      if(iacptl(jorptl(n)).ne.0)then
      k2=k2+1
      ifr(k2)=jorptl(n)
      ifrx(k2)=0
      iacptl(jorptl(n))=0
      endif
      endif
      endif
      k1=k1+1
      if(k1.le.k2)then
      n=ifr(k1)
      if(iacptl(n).ne.0)
     *call utstop('wspa  : should not happen (2)           ')
      goto12
      endif
 
46    continue
45    continue
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wspafi(ise,kmx)
c     analysis. space-time evolution.
c     ise=dcba: a=1:mesons ,a=2:baryons  (iwcent=1:a set to 1)
c               b=1:velocity, b=2:rapidity, b=3:density, b=4:en.density
c               b=5:mass/ptl, b=6:en*/ptl, b=7:pz2, b=8:px2, b=9:px2/pz2
c               b=0:strangeness/ptl
c               dc:time
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifhi,ifdt
      common/paro1/ptmx,gaumx,sigppi,core,fctrmx,neqmn,iaqu,qmust,cutmsq
     *,taunll,maxres,ptf,ptq,xcut,ioptq,irescl,ko1ko2,kentro,labsys
     *,ntrymx,delmss,pud,pspinl,pspinh,pispn,ncolmx,tensn,bmaxim,bminim
     *,cutmss,qvapc,qvatc,qsepc,qsetc,rstras,neqmx,taumx,nsttau,sigj
     *,pdiqua,parea,delrem,taumin,taumax,numtau,overlp,iopenu,iopent
     *,themas,amsiac,wproj,wtarg,iopbrk,ndecaw,radiac,radias,nclean
     *,ifrade,amprif,iojint,pth,pvalen,wtmini,wtstep,iwcent,phard,ioptf
      common/accum/imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl,amsac,ipage
     *,inoiac,ilamas
      parameter(kmaxmx=10000)
      common/canar/ar(kmaxmx,5)
      parameter (nspaz=51,nspat=9)
      common /cspati/dspa(nspaz,nspat,2),cspa(nspaz,nspat,2)
     *,espa(nspaz,nspat,2),pspa(nspaz,nspat,2),tspa(nspaz,nspat,2)
     *,aspa(nspaz,nspat,2),px2spa(nspaz,nspat,2),py2spa(nspaz,nspat,2)
     *,zspa(nspaz),sspa(nspaz,nspat,2),atra,dznull
      real vav(nspaz),rav(nspaz),pz2(nspaz)
 
      kmx=nspaz
      k=mod(ise,10)
      m=mod(ise/10,10)
      j=mod(ise/100,100)
      if(iwcent.eq.1)k=1
 
      do 44 i=1,nspaz
      pz2(i)=0
      vav(i)=0
      rav(i)=0
      if(dspa(i,j,k).eq.0.)goto44
      vav(i)=pspa(i,j,k)/espa(i,j,k)
      rav(i)=.5*alog((1+vav(i))/(1-vav(i)))
44    continue
 
      rewind ifdt
1     read(ifdt,end=2)nr,i,js,ks,ps,es
      if(js.eq.j.and.ks.eq.k)then
      p=-sinh(rav(i))*es+cosh(rav(i))*ps
      pz2(i)=pz2(i)+p**2
      endif
      goto1
2     continue
 
      do 45 i=1,nspaz
      ar(i,1)=zspa(i)
      if(dspa(i,j,k).eq.0.)goto45
      if(i.eq.1)dz=zspa(2)-zspa(1)
      if(i.gt.1)dz=zspa(i)-zspa(i-1)
      if(iwcent.eq.1)dz=dznull
      ano=1./nrevt/atra/dz
      gma=1./sqrt(1-vav(i)**2)
      den=(dspa(i,j,k)-cspa(i,j,k)*vav(i))*gma*ano
      ede=(espa(i,j,k)-2*pspa(i,j,k)*vav(i)+tspa(i,j,k)*vav(i)**2)*
     *gma**2*ano
      if(m.eq.1)ar(i,3)=vav(i)
      if(m.eq.2)ar(i,3)=rav(i)
      if(m.eq.3)ar(i,3)=den
      if(m.eq.4)ar(i,3)=ede
      if(m.eq.5)ar(i,3)=aspa(i,j,k)/dspa(i,j,k)
      if(m.eq.6)ar(i,3)=ede/den
      if(m.eq.7)ar(i,3)=pz2(i)/dspa(i,j,k)
      if(m.eq.8)ar(i,3)=px2spa(i,j,k)/dspa(i,j,k)
      if(m.eq.9.and.pz2(i).ne.0.)ar(i,3)=px2spa(i,j,k)/pz2(i)
      if(m.eq.0)ar(i,3)=sspa(i,j,k)/dspa(i,j,k)
45    continue
 
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wzzzz
c     calls arbitrary histogram routine; the call has to be inserted
c-----------------------------------------------------------------------
c     insert `call sr' here
      return
      end
 
c-----------------------------------------------------------------------
      subroutine wzzzzf(ise,kmx)
c     calls arbitrary normalization routine; the call has to be inserted
c-----------------------------------------------------------------------
c     insert `call sr' here
      return
      end

c-----------------------------------------------------------------------
      subroutine wzzzzi
c     calls arbitrary initializ. subroutine; the call has to be inserted
c-----------------------------------------------------------------------
c     insert `call sr' here
      return
      end

