      DOUBLE PRECISION Wmax,Wmin, Ymax,Ymin
C  Wmax, Wmin are the maximum and minimum photon energies (lab frame)
C  Ymax, YMIN are the maximum and minimum rapidities 
C  (Ymin only used for gamma p)
C  usually ymin=-ymax

      INTEGER numy,numw
C  numy and numw are the number of y and w bins respectively
C  for gammap, numy runs from ymin to ymax

      DOUBLE PRECISION	xinterfere
      INTEGER iinterfere

C  iiinterfere=0 (no interference) =1 (interference)
C  xinterfere is the interference fraction (0=none, 1 = full)

      DOUBLE PRECISION Eth, EgMax
C  Eth is the 
C  EgMax is the maximum photon energy (lab frame, for gP)

      DOUBLE PRECISION rdummy1
      INTEGER idummy1,idummy2
C  spare variables

      INTEGER iout,ievents,ip,iseed,gg_or_gP,ibreakup

      COMMON/inputp/iout,ievents,Wmax,Wmin,Ymax,Ymin,numw,numy,ip,iseed,
     & gg_or_gP, ibreakup, iinterfere,idummy1, xinterfere, rdummy1,
     &Eth,EgMax
C  iinterfere=0 (no interference) iinterference=1 (interference)
C  xinterfere = linear variable between 0 (no interference) and 1 (full
C  interference).  Only used when iinterference=1




