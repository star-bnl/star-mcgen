      Subroutine HemiCosm(Len,Rad,Off,MSpec,P_High,P_Low,D_Mode,
     +                    Theta,Phi,P,X,Y,Z,Charge,
     +                    End_Hits,Rate_PD,End_Nsph,Rate_Sph,
     +                    Version)

******-****************************************************************-******
*
*..Description: Simulates Cosmic Rays incident on Detector inside Hemi-sphere
*
*                 This Monte Carlo simulation models the incident Muon
*                 spectrum on a particle detector.
*
*                 A hemi-sphere placement method is used to be able to
*                   A)  generate rays of all Zenith angles (including 90 deg)
*                   B)  estimate the flux rate reliably
*
*                 The simulation expects:
*                   1) Detector mother volume length(cm)  [Len]
*                   2) Detector mother volume radius(cm)  [Rad]
*                   3) Detector mother volume offset(cm)  [Off]
*                   4) Momentum spectrum mode             [MSpec]
*                   5) Momentum Maximum cutoff(GeV/c)     [P_High]
*                   6) Momentum Minimum cutoff(GeV/c)     [P_Low]
*                   7) Data return mode                   [D_Mode]
*                 After the first(initialization) call, only the D_Mode
*                 parameter is relevant.
*
*                 The simulation returns:
*                   1) A momentum Vector(Rad and GeV/c)   [Theta Phi P]
*                   2) A position in 3-space(cm)          [X Y Z]
*                   3) A Muon charge(charge)              [Charge]
*                   4) Events on sphere, Rate on sphere,  [End_Nsph Rate_Sph
*                      Hits on detector, Rate on detector  End_Hits Rate_PD]
*                 The flux information is only returned after the final call,
*                 calculated from the available data concerning momentum and
*                 angular distribution, except for the rate on the sphere
*                 which is available already after the first call.
*                 The sources for each of the spectra are listed in
*                 the corresponding subroutines.
*
*                 The final arguments passed back to the main program are
*                 dependent upon the value of the D_Mode parameter, which sets
*                 what frame of reference you wish the data to be reported
*                 in.  Within the program, the Earth's frame of reference (_E
*                 suffix) is used.  Before the data is reported though, it is
*                 converted to the proper frame of reference depending
*                 on the value of the D_Mode Parameter.  D_Mode selects the
*                 data frame that the calling program wishes reported. The
*                 following options are available:
*                   D_Mode= 0: Null Call
*                      This mode does nothing at all.
*                   D_Mode= 1: Initialization of program
*                      The first call ever done to HemiCosm.  It sets-up
*                      all of the initial parameters. This mode is only
*                      called once.  It also reports some 'useful'
*                      information.
*                   D_Mode= 2: Do one simulation & report data in the detector
*                      frame.  This mode is selected every time a new event is
*                      required.  This will track a Muon from space to
*                      the cylinder and report the data back to the user
*                      in the detector's frame of reference(_D suffix.)
*                      Returns:[Theta Phi P X Y Z Charge]
*                   D_Mode=21 Report data in Earth frame
*                      This mode will not cause a new event to be created,
*                      but will instead report back the data from the
*                      current event in the Earth's frame of reference.
*                      Returns:[Theta Phi X Y Z]
*                   D_Mode=22 Report data in Shadow frame
*                      This mode will not cause a new event to be created,
*                      but will instead report back the data from the
*                      current event as if the particle were to have hit
*                      the X-Y Plane(Earth coordinate Frame).
*                      Returns:[X Y]
*                   D_Mode= 23: Do one simulation & report data in the detector
*                      frame.  This mode is similar to the D_Mode=2 call
*                      except, that in this mode the particle is only tracked
*                      to the hemisphere.  The particles position is returned
*                      in the detector's frame of reference.
*                      Returns:[Theta Phi P X Y Z Charge]
*                   D_Mode= 3 Final call of program
*                      This is the last call for the HemiCosm program.
*                      It is a final report of run data and also
*                      reports other information that can be useful to
*                      ensure that the program is running correctly.
*                      Returns:[Trial Rate_PD]
*                 The User can also choose which momentum spectrum they
*                 wish used by the setting of the P_Mode(MSpec) variable.  The
*                 P_Mode variable is passed into HemiCosm as MSpec and then
*                 stored under the P_Mode variable name.  The following
*                 momentum spectrum representations are present:
*                   P_Mode=0 Constant Spectrum = P_Min
*                   P_Mode=1 1/P^2 Spectrum
*                   P_Mode=2 Spectrum from measured values
*                   P_Mode=3 DarInv Spectrum
*                 The modes are explained in more detail in the flight
*                 subroutine.
*
*                 The Program works in the following outlined steps:
*                     i) Initialize the variables
*
*                    ii) Get Phi(Earth) angle for the momentum vector(Radian)
*                   iii) Get Theta angle(Earth) for the momentum vector(Rad)
*                    iv) Get P(momentum magnitude) of the particle(GeV/c)
*                     v) Get the particles charge(+1 or -1)
*
*                    vi) Calculate where the particle would hit the X-Y plane
*                   vii) Trace the impact point back to find where it would
*                         hit a hemisphere enclosing the particle detector
*                  viii) Convert the data to the Detectors frame of reference
*                    ix) Trace the particle to a Cylinder that just encloses
*                         the detector.
*
*                     x) Report the final data.
*
*..Optional program super-structure:
*    BegCosmFF*                        Optional FFREAD card definitions
*    /beghemicosm_c/ *                 Common block for interface to BEGET
*    BegHemiCosm*                      Optional interface to BEGET package
*      Begcosm*                          Optional filling of /STDHEP/
*
*..Program Structure:
*        HemiCosm(D_Mode=1)                Init phase:
*           Startup                           Detector, flux, printout
*              Shadow                            Area calculation variables
*        HemiCosm(D_Mode=2)                   Event phase:
*           Flight                               Sample phi, theta, p and q
*              GetRan*                              Random number, open (0, 1)
*              Getphi1                              Uniform phi sampling
*              Gettheta1                            Weighted cos(theta)**2
*                 CubeRoot*                            Solve cubic equation
*              GetP0                                P_mode=0:  flat
*              GetP1                                P_mode=1:  1/p**2
*              GetP2                                P_mode=2:  Measured flux
*                 GetP2ini                             First call:  integrate
*                    PhysDatSpec                          Physics Data flux
*                 IntrPol*                             Interpolate
*              GetP3                                P_mode=3:  Theory flux
*                 GetP3ini                             First call:  read data
*                    CosDarNN.dat* files                  Dar integral flux
*                 IntrPol*                             Interpolate
*              GetCharge                            Pick electric charge
*           Locate                               Place ray on detector
*              Shadow                               Hemisphere base area
*                 GetRan*                              Random number
*              Sphere                            Shadow to sphere impact
*              Frame                             Earth to detector rotation
*              Cylinder                          Sphere to cylinder impact
*        HemiCosm(D_Mode=3)                   Finish phase:
*           FinalInfo                            Report counters and rates
*              Shadow                               Report area calculation
*
*   ___
*   *
*    Code or data located in separate files:
*        BegCosmFF        in begcosmff.F
*        BegHemiCosm      in beghemicosm.F
*        /beghemicosm_c/  in beghemicosm_c.inc
*        BegCosm          in begcosm.F
*        GetRan           in getran.F
*        IntrPol          in intrpol.F
*        CubeRoot         in cuberoot.F
*        Dar flux         in cosdarNN.dat, where theta NN = 00, 05, ..., 90.
*
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Laboratory (LBL)
*
*..Copyright  : University of California, Berkeley, (LBL) 1994
*               All rights reserved
*
*..Credits    : Hemisphere method by
*                 J Frederic Kral and Michael P Bringle, LBL
*               Dar inverted integral flux spectrum data files by
*                 Henry Kendall, Massachussetts Institute of Technology (MIT)
*               Additional help from
*                 Henry Kendall and Amit Lath, MIT
*               Funding from
*                 U.S. Department of Energy
*
*..Contact    : J Frederic Kral, JFKral@LBL.Gov
*               
*..Versions   :
*   1.0         Bringle, August  5, 1994 - First development version.
*   2.0         Bringle, August 19, 1994 - Second development version.
*   2.1         Kral,    August 22, 1994 - Add Sphere N events and Rate.
*   2.2         Bringle, August 23, 1994 - Fix anomaly in last bin of DarInv.
*   2.3         Kral,    August 25, 1994 - First production version.
*   2.4         Kral,    Sept   14, 1994 - Correct Physics Data bibliography.
*
******-****************************************************************-******
      Implicit none
*..Version :
      Real*4  HemiVer           ! Program version
      Parameter (HemiVer = 2.4)
*..Input :
      Integer D_Mode           ! Program execution mode
      Integer MSpec            ! Momentum Spectrum Mode
      Real*4 Len               ! Mother Cylinder Length (cm)
      Real*4 Rad               ! Mother Cylinder Rad (cm)
      Real*4 Off               ! Mother Cylinder Offset from ground (cm)
      Real*4 P_Low,P_High      ! Momentum Spectrum Cutoffs (GeV/c)
*..Output:
      Integer Charge           ! Muon Charge (+1 or -1)
      Integer End_Hits         ! Number of particles thrown at the cylinder
      Integer End_Nsph         ! Number of hemi-sphere events
      Real*4 Theta,Phi         ! Momentum direction(radian)
      Real*4 P                 ! Momentum (GeV/c)
      Real*4 X,Y,Z             ! Particle location (cm)
      Real*4 Rate_Sph          ! Predicted Rate (Hz) on the hemisphere
      Real*4 Rate_PD           ! Actual Rate (Hz) of particles hitting
                               !   the detector (See Shadow)
      Real*4  Version          ! Program version (output)
*..Local :
      Integer P_Mode           ! Momentum Spectrum Mode
      Real*4 Det_Len           ! Mother Cylinder Length(cm)
      Real*4 Det_Rad           ! Mother Cylinder Radius(cm)
      Real*4 Det_Off           ! Mother Cylinder Offset from ground(cm)
      Real*4 P_Min,P_Max       ! Momentum Spectrum Cutoffs(GeV/c)

      Integer Flag             ! Reports if particle hit cylinder or not
      Integer Runs             ! Number of particles thrown at the cylinder
      Integer Trials           ! Number of particles hitting the cylinder
      Real*4 dT                ! Time step for each event(seconds)
      Real*4 HemiRad           ! Enclosing hemisphere radius(cm)
      Real*4 Theta_E,Phi_E     ! Momentum direction(Earth Frame)(Radians)

*..Save Variables
      Save dT,HemiRad          ! Save initialized parameters
      Save Theta_E,Phi_E       ! Save incase of a second mode call
      Save Runs,Trials         ! Save for use in final call
      Save Det_Len,Det_Rad     ! Save because are used on every call
      Save Det_Off,P_Mode      ! Save because are used on every call
      Save P_Min,P_Max         ! Save because are used on every call
*=====-================================================================-=======
* Standard create event and report it in the Detector Frame
      If (D_Mode.eq.2) then
        Flag=0
        Trials=Trials+1

        Do While (Flag.ne.1)    ! Loop to ensure a cylinder hit
          Runs=Runs+1

*         Get Particle's physical description
          Call Flight(P_Mode,P_Min,P_Max,Theta_E,Phi_E,P,Charge)

*         Get the particle's location in space
          Call Locate(D_Mode,HemiRad,Det_Len,Det_Rad,Det_Off,Theta_E,
     +                Phi_E,Theta,Phi,X,Y,Z,Flag)
        Enddo


* Reformat the data and report it in the Earth Frame
      Elseif (D_Mode.eq.21) then
        Call Locate(D_Mode,HemiRad,Det_Len,Det_Rad,Det_Off,Theta_E,
     +              Phi_E,Theta,Phi,X,Y,Z,Flag)

* Reformat the data and report it in the Shadow Frame
      Elseif (D_Mode.eq.22) then
        Call Locate(D_Mode,HemiRad,Det_Len,Det_Rad,Det_Off,Theta_E,
     +              Phi_E,Theta,Phi,X,Y,Z,Flag)

* Generate an event and use the hemisphere as the target
      Elseif (D_Mode.eq.23) then
        Trials=Trials+1           !Ratio should be 1:1 for this mode
        Runs=Runs+1

*       Get Particle's physical description
        Call Flight(P_Mode,P_Min,P_Max,Theta_E,Phi_E,P,Charge)

*       Get the particle's location in space
        Call Locate(D_Mode,HemiRad,Det_Len,Det_Rad,Det_Off,Theta_E,
     +              Phi_E,Theta,Phi,X,Y,Z,Flag)

* Initialize the program and swap variable names(First Call)
*   Set these values on startup so subsequent calls to HemiCosm
*   from routines other that the main one(begcosm),ie debug
*   calls, need not include this information.
      Elseif (D_Mode.eq.1) then
        Version = HemiVer       ! Program version
        P_Mode=MSpec            ! Momentum Spectrum Mode
        P_Min=P_Low             ! Low momentum Cutoff
        P_Max=P_High            ! High momentum Cutoff
        If (P_Mode.eq.0) P_Max=P_Min ! In case of constant spectrum
        Det_Rad=Rad             ! Mother Cylinder Radius(cm)
        Det_Len=Len             ! Mother Cylinder Length(cm)
        Det_Off=Off             ! Mother Cylinder Offset from ground(cm)

        Call StartUp(Det_Len,Det_Rad,Det_Off,P_Min,P_Max,P_Mode,
     +                   dT,HemiRad,Version)
        Rate_Sph=1.0/dT         ! Rate on hemi-sphere

* Print the final data of the program(Final Call)
      Elseif (D_Mode.eq.3) then
        Call FinishInfo(Trials,Runs,dT,Rate_PD)
        End_Hits=Trials         ! Number of hits on the detectors
        End_Nsph=Runs           ! Number of hemi-sphere events

* Do completely nothing
      Elseif (D_Mode.eq.0) then !Null Call

* Improper mode select
      Else
        Write(*,*)
        Write(*,*)"HemiCosm:  Error - aborting"
        Write(*,*)"Improper D_Mode selected:"
        Write(*,*)"  D_Mode= ",D_Mode
        STOP                    ! Improper D_Mode
      Endif

 999  Return
      End

      Subroutine StartUp(Det_Len,Det_Rad,Det_Off,P_Min,P_Max,P_Mode,
     +                   dT,HemiRad,Version)

******-****************************************************************-******
*
*..Description: HemiCosm initialization Routine
*                 This subroutine first calculates the smallest radius
*                 hemisphere that will enclose the detector. Next it uses
*                 this radius to calculate the incident flux rate on the
*                 hemisphere.  From this it calculates the 'average' time
*                 step to be assigned to each particle impact(event).  The
*                 subroutine also reports some data which may be useful,
*                 and warning messages as needed.  Calls to other routines
*                 are used to initialize them with this data.
*
*                 NOTE: The flux is a function of momentum.  If the limits
*                   of P_Min and P_Max are changed in then the total
*                   percentage of flux between 0.3 and 1000 GeV/c that this
*                   corresponds to, must be recalculated and substituted into
*                   the FluxRatio parameter.  If the values of P_Min and
*                   P_Max differ, significantly, from the values used to
*                   calculate the current flux a warning message will be
*                   displayed stating that the flux rate will be incorrect.
*                   The current value for the flux percentage was calculated
*                   using the data taken from:
*                       Physics Data 25-1 "Cosmic Rays on Earth" (1984)
*                        (ISBN 0344-8401) by O.C. Allkofer and P.K.F. Grieder.
*                        Table 2.2.1                 page 240
*                   Using the intensity from around 0.3 GeV/c and 1.0
*                   GeV/c, I extrapolated how much was present at each
*                   point and then divided to find the ratio of flux
*                   at 1.0 GeV/c to 0.3 GeV/c.  This number is not
*                   very accurate(good to about 10%) so a P_Min +- 0.2GeV/c
*                   and a P_Max above 40 GeV/c will still give an 'accurate'
*                   flux rate.  If the Momentum cutoff fall outside this
*                   range a warning message will be displayed.  The value
*                   for the rate of incident flux was taken from:
*                      Particle Properties Data Booklet  June 1992
*                      American Institute of Physics         P.129
*                   The hard component(>0.3 GeV/c) crossing a unit
*                   horizontal area from above was used.  Because of the
*                   shadow method(discussed later) a factor of 5/4 was used
*                   to give the correct rate of the hemisphere.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Integer P_Mode
      Real*4 Det_Len     ! Mother Cylinder Length(cm)
      Real*4 Det_Rad     ! Mother Cylinder Radius(cm)
      Real*4 Det_Off     ! Mother Cylinder Offset from ground(cm)
      Real*4 P_Min,P_Max ! Momentum Spectrum Cutoffs(GeV/c)
      Real*4 Version     ! Version of Program (filled in HemiCosm)
*..Output :
      Real*4 dT          ! Time Step for each hit(seconds)
      Real*4 HemiRad     ! Hemisphere radius (cm)
*..Local :
      Real*4 Flux        ! Incident flux rate (>0.3GeV/c) (Part/cm^2 sec)
      Real*4 Flux0       ! Flux rate greater than 1.0GeV/c (Part/cm^2 sec)
      Real*4 FluxRatio   ! Ratio of flux accepted to total(See Note)
      Real*4 Rate        ! Rate of particles hitting the hemisphere(Hz)
      Real*4 Geomfact           ! Geometry factor for hemi-sphere to base
      Parameter(Flux=0.013)
      Parameter(FluxRatio=0.7642)
      Parameter(Flux0=Flux*FluxRatio)
      Parameter (Geomfact = 5.0/4.0)

*..Constants
      Real*4 Pi             ! Constant
      Parameter (Pi=3.14159265358979323846264338328)
*=====-================================================================-=======
*Find minimum enclosing sphere radius
      HemiRad=Sqrt( (2.0*Det_Rad+Det_Off)**2 + (Det_Len/2.0)**2 )

*Get rate and time step for each particle hit
      Rate=(Flux0*Pi*HemiRad**2)*Geomfact !Rate=Flux times area(Hz)
      dT=1.0e0 / Rate                     !Time per particle(Seconds)

*Initialize Shadow Routine
      Call Shadow(1,Flux0,Rate,0.0,0.0,0.0,0.0,0.0)

*Report data to whoever is listening
 10   Format(A,F13.4,A)
 11   Format(A,I8)
 12   Format(A)
 13   Format(A,F3.1,A)
 14   Format(A,F13.2)

      Write(*,*)
      Write(*,12)" ***********************************"
      Write(*,12)" * HemiCosm - Cosmic Ray Generator *"
      Write(*,13)" *            Version ",Version,
     +                                    "          *"
      Write(*,12)" ***********************************"
      Write(*,*)
      Write(*,12)" HemiCosm/StartUp Data:"
      Write(*,10)"   Predicted rate of impacts: ",Rate,
     +                                          " (Hz) on Hemisphere"
      Write(*,14)"   Hemisphere geomet. factor: ",Geomfact
      Write(*,10)"   Rate before geometry     : ",Rate/Geomfact,
     +                                          " (Hz) on base"
      Write(*,10)"   Incident Flux(>0.3 GeV/c): ",Flux, " (s^-1 cm^-2)"
      Write(*,10)"   Flux Ratio (1.0/0.3GeV/c): ",FluxRatio
      Write(*,10)"   Incident Flux(>1.0 GeV/c): ",Flux0," (s^-1 cm^-2)"
      Write(*,10)"   Hemisphere radius        : ",HemiRad," (cm)"
      Write(*,10)"   Cylinder Radius          : ",Det_Rad," (cm)"
      Write(*,10)"   Cylinder Length          : ",Det_Len," (cm)"
      Write(*,10)"   Cylinder Offset          : ",Det_Off," (cm)"
      Write(*,10)"   Momentum Minimum         : ",P_Min," (GeV/c)"
      Write(*,10)"   Momentum Maximum         : ",P_Max," (GeV/c)"
      Write(*,11)"   Momentum Mode            : ",P_Mode

*Warning Message for Momentum out of range
*Note P_Mode=0 sets off this error(correctly) because P_Max is set equal to
*  P_Min upon start-up, and the two values will always fail the test
      If (P_Min.gt.1.2.or.P_Min.lt.0.8.or.P_Max.lt.40) then
        Write(*,*)"   The requested momentum spectrum is outside the"
        Write(*,*)"     limits of accuracy of the calculated rate.  "
        Write(*,*)"     The reported rate and elapsed time will not be"
        Write(*,*)"     correct."
        Write(*,*)
      Endif

      Write(*,*)

 999  Return
      End


      Subroutine Flight(P_Mode,P_Min,P_Max,Theta_E,Phi_E,P,Charge)

******-****************************************************************-*******
*
*..Description: This program calls the programs that generate the
*                 Phi, Theta, P, and Charge distributions.  Four random
*                 numbers are generated and then passed to the respective
*                 subroutines which return the parameters.  The individual
*                 spectrum are generated in the following ways:
*                   Phi  : 0 -> 2Pi with a flat distribution(Analytical
*                          formula)
*                   Theta: Pi/2->Pi with a Cos^2 distribution(Analytical
*                          formula)
*                   P    : P_Min->P_Max using a method determined by the
*                          P_mode parameter.
*                          The choices include:
*                           P_Mode=0 Constant Spectrum P=P_Min
*                           P_Mode=1 1/p^2 Distribution(Analytical)
*                           P_Mode=2 Measured Spectrum [Physics Data(1984)]
*                                    (Interpolation)
*                           P_Mode=3 DarInv calculated spectrum(Interpolation)
*                                    [Arnon Dar and Henry Kendall]
*                   Charge: +1 or -1 (Analytical formula)
*                 A more detailed description of each procedure and references
*                    can be found in the respective subroutine.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Integer P_Mode            ! Momentum Spectrum Mode(GeV/c)
      Real*4 P_Min,P_Max        ! Momentum Range(GeV/c)
*..Output:
      Integer Charge            ! Muon Charge
      Real*4 Theta_E,Phi_E      ! Angular orientation(Earth Frame)(Radians)
      Real*4 P                  ! Momentum (GeV/c)
*..Local:
      Integer L1                ! Do loop variable
      Integer Flag              ! Check Flag
      Integer Skip/1/           ! Random Number Index Advancer
      Real*4 Rnum(4)            ! Random Number variables
      Logical*4 First/.true./   ! First call
      Save First
*..Functs:
      Real*4 GetRan             ! Random number generator
      External GetRan           ! Random number alias function
*=====-================================================================-=======
*Get four random numbers for picking the parameters
      Do L1=1, 4
        RNum(L1)=GetRan(L1+Skip)
      Enddo

*Get Phi
      Call Getphi1(Rnum(1),Phi_E)

*Get Theta
      Call Gettheta1(Rnum(2),Theta_E)

*Get P
      If (First) then
        First = .false.
        Write (*, '(/, 1X, A, I3)')
     +    'HemiCosm/Flight:  Momentum spectrum selection', P_Mode
        Write (*, *) '  Pmode   Type'
        Write (*, *) '    0     Constant:  p = p_min'
        Write (*, *) '    1     1/p^2'
        Write (*, *) '    2     Measured Spectrum (Physics Data 1984)'
        Write (*, *) '    3     Theoretical Spectrum (Arnon Dar 1984)'
        Write (*, *) ' '
      Endif                     ! (First)
      Flag=0
      Do While (Flag.eq.0) !Ensure that it is in the proper range
        If (P_Mode.eq.1) then
          Call GetP1(P_Min,P_Max,Rnum(3),P)
        Elseif (P_Mode.eq.2) then
          Call GetP2(P_Min,P_Max,Rnum(3),P)
        Elseif (P_Mode.eq.3) then
          Call GetP3(P_Min,P_Max,Rnum(3),Theta_E,P)
        Elseif (P_Mode.eq.0) then
          Call GetP0(P_Min,P)
          Flag=1
        Endif

        If (P.ge.P_Min.and.P.le.P_Max) Flag=1 !Check to see that P is
                                              !  in proper range.
        Skip=Skip+1
        Rnum(3)=GetRan(L1+Skip) !Get a new random number in case a new
                                ! value is needed.
      Enddo

*Get Charge
      Call GetCharge(Rnum(4),Charge)

 999  Return
      End


      Subroutine Getphi1(Rnum,Phi_E)

******-****************************************************************-*******
*
*..Description: This program gets a randomly distributed phi angle in
*                 the Earth frame of reference.  It is assumed that there
*                 is no preference in the phi direction so the entire
*                 range 0 -> 2Pi has equal probability.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 Rnum       ! Random Number
*..Output:
      Real*4 Phi_E      ! Phi Angle(Earth)
*..Local:
      Real*4 Pi         !Constant
      Parameter (Pi=3.14159265358979323846264338328)
*=====-================================================================-=======
      Phi_E=Rnum*2.0d0*Pi

 999  Return
      End


      Subroutine Gettheta1(Rnum,Theta_E)

******-****************************************************************-*******
*
*..Description: This program solves the integral of the cross-section
*                 times the flux.  Given the current value for Rnum the
*                 subroutine will return a value of Theta picked randomly
*                 between 0 and Pi/2 radians according to the well known
*                 Cos^2 distribution. At the end we multiply Cos(theta)
*                 by -1 so it will be in the range of Pi/2 -> Pi radians
*                 to correspond to its downward trajectory.
*
*                 How It Is Done:
*
*                                         (    1)
*                 cross-section*Intensity=(1 + -) * x^2       x=Cos(Theta)
*                                         (    x)
*
*                  i) Integrate that over all x(0 -> 1) to get a
*                       normalized spectrum.
*
*                 ii) Next find an expression for the integral from 0 -> x
*                       till we get the enclosed area u(our random number)
*                       the equation is:
*                       x^3+(3/2)x^2+0-(5/2)u=0 (X=Cos(Theta))
*
*                iii) Finally solve that in terms of x and you have your
*                       answer (in terms of Cos(Theta) at least)
*
*                 The method used to solve the equation is Cardan's formula
*                 from 1545, and can be found in:
*                  The Universal Encyclopedia of Mathematics
*                  By James R Newman     ISBN(0-671-20348-7)  1964
*                  Pg. 197-198
*
*                 Note:
*                  If you change this distribution then you must also account
*                  for the change in the flux.  With the Cos^2 distribution
*                  the multiplication factor is 5/4.  The factor arises
*                  because of the fact our object is not a flat plane.  We
*                  must therefore count the additional particle which hit
*                  outside the 'footprint' of our hemisphere.  The following
*                  calculation shows how to find the(and how we got our)
*                  weight factor for any other theta distribution:
*
*                            / 1
*                            |   Area(x)*F(x) dx
*                           _/ 0
*                  Weight = --------------------
*                                 / 1
*                                 |   F(x) dx
*                                _/ 0
*
*                         X=Cos(Theta)
*                         Area(X)=Area of Shadow(cross-section) at given X
*                             (     1 )
*                            =(1 + ---)
*                             (     X )
*                         F(x)=Theta distribution function
*                            =Cos^2(Theta)=X^2
*
*
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 Rnum           !Random number input
*..Output:
      Real*4 Theta_E        !Theta trajectory(Earth Frame)
*..Local:
      Real*4 p,q,u,v        !Equation coefficients and simplifications
      Real*4 a,b,c,Phi      !Equation coefficients and simplifications
      Real*4 y1             !Substitution solution
      Real*4 S1             !Equation Solution
      Real*4 test           !Applicability to Cardan's Rule test
      Real*4 CosTheta       !Proper Solution to the equation
      Real*4 INV_3          !Pre Calculated number for speed
      Real*4 Pi             !Constant

      Parameter (Pi=3.14159265358979323846264338328)
      Parameter (INV_3=1.0d0/3.0d0)
*..Functs:
      Real*4   CubeRoot        !Takes the proper cube root of a neg number
      External CubeRoot
*=====-================================================================-=======
*Solving x^3+(3/2)x^2+0-(5/2)Rnum=0 (X=Cos(Theta))
      a= 1.5d0
      b= 0.0d0
      c=-2.5d0*Rnum

      p=(-(a*a)*INV_3 + b)*INV_3
      q=(2.0d0*(a*a*a)/27.0 - a*b*INV_3 + c)*0.5d0

*..Test is what needs to be positive for Cardan's rule to work
      test = q*q + p*p*p

      If (test .ge. 0.0d0.and.Rnum.gt.0.0d0) then
*       Solution satisfying q^2+p^3>=0 (Valid for Rnum=1 -> 0.2)
        u=CubeRoot(-q+Sqrt(test))
        v=CubeRoot(-q-Sqrt(test))
        y1=u + v

*       Only the 'real' root remains.  The other two are imaginary
*         or equal to zero so they can be ignored.
        S1=y1-a*INV_3
        CosTheta=S1
      Else
*       Solution satisfying p^3+q^3<=0(Valid for Rnum=0 -> 0.2)
        Phi=Acos(-q/Sqrt(-p*p*p))

*       The two other roots are discarded because they give a negative
*         result.  So only this one remains.
        y1= 2.0d0*Sqrt(-p)*Cos(Phi*INV_3)
        S1=y1-a*INV_3
        CosTheta=S1
      EndIf

      Theta_E=Acos(-CosTheta)   !Particle is traveling downwards

999   Return
      End


      Subroutine GetP0(P_Min,P)

******-****************************************************************-*******
*
*..Description: This program returns a momentum value equal to P_Min
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 P_Min         ! Momentum Value(GeV/c)
*..Output:
      Real*4 P             ! Momentum Value(GeV/c)
*..Local:
      Integer Flag         ! So will do initialization routine only once
*..Saved
      Save Flag
*=====-================================================================-=======
      P=P_Min

      If (Flag.eq.0) then
        Flag=1
*       Report Information to User
        Write(*,*)
        Write(*,*)"Momentum Spectrum Data:"
        Write(*,*)"   This simulation is using a momentum equal"
        Write(*,*)"   to P_Min."
        Write(*,*)
      Endif

 999  Return
      End


      Subroutine GetP1(P_Min,P_Max,RNum,P)

******-****************************************************************-*******
*
*..Description: This program accepts a random number and returns
*                 a momentum value associated with that number.
*               The spectrum used in this method is:
*                   Intensity=(Constant/p^2)
*               The limits on this equation are P_Min and P_Max.  The
*                 equation has been solve explicitly so that the solution
*                 is an analytical one.
*
*               How it is done:
*                i) Integrate the 1/p^2 from P_Min to P_Max to find the
*                   normalization constant.
*
*               ii) Using the normalization constant now integrate the
*                   constant/p^2 from P_Min to x to get the enclosed area(the
*                   random number) and solve for x
*
*                    /- x   constant                    (  1      1 )
*                    |      -------- dp = u => constant*(----- - ---)
*                   _/ P_Min  p^2                       (P_Min    x )
*
*               iii) The solution is easy and analytical:
*                        (  1         u    )-1
*                    x=P=(----- -  --------)
*                        (P_Min    constant)
*
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 RNum          ! Random Number
      Real*4 P_Min,P_Max   ! Momentum Cutoff(GeV/c)
*..Output:
      Real*4 P             ! Momentum (GeV/c)
*..Local:
      Integer Flag         ! So will do initialization routine only once

                           ! (These are inverse to speed calculation)
      Real*4 IPmax,IPmin   ! Inverse P maximum and minimum (GeV/c)
      Real*4 Iconstant     ! Inverse Normalization Constant
*..Saved
      Save Flag,IPmin,IPMax,Iconstant
*=====-================================================================-=======
*Initialization Routine
      If (Flag.eq.0) then
        Flag=1
        IPmin=1.0d0/P_Min         !Inverse PMin
        IPmax=1.0d0/P_Max         !Inverse PMax
        Iconstant=(IPmin - IPMax) !Inverse constant

* Report Information to User
        Write(*,*)
        Write(*,*)"Momentum Spectrum Data:"
        Write(*,*)"   This simulation is using a momentum spectrum"
        Write(*,*)"   based on the 1/p^2 differential spectrum curve."
        Write(*,*)"   The curve is normalized between the P_Max and"
        Write(*,*)"   P_Min limits."
        Write(*,*)"   The momentum spectrum is valid for any values of"
        Write(*,*)"   of P but diverges from the actual spectrum at"
        Write(*,*)"   low values of momentum(<1.0GeV/c)."
        Write(*,*)
      Endif

      P=1.0d0/(IPmin-RNum*Iconstant)    !Finds the momentum from Rnum

 999  Return
      End


      Subroutine GetP2(P_Min,P_Max,Rnum,P)

******-****************************************************************-*******
*
*..Description:  Physics Data spectrum
*                  This method uses measured flux values of the Muon
*                  spectrum vs momentum to calculate the momentum
*                  distribution.  The data was taken from:
*                    Physics Data 25-1 "Cosmic Rays on Earth" (1984)
*                     (ISBN 0344-8401)   O.C. Allkofer and P.K.F. Grieder
*                     Table 2.2.3              page 241
*                  The data points were fit to a curve and that curve is
*                  used in the calculation of the momentum spectrum.  The
*                  curve was integrated numerically(see corresponding
*                  subroutine) to get the integrated spectrum.  The array
*                  is searched and when the closest bins to the random
*                  number are found the actual data point is then found by
*                  cubic interpolation.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit None
*..Input:
      Real*4 Rnum               ! Random Number
      Real*4 P_Max,P_Min        ! Momentum Range(GeV/c)
*..Output:
      Real*4 P                  ! Momentum(GeV/c)
*..Local:
      Integer Bins              ! Number of divisions to break the
                                !   momentum spectrum into.
      Parameter (Bins=2000)

      Integer L1                ! Loop Variable
      Integer Flag              ! Initialization flag so only called once
      Integer Indx              ! Stores first Interpolation bin number
      Real*4 X4(4),Y4(4)        ! Store points for Interpolation
      Real*8 F(Bins)            ! Integrated Spectrum values
      Real*8 X_F(Bins)          ! P values for a given F(#)
*..Funct:
      Real*4 IntrPol            ! An interpolation routine
      External IntrPol
*..Saved:
      Save Flag                 ! Save so initialization routine only done once
      Save F,X_F                ! Save for use on later calls
*=====-================================================================-=======
*     Get the integrated Spectrum
      If (Flag.eq.0) then
        Call GetP2Ini(Bins,P_Min,P_Max,F,X_F)
        Flag=1
      Endif

* Take Random number and find which bins it lies between
      Do L1=3,Bins-4
        If (RNum.lt.F(L1)) then
          Indx=L1-3
          goto 100
        Endif
      Enddo

 100  Do L1=1,4
        X4(L1)=F(Indx+L1)
        Y4(L1)=X_F(Indx+L1)
      Enddo

      P=IntrPol(3,Rnum,X4,Y4) !Cubic Interpolation

      End


      Subroutine GetP2Ini(Bins,P_Min,P_Max,F,X_F)

******-****************************************************************-*******
*
*..Description: The program performs an integration of the input data
*                 to return a normalized spectrum.  The integration
*                 method used is the trapezoidal method with resolution
*                 enhancements to reduce the error.  For each 'bin' the
*                 curve is calculated at resolution number of points and
*                 that area is integrated and stored under the bin value.
*               The normalized  spectrum is returned in the F array in
*                 which the point F(L1) corresponds to the area under the
*                 curve from  X_F(1)->X_F(L1).  So F(1)=0 and F(Bins)=1
*                 The program also does the integration starting with the
*                 side that has the smallest Y values to help reduce
*                 rounding-off errors.  It was found that real*8 was
*                 necessary or the area divergent.
*
*               The function used in this method(PhysDatSpec) was derived
*                 from the data in:
*                    Physics Data 25-1 "Cosmic Rays on Earth" (1984)
*                     (ISBN 0344-8401)   O.C. Allkofer and P.K.F. Grieder
*                     Table 2.2.3              page 241
*                  and is a curve fit based upon the data in the table.  The
*                  curve is valid between 0.5 and 1000 GeV/c.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
*
******-****************************************************************-*******
      Implicit None
*..Input:
      Integer Bins              ! Number of divisions to break the
                                !   momentum spectrum into
*..Output:
      Real*8 F(Bins)            ! Integrated Spectrum values
      Real*8 X_F(Bins)          ! P values for a give F(#)
*..Local:
      Integer L1,L2             ! Loop variables
      Integer Res_Start,Res_End ! Resolution loop limits
      Integer Resolution        ! Number of steps between bins
      Real*4 P_Min,P_Max        ! Momentum Cutoff(GeV/c)
      Real*8 h1,h2              ! Height of trapezoid sides
      Real*8 dx                 ! Width of trapezoid
      Real*8 X                  ! Current X value
      Real*8 Section_Area       ! Area under the curve of the current section
      Real*8 Tot_Area           ! Total cumulative area under curve

*..Functions:
      Real*8 PhysDatSpec        ! Differential Momentum Spectrum(See above)
      External PhysDatSpec

*=====-================================================================-=======
*Report some information to the user
      Write(*,*)"Momentum Spectrum Data:"
      Write(*,*)"   This simulation is using a momentum spectrum"
      Write(*,*)"   based on the data in:"
      Write(*,*)"     O.C. Allkofer and P.K.F. Grieder"
      Write(*,*)"     Cosmic Rays on Earth"
      Write(*,*)"     Physics Data, No. 25-1 (1984) ISBN 0344-8401"
      Write(*,*)"     Table 2.2.3 on page 241"
      Write(*,*)"   and is a curve fit to the data in the table."
      Write(*,*)"   The curve is valid between 0.5 and 1000 GeV/c."

      If (P_Max.gt.1000.0.or.P_Min.lt.0.5) then
        Write(*,*)
        Write(*,*)"GetP2Ini:  Error - program terminating"
        Write(*,*)"   *************** WARNING ******************"
        Write(*,*)"   The Momentum limits you have chosen are"
        Write(*,*)"   outside the allowed range(see above.)"
        Write(*,10)"     P_Min allowed:    0.5 GeV/c   Yours: ",
     +            P_Min," GeV/c"
        Write(*,10)"     P_Max allowed: 1000.0 GeV/c   Yours: ",
     +            P_Max," GeV/c"
        Write(*,*)"   Program terminating..."
        Write(*,*)
        Stop    !Stopped in GetP2Ini
      Endif

      Write(*,*)

*Set-up boundary Values
      Resolution=1000     !Points to use in between bin values to find
                          !  the Area of the Curve

      Res_Start=Resolution-1
      Res_End  =0

*Find step Value
      dX=(1.0d0*(P_Max-P_Min)) / ((Bins-1)*resolution)

*Setup Last Bin(Should be zero by definition)
      F(Bins)=0.0d0
      X_F(Bins)=P_Max

*Ready Parameters to start on the Next to last bin
      L1=Bins-1
      X=P_Max+dX

*Calculate Areas
*     First get the total area under curve at each point(Backwards)
      Do While (L1.ge.1)
        Section_Area=0.0d0
        L1=L1-1

        Do L2=Res_Start, Res_End, -1
          X=X-dX
          h1=PhysDatSpec(X)
          h2=PhysDatSpec(X-dX)
          Section_Area=Section_Area+0.5d0*(h1+h2)*dX
        Enddo

        Tot_Area=Tot_Area+Section_Area
        F(L1+1)=Tot_Area
        X_F(L1+1)=X-dX
      Enddo

*     Second normalize by dividing each point by total area
      Do L1=1, Bins
        F(L1)=F(L1)/Tot_Area
      Enddo

*     Third Make sure the order is such that F(1)=0
      Do L1=1, Bins
        F(L1)=1.0d0-F(L1)
      Enddo

*Report Information to User
 10   Format(A,F13.4,A)
 11   Format(A,I8)
 12   Format(A,D21.8,A)

      Write(*,*)
      Write(*,*)"Curve Integration Package:"
      Write(*,*) " The total area of the chosen momentum"
      Write(*,*) " spectrum between the momentum limits is:"
      Write(*,12)"   Total Area: ",Tot_Area
      Write(*,*) " The width of each momentum bin:"
      Write(*,12)"   dX        : ",dX," (GeV/c)"
      Write(*,10)"   P-Minimum : ",P_Min," (GeV/c)"
      Write(*,10)"   P-Maximum : ",P_Max," (GeV/c)"
      Write(*,*) " The number of points in each bin used for"
      Write(*,*) " calculation of the area:"
      Write(*,11)"   Resolution: ",Resolution
      Write(*,*) " The number of momentum bins:"
      Write(*,11)"   Bins      : ",Bins
      Write(*,*)

      Return
      End


      Real*8 Function  PhysDatSpec(X)

******-****************************************************************-*******
*
*..Description: Differential Momentum Spectrum
*               This spectrum is a function that is used my the
*                 integration package to obtain an integrated spectrum.
*                 This integrated spectrum can then be used to find a
*                 properly distributed momentum value by choosing a random
*                 number between 0 and 1 and finding where it lies in the
*                 spectrum.
*
*               The function used in this subroutine was derived from
*                  the data in:
*                    Physics Data 25-1 "Cosmic Rays on Earth" (1984)
*                     (ISBN 0344-8401)   O.C. Allkofer and P.K.F. Grieder
*                     Table 2.2.3              page 241
*                  and is a curve fit based upon the data in the table.  The
*                  curve is valid between 0.5 and 1000 GeV/c.
*
*..Author     : Michael P Bringle
*..History    :
*
*
******-****************************************************************-*******
      Implicit None
*..Input:
        Real*8 X
*..Local:
      Real*8 a,b,c,d,e
*=====-================================================================-=======
      a=3.095047d-3
      b=2.187907d0
      c=0.205373d0
      d=0.462988d0
      e=0.17077d0

      PhysDatSpec=a*x**(-3.0d0 + (b/(c*x+d)**e))

      Return
      End


      Subroutine GetP3(P_Min,P_Max,Rnum,Theta_R,P)

******-****************************************************************-*******
*
*..Description: Dar Inv Spectrum
*                 This method uses data tabulated in bins of integrated
*                   momentum vs angle theta.  That data for this subroutine
*                   comes from a BASIC program(darinv3.exe) that was given
*                   to us by Prof. Henry W. Kendall at Massachusetts Institute
*                   of Technology.  It calculates the momentum spectrum vs
*                   zenith angle using Arnon Dar's data from:
*                     Cosmic Ray Muons at Gound Level and Deep Underground
*                     By Arnon Dar               1984
*                     TECHNION-PHYS-84-41
*                     Department of Physics
*                     Technion _ Israel Institute of Technology
*                     Haifa, Israel
*                   The program (darinv3.exe) states that if the zenith angle
*                   is chosen with a Cos^2 this will yield the total momentum
*                   energy distribution for the given zenith angle.
*                   The data is tabulated in 1000 bins in each of 19 files.
*                   The files are tabulated for every 5 degrees of zenith
*                   angle.  This program takes the four bins nearest to each
*                   random number that are in the four theta files nearest
*                   to the current Theta_D value.  It then separately
*                   interpolates each Theta curve to find the momentum value
*                   on that curve then interpolates the four Theta curves
*                   to find the momentum of the event.  Basically it is a
*                   two dimensional interpolation
*                 The valid range for momentum values is between 1 and
*                   1000 GeV/c.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*               Last-bin anomaly fixed, August 23, 1994.
*
******-****************************************************************-*******
      Implicit None
*..Input:
      Real*4 P_Min,P_Max        !Momentum cutoffs(GeV/c)
      Real*4 Rnum               !Random Number
      Real*4 Theta_R            !Theta_E in Radians
*..Output:
      Real*4 P                  !Momentum value(GeV/c)
*..Local:
      Integer Bins              !Number of bins in each file
      Parameter (Bins=1001)

      Integer L1,L2             !Do loop variables
      Integer Th_Indx,Rn_Indx   !Bin indexes for Theta and Rnum bins
      Integer Flag              !Initialization flag
      Real*4 Theta_D            !Theta_E in degrees
      Real*4 ISpc(Bins,0:18)    !Integrated Spectrum value for a given bin
                                !  number and theta index.
      Real*4 BNum(Bins,0:18)    !Bin number for random number indexing that
                                !  corresponds to ISpc(#1,#2)
      Real*4 X4(4),Y4(4)        !Interpolation points
      Real*4 T4(4),P4(4)        !Interpolation points

      Real*4 Rad_Deg            !Constant
      Parameter (Rad_Deg=180.0d0/3.14159265358979323846264338328)
*..Funct:
      Real*4 IntrPol            !An interpolation Routine
      External IntrPol
*..Saved:
      Save ISpc,BNum            !Save integrated spectrum
      Save Flag                 !Save so initialization routine done only once
*=====-================================================================-=======
*First Put Theta in degrees
      Theta_D=Acos(Abs(Cos(Theta_R)))*Rad_Deg !We want theta in degrees

*Second find the bin associated with our random number
      Rnum=Rnum*(Bins-1)+1  !We want Rnum 1 -> 1000

*Get the integrated Spectrum
      If (Flag.eq.0) then
        Call GetP3Ini(Bins,P_Min,P_Max,BNum,ISpc)
        Flag=1
      Endif

*Take Theta value and find which theta bins it lies between
      Th_Indx=Int(Theta_D/5.0d0)-2
      If (Th_Indx.lt.-1) Th_Indx=-1
      If (Th_Indx.gt.14) Th_Indx=14

      
*Take Random number and find which P bin it lies between taking into
*    account the theta bins just selected
      Do L1=1, 4
        Rn_Indx=Int(Rnum)-2
        If (Rn_Indx.lt.0) Rn_Indx=0
        If (Rn_Indx.gt.Bins-4) Rn_Indx=Bins-4
        If (RNum.gt.999.and.RNum.lt.1000) RN_Indx=RN_Indx-1
        Do L2=1,4
          X4(L2)=BNum(Rn_Indx+L2,Th_Indx+L1)
          Y4(L2)=ISPC(Rn_Indx+L2,Th_Indx+L1)
        Enddo

        P4(L1)=IntrPol(3,RNum,X4,Y4)
        T4(L1)=(Th_Indx+L1)*5.0
      Enddo

      P=IntrPol(3,Theta_D,T4,P4)

      Return
      End


      Subroutine GetP3Ini(Bins,P_Min,P_Max,BNum,ISpc)

******-****************************************************************-*******
*
*..Description: DarInv Momentum Spectrum
*                 This method uses data tabulated in bins of integrated
*                   momentum vs angle theta.  That data for this subroutine
*                   comes from a BASIC program(darinv3.exe) that was given
*                   to us by Prof. Henry W. Kendall at Massachusetts Institute
*                   of Technology.  It calculates the momentum spectrum vs
*                   zenith angle using Arnon Dar's data:
*                     Cosmic Ray Muons at Gound Level and Deep Underground
*                     By Arnon Dar               1984
*                     TECHNION-PHYS-84-41
*                     Department of Physics
*                     Technion _ Israel Institute of Technology
*                     Haifa, Israel
*                   The program (darinv3.exe) states that if the zenith angle
*                   is chossen with a Cos^2 this will yield the total momentum
*                   energy distribution for the given zenith angle.
*                   The data is read into arrays that are indexed by bin
*                   number and theta angle.  With the array have 1000 bins,
*                   which contain the integrated momentum spectrum, by 19
*                   theta numbers(in 5 degree steps.)
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History:
*    Bringle    8/5/94  - First Version.
*    Kral       8/8/94  - Changed file names from thetaxx.dat to cosmdarxx.dat.
*                       - Logical unit is 78 (for time being).
*    Bringle    8/18/94 - Changed file names from cosmdarxx.dat to cosdarxx.dat
*                         to allow for use on machines with 8 character
*                         filename limits(IBM PC's)
*
******-****************************************************************-*******
      Implicit None
*..Input:
      Integer Bins              !Number of bins
      Real*4  P_Min,P_Max       !Momentum Cutoff(GeV/c)
*..Output:
      Real*4 ISpc(Bins,0:18)    !Integrated Spectrum value for a given bin
                                !  number and theta index.
      Real*4 BNum(Bins,0:18)    !Bin number for random number indexing that
                                !  corresponds to ISpc(#1,#2)
*..Local:
      Integer Error             !Indicates Error on file read
      Integer L1,L2             !Do Loop variables
      Integer Iluncosm          !Logical unit

*..File Names (one for every 5 degrees in theta):
      Character*40 File_Name(0:18) /  'cosdar00.dat','cosdar05.dat',
     +  'cosdar10.dat','cosdar15.dat','cosdar20.dat','cosdar25.dat',
     +  'cosdar30.dat','cosdar35.dat','cosdar40.dat','cosdar45.dat',
     +  'cosdar50.dat','cosdar55.dat','cosdar60.dat','cosdar65.dat',
     +  'cosdar70.dat','cosdar75.dat','cosdar80.dat','cosdar85.dat',
     +  'cosdar90.dat'/
*=====-================================================================-=======
*Initialization Routine
*     Report Information to user
 10   Format(A,F13.4,A)
      Write(*,*)"Momentum Spectrum Data:"
      Write(*,*)"   The method uses look up tables to produce a"
      Write(*,*)"   momentum spectrum that has a zenith angle"
      Write(*,*)"   dependence.  The values are taken from DarInv,"
      Write(*,*)"   inverted integral flux from Arnon Dar,"
      Write(*,*)
     +  "   Cosmic Ray Muons at Ground Level and Deep Underground,"
      Write(*,*)"   TECHNION-PHYS-84-41, Haifa, Israel (1984)."
      Write(*,*)"   Flux is valid between 1 and 1000 GeV/c."

      If (P_Max.gt.1000.0.or.P_Min.lt.1.0) then
        Write(*,*)
        Write(*,*)"GetP3Ini:  Error - program terminating"
        Write(*,*)"   *************** WARNING ******************"
        Write(*,*)"   The Momentum limits you have chosen are"
        Write(*,*)"   outside the allowed range(see above.)"
        Write(*,10)"     P_Min allowed:    1.0 GeV/c   Yours: ",
     +            P_Min," GeV/c"
        Write(*,10)"     P_Max allowed: 1000.0 GeV/c   Yours: ",
     +            P_Max," GeV/c"
        Write(*,*)"   Program terminating..."
        Write(*,*)
        Stop  !Stopped in GetP3Ini - momentum limits
      Endif

*     Read in Spectrum
      Iluncosm = 78             ! Hardwired logical unit number
      Do L1=0, 18
        Open(Iluncosm,File=File_Name(L1),STATUS='OLD',IOSTAT=Error,
     +    BLANK='NULL')

        If(Error.ne.0) then
          Write(*,*)
          Write(*,*)"GetP3Ini:  Error - program terminating"
          Write(*,*)" Fatal error in reading in DarInv files:"
          Write(*,*)" File Error in ",File_Name(L1)
          Write(*,*)" Does file or link exist in working directory?"
          Write(*,*)" Program terminating..."
          Write(*,*)
          Stop  !Stopped in GetP3Ini - file read error
        Endif
        Do L2=1, Bins
          BNum(L2,L1)=L2        !Random number Bin
          Read(Iluncosm,*) ISPC(L2,L1) !Momentum Value(Integrated Spectrum)
        Enddo

        Close (Iluncosm,Status='KEEP')

      Enddo

      Write(*,*)

      Return
      End


      Subroutine GetCharge(Rnum,Charge)

******-****************************************************************-*******
*
*..Description: This program gets a randomly distributed Muon charge.
*                 Each Muon can have either a +1 or a -1 charge.
*                 Using the ratio of the two charges a simple selection
*                 method.  The ratio of + to - was taken from:
*                      Particle Properties Data Booklet  June 1992
*                      American Institute of Physics         P.129
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 Rnum               !Random number
*..Output:
      Integer Charge            !Muon Charge
*..Local
      Real*4 Ratio              !Ratio of + to Total[(-) + (+)]
      Parameter (Ratio=5.0/(4.0+5.0))
*=====-================================================================-=======
      If (Rnum.le.Ratio) Charge=1
      If (Rnum.gt.Ratio) Charge=-1

 999  Return
      End


      Subroutine Locate(D_Mode,HemiRad,Det_Len,Det_Rad,Det_Off
     +                  ,Theta_E,Phi_E,Theta,Phi,X,Y,Z,Flag)

******-****************************************************************-*******
*
*..Description: Incident Muon Spectrum Simulation
*                  This program calls the programs that place the
*                  particle on the cylinder enclosing the detector.
*                  If a mode 21 or 22 call is made the program just
*                  reformats the data in the desired manner.
*                  The program works by:
*                   i) Shadow:
*                      Taking an incoming particle with a given theta
*                      and phi angle and it places it on the X-Y plane
*                      in a region that would cause it to hit the hemisphere
*                      enclosing the detector.
*                  ii) Sphere:
*                      Next the program takes the particle and tracks it
*                      from the plane to where it would have hit on the
*                      hemisphere.
*                 iii) Frame:
*                      Translate Theta, Phi and the impact points to
*                      the detectors frame of reference.
*                  iv) Cylinder:
*                      Place the particle from the hemisphere onto the
*                      cylinder.  If the particle misses set flag=1 so
*                      a new event can be generated.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input :
      Integer D_Mode           ! Data output mode
      Real*4 Theta_E,Phi_E     ! Momentum magnitude and direction(E Frame)
      Real*4 HemiRad           ! Enclosing hemisphere radius(cm)
      Real*4 Det_Len           ! Mother Cylinder Length(cm)
      Real*4 Det_Rad           ! Mother Cylinder Radius(cm)
      Real*4 Det_Off           ! Mother Cylinder Offset from ground(cm)
*..Output:
      Integer Flag             ! Indicates if particle hit cylinder
      Real*4 Theta,Phi         ! Momentum magnitude and direction(radians)
      Real*4 X,Y,Z             ! Particle location(cm)
*..Local :
      Real*4 Theta_D,Phi_D     ! Momentum magnitude and direction(D Frame)
      Real*4 X_D,Y_D,Z_D       ! Particle location(Detector Frame)
      Real*4 X_S,Y_S           ! Where particle hits Shadow Plane
      Real*4 X_E,Y_E,Z_E       ! Particle Location on hemisphere(Muon Frame)
      Real*4 X_C,Y_C,Z_C       ! Particle Location on cylinder(Muon Frame)
*..Saved :
      Save X_E,Y_E,Z_E         ! Save incase mode 21 call
      Save X_S,Y_S             ! Save incase mode 22 call
*=====-================================================================-=======
* Standard create event and report it in the Detector Frame
      If(D_Mode.eq.2) then
*       Throw the particle at the Spherical Shadow
        Call Shadow(D_Mode,0.0,0.0,HemiRad,Theta_E,Phi_E,X_S,Y_S)

*       Find where the particle hit on the hemisphere
        Call Sphere(HemiRad,Theta_E,Phi_E,X_S,Y_S,X_E,Y_E,Z_E)

*       Transformation to the detector coordinate system
        Call Frame(Det_Off,Det_Rad,Theta_E,Phi_E,X_E,Y_E,Z_E,
     +              Theta_D,Phi_D,X_D,Y_D,Z_D)

*       Place onto mother volume cylinder
        Call Cylinder(Det_Len,Det_Rad,Theta_D,Phi_D,X_D,Y_D,Z_D,
     +                Flag,X_C,Y_C,Z_C)

*       Swap Cylinder impact data into return variables
        X=X_C
        Y=Y_C
        Z=Z_C
        Theta=Theta_D
        Phi  =Phi_D

*     Reformat the data and report it in the Earth Frame
      Elseif(D_Mode.eq.21) then
*       Swap Earth Frame data into return variables
        X=X_E
        Y=Y_E
        Z=Z_E
        Theta=Theta_E
        Phi  =Phi_E

*     Reformat the data and report it in the Earth Frame
      Elseif(D_Mode.eq.22) then
*       Swap Shadow Frame data into return variables
*       Transform back to shadow  coordinates(Rotate to proper Phi angle)
        X=Sqrt(X_S*X_S+Y_S*Y_S)*Cos(ATan2(Y_S,X_S)-Phi_E)
        Y=Sqrt(X_S*X_S+Y_S*Y_S)*Sin(ATan2(Y_S,X_S)-Phi_E)

* Standard create event and report it in the Detector Frame but
*    only track it to the hemisphere.
      Elseif(D_Mode.eq.23) then
*       Throw the particle at the Spherical Shadow
        Call Shadow(D_Mode,0.0,0.0,HemiRad,Theta_E,Phi_E,X_S,Y_S)

*       Find where the particle hit on the hemisphere
        Call Sphere(HemiRad,Theta_E,Phi_E,X_S,Y_S,X_E,Y_E,Z_E)

*       Transformation to the detector coordinate system
        Call Frame(Det_Off,Det_Rad,Theta_E,Phi_E,X_E,Y_E,Z_E,
     +              Theta_D,Phi_D,X_D,Y_D,Z_D)

*       Swap Hemisphere impact data into return variables
        X=X_D
        Y=Y_D
        Z=Z_D
        Theta=Theta_D
        Phi  =Phi_D
      Endif

 999  Return
      End


      Subroutine Shadow(D_Mode,Flux,Rate,HemiRad,
     +                  Theta_E,Phi_E,X_S,Y_S)

******-****************************************************************-*******
*
*..Description: Spherical Shadow Impact Point Routine
*                 This program determines where a particle with a given
*                   Phi and Theta will impact on a sphere by assigning
*                   equal probability to each area element of the sphere's
*                   'shadow' and then picking one of these elements.
*                  How it works:
*                    We can assume that for a given angle of theta, the
*                    flux over a flat plane will be equal at every point.
*                    From this one takes the incident zenith angle and
*                    determines what area on the plane corresponds to the
*                    exposed area on the hemisphere.  This is analogous
*                    to the shadow cast by shining a very distant light source
*                    on a hemisphere.  Two random coordinates are chosen until
*                    the points lie in the shadow.  This is all done assuming
*                    phi=0 for simplicity.  Then the final result is rotated
*                    to the proper phi angle.
*                  Because part of our shadow lies outside the hemisphere
*                    base we have to account for that in our flux rates.
*                    For differing theta distributions this will give
*                    differing effects(set GetTheta1.)  Shadow provides a
*                    test to see if the current weight factor gives a
*                    reasonable weight factor.  This is done by summing up
*                    the inverse area elements for each particle and at the
*                    end multiplying by the flux.  This value should be close
*                    to the predicted flux rate (Rate).  The ratio shows
*                    how close they are to each other.
*                  The D_Mode variable here tells the program which mode to
*                    execute.
*                    D_Mode=1 Initialize the variables
*                    D_MOde=3 Report the Data
*
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input  :
      Integer D_Mode         ! Used to select Hit/Miss report options
      Real*4 Theta_E,Phi_E   ! Muon Trajectory(radians)
      Real*4 HemiRad         ! Sphere Radius(cm)
      Real*4 Flux0,Flux      ! Flux rate greater than 1.0GeV/c (Part/cm^2 sec)
      Real*4 Rate_PD,Rate    ! Rate Predicted(Hz)
*..Output :
      Real*4 X_S,Y_S         ! Shadow impact point(cm)
*..Local  :
      Integer L1             ! Do loop variable
      Integer Flag           ! Within shadow flag
      Integer LooP_Max       ! Maximum number of trials to pick point
      Integer Tries,Hits     ! Hit/Miss shadow variables
      Real*4 Ratio,dRatio    ! Predicted Ratio and expected error
      Real*4 Rate_MC         ! Flux Rate calculated

      Real*4 GetRan          ! Random Number Generator Alias
      Real*4 ICosTh          ! Inverse Cosine Theta
      Real*4 X_Max,Y_Max     ! Box Enclosing the Shadow(cm)
      Real*4 X_hit,Y_Hit     ! Randomly chosen impact points(cm)
      Real*4 X_Chk           ! Shadow Boundary for Y_hit point(cm)
      Real*4 Pi              ! Constant

      Real*8 IArea           ! Weighted Rate Counter
      Save Tries,Hits,IArea

      Parameter (Pi=3.14159265358979323846264338328)
      Parameter(LooP_Max=1000) !Maximum Trials to pick point
*=====-================================================================-=======

*Counter Routines
 10   Format(A,F15.6,A)
 11   Format(A,I8)
      If (D_Mode.eq.2.or.D_Mode.eq.23) then     !Increment Counter
        Hits=Hits+1
        IArea=IArea+(1.0d0/(0.5*Pi*HemiRad**2*
     +              (1.0d0+1.0d0/Abs(Cos(Theta_E)))))

      Elseif (D_Mode.eq.1) then !Reset Counters
        Tries=0
        Hits =0
        IArea=0.0d0
        Flux0=Flux
        Rate_PD=Rate
        Goto 999

      Elseif (D_Mode.eq.3) then !Report Information
        Ratio=Real(Hits)/Real(Tries)
        dRatio=SQRT( (Ratio)*(1.0-Ratio)/Tries)
        Write(*,*)

        Write(*,11)"HemiCosm/Shadow Hit Statistics:"
        Write(*,11)"   Total Hits :",Hits
        Write(*,11)"   Total Tries:",Tries
        Write(*,10)"   Ratio      :",Ratio
        Write(*,10)"   dRatio     :",dRatio
        Write(*,10)"   Ratio Test :",Ratio/(Pi/4.0)
        Write(*,10)"   Ratio dTest:",dRatio/(Pi/4.0)
        Write(*,*)"Flux Testing:"
        Write(*,10)"   Sum 1/Area       : ",IArea," (cm^-2)"
        Write(*,10)"   Rate(Predicted  ): ",Rate_PD," (Hz)"
        Rate_MC=Flux0*Hits/IArea
        Write(*,10)"   Rate(Monte Carlo): ",Rate_MC," (Hz)"
        Write(*,10)"   Ratio(MC/Predict): ",Rate_MC/Rate_PD
        Write(*,*)"     The ratio would not be 1 if the Cos^2(Theta)"
        Write(*,*)"      distribution has been changed without"
        Write(*,*)"      updating the predicted rate"
        Write(*,*)
        Goto 999
      Endif

*Predefined Variables
      If (Theta_E.lt.(Pi/2.0+1.0e-8)) Theta_E=Pi/2.0+1.0e-8
      ICosTh=1.0d0/Abs(Cos(Theta_E))

*Determine Shadow Boundary(Enclosing Rectangle)
      X_Max=HemiRad*ICosTh+HemiRad !Shadow length+Radius
      Y_Max=2.0*HemiRad            !Diameter of hemisphere

      Do L1=1, LooP_Max   ! To prevent infinite loop
        Flag=1            ! Set flag to failure
        Tries=Tries+1     ! Increment number of tries

*       Pick a Point in the rectangle
        X_Hit=GetRan(L1)*X_Max-HemiRad       !-HemiRad to X_Max
        Y_Hit=(GetRan(L1+1000)-0.5)*Y_Max    !-(Y_Max/2) to (Ymax/2)

*Check to see if it lies in shadow(Note this assumes that the Sphere
*       is centered on (0,0,0)(Earth))

*       Extended Part of Shadow
        If (X_Hit.gt.0) then
*         Shadow boundary Equation for a given Y and Theta_E
          X_Chk=( SQRT(HemiRad*HemiRad-Y_Hit*Y_Hit)*ICosTh)

          If(Abs(X_Hit).le.X_Chk) then !If within shadow then keep it
            Flag=0              !Set flag to success
            Goto 111            !Exit loop
          Endif
        Endif

**Standard Circle part of shadow
        If (X_Hit.le.0) then
*         Shadow boundary Equation for a given Y and Theta_E
          X_Chk=( SQRT(HemiRad*HemiRad - Y_Hit*Y_Hit) )

          If(Abs(X_Hit).le.X_Chk) then !If within shadow then keep it
            Flag=0              !Set flag to success
            Goto 111            !Exit loop
          Endif
        Endif
      Enddo

 111  If (Flag.eq.1) then
        Write(*,*)
        Write(*,*)"Shadow:  Error - Shadow pick loops exceeded!"
        Write(*,*)"   The Program was configured to throw ",LooP_Max
        Write(*,*)"   trials at the 'impact shadow' of the figure."
        Write(*,*)"   This limit prevented an infinite loop.  The"
        Write(*,*)"   program has exceed this limit.  You can try"
        Write(*,*)"   increasing this limit,LooP_Max Variable in the"
        Write(*,*)"   shadow subroutine.  The value of the input"
        Write(*,*)"   data at the time of program termination was:"
        Write(*,10)"   Theta   : ",Theta_E
        Write(*,10)"   Phi     : ",Phi_E
        Write(*,10)"   HemiRad : ",HemiRad
        Write(*,10)"   X_S     : ",X_S
        Write(*,10)"   Y_S     : ",Y_S
        Write(*,*)"   Program terminating..."
        Write(*,*)
        Stop !Stopped in Shadow
      Endif

*Transform to proper coordinates(Rotate to proper Phi angle)
      X_S=Sqrt(X_Hit*X_Hit+Y_Hit*Y_Hit)*Cos(ATan2(Y_Hit,X_Hit)+Phi_E)
      Y_S=Sqrt(X_Hit*X_Hit+Y_Hit*Y_Hit)*Sin(ATan2(Y_Hit,X_Hit)+Phi_E)

999   Return
      End


      Subroutine Sphere(HemiRad,Theta_E,Phi_E,X_S,Y_S,X_E,Y_E,Z_E)

******-****************************************************************-*******
*
*..Description: This program determines where the particle hit on the
*                 hemisphere from where it hit on the shadow.  This is done
*                 by substituting the parametric equation for a line into
*                 the equation for a sphere and then solving for the
*                 parameter variable.  This variable is then put back into
*                 the equation for a line and the sphere coordinates result.
*                 Because we should be guaranteed a hit on the hemisphere
*                 a square root error should not result.  That would mean
*                 that the particle missed the sphere.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 Theta_E, Phi_E    !Incident Trajectory
      Real*4 X_S,Y_S           !Shadow Impact points
      Real*4 HemiRad           !Sphere Radius
*..Output:
      Real*4 X_E,Y_E,Z_E       !Sphere Impact Points
*..Local:
      Real*4 CosTh,SinTh       !Predefined variables(for speed)
      Real*4 CosPh,SinPh       !Predefined variables(for speed)
      Real*4 a,b,c             !Quadratic Equation Coefficients
      Real*4 t1,t2             !Roots of Quadratic Equation
      Real*4 X1,Y1,Z1          !Incoming Ray Solution(First Root)
      Real*4 X2,Y2,Z2          !Incoming Ray Solution(Second Root)
      Real*4 Dist1,Dist2       !Distance from impact to shadow plane
      Real*4 temp1,temp2       !Temporary variables

*=====-================================================================-=======
*Predefined Variables(for speed)
      CosTh=Cos(-Theta_E) !Here we are assuming that Theta 90->180 deg
      SinTh=Sin(-Theta_E) !  so this will give us a positive value.
      CosPh=Cos(Phi_E)
      SinPh=Sin(Phi_E)

*Find Coefficients of quadratic Equation
      a=1.0
      b=2.0*SinTh*(X_S*CosPh + Y_S*SinPh)
      c=X_S*X_S + Y_S*Y_S - HemiRad*HemiRad

*Solve Quadratic Equation
      temp1=1.0d0/(2.0*a)
      temp2=(b*b-4.0*a*c)
      if (temp2.lt.0) temp2=0.0!sometimes when subtracting the two large
                               !numbers we get a small neg so set it equal
                               !to zero in that case.
      temp2=sqrt(temp2)
      t1=((-b)+temp2)*temp1
      t2=((-b)-temp2)*temp1

*Solve for Points(Parametric Line Equations)
      X1= ( t1*SinTh*CosPh + X_S )
      Y1= ( t1*SinTh*SinPh + Y_S )
      Z1=-( t1*CosTh             )

      X2= ( t2*SinTh*CosPh + X_S )
      Y2= ( t2*SinTh*SinPh + Y_S )
      Z2=-( t2*CosTh             )

*Find Distance from Impact Point
      Dist1=Sqrt((X1-X_S)*(X1-X_S)+(Y1-Y_S)*(Y1-Y_S)+(Z1)*(Z1))
      Dist2=Sqrt((X2-X_S)*(X2-X_S)+(Y2-Y_S)*(Y2-Y_S)+(Z2)*(Z2))

*Assign Furthest one from shadow to coordinates(Z>=0 of course)
      If ( (Dist2.gt.Dist1.and.z2.ge.0).or.(Z1.lt.0)) then
        X_E=X2
        Y_E=Y2
        Z_E=Z2
      Else
        X_E=X1
        Y_E=Y1
        Z_E=Z1
      Endif

999   Return
      End


      Subroutine Frame(Det_Off,Det_Rad,Theta_E,Phi_E,X_E,Y_E,Z_E,
     +                 Theta_D,Phi_D,X_D,Y_D,Z_D)

******-****************************************************************-*******
*
*..Description: This program translates the Earth coordinate frame to the
*                 detector coordinate frame.  The frames differ in that the
*                 Earth frame is centered on the ground with the Z axis
*                 pointing 'up'.  The Detector frame however, has Y pointing
*                 up and is centered in the middle of the detector with
*                 the Z axis running along its axis of symmetry.  The
*                 geometry is as follows:
*
*               Earth                              Detector
*                +Z                                   +Y
*                 ^                                    ^
*                 |                                    |
*                 |                                    |
*                 |                                    |
*                 |                                    |
*                 |                                    |
*                 |                   ___   -----------*--------------->+Z
*                 |                    |                \
*                 |             Det_Rad+Det_Off          \
*                 |                    |                  \
*  ---------------*-------------->+X  ---                  \
*                  \                                        \
*                   \                                        \
*                    \                                        \
*                     \                                        \
*                      \                                       -X
*                       \
*                        \
*                        -Y
*
*       Phi  :  Between X and Y                            Between  X and Y
*       Theta:     Off of Z                                     Off of Z
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 Det_Rad         ! Mother Cylinder Radius(cm)
      Real*4 Det_Off         ! Mother Cylinder Offset from ground(cm)
      Real*4 Theta_E, Phi_E  ! Particle Trajectory(Earth)
      Real*4 X_E,Y_E,Z_E     ! Impact Position(Earth)
*..Output:
      Real*4 Theta_D, Phi_D  ! Particle Trajectory(Detector)
      Real*4 X_D, Y_D, Z_D   ! Particle Trajectory(Detector)
*..Local:
      Real*4 a,b,c           !Vector coefficients
      Real*4 Pi,Pi_O2        !Constants
      Parameter (Pi=3.14159265358979323846264338328)
      Parameter (Pi_O2=Pi/2.0d0)
*=====-================================================================-=======
*Translate position coordinates
      X_D=Y_E
      Y_D=Z_E+Det_Off-Det_Rad
      Z_D=X_E

*Translate trajectory vectors
*    Construct Unit vector in Muon frame
      a=Cos(Phi_E)*Sin(Theta_E)
      b=Sin(Phi_E)*Sin(Theta_E)
      c=           Cos(Theta_E)

*     Translate it over to detector frame
      Phi_D=Atan2(c,b)
      Theta_D=Atan2(Sqrt(b*b+c*c),a)

 999  Return
      End


      Subroutine Cylinder(Det_Len,Det_Rad,Theta_D,Phi_D,X_D,Y_D,Z_D,
     +                 Flag,X_C,Y_C,Z_C)

******-****************************************************************-*******
*
*..Description: Find where particle hits the cylinder.
*                The program traces the particle from a given point
*                 in space(the hemisphere) to the cylinder that just
*                 encloses the detector.  It does this in the following
*                 order:
*                   i) Check to see if it hit the end caps
*                  ii) Check to see if it hit the body
*                 The end caps are checked first because particle which
*                 are traveling nearly axial will give the body hit part
*                 a bad time.  The method for each process is described
*                 before each process.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Real*4 Det_Len           ! Mother Cylinder Length(cm)
      Real*4 Det_Rad           ! Mother Cylinder Radius(cm)
      Real*4 X_D,Y_D,Z_D       ! Hemisphere impact points(Detector)
      Real*4 Theta_D,Phi_D     ! Particle Trajectory(Detector)
*..Output:
      Integer Flag             ! Did hit occur(1=yes)
      Real*4 X_C,Y_C,Z_C       ! Cylinder Impact points(Detector)
*..Local:
      Real*4 tparm1,tparm2,t
      Real*4 a,b,c
      Real*4 Z_max,Z_min       ! Cylinder Extents
      Real*4 Z_test,Test_solution
      Real*4 X_Solution1,X_Solution2
      Real*4 Y_Solution1,Y_Solution2
      Real*4 Rad_1, Rad_2
      Real*4 Pi,Pi_O2           !Constants
      Real*4 CosTh,SinTh        !Stored for speed
      Real*4 SinPh,CosPh        !Stored for speed

      Parameter (Pi=3.14159265358979323846264338328)
      Parameter (Pi_O2=Pi/2.0d0)
*=====-================================================================-=======
*Setup variables
      CosTh=Cos(Theta_D)
      SinTh=Sin(Theta_D)
      CosPh=Cos(Phi_D)
      SinPh=Sin(Phi_D)
      Flag=0
      t=0
      Z_Max= Det_Len/2.0 !Z Maximum of cylinder
      Z_Min=-Det_Len/2.0 !Z Minimum of cylinder

*First:
*  Check to see if it hits the two end caps.  Checking this first will help
*  remove some problems that might occur if the particle is traveling
*  down the barrel.  To find if it hits the barrel set the Z coordinate to
*  the ends of the barrel.  Then solve for t_parm at this z.  Plug this into
*  the parameterized version of the line to find the X & Y coordinates.  Then
*  check to see if they reside within the radius of the end cap.  If it does
*  then check to see if Phi is at the proper orientation that this would be
*  the first place struck.
      If (CosTh.ne.0.0) then
        tparm1=(Z_max-Z_D)/cosTh
        tparm2=(Z_min-Z_D)/CosTh

        Rad_1=Sqrt( (tparm1*CosPh*SinTh+X_D)**2 +
     +              (tparm1*SinPh*SinTh+Y_D)**2 )

        Rad_2=Sqrt( (tparm2*CosPh*SinTh+X_D)**2 +
     +              (tparm2*SinPh*SinTh+Y_D)**2 )

        If (Rad_1.le.Det_Rad.and.CosTh.lt.0) then
          Flag=1
          t=tparm1
          X_C=t*CosPh*SinTh + X_D
          Y_C=t*SinPh*SinTh + Y_D
          Z_C=t*           CosTh + Z_D
        Endif

        If (Rad_2.le.Det_Rad.and.CosTh.gt.0) then
          Flag=1
          t=tparm2
          X_C=t*CosPh*SinTh + X_D
          Y_C=t*SinPh*SinTh + Y_D
          Z_C=t*           CosTh + Z_D
        Endif
      Endif

*Second:
*  Project the cylinder onto the X_Z plane and find where the particle
*  would hit this circle.  Then using these Y and Z points pick the one
*  with the highest Y(vertical) value because our particle will always
*  be slanting downwards.  Then using these X & Z points plug one back into
*  the equation for the line(parameterized) to fine the proper 't' value.
*  Both variables should give the same t value.  Then plug this into the
*  Parameterized equation for z and check to see if our cylinder walls enclose
*  this value for z.

* X=t_P * CosPh + X_D
* Y=t_P * SinPh + Y_D
* X^2 + Y^2= Det_Rad
      If (Flag.eq.0) then
        a=1.0
        b=2.0*(CosPh*X_D + SinPh*Y_D)
        c=(X_D*X_D)+(Y_D*Y_D)-Det_Rad*Det_Rad

        Test_Solution=b*b-4.0*a*c

        If(Test_solution.ge.0) then !does it even hit the circle

          tparm1=(-b + Sqrt(b*b - 4.0*a*c))/2.0*a
          tparm2=(-b - Sqrt(b*b - 4.0*a*c))/2.0*a

          If (SinPh.ne.0.0) then
            Y_Solution1=tparm1*SinPh+Y_D
            Y_Solution2=tparm2*SinPh+Y_D
            If (Y_Solution1.gt.Y_Solution2) then
              t=(Y_Solution1-Y_D)/(SinTh*SinPh)
            Else
              t=(Y_Solution2-Y_D)/(SinTh*SinPh)

            Endif

          Else
            X_solution1=tparm1*CosPh+X_D
            X_solution2=tparm2*CosPh+X_D

            If (X_solution1.gt.X_solution2.and.Phi_D.lt.Pi_O2) then
              t=(X_Solution1-X_D)/(SinTh*CosPh)
            Else
              t=(X_Solution2-X_D)/(SinTh*CosPh)
            Endif
          Endif

          Z_test=t*CosTh+Z_D

          If (Z_test.le.Z_Max.and.Z_test.ge.Z_Min.and.t.ge.0) then
            Flag=1              !Signal: Found it
            X_C=t*CosPh*SinTh + X_D
            Y_C=t*SinPh*SinTh + Y_D
            Z_C=t*CosTh            + Z_D
          Endif
        endif
      Endif

      If(t.lt.0) then           !Check that it is going the right way
        Flag=-1                 !Set to -1 if traveled the wrong way
      Endif

 999  Return
      End


      Subroutine FinishInfo(Trial,Runs,dT,Rate_PD)

******-****************************************************************-*******
*
*..Description: This program prints some data for the user.
*
*..Author     : Michael P Bringle
*               Lawrence Berkeley Lab
*
*..History    : First Version August 5, 1994
*               Second Version August 19, 1994
*               Selected rate, August 22, 1994
*               Formatting Xrate, August 25, 1994
*
******-****************************************************************-*******
      Implicit none
*..Input:
      Integer Trial             ! Number of Particles hitting the detector
      Integer Runs              ! Number of particles hitting the hemisphere
      Real*4 dT                 ! Time step for each particle
      Real*4 Rate_PD            ! Rate of particle hitting cylinder
*..Local:
      Real*4 Xrate              ! Multiplication factor to get selected rate

*=====-================================================================-=======
 10   Format(A,F13.4,A)
 11   Format(A,I8)
 12   Format(A,F18.7,A)
 13   Format(A,G18.7,A)

      Call Shadow(3,0.0,0.0,0.0,0.0,0.0,0.0,0.0)

      Rate_PD=(Real(Trial)/Real(Runs))/dT ! Rate on detector
      Xrate = 1.0/Real(Runs)/dT ! Multiplication factor to get selected rate

      Write(*,*)
      Write(*,*)"HemiCosm/FinishInfo Final Report:"
      Write(*,10)"   Effective Flux      : ",Rate_PD," Hz"
      Write(*,10)"   Hit/total Ratio     : ",Real(Trial)/Real(Runs)
      Write(*,10)"   Time Span           : ",dT*Runs," Seconds"
      Write(*,11)"   Total Particles     : ",Runs
      Write(*,11)"   Total Detector Hits : ",Trial
      Write(*,*)
      Write(*,*)
     +  "To find selected rate (where you define selection criteria),"
      Write(*,*)
     +  " Selected Rate =  [N_Selected/Total Particles]*Rate(Predicted)"
      If (Xrate .ge. 0.10) then
        Write(*,12)
     +    " In this run, Selected Rate = N_Selected * ",Xrate," Hz"
      Else
        Write(*,13)
     +    " In this run, Selected Rate = N_Selected * ",Xrate," Hz"
      Endif                     ! (Xrate)
      Write(*,*)

 999  Return
      End
