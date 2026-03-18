#include "StarUrQMD.h"
ClassImp(StarUrQMD);

#include "StarCallf77.h"
#include "StarGenerator/EVENT/StarGenPPEvent.h"
#include "StarGenerator/EVENT/StarGenEPEvent.h"
#include "StarGenerator/EVENT/StarGenAAEvent.h"
#include "StarGenerator/EVENT/StarGenParticle.h"

#include "StarGenerator/UTIL/StarRandom.h"
#include "StarGenerator/UTIL/StarParticleData.h"
#include <map>
#include <iostream>
#include <vector>
#include <string>

#define urqmd_init F77_NAME(urqmd_init,URQMD_INIT)
#define urqmd_make F77_NAME(urqmd_make,URQMD_MAKE)
#define urqmd_pdgid F77_NAME(urqmd_pdgid,URQMD_PDGID)

auto& pdb = StarParticleData::instance();

namespace {
  vector<std::string> tokenize(const std::string &s, const std::string &sep_chars=";")
  {
    std::string::size_type prev_pos = 0, pos = 0;
    std::vector<std::string> result;
    
    if ( s.length() > 1 ) {
      while ((pos = s.find_first_of(sep_chars, pos)) != std::string::npos) {
	result.push_back(  ( s.substr(prev_pos, pos - prev_pos))  );
	pos += 1;
	prev_pos = pos;
      }
      result.push_back( (s.substr(prev_pos)));
    }

    return result;
  }

}

// ---------------------------------------------------------------------------
/// Remap UrQMD's and Pythia's random number generator to StarRandom
extern "C" {
  Double_t ranfstar_( Int_t *idummy ){    return StarRandom::Instance().flat(); };
  Double_t pyrstar_ ( Int_t *idummy ){    return StarRandom::Instance().flat(); };  
  //  StarRandom &ranf_ = StarRandom::Instance();
  //  StarRandom &pyr_  = StarRandom::Instance();
  void urqmd_init();
  void urqmd_make();
};
Double_t rndm(){ return StarRandom::Instance().flat(); };
  
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
StarUrQMD::StarUrQMD( const Char_t *name ) : StarGenerator(name)
{


  //  assert( 2+2==5 ); // UrQMD is not ready for prime time

  /// Setting up a map between UrQMD's status codes and the HepMC status codes used in StarEvent
  //JFN 11/19/12 15:50- I can't find any documentation on UrQMD status codes (or even if the status information is stored) so we are going to do this in a very general way just to clean it up for later.
  for ( UInt_t i=0; i<200; i++)
    {
      mStatusCode[i+100] = StarGenParticle::kFinal;
    }
  //JFN 11/19/12 15:53- This next bit is for reference
  /*mStatusCode[0]   = StarGenParticle::kNull;
  mStatusCode[1]   = StarGenParticle::kFinal;
  mStatusCode[2]   = StarGenParticle::kDocumentation;
  mStatusCode[3]   = StarGenParticle::kDocumentation;*/

}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
extern "C" {
  int urqmd_pdgid( int* ityp, int* iso3 );
}


Int_t StarUrQMD::Init()
{

  if ( SAttr("FRAME") ) SetFrame( SAttr("FRAME"), DAttr("Ecms") );
  if ( SAttr("BLUE")  ) SetBlue ( SAttr("BLUE" ) );
  if ( SAttr("YELL")  ) SetYell ( SAttr("YELL" ) );


  // Proton mass:  // TODO:  Get from StarParticle DB
  Double_t ProtonMass  = 0.938272046;
  // Neutron mass:
  Double_t NeutronMass = 0.939565378;

  // Map typical species run at RHIC
  map<TString,Int_t> A, Z;
  A["p"]  =   1;    Z["p"]  =  1;
  A["n"]  =   1;    Z["n"]  =  0;
  A["d"]  =   2;    Z["d"]  =  1;
  A["He3"]=   3;    Z["He3"]=  2;
  A["Au"] = 197;    Z["Au"] = 79;
  A["Cu"] =  64;    Z["Cu"] = 29;
  A["U"]  = 238;    Z["U"]  = 92;

  A["proton"]   =1;    Z["proton"]   =1;
  A["neutron"]  =1;    Z["neutron"]  =0;
  A["deuteron"] =2;    Z["deuteron"] =1;
  A["e-"]       =0;    Z["e-"]       =0;
  A["electron"] =0;    Z["electron"] =0;
  A["e+"]       =0;    Z["e+"]       =0;
  A["positron"] =0;    Z["positron"] =0;


  TString myBlue = mBlue;
  TString myYell = mYell;

  stringstream Blue;
  stringstream Yell;

  Blue << A[myBlue] << " " << Z[myBlue];
  Yell << A[myYell] << " " << Z[myYell];

  InputParametersString["pro"] = Blue.str().c_str();
  InputParametersString["tar"] = Yell.str().c_str();

  stringstream impactParameters;
  impactParameters << mImpactMin << " " << mImpactMax;
  InputParametersString["IMP"] = impactParameters.str().c_str();

  //
  // Create a new event record for either pp or eo events
  //
  if ( ( mBlue == "proton" ) && ( mYell == "proton" ) )                 mEvent = new StarGenPPEvent();
  else                                                                  mEvent = new StarGenAAEvent();   

  /// Remapt the ROOT names to UrQMD names
  std::map< TString, string > particle;
  /// Ex: particle["ROOT name"]="UrQMD name";
  particle["proton"] = "P       ";
  particle["e-"]     = "E-      ";

  /// Set up frames
  if ( mFrame=="COM" )
    {
      InputParametersDouble["ecm"]=mRootS; //???
    }
  if ( mFrame=="3MOM" || mFrame=="4MOM" || mFrame=="5MOM" )
    {
      ///JFN 11/21/12 4:11pm- I believe this calculation of mRootS from the momentum is correct, but I wouldn't stake my life on it. Additionally, I think mRootS should be something that should be calculated by the StarGenerator framework (ie, I shouldn't have to do it)
      mRootS = ( sqrt(pow(((Z[myBlue]*ProtonMass)+((A[myBlue]-Z[myBlue])*NeutronMass)),2) + sqrt( pow(mBlueMomentum.Px(),2) + pow(mBlueMomentum.Py(),2) + pow(mBlueMomentum.Pz(),2))) + sqrt(pow(((Z[myYell]*ProtonMass)+((A[myYell]-Z[myYell])*NeutronMass)),2) + sqrt( pow(mYellMomentum.Px(),2) + pow(mYellMomentum.Py(),2) + pow(mYellMomentum.Pz(),2))));
      InputParametersDouble["ecm"]=mRootS;
    }

  std::vector<int> stablepdg = {
    111, // pi0
    211, // pi+/-
    221, // eta
    321, // K+/-
    310, // K short 
    130, // K long
    3122, // Lambda0
    3112, // Sigma -
    3222, // Sigma +
    3212, // Sigma 0
    3312, // Xi -
    3322, // Xi 0
    3334  // Omega -
  };

  // Set the pi0, pi+/- and eta states stable
  StableParticles.push_back("101");
  StableParticles.push_back("-101");
  StableParticles.push_back("102");
  StableParticles.push_back("-102");

  // Set the kaon states stable
  StableParticles.push_back("106");
  StableParticles.push_back("-106");

  // Set the lambda stable
  StableParticles.push_back("41");
  StableParticles.push_back("-41");

  // Sigma0 Sigma +/-
  StableParticles.push_back("54");
  StableParticles.push_back("-54");

  // Xi0 Xi+/-
  StableParticles.push_back("63");
  StableParticles.push_back("-63");

  // OMEGA- 3334
  StableParticles.push_back("69");
  StableParticles.push_back("-69");

  // Lambda_c
  //  StableParticles.push_back("-69");

  /// Initialize UrQMD:
  InitializeUrQMD();

  /* call */  urqmd_init();


  // Check that UrQMD codes are all available for simulation...
  std::vector<int> urqmdidvec = {

    //Neutron
    1, -1,  2112,  
    //Proton
       1,  1,  2212,
    //N*
       2, -1, 12112,        2,  1, 12212,
       3, -1,  1214,        3,  1,  2124, 
       4, -1, 22112,        4,  1, 22212,
       5, -1, 32112,        5,  1, 32212,
       6, -1,  2116,        6,  1,  2216,
       7, -1, 12116,        7,  1, 12216,
       8, -1, 21214,        8,  1, 22124,
       9, -1, 42112,        9,  1, 42212, 
      10, -1, 31214,       10,  1, 32124, 
      14, -1,  1218,       14,  1,  2128, 
      23, -1,  1218,       23,  1,  2128,
    //Delta
      24, -3,  1114,  24, -1,  2114,  24, 1,  2214,  24, 3,  2224,
      25, -3, 31114,  25, -1, 32114,  25, 1, 32214,  25, 3, 32224,
      26, -3,  1112,  26, -1,  1212,  26, 1,  2122,  26, 3,  2222,
      27, -3, 11114,  27, -1, 12114,  27, 1, 12214,  27, 3, 12224,
      28, -3, 11112,  28, -1, 11212,  28, 1, 12122,  28, 3, 12222,
      29, -3,  1116,  29, -1,  1216,  29, 1,  2126,  29, 3,  2226,
      30, -3, 21112,  30, -1, 21212,  30, 1, 22122,  30, 3, 22222,
      31, -3, 21114,  31, -1, 22114,  31, 1, 22214,  31, 3, 22224,
      32, -3, 11116,  32, -1, 11216,  32, 1, 12126,  32, 3, 12226,
      33, -3,  1118,  33, -1,  2118,  33, 1,  2218,  33, 3,  2228,
      40, -3,  1118,  40, -1,  2118,  40, 1,  2218,  40, 3,  2228,
//Lambda
      41,  0,  3122,
      42,  0, 13122,   
      43,  0,  3124,   
      44,  0, 23122,   
      45,  0, 33122,
      46,  0, 13124,
      47,  0, 43122,   
      48,  0, 53122,   
      49,  0,  3126,   
      50,  0, 13126,   
      51,  0, 23124,   
      52,  0,  3128,   
      53,  0, 23126,   
//Sigma
      54, -2,  3112,    54,  0,  3212,    54,  2,  3222,
      55, -2,  3114,    55,  0,  3214,    55,  2,  3224,
      56, -2, 13112,    56,  0, 13212,    56,  2, 13222,
      57, -2, 13114,    57,  0, 13214,    57,  2, 13224,
      58, -2, 23112,    58,  0, 23212,    58,  2, 23222,
      59, -2,  3116,    59,  0,  3216,    59,  2,  3226,
      60, -2, 13116,    60,  0, 13216,    60,  2, 13226,
      61, -2, 23114,    61,  0, 23214,    61,  2, 23224,
      62, -2,  3118,    62,  0,  3218,    62,  2,  3228,
//Xi
      63, -1,  3312,     63,  1,  3322,
      64, -1,  3314,     64,  1,  3324,
      66, -1, 13314,     66,  1, 13324,
//Omega
      69,  0,  3334,
// Lambda_C
      70,  0,  4122,
//gamma
     100,  0,    22, 
//pion
     101, -2,  -211,    101,  0,   111,    101,  2,   211,
    //eta
     102,  0,   221,
//omega
     103,  0,   223,
//rho
     104, -2,  -213,    104,  0,   113,    104,  2,   213,
//f0(980)
     105,  0, 10221,
//kaon
     106, -1,   311,    106,  1,   321,
//eta'
     107,  0,   331,
//k*(892)
     108, -1,   313,    108,  1,   323,
//phi
     109,  0,   333,
//k0*(1430)
     110, -1, 10311,    110,  1, 10321,
//a0(980)
     111, -2,-10211,    111,  0, 10111,    111,  2, 10211,
//f0(1370)
     112,  0, 20221,
//k1(1270)
     113, -1, 10313,    113,  1, 10323,
//a1(1260)
     114, -2,-20213,    114,  0, 20113,    114,  2, 20213,
//f1(1285)
     115,  0, 20223,
    //f1'(1510)
     116,  0, 40223,
//k2*(1430)
     117, -1,   315,    117,  1,   325,
//a2(1329)
     118, -2,  -215,    118,  0,   115,    118,  2,   215,
    //f2(1270)
     119,  0,   225,
//f2'(1525)
     120,  0,   335,
    //k1(1400)
     121, -1, 20313,    121,  1, 20323,
//b1
     122, -2,-10213,    122,  0, 10113,    122,  2, 10213,
//h1
     123,  0, 10223,
//K* (1410)
     125, -1, 30313,    125,  1, 30323,
//rho (1450)
     126, -2,-40213,    126,  0, 40113,    126,  2, 40213,
//omega (1420)
     127,  0, 50223,
//phi(1680)
     128,  0, 10333,
//k*(1680)
     129, -1, 40313,    129,  1, 40323,
//rho(1700)
     130, -2,-30213,    130,  0, 30113,    130,  2, 30213,
//omega(1600)
     131,  0, 60223,
//phi(1850)     
     132,  0,   337,
//D
     133,  -1,   421,    133,   1,   411,  
//D*
     134,  -1, 10421,    134,   1, 10411,  
//J/Psi
     135,  0, 443,
//Chi_c
     136,  0, 10441,
//Psi'
     137,  0, 100443,
    //Ds
     138,   0,   431,  
    //Ds*
     139,   0,   433 
  };

  int invalid = 0;
  for ( unsigned int idx=0;idx<urqmdidvec.size();idx+=3 ) {
    
    int id    = urqmdidvec[idx];
    int iso3  = urqmdidvec[idx+1];
    int pdgid = urqmdidvec[idx+2];

    auto* pdg = pdb.GetParticle( pdgid );
    if ( 0==pdg ) invalid++;
    
    LOG_INFO << "URQMD: " << id << " " 
	     << iso3 << " " 
	     << "PDG: " << pdgid << " "
	     << ((pdg)? pdg->GetName() : "unknown") << " "
	     << ((pdg)? pdg->TrackingCode() : 0) << " "
	     << endm;
    
  }

  if (invalid) {
    LOG_INFO << "Number of missing particle states=" << invalid << endm;
    assert(0);
  }

  return StMaker::Init();
}
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------

Int_t StarUrQMD::Generate()
{
  /// Generate an event
  //  GenerateEvent();

  /* call */ urqmd_make();


#if 1

  // Blue beam is a proton, running PP
  if ( ( mBlue == "proton" ) && ( mYell == "proton" ) )    FillPP( mEvent );
  // Otherwise, runnin EP
  else                                                     FillAA( mEvent );

#endif 

#if 1

  //Do Stuff with the particles
  mNumberOfParticles = sys().npart;

  LOG_INFO << "Number of particles=" << mNumberOfParticles << endm;

  //mNumberOfParticles = isys().ncoll(1); //JFN 11/28/12- this is wrong
  for ( Int_t idx=1; idx<=mNumberOfParticles; idx++ )
    {
      /* osc99 event
	i = idx
	id = pdgid(ityp(i), iso3(i))         // pdg id
	write(20,997) 
              uid(i),                        // urqmd id
              id,                            // pdgid
              0,                             // status (always 0)
     .        px(i)+ffermpx(i),              // px
              py(i)+ffermpy(i),              // py
              pz(i)+ffermpz(i),              // pz
     .        p0(i),                         // E
              fmass(i),                      // mass
x     .        frrx(i),                       // freezeout position
              frry(i), 
              frrz(i), 
              frr0(i)

       */
      int& ityp  = isys().ityp( idx );
      int& iso3  = isys().iso3( idx );
      int  id = urqmd_pdgid( &ityp, &iso3 );




      //      Int_t    id = isys().uid(idx); // or isys().itypd(idx). It isn't clear which is right


      // Not clear, but I believe that every particle is a final state particle
      Int_t    stat = StarGenParticle::kFinal; 
      // But... urqmd stores the mother particle(s) for this one.  
      Int_t    m1 = itdelay().ityptd(idx,1);
      Int_t    m2 = itdelay().ityptd(idx,2);
      // Daughter information is not stored... but in principle could be reconstructed here
      Int_t    d1 = 0; 
      Int_t    d2 = 0;
      // Momentum, energy and mass
      Double_t px = coor().px(idx);
      Double_t py = coor().py(idx);
      Double_t pz = coor().pz(idx);
      Double_t E  = coor().p0(idx);
      Double_t M  = coor().fmass(idx);
      // UrQMD does store the freezeout position... but this is well below our resolution as it is...
      Double_t vx = 0.0;
      Double_t vy = 0.0;
      Double_t vz = 0.0;
      Double_t vt = 0.0;

      //      LOG_INFO << "idx=" << idx << " id=" << id << " stat=" << stat << " m1=" << m1 << " m2=" << m2 << " d1=" << d1 << " d2=" << d2 << " px=" << px << " py=" << py << " pz=" << pz 
      //      	       << " vx=" << vx << " vy=" << vy << " vz=" << vz 
      //      	       << endm;


      // K0 to K0 short/long
      if ( id==311 || id==-311 ) {
	if ( rndm() > 0.5 ) id = 130  ; // K0 long
	else                id = 310  ; // K0 short
      }



      mEvent -> AddParticle( stat, id, m1, m2, d1, d2, px, py, pz, E, M, vx, vy, vz, vt );
    }


  mEvent->Print();

#endif

  return kStOK;
}
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
///JFN 11/18/12 13:24 - I think having a clear function is optional, and UrQMD doesn't have any explicit cleanup. Although, we could toss the un-needed output files here.
/*Int_t StarUrQMD::Clear()
{
  return kStOK;
}*/
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
void StarUrQMD::FillPP( StarGenEvent *event)
{

  // Fill the event-wise information
  StarGenPPEvent *myevent = (StarGenPPEvent *)event;
  myevent -> idBlue     = 0;
  myevent -> idYell     = 0;
  myevent -> process    = 0;
  /*myevent -> idBlue     = hwbeam().ipart1; // Int //JFN 11/26/12- "Blue beam ID". I don't know the ID convention
  myevent -> idYell     = hwbeam().ipart2;
  myevent -> process    = hwproc().iproc; //JFN 11/26/12- in principle there are process and subprocess ids because they get written in the headers of the output files, but I cant figure out where they are stored*/
  myevent -> subprocess = 0;

  myevent -> idParton1  = -999;
  myevent -> idParton2  = -999;
  myevent -> xParton1   = 0;
  myevent -> xParton2   = 0;
  myevent -> xPdf1      = -999;
  myevent -> xPdf2      = -999;
  myevent -> Q2fac      = -999;
  myevent -> Q2ren      = -999;
  myevent -> valence1   = 0;
  myevent -> valence2   = 0;

  myevent -> sHat       = 0;
  myevent -> tHat       = 0;
  myevent -> uHat       = 0;
  myevent -> ptHat      = -999;
  myevent -> thetaHat   = -999;
  myevent -> phiHat     = -999;
  
  myevent -> weight     = -999;

}


void StarUrQMD::FillAA( StarGenEvent *event)
{

  // Fill the event-wise information
  StarGenAAEvent *myevent = (StarGenAAEvent *)event;
  myevent -> idBlue     = 0;
  myevent -> idYell     = 0;
  myevent -> process    = 0;
  /*myevent -> idBlue     = hwbeam().ipart1; // Int //JFN 11/26/12- "Blue beam ID". I don't know the ID convention
  myevent -> idYell     = hwbeam().ipart2;
  myevent -> process    = hwproc().iproc; //JFN 11/26/12- in principle there are process and subprocess ids because they get written in the headers of the output files, but I cant figure out where they are stored*/
  myevent -> subprocess = 0;

  myevent -> idParton1  = -999;
  myevent -> idParton2  = -999;
  myevent -> xParton1   = 0;
  myevent -> xParton2   = 0;
  myevent -> xPdf1      = -999;
  myevent -> xPdf2      = -999;
  myevent -> Q2fac      = -999;
  myevent -> Q2ren      = -999;
  myevent -> valence1   = 0;
  myevent -> valence2   = 0;

  myevent -> sHat       = 0;
  myevent -> tHat       = 0;
  myevent -> uHat       = 0;
  myevent -> ptHat      = -999;
  myevent -> thetaHat   = -999;
  myevent -> phiHat     = -999;
  
  myevent -> weight     = -999;

}
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
///JFN 11/19/12 7:54am- to initialize UrQMD we first have to create input files for UrQMD to read and set enviormental variables, then we call UrQMD's init function
///JFN 11/19/12 7:37pm- In hindsight, I think I might rather have these functions defined in StarUrQMD.cxx.... maybe
///JFN 11/20/12 12:09pm- I have moved these two function (InitializeUrQMD and GenerateEvent) to StarUrQMD.cxx
void StarUrQMD::InitializeUrQMD()
  {
    /// Set the enviormental variables
    ///JFN 11/19/12 7:55am- Honestly, it pains me to set enviormental variables just to run this. I need to check if there is a way I can just set common block variables
    ///JFN 11/19/12 7:38pm- I would like to be able to just define the name on the input output files....
    ///JFN 11/21/12 5:19pm- I have found where the enviormental variables are loaded for UrQMD: line 359 of input.F. I see no easy way to spoof this behavior without modifying the source code, which I am trying to keep to an absolute minimum,
    //SetEnVars();

    /// Print the input file
    std::ofstream inputfile;
    inputfile.open( "UrQMD.in" );
    for(map<TString,Int_t>::iterator i = InputParametersInt.begin(); i != InputParametersInt.end(); i++)
      { 
        inputfile << i->first << " " << i->second << endl;
      }
    for(map<TString,Double_t>::iterator i = InputParametersDouble.begin(); i != InputParametersDouble.end(); i++)
      { 
        inputfile << i->first << " " << i->second << endl;
      }
    for(map<TString,TString>::iterator i = InputParametersString.begin(); i != InputParametersString.end(); i++)
      { 
        inputfile << i->first << " " << i->second << endl;
      }
    /// Set particles to not decay
    for(unsigned int i = 0; i< StableParticles.size(); i++)
      {
        inputfile << "stb " << StableParticles[i] << endl;
      }
    /// 
    if ( mFrame=="CMS" ) {
      inputfile << "ecm " << mRootS << "\n";
    }
    /// Define calculation time: total time span, interval at which output is written (both in fm/c)
    inputfile << "tim 200 200" << endl;
    /// This supresses the output files
    inputfile << "f13 \n #f14 \n f15 \n f16 \n #f17 \n #f18 \n f19 \n f20" << endl;
    /// This sets the number of events for UrQMD to run (but we aren't using this, so it doens't matter)
    inputfile << "nev 1000" << endl;
    /// This sets the random number generator seed (but we aren't using UrQMD's random number generator, so this doens't matter)
    inputfile << "rsd 111" << endl;

    
    if ( SAttr("SET" )  ) {
      for ( auto cmd : tokenize( SAttr("Set"), ";" ) ) {
	inputfile << cmd << "\n";
	LOG_INFO << "Set urqmd option:" << cmd << endm;
      }
    }


    /// This marks the end of the input file
    inputfile << "xxx and done" << endl;
    inputfile.close();

    /// Initialize UrQMD
    //    iurqmd();
  }
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
///JFN 11/19/12 7:57am- I think we can grab the event information right from the common blocks, so there may be no need to process anything, or for this function
///JFN 11/21/12 5:28pm- just for reference, the first calls for setting up output files in RunFunctin.F is at line 169. Final output is done at line 408.
/*void GenerateEvent()
  {
    /// Call the UrQMD function to make a new event
    genevt();

    /// Process the results
    //ProcessEvent();
  }*/

#if 0
#include "StarUrQMD.h"
ClassImp(StarUrQMD);

#include "StarCallf77.h"
#include "StarGenerator/EVENT/StarGenPPEvent.h"
#include "StarGenerator/EVENT/StarGenEPEvent.h"
#include "StarGenerator/EVENT/StarGenAAEvent.h"
#include "StarGenerator/EVENT/StarGenParticle.h"

#include "StarGenerator/UTIL/StarRandom.h"
#include <map>
#include <iostream>
#include <vector>
#include <string>

#define urqmd_init F77_NAME(urqmd_init,URQMD_INIT)
#define urqmd_make F77_NAME(urqmd_make,URQMD_MAKE)

namespace {
  vector<std::string> tokenize(const std::string &s, const std::string &sep_chars=";")
  {
    std::string::size_type prev_pos = 0, pos = 0;
    std::vector<std::string> result;
    
    if ( s.length() > 1 ) {
      while ((pos = s.find_first_of(sep_chars, pos)) != std::string::npos) {
	result.push_back(  ( s.substr(prev_pos, pos - prev_pos))  );
	pos += 1;
	prev_pos = pos;
      }
      result.push_back( (s.substr(prev_pos)));
    }

    return result;
  }

}

/// Remap UrQMD's and Pythia's random number generator to StarRandom
extern "C" {
  Double_t ranfstar_( Int_t *idummy ){    return StarRandom::Instance().flat(); };
  Double_t pyrstar_ ( Int_t *idummy ){    return StarRandom::Instance().flat(); };  
  //  StarRandom &ranf_ = StarRandom::Instance();
  //  StarRandom &pyr_  = StarRandom::Instance();
  void urqmd_init();
  void urqmd_make();
};

  


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
StarUrQMD::StarUrQMD( const Char_t *name ) : StarGenerator(name)
{


  //  assert( 2+2==5 ); // UrQMD is not ready for prime time

  /// Setting up a map between UrQMD's status codes and the HepMC status codes used in StarEvent
  //JFN 11/19/12 15:50- I can't find any documentation on UrQMD status codes (or even if the status information is stored) so we are going to do this in a very general way just to clean it up for later.
  for ( UInt_t i=0; i<200; i++)
    {
      mStatusCode[i+100] = StarGenParticle::kFinal;
    }
  //JFN 11/19/12 15:53- This next bit is for reference
  /*mStatusCode[0]   = StarGenParticle::kNull;
  mStatusCode[1]   = StarGenParticle::kFinal;
  mStatusCode[2]   = StarGenParticle::kDocumentation;
  mStatusCode[3]   = StarGenParticle::kDocumentation;*/

}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
Int_t StarUrQMD::Init()
{

  if ( SAttr("FRAME") ) SetFrame( SAttr("FRAME"), DAttr("Ecms") );
  if ( SAttr("BLUE")  ) SetBlue ( SAttr("BLUE" ) );
  if ( SAttr("YELL")  ) SetYell ( SAttr("YELL" ) );


  // Proton mass:  // TODO:  Get from StarParticle DB
  Double_t ProtonMass  = 0.938272046;
  // Neutron mass:
  Double_t NeutronMass = 0.939565378;

  // Map typical species run at RHIC
  map<TString,Int_t> A, Z;
  A["p"]  =   1;    Z["p"]  =  1;
  A["n"]  =   1;    Z["n"]  =  0;
  A["d"]  =   2;    Z["d"]  =  1;
  A["He3"]=   3;    Z["He3"]=  2;
  A["Au"] = 197;    Z["Au"] = 79;
  A["Cu"] =  64;    Z["Cu"] = 29;
  A["U"]  = 238;    Z["U"]  = 92;

  A["proton"]   =1;    Z["proton"]   =1;
  A["neutron"]  =1;    Z["neutron"]  =0;
  A["deuteron"] =2;    Z["deuteron"] =1;
  A["e-"]       =0;    Z["e-"]       =0;
  A["electron"] =0;    Z["electron"] =0;
  A["e+"]       =0;    Z["e+"]       =0;
  A["positron"] =0;    Z["positron"] =0;


  TString myBlue = mBlue;
  TString myYell = mYell;

  stringstream Blue;
  stringstream Yell;

  Blue << A[myBlue] << " " << Z[myBlue];
  Yell << A[myYell] << " " << Z[myYell];

  InputParametersString["pro"] = Blue.str().c_str();
  InputParametersString["tar"] = Yell.str().c_str();

  stringstream impactParameters;
  impactParameters << mImpactMin << " " << mImpactMax;
  InputParametersString["IMP"] = impactParameters.str().c_str();

  //
  // Create a new event record for either pp or eo events
  //
  if ( ( mBlue == "proton" ) && ( mYell == "proton" ) )                 mEvent = new StarGenPPEvent();
  else                                                                  mEvent = new StarGenAAEvent();   

  /// Remapt the ROOT names to UrQMD names
  std::map< TString, string > particle;
  /// Ex: particle["ROOT name"]="UrQMD name";
  particle["proton"] = "P       ";
  particle["e-"]     = "E-      ";

  /// Set up frames
  if ( mFrame=="COM" )
    {
      InputParametersDouble["ecm"]=mRootS; //???
    }
  if ( mFrame=="3MOM" || mFrame=="4MOM" || mFrame=="5MOM" )
    {
      ///JFN 11/21/12 4:11pm- I believe this calculation of mRootS from the momentum is correct, but I wouldn't stake my life on it. Additionally, I think mRootS should be something that should be calculated by the StarGenerator framework (ie, I shouldn't have to do it)
      mRootS = ( sqrt(pow(((Z[myBlue]*ProtonMass)+((A[myBlue]-Z[myBlue])*NeutronMass)),2) + sqrt( pow(mBlueMomentum.Px(),2) + pow(mBlueMomentum.Py(),2) + pow(mBlueMomentum.Pz(),2))) + sqrt(pow(((Z[myYell]*ProtonMass)+((A[myYell]-Z[myYell])*NeutronMass)),2) + sqrt( pow(mYellMomentum.Px(),2) + pow(mYellMomentum.Py(),2) + pow(mYellMomentum.Pz(),2))));
      InputParametersDouble["ecm"]=mRootS;
    }

  /// Set particles to not decay
  ///JFN 11/20/12 12:36pm- I have figured out how to set particles to not decay; they must be listed in the input file in the form "stb [itpy#]". I need to figure out the convention for ityp id numbers.
  ///JFN 11/20/12 12:46pm- See the UrQMD manual, page 12, tables 2 and 3 for particle IDs
  ///JFN 11/25/12 1:17pm- ityp doesn't seem to be completely specific for defingin particle type. Also, "antibaryons carry a negative sign", so some of these may need to be negative.
  // PI0 111
  // PI+ 211
  StableParticles.push_back("101");
  // ETA 221
  StableParticles.push_back("102");
  // K+ 321
  //$$$  StableParticles.push_back("");
  // K_SHORT 310
  //$$$  StableParticles.push_back("");
  // K_LONG 130
  //$$$  StableParticles.push_back("");
  // LAMBDA0 3122
  StableParticles.push_back("27");
  // SIGMA0 3212
  // SIGMA- 3112
  // SIGMA+ 3222
  StableParticles.push_back("40");
  // Xi- 3312
  // Xi0 3322
  StableParticles.push_back("49");
  // OMEGA- 3334
  StableParticles.push_back("55");

  /// Initialize UrQMD:
  InitializeUrQMD();



  /* call */  urqmd_init();


  return StMaker::Init();
}
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
Int_t StarUrQMD::Generate()
{
  /// Generate an event
  //  GenerateEvent();

  /* call */ urqmd_make();


#if 1

  // Blue beam is a proton, running PP
  if ( ( mBlue == "proton" ) && ( mYell == "proton" ) )    FillPP( mEvent );
  // Otherwise, runnin EP
  else                                                     FillAA( mEvent );

#endif 

#if 1

  //Do Stuff with the particles
  mNumberOfParticles = 1;
  //mNumberOfParticles = isys().ncoll(1); //JFN 11/28/12- this is wrong
  for ( Int_t idx=1; idx<=mNumberOfParticles; idx++ )
    {

      Int_t    id = isys().uid(idx); // or isys().itypd(idx). It isn't clear which is right
      Int_t    stat = StarGenParticle::kFinal; //JFN 11/25/12 12:28pm- for the moment I am setting every status to be a final state particle.
      //Int_t    stat = mStatusCode[ hepevt().isthep(idx) ]; */@
      //    if ( !stat ) {
      //	stat = StarGenParticle::kUnknown;
      //    }
      Int_t    m1 = itdelay().ityptd(idx,1);
      Int_t    m2 = itdelay().ityptd(idx,2);
      Int_t    d1 = 0; //JFN- I don't think daughter information is preserved
      Int_t    d2 = 0;
      Double_t px = coor().px(idx);
      Double_t py = coor().py(idx);
      Double_t pz = coor().pz(idx);
      Double_t E  = coor().p0(idx);
      Double_t M  = coor().fmass(idx);
      Double_t vx = px/M;
      Double_t vy = py/M;
      Double_t vz = pz/M;
      Double_t vt = sqrt(E*2/M); //E=m*v^2/2, v=sqrt(E*2/M)

      mEvent -> AddParticle( stat, id, m1, m2, d1, d2, px, py, pz, E, M, vx, vy, vz, vt );
    }


  mEvent->Print();

#endif

  return kStOK;
}
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
///JFN 11/18/12 13:24 - I think having a clear function is optional, and UrQMD doesn't have any explicit cleanup. Although, we could toss the un-needed output files here.
/*Int_t StarUrQMD::Clear()
{
  return kStOK;
}*/
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
void StarUrQMD::FillPP( StarGenEvent *event)
{

  // Fill the event-wise information
  StarGenPPEvent *myevent = (StarGenPPEvent *)event;
  myevent -> idBlue     = 0;
  myevent -> idYell     = 0;
  myevent -> process    = 0;
  /*myevent -> idBlue     = hwbeam().ipart1; // Int //JFN 11/26/12- "Blue beam ID". I don't know the ID convention
  myevent -> idYell     = hwbeam().ipart2;
  myevent -> process    = hwproc().iproc; //JFN 11/26/12- in principle there are process and subprocess ids because they get written in the headers of the output files, but I cant figure out where they are stored*/
  myevent -> subprocess = 0;

  myevent -> idParton1  = -999;
  myevent -> idParton2  = -999;
  myevent -> xParton1   = 0;
  myevent -> xParton2   = 0;
  myevent -> xPdf1      = -999;
  myevent -> xPdf2      = -999;
  myevent -> Q2fac      = -999;
  myevent -> Q2ren      = -999;
  myevent -> valence1   = 0;
  myevent -> valence2   = 0;

  myevent -> sHat       = 0;
  myevent -> tHat       = 0;
  myevent -> uHat       = 0;
  myevent -> ptHat      = -999;
  myevent -> thetaHat   = -999;
  myevent -> phiHat     = -999;
  
  myevent -> weight     = -999;

}


void StarUrQMD::FillAA( StarGenEvent *event)
{

  // Fill the event-wise information
  StarGenAAEvent *myevent = (StarGenAAEvent *)event;
  myevent -> idBlue     = 0;
  myevent -> idYell     = 0;
  myevent -> process    = 0;
  /*myevent -> idBlue     = hwbeam().ipart1; // Int //JFN 11/26/12- "Blue beam ID". I don't know the ID convention
  myevent -> idYell     = hwbeam().ipart2;
  myevent -> process    = hwproc().iproc; //JFN 11/26/12- in principle there are process and subprocess ids because they get written in the headers of the output files, but I cant figure out where they are stored*/
  myevent -> subprocess = 0;

  myevent -> idParton1  = -999;
  myevent -> idParton2  = -999;
  myevent -> xParton1   = 0;
  myevent -> xParton2   = 0;
  myevent -> xPdf1      = -999;
  myevent -> xPdf2      = -999;
  myevent -> Q2fac      = -999;
  myevent -> Q2ren      = -999;
  myevent -> valence1   = 0;
  myevent -> valence2   = 0;

  myevent -> sHat       = 0;
  myevent -> tHat       = 0;
  myevent -> uHat       = 0;
  myevent -> ptHat      = -999;
  myevent -> thetaHat   = -999;
  myevent -> phiHat     = -999;
  
  myevent -> weight     = -999;

}
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
///JFN 11/19/12 7:54am- to initialize UrQMD we first have to create input files for UrQMD to read and set enviormental variables, then we call UrQMD's init function
///JFN 11/19/12 7:37pm- In hindsight, I think I might rather have these functions defined in StarUrQMD.cxx.... maybe
///JFN 11/20/12 12:09pm- I have moved these two function (InitializeUrQMD and GenerateEvent) to StarUrQMD.cxx
void StarUrQMD::InitializeUrQMD()
  {
    /// Set the enviormental variables
    ///JFN 11/19/12 7:55am- Honestly, it pains me to set enviormental variables just to run this. I need to check if there is a way I can just set common block variables
    ///JFN 11/19/12 7:38pm- I would like to be able to just define the name on the input output files....
    ///JFN 11/21/12 5:19pm- I have found where the enviormental variables are loaded for UrQMD: line 359 of input.F. I see no easy way to spoof this behavior without modifying the source code, which I am trying to keep to an absolute minimum,
    //SetEnVars();

    /// Print the input file
    std::ofstream inputfile;
    inputfile.open( "UrQMD.in" );
    for(map<TString,Int_t>::iterator i = InputParametersInt.begin(); i != InputParametersInt.end(); i++)
      { 
        inputfile << i->first << " " << i->second << endl;
      }
    for(map<TString,Double_t>::iterator i = InputParametersDouble.begin(); i != InputParametersDouble.end(); i++)
      { 
        inputfile << i->first << " " << i->second << endl;
      }
    for(map<TString,TString>::iterator i = InputParametersString.begin(); i != InputParametersString.end(); i++)
      { 
        inputfile << i->first << " " << i->second << endl;
      }
    /// Set particles to not decay
    for(unsigned int i = 0; i< StableParticles.size(); i++)
      {
        inputfile << "stb " << StableParticles[i] << endl;
      }
    /// 
    if ( mFrame=="CMS" ) {
      inputfile << "ecm " << mRootS << "\n";
    }
    /// Define calculation time: total time span, interval at which output is written (both in fm/c)
    inputfile << "tim 200 200" << endl;
    /// This supresses the output files
    inputfile << "f13 \n #f14 \n f15 \n f16 \n #f17 \n #f18 \n f19 \n f20" << endl;
    /// This sets the number of events for UrQMD to run (but we aren't using this, so it doens't matter)
    inputfile << "nev 1000" << endl;
    /// This sets the random number generator seed (but we aren't using UrQMD's random number generator, so this doens't matter)
    inputfile << "rsd 111" << endl;

    
    if ( SAttr("SET" )  ) {
      for ( auto cmd : tokenize( SAttr("Set"), ";" ) ) {
	inputfile << cmd << "\n";
	LOG_INFO << "Set urqmd option:" << cmd << endm;
      }
    }


    /// This marks the end of the input file
    inputfile << "xxx and done" << endl;
    inputfile.close();

    /// Initialize UrQMD
    //    iurqmd();
  }
// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
///JFN 11/19/12 7:57am- I think we can grab the event information right from the common blocks, so there may be no need to process anything, or for this function
///JFN 11/21/12 5:28pm- just for reference, the first calls for setting up output files in RunFunctin.F is at line 169. Final output is done at line 408.
/*void GenerateEvent()
  {
    /// Call the UrQMD function to make a new event
    genevt();

    /// Process the results
    //ProcessEvent();
  }*/
#endif
