#include "StarUrQMDReader.h"

#include "StarGenerator/EVENT/StarGenAAEvent.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include "StarGenerator/UTIL/StarRandom.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

extern "C" {
  int pdgid_( int& ityp, int& iso3 );
};


//___________________________________________________________________________________
struct UrQMDTrack_t {
  double r0;
  double rx;
  double ry;
  double rz;
  double p0;
  double px;
  double py;
  double pz;
  double mass;
  int ityp;
  int iso3;
  int charge;
  int lcl_;
  int ncl_;
  int or_; 
};

UrQMDTrack_t read_particle ( std::ifstream& urqmd, std::string& line ) {
  std::getline( urqmd, line ); 
  std::cout << "      [particle] " << line << std::endl;
  UrQMDTrack_t track;
  std::istringstream ss(line);
  ss>>  track.r0 ;
  ss>>  track.rx ;
  ss>>  track.ry ;
  ss>>  track.rz ;
  ss>>  track.p0 ;
  ss>>  track.px ;
  ss>>  track.py ;
  ss>>  track.pz ;
  ss>>  track.mass ;
  ss>>  track.ityp;
  ss>>  track.iso3;
  ss>>  track.charge;
  
  //  track.lcl_ ;
  //  track.ncl_ ;
  //  track.or_ ;
  return track;
};

std::vector<UrQMDTrack_t> read_particle_vector( std::ifstream& urqmd, std::string& line ) {
  std::vector<UrQMDTrack_t> result;
  std::getline( urqmd, line );  
  std::cout << "    [nparticles] " << line << std::endl;
  std::istringstream ss(line);
  int ntracks = 0;
  int time=0;
  ss>>ntracks;
  ss>>time;
  //  std::cout << " ntracks=" << ntracks << " time=" << time << std::endl;
  std::getline( urqmd, line );  //std::cout << "RPV: " << line << std::endl;
  for ( int i=0;i<ntracks;i++ ) {
    UrQMDTrack_t track = read_particle( urqmd, line );
    //std::cout << track.px << " " << track.py << " " << track.pz << " " << track.ityp << " " << track.iso3 << " " << track.charge << std::endl;
    result.push_back(track);
  }
  return result;
}

std::vector<UrQMDTrack_t> read_event( std::ifstream& urqmd, std::string& line ) {
  std::cout << "[UrQMD] " << line <<   std::endl;
  std::vector<UrQMDTrack_t> result;
  while ( std::getline( urqmd, line ) ) {
    std::cout << "  [header] " << line << std::endl;
    // For now exit on the particle vector
    if ( line.find( "pvec:" ) != std::string::npos ) {
      result = read_particle_vector( urqmd, line );
      break;
    };
  };
  return result;
}
//___________________________________________________________________________________
//
//
std::ifstream urqmd;
//
//
//___________________________________________________________________________________
StarUrQMDReader::StarUrQMDReader( const char* name ) : StarGenerator(name) 
{
  // Default input file
  SetAttr("input","fort.14");
}
//___________________________________________________________________________________
int StarUrQMDReader::Init() {

  // Open the urqmd output file
  urqmd.open( SAttr("filename"),  std::ifstream::in);

  // Set the event
  mEvent = new StarGenAAEvent();   

  return StMaker::Init();
}
double rndm(){ return StarRandom::Instance().flat(); };
//___________________________________________________________________________________
int StarUrQMDReader::Generate() {

  StarParticleData& pdb = StarParticleData::instance();

  if ( urqmd.eof() ) return kStEOF;

  std::string line;
  std::vector<UrQMDTrack_t> tracks = read_event( urqmd, line );

  for ( auto& track : tracks ) {
    double x=track.rx;
    double y=track.ry;
    double z=track.rz;
    double px=track.px;
    double py=track.py;
    double pz=track.pz;
    double mass=track.mass;
    int    ityp=track.ityp;
    int    iso3=track.iso3;
    int    charge=track.charge;

    int pdgid = pdgid_( ityp, iso3 );    
    // K0 to K0 short/long
    if ( pdgid==311 || pdgid==-311 ) {
      if ( rndm() > 0.5 ) pdgid = 130  ; // K0 long
      else                pdgid = 310  ; // K0 short
    }

    // update mass based on particle database
    auto* pdg = pdb.GetParticle( pdgid );
    mass = pdg->Mass();

    double E2 = px*px + py*py + pz*pz + mass*mass;
    double E  = sqrt(E2);

    //    std::cout << "[pdgid] " << pdgid << " p=" << px << " " << py << " " << pz << std::endl;

    // o Add the particle as a final state particle
    // o This is a primary particle.  Mother and daughter are set to zero.  Geant responsble for decays.
    // o 
    mEvent -> AddParticle( StarGenParticle::kFinal, pdgid, 0, 0, 0, 0, px, py, pz, E, mass, 0., 0., 0., 0. );

  }

  mEvent->Print();


  return kStOK;
}
//___________________________________________________________________________________
void StarUrQMDReader::FillAA( StarGenEvent *event)
{

#if 0
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
#endif

}
