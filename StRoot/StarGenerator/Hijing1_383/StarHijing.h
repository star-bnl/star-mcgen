#ifndef __StarHijing_h__
#define __StarHijing_h__

#include "StarGenerator/BASE/StarGenerator.h"

#include "Hijing.h"
#include <map>
using namespace std;

class StarHijing : public StarGenerator
{
 public:
  StarHijing( const Char_t *name="Hijing" );
  ~StarHijing(){ /* nada */ };

  Int_t Init();
  Int_t Generate();

  /// Returns a reference to the hijing parameters
  HiParnt_t &hiparnt(){ return *address_of_hiparnt();}

  /// Returns a reference to the hijing main1 block
  HiMain1_t &himain1(){ return *address_of_himain1(); }
  /// Returns a refernece to the hijing main2 block
  HiMain2_t &himain2(){ return *address_of_himain2(); }

  /// Returns a reference to the ludat3 (pydat3) common block
  Ludat3_t &ludat3(){ return *address_of_ludat3(); }

 private:
 protected:
  ClassDef(StarHijing,1);

  void FillAA( StarGenEvent *event );

  map<Int_t, Int_t> mStatusCode;
  map<Int_t, Int_t> mParticleCode;

  // Count number of spectator protons and neutrons
  Int_t mNumberOfSpectatorProtons[2];
  Int_t mNumberOfSpectatorNeutrons[2];

  Int_t mNumberOfBeamProtons[2];
  Int_t mNumberOfBeamNeutrons[2];

  /// Given the event generator's native particle code, returns the corresponding
  /// PDG code
  /// @param code The event generator's native code
  Int_t pdgid(const Int_t &code);

};

#endif
