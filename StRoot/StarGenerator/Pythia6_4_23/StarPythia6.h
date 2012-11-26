#ifndef __StarPythia5_h__
#define __StarPythia6_h__

#include "StarGenerator/BASE/StarGenerator.h"
#include "Pythia6.h"
#include <map>
using namespace std;

class StarPythia6 : public StarGenerator
{
 public:
  StarPythia6( const Char_t *name="Pythia6" );
  ~StarPythia6(){ /* nada */ };

  Int_t Init();
  Int_t Generate();

  /// Returns a reference to the /PYJETS/ common block
  PyJets_t &pyjets(){ return *address_of_pyjets(); }
  /// Returns a reference to the /PYSUBS/ common block
  PySubs_t &pysubs(){ return *address_of_pysubs(); }
  /// Returns a reference to the /PYDAT3/ common block
  PyDat3_t &pydat3(){ return *address_of_pydat3(); }
  /// Returns a reference to the /PYPARS/ common block
  PyPars_t &pypars(){ return *address_of_pypars(); }

  /// Calls the pytune function
  void PyTune( Int_t tune );

 private:
 protected:
  ClassDef(StarPythia6,1);

  void FillPP( StarGenEvent *event );
  void FillEP( StarGenEvent *event );

  map<Int_t,Int_t> mStatusCode;

};

#endif
