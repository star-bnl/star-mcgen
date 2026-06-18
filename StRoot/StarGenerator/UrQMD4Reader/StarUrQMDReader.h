#ifndef __StarUrQMDReader__
#define __StarUrQMDReader__

#include "StarGenerator/BASE/StarGenerator.h"

#include <map>
#include <vector>

//________________________________________________________________________________________
class StarUrQMDReader : public StarGenerator {
public:

  StarUrQMDReader( const Char_t *name="UrQMD4Reader" );
  ~StarUrQMDReader(){ /* nothing */ };

  Int_t Init();
  Int_t Generate();
  void Clear( const Option_t *opts="" ){ /* nada */ };

  void FillAA( StarGenEvent *event );

  ClassDef(StarUrQMDReader,0);

};

#endif
