#ifndef __StarPythia8_h__
#define __StarPythia8_h__

#include "StarGenerator/BASE/StarGenerator.h"
#include "Pythia.h"

class StarGenPPEvent;
class StarGenEPEvent;

class StarPythia8 : public StarGenerator
{

 public:
  StarPythia8( const Char_t *name="Pythia8" );
  ~StarPythia8(){ /* nada */ };

  Int_t Init();
  Int_t Generate();

  void Set( const Char_t *s ){ mPythia -> readString(s); }

  virtual const char *GetCVS() const
  {static const char cvs[]="Tag $Name:  $ $Id: StarPythia8.h,v 1.1 2012/11/26 18:57:08 jwebb Exp $ built "__DATE__" "__TIME__ ; return cvs;}

 private:
 protected:

  Int_t InitCMS ( Int_t blue, Int_t yell );
  Int_t InitFIXT( Int_t blue, Int_t yell );
  Int_t Init3MOM( Int_t blue, Int_t yell );
  Int_t Init4MOM( Int_t blue, Int_t yell );
  Int_t Init5MOM( Int_t blue, Int_t yell );

  Pythia8::Pythia *mPythia;

  void FillPP( StarGenEvent *event );
  void FillEP( StarGenEvent *event );

  ClassDef(StarPythia8,1);

};

#endif
