#ifndef __StarPythia8_h__
#define __StarPythia8_h__

/*!
  \class StarPythaia8
  \author Jason C. Webb
  \brief Interface to Pythia 8

  StarPythia8 provides the user interface to the Pythia 8 event generator.
  To configure pythia8 for a specific run, users are encouraged to consult
  the Pythia 8 manual: http://home.thep.lu.se/~torbjorn/pythia81html/Welcome.html.
  In particular the section on the settings scheme http://home.thep.lu.se/~torbjorn/pythia81html/SettingsScheme.html.

  The StarPythia8::Set() method passes configuration strings to the instance
  of pythia, as described under the "Operation" heading.

  Configuration of the beam parameters (energy, frame, species, etc...)
  is through methods defined on the StarGenerator base class.


 */

#include "StarGenerator/BASE/StarGenerator.h"
#include "Pythia.h"

class StarGenPPEvent;
class StarGenEPEvent;

class StarPythia8 : public StarGenerator
{

 public:
  StarPythia8( const Char_t *name="Pythia8" );
  ~StarPythia8(){ /* nada */ };

  /// Initialize the event generator
  Int_t Init();
  /// Generate one event
  Int_t Generate();

  /// Pass a string to Pythia8::Pythia::readString(), for user configuration.
  void Set( const Char_t *s ){ mPythia -> readString(s); }

  virtual const char *GetCVS() const
  {static const char cvs[]="Tag $Name:  $ $Id: StarPythia8.h,v 1.2 2012/12/05 21:50:57 jwebb Exp $ built "__DATE__" "__TIME__ ; return cvs;}

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
