// Example STAR macro for running the UrQMD 4 bridge through StarPrimaryMaker.
//
// Usage:
//   root4star starsim.urqmd4.C
//
// Note:
//   The checked-in UrQMD4 bridge currently contains a diagnostic assert in
//   StarUrQMD::Init(). This macro documents the intended STAR-side setup.

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarUrQMD;
StarUrQMD *urqmd = 0;

// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  for ( Int_t i=0; i<n; ++i ) {
    chain->Clear();
    chain->Make();
    _primary -> event() -> Print();
  }
}
// ----------------------------------------------------------------------------
Int_t loadUrQMD4()
{
  Int_t ierr = gSystem->Load( "libUrQMD4_0_0.so" );
  if ( ierr < 0 ) ierr = gSystem->Load( "UrQMD4_0_0.so" );
  return ierr;
}
// ----------------------------------------------------------------------------
void UrQMD4( Double_t ecm=200.0,
             const Char_t *blue="Au",
             const Char_t *yell="Au",
             Double_t bmin=0.0,
             Double_t bmax=30.0,
             const Char_t *urqmdSet="" )
{
  urqmd = new StarUrQMD("urqmd4");
  urqmd->SetTitle("UrQMD 4.0.0");

  urqmd->SetFrame("CMS", ecm);
  urqmd->SetBlue( blue );
  urqmd->SetYell( yell );
  urqmd->SetImpact( bmin, bmax );

  if ( urqmdSet && TString(urqmdSet).Length() ) {
    // The bridge checks both SET and Set inconsistently, so set both.
    urqmd->SetAttr("SET", urqmdSet);
    urqmd->SetAttr("Set", urqmdSet);
  }

  _primary -> AddGenerator( urqmd );
}
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=10,
              Int_t rngSeed=1234,
              const Char_t *tag="y2018",
              Double_t ecm=200.0,
              const Char_t *blue="Au",
              const Char_t *yell="Au",
              Double_t bmin=0.0,
              Double_t bmax=30.0,
              const Char_t *urqmdSet="" )
{
  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = Form("%s geant gstar usexgeom agml ", tag);
    bfc( 0, simple );
  }

  gSystem->Load( "libVMC.so" );
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  gSystem->Load( "libMathMore.so" );
  gSystem->Load( "xgeometry.so" );

  if ( loadUrQMD4() < 0 ) {
    std::cout << "Unable to load the UrQMD4 generator library." << std::endl;
    std::cout << "Tried libUrQMD4_0_0.so and UrQMD4_0_0.so." << std::endl;
    return;
  }

  StarRandom::seed( rngSeed );
  StarRandom::capture();

  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "urqmd4.starsim.root" );
    _primary -> SetVertex( 0.0, 0.0, 0.0 );
    _primary -> SetSigma ( 0.1, 0.1, 30.0 );
    chain -> AddBefore( "geant", _primary );
  }

  UrQMD4( ecm, blue, yell, bmin, bmax, urqmdSet );

  _primary->SetPtRange  ( 0.0, -1.0 );
  _primary->SetEtaRange ( -5.0, +5.0 );
  _primary->SetPhiRange ( 0.0, TMath::TwoPi() );

  _primary -> Init();

  command("gkine -4 0");
  command("gfile o urqmd4.starsim.fzd");

  trig( nevents );

  chain->Finish();

  if ( gROOT->IsBatch() ) {
    command("call agexit");
  }
  else {
    std::cout << "Interactive mode. Please call AgExit() before .q" << std::endl;
  }
}
// ----------------------------------------------------------------------------
void AgExit() { command("call agexit"); }
