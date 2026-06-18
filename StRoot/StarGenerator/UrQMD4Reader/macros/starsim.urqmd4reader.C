// Example STAR macro for running the UrQMD 4 bridge through StarPrimaryMaker.
//
// Usage:
//   root4star starsim.urqmd4reader.C
//
// Note:
// Requires the output of the urqmd event generator.  By default it looks for 

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarUrQMDReader;
StarUrQMDReader *urqmd = 0;

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
  Int_t ierr = gSystem->Load( "libUrQMD4Reader.so" );
  if ( ierr < 0 ) ierr = gSystem->Load( "UrQMD4Reader.so" );
  return ierr;
}
// ----------------------------------------------------------------------------
void UrQMD4( const char* filename="urqmd-4.0/fort.14" ) 
{
  urqmd = new StarUrQMDReader("urqmd4");
  //  urqmd->SetTitle("UrQMDReader 4.0.0");
  urqmd->SetAttr("filename",filename);
  _primary -> AddGenerator( urqmd );
}
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=1,
              Int_t rngSeed=1234,
              const Char_t *tag="y2018",
              const Char_t *filename="urqmd-4.0/fort.14" )
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

  UrQMD4("urqmd-4.0/fort.14");

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
