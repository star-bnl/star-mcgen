#ifndef __Hijing_h__
#define __Hijing_h__

#include "StarCallf77.h"
#include <string>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//
// Wrappers for FORTRAN routines
//
//////////////////////////////////////////////////////////////////////////////////////////

/// Configuration of HIJING
void Hijset( float   efrm, string frame, string proj, string targ, int ap, int zp, int at, int zt );
/// Generation of one event in hijing
void Hijing( string frame, float bmin, float bmax );

/// Lookup particle mass given jetset code
float Ulmass( int jetsetid );


//////////////////////////////////////////////////////////////////////////////////////////
//
// Interfaces to COMMON blocks
//
//////////////////////////////////////////////////////////////////////////////////////////

#define address_of_himain1 F77_NAME( address_of_himain1, ADDRESS_OF_HIMAIN1 )
struct HiMain1_t {
  /* memory layout */
  int   natt;
  float eatt;
  int   jatt;
  int nwounded_yell;
  int nwounded_blue;
  int n0, n01, n10, n11;  
  /* access fuctions */
};
extern "C" HiMain1_t *address_of_himain1();

#define address_of_himain2 F77_NAME( address_of_himain2, ADDRESS_OF_HIMAIN2 )
struct HiMain2_t {
  /* memory layout */
  int   _katt[4][130000];
  float _patt[4][130000];
  float _vatt[4][130000];
  /* access functions */ 
  int   &katt(int i, int j){ return _katt[j-1][i-1]; }
  float &patt(int i, int j){ return _patt[j-1][i-1]; }
  float &vatt(int i, int j){ return _vatt[j-1][i-1]; }
};
extern "C" HiMain2_t *address_of_himain2();

#define address_of_hiparnt F77_NAME( address_of_hiparnt, ADDRESS_OF_HIPARNT )
struct HiParnt_t
{
  /* memory layout */
  float _hipr1[100]; 
  int   _ihpr2[50];  
  float _hint1[100]; 
  int   _ihnt2[50];  
  /* access functions */
  float &hipr1(int i){ return _hipr1[i-1]; }
  int   &ihpr2(int i){ return _ihpr2[i-1]; }
  float &hint1(int i){ return _hint1[i-1]; }
  int   &ihnt2(int i){ return _ihnt2[i-1]; }

};
extern "C" HiParnt_t *address_of_hiparnt();

//  COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)  
#define address_of_ludat3 F77_NAME( address_of_ludat3, ADDRESS_OF_LUDAT3 )
struct Ludat3_t {
  /* memory layout */
  int _mdcy[3][500];
  int _mdme[2][2000];
  float _brat[2000];
  int _kfdp[5][2000];
  /* access functions */
  int &mdcy(int i, int j){ return _mdcy[j-1][i-1]; }
  int &mdme(int i, int j){ return _mdme[j-1][i-1]; }
  float &brat( int i ){ return _brat[i-1]; }
  int &kfdp( int i, int j ){ return _kfdp[j-1][i-1]; }
};

extern "C" Ludat3_t *address_of_ludat3();



#endif
