#ifndef __UrQMD4_0_0__
#define __UrQMD4_0_0__

#include "StarCallf77.h"
#include <string>
using namespace std;

#include "TObject.h"

//
// Interface to common Blocks
//

#define energies F77_NAME(energies,ENERGIES)
struct ENERGIES_t {
      /*real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau*/
  Double_t  Ekinbar;
  Double_t  Ekinmes;
  Double_t  ESky2;
  Double_t  ESky3;
  Double_t  EYuk;
  Double_t  ECb;
  Double_t  EPau;
};
extern "C" ENERGIES_t *address_of_energies();

#define sys F77_NAME(sys,SYS)
struct SYS_t {
    /*
      integer npart, nbar, nmes, ctag, nsteps, uid_cnt, ranseed,
     &        event, Ap, At, Zp, Zt, eos, dectag, NHardRes, NSoftRes,
     &        NDecRes, NElColl, NBlColl
      logical success
      integer npartcoal, nclus
      common /sys/ npart, nbar, nmes, ctag, nsteps, uid_cnt,
     &             ranseed, event, Ap, At, Zp, Zt, eos, dectag,
     &             NHardRes, NSoftRes, NDecRes, NElColl, NBlColl,
     &             success, npartcoal, nclus
    */

  Int_t  npart;
  Int_t  nbar;
  Int_t  nmes;
  Int_t  ctag;
  Int_t  nsteps;
  Int_t  uid_cnt;
  Int_t  ranseed;
  Int_t  event;
  Int_t  Ap;
  Int_t  At;
  Int_t  Zp;
  Int_t  Zt;
  Int_t  eos;
  Int_t  dectag;
  Int_t  NHardRes;
  Int_t  NSoftRes;
  Int_t  NDecRes;
  Int_t  NElColl;
  Int_t  NBlColl; 
  Int_t  success;  // believe that logical maps onto integer
  Int_t npartcoal;
  Int_t nclus;
};
extern "C" SYS_t *address_of_sys();

#define rsys F77_NAME(rsys,RSYS)
struct RSYS_t {
    /*real*8  time,  acttime, bdist, ebeam, bimp, bmin, ecm
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm*/
  Double_t  time;
  Double_t  acttime;
  Double_t  bdist;
  Double_t  bimp;
  Double_t  bmin;
  Double_t  ebeam;
  Double_t  ecm;
};
extern "C" RSYS_t *address_of_rsys();

#define cuts F77_NAME(cuts,CUTS)
struct CUTS_t {
    /*real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww*/
  Double_t  cutmax;
  Double_t  cutPau;
  Double_t  cutCb;
  Double_t  cutYuk;
  Double_t  cutSky;
  Double_t  cutdww;
};
extern "C" CUTS_t *address_of_cuts();

#define spdata F77_NAME(spdata,SPDATA)
struct SPDATA_t {
   /* parameter (nspl = 500)  ! dimension of spline arrays
      real*8 spx(nspl), spPauy(nspl), outPau(nspl),
     &                spCby(nspl),  outCb(nspl),
     &                spYuky(nspl), outYuk(nspl),
     &                spSkyy(nspl), outSky(nspl),
     &                spdwwy(nspl), outdww(nspl)
      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     &                     spYuky, outYuk, spSkyy, outSky,
     &                     spdwwy, outdww*/
  Double_t  _spx[500];
  Double_t  _spPauy[500];
  Double_t  _outPau[500];
  Double_t  _spCby[500];
  Double_t  _outCb[500];
  Double_t  _spYuky[500];
  Double_t  _outYuk[500];
  Double_t  _spSkyy[500];
  Double_t  _outSky[500];
  Double_t  _spdwwy[500];
  Double_t  _outdww[500];
  Double_t  &spx( Int_t i ){return _spx[i-1]; }
  Double_t  &spPauy( Int_t i ){return _spPauy[i-1]; }
  Double_t  &outPau( Int_t i ){return _outPau[i-1]; }
  Double_t  &spCby( Int_t i ){return _spCby[i-1]; }
  Double_t  &outCb( Int_t i ){return _outCb[i-1]; }
  Double_t  &spYuky( Int_t i ){return _spYuky[i-1]; }
  Double_t  &outYuk( Int_t i ){return _outYuk[i-1]; }
  Double_t  &spSkyy( Int_t i ){return _spSkyy[i-1]; }
  Double_t  &outSky( Int_t i ){return _outSky[i-1]; }
  Double_t  &spdwwy( Int_t i ){return _spdwwy[i-1]; }
  Double_t  &outdww( Int_t i ){return _outdww[i-1]; }
};
extern "C" SPDATA_t *address_of_spdata();

#define isys F77_NAME(isys,ISYS)
struct ISYS_t {
   /* parameter (nmax = 100000) ! maximum number of particles
      integer spin(nmax), ncoll(nmax), charge(nmax),
     &        ityp(nmax), lstcoll(nmax), iso3(nmax), origin(nmax), uid(nmax)
      common/isys/ spin, ncoll, charge, ityp, lstcoll, iso3, origin, uid */
  Int_t  _spin[100000];
  Int_t  _ncoll[100000];
  Int_t  _charge[100000];
  Int_t  _ityp[100000];
  Int_t  _lstcoll[100000];
  Int_t  _iso3[100000];
  Int_t  _origin[100000];
  //Int_t  _strid[100000];
  Int_t  _uid[100000];
  Int_t  &spin( Int_t i ){return _spin[i-1]; }
  Int_t  &ncoll( Int_t i ){return _ncoll[i-1]; }
  Int_t  &charge( Int_t i ){return _charge[i-1]; }
  Int_t  &ityp( Int_t i ){return _ityp[i-1]; }
  Int_t  &lstcoll( Int_t i ){return _lstcoll[i-1]; }
  Int_t  &iso3( Int_t i ){return _iso3[i-1]; }
  Int_t  &origin( Int_t i ){return _origin[i-1]; }
  //  Int_t  &strid( Int_t i ){return _strid[i-1]; }
  Int_t  &uid( Int_t i ){return _uid[i-1]; }
};
extern "C" ISYS_t *address_of_isys();

#define coor F77_NAME(coor,COOR)
struct COOR_t {
   /* parameter (nmax = 100000) ! maximum number of particles
      real*8
     &     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     &     p0(nmax), px(nmax), py(nmax), pz(nmax),
     &     fmass(nmax), rww(nmax), dectime(nmax)
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime*/
  Double_t  _r0[100000];
  Double_t  _rx[100000];
  Double_t  _ry[100000];
  Double_t  _rz[100000];
  Double_t  _p0[100000];
  Double_t  _px[100000];
  Double_t  _py[100000];
  Double_t  _pz[100000];
  Double_t  _fmass[100000];
  Double_t  _rww[100000];
  Double_t  _dectime[100000];
  Double_t  &r0( Int_t i ){return _r0[i-1]; }
  Double_t  &rx( Int_t i ){return _rx[i-1]; }
  Double_t  &ry( Int_t i ){return _ry[i-1]; }
  Double_t  &rz( Int_t i ){return _rz[i-1]; }
  Double_t  &p0( Int_t i ){return _p0[i-1]; }
  Double_t  &px( Int_t i ){return _px[i-1]; }
  Double_t  &py( Int_t i ){return _py[i-1]; }
  Double_t  &pz( Int_t i ){return _pz[i-1]; }
  Double_t  &fmass( Int_t i ){return _fmass[i-1]; }
  Double_t  &rww( Int_t i ){return _rww[i-1]; }
  Double_t  &dectime( Int_t i ){return _dectime[i-1]; }
};
extern "C" COOR_t *address_of_coor();

#define frag F77_NAME(frag,FRAG)
struct FRAG_t {
   /* parameter (nmax = 100000) ! maximum number of particles
      real*8 tform(nmax), xtotfac(nmax)
      common /frag/ tform, xtotfac*/
  Double_t  _tform[100000];
  Double_t  _xtotfac[100000];
  Double_t  &tform( Int_t i ){return _tform[i-1]; }
  Double_t  &xtotfac( Int_t i ){return _xtotfac[i-1]; }
};
extern "C" FRAG_t *address_of_frag();

#define aios F77_NAME(aios,AIOS)
struct AIOS_t {
   /* parameter (nmax = 100000) ! maximum number of particles
      real*8 airx(nmax), airy(nmax), airz(nmax),
     &     aipx(nmax), aipy(nmax), aipz(nmax),
     &     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     &     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4)
      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     &              aorx, aory, aorz, aopx, aopy, aopz*/
  Double_t  _airx[100000];
  Double_t  _airy[100000];
  Double_t  _airz[100000];
  Double_t  _aipx[100000];
  Double_t  _aipy[100000];
  Double_t  _aipz[100000];
  Double_t  _aorx[4][100000];
  Double_t  _aory[4][100000];
  Double_t  _aorz[4][100000];
  Double_t  _aopx[4][100000];
  Double_t  _aopy[4][100000];
  Double_t  _aopz[4][100000];
  Double_t  &airx( Int_t i ){return _airx[i-1]; }
  Double_t  &airy( Int_t i ){return _airy[i-1]; }
  Double_t  &airz( Int_t i ){return _airz[i-1]; }
  Double_t  &aipx( Int_t i ){return _aipx[i-1]; }
  Double_t  &aipy( Int_t i ){return _aipy[i-1]; }
  Double_t  &aipz( Int_t i ){return _aipz[i-1]; }
  Double_t  &aorx( Int_t i, Int_t j ){return _aorx[j-1][i-1]; }
  Double_t  &aory( Int_t i, Int_t j ){return _aory[j-1][i-1]; }
  Double_t  &aorz( Int_t i, Int_t j ){return _aorz[j-1][i-1]; }
  Double_t  &aopx( Int_t i, Int_t j ){return _aopx[j-1][i-1]; }
  Double_t  &aopy( Int_t i, Int_t j ){return _aopy[j-1][i-1]; }
  Double_t  &aopz( Int_t i, Int_t j ){return _aopz[j-1][i-1]; }
};
extern "C" AIOS_t *address_of_aios();

#define pots F77_NAME(pots,POTS)
struct POTS_t {
    /*real*8
     &     gw, sgw, delr, fdel, dt,da, db, Cb0, Yuk0, Pau0, Sky20,
     &     Sky30, gamSky, gamYuk, drPau, dpPau, dtimestep
      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky,
     &              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     &              dt,da, db,dtimestep*/
  Double_t  Cb0;
  Double_t  Yuk0;
  Double_t  Pau0;
  Double_t  Sky20;
  Double_t  Sky30;
  Double_t  gamSky;
  Double_t  gamYuk;
  Double_t  drPau;
  Double_t  dpPau;
  Double_t  gw;
  Double_t  sgw;
  Double_t  delr;
  Double_t  fdel;
  Double_t  dt;
  Double_t  da;
  Double_t  db;
  Double_t  dtimestep;
};
extern "C" POTS_t *address_of_pots();

#define scoor F77_NAME(scoor,SCOOR)
struct SCOOR_t {
   /* parameter(smax=500)  ! maximum number of spectators
      real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     &         p0s(smax), pxs(smax), pys(smax), pzs(smax),
     &         sfmass(smax)
      common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass*/
  Double_t  _r0s[500];
  Double_t  _rxs[500];
  Double_t  _rys[500];
  Double_t  _rzs[500];
  Double_t  _p0s[500];
  Double_t  _pxs[500];
  Double_t  _pys[500];
  Double_t  _pzs[500];
  Double_t  _sfmass[500];
  Double_t  &r0s( Int_t i ){return _r0s[i-1]; }
  Double_t  &rxs( Int_t i ){return _rxs[i-1]; }
  Double_t  &rys( Int_t i ){return _rys[i-1]; }
  Double_t  &rzs( Int_t i ){return _rzs[i-1]; }
  Double_t  &p0s( Int_t i ){return _p0s[i-1]; }
  Double_t  &pxs( Int_t i ){return _pxs[i-1]; }
  Double_t  &pys( Int_t i ){return _pys[i-1]; }
  Double_t  &pzs( Int_t i ){return _pzs[i-1]; }
  Double_t  &sfmass( Int_t i ){return _sfmass[i-1]; }
};
extern "C" SCOOR_t *address_of_scoor();

#define sisys F77_NAME(sisys,SISYS)
struct SISYS_t {
   /* parameter(smax=500)  ! maximum number of spectators
      integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     &          suid(smax)
      common /sisys/ sspin, scharge, sityp, siso3, suid*/
  Int_t  _sspin[500];
  Int_t  _scharge[500];
  Int_t  _sityp[500];
  Int_t  _siso3[500];
  Int_t  _suid[500];
  Int_t  &sspin( Int_t i ){return _sspin[i-1]; }
  Int_t  &scharge( Int_t i ){return _scharge[i-1]; }
  Int_t  &sityp( Int_t i ){return _sityp[i-1]; }
  Int_t  &siso3( Int_t i ){return _siso3[i-1]; }
  Int_t  &suid( Int_t i ){return _suid[i-1]; }
};
extern "C" SISYS_t *address_of_sisys();

#define ssys F77_NAME(ssys,SSYS)
struct SSYS_t {
    /*integer nspec
      common /ssys/ nspec*/
  Int_t nspec;
};
extern "C" SSYS_t *address_of_ssys();

#define rtdelay F77_NAME(rtdelay,RTDELAY)
struct RTDELAY_t {
   /* parameter (nmax = 100000) ! maximum number of particles
      real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     &         fmasstd(2,nmax)
      common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd*/
  Double_t  _p0td[100000][2];
  Double_t  _pxtd[100000][2];
  Double_t  _pytd[100000][2];
  Double_t  _pztd[100000][2];
  Double_t  _fmasstd[100000][2];
  Double_t  &p0td( Int_t i, Int_t j ){return _p0td[j-1][i-1]; }
  Double_t  &pxtd( Int_t i, Int_t j ){return _pxtd[j-1][i-1]; }
  Double_t  &pytd( Int_t i, Int_t j ){return _pytd[j-1][i-1]; }
  Double_t  &pztd( Int_t i, Int_t j ){return _pztd[j-1][i-1]; }
  Double_t  &fmasstd( Int_t i, Int_t j ){return _fmasstd[j-1][i-1]; }
};
extern "C" RTDELAY_t *address_of_rtdelay();

#define itdelay F77_NAME(itdelay,ITDELAY)
struct ITDELAY_t {
   /* parameter (nmax = 100000) ! maximum number of particles
      integer ityptd(2,nmax),iso3td(2,nmax)
      common /itdelay/ityptd,iso3td*/
  Int_t  _ityptd[100000][2];
  Int_t  _iso3td[100000][2];
  Int_t  &ityptd( Int_t i, Int_t j ){return _ityptd[j-1][i-1]; }
  Int_t  &iso3td( Int_t i, Int_t j ){return _iso3td[j-1][i-1]; }
};
extern "C" ITDELAY_t *address_of_itdelay();

#define svinfo F77_NAME(svinfo,SVINFO)
struct SVINFO_t {
    /*integer itypt(2),uidt(2),origint(2),iso3t(2)
      common /svinfo/itypt,uidt,origint,iso3t*/
  Int_t  _itypt[2];
  Int_t  _uidt[2];
  Int_t  _origint[2];
  Int_t  _iso3t[2];
  Int_t  &itypt( Int_t i ){return _itypt[i-1]; }
  Int_t  &uidt( Int_t i ){return _uidt[i-1]; }
  Int_t  &origint( Int_t i ){return _origint[i-1]; }
  Int_t  &iso3t( Int_t i ){return _iso3t[i-1]; }
};
extern "C" SVINFO_t *address_of_svinfo();

#define ffermi F77_NAME(ffermi,FFERMI)
struct FFERMI_t {
  /* parameter (nmax = 100000) ! maximum number of particles
      real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
      common /ffermi/ ffermpx, ffermpy, ffermpz*/
  Double_t  _ffermpx[100000];
  Double_t  _ffermpy[100000];
  Double_t  _ffermpz[100000];
  Double_t  &ffermpx( Int_t i ){return _ffermpx[i-1]; }
  Double_t  &ffermpy( Int_t i ){return _ffermpy[i-1]; }
  Double_t  &ffermpz( Int_t i ){return _ffermpz[i-1]; }
};
extern "C" FFERMI_t *address_of_ffermi();

#define peq F77_NAME(peq,PEQ)
struct PEQ_t {
    /*real*8 peq1, peq2
      common /peq/ peq1,peq2*/
  Double_t  peq1;
  Double_t  peq2;
};
extern "C" PEQ_t *address_of_peq();

void iurqmd( void );
void genevt( void );

/// Global Parameters which are needed in UrQMD.* and StarUrQMD.*
///JFN 11/20/12 12:05pm- I am moving these to StarUrQMD.h
/*std::vector< vector <string> > InputParameters;
map<TString,Int_t> InputParametersInt;
map<TString,Double_t> InputParametersDouble;
map<TString,TString> InputParametersString;*/

/*
 * ============================================================================
 * MISSING COMMON BLOCKS — UrQMD 4.0.0
 * ============================================================================
 * The following Fortran COMMON blocks exist in the UrQMD 4.0.0 source but have
 * no corresponding C++ struct in this header.  C++ struct definitions and
 * address_of_*() accessor declarations should be added here as they are needed
 * by the STAR bridge code.
 *
 * NOTE on array index transposition:
 *   Fortran arrays are column-major.  A Fortran array declared A(m,n) must be
 *   declared in C++ as _A[n][m], with accessor A(i,j) = _A[j-1][i-1].
 *
 * NOTE on Fortran LOGICAL:
 *   Fortran LOGICAL maps to Int_t (4 bytes) in C++.
 *
 * NOTE on blank (unnamed) COMMON:
 *   The Fortran blank common "common cross" (without /name/ delimiters) has no
 *   portable C-linkage symbol.  Accessing it from C++ requires a thin Fortran
 *   wrapper function returning the address of the array.
 *
 * NOTE on C++ reserved keyword 'try':
 *   The Fortran variable 'try' in /outcom/ cannot be used as a C++ method name.
 *   Use the name try_() for the accessor.
 *
 * ----------------------------------------------------------------------------
 * From coms.inc / coms90.inc  (new common blocks added in UrQMD 4)
 * ----------------------------------------------------------------------------
 *
 *   /eccA/ — eccentricity observables
 *     real*8 eccentricity, num_part, proton_part, neutron_part,
 *    &       eps1, eps2, eps3, eps4
 *     common /eccA/ eccentricity, num_part, proton_part, neutron_part,
 *    &              eps1, eps2, eps3, eps4
 *
 *   /comseed/ — random-seed first-use flag
 *     logical firstseed
 *     common /comseed/ firstseed
 *
 *   /logic/ — per-particle potential-type flags
 *     logical lsct(nmax), logSky, logYuk, logCb, logPau
 *     common /logic/ lsct, logSky, logYuk, logCb, logPau
 *
 *   /mdprop/ — molecular-dynamics propagator temporary coordinate arrays
 *     real*8 r0_t(nmax), rx_t(nmax), ry_t(nmax), rz_t(nmax)
 *     common /mdprop/ r0_t, rx_t, ry_t, rz_t
 *
 *   /cluster/ — nuclear coalescence cluster data  (smax=500)
 *     real*8    r0cl(500), rxcl(500), rycl(500), rzcl(500)
 *     integer   lstcollcl(500), ncollcl(500)
 *     real*8    p0cl(500), pxcl(500), pycl(500), pzcl(500), mcl(500)
 *     integer   chgcl(500), isocl(500), itypcl(500), idcl(500)
 *     common /cluster/ r0cl, rxcl, rycl, rzcl, lstcollcl,
 *    &                 ncollcl, p0cl, pxcl, pycl, pzcl, mcl,
 *    &                 chgcl, isocl, itypcl, idcl
 *
 *   blank common — ccbar-type crossing counter
 *     integer cross(60)
 *     common cross
 *
 * ----------------------------------------------------------------------------
 * From inputs.inc / inputs90.inc
 * ----------------------------------------------------------------------------
 *
 *   /inputs/ — primary event-configuration integers
 *     integer nevents, spityp(2), prspflg, trspflg, spiso3(2),
 *    &        outsteps, bflag, srtflag, efuncflag, nsrt, firstev, npb
 *     common /inputs/ nevents, spityp, prspflg, trspflg,
 *    &                spiso3, outsteps, bflag, srtflag, efuncflag,
 *    &                nsrt, firstev, npb
 *
 *   /input2/ — beam kinematics
 *     real*8 srtmin, srtmax, pbeam, betann, betatar, betapro, pbmin, pbmax
 *     common /input2/ srtmin, srtmax, pbeam, betann, betatar, betapro,
 *    &                pbmin, pbmax
 *
 *   /ProTarInts/ — projectile/target integer particle data
 *     (AAmax=300, MaxNTest=20 → dim = AAmax*MaxNTest = 6000)
 *     integer PT_iso3(6000,2), PT_ityp(6000,2), PT_spin(6000,2),
 *    &        PT_charge(6000,2), PT_AA(2), PT_uid(6000,2)
 *     common /ProTarInts/ PT_iso3, PT_ityp, PT_spin, PT_charge,
 *    &                    PT_AA, PT_uid
 *     C++ layout: PT_iso3(6000,2) → _PT_iso3[2][6000], accessor [j-1][i-1]
 *
 *   /ProTarReals/ — projectile/target real particle data
 *     real*8 PT_r0(6000,2), PT_rx(6000,2), PT_ry(6000,2), PT_rz(6000,2),
 *    &       PT_fmass(6000,2), PT_dectime(6000,2),
 *    &       PT_p0(6000,2), PT_px(6000,2), PT_py(6000,2), PT_pz(6000,2),
 *    &       PT_rho(6000,2), PT_pmax(6000,2)
 *     common /ProTarReals/ PT_r0, PT_rx, PT_ry, PT_rz, PT_fmass, PT_dectime,
 *    &                     PT_p0, PT_px, PT_py, PT_pz, PT_rho, PT_pmax
 *     C++ layout: all (6000,2) → [2][6000], accessor [j-1][i-1]
 *
 * ----------------------------------------------------------------------------
 * From options.inc / options90.inc
 * ----------------------------------------------------------------------------
 *
 *   /options/ — simulation option arrays
 *     integer CTOption(400)
 *     real*8  CTParam(400)
 *     common /options/ CTOption, CTParam
 *
 *   /optstrings/ — option documentation character arrays
 *     character ctodc(400)*2, ctpdc(400)*2
 *     common /optstrings/ ctodc, ctpdc
 *     C++ type: char _ctodc[400][2], char _ctpdc[400][2]  (not null-terminated)
 *
 *   /loptions/ — logical option flags
 *     logical fixedseed, bf13, bf14, bf15, bf16, bf19, bf20
 *     common /loptions/ fixedseed, bf13, bf14, bf15, bf16, bf19, bf20
 *
 *   /stables/ — stable-particle configuration
 *     integer nstable, stabvec(20)
 *     common /stables/ nstable, stabvec
 *
 * ----------------------------------------------------------------------------
 * From newpart.inc / newpart90.inc
 * ----------------------------------------------------------------------------
 *
 *   /inewpart/ — new-particle integer data per collision  (mprt=1000, oprt=2)
 *     integer itypnew(mprt), i3new(mprt), itot(mprt), inew(mprt),
 *    &        nexit, iline, pslot(oprt), nstring1, nstring2,
 *    &        itypold(oprt), iso3old(oprt)
 *     common /inewpart/ itypnew, i3new, itot, inew, nexit, iline,
 *    &                  pslot, nstring1, nstring2, itypold, iso3old
 *
 *   /rnewpart/ — new-particle real data per collision
 *     real*8 pnew(5,mprt), xnew(4,mprt), betax, betay, betaz,
 *    &       pold(5,oprt), p0nn, pxnn, pynn, pznn, pnn,
 *    &       mstring(oprt), pnnout, xtotfacold(oprt)
 *     common /rnewpart/ pnew, xnew, betax, betay, betaz, pold,
 *    &                  p0nn, pxnn, pynn, pznn, pnn, mstring,
 *    &                  pnnout, xtotfacold
 *     C++ layout: pnew(5,1000) → _pnew[1000][5], accessor [j-1][i-1]
 *                 xnew(4,1000) → _xnew[1000][4], accessor [j-1][i-1]
 *                 pold(5,2)    → _pold[2][5],     accessor [j-1][i-1]
 *
 *   /fnewpart/ — leading-hadron formation factor  (mprt=1000)
 *     real*8 leadfac(mprt)
 *     common /fnewpart/ leadfac
 *
 * ----------------------------------------------------------------------------
 * From comhis.inc  (entirely new file in UrQMD 4)
 * ----------------------------------------------------------------------------
 *
 *   /historyhead/ — collision-history header arrays  (wwmax=20000)
 *     integer   hisnin(wwmax), hisnexit(wwmax), hisiline(wwmax),
 *    &           hisctag(wwmax)
 *     real*8    hisacttime(wwmax), hissqrts(wwmax), hisstot(wwmax),
 *    &           hissigpart(wwmax), hiscdens(wwmax)
 *     common /historyhead/ hisnin, hisnexit, hisiline, hisctag,
 *    &                     hisacttime, hissqrts, hisstot, hissigpart,
 *    &                     hiscdens
 *
 *   /historyin/ — ingoing-particle collision history  (wwmax=20000, inmax=10)
 *     integer   INind(wwmax,inmax)
 *     real*8    INr0(wwmax,inmax), INrx(wwmax,inmax), INry(wwmax,inmax),
 *    &           INrz(wwmax,inmax), INp0(wwmax,inmax), INpx(wwmax,inmax),
 *    &           INpy(wwmax,inmax), INpz(wwmax,inmax), INmass(wwmax,inmax)
 *     integer   INityp(wwmax,inmax), INiso3(wwmax,inmax), INch(wwmax,inmax),
 *    &           INlcoll(wwmax,inmax), INcoll(wwmax,inmax),
 *    &           INistr(wwmax,inmax), INorigin(wwmax,inmax)
 *     common /historyin/ INind, INr0, INrx, INry, INrz,
 *    &                   INp0, INpx, INpy, INpz, INmass, INityp,
 *    &                   INiso3, INch, INlcoll, INcoll, INistr, INorigin
 *     C++ layout: all (20000,10) → [10][20000], accessor [j-1][i-1]
 *
 *   /historyout/ — outgoing-particle collision history  (wwmax=20000, outmax=50)
 *     Same structure as /historyin/ but prefixed OUT and with outmax=50.
 *     integer   OUTind(wwmax,outmax), OUTityp, OUTiso3, OUTch,
 *    &           OUTlcoll, OUTcoll, OUTistr, OUTorigin  (all (wwmax,outmax))
 *     real*8    OUTr0, OUTrx, OUTry, OUTrz, OUTp0, OUTpx, OUTpy, OUTpz,
 *    &           OUTmass  (all (wwmax,outmax))
 *     common /historyout/ OUTind, OUTr0, OUTrx, OUTry, OUTrz,
 *    &                    OUTp0, OUTpx, OUTpy, OUTpz, OUTmass, OUTityp,
 *    &                    OUTiso3, OUTch, OUTlcoll, OUTcoll, OUTistr, OUTorigin
 *     C++ layout: all (20000,50) → [50][20000], accessor [j-1][i-1]
 *
 *   /historyrec/ — reconstructable-resonance record  (resomax=40000)
 *     integer   resctr, www
 *     real*8    RECr0(resomax), RECrx(resomax), RECry(resomax),
 *    &           RECrz(resomax), RECp0(resomax), RECpx(resomax),
 *    &           RECpy(resomax), RECpz(resomax), RECmass(resomax)
 *     integer   RECityp(resomax), RECiso3(resomax), RECch(resomax),
 *    &           REClcoll(resomax), RECcoll(resomax), RECorigin(resomax)
 *     integer   RECind(resomax,0:4)    ! NOTE: second index is 0-based
 *     real*8    ORIr0(resomax), ORIrx(resomax), ORIry(resomax),
 *    &           ORIrz(resomax), ORIp0(resomax), ORIpx(resomax),
 *    &           ORIpy(resomax), ORIpz(resomax), ORImass(resomax)
 *     integer   ORIityp(resomax), ORIiso3(resomax), ORIch(resomax),
 *    &           ORIlcoll(resomax), ORIcoll(resomax), ORIorigin(resomax)
 *     common /historyrec/ resctr, www, RECr0, RECrx, RECry, RECrz,
 *    &                    RECp0, RECpx, RECpy, RECpz, RECmass,
 *    &                    RECityp, RECiso3, RECch, REClcoll, RECcoll,
 *    &                    RECorigin, RECind,
 *    &                    ORIr0, ORIrx, ORIry, ORIrz, ORIp0, ORIpx,
 *    &                    ORIpy, ORIpz, ORImass, ORIityp, ORIiso3, ORIch,
 *    &                    ORIlcoll, ORIcoll, ORIorigin
 *     C++ layout: RECind(resomax,0:4) → _RECind[5][40000]
 *                 accessor: RECind(i,j) = _RECind[j][i-1]  (j already 0-based)
 *
 * ----------------------------------------------------------------------------
 * From comstr.inc
 * ----------------------------------------------------------------------------
 *
 *   /FRGSPA/ — scalar string-fragmentation parameters  (njspin=8)
 *     real*8 PJSPNS, PMIX1S(3,njspin), PMIX2S(3,njspin), PBARS, PARQLS, PARRS
 *     COMMON /FRGSPA/ PJSPNS, PMIX1S, PMIX2S, PBARS, PARQLS, PARRS
 *     C++ layout: PMIX1S(3,8) → _PMIX1S[8][3], accessor [j-1][i-1]
 *
 *   /FRGCPA/ — charm string-fragmentation parameters  (njspin=8)
 *     real*8 PJSPNC, PMIX1C(3,njspin), PMIX2C(3,njspin), PBARC
 *     COMMON /FRGCPA/ PJSPNC, PMIX1C, PMIX2C, PBARC
 *     C++ layout: PMIX1C(3,8) → _PMIX1C[8][3], accessor [j-1][i-1]
 *
 *   /coparm/ — combinatorial-background fragmentation probabilities (njspin=8)
 *     real*8 parm(njspin)
 *     COMMON /coparm/ parm
 *
 *   /CONST/ — runtime value of pi  (distinct from the parameter pi in coms.inc)
 *     real*8 PI
 *     COMMON /CONST/ PI
 *
 * ----------------------------------------------------------------------------
 * From freezeout.inc
 * ----------------------------------------------------------------------------
 *
 *   /frcoor/ — freeze-out space-time coordinates  (nmax=100000)
 *     real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
 *    &       frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)
 *     common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz
 *
 * ----------------------------------------------------------------------------
 * From outcom.inc
 * ----------------------------------------------------------------------------
 *
 *   /outcom/ — per-collision output buffer  (3-particle slots)
 *     real*8 tsqrts, tstot, tsigpart
 *     real*8 tr0(3), trx(3), try(3), trz(3), ttform(3), txtotfac(3),
 *    &       tp0(3), tpx(3), tpy(3), tpz(3), tm(3)
 *     integer tind(3), tityp(3), tiso3(3), tcoll(3), tstrange(3),
 *    &        tlcoll(3), tcharge(3), torigin(3), tuid(3)
 *     common /outcom/ tsqrts, tstot, tsigpart,
 *    &                tr0, trx, try, trz, ttform, txtotfac,
 *    &                tp0, tpx, tpy, tpz, tm,
 *    &                tind, tityp, tiso3, tcoll, tstrange,
 *    &                tlcoll, tcharge, torigin, tuid
 *     WARNING: 'try' is a C++ reserved keyword — the accessor must be try_()
 *
 * ============================================================================
 */

#endif
