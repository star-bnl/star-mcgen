# STAR UrQMD 4.0 bridge

This directory contains the STAR bridge for the UrQMD 4.0 event generator. The bridge is centered on `StarUrQMD.cxx`, which adapts the FORTRAN generator to the `StarGenerator` C++ framework used by STAR simulation and embedding jobs.

At a high level, the bridge works in three layers:

1. `StarUrQMD.cxx` / `StarUrQMD.h` steer the run from C++, create the STAR event record, and expose STAR-style configuration such as beam species, frame, and user-supplied `Set` options.
2. `UrQMD.h` / `UrQMD.cxx` provide C++ overlays for UrQMD FORTRAN common blocks, plus symbol binding for FORTRAN entry points.
3. The FORTRAN side (`urqmd.F90`, `input.F`, `ityp2pdg.F`, `urqmd2pdg.F90`, and the rest of the generator) performs initialization, time evolution, decays, and particle bookkeeping.

## Main bridge files

The files most relevant to the STAR interface are:

- `StarUrQMD.h`  
  Declares the `StarUrQMD` class and exposes accessor methods for each FORTRAN common block, for example `sys()`, `coor()`, `isys()`, and `itdelay()`.

- `StarUrQMD.cxx`  
  Implements initialization, input-file generation, event generation, event-record filling, random-number remapping, and UrQMD-to-PDG conversion calls.

- `UrQMD.h`  
  Declares C++ `struct`s whose memory layout mirrors the FORTRAN common blocks. It also provides 1-based access helpers so the C++ code can read FORTRAN arrays without rewriting UrQMD internals.

- `UrQMD.cxx`  
  Binds the common-block symbols and FORTRAN routines with `F77_NAME(...)`.

- `urqmd.F90`  
  Exposes the C-callable entry points `urqmd_init_` and `urqmd_make_`, which initialize UrQMD and generate one event, respectively.

- `input.F`  
  Reads `UrQMD.in` and fills UrQMD control state, including the `stb` commands used to keep selected internal particle species stable.

- `ityp2pdg.F` and `urqmd2pdg.F90`  
  Convert UrQMD internal particle identifiers `(ityp, iso3)` into PDG IDs.

## How FORTRAN common blocks are accessed from C++

The bridge uses the standard STAR pattern for older FORTRAN generators: it overlays the FORTRAN common blocks with C++ `struct`s and reads them directly after each generator call.

`UrQMD.h` declares a `struct` for each common block. A few important examples are:

- `SYS_t` for `/sys/`  
  Run- and event-level counters such as `npart`, `event`, `Ap`, `At`, `Zp`, `Zt`, and `success`.

- `COOR_t` for `/coor/`  
  Particle four-vectors and coordinates such as `p0`, `px`, `py`, `pz`, `fmass`, `r0`, `rx`, `ry`, `rz`.

- `ISYS_t` for `/isys/`  
  Per-particle identity and bookkeeping such as `ityp`, `iso3`, `charge`, `uid`, `origin`, and `ncoll`.

- `ITDELAY_t` and `RTDELAY_t`  
  Delayed-decay parent information used by the bridge when it fills mother fields in the STAR event record.

Each common block has an `extern "C"` accessor such as `address_of_sys()` or `address_of_coor()`. `StarUrQMD.h` wraps those accessors in convenience methods:

```c++
SYS_t  &sys()  { return *address_of_sys();  }
COOR_t &coor() { return *address_of_coor(); }
ISYS_t &isys() { return *address_of_isys(); }
```

The array accessors inside the structs intentionally use FORTRAN-style indexing. For example:

```c++
Int_t    &ityp( Int_t i )   { return _ityp[i-1]; }
Double_t &px  ( Int_t i )   { return _px[i-1];   }
```

This lets the bridge write loops like `for (Int_t idx=1; idx<=sys().npart; ++idx)` and read UrQMD storage directly without index translation bugs.

`UrQMD.cxx` binds the underlying FORTRAN symbols with `F77_NAME(...)`, so the same bridge can tolerate compiler-specific name mangling.

## How the simulation workflow is steered

### 1. STAR-side configuration

`StarUrQMD::Init()` reads standard `StarGenerator` maker attributes:

- `FRAME`
- `BLUE`
- `YELL`
- `Ecms`
- optional `Set`

Beam species are converted into UrQMD projectile and target definitions through the local `A` and `Z` maps. The bridge writes those values into the string parameters:

- `pro` for projectile `(A, Z)`
- `tar` for target `(A, Z)`
- `IMP` for impact-parameter range

For CMS running, the bridge writes `ecm`. For momentum-defined frames, it estimates `mRootS` from the registered beam momenta and nucleon masses, then also writes `ecm`.

### 2. Input-file generation

UrQMD is not configured by direct C++ setter calls. Instead, `StarUrQMD::InitializeUrQMD()` writes a text file named `UrQMD.in`, which `input.F` later parses.

The generated file contains:

- integer, floating-point, and string parameters from the `InputParameters*` maps
- one `stb <ityp>` line per particle species that should be kept stable
- `tim 200 200`
- output-file flags (`f13`, `f15`, `f16`, `f19`, `f20`, etc.) to suppress or enable UrQMD text output
- `nev 1000`
- `rsd 111`
- any additional commands passed through the STAR maker `Set` attribute
- the end marker `xxx and done`

On the FORTRAN side, `input.F` opens `UrQMD.in` on unit 9 and parses those flags. The `stb` records are stored in the common block:

- `nstable`
- `stabvec(maxstables)`

### 3. Generator initialization

After writing `UrQMD.in`, the C++ bridge calls:

```c++
urqmd_init();
```

This enters `urqmd.F90`, where `urqmd_init()` simply calls `uinit(0)`. That initializes UrQMD's technical state, reads the input, and prepares the generator for event production.

### 4. Per-event generation

`StarUrQMD::Generate()` drives one event by calling:

```c++
urqmd_make();
```

`urqmd_make()` in `urqmd.F90` performs the full UrQMD event evolution:

- initializes event-level state
- seeds the generator if needed
- calls `init` and `init_eccentricity`
- iterates over time steps
- loads and updates the collision table
- propagates particles
- performs scatterings and decays
- optionally forces final decays for unstable particles at the end of the event

When control returns to C++, the particle record is already present in the FORTRAN common blocks. `StarUrQMD::Generate()` then:

1. creates or fills the STAR event header with `FillPP()` or `FillAA()`
2. reads the particle count from `sys().npart`
3. loops over particles `idx = 1 .. npart`
4. reads identity from `isys()`
5. reads kinematics from `coor()`
6. reads mother information from `itdelay()`
7. converts the UrQMD particle code to a PDG code
8. appends a `StarGenParticle` to the STAR event record with `mEvent->AddParticle(...)`

The bridge marks all produced particles as `StarGenParticle::kFinal`, because it does not reconstruct a richer UrQMD-specific status scheme.

## Random-number handling

`StarUrQMD.cxx` remaps UrQMD's random-number entry points to `StarRandom`:

```c++
Double_t ranfstar_( Int_t *idummy ){ return StarRandom::Instance().flat(); }
Double_t pyrstar_ ( Int_t *idummy ){ return StarRandom::Instance().flat(); }
```

This matches the pattern used by other STAR generator bridges and keeps random-number control centralized in the framework.

## Particle ID handling

### UrQMD internal IDs

UrQMD does not store PDG IDs in its event record. Instead, each particle is described by:

- `ityp` — the UrQMD particle type
- `iso3` — the third component of isospin in UrQMD's convention

Those values live in the `/isys/` common block and are read in C++ as:

```c++
int& ityp = isys().ityp(idx);
int& iso3 = isys().iso3(idx);
```

### Conversion to PDG

The bridge converts every UrQMD particle to the STAR-standard PDG code by calling:

```c++
int id = urqmd_pdgid(&ityp, &iso3);
```

`urqmd_pdgid` is a thin C-binding wrapper in `urqmd2pdg.F90` around the FORTRAN function `pdgid()` implemented in `ityp2pdg.F`.

`pdgid()` uses a static lookup table of `(ityp, iso3, pdgid)` triples. A few examples are:

- `(1, -1) -> 2112` for the neutron
- `(1,  1) -> 2212` for the proton
- `(101, 0) -> 111` for `pi0`
- `(101, 2) -> 211` for `pi+`
- `(106, -1) -> 311` for neutral kaons
- `(106,  1) -> 321` for charged kaons
- `(54, -2/0/2) -> 3112/3212/3222` for the sigma triplet

The FORTRAN converter also handles:

- antiparticles through negative `ityp`
- the corresponding `iso3` sign flip for anti-states
- PYTHIA-style embedded IDs when `abs(ityp) > 1000`, in which case the offset `1000` is stripped directly

### Neutral kaons in the C++ bridge

UrQMD generates `K0` / `K0bar` at the internal-ID stage. After conversion, the bridge checks for PDG `311` or `-311` and replaces the particle with a 50/50 choice of:

- `130` (`K0L`)
- `310` (`K0S`)

This is the same convention used in other STAR bridges such as `StarHijing.cxx`, and it is how the STAR event record receives the neutral kaon state expected downstream.

## Stable-particle handling

### The STAR convention

The standard stable-particle list used in the STAR bridges is visible in `StarPythia6.cxx`, `StarPythia8.cxx`, and `StarHijing.cxx`. Those bridges disable decays for the following PDG states so that GEANT can decide whether and where to decay them:

- `111` `pi0`
- `211` `pi+/-`
- `221` `eta`
- `321` `K+/-`
- `310` `K0S`
- `130` `K0L`
- `3122` `Lambda0`
- `3112` `Sigma-`
- `3222` `Sigma+`
- `3212` `Sigma0`
- `3312` `Xi-`
- `3322` `Xi0`
- `3334` `Omega-`

### How UrQMD implements that list

`StarUrQMD::Init()` contains an unused `stablepdg` vector that documents exactly that same 13-particle STAR list. The actual control path, however, is the `StableParticles` string vector written to `UrQMD.in` through `stb` commands.

The bridge currently writes these UrQMD internal IDs:

- `101`, `-101` for the pion family
- `102`, `-102` for eta
- `106`, `-106` for kaons
- `41`, `-41` for lambda
- `54`, `-54` for sigma
- `63`, `-63` for xi
- `69`, `-69` for omega

This matches UrQMD's own species organization:

- `101` covers the full pion triplet through `iso3 = -2, 0, +2`
- `106` covers both neutral and charged kaons through `iso3 = -1, +1`
- `54` covers `Sigma-`, `Sigma0`, and `Sigma+`
- `63` covers `Xi-` and `Xi0`

On the FORTRAN side, `input.F` stores each `stb` entry in `stabvec`, and both `coload.F` and `urqmd.F90` test stability with:

```fortran
if (ityp(i).eq.stabvec(stidx)) then
   isstable = .true.
endif
```

That comparison is done on `ityp` only, not on `(ityp, iso3)`. For UrQMD this is exactly what is needed, because a single `ityp` value typically represents an isospin multiplet and `iso3` selects the member.

### Verification result

The bridge does implement the standard STAR stable-particle policy, with one UrQMD-specific nuance:

- the STAR-side intent matches the standard 13-particle list exactly
- the UrQMD side enforces stability by internal family code (`ityp`) rather than by PDG ID
- neutral kaons are kept stable as UrQMD kaons and only converted to `K0S`/`K0L` afterward in C++

So, in practice, the bridge preserves the standard STAR stable set. It does not set `310` and `130` directly inside UrQMD because those are not native UrQMD particle species.

Because UrQMD compares only `ityp(i)` against `stabvec`, the stability setting applies at the UrQMD family-code level, not at the `(ityp, iso3)` level. The signed entries therefore matter for species that can appear with negative `ityp` as anti-states, while the positive family IDs define the particle families that STAR wants to preserve for transport and later decay handling.

## Event-record details

For each final UrQMD particle, the bridge fills:

- status: always `kFinal`
- PDG ID: via `urqmd_pdgid()` and the neutral-kaon post-processing
- mothers: from `itdelay().ityptd(idx,1/2)`
- daughters: left as `0`
- momentum and energy: from `/coor/`
- mass: from `/coor/`
- production vertex: set to zero in the current bridge

UrQMD does carry more detailed transport history and coordinate information than the bridge exports. In particular, the bridge currently ignores freeze-out coordinates even though comments in `StarUrQMD.cxx` note that UrQMD stores them.

## Comparison with other STAR bridges

Relative to other STAR generator adapters in `StRoot/StarGenerator`:

- the attribute-driven steering (`FRAME`, `BLUE`, `YELL`, `Set`) is consistent with `StarPythia6`, `StarPythia8`, and `StarHijing`
- the stable-particle policy is consistent with `StarPythia6`, `StarPythia8`, and `StarHijing`
- the random-number remapping follows the same STAR pattern used by other FORTRAN-backed bridges
- the biggest difference is that UrQMD is driven through a generated `UrQMD.in` file and common-block overlays instead of a direct C++ API or HEPEVT-style event record

## Current caveat in the checked-in bridge

Earlier development versions of this bridge used a diagnostic loop in `StarUrQMD::Init()` to print the UrQMD-to-PDG lookup coverage through `StarParticleData` and then exited early. That diagnostic `exit(0)` has been removed from the checked-in `StarUrQMD.cxx`, so the bridge now proceeds from initialization into ordinary event generation as described in the rest of this document.

If you reintroduce similar diagnostics in local development, ensure they are either disabled or properly guarded before using the bridge in production simulation or embedding jobs.
