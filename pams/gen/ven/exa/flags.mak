ifeq ($(STAR_ARCH),sun4)
ARFLAGS=rl
CFLAGS=-g -sb -D SUN
FFLAGS=-g -xl -e -sb -C -time -v
RANLIB=ranlib
endif
ifeq ($(STAR_ARCH),hpux)
ARFLAGS=r
CFLAGS=-c -D HPUX
FFLAGS=-C +e +es +ppu +E0 -K -g
RANLIB=ranlib
endif
ifeq ($(STAR_ARCH),aix)
ARFLAGS=srv 
# -C :  array bound checking
# -qflttrap=ov:und:zero:inv:en   floating point exeception handling
# -Pk -Wp,-tr=NV,-INTL,-SV=A,-CTYPLSS   Pk precompiler stuff
FFLAGS=-g -C -w -qextname -qextchk -qsave -qfixed=132 -qctyplss -qintlog
FPPFLAGS=-DAIX
CCFLAGS=-c -DAIX 
RANLIB=ranlib
endif
ifeq ($(STAR_ARCH),irix)
ARFLAGS=slrv
FFLAGS=-extend_source -C -strictIEEE -g -static -v
CCFLAGS=-c -DIRIX
RANLIB=ar ts
endif













