#!/bin/csh
#
#    20 proton+proton events 
#    ===================================
setenv Z_TAPE "`pwd`/test_pp.ven"
#
# make the optns dat for venus...
#
echo "'had'">optns.dat
echo "1 1  -200. 1 1 1 1 3">>optns.dat
echo "'irescl' 0">>optns.dat
echo "'bmaxim' 1.00">>optns.dat
echo "'ish' 0">>optns.dat
echo "'stp'">>optns.dat
#
time $STAR_LIB/star/$STAR_LEVEL/sim/bin/$STAR_ARCH/venus307 <<eof
eof
#
unsetenv Z_TAPE 
if ("`alias rm `" == "rm -i") then
   unalias rm
   rm check.dat histo.dat optns.dat
   alias rm 'rm -i'
else
   rm check.dat histo.dat optns.dat
endif

