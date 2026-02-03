#!/bin/sh

# fix version number to the content of file version:
version=$(head -n 1 version)
number=$(tail -n 1 version)
# rename "aa- UrQMD..."-file:
mv aa-* "aa- UrQMD $version -"
# next: README.
perl -pi -e "s/UrQMD version .*$/UrQMD version $version/" README
# blockres.f:
perl -pi -e "s{UrQMD: Version \S+ \([\d/]+\)}{UrQMD: Version $version ($number)}" blockres.f
# coms.f (three lines need to be changed):
perl -pi -e "s{version\s+=\s+\d+}{version = $number}" coms.f
perl -pi -e "s{laires\s+=\s+\d+}{laires  = $number}" coms.f
perl -pi -e "s{versiontxt\s+=\s+\'[^']*'}{versiontxt = '$version'}" coms.f
# user guide (needs to be recompiled):
cd doc/
perl -pi -e "s/\\\\newcommand\{\\\\uversion.+\$/\\\\newcommand{\\\\uversion}{\\\\textbf{$version}}/" urqmd-user.tex
make ps
make pdf
cd ..

tar czvvf urqmd-$version.tar.gz --transform "s,^,urqmd-$version/," \
        version \
        Copyright \
        README \
        ChangeLog \
        doc/urqmd-user.ps \
        doc/urqmd-user.pdf \
        *.f *.f90 \
        lhc.patch \
        inputfile.example \
        GNUmakefile \
        mk/ \
        maketables \
        runqmd.bash \
        eosfiles/

exit;
