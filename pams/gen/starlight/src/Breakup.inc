      REAL*8 BMin,ProbTot(1000)
      Common/Breakup/Bmin,ProbTot
C Bmin is in units of fermi
C ProbTot is computed once per run, 
C depending on your choice of ibeakup, and then passed to integrand() 
C Choice dim(PTot)=1000 is hard-coded into PofB.f, vladimir.f and flux.f  

