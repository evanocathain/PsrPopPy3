
# This plots scattering relation at 190 MHz from Bhat et al. (2004) a la psrpoppy
# and the Krishnakumar et al. 2015 relation a la Vanessa's code

ramesh(x) = 10**(-6.46 + 0.154*log(x)/log(10) + 1.07*(log(x)/log(10))*(log(x)/log(10)) - 3.86*log(freq)/log(10))
kk(x)=1000*(327.0/190)**4.4*(3.6e-9*x**2.2*(1+1.94e-3*x**2))

set xlabel "DM (pc/cc)"
set ylabel "{/Symbol t}_s (ms)"
plot [1:200]kk(x), ramesh(x)
 

