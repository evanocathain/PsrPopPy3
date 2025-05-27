f(x)=(1.0/(1.0+0.0473*x**(-0.627)))
set ylabel "FFT efficiency factor"
set xlabel "Duty cycle, {/Symbol d}=W/P"
set logscale x
plot [0.0001:0.5][0:1]f(x) notitle

set terminal postscript enhanced color solid
set output "fft_factor.ps"
replot
