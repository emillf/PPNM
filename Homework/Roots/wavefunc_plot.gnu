set terminal svg size 500,500 background 'white' 
set output 'Wavefunc.svg'
set title 'Plot of odesolver wavefunction at E_0 compared to exact wavefunction'
set xlabel 'r [a_0]'
set ylabel 'f'
set grid
set xrange [0:8]
plot 'wavefunc.dat' using 1:2 with lines title 'exact wavefunction r*exp(-r)', \
     'wavefunc.dat' using 3:4 with points pointtype 6 title 'Numerically calculated wavefunction', \
