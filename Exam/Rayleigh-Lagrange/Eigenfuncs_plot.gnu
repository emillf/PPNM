set terminal svg size 1000,800 background 'white'
set output 'Eigenfuncs_plot.svg'
set title 'Calculated lowest energy reduced radial hydrogen eigenfunctions vs real eigenfunctions'
set ylabel 'u(r)'
set xlabel 'r'
set grid
set xrange [0:20]
set multiplot layout 3,1 rowsfirst
plot 'Eigenfuncs.dat' using 1:2 with lines title 'calculated u_{10}(r)', 'Eigenfuncs.dat' using 1:5 with lines lw 2 dashtype 2 title 'real u_{10}(r)'
plot 'Eigenfuncs.dat' using 1:3 with lines title 'calculated u_{20}(r)', 'Eigenfuncs.dat' using 1:6 with lines lw 2 dashtype 2 title ' real u_{20}(r)'
plot 'Eigenfuncs.dat' using 1:4 with lines title 'calculated u_{21}(r)', 'Eigenfuncs.dat' using 1:7 with lines lw 2 dashtype 2 title ' real u_{21}(r)'
