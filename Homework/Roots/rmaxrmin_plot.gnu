set terminal svg size 1200,500 background 'white' 
set output 'rmax_rmin_convergence.svg'
set title 'Convergence of E_0 as a function of rmax or rmin with acc = eps = 0.03'
set xlabel 'rmin [a_0]'
set ylabel 'E_0 [Hartree]'
set grid
set xrange [*:*]
set multiplot layout 1,2 rowsfirst ;\
plot \
'Convrminrmax.dat' using 3:4 with points pointtype 6 title 'fixed rmax = 5.0' ;\
set xlabel 'rmax [a_0]' ;\
plot \
'Convrminrmax.dat' using 1:2 with points pointtype 6 title 'fixed rmin = 0.1' ;\
unset multiplot
