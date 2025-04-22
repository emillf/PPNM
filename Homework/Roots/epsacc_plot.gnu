set terminal svg size 1200,500 background 'white' 
set output 'acc_eps_convergence.svg'
set title 'Convergence of E as a function of acc or eps with (rmin,rmax)=(0.3,5.0)'
set xlabel 'acc/eps'
set ylabel 'E_0 [Hartree]'
set grid
set xrange [0:0.07]
plot 'Convacceps.dat' using 1:2 with points pointtype 2 title 'fixed eps = 0.03', \
     'Convacceps.dat' using 3:4 with points pointtype 6 title 'fixed acc = 0.03', \
