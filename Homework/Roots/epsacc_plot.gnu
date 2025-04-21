set terminal svg size 1200,500 background 'white' 
set output 'acc_eps_convergence.svg'
set title 'Convergence of E as a function of acc or eps'
set xlabel 'acc/eps'
set ylabel 'E'
set grid
set xrange [0:0.07]
plot 'Convacceps.dat' using 1:2 with lines title 'fixed eps = 0.05', \
     'Convacceps.dat' using 3:4 with lines title 'fixed acc = 0.05', \
