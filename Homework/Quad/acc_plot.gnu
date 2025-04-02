set terminal png
set output 'acc_convergence.png'
set title 'Difference between erf(1) and tabulated value as a function of acc'
set xlabel 'log_{10}(acc)'
set ylabel 'log_{10}(|erf(x)-exact value|)'
set grid
set key bottom right
set xrange [*:*]
set yrange [*:*]
set format x "%2.0t{×}10^{%L}"
set format y "%2.0t{×}10^{%L}"
set logscale x 10
set logscale y 10
plot 'Vals2.dat' using 1:2 with lines title '|erf(x)-exact value|', \
