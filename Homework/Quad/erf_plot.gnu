set terminal png
set output 'erf(x)_comparison.png'
set title 'Comparison of erf(x) with tabulated values'
set xlabel 'x'
set ylabel 'y'
set grid
set key top right
set xrange [0:1]
set yrange [*:*]
plot 'Vals1.dat' using 1:2 with lines title 'Tabulated values of erf(x)', \
     'Vals1.dat' using 1:3 with lines title 'erf(x)' \
