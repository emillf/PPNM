set terminal png
set output 'erf(x)_comparison.png'
set title 'Comparison of erf(x) with tabulated values'
set xlabel 'x'
set ylabel 'erf(x)'
set grid
set key bottom right
set xrange [-3.5:3.5]
set yrange [-1.05:1.05]
plot 'Tabvals.txt' using 1:2 with points pointtype 6 title 'Tabulated values of erf(x)', \
     'Vals1.dat' using 1:2 with lines title 'erf(x)' \
