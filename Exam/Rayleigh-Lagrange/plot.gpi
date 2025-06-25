set terminal svg size 800,600 background 'white'
set output 'Time_plot.svg'
set title 'Time to find eigenvalue and eigenvector of n x n matrix using Lagrange multipliers'
set xlabel 'n (logscale)'
set ylabel 't [s] (logscale)'
set grid
set logscale x
set logscale y
set xrange [40:400]
plot 'Times.dat' using 1:2 with lines title 'measured time in seconds'
