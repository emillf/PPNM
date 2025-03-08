set terminal png
set output 'rmax_convergence.png'
set title 'Convergence of Ground State Energy with r_{max}'
set xlabel 'r_{max}'
set ylabel 'ε₀'
set grid
set key bottom right
set yrange [*:0]
set label 'Δr = ' at graph 0.7, 0.8
plot 'rmax_convergence.dat' using 1:2 with linespoints pt 7 title 'Ground State Energy', \
     -0.5 with lines lt 2 title 'Exact Value'
