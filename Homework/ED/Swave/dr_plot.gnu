set terminal png
set output 'dr_convergence.png'
set title 'Convergence of Ground State Energy with Δr'
set xlabel 'Δr'
set ylabel 'ε₀'
set grid
set key top right
set xrange [0:*]
set yrange [-0.6:0.4]
set label 'r_{max} = 10' at graph 0.7, 0.8
plot 'dr_convergence.dat' using 1:2 with linespoints pt 7 title 'Ground State Energy', \
     -0.5 with lines lt 2 title 'Exact Value'
