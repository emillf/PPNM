set terminal png
set output 'Linear_splines.png'
set title 'Linear spline interpolation and integration of cos(x)'
set xlabel 'x'
set ylabel 'y'
set grid
set key top right
set xrange [0:8]
set yrange [-1.5:1.5]
plot 'Vals.dat' using 1:2 with lines title 'Linear spline of cos(x)', \
     'Vals.dat' using 3:4 with lines title 'Linear spline integration of cos(x)', \
