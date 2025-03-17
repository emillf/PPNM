set terminal png
set output 'Quadratic_splines.png'
set title 'Quadratic spline interpolation, integration, and derivative of Sin(x)'
set xlabel 'x'
set ylabel 'y'
set grid
set key top right
set xrange [0:8]
set yrange [-1.5:1.5]
plot 'Vals.dat' using 1:2 with lines title 'Quadratic spline of Sin(x)', \
     'Vals.dat' using 1:3 with lines title 'Quadratic spline integration of Sin(x) (-Cos(x))', \
