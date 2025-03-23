set terminal png
set output 'Test_differential_equation.png'
set title 'Numerical Solution of d^2u/dx^2=-u with yinit=(0,1)'
set xlabel 'x'
set ylabel 'y'
set grid
set key top right
set xrange [*:*]
set yrange [*:*]
plot 'Test.dat' using 1:2 with lines title 'u', \
