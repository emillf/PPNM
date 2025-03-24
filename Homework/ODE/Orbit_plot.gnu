set terminal png
set output 'Planetary_precession.png'
set title 'Numerical calculation of planetary precession per general relativity'
set xlabel 'x'
set ylabel 'y'
set grid
set key top right
set xrange [-2.5:2.5]
set yrange [-2.5:2.5]
plot 'Orbit.dat' using 1:2 with lines title 'ε = 0, Circular', \
     'Orbit.dat' using 3:4 with lines title 'ε = 0, Elliptical', \
     'Orbit.dat' using 5:6 with lines title 'ε = 0.01, Elliptical'
