set terminal png
set output 'ThX_Activity.png'
set title 'Least squares fit of measured activity of ThX'
set xlabel 't(days)'
set ylabel 'Relative activity'
set grid
set key top right
set xrange [0:14]
set yrange [0:*]
plot 'Vals.dat' using 1:2 with lines title 'Fitted curve', \
     'Vals.dat' using 3:4:5 with yerrorbars title 'Experimental data with errorbars', \
