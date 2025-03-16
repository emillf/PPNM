set terminal png
set output 'Other_decay_functions_with_errors.png'
set title 'Other decay functions'
set xlabel 't(days)'
set ylabel 'Relative activity'
set grid
set key top right
set xrange [0:14]
set yrange [0:*]
plot 'OtherVals.dat' using 1:2 with lines title '(λ+δ,a+δ)', \
     'OtherVals.dat' using 1:3 with lines title '(λ-δ,a-δ)', \
     'OtherVals.dat' using 1:4 with lines title '(λ+δ,a-δ)', \
     'OtherVals.dat' using 1:5 with lines title '(λ-δ,a+δ)', \
