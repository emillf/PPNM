set terminal svg size 1200,500 background 'white' 
set output 'errscal.svg'
set title 'Comparison of error scaling of quasimc and plainmc when calculating area of unit circle '
set xlabel 'Npoints'
set ylabel 'Error in area'
set grid
set xrange [0:5001]
set multiplot layout 1,2 rowsfirst
plot 'Vals2.dat' using 1:2 with points pointtype 6 title 'Estimated error plainmc', \
     'Vals2.dat' using 1:3 with points pointtype 2 title 'Exact error plainmc', \
     'Vals2.dat' using 1:4 with lines title 'Error scaling plainmc'
plot 'Vals2.dat' using 1:5 with points pointtype 6 title 'Estimated error quasimc', \
     'Vals2.dat' using 1:6 with points pointtype 2 title 'Exact error quasimc'
unset multiplot
