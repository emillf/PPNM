set terminal svg size 1200,500 background 'white' 
set output 'unit_circle.svg'
set title 'Area of unit circle with error estimate as function of Npoints in plainmc'
set xlabel 'Npoints'
set ylabel 'Area'
set grid
set multiplot layout 1,2 rowsfirst
plot 'Vals1.dat' using 1:2 with points pointtype 6 title 'Estimated error', \
     'Vals1.dat' using 1:3 with points pointtype 2 title 'Exact error', \
     'Vals1.dat' using 1:4 with lines title '1/Sqrt(N)'
plot 'Vals1.dat' using 1:5 with points pointtype 6 title 'Calculated area', \
     'Vals1.dat' using 1:6 with lines title 'Real area (Pi)'
unset multiplot
