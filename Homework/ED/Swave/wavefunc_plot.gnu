set terminal png
set output 'Wavefuncs.png'
set title 'Hydrogen eigenfunctions'
set xlabel 'r'
set ylabel 'probability density'
set grid
set key bottom right
set yrange [*:*]
plot 'WaveFuncs.dat' using 1:2 with linespoints pt 7 title 'n=0', \
 'WaveFuncs.dat' using 1:3 with linespoints pt 7 title 'n=1', \
