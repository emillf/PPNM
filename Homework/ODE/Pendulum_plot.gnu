set terminal png
set output 'Pendulum_with_friction.png'
set title 'Numerical Solution to pendulum with friction differential equation from scipy odeint'
set xlabel 't'
set ylabel 'theta'
set grid
set key top right
set xrange [*:*]
set yrange [*:*]
plot 'Pend.dat' using 1:2 with lines title '{\Symbol theta}(t)', \
