libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
genlist.dll : $(libs)/genlist.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll genlist.dll
	$(MKEXE)
plot: Test_plot Pendulum_plot Orbit_plot
	
Test_plot : Test.dat Test_plot.gnu
	gnuplot Test_plot.gnu
Pendulum_plot: Pend.dat Pendulum_plot.gnu
	gnuplot Pendulum_plot.gnu
Orbit_plot : Orbit.dat Orbit_plot.gnu
	gnuplot Orbit_plot.gnu
.PHONY: Test_plot.gnu
Test_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'Test_differential_equation.png'" >> $@
	@echo "set title 'Numerical Solution of d^2u/dx^2=-u with yinit=(0,1)'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'y'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [*:*]" >> $@
	@echo "set yrange [*:*]" >> $@
	@echo "plot 'Test.dat' using 1:2 with points title 'u', \\" >> $@
	@echo  "    'Test.dat' using 1:3 with lines lw 2 dashtype 2 title 'Sin(x)', \\" >> $@
Pendulum_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'Pendulum_with_friction.png'" >> $@
	@echo "set title 'Numerical Solution to pendulum with friction differential equation from scipy odeint'" >> $@
	@echo "set xlabel 't'" >> $@
	@echo "set ylabel 'theta'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [*:*]" >> $@
	@echo "set yrange [*:*]" >> $@
	@echo "plot 'Pend.dat' using 1:2 with lines title '{\Symbol theta}(t)', \\" >> $@
Orbit_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'Planetary_precession.png'" >> $@
	@echo "set title 'Numerical calculation of planetary precession per general relativity'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'y'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [-2.5:2.5]" >> $@
	@echo "set yrange [-2.5:2.5]" >> $@
	@echo "plot 'Orbit.dat' using 1:2 with lines title 'ε = 0, Circular', \\" >> $@
	@echo "     'Orbit.dat' using 3:4 with lines title 'ε = 0, Elliptical', \\" >> $@
	@echo "     'Orbit.dat' using 5:6 with lines title 'ε = 0.01, Elliptical'" >> $@
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
