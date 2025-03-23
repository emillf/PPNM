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
plot: Test_plot Pendulum_plot
	
Test_plot : Test.dat Test_plot.gnu
	gnuplot Test_plot.gnu
Pendulum_plot: Pend.dat Pendulum_plot.gnu
	gnuplot Pendulum_plot.gnu
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
	@echo "plot 'Test.dat' using 1:2 with lines title 'u', \\" >> $@
Pendulum_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'Pendulum_with_friction.png'" >> $@
	@echo "set title 'Numerical Solution to pendulum with friction differential equation from scipy odeint'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'y'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [*:*]" >> $@
	@echo "set yrange [*:*]" >> $@
	@echo "plot 'Pend.dat' using 1:2 with lines title 'theta(t)', \\" >> $@
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
