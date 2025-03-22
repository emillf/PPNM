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
main.exe : main.cs matrix.dll
	$(MKEXE)
plot: QuadraticSpline_plot
#	
QuadraticSpline_plot : Vals.dat QuadraticSpline_plot.gnu
	gnuplot QuadraticSpline_plot.gnu
.PHONY: QuadraticSpline_plot.gnu
QuadraticSpline_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'Cubic_splines.png'" >> $@
	@echo "set title 'Quadratic spline interpolation, integration, and derivative of Sin(x)'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'y'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [0:8]" >> $@
	@echo "set yrange [-1.5:1.5]" >> $@
	@echo "plot 'Vals.dat' using 1:2 with lines title 'Cubic spline of Sin(x)', \\" >> $@
	@echo "     'Vals.dat' using 1:3 with lines title 'Cubic spline integration of Sin(x) (-Cos(x))', \\" >> $@
	@echo "     'Vals.dat' using 1:2 smooth csplines title 'Gnuplots own quadratic spline implimentation', \\" >> $@
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
