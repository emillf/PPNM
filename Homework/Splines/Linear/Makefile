libs = ../../matlib
	
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
plot: LinearSpline_plot
#	
LinearSpline_plot : Vals.dat LinearSpline_plot.gnu
	gnuplot LinearSpline_plot.gnu
.PHONY: LinearSpline_plot.gnu
LinearSpline_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'Linear_splines.png'" >> $@
	@echo "set title 'Linear spline interpolation and integration of cos(x)'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'y'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [0:8]" >> $@
	@echo "set yrange [-1.5:1.5]" >> $@
	@echo "plot 'Vals.dat' using 1:2 with lines title 'Linear spline of cos(x)', \\" >> $@
	@echo "     'Vals.dat' using 3:4 with lines title 'Linear spline integration of cos(x)', \\" >> $@
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
