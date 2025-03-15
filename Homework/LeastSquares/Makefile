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
plot: Activity_plot
	
Activity_plot : Vals.dat Activity_plot.gnu
	gnuplot Activity_plot.gnu
.PHONY: Activity_plot.gnu
Activity_plot.gnu:
	@echo "set terminal png" > Activity_plot.gnu
	@echo "set output 'ThX_Activity.png'" >> Activity_plot.gnu
	@echo "set title 'Least squares fit of measured activity of ThX'" >> Activity_plot.gnu
	@echo "set xlabel 't(days)'" >> Activity_plot.gnu
	@echo "set ylabel 'Relative activity'" >> Activity_plot.gnu
	@echo "set grid" >> Activity_plot.gnu
	@echo "set key top right" >> Activity_plot.gnu
	@echo "set xrange [0:14]" >> Activity_plot.gnu
	@echo "set yrange [0:*]" >> Activity_plot.gnu
	@echo "plot 'Vals.dat' using 1:2 with lines title 'Fitted curve', \\" >> Activity_plot.gnu
	@echo "     'Vals.dat' using 3:4:5 with yerrorbars title 'Experimental data with errorbars', \\" >> Activity_plot.gnu
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
