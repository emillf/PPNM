libs = ../matlib
CSS = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
MKLIBS = mcs -target:library -out:$@ $(CSS)
default: Out.txt Time_plot.svg
	
Out.txt Times.dat: main.exe
	mono main.exe > Out.txt
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll
	$(MKEXE)
Time_plot.svg: Times.dat
	@echo "set terminal svg size 800,600 background 'white'" > plot.gpi
	@echo "set output 'Time_plot.svg'" >> plot.gpi
	@echo "set title 'Time to find eigenvalue and eigenvector of n x n matrix using Lagrange multipliers'" >> plot.gpi
	@echo "set xlabel 'n (logscale)'" >> plot.gpi
	@echo "set ylabel 't [s] (logscale)'" >> plot.gpi
	@echo "set grid" >> plot.gpi
	@echo "set logscale x" >> plot.gpi
	@echo "set logscale y" >> plot.gpi
	@echo "set xrange [40:400]" >> plot.gpi
	@echo "plot 'Times.dat' using 1:2 with lines title 'measured time in seconds'">> plot.gpi
	gnuplot plot.gpi
	@cat plot.gpi >> log.gpi
	@rm -f plot.gpi
clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr plot.gpi
