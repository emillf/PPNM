libs = ../matlib
CSS = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
MKLIBS = mcs -target:library -out:$@ $(CSS)
default: Out.txt plot
	
plot: Eigenfuncs_plot Times_plot
	
Out.txt Times.dat Eigenfuncs.dat : main.exe
	mono main.exe > Out.txt
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll
	$(MKEXE)
Eigenfuncs_plot: Eigenfuncs.dat Eigenfuncs_plot.gnu
	gnuplot Eigenfuncs_plot.gnu
Times_plot : Times.dat Times_plot.gnu
	gnuplot Times_plot.gnu
.PHONY: Eigenfuncs_plot.gnu 
Eigenfuncs_plot.gnu:
	@echo "set terminal svg size 1000,800 background 'white'" > $@
	@echo "set output 'Eigenfuncs_plot.svg'" >> $@
	@echo "set title 'Calculated lowest energy reduced radial hydrogen eigenfunctions vs real eigenfunctions'" >> $@
	@echo "set ylabel 'u(r)'" >> $@
	@echo "set xlabel 'r'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [0:15]" >> $@
	@echo "set multiplot layout 3,1 rowsfirst" >> $@
	@echo "plot 'Eigenfuncs.dat' using 1:2 with lines title 'calculated u_{10}(r)', 'Eigenfuncs.dat' using 1:5 with lines lw 2 dashtype 2 title 'real u_{10}(r)'">>$@
	@echo "plot 'Eigenfuncs.dat' using 1:3 with lines title 'calculated u_{20}(r)', 'Eigenfuncs.dat' using 1:6 with lines lw 2 dashtype 2 title ' real u_{20}(r)'">>$@
	@echo "plot 'Eigenfuncs.dat' using 1:4 with lines title 'calculated u_{21}(r)', 'Eigenfuncs.dat' using 1:7 with lines lw 2 dashtype 2 title ' real u_{21}(r)'">>$@
.PHONY : Times_plot.gnu
Times_plot.gnu:
	@echo "set terminal svg size 800,600 background 'white'" > $@
	@echo "set output 'Time_plot.svg'" >> $@
	@echo "set title 'Time to find eigenvalue and eigenvector of n x n matrix using Lagrange multipliers'" >> $@
	@echo "set xlabel 'n (logscale)'" >> $@
	@echo "set ylabel 't [s] (logscale)'" >> $@
	@echo "set grid" >> $@
	@echo "set logscale x" >> $@
	@echo "set logscale y" >> $@
	@echo "set xrange [40:400]" >> $@
	@echo "plot 'Times.dat' using 1:2 with lines title 'measured time in seconds'">> $@
clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr log.gpi log1.gpi
