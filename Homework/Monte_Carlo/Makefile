libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
plot : circle_plot
	
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
genlist.dll : $(libs)/genlist.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll genlist.dll
	$(MKEXE)
circle_plot : Vals1.dat circle_plot.gnu
	gnuplot circle_plot.gnu
#acc_plot: Vals2.dat acc_plot.gnu
#	gnuplot acc_plot.gnu
.PHONY: circle_plot.gnu
circle_plot.gnu:
	@echo "set terminal svg size 1200,500 background 'white' " > $@
	@echo "set output 'unit_circle.svg'" >> $@
	@echo "set title 'Area of unit circle with error estimate as function of Npoints in plainmc'" >> $@
	@echo "set xlabel 'Npoints'" >> $@
	@echo "set ylabel 'Area'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [-1:25001]" >> $@
	@echo "set multiplot layout 1,2 rowsfirst" >>$@
	@echo "plot 'Vals1.dat' using 1:2 with points pointtype 6 title 'Estimated error', \\" >> $@
	@echo "     'Vals1.dat' using 1:3 with points pointtype 2 title 'Exact error', \\" >> $@
	@echo "     'Vals1.dat' using 1:4 with lines title '1/Sqrt(N)'" >> $@
	@echo "plot 'Vals1.dat' using 1:5 with points pointtype 6 title 'Calculated area', \\" >> $@
	@echo "     'Vals1.dat' using 1:6 with lines title 'Real area (Pi)'" >> $@
	@echo "unset multiplot" >> $@
#.PHONY: acc_plot.gnu
#acc_plot.gnu:
#	@echo "set terminal png" > $@
#	@echo "set output 'acc_convergence.png'" >> $@
#	@echo "set title 'Difference between erf(1) and tabulated value as a function of acc'" >> $@
#	@echo "set xlabel 'log_{10}(acc)'" >> $@
#	@echo "set ylabel 'log_{10}(|erf(x)-exact value|)'" >> $@
#	@echo "set grid" >> $@
#	@echo "set key bottom right" >> $@
#	@echo "set xrange [*:*]" >> $@
#	@echo "set yrange [*:*]" >> $@
#	@echo 'set format x "%2.0t{×}10^{%L}"' >> $@
#	@echo 'set format y "%2.0t{×}10^{%L}"' >> $@
#	@echo "set logscale x 10" >> $@
#	@echo "set logscale y 10" >> $@
#	@echo "plot 'Vals2.dat' using 1:2 with lines title '|erf(x)-exact value|', \\" >> $@
clean:
	$(RM) *.exe *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
