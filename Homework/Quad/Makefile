libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
plot : erf_plot acc_plot
	
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
genlist.dll : $(libs)/genlist.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll genlist.dll
	$(MKEXE)
erf_plot : Vals1.dat erf_plot.gnu
	gnuplot erf_plot.gnu
acc_plot: Vals2.dat acc_plot.gnu
	gnuplot acc_plot.gnu
.PHONY: erf_plot.gnu
erf_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'erf(x)_comparison.png'" >> $@
	@echo "set title 'Comparison of erf(x) with tabulated values'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'erf(x)'" >> $@
	@echo "set grid" >> $@
	@echo "set key bottom right" >> $@
	@echo "set xrange [-3.5:3.5]" >> $@
	@echo "set yrange [-1.05:1.05]" >> $@
	@echo "plot 'Tabvals.txt' using 1:2 with points pointtype 6 title 'Tabulated values of erf(x)', \\" >> $@
	@echo "     'Vals1.dat' using 1:2 with lines title 'erf(x)' \\" >> $@
.PHONY: acc_plot.gnu
acc_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'acc_convergence.png'" >> $@
	@echo "set title 'Difference between erf(1) and tabulated value as a function of acc'" >> $@
	@echo "set xlabel 'log_{10}(acc)'" >> $@
	@echo "set ylabel 'log_{10}(|erf(x)-exact value|)'" >> $@
	@echo "set grid" >> $@
	@echo "set key bottom right" >> $@
	@echo "set xrange [*:*]" >> $@
	@echo "set yrange [*:*]" >> $@
	@echo 'set format x "%2.0t{×}10^{%L}"' >> $@
	@echo 'set format y "%2.0t{×}10^{%L}"' >> $@
	@echo "set logscale x 10" >> $@
	@echo "set logscale y 10" >> $@
	@echo "plot 'Vals2.dat' using 1:2 with lines title '|erf(x)-exact value|', \\" >> $@
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
