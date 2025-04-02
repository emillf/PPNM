libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
plot : erf_plot
	
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
.PHONY: erf_plot.gnu
erf_plot.gnu:
	@echo "set terminal png" > $@
	@echo "set output 'erf(x)_comparison.png'" >> $@
	@echo "set title 'Comparison of erf(x) with tabulated values'" >> $@
	@echo "set xlabel 'x'" >> $@
	@echo "set ylabel 'y'" >> $@
	@echo "set grid" >> $@
	@echo "set key top right" >> $@
	@echo "set xrange [0:1]" >> $@
	@echo "set yrange [*:*]" >> $@
	@echo "plot 'Vals1.dat' using 1:2 with lines title 'Tabulated values of erf(x)', \\" >> $@
	@echo "     'Vals1.dat' using 1:3 with lines title 'erf(x)' \\" >> $@
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
