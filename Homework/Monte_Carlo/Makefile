libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
plot : circle_plot errscal_plot
	
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
errscal_plot: Vals2.dat errscal_plot.gnu
	gnuplot errscal_plot.gnu
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
.PHONY: errscal_plot.gnu
errscal_plot.gnu:
	@echo "set terminal svg size 1200,500 background 'white' " > $@
	@echo "set output 'errscal.svg'" >> $@
	@echo "set title 'Comparison of error scaling of quasimc and plainmc when calculating area of unit circle '" >> $@
	@echo "set xlabel 'Npoints'" >> $@
	@echo "set ylabel 'Error in area'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [-1:25001]" >> $@
	@echo "set multiplot layout 1,2 rowsfirst" >>$@
	@echo "plot 'Vals2.dat' using 1:2 with points pointtype 6 title 'Estimated error plainmc', \\" >> $@
	@echo "     'Vals2.dat' using 1:3 with points pointtype 2 title 'Exact error plainmc', \\" >> $@
	@echo "     'Vals2.dat' using 1:4 with lines title 'Error scaling plainmc'" >> $@
	@echo "plot 'Vals2.dat' using 1:5 with points pointtype 6 title 'Estimated error quasimc', \\" >> $@
	@echo "     'Vals2.dat' using 1:6 with points pointtype 2 title 'Exact error quasimc'" >> $@
	@echo "unset multiplot" >> $@
clean:
	$(RM) *.exe *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
