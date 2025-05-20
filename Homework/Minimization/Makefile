libs = ../matlib
CSS = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
MKLIBS = mcs -target:library -out:$@ $(CSS)

default: Out.txt Higgsplot.svg

Out.txt : main.exe
	mono main.exe < higgs.data.txt 1> Out.txt 2> log

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)

main.exe : main.cs matrix.dll
	$(MKEXE)

Higgsplot.svg : Results.dat higgs.data.txt
	@echo "set terminal svg size 500,500 background 'white'" > plot.gpi
	@echo "set out '$@'" >> plot.gpi
	@echo "set title 'Plot of deviation-minimized Breit-Wigner function with experimental data'" >> plot.gpi
	@echo "set xlabel 'Energy [GeV]'" >> plot.gpi
	@echo "set ylabel 'signal [Certain units]'" >> plot.gpi
	@echo "set grid" >> plot.gpi
	@echo "set xrange [100:160]" >> plot.gpi
	@echo "plot 'Results.dat' with lines lt 8 lw 2 dashtype 2 title 'Breit-Wigner function', \\" >> plot.gpi
	@echo "     'higgs.data.txt' with errorbars pointtype 6 lc rgb 'red' title 'Experimental data'" >> plot.gpi
	gnuplot plot.gpi
	@cat plot.gpi>>log.gpi
	@rm -f plot.gpi

clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr plot.gpi
