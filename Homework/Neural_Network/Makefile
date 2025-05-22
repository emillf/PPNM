libs = ../matlib
CSS = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
MKLIBS = mcs -target:library -out:$@ $(CSS)

default: Out.txt NN_plot.svg

Out.txt Vals.dat: main.exe
	mono main.exe > Out.txt

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs minimizer.cs integration.cs
	$(MKLIBS)

main.exe : main.cs matrix.dll
	$(MKEXE)
NN_plot.svg: Vals.dat
	@echo "set terminal svg size 1000,800 background 'white'" > plot.gpi
	@echo "set output 'NN_plot.svg'" >> plot.gpi
	@echo "set title 'Derivatives of f(x)=Cos(5x-1)*Exp(-x*x) and the neural network response NN(x)'" >> plot.gpi
	@echo "set xlabel 'x'" >> plot.gpi
	@echo "set ylabel 'y'" >> plot.gpi
	@echo "set grid" >> plot.gpi
	@echo "set xrange [-1:1]" >> plot.gpi
	@echo "set multiplot layout 2,2 rowsfirst" >> plot.gpi
	@echo "plot 'Vals.dat' using 1:2 with lines title 'f(x)', 'Vals.dat' using 1:3 with lines lw 2 dashtype 2 title 'NN(x)'" >> plot.gpi
	@echo "plot 'Vals.dat' using 1:4 with lines title 'first derivative of f(x)', 'Vals.dat' using 1:5 with lines lw 2 dashtype 2 title 'first derivative of NN(x)'" >> plot.gpi
	@echo "plot 'Vals.dat' using 1:6 with lines title 'second derivative of f(x)', 'Vals.dat' using 1:7 with lines lw 2 dashtype 2 title 'second derivative of NN(x)'" >> plot.gpi
	@echo "plot 'Vals.dat' using 1:8 with lines title 'antiderivative of f(x)', 'Vals.dat' using 1:9 with lines lw 2 dashtype 2 title 'antiderivative of NN(x)'" >> plot.gpi
	@echo "unset multiplot" >> plot.gpi
	gnuplot plot.gpi
	@cat plot.gpi >> log.gpi
	@rm -f plot.gpi

clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr plot.gpi
