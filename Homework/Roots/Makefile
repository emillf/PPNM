libs = ../matlib
CSS = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
MKLIBS = mcs -target:library -out:$@ $(CSS)
DATS = wavefunc.dat Convacceps.dat Convrminrmax.dat

default: Out.txt Wavefunc.svg acc_eps_convergence.svg rmax_rmin_convergence.svg

Out.txt wavefunc.dat Convacceps.dat Convrminrmax.dat : main.exe
	mono $< >Out.txt

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs ode.cs
	$(MKLIBS)

main.exe : main.cs matrix.dll
	$(MKEXE)

Wavefunc.svg: wavefunc.dat
	@echo "set terminal svg size 500,500 background 'white'" > plot.gpi
	@echo "set out '$@'" >> plot.gpi
	@echo "set title 'Plot of odesolver wavefunction at E_0 compared to exact wavefunction'" >> plot.gpi
	@echo "set xlabel 'r [a_0]'" >> plot.gpi
	@echo "set ylabel 'f'" >> plot.gpi
	@echo "set grid" >> plot.gpi
	@echo "set xrange [0:8]" >> plot.gpi
	@echo "plot '$<' using 1:2 with lines title 'exact wavefunction r*exp(-r)', \\" >> plot.gpi
	@echo "     '$<' using 3:4 with points pointtype 6 title 'Numerically calculated wavefunction'" >> plot.gpi
	gnuplot plot.gpi
	@cat plot.gpi>>log.gpi
	@rm -f plot.gpi

acc_eps_convergence.svg: Convacceps.dat Makefile
	@echo "set terminal svg size 1200,500 background 'white'" > plot.gpi
	@echo "set out '$@'" >> plot.gpi
	@echo "set title 'Convergence of E as a function of acc or eps with (rmin,rmax)=(0.05,8.0)'" >> plot.gpi
	@echo "set ylabel 'E_0 [Hartree]'" >> plot.gpi
	@echo "set grid" >> plot.gpi
	@echo "set xrange [*:*]" >> plot.gpi
	@echo "set xlabel 'acc'" >> plot.gpi
	@echo "set multiplot layout 1,2 rowsfirst" >> plot.gpi
	@echo "plot '$<' using 1:2 with points pointtype 2 title 'fixed eps = 0.01'" >> plot.gpi
	@echo "set xlabel 'eps'" >> plot.gpi
	@echo "plot '$<' using 3:4 with points pointtype 6 title 'fixed acc = 0.01'" >> plot.gpi
	@echo "unset multiplot" >> plot.gpi
	gnuplot plot.gpi
	@cat plot.gpi>>log.gpi
	@rm -f plot.gpi

rmax_rmin_convergence.svg: Convrminrmax.dat Makefile
	@echo "set terminal svg size 1000,500 background 'white'" > plot.gpi
	@echo "set out '$@'" >> plot.gpi
	@echo "set title 'Convergence of E_0 as a function of rmax or rmin with acc = eps = 0.03'" >> plot.gpi
	@echo "set xlabel 'rmin [a_0]'" >> plot.gpi
	@echo "set ylabel 'E_0 [Hartree]'" >> plot.gpi
	@echo "set grid" >> plot.gpi
	@echo "set xrange [*:*]" >> plot.gpi
	@echo "set multiplot layout 1,2 rowsfirst" >> plot.gpi
	@echo "plot '$<' using 3:4 with points pointtype 6 title 'fixed rmax = 5.0'" >> plot.gpi
	@echo "set xlabel 'rmax [a_0]'" >> plot.gpi
	@echo "plot [][:0] '$<' using 1:2 with points pointtype 6 title 'fixed rmin = 0.1'" >> plot.gpi
	@echo "unset multiplot" >> plot.gpi
	gnuplot plot.gpi
	@cat plot.gpi>>log.gpi
	@rm -f plot.gpi

clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr plot.gpi
