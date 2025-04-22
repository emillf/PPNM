libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
plot : epsacc_plot rmaxrmin_plot wavefunc_plot
	
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs ode.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll
	$(MKEXE)
epsacc_plot : Convacceps.dat epsacc_plot.gnu
	gnuplot epsacc_plot.gnu
rmaxrmin_plot : Convrminrmax.dat rmaxrmin_plot.gnu
	gnuplot rmaxrmin_plot.gnu
wavefunc_plot : wavefunc.dat wavefunc_plot.gnu
	gnuplot wavefunc_plot.gnu
wavefunc_plot.gnu:
	@echo "set terminal svg size 1200,500 background 'white' " > $@
	@echo "set output 'Wavefunc.svg'" >> $@
	@echo "set title 'Plot of odesolver wavefunction at E_0 compared to exact wavefunction'" >> $@
	@echo "set xlabel 'r [a_0]'" >> $@
	@echo "set ylabel 'f'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [*:*]" >> $@
	@echo "plot 'wavefunc.dat' using 1:2 with lines title 'exact wavefunction r*exp(-r)', \\" >> $@
	@echo "     'wavefunc.dat' using 3:4 with points pointtype 6 title 'Numerically calculated wavefunction', \\" >> $@
epsacc_plot.gnu:
	@echo "set terminal svg size 1200,500 background 'white' " > $@
	@echo "set output 'acc_eps_convergence.svg'" >> $@
	@echo "set title 'Convergence of E as a function of acc or eps with (rmin,rmax)=(0.3,5.0)'" >> $@
	@echo "set xlabel 'acc/eps'" >> $@
	@echo "set ylabel 'E_0 [Hartree]'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [0:0.07]" >> $@
	@echo "plot 'Convacceps.dat' using 1:2 with points pointtype 2 title 'fixed eps = 0.03', \\" >> $@
	@echo "     'Convacceps.dat' using 3:4 with points pointtype 6 title 'fixed acc = 0.03', \\" >> $@
rmaxrmin_plot.gnu:
	@echo "set terminal svg size 1200,500 background 'white' " > $@
	@echo "set output 'rmax_rmin_convergence.svg'" >> $@
	@echo "set title 'Convergence of E_0 as a function of rmax or rmin with acc = eps = 0.03'" >> $@
	@echo "set xlabel 'rmin [a_0]'" >> $@
	@echo "set ylabel 'E_0 [Hartree]'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [*:*]" >> $@
	@echo "set multiplot layout 1,2 rowsfirst ;\\" >>$@
	@echo "plot \\" >> $@ 
	@echo "'Convrminrmax.dat' using 3:4 with points pointtype 6 title 'fixed rmax = 5.0' ;\\" >> $@
	@echo "set xlabel 'rmax [a_0]' ;\\" >> $@
	@echo "plot \\" >> $@
	@echo "'Convrminrmax.dat' using 1:2 with points pointtype 6 title 'fixed rmin = 0.1' ;\\" >> $@
	@echo "unset multiplot" >> $@
clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
