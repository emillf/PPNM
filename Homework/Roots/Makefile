libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt plot
	
plot : epsacc_plot
	
gnus : epsacc_plot.gnu
	
epsacc_plot : Convacceps.dat epsacc_plot.gnu
	gnuplot epsacc_plot.gnu
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs ode.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll
	$(MKEXE)
epsacc_plot.gnu:
	@echo "set terminal svg size 1200,500 background 'white' " > $@
	@echo "set output 'acc_eps_convergence.svg'" >> $@
	@echo "set title 'Convergence of E as a function of acc or eps'" >> $@
	@echo "set xlabel 'acc/eps'" >> $@
	@echo "set ylabel 'E'" >> $@
	@echo "set grid" >> $@
	@echo "set xrange [0:0.07]" >> $@
	@echo "plot 'Convacceps.dat' using 1:2 with lines title 'fixed eps = 0.05', \\" >> $@
	@echo "     'Convacceps.dat' using 3:4 with lines title 'fixed acc = 0.05', \\" >> $@
clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
