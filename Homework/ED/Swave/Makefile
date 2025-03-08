libs = ../../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default : Hyd.exe run_tests plot
	
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
Hyd.exe : Hyd.cs matrix.dll
	$(MKEXE)
run_tests: dr_test rmax_test
	
dr_test:
	mono Hyd.exe -mode dr -rmax 10
rmax_test:
	mono Hyd.exe -mode rmax -dr 0.3
plot: dr_plot rmax_plot
	
dr_plot : dr_convergence.dat dr_plot.gnu
	gnuplot dr_plot.gnu
rmax_plot: rmax_convergence.dat rmax_plot.gnu
	gnuplot rmax_plot.gnu
dr_plot.gnu:
	@echo "Creating dr plot script..."
	@echo "set terminal png" > dr_plot.gnu
	@echo "set output 'dr_convergence.png'" >> dr_plot.gnu
	@echo "set title 'Convergence of Ground State Energy with Δr'" >> dr_plot.gnu
	@echo "set xlabel 'Δr'" >> dr_plot.gnu
	@echo "set ylabel 'ε₀'" >> dr_plot.gnu
	@echo "set grid" >> dr_plot.gnu
	@echo "set key top right" >> dr_plot.gnu
	@echo "set xrange [0:*]" >> dr_plot.gnu
	@echo "set yrange [*:0]" >> dr_plot.gnu  # Negative energy values
	@echo "set label 'r_{max} = $(RMAX)' at graph 0.7, 0.8" >> dr_plot.gnu
	@echo "plot 'dr_convergence.dat' using 1:2 with linespoints pt 7 title 'Ground State Energy', \\" >> dr_plot.gnu
	@echo "     -0.5 with lines lt 2 title 'Exact Value'" >> dr_plot.gnu
rmax_plot.gnu:
	@echo "Creating rmax plot script..."
	@echo "set terminal png" > rmax_plot.gnu
	@echo "set output 'rmax_convergence.png'" >> rmax_plot.gnu
	@echo "set title 'Convergence of Ground State Energy with r_{max}'" >> rmax_plot.gnu
	@echo "set xlabel 'r_{max}'" >> rmax_plot.gnu
	@echo "set ylabel 'ε₀'" >> rmax_plot.gnu
	@echo "set grid" >> rmax_plot.gnu
	@echo "set key bottom right" >> rmax_plot.gnu
	@echo "set yrange [*:0]" >> rmax_plot.gnu  # Negative energy values
	@echo "set label 'Δr = $(DR)' at graph 0.7, 0.8" >> rmax_plot.gnu
	@echo "plot 'rmax_convergence.dat' using 1:2 with linespoints pt 7 title 'Ground State Energy', \\" >> rmax_plot.gnu
	@echo "     -0.5 with lines lt 2 title 'Exact Value'" >> rmax_plot.gnu
clean:
	$(RM) *.dat *.png *.gnu *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
