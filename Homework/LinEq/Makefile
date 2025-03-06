libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: Out.txt out.n3fit.svg
	
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll
	$(MKEXE)
out.times.data : time.exe
	>$@
	for N in $$(seq 100 20 200); do \
		time --format "$$N %e" --output $@ --append \
	mono $< -size:$$N 1>out 2>err ;\
	done
time.exe : time.cs matrix.dll
	$(MKEXE)
out.n3fit.svg : out.times.data
	echo '\
		set terminal svg background "white";\
		set output "out.n3fit.svg";\
		set xlabel "N";\
		set ylabel "Elapsed time (s)";\
		set grid;\
		set title "Elapsed time as a function of N";\
		f(x) = A*x**3 + B;\
		fit f(x) "out.times.data" via A,B;\
		plot "out.times.data" with points pointtype 6 title "Data", \
			f(x) title "A*N^3 + C";' > log.gpi
	ls -l log.gpi  # Debug: Check if log.gpi exists
	cat log.gpi    # Debug: Show log.gpi content
	gnuplot log.gpi
clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
