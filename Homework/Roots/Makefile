libs = ../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: main.exe Out.txt #plot
	
#plot : erf_plot acc_plot
	
Out.txt : main.exe
	mono $< >$@
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
ode.dll : ode.cs
	$(MKLIBS)
genlist.dll : $(libs)/genlist.cs
	$(MKLIBS)
main.exe : main.cs matrix.dll genlist.dll ode.dll
	$(MKEXE)
clean:
	$(RM) *.exe *.dll *.png *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr
