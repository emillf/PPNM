libs = ../../matlib
	
CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)
Hyd.exe : Hyd.cs matrix.dll
	$(MKEXE)
clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
