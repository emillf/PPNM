CSS = $(filter %.cs,$^)
	 
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
	
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
	
MKLIBS = mcs -target:library -out:$@ $(CSS)
	
default: Out.txt
	
Out.txt : main.exe
	mono $< >$@
vec.dll : vec.cs
	$(MKLIBS)
main.exe : main.cs vec.dll
	$(MKEXE)
clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
