libs = ../matlib
CSS = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CSS)
MKLIBS = mcs -target:library -out:$@ $(CSS)

default: Out.txt

Out.txt : main.exe
	mono main.exe < higgs.data.txt 1> Out.txt 2> log

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs
	$(MKLIBS)

main.exe : main.cs matrix.dll
	$(MKEXE)

clean:
	$(RM) *.exe *.dat *.dll *.png *.svg *.gnu *.log [Oo]ut* [Ll]og* [Ee]rr plot.gpi
