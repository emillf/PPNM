MKEXE = mcs -target:exe -out:$@ GenlistGenerator.cs
	
default: Out.txt
	
main.exe : GenlistGenerator.cs
	$(MKEXE)
	
Out.txt : main.exe input.txt
	cat 'input.txt'| ./main.exe > Out.txt
	
clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
