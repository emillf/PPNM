main:plot.ps
	
vals.txt : errfunc.exe
	mono errfunc.exe > vals.txt
errfunc.exe : errfunc.cs
	mcs -target:exe -out:errfunc.exe errfunc.cs
plot.ps : vals.txt
	graph -T ps < vals.txt > plot.ps
clean:                           
	rm -f vals.txt errfunc.exe plot.ps # remove secondary files





