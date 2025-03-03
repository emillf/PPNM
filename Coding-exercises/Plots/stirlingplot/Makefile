main: plot.ps
	
vals.txt : stirling.exe
	mono stirling.exe > vals.txt
stirling.exe : stirling.cs
	mcs -target:exe -out:stirling.exe stirling.cs
plot.ps : vals.txt
	graph -T ps -X "x" -Y "gamma(x)" < vals.txt > plot.ps
clean:
	rm -f vals.txt stirling.exe plot.ps
