Out.txt :maxminint.exe
	mono maxminint.exe > Out.txt
maxminint.exe : maxminint.cs machep.dll MathLibrary.dll
	mcs -target:exe -out:maxminint.exe maxminint.cs -reference:MathLibrary.dll -reference:machep.dll
machep.dll: machep.cs
	mcs -target:library -out:machep.dll machep.cs
MathLibrary.dll: MathLibrary.cs
	mcs -target:library -out:MathLibrary.dll MathLibrary.cs
clean:                           
	rm -f Out.txt maxminint.exe  # remove secondary files





