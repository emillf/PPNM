Outln.txt: mainln.exe
	mono mainln.exe > Outln.txt
sfunsln.dll : sfunsln.cs
	mcs -target:library -out:sfunsln.dll sfunsln.cs
mainln.exe : mainln.cs sfunsln.dll
	mcs -target:exe -out:mainln.exe -reference:sfunsln.dll mainln.cs
