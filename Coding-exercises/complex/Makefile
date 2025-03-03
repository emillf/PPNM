LIBS := approx.dll
CSS := complexcalc.cs approx.cs

Out.txt :complex.exe
	mono complex.exe > Out.txt
complex.exe : complexcalc.cs $(LIBS)
	mcs -target:exe -out:complex.exe complexcalc.cs $(addprefix -reference:,$(LIBS)) -r:System.Numerics.dll
approx.dll: approx.cs
	mcs -target:library -out:approx.dll approx.cs
clean:                           
	rm -f Out.txt complex.exe  # remove secondary files





