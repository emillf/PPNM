N = 1e8
TIME = time --portability --append --output $@
Out.txt : MultiThreading.exe Makefile
	>$@
	$(TIME) mono $< -terms:$N -threads:1 >>$@
	$(TIME) mono $< -terms:$N -threads:2 >>$@
	$(TIME) mono $< -terms:$N -threads:3 >>$@
	$(TIME) mono $< -terms:$N -threads:4 >>$@
MultiThreading.exe : MultiThreading.cs
	mcs -target:exe -out:MultiThreading.exe MultiThreading.cs
clean:
	rm -f out.txt MultiThreading.exe
