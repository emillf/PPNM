N = 1e8
TIME = time --portability --append --output $@
Out.txt : MultiThreading.exe Makefile
	>$@
	$(TIME) mono $< -terms:$N >>$@
MultiThreading.exe : MultiThreading.cs
	mcs -target:exe -out:MultiThreading.exe MultiThreading.cs
clean:
	rm -f out.txt MultiThreading.exe
