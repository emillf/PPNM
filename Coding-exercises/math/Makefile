Out.txt : funcvals.exe              # Out.txt depends on funcvals.exe
	mono funcvals.exe > Out.txt # run funcvals.exe, send output to Out.txt
funcvals.exe : mathexercise.cs             # funcvals.exe depends on mathexercise.cs
	mcs -target:exe -out:funcvals.exe mathexercise.cs # compile mathexercise.cs, save bytecode in funcvals.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt funcvals.exe  # remove secondary files

