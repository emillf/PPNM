StandardInput.exe : StandardInput.cs
	mono StandardInput.exe # run hello.exe, send output to Out.txt
	mcs -target:exe -out:StandardInput.exe StandardInput.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f StandardInput.exe  # remove secondary files
