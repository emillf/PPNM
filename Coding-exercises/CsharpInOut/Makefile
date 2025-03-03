cmdline.exe: cmdline.cs
	mcs -target:exe -out:cmdline.exe cmdline.cs
.PHONY: Out.txt  # Forces this rule to always execute
Out.txt: cmdline.exe
	mono cmdline.exe -numbers:1,2,3,4,5 > Out.txt
clean:  
	rm -f cmdline.exe Out.txt
