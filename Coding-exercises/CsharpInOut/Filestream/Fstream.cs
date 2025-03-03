using System;
using System.IO;
using static System.Console;
class Program
{
	public static int Main(string[] args)
	{
		string infile=null,outfile=null;
		foreach(var arg in args){
			var words=arg.Split(':');
			if(words[0]=="-input")infile=words[1];
			if(words[0]=="-output")outfile=words[1];
			}
		if( infile==null || outfile==null) {
			Error.WriteLine("wrong filename argument");
			return 1;
			}
		var instream =new StreamReader(infile);
		var outstream =new StreamWriter(outfile,append:false);
		for(string line=instream.ReadLine();line!=null;line=instream.ReadLine()){
			double x=double.Parse(line);
			outstream.WriteLine($"{x} {Math.Sin(x)} {Math.Cos(x)}");
        		}
		instream.Close();
		outstream.Close();
		return 0;
	}
}
