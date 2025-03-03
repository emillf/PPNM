using System;
using static System.Console;
class program {
	static public void Main(string[] args)
	{
		char[] split_delimiters = {' ','\t','\n'};
		var split_options = StringSplitOptions.RemoveEmptyEntries;
		for( string line = ReadLine(); line != null; line = ReadLine() )
		{
			var numbers = line.Split(split_delimiters,split_options);
			foreach(var number in numbers)
			{
				double x = double.Parse(number);
				Error.WriteLine($"{x} {Math.Sin(x)} {Math.Cos(x)}");
                	}
        	}
	}
}
