using genlist;
using System;
using static System.Console;
class FileReader{
	static void Main(){
	var list = new genlist<double[]>;
	char[] delimiters = {' ','\t'};
	var options = StringSplitOptions.RemoveEmptyEntries;
	for(string line = ReadLine(); line!=null; line = ReadLine()){
		var words = line.Split(delimiters,options);
		int n = words.Length;
		var numbers = new double[n];
		for(int i=0;i<n;i++){
			numbers[i]=double.Parse(words[i]);
			list.add(numbers);
			foreach(var number in numbers) Write($"{number : 0.00e+00; -0.00e+00}");
			Writeline();
		}
	WriteLine("\nimplicit cast to array...");
	var array=(double[][])list;
	WriteLine("printing array as: "+array);
	foreach(var i in array){
		foreach(var j in i)Write($"{j}\t");
		WriteLine();
	}
}//Main
}//main

}
