using System;
using MathLibrary;
class minmaxint {
static int Main() {
	double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
	double d2 = 8*0.1;
	bool result = MathUtils.Approx(d1,d2);
	Machep machep = new Machep();
	int maxintreal=int.MaxValue;
	int minintreal=int.MinValue;
	int i=1;
	while(i+1>i) {i++;}
	int j=1;
	while (j-1<j) {j++;}
	Console.WriteLine("my max int = {0}\n",i);
	Console.WriteLine($"actual max int is {maxintreal}");
	Console.WriteLine("my min int = {0}\n",j);
	Console.WriteLine($"actual min int is {minintreal}");
	Console.WriteLine($"Machine epsilon for a float is {machep.epvalfloat}");
	Console.WriteLine($"Machine epsilon for a double is {machep.epvaldouble}");
	Console.WriteLine($"Machine epsilon should be for a float {machep.epvalaprxfloat}");
	Console.WriteLine($"Machine epsilon should be for a double {machep.epvalaprxdouble}");
	Console.WriteLine($"Machine epsilon should be for a double {machep.epvalaprxdouble}");
	Console.WriteLine($"machep.a==machep.b ? {machep.a==machep.b}\n");
	Console.WriteLine($"machep.a>1  ? {machep.a>1}\n");
	Console.WriteLine($"machep.b>1  ? {machep.b>1}\n");
	Console.WriteLine($"If no approximation is used d1==d2 ? {d1==d2}\n");
	Console.WriteLine($"If approximation is used d1==d2 ? {result}\n");
	return 0;
	}
}
