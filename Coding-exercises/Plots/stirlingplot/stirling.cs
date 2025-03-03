using static System.Math;
using System;
class Program{
	public static double sgamma(double x){
		if(x<0)return PI/Sin(PI*x)/sgamma(1-x);
		if(x<9)return sgamma(x+1)/x;
		double lnsgamma=Log(2*PI)/2+(x-0.5)*Log(x)-x
    		+(1.0/12)/x-(1.0/360)/(x*x*x)+(1.0/1260)/(x*x*x*x*x);
		return Exp(lnsgamma);
		}
	static void Main(){
	double i=1;
	double n=100;
	double imin=1;
	double imax=4;
	while(i<imax)
		{
		Console.WriteLine($"{i}");
		Console.WriteLine($"{sgamma(i)}");
		i=i+(imax-imin)/n;
		}
	}
}
