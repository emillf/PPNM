using System;
using static System.Console;
using MathLibrary;
using System.Numerics;
public class Program {
	public static void Main(string[] args){
		Complex csqrti = Complex.Sqrt(Complex.ImaginaryOne);
        	Complex csqrtneg1 = Complex.Sqrt(new Complex(-1,0));
        	Complex cexpi = Complex.Exp(Complex.ImaginaryOne);
        	Complex cexpipi = Complex.Exp(Complex.ImaginaryOne*Math.PI);
        	Complex ipowi = Complex.Pow(Complex.ImaginaryOne,Complex.ImaginaryOne);
        	Complex lni = Complex.Log(Complex.ImaginaryOne);
        	Complex sinipi = Complex.Sin(Complex.ImaginaryOne*Math.PI);
		bool valcheck = MathUtils.Approx(Math.Pow(Math.E,- Math.PI / 2),ipowi.Real);
		WriteLine("Square root of i: " + csqrti );
        	WriteLine("Square root of -1: " + csqrtneg1 );
        	WriteLine("Exponential of i: " + cexpi );
        	WriteLine("Exponential of i*pi: " + cexpipi );
        	WriteLine("i raised to the power of i: " + ipowi );
        	WriteLine("Natural log of i: " + lni );
        	WriteLine("Sin of i*pi: " + sinipi );
		WriteLine("i^i manually calculated = i^i in complex? " +valcheck);
	}
}
