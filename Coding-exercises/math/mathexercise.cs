using System;
class mathvals
{
	public static double sqrt2=Math.Sqrt(2.0);
	public static double fifthroot2=Math.Pow(2.0,0.20);
        public static double etothepowerpi=Math.Pow(Math.E,Math.PI);
        public static double pitothepowere=Math.Pow(Math.PI,Math.E);

	static int Main(){
        	Console.WriteLine($"Sqrt(2)={sqrt2}, 2^(1/5) = {fifthroot2}, e^pi={etothepowerpi}, pi^e={pitothepowere}");
        return 0;

	}
}
