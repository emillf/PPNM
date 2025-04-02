using static System.Math;
using System;
using System.IO;
using static System.Console;
class Program{
	static (double,int,double) integrate(Func<double,double> f, double a, double b, double acc=0.001, double eps=0.001,
	double f2=double.NaN, double f3=double.NaN){
		if(a == double.NegativeInfinity && b == double.PositiveInfinity){
			Func<double,double> newfunc = x => f(x/(1-x*x))*(1+x*x)/Pow(1-x*x,2);
			return CCintegrate(newfunc,-1,1,acc,eps);
			}
		else if(a == double.NegativeInfinity){
			Func<double,double> newfunc = x => f(b-(1-x)/x)/(x*x);
			return CCintegrate(newfunc,0,1,acc,eps);
			}
		else if(b == double.PositiveInfinity){
			Func<double,double> newfunc = x => f(a+(1-x)/x)/(x*x);
			return CCintegrate(newfunc,0,1,acc,eps);
			}
		int ncalls = 0;
		double h = b-a;
		if (double.IsNaN(f2)){
			f2 =f(a+2*h/6);
			f3=f(a+4*h/6);
			ncalls+=1;
			}// first call no points to reuse
		double f1=f(a+h/6);
		double f4=f(a+5*h/6); ncalls+=2;
		double Q =  (2*f1+f2+f3+2*f4)/6*h; //higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
		double err = Abs(Q-q);
		if (err <= acc+eps*Abs(Q)) return (Q,ncalls,err);
		else{
			(double int1,int ncalls1,double err1) = integrate(f,a,(a+b)/2,acc/Sqrt(2),eps,f1,f2);
			(double int2,int ncalls2, double err2) = integrate(f,(a+b)/2,b,acc/Sqrt(2),eps,f3,f4);
			double Int = int1 + int2;
			return (Int,ncalls+ncalls1+ncalls2,Sqrt(err1*err1+err2*err2));
			}
		}
	static (double,int,double) CCintegrate(Func<double,double> f, double a, double b, double acc=0.001, double eps=0.001){
		Func<double,double> newfunc = theta => f((a+b)/2+(b-a)/2*Cos(theta) )*Sin(theta)*(b-a)/2;
		return integrate(newfunc,0,PI,acc,eps);
		}
	static double erf(double z, double acc = 0.001 , double eps = 0.001){
		Func<double,double> exp = x=> Exp(-x*x);
		Func<double,double> expvar = x => Exp(-Pow((z+(1-x)/x),2))/x/x;
		if(z<0) return -erf(-z);
		else{
			if(Abs(z)<=1) return 2.0/Sqrt(PI)*integrate(exp,0,z,acc,eps).Item1;
			if(z>1) return 1.0-2.0/Sqrt(PI)*integrate(expvar,0,1,acc,eps).Item1;
			}
		return double.NaN;
		}
	static int Main(){
		Func<double,double> sqrt = x => Sqrt(x);
		Func<double,double> recipsqrt = x=> 1/Sqrt(x);
		Func<double,double> trigstuff = x=> 1/Sqrt(1-x*x);
		Func<double,double> lnsqrt = x => Log(x)/Sqrt(x);
		WriteLine("Lets check if our recursive open quadrature integrator works \n");
		WriteLine($"We use a default absolute and relative accuracy of 0.001. All integrals are from 0 to 1\n");
		WriteLine($"∫ dx √(x) = {integrate(sqrt,0,1).Item1} is within accuracy limits of 2/3? {matrix.approx(0.66666666666666,integrate(sqrt,0,1).Item1,0.001,0.001)}\n");
 		WriteLine($"∫ dx 1/√(x) = {integrate(recipsqrt,0,1).Item1} is within accuracy limits of 2? {matrix.approx(2,integrate(recipsqrt,0,1).Item1,0.001,0.001)}\n");
		WriteLine($"∫ dx √(1-x²) = {integrate(trigstuff,0,1).Item1} is within accuracy limits of π/2? {matrix.approx(PI/2,integrate(trigstuff,0,1).Item1,0.001,0.001)}\n");
		WriteLine($"∫ dx ln(x)/√(x) = {integrate(lnsqrt,0,1).Item1} is within accuracy limits of -4? {matrix.approx(-4,integrate(lnsqrt,0,1).Item1,0.001,0.001)}\n");
		WriteLine($"\n Part b) \n \n Lets compare our integrator with Clenshaw Curtis variable transformation and our normal integrator\n");
                WriteLine($"∫ dx 1/√(x) = {CCintegrate(recipsqrt,0,1).Item1} with calls {CCintegrate(recipsqrt,0,1).Item2} compared to {integrate(recipsqrt,0,1).Item2} from our normal integrator\n");
                WriteLine($"∫ dx ln(x)/√(x) = {CCintegrate(lnsqrt,0,1).Item1} with calls {CCintegrate(lnsqrt,0,1).Item2} compared to {integrate(lnsqrt,0,1).Item2} from our normal integrator\n");
                WriteLine($"So much fewer evaluations using the Clenshaw Curtis method\n");
                WriteLine($"Scipy needed 213(315) calls for the first(second) integral\n");
		WriteLine($"\n Now lets test our integrator on some converging infinite limit integrals\n");
		Func<double,double> Gaussian = x => Exp(-x*x);
		Func<double,double> recipx2 = x => 1/(x*x);
		Func<double,double> negexp = x => Exp(-x);
		var Gaussint = integrate(Gaussian,double.NegativeInfinity,double.PositiveInfinity);
		var recipx2int = integrate(recipx2,double.NegativeInfinity, -1.0);
		var negexpint = integrate(negexp, 0, double.PositiveInfinity);
		WriteLine($"(from -∞ to ∞) ∫ dx Exp(-x²) = {Gaussint.Item1} which is approximately √π={Sqrt(PI)} ({Gaussint.Item2} calls)\n");
                WriteLine($"(from -∞ to 1) ∫ dx 1/x² = {recipx2int.Item1} which is approximately 1 ({recipx2int.Item2} calls)\n");
                WriteLine($"(from 0 to ∞) ∫ dx Exp(-x) = {negexpint.Item1} which is approximately 1 ({negexpint.Item2} calls)\n");
		WriteLine($"\n Part c) \n");
		WriteLine($"We check some difficult integrals and their error");
                WriteLine($"First lets check the (in)famous Cleo integral see ( https://en.wikipedia.org/wiki/Cleo_(mathematician) )\n");
		Func<double,double> cleofunc = x => 1/x*Sqrt((1+x)/(1-x))*Log((2*x*x+2*x+1)/(2*x*x-2*x+1));
		Func<double,double> upperfunc = x => 1/(4+x*x);
		Func<double,double> doublelim = x => Exp(x)/(1+Exp(2*x));
		var cleoint = integrate(cleofunc,-1,1);
		var doublelimint = integrate(doublelim,double.NegativeInfinity,double.PositiveInfinity);
		var upperfuncint = integrate(upperfunc,0,double.PositiveInfinity);
                WriteLine($"(from -1 to 1) ∫ dx 1/x*Sqrt((1+x)/(1-x))*Log((2*x*x+2*x+1)/(2*x*x-2*x+1))={cleoint.Item1}+-{cleoint.Item3}. \n \n The true error is {Abs(cleoint.Item1-4*PI*Atan(1/Sqrt((1+Sqrt(5))/2)))} compared to 4π*arccot(√ϕ)\n");
                WriteLine($"(from 0 to ∞) ∫ dx 1/(4+x²) = {upperfuncint.Item1}+-{upperfuncint.Item3} The true error is {Abs(upperfuncint.Item1-PI/4)} compared to π/4 \n");
                WriteLine($"(from -∞ to ∞) ∫ dx exp(x)/(exp(2x)+1) = {doublelimint.Item1}+-{doublelimint.Item3} The true error is {Abs(doublelimint.Item1-PI/2)} compared to π/2");
		WriteLine($"\n So the error estimation is fairly high quality it seems");
		using (StreamWriter writer = new StreamWriter("Vals1.dat")){
                writer.WriteLine("# x1 y1");
		double xinterval = 3.5-(-3.5);
		int npoints = 100;
               	for(int i=0;i<npoints;i++){
			double x1 = -3.5+xinterval/100*i;
			double y1 = erf(x1);
                        writer.WriteLine($"{x1} {y1}");
                        }
		}
                using (StreamWriter writer = new StreamWriter("Vals2.dat")){
                writer.WriteLine("# x1 y1");
                int npoints = 10;
		double startacc = 0.1;
                for(int i=1;i<npoints;i++){
			double x1 = Pow(startacc,i);
                        double apxval = erf(1,x1,0);
			double y1 = Abs(0.84270079294971486934-apxval);
                        writer.WriteLine($"{x1} {y1}");
                        }
                }
		return 0;
	}//main
}//Program
