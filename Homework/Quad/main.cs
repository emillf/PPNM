using static System.Math;
using System;
using System.IO;
using static System.Console;
class Program{
	static double integrate(Func<double,double> f, double a, double b, double acc=0.001, double eps=0.001,
	double f2=double.NaN, double f3=double.NaN){
		double h = b-a;
		if (double.IsNaN(f2)){
			f2 =f(a+2*h/6);
			f3=f(a+4*h/6);
			}// first call no points to reuse
		double f1=f(a+h/6);
		double f4=f(a+5*h/6);
		double Q =  (2*f1+f2+f3+2*f4)/6*h; //higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
		double err = Abs(Q-q);
		if (err <= acc+eps*Abs(Q)) return Q;
			else return integrate(f,a,(a+b)/2,acc/Sqrt(2),eps,f1,f2)+
            		integrate(f,(a+b)/2,b,acc/Sqrt(2),eps,f3,f4);
		}
	static double erf(double z, double acc = 0.001 , double eps = 0.001){
		Func<double,double> exp = x=> Exp(-x*x);
		Func<double,double> expvar = x => Exp(-Pow((z+(1-x)/x),2))/x/x;
		if(z<0) return -erf(-z);
		else{
			if(Abs(z)<=1) return 2.0/Sqrt(PI)*integrate(exp,0,z,acc,eps);
			if(z>1) return 1.0-2.0/Sqrt(PI)*integrate(expvar,0,1,acc,eps);
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
		WriteLine($"∫ dx √(x) = {integrate(sqrt,0,1)} is within accuracy limits of 2/3? {matrix.approx(0.66666666666666,integrate(sqrt,0,1),0.001,0.001)}\n");
 		WriteLine($"∫ dx 1/√(x) = {integrate(recipsqrt,0,1)} is within accuracy limits of 2? {matrix.approx(2,integrate(recipsqrt,0,1),0.001,0.001)}\n");
		WriteLine($"∫ dx √(1-x²) = {integrate(trigstuff,0,1)} is within accuracy limits of π/2? {matrix.approx(PI/2,integrate(trigstuff,0,1),0.001,0.001)}\n");
		WriteLine($"∫ dx ln(x)/√(x) = {integrate(lnsqrt,0,1)} is within accuracy limits of -4? {matrix.approx(-4,integrate(lnsqrt,0,1),0.001,0.001)}\n");
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
