using static System.Math;
using System;
using System.IO;
using static System.Console;
public partial class integration{
	public static (double,int,double) integrate(Func<double,double> f, double a, double b, double acc=0.001, double eps=0.001,
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
	public static (double,int,double) CCintegrate(Func<double,double> f, double a, double b, double acc=0.001, double eps=0.001){
		Func<double,double> newfunc = theta => f((a+b)/2+(b-a)/2*Cos(theta) )*Sin(theta)*(b-a)/2;
		return integrate(newfunc,0,PI,acc,eps);
		}
	}//Program
