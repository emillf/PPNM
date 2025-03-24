using static System.Math;
using System.IO;
using System;
public class Program{
	public static (vector,vector) rkstep12(
	Func<double,vector,vector> f,/* the f from dy/dx=f(x,y) */
	double x,                    /* the current value of the variable */
	vector y,                    /* the current value y(x) of the sought function */
	double h                     /* the step to be taken */
	)
	{
	vector k0 = f(x,y);              /* embedded lower order formula (Euler) */
	vector k1 = f(x+h/2,y+k0*(h/2)); /* higher order formula (midpoint) */
	vector yh = y+k1*h;              /* y(x+h) estimate */
	vector δy = (k1-k0)*h;           /* error estimate */
	return (yh,δy);
	}
	public static (genlist<double>,genlist<vector>) driver(
	Func<double,vector,vector> F,/* the f from dy/dx=f(x,y) */
	(double,double) interval,    /* (initial-point,final-point) */
	vector yinit,                /* y(initial-point) */
	double h=0.125,              /* initial step-size */
	double acc=0.01,             /* absolute accuracy goal */
	double eps=0.01, double adj=1.5,              /* relative accuracy goal */
	bool writeh=false
	){
	var (a,b)=interval; double x=a; vector y=yinit.copy();
	var xlist=new genlist<double>(); xlist.add(x);
	var ylist=new genlist<vector>(); ylist.add(y);
	double hmax =(b-a)/500;
	do{
		if(x>=b) return (xlist,ylist); /* job done */
		if(x+h>b) h=b-x;               /* last step should end at b */
		var (yh,δy) = rkstep12(F,x,y,h);
		double tol = (acc+eps*yh.norm()) * Sqrt(h/(b-a));
		double err = δy.norm();
		if(err<=tol){ // accept step
			x+=h; y=yh;
			xlist.add(x);
			ylist.add(y);
			}
		if(err>0){
			double factor = Min( Pow(tol/err,0.25)*0.95 , adj); // readjust stepsize
			h=Min(h*factor,hmax); //enforce minimum step length
			}
		else h=Min(h*2,hmax); //enforce minimum step length
		if(writeh==true) System.Console.WriteLine($"{h}");
		}while(true);
	}//driver
	public static int Main(){
		Func<double, vector, vector> F = delegate(double x, vector y){
			return new vector(y[1],-y[0]);
			};
		vector ystart = new vector(0,1);
		var (xlist,ylist) = driver(F,(0,10.0),ystart);
		using (StreamWriter writer = new StreamWriter("Test.dat")){
        		writer.WriteLine("# xs ys");
                	for(int i=0;i<xlist.size;i++){
                        	double xs = xlist[i];
				double ys = ylist[i][0];
                                writer.WriteLine($"{xs} {ys}");
                                		}
                        	}
		Func<double, vector, vector> Pend = delegate(double x, vector y){
			var theta = y[0];
			var omega = y[1];
			return new vector(omega,-0.25*omega-5.0*Sin(theta));
			};
		vector pendy0 = new vector(PI-0.1,0.0);
		var (pendxlist,pendylist) = driver(Pend,(0,10.0),pendy0);
		using (StreamWriter writer = new StreamWriter("Pend.dat")){
                        writer.WriteLine("# xs ys");
                        for(int i=0;i<pendxlist.size;i++){
                                double xs = pendxlist[i];
                                double ys = pendylist[i][0];
                                writer.WriteLine($"{xs} {ys}");
                                                }
                        }
		double eps1= 0;
		var yinit1=new vector(1,0);
		double eps2 = 0;
		var yinit2 = new vector(1,-0.5);
		double eps3 = 0.01;
		var yinit3 = new vector(1,-0.5);
		Func<double, vector, vector> Orbiteps1 =delegate(double x, vector y){
			return new vector(y[1],1-y[0]+eps1*y[0]*y[0]);
		};
                Func<double, vector, vector> Orbiteps2 =delegate(double x, vector y){
                        return new vector(y[1],1-y[0]+eps2*y[0]*y[0]);
                };
                Func<double, vector, vector> Orbiteps3 =delegate(double x, vector y){
                        return new vector(y[1],1-y[0]+eps3*y[0]*y[0]);
                };
		var (Orbiteps1xlist,Orbiteps1ylist) = driver(Orbiteps1,(0,10*2*PI),yinit1,eps:1e-6,h:10*2*PI/500,acc:1e-6,writeh:true);
                var (Orbiteps2xlist,Orbiteps2ylist) = driver(Orbiteps2,(0,10*2*PI),yinit2);
                var (Orbiteps3xlist,Orbiteps3ylist) = driver(Orbiteps3,(0,10*2*PI),yinit3);
                using (StreamWriter writer = new StreamWriter("Orbit.dat")){
                        writer.WriteLine("# x1 y1 x2 y2 x3 y3");
                        for(int i=0;i<Orbiteps1xlist.size;i++){
                                double x1t = Orbiteps1xlist[i];
                                double y1t = Orbiteps1ylist[i][0];
				//double x1r = 1/y1t*Cos(x1t);
				//double y1r = 1/y1t*Sin(x1t);
				double x1r = y1t != 0 ? Cos(x1t)/y1t : double.NaN;
            			double y1r = y1t != 0 ? Sin(x1t)/y1t : double.NaN;
				writer.WriteLine($"{x1r} {y1r} {double.NaN} {double.NaN} {double.NaN} {double.NaN}");
				}
                        for(int i=0;i<Orbiteps3xlist.size;i++){
                                double x2t = Orbiteps2xlist[i];
                                double y2t = Orbiteps2ylist[i][0];
                                double x2r = 1/y2t*Cos(x2t);
                                double y2r = 1/y2t*Sin(x2t);
                                writer.WriteLine($"{double.NaN} {double.NaN} {x2r} {y2r} {double.NaN} {double.NaN}");
                                }
                        for(int i=0;i<Orbiteps3xlist.size;i++){
                                double x3t = Orbiteps3xlist[i];
                                double y3t = Orbiteps3ylist[i][0];
                                double x3r = 1/y3t*Cos(x3t);
                                double y3r = 1/y3t*Sin(x3t);
                                writer.WriteLine($"{double.NaN} {double.NaN} {double.NaN} {double.NaN} {x3r} {y3r}");
                                }
                        }
		return 0;
		}//Main
	}//Program
