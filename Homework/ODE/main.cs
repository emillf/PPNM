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
	double eps=0.01              /* relative accuracy goal */
	){
	var (a,b)=interval; double x=a; vector y=yinit.copy();
	var xlist=new genlist<double>(); xlist.add(x);
	var ylist=new genlist<vector>(); ylist.add(y);
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
		if(err>0) h *= Min( Pow(tol/err,0.25)*0.95 , 2); // readjust stepsize
		else h*=2;
		}while(true);
	}//driver
	public static Main(){
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
                        	}

