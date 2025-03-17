using static System.Console;
using static System.Math;
using System;
using System.IO;
public class qspline{
	vector x,y,b,c;
	public qspline(vector xs, vector ys){
		x = xs.copy();
		y = ys.copy();
		b = new vector(x.size-1);
		c = new vector(x.size-1);
		//Build c with forwards then backwards recursion.
		c[0] = 0;
		for(int i=0;i<c.size-1;i++){
			double dxi = x[i+1]-x[i];
			double dxip1 = x[i+2]-x[i+1];

			double dyi = y[i+1]-y[i];
			double dyip1 = y[i+2]-y[i+1];

			double pi =dyi/dxi;
			double pip1 = dyip1/dxip1;

			c[i+1] = 1/dxip1*(pip1-pi-c[i]*dxi);
			}
		c[c.size-1]/=2.0;
		for(int i=c.size-2;i>1;i--){
			double dxi = x[i+1]-x[i];
                        double dxip1 = x[i+2]-x[i+1];

                        double dyi = y[i+1]-y[i];
                        double dyip1 = y[i+2]-y[i+1];

                        double pi =dyi/dxi;
                        double pip1 = dyip1/dxip1;

			c[i]=1/dxi*(pip1-pi-c[i+1]*dxip1);
	                }
		//Build b from resulting c
		for(int i=0;i<b.size;i++){
			double dxi = x[i+1]-x[i];
			double dyi = y[i+1]-y[i];
			double pi = dyi/dxi;
			b[i] = pi-c[i]*dxi;
			}
		}

		public static int binsearch(vector x, double z){
			if( z<x[0] || z>x[x.size-1] ) throw new Exception("binsearch bad z");
			int i=0, j=x.size-1;
			while(j-i>1){
				int mid=(i+j)/2;
				if(z>x[mid]) i=mid; else j=mid;
				}
		return i;
			}

		public double evaluate(double z){
			int i = binsearch(x,z);
			double siz = y[i]+b[i]*(z-x[i])+c[i]*Pow((z-x[i]),2.0);
			return siz;
			}
		public double derivative(double z){
			int i = binsearch(x,z);
			double dsiz = b[i]+2*c[i]*(z-x[i]);
			return dsiz;
			}
		public double integral(double z, double startval){
			int i = binsearch(x,z);
			double Intval = 0;
			Intval+=startval;
			for(int j=0;j<i;j++){
				Intval+= (y[j]*(x[j+1]-x[j])+b[j]*Pow((x[j+1]-x[j]),2)/2+c[j]*Pow((x[j+1]-x[j]),3)/3);
				}
			double Vallast = (y[i]*(z-x[i])+b[i]*Pow((z-x[i]),2)/2+c[i]*Pow((z-x[i]),3)/3);
			Intval+=Vallast;
			return Intval;
			}
		}
public class Program{
	static int Main(){
		double[] xi = new double[9];
        	double[] yi = new double[9];
        	for(int i=0;i<9;i++){
                	xi[i]=(double)i;
                	yi[i]=Sin((double)i);
                	}
		double npoints = 50;
		var myspline = new qspline(xi,yi);
	using (StreamWriter writer = new StreamWriter("Vals.dat")){
        writer.WriteLine("# xs ys ysinteg");
                for(int i=0;i<xi.Length-1;i++){
			for(int j=0;j<npoints;j++){
                                double zs = xi[i] + (xi[i+1]-xi[i])/npoints*j;
				double xs = xi[i]+(xi[i+1]-xi[i])/npoints*j;
                                double ys = myspline.evaluate(zs);
                                double ysinteg = myspline.integral(zs,-1.0);
				double ysder = myspline.derivative(zs);
                                writer.WriteLine($"{xs} {ys} {ysinteg} {ysder}");
                                }
                        }
			}
		return 0;
		}
	}
