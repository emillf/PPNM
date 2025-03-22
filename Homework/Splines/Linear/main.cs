using static System.Console;
using static System.Math;
using System;
using System.IO;
public class Program{
	public static int binsearch(double[] x, double z){/* locates the interval for z by bisection */ 
		if( z<x[0] || z>x[x.Length-1] ) throw new Exception("binsearch bad z");
		int i=0, j=x.Length-1;
		while(j-i>1){
			int mid=(i+j)/2;
			if(z>x[mid]) i=mid; else j=mid;
			}
		return i;
		}
	 public static double linterp(double[] x, double[] y, double z){
       		int i=binsearch(x,z);
       		double dx=x[i+1]-x[i]; if(!(dx>0)) throw new Exception("x's are not arranged in increasing order");
       		double dy=y[i+1]-y[i];
       		return y[i]+dy/dx*(z-x[i]);
       		}
	 public static double linterpInteg(double[] x ,double[] y, double z){
		int xmaxi=binsearch(x,z);
		double Integval = 0;
		for(int i=0;i<xmaxi;i++){
			double dx = x[i+1]-x[i];
			Integval+=(y[i] + y[i + 1]) * dx / 2;
			}
		double dx_last = z-x[xmaxi];
		Integval+= (y[xmaxi] + linterp(x, y, z)) * dx_last / 2;
		return Integval;
		}
	static int Main(){
		double[] xi = new double[9];
        	double[] yi = new double[9];
        	for(int i=0;i<9;i++){
                	xi[i]=(double)i;
                	yi[i]=Cos((double)i);
                	}
	using (StreamWriter writer = new StreamWriter("Vals.dat")){
        writer.WriteLine("# xs ys xsinteg ysinteg");
		writer.WriteLine($"{xi[0]} {yi[0]} {0} {0}"); // We cant integrate these values unless we change xi so we just hardcode them
                for(int i=0;i<xi.Length-1;i++){
                                double zs = xi[i+1];
				double xsinteg = xi[i+1];
                                double ys = linterp(xi,yi,zs);
                                double ysinteg = linterpInteg(xi,yi,zs);
                                writer.WriteLine($"{xsinteg} {ys} {xsinteg} {ysinteg}");
                                }
                        }
		return 0;
		}
	}
