using static System.Console;
using static System.Math;
using System;
using System.IO;
public class cspline{
	vector x,y,b,B,d,c;
	matrix Bmatrix;
	public cspline(vector xs, vector ys){
		x = xs.copy();
		y = ys.copy();
		b = new vector(x.size);
		B = new vector(x.size);
		c = new vector(x.size-1);
		d = new vector(x.size-1);
		//Build b from Bmatrix and B with gauss elimination and back substitution
		//First build Bmatrix and B
		Bmatrix = new matrix(x.size,x.size);
		Bmatrix[0,0] = 2;
		Bmatrix[x.size-1,x.size-1] = 2;
                Bmatrix[1,0]=1;
                Bmatrix[x.size-2,x.size-1]=1;
		B[0]=3*(y[1]-y[0]);
		B[x.size-1]=3*(y[x.size-1]-y[x.size-2]);
		for(int i=0;i<x.size-3;i++){
			double hi=x[i+1]-x[i];
			double hip1=x[i+2]-x[i+1];
			double pi =y[i+1]-y[i];
			double pip1 =y[i+2]-y[i+1];
			Bmatrix[i+1,i+1]=2*(hi/hip1)+2;
                        Bmatrix[i+2,i+1]=hi/hip1;
                        Bmatrix[i,i+1]=1;
			B[i+1]=3*(pi+pip1*(hi/hip1));
			}
		// Now do Gauss elimination
		for(int i=1;i<x.size-1;i++){
			Bmatrix[i,i]=Bmatrix[i,i]-Bmatrix[i,i-1]/Bmatrix[i-1,i-1];
			B[i]=B[i]-B[i-1]/Bmatrix[i-1,i-1];
			}
		// Now back substitution
		b[x.size-1]=B[x.size-1]/Bmatrix[x.size-1,x.size-1];
		for(int i=x.size-2;i>0;i--){
			b[i]=(B[i]-Bmatrix[i+1,i]*b[i+1])/Bmatrix[i,i];
			}
		// Finally construct c and d
		for(int i=0;i<c.size-1;i++){
			double hi=x[i+1]-x[i];
			double pi=y[i+1]-y[i];
			c[i]=(-2*b[i]-b[i+1]+3*pi)/hi;
			d[i]=(b[i]+b[i+1]-2*pi)/Pow(hi,2.0);
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
		double siz = y[i]+b[i]*(z-x[i])+c[i]*Pow((z-x[i]),2.0)+d[i]*Pow(z-x[i],3.0);
		return siz;
		}
	public double derivative(double z){
		int i = binsearch(x,z);
		double dsiz = b[i]+2*c[i]*(z-x[i])+3*d[i]*Pow(z-x[i],2);
		return dsiz;
		}
	public double integral(double z, double startval){
		int i = binsearch(x,z);
		double Intval = 0;
		Intval+=startval;
		for(int j=0;j<i;j++){
			Intval+= (y[j]*(x[j+1]-x[j])+b[j]*Pow((x[j+1]-x[j]),2)/2+c[j]*Pow((x[j+1]-x[j]),3)/3)+d[j]*Pow((x[j+1]-x[j]),4)/4.0;
			}
		double Vallast = (y[i]*(z-x[i])+b[i]*Pow((z-x[i]),2)/2+c[i]*Pow((z-x[i]),3)/3+d[i]*Pow((z-x[i]),4)/4);
		Intval+=Vallast;
		return Intval;
		}
	}
public class Program{
	static int Main(){
		vector xi = new vector(9);
        	vector yi = new vector(9);
        	for(int i=0;i<9;i++){
                	xi[i]=(double)i;
                	yi[i]=Sin((double)i);//*Exp(-Pow((double)i/3.0,2.0));
                	}
		double npoints = 50;
		var myspline = new cspline(xi,yi);
	using (StreamWriter writer = new StreamWriter("Vals.dat")){
        writer.WriteLine("# xs ys ysinteg");
                for(int i=0;i<xi.size-1;i++){
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
