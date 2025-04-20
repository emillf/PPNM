using System;
using System.IO;
using static System.Math;
using static System.Console;
public class Program{
	static int Main(){
		public static matrix jacobian(Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
			if(dx == null) dx = x.map(xi => Abs(xi)*Pow(2,-26));
			if(fx == null) fx = f(x);
			matrix J=new matrix(x.size);
			for(int j=0;j < x.size;j++){
				x[j]+=dx[j];
				vector df=f(x)-fx;
				for(int i=0;i < x.size;i++) J[i,j]=df[i]/dx[j];
				x[j]-=dx[j];
				}
			return J;
			}
		static vector newton(
		Func<vector,vector>f /* the function to find the root of */
		,vector start        /* the start point */
		,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
		,vector δx=null      /* optional δx-vector for calculation of jacobian */
		){
			vector x=start.copy();
			vector fx=f(x),z,fz;
			do{ /* Newton's iterations */
				if(fx.norm() < acc) break; /* job done */
				matrix J=jacobian(f,x,fx,δx);
				var QRJ = matrix.QR.decomp(J);
				vector Dx = matrix.QR.solve(QRJ.Item1,QRJ.Item2,-fx); /* Newton's step */
				double λ=1;
				do{ /* linesearch */
					z=x+λ*Dx;
					fz=f(z);
					if( fz.norm() < (1-λ/2)*fx.norm() ) break;
					if( λ < λmin ) break;
					λ/=2;
					}while(true);
				x=z; fx=fz;
				}while(true);
			return x;
			}
		Func<vector,vector> Rosengrad = vector x => new vector y(-2*(1-x[0])-100*2*(x[1]-x[0]*x[0])*2*x[0],2*100*(x[1]-x[0]*x[0]));
		Func<vector, vector> Himmelgrad = vector x => new vector y(2*(x[0]*x[0]+x[1]-11)*2*x[0]+2*(x[0]+x[1]-7),2*(x[0]*x[0]+x[1]-11)+2*(x[0]+x[1]*x[1]-7)*2*x[1]);
		var eksHimmel = newton(Himmelgrad);
		var eksRosen = newton(Rosengrad);
		eksHimmel.Print();
		eksRosen.Print();
		WriteLine($"{}\n");
		}//Main
	}//Program
