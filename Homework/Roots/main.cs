using System;
using System.IO;
using static System.Math;
using static System.Console;
public class Program{
		static matrix jacobian(Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
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
				double λmin =Pow(10,-5);
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
		static int Main(){
		var xstarthimmel = new vector(10,10);
		Func<vector,vector> test = x => new vector(3*x[0]*x[0],3*x[1]*x[1]);
		Func<vector,vector> Rosengrad = x => new vector(-2*(1-x[0])-100*2*(x[1]-x[0]*x[0])*2*x[0],2*100*(x[1]-x[0]*x[0]));
		Func<vector, vector> Himmelgrad = x => new vector(2*(x[0]*x[0]+x[1]-11)*2*x[0]+2*(x[0]+x[1]-7),2*(x[0]*x[0]+x[1]-11)+2*(x[0]+x[1]*x[1]-7)*2*x[1]);
		var eksHimmel = newton(Himmelgrad,xstarthimmel);
		var testres = newton(test,xstarthimmel);
		var eksRosen = newton(Rosengrad,xstarthimmel);
		eksHimmel.print();
		testres.print();
		eksRosen.print();
		WriteLine($"{1}\n");
		return 1;
		}//Main
	}//Program
