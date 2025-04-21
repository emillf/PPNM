using System;
using System.IO;
using static System.Math;
using static System.Console;
using System.Collections.Generic;
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
				double λ=1.0;
				double λmin =Pow(10,-2);
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
        public static Func<vector,vector> Create_M_E(double acc=0.01,double eps=0.01,double rmax=8.0,double rmin=0.05){
                        return (vector E) => {
                                var FEdifinit = new vector(rmin-rmin*rmin,1-2*rmin);
                                Func<double,vector,vector> FEdif = delegate(double x, vector y){
                                        // We rewrite the second order differential equation to 2 first order ones by introducing u0=f' a>                                        // This gives the output equation in the form we want
                                        var u0=y[0];
                                        var u1=y[1];
                                        return new vector(u1,-2*(E[0]*u0+1/x*u0));
                                        };
                                var F_E = ODE.driver(FEdif,(rmin,rmax),FEdifinit,acc,eps);
                                vector M_Eres= new vector(F_E.Item2[F_E.Item2.Count-1][0]);
                               return M_Eres;
                               };
                       }
	static int Main(){
		WriteLine($"Part A)\n");
		WriteLine("The analytically calculated gradients of Rosenbrocks and Himmelblaus functions are found to be:\n");
		WriteLine("∇f(x,y) = (2*(-1+x+200*x^3 - 200*x*y),200*(-x^2+y)) and ∇f(x,y) = (2*(2*x*(x^2+y-11)+x+y^2-7),2*(x^2+2*y*(x+y^2-7)+y-11)) respectively.\n");
                WriteLine($"We find the extremum of Rosenbrocks valley function using the initial guess (xi,yi)=(10,10):\n");
		List<vector> xstartHimmellist = new List<vector>();
		List<vector> resHimmellist = new List<vector>();
                Func<vector,vector> test = x => new vector(3*x[0]*x[0],3*x[1]*x[1]);
                Func<vector,vector> Rosengrad = x => new vector(-2*(1-x[0])-100*2*(x[1]-x[0]*x[0])*2*x[0],2*100*(x[1]-x[0]*x[0]));
                Func<vector, vector> Himmelgrad = x => new vector(2*(x[0]*x[0]+x[1]-11)*2*x[0]+2*(x[0]+x[1]*x[1]-7),2*(x[0]*x[0]+x[1]-11)+2*(x[0]+x[1]*x[1]-7)*2*x[1]);
		for(int i=0;i<2;i++){
			vector xstart1 = new vector(Pow(-1,i)*4,-4);
			vector xstart2 = new vector(Pow(-1,i)*4,4);
			xstartHimmellist.Add(xstart1);
			xstartHimmellist.Add(xstart2);
			resHimmellist.Add(newton(Himmelgrad,xstart1));
			resHimmellist.Add(newton(Himmelgrad,xstart2));
			}
		var xstartrosen = new vector(10,10);
		var eksRosen = newton(Rosengrad,xstartrosen);
		eksRosen.print("The result is found to be: \n");
		WriteLine("In line with the expected result. \nWe find the minima of Himmelblaus function using the following initial guesses:");
		for(int i=0;i<4;i++){
			xstartHimmellist[i].print("\ninitial guess:");
			resHimmellist[i].print("gives result:");
			}
		WriteLine($"\nThis is also in correspondance with expected values\n");
		WriteLine("Part B)\n");
		Func<vector,vector> M_E = Create_M_E();
		vector Einit = new vector(-10.0);
		newton(M_E,Einit).print("With initial guess E=-10 and rmin=0.05 rmax=8.0 we get the result E_0   =");
		WriteLine("\nagain in correspondence with the expected result \n");
		WriteLine("The plots contain an investigation of convergence of these results\n");
		using (StreamWriter writer = new StreamWriter("Convacceps.dat")){
			writer.WriteLine("# acct Eacct epst Eepst");
			double accmax = 1.0;
			double epsmax = 1.0;
			double epsfix = 0.05;
			double accfix = 0.05;
			int Nmax = 100;
			vector Einitnew = new vector(-1.0);
			for(int N=1;N<Nmax;N++){
				double acct=accmax/Nmax*N;
				double epst=epsmax/Nmax*N;
				var M_Eepst = Create_M_E(accfix,epst);
				var M_Eacct = Create_M_E(acct,epsfix);
				var Eepst = newton(M_Eepst,Einitnew)[0];
				var Eacct = newton(M_Eacct,Einitnew)[0];
				writer.WriteLine($"{acct}{Eacct}{epst}{Eepst}");
				}
			}
		return 1;
		}//Main
	}//Program
