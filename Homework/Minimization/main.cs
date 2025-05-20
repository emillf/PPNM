using System;
using static System.Console;
using static System.Math;
using System.IO;
using System.Collections.Generic;
public class Program{
	static vector gradient(Func<vector,double> φ,vector x){
		var φx = φ(x);
    		vector gφ = new vector(x.size);
    		for(int i = 0;i<x.size;i++){
        		double dxi = (1+Abs(x[i]))*Pow(2,-26);
        		x[i]+=dxi;
			var φxny=φ(x);
        		gφ[i]=(φxny-φx)/dxi;
        		x[i]-=dxi;
			}
    		return gφ;
		}
        static (vector,matrix) gradientHessianCD(Func<vector,double> φ,vector x){
                var φx = φ(x);
                vector gφ = new vector(x.size);
		matrix H = new matrix(x.size,x.size);
		vector dxlist = new vector(x.size);
		vector φxuplist = new vector(x.size);
		vector φxdownlist = new vector(x.size);
		for(int i = 0;i<x.size;i++){
			double dxi = (1+Abs(x[i]))*Pow(2,-26);
			dxlist[i]=dxi;
			vector xtempup = x.copy();
			vector xtempdown = x.copy();
			xtempup[i]=x[i]+dxi;
			xtempdown[i]=x[i]-dxi;
			φxuplist[i]= φ(xtempup);
			φxdownlist[i]= φ(xtempdown);
			gφ[i]=(φxuplist[i]-φxdownlist[i])/(2*dxi);
			}
		for(int i = 0;i<x.size;i++){
			double dxi = dxlist[i];
			H[i,i] = (φxuplist[i]-2*φx+φxdownlist[i])/(dxi*dxi);
			}
		for(int i = 0;i<x.size;i++){
			for(int j = i+1;j<x.size;j++){
				var dxi = dxlist[i];
				var dxj = dxlist[j];
				var xupup=x.copy();
				var xdowndown=x.copy();
				xupup[i]=x[i]+dxi;
				xupup[j]=x[j]+dxj;
                                xdowndown[i]=x[i]-dxi;
                                xdowndown[j]=x[j]-dxj;
				H[i,j]= (φ(xupup)-φxuplist[i]-φxuplist[j]+2*φx-φxdownlist[i]-φxdownlist[j]+φ(xdowndown))/(2*dxi*dxj);
				H[j,i]=H[i,j];
                        	}
			}
                return (gφ,H);
                }
	static matrix hessian(Func<vector,double> φ,vector x){
		matrix H = new matrix(x.size,x.size);
		vector gφx = gradient(φ,x);
		for(int j = 0;j<x.size; j++){
			double dxj = (1+Abs(x[j]))*Pow(2,-26);
			x[j] +=dxj;
			vector gφx_d = gradient(φ, x);
			var dgφ = gφx_d-gφx;
			H[j,j]=dgφ[j]/dxj;
			for(int i=j+1;i<x.size;i++){
				H[i,j]=dgφ[i]/dxj;
				H[j,i]=H[i,j];
				}
			x[j] -= dxj;
			}
		return H;
		}
        static vector newton(Func<vector, double> φ,vector x,double acc=1e-3, bool show_steps=false){
		var stopwatch = System.Diagnostics.Stopwatch.StartNew();
		int steps=0;
		int maxsteps = 1000;
                do{                   // Newton iterations
                        steps++;
			vector g = gradient(φ,x);
                        if(g.norm() < acc)break;   // job done
                        matrix H = hessian(φ,x);
                        var (Q,R) = matrix.QR.decomp(H);
                        vector dx = matrix.QR.solve(Q,R,-g);
                        double λ = 1.0 ;
			double φx = φ(x);
			double λmin = 1.0/1024;
                	do{
				vector xnew = x + λ * dx;
                        	if(φ(xnew) < φx){
					x=xnew;
					break; // good step
					}
				if(λ<λmin){
					x=xnew;
					break; //accept anyway
					}
                        	λ /= 2;
                        	}while(true);
			}while(steps<maxsteps);
		stopwatch.Stop();
		var time = stopwatch.Elapsed.TotalMilliseconds;
		if(show_steps) WriteLine($"Newtons method completed in {steps} steps and  {time} ms \n");
                return x;
                }
        static vector newtonCD(Func<vector, double> φ,vector x,double acc=1e-3, bool show_steps=false){
                var stopwatch = System.Diagnostics.Stopwatch.StartNew();
		int steps=0;
                int maxsteps = 100000;
                do{                   // Newton iterations
                        steps++;
			var (g,H) = gradientHessianCD(φ,x);
                        if(g.norm() < acc)break;   // job done
                        var (Q,R) = matrix.QR.decomp(H);
                        vector dx = matrix.QR.solve(Q,R,-g);
                        double λ = 1.0 ;
                        double φx = φ(x);
                        double λmin = 1.0/1024;
                        do{
                                vector xnew = x + λ * dx;
                                if(φ(xnew) < φx){
					x=xnew;
                                        break; // good step
                                        }
                                if(λ<λmin){
					x=xnew;
                                        break; //accept anyway
                                        }
                                λ /= 2;
                                }while(true);
                        }while(steps<maxsteps);
                stopwatch.Stop();
                var time = stopwatch.Elapsed.TotalMilliseconds;
                if(show_steps) WriteLine($"Newtons method completed in {steps} steps and {time} ms\n");
                return x;
                }
        public static Func<vector,double> Create_Devfunc(List<double> E,List<double> sig,List<double> err){
			Func<vector,double,double> BreitWigner = (y,En) => y[0]/(Pow(En-y[1],2)+y[2]*y[2]/4);
                        return (vector x) => {
				double D =0;
                        	for(int i =0;i<E.Count;i++){
					D += Pow((BreitWigner(x,E[i])-sig[i])/err[i],2);
					}
				return D;
				};
			}
static int Main(string[] filecontent){
	var energy = new List<double>();
	var signal = new List<double>();
	var error  = new List<double>();
	var separators = new char[] {' ','\t'};
	var options = StringSplitOptions.RemoveEmptyEntries;
	do{
        	string line=Console.In.ReadLine();
        	if(line==null)break;
        	string[] words=line.Split(separators,options);
        	energy.Add(double.Parse(words[0]));
        	signal.Add(double.Parse(words[1]));
        	error .Add(double.Parse(words[2]));
		}while(true);
	WriteLine("Part A) \n I have implemented newtons minimization method modified so that it exploits hessian symmetry \n");
	Func<vector, double> Rosenfunc = x => Pow(1-x[0],2)+100*Pow(x[1]-x[0]*x[0],2);
        Func<vector,double> Himmelfunc = x => Pow(x[0]*x[0]+x[1]-11,2)+Pow(x[1]*x[1]+x[0]-7,2);
	WriteLine($"We find a minimum of Rosenbrocks valley function using the initial guess (xi,yi)=(10,10):\n");
	var xstartrosen = new vector(10,10);
        var eksRosen = newton(Rosenfunc,xstartrosen,1e-3,true);
        eksRosen.print("The result is found to be: \n");
        WriteLine($"\n We find a minimum of Himmelblaus function using the initial guess (xi,yi)=(4,4):\n");
        var xstarthimmel = new vector(4,4);
        var ekshimmel = newton(Himmelfunc,xstarthimmel,1e-3,true);
        ekshimmel.print("The result is found to be: \n");
	WriteLine("\n Both of these minima are in accordance with the minima in the wiki \n");
	WriteLine("\n Part B) \n We try to find the mass of the Higgs boson by minimizing the deviation function \n");
	var Deviationfunc = Create_Devfunc(energy,signal,error);
	vector initvals = new vector(8,126,2);
	var Higgslist = newton(Deviationfunc,initvals,1e-3,true);
	var (A,m,Γ) = (Higgslist[0],Higgslist[1],Higgslist[2]);
	Higgslist.print("Our result with initial values (A,m,Γ) is: \n");
	WriteLine($"So this gives a result of {Higgslist[1]:F2}GeV which is the expected value");
	Func<double,double> BWres = x => A/(Pow(x-m,2)+Γ*Γ/4.0);
	using (StreamWriter writer = new StreamWriter("Results.dat")){
        	writer.WriteLine("# Energy[GeV] signal[Certain units]");
		double Einit = energy[0];
		double Efinal = energy[energy.Count-1];
		double Einterval = Efinal-Einit;
		double Nmax = 200;
		double Eslice = Einterval/Nmax;
                for(int N=1;N<Nmax;N++){
                	double Et = Einit+Eslice*N;
			double sigt = BWres(Et);
                        writer.WriteLine($"{Et} {sigt}");
                        }
                }
	WriteLine($"\n \n Part C) \n Now lets test our central difference method, where I have also exploited hessian symmetry \n");
        WriteLine($"Using the same initial values as for the forward difference approach\n");
	WriteLine("We get the following minima for Rosenbrocks and Himmelblaus functions respectively\n Rosen time and steps:");
	var eksRosenCD = newtonCD(Rosenfunc,xstartrosen,1e-3,true);
	eksRosenCD.print("Rosenbrocks minimum: \n");
	WriteLine("Himmelblau time and steps:");
	var eksHimmelCD = newtonCD(Himmelfunc,xstarthimmel,1e-3,true);
	eksHimmelCD.print("Himmelblaus minimum: \n");
	WriteLine("So with these initial guesses the method is faster with the same steps.\nAditionally, if we take (xi,yi)=(0,1)");
	WriteLine("which is a much more difficult initial value. \n The forward difference method gives us: \n");
	vector xstartbad= new vector(0,1);
	var eksHimmel1 = newton(Himmelfunc,xstartbad,1e-3,true);
	eksHimmel1.print("Himmelblaus minimum\n");
	WriteLine("\nwhilst the central difference method gives us \n");
	var eksHimmel1CD = newtonCD(Himmelfunc,xstartbad,1e-3,true);
	eksHimmel1CD.print("Himmelblaus minimum\n");
	WriteLine("\n so it depends on the specific use case \n");
	return 0;
	}//Main
	}//Program
