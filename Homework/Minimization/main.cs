using System;
using static System.Console;
using static System.Math;
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
	static matrix hessian(Func<vector,double> φ,vector x){
		matrix H = new matrix(x.size,x.size);
		vector gφx = gradient(φ,x);
		for(int j = 0;j<x.size; j++){
			double dxj = (1+Abs(x[j]))*Pow(2,-26);
			x[j] +=dxj;
			vector gφx_d = gradient(φ, x);
			var dgφ = gφx_d-gφx;
			for(int i=0;i<x.size;i++){
				H[i,j]=dgφ[i]/dxj;
				}
			x[j] -= dxj;
			}
		return H;
		}
        static vector newton(Func<vector, double> φ,vector x,double acc=1e-3, bool show_steps=false){
		int steps=0;
                while(true){                   // Newton iterations
                        vector g = gradient(φ,x);
                        if(g.norm() < acc)break;   // job done
                        matrix H = hessian(φ,x);
                        var (Q,R) = matrix.QR.decomp(H);
                        vector dx = matrix.QR.solve(Q,R,-g);
                        double λ = 1.0 ;
			double φx = φ(x);
			bool step_accepted=false;
                	while(λ >= 1.0/1024){
				vector xnew = x + λ * dx;
                        	if(φ(xnew) < φx){
					x=xnew;
                                        step_accepted=true;
					break; // good step
					}
                        	λ /= 2;
                        	}
			steps++;
			if(steps>1000){
				WriteLine("Convergence criterion could not be reached");
				break;
				}
            		if(!step_accepted){
               			WriteLine("Line search failed to make progress");
                		break;
            			}
			}
		if(show_steps) WriteLine($"Newtons method completed in {steps} steps \n");
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
	WriteLine("Part A) \n");
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
	WriteLine("\n Part B) \n");
	var Deviationfunc = Create_Devfunc(energy,signal,error);
	vector initvals = new vector(1.0,1.0,1.0);
	var Higgslist = newton(Deviationfunc,initvals,1e-3,true);
	Higgslist.print();
	return 0;
	}//Main
	}//Program
