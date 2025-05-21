using System;
using static System.Math;
using static System.Console;
using System.Collections.Generic;
using System.IO;

public class Program{
public class ann{
	public int n; /* number of hidden neurons */
	public Func<double,double> f = x => x*Exp(-x*x); /* activation function */
	public Func<double,double> df = x => Exp(-x*x)-2*x*x*Exp(-x*x);
	public Func<double, double> d2f = x => 2.0*Exp(-x*x)*x*(2.0*x*x-3.0);
	public Func<double, double> adf = x => -Exp(-x*x)/2.0;
	public vector p; /* network parameters */
	public ann(int n){
		this.n=n;
		p = new vector(3*n);
		for(int i=0;i<p.size;i+=3){
			p[i]=1.0;
			p[i+1]=1.0;
			p[i+2]=1.0;
			}
		}
	public double response(double x,vector p = null){
	/* return the response of the network to the input signal x */
		double result=0;
		if(p==null) p = this.p;
		for(int i=0;i<p.size;i+=3){
			result+=f((x-p[i])/p[i+1])*p[i+2];
			}
		return result;
		}
	public double response2d(double x, vector p = null){
                double result=0;
                if(p==null) p = this.p;
                for(int i=0;i<p.size;i+=3){
                        result+=d2f((x-p[i])/p[i+1])*(p[i+2]/(p[i+1]*p[i+1]));
                        }
                return result;
                }
        public double responsead(double x, vector p = null){
                double result=0;
                if(p==null) p = this.p;
                for(int i=0;i<p.size;i+=3){
                        result += p[i+1]*p[i+2]*adf((x-p[i])/p[i+1]);
                        }
                return result;
                }

	public Func<vector,vector> GradientCostFunc1(vector x, vector y){
		// Return a function expressing analytic gradient of the first cost function
		Func<vector,vector> Gradient = delegate(vector param){
			vector g = new vector(3*n);
			for(int i=0;i<param.size;i+=3){
				double ai=param[i];
				double bi=param[i+1];
				double wi=param[i+2];
                    		if(Abs(bi) < 1e-10) {
                        		bi = (bi >= 0) ? 1e-10 : -1e-10;
                    			}
				double gradai =0;
				double gradbi =0;
				double gradwi =0;
				for(int k=0; k<x.size;k++){
					double fval=f((x[k]-ai)/bi);
					double dfval= df((x[k]-ai)/bi);
					double rfactor=(1.0/x.size)*2.0*(response(x[k],param)-y[k]);
                            		gradai += (-1.0/bi) * dfval * wi * rfactor;
                            		gradbi += (ai-x[k])/(bi*bi) * dfval * wi * rfactor;
                            		gradwi += fval * rfactor;
                        		}
                    		g[i] = gradai;
                    		g[i+1] = gradbi;
                    		g[i+2] = gradwi;
                		}
			return g;
			};
		return Gradient;
		}
	public void trainnumeric(vector x,vector y){
      /* train the network to interpolate the given table {x,y} */
		Func<vector,double> cost_func = delegate(vector param){
			double result = 0.0;
			for(int i=0;i<x.size;i++){
				result+= Pow(response(x[i],param)-y[i],2.0);
				}
			double output = result/x.size;
			return output;
			};
		p = minimizer.newton(cost_func,p);
		}
	public void trainanalytic1(vector x,vector y){
                Func<vector,double> cost_func = delegate(vector paramt){
                        double result = 0.0;
                        for(int i=0;i<x.size;i++){
                                result+= Pow(response(x[i],paramt)-y[i],2.0);
                                }
                        double output = result/x.size;
                        return output;
                        };
		var Gradient = GradientCostFunc1(x,y);
		p = minimizer.newtonAnalytic(cost_func,Gradient,p);
		}

	}//Ann class
static int Main(){
	int neurons = 3;
	int npoints = 100;
	WriteLine($"Part A)\n\nWe make a neural network of {neurons} neurons and train it on the function Cos(5*x-1)*Exp(-x*x)\n");
	WriteLine($"Training data will be {npoints} random x values on the interval [-1,1] and their corresponding y values from the training function\n");
	Func<double,double> gx = x => Cos(5*x-1)*Exp(-x*x);
	Func<double,double> dgx = x =>  Exp(-x*x) * (-5*Sin(5*x - 1) - 2*x*Cos(5*x - 1));
	Func<double,double> d2gx = x => Exp(-x*x) * ((4*x*x - 27)*Cos(5*x-1)+20*x*Sin(5*x-1));
	vector xs = new vector(npoints);
	vector ys = new vector(npoints);
	vector ysd = new vector(npoints);
	vector ys2d = new vector(npoints);
	vector ysad = new vector(npoints);
	for(int i=0;i<xs.size;i++){
		xs[i]=-1+i*2/npoints;
		ys[i]=gx(xs[i]);
		ysd[i]=dgx(xs[i]);
		ys2d[i]=d2gx(xs[i]);
		if(xs[i]<0.0){
			ysad[i]=-integration.integrate(gx,0.0,xs[i]).Item1;
			}
		else{
			ysad[i]=integration.integrate(gx,0.0,xs[i]).Item1;
			}
		}
	var nn = new ann(neurons);
	var nnum = new ann(neurons);
	nnum.trainnumeric(xs,ys);
	nn.trainanalytic1(xs,ys);
	nn.p.print("The final parameters are: \n ");
	nnum.p.print("The final numerical parameters are: \n ");
	WriteLine($"numerical result at {xs[20]} is {nnum.response(xs[20])} analytic is {nn.response(xs[20])} real result is {gx(xs[20])}");
	return 0;
		}//Main
	}//Program
