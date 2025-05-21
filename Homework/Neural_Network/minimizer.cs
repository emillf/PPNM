using System;
using static System.Console;
using static System.Math;
using System.IO;
using System.Collections.Generic;
public partial class minimizer{
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
			H[j,j]=dgφ[j]/dxj;
			for(int i=j+1;i<x.size;i++){
				H[i,j]=dgφ[i]/dxj;
				H[j,i]=H[i,j];
				}
			x[j] -= dxj;
			}
		return H;
		}
        public static vector newton(Func<vector, double> φ,vector x,double acc=1e-3, int maxsteps =1000){
		int steps=0;
                do{                   // Newton iterations
                        steps++;
			vector g = gradient(φ,x);
                        if(g.norm() < acc)break;   // job done
                        matrix H = hessian(φ,x);
                        var (Q,R) = matrix.QR.decomp(H);
                        vector dx = matrix.QR.solve(Q,R,-g);
                        double λ = 1.0 ;
			double φx = φ(x);
			double λmin = 1.0/64.0;
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
                return x;
                }
        public static vector newtonAnalytic(Func<vector, double> φ,Func<vector,vector> grad, vector x,double acc=1e-3, int maxsteps=1000){
                int steps=0;
                do{                   // Newton iterations
                        steps++;
                        vector g = grad(x);
			WriteLine($"Step {steps}, ||grad|| = {g.norm()}");
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
                return x;
                }
	}//Class
