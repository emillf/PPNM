using System;
using static System.Math;
public static class jacobi{
	public static void timesJ(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int i=0;i<A.size1;i++){
			double aip = A[i,p];
			double aiq = A[i,q];
			A[i,p]=c*aip-s*aiq;
			A[i,q]=s*aip+c*aiq;
			}
		}
	public static void Jtimes(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int j =0;j<A.size2;j++){
			double apj = A[p,j];
			double aqj = A[q,j];
			A[p,j]=c*apj-s*aqj;
			A[q,j]=s*apj+c*aqj;
			}
		}

	public static (matrix M,vector w, matrix V) cyclic(matrix M){
		matrix A=M.copy();
		matrix V=matrix.id(M.size1);
		vector w=new vector(M.size1);
		bool changed;
		double tolerance = Pow(10,-18);
		int maxIterations = 100;
		int iteration = 0;
		do{
			double[] diagOld = new double[A.size1];
			for (int i = 0; i < A.size1; i++){
				diagOld[i] = A[i, i];
        			}
			changed = false;
			for(int p = 0;p<A.size1-1;p++){
				for(int q=p+1;q<A.size1;q++){
					double apq=A[p,q];
					double aqq=A[q,q];
					double app=A[p,p];
					if(Abs(apq)>tolerance)
						{
						double theta=0.5*Atan2(2*apq,aqq-app);
						timesJ(A,p,q,theta);
						Jtimes(A,p,q,theta);
						timesJ(V,p,q,theta);
						A[p,q]=0;
						A[q,p]=0;
						}
					}
				}
			for (int i = 0; i < A.size1; i++){
				if (Abs(A[i, i]-diagOld[i])>tolerance)
					{
                			changed = true;
                			break;
            				}
        			}
		iteration++;
		}while(changed && iteration<maxIterations);
		for(int i = 0 ;i<A.size1;i++){
			w[i]=A[i,i];
			}
		return (A,w,V);
		}
	}
class Program{
	static void Main(){
		var rnd = new Random(299);
		var A = matrix.symrandom(5,rnd);
		A.print("A random symmetric matrix A");
		Console.WriteLine($"Lets do eigenvalue decomposition");
		var tups = jacobi.cyclic(A);
		var Adiag=tups.Item1;
		Adiag.print("A after using the algorithm:");
		var w=tups.Item2;
		var V=tups.Item3;
		var VT=V.T;
		var VTV=VT*V;
		var VVT=V*VT;
		var VTAV =V.T*A*V;
		var VDVT=V*Adiag*VT;
		VTV.print("V^(T)*V =");
		VVT.print("V*V^(T) =");
		w.print("w = ");
		VDVT.print("VDV^(T)=");
		VTAV.print("V^TAV=");
		bool check1 = VDVT.approx(A);
		bool check2 = VTAV.approx(Adiag);
		Console.WriteLine($"VDV^(T) = A? {check1}");
		Console.WriteLine($"V^(T)AV = D? {check2}");
		}
	}

